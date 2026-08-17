#!/usr/bin/env python3
"""Check that the things BacFlux's comments name actually exist.

WHY THIS EXISTS
    Comments rot silently. A rule gets renamed, a stage directory is renumbered, a
    constant moves to another file — and every comment pointing at the old name
    keeps sitting there looking authoritative. Nothing fails, no test goes red, and
    the reader who trusts it loses an afternoon.

    A large comment rewrite of this repo (2026-08-14) turned up dozens of exactly
    that: a cited `rule mobilome_name_transposons` that was never built, v1 stage
    paths (`05.annotation` where v2 has `04.annotation`), a count that said four
    where the code had five. Every one of them was mechanically checkable. This
    script does that check, so the next batch gets caught in seconds instead of by
    someone reading the file two years later.

WHAT IT CHECKS
    Four kinds of claim, chosen because each is unambiguous — either the thing is
    there or it is not, with no judgement involved:

      rules    every "rule <name>" mentioned in a comment must be defined somewhere
               under workflow/rules/
      paths    every repo-relative path in a comment (workflow/..., config/...,
               docs/..., miscellaneous/...) must exist on disk
      globals  every ALL_CAPS constant named in a comment must be assigned
               somewhere in the .smk tree
      config   every dotted config key (mobilome.run, parameters.eggnog.dbmem)
               must exist in config/config.yaml

    It also reports, separately and as information rather than an error, any comment
    living INSIDE a shell: block. Those are part of the script Snakemake hashes, so
    rewording one makes Snakemake re-run the rule. Harmless on a database download,
    expensive on GTDB-Tk — worth knowing before you edit one.

WHAT IT DOES NOT CHECK
    Whether an explanation is TRUE. A comment saying "we filter here because X" is
    beyond any script; only a reader who knows the biology can judge it. This tool
    clears away the mechanical errors so that human attention goes where it is the
    only thing that works.

RUN
    python miscellaneous/check_comment_references.py           # from the repo root
    python miscellaneous/check_comment_references.py --quiet   # only the failures

    Exit code 0 = every reference resolves. 1 = at least one does not.
"""

import argparse
import difflib
import os
import re
import sys


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Files whose comments we read. Env YAMLs are conda pins with no prose.
SOURCE_SUFFIXES = (".smk", ".py", ".sh", ".md")
SOURCE_NAMES = ("Snakefile",)
SKIP_DIRS = {".git", ".snakemake", "__pycache__", "ref"}

# ── What a claim looks like ──────────────────────────────────────────────────
# "rule foo_bar", "rules foo_bar and baz_qux". The UNDERSCORE is required, and that
# is a deliberate trade. Without it the pattern matches ordinary English — "the rule
# always fires", "the rule asks for" — and produces so much noise that a real finding
# is invisible, which is how checkers end up ignored. The cost is that a reference to
# a single-word rule (isescan, conjscan, multiqc, annotation) is not checked. Most
# BacFlux rule names are snake_case, and every stale reference the 2026-08-14 comment
# pass turned up was of that shape, so this catches the class that actually occurs.
RULE_REFERENCE = re.compile(r"\brules?\s+([a-z][a-z0-9]*(?:_[a-z0-9]+)+)\b")

# A repo-relative path. Anchored on the four real top-level directories so that
# "AMR/mapping" or "identity/coverage" in prose is not mistaken for a file.
PATH_REFERENCE = re.compile(
    r"\b((?:workflow|config|docs|miscellaneous)/[A-Za-z0-9_./-]*[A-Za-z0-9_])"
)

# An ALL_CAPS global. Three characters minimum, and it must contain an underscore
# or be long, so that shouted prose (NOT, ONLY, THIS) is not treated as a name.
GLOBAL_REFERENCE = re.compile(r"\b([A-Z][A-Z0-9]*(?:_[A-Z0-9]+)+)\b")

# A dotted config key, all lowercase: mobilome.run, parameters.long_read_qc.min_length
# A snake_case identifier written either in backticks or as the first line of a
# Mermaid node label — the two places BacFlux's documentation names a rule.
# How close a name must be to a real rule before it is called a probable typo.
# 0.70 is empirical: 'check_medaka_check' scores 0.72 against 'check_medaka_model',
# the actual mistake this check exists for. Raising it past 0.75 misses that; lowering
# it much below 0.70 starts pairing unrelated names that merely share a prefix.
NEAR_MISS_CUTOFF = 0.70

SNAKE_NAME = re.compile(r"`([a-z][a-z0-9]*(?:_[a-z0-9]+)+)`"
                        r"|\[([a-z][a-z0-9]*(?:_[a-z0-9]+)+)(?:<br/>|\])")

CONFIG_REFERENCE = re.compile(r"\b([a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)+)\b")

# A comment that is deliberately describing what something USED to be called is
# correct precisely because the name no longer resolves. BacFlux does this
# constantly — "(v1 message: ...)", "v1 rule was called map_qc", "these used to be
# flat" — and without this the checker flags the repo's own change log as broken.
# Any comment line carrying one of these markers is skipped entirely.
HISTORICAL_MARKERS = re.compile(
    r"\b(v1|used to|was called|were called|renamed|formerly|old name|previously|"
    r"superseded|no longer|used to be|predates|pre-dates|stale|"
    r"that path|never built|does not exist|was never)\b", re.IGNORECASE
)

# Snakemake's own directive namespaces. In a comment, "input.contigs" is nearly
# always a reference to a rule's own input: block, not to config["input"]["contigs"]
# — but config really does have an `input:` section, so the two collide. Resolved by
# only flagging such a key when the FULL dotted path is a real config key that has
# been mistyped; a first segment that is a Snakemake directive is otherwise ignored.
SNAKEMAKE_DIRECTIVES = {"input", "output", "params", "wildcards", "log",
                        "resources", "threads", "config", "rules"}

# Dotted things that are not config keys: filenames, module paths, tool invocations.
CONFIG_SKIP_SUFFIXES = (
    ".py", ".sh", ".smk", ".tsv", ".txt", ".fasta", ".fa", ".fna", ".faa", ".gff",
    ".gff3", ".gbff", ".json", ".yaml", ".yml", ".log", ".html", ".zip", ".gz",
    ".bz2", ".md", ".fastq", ".fq", ".bam", ".sam", ".ms", ".gbk", ".fas", ".xz",
)


def source_files():
    """Every file in the repo whose comments this tool should read."""
    found = []
    for root, dirs, files in os.walk(REPO_ROOT):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for name in files:
            if name.endswith(SOURCE_SUFFIXES) or name in SOURCE_NAMES:
                found.append(os.path.join(root, name))
    return sorted(found)


def comment_lines(path):
    """Yield (line_number, text, inside_shell_block) for every comment in a file.

    `inside_shell_block` tracks whether the comment sits within a Snakemake
    shell: body, because those count as code (see the module docstring). The
    tracker is deliberately simple — it opens on a `shell:` line and closes on the
    next line that is a triple quote at lower-or-equal indentation. Snakemake rule
    bodies are regular enough for that to hold.
    """
    inside_shell = False
    shell_indent = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for number, raw in enumerate(handle, start=1):
            stripped = raw.strip()
            indent = len(raw) - len(raw.lstrip())

            if not inside_shell and stripped.startswith("shell:"):
                inside_shell = True
                shell_indent = indent
                continue
            if inside_shell and stripped.startswith('"""') and indent <= shell_indent + 4:
                # The opening triple quote also matches; only close on a later one.
                if number > 0 and stripped == '"""':
                    inside_shell = False
                continue

            if path.endswith(".md"):
                # A Markdown file is prose end to end: there is nothing to strip and
                # no shell block to be inside.
                yield number, raw.rstrip("\n"), False
            elif stripped.startswith("#"):
                yield number, stripped.lstrip("#").strip(), inside_shell


def collect_defined_rules():
    """Every rule name defined anywhere under workflow/."""
    names = set()
    pattern = re.compile(r"^\s*rule\s+([a-z_][a-z0-9_]*)\s*:")
    for path in source_files():
        if not (path.endswith(".smk") or os.path.basename(path) in SOURCE_NAMES):
            continue
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                match = pattern.match(line)
                if match:
                    names.add(match.group(1))
    return names


def collect_defined_globals():
    """Every ALL_CAPS name assigned anywhere in the .smk tree or the scripts."""
    names = set()
    pattern = re.compile(r"^\s*([A-Z][A-Z0-9_]*)\s*=")
    for path in source_files():
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                match = pattern.match(line)
                if match:
                    names.add(match.group(1))
    return names


def collect_config_keys():
    """Every dotted key path present in config/config.yaml.

    Parsed with PyYAML when available, and by indentation otherwise, so the check
    still works in a bare interpreter with no third-party packages.
    """
    config_path = os.path.join(REPO_ROOT, "config", "config.yaml")
    if not os.path.exists(config_path):
        return set()
    keys = set()
    try:
        import yaml
        with open(config_path, encoding="utf-8") as handle:
            data = yaml.safe_load(handle)

        def walk(node, prefix):
            if not isinstance(node, dict):
                return
            for key, value in node.items():
                dotted = f"{prefix}.{key}" if prefix else str(key)
                keys.add(dotted)
                walk(value, dotted)

        walk(data, "")
    except Exception:                                           # noqa: BLE001
        stack = []
        with open(config_path, encoding="utf-8") as handle:
            for line in handle:
                match = re.match(r"^(\s*)([a-z][a-z0-9_]*)\s*:", line)
                if not match:
                    continue
                depth = len(match.group(1)) // 2
                stack = stack[:depth] + [match.group(2)]
                keys.add(".".join(stack))
    return keys


def check():
    rules = collect_defined_rules()
    globals_defined = collect_defined_globals()
    config_keys = collect_config_keys()

    # Every snake_case name that legitimately is NOT a rule: config keys, output
    # column names, filenames, constants. Without this the near-miss check would
    # flag things like `amr_gene_family` for merely resembling a rule.
    defined_elsewhere = set()
    for _p in source_files():
        try:
            _t = open(_p, encoding="utf-8", errors="replace").read()
        except Exception:                                   # noqa: BLE001
            continue
        defined_elsewhere.update(re.findall(r"^\s*([a-z][a-z0-9_]+)\s*[:=]", _t, re.M))
    defined_elsewhere.update(k.split(".")[-1] for k in config_keys)

    problems = []
    advisories = []
    shell_comments = []

    for path in source_files():
        relative = os.path.relpath(path, REPO_ROOT)
        if os.path.abspath(path) == os.path.abspath(__file__):
            continue                        # this file's own docstring examples
        in_historical_block = False
        previous_number = -10
        for number, text, in_shell in comment_lines(path):
            if number != previous_number + 1:
                in_historical_block = False     # a gap ends the block
            previous_number = number
            if HISTORICAL_MARKERS.search(text):
                in_historical_block = True
            if in_historical_block:
                continue                    # deliberately describing a former name
            where = f"{relative}:{number}"
            if in_shell:
                shell_comments.append(where)

            for name in RULE_REFERENCE.findall(text):
                if name not in rules and name not in {"all", "above", "below", "that",
                                                      "the", "and", "in", "for", "is",
                                                      "it", "this", "which", "with"}:
                    problems.append(f"{where}: no rule named '{name}'")

            for candidate in PATH_REFERENCE.findall(text):
                if "{" in candidate or "*" in candidate:
                    continue                                    # a template, not a path
                if candidate.endswith("envs/x.yaml"):
                    continue                                    # the stand-in used to
                                                                # explain how the
                                                                # ../../envs/ hop works
                if not os.path.exists(os.path.join(REPO_ROOT, candidate)):
                    problems.append(f"{where}: path does not exist: {candidate}")

            for match in SNAKE_NAME.findall(text):
                # Two alternatives in the pattern, so findall yields a pair; exactly
                # one half is ever populated.
                bare = match[0] or match[1]
                if bare in rules or bare in defined_elsewhere:
                    continue
                near = difflib.get_close_matches(bare, sorted(rules), n=1, cutoff=NEAR_MISS_CUTOFF)
                if near:
                    # ADVISORY, not a failure. Measured on this repo the ratio is about
                    # nine false alarms to one real find: output column names
                    # (plasmid_score, same_replicon) and config keys are legitimately
                    # snake_case and resemble rule names closely. As a gate it would be
                    # switched off within a week; as a list you read after editing the
                    # documentation it is worth having, because it does catch the real
                    # thing — 'check_medaka_check' for 'check_medaka_model', which no
                    # exact-match test can see.
                    advisories.append(
                        f"{where}: '{bare}' is not a rule; did you mean '{near[0]}'?")

            for name in GLOBAL_REFERENCE.findall(text):
                if name in globals_defined:
                    continue
                advisories.append(f"{where}: '{name}' is not defined in this repo")

            for key in CONFIG_REFERENCE.findall(text):
                if key.endswith(CONFIG_SKIP_SUFFIXES) or "/" in key:
                    continue
                if key.split(".")[0] in SNAKEMAKE_DIRECTIVES:
                    continue                                    # a rule's own input:/params:
                if key.split(".")[0] not in {k.split(".")[0] for k in config_keys}:
                    continue                                    # not a config path at all
                if key not in config_keys:
                    problems.append(f"{where}: no config key '{key}'")

    return problems, advisories, shell_comments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true",
                        help="print only the failures, not the summary counts")
    parser.add_argument("--advisories", action="store_true",
                        help="also list undefined ALL_CAPS names (noisy by nature)")
    args = parser.parse_args()

    problems, advisories, shell_comments = check()

    for problem in problems:
        print(problem)
    if args.advisories:
        for advisory in advisories:
            print(advisory)

    if not args.quiet:
        print()
        print(f"{len(advisories)} ALL_CAPS name(s) in comments are not defined in this")
        print("  repo. Most are accessions, locus tags, HMM profile names or another")
        print("  tool's constants, which is fine — listed with --advisories, not failed.")
        print()
        print(f"{len(shell_comments)} comment line(s) sit inside a shell: block.")
        print("  Those are part of the script Snakemake hashes: rewording one makes")
        print("  the rule re-run. Cheap on a download, expensive on GTDB-Tk.")
        print()
        if problems:
            print(f"FAIL: {len(problems)} reference(s) in comments do not resolve.")
        else:
            print("OK: every rule, path, global and config key named in a comment exists.")

    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
