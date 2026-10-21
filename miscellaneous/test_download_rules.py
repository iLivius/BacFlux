"""Every rule that downloads must run where Snakemake runs, with pinned download tools.

What this test protects
-----------------------
Two things go wrong with downloads when BacFlux leaves the computer it was written on.

1. Cluster compute nodes often have no internet. On a cluster, Snakemake sends every
   job to a compute node, except the rules marked as local rules (`localrule: True`),
   which it runs itself, on the computer it was started on. So every rule that
   fetches something from the internet must be a local rule.
2. wget, curl, tar, gzip, bzip2 and unzip are not on every Linux computer (a minimal
   Debian or Ubuntu system has no wget, a minimal Red Hat-family system no tar). A
   rule that runs them must get them from a pinned conda environment: its pin list
   (<name>.linux-64.pin.txt, the exact package list Snakemake builds the environment
   from) must hold each program the shell runs. tar also needs gzip and bzip2, which
   it calls to unpack .tar.gz and .tar.bz2 archives.

So this test reads every rule of every .smk file under workflow/rules/, with the rule
reader of test_python_environments.py (the same folder), and requires:
  * for each rule whose shell runs one of those programs: a conda environment whose
    pin list holds the program's package;
  * for each rule that downloads (its shell runs wget or curl, or it is one of the
    rules in OTHER_DOWNLOAD_RULES, which download through a tool of their own): that
    it is a local rule.
A rule that breaks this is reported with its name, its file and the reason.

It also checks rule download_all, which a user names on the command line to fetch
everything a run downloads before the jobs go somewhere without internet. Its input
list is download_targets() in workflow/rules/shared/00_common.smk. The test requires
that the list names an output of every download rule found above (apart from the
Medaka model check, see NOT_IN_DOWNLOAD_ALL), and that download_all is a local rule
whose input is download_targets(). Whether each name sits under the right switch
(PhiX only with Illumina reads, and so on) is checked by dry runs, which this test,
reading the files as text, cannot do.

Run with:  python -m pytest miscellaneous -q
"""

import os
import re

# The repository root, the rule folder, the rule reader and the conda-line reader of the
# Python-environment test, which reads the rule files the same way (as plain text,
# without Snakemake).
from test_python_environments import REPO_ROOT, RULES_DIR, conda_path, read_rules


# The rules that download without running wget or curl themselves, by name:
#   secondary_metabolites_db  download-antismash-databases, a command of antiSMASH
#   virsorter2_db             `virsorter setup`, which runs wget itself
#   conjscan_models           msf_data (MacSyFinder's model installer)
#   check_medaka_model        medaka_model_check.py, which makes Medaka store the
#                             weights of the model it uses
OTHER_DOWNLOAD_RULES = {"secondary_metabolites_db", "virsorter2_db", "conjscan_models",
                        "check_medaka_model"}

# Programs a rule runs that its own shell text does not show. `virsorter setup` (rule
# virsorter2_db) runs wget, gunzip and tar from VirSorter2's packaged setup rules.
HIDDEN_PROGRAMS = {"virsorter2_db": ["wget", "tar", "gzip"]}

# The programs, each used as a command word: not inside a longer name such as libcurl,
# card_tarball or .tar.gz, and not as part of a path (envs/download.yaml).
#   (?<![\w./$-])  the character in front is not a letter, digit, _, ., /, $ or -
#   (?![\w./-])    nor is the character after it
PROGRAMS = ("wget", "curl", "tar", "gzip", "gunzip", "bzip2", "unzip")
PROGRAM_PATTERNS = {
    name: re.compile(r"(?<![\w./$-])" + name + r"(?![\w./-])")
    for name in PROGRAMS
}

# The conda package that brings each program, where the names differ.
PACKAGE_OF = {"gunzip": "gzip"}

# The shell line of a rule, and the local-rule flag.
SHELL_LINE = re.compile(r"^\s*shell:")
LOCALRULE_FLAG = re.compile(r"^\s*localrule:\s*True\s*$")


def all_rule_files():
    """Every .smk file under workflow/rules/, sorted, as full paths."""
    found = []
    for folder, _subfolders, names in os.walk(RULES_DIR):
        for name in names:
            if name.endswith(".smk"):
                found.append(os.path.join(folder, name))
    return sorted(found)


def shell_lines(rule_lines):
    """The lines of one rule from its `shell:` line to its end, comment lines left out."""
    lines = []
    in_shell = False
    for line in rule_lines:
        if SHELL_LINE.match(line):
            in_shell = True
        if in_shell and not line.strip().startswith("#"):
            lines.append(line)
    return lines


def is_local(rule_lines):
    """True when the rule carries `localrule: True`."""
    return any(LOCALRULE_FLAG.match(line) for line in rule_lines)


def packages_needed(rule_name, rule_lines):
    """The conda packages the rule's shell needs, as a sorted list of names."""
    shell_text = "\n".join(shell_lines(rule_lines))
    programs = [name for name, pattern in PROGRAM_PATTERNS.items() if pattern.search(shell_text)]
    programs += HIDDEN_PROGRAMS.get(rule_name, [])
    packages = {PACKAGE_OF.get(program, program) for program in programs}
    # tar unpacks .tar.gz and .tar.bz2 archives by calling gzip and bzip2.
    if "tar" in packages:
        packages.update({"gzip", "bzip2"})
    return sorted(packages)


def downloads(rule_name, rule_lines):
    """True when the rule fetches something from the internet."""
    packages = packages_needed(rule_name, rule_lines)
    return "wget" in packages or "curl" in packages or rule_name in OTHER_DOWNLOAD_RULES


def pin_list_text(smk_path, environment):
    """The text of the pin list next to the rule's environment file, or None if absent.

    The environment path is written relative to the rule file, for example
    ../../envs/download.yaml; its pin list is ../../envs/download.linux-64.pin.txt."""
    yaml_path = os.path.normpath(os.path.join(os.path.dirname(smk_path), environment))
    pin_path = yaml_path[: -len(".yaml")] + ".linux-64.pin.txt"
    if not os.path.isfile(pin_path):
        return None
    with open(pin_path, encoding="utf-8") as handle:
        return handle.read()


def pin_holds(pin_text, package):
    """True when the pin list holds the conda package of that name (for example
    .../linux-64/wget-1.25.0-h653f8fd_1.conda for wget)."""
    return re.search(r"/" + package + r"-\d[^/\s]*\.(conda|tar\.bz2)\s*$", pin_text, flags=re.M) is not None


def problems_in(smk_path):
    """Every rule of one file that breaks the rules above, as readable sentences."""
    problems = []
    relative = os.path.relpath(smk_path, os.path.dirname(RULES_DIR))
    for rule_name, rule_lines in read_rules(smk_path):
        packages = packages_needed(rule_name, rule_lines)
        if packages:
            environment = conda_path(rule_lines)
            if environment is None:
                problems.append(
                    f"rule {rule_name} ({relative}): its shell runs {', '.join(packages)} but it has "
                    "no conda environment, so these come from the computer, which may lack them"
                )
            else:
                pin_text = pin_list_text(smk_path, environment)
                if pin_text is None:
                    problems.append(f"rule {rule_name} ({relative}): {environment} has no pin list next to it")
                else:
                    for package in packages:
                        if not pin_holds(pin_text, package):
                            problems.append(
                                f"rule {rule_name} ({relative}): it needs {package}, but the pin list "
                                f"of {environment} holds no {package} package"
                            )
        if downloads(rule_name, rule_lines) and not is_local(rule_lines):
            problems.append(
                f"rule {rule_name} ({relative}): it downloads from the internet but is not a local rule "
                "(localrule: True), so on a cluster it would run on a compute node, which may have no internet"
            )
    return problems


def local_download_rules():
    """{rule name: its rule lines} for every rule that downloads and is local. A rule
    defined in two files (download_phix) appears once."""
    found = {}
    for smk_path in all_rule_files():
        for rule_name, rule_lines in read_rules(smk_path):
            if downloads(rule_name, rule_lines) and is_local(rule_lines):
                found[rule_name] = rule_lines
    return found


def test_download_rules_are_local_and_pinned():
    problems = []
    for smk_path in all_rule_files():
        problems.extend(problems_in(smk_path))
    assert problems == [], "\n".join(problems)


def test_the_known_download_rules_are_found():
    # A guard on the guard: if the reader stopped finding the rules (a renamed rule, a
    # changed layout), the test above would pass on nothing. These rules download in
    # BacFlux and must be seen as downloading and local.
    expected = {"download_phix", "download_amr_db", "cazyme_db_download",
                "secondary_metabolites_db", "virsorter2_db", "genomad_db", "checkv_db",
                "conjscan_models", "icescan_models", "tncentral_db", "iceberg_db",
                "isosdb_db", "check_medaka_model"}
    found = set(local_download_rules())
    assert expected <= found, sorted(expected - found)


def test_download_phix_is_local_in_both_read_modes():
    # download_phix is written twice, once per short-read mode; both copies must
    # be local.
    seen_in = []
    for smk_path in all_rule_files():
        for rule_name, rule_lines in read_rules(smk_path):
            if rule_name == "download_phix" and is_local(rule_lines):
                seen_in.append(os.path.basename(os.path.dirname(smk_path)))
    assert sorted(seen_in) == ["hybrid", "illumina"], seen_in


# ── rule download_all ──

COMMON_FILE = os.path.join(RULES_DIR, "shared", "00_common.smk")
SNAKEFILE = os.path.join(REPO_ROOT, "workflow", "Snakefile")

# A name of a constant as the rule files write it: capital letters, digits and _, with
# an optional leading _ (PHIX_FASTA, CARD_DB_DIR).
CONSTANT_NAME = re.compile(r"\b_?[A-Z][A-Z0-9_]*\b")

# The output line of a rule; an output may follow on the same line.
OUTPUT_LINE = re.compile(r"^(\s*)output:\s*(.*)$")

# Download rules download_all deliberately leaves out, with the reason. The Medaka
# model check reads the filtered long reads, so it cannot run before the reads are
# filtered; download_targets() says so in its comment.
NOT_IN_DOWNLOAD_ALL = {"check_medaka_model"}


def download_targets_names():
    """The constant names used in the body of download_targets() in 00_common.smk.

    The body runs from the `def download_targets():` line to the next line of code at
    the left margin (or the end of the file); comment lines are left out."""
    with open(COMMON_FILE, encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    names = set()
    inside = False
    for line in lines:
        stripped = line.strip()
        if line.startswith("def download_targets():"):
            inside = True
            continue
        if not inside or stripped == "" or stripped.startswith("#"):
            continue
        # A line of code at the left margin: the function has ended.
        if not line[0].isspace():
            break
        names.update(CONSTANT_NAME.findall(line))
    return names


def output_names(rule_lines):
    """The constant names on the output lines of one rule (from read_rules).

    The outputs run from the `output:` line to the next line indented no deeper than it
    (the next directive, such as params: or conda:)."""
    names = set()
    output_indent = None
    for line in rule_lines:
        stripped = line.strip()
        if output_indent is None:
            match = OUTPUT_LINE.match(line)
            if match:
                output_indent = len(match.group(1))
                names.update(CONSTANT_NAME.findall(match.group(2)))
            continue
        if stripped == "" or stripped.startswith("#"):
            continue
        indent = len(line) - len(line.lstrip())
        if indent <= output_indent:
            break
        names.update(CONSTANT_NAME.findall(line))
    return names


def test_download_all_asks_for_every_download_rule():
    asked_for = download_targets_names()
    assert asked_for, "download_targets() was not found in 00_common.smk, or names nothing"
    problems = []
    # One output of a rule is enough: asking for it makes Snakemake run the rule.
    for rule_name, rule_lines in sorted(local_download_rules().items()):
        if rule_name in NOT_IN_DOWNLOAD_ALL:
            continue
        outputs = output_names(rule_lines)
        if not outputs & asked_for:
            problems.append(
                f"rule {rule_name}: download_targets() in 00_common.smk names none of "
                f"{sorted(outputs)}, so snakemake download_all would not run it"
            )
    assert problems == [], "\n".join(problems)


def test_download_all_is_a_local_rule_asking_for_download_targets():
    rules = read_rules(SNAKEFILE)
    names = [name for name, _lines in rules]
    assert "download_all" in names, "workflow/Snakefile defines no rule download_all"
    # rule all must stay the first rule: Snakemake runs it when no target is named.
    assert names[0] == "all", f"the first rule of the Snakefile is {names[0]}, not all"
    rule_lines = dict(rules)["download_all"]
    assert is_local(rule_lines), "rule download_all is not a local rule"
    assert any(line.strip() == "download_targets()" for line in rule_lines), (
        "rule download_all does not take download_targets() as its input"
    )


def test_the_pattern_tells_commands_from_other_words():
    # Shell words that run a program must be recognised.
    assert PROGRAM_PATTERNS["tar"].search("tar -xzf archive.tar.gz")
    assert PROGRAM_PATTERNS["wget"].search('wget -O "$TAR" "$URL"')
    assert PROGRAM_PATTERNS["gunzip"].search("gunzip -c file.gz")
    # Names and paths that only contain the word must not count.
    assert not PROGRAM_PATTERNS["tar"].search("card_tarball = temp(CARD_TARBALL)")
    assert not PROGRAM_PATTERNS["tar"].search("x.tar.bz2")
    assert not PROGRAM_PATTERNS["tar"].search('echo "$tar"')
    assert not PROGRAM_PATTERNS["curl"].search("libcurl-8.21.0")
