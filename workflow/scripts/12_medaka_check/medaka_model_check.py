#!/usr/bin/env python3
"""Validate the configured Medaka consensus model, or resolve one from the
reads, before the assembler starts — and when the configured model is wrong,
print the closest available ones so the fix is obvious.

Medaka polishing is one of the LAST steps of a long-read run, because it needs
the finished assembly. A model this Medaka version does not have — a typo, or
one dropped in a newer release — therefore used to show up only after Flye had
already spent an hour or more. rule check_medaka_model
(shared/12_medaka_check.smk) runs this script right after read filtering and
gates the assembler on its result, so an unusable model kills the run in
seconds.

Arguments, all handed over by that rule:
  --model    the explicit model name the user set in
             parameters.{nanopore|hybrid}.medaka_model. Empty means auto mode:
             infer the model from the reads.
  --reads    FILT_LONG, the filtlong-filtered ONT reads Medaka will polish
             against. Auto-inference reads the basecaller tag out of their
             FASTQ headers, so the model resolved here is the one the polish
             step would have inferred for itself.
  --fallback-to-auto   the opt-in middle option (config key
             medaka_model_fallback_auto, default off): when an EXPLICIT model
             is invalid, try auto-inference instead of stopping, and use what
             it gives back with a loud warning.
  --out      the validated or resolved model, written as one line. rule
             long_read_consensus (nanopore/30_polish.smk, and its twin in
             hybrid/40_ont_assembly.smk) cats it back and hands it to
             `medaka_consensus -m`, so the model is resolved once, not twice,
             and the two steps cannot disagree.

Exit 0 and write --out on success; print a filtered suggestion table to stderr
and exit 1 on failure.

Auto mode asks for the BACTERIAL consensus model (`--auto_model
consensus_bacteria` below), which corrects the systematic errors bacterial DNA
methylation leaves in canonical basecalls. That is a different job from
Dorado's modified-base calling, and it is still the right model for reads
basecalled without one. It does assume native, unamplified DNA: on an
amplified library pin the matching standard model explicitly instead — see
docs/methods_medaka_model.md.

Downloading the weights here
----------------------------
Medaka ships the weights of only a few models; for any other it downloads them
from GitHub the first time they are needed. Without this step that would
happen in rule long_read_consensus, at polishing time, which on a cluster runs
on a compute node, and compute nodes often have no internet. So the weights are
fetched here instead: in auto mode, working out the model already stores them;
for a named model, download_model() asks Medaka to store them. Rule
check_medaka_model is a local rule, so it runs on the computer that runs
Snakemake, which normally has internet. Medaka keeps the weights inside its
own conda environment (or in ~/.medaka/data when that folder cannot be
written), and the polishing job reads them from there. When the weights
cannot be fetched, the check stops the run with Medaka's own message, before
the assembler starts.

Runs inside the Medaka conda env, so `medaka` is on PATH. Only
list_available_models(), resolve_from_reads() and download_model() call the tool; the parsing,
tokenising and suggesting are kept free of it so they can be unit-tested with
no Medaka installed (test_medaka_model_check.py).
"""

import argparse
import os
import re
import subprocess
import sys


# ── Medaka calls (the only functions that touch the tool) ────────────────────

def list_available_models():
    """Return the model names this Medaka build knows, read from
    `medaka tools list_models` (one 'Available: m1, m2, ...' line, plus Default
    lines that are ignored). Asking the installed tool instead of carrying a
    hard-coded list is the whole point: the valid set is version-specific, which
    is exactly how a config model that worked last year stops working."""
    proc = subprocess.run(
        ["medaka", "tools", "list_models"],
        capture_output=True, text=True,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("Available:"):
            names = line[len("Available:"):].split(",")
            return [n.strip() for n in names if n.strip()]
    # No "Available:" line means the tool itself failed (env broken, etc.).
    raise RuntimeError(
        "Could not read the Medaka model list ('medaka tools list_models' gave no "
        "'Available:' line). stderr:\n" + proc.stderr
    )


def resolve_from_reads(reads):
    """Ask Medaka to infer the consensus model from the basecaller tag in the
    reads, requesting the bacterial variant. Working the model out also makes
    Medaka store its weights.

    Returns a pair (model, problem):
      model    what Medaka prints, a PATH to the model file rather than a bare
               name (a valid `-m` argument all the same), or None when Medaka
               could not give one;
      problem  None on success, otherwise the last lines of Medaka's own message,
               which say why: no basecaller tag in the reads, a tag this Medaka
               does not know, a model the bacterial variant does not cover, or
               weights that could not be downloaded."""
    proc = subprocess.run(
        ["medaka", "tools", "resolve_model", "--auto_model", "consensus_bacteria", reads],
        capture_output=True, text=True,
    )
    model = proc.stdout.strip()
    if proc.returncode == 0 and model:
        return model, None
    return None, last_lines(proc.stderr)


def download_model(model):
    """Make Medaka store the weights of the NAMED model on this computer, so that
    the polishing step finds them without the internet (see "Downloading the
    weights here" at the top of this file). Runs
    `medaka tools download_models --models <model>`, which downloads the weights
    from GitHub only when they are not stored yet. The --models list matters:
    without it Medaka downloads every model it knows.

    Returns None on success, or the last lines of Medaka's message when the
    download failed (for example: no internet, or no folder Medaka may write to)."""
    proc = subprocess.run(
        ["medaka", "tools", "download_models", "--models", model],
        capture_output=True, text=True,
    )
    if proc.returncode == 0:
        return None
    return last_lines(proc.stderr)


def last_lines(text, n=3):
    """The last n non-empty lines of a tool's error output. Medaka ends a failure
    with a Python traceback whose last line says what went wrong."""
    lines = [line for line in text.strip().splitlines() if line.strip()]
    return "\n".join(lines[-n:]) if lines else "(Medaka gave no message)"


# ── Reading a model name: pure logic, unit-tested without Medaka ─────────────

# The axes a user actually scans for when picking a model, in the words the
# Medaka names use. Everything else in a name (pore e82, speed 400bps, the
# neural-net suffixes rl_lstm384_dwells, guppy/version tag) is kept but not used
# for matching.
_DEVICES = {"min": "MinION", "prom": "PromethION"}
_ACCURACIES = {"fast", "hac", "sup", "high"}


def tokenize(model):
    """Split a Medaka model name into the axes a user picks by: flowcell,
    device, accuracy, version (any of them None when the name omits it) plus
    is_variant. Names are underscore-joined and irregular — some carry a device
    token, some do not, versions are either gNNN or vX.Y.Z — so every token is
    classified on its own rather than read from a fixed position."""
    parts = model.split("_")
    out = {"flowcell": None, "device": None, "accuracy": None,
           "version": None, "is_variant": False}
    for p in parts:
        if re.fullmatch(r"r\d+", p):
            out["flowcell"] = p
        elif p in _DEVICES:
            out["device"] = p
        elif p in _ACCURACIES:
            out["accuracy"] = p
        elif re.fullmatch(r"g\d+", p) or re.fullmatch(r"v\d[\w.]*", p):
            # First version-like token wins (g507, v5.2.0); later _rl_lstm... etc.
            # are network variants we keep in the name but do not surface.
            if out["version"] is None:
                out["version"] = p
        elif p in {"variant", "snp"}:
            out["is_variant"] = True
    return out


def consensus_models(available):
    """The subset usable for medaka_consensus polishing: drop the variant/snp
    models (they are for variant calling and would be the wrong choice here)."""
    return [m for m in available if not tokenize(m)["is_variant"]]


def _fmt_table(models):
    """A compact, aligned table with the columns a user picks by. Device tokens
    are spelled out (min -> MinION) since that is what people recognise."""
    rows = []
    for m in sorted(models):
        t = tokenize(m)
        rows.append((
            m,
            t["flowcell"] or "-",
            _DEVICES.get(t["device"], "-"),
            t["accuracy"] or "-",
            t["version"] or "-",
        ))
    headers = ("model", "flowcell", "device", "accuracy", "version")
    widths = [max(len(headers[i]), *(len(r[i]) for r in rows)) for i in range(5)] \
        if rows else [len(h) for h in headers]
    line = lambda cols: "  ".join(c.ljust(widths[i]) for i, c in enumerate(cols))
    return "\n".join([line(headers)] + [line(r) for r in rows])


def suggest(bad_model, available):
    """Build the suggestion text printed when an explicit model is not available.
    Narrows the consensus models to those sharing whatever axes the bad name did
    specify — a typo'd 'r941_min_hac_g508' narrows to the r941 / MinION / hac
    rows, where the real 'r941_min_hac_g507' is the obvious pick — and falls back
    to the whole consensus table when nothing in the name is recognisable."""
    pool = consensus_models(available)
    want = tokenize(bad_model)

    def matches(m):
        t = tokenize(m)
        for axis in ("flowcell", "device", "accuracy"):
            # Only constrain on axes the user's (bad) name actually specified.
            if want[axis] is not None and t[axis] != want[axis]:
                return False
        return True

    narrowed = [m for m in pool if matches(m)] if any(
        want[a] is not None for a in ("flowcell", "device", "accuracy")
    ) else []

    if narrowed:
        return (
            f"Models matching your flowcell/device/accuracy "
            f"(flowcell={want['flowcell']}, device={_DEVICES.get(want['device'], want['device'])}, "
            f"accuracy={want['accuracy']}):\n\n" + _fmt_table(narrowed)
        )
    return (
        "Could not match your model name to a flowcell/device/accuracy. "
        "All consensus models available in this Medaka version:\n\n" + _fmt_table(pool)
    )


def auto_failure_message(available, medaka_said):
    """Guidance when AUTO inference fails. Medaka's own last lines say which of
    the causes it was, so they are shown first; then the user is pointed at
    setting an explicit model, with the menu."""
    return (
        "Medaka could not work out a model from the reads. Medaka said:\n"
        f"{medaka_said}\n\n"
        "The usual causes: the FASTQ headers carry no basecaller tag, or one this "
        "Medaka version does not know; the reads come from a basecaller model the "
        "bacterial model does not cover; or the model's weights could not be "
        "downloaded, because this computer cannot reach GitHub. In the last case, "
        "run Snakemake on a computer that can. Otherwise, set "
        "parameters.<mode>.medaka_model to one of the consensus models below (pick "
        "by your flowcell / device / accuracy), or to FALSE to skip Medaka.\n\n"
        + _fmt_table(consensus_models(available))
    )


# ── Decide the model: auto, explicit, or explicit-then-auto fallback ─────────

# Three ways in, tried in this order: an empty --model means auto, so the reads
# decide; a non-empty one is honoured when it names a model this Medaka build has,
# or a model file already on disk; and only with --fallback-to-auto does an
# invalid explicit name fall back to auto rather than stopping the run. An
# explicit choice is otherwise never swapped for a guess — it is honoured or
# reported, which is why the default is off.
def main(argv=None):
    ap = argparse.ArgumentParser(description="Validate/resolve the Medaka model early.")
    ap.add_argument("--model", default="", help="explicit model, or empty for auto")
    ap.add_argument("--reads", required=True, help="filtlong FASTQ for auto-inference")
    ap.add_argument("--fallback-to-auto", default="false",
                    help="true: on an invalid explicit model, fall back to auto")
    ap.add_argument("--out", required=True, help="write the resolved model name here")
    args = ap.parse_args(argv)

    model = args.model.strip()
    fallback = args.fallback_to_auto.strip().lower() in {"true", "1", "yes", "on"}

    def succeed(resolved, note=None):
        if note:
            sys.stderr.write(note + "\n")
        with open(args.out, "w") as fh:
            fh.write(resolved + "\n")
        print(f"Medaka model resolved to: {resolved}")
        return 0

    def fail(message):
        sys.stderr.write("ERROR: " + message + "\n")
        return 1

    available = list_available_models()

    # ── AUTO mode ────────────────────────────────────────────────────────────
    # resolve_model returns a PATH to the model file (e.g. .../r1041_e82_400bps_
    # bacterial_methylation_model_pt.tar.gz), NOT a bare name, and only on
    # success. So a returned model is already a valid, usable `-m` argument;
    # trust it exactly as v1 did. Do NOT cross-check it against the bare-name
    # list (a path can never match, which would fail every auto run).
    if not model:
        resolved, problem = resolve_from_reads(args.reads)
        if resolved:
            return succeed(resolved)
        return fail(auto_failure_message(available, problem))

    # ── EXPLICIT model ───────────────────────────────────────────────────────
    # A path to a local model file on disk is accepted as-is: Medaka's -m accepts
    # a file, and only NAMED models can be checked against the version's list.
    # (Preserves v1's `[ ! -e "$model" ]` escape hatch.)
    if os.path.exists(model):
        return succeed(model)
    if model in available:
        # The name is right; now make sure its weights are on this computer, so
        # that the polishing step needs no internet.
        download_problem = download_model(model)
        if download_problem is not None:
            return fail(
                f"Medaka model '{model}' exists, but its weights are not stored on "
                f"this computer and could not be downloaded (Medaka downloads them "
                f"from GitHub). Start Snakemake on a computer that can reach the "
                f"internet (on a cluster this check runs where Snakemake runs, not on "
                f"a compute node), or set parameters.<mode>.medaka_model to false to "
                f"skip Medaka. Medaka said:\n{download_problem}"
            )
        return succeed(model)

    # Invalid explicit model. The opt-in middle option: try auto before failing.
    # As in AUTO mode above, resolve_from_reads returns a usable path on success.
    if fallback:
        resolved, problem = resolve_from_reads(args.reads)
        if resolved:
            return succeed(
                resolved,
                note=(
                    f"WARNING: configured medaka_model '{model}' is not available in "
                    f"this Medaka version; fell back to the auto-inferred model "
                    f"'{resolved}' (parameters.<mode>.medaka_model_fallback_auto is "
                    f"on). Set medaka_model to auto to silence this."
                ),
            )
        # Fall-through: fallback was requested but auto could not help either.
        return fail(
            f"Configured medaka_model '{model}' is not available, and fallback "
            f"auto-inference could not infer a usable model from the reads. "
            f"Medaka said:\n{problem}\n\n"
            + suggest(model, available)
        )

    return fail(
        f"Configured medaka_model '{model}' is not available in this Medaka version.\n\n"
        + suggest(model, available)
    )


if __name__ == "__main__":
    sys.exit(main())
