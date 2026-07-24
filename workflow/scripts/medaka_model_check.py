#!/usr/bin/env python3
"""Validate (or auto-infer) the Medaka consensus model BEFORE the expensive
assembly runs — and, when a configured model is wrong, suggest the closest
available ones so the fix is obvious.

WHY THIS EXISTS
    Medaka polishing is one of the LAST steps of a long-read run (it needs the
    finished assembly). So a bad model — a typo, or one dropped in a newer Medaka
    version — used to surface only after Flye had already spent an hour or more.
    rule check_medaka_model runs THIS script right after read filtering and gates
    the assembler on it, so an unusable model kills the run in seconds, not hours.

WHERE IT RUNS
    Inside the Medaka conda env (rule check_medaka_model), so `medaka` is on PATH.
    The pure string logic (parsing/tokenising/suggesting) is deliberately kept
    free of any Medaka call so it can be unit-tested without the tool.

INPUT / OUTPUT
    --model    the configured model. Non-empty = an explicit name the user set;
               empty = auto mode ("infer from the reads").
    --reads    the (filtlong) FASTQ Medaka will polish against; the basecaller tag
               in its headers is what auto-inference reads.
    --fallback-to-auto   the opt-in middle option: when an EXPLICIT model is
               invalid, instead of failing, try auto-inference from the reads and
               use that if it works (with a loud warning). Default off.
    --out      on success, the single validated/resolved model NAME is written
               here (one line). rule long_read_consensus reads it back, so the
               model is resolved once, not twice.

    Exit 0 and write --out on success; print a filtered suggestion table to
    stderr and exit 1 on failure.
"""

import argparse
import os
import re
import subprocess
import sys


# ── Medaka calls (the only functions that touch the tool) ────────────────────

def list_available_models():
    """Return the list of model names Medaka knows about, from
    `medaka tools list_models`. That command prints one line
    'Available: m1, m2, ...' (plus Default lines we ignore). Parsing the live
    tool — not a hard-coded list — is the whole point: the valid set is
    version-specific, which is exactly why a stale config model breaks."""
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
    reads. Returns the resolved model name, or None if Medaka cannot infer one
    (no tag, or a tag its lookup does not recognise). Never raises for the
    ordinary 'could not infer' case — that is a result, not an error."""
    proc = subprocess.run(
        ["medaka", "tools", "resolve_model", "--auto_model", "consensus_bacteria", reads],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        return None
    model = proc.stdout.strip()
    return model or None


# ── Pure logic (unit-tested; no Medaka needed) ───────────────────────────────

# The axes a user actually scans for when picking a model, in the words the
# Medaka names use. Everything else in a name (pore e82, speed 400bps, the
# neural-net suffixes rl_lstm384_dwells, guppy/version tag) is kept but not used
# for matching.
_DEVICES = {"min": "MinION", "prom": "PromethION"}
_ACCURACIES = {"fast", "hac", "sup", "high"}


def tokenize(model):
    """Split a Medaka model name into the axes we care about. Names are
    underscore-joined and irregular (some carry a device token, some do not;
    versions are either gNNN or vX.Y.Z), so classify token-by-token rather than
    assume fixed positions. Returns a dict with flowcell/device/accuracy/version
    (any of which may be None) plus is_variant (snp/variant models are for
    variant CALLING, not consensus polishing)."""
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
    """Build the suggestion text for an invalid explicit model. If the bad name
    has recognisable axes, narrow the (consensus) models to those sharing them —
    e.g. a typo'd version 'r941_min_hac_g508' narrows to the r941/MinION/hac
    rows, where the real 'r941_min_hac_g507' is the obvious pick. If nothing in
    the name is recognisable, fall back to the full consensus table grouped the
    same way."""
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


def auto_failure_message(available):
    """Guidance when AUTO inference fails — the reads carry no basecaller tag, or
    one Medaka's lookup does not know. Point the user at setting an explicit
    model and show the menu."""
    return (
        "Medaka could not infer a model from the reads: their FASTQ headers carry "
        "no basecaller tag, or one this Medaka version does not recognise. Set "
        "parameters.<mode>.medaka_model to one of the consensus models below (pick "
        "by your flowcell / device / accuracy), or to FALSE to skip Medaka.\n\n"
        + _fmt_table(consensus_models(available))
    )


# ── Orchestration ────────────────────────────────────────────────────────────

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
    if not model:
        resolved = resolve_from_reads(args.reads)
        if resolved and resolved in available:
            return succeed(resolved)
        if resolved:
            # Medaka named a model its own list does not contain — treat as a
            # failure rather than trust an unusable name downstream.
            return fail(
                f"Auto-inference returned '{resolved}', which is not in this Medaka "
                f"version's model list.\n\n" + auto_failure_message(available)
            )
        return fail(auto_failure_message(available))

    # ── EXPLICIT model ───────────────────────────────────────────────────────
    # A path to a local model file on disk is accepted as-is: Medaka's -m accepts
    # a file, and only NAMED models can be checked against the version's list.
    # (Preserves v1's `[ ! -e "$model" ]` escape hatch.)
    if os.path.exists(model):
        return succeed(model)
    if model in available:
        return succeed(model)

    # Invalid explicit model. The opt-in middle option: try auto before failing.
    if fallback:
        resolved = resolve_from_reads(args.reads)
        if resolved and resolved in available:
            return succeed(
                resolved,
                note=(
                    f"WARNING: configured medaka_model '{model}' is not available in "
                    f"this Medaka version; fell back to the auto-inferred '{resolved}' "
                    f"(parameters.<mode>.medaka_model_fallback_auto is on). Set "
                    f"medaka_model to '{resolved}' or to auto to silence this."
                ),
            )
        # Fall-through: fallback was requested but auto could not help either.
        return fail(
            f"Configured medaka_model '{model}' is not available, and fallback "
            f"auto-inference could not infer a usable model from the reads.\n\n"
            + suggest(model, available)
        )

    return fail(
        f"Configured medaka_model '{model}' is not available in this Medaka version.\n\n"
        + suggest(model, available)
    )


if __name__ == "__main__":
    sys.exit(main())
