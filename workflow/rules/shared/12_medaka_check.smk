# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — early Medaka model check (rules/shared/12_medaka_check.smk)
#
# Long-read modes only, and only when Medaka is enabled. Medaka polishing is one
# of the LAST steps of a run (it needs the finished assembly), so a bad model — a
# typo, or one removed in a newer Medaka — used to fail only after Flye had
# already spent an hour+ assembling. This rule validates the model right after
# read filtering and gates the assembler on it (see the MEDAKA_MODEL_RESOLVED
# input added to ont_assembly in nanopore/ and hybrid/), turning "fail after
# hours" into "fail in seconds". It also RESOLVES the model once (in auto mode)
# and writes the name for long_read_consensus to reuse, so resolution never runs
# twice.
#
# All the real logic — validating an explicit model against the version's model
# list, auto-inferring from the reads, the opt-in fallback, and the flowcell/
# device/accuracy suggestion table on failure — lives in
# scripts/medaka_model_check.py (unit-tested; see its test file).
#
# Guarded on HAS_LONG_READS and USE_MEDAKA, so the rule exists for nanopore/hybrid
# with Medaka on and never for illumina/contigs or a Medaka-off run. Both those
# names come from 00_common, as do FILT_LONG, MEDAKA_MODEL, MEDAKA_MODEL_RESOLVED,
# MEDAKA_MODEL_FALLBACK_AUTO, MEDAKA_CHECK_SCRIPT and LOGS.
# ─────────────────────────────────────────────────────────────────────────────
if HAS_LONG_READS and USE_MEDAKA:

    # ── Rule: check_medaka_model — validate/resolve the model before assembly ─
    # Takes in: reads = FILT_LONG (the filtered reads Medaka will polish; their
    #           basecaller tag is what auto-inference reads — depending on
    #           FILT_LONG rather than raw reads keeps the resolved model identical
    #           to what the polish step would itself have inferred).
    # Does: run medaka_model_check.py — for an explicit model, confirm it is in
    #       'medaka tools list_models' (or is a local model file); for auto, infer
    #       from the reads; on failure, print a filtered suggestion table and exit
    #       non-zero so the run stops here, before assembly.
    # Produces: MEDAKA_MODEL_RESOLVED — a one-line file holding the confirmed/
    #           resolved model NAME. A plain file (not a directory()), so nothing
    #           is ever wiped; consumed by ont_assembly (as a gate) and by
    #           long_read_consensus (which reads the name back).
    rule check_medaka_model:
        input:
            reads = FILT_LONG,
        output:
            resolved = MEDAKA_MODEL_RESOLVED,
        params:
            # "" here means auto mode; a non-empty value is the explicit model.
            model = MEDAKA_MODEL if MEDAKA_MODEL is not None else "",
            fallback = "true" if MEDAKA_MODEL_FALLBACK_AUTO else "false",
            script = MEDAKA_CHECK_SCRIPT,
        conda:
            "../../envs/medaka.yaml"
        log:
            LOGS + "/check_medaka_model_{sample}.log"
        # Runs before the assembler; give it high priority so the scheduler picks
        # it up promptly and a bad model surfaces as early as possible.
        priority: 10
        shell:
            # On failure, echo the log (which holds the suggestion table — the
            # whole point of this rule) to stderr so it lands on the Snakemake
            # console, not only in the per-sample log file. Mirrors v1's
            # `cat {log} >&2` in every error branch.
            """
            python {params.script} \
              --model "{params.model}" \
              --reads {input.reads} \
              --fallback-to-auto {params.fallback} \
              --out {output.resolved} > {log} 2>&1 || {{ cat {log} >&2; exit 1; }}
            """
