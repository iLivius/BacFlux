# BacFlux v2.0.0 — early Medaka model check: validate the consensus model, or in
# auto mode resolve it, before the assembler starts, so a bad model costs seconds
# instead of hours.
#
# Medaka polishing is one of the LAST steps of a long-read run, because it needs
# the finished assembly. A bad model (a typo, or one dropped in a newer Medaka)
# therefore used to surface only after Flye had already spent an hour or more.
# check_medaka_model runs right after read filtering and the assembler depends on
# its output, so an unusable model kills the run before any assembly happens — see
# the medaka_ok input on rule ont_assembly in nanopore/20_assembly.smk and
# hybrid/40_ont_assembly.smk. It also resolves the model ONCE and writes the name
# down, so long_read_consensus reads it back rather than inferring a second time
# and possibly disagreeing.
#
# Takes in: FILT_LONG — the filtlong-filtered ONT reads Medaka will polish
#           against. Auto-inference reads the basecaller tag out of their FASTQ
#           headers, so pointing this at FILT_LONG rather than the raw input
#           gives exactly the model the polish step would infer for itself.
# Does:     runs scripts/12_medaka_check/medaka_model_check.py — for an explicit model, confirm it
#           appears in `medaka tools list_models` (or is a local model file); for
#           auto, infer it from the reads; on failure, print a filtered suggestion
#           table and exit non-zero so the run stops here, before assembly.
# Produces: MEDAKA_MODEL_RESOLVED — a one-line file holding the confirmed or
#           resolved model NAME. A plain file, deliberately not a directory(),
#           because Snakemake wipes a directory() output before re-running its
#           rule.
# Consumed by: ont_assembly as a gate only — it never reads the contents — and by
#              long_read_consensus, which does read the name back
#              (nanopore/30_polish.smk, hybrid/40_ont_assembly.smk).
#
# All the real logic — validating an explicit model against that Medaka version's
# own model list, auto-inferring from the reads, the opt-in fallback, and the
# flowcell/device/accuracy suggestion table printed on failure — lives in
# scripts/12_medaka_check/medaka_model_check.py, which is unit-tested without Medaka installed.
#
# Guarded on HAS_LONG_READS and USE_MEDAKA, so the rule exists for nanopore and
# hybrid with Medaka on, and never for illumina, contigs, or a Medaka-off run.
# Both flags come from 00_common.smk, as do FILT_LONG, MEDAKA_MODEL,
# MEDAKA_MODEL_RESOLVED, MEDAKA_MODEL_FALLBACK_AUTO, MEDAKA_CHECK_SCRIPT and LOGS.
if HAS_LONG_READS and USE_MEDAKA:

    rule check_medaka_model:
        input:
            reads = FILT_LONG,
        output:
            resolved = MEDAKA_MODEL_RESOLVED,
        params:
            # An empty string means auto mode; a non-empty value is the explicit
            # model the user set in parameters.{nanopore|hybrid}.medaka_model.
            model = MEDAKA_MODEL if MEDAKA_MODEL is not None else "",
            # parameters.{nanopore|hybrid}.medaka_model_fallback_auto, default
            # false: the opt-in middle option where an INVALID explicit model
            # falls back to auto-inference (loudly) instead of stopping the run.
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
