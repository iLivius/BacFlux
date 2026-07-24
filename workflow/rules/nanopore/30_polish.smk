# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — nanopore front end, polishing (rules/nanopore/30_polish.smk)
#
# The biology: ONT reads are excellent at telling you the STRUCTURE of a genome
# and still imperfect at telling you the exact base, especially in homopolymers.
# Medaka re-reads the raw signal-derived reads against the draft assembly with a
# neural network trained on the same basecaller, and rewrites the consensus. This
# is the step that turns a structurally-correct ONT assembly into one whose gene
# calls can be trusted.
#
# Data flow through this module:
#
#   DECONTAM_CONTIGS ──┐
#   (shared/10_decontam)│
#   FILT_LONG ──────────┼──► long_read_consensus (Medaka) ──► MEDAKA_CONSENSUS ─┐
#   raw ONT FASTQ ──────┘        [only when USE_MEDAKA]                          │
#                                                                                ▼
#                                                              finalize_contigs (cp)
#                                                                                │
#                                                                                ▼
#                                                                       FINAL_CONTIGS
#                                                                (the D2 hand-off:
#                                                                 QC, taxonomy,
#                                                                 annotation, AMR,
#                                                                 plasmids, phages)
#
# D3 ORDERING: Medaka polishes the DECONTAMINATED assembly, not the raw one. That
# is v1 nanopore behaviour and it is deliberate — polishing a contaminant contig
# with this isolate's reads would be wasted effort at best and would smear real
# differences at worst.
#
# WHY finalize_contigs EXISTS AT ALL: whether Medaka runs is a config choice, so
# the last file in this mode's chain is either the Medaka consensus or the
# decontaminated assembly. `cp`ing whichever it is to ONE canonical name means
# every downstream module can depend on FINAL_CONTIGS and never ask the question.
#
# v1 -> v2 changes in this file:
#  * v1 nanopore had a PAIR of anonymous rules (`rule:` with no name), selected by
#    an if/else on the config. v2 has ONE named rule under `if USE_MEDAKA:`, with
#    the explicit-model and auto-inference paths as a bash if/else inside it —
#    i.e. the hybrid v1 form, which is easier to read and easier to log about.
#  * The model pre-flight check (validate the model name BEFORE starting) is kept
#    from v1 nanopore; the post-failure hint text is kept from v1 hybrid.
#  * Auto-inference now reads the RAW ONT file rather than the filtlong output,
#    because the basecaller tag lives in the original read headers and filtlong
#    is not guaranteed to keep a read that carries it. This is v1 hybrid's choice,
#    adopted here too.
#  * v1's `final_contigs(wc)` input FUNCTION (which re-read the config at DAG
#    build time) is replaced by the parse-time constant FINALIZE_SOURCE.
#
# Everything referenced here comes from 00_common.smk: USE_MEDAKA, MEDAKA_MODEL,
# MEDAKA_INPUT (= DECONTAM_CONTIGS in this mode), MEDAKA_DIR, MEDAKA_CONSENSUS,
# FINALIZE_SOURCE, FINAL_CONTIGS, FILT_LONG, NANOPORE_DIR, ONT, LOGS, capped_cpus.
# ─────────────────────────────────────────────────────────────────────────────


# The rule below exists ONLY when Medaka is enabled. USE_MEDAKA is resolved once
# at parse time in 00_common section 8: it is False only when the user explicitly
# set parameters.nanopore.medaka_model to a false-like value. A missing or "auto"
# value means "run Medaka and work the model out yourself", NOT "skip".
if USE_MEDAKA:

    # ── Rule: long_read_consensus — ONT consensus polishing (Medaka) ─────────
    # Takes in:
    #   reads    = FILT_LONG        the same reads Flye assembled
    #   contigs  = MEDAKA_INPUT     which in nanopore mode is DECONTAM_CONTIGS,
    #                               i.e. the assembly AFTER the contamination
    #                               screen (D3)
    #
    # Model auto-inference reads the basecaller tag from FILT_LONG (the same reads
    # being polished), exactly as v1 BacFluxL did. An earlier draft pointed it at
    # the RAW ONT file instead: that is both a fidelity change and a new failure
    # mode, because the raw input may be gzipped (BacFluxL accepts fastq.gz/fq.gz)
    # while the filtlong output is always plain uncompressed FASTQ.
    # Produces: MEDAKA_DIR and MEDAKA_CONSENSUS (consensus.fasta, the name Medaka
    #           chooses itself).
    # Consumed by: finalize_contigs.
    #
    # params.model is the explicit model name from the config, or "" when the
    # user asked for automatic inference — that empty string is what the bash
    # if/else below branches on.
    #
    # (v1 message: "--- Medaka: Improve contig consensus with long reads. ---")
    rule long_read_consensus:
        input:
            reads = FILT_LONG,
            contigs = MEDAKA_INPUT,
            # The model NAME, already validated (explicit) or inferred (auto) by
            # check_medaka_model (shared/12_medaka_check.smk), which also gated the
            # assembler. Reading it here means the model is resolved once, and any
            # bad-model failure already happened before assembly — so this rule
            # trusts the name unconditionally and stays a plain Medaka call.
            model = MEDAKA_MODEL_RESOLVED,
        output:
            consensus_dir = directory(MEDAKA_DIR),
            consensus_contigs = MEDAKA_CONSENSUS,
        conda:
            "../../envs/medaka.yaml"
        threads: capped_cpus(24)
        log:
            LOGS + "/long_read_consensus_{sample}.log"
        priority: 9
        shell:
            # check_medaka_model already validated the NAME (or resolved auto), so
            # a failure here is a RUNTIME one — model weights not installed (no
            # network to fetch them), or a legacy model broken under Medaka v2.
            # The hint keeps that case actionable instead of a bare stack trace.
            """
            model=$(cat {input.model})
            echo "Polishing with Medaka model: $model" > {log}
            medaka_consensus \
              -i {input.reads} \
              -d {input.contigs} \
              -t {threads} \
              -m "$model" \
              -o {output.consensus_dir} >> {log} 2>&1 || {{
                cat {log} >&2
                echo "" >&2
                echo "Medaka failed at polishing with the pre-validated model '$model'." >&2
                echo "The name was accepted by check_medaka_model, so this is a RUNTIME failure: the model weights may not be installed locally (no network to fetch them), or a legacy model may be broken under Medaka v2." >&2
                echo "Set 'parameters.nanopore.medaka_model' to a supported model, or to FALSE to skip Medaka." >&2
                exit 1
              }}
            """


# ── Rule: finalize_contigs — publish the delivered genome under one name ─────
# Takes in: FINALIZE_SOURCE, a PARSE-TIME constant from 00_common that is
#           MEDAKA_CONSENSUS when Medaka is enabled and DECONTAM_CONTIGS when it
#           is not. Because it is resolved at parse time, the DAG is fixed before
#           the run starts and `snakemake -n` shows the real chain.
# Does:     copies the file and records, in the log, WHICH file it copied. That
#           one line is the provenance of the delivered genome.
# Produces: FINAL_CONTIGS — the single canonical hand-off (D2) that every shared
#           downstream module consumes.
# Consumed by: shared/15_replicons.smk, 20_qc, 30_taxonomy, 40_annotation,
#              50_amr, 60_plasmid, 70_phage, and blast_final_contigs in
#              shared/10_decontam.smk.
#
# conda: NONE — cp and echo only, as in v1.
rule finalize_contigs:
    input:
        contigs = FINALIZE_SOURCE,
    output:
        final = FINAL_CONTIGS,
    log:
        LOGS + "/finalize_contigs_{sample}.log"
    priority: 9
    shell:
        """
        mkdir -p $(dirname {output.final})

        echo "Final contigs source: {input.contigs}" | tee {log}

        cp {input.contigs} {output.final}
        """
