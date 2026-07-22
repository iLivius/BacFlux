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
        output:
            consensus_dir = directory(MEDAKA_DIR),
            consensus_contigs = MEDAKA_CONSENSUS,
        params:
            model = MEDAKA_MODEL if MEDAKA_MODEL is not None else "",
        conda:
            "../../envs/medaka.yaml"
        resources:
            cpus = capped_cpus(24)
        log:
            LOGS + "/long_read_consensus_{sample}.log"
        priority: 9
        shell:
            # A NOTE ON BACKSLASHES BELOW: this shell block is an ordinary Python
            # string, so a backslash-n written once would become a real newline
            # before bash ever sees it. Where bash itself needs the two characters
            # \n (in printf and tr), they are written doubled: \\n.
            """
            if [ -n "{params.model}" ]; then

              # ── The user named a model explicitly ─────────────────────────
              model="{params.model}"

              # PRE-FLIGHT CHECK (kept from v1 nanopore): a wrong model name is
              # the most common Medaka failure, and without this check it only
              # shows up after Medaka has already loaded the reads. Skip the
              # check when the value points at a local model file on disk.
              if [ ! -e "$model" ]; then
                models=$(medaka tools list_models 2>&1) || {{
                  {{
                    echo "ERROR: Unable to query Medaka models in the current conda environment."
                    echo "  medaka_model: $model"
                    echo "  reason: 'medaka tools list_models' failed."
                    echo ""
                    echo "$models"
                  }} > {log}
                  cat {log} >&2
                  exit 1
                }}
                if ! printf '%s\\n' "$models" | sed -n 's/^Available: //p' | tr ',' '\\n' | sed 's/^ *//; s/ *$//' | grep -Fxq "$model"; then
                  {{
                    echo "ERROR: Invalid Medaka model configured for consensus polishing."
                    echo "  medaka_model: $model"
                    echo "  reason: this model is not available in the current Medaka environment."
                    echo "  next steps: choose a model from 'medaka tools list_models', set 'parameters.nanopore.medaka_model' to auto to infer it from the FASTQ headers, set it to FALSE to skip Medaka, or use a Medaka 1.x environment if this exact legacy model is required."
                  }} > {log}
                  cat {log} >&2
                  exit 1
                fi
              fi

              medaka_consensus \
                -i {input.reads} \
                -d {input.contigs} \
                -t {resources.cpus} \
                -m "$model" \
                -o {output.consensus_dir} > {log} 2>&1 || {{
                  cat {log}
                  echo "" >&2
                  echo "Medaka failed while using the explicit model setting '$model'." >&2
                  echo "If this is an older ONT model, it may be deprecated in Medaka v2 or its model weights may not be installed locally." >&2
                  echo "For Medaka v2, prefer a supported model such as 'r941_min_fast_g507' over deprecated legacy names like 'r941_min_fast_g303'." >&2
                  echo "Otherwise set 'parameters.nanopore.medaka_model' to FALSE to skip Medaka." >&2
                  exit 1
                }}

            else

              # ── Infer the model from the read headers ─────────────────────
              # ONT basecallers stamp the model into the FASTQ header, and
              # filtlong preserves headers on the reads it keeps — so the filtered
              # file carries the tag just as v1 BacFluxL relied on. Using the same
              # reads being polished also keeps this immune to the raw input's
              # compression (the raw file may be .gz; this one never is).
              resolved_model=$(medaka tools resolve_model --auto_model consensus_bacteria {input.reads} 2> {log}) || {{
                cat {log}
                echo "" >&2
                echo "Medaka could not auto-infer a consensus model from {input.reads}." >&2
                echo "This usually means the ONT FASTQ headers do not contain exactly one basecaller model reference." >&2
                echo "Set 'parameters.nanopore.medaka_model' in the config to auto, an explicit Medaka model name, or FALSE to skip Medaka." >&2
                exit 1
              }}
              echo "Resolved Medaka model: $resolved_model" >> {log}

              medaka_consensus \
                -i {input.reads} \
                -d {input.contigs} \
                -t {resources.cpus} \
                -m "$resolved_model" \
                -o {output.consensus_dir} >> {log} 2>&1

            fi
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
