# BacFlux v2.0.0 — nanopore front end, polishing and the delivered genome.
#
# ONT reads are excellent at telling you the STRUCTURE of a genome and still
# imperfect at telling you the exact base, especially in homopolymers. Medaka
# re-reads the ONT reads against the draft assembly with a neural network trained
# on the same basecaller and rewrites the consensus. This is the step that turns a
# structurally-correct ONT assembly into one whose gene calls can be trusted.
#
# Stage chain: DECONTAM_CONTIGS + FILT_LONG → long_read_consensus →
# MEDAKA_CONSENSUS → finalize_contigs → FINAL_CONTIGS.
#
# long_read_consensus : Medaka. Exists only when Medaka is enabled, and polishes
#                       the DECONTAMINATED assembly.
# finalize_contigs    : copies whichever file ended the chain to one canonical
#                       name and logs which file that was.
#
# D3 ORDERING: Medaka polishes the DECONTAMINATED assembly, not the raw one. That
# is v1 nanopore behaviour and it is deliberate — polishing a contaminant contig
# with this isolate's reads would be wasted effort at best and would smear real
# differences at worst.
#
# WHY finalize_contigs EXISTS AT ALL: whether Medaka runs is a config choice, so
# the last file in this mode's chain is either the Medaka consensus or the
# decontaminated assembly. Copying whichever it is to ONE canonical name means
# every downstream module can depend on FINAL_CONTIGS and never ask the question.
#
# v1 → v2 changes in this file:
#  * v1 nanopore had a PAIR of anonymous rules (`rule:` with no name), selected by
#    an if/else on the config. v2 has ONE named rule under `if USE_MEDAKA:`.
#  * The model pre-flight check that v1 ran inside the polish rule has MOVED OUT
#    to check_medaka_model (shared/12_medaka_check.smk), which validates an
#    explicit model or infers one right after read filtering, gates the assembler
#    on the result, and writes the resolved name to a file. A bad model now fails
#    in seconds instead of after Flye, and this rule is a plain Medaka call that
#    reads the name back. The post-failure hint text in its shell is v1 hybrid's,
#    kept word for word.
#  * v1's `final_contigs(wc)` input FUNCTION (which re-read the config at DAG
#    build time) is replaced by the parse-time constant FINALIZE_SOURCE.
#
# Everything referenced here comes from 00_common.smk: USE_MEDAKA,
# MEDAKA_MODEL_RESOLVED, MEDAKA_INPUT (= DECONTAM_CONTIGS in this mode),
# MEDAKA_DIR, MEDAKA_CONSENSUS, FINALIZE_SOURCE, FINAL_CONTIGS, FILT_LONG, LOGS,
# capped_cpus.


# The rule below exists ONLY when Medaka is enabled — with Medaka off it is not in
# the DAG at all and FINALIZE_SOURCE points straight at DECONTAM_CONTIGS.
# USE_MEDAKA is resolved once at parse time in 00_common.smk section 8: it is
# False only when the user explicitly set parameters.nanopore.medaka_model to a
# false-like value. A missing or "auto" value means "run Medaka and work the model
# out yourself", NOT "skip".
if USE_MEDAKA:

    # ── Consensus polishing (Medaka) ──
    # Takes in:
    #   reads   = FILT_LONG             the same reads Flye assembled, from rule
    #                                   filter_long_reads (nanopore/10_reads.smk).
    #   contigs = MEDAKA_INPUT          in nanopore mode that is DECONTAM_CONTIGS,
    #                                   the assembly AFTER the contamination screen
    #                                   (D3), written by rule select_contigs
    #                                   (shared/10_decontam.smk).
    #   model   = MEDAKA_MODEL_RESOLVED the one-line file holding the model NAME,
    #                                   from check_medaka_model
    #                                   (shared/12_medaka_check.smk) — see the
    #                                   inline note on the input below.
    # Does:     medaka_consensus aligns those reads back to that draft and rewrites
    #           the consensus with the network trained on the matching basecaller.
    # Produces:
    #   consensus_dir     = MEDAKA_DIR, 02.assembly/{sample}/medaka/
    #   consensus_contigs = MEDAKA_CONSENSUS, consensus.fasta — the filename is
    #                       Medaka's choice, not ours
    # Consumed by: finalize_contigs below, which reaches it through FINALIZE_SOURCE.
    #
    # WHY the model is inferred from FILT_LONG and not from the raw ONT file (that
    # happens in check_medaka_model, not here): the raw input may be gzipped, while
    # the filtlong output is always plain uncompressed FASTQ, and the basecaller tag
    # in the read headers is the same in both. v1 BacFluxL did it this way too.
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


# ──────────────────────── Final contigs hand-off ───────────────
# Takes in: contigs = FINALIZE_SOURCE, a PARSE-TIME constant from 00_common.smk:
#           MEDAKA_CONSENSUS (from long_read_consensus above) when Medaka is
#           enabled, DECONTAM_CONTIGS (from select_contigs, shared/10_decontam.smk)
#           when it is not. Because it is resolved before the run starts, the DAG
#           is fixed and `snakemake -n` shows the real chain.
# Does:     copies that file to one canonical name, and logs WHICH file it copied
#           — that one line is the provenance of the delivered genome.
# Produces: final = FINAL_CONTIGS, 02.assembly/{sample}/contigs_final.fasta, the
#           single canonical hand-off (D2).
# Consumed by: build_replicons (shared/15_replicons.smk); stage_qc_genomes
#              (shared/20_qc.smk), which is also how taxonomic_assignment
#              (shared/30_taxonomy.smk) gets it; annotation
#              (shared/40_annotation.smk); amr_contigs (shared/50_amr.smk);
#              plasmid_search (shared/60_plasmid.smk); the phage caller
#              (viral_identification_virsorter2, or genomad_end_to_end when
#              PHAGE_CALLER is genomad — shared/70_phage.smk);
#              blast_final_contigs (shared/10_decontam.smk); and, when the
#              mobilome module is on, amrfinderplus, isescan and six more rules
#              in shared/80_mobilome.smk.
#
# conda: NONE — cp and echo only, as in v1.
rule finalize_contigs:
    input:
        contigs = FINALIZE_SOURCE,
    output:
        final = FINAL_CONTIGS,
    log:
        LOGS + "/finalize_contigs_{sample}.log"
    shell:
        """
        mkdir -p $(dirname {output.final})

        echo "Final contigs source: {input.contigs}" | tee {log}

        cp {input.contigs} {output.final}
        """
