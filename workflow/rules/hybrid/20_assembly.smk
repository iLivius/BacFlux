# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — hybrid front end, Illumina draft assembly
# (rules/hybrid/20_assembly.smk)
#
# In hybrid mode the Illumina assembly is NOT the deliverable. It has two jobs:
#   1. it is what the contamination screen actually screens, and the reads that
#      map to the CLEAN version of it become the read set that filters the ONT
#      reads (30_ont_reads.smk). So the Illumina leg IS the ONT leg's
#      decontamination filter — the ONT leg has no screen of its own.
#   2. it is kept as the comparator genome: QC and taxonomy are run on both the
#      Illumina draft and the delivered ONT genome, and Snippy compares them.
#
# NOTE: both rules below are IDENTICAL to rules/illumina/20_assembly.smk — keep
# them in sync. Duplicated deliberately (see the banner in hybrid/10_reads.smk).
#
# Data flow through this module:
#
#   TRIM_R1 / TRIM_R2 ──► illumina_assembly (SPAdes) ──► SPADES_CONTIGS
#   (10_reads.smk)                                            │
#                                                             ▼
#                                                      filter_contigs
#                                                             │
#                                                             ▼
#                                                       DRAFT_CONTIGS
#                                                             │
#                              shared/10_decontam.smk (BLAST + BlobTools + selector)
#                                                             │
#                                                             ▼
#                                    DECONTAM_CONTIGS = contaminants/contigs_sel.fasta
#                             ┌───────────────┬───────────────┴──────────────┐
#                             ▼               ▼                              ▼
#                index_selected_contigs   Snippy --ref            QC comparator genome
#                (30_ont_reads.smk)       (50_polish.smk)         "{sample}_illumina"
#                                                                 (shared/20_qc.smk)
#
# THIS MODULE MUST NOT WRITE FINAL_CONTIGS. In hybrid mode FINAL_CONTIGS is the
# ONT+Polypolish genome, written by short_read_correction in 50_polish.smk.
#
# Everything referenced here comes from 00_common.smk.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: illumina_assembly — de novo short-read assembly (SPAdes) ───────────
# NOTE: identical to rules/illumina/20_assembly.smk — keep in sync.
#
# Takes in: TRIM_R1 / TRIM_R2 (fastp-trimmed, PhiX-free).
# Does:     SPAdes --isolate over k = 21..127. OMP_NUM_THREADS is exported so the
#           OpenMP sections obey the same budget as -t and do not grab every core.
# Produces: SPADES_DIR + contigs.fasta inside it.
# Consumed by: filter_contigs, below.
#
# LAYOUT LANDMINE: declare directory(SPADES_DIR), never the sample directory. In
# the v2 layout 02.assembly/{sample}/ also holds contaminants/, eval/, flye/,
# fix_start/, medaka/, polypolish/ and snps/, all written by other rules, and a
# directory() output is deleted and recreated whenever the rule re-runs.
rule illumina_assembly:
    input:
        r1 = TRIM_R1,
        r2 = TRIM_R2,
    output:
        dir = directory(SPADES_DIR),
        contigs = SPADES_CONTIGS,
    conda:
        "../../envs/spades.yaml"
    threads: CPUS
    resources:
        ram = RAM
    log:
        LOGS + "/illumina_assembly_{sample}.log"
    priority: 10
    shell:
        """
        OMP_NUM_THREADS={threads} \
        spades.py -k 21,33,55,77,99,127 --isolate \
          --pe1-1 {input.r1} \
          --pe1-2 {input.r2} \
          -o {output.dir} \
          -t {threads} \
          -m {resources.ram} > {log} 2>&1
        """


# ── Rule: filter_contigs — drop short and low-coverage contigs ───────────────
# NOTE: identical to rules/illumina/20_assembly.smk — keep in sync.
#
# Takes in: SPADES_CONTIGS.
# Does:     FASTA_LIN_CMD puts each sequence on one line; FASTA_SEL_CMD then keeps
#           records whose SPAdes header says length >= 500 (field 4 when the
#           header is split on "_") and coverage >= 2.0 (field 6). Both awk
#           programs live in 00_common and are passed with {NAME:q} so Snakemake
#           quotes them for the shell.
# Produces: DRAFT_CONTIGS — the input to the shared contamination screen.
# Consumed by: index_contigs, map_contigs, blast_contigs, blob_json,
#              select_contigs (all shared/10_decontam.smk).
#
# conda: NONE — cat and awk only.
rule filter_contigs:
    input:
        contigs = SPADES_CONTIGS,
    output:
        contigs = DRAFT_CONTIGS,
    priority: 9
    shell:
        """
        cat {input.contigs} | \
        awk {FASTA_LIN_CMD:q} | \
        awk -F"_" {FASTA_SEL_CMD:q} > {output.contigs}
        """
