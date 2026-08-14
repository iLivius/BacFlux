# Hybrid front end, Illumina draft assembly: SPAdes over the trimmed pairs, then
# a length and coverage filter.
#
# In hybrid mode the Illumina assembly is NOT the deliverable. It has two jobs.
# First, it is what the contamination screen actually screens, and the reads that
# map back to the CLEAN version of it become the read set that filters the ONT
# reads (hybrid/30_ont_reads.smk) — so the Illumina leg IS the ONT leg's
# decontamination filter, and the ONT leg has no screen of its own. Second, it is
# kept as the comparator genome: QC and taxonomy run over both the Illumina draft
# and the delivered ONT genome, and Snippy compares them.
#
# Chain: TRIM_R1/TRIM_R2 (hybrid/10_reads.smk) → illumina_assembly →
# SPADES_CONTIGS → filter_contigs → DRAFT_CONTIGS → the shared screen (BLAST +
# BlobTools + selector, shared/10_decontam.smk) → DECONTAM_CONTIGS
# (contaminants/contigs_sel.fasta). From there DECONTAM_CONTIGS is read by
# index_selected_contigs (hybrid/30_ont_reads.smk), by Snippy's --ref in
# hybrid/50_polish.smk, and staged as the "{sample}_illumina" comparator genome
# by shared/20_qc.smk.
#
# This module must NOT write FINAL_CONTIGS. In hybrid mode FINAL_CONTIGS is the
# ONT+Polypolish genome, written by short_read_correction in hybrid/50_polish.smk.
#
# Both rules are identical to rules/illumina/20_assembly.smk — keep them in sync.
# Duplicated deliberately, for the reasons given in the header of
# hybrid/10_reads.smk.
#
# illumina_assembly : SPAdes --isolate over k = 21..127.
# filter_contigs    : drop contigs under 500 bp or under 2x coverage, using the
#                     length and coverage SPAdes already wrote into each header.
#
# Every constant used here comes from shared/00_common.smk.


# ──────────────────────── Draft assembly (SPAdes) ──────────────
# Takes in: r1 / r2 = TRIM_R1 / TRIM_R2, the fastp-trimmed, PhiX-free pairs
#           written by rule trim_adapters (hybrid/10_reads.smk).
# Does:     SPAdes --isolate over k = 21,33,55,77,99,127. OMP_NUM_THREADS is
#           exported so the OpenMP sections obey the same budget as -t instead of
#           grabbing every core on the machine.
# Produces:
#   dir     = SPADES_DIR, SPAdes' own output tree
#   contigs = SPADES_CONTIGS, the contigs.fasta inside it
# Consumed by: filter_contigs, below, which reads SPADES_CONTIGS. Nothing reads
#              SPADES_DIR — it is declared so Snakemake owns SPAdes' whole tree.
#
# Identical to rules/illumina/20_assembly.smk — keep in sync.
#
# Declare directory(SPADES_DIR), never the sample directory. Under the v2 layout
# 02.assembly/{sample}/ also holds contaminants/, eval/, flye/, fix_start/,
# medaka/, polypolish/ and snps/, all written by other rules — and Snakemake
# deletes and recreates a directory() output whenever its rule re-runs, so
# declaring the parent would wipe those out.
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


# ──────────────────────── Length and coverage filter ───────────
# Drop the short, thin contigs SPAdes leaves behind before anything downstream
# has to reason about them.
#
# Takes in: contigs = SPADES_CONTIGS, from rule illumina_assembly above.
# Does:     FASTA_LIN_CMD puts each sequence on one line, then FASTA_SEL_CMD
#           keeps records whose SPAdes header says length >= 500 (field 4 when
#           the header is split on "_") and coverage >= 2.0 (field 6). Both awk
#           programs live in shared/00_common.smk and are passed with {NAME:q} so
#           Snakemake quotes them for the shell.
# Produces: contigs = DRAFT_CONTIGS, the filtered draft and the input to the
#           shared contamination screen.
# Consumed by: index_contigs, blast_contigs, blob_json and select_contigs, all in
#              shared/10_decontam.smk. map_contigs reads the Bowtie2 index that
#              index_contigs builds from this file, not the FASTA itself.
#
# Identical to rules/illumina/20_assembly.smk — keep in sync.
#
# No conda environment: cat and awk only.
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
