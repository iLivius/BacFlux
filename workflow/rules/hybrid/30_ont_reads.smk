# The module that makes hybrid mode hybrid. Everything here is hybrid-only; none
# of these rules exists in the illumina or nanopore front ends.
#
# The idea, preserved exactly from v1 BacFluxL+: the ONT leg never gets a
# contamination screen of its own. Instead the Illumina reads that map to the
# ALREADY DECONTAMINATED Illumina assembly are pulled back out, and those reads
# are handed to filtlong as the reference for scoring ONT reads. An ONT read that
# does not look like the clean Illumina genome scores badly and is dropped. So
# the Illumina leg IS the ONT leg's decontamination filter.
#
# Chain: DECONTAM_CONTIGS (shared/10_decontam.smk) → index_selected_contigs →
# map_sel_contigs, which also takes TRIM_R1/TRIM_R2 (hybrid/10_reads.smk) and
# emits SEL_R1/SEL_R2, the decontaminated read pairs. SEL_R1/SEL_R2 then feed
# filter_long_reads (with the raw ONT FASTQ) and short_read_correction
# (hybrid/50_polish.smk). filter_long_reads emits FILT_LONG →
# filtered_long_read_qc and ont_assembly (hybrid/40_ont_assembly.smk).
#
# There are TWO Bowtie2 indexes in play and they are not the same index.
# shared/10_decontam.smk builds one of the DRAFT assembly, at
# DECONTAM_DIR + "/{sample}_contigs" (its module-local _BT2_PREFIX). This module
# builds a second one, of the DECONTAMINATED assembly, at
# DECONTAM_DIR + "/{sample}_contigs_sel". Different prefixes, different purposes,
# no collision — v1 used exactly this pair.
#
# index_selected_contigs : bowtie2-build over the decontaminated assembly.
# map_sel_contigs        : map the trimmed pairs to it and keep only concordant
#                          pairs — that read set is the decontaminated Illumina
#                          data, and it is not temp().
# filter_long_reads      : filtlong, scoring ONT reads against those short reads.
# filtered_long_read_qc  : NanoPlot profile of what survived.
#
# Every constant used here comes from shared/00_common.smk: DECONTAM_CONTIGS,
# DECONTAM_DIR, TRIM_R1/TRIM_R2, SEL_R1/SEL_R2, NANOPORE_DIR, ONT, FILT_LONG,
# NANOPLOT_FILT_DIR, LOGS, CPUS, capped_cpus.


# Prefix of the six Bowtie2 index files for the DECONTAMINATED Illumina assembly.
# Module-local: the producer and the only consumer are both in this file, so it
# is not a cross-module contract and does not belong in 00_common (the same
# reasoning as _BT2_PREFIX in shared/10_decontam.smk). Built off DECONTAM_DIR so
# no stage number is re-derived.
_BT2_SEL_PREFIX = DECONTAM_DIR + "/{sample}_contigs_sel"


# ────────────────────── Clean Illumina read set (Bowtie2) ──────
# Takes in: contigs = DECONTAM_CONTIGS, written by select_contigs in
#           shared/10_decontam.smk — in hybrid mode that file is
#           contaminants/contigs_sel.fasta, not FINAL_CONTIGS (the per-mode table
#           in that module's header spells this out).
# Does:     bowtie2-build, under the module-local prefix _BT2_SEL_PREFIX above.
# Produces: six temp() index files — .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2
#           and .rev.2.bt2.
# Consumed by: map_sel_contigs, below.
#
# (v1 message: "--- Bowtie2: Build selected contig db. ---")
rule index_selected_contigs:
    input:
        contigs = DECONTAM_CONTIGS,
    output:
        idx1  = temp(_BT2_SEL_PREFIX + ".1.bt2"),
        idx2  = temp(_BT2_SEL_PREFIX + ".2.bt2"),
        idx3  = temp(_BT2_SEL_PREFIX + ".3.bt2"),
        idx4  = temp(_BT2_SEL_PREFIX + ".4.bt2"),
        ridx1 = temp(_BT2_SEL_PREFIX + ".rev.1.bt2"),
        ridx2 = temp(_BT2_SEL_PREFIX + ".rev.2.bt2"),
    params:
        basename = _BT2_SEL_PREFIX,
    conda:
        "../../envs/bowtie.yaml"
    log:
        LOGS + "/index_selected_contigs_{sample}.log"
    shell:
        """
        bowtie2-build \
          -f {input.contigs} \
          {params.basename} > {log} 2>&1
        """


# ─────────────────── Extracting the clean pairs (Bowtie2) ──────
# Pull the decontaminated Illumina read set back out of the clean assembly. A
# read that aligns as a proper pair to that assembly is, by construction, a read
# from the organism the screen decided to keep, so extracting those reads gives a
# contamination-filtered read set without ever classifying a read directly.
#
# Takes in: the six index files from index_selected_contigs, which is what
#           creates the DAG edge to it, plus r1 / r2 = TRIM_R1 / TRIM_R2, the
#           fastp-trimmed pairs from trim_adapters (hybrid/10_reads.smk).
# Does:     bowtie2 with --no-unal, so the unaligned reads are never even
#           written, then `samtools fastq -f 0x2`, which keeps only records
#           carrying SAM flag 0x2, "read mapped in a proper pair" — concordant
#           pairs only.
# Produces:
#   sam    = {sample}_map.sam under DECONTAM_DIR, the alignment itself. temp().
#   r1, r2 = SEL_R1 / SEL_R2, the decontaminated Illumina pairs. Deliberately NOT
#            temp(), for two reasons both inherited from v1: the decontaminated
#            read set is a useful deliverable on its own, and if it were deleted,
#            re-running Polypolish alone would drag SPAdes and the whole
#            contamination screen back through the DAG to recreate it.
# Consumed by: filter_long_reads below, and short_read_correction
#              (hybrid/50_polish.smk).
#
# (v1 message: "--- Bowtie2: Map trimmed reads against selected contigs. ---")
rule map_sel_contigs:
    input:
        idx1  = _BT2_SEL_PREFIX + ".1.bt2",
        idx2  = _BT2_SEL_PREFIX + ".2.bt2",
        idx3  = _BT2_SEL_PREFIX + ".3.bt2",
        idx4  = _BT2_SEL_PREFIX + ".4.bt2",
        ridx1 = _BT2_SEL_PREFIX + ".rev.1.bt2",
        ridx2 = _BT2_SEL_PREFIX + ".rev.2.bt2",
        r1 = TRIM_R1,
        r2 = TRIM_R2,
    output:
        sam = temp(DECONTAM_DIR + "/{sample}_map.sam"),
        r1 = SEL_R1,
        r2 = SEL_R2,
    params:
        db = _BT2_SEL_PREFIX,
    conda:
        "../../envs/bowtie.yaml"
    threads: CPUS
    log:
        LOGS + "/map_sel_contigs_{sample}.log"
    shell:
        """
        bowtie2 \
          -x {params.db} \
          -1 {input.r1} \
          -2 {input.r2} \
          -p {threads} \
          -t \
          --no-unal \
          -S {output.sam} > {log} 2>&1

        samtools fastq \
          -f 0x2 \
          -1 {output.r1} \
          -2 {output.r2} \
          {output.sam} >> {log} 2>&1
        """


# ──────────────────────── Long-read filtering (filtlong) ───────
# filtlong scores every ONT read and keeps the best of them. Given -1/-2 it stops
# scoring by quality alone and asks instead "how much of this long read is
# supported by k-mers from the trusted short reads?" — which is simultaneously a
# quality filter and a decontamination filter, because a long read from a
# contaminant has no short-read support once the short reads themselves have been
# decontaminated.
#
# Takes in: long   = the raw ONT FASTQ, ONT under NANOPORE_DIR, i.e.
#                    {sample}_ont.fastq[.gz] in config input.nanopore_dir
#           r1, r2 = SEL_R1 / SEL_R2, the decontaminated Illumina pairs from
#                    map_sel_contigs above
# Does:     filtlong, with a flag set that differs from nanopore mode on purpose:
#   --min_length      FILTLONG_MIN_LENGTH — drop reads too short to help
#                     contiguity. Set by parameters.long_read_qc.min_length,
#                     default 1000.
#   --trim            cut low-support read ENDS rather than discard the read.
#   --split 1000      break a read at any unsupported stretch of 1000 bp, which
#                     is how a chimeric ONT read gets separated. The one number
#                     still hard-coded in this rule.
#   --keep_percent    FILTLONG_KEEP_PERCENT — how much of the remaining sequence
#                     to keep, by score. Set by
#                     parameters.long_read_qc.keep_percent, default 95.
#   --length_weight   FILTLONG_LENGTH_WEIGHT — how much length dominates the
#                     ranking. Set by parameters.long_read_qc.length_weight,
#                     default 1, which is filtlong's own default.
# Produces: filt_long = FILT_LONG, the ONT reads that survived.
# Consumed by: filtered_long_read_qc below, ont_assembly and long_read_consensus
#              (hybrid/40_ont_assembly.smk), and check_medaka_model
#              (shared/12_medaka_check.smk).
#
# Nanopore mode instead uses --target_bases (a plain coverage cap) and no short
# reads, because it has none.
#
# map_contigs in shared/10_decontam.smk also reads FILT_LONG, but only in
# nanopore mode: in HYBRID mode that rule takes the short-read branch, so it
# never touches this file.
#
# (v1 message: "--- Filtlong: Filter and short-read correct long reads. ---")
rule filter_long_reads:
    input:
        long = os.path.join(NANOPORE_DIR, ONT),
        r1 = SEL_R1,
        r2 = SEL_R2,
    output:
        filt_long = FILT_LONG,
    params:
        min_length = FILTLONG_MIN_LENGTH,
        split = 1000,
        keep_percent = FILTLONG_KEEP_PERCENT,
        # keep_percent (the line above) and length_weight were both hard-coded
        # here once, at 90 and 10, and together they destroyed the 5,596 bp Col2
        # plasmid of K. pneumoniae TUM24772: 93 of the 603 ONT reads covering it
        # survived. keep_percent is the one to reach for — at 95, the value that
        # now ships, 500 survive even with length_weight left at 10 — but both
        # keys matter: at the old 90, dropping length_weight to 1 on its own also
        # recovers most of them (413). length_weight 1 additionally protects a
        # size class keep_percent cannot: at 10 the 1-3 kb read band is emptied
        # whatever keep_percent says. Both are config keys now, defaulting to 95
        # and 1. The full 2x2, the read-length bands, the strain's accessions and
        # what Ryan Wick's Feb 2026 benchmark does and does not say about this all
        # sit next to FILTLONG_LENGTH_WEIGHT in shared/00_common.smk.
        length_weight = FILTLONG_LENGTH_WEIGHT,
    conda:
        "../../envs/filtlong.yaml"
    log:
        LOGS + "/filter_long_reads_{sample}.log"
    shell:
        """
        # The selected Illumina reads define the target genome, so this single
        # step does the enrichment AND the decontamination for the ONT leg; there
        # is deliberately no separate ONT screening block after Flye/dnaapler.
        filtlong \
          -1 {input.r1} \
          -2 {input.r2} \
          --min_length {params.min_length} \
          --trim \
          --split {params.split} \
          --keep_percent {params.keep_percent} \
          --length_weight {params.length_weight} \
          {input.long} > {output.filt_long} 2>{log}
        """


# ──────────────────── Filtered ONT read profile (NanoPlot) ─────
# Takes in: fastq = FILT_LONG, from filter_long_reads above.
# Does:     the same NanoPlot run as the raw profile in hybrid/10_reads.smk, this
#           time over the reads that survived filtering.
# Produces: nanoplot_filt_dir = NANOPLOT_FILT_DIR, a directory of plots plus
#           {sample}_NanoStats.txt.
# Consumed by: multiqc (shared/90_report.smk). Read side by side with
#              NANOPLOT_RAW_DIR from hybrid/10_reads.smk it shows exactly what
#              the short-read-guided filter removed.
#
# Identical to rules/nanopore/10_reads.smk. Keep the two in step: a change to
# the filtlong invocation here belongs there too.
#
# `--prefix "{sample}_"` is load-bearing for the same silent-failure reason given
# in hybrid/10_reads.smk: the MultiQC relabel regex in shared/90_report.smk needs
# the sample name inside the NanoStats filename.
rule filtered_long_read_qc:
    input:
        fastq = FILT_LONG,
    output:
        nanoplot_filt_dir = directory(NANOPLOT_FILT_DIR),
    conda:
        "../../envs/nanoplot.yaml"
    threads: capped_cpus(8)
    log:
        LOGS + "/filtered_long_read_qc_{sample}.log"
    shell:
        """
        NanoPlot \
          --fastq {input.fastq} \
          --threads {threads} \
          --loglength \
          --prefix "{wildcards.sample}_" \
          --outdir {output.nanoplot_filt_dir} > {log} 2>&1
        """
