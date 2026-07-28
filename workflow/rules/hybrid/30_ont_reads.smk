# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — hybrid front end, the coupling point
# (rules/hybrid/30_ont_reads.smk)
#
# THIS IS THE MODULE THAT MAKES HYBRID MODE HYBRID. Everything here is
# hybrid-only; none of these rules exists in the illumina or nanopore front ends.
#
# The idea, preserved exactly from v1 BacFluxL+: the ONT leg never gets its own
# contamination screen. Instead, the Illumina reads that map to the ALREADY
# DECONTAMINATED Illumina assembly are extracted, and those reads are handed to
# filtlong as the reference for scoring ONT reads. An ONT read that does not look
# like the clean Illumina genome scores badly and is dropped. So the Illumina leg
# IS the ONT leg's decontamination filter.
#
# Data flow through this module:
#
#   DECONTAM_CONTIGS ──► index_selected_contigs ──► (6 bowtie2 index files)
#   (shared/10_decontam)                                     │
#   TRIM_R1 / TRIM_R2 ───────────────────────────────────────┴──► map_sel_contigs
#   (10_reads.smk)                                                       │
#                                                    (samtools fastq -f 0x2)
#                                                                        ▼
#                                                             SEL_R1 / SEL_R2
#                                                     the decontaminated read pairs
#                                          ┌──────────────────────┴──────────────┐
#                                          ▼                                     ▼
#   raw ONT FASTQ ──────────────► filter_long_reads (filtlong -1 -2)    short_read_
#   (config input.nanopore_dir)              │                          correction
#                                            ▼                          (50_polish)
#                                        FILT_LONG
#                                            │
#                        ┌───────────────────┴────────────────┐
#                        ▼                                    ▼
#             filtered_long_read_qc                    ont_assembly (Flye)
#             -> NANOPLOT_FILT_DIR -> multiqc          (40_ont_assembly.smk)
#
# NOTE ON THE TWO BOWTIE2 INDEXES: shared/10_decontam.smk builds an index of the
# DRAFT assembly at DECONTAM_DIR + "/{sample}_contigs" (its module-local
# _BT2_PREFIX). This module builds a SECOND index, of the DECONTAMINATED assembly,
# at DECONTAM_DIR + "/{sample}_contigs_sel". Different prefixes, different
# purposes, no collision — v1 used exactly this pair.
#
# Everything referenced here comes from 00_common.smk: DECONTAM_CONTIGS,
# DECONTAM_DIR, TRIM_R1/TRIM_R2, SEL_R1/SEL_R2, NANOPORE_DIR, ONT, FILT_LONG,
# NANOPLOT_FILT_DIR, LOGS, CPUS, capped_cpus.
# ─────────────────────────────────────────────────────────────────────────────


# Prefix of the six Bowtie2 index files for the DECONTAMINATED Illumina assembly.
# Module-local: both the producer and the only consumer are in this file, so it is
# not a cross-module contract and does not belong in 00_common (same reasoning as
# _BT2_PREFIX in shared/10_decontam.smk). Built off DECONTAM_DIR so no stage
# number is re-derived.
_BT2_SEL_PREFIX = DECONTAM_DIR + "/{sample}_contigs_sel"


# ── Rule: index_selected_contigs — Bowtie2 index of the CLEAN assembly ───────
# Takes in: DECONTAM_CONTIGS, written by select_contigs in shared/10_decontam.smk
#           (in hybrid mode that file is contaminants/contigs_sel.fasta, not
#           FINAL_CONTIGS — see the per-mode table in that module's banner).
# Does:     bowtie2-build.
# Produces: six temp() index files.
# Consumed by: map_sel_contigs.
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
    priority: 6
    shell:
        """
        bowtie2-build \
          -f {input.contigs} \
          {params.basename} > {log} 2>&1
        """


# ── Rule: map_sel_contigs — extract the decontaminated Illumina read set ─────
# Biology: a read that aligns as a proper pair to the clean assembly is, by
# construction, a read from the organism we decided to keep. Pulling those reads
# back out gives a contamination-filtered read set without ever having to
# classify a read directly.
#
# Takes in: the six index files (the DAG edge to index_selected_contigs) and the
#           fastp-trimmed pairs TRIM_R1/TRIM_R2.
# Does:     bowtie2 with --no-unal (do not even write the unaligned reads), then
#           `samtools fastq -f 0x2`, where 0x2 is the SAM flag "read mapped in a
#           proper pair" — so only concordant pairs survive.
# Produces: SEL_R1 / SEL_R2, deliberately NOT temp(). Two reasons, both from v1:
#           the decontaminated read set is useful on its own, and if it were
#           deleted, re-running Polypolish alone would force SPAdes and the whole
#           screen to be re-run to recreate it.
# Consumed by: filter_long_reads (below) and short_read_correction (50_polish.smk).
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
    priority: 6
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


# ── Rule: filter_long_reads — score ONT reads against the clean Illumina reads ─
# Biology: filtlong scores every ONT read and keeps the best. Given -1/-2 it stops
# scoring by quality alone and asks "how much of this long read is supported by
# k-mers from the trusted short reads?" — which is simultaneously a quality filter
# and a decontamination filter, because a long read from a contaminant has no
# short-read support once the short reads have been decontaminated.
#
# The flag set differs from nanopore mode ON PURPOSE (both preserved from v1):
#   --min_length 1000   drop reads too short to help contiguity
#   --trim              cut low-support read ends rather than discarding the read
#   --split 1000        break a read at any unsupported stretch of 1000 bp, which
#                       is how a chimeric ONT read gets separated
#   --keep_percent 90   discard the worst 10% of the remaining bases
#   --length_weight 10  weight length heavily, since contiguity is the reason for
#                       sequencing ONT at all
# Nanopore mode instead uses --target_bases (a coverage cap) and no short reads,
# because it has none.
#
# Takes in: the raw ONT FASTQ plus SEL_R1/SEL_R2.
# Produces: FILT_LONG.
# Consumed by: filtered_long_read_qc (below), ont_assembly and
#              long_read_consensus (40_ont_assembly.smk), and map_contigs in
#              shared/10_decontam.smk — note that in HYBRID mode map_contigs takes
#              the SHORT-read branch, so FILT_LONG is not used there; only in
#              nanopore mode does that rule read it.
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
        # Was hard-coded to 10 - ten times filtlong's default - which destroyed
        # 85% of the ONT reads belonging to a 5.6 kb plasmid on the KPNIH-class
        # clinical isolates. See the note on FILTLONG_LENGTH_WEIGHT in
        # 00_common.smk for the measurement.
        length_weight = FILTLONG_LENGTH_WEIGHT,
    conda:
        "../../envs/filtlong.yaml"
    log:
        LOGS + "/filter_long_reads_{sample}.log"
    priority: 6
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


# ── Rule: filtered_long_read_qc — ONT profile AFTER filtering ────────────────
# NOTE: identical to rules/nanopore/10_reads.smk — keep in sync.
#
# Takes in: FILT_LONG. Produces: NANOPLOT_FILT_DIR. Consumed by: multiqc.
# Side by side with NANOPLOT_RAW_DIR (10_reads.smk) it shows exactly what the
# short-read-guided filter removed.
#
# --prefix "{sample}_" IS LOAD-BEARING — see the note in hybrid/10_reads.smk.
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
    priority: 6
    shell:
        """
        NanoPlot \
          --fastq {input.fastq} \
          --threads {threads} \
          --loglength \
          --prefix "{wildcards.sample}_" \
          --outdir {output.nanoplot_filt_dir} > {log} 2>&1
        """
