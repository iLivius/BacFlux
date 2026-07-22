# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — hybrid front end, short-read prep + raw ONT QC
# (rules/hybrid/10_reads.smk)
#
# Hybrid mode sequences the same isolate twice: Illumina for base accuracy, ONT
# for contiguity. This module prepares the Illumina side exactly as illumina mode
# does, and takes a first look at the raw ONT reads. The two legs meet later, in
# 30_ont_reads.smk.
#
# NOTE: the four Illumina rules below are IDENTICAL to rules/illumina/10_reads.smk
# — keep them in sync. They are duplicated rather than shared because only one
# mode's rule directory is ever included, so duplicate rule NAMES are safe, and
# because a reader working in hybrid mode should be able to read this mode's whole
# front end in one place. Both copies reference the same 00_common constants, so
# any drift is bounded to tool flags, never to paths.
#
# Data flow through this module:
#
#   NCBI ─► download_phix ─► build_phix ─► (6 bowtie2 index files)
#                                                   │
#   {sample}_R1/R2.fastq[.gz] ───────────────────────┴─► map_phix
#   (config input.illumina_dir)                            │
#                                                          ▼
#                                          {sample}.1/.2.fastq (PhiX-free)
#                                                          │
#                                                          ▼
#                                                  trim_adapters (fastp)
#                                       ┌──────────────────┴──────────────┐
#                                       ▼                                 ▼
#                              TRIM_R1 / TRIM_R2               FASTP_JSON / HTML
#                                       │                          (multiqc)
#            ┌──────────────────────────┼──────────────────────────┐
#            ▼                          ▼                          ▼
#   illumina_assembly (20)   map_contigs (shared/10)   map_sel_contigs (30) +
#                                                      map_amr_db (shared/50)
#
#   {sample}_ont.fastq[.gz] ─► raw_long_read_qc ─► NANOPLOT_RAW_DIR ─► multiqc
#   (config input.nanopore_dir)
#
# OWNERSHIP NOTE: `index_contigs` and `map_contigs` (trimmed reads onto the
# ILLUMINA draft, for BlobTools coverage) belong to shared/10_decontam.smk, not
# here. The separate `index_selected_contigs` / `map_sel_contigs` pair in
# 30_ont_reads.smk is a DIFFERENT job (reads onto the DECONTAMINATED contigs, to
# extract the clean read set) and is hybrid-only.
#
# Everything referenced here comes from 00_common.smk.
# ─────────────────────────────────────────────────────────────────────────────


# Where this mode's per-sample Illumina files live. Derived from TRIM_R1 rather
# than retyped — same anti-drift trick as SPADES_DIR in 00_common.
_READS_DIR = os.path.dirname(TRIM_R1)


# ── Rule: download_phix — fetch the PhiX control genome ──────────────────────
# NOTE: identical to rules/illumina/10_reads.smk — keep in sync.
#
# PhiX is the spike-in control present in essentially every Illumina lane. Its
# reads are real sequence from a different genome, so they must go before
# assembly. Runs ONCE per run (no {sample} in the path).
# Produces: PHIX_FASTA (temp). Consumed by: build_phix.
# conda: NONE — system wget, v1 parity (same deferred decision as
# cazyme_db_download in shared/40_annotation.smk).
rule download_phix:
    output:
        phix = temp(PHIX_FASTA),
    params:
        link = PHIX_LINK,
    log:
        LOGS + "/download_phix.log"
    priority: 10
    shell:
        """
        wget {params.link} -O {output.phix} > {log} 2>&1
        """


# ── Rule: build_phix — index the PhiX genome for Bowtie2 ─────────────────────
# NOTE: identical to rules/illumina/10_reads.smk — keep in sync.
# Takes in: PHIX_FASTA. Produces: six temp() index files. Consumed by: map_phix.
rule build_phix:
    input:
        phix = PHIX_FASTA,
    output:
        idx1  = temp(PHIX_BT2_PREFIX + ".1.bt2"),
        idx2  = temp(PHIX_BT2_PREFIX + ".2.bt2"),
        idx3  = temp(PHIX_BT2_PREFIX + ".3.bt2"),
        idx4  = temp(PHIX_BT2_PREFIX + ".4.bt2"),
        ridx1 = temp(PHIX_BT2_PREFIX + ".rev.1.bt2"),
        ridx2 = temp(PHIX_BT2_PREFIX + ".rev.2.bt2"),
    params:
        basename = PHIX_BT2_PREFIX,
    conda:
        "../../envs/bowtie.yaml"
    log:
        LOGS + "/build_phix.log"
    priority: 10
    shell:
        """
        bowtie2-build \
          {input.phix} \
          {params.basename} > {log} 2>&1
        """


# ── Rule: map_phix — throw away the reads that are PhiX ──────────────────────
# NOTE: identical to rules/illumina/10_reads.smk — keep in sync.
#
# Takes in: the six index files + this sample's RAW Illumina pair.
# Produces: {sample}.1.fastq / {sample}.2.fastq — the pairs that did NOT align
#           concordantly to PhiX, written by bowtie2's --un-conc. The SAM of the
#           reads that DID align is temp() and discarded.
# Consumed by: trim_adapters.
#
# --un-conc NAMING is load-bearing: from the basename "{sample}.fastq" bowtie2
# builds "{sample}.1.fastq" and "{sample}.2.fastq" by inserting .1/.2 before the
# final extension. Change params.basename and the declared outputs stop existing.
rule map_phix:
    input:
        idx1  = PHIX_BT2_PREFIX + ".1.bt2",
        idx2  = PHIX_BT2_PREFIX + ".2.bt2",
        idx3  = PHIX_BT2_PREFIX + ".3.bt2",
        idx4  = PHIX_BT2_PREFIX + ".4.bt2",
        ridx1 = PHIX_BT2_PREFIX + ".rev.1.bt2",
        ridx2 = PHIX_BT2_PREFIX + ".rev.2.bt2",
        r1 = os.path.join(ILLUMINA_DIR, R1),
        r2 = os.path.join(ILLUMINA_DIR, R2),
    output:
        sam = temp(_READS_DIR + "/{sample}_contam.sam"),
        r1  = temp(_READS_DIR + "/{sample}.1.fastq"),
        r2  = temp(_READS_DIR + "/{sample}.2.fastq"),
    params:
        db = PHIX_BT2_PREFIX,
        basename = _READS_DIR + "/{sample}.fastq",
    conda:
        "../../envs/bowtie.yaml"
    threads: CPUS
    log:
        LOGS + "/map_phix_{sample}.log"
    priority: 10
    shell:
        """
        bowtie2 \
          -x {params.db} \
          -1 {input.r1} -2 {input.r2} \
          --threads {threads} \
          --un-conc {params.basename} \
          -S {output.sam} \
          --local \
          > {log} 2>&1
        """


# ── Rule: trim_adapters — adapter removal and quality trimming (fastp) ───────
# NOTE: identical to rules/illumina/10_reads.smk — keep in sync.
#
# Takes in: the PhiX-free pairs from map_phix.
# Produces: TRIM_R1 / TRIM_R2 (temp — kept until the last consumer finishes) plus
#           the fastp JSON (MultiQC) and HTML (human) reports.
# Consumed by: illumina_assembly (20_assembly.smk), map_contigs
#              (shared/10_decontam.smk), map_sel_contigs (30_ont_reads.smk) and
#              map_amr_db (shared/50_amr.smk).
#
# v1 -> v2: hybrid v1 asked for min(CPUS, 24) here and BacFlux v1 for a hard-coded
# 16. Both short-read front ends now use capped_cpus(16): fastp is I/O-bound past
# roughly that point, and one value across the two modes is one thing less to
# reconcile. The hybrid env pin also moves from fastp 1.1.0 to the tree's existing
# fastp 1.0.1, the only env divergence the migration found.
rule trim_adapters:
    input:
        r1 = _READS_DIR + "/{sample}.1.fastq",
        r2 = _READS_DIR + "/{sample}.2.fastq",
    output:
        r1   = temp(TRIM_R1),
        r2   = temp(TRIM_R2),
        html = FASTP_HTML,
        json = FASTP_JSON,
    conda:
        "../../envs/fastp.yaml"
    threads: capped_cpus(16)
    log:
        LOGS + "/trim_adapters_{sample}.log"
    priority: 10
    shell:
        """
        fastp \
          --detect_adapter_for_pe \
          --length_required 100 \
          --cut_front \
          --cut_right \
          --thread {threads} \
          --verbose \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          -j {output.json} \
          -h {output.html} > {log} 2>&1
        """


# ── Rule: raw_long_read_qc — ONT length/quality profile BEFORE filtering ─────
# NOTE: identical to rules/nanopore/10_reads.smk — keep in sync.
#
# Takes in: this sample's raw ONT FASTQ.
# Produces: NANOPLOT_RAW_DIR (plots + {sample}_NanoStats.txt).
# Consumed by: multiqc (shared/90_report.smk).
#
# --prefix "{sample}_" IS LOAD-BEARING: shared/90_report.smk relabels this
# directory with a regex whose backreference requires the MultiQC sample name
# inside it to equal {sample}, which only holds if NanoPlot writes
# {sample}_NanoStats.txt. Without the prefix NanoPlot vanishes from the report,
# silently.
#
# The matching "after filtering" run lives in 30_ont_reads.smk, because in hybrid
# mode filtering cannot happen until the Illumina leg has produced its
# decontaminated reads.
rule raw_long_read_qc:
    input:
        fastq = os.path.join(NANOPORE_DIR, ONT),
    output:
        nanoplot_raw_dir = directory(NANOPLOT_RAW_DIR),
    conda:
        "../../envs/nanoplot.yaml"
    threads: capped_cpus(8)
    log:
        LOGS + "/raw_long_read_qc_{sample}.log"
    priority: 10
    shell:
        """
        NanoPlot \
          --fastq {input.fastq} \
          --threads {threads} \
          --loglength \
          --prefix "{wildcards.sample}_" \
          --outdir {output.nanoplot_raw_dir} > {log} 2>&1
        """
