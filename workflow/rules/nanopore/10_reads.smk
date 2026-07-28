# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — nanopore front end, read preparation (rules/nanopore/10_reads.smk)
#
# The biology: an ONT run produces reads of wildly varying length and quality,
# and far more data than a bacterial isolate needs. Assembling everything is
# slower AND worse: short, poor reads add errors without adding contiguity. So we
# keep the longest, best subset, and we look at the length/quality distribution
# before and after so the choice is visible in the report.
#
# Data flow through this module:
#
#   {sample}_ont.fastq[.gz] ──┬──► raw_long_read_qc ──► NANOPLOT_RAW_DIR ──┐
#   (raw input, config          │                                          │
#    input.nanopore_dir)        └──► filter_long_reads (filtlong)          │
#                                              │                           │
#                                              ▼                           │
#                                          FILT_LONG                       │
#                                    ┌─────────┼─────────┐                 ▼
#                                    ▼         ▼         ▼             multiqc
#                          filtered_long_    Flye    map_contigs      (shared/90)
#                            read_qc      (20_asm)  (shared/10_decontam)
#                                    │
#                                    ▼
#                           NANOPLOT_FILT_DIR ──────────────────────────► multiqc
#
# There is NO PhiX step here: PhiX is an Illumina spike-in and does not exist in
# an ONT library.
#
# Everything referenced here comes from 00_common.smk: NANOPORE_DIR, ONT,
# FILT_LONG, NANOPLOT_RAW_DIR, NANOPLOT_FILT_DIR, LOGS, capped_cpus.
#
# conda: paths resolve relative to THIS file (workflow/rules/nanopore/), so
# "../../envs/x.yaml" climbs nanopore/ -> rules/ -> workflow/ -> workflow/envs/.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: raw_long_read_qc — length/quality profile BEFORE filtering ─────────
# Takes in: this sample's raw ONT FASTQ from config input.nanopore_dir.
# Does:     NanoPlot, which summarises read length, quality and yield.
# Produces: NANOPLOT_RAW_DIR (a DIRECTORY of plots plus {sample}_NanoStats.txt).
# Consumed by: multiqc (shared/90_report.smk), and nothing else.
#
# --prefix "{sample}_" IS LOAD-BEARING, not cosmetic. MultiQC names a sample
# after the file it read, and shared/90_report.smk rewrites those names with the
# regex '^01\\.reads \\| ([^|]+) \\| ont \\| raw_qc \\| \\1$' — the backreference
# requires the name INSIDE the directory to equal the sample name, which only
# happens when NanoPlot writes {sample}_NanoStats.txt. Drop the prefix and
# NanoPlot silently disappears from the report.
#
# Rule NAME: v1 BacFluxL called this `raw_read_qc`; we adopt BacFluxL+'s
# `raw_long_read_qc` so the nanopore and hybrid front ends read identically.
#
# (v1 message: "--- NanoPlot: Raw long-read QC. ---")
rule raw_long_read_qc:
    input:
        fastq = os.path.join(NANOPORE_DIR, ONT),
    output:
        nanoplot_raw_dir = directory(NANOPLOT_RAW_DIR),
    conda:
        "../../envs/nanoplot.yaml"
        # NanoPlot is plotting, not aligning; it stops scaling early.
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


# ── Rule: filter_long_reads — keep the long, high-quality subset (filtlong) ──
# Biology / the three thresholds, all preserved from v1 BacFluxL:
#   --min_length 1000     reads under 1 kb bring ONT's error rate without ONT's
#                         main benefit (spanning repeats), so they are dropped.
#   --keep_percent 90     discard the worst-scoring 10% of the remaining bases.
#   --target_bases 5e8    stop at 500 Mbp, which is roughly 100x for a typical
#                         5 Mbp bacterial genome. More coverage than that costs
#                         Flye time without improving the assembly.
# Note this is the NANOPORE-mode flag set. Hybrid mode uses a different one
# (--trim --split --length_weight, and no --target_bases) because there the
# selected Illumina reads, not a coverage cap, define what is worth keeping.
#
# Takes in: the raw ONT FASTQ.
# Produces: FILT_LONG.
# Consumed by: filtered_long_read_qc (below), ont_assembly (20_assembly.smk),
#              long_read_consensus (30_polish.smk) and map_contigs
#              (shared/10_decontam.smk, the ONT coverage track for BlobTools).
#
# Rule NAME: v1 BacFluxL called this `filter_reads`; renamed to BacFluxL+'s
# `filter_long_reads` for the same read-alike reason as above.
#
# (v1 message: "--- Filtlong: Filter long reads. ---")
rule filter_long_reads:
    input:
        long = os.path.join(NANOPORE_DIR, ONT),
    output:
        filt_long = FILT_LONG,
    params:
        min_length = FILTLONG_MIN_LENGTH,
        keep_percent = FILTLONG_KEEP_PERCENT,
        target_bases = 500000000,
    conda:
        "../../envs/filtlong.yaml"
    log:
        LOGS + "/filter_long_reads_{sample}.log"
    priority: 10
    shell:
        """
        filtlong \
          --min_length {params.min_length} \
          --keep_percent {params.keep_percent} \
          --target_bases {params.target_bases} \
          {input.long} > {output.filt_long} 2>{log}
        """


# ── Rule: filtered_long_read_qc — the same profile AFTER filtering ───────────
# Takes in: FILT_LONG.
# Produces: NANOPLOT_FILT_DIR.
# Consumed by: multiqc only. Side by side with the raw plot it shows exactly what
#              filtlong removed, which is the honest way to report a filter.
#
# Same --prefix rule as raw_long_read_qc: the report regex depends on it.
#
# (v1 message: "--- NanoPlot: Filtered long-read QC. ---")
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
    priority: 10
    shell:
        """
        NanoPlot \
          --fastq {input.fastq} \
          --threads {threads} \
          --loglength \
          --prefix "{wildcards.sample}_" \
          --outdir {output.nanoplot_filt_dir} > {log} 2>&1
        """
