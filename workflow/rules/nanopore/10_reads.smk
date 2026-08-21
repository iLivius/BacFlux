# BacFlux v2.0.0 — nanopore front end, read preparation.
#
# An ONT run produces reads of wildly varying length and quality, and far more
# data than a bacterial isolate needs. Assembling everything is slower AND worse:
# short, poor reads add errors without adding contiguity. So filtlong keeps the
# longest, best subset, and NanoPlot runs on both sides of it so the choice is
# visible in the report rather than hidden.
#
# Stage chain: raw {sample}_ont → filter_long_reads → FILT_LONG, with a NanoPlot
# profile hanging off either side.
#
# raw_long_read_qc      : NanoPlot over the raw FASTQ → NANOPLOT_RAW_DIR.
# filter_long_reads     : filtlong length/quality selection → FILT_LONG.
# filtered_long_read_qc : the same NanoPlot over FILT_LONG → NANOPLOT_FILT_DIR.
#
# There is NO PhiX step here: PhiX is an Illumina spike-in and does not exist in
# an ONT library.
#
# Everything referenced here comes from 00_common.smk: NANOPORE_DIR, ONT,
# FILT_LONG, NANOPLOT_RAW_DIR, NANOPLOT_FILT_DIR, FILTLONG_MIN_LENGTH,
# FILTLONG_KEEP_PERCENT, LOGS, capped_cpus.
#
# conda: paths resolve relative to THIS file (workflow/rules/nanopore/), so
# "../../envs/x.yaml" climbs nanopore/ → rules/ → workflow/ → workflow/envs/.


# ──────────────────────── Raw read profile (NanoPlot) ──────────
# Takes in: fastq = os.path.join(NANOPORE_DIR, ONT), this sample's raw ONT FASTQ.
#           No rule produces it — NANOPORE_DIR is config input.nanopore_dir and
#           ONT is the {sample}_ont.<ext> filename pattern, both from 00_common.smk.
# Does:     NanoPlot summarises read length, quality and yield for that FASTQ.
#           --loglength adds log-scaled length plots, because ONT read lengths span
#           orders of magnitude and collapse into one spike on a linear axis.
# Produces: nanoplot_raw_dir = NANOPLOT_RAW_DIR, 01.reads/{sample}/ont/raw_qc — a
#           DIRECTORY of plots plus {sample}_NanoStats.txt.
# Consumed by: multiqc (shared/90_report.smk), and nothing else.
#
# --prefix "{sample}_" IS LOAD-BEARING, not cosmetic. MultiQC names a sample after
# the file it read, and shared/90_report.smk rewrites those names with the regex
#   '^01\.reads \| ([^|]+) \| ont \| raw_qc \| \1$'
# whose backreference requires the name INSIDE the directory to equal the sample
# name — which only happens when NanoPlot writes {sample}_NanoStats.txt. Drop the
# prefix and NanoPlot silently disappears from the report.
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
    shell:
        """
        # NanoPlot can write its stats file before --outdir exists under load.
        mkdir -p {output.nanoplot_raw_dir}

        NanoPlot \
          --fastq {input.fastq} \
          --threads {threads} \
          --loglength \
          --prefix "{wildcards.sample}_" \
          --outdir {output.nanoplot_raw_dir} > {log} 2>&1
        """


# ──────────────────────── Read selection (filtlong) ────────────
# Takes in: long = os.path.join(NANOPORE_DIR, ONT), the same raw ONT FASTQ
#           raw_long_read_qc profiled above — filtlong reads the reads, not the
#           plots, so the two rules are independent of each other.
# Does:     filtlong scores every raw ONT read and keeps the longest, best-scoring
#           subset. Three knobs, and only one of them is fixed in this rule:
#             --min_length     parameters.long_read_qc.min_length, default 1000.
#                              Reads under ~1 kb bring ONT's error rate without
#                              ONT's main benefit, which is spanning repeats.
#             --keep_percent   parameters.long_read_qc.keep_percent, default 95 —
#                              keep the best N% of the remaining BASES. It was 90;
#                              the long_read_qc block in config.yaml explains why
#                              being less selective is safer for small plasmids.
#             --target_bases   hard-coded 5e8 here: stop at 500 Mbp, roughly 100x
#                              for a 5 Mbp bacterial genome. More coverage than
#                              that costs Flye time without improving the assembly.
# Produces: filt_long = FILT_LONG, 01.reads/{sample}/ont/{sample}_filt.fastq — the
#           kept reads, and the only ONT reads anything downstream ever sees.
# Consumed by: filtered_long_read_qc (below), check_medaka_model
#              (shared/12_medaka_check.smk, which infers the Medaka model from
#              these reads), ont_assembly (nanopore/20_assembly.smk),
#              long_read_consensus (nanopore/30_polish.smk) and map_contigs
#              (shared/10_decontam.smk, the ONT coverage track for BlobTools).
#
# GOTCHA: setting parameters.long_read_qc.length_weight has NO effect in this
# mode. There is no --length_weight below, so filtlong uses its own default of 1.
# That key reaches filtlong in hybrid mode only.
#
# keep_percent, above, DOES apply here, and it is the key to reach for if a small
# plasmid goes missing. On K. pneumoniae TUM24772 the reads covering a 5,596 bp
# Col2 plasmid went from 93 of 603 at the old keep_percent 90 to 500 at the 95
# that now ships, without touching length_weight. Note that recovering the reads
# did NOT recover the plasmid: nothing we tried assembled it. The full 2x2 of
# read counts is in config.yaml, next to the long_read_qc block.
#
# Hybrid's whole flag set differs (--trim --split --length_weight, and no
# --target_bases) because there the decontaminated Illumina reads, not a coverage
# cap, define what is worth keeping — see hybrid/30_ont_reads.smk.
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
    shell:
        """
        filtlong \
          --min_length {params.min_length} \
          --keep_percent {params.keep_percent} \
          --target_bases {params.target_bases} \
          {input.long} > {output.filt_long} 2>{log}
        """


# ──────────────────── Filtered read profile (NanoPlot) ─────────
# Takes in: fastq = FILT_LONG, the kept reads from rule filter_long_reads above.
# Does:     the same NanoPlot profile again, same flags, this time over the subset
#           that survived the filter.
# Produces: nanoplot_filt_dir = NANOPLOT_FILT_DIR, 01.reads/{sample}/ont/filt_qc —
#           the same directory of plots plus {sample}_NanoStats.txt.
# Consumed by: multiqc (shared/90_report.smk) only. Side by side with the raw plot
#              it shows exactly what filtlong removed, which is the honest way to
#              report a filter.
#
# Same --prefix rule as raw_long_read_qc — the report regex depends on it, here on
# its filt_qc branch.
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
    shell:
        """
        # NanoPlot can write its stats file before --outdir exists under load.
        mkdir -p {output.nanoplot_filt_dir}

        NanoPlot \
          --fastq {input.fastq} \
          --threads {threads} \
          --loglength \
          --prefix "{wildcards.sample}_" \
          --outdir {output.nanoplot_filt_dir} > {log} 2>&1
        """
