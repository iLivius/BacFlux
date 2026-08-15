# Hybrid front end, short-read side: PhiX removal, adapter trimming, and a first
# look at the raw ONT reads. Hybrid sequences one isolate twice — Illumina for
# base accuracy, ONT for contiguity — and the two legs only meet in
# hybrid/30_ont_reads.smk, where the cleaned Illumina reads become the filter that
# scores the ONT reads.
#
# Illumina chain: raw {sample}_R1/R2.fastq[.gz] (config input.illumina_dir) →
# map_phix → trim_adapters → TRIM_R1/TRIM_R2, read by illumina_assembly
# (hybrid/20_assembly.smk), map_contigs (shared/10_decontam.smk), map_sel_contigs
# (hybrid/30_ont_reads.smk) and map_amr_db (shared/50_amr.smk). The fastp JSON
# goes to multiqc (shared/90_report.smk); the HTML twin is for a human only.
# ONT chain: raw {sample}_ont.fastq[.gz] (config input.nanopore_dir) →
# raw_long_read_qc → NANOPLOT_RAW_DIR → multiqc.
#
# The four Illumina rules are identical to rules/illumina/10_reads.smk and
# raw_long_read_qc is identical to rules/nanopore/10_reads.smk — keep them in
# sync. They are duplicated rather than shared because only one mode's rule
# directory is ever included, so a duplicate rule NAME can never collide, and
# because someone working in hybrid mode should find the whole front end in one
# place. Both copies read the same constants from shared/00_common.smk, so drift
# is bounded to tool flags and can never reach a path.
#
# index_contigs and map_contigs — the trimmed reads mapped onto the ILLUMINA
# draft, for the BlobTools coverage track — belong to shared/10_decontam.smk, not
# here. The index_selected_contigs / map_sel_contigs pair in
# hybrid/30_ont_reads.smk is a DIFFERENT job: reads mapped onto the
# DECONTAMINATED contigs, to pull the clean read set back out. That pair is
# hybrid-only.
#
# download_phix    : fetch the PhiX control genome, once per run.
# build_phix       : bowtie2-build an index of it.
# map_phix         : map each raw pair against PhiX and keep what does NOT align.
# trim_adapters    : fastp adapter removal and quality trimming.
# raw_long_read_qc : NanoPlot length/quality profile of the ONT reads BEFORE
#                    filtering. The matching "after filtering" run lives in
#                    hybrid/30_ont_reads.smk, because in hybrid mode filtering
#                    cannot start until the Illumina leg has produced its
#                    decontaminated reads.
#
# Every constant used here comes from shared/00_common.smk.


# Where this mode's per-sample Illumina working files live. Derived from TRIM_R1
# rather than retyped, so a change to the stage layout cannot leave this behind —
# the same anti-drift trick as SPADES_DIR in shared/00_common.smk.
_READS_DIR = os.path.dirname(TRIM_R1)


# ──────────────────────── PhiX removal (Bowtie2) ───────────────
# PhiX is the spike-in control present in essentially every Illumina lane. Its
# reads are real sequence from a different genome, so they have to go before
# assembly or they become contigs.
#
# Takes in: nothing on disk — params.link = PHIX_LINK, the download URL resolved
#           once from config links.phix_link in shared/00_common.smk.
# Does:     wget the PhiX genome. Runs once per RUN, not once per sample: there
#           is no {sample} anywhere in the path.
# Produces: phix = PHIX_FASTA, the gzipped PhiX genome. temp(), so it goes as
#           soon as the index below has been built from it.
# Consumed by: build_phix.
#
# Identical to rules/illumina/10_reads.smk — keep in sync.
# No conda environment: system wget, kept for v1 parity (the same deferred
# decision as cazyme_db_download in shared/40_annotation.smk).
rule download_phix:
    output:
        phix = temp(PHIX_FASTA),
    params:
        link = PHIX_LINK,
    log:
        LOGS + "/download_phix.log"
    shell:
        """
        wget {params.link} -O {output.phix} > {log} 2>&1
        """


# ──────────────────────── PhiX index (Bowtie2) ─────────────────
# Takes in: phix = PHIX_FASTA, from rule download_phix.
# Does:     bowtie2-build, writing under the prefix PHIX_BT2_PREFIX.
# Produces: six temp() index files — .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2
#           and .rev.2.bt2 — declared individually so Snakemake tracks all six.
# Consumed by: map_phix, and nothing else.
#
# Identical to rules/illumina/10_reads.smk — keep in sync.
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
    shell:
        """
        bowtie2-build \
          {input.phix} \
          {params.basename} > {log} 2>&1
        """


# ──────────────────────── PhiX screening (Bowtie2) ─────────────
# Takes in: the six PhiX index files from build_phix — taking them as inputs is
#           what creates the DAG edge — plus this sample's RAW Illumina pair, R1
#           and R2 under ILLUMINA_DIR, i.e. {sample}_R1/R2.fastq[.gz] in config
#           input.illumina_dir.
# Does:     bowtie2 --local against the PhiX index, keeping the pairs that did
#           NOT align concordantly. Those are the reads that are not PhiX.
# Produces:
#   sam = {sample}_contam.sam, the alignments of the reads that DID match PhiX.
#         temp(), and thrown away.
#   r1  = {sample}.1.fastq, and
#   r2  = {sample}.2.fastq — the PhiX-free pair, both temp(), written by
#         bowtie2's --un-conc.
# Consumed by: trim_adapters.
#
# Identical to rules/illumina/10_reads.smk — keep in sync.
#
# bowtie2's `--un-conc` does not write the filename it is given. From the
# basename "{sample}.fastq" it builds "{sample}.1.fastq" and "{sample}.2.fastq"
# by inserting .1 / .2 before the final extension. So params.basename is spelled
# to make those two names come out equal to the declared outputs — edit it and
# the outputs Snakemake is waiting for are simply never created.
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


# ──────────────────────── Adapter trimming (fastp) ─────────────
# Takes in: r1 / r2 = the PhiX-free pairs {sample}.1.fastq / {sample}.2.fastq,
#           from rule map_phix.
# Does:     fastp adapter removal and quality trimming.
# Produces:
#   r1, r2 = TRIM_R1 / TRIM_R2, the trimmed pair. temp(), so they survive only
#            until the last consumer listed below has finished with them.
#   json   = FASTP_JSON, the machine-readable fastp report
#   html   = FASTP_HTML, its human-readable twin
# Consumed by: the trimmed pair goes to illumina_assembly
#              (hybrid/20_assembly.smk), map_contigs (shared/10_decontam.smk),
#              map_sel_contigs (hybrid/30_ont_reads.smk), map_amr_db
#              (shared/50_amr.smk) and, when mobilome.run is true, assembly_depth
#              and isosdb_map (shared/80_mobilome.smk); FASTP_JSON goes to
#              multiqc (shared/90_report.smk). Nothing reads FASTP_HTML — it is
#              there for a human.
#
# Identical to rules/illumina/10_reads.smk — keep in sync.
#
# v1 → v2 on the thread count: hybrid v1 asked for min(CPUS, 24) here and BacFlux
# v1 for a hard-coded 16. Both short-read front ends now use capped_cpus(16),
# because fastp is I/O-bound past roughly that point and one value across the two
# modes is one thing less to reconcile. The hybrid env pin also moved from fastp
# 1.1.0 to the tree's existing fastp 1.0.1, the only env divergence the migration
# turned up.
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


# ──────────────────────── Raw ONT read profile (NanoPlot) ──────
# Takes in: fastq = this sample's raw ONT FASTQ, ONT under NANOPORE_DIR, i.e.
#           {sample}_ont.fastq[.gz] in config input.nanopore_dir. Nothing has
#           been filtered out of it yet.
# Does:     NanoPlot, with --loglength so the read-length distribution is plotted
#           on a log axis — ONT lengths span orders of magnitude.
# Produces: nanoplot_raw_dir = NANOPLOT_RAW_DIR, a directory of plots plus
#           {sample}_NanoStats.txt.
# Consumed by: multiqc (shared/90_report.smk).
#
# Identical to rules/nanopore/10_reads.smk — keep in sync.
#
# `--prefix "{sample}_"` is load-bearing. Without it NanoPlot writes a bare
# NanoStats.txt, and shared/90_report.smk relabels this directory with a regex
# whose backreference needs the MultiQC sample name to equal {sample} — which
# only holds when the stats file is called {sample}_NanoStats.txt. Drop the
# prefix and NanoPlot disappears from the report silently, with no error
# anywhere.
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
    shell:
        """
        NanoPlot \
          --fastq {input.fastq} \
          --threads {threads} \
          --loglength \
          --prefix "{wildcards.sample}_" \
          --outdir {output.nanoplot_raw_dir} > {log} 2>&1
        """
