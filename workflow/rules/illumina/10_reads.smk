# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — illumina front end, read preparation (rules/illumina/10_reads.smk)
#
# The biology: an Illumina run is not just the isolate's DNA. Every Illumina lane
# carries a PhiX spike-in (a small bacteriophage genome added as a sequencing
# control), and every read may still carry sequencing-adapter tails and a
# low-quality 3' end. Both would be assembled into junk contigs, so we remove
# the PhiX reads first and then quality/adapter-trim what is left.
#
# Data flow through this module:
#
#   NCBI ──► download_phix ──► build_phix ──► (6 bowtie2 index files)
#                                                     │
#   {sample}_R1/R2.fastq[.gz] ─────────────────────────┴──► map_phix
#   (raw input, config input.illumina_dir)                     │
#                                                              ▼
#                                        {sample}.1.fastq / {sample}.2.fastq
#                                        (the reads that did NOT map to PhiX)
#                                                              │
#                                                              ▼
#                                                        trim_adapters (fastp)
#                                                              │
#                        ┌─────────────────────────────────────┤
#                        ▼                                     ▼
#                 TRIM_R1 / TRIM_R2                    FASTP_JSON / FASTP_HTML
#                        │                             (MultiQC / human reading)
#     ┌──────────────────┼──────────────────┐
#     ▼                  ▼                  ▼
# illumina_assembly  map_contigs        map_amr_db
# (20_assembly)      (shared/10_decontam) (shared/50_amr, CARD leg)
#
# OWNERSHIP NOTE — read before adding a rule here. `index_contigs` and
# `map_contigs` (trimmed reads back onto the draft assembly, for BlobTools
# coverage) are NOT in this module: they are owned by shared/10_decontam.smk,
# which defines them inside its `if HAS_SHORT_READS:` branch. This front end's
# only obligation to that module is to PRODUCE TRIM_R1 and TRIM_R2.
#
# Every path comes from 00_common.smk (TRIM_R1, TRIM_R2, FASTP_JSON, FASTP_HTML,
# PHIX_LINK, PHIX_FASTA, PHIX_BT2_PREFIX, ILLUMINA_DIR, R1, R2, LOGS, CPUS,
# capped_cpus). Nothing here re-derives a stage number or an output root.
#
# conda: paths resolve relative to THIS file (workflow/rules/illumina/), so
# "../../envs/x.yaml" climbs illumina/ -> rules/ -> workflow/ -> workflow/envs/.
#
# Resource convention: cpu-bound rules request `resources: cpus = ...` and use
# {resources.cpus} in the shell. We do NOT use Snakemake's `threads:` keyword,
# which is why `--resources cpus=N` is mandatory on the command line.
# ─────────────────────────────────────────────────────────────────────────────


# Where this mode's per-sample read files live. DERIVED from TRIM_R1 rather than
# retyped, the same anti-drift trick 00_common uses for BLOB_COV and SPADES_DIR:
# these intermediates must sit beside the trimmed reads, and if TRIM_R1 ever
# moves they follow automatically. Module-local (only the rules below use it), so
# it does not belong in 00_common.
_READS_DIR = os.path.dirname(TRIM_R1)


# ── Rule: download_phix — fetch the PhiX control genome ──────────────────────
# Takes in:  nothing (the URL is config links.phix_link, validated in 00_common).
# Does:      one wget into 01.reads/phix/.
# Produces:  PHIX_FASTA, temp() — it is only needed until the index is built.
# Consumed by: build_phix.
#
# Runs ONCE per run, not once per sample: there is no {sample} in the path.
#
# conda: NONE — deliberately inherited from v1, which used the system `wget`.
# Same deferred decision as cazyme_db_download in shared/40_annotation.smk: if we
# ever give the download rules an env, all of them should get it together.
#
# (v1 message: "--- Download PhiX genome from NCBI. ---")
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
# Takes in:  PHIX_FASTA from download_phix.
# Does:      bowtie2-build, which writes six binary index files sharing a prefix.
# Produces:  the six index files, all temp() — cheaper to rebuild than to keep,
#            and useless once every sample has been screened.
# Consumed by: map_phix.
#
# (v1 message: "--- Bowtie2: Build PhiX genome db. ---")
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
        # bowtie2-build takes the shared PREFIX of those six files, not a filename.
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
# Biology: the spike-in is a real, different genome. Leaving its reads in would
# put a ~5.4 kb phage contig in the assembly and skew the coverage statistics
# BlobTools later uses to separate organisms.
#
# Takes in: the six PhiX index files (the DAG edge to build_phix) and this
#           sample's RAW read pair from config input.illumina_dir.
# Does:     bowtie2 in --local mode against PhiX. The trick is --un-conc, which
#           writes the pairs that did NOT align concordantly to PhiX — i.e. the
#           reads we want to keep. The SAM of the reads that DID align is
#           produced only because bowtie2 must write one somewhere; it is temp().
# Produces: {sample}.1.fastq / {sample}.2.fastq (PhiX-free pairs, temp()).
# Consumed by: trim_adapters.
#
# --un-conc NAMING, preserved verbatim from v1 and load-bearing: given the
# basename "{sample}.fastq", bowtie2 inserts ".1"/".2" BEFORE the final extension
# and writes "{sample}.1.fastq" and "{sample}.2.fastq". Change params.basename and
# the two declared outputs stop existing, with no error from bowtie2 itself.
#
# (v1 message: "--- Bowtie2: Map reads against PhiX genome db. ---")
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
        # See the --un-conc note above: bowtie2 derives the two output names from
        # this one basename.
        basename = _READS_DIR + "/{sample}.fastq",
    conda:
        "../../envs/bowtie.yaml"
    resources:
        cpus = CPUS
    log:
        LOGS + "/map_phix_{sample}.log"
    priority: 10
    shell:
        """
        bowtie2 \
          -x {params.db} \
          -1 {input.r1} -2 {input.r2} \
          --threads {resources.cpus} \
          --un-conc {params.basename} \
          -S {output.sam} \
          --local \
          > {log} 2>&1
        """


# ── Rule: trim_adapters — adapter removal and quality trimming (fastp) ───────
# Biology: adapter read-through and low-quality tails create false k-mers, which
# SPAdes turns into short spurious contigs and mis-assemblies. fastp detects the
# adapter from the pairing itself (--detect_adapter_for_pe), trims a sliding
# window from both ends (--cut_front --cut_right) and drops anything shorter than
# 100 bp, which is the length below which a read stops helping a 127-mer assembly.
#
# Takes in: the PhiX-free pairs from map_phix.
# Produces: TRIM_R1 / TRIM_R2 (temp() — Snakemake keeps them until the LAST
#           consumer is done, i.e. SPAdes, map_contigs and the CARD leg), plus
#           the fastp JSON (MultiQC) and HTML (for a human) reports, both kept.
# Consumed by: illumina_assembly (20_assembly.smk), map_contigs
#              (shared/10_decontam.smk), map_amr_db (shared/50_amr.smk).
#
# v1 -> v2: `resources: cpus` was a hard-coded 16, which asked for 16 threads even
# on an 8-core machine. capped_cpus(16) requests min(CPUS, 16) instead — the same
# ceiling, but never more than the run was given.
#
# (v1 message: "--- Fastp: Remove adapters and quality filter. ---")
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
    resources:
        cpus = capped_cpus(16)
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
          --thread {resources.cpus} \
          --verbose \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          -j {output.json} \
          -h {output.html} > {log} 2>&1
        """
