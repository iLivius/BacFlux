# BacFlux v2.0.0 — illumina front end, read preparation.
#
# An Illumina run is not just the isolate's DNA. Every lane carries a PhiX
# spike-in (a small bacteriophage genome added as a sequencing control), and every
# read may still carry an adapter tail and a low-quality 3' end. Both would be
# assembled into junk contigs, so the PhiX reads go first and what survives is
# trimmed.
#
# Stage chain: raw {sample}_R1/R2 → map_phix → trim_adapters → TRIM_R1/TRIM_R2.
# The PhiX reference is fetched and indexed ONCE per run, not per sample:
# download_phix → build_phix → six Bowtie2 index files shared by every sample.
#
# download_phix  : wget the PhiX genome named by links.phix_link.
# build_phix     : bowtie2-build over it — six index files sharing one prefix.
# map_phix       : bowtie2 --local per sample, keeping the pairs that do NOT
#                  align to PhiX.
# trim_adapters  : fastp adapter detection, sliding-window quality trimming and a
#                  100 bp minimum length. Writes TRIM_R1/TRIM_R2 plus the fastp
#                  JSON (read by MultiQC) and HTML (read by a human) reports.
#
# OWNERSHIP — read before adding a rule here. index_contigs and map_contigs (the
# trimmed reads mapped back onto the draft assembly, giving BlobTools its coverage
# track) are NOT in this module: they live in the `if HAS_SHORT_READS:` branch of
# shared/10_decontam.smk. This front end's only obligation to that module is to
# PRODUCE TRIM_R1 and TRIM_R2.
#
# Every path comes from 00_common.smk (TRIM_R1, TRIM_R2, FASTP_JSON, FASTP_HTML,
# PHIX_LINK, PHIX_FASTA, PHIX_BT2_PREFIX, ILLUMINA_DIR, R1, R2, LOGS, CPUS,
# capped_cpus). Nothing here re-derives a stage number or an output root.
#
# conda: paths resolve relative to THIS file (workflow/rules/illumina/), so
# "../../envs/x.yaml" climbs illumina/ → rules/ → workflow/ → workflow/envs/.
#
# Resource convention: cpu-bound rules declare Snakemake's BUILT-IN `threads:` —
# CPUS, or capped_cpus(N) where a tool stops scaling past N — and refer to
# {threads} in the shell. Using the built-in keyword rather than a `resources:
# cpus` of our own is what makes `--cores N` actually enforce the limit, so a
# plain `snakemake --cores N` is safe on its own, with no extra --resources flag.


# Where this mode's per-sample read files live. DERIVED from TRIM_R1 rather than
# retyped: these intermediates must sit beside the trimmed reads, and if TRIM_R1
# ever moves they follow automatically. The same anti-drift trick 00_common.smk
# uses for BLOB_COV. Module-local — only the rules below use it — so it does not
# belong in 00_common.
_READS_DIR = os.path.dirname(TRIM_R1)


# ──────────────────────── PhiX reference (NCBI) ────────────────
# Takes in: nothing on disk. The address is PHIX_LINK, which 00_common.smk reads
#           from the config key links.phix_link and validates at parse time.
# Does:     one wget of the PhiX control genome into 01.reads/phix/.
# Produces: phix = PHIX_FASTA (01.reads/phix/phix.fna.gz), temp() — the FASTA is
#           wanted only until build_phix has indexed it.
# Consumed by: build_phix, below.
#
# Runs once per RUN, not once per sample — there is no {sample} in the path, so
# every sample in the batch screens against this one download.
#
# conda: NONE — deliberately inherited from v1, which used the system `wget`. Same
# deferred decision as cazyme_db_download in shared/40_annotation.smk: if the
# download rules ever get an env, they should all get one together.
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


# ──────────────────────── PhiX index (Bowtie2) ─────────────────
# Takes in: phix = PHIX_FASTA, the genome fetched by rule download_phix.
# Does:     bowtie2-build over it, once per run. The tool is handed the shared
#           PREFIX the index files will have (params.basename = PHIX_BT2_PREFIX),
#           not a filename, and works out the six names itself.
# Produces:
#   idx1, idx2, idx3, idx4 = PHIX_BT2_PREFIX + ".1.bt2" through ".4.bt2", the
#                            forward index plus bowtie2's own packed copy of the
#                            reference sequence, all temp()
#   ridx1, ridx2           = PHIX_BT2_PREFIX + ".rev.1.bt2" and ".rev.2.bt2", the
#                            mirror index of the same genome, also temp()
# Consumed by: map_phix, below.
#
# All six are temp() because they are cheaper to rebuild than to keep, and useless
# once every sample has been screened.
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


# ──────────────────────── PhiX screening (Bowtie2) ─────────────
# Throws away the reads that are PhiX. The spike-in is a real, different genome:
# leaving its reads in puts a ~5.4 kb phage contig in the assembly and skews the
# coverage statistics BlobTools later uses to separate organisms.
#
# Takes in:
#   idx1–idx4,   = the six Bowtie2 index files from rule build_phix. They are
#   ridx1, ridx2   listed only to create the DAG edge to that rule; bowtie2 itself
#                  is handed their shared prefix in params.db.
#   r1, r2       = this sample's RAW read pair, os.path.join(ILLUMINA_DIR, R1) and
#                  os.path.join(ILLUMINA_DIR, R2) — the files found in the config
#                  key input.illumina_dir, untouched by any earlier rule.
# Does:     bowtie2 in --local mode against the PhiX index. --un-conc is what does
#           the work: it writes out the pairs that did NOT align concordantly to
#           PhiX, i.e. the reads we keep.
# Produces:
#   sam = {sample}_contam.sam, the alignments of the reads that DID match PhiX.
#         temp(), and it exists only because bowtie2 must write a SAM somewhere.
#   r1  = {sample}.1.fastq, temp() — the PhiX-free forward reads
#   r2  = {sample}.2.fastq, temp() — the PhiX-free reverse reads
# Consumed by: trim_adapters, below, which takes r1 and r2. Nothing reads the SAM.
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


# ──────────────────── Adapter and quality trim (fastp) ─────────
# Adapter read-through and low-quality tails create false k-mers, which SPAdes
# turns into short spurious contigs and mis-assemblies. This is where both are cut
# off, before the reads ever reach the assembler.
#
# Takes in: r1, r2 = {sample}.1.fastq and {sample}.2.fastq, the PhiX-free pairs
#           from rule map_phix.
# Does:     fastp adapter removal and quality trimming. --detect_adapter_for_pe
#           infers the adapter from the pairing itself, --cut_front and --cut_right
#           trim a sliding quality window from each end, and --length_required 100
#           drops whatever is left of a read shorter than 100 bp — the length below
#           which a read stops helping a 127-mer assembly. The 100 is the v1 value,
#           hard-coded here rather than exposed in the config.
# Produces:
#   r1   = TRIM_R1, temp() — the trimmed, PhiX-free forward reads
#   r2   = TRIM_R2, temp() — the trimmed, PhiX-free reverse reads
#   html = FASTP_HTML, fastp's own report, kept for a human to open
#   json = FASTP_JSON, the same numbers in parsable form, kept for MultiQC
# Consumed by: TRIM_R1/TRIM_R2 by illumina_assembly (illumina/20_assembly.smk),
#              map_contigs (shared/10_decontam.smk), map_amr_db — the CARD
#              read-mapping leg — (shared/50_amr.smk) and, when the IS copy-number
#              leg is on (MOBILOME_COPY_NUMBER — mobilome.run AND a configured
#              ISOSDB source), assembly_depth and isosdb_map (80_mobilome.smk).
#              FASTP_JSON by MultiQC (shared/90_report.smk). FASTP_HTML by the
#              user only: no rule reads it.
#
# Because TRIM_R1/TRIM_R2 are temp(), Snakemake keeps them on disk until the LAST
# of those consumers is done and then deletes them. The two fastp reports are kept.
#
# v1 → v2: `resources: cpus` was a hard-coded 16, which asked for 16 threads even
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
