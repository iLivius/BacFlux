# Stage 20 assembly QC: four questions about the finished genome, one tool each.
#   QUAST    — is the ASSEMBLY any good? Contig count, N50, total length, GC.
#   CheckM   — is the GENOME complete, and is it one organism? Single-copy marker
#              genes, reported as completeness % and contamination %.
#   Qualimap — did the READS map back sensibly? Depth, evenness, insert size.
#   BLAST    — does the GENOME still carry sequencing adapters or barcodes? The
#              adapter screen, at the end of this file.
# QUAST, CheckM and the adapter screen need a genome; Qualimap needs the BAM
# shared/10_decontam.smk already made for BlobTools.
#
# stage_qc_genomes               : copies this sample's 1-or-2 assemblies into
#                                  one temporary directory under controlled
#                                  names, and writes the table saying which is
#                                  which.
# genome_assembly_evaluation     : QUAST over that directory.
# completeness_and_contamination : CheckM lineage_wf then qa, same directory.
# map_evaluation                 : Qualimap on the decontamination BAM. Only in
#                                  modes that have reads — see its banner below.
# univec_illumina                : once per run, downloads NCBI's UniVec and keeps
#                                  its Illumina adapter records.
# adapter_search_input           : the files the adapter search reads: every
#                                  adapter in one FASTA, and the genome with its
#                                  circular contigs extended across their start.
# blast_adapters                 : BLAST of those adapters against the genome.
# screen_adapters                : the hits turned into adapter loci, with a
#                                  reason for each, and the sample's row of the
#                                  Adapters table of the MultiQC report.
#
# GTDB-Tk reads the same staged directory from shared/30_taxonomy.smk, and all
# five results end up in MultiQC (shared/90_report.smk).
#
# The hybrid complication, and why stage_qc_genomes exists. Three modes deliver
# one genome per sample; hybrid delivers the ONT+Polypolish genome AND keeps the
# decontaminated Illumina draft, and v1 BacFluxL+ deliberately QC'd and classified
# both so the two could be compared. Rather than duplicate rule bodies or invent a
# second wildcard, one small rule stages the 1-or-2 assemblies and QUAST, CheckM
# and GTDB-Tk all read that single directory. Names, roles and labels come from
# the QC_GENOMES list in 00_common.smk — one source for the staged filenames, the
# QUAST labels, the CheckM/GTDB-Tk bin ids and the MultiQC relabelling.
#
# Two rejected alternatives, recorded so they are not re-litigated:
#   * fan out over a {genome} wildcard — this would run GTDB-Tk twice per hybrid
#     sample. classify_wf loads the reference tree and skani sketches on every
#     invocation, so that doubles the single most expensive step in the pipeline;
#     CheckM's lineage_wf likewise shares its HMM/pplacer setup across bins in one
#     run. It would also add a wildcard and a directory level in the three
#     single-genome modes for no gain.
#   * `if IS_HYBRID:` with two near-duplicate rule bodies — duplicated shell is
#     what this unification exists to delete, and duplicates drift apart.
#
# Everything referenced here is defined once in 00_common.smk and never re-derived:
# qc_genome_fastas, qc_stage_commands, qc_genome_table_text, QC_GENOMES_DIR,
# QC_GENOME_TABLE, QUAST_DIR, CHECKM_DIR, CHECKM_STATS, CHECKM_LINEAGE,
# QUALIMAP_DIR, DECONTAM_BAM, LOGS, RAM, capped_cpus, HAS_READS; for the adapter
# screen FINAL_CONTIGS, BAKTA_REPLICON_INPUT, UNIVEC_LINK, UNIVEC_SHA256,
# UNIVEC_ILLUMINA, UNIVEC_PROVENANCE, ADAPTER_LIST_FILES, ADAPTER_SCREEN_TRIMMED,
# ADAPTER_SCREEN_CHOICE, ADAPTER_SCREEN_CHECK, ADAPTER_QUERY, ADAPTER_SUBJECT,
# ADAPTER_BLAST, ADAPTER_SCREEN, ADAPTER_SCREEN_MQC and ADAPTER_SCREEN_SCRIPT.
#
# conda: env paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/x.yaml" climbs shared/ → rules/ → workflow/ → workflow/envs/x.yaml.


# ──────────────────────── Staging the genomes to QC ────────────
# Not a biology step: bookkeeping that keeps the three tools below simple.
#
# Takes in: an input FUNCTION, because the NUMBER of files varies by mode —
#           qc_genome_fastas() returns FINAL_CONTIGS in illumina/nanopore/contigs,
#           and [DECONTAM_CONTIGS (the Illumina draft), FINAL_CONTIGS (the ONT
#           genome)] in hybrid. A templated string cannot express "1 or 2 files".
# Does:     copies each into eval/genomes/ under its bin id — {sample}.fasta, or
#           {sample}_illumina.fasta + {sample}_ont.fasta in hybrid. Those
#           basenames become the ids CheckM, GTDB-Tk and QUAST report under, which
#           is why they come from the shared QC_GENOMES list rather than from here.
# Produces:
#   genomes_dir  = a temp() DIRECTORY, and temp() is the whole point: Snakemake
#                  owns the cleanup and removes it once the last consumer (QUAST,
#                  CheckM, GTDB-Tk) is done. v1 had no such rule — CheckM copied
#                  the genomes into its own output directory and the taxonomy rule
#                  ended with `rm -rf {input.checkm_dir}/*.fasta`, one rule
#                  deleting files inside another rule's declared output. Re-running
#                  taxonomy on its own then found an empty --genome_dir and GTDB-Tk
#                  classified nothing, without a useful error. No rule deletes
#                  anything now.
#   genome_table = a kept TSV — bin_id, role, technology, source_path. It says in
#                  writing which of the two hybrid rows is the DELIVERED genome
#                  ("primary") and which is the comparator, so a reader never has
#                  to infer it from a filename suffix. No rule reads it; it is a
#                  terminal output that all_targets() in 00_common.smk asks for so
#                  that it is always written, and it follows the project's
#                  audit-file convention.
# Consumed by: genome_assembly_evaluation and completeness_and_contamination in
#              this file, taxonomic_assignment in shared/30_taxonomy.smk.
#
# Runtime to expect (a v1 hazard made visible, not fixed): a hybrid sample does
# twice the CheckM and GTDB-Tk work inside ONE job, under one thread reservation.
# That was already true in v1; it is written down here and in genome_table so the
# wall-clock is not a surprise.
#
# conda: NONE — mkdir and cp only, from the launch environment.
rule stage_qc_genomes:
    input:
        genomes = qc_genome_fastas,
    output:
        genomes_dir = temp(directory(QC_GENOMES_DIR)),
        genome_table = QC_GENOME_TABLE,
    params:
        # Both are built from QC_GENOMES while Snakemake works out the job, before
        # anything runs, so a dry run and the log show the literal `cp` lines and
        # the literal table content instead of a loop over hidden bash arrays.
        copy_cmds = qc_stage_commands,
        table_text = qc_genome_table_text,
    log:
        LOGS + "/stage_qc_genomes_{sample}.log"
    shell:
        # {params.copy_cmds} expands to one `cp` line per genome. Snakemake formats
        # the shell TEMPLATE only and never re-scans a substituted value for {}
        # placeholders, so generated multi-line commands are safe here. The heredoc
        # is quoted ('EOF'), so the table text lands in the file exactly as built —
        # no shell expansion, no printf format-string hazard.
        """
        mkdir -p {output.genomes_dir}

        {params.copy_cmds}

        cat > {output.genome_table} << 'EOF'
{params.table_text}
EOF

        ls -l {output.genomes_dir} > {log} 2>&1
        """


# ──────────────────────── Assembly metrics (QUAST) ─────────────
# QUAST reports the structural quality of the assembly — how many contigs it took
# to represent the genome, N50, largest contig, total length, GC. It says nothing
# about biological completeness, which is CheckM's job below; it says how
# fragmented the assembly is, and fragmentation limits everything downstream that
# depends on gene context.
#
# Takes in: the staged genomes DIRECTORY. The shell globs *.fasta inside it, so a
#           hybrid sample gets both assemblies side by side in ONE report —
#           exactly the comparison a hybrid run is for.
# Produces: 02.assembly/{sample}/eval/quast/ (a directory; QUAST names the files).
# Consumed by: MultiQC (shared/90_report.smk).
#
# Two v1→v2 additions, both purely additive — no existing output changes:
#   1. nanopore had NO QUAST rule at all in v1. It runs here, on the quast.yaml env
#      that already exists in the v2 tree, so no new dependency.
#   2. hybrid QUAST'd only the Illumina draft in v1 and never evaluated the genome
#      it actually delivered. Both are evaluated now.
# One consequence to know about: QUAST labels an assembly by the input file's
# basename, so the label changed from v1's `contigs_sel` to `{sample}` (or
# {sample}_illumina / {sample}_ont). The MultiQC rename patterns are generated from
# the same QC_GENOMES list and already match — see multiqc_replace_block() in
# shared/90_report.smk.
#
# (v1 message: "--- QUAST: Genome assembly evaluation. ---")
rule genome_assembly_evaluation:
    input:
        genomes_dir = QC_GENOMES_DIR,
    output:
        quast_dir = directory(QUAST_DIR),
    conda:
        "../../envs/quast.yaml"
    threads: capped_cpus(24)
    log:
        LOGS + "/assembly_evaluation_{sample}.log"
    shell:
        # The shell expands the *.fasta glob, not Snakemake; the files are
        # guaranteed to be there because genomes_dir is an input of this rule.
        # QUAST is a report: nothing downstream reads its output and MultiQC only
        # searches the directory, so a QUAST crash must not cost the batch. The output
        # directory is created first, the tool is allowed to fail, and a note goes to
        # the log - the same three-part guard rule isescan uses.
        """
        mkdir -p {output.quast_dir}
        quast \
          {input.genomes_dir}/*.fasta \
          -o {output.quast_dir} \
          --no-icarus \
          -t {threads} > {log} 2>&1 || {{
            echo "" >> {log}
            echo "NOTE: QUAST exited non-zero. This is a QC report, not an input to any other step, so the run continues; the assembly-metrics table will be missing from MultiQC for this sample. Check the log above." >> {log}
          }}
        """


# ────────────── Completeness and contamination (CheckM) ────────
# CheckM places the genome in its reference tree, picks the set of single-copy
# marker genes for that lineage, and counts them. Markers missing means an
# incomplete assembly; markers present in more copies than expected means either
# contamination — more than one organism in the bin — or a genuine duplication.
# This is the standard pass/fail gate for an isolate genome, and it is why
# decontamination (shared/10_decontam.smk) has to run first: leave one contaminant
# contig in and it shows up here as contamination of the isolate.
#
# Takes in: the staged genomes DIRECTORY (1 genome, or 2 in hybrid — CheckM takes
#           several "bins" in one run and reports a row for each).
# Does:     lineage_wf to place the genome and count markers, then qa -o 2
#           --tab_table to write the extended per-bin statistics as a parsable TSV.
# Produces:
#   checkm_dir     = 02.assembly/{sample}/eval/checkm/, CheckM's own working and
#                    result tree
#   checkm_stats   = {sample}_checkm_stats.tsv, the table MultiQC reads
#   checkm_lineage = lineage.ms, the marker-set file lineage_wf writes and qa reads
#                    back. The name is CheckM's, not ours.
# Consumed by: MultiQC (shared/90_report.smk).
#
# v1→v2 fix: v1 ran `checkm lineage_wf -x fasta <dir> <dir>` with the SAME
# directory as input and output, which is the whole reason the staged FASTAs had
# to be deleted afterwards by another rule. Input and output are separate now
# ({input.genomes_dir} → {output.checkm_dir}), and that is what lets Snakemake own
# the staged directory as temp().
#
# The stats file gained a {sample}_ prefix (v1 wrote a bare checkm_stats.tsv).
# Under the v2 layout every sample's containing directory is called "checkm", so
# the report's staging loop can no longer recover the sample name from the
# directory and takes it from the filename instead — see the checkm_stats loop in
# shared/90_report.smk.
#
# (v1 message: "--- CheckM: Assessment of genome completeness and contamination. ---";
#  FastaFlux spelled it "completenness" in its message and log name — typo fixed.)
rule completeness_and_contamination:
    input:
        genomes_dir = QC_GENOMES_DIR,
    output:
        checkm_dir = directory(CHECKM_DIR),
        checkm_stats = CHECKM_STATS,
        checkm_lineage = CHECKM_LINEAGE,
    conda:
        "../../envs/checkm.yaml"
    threads: capped_cpus(24)
    resources:
        # Cluster requests (00_common.smk §3a). A single computer reads mem_mb
        # only when the launch line passes --resources mem_mb=N.
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    log:
        LOGS + "/completeness_and_contamination_{sample}.log"
    shell:
        # `checkm qa` takes the lineage_wf OUTPUT directory as its second
        # positional argument — the one holding bins/ and storage/ — and not the
        # directory of input FASTAs. {output.checkm_dir} is right on that line.
        """
        mkdir -p {output.checkm_dir}

        checkm lineage_wf \
          -t {threads} \
          -x fasta {input.genomes_dir} \
          {output.checkm_dir} > {log} 2>&1

        checkm qa \
          -o 2 \
          -t {threads} \
          --tab_table \
          -f {output.checkm_stats} {output.checkm_lineage} {output.checkm_dir} >> {log} 2>&1
        """


# ──────────────────────── Mapping QC (reads-only modes) ────────
# Gated on HAS_READS — the capability flag rather than the mode name — so
# in contigs mode the rule below is never defined and QUALIMAP_DIR is never
# requested by all_targets() in 00_common.smk. Contigs mode does have a BAM, but
# it is the contigs aligned against THEMSELVES, made only to satisfy
# `blobtools create -b`; a mapping-quality report on a self-alignment says
# nothing. v1 FastaFlux had no Qualimap rule and v2 keeps that absence — a
# mode-specific absence, not a gap.
if HAS_READS:

    # ── Alignment quality report (Qualimap) ──
    # Summarises how this sample's reads sit on its own assembly — mean depth and
    # how even it is, GC bias, mapping rate, and insert size for paired reads.
    # Uneven or unexpectedly low coverage is an early warning of a mixed culture,
    # a mis-assembly, or simply not enough data.
    #
    # Takes in: DECONTAM_BAM, the alignment shared/10_decontam.smk made for
    #           BlobTools (bowtie2 in the short-read modes, minimap2 map-ont in
    #           nanopore). That BAM is temp(), so taking a DAG edge on it here is
    #           what keeps it alive until this rule has run.
    # Produces: 02.assembly/{sample}/eval/qualimap/, an HTML report directory.
    # Consumed by: MultiQC (shared/90_report.smk). In hybrid the panel is labelled
    #           "mapping Illumina QC", because the BAM is Illumina reads on the
    #           pre-decontamination SPAdes draft and not coverage of the delivered
    #           ONT genome — see multiqc_replace_block() in shared/90_report.smk.
    #
    # v1→v2 unification: illumina passed neither --java-mem-size nor -nt, so it ran
    # single-threaded on Qualimap's built-in 1,200 MB heap and ran out of memory on
    # large BAMs; nanopore hard-coded 24G; hybrid passed a config-aware heap plus
    # -nt. v2 uses the hybrid form everywhere. The RESULTS are unchanged — this
    # only makes the rule faster and stops it dying on heap.
    #
    # java_mem is the most memory Java may use (its heap), in GB: up to 64, and
    # always 4 GB below ram_gb so that Java's own overhead fits too
    # (QUALIMAP_HEAP_GB in 00_common.smk). mem_mb is what a cluster job asks for:
    # the heap plus those 4 GB. On a single computer mem_mb is ignored, and
    # Snakemake schedules on java_mem only if the launch line passes
    # `--resources java_mem=N`, which then also lowers each heap to N.
    #
    # (v1 message: "--- Qualimap: Mapping evaluation. ---"; the nanopore v1 rule
    #  was named map_qc — v2 uses map_evaluation in every mode that has reads.)
    rule map_evaluation:
        input:
            bam = DECONTAM_BAM,
        output:
            qualimap_dir = directory(QUALIMAP_DIR),
        conda:
            "../../envs/qualimap.yaml"
        threads: capped_cpus(24)
        resources:
            java_mem = QUALIMAP_HEAP_GB,
            # Cluster requests (00_common.smk §3a). A single computer reads mem_mb
            # only when the launch line passes --resources mem_mb=N.
            mem_mb = QUALIMAP_JOB_MB,
            runtime = LONG_RUNTIME_MINUTES,
        log:
            LOGS + "/map_evaluation_{sample}.log"
        shell:
            # Qualimap reports on the read mapping; nothing downstream reads it. Same guard
            # as QUAST: make the output directory, tolerate the tool, note the log.
            """
            # Qualimap draws its charts with Java's graphics library. If DISPLAY is
            # set -- an SSH session with X11 forwarding, say -- Java tries to open
            # that screen and dies with "AWTError: Can't connect to X11 window
            # server" after the analysis has finished, at the report-writing step.
            # Nothing here needs a screen, so DISPLAY is removed.
            unset DISPLAY

            # Qualimap's start script uses JAVA_OPTS, when it is set, as the whole
            # list of Java options and then ignores --java-mem-size, so Java would
            # run on its own default memory instead of the one asked for below. A
            # JAVA_OPTS left by a login script or a cluster module does exactly
            # that, so it is removed for this step.
            unset JAVA_OPTS

            # Tell Java there is no screen. Java reads JAVA_TOOL_OPTIONS itself
            # and adds it to the start script's options, so the memory from
            # --java-mem-size is kept, and so is any value the user already set.
            # Java notes it in the log with one line, "Picked up
            # JAVA_TOOL_OPTIONS: ...".
            export JAVA_TOOL_OPTIONS="-Djava.awt.headless=true ${{JAVA_TOOL_OPTIONS:-}}"

            mkdir -p {output.qualimap_dir}
            qualimap bamqc \
              -bam {input.bam} \
              --java-mem-size={resources.java_mem}G \
              -nt {threads} \
              -outdir {output.qualimap_dir} \
              -outformat html > {log} 2>&1 || {{
                echo "" >> {log}
                echo "NOTE: Qualimap exited non-zero. This is a QC report, not an input to any other step, so the run continues; the mapping-quality report will be missing from MultiQC for this sample. Check the log above." >> {log}
              }}
            """


# ═══════════════════════ The adapter screen (every mode) ═══════════════════════
# Does the finished genome still carry sequencing adapters or barcodes? Left in the
# reads, they end up in the assembly where many reads end at the same place: at the
# ends of linear contigs. NCBI screens genome submissions for contamination, and its
# FCS-adaptor, provided to clean genomes before submission to GenBank, trims adaptor
# found within 100 bp of a sequence end. BacFlux trims ONT reads before
# assembly (shared/05_ont_adapters.smk); this screen checks contigs_final.fasta in
# every mode, including assemblies made elsewhere (contigs mode) and Illumina
# assemblies. It reports and removes nothing. In the hybrid genome of K. pneumoniae
# TUM24772 made without trimming, it finds the 112 bp of rapid-kit construct that end
# the linear contig_2.
#
# Chain, per sample:
#   contigs_final.fasta (+ the Bakta replicon table in nanopore/hybrid)
#   + the adapter lists (shipped ONT, PacBio and MGI records; UniVec's Illumina
#     records from univec_illumina; the sequences trimmed from this sample's reads)
#     --adapter_search_input--> query FASTA (kept) + the genome with circular contigs
#                               extended (temporary)
#     --blast_adapters-------> the BLAST hits (kept)
#     --screen_adapters------> the adapter loci (kept) + the sample's MultiQC row (kept)
# The method and its thresholds are explained in the opening text of
# scripts/20_qc/adapter_screen.py and in docs/methods_ont_adapters.md.
#
# Why the genome as delivered, and not Flye's draft: it is the genome every later step
# annotates and the one a user submits, it exists in all four modes, and dnaapler
# rotates only circular contigs, whose ends are no ends at all. A linear contig keeps
# Flye's ends through dnaapler, Medaka and Polypolish.


# ──────────────── The Illumina adapters (UniVec), once per run ────────────────
# Takes in: nothing from other rules; links.univec_link (UNIVEC_LINK, by default
#           build 10.1 in NCBI's UniVec archive) and links.univec_sha256
#           (UNIVEC_SHA256).
# Does:     downloads UniVec (about 1.7 MB) with wget, checks its SHA-256 checksum,
#           keeps the records whose description names Illumina (151 in build 10.1:
#           TruSeq, Nextera and other adapters and primers) and deletes the rest. A
#           file that fails the check is renamed with a .failed suffix and kept, and
#           the run stops. awk renames each kept record from UniVec's
#           gnl|uv|NGB00747.1:1-70 to Illumina_UniVec_NGB00747.1:1-70, so that the
#           screen's tables say at a glance that a hit is an Illumina adapter: the
#           variable keep is set at every header line (starting with >) and every line
#           is printed while it is true; [|] stands for a literal |, and the doubled
#           braces are Snakemake's way of writing single braces in a shell command.
#           grep -c ends with status 1 when it counts nothing, which the strict shell
#           would take for a failure; "|| true" keeps the count of 0, and a NOTE then
#           says that the screen searches the shipped lists only. When the checksum
#           does not match, the message gives the checksum of the file received, to
#           put into links.univec_sha256 if that file is the one meant.
# Produces: UNIVEC_ILLUMINA, the Illumina records (FASTA), and UNIVEC_PROVENANCE, one
#           item per line: the address, the date in UTC, the checksum and the number
#           of records. Both in 02.assembly/adapter_lists/, kept, and shared by every
#           sample of the output folder.
# Consumed by: adapter_search_input, for every sample.
#
# Why downloaded and not shipped: UniVec notes on every Illumina record that the
# oligonucleotide sequence is copyright Illumina, Inc. The ONT, PacBio and MGI records
# are shipped in workflow/data/adapters/ (its README names their sources).
#
# Where it runs: on the computer that runs Snakemake (localrule: True), because it
# downloads; `snakemake download_all` fetches it too. wget comes from
# envs/download.yaml; sha256sum, awk, grep and date from the computer itself.
rule univec_illumina:
    localrule: True
    output:
        fasta = UNIVEC_ILLUMINA,
        provenance = UNIVEC_PROVENANCE,
    params:
        url = UNIVEC_LINK,
        sha256 = UNIVEC_SHA256,
        download = UNIVEC_ILLUMINA + ".download",
    conda:
        "../../envs/download.yaml"
    log:
        LOGS + "/univec_illumina.log"
    shell:
        """
        mkdir -p "$(dirname {output.fasta:q})"

        echo "downloading {params.url}" > {log}
        wget -O {params.download:q} {params.url:q} >> {log} 2>&1

        echo "checking the file against links.univec_sha256 ({params.sha256})" >> {log}
        if ! printf '%s  %s\\n' {params.sha256:q} {params.download:q} | sha256sum -c - >> {log} 2>&1; then
            got=$(sha256sum {params.download:q} | cut -d ' ' -f 1)
            mv {params.download:q} {params.download:q}.failed
            echo "[BacFlux] univec_illumina: the file downloaded from {params.url} has the checksum $got, not the one in links.univec_sha256, so it was NOT used (the file is kept as {params.download}.failed; details in {log}). If this is the file you meant to use, put $got into links.univec_sha256." | tee -a {log} >&2
            exit 1
        fi

        awk '/^>/ {{ keep = /Illumina/; if (keep) sub(/^>gnl[|]uv[|]/, ">Illumina_UniVec_") }} keep' {params.download:q} > {output.fasta:q}
        rm {params.download:q}
        records=$(grep -c '^>' {output.fasta:q} || true)
        if [ "$records" = 0 ]; then
            echo "[BacFlux] NOTE univec_illumina: no record of {params.url} names Illumina, so the adapter screen searches for the shipped adapters (ONT, PacBio, MGI) only." | tee -a {log} >&2
        fi

        echo "source_url: {params.url}" > {output.provenance:q}
        echo "downloaded_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> {output.provenance:q}
        echo "sha256: {params.sha256}" >> {output.provenance:q}
        echo "illumina_records: $records" >> {output.provenance:q}
        echo "done: $records Illumina records" >> {log}
        """


# ──────────────── What to search, and where (adapter_search_input) ────────────────
# Takes in: contigs = FINAL_CONTIGS, the finished genome; replicons =
#           BAKTA_REPLICON_INPUT, the Bakta replicon table (build_replicons,
#           shared/15_replicons.smk) in nanopore and hybrid mode and an empty list in
#           the other two, where the topology is not known; illumina =
#           UNIVEC_ILLUMINA from univec_illumina; trimmed = ADAPTER_SCREEN_TRIMMED,
#           the sequences fastplong trimmed from this sample's reads
#           (ONT_ADAPTER_SEQUENCE, rule find_ont_adapter) when BacFlux trimmed them,
#           otherwise an empty list. The shipped lists (ADAPTER_LIST_FILES:
#           workflow/data/adapters/ont_chtd.fasta and univec_pacbio_mgi.fasta) go in
#           as a param, not an input: a file in the repository gets a new modification
#           time at every checkout, and as an input that alone would re-run every
#           sample's screen and the MultiQC report. The cost: a change to a shipped
#           list does not re-run the screen by itself.
# Does:     python adapter_screen.py prepare. It joins every adapter record into one
#           query FASTA, and writes the genome with each CIRCULAR contig extended by
#           its own first 200 bp: a circular contig's ends are an arbitrary cut, and
#           an adapter lying across it would otherwise be seen as two short pieces,
#           or not at all. An empty list renders as nothing in the command, and the
#           script reads the option without a value as "not given".
# Produces: query = ADAPTER_QUERY (kept: screen_adapters reads it to know where the
#           Mu end sits on each query, and it records what was searched for);
#           subject = ADAPTER_SUBJECT, temporary, deleted after BLAST.
# Consumed by: blast_adapters; screen_adapters (the query).
rule adapter_search_input:
    input:
        contigs = FINAL_CONTIGS,
        replicons = BAKTA_REPLICON_INPUT,
        illumina = UNIVEC_ILLUMINA,
        trimmed = ADAPTER_SCREEN_TRIMMED,
    output:
        query = ADAPTER_QUERY,
        subject = temp(ADAPTER_SUBJECT),
    params:
        script = ADAPTER_SCREEN_SCRIPT,
        lists = ADAPTER_LIST_FILES,
    conda:
        "../../envs/python.yaml"
    log:
        LOGS + "/adapter_search_input_{sample}.log"
    shell:
        """
        python {params.script:q} prepare \
          --contigs {input.contigs:q} \
          --replicons {input.replicons:q} \
          --lists {params.lists:q} {input.illumina:q} \
          --trimmed {input.trimmed:q} \
          --out-query {output.query:q} \
          --out-subject {output.subject:q} > {log} 2>&1
        """


# ──────────────── The search (blast_adapters) ────────────────
# Takes in: query and subject from adapter_search_input.
# Does:     blastn with the adapters as queries and the genome as subject. -task
#           blastn-short is BLAST's setting for queries shorter than 50 bp (word size
#           7), and -evalue 1 keeps short matches, which screen_adapters then judges
#           by length, identity and position. -subject reads the FASTA directly, so no
#           BLAST database is built; on a bacterial genome the search takes seconds.
#           BLAST names every contig and adapter exactly as written in the FASTA
#           headers, which is what screen_adapters joins on.
# Produces: hits = ADAPTER_BLAST, the hits in BLAST's tabular format, columns named
#           in the -outfmt below. Kept: every hit, also those that do not pass the
#           screen's rules, so the locus table can be checked against it.
# Consumed by: screen_adapters.
rule blast_adapters:
    input:
        query = ADAPTER_QUERY,
        subject = ADAPTER_SUBJECT,
    output:
        hits = ADAPTER_BLAST,
    conda:
        "../../envs/blast.yaml"
    log:
        LOGS + "/blast_adapters_{sample}.log"
    shell:
        """
        blastn -task blastn-short \
          -query {input.query:q} \
          -subject {input.subject:q} \
          -evalue 1 \
          -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen" \
          -out {output.hits:q} 2> {log}
        """


# ──────────────── The verdict (screen_adapters) ────────────────
# Takes in: hits = ADAPTER_BLAST; contigs = FINAL_CONTIGS (for the contig lengths);
#           query = ADAPTER_QUERY; replicons = BAKTA_REPLICON_INPUT (the topology, or
#           an empty list); and, for a sample whose ONT reads BacFlux trimmed, choice
#           and check = the files of rules find_ont_adapter and check_ont_trimming
#           (ADAPTER_SCREEN_CHOICE, ADAPTER_SCREEN_CHECK; empty lists otherwise).
# Does:     python adapter_screen.py report. A hit counts at 90 % identity or more,
#           from 20 bp within 100 bp of the end of a linear contig, from 28 bp (a whole
#           ONT adapter) inside a contig; a circular contig has no ends. Overlapping
#           hits form one locus. A locus made only of the phage Mu end that follows
#           every rapid barcode is listed but not counted when it lies inside a
#           contig, because a Mu prophage holds the same 50 bp. Every locus gets a
#           counted (yes/no) column and a reason. The script's opening text gives the
#           reasons for every rule. Its summary line is copied from the log to the
#           console.
# Produces: loci = ADAPTER_SCREEN, one row per locus (a header line alone when there
#           is none); mqc = ADAPTER_SCREEN_MQC, the sample's row of the Adapters table
#           of the MultiQC report. Both kept.
# Consumed by: multiqc (shared/90_report.smk) reads the MultiQC row; _downstream_targets()
#              in 00_common.smk asks for the locus table, which nothing else reads.
rule screen_adapters:
    input:
        hits = ADAPTER_BLAST,
        contigs = FINAL_CONTIGS,
        query = ADAPTER_QUERY,
        replicons = BAKTA_REPLICON_INPUT,
        choice = ADAPTER_SCREEN_CHOICE,
        check = ADAPTER_SCREEN_CHECK,
    output:
        loci = ADAPTER_SCREEN,
        mqc = ADAPTER_SCREEN_MQC,
    params:
        script = ADAPTER_SCREEN_SCRIPT,
    conda:
        "../../envs/python.yaml"
    log:
        LOGS + "/screen_adapters_{sample}.log"
    shell:
        """
        python {params.script:q} report \
          --blast {input.hits:q} \
          --contigs {input.contigs:q} \
          --replicons {input.replicons:q} \
          --query {input.query:q} \
          --sample {wildcards.sample:q} \
          --read-choice {input.choice:q} \
          --read-check {input.check:q} \
          --out {output.loci:q} \
          --out-mqc {output.mqc:q} > {log} 2>&1

        grep -F '[BacFlux]' {log} >&2 || true
        """
