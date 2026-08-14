# Stage 20 assembly QC: three questions about the finished genome, one tool each.
#   QUAST    — is the ASSEMBLY any good? Contig count, N50, total length, GC.
#   CheckM   — is the GENOME complete, and is it one organism? Single-copy marker
#              genes, reported as completeness % and contamination %.
#   Qualimap — did the READS map back sensibly? Depth, evenness, insert size.
# The first two need a genome; the third needs the BAM shared/10_decontam.smk
# already made for BlobTools.
#
# stage_qc_genomes               : copies this sample's 1-or-2 assemblies into
#                                  one temporary directory under controlled
#                                  names, and writes the table saying which is
#                                  which.
# genome_assembly_evaluation     : QUAST over that directory.
# completeness_and_contamination : CheckM lineage_wf then qa, same directory.
# map_evaluation                 : Qualimap on the decontamination BAM. Only in
#                                  modes that have reads — see its banner below.
#
# GTDB-Tk reads the same staged directory from shared/30_taxonomy.smk, and all
# four results end up in MultiQC (shared/90_report.smk).
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
# QUALIMAP_DIR, DECONTAM_BAM, LOGS, RAM, capped_cpus, HAS_READS.
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
    priority: 6
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
    priority: 5
    shell:
        # The shell expands the *.fasta glob, not Snakemake; the files are
        # guaranteed to be there because genomes_dir is an input of this rule.
        """
        quast \
          {input.genomes_dir}/*.fasta \
          -o {output.quast_dir} \
          --no-icarus \
          -t {threads} > {log} 2>&1
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
# Under the D1 layout every sample's containing directory is called "checkm", so
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
    log:
        LOGS + "/completeness_and_contamination_{sample}.log"
    priority: 5
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
# Gated on HAS_READS — the capability flag rather than the mode name, per D7 — so
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
    # single-threaded on the JVM's default heap and ran out of memory on large
    # BAMs; nanopore hard-coded 24G; hybrid used the config-aware min(RAM, 64) plus
    # -nt. v2 uses the hybrid form everywhere. The RESULTS are unchanged — this
    # only makes the rule faster and stops it dying on heap.
    #
    # java_mem is a GIGABYTE figure, not a core count, which is why it sits in
    # `resources:` and not `threads:`. It is handed straight to the JVM, and
    # Snakemake schedules against it only if the launch line passes
    # `--resources java_mem=N`.
    #
    # (v1 message: "--- Qualimap: Mapping evaluation. ---"; the nanopore v1 rule
    #  was named map_qc — D5 unifies on map_evaluation.)
    rule map_evaluation:
        input:
            bam = DECONTAM_BAM,
        output:
            qualimap_dir = directory(QUALIMAP_DIR),
        conda:
            "../../envs/qualimap.yaml"
        threads: capped_cpus(24)
        resources:
            java_mem = min(RAM, 64)
        log:
            LOGS + "/map_evaluation_{sample}.log"
        priority: 5
        shell:
            """
            qualimap bamqc \
              -bam {input.bam} \
              --java-mem-size={resources.java_mem}G \
              -nt {threads} \
              -outdir {output.qualimap_dir} \
              -outformat html > {log} 2>&1
            """
