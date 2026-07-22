# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — Stage 20 assembly-QC module (rules/shared/20_qc.smk)
#
# Three questions about the finished genome, answered by three tools:
#   QUAST    — is the ASSEMBLY any good?  (contig count, N50, total length, GC)
#   CheckM   — is the GENOME complete, and is it one organism?  (single-copy
#              marker genes: completeness % and contamination %)
#   Qualimap — did the READS map back sensibly?  (depth, evenness, insert size)
#
# The first two need a genome; the third needs the BAM that shared/10_decontam.smk
# already produced for BlobTools.
#
# THE HYBRID COMPLICATION, and how this module deals with it:
# three modes deliver one genome per sample, but hybrid delivers the ONT+Polypolish
# genome AND keeps the decontaminated Illumina draft, and v1 BacFluxL+ deliberately
# QC'd and classified both so they could be compared. Rather than duplicating rule
# bodies or inventing a second wildcard, ONE small rule (stage_qc_genomes) copies
# this sample's 1-or-2 assemblies into a single temporary directory under
# controlled names, and QUAST, CheckM and GTDB-Tk (shared/30_taxonomy.smk) all
# read that one directory. The names, roles and labels come from the QC_GENOMES
# list in 00_common — one source of truth for the staged filenames, the QUAST
# labels, the CheckM/GTDB-Tk bin ids and the MultiQC relabelling.
#
# Two rejected alternatives, recorded so they are not re-litigated:
#   * fan out over a {genome} wildcard — would run GTDB-Tk twice per hybrid
#     sample. classify_wf loads the reference tree and skani sketches on every
#     invocation, so that doubles the single most expensive step in the pipeline;
#     CheckM's lineage_wf likewise shares its HMM/pplacer setup across bins in one
#     run. It would also need a new wildcard and an extra directory level in the
#     three single-genome modes for no gain.
#   * `if IS_HYBRID:` with two near-duplicate rule bodies — duplicated shell is
#     exactly what this unification exists to delete, and duplicates drift.
#
# Data flow (top to bottom):
#
#   FINAL_CONTIGS ────────┐
#   (+ DECONTAM_CONTIGS   ├─► stage_qc_genomes ─► eval/genomes/  (temp)
#      in hybrid only)    │                    └─► {sample}_qc_genomes.tsv (kept)
#                         │                            │
#            ┌────────────┴────────────┬───────────────┴──────────┐
#            ▼                         ▼                          ▼
#   genome_assembly_evaluation  completeness_and_       taxonomic_assignment
#           (QUAST)              contamination (CheckM)  (GTDB-Tk, 30_taxonomy)
#            │                         │                          │
#            └──────────────► MultiQC (shared/90_report.smk) ◄─────┘
#
#   DECONTAM_BAM ─► map_evaluation (Qualimap) ─► MultiQC       [modes with reads]
#
# Everything referenced here is defined once in 00_common.smk and never
# re-derived: QC_GENOMES, qc_genome_fastas, qc_stage_commands,
# qc_genome_table_text, QC_GENOMES_DIR, QC_GENOME_TABLE, QUAST_DIR, CHECKM_DIR,
# CHECKM_STATS, CHECKM_LINEAGE, QUALIMAP_DIR, DECONTAM_BAM, LOGS, RAM,
# capped_cpus, HAS_READS.
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/x.yaml" climbs shared/ -> rules/ -> workflow/ -> workflow/envs/x.yaml.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: stage_qc_genomes — put this sample's genome(s) in one place ─────────
# Not a biology step: pure bookkeeping that makes the three tools below simple.
#
# Takes in: an input FUNCTION, because the NUMBER of files varies by mode —
#           qc_genome_fastas() returns FINAL_CONTIGS in illumina/nanopore/contigs,
#           and [DECONTAM_CONTIGS (Illumina draft), FINAL_CONTIGS (ONT genome)] in
#           hybrid. A plain templated string cannot express "1 or 2 files".
# Does:     copies each one into eval/genomes/ named after its bin id —
#           {sample}.fasta, or {sample}_illumina.fasta + {sample}_ont.fasta in
#           hybrid. Those basenames become the ids CheckM, GTDB-Tk and QUAST
#           report under, which is why they come from the shared QC_GENOMES list.
# Produces:
#   genomes_dir  = a temp() DIRECTORY. temp() is the whole point: Snakemake OWNS
#                  the cleanup and removes it once the last consumer (QUAST,
#                  CheckM, GTDB-Tk) is finished. In v1 there was no such rule —
#                  CheckM copied the genomes into its own output directory and the
#                  taxonomy rule ended with `rm -rf {input.checkm_dir}/*.fasta`,
#                  i.e. one rule deleting files inside another rule's declared
#                  output. Re-running taxonomy alone then found an empty
#                  --genome_dir and GTDB-Tk classified nothing, with no useful
#                  error. No rule deletes anything now.
#   genome_table = a kept TSV: bin_id, role, technology, source_path. It states
#                  which of the two hybrid rows is the DELIVERED genome
#                  ("primary") and which is the comparator, follows this project's
#                  audit-file convention, and gives the later mobilome module an
#                  unambiguous key instead of a filename-suffix guess.
# Consumed by: genome_assembly_evaluation, completeness_and_contamination (this
#              file) and taxonomic_assignment (shared/30_taxonomy.smk).
#
# SCHEDULING NOTE (v1 hazard, now visible rather than fixed): a hybrid sample does
# twice the CheckM/GTDB-Tk work inside ONE job, under ONE cpus reservation. That
# was already true in v1; it is stated here and in genome_table so nobody is
# surprised by the runtime.
#
# conda: NONE — mkdir and cp only, from the launch environment.
rule stage_qc_genomes:
    input:
        genomes = qc_genome_fastas,
    output:
        genomes_dir = temp(directory(QC_GENOMES_DIR)),
        genome_table = QC_GENOME_TABLE,
    params:
        # Both are generated at PARSE time from QC_GENOMES, so the dry run and the
        # log show the literal `cp` commands and the literal table content rather
        # than a loop over hidden bash arrays.
        copy_cmds = qc_stage_commands,
        table_text = qc_genome_table_text,
    log:
        LOGS + "/stage_qc_genomes_{sample}.log"
    priority: 6
    shell:
        # {params.copy_cmds} expands to one `cp` line per genome. Snakemake
        # formats the shell TEMPLATE only — it does not re-scan substituted values
        # for {} placeholders — so multi-line generated commands are safe here.
        # The heredoc is quoted ('EOF'), so the table text is written literally,
        # with no shell expansion and no printf format-string hazards.
        """
        mkdir -p {output.genomes_dir}

        {params.copy_cmds}

        cat > {output.genome_table} << 'EOF'
{params.table_text}
EOF

        ls -l {output.genomes_dir} > {log} 2>&1
        """


# ── Rule: genome_assembly_evaluation — assembly contiguity metrics (QUAST) ────
# Biology: QUAST reports the structural quality of the assembly — how many contigs
# it took to represent the genome, N50, largest contig, total length and GC. It
# says nothing about biological completeness (that is CheckM's job); it tells you
# how fragmented the assembly is, which is what limits everything that depends on
# gene context.
#
# Takes in: the staged genomes DIRECTORY. The shell globs *.fasta inside it, so
#           hybrid gets both assemblies side by side in ONE QUAST report — exactly
#           the comparison you want from a hybrid run.
# Produces: 02.assembly/{sample}/eval/quast/ (a directory; QUAST names the files).
# Consumed by: MultiQC (shared/90_report.smk).
#
# Two deliberate v1->v2 additions, both purely additive (no existing output
# changes), both to be listed in the changelog:
#   1. nanopore had NO QUAST rule at all in v1. It runs here, using the
#      quast.yaml env that already exists in the v2 tree — no new dependency.
#   2. hybrid QUAST'd only the Illumina draft in v1; the delivered ONT genome was
#      never evaluated. It now evaluates both.
# Consequence to be aware of: QUAST's assembly label is the input file's basename,
# so it changes from v1's `contigs_sel` to `{sample}` (or {sample}_illumina /
# {sample}_ont). The MultiQC rename patterns in 90_report.smk match the new form.
#
# (v1 message: "--- QUAST: Genome assembly evaluation. ---")
rule genome_assembly_evaluation:
    input:
        genomes_dir = QC_GENOMES_DIR,
    output:
        quast_dir = directory(QUAST_DIR),
    conda:
        "../../envs/quast.yaml"
    resources:
        cpus = capped_cpus(24)
    log:
        LOGS + "/assembly_evaluation_{sample}.log"
    priority: 5
    shell:
        # The *.fasta glob is expanded by the shell, not by Snakemake; the files
        # are guaranteed to exist because genomes_dir is an input of this rule.
        """
        quast \
          {input.genomes_dir}/*.fasta \
          -o {output.quast_dir} \
          --no-icarus \
          -t {resources.cpus} > {log} 2>&1
        """


# ── Rule: completeness_and_contamination — is this one complete genome? (CheckM) ─
# Biology: CheckM places the genome in the reference tree, picks the lineage-
# specific set of single-copy marker genes, and counts them. Markers missing =>
# incomplete assembly; markers present in more copies than expected => either
# contamination (more than one organism in the bin) or a genuine duplication.
# This is the standard pass/fail gate for an isolate genome, and it is exactly why
# decontamination (shared/10_decontam.smk) has to happen first: a single
# contaminant contig set inflates the contamination figure.
#
# Takes in: the staged genomes DIRECTORY (1 genome, or 2 in hybrid — CheckM
#           handles multiple "bins" in one run and reports one row each).
# Does:     lineage_wf (place + count markers), then qa -o 2 --tab_table to write
#           the extended per-bin statistics as a parsable TSV.
# Produces:
#   checkm_dir     = 02.assembly/{sample}/eval/checkm/ (CheckM's own working +
#                    result tree)
#   checkm_stats   = {sample}_checkm_stats.tsv — the table MultiQC reads
#   checkm_lineage = lineage.ms — CheckM's marker-set file, produced by lineage_wf
#                    and consumed by qa (the name is CheckM's, not ours)
# Consumed by: MultiQC (shared/90_report.smk).
#
# Deliberate v1->v2 fix: v1 ran `checkm lineage_wf -x fasta <dir> <dir>` with the
# SAME directory as input and output — which is the only reason the staged FASTAs
# had to be deleted afterwards by another rule. Input and output are separate now
# ({input.genomes_dir} -> {output.checkm_dir}), which is what lets Snakemake own
# the staged directory as temp().
#
# Filename note: the stats file gained a {sample}_ prefix (v1: bare
# checkm_stats.tsv). Under the D1 layout the containing directory is called
# "checkm" in every sample, so the report's staging loop can no longer recover the
# sample name from the directory — it takes it from the filename instead.
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
    resources:
        cpus = capped_cpus(24)
    log:
        LOGS + "/completeness_and_contamination_{sample}.log"
    priority: 5
    shell:
        # `checkm qa`'s second positional argument is the lineage_wf OUTPUT
        # directory (it holds bins/ and storage/), not the directory of input
        # FASTAs — so {output.checkm_dir} is correct on that line.
        """
        mkdir -p {output.checkm_dir}

        checkm lineage_wf \
          -t {resources.cpus} \
          -x fasta {input.genomes_dir} \
          {output.checkm_dir} > {log} 2>&1

        checkm qa \
          -o 2 \
          -t {resources.cpus} \
          --tab_table \
          -f {output.checkm_stats} {output.checkm_lineage} {output.checkm_dir} >> {log} 2>&1
        """


# ── Mapping QC — only where real reads were mapped ────────────────────────────
# Gated on HAS_READS (the capability flag, not the mode name, per D7). Contigs
# mode has a BAM too, but it is the contigs aligned against THEMSELVES, made only
# to satisfy `blobtools create -b`; a mapping-quality report on a self-alignment
# is meaningless, so v1 FastaFlux had no Qualimap rule and v2 keeps that absence.
# It is a legitimate mode-specific absence, not a gap.
if HAS_READS:

    # ── Rule: map_evaluation — read-alignment quality report (Qualimap) ──────
    # Biology: summarises how this sample's reads sit on its own assembly — mean
    # depth and its evenness, GC bias, mapping rate, insert size (paired reads).
    # Uneven or unexpectedly low coverage is an early warning of a mixed culture,
    # a mis-assembly, or simply not enough data.
    #
    # Takes in: DECONTAM_BAM — the alignment shared/10_decontam.smk made for
    #           BlobTools (bowtie2 in the short-read modes, minimap2 map-ont in
    #           nanopore). The BAM is temp(); depending on it here keeps it alive
    #           until this rule has run.
    # Produces: 02.assembly/{sample}/eval/qualimap/ (an HTML report directory).
    # Consumed by: MultiQC (shared/90_report.smk).
    #
    # Deliberate v1->v2 unification: v1 illumina passed neither --java-mem-size nor
    # -nt (single-threaded, JVM default heap, which OOMs on large BAMs); nanopore
    # hard-coded 24G; hybrid used the config-aware min(RAM, 64) plus -nt. v2 uses
    # the hybrid form everywhere. The RESULTS are unchanged — this only makes the
    # rule faster and stops it running out of heap.
    #
    # Resource note: java_mem is a GIGABYTE figure, not a CPU count. It is passed
    # straight through to the JVM and is NOT scheduled against, so it sits outside
    # this project's `--resources cpus=N` convention.
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
        resources:
            cpus = capped_cpus(24),
            java_mem = min(RAM, 64)
        log:
            LOGS + "/map_evaluation_{sample}.log"
        priority: 5
        shell:
            """
            qualimap bamqc \
              -bam {input.bam} \
              --java-mem-size={resources.java_mem}G \
              -nt {resources.cpus} \
              -outdir {output.qualimap_dir} \
              -outformat html > {log} 2>&1
            """
