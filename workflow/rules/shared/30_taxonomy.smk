# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — Stage 30 taxonomy module (rules/shared/30_taxonomy.smk)
#
# One question, one tool: WHAT IS THIS ISOLATE? GTDB-Tk answers it against the
# Genome Taxonomy Database, which classifies genomes by genome-wide relatedness
# (concatenated marker proteins placed on a reference tree, plus average
# nucleotide identity to the closest reference) rather than by 16S rRNA alone.
# That matters for a food/feed-safety workflow: the intrinsic-versus-acquired
# reading of any AMR gene found later depends on knowing the species correctly.
#
# GTDB names are not always NCBI names. That is expected and is a GTDB property,
# not a bug — note it before comparing this output with an NCBI-based call.
#
# Data flow:
#
#   eval/genomes/  ─► taxonomic_assignment (GTDB-Tk classify_wf) ─► 03.taxonomy/{sample}/
#   (from stage_qc_genomes,                                            │
#    shared/20_qc.smk)                                                 ▼
#                                                        MultiQC (shared/90_report.smk)
#
# Hybrid: the staged directory holds TWO genomes ({sample}_illumina and
# {sample}_ont), so one classify_wf run emits a summary with two rows keyed by
# those staged basenames. They do not need separate output files — the report
# module rewrites the two keys into "taxonomy Illumina | {sample}" and
# "taxonomy ONT | {sample}", and {sample}_qc_genomes.tsv records which one is the
# delivered ("primary") genome. Both the staged names and the report labels come
# from the single QC_GENOMES list in 00_common.
#
# Everything referenced here is defined once in 00_common.smk and never
# re-derived: QC_GENOMES_DIR, GTDBTK_DIR, GTDBTKDB, LOGS, capped_cpus.
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/gtdbtk.yaml" climbs shared/ -> rules/ -> workflow/ -> workflow/envs/.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: taxonomic_assignment — genome-based classification (GTDB-Tk) ────────
# Biology: classify_wf identifies the ~120 bacterial marker proteins in the
# genome, aligns and places them on the GTDB reference tree, then refines the call
# with skani ANI comparisons against the closest reference genomes. The result is
# a full lineage from domain to species, with the ANI evidence behind it.
#
# Takes in: genomes_dir = QC_GENOMES_DIR — the temporary staging directory that
#           stage_qc_genomes (shared/20_qc.smk) filled with this sample's 1-or-2
#           assemblies, named after their bin ids. NOTE this is the STAGED
#           directory, NOT CheckM's output directory: in v1 GTDB-Tk read CheckM's
#           directory, so its genome set was a side effect of another rule's shell
#           block and the two rules were needlessly serialised.
# Does:     one classify_wf run per sample, over every .fasta in that directory.
# Produces: 03.taxonomy/{sample}/ (a directory; GTDB-Tk names the files inside,
#           notably classify/gtdbtk.bac120.summary.tsv). This is exactly what
#           _downstream_targets() in 00_common requests.
# Consumed by: MultiQC (shared/90_report.smk), which reads the summary TSV.
#
# Two deliberate v1->v2 changes, both to be listed in the changelog:
#
#  1. The trailing `rm -rf {input.checkm_dir}/*.fasta` is DELETED. It removed
#     files inside another rule's declared output, which broke re-entrancy:
#     re-running taxonomy on its own found an empty --genome_dir and GTDB-Tk
#     classified nothing without a useful error. Cleanup is Snakemake's job now,
#     through the temp() on the staged directory. This bug was present in all four
#     v1 workflows; it is fixed once, here.
#
#  2. --pplacer_cpus is unified at the same cap as --cpus (24). v1 used
#     min(CPUS, 64) in illumina/contigs and min(CPUS, 24) in nanopore/hybrid — an
#     accident rather than a decision. pplacer's MEMORY scales with its thread
#     count (tens of GB at high counts), so the 64 cap was the riskier of the two.
#     With both at 24 the separate cpus_p resource key is no longer needed and
#     {threads} drives both flags: one knob. Thread count affects speed and
#     memory only, not the classification.
#
# VERIFIED 2026-07-27, moving from R226 to R232: GTDB-Tk's own
# COMPATIBLE_REF_DATA_VERSIONS hard-pins ONE compatible reference release per
# tool version - 2.6.1 only accepted r220/r226, so R232 needs gtdbtk>=2.7.0
# (see workflow/envs/gtdbtk.yaml). That version bump also changes how the skani
# reference is laid out: R226's release directory ships a raw skani reference
# and GTDB-Tk built its OWN sketch CACHE the first time it ran (hence the
# --skani_sketch_dir flag and the mkdir this rule used to do for it, 2.6.1
# only). From R232 onward the release directory ships that sketch PRE-BUILT
# (GTDB-Tk 2.7's release notes: "reduces the database storage footprint from
# 198 GB down to 98 GB"), classify_wf reads it straight out of
# {gtdbtk_db}/skani/ via GTDBTK_DATA_PATH, and --skani_sketch_dir no longer
# exists in the CLI at all (confirmed against `gtdbtk classify_wf --help`) - so
# there is nothing left for this rule to build or point at separately.
#
# (v1 message: "--- GTDB-Tk: Taxonomic assignment. ---")
rule taxonomic_assignment:
    input:
        genomes_dir = QC_GENOMES_DIR,
        # SCHEDULING EDGE, not a data dependency: GTDB-Tk reads only genomes_dir.
        # In all four v1 workflows taxonomic_assignment took CheckM's output as its
        # input, so GTDB-Tk was structurally guaranteed to run AFTER CheckM for a
        # given sample and the two never overlapped. Both tools run pplacer and can
        # each use tens of GB; on a large box (e.g. --cores 64)
        # dropping that edge would let them run concurrently on the same sample and
        # risk an OOM that v1 could not produce. Keeping this input preserves v1's
        # ordering at zero cost.
        checkm_stats = CHECKM_STATS,
    output:
        gtdbtk_dir = directory(GTDBTK_DIR),
    params:
        gtdbtk_db = GTDBTKDB,
    conda:
        "../../envs/gtdbtk.yaml"
    threads: capped_cpus(24)
    log:
        LOGS + "/taxonomic_assignment_{sample}.log"
    priority: 5
    shell:
        # GTDBTK_DATA_PATH is how GTDB-Tk finds its reference release; :q quotes
        # the database path in case it contains spaces.
        """
        mkdir -p {output.gtdbtk_dir}

        GTDBTK_DATA_PATH={params.gtdbtk_db:q} \
        gtdbtk classify_wf \
          -x fasta \
          --genome_dir {input.genomes_dir} \
          --out_dir {output.gtdbtk_dir} \
          --cpus {threads} \
          --pplacer_cpus {threads} > {log} 2>&1
        """
