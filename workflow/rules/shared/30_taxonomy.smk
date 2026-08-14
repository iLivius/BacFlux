# Stage 03 — what is this isolate? One question, one tool: GTDB-Tk.
#
# classify_wf finds the ~120 bacterial marker proteins in the genome, aligns and
# places them on the Genome Taxonomy Database reference tree, then refines the
# call with skani ANI comparisons against the closest reference genomes. The
# answer is a full lineage from domain to species with the ANI evidence behind
# it — a genome-wide call rather than a 16S rRNA one. That matters for a
# food/feed-safety workflow: whether an AMR gene found later reads as intrinsic
# or acquired depends on having the species right in the first place.
#
# GTDB names are not always NCBI names — GTDB splits genera (Pseudomonas_E) and
# uses placeholder species (sp024807945). That is a property of GTDB, not a
# fault, and it is why rule amrfinder_organism in shared/80_mobilome.smk exists
# to translate between the two vocabularies.
#
# Data flow:
#
#   eval/genomes/ ─► taxonomic_assignment (classify_wf) ─► 03.taxonomy/{sample}/
#   (from stage_qc_genomes,                                          │
#    shared/20_qc.smk)                                               │
#                          MultiQC (shared/90_report.smk) ◄──────────┤
#                  amrfinder_organism (shared/80_mobilome.smk) ◄─────┘
#
# Takes in: genomes_dir = QC_GENOMES_DIR, the temporary staging directory that
#           stage_qc_genomes (shared/20_qc.smk) filled with this sample's one —
#           or, in hybrid, two — assemblies, named after their bin ids. This is
#           the STAGED directory, NOT CheckM's output directory: in v1 GTDB-Tk
#           read CheckM's directory, so its genome set was a side effect of
#           another rule's shell block and the two rules were needlessly
#           serialised.
# Does:     one classify_wf run per sample, over every .fasta in that directory.
# Produces: 03.taxonomy/{sample}/ — a directory, because GTDB-Tk names the files
#           inside it, notably classify/gtdbtk.bac120.summary.tsv. This is
#           exactly what _downstream_targets() in 00_common.smk requests.
# Consumed by: MultiQC (shared/90_report.smk), which reads the summary TSV, and
#           — when mobilome.run is true — amrfinder_organism
#           (shared/80_mobilome.smk), where the species call is what unlocks
#           AMRFinderPlus point-mutation detection. An unmapped species there
#           degrades silently, which is the common case for environmental
#           isolates.
#
# Hybrid runs ONE classify_wf over TWO genomes ({sample}_illumina and
# {sample}_ont), so the summary comes back with two rows keyed by those staged
# basenames. They need no separate output files: the report module rewrites the
# keys into "taxonomy Illumina | {sample}" and "taxonomy ONT | {sample}", and
# {sample}_qc_genomes.tsv records which of the two is the delivered ("primary")
# genome. Staged names and report labels both come from the single QC_GENOMES
# list in 00_common.smk.
#
# Two deliberate v1→v2 changes, both still to be written up in the changelog:
#
#  1. The trailing `rm -rf {input.checkm_dir}/*.fasta` is DELETED. It removed
#     files inside another rule's declared output, which broke re-entrancy:
#     re-running taxonomy on its own found an empty --genome_dir and GTDB-Tk
#     classified nothing without a useful error. Cleanup is Snakemake's job now,
#     through the temp() on the staged directory. The bug was present in all four
#     v1 workflows; it is fixed once, here.
#
#  2. --pplacer_cpus is unified at the same cap as --cpus (24). v1 used
#     min(CPUS, 64) in illumina/contigs and min(CPUS, 24) in nanopore/hybrid — an
#     accident rather than a decision. pplacer's MEMORY scales with its thread
#     count (tens of GB at high counts), so the 64 cap was the riskier of the
#     two. With both at 24 the separate cpus_p resource key is gone and {threads}
#     drives both flags: one knob. Thread count changes speed and memory only,
#     never the classification.
#
# Verified 2026-07-27 while moving from GTDB R226 to R232: GTDB-Tk's own
# COMPATIBLE_REF_DATA_VERSIONS hard-pins ONE compatible reference release per
# tool version — 2.6.1 accepted only r220/r226, so R232 needs gtdbtk>=2.7.0
# (pinned at 2.7.2 in workflow/envs/gtdbtk.yaml). That version bump also changes
# how the skani reference is laid out. Under R226 the release directory shipped a
# raw skani reference and GTDB-Tk built its OWN sketch cache the first time it
# ran, which is why this rule used to pass --skani_sketch_dir and mkdir a
# directory for it (2.6.1 only). From R232 the release ships that sketch
# PRE-BUILT — GTDB-Tk 2.7's release notes: "reduces the database storage
# footprint from 198 GB down to 98 GB" — classify_wf reads it straight out of
# {gtdbtk_db}/skani/ via GTDBTK_DATA_PATH, and --skani_sketch_dir no longer
# exists in the CLI at all (confirmed against `gtdbtk classify_wf --help`). So
# there is nothing left for this rule to build or point at separately.
#
# Everything referenced here is defined once in 00_common.smk and never
# re-derived: QC_GENOMES_DIR, CHECKM_STATS, GTDBTK_DIR, GTDBTKDB, LOGS,
# capped_cpus.
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/gtdbtk.yaml" climbs shared/ → rules/ → workflow/ → workflow/envs/gtdbtk.yaml.
#
# (v1 message: "--- GTDB-Tk: Taxonomic assignment. ---")
rule taxonomic_assignment:
    input:
        genomes_dir = QC_GENOMES_DIR,
        # A scheduling edge, not a data dependency: GTDB-Tk reads only
        # genomes_dir. In all four v1 workflows taxonomic_assignment took CheckM's
        # output as its input, so GTDB-Tk was structurally guaranteed to run AFTER
        # CheckM for a given sample and the two never overlapped. Both tools run
        # pplacer and each can use tens of GB; on a large box (say --cores 64)
        # dropping this edge would let them run concurrently on the same sample
        # and risk an out-of-memory kill that v1 could not produce. Keeping the
        # input preserves v1's ordering at zero cost.
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
