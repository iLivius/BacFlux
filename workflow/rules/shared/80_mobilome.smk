# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — Stage 08 mobilome / AMR-mobility module
# (rules/shared/80_mobilome.smk)
#
# THE QUESTION THIS MODULE ANSWERS, once per sample:
#
#     For each AMR gene we found, is it embedded in a mobile genetic element,
#     and if so, how transferable is that element?
#
# That is the evidence behind the distinction regulators actually care about:
# INTRINSIC resistance (species-wide, chromosomal, not transferable) versus
# ACQUIRED resistance (arrived horizontally on a plasmid or transposon, and may
# move again). Full design and the reasoning behind every tool choice is in
# docs/mobilome_module_SPEC.md; the empirically verified tool contracts are in
# docs/mobilome_wpA_ground_truth.md.
#
# OPT-IN. Everything here is gated on `config.mobilome.run` (default false), so a
# normal BacFlux run never builds these envs or executes these rules.
#
# ── Data flow ────────────────────────────────────────────────────────────────
#
#   contigs_final.fasta ─┬─► contig_lengths ──────────────────┐
#                        │                                     │
#                        ├─► isescan ──► is_table ─────────────┤   (WP-C)
#                        │   (IS elements on the genome)       │
#                        │                                     │
#   bakta/{sample}.faa ──┼─► amrfinderplus ──► amr TSV ────────┤   (WP-A)
#   bakta/{sample}.gff3 ─┘        ▲                            │
#                                 │                            ├─► amr_mge_colocalisation
#   03.taxonomy/{sample} ─► amrfinder_organism ────────────────┤        (WP-D)
#              (GTDB-Tk)   (unlocks point mutations)           │           │
#                                                              │           ▼
#   bakta/{sample}.faa ──► conjscan ──► conjugation table ─────┤   {sample}_amr_mobility.tsv
#              (relaxase / T4CP / T4SS = can it self-transmit?) │   + _amr_mobility_audit.tsv
#                                                              │
#   06.plasmids/…/platon ─► replicon calls (chromosome|plasmid)┘
#
# ── Two things a reader must keep in mind about the OUTPUT ───────────────────
#
# 1. ON A FRAGMENTED (short-read) ASSEMBLY THE IS COUNT IS A FLOOR, NOT A COUNT.
#    IS elements are the single biggest cause of contig breaks, because multiple
#    identical copies collapse in the assembly graph. So the very element you are
#    hunting is often what broke the contig you are looking at. Every call
#    therefore carries a contig-boundary flag and a confidence tier, and the IS
#    summary reports what FRACTION of calls sit at a contig end. Read those.
#
# 2. PUBLISHED IS-DETECTION FDR IS 8–24% EVEN ON CURATED DATA. The output is
#    deliberately tiered evidence, never a bare count. And the language is
#    "PREDICTED self-transmissible" — the confirmatory experiment is a filter or
#    broth mating assay, not bioinformatics.
#
# Everything referenced here is defined once in 00_common.smk: MOBILOME_RUN,
# MOBILOME_DIR and the per-file constants, MOBILOME_MAX_COMPOSITE_SPAN,
# MOBILOME_BOUNDARY_BP, the three script paths, FINAL_CONTIGS, DIR_ANNOTATION,
# DIR_TAXONOMY, PLATON_TABLE, LOGS, BAKTADB, capped_cpus.
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/x.yaml" climbs shared/ -> rules/ -> workflow/ -> workflow/envs/x.yaml.
# ─────────────────────────────────────────────────────────────────────────────
if MOBILOME_RUN:

    # ── Rule: contig_lengths — how long is every contig? ─────────────────────
    # Biology: nothing on its own — but every downstream honesty check needs it.
    # "Is this AMR gene 200 bp from the end of its contig?" is only answerable
    # against the contig's length, and that question is what separates "there is
    # no IS next to this gene" from "the contig ended before we could tell".
    #
    # Takes in: FINAL_CONTIGS, the one delivered genome (D2 hand-off), whichever
    #           mode produced it.
    # Does: a stdlib FASTA scan (no biopython) writing contig + length.
    # Produces: {sample}_contig_lengths.tsv.
    # Consumed by: isescan_table and amr_mge_colocalisation.
    #
    # conda: NONE — plain Python from the launch environment, like the other
    # small helper steps in this workflow.
    rule contig_lengths:
        input:
            contigs = FINAL_CONTIGS,
        output:
            lengths = CONTIG_LENGTHS,
        params:
            script = ISESCAN_TABLE_SCRIPT,
        log:
            LOGS + "/mobilome_contig_lengths_{sample}.log"
        priority: 3
        shell:
            """
            python {params.script} contig-lengths \
              --contigs {input.contigs} \
              --out {output.lengths} > {log} 2>&1
            """

    # ── Rule: amrfinder_organism — can we ask for point mutations? ───────────
    # Biology: AMRFinderPlus can additionally report RESISTANCE-CONFERRING POINT
    # MUTATIONS, but only when told which organism it is looking at, because those
    # mutations are curated per species. Point mutations are precisely the
    # INTRINSIC, non-transferable category that ABRicate structurally cannot see —
    # so this small step is what lets the module distinguish "resistant because of
    # a chromosomal mutation" from "resistant because it acquired a gene".
    #
    # The catch: AMRFinderPlus curates only 31 organisms under NCBI names, and
    # GTDB-Tk speaks GTDB names, which split genera (Pseudomonas_E) and use
    # placeholder species (sp024807945). So a careful mapping is needed, and for
    # most environmental isolates there is legitimately NO match — verified: none
    # of the six isolates screened in this project map to a curated organism.
    # A no-match is therefore the COMMON case and must be silent and graceful.
    #
    # Takes in: the GTDB-Tk output directory for this sample (03.taxonomy).
    # Does: map the GTDB classification to a curated AMRFinderPlus --organism, or
    #       to nothing at all.
    # Produces: a one-line file holding the organism name (or an EMPTY file), plus
    #           an audit TSV recording the decision and its reason.
    # Consumed by: amrfinderplus (below).
    rule amrfinder_organism:
        input:
            gtdbtk_dir = GTDBTK_DIR,
        output:
            organism = AMRFINDER_ORGANISM,
            audit = AMRFINDER_ORGANISM_AUDIT,
        params:
            script = ORGANISM_SCRIPT,
        log:
            LOGS + "/mobilome_amrfinder_organism_{sample}.log"
        priority: 3
        shell:
            """
            python {params.script} \
              --gtdbtk-dir {input.gtdbtk_dir} \
              --sample {wildcards.sample} \
              --out-organism {output.organism} \
              --out-audit {output.audit} > {log} 2>&1
            """

    # ── Rule: amrfinderplus — the structured AMR input (WP-A) ────────────────
    # Biology: Bakta already runs AMRFinderPlus internally, but surfaces only the
    # gene name and product. The full report adds five things this module needs:
    #   1. the Method column (EXACTX / BLASTX / PARTIALX / HMM / POINTX …) — a
    #      free confidence tier, straight from the tool;
    #   2. PARTIAL* methods, which flag a hit truncated at a contig end — the
    #      honest signal for short-read fragmentation;
    #   3. element type/subtype (AMR / STRESS / VIRULENCE);
    #   4. drug class and subclass;
    #   5. point mutations, when --organism is known (see the rule above).
    #
    # This costs NO new conda env and NO new database: amrfinder ships inside the
    # Bakta env, and its database ships inside the Bakta database. Verified by
    # running it for real (see docs/mobilome_wpA_ground_truth.md): the 22-column
    # schema there is the contract the co-localisation parser keys on.
    #
    # ABRicate is NOT replaced. The three AMR legs are complementary and stay:
    #   BBMap->CARD (reads)   — immune to assembly collapse
    #   ABRicate (contigs)    — multi-database breadth + the EFSA-threshold report
    #   AMRFinderPlus         — structured, coordinate-bearing input for THIS module
    # EFSA thresholds stay on the ABRicate leg only; imposing a blanket 80/70 cut
    # here would fight AMRFinderPlus's curated per-gene cutoffs.
    #
    # Takes in: Bakta's proteins (.faa) and gene coordinates (.gff3), the genome
    #           itself, and the organism file from the rule above.
    # Does: amrfinder with --annotation_format bakta (verified accepted in 4.2.7),
    #       --plus for the STRESS/VIRULENCE categories, and --organism only when
    #       one was matched. --mutation_all is written only in that same case,
    #       because point mutations are meaningless without an organism.
    # Produces: the AMR TSV (and a mutations TSV, empty when no organism matched).
    # Consumed by: amr_mge_colocalisation.
    rule amrfinderplus:
        input:
            bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
            contigs = FINAL_CONTIGS,
            organism = AMRFINDER_ORGANISM,
        output:
            tsv = AMRFINDER_TSV,
            mutations = AMRFINDER_MUTATIONS,
        params:
            # The versioned symlink, NOT the parent: the parent holds versioned
            # subfolders and the built BLAST index lives inside them, so pointing
            # at the parent fails with "BLAST database ... was not found".
            db = os.path.join(BAKTADB, "amrfinderplus-db", "latest"),
        conda:
            "../../envs/bakta.yaml"
        threads: capped_cpus(8)
        log:
            LOGS + "/mobilome_amrfinderplus_{sample}.log"
        priority: 3
        shell:
            # The organism file is EMPTY whenever no curated organism matched,
            # which is the common case for environmental isolates. An empty file
            # means we simply omit --organism (and --mutation_all with it), and
            # still write a well-formed empty mutations file so the output set is
            # the same either way and nothing downstream has to special-case it.
            """
            organism=$(cat {input.organism} 2>/dev/null | tr -d '[:space:]')

            if [ -n "$organism" ]; then
                echo "Using AMRFinderPlus --organism $organism (point mutations enabled)." > {log}
                amrfinder \
                  -p {input.bakta_dir}/{wildcards.sample}.faa \
                  -n {input.contigs} \
                  -g {input.bakta_dir}/{wildcards.sample}.gff3 \
                  --annotation_format bakta \
                  --plus \
                  --organism "$organism" \
                  --mutation_all {output.mutations} \
                  -d {params.db} \
                  --threads {threads} \
                  -o {output.tsv} >> {log} 2>&1
            else
                echo "No curated AMRFinderPlus organism for this sample; running without --organism (point mutations not available). This is normal for environmental isolates." > {log}
                amrfinder \
                  -p {input.bakta_dir}/{wildcards.sample}.faa \
                  -n {input.contigs} \
                  -g {input.bakta_dir}/{wildcards.sample}.gff3 \
                  --annotation_format bakta \
                  --plus \
                  -d {params.db} \
                  --threads {threads} \
                  -o {output.tsv} >> {log} 2>&1
                # Keep the output set identical on both paths: a header-only file
                # says "we looked and there is nothing", which is not the same as
                # a missing file, and every downstream reader can stay simple.
                printf 'No --organism was matched for this sample, so AMRFinderPlus point-mutation screening was not available.\n' > {output.mutations}
            fi
            """
