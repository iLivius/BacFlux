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
            python {params.script} \
              --genome-fasta {input.contigs} \
              --out-contig-lengths {output.lengths} > {log} 2>&1
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

    # ── Rule: isescan — find insertion sequences on the genome (WP-C) ────────
    # Biology: IS elements are the workhorse of the bacterial mobilome — small,
    # self-mobile, and present in many copies. They matter here for two reasons:
    # an IS beside an AMR gene can supply a promoter (raising expression without
    # moving anything), and two copies of the SAME IS bracketing a gene make a
    # composite transposon that can move the gene within the cell.
    #
    # Takes in: FINAL_CONTIGS.
    # Does: ISEScan with its bundled profile HMMs — no external database.
    #       DELIBERATELY WITHOUT --removeShortIS: that flag drops partial copies,
    #       and on a fragmented assembly a "partial" IS is usually a real IS the
    #       contig ran out on. We keep them and tier them honestly instead.
    # Produces: 07/08 ISESCAN_DIR — declared a DIRECTORY because ISEScan writes a
    #       tree, not one file (see the two gotchas below).
    # Consumed by: isescan_table.
    #
    # THREE VERIFIED GOTCHAS, each handled in the shell below:
    #  1. ISEScan writes NO FILES AT ALL when it finds no IS. A genome with no
    #     detectable IS is a perfectly ordinary result, so the rule creates the
    #     output directory itself and never depends on the tool having written
    #     into it. Without this the run dies on a MissingOutputException for a
    #     biologically normal sample.
    #  2. It CACHES: if a previous run's proteome/ and hmm/ sub-directories are
    #     present it silently reuses them, so a re-run on changed contigs would
    #     report the OLD genome's IS. The directory is cleared first.
    #  3. It must run with its own env's python (it loads libssw.so relative to
    #     the interpreter), so the command calls `isescan.py` BY NAME and lets
    #     the activated conda env resolve it — never by absolute path.
    rule isescan:
        input:
            contigs = FINAL_CONTIGS,
        output:
            isescan_dir = directory(ISESCAN_DIR),
        conda:
            "../../envs/isescan.yaml"
        threads: capped_cpus(8)
        log:
            LOGS + "/mobilome_isescan_{sample}.log"
        priority: 3
        shell:
            # `|| true` on the tool itself, then create the directory
            # unconditionally: "no IS found" exits non-zero on some inputs and is
            # not an error. isescan_table (next rule) decides what the result
            # means, and writes a well-formed empty table when there is nothing.
            """
            rm -rf {output.isescan_dir}
            mkdir -p {output.isescan_dir}

            isescan.py \
              --seqfile {input.contigs} \
              --output {output.isescan_dir} \
              --nthread {threads} > {log} 2>&1 || {{
                echo "" >> {log}
                echo "NOTE: ISEScan exited non-zero. On a genome with no detectable IS it writes no output at all, which is a normal result, so the module continues and the IS table will be empty. Check the log above if you expected IS elements." >> {log}
              }}
            """

    # ── Rule: isescan_table — tidy the IS calls and measure the honesty ──────
    # Takes in: the ISEScan output directory and the contig lengths.
    # Does: normalise ISEScan's 24-column table to one tidy row per IS, and
    #       compute the QC that keeps this module truthful — how many IS calls
    #       sit within MOBILOME_BOUNDARY_BP of a contig end, and what fraction of
    #       the total that is.
    # Produces: IS_TABLE (rows), IS_SUMMARY (the QC line), IS_AUDIT (dropped rows
    #       with a reason).
    # Consumed by: amr_mge_colocalisation; IS_SUMMARY is also a rule-all target,
    #       so the QC is always produced, never optional.
    #
    # WHY THE QC MATTERS: IS elements are the main cause of contig breaks, so a
    # high boundary fraction means the assembly fragmented exactly where the
    # elements are, and the located count is a FLOOR rather than a count.
    rule isescan_table:
        input:
            isescan_dir = ISESCAN_DIR,
            lengths = CONTIG_LENGTHS,
        output:
            table = IS_TABLE,
            summary = IS_SUMMARY,
            audit = IS_AUDIT,
        params:
            script = ISESCAN_TABLE_SCRIPT,
            boundary_bp = MOBILOME_BOUNDARY_BP,
        log:
            LOGS + "/mobilome_isescan_table_{sample}.log"
        priority: 3
        shell:
            """
            python {params.script} \
              --sample {wildcards.sample} \
              --isescan-out {input.isescan_dir} \
              --contig-lengths {input.lengths} \
              --boundary-bp {params.boundary_bp} \
              --out-table {output.table} \
              --out-summary {output.summary} \
              --out-audit {output.audit} > {log} 2>&1
            """

    # ── Rule: conjscan_models — fetch the CONJscan model package once ────────
    # Biology: CONJscan is a set of profile models describing the machinery a
    # cell needs to conjugate — a relaxase (nicks the DNA), a coupling protein,
    # and a type IV secretion system (the mating apparatus). Which of those are
    # present is what separates "can be moved by a helper" from "moves itself".
    #
    # ⚠ LICENCE (verified from the package's own metadata.yml): the CONJscan
    # MODELS are CC BY-NC-SA 4.0 (Institut Pasteur / CNRS) — academic and
    # non-commercial use only. BacFlux never vendors them: they are fetched here,
    # at the user's request, under the user's own agreement with the licensor,
    # exactly as this workflow already treats bakta_db, blast_db and the rest.
    # BacFlux's own MIT licence is unaffected. Turning the mobilome module on
    # means accepting that dependency — it is stated in the README and printed
    # at parse time.
    #
    # Takes in: nothing (a network fetch).
    # Produces: 08.mobilome/conjscan_models/ — one shared copy for all samples.
    # Consumed by: conjscan.
    #
    # GOTCHA: `macsydata` is DEPRECATED in MacSyFinder 2.1.6 and prints a rename
    # warning; the current command is `msf_data`. The shell prefers msf_data and
    # falls back, so the rule works across both spellings.
    rule conjscan_models:
        output:
            models = directory(CONJSCAN_MODELS_DIR),
        conda:
            "../../envs/macsyfinder.yaml"
        log:
            LOGS + "/mobilome_conjscan_models.log"
        priority: 4
        shell:
            """
            mkdir -p {output.models}
            {{
              echo "Installing the CONJScan model package."
              echo "NOTE: these models are licensed CC BY-NC-SA 4.0 (Institut Pasteur/CNRS) - academic / non-commercial use only. They are downloaded here at your request and are never redistributed by BacFlux."
            }} > {log}

            if command -v msf_data >/dev/null 2>&1; then
                msf_data install --target {output.models} CONJScan >> {log} 2>&1
            else
                macsydata install --target {output.models} CONJScan >> {log} 2>&1
            fi
            """

    # ── Rule: conjscan — conjugation machinery on THIS genome ────────────────
    # Takes in: Bakta's proteins (.faa) and the CONJscan models.
    # Does: MacSyFinder with the CONJScan/Chromosome model set. Chromosome, not
    #       Plasmids, because that is the question nothing else in BacFlux can
    #       answer: Platon already reports plasmid mobility from its own table,
    #       but it SKIPS any contig over 500 kb, so it never looks at the
    #       chromosome. Conjugation machinery on a chromosome means an ICE — the
    #       case that breaks the naive "chromosomal, therefore not transferable"
    #       assumption.
    # Produces: 08.mobilome/{sample}/conjscan/ (best_solution.tsv and friends).
    # Consumed by: conjscan_ice.
    #
    # --db-type ordered_replicon is right for one genome's proteome, and Bakta
    # writes its .faa in genome order. CAVEAT recorded in the ground-truth doc:
    # on a fragmented assembly MacSyFinder treats the whole proteome as one
    # pseudo-replicon and can cluster hits ACROSS contigs; conjscan_ice therefore
    # flags any system whose hits span contigs and caps it at low confidence.
    #
    # Most isolates carry no conjugative system at all, so a run finding nothing
    # is the common case and must not fail the pipeline.
    rule conjscan:
        input:
            bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
            models = CONJSCAN_MODELS_DIR,
        output:
            conjscan_dir = directory(CONJSCAN_DIR),
        conda:
            "../../envs/macsyfinder.yaml"
        threads: capped_cpus(8)
        log:
            LOGS + "/mobilome_conjscan_{sample}.log"
        priority: 3
        shell:
            """
            rm -rf {output.conjscan_dir}
            mkdir -p {output.conjscan_dir}

            macsyfinder \
              --models CONJScan/Chromosome all \
              --sequence-db {input.bakta_dir}/{wildcards.sample}.faa \
              --db-type ordered_replicon \
              --models-dir {input.models} \
              --out-dir {output.conjscan_dir} \
              --worker {threads} \
              --force > {log} 2>&1 || {{
                echo "" >> {log}
                echo "NOTE: MacSyFinder exited non-zero. A genome with no conjugative system is the common case for environmental isolates; the module continues and the ICE/IME table will be empty." >> {log}
              }}
            """

    # ── Rule: conjscan_ice — turn machinery hits into ICE / IME candidates ───
    # Takes in: CONJscan's best_solution.tsv, Bakta's GFF3 (for the genomic
    #           coordinates of each protein hit, for integrase genes found by
    #           product regex, AND for the tRNAs the att search anchors on), the
    #           genome itself, the IS table, and the contig lengths.
    # Does: spec §8 Phases 0-4 — collect anchors (relaxase, coupling protein,
    #       T4SS, integrase), cluster them on one contig, and classify:
    #         integrase + relaxase + T4SS -> ICE  (predicted self-transmissible)
    #         integrase + relaxase        -> IME  (mobilisable, needs a helper)
    #         integrase only              -> passive island
    #         relaxase + T4SS, no integrase -> conjugative region, NOT an ICE
    #       then Phase 3, the att-site search: look for the attL/attR direct
    #       repeats that mark where the element really starts and stops, and
    #       widen the interval to them when they are found.
    #
    # WHY THE att SEARCH RUNS IN EVERY MODE, not just long-read. The spec scoped
    # Phase 3 to BacFluxL because of short-read FRAGMENTATION, and that reasoning
    # is about assembly contiguity, not about the sequencer: a closed genome
    # arriving through `contigs` mode has exactly the flanking sequence the search
    # needs. So it is attempted always and degrades honestly — when the flanks are
    # missing, or the element runs off the end of a contig, nothing is found and
    # boundary_method stays 'none' with the machinery span reported unchanged.
    # The existing spans_contigs / at_contig_boundary flags already cap confidence
    # for exactly those cases.
    #
    # The IS table is an input because it is MASKED OUT before the de novo half of
    # the search: insertion sequences carry terminal repeats and duplicate target
    # DNA when they transpose, so an IS-rich neighbourhood is full of direct
    # repeats that have nothing to do with ICE integration. The spec names this as
    # the most likely way to get Phase 3 wrong.
    #
    # Produces: the ICE/IME element table + its audit.
    # Consumed by: amr_mge_colocalisation, as a SECOND element source alongside
    #              the IS table.
    rule conjscan_ice:
        input:
            conjscan_dir = CONJSCAN_DIR,
            bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
            lengths = CONTIG_LENGTHS,
            genome = FINAL_CONTIGS,
            is_table = IS_TABLE,
        output:
            table = ICE_TABLE,
            audit = ICE_AUDIT,
        params:
            script = CONJSCAN_ICE_SCRIPT,
            # Empty unless the user opted in to the strict spec §8 Phase 6 rule
            # (mobilome.require_trna_boundary_for_high). Off by default so that a
            # well-evidenced ICE is not capped at medium merely because a
            # short-read assembly could not show us where it ends.
            strict_boundary = ("--require-trna-boundary-for-high"
                               if MOBILOME_REQUIRE_TRNA_BOUNDARY else ""),
        log:
            LOGS + "/mobilome_conjscan_ice_{sample}.log"
        priority: 3
        shell:
            # best_solution.tsv is absent when MacSyFinder found nothing; the
            # script treats a missing file as "no systems" and still writes a
            # well-formed empty table, so no guard is needed here.
            """
            python {params.script} \
              --sample {wildcards.sample} \
              --conjscan-tsv {input.conjscan_dir}/best_solution.tsv \
              --bakta-gff {input.bakta_dir}/{wildcards.sample}.gff3 \
              --contig-lengths {input.lengths} \
              --genome {input.genome} \
              --is-table {input.is_table} \
              {params.strict_boundary} \
              --out-table {output.table} \
              --out-audit {output.audit} > {log} 2>&1
            """

    # ══ The TnCentral naming layer (ladder tier 4) ═══════════════════════════
    # Only exists when the user configured a TnCentral source. Without it the
    # module behaves exactly as before and tier 4 stays unreachable.
    if MOBILOME_NAME_ELEMENTS:

        # ── Rule: tncentral_db — fetch and index the curated transposon set ──
        # Biology: TnCentral catalogues transposons and integrons that people have
        #          characterised and named. Matching one is qualitatively different
        #          from inferring a composite from two IS copies: the architecture
        #          is already known, so it gets awkward cases like IS26 right for
        #          free (spec §7).
        # Takes in: nothing from the workflow — a URL from the config, or a
        #           directory the user already holds.
        # Does: download, verify, unpack, and rebuild the BLAST database.
        #       The shipped index is BLAST v4; we dump to FASTA and rebuild as v5
        #       so it works with current blast+ (spec §5.3), and record what was
        #       actually fetched in PROVENANCE.txt because the endpoint is
        #       unversioned and otherwise there is no way to say later WHICH
        #       release a result came from.
        # Produces: TNCENTRAL_BLAST_DB (a prefix) + PROVENANCE.txt.
        # Consumed by: tncentral_blast.
        rule tncentral_db:
            output:
                db_dir = directory(TNCENTRAL_DB_DIR),
            params:
                url = TNCENTRAL_URL,
                sha256 = TNCENTRAL_SHA256,
                local_dir = TNCENTRAL_LOCAL,
                # TnCentral's server 403s the default curl user-agent. This is a
                # bot block rather than a wall, so we identify as a browser. If
                # they tighten it the rule fails loudly - it never leaves an empty
                # database behind that would silently produce zero named elements.
                user_agent = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                              "(KHTML, like Gecko) Chrome/120.0 Safari/537.36"),
            conda:
                "../../envs/tncentral.yaml"
            log:
                LOGS + "/mobilome_tncentral_db.log"
            priority: 4
            shell:
                """
                exec > {log} 2>&1
                set -euo pipefail
                mkdir -p {output.db_dir}

                if [ -n "{params.local_dir}" ]; then
                    echo "Using the local TnCentral directory '{params.local_dir}'. Nothing will be downloaded."
                    cp "{params.local_dir}"/tncentral.fa {output.db_dir}/tncentral.fa
                    SOURCE="local:{params.local_dir}"
                else
                    echo "Fetching TnCentral from {params.url}"
                    curl -sSL -A "{params.user_agent}" -o {output.db_dir}/tncentral.zip "{params.url}"

                    # A bot block or an error page returns HTML with a 200, which
                    # would otherwise unzip-fail confusingly or, worse, produce an
                    # empty database and zero named elements with no complaint.
                    if ! unzip -t {output.db_dir}/tncentral.zip > /dev/null 2>&1; then
                        echo "ERROR: the download is not a ZIP archive. The server may have" >&2
                        echo "       blocked the request or changed its API. First bytes:" >&2
                        head -c 200 {output.db_dir}/tncentral.zip >&2
                        exit 1
                    fi

                    OBSERVED=$(sha256sum {output.db_dir}/tncentral.zip | cut -d' ' -f1)
                    if [ -n "{params.sha256}" ] && [ "$OBSERVED" != "{params.sha256}" ]; then
                        echo "ERROR: TnCentral checksum mismatch." >&2
                        echo "       expected {params.sha256}" >&2
                        echo "       observed $OBSERVED" >&2
                        echo "       The endpoint is unversioned, so this means the upstream" >&2
                        echo "       data changed. Update mobilome.tncentral.sha256 once you" >&2
                        echo "       have decided the new release is the one you want." >&2
                        exit 1
                    fi

                    unzip -o -q -j {output.db_dir}/tncentral.zip -d {output.db_dir}
                    SOURCE="{params.url} (sha256 $OBSERVED)"
                    rm -f {output.db_dir}/tncentral.zip
                fi

                # The archive ships a BLAST v4 index. Rebuild as v5 so it works
                # with current blast+ (spec §5.3), and drop the shipped index files
                # so there is no chance of the old one being picked up instead.
                rm -f {output.db_dir}/tncentral.fa.n*
                makeblastdb -in {output.db_dir}/tncentral.fa -dbtype nucl \
                  -out {output.db_dir}/tncentral_v5 -blastdb_version 5

                N_SEQ=$(grep -c '^>' {output.db_dir}/tncentral.fa)
                {{
                  echo "source:      $SOURCE"
                  echo "fetched:     $(date -u +%Y-%m-%dT%H:%M:%SZ)"
                  echo "fasta_sha256: $(sha256sum {output.db_dir}/tncentral.fa | cut -d' ' -f1)"
                  echo "sequences:   $N_SEQ"
                  echo ""
                  echo "The TnCentral download endpoint is UNVERSIONED, so this file is the"
                  echo "only record of which release was used. Quote it in a methods section."
                  echo "TnCentral carries an 'All Rights Reserved' notice; BacFlux ships no"
                  echo "TnCentral data, only this URL. See config mobilome.tncentral."
                }} > {output.db_dir}/PROVENANCE.txt

                echo "TnCentral ready: $N_SEQ sequences."
                """

        # ── Rule: tncentral_blast — where are the curated elements? ──────────
        # Takes in: this sample's contigs + the v5 database above.
        # Does: blastn the whole assembly against the curated set. The whole
        #       assembly rather than the ISEScan calls, because a transposon is
        #       bigger than the IS that bounds it and would be clipped otherwise.
        # Produces: TNCENTRAL_BLAST_HITS, tabular, with the exact columns
        #           name_transposons.py expects (BLAST_COLUMNS there).
        rule tncentral_blast:
            input:
                contigs = FINAL_CONTIGS,
                db_dir = TNCENTRAL_DB_DIR,
            output:
                hits = TNCENTRAL_BLAST_HITS,
            params:
                db = lambda w, input: os.path.join(input.db_dir, "tncentral_v5"),
                # Loose enough that the naming thresholds in the script - not the
                # search - decide what counts, so every near miss is auditable.
                evalue = "1e-20",
            conda:
                "../../envs/tncentral.yaml"
            threads: capped_cpus(8)
            log:
                LOGS + "/mobilome_tncentral_blast_{sample}.log"
            priority: 3
            shell:
                """
                blastn \
                  -query {input.contigs} \
                  -db {params.db} \
                  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore slen qlen" \
                  -evalue {params.evalue} \
                  -num_threads {threads} \
                  -out {output.hits} > {log} 2>&1
                """

        # ── Rule: name_elements — turn BLAST hits into named elements ────────
        # Does: merge HSPs into element COPIES (separate copies of one transposon
        #       must not be joined - see the script), apply the naming thresholds,
        #       drop plain IS entries that ISEScan already covers, and write the
        #       result in the element-table shape colocalise.py consumes.
        # Produces: NAMED_ELEMENTS_TABLE + its discard audit.
        # Consumed by: amr_mge_colocalisation, as a third --is-table.
        rule name_elements:
            input:
                hits = TNCENTRAL_BLAST_HITS,
            output:
                table = NAMED_ELEMENTS_TABLE,
                audit = NAMED_ELEMENTS_AUDIT,
            params:
                script = NAME_TRANSPOSONS_SCRIPT,
                min_identity = TNCENTRAL_MIN_IDENTITY,
                min_coverage = TNCENTRAL_MIN_COVERAGE,
            log:
                LOGS + "/mobilome_name_elements_{sample}.log"
            priority: 3
            shell:
                """
                python {params.script} \
                  --sample {wildcards.sample} \
                  --blast {input.hits} \
                  --min-identity {params.min_identity} \
                  --min-coverage {params.min_coverage} \
                  --out-table {output.table} \
                  --out-audit {output.audit} > {log} 2>&1
                """

    # ── Rule: mobilome_replicons — is each contig chromosome or plasmid? ─────
    # Takes in: the Platon directory this sample's plasmid stage already produced,
    #           plus the genome (so contigs Platon skipped are still listed).
    # Does: read Platon's chromosome/plasmid split and its own # Conjugation /
    #       # Mobilization / # OriT counts, and turn them into a per-contig call
    #       plus a plasmid mobility class.
    # Produces: MOBILOME_REPLICONS.
    # Consumed by: amr_mge_colocalisation, where it decides ladder tiers 5 and 6
    #              for anything sitting on a plasmid.
    rule mobilome_replicons:
        input:
            platon_dir = PLATON_DIR,
            contigs = FINAL_CONTIGS,
            # geNomad's plasmid calls, as a SECOND OPINION on Platon's.
            #
            # Present only when the user opted in to geNomad: rule
            # plasmid_concordance lives under `if PHAGE_CALLER == "genomad"` in
            # 60_plasmid.smk, because geNomad is academic/non-commercial licensed
            # and so cannot be BacFlux's default (spec section 3.1). Unpacking a
            # dict keeps the input list valid in both configurations.
            #
            # Why it matters here: Platon decides tiers 5 and 6 on its own, and a
            # plasmid it misses becomes "chromosomal, intrinsic candidate" for
            # every AMR gene on it - the one error direction that HIDES
            # transferability. The concordance table already knew better; until
            # now nothing read it.
            **({"concordance": PLASMID_CONCORDANCE} if PHAGE_CALLER == "genomad" else {}),
        output:
            replicons = MOBILOME_REPLICONS,
        params:
            script = REPLICONS_MOBILOME_SCRIPT,
            prefix = GENOMAD_PREFIX,
            # Empty string when geNomad was not run, so the script simply behaves
            # as it did before (Platon alone).
            genomad_flag = lambda w: (
                "--genomad-concordance " + PLASMID_CONCORDANCE.format(sample=w.sample)
                if PHAGE_CALLER == "genomad" else ""
            ),
        log:
            LOGS + "/mobilome_replicons_{sample}.log"
        priority: 3
        shell:
            """
            python {params.script} \
              --sample {wildcards.sample} \
              --platon-dir {input.platon_dir} \
              --prefix {params.prefix} \
              --contigs {input.contigs} \
              {params.genomad_flag} \
              --out {output.replicons} > {log} 2>&1
            """

    # ── Rule: amr_mge_colocalisation — the deliverable (WP-D) ────────────────
    # Biology: this is where the module answers its question. For every AMR gene
    # AMRFinderPlus found, look at what mobile elements sit around it and decide
    # where it lands on the mobility ladder:
    #   1 chromosomal, nothing nearby      -> intrinsic candidate
    #   2 IS adjacent, pointing at it      -> expression change, NOT mobilisation
    #   3 between two copies of one IS     -> composite transposon, moves in-cell
    #   4 in a named transposon / integron -> mobilisable, named architecture
    #   5 on a mobilisable plasmid         -> transferable with a helper
    #   6 in an ICE, or on a conjugative plasmid -> PREDICTED self-transmissible
    # An IS sitting INSIDE an AMR gene is reported separately as likely
    # inactivation — it must not be counted as mobilisation.
    #
    # Takes in: the AMR calls, BOTH element sources (IS and ICE/IME — --is-table
    #           is repeatable), the replicon calls, and the contig lengths.
    # Produces: MOBILITY_TABLE (the deliverable, one row per AMR gene) and
    #           MOBILITY_AUDIT (why every gene without context got none).
    # Consumed by: the user. This is the module's terminal product.
    rule amr_mge_colocalisation:
        input:
            amrfinder = AMRFINDER_TSV,
            is_table = IS_TABLE,
            ice_table = ICE_TABLE,
            # Curated transposons/integrons, present only when a TnCentral source
            # was configured. Unpacking a dict keeps the input list valid either
            # way; without it tier 4 is simply never awarded.
            **({"named_table": NAMED_ELEMENTS_TABLE} if MOBILOME_NAME_ELEMENTS else {}),
            replicons = MOBILOME_REPLICONS,
            lengths = CONTIG_LENGTHS,
        output:
            table = MOBILITY_TABLE,
            audit = MOBILITY_AUDIT,
        params:
            script = COLOCALISE_SCRIPT,
            max_span = MOBILOME_MAX_COMPOSITE_SPAN,
            # Third element table, empty unless the TnCentral naming layer is on.
            # --is-table is `action="append"`, so extra tables just add elements.
            named_flag = lambda w: (
                "--is-table " + NAMED_ELEMENTS_TABLE.format(sample=w.sample)
                if MOBILOME_NAME_ELEMENTS else ""
            ),
        log:
            LOGS + "/mobilome_colocalisation_{sample}.log"
        priority: 3
        shell:
            """
            python {params.script} \
              --sample {wildcards.sample} \
              --amrfinder {input.amrfinder} \
              --is-table {input.is_table} \
              --is-table {input.ice_table} \
              {params.named_flag} \
              --replicons {input.replicons} \
              --contig-lengths {input.lengths} \
              --max-composite-span {params.max_span} \
              --out-table {output.table} \
              --out-audit {output.audit} > {log} 2>&1
            """
