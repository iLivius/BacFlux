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
# (WP-A, WP-C, WP-D below are the spec's "work package" labels — the order the
# module was built in, kept because docs/mobilome_module_SPEC.md refers to them:
# WP-A = surface the AMR calls, WP-C = find the IS elements, WP-D = intersect the
# two and award a mobility tier. WP-B is the database layer, WP-E the ICE work.)
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
#   bakta/{sample}.faa ──┬─► conjscan ──► conjugation table ───┤   {sample}_amr_mobility.tsv
#             (relaxase / T4CP / T4SS = can it self-transmit?) │   + _amr_mobility_audit.tsv
#                        └─► icescan ──► integrase + IME/AICE ─┤
#                           (OPTIONAL second model set; both   │
#                            feed conjscan_ice, which merges   │
#                            them — see mobilome.icescan.run)  │
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
# WHERE THE NAMES IN THIS FILE COME FROM. Every constant used below — the on/off
# switches (MOBILOME_RUN, MOBILOME_NAME_ELEMENTS, MOBILOME_NAME_ICE,
# MOBILOME_COPY_NUMBER), every output path, all eight mobilome script paths, the
# database URLs and thresholds, and the shared workflow names (FINAL_CONTIGS,
# DIR_ANNOTATION, GTDBTK_DIR, PLATON_DIR, LOGS, BAKTADB, RAM, TRIM_R1/TRIM_R2,
# PHAGE_CALLER, capped_cpus) — is defined once in the mobilome block of
# 00_common.smk. Read that block before this file; nothing here re-derives a path
# or re-reads the config, so producer and consumer can never drift apart. The
# only exceptions are declared just below: MOBILOME_COVERAGE_PROFILE and the four
# ICEscan constants, which are read here because nothing outside this file uses
# them.
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/x.yaml" climbs shared/ -> rules/ -> workflow/ -> workflow/envs/x.yaml.
# ─────────────────────────────────────────────────────────────────────────────

# ── How much of an HMM profile a protein must match (--coverage-profile) ─────
# Lives in this file rather than 00_common.smk for the same reason as the
# ICEscan constants below: it is read by exactly the two MacSyFinder rules in
# this file and by nothing else, so it reads best next to them.
#
# WHAT THE NUMBER IS. MacSyFinder finds machinery proteins with profile HMMs —
# statistical descriptions of what a relaxase, a coupling protein or a VirB4
# looks like across many species. A protein can score well against a profile
# while aligning to only part of it, which happens for two very different
# reasons: the protein really is a fragment or a decayed remnant, or it is a
# genuine full-length member of a DIVERGENT family that only shares the
# catalytic core. This threshold is the fraction of the PROFILE's length the
# alignment must cover before the hit is kept. It is not identity and it is not
# an E-value; a hit can be beyond doubt statistically and still be dropped here.
#
# WHY 0.5 IS THE DEFAULT. It is MacSyFinder's own default, and it is the value
# every scored result in docs/ was measured at — the ICE pilot, the IME pilot,
# the negative controls and the head-to-head against the EBI pipeline. Keeping
# it the default means the shipped configuration is the validated one.
#
# WHAT LOWERING IT BUYS AND COSTS. Measured end to end at 0.5 / 0.4 / 0.3 over
# the three benchmark sets — 18 curated ICEs, 12 curated IMEs, and 12 genomes
# with no curated element (32.6 Mb). With the ICEscan layer on:
#
#                curated ICEs   curated IMEs   calls on the negative set
#     0.5         15 of 18        5 of 12       5   (2 without ICEscan)
#     0.4         15 of 18        6 of 12       7   (4 without ICEscan)
#     0.3         15 of 18        6 of 12       9   (5 without ICEscan)
#
# One more curated element out of thirty, at double the calls on genomes that
# should have none. The one gained is a 23 kb IME in Faecalibacterium duncaniae,
# recovered in full and correctly classed — a real detection, not a scoring
# artefact — and it arrives at 0.4 with nothing further at 0.3.
#
# The extra negative-control calls at 0.4 are all `cime_or_island`, mobility
# `passive`, confidence `low`, evidence_level `profile_hits_only`: an integrase
# and a T4SS-like protein near each other, claiming nothing and carrying no
# mobility tier. The tiering absorbs them, which is its job. At 0.3 that stops
# holding — a 21.6 kb IME appears in Staphylococcus aureus N315 at medium
# confidence with an actual mobility claim, unverified either way.
#
# WHY THE GAIN IS SO SMALL WHEN THE RAW EVIDENCE MOVES A LOT. 0.4 adds 5% more
# MacSyFinder hits and changes the hit table in 13 of 40 benchmark genomes; 0.3
# adds 12% and changes 19 of 40. Almost none of it reaches the element table,
# because conjscan_ice needs an integrase within 50 kb of conjugation machinery
# before it seeds anything. Profile coverage sits upstream of a stronger
# constraint, and CO-LOCALISATION is what actually binds. Concretely: the
# pilot's one Streptococcus salivarius IME is still missed at 0.3, and its audit
# TSV explains it — the MOBT relaxases were already found at 0.5 and sit 22 kb
# and 518 kb from the integrase marking the element.
#
# The IME ceiling study says 73 of 395 curated IMEs carry a relaxase hit that is
# clean on E-value and fails only this rule (60 of them T4SS_MOBT at ~0.32
# coverage), so 0.4 is a reasonable thing to try on that biology — with the audit
# TSV as the check on what changed. Nothing measured supports 0.3.
# (73, not the 68 this comment carried until 2026-07-31: the 68 was stale in the
# script that builds the table, not in the table. 279 + 73 + 43 = 395. See
# docs/mobilome_tuning_guide.md §3.2.)
#
# The EBI mobilome-annotation-pipeline runs 0.3. That is evidence about a
# metagenome pipeline's priorities — recall over precision, because a MAG's
# proteins are fragmentary anyway — not about what a single-isolate workflow
# reporting a regulatory-facing mobility tier should do.
#
# BOTH SEARCHES GET THE SAME VALUE, deliberately. conjscan and icescan hit
# tables are UNIONED by the caller, so running them at different stringencies
# would mean an element's class depended on which model set happened to be more
# permissive, which nothing downstream could untangle.
MOBILOME_COVERAGE_PROFILE = float(
    (config.get("mobilome") or {}).get("coverage_profile", 0.5))

# A fraction outside (0, 1] is a typo (a percentage, most likely). MacSyFinder
# would accept 30 and then silently find nothing at all, which looks exactly
# like a genome with no conjugative system - so fail here instead.
if MOBILOME_RUN and not (0.0 < MOBILOME_COVERAGE_PROFILE <= 1.0):
    sys.exit(
        "[BacFlux] mobilome.coverage_profile must be a fraction greater than 0 "
        f"and at most 1.0, but it is {MOBILOME_COVERAGE_PROFILE}. It is the "
        "fraction of an HMM profile a protein must align to (0.5 = half), not a "
        "percentage."
    )

# ── The ICEscan model set: switch, paths, and why they live here ─────────────
# Every other mobilome path constant is declared in 00_common.smk. These four
# stay in this file on purpose, so that the whole ICEscan layer — the config
# switch, the two output paths and the two rules that write them — can be read
# in one place. Move them into 00_common.smk if you prefer the convention;
# nothing else refers to them.
#
# WHAT THIS LAYER IS. ICEscan is a SECOND MacSyFinder model package, run over the
# same Bakta proteins as CONJscan and merged with it. It is a fork of CONJScan by
# the same Institut Pasteur authors, one minor version behind the release BacFlux
# installs, so it is added ALONGSIDE and never instead of it (a swap loses the MOB
# relaxase models, the decayed-machinery models and the plasmid set). What it adds
# is the IME and AICE classes, which CONJScan 2.1.0 has no model for at all, plus
# integrase profiles. What it does NOT add is boundaries: MacSyFinder reports gene
# ORDINALS, not base pairs, so every coordinate still comes from our own Bakta
# GFF3 join, our own clustering and our own att search. The reasoning in full is
# in the config comments next to mobilome.icescan.
#
# OFF unless the user both switches it on AND gives a source, because turning it
# on downloads a CC BY-NC-SA (non-commercial) model set.
_icescan_cfg = ((config.get("mobilome") or {}).get("icescan") or {})
ICESCAN_ENABLED = _config_bool(_icescan_cfg.get("run"), False)
ICESCAN_URL = str(_icescan_cfg.get("url") or "").strip()
ICESCAN_SHA256 = str(_icescan_cfg.get("sha256") or "").strip()
# A models directory the user already holds (one CONTAINING an ICEscan/ folder)
# beats downloading, exactly as directories.* beats links.* everywhere else.
ICESCAN_LOCAL = str(_icescan_cfg.get("dir") or "").strip()

# Switched on but with nowhere to fetch from is a configuration mistake, not a
# reason to quietly run without it: the user would get today's results back and
# no indication that the layer they asked for never existed.
if MOBILOME_RUN and ICESCAN_ENABLED and not (ICESCAN_URL or ICESCAN_LOCAL):
    sys.exit(
        "[BacFlux] mobilome.icescan.run is true but neither mobilome.icescan.url "
        "nor mobilome.icescan.dir is set, so there is no model package to use. "
        "Set the url (see config/config_v2.yaml for the ICEfinder2 bundle it "
        "comes in), point dir at a MacSyFinder models directory that already "
        "holds ICEscan/, or set run: false."
    )

MOBILOME_ICESCAN = MOBILOME_RUN and ICESCAN_ENABLED

# One shared copy of the models for all samples, next to conjscan_models.
ICESCAN_MODELS_DIR = DIR_MOBILOME + "/icescan_models"   # a DIRECTORY (rule icescan_models)
ICESCAN_DIR = MOBILOME_DIR + "/icescan"                 # a DIRECTORY (rule icescan)

# ── What to expect on a DRAFT assembly (printed once, at parse time) ─────────
# The same pattern as every other banner in this workflow: a plain print() next
# to the constants it describes, so it lands in the log header before the DAG is
# built (see the MODE and phage-caller banners in 00_common.smk).
#
# WHY IT EXISTS. Every validation this module had was on CLOSED genomes, where
# spans_contigs was TRUE on 0 of 63 calls - so its fragmentation guards had
# never been measured on the input BacFlux actually gets. They have now been:
# 40 benchmark genomes were cut to ~150 kb, ~50 kb and ~20 kb N50, all 120
# assemblies re-annotated and re-run end to end. The numbers below are from that
# run, not from an estimate. Kept to one short paragraph on purpose - a wall of
# text at every run gets skipped, and this one has to be read.
if MOBILOME_RUN:

    print(
        "ICE/IME calling was validated on CLOSED genomes; measured on drafts it "
        "degrades honestly. At ~50 kb N50 the CLASS still holds (ICE 15/18) but "
        "the EXTENT does not (median 0.34x the true length) and high-confidence "
        "calls fall from 24% to 10%. So read mge_class + confidence next to "
        "spans_contigs and at_contig_boundary, and read boundary_method before "
        "start/end: 'none' means the interval is only the machinery span, a "
        "floor. Low-confidence calls on a draft are expected, not a fault. Each "
        "sample's contig count and N50 head its _ice_discarded.tsv; details in "
        "docs/mobilome_draft_assemblies.md."
    )

    # ── Rule: contig_lengths — how long is every contig? ─────────────────────
    # Biology: nothing on its own — but every downstream honesty check needs it.
    # "Is this AMR gene 200 bp from the end of its contig?" is only answerable
    # against the contig's length, and that question is what separates "there is
    # no IS next to this gene" from "the contig ended before we could tell".
    #
    # Takes in: FINAL_CONTIGS — the single delivered genome, contigs_final.fasta.
    #           Every mode's front end ends by writing that one file (the
    #           workflow calls this the "D2 hand-off"), so no rule from stage 03
    #           onwards has to know which assembler produced it.
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
    #
    # --coverage-profile is written out even though the default value BacFlux
    # ships (0.5) is also MacSyFinder's own default. Two reasons: the number is
    # now a config key a user may change, and it must be identical here and in
    # rule icescan because the two hit tables get merged. Leaving it implicit
    # would mean a MacSyFinder release could move it under us and only one of the
    # two searches would notice. See MOBILOME_COVERAGE_PROFILE at the top of this
    # file for what lowering it buys and costs.
    rule conjscan:
        input:
            bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
            models = CONJSCAN_MODELS_DIR,
        output:
            conjscan_dir = directory(CONJSCAN_DIR),
        params:
            coverage = MOBILOME_COVERAGE_PROFILE,
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
              --coverage-profile {params.coverage} \
              --out-dir {output.conjscan_dir} \
              --worker {threads} \
              --force > {log} 2>&1 || {{
                echo "" >> {log}
                echo "NOTE: MacSyFinder exited non-zero. A genome with no conjugative system is the common case for environmental isolates; the module continues and the ICE/IME table will be empty." >> {log}
              }}
            """

    # ══ The ICEscan model set (a second machinery search) ════════════════════
    # Only exists when the user opted in. Without it the module behaves exactly
    # as it did before: CONJscan alone, no IME class, no AICE class.
    if MOBILOME_ICESCAN:

        # ── Rule: icescan_models — fetch the ICEscan model package once ──────
        # Biology: the same kind of package as CONJscan — profile models of the
        # proteins an element needs in order to move — but covering two classes
        # CONJScan 2.1.0 does not model at all:
        #   IME  = integrates into the chromosome and carries a relaxase (the
        #          enzyme that nicks the DNA to start transfer) but NO mating
        #          apparatus, so it moves only by borrowing one from a
        #          conjugative element in the same cell (mobility tier 5);
        #   AICE = the actinomycete elements that do not conjugate at all and
        #          instead push single-stranded DNA between hyphae with an
        #          FtsK/SpoIIIE translocase.
        # It also carries integrase profiles CONJscan lacks. Only some of those
        # are trusted as element integrases — that judgement lives in the caller
        # (conjscan_to_ice.py), not here.
        #
        # ⚠ LICENCE (from the package's own LICENSE file): CC BY-NC-SA 4.0,
        # Institut Pasteur / CNRS — academic and non-commercial use only, the
        # same terms as the CONJscan models. BacFlux never vendors either: they
        # are fetched here, at the user's request, under the user's own agreement
        # with the licensor. BacFlux's own MIT licence is unaffected.
        #
        # Takes in: nothing from the workflow — a URL from the config, or a
        #           models directory the user already holds.
        # Does: fetch ICEfinder2's database bundle and extract ONLY the
        #       macsydata/ICEscan package from it.
        # Produces: 08.mobilome/icescan_models/ICEscan/ — one shared copy for all
        #           samples — plus PROVENANCE.txt.
        # Consumed by: icescan.
        #
        # Why not `msf_data install` as conjscan_models does: ICEscan is not
        # published in the macsy-models registry that command reads. It is
        # distributed inside ICEfinder2's database bundle, so it has to be
        # downloaded and unpacked by hand.
        #
        # conda: NONE — wget and tar from the launch environment, the same
        # arrangement as rule download_amr_db in 50_amr.smk.
        rule icescan_models:
            output:
                models = directory(ICESCAN_MODELS_DIR),
            params:
                url = ICESCAN_URL,
                sha256 = ICESCAN_SHA256,
                local_dir = ICESCAN_LOCAL,
            log:
                LOGS + "/mobilome_icescan_models.log"
            priority: 4
            shell:
                """
                exec > {log} 2>&1
                set -euo pipefail
                mkdir -p {output.models}

                echo "Installing the ICEscan model package."
                echo "NOTE: these models are licensed CC BY-NC-SA 4.0 (Institut Pasteur/CNRS) - academic / non-commercial use only. They are downloaded here at your request and are never redistributed by BacFlux."

                if [ -n "{params.local_dir}" ]; then
                    echo "Using the local models directory '{params.local_dir}'. Nothing will be downloaded."
                    cp -r "{params.local_dir}"/ICEscan {output.models}/ICEscan
                    SOURCE="local:{params.local_dir}"
                else
                    echo "Fetching ICEfinder2's database bundle from {params.url}"
                    # Use the https:// spelling of this path: the ftp:// one
                    # times out from behind many firewalls, including ours.
                    # -nv (not -q) so that a failed transfer says why in the log
                    # rather than leaving a bare non-zero exit to explain itself.
                    wget -nv -O {output.models}/icf2_dbs.tar.gz "{params.url}"

                    # A redirect to an HTML error page arrives with a 200 and
                    # would leave an empty models directory behind. MacSyFinder
                    # would then run happily and report no systems at all, which
                    # is a WRONG answer that looks exactly like a real one, so
                    # fail here instead.
                    if ! tar -tzf {output.models}/icf2_dbs.tar.gz > /dev/null 2>&1; then
                        echo "ERROR: the download is not a gzipped tar archive. The server may" >&2
                        echo "       have moved the file or returned an error page. First bytes:" >&2
                        head -c 200 {output.models}/icf2_dbs.tar.gz >&2
                        exit 1
                    fi

                    OBSERVED=$(sha256sum {output.models}/icf2_dbs.tar.gz | cut -d' ' -f1)
                    if [ -n "{params.sha256}" ] && [ "$OBSERVED" != "{params.sha256}" ]; then
                        echo "ERROR: ICEscan bundle checksum mismatch." >&2
                        echo "       expected {params.sha256}" >&2
                        echo "       observed $OBSERVED" >&2
                        echo "       The URL carries no version, so this means the upstream file" >&2
                        echo "       changed. Update mobilome.icescan.sha256 once you have decided" >&2
                        echo "       the new release is the one you want." >&2
                        exit 1
                    fi

                    # The bundle also ships ICEfinder2's own HMM sets and a
                    # UniProt BLAST index - about 60 MB this module never reads.
                    # Extract only the MacSyFinder package. --strip-components=2
                    # removes the leading 'icf2_dbs/macsydata/', so the package
                    # lands as <models>/ICEscan, which is the layout
                    # --models-dir expects.
                    tar -xzf {output.models}/icf2_dbs.tar.gz -C {output.models} \
                      --strip-components=2 icf2_dbs/macsydata/ICEscan
                    rm -f {output.models}/icf2_dbs.tar.gz
                    SOURCE="{params.url} (sha256 $OBSERVED)"
                fi

                # A MacSyFinder package is definitions (which genes make a
                # system) plus profiles (the HMMs that find those genes). With
                # either missing the search still runs and still finds nothing,
                # so check for both rather than trust the archive's shape.
                if [ ! -s {output.models}/ICEscan/definitions/Chromosome/IME.xml ]; then
                    echo "ERROR: the ICEscan package has no definitions/Chromosome/IME.xml." >&2
                    echo "       Either the archive layout changed or the local directory does" >&2
                    echo "       not hold a MacSyFinder model package." >&2
                    exit 1
                fi
                # The file existing is not enough. ICEscan ships TWO IME
                # definitions and they disagree about which relaxases count:
                #
                #   definitions/Chromosome/IME.xml  18 relaxase families,
                #       including the eight Relaxase_* profiles ICEscan adds on
                #       top of CONJScan (the Gram-positive and IME-specific ones);
                #   IME_type.xml, at the package root   11 families, and NOT ONE
                #       of those eight.
                #
                # MacSyFinder only ever reads the definitions/ tree - it lists
                # <package>/definitions and recurses from there - so the root-level
                # file is dead weight and we load the wide one. Two things confirm
                # that rather than assume it: every model_fqn MacSyFinder writes is
                # 'ICEscan/Chromosome/...', never 'ICEscan/IME_type', and the root
                # AICE_type.xml names HMMs (AICE_rep1, AICE_tra) that the package
                # does not even ship, so loading it would fail outright.
                #
                # We still check, because the download URL carries no version.
                # Measured over the 395 curated ICEberg IMEs, the narrow set would
                # cost 41 elements (10.4%) that the CONJScan leg does not rescue -
                # mostly Streptococcus salivarius, found through
                # Relaxase_firmi_Rep_2 - and would gain nothing, since its one
                # exclusive family (T4SS_MOBL) matches none of the 395. That loss
                # would surface as a thinner IME table, not as an error, so assert
                # the shape of the definition instead of trusting the archive.
                for PROFILE in Relaxase_firmi_Rep_2 Relaxase_PHA_IME_A1 Relaxase_profile_MOBT; do
                    if ! grep -q "$PROFILE" {output.models}/ICEscan/definitions/Chromosome/IME.xml; then
                        echo "ERROR: definitions/Chromosome/IME.xml does not list $PROFILE." >&2
                        echo "       This looks like the NARROW IME definition (the 11-family" >&2
                        echo "       one ICEscan also ships as IME_type.xml). Running it would" >&2
                        echo "       silently drop roughly 10% of detectable IMEs, mostly in" >&2
                        echo "       Gram-positives, with no other sign that anything changed." >&2
                        echo "       The upstream package has probably been reorganised; check" >&2
                        echo "       it before updating mobilome.icescan.sha256." >&2
                        exit 1
                    fi
                done

                if [ ! -d {output.models}/ICEscan/profiles ]; then
                    echo "ERROR: the ICEscan package has no profiles/ directory." >&2
                    exit 1
                fi
                N_PROFILES=$(find {output.models}/ICEscan/profiles -name '*.hmm' | wc -l)
                if [ "$N_PROFILES" -lt 1 ]; then
                    echo "ERROR: no HMM profiles in {output.models}/ICEscan/profiles." >&2
                    exit 1
                fi

                # metadata.yml carries the package's own version line, which is
                # the only version string in the whole download.
                VERSION=$(grep '^vers:' {output.models}/ICEscan/metadata.yml | cut -d' ' -f2)
                {{
                  echo "source:       $SOURCE"
                  echo "fetched:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
                  echo "package_vers: $VERSION"
                  echo "profiles:     $N_PROFILES"
                  echo ""
                  echo "The ICEscan download URL is UNVERSIONED, so this file is the only"
                  echo "record of which release was used. Quote it in a methods section."
                  echo ""
                  echo "ICEscan is a fork of CONJScan by the same Institut Pasteur authors,"
                  echo "distributed inside ICEfinder2's database bundle and licensed"
                  echo "CC BY-NC-SA 4.0 - academic / non-commercial use only. BacFlux ships"
                  echo "no models, only this URL. Cite Coluzzi et al. 2022 and Wang et al."
                  echo "2024 (ICEfinder2)."
                }} > {output.models}/PROVENANCE.txt

                echo "ICEscan ready: version $VERSION, $N_PROFILES profiles."
                """

        # ── Rule: icescan — the second machinery search on THIS genome ───────
        # Takes in: the same Bakta proteins (.faa) rule conjscan reads, and the
        #           ICEscan models above.
        # Does: MacSyFinder again, with the ICEscan package. `--models ICEscan
        #       all` (not ICEscan/Chromosome, the form conjscan uses) because
        #       ICEscan ships only a Chromosome set - the plasmid models were
        #       removed in the fork - so `all` is the whole package.
        # Produces: 08.mobilome/{sample}/icescan/ (best_solution.tsv and friends).
        # Consumed by: conjscan_ice, which UNIONS this table with CONJscan's.
        #
        # --coverage-profile comes from the SAME config key rule conjscan reads
        # (mobilome.coverage_profile, default 0.5). That is not tidiness: the two
        # searches are merged by the caller, so running them at different
        # stringencies would make an element's class depend on which model set was
        # more permissive. See MOBILOME_COVERAGE_PROFILE at the top of this file.
        #
        # Same tolerance of a non-zero exit as rule conjscan: a genome with no
        # detectable IME or AICE is the ordinary result, not an error.
        rule icescan:
            input:
                bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
                models = ICESCAN_MODELS_DIR,
            output:
                icescan_dir = directory(ICESCAN_DIR),
            params:
                coverage = MOBILOME_COVERAGE_PROFILE,
            conda:
                "../../envs/macsyfinder.yaml"
            threads: capped_cpus(8)
            log:
                LOGS + "/mobilome_icescan_{sample}.log"
            priority: 3
            shell:
                """
                rm -rf {output.icescan_dir}
                mkdir -p {output.icescan_dir}

                macsyfinder \
                  --models ICEscan all \
                  --sequence-db {input.bakta_dir}/{wildcards.sample}.faa \
                  --db-type ordered_replicon \
                  --models-dir {input.models} \
                  --coverage-profile {params.coverage} \
                  --out-dir {output.icescan_dir} \
                  --worker {threads} \
                  --force > {log} 2>&1 || {{
                    echo "" >> {log}
                    echo "NOTE: MacSyFinder exited non-zero on the ICEscan model set. A genome with no integrative element is the common case; the module continues, and the CONJscan search is unaffected." >> {log}
                  }}
                """

    # ── Rule: conjscan_ice — turn machinery hits into ICE / IME candidates ───
    # Takes in, and which rule produced each:
    #   * CONJscan's best_solution.tsv and its hmmer_results/ directory  (conjscan)
    #   * ICEscan's best_solution.tsv, ONLY when that layer is on        (icescan)
    #   * Bakta's GFF3 — the genomic coordinates of every protein hit, the
    #     integrase genes found by product regex, and the tRNAs the att search
    #     anchors on                                                    (annotation)
    #   * the genome itself, for the att-site sequence search            (D2 hand-off)
    #   * the IS table, masked out of the att search (see below)         (isescan_table)
    #   * the contig lengths, for the contig-edge flags                  (contig_lengths)
    #   * chromosome-or-plasmid, per contig                       (mobilome_replicons)
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
    # The IS table is an input because every IS interval is MASKED OUT of the
    # flanking sequence before the repeat search runs — all of it, not one half:
    # insertion sequences carry terminal repeats and duplicate their target DNA
    # when they transpose, so an IS-rich neighbourhood is full of direct repeats
    # that have nothing to do with ICE integration. The spec names this as the
    # most likely way to get Phase 3 wrong.
    #
    # Produces: the ICE/IME element table + its audit.
    # Consumed by: amr_mge_colocalisation, as a SECOND element source alongside
    #              the IS table.
    rule conjscan_ice:
        input:
            conjscan_dir = CONJSCAN_DIR,
            # The optional second machinery search. Present only when the user
            # opted in to the ICEscan model set; unpacking a dict keeps the input
            # list valid either way. The script UNIONS the two hit tables and
            # takes from ICEscan only its integrase hits and its IME/AICE
            # classes - never its spans (MacSyFinder reports gene ordinals, not
            # base pairs) and never its quorum for the mating apparatus.
            **({"icescan_dir": ICESCAN_DIR} if MOBILOME_ICESCAN else {}),
            bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
            lengths = CONTIG_LENGTHS,
            genome = FINAL_CONTIGS,
            is_table = IS_TABLE,
            # Which contigs are plasmids. An ICE is by definition an element that
            # integrates into the CHROMOSOME (spec section 2.4), so conjugation
            # machinery sitting on a plasmid is a conjugative plasmid, not an ICE.
            # Without this the classifier cannot tell the two apart and reports
            # every self-transmissible plasmid as a predicted ICE.
            replicons = MOBILOME_REPLICONS,
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
            # Empty unless the ICEscan layer is on, so a user who left it off
            # gets exactly the behaviour this rule had before it existed.
            #
            # Only --icescan-tsv is passed, with no --icescan-hmmer-dir beside
            # it. That is not an oversight: when the flag is omitted the script
            # looks for hmmer_results/ next to the best_solution.tsv it was
            # given, which is exactly the path we would have written. The
            # CONJscan side passes both only because it was wired first.
            icescan_flag = lambda w: (
                "--icescan-tsv " + ICESCAN_DIR.format(sample=w.sample) + "/best_solution.tsv"
                if MOBILOME_ICESCAN else ""
            ),
        log:
            LOGS + "/mobilome_conjscan_ice_{sample}.log"
        priority: 3
        shell:
            # best_solution.tsv is absent when MacSyFinder found nothing; the
            # script treats a missing file as "no systems" and still writes a
            # well-formed empty table, so no guard is needed here - and that
            # holds for the ICEscan table too.
            """
            python {params.script} \
              --sample {wildcards.sample} \
              --conjscan-tsv {input.conjscan_dir}/best_solution.tsv \
              --conjscan-hmmer-dir {input.conjscan_dir}/hmmer_results \
              {params.icescan_flag} \
              --bakta-gff {input.bakta_dir}/{wildcards.sample}.gff3 \
              --contig-lengths {input.lengths} \
              --genome {input.genome} \
              --is-table {input.is_table} \
              --replicons {input.replicons} \
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
        # Produces: the directory TNCENTRAL_DB_DIR, holding tncentral.fa, the
        #           BLAST index files under the prefix tncentral_v5, and
        #           PROVENANCE.txt. Those filenames are written literally in the
        #           shell block below because they are shell-side paths; the only
        #           one another rule needs is the tncentral_v5 prefix, which
        #           tncentral_blast rebuilds from input.db_dir.
        # Consumed by: tncentral_blast. Fetched once and shared by all samples.
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

                # REPAIR MALFORMED RECORDS BEFORE INDEXING. The upstream FASTA has
                # deflines glued onto the END of a sequence line instead of
                # starting their own, e.g.
                #     ...gtgcagccgtcttctgaaaacgaca>In1223-KX784502
                # In the release checked (2025-05-16) 21 of 533 records were like
                # this, and the damage runs BOTH ways:
                #   * those 21 elements are invisible to makeblastdb - among them
                #     Tn7 itself and eleven integrons, the very class tier 4 exists
                #     to name;
                #   * and the 16 records they were glued to became CHIMERIC,
                #     absorbing the defline text plus the next element's sequence.
                #     In_Tn6162 measured 41,492 bp instead of 8,911 - 4.7x its real
                #     length. Since the naming cascade tests coverage as
                #     alignment/slen, an inflated slen makes those elements almost
                #     impossible to name, silently.
                # Splitting on '>' is safe here because these deflines carry no
                # free-text description in which a '>' could legitimately appear.
                awk '{{ if (substr($0,1,1) != ">") gsub(/>/, "\\n>"); print }}' \
                  {output.db_dir}/tncentral.fa > {output.db_dir}/tncentral.repaired.fa
                mv {output.db_dir}/tncentral.repaired.fa {output.db_dir}/tncentral.fa

                # Every '>' must now begin a line. If not, the file has a shape we
                # did not anticipate and indexing it would silently lose or merge
                # records - fail instead of producing a quietly wrong database.
                N_SEQ=$(grep -c '^>' {output.db_dir}/tncentral.fa)
                N_MARK=$(grep -o '>' {output.db_dir}/tncentral.fa | wc -l)
                if [ "$N_SEQ" != "$N_MARK" ]; then
                    echo "ERROR: $N_MARK '>' characters but only $N_SEQ deflines start a line." >&2
                    echo "       The FASTA is malformed in a way this rule does not handle." >&2
                    exit 1
                fi

                # The archive ships a BLAST v4 index. Rebuild as v5 so it works
                # with current blast+ (spec §5.3), and drop the shipped index files
                # so there is no chance of the old one being picked up instead.
                # NOTE the shipped index was built from the UNREPAIRED FASTA and is
                # therefore missing those 21 elements - another reason to rebuild.
                rm -f {output.db_dir}/tncentral.fa.n*
                makeblastdb -in {output.db_dir}/tncentral.fa -dbtype nucl \
                  -out {output.db_dir}/tncentral_v5 -blastdb_version 5
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
        # Biology: a curated transposon or integron is a named architecture
        # somebody characterised and deposited. Finding one on this genome is a
        # much stronger statement than inferring a composite from two IS copies,
        # which is why a hit here is what makes ladder tier 4 reachable.
        #
        # Takes in: this sample's contigs (the D2 hand-off) + the v5 BLAST
        #           database from rule tncentral_db above.
        # Does: blastn the whole assembly against the curated set. The whole
        #       assembly rather than the ISEScan calls, because a transposon is
        #       bigger than the IS that bounds it and would be clipped otherwise.
        # Produces: TNCENTRAL_BLAST_HITS, tabular, with the exact columns
        #           name_transposons.py expects (BLAST_COLUMNS there).
        # Consumed by: name_elements, which applies the naming thresholds. This
        #           rule deliberately applies none of them — see the evalue note.
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
        # FOUR NAMES, ONE STEP, so nobody has to guess: the rule is
        # `name_elements`, the script it runs is `name_transposons.py` (reached
        # through the constant NAME_TRANSPOSONS_SCRIPT), the table it writes is
        # NAMED_ELEMENTS_TABLE, and its log is mobilome_name_elements_*.log. The
        # names have not been unified because renaming the rule or the script
        # would change the job names and the log filenames on disk for no gain.
        #
        # Biology: a BLAST hit is not yet an element. One curated transposon
        # usually matches as several high-scoring segments, and the same
        # transposon may be present in two separate copies on the genome — those
        # must stay two elements, not become one enormous one.
        #
        # Takes in: TNCENTRAL_BLAST_HITS from tncentral_blast above.
        # Does: merge HSPs into element COPIES (separate copies of one transposon
        #       must not be joined - see the script), apply the naming thresholds,
        #       drop plain IS entries that ISEScan already covers, and write the
        #       result in the element-table shape colocalise.py consumes.
        # Produces: NAMED_ELEMENTS_TABLE + its discard audit, which records a
        #       reason for every hit refused a name (identity, coverage, or being
        #       a plain IS).
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

    # ══ The ICEberg naming layer (which ICE is it?) ══════════════════════════
    # Only exists when the user configured an ICEberg source. This layer adds NO
    # elements and can change NO tier: conjscan_to_ice.py decides what is an ICE,
    # and this only says which one.
    if MOBILOME_NAME_ICE:

        # ── Rule: iceberg_db — fetch and index the curated ICE catalogue ─────
        # Takes in: nothing from the workflow — URLs from the config, or a
        #           directory the user already holds.
        # Does: fetch each .fas, concatenate, index as a v5 BLAST database, and
        #       record provenance (ICEberg is versioned - 3.0, June 2023 - but the
        #       download URLs are not, so the fetch date is worth keeping).
        # Produces: the directory ICEBERG_DB_DIR, holding iceberg.fa, the BLAST
        #           index under the prefix iceberg_v5, and PROVENANCE.txt — the
        #           same arrangement as tncentral_db, and for the same reason:
        #           those are shell-side filenames, and iceberg_blast rebuilds the
        #           one prefix it needs from input.db_dir.
        # Consumed by: iceberg_blast. Fetched once and shared by all samples.
        #
        # NO CHECKSUM KEY, unlike the ICEscan and TnCentral layers. Not an
        # omission: ICEberg serves two files rather than one archive, and the
        # larger is ~100 MB over a slow link that the rule has to be allowed to
        # resume — so a single pinned digest is not the right instrument here.
        # PROVENANCE.txt records the observed sha256 of the concatenated FASTA
        # and the fetch date, which is what a methods section needs.
        rule iceberg_db:
            output:
                db_dir = directory(ICEBERG_DB_DIR),
            params:
                urls = " ".join(ICEBERG_URLS),
                local_dir = ICEBERG_LOCAL,
                user_agent = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                              "(KHTML, like Gecko) Chrome/120.0 Safari/537.36"),
            conda:
                "../../envs/tncentral.yaml"
            log:
                LOGS + "/mobilome_iceberg_db.log"
            priority: 4
            shell:
                """
                exec > {log} 2>&1
                set -euo pipefail
                mkdir -p {output.db_dir}

                if [ -n "{params.local_dir}" ]; then
                    echo "Using the local ICEberg directory '{params.local_dir}'. Nothing will be downloaded."
                    cat "{params.local_dir}"/*.fas > {output.db_dir}/iceberg.fa
                    SOURCE="local:{params.local_dir}"
                else
                    : > {output.db_dir}/iceberg.fa
                    for url in {params.urls}; do
                        echo "Fetching $url"
                        # ICE_seq_all.fas is ~100 MB and the server is slow, so allow
                        # resuming rather than restarting a part-finished transfer.
                        curl -sSL -C - -A "{params.user_agent}" -o {output.db_dir}/part.fas "$url"
                        # Reject an HTML error page, which would otherwise be
                        # indexed as an empty database that silently names nothing.
                        # Test that a defline appears near the START of the file
                        # rather than that byte 0 is '>': the real ICE_seq_all.fas
                        # begins with a newline, and the stricter test rejected it,
                        # which meant this layer could never build at all.
                        if ! head -c 4096 {output.db_dir}/part.fas | grep -q '^>'; then
                            echo "ERROR: $url did not return FASTA. First bytes:" >&2
                            head -c 200 {output.db_dir}/part.fas >&2
                            exit 1
                        fi
                        cat {output.db_dir}/part.fas >> {output.db_dir}/iceberg.fa
                        rm -f {output.db_dir}/part.fas
                    done
                    SOURCE="{params.urls}"
                fi

                # Same defect as TnCentral, and repaired the same way: the ICEberg
                # release checked (2023-06-01) had one defline glued onto the end
                # of a sequence line, losing an IME and leaving the record before
                # it chimeric. One in 1,774 is rarer than TnCentral's 21 in 533,
                # but a silently merged reference is exactly as wrong.
                awk '{{ if (substr($0,1,1) != ">") gsub(/>/, "\\n>"); print }}' \
                  {output.db_dir}/iceberg.fa > {output.db_dir}/iceberg.repaired.fa
                mv {output.db_dir}/iceberg.repaired.fa {output.db_dir}/iceberg.fa

                N_SEQ=$(grep -c '^>' {output.db_dir}/iceberg.fa)
                N_MARK=$(grep -o '>' {output.db_dir}/iceberg.fa | wc -l)
                if [ "$N_SEQ" != "$N_MARK" ]; then
                    echo "ERROR: $N_MARK '>' characters but only $N_SEQ deflines start a line." >&2
                    exit 1
                fi

                makeblastdb -in {output.db_dir}/iceberg.fa -dbtype nucl \
                  -out {output.db_dir}/iceberg_v5 -blastdb_version 5
                {{
                  echo "source:      $SOURCE"
                  echo "fetched:     $(date -u +%Y-%m-%dT%H:%M:%SZ)"
                  echo "fasta_sha256: $(sha256sum {output.db_dir}/iceberg.fa | cut -d' ' -f1)"
                  echo "sequences:   $N_SEQ"
                  echo ""
                  echo "ICEberg 3.0 (released June 2023) publishes no licence or terms of"
                  echo "use; its pages carry only 'Copyright (c) 2023 All Rights Reserved by"
                  echo "Microbial Bioinformatics Group in MML, SJTU.'  BacFlux ships no"
                  echo "ICEberg data, only these URLs. Cite Wang et al. 2024, NAR."
                }} > {output.db_dir}/PROVENANCE.txt

                echo "ICEberg ready: $N_SEQ sequences."
                """

        # ── Rule: iceberg_blast — where do curated ICEs match this genome? ───
        # Biology: ICEberg curates ICEs and IMEs that have been described and
        # named in the literature. Matching one turns "a predicted
        # self-transmissible element" into a name you can look up.
        #
        # Takes in: this sample's contigs (the D2 hand-off) + the v5 database from
        #           rule iceberg_db above.
        # Does: blastn the WHOLE genome, not our ICE intervals, so a curated
        #       element that OVERHANGS our interval still shows up — which is how
        #       the naming step can report that our boundaries fall short.
        # Produces: ICEBERG_BLAST_HITS, the same tabular column set as the
        #       TnCentral search, which is what name_ice_elements parses.
        # Consumed by: name_ice_elements.
        rule iceberg_blast:
            input:
                contigs = FINAL_CONTIGS,
                db_dir = ICEBERG_DB_DIR,
            output:
                hits = ICEBERG_BLAST_HITS,
            params:
                db = lambda w, input: os.path.join(input.db_dir, "iceberg_v5"),
                evalue = "1e-50",
                # ICEs of one species are near-identical across strains, so a
                # single element matches many entries. Keep enough of them that
                # the naming step can report how ambiguous the name really is.
                max_targets = 50,
            conda:
                "../../envs/tncentral.yaml"
            threads: capped_cpus(8)
            log:
                LOGS + "/mobilome_iceberg_blast_{sample}.log"
            priority: 3
            shell:
                """
                blastn \
                  -query {input.contigs} \
                  -db {params.db} \
                  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore slen qlen" \
                  -evalue {params.evalue} \
                  -max_target_seqs {params.max_targets} \
                  -num_threads {threads} \
                  -out {output.hits} > {log} 2>&1
                """

        # ── Rule: name_ice_elements — put the curated name on the candidate ──
        # Biology: this LABELS, it never DECIDES. conjscan_to_ice.py has already
        # said what is an ICE and what class it is; all this step adds is which
        # published element it looks like, and how much of that element our
        # interval actually covers. It can change no gene's mobility tier.
        #
        # Takes in: the ICE/IME candidates (rule conjscan_ice) and the ICEberg
        #           BLAST hits (rule iceberg_blast).
        # Does: match each candidate interval to overlapping curated elements,
        #       above the identity and overlap floors from the config. A name is
        #       suffixed "-like" when less than 80% of the reference is present.
        # Produces: ICE_TABLE_NAMED — the same rows and the same columns as
        #       ICE_TABLE, with mge_name filled in — plus a naming audit that
        #       records what each candidate matched and what was refused.
        # Consumed by: amr_mge_colocalisation, which reads this table INSTEAD of
        #       ICE_TABLE whenever this layer is on. That swap is decided once, in
        #       00_common.smk, by ICE_TABLE_FOR_COLOCALISE.
        #
        # Its log is mobilome_name_ice_{sample}.log — shortened, where every other
        # rule's log is "mobilome_" plus the rule name in full. Left alone because
        # renaming it changes filenames on disk for no benefit.
        rule name_ice_elements:
            input:
                ice_table = ICE_TABLE,
                hits = ICEBERG_BLAST_HITS,
            output:
                table = ICE_TABLE_NAMED,
                audit = ICE_NAMING_AUDIT,
            params:
                script = NAME_ICE_SCRIPT,
                min_identity = ICEBERG_MIN_IDENTITY,
                min_overlap = ICEBERG_MIN_OVERLAP,
            log:
                LOGS + "/mobilome_name_ice_{sample}.log"
            priority: 3
            shell:
                """
                python {params.script} \
                  --sample {wildcards.sample} \
                  --ice-table {input.ice_table} \
                  --blast {input.hits} \
                  --min-identity {params.min_identity} \
                  --min-overlap-fraction {params.min_overlap} \
                  --out-table {output.table} \
                  --out-audit {output.audit} > {log} 2>&1
                """

    # ══ How many IS copies did the assembly lose? (spec WP-C) ═══════════════
    # Short-read modes only - it needs reads - and only when an ISOSDB source is
    # configured. Changes no AMR gene's tier: this is a quality metric on the IS
    # inventory, quantifying the collapse the module warns about everywhere else.
    #
    # THIS LEG NEEDS TWO FILES, and the switch that turns it on only looks for
    # one. MOBILOME_COPY_NUMBER (00_common.smk) is satisfied by isosdb.fasta_url
    # alone, but rule isosdb_db also downloads the IS family map. With the family
    # URL left empty the download runs `curl -o ... ""`, which dies with a bare
    # "curl: (3) URL using bad/illegal format" three quarters of the way through a
    # run and never reaches the friendly message written for a missing family map.
    # So say it here, at parse time, in the same way and for the same reason as
    # the ICEscan check near the top of this file.
    if MOBILOME_COPY_NUMBER and not ISOSDB_LOCAL and not ISOSDB_FAMILY_URL:
        sys.exit(
            "[BacFlux] mobilome.isosdb.fasta_url is set but "
            "mobilome.isosdb.family_map_url is empty, so the IS family map "
            "cannot be downloaded and every element would be reported as "
            "'unassigned'. Set family_map_url (see config/config_v2.yaml for the "
            "pseudoR link it sits beside), or point mobilome.isosdb.dir at a "
            "directory already holding BOTH ISOSDB.V3.fna and IS_fam_annot.txt, "
            "or clear fasta_url to leave the IS copy-number leg switched off."
        )

    if MOBILOME_COPY_NUMBER:

        # ── Rule: isosdb_db — fetch the openly licensed IS sequence set ──────
        # Biology: ISOSDB is a catalogue of IS nucleotide sequences, dereplicated
        # at 95% identity, plus a map from each entry to its IS family. Reads are
        # mapped against it in the next rules; the family map is what lets the
        # answer be reported per FAMILY, which is the robust number.
        #
        # ISOSDB lives in the pseudoR repository under the MIT licence, which
        # makes it the one mobilome database with no redistribution question. It
        # is still fetched rather than vendored, to keep one rule for all of them.
        # For the same reason it carries no sha256 config key: the URLs are pinned
        # GitHub raw paths under version control at the source, so PROVENANCE.txt
        # recording the observed checksum and fetch date is enough.
        #
        # Takes in: nothing from the workflow — two URLs from the config, or a
        #           directory the user already holds.
        # Does: download and unzip the sequence set, fetch the family map beside
        #       it, and refuse to continue if the family map did not arrive.
        # Produces: the directory ISOSDB_DB_DIR, holding ISOSDB.V3.fna,
        #           IS_fam_annot.txt and PROVENANCE.txt.
        # Consumed by: isosdb_map (the sequences) and is_copy_number (the family
        #           map). Fetched once and shared by all samples.
        rule isosdb_db:
            output:
                db_dir = directory(ISOSDB_DB_DIR),
            params:
                fasta_url = ISOSDB_FASTA_URL,
                family_url = ISOSDB_FAMILY_URL,
                local_dir = ISOSDB_LOCAL,
            conda:
                "../../envs/tncentral.yaml"
            log:
                LOGS + "/mobilome_isosdb_db.log"
            priority: 4
            shell:
                """
                exec > {log} 2>&1
                set -euo pipefail
                mkdir -p {output.db_dir}

                if [ -n "{params.local_dir}" ]; then
                    echo "Using the local ISOSDB directory '{params.local_dir}'."
                    cp "{params.local_dir}"/ISOSDB.V3.fna {output.db_dir}/ISOSDB.V3.fna
                    cp "{params.local_dir}"/IS_fam_annot.txt {output.db_dir}/IS_fam_annot.txt
                    SOURCE="local:{params.local_dir}"
                else
                    curl -sSL -o {output.db_dir}/isosdb.zip "{params.fasta_url}"
                    if ! unzip -t {output.db_dir}/isosdb.zip > /dev/null 2>&1; then
                        echo "ERROR: the ISOSDB download is not a ZIP archive. First bytes:" >&2
                        head -c 200 {output.db_dir}/isosdb.zip >&2
                        exit 1
                    fi
                    unzip -o -q -j {output.db_dir}/isosdb.zip -d {output.db_dir}
                    rm -f {output.db_dir}/isosdb.zip
                    curl -sSL -o {output.db_dir}/IS_fam_annot.txt "{params.family_url}"
                    SOURCE="{params.fasta_url}"
                fi

                # An IS family map that did not download leaves every element
                # "unassigned", which looks like a real result rather than a
                # missing file. Fail instead.
                if [ ! -s {output.db_dir}/IS_fam_annot.txt ]; then
                    echo "ERROR: IS_fam_annot.txt is missing or empty; families could not be" >&2
                    echo "       assigned and the summary would silently be meaningless." >&2
                    exit 1
                fi

                N_SEQ=$(grep -c '^>' {output.db_dir}/ISOSDB.V3.fna)
                {{
                  echo "source:      $SOURCE"
                  echo "fetched:     $(date -u +%Y-%m-%dT%H:%M:%SZ)"
                  echo "fasta_sha256: $(sha256sum {output.db_dir}/ISOSDB.V3.fna | cut -d' ' -f1)"
                  echo "sequences:   $N_SEQ"
                  echo ""
                  echo "ISOSDB is distributed in the pseudoR repository under the MIT licence."
                  echo "Cite Kirsch et al. 2024, Cell Host & Microbe."
                }} > {output.db_dir}/PROVENANCE.txt

                echo "ISOSDB ready: $N_SEQ sequences."
                """

        # ── Rule: assembly_depth — what does single-copy look like? ──────────
        # Biology: THE DENOMINATOR. Mapping the same reads to the sample's OWN
        # assembly gives the depth of ordinary single-copy sequence. Without it a
        # raw IS depth means nothing, because it scales with how deeply the
        # sample happened to be sequenced — 200x over an IS is unremarkable in a
        # 200x library and means four copies in a 50x one.
        #
        # Takes in: the sample's contigs (the D2 hand-off) and its trimmed
        #           Illumina reads (rule trim_adapters, in the illumina/hybrid
        #           front end).
        # Does: BBMap the reads back onto the assembly. ambiguous=best and
        #       secondary=f match rule isosdb_map exactly, so the two depths are
        #       measured the same way and their ratio means something.
        # Produces: ASSEMBLY_COVSTATS, BBMap's per-contig coverage table. The
        #           BBMap index is a temp() directory and is deleted after.
        # Consumed by: is_copy_number.
        rule assembly_depth:
            input:
                contigs = FINAL_CONTIGS,
                r1 = TRIM_R1,
                r2 = TRIM_R2,
            output:
                covstats = ASSEMBLY_COVSTATS,
                ref_dir = temp(directory(MOBILOME_DIR + "/assembly_depth_ref")),
            params:
                max_ram = min(RAM, 32),
            conda:
                "../../envs/bbmap.yaml"
            threads: capped_cpus(16)
            log:
                LOGS + "/mobilome_assembly_depth_{sample}.log"
            priority: 3
            shell:
                """
                bbmap.sh \
                  -in={input.r1} -in2={input.r2} \
                  ref={input.contigs} path={output.ref_dir} \
                  -Xmx{params.max_ram}g threads={threads} \
                  ambiguous=best secondary=f \
                  covstats={output.covstats} > {log} 2>&1
                """

        # ── Rule: isosdb_map — how deep are the IS elements? ─────────────────
        # Biology: THE NUMERATOR. Reads are immune to assembly collapse — every
        # copy of an IS contributes its own reads whether or not the assembler
        # kept them apart — so an IS present in five copies attracts about five
        # times the depth of single-copy sequence.
        #
        # Takes in: the ISOSDB sequences (rule isosdb_db) and the SAME trimmed
        #           reads rule assembly_depth used.
        # Does: BBMap with settings identical to assembly_depth. ambiguous=best on
        #       purpose: ISOSDB is dereplicated at 95% but families remain
        #       similar, and letting one read count for every near-identical entry
        #       would multiply the totals. The cost is that the split BETWEEN
        #       near-identical entries is arbitrary, which is why the next rule
        #       reports the FAMILY sum as the headline.
        # Produces: ISOSDB_COVSTATS, per-database-entry coverage.
        # Consumed by: is_copy_number.
        rule isosdb_map:
            input:
                db_dir = ISOSDB_DB_DIR,
                r1 = TRIM_R1,
                r2 = TRIM_R2,
            output:
                covstats = ISOSDB_COVSTATS,
                ref_dir = temp(directory(MOBILOME_DIR + "/isosdb_ref")),
            params:
                max_ram = min(RAM, 32),
                fasta = lambda w, input: os.path.join(input.db_dir, "ISOSDB.V3.fna"),
            conda:
                "../../envs/bbmap.yaml"
            threads: capped_cpus(16)
            log:
                LOGS + "/mobilome_isosdb_map_{sample}.log"
            priority: 3
            shell:
                """
                bbmap.sh \
                  -in={input.r1} -in2={input.r2} \
                  ref={params.fasta} path={output.ref_dir} \
                  -Xmx{params.max_ram}g threads={threads} \
                  ambiguous=best secondary=f \
                  covstats={output.covstats} > {log} 2>&1
                """

        # ── Rule: is_copy_number — located vs implied ────────────────────────
        # Biology: the point of the whole leg. Divide the IS depth by the genome
        # depth and you get a copy number that the assembly could not collapse:
        #
        #     copy number  ~=  depth over the IS  /  depth over the genome
        #
        # Put that next to how many copies ISEScan actually LOCATED on the contigs
        # and the difference is the assembler's collapse, measured rather than
        # warned about.
        #
        # Takes in: both coverage tables (isosdb_map and assembly_depth), the
        #           family map (isosdb_db), and the located IS calls
        #           (isescan_table).
        # Does: convert depths to copy numbers, discard entries covered over less
        #       than --min-covered-percent of their length or sitting below
        #       --min-copy-number (each with a reason in the audit), and report
        #       both per entry and summed per IS family.
        # Produces: IS_COPY_NUMBER + IS_COPY_NUMBER_AUDIT.
        # Consumed by: NOBODY downstream — deliberately. This is a workflow target
        #       in its own right, read by a person. It changes no AMR gene's tier
        #       and must not: it says how many copies exist, never WHERE they are,
        #       so it cannot place a gene inside anything.
        rule is_copy_number:
            input:
                isosdb_covstats = ISOSDB_COVSTATS,
                assembly_covstats = ASSEMBLY_COVSTATS,
                db_dir = ISOSDB_DB_DIR,
                is_table = IS_TABLE,
            output:
                table = IS_COPY_NUMBER,
                audit = IS_COPY_NUMBER_AUDIT,
            params:
                script = ISOSDB_COPY_SCRIPT,
                family_map = lambda w, input: os.path.join(input.db_dir, "IS_fam_annot.txt"),
                min_covered = ISOSDB_MIN_COVERED,
                min_copies = ISOSDB_MIN_COPIES,
            log:
                LOGS + "/mobilome_is_copy_number_{sample}.log"
            priority: 3
            shell:
                """
                python {params.script} \
                  --sample {wildcards.sample} \
                  --isosdb-covstats {input.isosdb_covstats} \
                  --assembly-covstats {input.assembly_covstats} \
                  --family-map {params.family_map} \
                  --is-table {input.is_table} \
                  --min-covered-percent {params.min_covered} \
                  --min-copy-number {params.min_copies} \
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
    #              for anything sitting on a plasmid, and by conjscan_ice, which
    #              needs it to tell a conjugative PLASMID from an ICE.
    #
    # NO SEPARATE AUDIT FILE, and that is deliberate — it is the one mobilome
    # rule without one. The project rule is that every decision states a reason,
    # not that every reason lives in its own file: this step drops no rows (one
    # row per contig, always), so there is nothing to explain the absence of.
    # The reasons ride along on the table itself, in replicon_call_source and
    # mobility_evidence.
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
            # GENOMAD_PREFIX reads oddly here, since this is a Platon parser. It
            # is not a geNomad-only constant: it is the basename of the input
            # genome with the extension stripped (the literal "contigs_final"),
            # and BOTH tools name every output file after their input's basename.
            # One constant therefore names both tools' per-sample files, which is
            # why 00_common.smk defines it once and reuses it deliberately.
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
    # Takes in, and which rule produced each:
    #   * the AMR calls with their coordinates                    (amrfinderplus)
    #   * the IS elements                                         (isescan_table)
    #   * the ICE/IME candidates — the NAMED copy when the ICEberg layer is on,
    #     the raw one otherwise (ICE_TABLE_FOR_COLOCALISE decides)
    #                                    (name_ice_elements or conjscan_ice)
    #   * curated transposons and integrons, only when that layer is configured
    #                                                             (name_elements)
    #   * chromosome-or-plasmid per contig                  (mobilome_replicons)
    #   * the contig lengths, for the contig-edge honesty flags  (contig_lengths)
    #   --is-table is repeatable, which is how the three element sources arrive
    #   as one pool of elements rather than three special cases in the script.
    # Does: for every AMR gene, find the elements around it, apply the ladder
    #       above top-down (highest rung that the evidence supports wins), and
    #       record what capped the confidence.
    # Produces: MOBILITY_TABLE (the deliverable, one row per AMR gene) and
    #           MOBILITY_AUDIT (why every gene without context got none).
    # Consumed by: the user. This is the module's terminal product.
    #
    # Its log is mobilome_colocalisation_{sample}.log rather than the rule name
    # in full; left alone, like name_ice_elements's, because renaming it only
    # moves files on disk.
    rule amr_mge_colocalisation:
        input:
            amrfinder = AMRFINDER_TSV,
            is_table = IS_TABLE,
            ice_table = ICE_TABLE_FOR_COLOCALISE,
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
