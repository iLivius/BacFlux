# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — Stage 07 phage module (rules/shared/70_phage.smk)  [D8]
#
# This module finds viruses / prophages in every finished genome and grades their
# quality. Like the rest of the shared tail it consumes only FINAL_CONTIGS (D2),
# so it is identical in all four modes.
#
# Two callers, one QC step:
#   - VirSorter2 (DEFAULT) — permissively licensed (GPLv2, commercial use OK), so
#     it is the out-of-the-box caller. Its conda env is fragile (see the fix note
#     at its rules), so its two rules are defined behind an
#     `if PHAGE_CALLER == "virsorter2":` guard.
#   - geNomad (OPT-IN) — actively maintained and does viruses AND plasmids in one
#     end-to-end run, BUT is licensed ACADEMIC / NON-COMMERCIAL-USE-ONLY (Berkeley
#     Lab). BacFlux is MIT and must not force a non-commercial restriction on
#     users, so geNomad is opt-in (config.phage.caller: genomad), NEVER the
#     default: its rules are defined behind an `if PHAGE_CALLER == "genomad":`
#     guard, so on a default run geNomad never builds its env or runs. When opted
#     in, that same run also feeds the plasmid module's D9 concordance
#     (60_plasmid.smk) — see PHAGE_CALLER in 00_common.smk.
#   - CheckV (ALWAYS) — neither caller reports completeness/contamination, so
#     CheckV runs unconditionally on whichever caller's virus FASTA was produced.
#
# Exactly ONE caller's rules are defined per run (the guards are mutually
# exclusive), so there is never an ambiguity over who produces the virus FASTA.
#
# Data flow (top to bottom):
#
#   contigs_final.fasta ─(caller==genomad)─► genomad_end_to_end ─► 07.phages/genomad/{sample}/
#          │                       │  (…_summary/contigs_final_virus.fna  → CheckV)
#          │                       └  (…_summary/contigs_final_plasmid_summary.tsv
#          │                          → consumed by 60_plasmid.smk, the D9 concordance)
#          │
#          └(caller==virsorter2, default)─► viral_identification_virsorter2 ─►
#                                          07.phages/virsorter/{sample}/final-viral-combined.fa
#                                                        │
#          whichever caller ran ─────────────────────────┴─► viral_quality (CheckV)
#                                          ─► 07.phages/checkv/{sample}/quality_summary.tsv
#
# Databases: separate one-off downloads (checkv_db always; genomad_db OR
# virsorter2_db depending on the caller), each in its OWN conda env — replacing
# v1's single fused `viral_db` rule, so the default path never builds geNomad's or
# VirSorter2's env unnecessarily.
#
# CROSS-STAGE NOTE: when geNomad is opted in, 60_plasmid.smk (stage 06) reads a file
# that physically lives under THIS stage's 07.phages/genomad/{sample}/ directory.
# That backwards-numbered 06←07 edge is safe: stage numbers here are organisational
# only (D1 groups by tool family); Snakemake orders work by the input/output DAG.
#
# Everything referenced here is defined once in 00_common.smk (never re-derived):
# FINAL_CONTIGS, DIR_PHAGES, PHAGE_CALLER, GENOMAD_DB_DIR, GENOMAD_DIR,
# GENOMAD_PREFIX, VS2_DB_DIR, VS2_DIR, CHECKV_DB_DIR, CHECKV_LINK, CHECKV_DB_ID,
# LOGS, capped_cpus.
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/x.yaml" climbs shared/ -> rules/ -> workflow/ -> workflow/envs/x.yaml.
#
# Resource convention: cpu-bound rules declare Snakemake's built-in
# `threads: capped_cpus(N)` and refer to `{threads}` in the shell. Using the
# BUILT-IN keyword (rather than a custom `resources: cpus`) is what makes
# `--cores N` actually enforce the limit, so a plain `snakemake --cores N` is
# safe on its own and no extra `--resources` flag is needed.
# ─────────────────────────────────────────────────────────────────────────────


# ── Parse-time caller selection for CheckV ───────────────────────────────────
# CheckV takes the SAME command whichever caller ran; only the input virus FASTA
# differs. Resolve that choice ONCE here (PHAGE_CALLER comes from 00_common), so
# the CheckV rule below stays caller-agnostic. CHECKV_CALLER_DIR is the {sample}-
# templated DIRECTORY output of the chosen caller; CHECKV_VIRAL_REL is the fixed
# relative path of the virus FASTA inside it (NON-templated — no {sample} — so it
# is safe to pass as a plain params string, which Snakemake does not expand).
if PHAGE_CALLER == "genomad":
    CHECKV_CALLER_DIR = GENOMAD_DIR
    CHECKV_VIRAL_REL = GENOMAD_PREFIX + "_summary/" + GENOMAD_PREFIX + "_virus.fna"
else:  # virsorter2 (default)
    CHECKV_CALLER_DIR = VS2_DIR
    CHECKV_VIRAL_REL = "final-viral-combined.fa"


# ── geNomad OPT-IN path (only defined when PHAGE_CALLER == "genomad") ─────────
# Wrapping both geNomad rules in this guard means that on the default (VirSorter2)
# path they are never defined — geNomad's env is never built and it never runs, so
# a non-academic / commercial user is never forced to execute a non-commercial
# tool. Opting in (config.phage.caller: genomad) is the user's own informed choice.
if PHAGE_CALLER == "genomad":

    # ── Rule: genomad_db — one-off geNomad reference download ────────────────
    # Biology: geNomad classifies sequences with marker-gene profiles bundled in a
    # versioned database it must fetch once (~1.5 GB extracted). geNomad's own
    # downloader pulls the Zenodo tarball and extracts it.
    #
    # Takes in: nothing (pure download).
    # Does: `genomad download-database <parent>` — geNomad always creates the
    #       subfolder "genomad_db" inside <parent>, so we point <parent> at
    #       DIR_PHAGES and declare the output as DIR_PHAGES/genomad_db.
    # Produces: 07.phages/genomad_db/ — the shared geNomad DB directory.
    # Consumed by: genomad_end_to_end (every sample waits on this DB).
    rule genomad_db:
        output:
            genomad_db = directory(GENOMAD_DB_DIR),
        params:
            # geNomad appends "genomad_db" to this parent path; see GENOMAD_DB_DIR.
            parent = DIR_PHAGES,
        conda:
            "../../envs/genomad.yaml"
        log:
            LOGS + "/genomad_db.log"
        priority: 9
        shell:
            """
            genomad download-database {params.parent} > {log} 2>&1
            """

    # ── Rule: genomad_end_to_end — virus + plasmid calling in one run ────────
    # Biology: geNomad scans the finished genome and, in a single end-to-end run,
    # reports which contigs (or contig regions) are viral and which are plasmids.
    # When opted in, its virus output feeds CheckV here, and its plasmid summary
    # feeds the D9 concordance in 60_plasmid.smk.
    #
    # Takes in:
    #   contigs    = FINAL_CONTIGS — the finished, decontaminated assembly (D4: all
    #                modes scan decontaminated contigs; v1 illumina scanned the
    #                pre-decontam contigs_filt.fasta — flagged in the changelog).
    #   genomad_db = the geNomad DB directory (from genomad_db).
    # Does: `genomad end-to-end` with positional args INPUT OUTPUT DATABASE. Default
    #       presets (neither --conservative nor --relaxed) match the standard
    #       geNomad→CheckV combo; --cleanup deletes intermediates for disk hygiene.
    #       (RAM lever if ever needed on a huge input: `--splits N` caps peak memory
    #       at a speed cost — not needed for isolate-sized genomes.)
    # Produces: 07.phages/genomad/{sample}/ — declared a DIRECTORY so a failed/rerun
    #       job is wiped clean, stopping geNomad from resuming on stale
    #       intermediates. geNomad names every output after the input basename
    #       (contigs_final), so the key files land at:
    #         …/contigs_final_summary/contigs_final_virus.fna          (→ CheckV)
    #         …/contigs_final_summary/contigs_final_plasmid_summary.tsv (→ 60_plasmid)
    # Consumed by: viral_quality (CheckV) and plasmid_concordance (60_plasmid.smk) —
    #              each depends on this DIRECTORY and reaches inside for its file,
    #              keeping the inner-file edge off the DAG (as 40_annotation does for
    #              the Bakta directory).
    # VERIFY on the first real geNomad run: that end-to-end runs fresh into the
    # Snakemake-managed output dir (no unwanted resume); if it resumes, add
    # `--restart`. Output filenames follow the basename pattern above — confirm.
    rule genomad_end_to_end:
        input:
            contigs = FINAL_CONTIGS,
            genomad_db = GENOMAD_DB_DIR,
        output:
            genomad_dir = directory(GENOMAD_DIR),
        conda:
            "../../envs/genomad.yaml"
        threads: capped_cpus(24)
        log:
            LOGS + "/genomad_{sample}.log"
        priority: 8
        shell:
            """
            genomad end-to-end \
              --cleanup \
              --threads {threads} \
              {input.contigs} \
              {output.genomad_dir} \
              {input.genomad_db} > {log} 2>&1
            """


# ── VirSorter2 DEFAULT path (only defined when PHAGE_CALLER == "virsorter2") ──
# VS2 2.2.4 (Jan 2023) is itself a Snakemake workflow, and it manages its OWN
# nested conda env at runtime. That nesting is what crashed every v1.3.1 phage run
# (an ancient transitive `mamba` against a modern `conda` -> "No module named
# 'conda._vendor.auxlib'", the `virsorter_deps_env` saga).
#
# The fix is to switch the nesting OFF entirely and make our single env carry
# everything VS2 needs. Three parts, applied below and in envs/virsorter.yaml:
#   1. envs/virsorter.yaml installs VS2's OWN internal dependency list (copied
#      from the `envs/vs2.yaml` that ships inside the virsorter package) next to
#      virsorter itself, so the one env is self-sufficient.
#   2. `virsorter setup --skip-deps-install` downloads only the DB and does NOT
#      build VS2's per-rule nested dependency envs.
#   3. `virsorter run --use-conda-off`, plus exporting the env's own bin onto PATH,
#      makes VS2 subprocesses resolve OUR in-env binaries instead of nested envs.
#
# VERIFIED 2026-07-22 on the v2 illumina validation run. The check this comment
# used to ask for has been done, and the answer was NOT the comfortable one: the
# bioconda `virsorter=2.2.4` package does NOT pin the runtime tool closure. With
# only `virsorter` in the env, screed, hmmer, prodigal, last, pandas,
# scikit-learn, numpy, seaborn, imbalanced-learn and ncbi-genome-download were ALL
# absent, and the run died at the first internal rule on `No module named
# 'screed'`. Hence part 1 above. If VS2 is ever unpinned from 2.2.4, re-read its
# packaged envs/vs2.yaml and re-sync envs/virsorter.yaml against it.
if PHAGE_CALLER == "virsorter2":

    # ── Rule: virsorter2_db — one-off VirSorter2 reference download ──────────
    # Biology: VirSorter2 scores contigs against curated viral HMM groups from a
    # database it fetches once.
    # Takes in: nothing (pure download).
    # Does: `virsorter setup` to download the DB; --skip-deps-install skips the
    #       fragile nested-env build (fix part 2 above).
    # Produces: 07.phages/vs2_db/ (VS2_DB_DIR).
    # Consumed by: viral_identification_virsorter2.
    rule virsorter2_db:
        output:
            vs2_db = directory(VS2_DB_DIR),
        conda:
            "../../envs/virsorter.yaml"
        threads: capped_cpus(4)
        log:
            LOGS + "/virsorter2_db.log"
        priority: 9
        shell:
            """
            virsorter setup \
              -d {output.vs2_db} \
              -j {threads} \
              --skip-deps-install > {log} 2>&1
            """

    # ── Rule: viral_identification_virsorter2 — VS2 virus calling ────────────
    # Biology: identify phages / prophages on the finished genome. Faithful port
    # of v1's `viral_identification` (VS2 half) onto the 00_common API.
    # Takes in: contigs = FINAL_CONTIGS (D4 change: v1 illumina used the
    #           pre-decontam contigs_filt.fasta), vs2_db = the VS2 DB directory.
    # Does: `virsorter run ... all`, over the same viral groups and min-score as
    #       v1, with --use-conda-off and a PATH preamble (fix part 3 above).
    # Produces: 07.phages/virsorter/{sample}/ — key file final-viral-combined.fa.
    # Consumed by: viral_quality (CheckV).
    rule viral_identification_virsorter2:
        input:
            contigs = FINAL_CONTIGS,
            vs2_db = VS2_DB_DIR,
        output:
            vs2_dir = directory(VS2_DIR),
        params:
            viral_groups = "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae",
            min_score = 0.5,
        conda:
            "../../envs/virsorter.yaml"
        threads: capped_cpus(24)
        log:
            LOGS + "/viral_identification_{sample}.log"
        priority: 8
        shell:
            # Put the env's own bin first so VS2 subprocesses resolve in-env tools.
            """
            export PATH="$CONDA_PREFIX/bin:$PATH"
            virsorter run \
              -i {input.contigs} \
              -w {output.vs2_dir} \
              -d {input.vs2_db} \
              --keep-original-seq \
              --include-groups {params.viral_groups} \
              --min-score {params.min_score} \
              --use-conda-off \
              -j {threads} \
              all > {log} 2>&1
            """


# ── Rule: checkv_db — one-off CheckV reference download ──────────────────────
# Biology: CheckV grades viral genome completeness/contamination against a
# reference database. Faithful port of v1 `viral_db` (CheckV half) onto its OWN
# env (envs/checkv.yaml), so the default path never builds another tool's env just
# to get CheckV. CheckV is permissively licensed (LBNL BSD), commercial use OK.
#
# Takes in: nothing (pure download).
# Does: if CHECKV_LINK is empty, let CheckV download its own default DB; otherwise
#       wget the given .tar.gz, extract it, and (re)build the diamond DB. Both the
#       link and the derived folder id (CHECKV_DB_ID) come from 00_common.
# Produces: 07.phages/checkv_db/ (CHECKV_DB_DIR).
# Consumed by: viral_quality.
# DEFINED ONLY when BacFlux is the one providing the database. If the user pointed
# directories.checkv_db at a copy they already hold, this rule must NOT exist.
#
# That is a safety requirement, not tidiness. The output below is a `directory()`,
# and Snakemake DELETES a directory output before re-running its rule. If this rule
# were defined with CHECKV_DB_DIR pointing at a shared database on a NAS, then any
# trigger to re-run it — a changed env file, a --forcerun, an interrupted job —
# would wipe that shared database for everyone using it. Leaving the rule undefined
# means Snakemake treats the path as a plain existing input it may only read.
if not CHECKVDB:

    rule checkv_db:
        output:
            checkv_db = directory(CHECKV_DB_DIR),
        params:
            checkv_link = CHECKV_LINK,
            db_id = CHECKV_DB_ID,
            tries = 5,
        conda:
            "../../envs/checkv.yaml"
        log:
            LOGS + "/checkv_db.log"
        priority: 9
        shell:
            """
            if [ -z "{params.checkv_link}" ]; then
                checkv download_database {output.checkv_db} > {log} 2>&1
            else
                wget --tries={params.tries} -c {params.checkv_link} -P {output.checkv_db} > {log} 2>&1
                tar -xzvf {output.checkv_db}/{params.db_id}.tar.gz -C {output.checkv_db} >> {log} 2>&1
                diamond makedb \
                  --in {output.checkv_db}/{params.db_id}/genome_db/checkv_reps.faa \
                  --db {output.checkv_db}/{params.db_id}/genome_db/checkv_reps >> {log} 2>&1
            fi
            """


# ── Rule: viral_quality — completeness/contamination of virus calls (CheckV) ──
# Biology: CheckV estimates how complete each predicted viral sequence is and
# flags host contamination. It is caller-agnostic: the SAME `checkv end_to_end`
# command grades geNomad's virus FASTA or VirSorter2's, so this single rule serves
# both paths (the input FASTA was chosen once at parse time, above).
#
# Takes in:
#   caller_dir = CHECKV_CALLER_DIR — the chosen caller's {sample} output DIRECTORY.
#                The shell reaches inside for the virus FASTA at the fixed relative
#                path CHECKV_VIRAL_REL (params.viral_rel).
#   checkv_db  = the CheckV DB directory (from checkv_db).
# Does: resolve the actual DB sub-directory at runtime (ported verbatim from v1),
#       then run `checkv end_to_end`.
# Produces: 07.phages/checkv/{sample}/ — key file quality_summary.tsv. This is the
#       rule-all leaf that "pulls in whichever caller ran".
# Consumed by: the report / the user (a terminal phage product).
#
# GRACEFUL EMPTY-INPUT HANDLING (BacFlux "degrade, don't hard-fail" convention):
# a virus-free genome is common, and both callers can emit an EMPTY (or absent)
# virus FASTA. `checkv end_to_end` errors on an empty input, which would fail the
# whole sample. So we test the FASTA first: if it has no sequences, we write a
# header-only quality_summary.tsv and skip CheckV, rather than crash the run.
rule viral_quality:
    input:
        caller_dir = CHECKV_CALLER_DIR,
        checkv_db = CHECKV_DB_DIR,
    output:
        checkv_dir = directory(DIR_PHAGES + "/checkv/{sample}"),
    params:
        viral_rel = CHECKV_VIRAL_REL,
    conda:
        "../../envs/checkv.yaml"
    threads: capped_cpus(24)
    log:
        LOGS + "/viral_quality_{sample}.log"
    priority: 7
    shell:
        # Runtime DB-dir resolution (verbatim from v1): find the single
        # genome_db/checkv_reps.faa, assert exactly one, then take its
        # grandparent as the DB directory CheckV expects.
        """
        mkdir -p {output.checkv_dir}
        viral_fasta="{input.caller_dir}/{params.viral_rel}"

        # Empty/absent virus FASTA => no viral sequences to grade. Write a
        # header-only quality_summary.tsv and skip CheckV so the run does not fail.
        if [ ! -s "$viral_fasta" ] || ! grep -q '>' "$viral_fasta"; then
            echo "No viral sequences called for {wildcards.sample}; writing empty CheckV summary." > {log}
            printf "contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\twarnings\n" \
              > {output.checkv_dir}/quality_summary.tsv
        else
            checkv_rep_files=$(find {input.checkv_db} -type f -path "*/genome_db/checkv_reps.faa" | sort)
            checkv_rep_file=$(printf "%s\n" "$checkv_rep_files" | head -n 1)
            checkv_rep_count=$(printf "%s\n" "$checkv_rep_files" | sed '/^$/d' | wc -l)
            if [ "$checkv_rep_count" -ne 1 ]; then
                echo "Expected exactly one CheckV database in {input.checkv_db}, found $checkv_rep_count candidate(s)." > {log}
                exit 1
            fi
            checkv_db_dir=$(dirname "$(dirname "$checkv_rep_file")")

            checkv end_to_end \
              "$viral_fasta" \
              {output.checkv_dir} \
              -t {threads} \
              -d "$checkv_db_dir" > {log} 2>&1
        fi
        """
