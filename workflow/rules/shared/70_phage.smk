# BacFlux v2.0.0 — virus and prophage detection on the finished genome (stage
# 07.phages), followed by CheckV quality grading. Like the rest of the shared tail
# this module reads only FINAL_CONTIGS (the D2 hand-off), so it behaves identically
# in all four modes.
#
# Two callers are wired in and exactly one of them is defined per run, chosen by
# config phage.caller (resolved into PHAGE_CALLER in 00_common.smk; the choice
# itself is decision D8 in docs/unification_migration_plan.md):
#   virsorter2 (DEFAULT) — permissively licensed (GPLv2, commercial use OK), so it
#                          is what an out-of-the-box run gets. Its conda env is the
#                          fragile part — the whole story is at the VirSorter2
#                          caller section of this file.
#   genomad    (OPT-IN)  — actively maintained, and calls viruses AND plasmids in a
#                          single run, but Berkeley Lab licenses it for academic /
#                          non-commercial use only. BacFlux is MIT and must not
#                          push that restriction onto its users, so geNomad is
#                          never the default — see README Licensing. Opting in also
#                          switches on the plasmid concordance in
#                          shared/60_plasmid.smk (decision D9).
#   CheckV     (ALWAYS)  — neither caller reports completeness or contamination, so
#                          CheckV 1.0.3 grades whatever the caller produced.
#
# The callers sit behind mutually exclusive `if PHAGE_CALLER == ...` guards, so only
# one set of rules is ever defined and there is never a question about who produced
# the virus FASTA.
#
#   contigs_final.fasta
#     ├─(genomad)─► genomad_end_to_end ─► 07.phages/genomad/{sample}/
#     │                …_summary/contigs_final_virus.fna ─────────────┐
#     │                …_summary/…_plasmid_summary.tsv → 60_plasmid.smk
#     │                                                               │
#     └─(virsorter2, default)─► viral_identification_virsorter2 ──────┤
#                07.phages/virsorter/{sample}/final-viral-combined.fa │
#                                                                     ▼
#          viral_quality (CheckV) ─► 07.phages/checkv/{sample}/quality_summary.tsv
#
# genomad_db                      : download geNomad's marker database once, from
#                                   the configured mirror or geNomad's own
#                                   downloader.
# genomad_db_local                : symlink view of a geNomad database you already
#                                   hold.
# genomad_end_to_end              : virus + plasmid calling in one run, per sample.
# virsorter2_db                   : download VirSorter2's reference database once.
# virsorter2_db_local             : symlink view of a VirSorter2 database you
#                                   already hold.
# viral_identification_virsorter2 : VirSorter2 virus calling, per sample.
# checkv_db_local                 : local view of a CheckV database you already
#                                   hold, with the DIAMOND index rebuilt here.
# checkv_db                       : download, checksum and index CheckV's database.
# viral_quality                   : CheckV completeness/contamination on whichever
#                                   caller ran.
#
# Every database therefore comes in a matched pair of rules — one that downloads
# it, one that builds a view of a copy the user already holds — and exactly one of
# each pair is defined, gated on the matching directories.* config key. v1 fused
# all of this into a single `viral_db` rule, which meant the default path built
# another tool's conda env just to get CheckV.
#
# One cross-stage edge: when geNomad is opted in, shared/60_plasmid.smk (stage 06)
# reads a file that physically lives under THIS stage's 07.phages/genomad/{sample}/
# directory. Reading "backwards" from 06 into 07 is safe: the stage numbers group
# tools for the reader (D1), and Snakemake orders work by the input/output DAG.
#
# Defined once in 00_common.smk and never re-derived here: FINAL_CONTIGS,
# DIR_PHAGES, PHAGE_CALLER, GENOMAD_DB_DIR, GENOMAD_DIR, GENOMAD_PREFIX,
# GENOMAD_LINK, GENOMAD_MD5, VS2_DB_DIR, VS2_DIR, CHECKV_DB_DIR, CHECKV_LINK,
# CHECKV_SHA_URL, CHECKV_DB_ID, LOGS, capped_cpus.
#
# conda: env paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/x.yaml" climbs shared/ → rules/ → workflow/ → workflow/envs/x.yaml.
#
# Resource convention: cpu-bound rules declare Snakemake's built-in
# `threads: capped_cpus(N)` and refer to `{threads}` in the shell. Using the
# BUILT-IN keyword (rather than a custom `resources: cpus`) is what makes
# `--cores N` actually enforce the limit, so a plain `snakemake --cores N` is safe
# on its own and no extra `--resources` flag is needed.


# ─────────────────── CheckV input selection ────────────────────
# CheckV takes the SAME command whichever caller ran; only the input virus FASTA
# differs. That choice is resolved ONCE here, at parse time (before any job runs),
# so the CheckV rule itself stays caller-agnostic. CHECKV_CALLER_DIR is the
# {sample}-templated DIRECTORY output of the chosen caller; CHECKV_VIRAL_REL is the
# fixed relative path of the virus FASTA inside it — NON-templated, no {sample}, so
# it is safe to hand to the rule as a plain params string, which Snakemake does not
# expand.
if PHAGE_CALLER == "genomad":
    CHECKV_CALLER_DIR = GENOMAD_DIR
    CHECKV_VIRAL_REL = GENOMAD_PREFIX + "_summary/" + GENOMAD_PREFIX + "_virus.fna"
else:  # virsorter2 (default)
    CHECKV_CALLER_DIR = VS2_DIR
    CHECKV_VIRAL_REL = "final-viral-combined.fa"


# ─────────────────── geNomad caller (opt-in) ───────────────────
# Wrapping the geNomad rules in this guard is the whole mechanism behind "opt-in":
# on the default (VirSorter2) path they are never defined, geNomad's env is never
# built and it never runs, so a commercial user is never made to execute a
# non-commercial tool. Setting config phage.caller: genomad is that user's own
# informed choice, and it also turns on the plasmid concordance in
# shared/60_plasmid.smk.
if PHAGE_CALLER == "genomad":

    # ── geNomad database, downloaded here ──
    # geNomad (v1.12.0, pinned in envs/genomad.yaml) classifies sequences against
    # marker-gene profiles bundled in a versioned database it must fetch once
    # (~1.5 GB extracted). Consumed by genomad_end_to_end, which every sample waits
    # on.
    #
    # Two ways to get it, chosen by whether links.genomad_link is set:
    #
    #   links.genomad_link SET (the shipped default) — we fetch the archive
    #     ourselves from the mirror. geNomad's authors publish the same database on
    #     Zenodo and link to it from their own README, so this is the same data,
    #     just from a host that is actually up. The archive expands to a top-level
    #     genomad_db/ directory, which is exactly GENOMAD_DB_DIR, so it is
    #     extracted into DIR_PHAGES (verified against the real archive).
    #
    #   links.genomad_link EMPTY — fall back to `genomad download-database`, whose
    #     URL is hard-coded to portal.nersc.gov inside the package. That host is
    #     frequently unreachable and there is no --url option, so if this branch
    #     fails with "No route to host", set links.genomad_link (or
    #     directories.genomad_db) rather than retrying.
    #
    # VERIFIED 2026-07-24: `genomad download-database DESTINATION` does create the
    # genomad_db subfolder inside DESTINATION, so pointing it at DIR_PHAGES and
    # declaring DIR_PHAGES/genomad_db as the output is right.
    #
    # Integrity: Zenodo publishes an MD5 per file rather than the .sha256 sidecar
    # the CheckV and dbCAN mirrors carry, so the expected hash comes from
    # links.genomad_md5. A mismatch is fatal (same discipline as checkv_db and
    # cazyme_db_download — a truncated archive must never reach `tar`). An empty
    # hash downloads unverified and SAYS so in the log.
    #
    # Defined only when BacFlux is the one downloading — mutually exclusive with
    # genomad_db_local, same reasoning as checkv_db / virsorter2_db.
    if not GENOMADDB:

        rule genomad_db:
            output:
                genomad_db = directory(GENOMAD_DB_DIR),
            params:
                # geNomad appends "genomad_db" to this parent path; see
                # GENOMAD_DB_DIR. The Zenodo archive also expands to genomad_db/,
                # so both branches land in the same place.
                parent = DIR_PHAGES,
                link = GENOMAD_LINK,
                md5 = GENOMAD_MD5,
                tries = 5,
            conda:
                "../../envs/genomad.yaml"
            log:
                LOGS + "/genomad_db.log"
            shell:
                """
                if [ -z "{params.link}" ]; then
                    echo "No links.genomad_link set; using geNomad's own downloader (portal.nersc.gov)." > {log}
                    genomad download-database {params.parent} >> {log} 2>&1
                else
                    mkdir -p {params.parent}
                    TAR="{params.parent}/$(basename '{params.link}')"
                    echo "Fetching the geNomad database from the configured mirror:" > {log}
                    echo "  {params.link}" >> {log}
                    wget --tries={params.tries} -c "{params.link}" -O "$TAR" >> {log} 2>&1

                    if [ -n "{params.md5}" ]; then
                        actual="$(md5sum "$TAR" | awk '{{print $1}}')"
                        if [ "{params.md5}" != "$actual" ]; then
                            echo "ERROR: checksum mismatch for $TAR" >> {log}
                            echo "  expected (links.genomad_md5): {params.md5}" >> {log}
                            echo "  actual:                       $actual" >> {log}
                            echo "If you changed links.genomad_link, update links.genomad_md5" >> {log}
                            echo "to match the new file (Zenodo shows the MD5 next to it), or" >> {log}
                            echo "clear it to download without verification." >> {log}
                            rm -f "$TAR"
                            exit 1
                        fi
                        echo "Checksum OK ({params.md5})." >> {log}
                    else
                        echo "links.genomad_md5 is empty: extracting WITHOUT verifying the download." >> {log}
                    fi

                    # The archive's own top-level directory is genomad_db/, so this
                    # produces {output.genomad_db}. --no-same-permissions keeps a
                    # restrictive archive from producing a database the workflow
                    # cannot read back (geNomad reads every file, including
                    # version.txt).
                    tar -xzf "$TAR" -C {params.parent} --no-same-permissions >> {log} 2>&1
                    rm -f "$TAR"
                    chmod -R a+rX {output.genomad_db} >> {log} 2>&1 || true
                    echo "geNomad database ready at {output.genomad_db}." >> {log}
                fi
                """


    # ── geNomad database, a copy you already hold ──
    # Defined only when directories.genomad_db is set. Symlinks the database files
    # into BacFlux's own directory rather than reading the user's path directly,
    # because a directory() output is WIPED before its rule reruns — pointing that
    # at a shared database would delete it for everyone. Same shape as
    # virsorter2_db_local.
    #
    # No index is rebuilt (unlike checkv_db_local): geNomad's files are MMseqs2
    # databases and plain tables, which MMseqs2 reads through symlinks without
    # complaint, and no geNomad step writes into the database directory.
    #
    # A database that is present can still fail in two ways, and 00_common.smk
    # catches both at parse time rather than an hour into the run: a copy that is
    # only partly readable (11 of 27 files mode 0640, seen for real here), and a
    # copy too old for the installed geNomad (which parses its marker metadata
    # positionally and dies on "invalid literal for int()").
    #
    # Produces 07.phages/genomad_db/ — the same path the download rule produces, so
    # genomad_end_to_end is identical either way.
    if GENOMADDB:

        rule genomad_db_local:
            input:
                src = GENOMADDB,
            output:
                genomad_db = directory(GENOMAD_DB_DIR),
            log:
                LOGS + "/genomad_db_local.log"
            shell:
                # The symlink TARGET must be absolute. A relative
                # directories.genomad_db would put a relative target inside the view
                # directory, where it would resolve against the VIEW's location
                # instead of the launch directory — i.e. every link dangles, and
                # geNomad fails on a database that is actually there. `cd ... && pwd`
                # resolves it before any link is made.
                """
                mkdir -p {output.genomad_db}
                src_abs=$(cd "{input.src}" && pwd)
                {{
                  echo "Building a local geNomad database view"
                  echo "  source (read-only): $src_abs"
                  echo "  view:               {output.genomad_db}"
                }} > {log}
                for f in "$src_abs"/*; do
                    ln -sfn "$f" "{output.genomad_db}/$(basename "$f")"
                done
                """


    # ── Virus and plasmid calling in one run ──
    # geNomad scans the finished genome and reports, in a single end-to-end run,
    # which contigs (or contig regions) are viral and which are plasmids. The virus
    # output feeds CheckV here; the plasmid summary feeds the concordance in
    # shared/60_plasmid.smk. Default presets (neither --conservative nor --relaxed)
    # match the standard geNomad→CheckV combination, and --cleanup deletes the
    # intermediates. If a huge input ever needs a RAM lever, `--splits N` caps peak
    # memory at a speed cost; isolate-sized genomes do not need it.
    #
    # Input is FINAL_CONTIGS, the decontaminated assembly (D4: all four modes scan
    # decontaminated contigs, where v1 illumina scanned the pre-decontam
    # contigs_filt.fasta — flagged in the changelog).
    #
    # The output is declared a DIRECTORY so a failed or rerun job is wiped clean and
    # geNomad cannot resume on stale intermediates. geNomad names every output after
    # the input basename (contigs_final), so the two files that matter land at
    #   …/contigs_final_summary/contigs_final_virus.fna           → viral_quality
    #   …/contigs_final_summary/contigs_final_plasmid_summary.tsv → 60_plasmid.smk
    # Both consumers depend on this DIRECTORY and reach inside for their file, which
    # keeps the inner-file edge off the DAG (as 40_annotation.smk does for Bakta).
    #
    # The CLI contract — positional INPUT OUTPUT DATABASE, that summary layout, and
    # a clean run into the directory Snakemake pre-creates — was checked against the
    # real tool on the first geNomad run, 2026-07-25, and matches.
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
        shell:
            """
            genomad end-to-end \
              --cleanup \
              --threads {threads} \
              {input.contigs} \
              {output.genomad_dir} \
              {input.genomad_db} > {log} 2>&1
            """


# ───────────────── VirSorter2 caller (default) ─────────────────
# VirSorter2 2.2.4 (Jan 2023) is itself a Snakemake workflow, and it manages its
# OWN nested conda envs at runtime. That nesting is what crashed every v1.3.1 phage
# run: an ancient transitive `mamba` against a modern `conda`, failing with
# "No module named 'conda._vendor.auxlib'" — the `virsorter_deps_env` saga.
#
# The fix is to switch the nesting OFF entirely and make our single env carry
# everything VS2 needs. Three parts, applied in the rules that follow and in
# envs/virsorter.yaml:
#   1. envs/virsorter.yaml installs VS2's OWN internal dependency list (copied from
#      the envs/vs2.yaml that ships inside the virsorter package) next to virsorter
#      itself, so the one env is self-sufficient.
#   2. `virsorter setup --skip-deps-install` downloads only the database and does
#      NOT build VS2's per-rule nested dependency envs.
#   3. `virsorter run --use-conda-off`, plus exporting the env's own bin onto PATH,
#      makes VS2 subprocesses resolve OUR in-env binaries instead of nested envs.
#
# VERIFIED 2026-07-22 on the v2 illumina validation run, and the answer was NOT the
# comfortable one: the bioconda `virsorter=2.2.4` package does NOT pin the runtime
# tool closure. With only `virsorter` in the env, screed, hmmer, prodigal, last,
# pandas, scikit-learn, numpy, seaborn, imbalanced-learn and ncbi-genome-download
# were ALL absent, and the run died at the first internal rule on "No module named
# 'screed'". Hence part 1 above. If VS2 is ever unpinned from 2.2.4, re-read its
# packaged envs/vs2.yaml and re-sync envs/virsorter.yaml against it.
if PHAGE_CALLER == "virsorter2":

    # ── VirSorter2 database, downloaded here ──
    # VirSorter2 scores contigs against curated viral HMM groups held in a database
    # it fetches once (~10 GB). `virsorter setup` does the download;
    # --skip-deps-install keeps it from building the fragile nested envs (fix part 2
    # above). Produces 07.phages/vs2_db/, consumed by
    # viral_identification_virsorter2.
    #
    # Defined only when BacFlux is the one downloading — mutually exclusive with
    # virsorter2_db_local, same safety reasoning as checkv_db / checkv_db_local: a
    # directory() output is wiped before its rule reruns, so it must never point at
    # a database someone else shares.
    if not VS2DB:

        rule virsorter2_db:
            output:
                vs2_db = directory(VS2_DB_DIR),
            conda:
                "../../envs/virsorter.yaml"
            threads: capped_cpus(4)
            log:
                LOGS + "/virsorter2_db.log"
            shell:
                """
                virsorter setup \
                  -d {output.vs2_db} \
                  -j {threads} \
                  --skip-deps-install > {log} 2>&1
                """


    # ── VirSorter2 database, a copy you already hold ──
    # Defined only when directories.vs2_db is set. Symlinks what `virsorter setup`
    # produced (hmm/, group/, rbs/, Done_all_setup) into BacFlux's own directory
    # rather than reading the user's path directly, so the wipe-on-rerun hazard
    # described at virsorter2_db can never reach it. Nothing is rebuilt here,
    # unlike checkv_db_local — no cross-build incompatibility has turned up for
    # VS2's HMM files, so this is a plain, cheap view.
    #
    # Produces 07.phages/vs2_db/ — the same path the download rule produces, so
    # viral_identification_virsorter2 is identical either way.
    if VS2DB:

        rule virsorter2_db_local:
            input:
                src = VS2DB,
            output:
                vs2_db = directory(VS2_DB_DIR),
            log:
                LOGS + "/virsorter2_db_local.log"
            shell:
                # Absolute target, for the same reason as in genomad_db_local: a
                # relative directories.vs2_db would produce links that resolve
                # against the view directory and therefore dangle.
                """
                mkdir -p {output.vs2_db}
                src_abs=$(cd "{input.src}" && pwd)
                {{
                  echo "Building a local VirSorter2 database view"
                  echo "  source (read-only): $src_abs"
                  echo "  view:               {output.vs2_db}"
                }} > {log}
                for f in "$src_abs"/*; do
                    ln -sfn "$f" "{output.vs2_db}/$(basename "$f")"
                done
                """


    # ── VirSorter2 virus calling ──
    # Identify phages and prophages on the finished genome — a faithful port of v1's
    # `viral_identification` (VS2 half) onto the 00_common.smk paths, with the same
    # viral groups and the same score cutoff. Input is FINAL_CONTIGS (D4 change: v1
    # illumina used the pre-decontam contigs_filt.fasta) plus the VS2 database.
    #
    # min_score 0.5 and the five viral groups come from the VirSorter2 SOP the
    # workflow follows (protocols.io, linked from README 07.phages): a deliberately
    # loose cutoff, taken for maximal sensitivity — CheckV downstream is what grades
    # the result. --keep-original-seq preserves the original sequence of
    # circular and near-fully-viral contigs instead of VS2's trimmed version, so
    # CheckV sees the real contig.
    #
    # Produces 07.phages/virsorter/{sample}/ — key file final-viral-combined.fa,
    # consumed by viral_quality.
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
        shell:
            # Put the env's own bin first so VS2 subprocesses resolve in-env tools
            # (fix part 3 above); without it they hunt for the nested envs that
            # --use-conda-off just told VS2 not to build.
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


# ───────────── CheckV database, a copy you already hold ────────
# Defined only when directories.checkv_db is set. BacFlux does not simply point
# CheckV at that path, for a reason worth remembering: the DIAMOND index inside a
# shared database was built by whatever DIAMOND that site had, DIAMOND's database
# format is versioned, and CheckV then fails deep into the completeness stage with
# "DIAMOND task failed". Shared databases are also usually read-only to the person
# running the workflow. (The longer version of this note is in 00_common.smk, where
# directories.checkv_db is resolved.)
#
# So the versioned folder is recreated locally, every large file symlinked (no
# gigabytes copied), and the DIAMOND index rebuilt with THIS workflow's DIAMOND so
# it is guaranteed compatible. genome_db/ must be a REAL directory, because the new
# index is written into it; hmm_db/ can be one symlink for the whole tree because
# nothing writes there. Disk cost: the index only (~950 MB), against ~6.4 GB for a
# full downloaded copy.
#
# Produces 07.phages/checkv_db/ — the same path the download rule produces, so
# viral_quality is identical either way.
if CHECKVDB:

    rule checkv_db_local:
        input:
            src = CHECKVDB,
        output:
            checkv_db = directory(CHECKV_DB_DIR),
        params:
            db_id = CHECKV_DB_ID,
        conda:
            "../../envs/checkv.yaml"
        threads: capped_cpus(8)
        log:
            LOGS + "/checkv_db_local.log"
        shell:
            """
            # Locate the versioned DB folder inside the user's directory the same
            # way viral_quality does, so both agree on what "the database" is.
            # -L for the same reason as in viral_quality: the user's own database
            # may itself be reached through symlinks, and a plain `find -type f`
            # would silently find nothing.
            src_db=$(dirname "$(dirname "$(find -L {input.src} -type f -path '*/genome_db/checkv_reps.faa' | sort | head -n 1)")")
            dst="{output.checkv_db}/{params.db_id}"

            {{
              echo "Building a local CheckV view"
              echo "  source (read-only): $src_db"
              echo "  view:               $dst"
            }} > {log}

            mkdir -p "$dst/genome_db"

            # Symlink every genome_db file EXCEPT the index: that one we rebuild,
            # because the source copy may be in a DIAMOND format this environment's
            # DIAMOND cannot read.
            for f in "$src_db"/genome_db/*; do
                case "$(basename "$f")" in
                    checkv_reps.dmnd) continue ;;
                esac
                ln -sfn "$f" "$dst/genome_db/$(basename "$f")"
            done

            # Nothing writes into hmm_db/, so one symlink for the whole tree.
            ln -sfn "$src_db/hmm_db" "$dst/hmm_db"
            [ -e "$src_db/README.txt" ] && ln -sfn "$src_db/README.txt" "$dst/README.txt"

            echo "Building DIAMOND index with $(diamond --version 2>&1 | head -n1)" >> {log}
            diamond makedb \
              --in "$dst/genome_db/checkv_reps.faa" \
              --db "$dst/genome_db/checkv_reps" \
              --threads {threads} >> {log} 2>&1
            """


# ─────────────── CheckV database, downloaded here ──────────────
# CheckV grades viral genome completeness and contamination against a reference
# database. Defined only when directories.checkv_db is NOT set, which makes it
# mutually exclusive with checkv_db_local — and that exclusivity is a safety
# requirement, not tidiness. The output is a `directory()`, and Snakemake DELETES a
# directory output before re-running its rule, so both rules write to BacFlux's own
# 07.phages/checkv_db and never to the user's directory. If a rule's output ever
# pointed at a shared database, any re-run trigger (a changed env file, a
# --forcerun, an interrupted job) would wipe it for everyone using it.
#
# With links.checkv_link empty, CheckV downloads its own default database and is
# responsible for its own integrity. With a link set, the .tar.gz AND its .sha256
# are fetched, the hash must match or the rule hard-fails (same discipline as
# dbCAN's cazyme_db_download — a truncated or tampered download must never be
# silently extracted), and the archive is then extracted and its DIAMOND index
# built. The link, the .sha256 URL and the derived folder id (CHECKV_DB_ID) are all
# resolved in 00_common.smk.
#
# CheckV is permissively licensed (LBNL BSD, commercial use OK), so running it on
# every path imposes nothing on the user — unlike geNomad.
#
# The DIAMOND index build is the one real CPU cost here, and capped_cpus(8) matches
# checkv_db_local, which does the identical makedb step.
if not CHECKVDB:

    rule checkv_db:
        output:
            checkv_db = directory(CHECKV_DB_DIR),
        params:
            checkv_link = CHECKV_LINK,
            sha_url = CHECKV_SHA_URL,
            db_id = CHECKV_DB_ID,
            tries = 5,
        conda:
            "../../envs/checkv.yaml"
        threads: capped_cpus(8)
        log:
            LOGS + "/checkv_db.log"
        shell:
            """
            if [ -z "{params.checkv_link}" ]; then
                checkv download_database {output.checkv_db} > {log} 2>&1
            else
                TAR="{output.checkv_db}/{params.db_id}.tar.gz"
                SHA="{output.checkv_db}/{params.db_id}.sha256"

                wget --tries={params.tries} -c {params.checkv_link} -P {output.checkv_db} > {log} 2>&1
                wget --tries={params.tries} -c {params.sha_url} -O "$SHA" >> {log} 2>&1

                # Hard-fail on a mismatch, exactly like dbCAN's cazyme_db_download:
                # a truncated or tampered archive must never reach `tar`/`diamond`.
                expected="$(awk 'NR==1{{print $1}}' "$SHA")"
                actual="$(sha256sum "$TAR" | awk '{{print $1}}')"
                if [ "$expected" != "$actual" ]; then
                    echo "ERROR: checksum mismatch for $TAR" >> {log}
                    echo "  expected: $expected" >> {log}
                    echo "  actual:   $actual" >> {log}
                    exit 1
                fi

                tar -xzvf "$TAR" -C {output.checkv_db} >> {log} 2>&1
                diamond makedb \
                  --in {output.checkv_db}/{params.db_id}/genome_db/checkv_reps.faa \
                  --db {output.checkv_db}/{params.db_id}/genome_db/checkv_reps \
                  --threads {threads} >> {log} 2>&1
            fi
            """


# ─────────────────── Viral quality (CheckV) ────────────────────
# CheckV estimates how complete each predicted viral sequence is and flags host
# contamination in it. The same `checkv end_to_end` command grades geNomad's virus
# FASTA or VirSorter2's, so one rule serves both paths: the caller's output
# DIRECTORY arrives as CHECKV_CALLER_DIR and the shell reaches inside for the virus
# FASTA at CHECKV_VIRAL_REL, both fixed at parse time in the CheckV input selection
# block at the top of this file.
#
# Produces 07.phages/checkv/{sample}/ — key file quality_summary.tsv. This is the
# last phage step, and the file rule all asks for: asking for it is what pulls in
# whichever caller ran.
#
# GRACEFUL EMPTY-INPUT HANDLING (BacFlux "degrade, don't hard-fail" convention): a
# virus-free genome is common, and both callers can emit an empty or absent virus
# FASTA. `checkv end_to_end` errors on an empty input, which would fail the whole
# sample. So we test the FASTA first and, if it holds no sequences, write a
# header-only quality_summary.tsv and skip CheckV instead of crashing the run.
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
    shell:
        # The database directory is resolved at runtime rather than being passed in
        # (ported verbatim from v1): find the single genome_db/checkv_reps.faa,
        # insist there is exactly one, and take its grandparent as the directory
        # CheckV expects. That is what lets checkv_db and checkv_db_local produce
        # differently-named versioned folders without this rule knowing.
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
            # -L (follow symlinks) is REQUIRED, not cosmetic. When the database is
            # a local view of a user-provided copy (rule checkv_db_local), every
            # file in it is a symlink, and a plain `find -type f` reports those as
            # type l and matches NOTHING — the resolver then reports "found 0
            # candidate(s)" for a database that is perfectly fine. With -L, find
            # follows the link and tests the TARGET's type, so real files and
            # symlinked files both match.
            checkv_rep_files=$(find -L {input.checkv_db} -type f -path "*/genome_db/checkv_reps.faa" | sort)
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
              -d "$checkv_db_dir" > {log} 2>&1 || {{
                echo "" >> {log}
                echo "NOTE: CheckV exited non-zero. It grades prophage predictions and nothing downstream reads the grades, so the run continues with a header-only quality_summary.tsv - the same file an isolate with no prophage gets. Check the log above." >> {log}
                printf "contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\twarnings\n" \
                  > {output.checkv_dir}/quality_summary.tsv
              }}
        fi
        """
