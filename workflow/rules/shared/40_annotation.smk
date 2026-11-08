# Stage 04 annotation: the four tools that turn a finished genome into gene calls
# and functional labels, plus the two database fetches they need, and two opt-in
# layers on top of Bakta, with the databases they need: structure-based names
# from Baktfold, and the EcoFoldDB trait profile. Shared by all four entry points
# — every front end has already converged on one canonical assembly
# (FINAL_CONTIGS) and one kept-genus file, so nothing here branches on the
# sequencing technology.
#
# Chain: contigs_final.fasta + {sample}_kept_genus.txt → Bakta → a DIRECTORY of
# {sample}.faa/.gff3/.gbff/.tsv/.json, which the three analysis rules, and the
# opt-in layers when they are on, each reach into for the files they need.
#
# annotation                     : Bakta. Gene calling plus structural and
#                                  functional annotation. Everything below it
#                                  depends on its output directory.
# functional_annotation          : eggNOG-mapper over Bakta's .faa — orthologous
#                                  group per protein, and the COG/GO/KEGG/EC
#                                  labels that group carries.
# secondary_metabolites_db       : fetches the antiSMASH reference database once.
# secondary_metabolites_db_local : same directory, built from a copy the user
#                                  already holds. Mutually exclusive with the
#                                  rule above.
# secondary_metabolites_analysis : antiSMASH over Bakta's .gbff — biosynthetic
#                                  gene clusters.
# cazyme_db_download             : fetches and checksum-verifies dbCAN.
# cazyme_db_local                : same directory from a user-held copy.
#                                  Mutually exclusive with the rule above.
# cazyme_gene_cluster            : run_dbcan over Bakta's .faa + .gff3 — CAZymes,
#                                  the clusters they form, and each cluster's
#                                  predicted substrate.
# structure_annotation_db        : fetches Baktfold's database once (opt-in,
#                                  baktfold.run or traits.run: the trait profile
#                                  uses its ProstT5 model).
# structure_annotation_db_local  : same directory from a user-held copy.
#                                  Mutually exclusive with the rule above.
# structure_annotation           : Baktfold over Bakta's .json — names for the
#                                  proteins Bakta called hypothetical, found by
#                                  predicted structure (opt-in, baktfold.run).
# trait_db                       : fetches and MD5-checks EcoFoldDB once (opt-in,
#                                  traits.run).
# trait_db_local                 : same directory from a user-held copy.
#                                  Mutually exclusive with the rule above.
# trait_structure_query          : ProstT5 3Di strings for every protein of up
#                                  to 4000 aa in Bakta's .faa, and a Foldseek
#                                  query database from them (opt-in,
#                                  traits.run).
# trait_structure_search         : Foldseek search of that database against
#                                  EcoFoldDB (opt-in, traits.run).
# trait_profile                  : the trait calls, from the structure hits and
#                                  the Bakta and Baktfold names, with antiSMASH's
#                                  region types reported beside them (opt-in,
#                                  traits.run).
#
# Every name used here — DIR_ANNOTATION, FINAL_CONTIGS, KEPT_GENUS,
# BAKTA_REPLICON_INPUT, BAKTADB, DMNDDB, EGGNOG_DBMEM, EGGNOG_DBMEM_GB,
# ANTISMASHDB, ANTISMASH_DB_DIR, DBCANDB, DBCAN_DB_DIR, DBCAN_SENTINEL,
# DBCAN_LINK, DBCAN_SHA_URL, BAKTFOLD_RUN, BAKTFOLD_FAST, BAKTFOLD_DB_NEEDED,
# BAKTFOLDDB, BAKTFOLD_DB_DIR, TRAITS_RUN, TRAITS_MAX_GENE_GAP,
# TRAITS_MAX_PROTEIN_LENGTH, TRAITS_FOLDSEEK_SENSITIVITY, ECOFOLDDB,
# ECOFOLDDB_DB_DIR, ECOFOLDDB_LINK, ECOFOLDDB_MD5, DIR_TRAITS,
# TRAIT_STRUCTURE_DIR, TRAIT_HITS, TRAIT_GENES, TRAIT_PROFILE,
# TRAIT_PROFILE_AUDIT, TRAIT_PROFILE_SCRIPT, TRAIT_RULES_TSV, TRAIT_GENES_TSV,
# HEAVY_MEM_MB, LONG_RUNTIME_MINUTES, LOGS, capped_cpus,
# bakta_locus_tag — is defined once
# in 00_common.smk and inherited through the Snakefile include order. No path,
# stage number or database version is re-derived here.
#
# conda: env paths are written relative to THIS file. Snakemake resolves a rule's
# conda env against the file that DEFINES the rule, and this one sits at
# workflow/rules/shared/, so "../../envs/x.yaml" climbs shared/ → rules/ →
# workflow/ and lands on workflow/envs/x.yaml (the single shared env copy).
#
# Resources: cpu-bound rules declare Snakemake's built-in `threads: capped_cpus(N)`
# and use {threads} in the shell. The BUILT-IN keyword — rather than a custom
# `resources: cpus` — is what makes `--cores N` enforce the limit, so a plain
# `snakemake --cores N` is safe on its own. The one exception is eggNOG's mem_gb,
# a gigabyte figure rather than a core count, which Snakemake schedules against
# only when the launch line also passes `--resources mem_gb=N`.


# ──────────────────────── Genome annotation (Bakta) ────────────
# Bakta calls genes (CDS, tRNA, rRNA, ncRNA, CRISPR arrays, …) on the finished
# genome and writes the whole standard annotation set. It is also given the
# isolate's genus, which it writes as the organism name in its GenBank, EMBL and
# GFF3 files and on its genome plot. The genus does not change which genes are
# called or how they are annotated.
#
# Takes in: FINAL_CONTIGS, the one decontaminated assembly every front end
#           produces (02.assembly/{sample}/contigs_final.fasta); KEPT_GENUS, the
#           one-line genus file select_contigs wrote in shared/10_decontam.smk;
#           and, in the long-read modes only, the replicon table from
#           build_replicons (shared/15_replicons.smk).
# Produces: 04.annotation/bakta/{sample}/, declared as a DIRECTORY (faithful to
#           v1) because Bakta names the files inside it — {sample}.faa, .gff3,
#           .gbff, .tsv and more. Every downstream rule therefore takes its DAG
#           edge on the directory and reaches inside for the file it wants.
# Consumed by: functional_annotation (.faa), secondary_metabolites_analysis
#              (.gbff), cazyme_gene_cluster (.faa + .gff3) and the report module;
#              when switched on, also structure_annotation (.json),
#              trait_structure_query (.faa) and trait_profile (.tsv + .faa).
#
# (v1 named this `annotation` here and `accurate_annotation` in FastaFlux; v2
#  uses `annotation` in every mode.
#  v1 message: "--- Bakta: Genome annotation. ---")
rule annotation:
    input:
        # The finished genome — one canonical hand-off from any front end.
        contigs = FINAL_CONTIGS,
        # The genus to name the genome after, a CROSS-STAGE edge: select_contigs
        # writes it in shared/10_decontam.smk as the genus holding most of the
        # KEPT contigs. Consumer and producer both name the single KEPT_GENUS
        # constant in 00_common.smk, so the two can never drift onto different
        # paths.
        genus = KEPT_GENUS,
        # Replicon table, LONG-READ MODES ONLY. It tells Bakta which contigs are
        # circular, and that matters because Pyrodigal is then allowed to call a
        # gene that crosses the start of a circular sequence whole, for example
        # on a replicon dnaapler could not rotate. Built by build_replicons in
        # shared/15_replicons.smk from Flye's circularity call plus dnaapler's
        # start-gene marker.
        #
        # In illumina and contigs mode BAKTA_REPLICON_INPUT is an EMPTY LIST, so
        # there is no DAG edge and no producer is needed; Snakemake renders it as
        # an empty string in the shell, the `[ -s ]` test below is false, and the
        # flag is simply absent. One rule body, valid in all four modes.
        replicons = BAKTA_REPLICON_INPUT,
    output:
        bakta_dir = directory(DIR_ANNOTATION + "/bakta/{sample}"),
    params:
        bakta_db = BAKTADB,
        # Sample names can hold characters Bakta rejects in a locus tag;
        # bakta_locus_tag() in 00_common.smk trims one down to a safe 24-character
        # prefix. The sample name the user sees is untouched.
        locus_tag = lambda wc: bakta_locus_tag(wc.sample),
    conda:
        "../../envs/bakta.yaml"
    # Bakta gains little past ~24 threads; capped_cpus takes the lower of that and
    # resources.threads.
    threads: capped_cpus(24)
    resources:
        # Cluster requests (00_common.smk §3a). A single computer reads mem_mb
        # only when the launch line passes --resources mem_mb=N.
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    log:
        LOGS + "/annotation_{sample}.log"
    shell:
        # The organism name for Bakta: the genus in the one-line KEPT_GENUS file,
        # i.e. the genus holding most of the contigs decontamination kept. The
        # file is empty when no kept contig has a genus (for example, every kept
        # contig is one BLAST could not place); Bakta then runs without a genus,
        # so the sample is still annotated. Either way the first line of the log
        # says which genus was used.
        """
        genus=$(head -n 1 {input.genus})

        if [ -n "$genus" ]; then
            taxon_args="--genus $genus --species sp."
            echo "Annotating {wildcards.sample} as genus: $genus" > {log}
        else
            taxon_args=""
            echo "No genus for {wildcards.sample}: no kept contig has a genus, so Bakta runs without one." > {log}
        fi

        # Replicon table, same idiom as the genus hint above. The -s test ("exists
        # and is not empty") does double duty: it is false in the short-read modes,
        # where the input is bound to an empty list and renders as "", AND it is
        # the second guard against handing Bakta an EMPTY file, which it treats as
        # a fatal format error rather than as "no information" — the first is in
        # build_bakta_replicons.py, which warns when it writes one.
        replicon_args=""
        if [ -s "{input.replicons}" ]; then
            replicon_args="--replicons {input.replicons}"
            echo "Using replicon table: {input.replicons}" >> {log}
        fi

        bakta \
          --db {params.bakta_db} \
          --verbose \
          $taxon_args \
          $replicon_args \
          --strain {wildcards.sample} \
          --translation-table 11 \
          --min-contig-length 500 \
          --locus-tag {params.locus_tag} \
          --prefix {wildcards.sample} \
          --keep-contig-headers \
          --output {output.bakta_dir} \
          --threads {threads} \
          --force {input.contigs} >> {log} 2>&1
        """


# ──────────────── Functional annotation (eggNOG-mapper) ────────
# eggNOG-mapper assigns each protein Bakta predicted to an orthologous group and
# attaches what that group is known to do — COG category, GO terms, KEGG KO and
# pathway, EC number. It is the "what do these genes DO" layer on top of Bakta's
# "where the genes are".
#
# Takes in: the Bakta output DIRECTORY as the DAG edge; the shell reads the
#           protein FASTA {sample}.faa from inside it.
# Does:     a DIAMOND-mode search of those proteins against the eggNOG database
#           (DMNDDB, the user-supplied directories.eggnog_db).
# Produces: 04.annotation/eggnog/{sample}/ with the {sample}.emapper.* tables.
# Consumed by: nobody — you. A terminal product: rule all asks for it, and
#              MultiQC does NOT read it (multiqc_qc_inputs() in
#              shared/90_report.smk aggregates fastp, NanoPlot, Qualimap,
#              QUAST, CheckM, GTDB-Tk and Bakta, plus the trait profile when
#              traits.run is on, and nothing else).
#
# This is the slow tail of a BacFlux run. parameters.eggnog.dbmem loads the 39 GB
# eggnog.db into memory instead of reading it off disk, at ~42 GB per concurrent
# job; it is off by default — see the params and resources blocks below.
#
# (v1 message: "--- EggNOG: Functional annotation. ---")
rule functional_annotation:
    input:
        bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
    output:
        # emapper's scratch space, wiped by Snakemake when the rule ends. It is a
        # SIBLING of the kept output (…/{sample}_tmp), never a subdirectory of it,
        # so temp() cleanup never reaches inside a directory output that survives.
        # (00_common.smk applies the same house rule to CARD_TARBALL.)
        temp_dir = temp(directory(DIR_ANNOTATION + "/eggnog/{sample}_tmp")),
        eggnog_dir = directory(DIR_ANNOTATION + "/eggnog/{sample}"),
    params:
        dmnd_db = DMNDDB,
        # "--dbmem" when parameters.eggnog.dbmem is on, empty string otherwise.
        # Resolved once in 00_common.smk, which also refuses at parse time to turn
        # it on when resources.ram_gb cannot fit one job; see EGGNOG_DBMEM there.
        dbmem_flag = "--dbmem" if EGGNOG_DBMEM else "",
    conda:
        "../../envs/eggnog-mapper.yaml"
    threads: capped_cpus(24)
    resources:
        # RAM this job needs, in GB: 0 unless --dbmem is on, then the ~42 GB that
        # holds the 39 GB eggnog.db in memory. Like Qualimap's java_mem, mem_gb is
        # a GIGABYTE figure and not a core count, and Snakemake schedules against
        # it only when the launch line passes `--resources mem_gb=N` — 00_common
        # prints that exact flag, with the user's own ram_gb filled in, whenever
        # --dbmem is switched on.
        mem_gb = EGGNOG_DBMEM_GB if EGGNOG_DBMEM else 0,
        # Cluster requests (00_common.smk §3a). A single computer reads mem_mb
        # only when the launch line passes --resources mem_mb=N.
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    log:
        LOGS + "/functional_annotation_{sample}.log"
    shell:
        # emapper.py expects both directories to exist already, hence the mkdir.
        # {params.dbmem_flag} is an empty string on the default path and adds
        # nothing to the command line.
        """
        mkdir -p {output.temp_dir} {output.eggnog_dir}

        emapper.py \
          -i {input.bakta_dir}/{wildcards.sample}.faa \
          --output_dir {output.eggnog_dir} \
          --cpu {threads} \
          -m diamond \
          --data_dir {params.dmnd_db} \
          {params.dbmem_flag} \
          --output {wildcards.sample} \
          --temp_dir {output.temp_dir} \
          --override > {log} 2>&1
        """


# ──────────────────────── antiSMASH reference database ─────────
# antiSMASH recognises biosynthetic gene clusters by matching curated cluster and
# domain models, which live in a reference database it has to fetch once. A pure
# download: no per-sample input, and every sample's antiSMASH run waits on it.
#
# Defined ONLY when directories.antismash_db is empty. Setting that key defines
# secondary_metabolites_db_local below instead, and the two are never both in the
# DAG. Keeping them exclusive is a safety requirement rather than tidiness: the
# output is a `directory()`, Snakemake DELETES a directory output before re-running
# its rule, and if that path were ever the user's own database, any re-run trigger
# — a changed env file, a --forcerun, an interrupted job — would wipe it for
# everyone sharing it. Both rules therefore write only into BacFlux's own
# 04.annotation/antismash/databases. Same reasoning as checkv_db /
# checkv_db_local in shared/70_phage.smk.
#
# Path collision that is not one: ANTISMASH_DB_DIR sits at antismash/databases,
# right beside the antismash/{sample} directories the analysis rule writes. The
# project-wide `wildcard_constraints: sample=...` in 00_common.smk pins {sample} to
# the sample names actually discovered, so "databases" can never match it and this
# path can only ever be produced by this rule. (v1 leaned on Snakemake's
# concrete-beats-wildcard tie-break; the constraint makes it explicit.)
#
# (v1 message: "--- antiSMASH: database download. ---")
if not ANTISMASHDB:

    # Where it runs: on the computer that runs Snakemake (localrule: True), so on
    # a cluster the download never lands on a compute node without internet.
    rule secondary_metabolites_db:
        localrule: True
        output:
            antismash_db = directory(ANTISMASH_DB_DIR),
        conda:
            "../../envs/antismash.yaml"
        log:
            LOGS + "/secondary_metabolites_database.log"
        shell:
            """
            download-antismash-databases \
              --database-dir {output.antismash_db} > {log} 2>&1
            """


# Defined ONLY when directories.antismash_db is set, and then in place of the
# download rule above. Symlinks the database's top-level entries (clusterblast/,
# pfam/, ...) into BacFlux's own directory rather than pointing antiSMASH at the
# user's path, for the directory()-wipe reason spelled out above. Nothing is
# rebuilt: unlike CheckV's DIAMOND index, no cross-build incompatibility has shown
# up in antiSMASH's database — this exact database was reused as-is, by plain copy,
# across two real screening batches with zero errors.
if ANTISMASHDB:

    rule secondary_metabolites_db_local:
        input:
            src = ANTISMASHDB,
        output:
            antismash_db = directory(ANTISMASH_DB_DIR),
        log:
            LOGS + "/secondary_metabolites_database_local.log"
        shell:
            """
            mkdir -p {output.antismash_db}
            {{
              echo "Building a local antiSMASH database view"
              echo "  source (read-only): {input.src}"
              echo "  view:               {output.antismash_db}"
            }} > {log}
            for f in "{input.src}"/*; do
                ln -sfn "$f" "{output.antismash_db}/$(basename "$f")"
            done
            """


# ──────────────────── Secondary metabolites (antiSMASH) ────────
# Scans one genome for biosynthetic gene clusters — antibiotics, siderophores and
# the rest of the secondary metabolism — and reports each cluster with its type
# and its closest known relative.
#
# Takes in: two DAG edges, the antiSMASH database directory from whichever of the
#           two rules above is defined, and the Bakta output DIRECTORY. The shell
#           reads Bakta's annotated GenBank file {sample}.gbff from inside it.
# Does:     runs antiSMASH with --genefinding-tool none, meaning DO NOT re-predict
#           genes: reuse the calls Bakta already made, which travel inside the
#           .gbff. That is why this rule consumes .gbff and not the bare .fna, and
#           it keeps the cluster coordinates on the same gene set as every other
#           annotation output.
# Produces: 04.annotation/antismash/{sample}/.
# Consumed by: you; and, when traits.run is on, trait_profile, which reads its
#              {sample}.json for the region types it reports beside a trait call
#              (they never change the call). Rule all asks for it, and
#              MultiQC does NOT read it (multiqc_qc_inputs() in
#              shared/90_report.smk aggregates fastp, NanoPlot, Qualimap,
#              QUAST, CheckM, GTDB-Tk and Bakta, plus the trait profile when
#              traits.run is on, and nothing else).
#
# (v1 message: "--- antiSMASH: secondary metabolite annotation. ---")
rule secondary_metabolites_analysis:
    input:
        antismash_db = ANTISMASH_DB_DIR,
        bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
    output:
        antismash_dir = directory(DIR_ANNOTATION + "/antismash/{sample}"),
    params:
        # Rule-local literals rather than paths or config, so they stay inline.
        taxon = 'bacteria',
        genefinding_tool = 'none',
    conda:
        "../../envs/antismash.yaml"
    threads: capped_cpus(24)
    resources:
        # Cluster requests (00_common.smk §3a). A single computer reads mem_mb
        # only when the launch line passes --resources mem_mb=N.
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    log:
        LOGS + "/secondary_metabolites_{sample}.log"
    shell:
        """
        antismash \
          --output-dir {output.antismash_dir} \
          --output-basename {wildcards.sample} \
          --databases {input.antismash_db} \
          --taxon {params.taxon} \
          --genefinding-tool {params.genefinding_tool} \
          --cpus {threads} \
          {input.bakta_dir}/{wildcards.sample}.gbff > {log} 2>&1
        """


# ──────────────────────── dbCAN reference database ─────────────
# dbCAN recognises carbohydrate-active enzymes against HMM, DIAMOND and dbCAN-sub
# reference sets shipped together in one versioned tarball. Downloads that
# tarball and — the point of the rule — VERIFIES it before anything trusts it.
#
# Does: downloads the .tar.gz and its .sha256; reads the expected hash out of the
#       .sha256 (first field of line 1); computes the tarball's actual hash; and
#       hard-fails on `test "$expected" = "$actual"` if the two differ, so a
#       truncated or tampered download can never be extracted quietly. Only on a
#       match does it extract (flattening the top-level directory) and write the
#       verified hash into the sentinel.
# Produces:
#   dbcan_db       = DBCAN_DB_DIR, the extracted database. Its version folder name
#                    is DERIVED FROM links.dbcan_link in 00_common.smk
#                    (DBCAN_DB_ID), so the folder can never disagree with the
#                    configured URL — v1 hard-coded "dbcan_db_v5.1.2" and would
#                    have mismatched silently against any other link. It sits
#                    beside the per-sample dbcan/{sample} directories and cannot
#                    be confused with one, for the wildcard_constraints reason
#                    given under secondary_metabolites_db.
#   dbcan_verified = DBCAN_SENTINEL, a marker file holding the verified checksum.
#                    Its existence is the "present AND integrity-checked" gate.
# Consumed by: cazyme_gene_cluster, which takes its DAG edge on the SENTINEL
#              rather than the directory, so dbCAN cannot start against a database
#              that failed its check.
#
# Defined ONLY when directories.dbcan_db is empty; setting that key defines
# cazyme_db_local below instead. Mutually exclusive for the directory()-wipe
# reason given under secondary_metabolites_db.
#
# Environment: envs/download.yaml, for wget and tar (with the gzip tar calls);
# sha256sum and awk are part of every Linux base system.
# Where it runs: on the computer that runs Snakemake (localrule: True), so on a
# cluster the download never lands on a compute node without internet.
#
# The doubled braces in the shell are not a typo: Snakemake reads {…} as one of
# its own placeholders, so a literal brace for awk has to be written {{…}} to
# reach the shell intact.
if not DBCANDB:

    rule cazyme_db_download:
        localrule: True
        output:
            dbcan_db = directory(DBCAN_DB_DIR),
            dbcan_verified = DBCAN_SENTINEL,
        params:
            dbcan_db_url = DBCAN_LINK,
            sha_url = DBCAN_SHA_URL,
        conda:
            "../../envs/download.yaml"
        log:
            LOGS + "/cazyme_db_download.log"
        shell:
            """
            mkdir -p "{output.dbcan_db}"

            TAR="{output.dbcan_db}/$(basename "{params.dbcan_db_url}")"
            SHA="{output.dbcan_db}/$(basename "{params.sha_url}")"

            # --progress=dot:giga prints one dot per megabyte, so the log gets a
            # few dozen lines for the whole database instead of tens of thousands.
            wget --progress=dot:giga -O "$TAR" "{params.dbcan_db_url}" > "{log}" 2>&1
            wget --progress=dot:giga -O "$SHA" "{params.sha_url}" >> "{log}" 2>&1

            # verify checksum
            expected="$(awk 'NR==1{{print $1}}' "$SHA")"
            actual="$(sha256sum "$TAR" | awk '{{print $1}}')"
            test "$expected" = "$actual"

            # extract (flatten top-level directory)
            tar -xzf "$TAR" -C "{output.dbcan_db}" --strip-components=1 >> "{log}" 2>&1

            # create sentinel containing verified checksum
            echo "$actual" > "{output.dbcan_verified}"
            """


# Defined ONLY when directories.dbcan_db is set, and then in place of the download
# rule above. Symlinks the database's top-level files (dbCAN.hmm, CAZy.dmnd, ...)
# into BacFlux's own directory for the directory()-wipe reason given under
# secondary_metabolites_db. Nothing is rebuilt: no cross-build incompatibility has
# shown up in dbCAN's HMM/DIAMOND files, and this exact database was reused as-is,
# by plain copy, across two real screening batches with zero errors.
#
# The sentinel it writes says in words that this copy is UNVERIFIED, as opposed to
# cazyme_db_download's, which holds a checksum this run actually computed.
# cazyme_gene_cluster tests only that the sentinel exists, never its contents, so
# a user-supplied database opens the same gate — the integrity claim in that rule
# holds for the download path alone.
if DBCANDB:

    rule cazyme_db_local:
        input:
            src = DBCANDB,
        output:
            dbcan_db = directory(DBCAN_DB_DIR),
            dbcan_verified = DBCAN_SENTINEL,
        log:
            LOGS + "/cazyme_db_local.log"
        shell:
            """
            mkdir -p "{output.dbcan_db}"
            {{
              echo "Building a local dbCAN database view"
              echo "  source (read-only, unverified): {input.src}"
              echo "  view:                           {output.dbcan_db}"
            }} > {log}
            for f in "{input.src}"/*; do
                ln -sfn "$f" "{output.dbcan_db}/$(basename "$f")"
            done
            echo "local copy, not independently checksummed" > {output.dbcan_verified}
            """


# ──────────────────── CAZymes and gene clusters (dbCAN) ────────
# Finds the carbohydrate-active enzymes, groups neighbouring ones into CAZyme Gene
# Clusters, and predicts the substrate each cluster acts on — in short, which
# sugars and polysaccharides this genome can process.
#
# Takes in: two DAG edges —
#   bakta_dir      = the Bakta output DIRECTORY; the shell reads BOTH {sample}.faa
#                    (proteins) and {sample}.gff3 (gene coordinates) from inside.
#   dbcan_verified = DBCAN_SENTINEL, written by whichever of cazyme_db_download /
#                    cazyme_db_local is defined. Only its existence is tested, so
#                    it means "a database is in place", and additionally "its
#                    checksum matched" on the download path.
# Does: four run_dbcan sub-commands in sequence, all writing into ONE directory —
#         1. CAZyme_annotation    HMM + DIAMOND + dbCAN-sub over the proteins
#         2. gff_process          put those calls back on the genome via the GFF
#         3. cgc_finder           group adjacent CAZymes into clusters
#         4. substrate_prediction predict what each cluster acts on
#       Step 2 is passed `--gff_type prodigal` for a Bakta GFF3, carried over
#       unchanged from v1; Bakta calls its CDS with Pyrodigal.
# Produces: 04.annotation/dbcan/{sample}/.
# Consumed by: nobody — you. A terminal product: rule all asks for it, and
#              MultiQC does NOT read it (multiqc_qc_inputs() in
#              shared/90_report.smk aggregates fastp, NanoPlot, Qualimap,
#              QUAST, CheckM, GTDB-Tk and Bakta, plus the trait profile when
#              traits.run is on, and nothing else).
#
# run_dbcan is pinned to 5.1.2 (workflow/envs/dbcan.yaml) against a version-pinned
# Zenodo copy of its database, and the pin is load-bearing: 5.2.x changed both the
# subcommand arguments and the output schema, and on the same input 5.2.9 called
# 205 CAZyme genes where 5.1.2 called 315 — see README, links.dbcan_link.
#
# (v1 message: "--- dbCAN: Finding CAZyme gene clusters. ---")
rule cazyme_gene_cluster:
    input:
        bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
        dbcan_verified = DBCAN_SENTINEL,
    output:
        dbcan_dir = directory(DIR_ANNOTATION + "/dbcan/{sample}"),
    params:
        dbcan_db = DBCAN_DB_DIR,
    conda:
        "../../envs/dbcan.yaml"
    threads: capped_cpus(24)
    log:
        LOGS + "/cazyme_{sample}.log"
    shell:
        """
        # CAZyme annotation of protein sequences
        run_dbcan CAZyme_annotation \
          --input_raw_data {input.bakta_dir}/{wildcards.sample}.faa \
          --output_dir {output.dbcan_dir} \
          --db_dir {params.dbcan_db} \
          --mode protein \
          --threads {threads} \
          --methods hmm \
          --methods diamond \
          --methods dbCANsub > {log} 2>&1

        # CAZyme Gene Cluster (CGC) Annotation
        run_dbcan gff_process \
          --input_gff {input.bakta_dir}/{wildcards.sample}.gff3 \
          --output_dir {output.dbcan_dir} \
          --db_dir {params.dbcan_db} \
          --gff_type prodigal \
          --threads {threads} >> {log} 2>&1

        # CAZyme Gene Cluster (CGC) Identification
        run_dbcan cgc_finder \
          --output_dir {output.dbcan_dir} >> {log} 2>&1

        # CGC Substrate Prediction
        run_dbcan substrate_prediction \
          --output_dir {output.dbcan_dir} \
          --db_dir {params.dbcan_db} >> {log} 2>&1
        """


# ──────────── Structure-based names (Baktfold), opt-in ────────────
# Bakta names a protein by sequence similarity. When nothing in its databases is
# similar enough, it writes "hypothetical protein". Protein structure is more
# conserved than sequence, so Baktfold looks again at exactly those proteins by
# predicted 3D structure: the ProstT5 protein language model turns each sequence
# into a string of Foldseek's 3Di structural alphabet, and Foldseek searches
# SwissProt, AlphaFold clusters, PDB and CATH with it. 3Di has one letter per
# residue, describing how that residue sits against the residue closest to it in
# space in the folded protein. ProstT5 predicts these letters from the sequence
# alone, with no full structure prediction, and Foldseek searches them like a
# sequence. A hit gives the protein a name; the database and the alignment
# behind it are reported beside it, because a structure match is weaker evidence
# than a sequence match.
#
# Rule structure_annotation exists only when baktfold.run is true. The two
# database rules also exist when only traits.run is on (BAKTFOLD_DB_NEEDED),
# because the trait profile below makes its 3Di strings with the ProstT5 model
# kept in this database. Like antiSMASH above, the database comes from exactly
# one of two rules, never both: a download when directories.baktfold_db is
# empty, or a symlink view of the user's own copy when it is set. Both write only
# into BacFlux's own 04.annotation/baktfold/databases, for the directory()-wipe
# reason given at secondary_metabolites_db.

if BAKTFOLD_DB_NEEDED and not BAKTFOLDDB:

    # ── structure_annotation_db — Baktfold's database, downloaded ──
    # Fetches Baktfold's database and the ProstT5 model (about 22 GB) once per
    # output folder, with `baktfold install`.
    # Takes in: nothing from the workflow.
    # Produces: 04.annotation/baktfold/databases: the Foldseek databases of
    #           SwissProt, AlphaFold clusters, PDB and CATH, Baktfold's own SQLite
    #           database and the ProstT5 model.
    # Consumed by: structure_annotation (all of it, when baktfold.run is on) and
    #              trait_structure_query (the ProstT5 model, when traits.run is
    #              on).
    #
    # Defined only when directories.baktfold_db is empty; setting that key
    # defines structure_annotation_db_local below instead.
    # Environment: envs/baktfold.yaml, for `baktfold install`.
    # Where it runs: on the computer that runs Snakemake (localrule: True), so on
    # a cluster the download never lands on a compute node without internet;
    # download_all asks for it too.
    rule structure_annotation_db:
        localrule: True
        output:
            baktfold_db = directory(BAKTFOLD_DB_DIR),
        conda:
            "../../envs/baktfold.yaml"
        log:
            LOGS + "/structure_annotation_database.log"
        shell:
            # --threads: Baktfold uses it to make the download faster.
            """
            baktfold install --database {output.baktfold_db} --threads 4 > {log} 2>&1
            """


if BAKTFOLD_DB_NEEDED and BAKTFOLDDB:

    # ── structure_annotation_db_local — Baktfold's database, from your copy ──
    # Defined only when directories.baktfold_db is set, in place of the download
    # above. Links each top-level entry of the user's database (the Foldseek
    # databases and the ProstT5 model folder) into BacFlux's own directory, so no
    # rule ever writes to, or deletes, the user's copy.
    # Takes in: the user's database (directories.baktfold_db), checked at
    #           start-up to hold swissprot.dbtype (00_common.smk).
    # Produces: 04.annotation/baktfold/databases, the same path the download
    #           rule writes, so the rules that read it are the same either way.
    # Consumed by: structure_annotation and trait_structure_query, as above.
    # No environment: it uses only echo, mkdir, ln and basename, which every
    # Linux base system has.
    rule structure_annotation_db_local:
        input:
            src = BAKTFOLDDB,
        output:
            baktfold_db = directory(BAKTFOLD_DB_DIR),
        log:
            LOGS + "/structure_annotation_database_local.log"
        shell:
            """
            mkdir -p {output.baktfold_db}
            {{
              echo "Building a local Baktfold database view"
              echo "  source (read-only): {input.src}"
              echo "  view:               {output.baktfold_db}"
            }} > {log}
            for f in "{input.src}"/*; do
                ln -sfn "$f" "{output.baktfold_db}/$(basename "$f")"
            done
            """


if BAKTFOLD_RUN:

    # ── structure_annotation — names for Bakta's hypothetical proteins ──
    # Takes in: the Bakta output DIRECTORY as the DAG edge (the shell reads
    #           Bakta's {sample}.json from inside it: Baktfold needs the full
    #           annotation, not just the proteins, to write Bakta-format files
    #           back), and the Baktfold database directory from whichever of the
    #           two rules above is defined.
    # Does:     baktfold run on the CPU. Only the proteins Bakta called
    #           hypothetical are searched; every other annotation is copied from
    #           Bakta unchanged. --fast (baktfold.fast) skips the AlphaFold-
    #           clusters search.
    # Produces: 04.annotation/baktfold/{sample}/: the same file set as Bakta
    #           ({sample}.tsv, .gff3, .gbff, .json, .faa, ...) with Baktfold's
    #           names filled in, plus {sample}.inference.tsv, one row per
    #           searched protein with its top hit in each database searched
    #           (four; three with baktfold.fast), and one
    #           baktfold_<database>_tophit.tsv per database searched, with the
    #           alignment details (these names always start with baktfold_, not
    #           with the sample name --prefix gives the other files).
    #           A directory() output, as for Bakta: Baktfold names the files
    #           inside it.
    # Consumed by: you, and trait_profile when traits.run is on too ({sample}.tsv:
    #              Baktfold's names count there as evidence for a gene). Bakta's
    #              annotation stays the one eggNOG, antiSMASH, dbCAN, the AMR stage
    #              and the mobilome module read, so switching this on changes none
    #              of their results.
    #
    # Environment: envs/baktfold.yaml, Baktfold 0.3.0, which brings ProstT5 (run
    # by PyTorch on the CPU) and Foldseek. Threads: up to 24, the cap the other
    # annotation rules use; the ProstT5 step and each Foldseek search use them in
    # turn.
    rule structure_annotation:
        input:
            bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
            baktfold_db = BAKTFOLD_DB_DIR,
        output:
            baktfold_dir = directory(DIR_ANNOTATION + "/baktfold/{sample}"),
        params:
            # "--fast" when baktfold.fast is on, an empty string otherwise.
            fast_flag = "--fast" if BAKTFOLD_FAST else "",
        conda:
            "../../envs/baktfold.yaml"
        threads: capped_cpus(24)
        resources:
            # Cluster requests (00_common.smk §3a). Baktfold needs well under
            # 16 GB on a bacterial genome (about 8 GB at its peak), so it asks
            # for no more; a smaller ram_gb wins. The time limit is the one for
            # long jobs (runtime_hours), as for Bakta.
            mem_mb = min(HEAVY_MEM_MB, 16 * 1024),
            runtime = LONG_RUNTIME_MINUTES,
        log:
            LOGS + "/structure_annotation_{sample}.log"
        shell:
            # The ProstT5 step runs in PyTorch, which by default starts a thread
            # per core of the computer; OMP_NUM_THREADS and MKL_NUM_THREADS hold it
            # to {threads}, the same budget Foldseek gets through --threads. The
            # ProstT5 model is already in the database folder: HF_HUB_OFFLINE=1
            # stops the job asking the Hugging Face hub for updates, so it runs on
            # a compute node without internet.
            """
            export OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} HF_HUB_OFFLINE=1
            baktfold run \
              --input {input.bakta_dir}/{wildcards.sample}.json \
              --output {output.baktfold_dir} \
              --prefix {wildcards.sample} \
              --database {input.baktfold_db} \
              --threads {threads} \
              --cpu {params.fast_flag} \
              --force > {log} 2>&1
            """


# ─────── Trait profile from predicted structure (EcoFoldDB), opt-in ───────
# What it is for: a profile of the traits that make a strain of interest for
# agriculture or bioremediation — nitrogen fixation, phosphate solubilisation,
# siderophores, plant hormones, stress protectants, breakdown of aromatic
# pollutants and more: one trait per row of trait_rules_v1.tsv. Each trait is
# called from the set of genes that table lays down: for most traits several
# genes that must all be present, often side by side as in an operon; for a few,
# one marker gene. trait_genes_v1.tsv, the gene dictionary, says how each gene
# is recognised.
#
# A gene counts on a protein when a name says so (Bakta's, or Baktfold's for a
# protein Bakta could only call hypothetical, when baktfold.run is on) or, for
# the genes the gene dictionary trusts on structure alone (its s_alone column),
# when the protein's predicted structure matches one of the gene's EcoFoldDB
# entries.
# EcoFoldDB (Ghaly et al. 2025) is a curated database of protein structures with
# known ecological functions. ProstT5, a protein language model, turns every
# protein of up to 4000 amino acids into a string of Foldseek's 3Di structural
# alphabet (defined in the Baktfold section above), and Foldseek searches
# EcoFoldDB with those strings.
#
# The settings are those of EcoFoldDB's authors: that length limit, Foldseek's
# sensitivity -s at 8, and their filter on each protein's best hit, applied by
# trait_profile.py. The best hit must be one of the database's catalogue
# entries, not one of the background structures it holds beside them, with an
# e-value (the number of hits this good expected by chance in a search of this
# size) below 1e-10 and an alignment covering more than 80% of both the protein
# (query coverage) and the EcoFoldDB structure (target coverage). -s trades
# speed for reach: a higher value finds more distant matches but runs longer;
# Foldseek's own default is 9.5.
#
# One step differs from the authors' script: the 3Di strings come from Baktfold's
# ProstT5 (baktfold proteins-predict), not from the ProstT5 built into
# Foldseek that their script uses. Baktfold's ProstT5 comes with
# envs/baktfold.yaml and with Baktfold's database, which BacFlux already uses
# for the Baktfold layer, so no new environment or model is needed. The price
# is Baktfold's database (about 22 GB), which holds the model and is needed
# even when baktfold.run is off.
#
# All five rules exist only when traits.run is true:
#   trait_db / trait_db_local : EcoFoldDB, downloaded once per output folder or
#                               linked from the user's copy (never both).
#   trait_structure_query     : the 3Di strings and a Foldseek query database
#                               made from them; the slow step.
#   trait_structure_search    : the Foldseek search against EcoFoldDB.
#   trait_profile             : the trait calls.
# Rule traits_summary in shared/90_report.smk then sets every sample's calls side
# by side. The ProstT5 model sits in Baktfold's database, which is why the two
# Baktfold database rules above also exist when only traits.run is on.
#
# Why the 3Di step and the search are separate rules: the 3Di step takes about
# 1 h 50 min per genome on 16 cores and needs only Bakta's proteins and the
# ProstT5 model, never EcoFoldDB. So a new copy of EcoFoldDB (a re-download, a
# changed link, another directories.ecofolddb_db) re-runs only the search, about
# 2 minutes per genome, and the calls after it.
#
# Environments: each rule names its own below. envs/baktfold.yaml,
# envs/download.yaml and envs/python.yaml each list in their header the rules
# that use them, and those lists leave these rules out on purpose: Snakemake
# re-builds an environment, and re-runs every rule that uses it in every finished
# output folder, whenever a byte of its file changes, comments included.
#
# Outputs: 04.annotation/traits/{sample}/ and 04.annotation/traits/databases.
# Apart from the MultiQC section built from the summary, nothing else in the
# workflow reads them, so switching the layer on adds jobs and changes no other
# result. User documentation: docs/analysis/traits.md.

if TRAITS_RUN and not ECOFOLDDB:

    # ── trait_db — EcoFoldDB, downloaded and checked ──
    # Takes in: nothing from the workflow; links.ecofolddb_link and
    #           links.ecofolddb_md5 from the config.
    # Does:     downloads the EcoFoldDB archive; checks its MD5 when one is set
    #           (a mismatch deletes the archive and stops the job, so a truncated
    #           download never reaches tar) and says in the log when none is;
    #           unpacks it; removes the archive; and stops unless the two files
    #           every later step needs are in place.
    # Produces: 04.annotation/traits/databases (about 810 MB): the Foldseek
    #           database EcoFoldDB and EcoFoldDB_descriptions.txt, the table that
    #           names the gene and pathway behind each entry.
    # Consumed by: trait_structure_search (the database) and trait_profile (the
    #              descriptions table).
    #
    # Defined only when directories.ecofolddb_db is empty; setting that key
    # defines trait_db_local below instead, for the directory()-wipe reason given
    # at secondary_metabolites_db. Modelled on rule genomad_db in
    # shared/70_phage.smk, which fetches a Zenodo archive the same way, minus its
    # fall-back to the tool's own downloader: EcoFoldDB has none.
    #
    # Environment: envs/download.yaml, for wget and tar (with the gzip tar calls);
    # md5sum, awk and ls are part of every Linux base system.
    # Where it runs: on the computer that runs Snakemake (localrule: True), so on a
    # cluster the download never lands on a compute node without internet;
    # download_all asks for it too.
    rule trait_db:
        localrule: True
        output:
            ecofolddb = directory(ECOFOLDDB_DB_DIR),
        params:
            link = ECOFOLDDB_LINK,
            md5 = ECOFOLDDB_MD5,
            # How many times wget tries the download before it gives up, the
            # same number rule genomad_db uses.
            tries = 5,
        conda:
            "../../envs/download.yaml"
        log:
            LOGS + "/trait_db.log"
        shell:
            # The archive is saved inside the output folder, so if the job fails,
            # Snakemake's removal of that folder takes the archive with it.
            # --progress=dot:giga prints one dot per megabyte, 32 to a line, so
            # the log gets about twenty lines instead of tens of thousands.
            # The archive holds one top-level folder, EcoFoldDB_v2.0/, with no
            # leading ./ in its member names; --strip-components=1 drops that
            # folder, so the files land directly in {output.ecofolddb}.
            # --no-same-permissions applies your umask to the modes stored in the
            # archive. tar does that anyway for an ordinary user; the option
            # matters only when the job runs as root (in a container, for
            # example), where tar would otherwise keep the archive's modes as
            # they are. With any other layout the two files end up somewhere
            # else, or are not unpacked at all; the check at the end then stops
            # the job and lists in the log what was unpacked.
            """
            mkdir -p {output.ecofolddb}
            TAR="{output.ecofolddb}/$(basename '{params.link}')"
            echo "Fetching the EcoFoldDB database:" > {log}
            echo "  {params.link}" >> {log}
            wget --tries={params.tries} --progress=dot:giga -c "{params.link}" -O "$TAR" >> {log} 2>&1

            if [ -n "{params.md5}" ]; then
                actual="$(md5sum "$TAR" | awk '{{print $1}}')"
                if [ "{params.md5}" != "$actual" ]; then
                    echo "ERROR: checksum mismatch for $TAR" >> {log}
                    echo "  expected (links.ecofolddb_md5): {params.md5}" >> {log}
                    echo "  actual:                         $actual" >> {log}
                    echo "If you changed links.ecofolddb_link, update links.ecofolddb_md5" >> {log}
                    echo "to match the new file (Zenodo shows the MD5 next to it), or" >> {log}
                    echo "clear it to download without verification." >> {log}
                    rm -f "$TAR"
                    exit 1
                fi
                echo "Checksum OK ({params.md5})." >> {log}
            else
                echo "links.ecofolddb_md5 is empty: extracting WITHOUT verifying the download." >> {log}
            fi

            tar -xzf "$TAR" -C {output.ecofolddb} --strip-components=1 --no-same-permissions >> {log} 2>&1
            rm -f "$TAR"

            for needed in EcoFoldDB.dbtype EcoFoldDB_descriptions.txt; do
                if [ ! -f "{output.ecofolddb}/$needed" ]; then
                    echo "ERROR: $needed is not at the top of {output.ecofolddb} after unpacking." >> {log}
                    echo "The archive does not have the layout of EcoFoldDB_v2.0.tar.gz. Unpacked:" >> {log}
                    ls -lR {output.ecofolddb} >> {log} 2>&1 || true
                    exit 1
                fi
            done
            echo "EcoFoldDB ready at {output.ecofolddb}." >> {log}
            """


if TRAITS_RUN and ECOFOLDDB:

    # ── trait_db_local — EcoFoldDB, from your copy ──
    # Defined only when directories.ecofolddb_db is set, in place of the download
    # above. Links each top-level entry of the user's EcoFoldDB into BacFlux's own
    # 04.annotation/traits/databases, so no rule ever writes to, or deletes, the
    # user's copy. Its entry names are checked against the gene dictionary
    # (trait_genes_v1.tsv) when the workflow starts (00_common.smk).
    # Takes in: the user's EcoFoldDB folder (directories.ecofolddb_db).
    # Produces: 04.annotation/traits/databases, the same path the download rule
    #           writes, so the rules that read it are the same either way.
    # Consumed by: trait_structure_search and trait_profile, as above.
    # No environment: it uses only shell built-ins (cd, pwd, echo) and mkdir, ln
    # and basename, which every Linux base system has.
    rule trait_db_local:
        input:
            src = ECOFOLDDB,
        output:
            ecofolddb = directory(ECOFOLDDB_DB_DIR),
        log:
            LOGS + "/trait_db_local.log"
        shell:
            # The link TARGET must be absolute: a relative directories.ecofolddb_db
            # would otherwise resolve against the view folder instead of the
            # launch folder, and every link would dangle. `cd ... && pwd` resolves
            # it first, as rule genomad_db_local in shared/70_phage.smk does.
            """
            mkdir -p {output.ecofolddb}
            src_abs=$(cd "{input.src}" && pwd)
            {{
              echo "Building a local EcoFoldDB database view"
              echo "  source (read-only): $src_abs"
              echo "  view:               {output.ecofolddb}"
            }} > {log}
            for f in "$src_abs"/*; do
                ln -sfn "$f" "{output.ecofolddb}/$(basename "$f")"
            done
            """


if TRAITS_RUN:

    # ── trait_structure_query — a 3Di string for every protein (the slow step) ──
    # Takes in: the Bakta output DIRECTORY as the DAG edge (the shell reads
    #           {sample}.faa, every protein of the genome, from inside it), and
    #           the Baktfold database directory, for its ProstT5 model, from
    #           whichever Baktfold database rule above is defined. Never
    #           EcoFoldDB: see the section comment.
    # Does:     1. keeps the proteins of at most TRAITS_MAX_PROTEIN_LENGTH (4000)
    #              amino acids, as EcoFoldDB's own method does; trait_profile
    #              lists every longer protein as not searched;
    #           2. turns each kept protein into a 3Di string, Foldseek's
    #              structural alphabet, with ProstT5 on the CPU (baktfold
    #              proteins-predict). Baktfold writes a residue as X, masked
    #              for the search, when ProstT5 is less than 25% sure of its
    #              letter (Baktfold's default);
    #           3. builds a Foldseek query database from the amino-acid and 3Di
    #              strings (baktfold createdb).
    # Produces: 04.annotation/traits/{sample}/structure/: {sample}_searched.faa
    #           (the kept proteins), prostt5/ (the 3Di strings, with ProstT5's
    #           confidence for each residue) and querydb/ (the query database,
    #           named {sample}, and Baktfold's logs); about 25 MB per genome. A
    #           directory() output, because Baktfold names the files it writes
    #           there. It is kept after the run, so a later re-run of the search
    #           starts from it instead of repeating this step.
    # Consumed by: trait_structure_search.
    #
    # Environment: envs/baktfold.yaml, for Baktfold's ProstT5 and its Foldseek.
    # Threads, memory and offline settings as rule structure_annotation, which
    # runs ProstT5 too, but only on Bakta's hypothetical proteins, a small share
    # of the genome. Here every protein of up to 4000 aa goes through ProstT5, so
    # this is the slow step: about 1 h 50 min on 16 cores for a genome of about
    # 6,000 proteins, with a peak of about 10 GB of memory. The time grows with
    # the number of proteins, and with other jobs on the same computer.
    rule trait_structure_query:
        input:
            bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
            baktfold_db = BAKTFOLD_DB_DIR,
        output:
            query_dir = directory(TRAIT_STRUCTURE_DIR),
        params:
            # 4000 amino acids, EcoFoldDB's limit; the awk filter in the shell
            # applies it.
            max_length = TRAITS_MAX_PROTEIN_LENGTH,
        conda:
            "../../envs/baktfold.yaml"
        threads: capped_cpus(24)
        resources:
            # Cluster requests (00_common.smk §3a). ProstT5 needs under 16 GB
            # on a bacterial genome (about 10 GB at its peak), so the job asks
            # for no more; a smaller ram_gb wins. The time limit is the one for
            # long jobs (runtime_hours): this is the slow step.
            mem_mb = min(HEAVY_MEM_MB, 16 * 1024),
            runtime = LONG_RUNTIME_MINUTES,
        log:
            LOGS + "/trait_structure_query_{sample}.log"
        shell:
            # The awk program reads Bakta's protein FASTA one record at a time:
            # it joins a record's sequence lines first, so it works whether or
            # not the FASTA is wrapped, and writes the record (sequence on one
            # line) only when the sequence is at most {params.max_length} letters
            # long. Its braces are doubled because Snakemake would otherwise read
            # {...} as one of its own placeholders.
            #
            # PyTorch, which runs ProstT5, by default starts a thread per core of
            # the computer: OMP_NUM_THREADS and MKL_NUM_THREADS hold it to
            # {threads}. The model is already in the database folder:
            # HF_HUB_OFFLINE=1 stops the job asking the Hugging Face hub for
            # updates, so it runs on a compute node without internet. -f
            # (overwrite) matters only when these commands are re-run by hand in
            # the same folder: Snakemake deletes structure/ before every attempt.
            """
            mkdir -p {output.query_dir}
            export OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} HF_HUB_OFFLINE=1

            awk -v max_len={params.max_length} '
                /^>/ {{
                    if (header != "" && length(seq) <= max_len) {{ print header; print seq }}
                    header = $0
                    seq = ""
                    next
                }}
                {{ seq = seq $0 }}
                END {{ if (header != "" && length(seq) <= max_len) {{ print header; print seq }} }}
            ' "{input.bakta_dir}/{wildcards.sample}.faa" > "{output.query_dir}/{wildcards.sample}_searched.faa"
            echo "Proteins of at most {params.max_length} aa, searched by structure: $(grep -c '^>' "{output.query_dir}/{wildcards.sample}_searched.faa")" > {log}

            baktfold proteins-predict \
              -i "{output.query_dir}/{wildcards.sample}_searched.faa" \
              -o "{output.query_dir}/prostt5" \
              -p {wildcards.sample} \
              -d {input.baktfold_db:q} \
              -t {threads} \
              --cpu -f >> {log} 2>&1

            baktfold createdb \
              --fasta-aa "{output.query_dir}/{wildcards.sample}_searched.faa" \
              --fasta-3di "{output.query_dir}/prostt5/{wildcards.sample}_3di.fasta" \
              -o "{output.query_dir}/querydb" \
              -p {wildcards.sample} \
              -t {threads} \
              -f >> {log} 2>&1
            """

    # ── trait_structure_search — the Foldseek search against EcoFoldDB ──
    # Takes in: the query folder of trait_structure_query (the shell names the
    #           query database querydb/{sample} inside it), and EcoFoldDB, from
    #           whichever of trait_db and trait_db_local is defined.
    # Does:     a Foldseek search of every searched protein against EcoFoldDB at
    #           sensitivity -s TRAITS_FOLDSEEK_SENSITIVITY (8, EcoFoldDB's own
    #           setting, a little below Foldseek's default of 9.5), then
    #           a conversion of the results to a table of query, target, e-value,
    #           query coverage and target coverage. Foldseek keeps its default
    #           limits (e-value 10, up to 1000 hits per protein), so the table
    #           holds far more than the hits that count. trait_profile needs
    #           them all: it takes each protein's best hit across the whole of
    #           EcoFoldDB, background structures included, before it applies
    #           EcoFoldDB's filter. For a best hit that passes, it also looks at
    #           the best hit to a different catalogue gene (another entry of the
    #           same gene does not count), or to a background structure: when
    #           that hit's e-value is at most ten times the best hit's, the
    #           borderline column of {sample}_trait_genes.tsv flags a near tie.
    # Produces: {sample}_ecofolddb_hits.tsv, about 3 million rows (170 MB) per
    #           genome. A temp() file: only trait_profile reads it, so Snakemake
    #           deletes it once the calls are made. A later re-run of the calls
    #           first repeats this search (about 2 minutes) from the query folder,
    #           which is kept; the 3Di step is not repeated.
    # Consumed by: trait_profile.
    #
    # Environment: envs/baktfold.yaml, for the Foldseek Baktfold uses, so the
    # query database and the search come from one Foldseek version.
    # Threads: up to 24, as for the 3Di step. No cluster requests (mem_mb,
    # runtime): on 16 cores the search takes about 2 minutes and about 2 GB of
    # memory.
    rule trait_structure_search:
        input:
            query_dir = TRAIT_STRUCTURE_DIR,
            ecofolddb = ECOFOLDDB_DB_DIR,
        output:
            hits = temp(TRAIT_HITS),
            # Foldseek's scratch space and its results database (about 150 MB
            # before conversion), deleted by Snakemake when the job ends, like
            # eggNOG's. A sibling of the kept outputs, never inside one of them.
            work_dir = temp(directory(DIR_TRAITS + "/{sample}/search_tmp")),
        params:
            # Foldseek's -s, 8 as in EcoFoldDB's own method.
            sensitivity = TRAITS_FOLDSEEK_SENSITIVITY,
        conda:
            "../../envs/baktfold.yaml"
        threads: capped_cpus(24)
        log:
            LOGS + "/trait_structure_search_{sample}.log"
        shell:
            # The results database is removed as soon as the table is written,
            # so with --notemp, which keeps temp() outputs, search_tmp/ holds
            # only Foldseek's small scratch folder (the hit table is kept as
            # well).
            """
            mkdir -p {output.work_dir}
            foldseek search \
              "{input.query_dir}/querydb/{wildcards.sample}" \
              "{input.ecofolddb}/EcoFoldDB" \
              "{output.work_dir}/results" \
              "{output.work_dir}/tmp" \
              -s {params.sensitivity} \
              --threads {threads} \
              --remove-tmp-files 1 > {log} 2>&1

            foldseek convertalis \
              "{input.query_dir}/querydb/{wildcards.sample}" \
              "{input.ecofolddb}/EcoFoldDB" \
              "{output.work_dir}/results" \
              {output.hits:q} \
              --format-output query,target,evalue,qcov,tcov \
              --threads {threads} >> {log} 2>&1

            rm -f "{output.work_dir}"/results*
            """

    # ── trait_profile — the trait calls, with the evidence behind each ──
    # Takes in: the Bakta output DIRECTORY (the shell reads {sample}.tsv, every
    #           feature with its position, gene symbol, product and KEGG and EC
    #           cross-references, and {sample}.faa, for each protein's length);
    #           the hit table of trait_structure_search; EcoFoldDB's descriptions
    #           table, which names the gene and pathway behind each hit; the two
    #           trait tables; the antiSMASH output DIRECTORY ({sample}.json, the
    #           region types); and, only when baktfold.run is on, the Baktfold
    #           output DIRECTORY ({sample}.tsv, the names Baktfold gave Bakta's
    #           hypothetical proteins).
    # Does:     trait_profile.py checks the two trait tables, keeps the structure
    #           hits that pass EcoFoldDB's published filter, decides which gene
    #           counts on which protein, and calls each trait from its gene sets,
    #           checking that the genes a trait needs side by side form a chain
    #           on one contig, each at most traits.max_gene_gap genes from the
    #           one before it. Its comments and docs/analysis/traits.md give the
    #           full logic and every reason word it writes.
    # Produces: three tables and their provenance in 04.annotation/traits/{sample}/:
    #             {sample}_trait_profile.tsv        one row per trait: called,
    #                                               partial or not called, with
    #                                               the genes behind it;
    #             {sample}_trait_genes.tsv          one row per protein and gene
    #                                               with any evidence, and whether
    #                                               and why it counted;
    #             {sample}_trait_profile_audit.tsv  one row per trait, with the
    #                                               reason for its decision.
    #             {sample}_trait_provenance.json    SHA-256 fingerprints of the
    #                                               two catalogue tables, the
    #                                               EcoFoldDB catalogue and the
    #                                               three tables above, plus the
    #                                               calling parameters, so the
    #                                               traits workbook can check that
    #                                               the catalogue it quotes is the
    #                                               one these calls were made with.
    # Consumed by: traits_summary in shared/90_report.smk (the profiles and the
    #              provenance records); traits_workbook there (all four files);
    #              the user, for the three tables.
    #
    # The two trait tables are inputs on purpose. The calls must always match
    # the tables they were made from, so editing either table re-makes them
    # (seconds per sample), and with them the summary and the MultiQC report.
    # Because the hit table is temp(), such a re-run first repeats the 2-minute
    # search; the 3Di step is not repeated. A fresh clone of BacFlux gives the
    # tables new file dates, which re-makes the calls the same way in an output
    # folder made before. The script itself is a parameter, as for every helper
    # script, so editing trait_profile.py re-makes nothing by itself.
    #
    # Environment: envs/python.yaml, the pinned Python of the stdlib helper rules.
    # One thread and no cluster requests: the script reads the hit table line by
    # line, keeping only each protein's best hit and its closest rival, and takes
    # seconds.
    rule trait_profile:
        input:
            bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
            hits = TRAIT_HITS,
            ecofolddb = ECOFOLDDB_DB_DIR,
            rules = TRAIT_RULES_TSV,
            genes = TRAIT_GENES_TSV,
            antismash_dir = DIR_ANNOTATION + "/antismash/{sample}",
            # Baktfold's annotation. The ** line (Python's way to add the
            # entries of a dictionary here) adds the input baktfold_dir only
            # when baktfold.run is on. With it off, the rule has no such input,
            # so it neither asks for nor waits on a Baktfold job.
            **({"baktfold_dir": DIR_ANNOTATION + "/baktfold/{sample}"} if BAKTFOLD_RUN else {}),
        output:
            genes = TRAIT_GENES,
            profile = TRAIT_PROFILE,
            audit = TRAIT_PROFILE_AUDIT,
            provenance = TRAIT_PROVENANCE,
        params:
            script = TRAIT_PROFILE_SCRIPT,
            # traits.max_gene_gap, for the side-by-side test.
            max_gene_gap = TRAITS_MAX_GENE_GAP,
            # The length limit rule trait_structure_query applied, so the script
            # can list each longer protein as not searched, rather than as a
            # protein the search found nothing for.
            max_protein_length = TRAITS_MAX_PROTEIN_LENGTH,
            # "--baktfold-tsv <path>" when baktfold.run is on, an empty string
            # otherwise, so the option is absent from the command line. Written
            # as a function of the job's wildcards (w.sample is the sample
            # name), so the path is built for each sample, in the same form as
            # icescan_flag in shared/80_mobilome.smk.
            baktfold_flag = lambda w: (
                "--baktfold-tsv " + DIR_ANNOTATION + "/baktfold/" + w.sample + "/" + w.sample + ".tsv"
                if BAKTFOLD_RUN else ""
            ),
        conda:
            "../../envs/python.yaml"
        log:
            LOGS + "/trait_profile_{sample}.log"
        shell:
            """
            python {params.script:q} \
              --sample {wildcards.sample:q} \
              --bakta-tsv "{input.bakta_dir}/{wildcards.sample}.tsv" \
              --bakta-faa "{input.bakta_dir}/{wildcards.sample}.faa" \
              --hits {input.hits:q} \
              --catalogue "{input.ecofolddb}/EcoFoldDB_descriptions.txt" \
              --rules {input.rules:q} \
              --genes {input.genes:q} \
              --antismash-json "{input.antismash_dir}/{wildcards.sample}.json" \
              {params.baktfold_flag} \
              --max-gene-gap {params.max_gene_gap} \
              --max-protein-length {params.max_protein_length} \
              --out-genes {output.genes:q} \
              --out-profile {output.profile:q} \
              --out-audit {output.audit:q} \
              --out-provenance {output.provenance:q} > {log} 2>&1
            """
