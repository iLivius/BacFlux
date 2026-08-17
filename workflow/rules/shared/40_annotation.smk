# Stage 04 annotation: the four tools that turn a finished genome into gene calls
# and functional labels, plus the two database fetches they need. Shared by all
# four entry points — every front end has already converged on one canonical
# assembly (FINAL_CONTIGS, D2) and one genus-composition file, so nothing here
# branches on the sequencing technology.
#
# Chain: contigs_final.fasta + {sample}_composition.txt → Bakta → a DIRECTORY of
# {sample}.faa/.gff3/.gbff/.tsv, which the three analysis rules each reach into
# for the one file they need.
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
#
# Every name used here — DIR_ANNOTATION, FINAL_CONTIGS, COMPOSITION,
# BAKTA_REPLICON_INPUT, BAKTADB, DMNDDB, EGGNOG_DBMEM, EGGNOG_DBMEM_GB,
# ANTISMASHDB, ANTISMASH_DB_DIR, DBCANDB, DBCAN_DB_DIR, DBCAN_SENTINEL,
# DBCAN_LINK, DBCAN_SHA_URL, LOGS, capped_cpus, bakta_locus_tag — is defined once
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
# genome and writes the whole standard annotation set. Telling it which genus the
# isolate belongs to sharpens the call, so the shell derives one from the
# decontamination composition table. Do not trust that genus without reading the
# shell preface: the line it picks is NOT the most abundant genus.
#
# Takes in: FINAL_CONTIGS, the one decontaminated assembly every front end
#           produces (02.assembly/{sample}/contigs_final.fasta); COMPOSITION, the
#           per-genus abundance table select_contigs wrote in
#           shared/10_decontam.smk; and, in the long-read modes only, the replicon
#           table from build_replicons (shared/15_replicons.smk).
# Produces: 04.annotation/bakta/{sample}/, declared as a DIRECTORY (faithful to
#           v1) because Bakta names the files inside it — {sample}.faa, .gff3,
#           .gbff, .tsv and more. Every downstream rule therefore takes its DAG
#           edge on the directory and reaches inside for the file it wants.
# Consumed by: functional_annotation (.faa), secondary_metabolites_analysis
#              (.gbff), cazyme_gene_cluster (.faa + .gff3) and the report module.
#
# (v1 named this `annotation` here and `accurate_annotation` in FastaFlux; D5
#  unifies on `annotation`. v1 message: "--- Bakta: Genome annotation. ---")
rule annotation:
    input:
        # The finished genome — one canonical hand-off from any front end (D2).
        contigs = FINAL_CONTIGS,
        # Per-genus abundance table, a CROSS-STAGE edge: select_contigs writes it
        # in shared/10_decontam.smk. Consumer and producer both name the single
        # COMPOSITION constant in 00_common.smk (alongside FINAL_CONTIGS), so the
        # two can never drift onto different paths.
        abund = COMPOSITION,
        # Replicon table, LONG-READ MODES ONLY. It tells Bakta which contigs are
        # circular, and that matters because Pyrodigal is then allowed to call
        # genes across the origin of a closed replicon — typically a handful of
        # genes at position 1 of a chromosome, often including dnaA itself. Built
        # by build_replicons in shared/15_replicons.smk from Flye's circularity
        # call plus dnaapler's start-gene marker.
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
    log:
        LOGS + "/annotation_{sample}.log"
    shell:
        # Genus hint for Bakta: drop the "no-hit" line from the composition table,
        # sort what is left, and pass the surviving genus name as --genus.
        #
        # Read that sort carefully: it does NOT rank by abundance any more. It did
        # while the table held one bare number per genus ("Bacillus: 0.30").
        # write_composition() in scripts/10_decontam/select_contigs_by_taxonomy.py now writes
        # two labelled figures instead — "Bacillus: bases 0.29  contigs 0.30" — so
        # the key after the ':' starts with the word "bases", `-n` scores every
        # line 0, and GNU sort breaks the resulting all-way tie on the whole line,
        # reversed. What survives is the alphabetically LAST genus. Worked through
        # on a four-genus composition file — 0.40 of the DNA on the first genus,
        # then 0.29, 0.16 and 0.11 — where the name handed to Bakta is the third
        # of them, not the genus holding most of the genome.
        #
        # It stays a hint, not a filter: a wrong genus costs annotation accuracy,
        # never contigs or reads, and nothing in the run announces it.
        #
        # Two deliberate v2 changes over the plain v1 illumina form:
        #
        # 1. Skip "no-hit". The composition table counts EVERY contig in the
        #    BlobTools table, and contigs with no informative BLAST assignment are
        #    counted under the literal genus "no-hit". Without the grep, a
        #    poorly-placed or heavily-contaminated sample can put "no-hit" on the
        #    winning line and Bakta runs with `--genus no-hit`.
        # 2. Run Bakta EXACTLY ONCE, with the genus hint only if one was found.
        #    v1 wrapped the bakta call in a `for` loop over that single genus, so
        #    if no genus survived the filter the loop body never ran and Bakta was
        #    never invoked — the rule then failed on missing output. Here the genus
        #    is resolved first and the hint added conditionally, so a sample with
        #    no usable genus is still annotated, just without the hint. (This
        #    matches v1 BacFluxL+, which likewise passed no --genus/--species when
        #    it had no clear winner. Reconciling all four modes onto L+'s stricter
        #    kept-contigs-only rule is logged as a follow-up.)
        """
        genus=$(grep -v '^no-hit:' {input.abund} | sort -t':' -k2 -nr | cut -d':' -f1 | sed -n '1p')

        if [ -n "$genus" ]; then
            taxon_args="--genus $genus --species sp."
            echo "Annotating {wildcards.sample} with genus hint: $genus" > {log}
        else
            taxon_args=""
            echo "No usable genus for {wildcards.sample} (composition empty or all no-hit); annotating without a genus hint." > {log}
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
#              QUAST, CheckM, GTDB-Tk and Bakta, and nothing else).
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

    rule secondary_metabolites_db:
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
# Consumed by: nobody — you. A terminal product: rule all asks for it, and
#              MultiQC does NOT read it (multiqc_qc_inputs() in
#              shared/90_report.smk aggregates fastp, NanoPlot, Qualimap,
#              QUAST, CheckM, GTDB-Tk and Bakta, and nothing else).
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
# conda: NONE, inherited from v1. It runs in the environment Snakemake was
# launched from and needs wget, tar, sha256sum and awk to be available there.
# Adding a small wget/coreutils env would be new behaviour, so it is left alone.
#
# The doubled braces in the shell are not a typo: Snakemake reads {…} as one of
# its own placeholders, so a literal brace for awk has to be written {{…}} to
# reach the shell intact.
if not DBCANDB:

    rule cazyme_db_download:
        output:
            dbcan_db = directory(DBCAN_DB_DIR),
            dbcan_verified = DBCAN_SENTINEL,
        params:
            dbcan_db_url = DBCAN_LINK,
            sha_url = DBCAN_SHA_URL,
        log:
            LOGS + "/cazyme_db_download.log"
        shell:
            """
            mkdir -p "{output.dbcan_db}"

            TAR="{output.dbcan_db}/$(basename "{params.dbcan_db_url}")"
            SHA="{output.dbcan_db}/$(basename "{params.sha_url}")"

            wget -O "$TAR" "{params.dbcan_db_url}" > "{log}" 2>&1
            wget -O "$SHA" "{params.sha_url}" >> "{log}" 2>&1

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
#              QUAST, CheckM, GTDB-Tk and Bakta, and nothing else).
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
