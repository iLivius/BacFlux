# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — Stage 04 annotation module (rules/shared/40_annotation.smk)
#
# This module runs the four annotation tools on every sample's finished genome,
# regardless of which front end (illumina / nanopore / hybrid / contigs) built
# that genome. It lives in the SHARED tail of the pipeline: by this point every
# mode has converged on the one canonical assembly file (FINAL_CONTIGS, D2) and
# on the genus-composition file, so the rules here never need to know the input
# technology — the paths and behaviour are identical in all four modes.
#
# Data flow through this module (top to bottom):
#
#   contigs_final.fasta ─┐
#                        ├─► annotation (Bakta) ─► 04.annotation/bakta/{sample}/
#   {sample}_composition ┘        │   (a DIRECTORY of {sample}.faa/.gff3/.gbff/...)
#                                 │
#                 ┌───────────────┼───────────────────────┐
#                 ▼               ▼                        ▼
#   functional_annotation  secondary_metabolites_analysis  cazyme_gene_cluster
#      (eggNOG-mapper)         (antiSMASH)                   (dbCAN)
#      reads .faa              reads .gbff                   reads .faa + .gff3
#
# Two tools need a reference database fetched first; those downloads are their
# own rules (secondary_metabolites_db, cazyme_db_download) and enter the DAG only
# because the analysis rules depend on them.
#
# Everything this module references — DIR_ANNOTATION, DIR_ASSEMBLY, FINAL_CONTIGS,
# LOGS, BAKTADB, DMNDDB, DBCAN_DB_DIR, DBCAN_SENTINEL, DBCAN_LINK, DBCAN_SHA_URL,
# capped_cpus, bakta_locus_tag — is defined once in 00_common.smk and inherited
# through the Snakefile include order. No path, stage number, or database version
# is re-derived here.
#
# conda: paths are written relative to THIS .smk file. Snakemake resolves a
# rule's conda env relative to the file that DEFINES the rule, and this file sits
# at workflow/rules/shared/, so "../../envs/x.yaml" climbs shared/ -> rules/ ->
# workflow/ and lands on workflow/envs/x.yaml (the single shared env copy).
#
# Resource convention: cpu-bound rules request `resources: cpus = capped_cpus(N)`
# and refer to `{resources.cpus}` in the shell — the same mechanism v1 used
# (`min(CPUS, N)`), now going through the 00_common helper. We do NOT switch to
# Snakemake's `threads:` keyword.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: annotation — whole-genome structural + functional annotation (Bakta) ──
# Biology: Bakta calls genes (CDS, tRNA, rRNA, ncRNA, CRISPR arrays, …) on the
# finished genome and writes a rich, standardised annotation set. Giving it the
# genus this isolate most likely belongs to sharpens the annotation, so we feed
# it the top genus from the decontamination composition table.
#
# Takes in:
#   contigs = FINAL_CONTIGS — the one decontaminated, finished assembly that every
#             mode's front end produces (02.assembly/{sample}/contigs_final.fasta).
#   abund   = the "Genus:percentage" composition table written by the
#             decontamination selector. The shell picks the single most abundant
#             genus from it to pass to Bakta's --genus (see the shell notes).
#
# Does: for the top genus in the composition file, run Bakta once over the
#       contigs. The for-loop is only a portable way to capture that one genus
#       name into $i; it iterates exactly once (sed -n '1p' keeps just line 1).
#
# Produces: 04.annotation/bakta/{sample}/ — declared as a DIRECTORY, because Bakta
#           names the files inside it ({sample}.faa/.gff3/.gbff/.tsv/…). We keep
#           only the directory as the output (faithful to v1): the DAG edge to
#           every downstream rule is therefore "depends on the Bakta directory",
#           and each downstream rule reaches inside for the specific file it needs.
#
# Consumed by: functional_annotation (.faa), secondary_metabolites_analysis
#              (.gbff), cazyme_gene_cluster (.faa + .gff3), and the report module.
#
# (v1 rule name was `annotation` here and `accurate_annotation` in FastaFlux; D5
#  unifies on `annotation`. v1 message: "--- Bakta: Genome annotation. ---")
rule annotation:
    input:
        # The finished genome — single canonical hand-off from any front end (D2).
        contigs = FINAL_CONTIGS,
        # Genus:abundance table from the decontamination selector. This is a
        # CROSS-STAGE edge: the file is PRODUCED by `select_contigs` in the future
        # shared/10_decontam.smk. Both this consumer and that producer reference the
        # single COMPOSITION constant defined in 00_common (alongside FINAL_CONTIGS),
        # so the two can never drift onto different paths — the same single-source
        # pattern the assembly hand-off uses.
        abund = COMPOSITION,
    output:
        # Whole directory, not per-file outputs (faithful to v1); Bakta writes
        # {sample}.faa/.gff3/.gbff/.tsv/... inside it.
        bakta_dir = directory(DIR_ANNOTATION + "/bakta/{sample}"),
    params:
        bakta_db = BAKTADB,
        # Short, filesystem-safe locus-tag prefix derived from the sample name.
        locus_tag = lambda wc: bakta_locus_tag(wc.sample),
    conda:
        "../../envs/bakta.yaml"
    resources:
        # Bakta gains little past ~24 threads; cap via the shared helper.
        cpus = capped_cpus(24)
    log:
        LOGS + "/annotation_{sample}.log"
    priority: 5
    shell:
        # Genus selection preserved byte-for-byte from v1: sort the composition
        # lines by the numeric abundance after the ':' (descending), take the
        # genus name before the ':', keep only the first (most abundant) line.
        """
        for i in $(cat {input.abund} | sort -t':' -k2 -nr | cut -d':' -f1 | sed -n '1p'); do \
        bakta \
          --db {params.bakta_db} \
          --verbose \
          --genus $i \
          --species sp. \
          --strain {wildcards.sample} \
          --translation-table 11 \
          --min-contig-length 500 \
          --locus-tag {params.locus_tag} \
          --prefix {wildcards.sample} \
          --keep-contig-headers \
          --output {output.bakta_dir} \
          --threads {resources.cpus} \
          --force {input.contigs}; \
        done > {log} 2>&1
        """


# ── Rule: functional_annotation — orthology / functional labels (eggNOG-mapper) ─
# Biology: eggNOG-mapper takes the proteins Bakta predicted and assigns each one
# to an orthologous group, then attaches the functional annotation that group
# carries (COG category, GO terms, KEGG KO/pathway, EC number, …). This is the
# "what do these genes DO" layer on top of Bakta's "where the genes are".
#
# Takes in: the Bakta output DIRECTORY (the DAG edge). In the shell it reads the
#           protein FASTA {sample}.faa from inside that directory.
# Does: a DIAMOND-mode search of those proteins against the eggNOG DB (DMNDDB).
# Produces: 04.annotation/eggnog/{sample}/ with the {sample}.emapper.* tables. A
#           scratch eggnog_tmp/ subdirectory is declared temp() so Snakemake wipes
#           it once the rule finishes.
# Consumed by: the report module (a terminal annotation product).
#
# (v1 message: "--- EggNOG: Functional annotation. ---")
rule functional_annotation:
    input:
        # Depends on the Bakta directory; the .faa is read from within it (shell).
        bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
    output:
        # emapper scratch space, deleted automatically (temp) when the rule ends.
        # Kept as a SIBLING of the persisted output (…/{sample}_tmp), not nested
        # inside it, so temp() cleanup never has to remove a subdirectory of a
        # directory output we keep — which avoids surprising cleanup behaviour.
        temp_dir = temp(directory(DIR_ANNOTATION + "/eggnog/{sample}_tmp")),
        eggnog_dir = directory(DIR_ANNOTATION + "/eggnog/{sample}"),
    params:
        dmnd_db = DMNDDB,
    conda:
        "../../envs/eggnog-mapper.yaml"
    resources:
        cpus = capped_cpus(24)
    log:
        LOGS + "/functional_annotation_{sample}.log"
    priority: 4
    shell:
        # Both output dirs are created up front (emapper expects them to exist).
        """
        mkdir -p {output.temp_dir} {output.eggnog_dir}

        emapper.py \
          -i {input.bakta_dir}/{wildcards.sample}.faa \
          --output_dir {output.eggnog_dir} \
          --cpu {resources.cpus} \
          -m diamond \
          --data_dir {params.dmnd_db} \
          --output {wildcards.sample} \
          --temp_dir {output.temp_dir} \
          --override > {log} 2>&1
        """


# ── Rule: secondary_metabolites_db — one-off antiSMASH reference download ───────
# Biology: antiSMASH detects biosynthetic gene clusters (BGCs) by matching known
# cluster/domain models, which live in a reference database it must fetch once.
# This rule is a pure download with no per-sample input.
#
# Takes in: nothing.
# Produces: 04.annotation/antismash/databases/ — the shared antiSMASH DB dir.
# Consumed by: secondary_metabolites_analysis (every sample waits on this DB).
#
# PATH NOTE: this DB directory (ANTISMASH_DB_DIR) lives at antismash/databases,
# right beside the antismash/{sample} directories the analysis rule writes. The
# project-wide `wildcard_constraints: sample=...` in 00_common pins {sample} to
# the real discovered sample names, so "databases" can never match the analysis
# rule's {sample} — this shared DB path can only ever be produced by THIS rule.
# (v1 relied on Snakemake's concrete-beats-wildcard tie-break instead; the
# constraint makes it explicit and collision-proof.)
#
# (v1 message: "--- antiSMASH: database download. ---")
rule secondary_metabolites_db:
    output:
        antismash_db = directory(ANTISMASH_DB_DIR),
    conda:
        "../../envs/antismash.yaml"
    log:
        LOGS + "/secondary_metabolites_database.log"
    priority: 4
    shell:
        """
        download-antismash-databases \
          --database-dir {output.antismash_db} > {log} 2>&1
        """


# ── Rule: secondary_metabolites_analysis — BGC detection per sample (antiSMASH) ─
# Biology: scans one genome for biosynthetic gene clusters (antibiotics,
# siderophores, and other secondary metabolites) and reports each cluster with
# its type and closest known cluster.
#
# Takes in: two DAG edges — the antiSMASH DB directory (from the rule above) and
#           the Bakta output DIRECTORY. In the shell it reads Bakta's annotated
#           GenBank file {sample}.gbff from inside the Bakta directory.
# Does: run antiSMASH with --genefinding-tool none, i.e. DO NOT re-predict genes
#       — reuse Bakta's existing gene calls carried in the .gbff. That is exactly
#       why this rule consumes .gbff (annotated GenBank) and not the bare .fna.
# Produces: 04.annotation/antismash/{sample}/ — the per-sample antiSMASH results.
# Consumed by: the report module (a terminal annotation product).
#
# (v1 message: "--- antiSMASH: secondary metabolite annotation. ---")
rule secondary_metabolites_analysis:
    input:
        antismash_db = ANTISMASH_DB_DIR,
        bakta_dir = DIR_ANNOTATION + "/bakta/{sample}",
    output:
        antismash_dir = directory(DIR_ANNOTATION + "/antismash/{sample}"),
    params:
        # Rule-local literals (not paths/config), so they stay inline here.
        taxon = 'bacteria',
        genefinding_tool = 'none',
    conda:
        "../../envs/antismash.yaml"
    resources:
        cpus = capped_cpus(24)
    log:
        LOGS + "/secondary_metabolites_{sample}.log"
    priority: 4
    shell:
        """
        antismash \
          --output-dir {output.antismash_dir} \
          --output-basename {wildcards.sample} \
          --databases {input.antismash_db} \
          --taxon {params.taxon} \
          --genefinding-tool {params.genefinding_tool} \
          --cpus {resources.cpus} \
          {input.bakta_dir}/{wildcards.sample}.gbff > {log} 2>&1
        """


# ── Rule: cazyme_db_download — fetch + checksum-verify the dbCAN database ────────
# Biology: dbCAN annotates carbohydrate-active enzymes (CAZymes) against HMM /
# DIAMOND / dbCAN-sub reference sets bundled in a versioned tarball. This rule
# downloads that tarball and, crucially, VERIFIES it before trusting it.
#
# Takes in: nothing (pure download).
# Does: download the .tar.gz and its .sha256; extract the expected hash from the
#       .sha256 (first field of line 1); compute the tarball's actual hash;
#       HARD-FAIL via `test "$expected" = "$actual"` if they differ, so a
#       truncated or tampered download can never be silently extracted; only on a
#       match, extract (flattening the top-level dir) and write the verified hash
#       into the sentinel file.
# Produces:
#   dbcan_db       = DBCAN_DB_DIR — the extracted database directory. Its version
#                    folder name is DERIVED FROM THE LINK in 00_common (DBCAN_DB_ID),
#                    so the folder can never disagree with the configured URL.
#   dbcan_verified = DBCAN_SENTINEL — a tiny marker holding the verified checksum;
#                    its existence is the "DB present AND integrity-checked" gate.
# Consumed by: cazyme_gene_cluster, which depends on the SENTINEL (not the whole
#              dir), so dbCAN never starts against an unverified database.
#
# conda: NONE — inherited from v1. This rule runs in the environment Snakemake was
# launched from and needs wget, tar, sha256sum and awk to be available there.
# Adding a small wget/coreutils env would be NEW behaviour, so it is left as-is.
#
# All four values used here (DBCAN_DB_DIR, DBCAN_SENTINEL, DBCAN_LINK,
# DBCAN_SHA_URL) come from 00_common — v1 hard-coded the folder name and derived
# the .sha256 URL inline with .replace(); neither is repeated here.
#
# Note on the doubled braces in the shell: Snakemake treats {…} as a placeholder,
# so a literal brace for awk must be written as {{…}} to survive to the shell.
rule cazyme_db_download:
    output:
        dbcan_db = directory(DBCAN_DB_DIR),
        dbcan_verified = DBCAN_SENTINEL,
    params:
        dbcan_db_url = DBCAN_LINK,
        sha_url = DBCAN_SHA_URL,
    log:
        LOGS + "/cazyme_db_download.log"
    priority: 4
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


# ── Rule: cazyme_gene_cluster — CAZyme + gene-cluster annotation (run_dbcan) ─────
# Biology: identifies carbohydrate-active enzymes, groups neighbouring ones into
# CAZyme Gene Clusters (CGCs), and predicts the substrate each cluster likely
# acts on — i.e. what sugars / polysaccharides this genome can process.
#
# Takes in: two DAG edges —
#   bakta_dir      = the Bakta output DIRECTORY; the shell reads BOTH {sample}.faa
#                    (proteins) and {sample}.gff3 (gene coordinates) from inside it.
#   dbcan_verified = DBCAN_SENTINEL, the integrity gate from cazyme_db_download;
#                    depending on the sentinel guarantees the database was
#                    checksum-verified before dbCAN runs.
# Does: four sequential run_dbcan sub-commands into ONE output dir, in order —
#         1. CAZyme_annotation    (HMM + DIAMOND + dbCAN-sub over the proteins)
#         2. gff_process          (map calls back onto the genome via the GFF)
#         3. cgc_finder           (cluster adjacent CAZymes into CGCs)
#         4. substrate_prediction (predict each CGC's substrate)
# Produces: 04.annotation/dbcan/{sample}/ — the per-sample CAZyme results.
# Consumed by: the report module (a terminal annotation product).
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
    resources:
        cpus = capped_cpus(24)
    log:
        LOGS + "/cazyme_{sample}.log"
    priority: 4
    shell:
        """
        # CAZyme annotation of protein sequences
        run_dbcan CAZyme_annotation \
          --input_raw_data {input.bakta_dir}/{wildcards.sample}.faa \
          --output_dir {output.dbcan_dir} \
          --db_dir {params.dbcan_db} \
          --mode protein \
          --threads {resources.cpus} \
          --methods hmm \
          --methods diamond \
          --methods dbCANsub > {log} 2>&1

        # CAZyme Gene Cluster (CGC) Annotation
        run_dbcan gff_process \
          --input_gff {input.bakta_dir}/{wildcards.sample}.gff3 \
          --output_dir {output.dbcan_dir} \
          --db_dir {params.dbcan_db} \
          --gff_type prodigal \
          --threads {resources.cpus} >> {log} 2>&1

        # CAZyme Gene Cluster (CGC) Identification
        run_dbcan cgc_finder \
          --output_dir {output.dbcan_dir} >> {log} 2>&1

        # CGC Substrate Prediction
        run_dbcan substrate_prediction \
          --output_dir {output.dbcan_dir} \
          --db_dir {params.dbcan_db} >> {log} 2>&1
        """
