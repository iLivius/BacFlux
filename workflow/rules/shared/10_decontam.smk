# BacFlux v2.0.0 — contamination screening of the draft assembly: reads mapped
# back for coverage, megablast against NCBI nt for taxonomy, BlobTools to join
# the two, and a Python selector that applies the user's keep/drop policy and
# records every decision.
#
# An isolate assembly is not guaranteed to be one organism. Culture contaminants,
# index hopping, carry-over from a neighbouring library and adapter/vector debris
# all arrive as extra contigs. Asking two questions per contig — what does it look
# like taxonomically, and how deeply is it covered by this sample's own reads —
# separates the isolate from everything else. Doing it BEFORE annotation matters:
# contaminant contigs inflate CheckM's contamination estimate, pull GTDB-Tk off
# the right lineage, and pollute every downstream annotation and AMR call.
#
# Chain: DRAFT_CONTIGS → map_contigs (coverage) + blast_contigs (taxonomy) →
# blob_json → blob_table → select_contigs → DECONTAM_CONTIGS + the audit files.
#
# index_contigs      : bowtie2-build over the draft. Short-read modes only, and
#                      only because map_contigs needs an index to map against.
# map_contigs        : this sample's own reads back onto its own draft, giving the
#                      per-contig depth track. Three alternative rule bodies, one
#                      per read type — see the section prose below.
# blast_contigs      : megablast of every draft contig against NCBI nt.
# blast_final_contigs: the same screen re-run on the DELIVERED genome. Long-read
#                      modes only, and it exists for the plasmid check in
#                      shared/60_plasmid.smk, not for decontamination.
# blob_json          : BlobTools joins coverage and taxonomy into one database.
# blob_table         : collapses each contig's many BLAST hits into one call per
#                      taxonomic rank.
# select_contigs     : applies the decontamination policy, writes the kept contigs
#                      and the per-contig audit trail.
#
# Where decontamination sits per mode (decision D3 — v1 order preserved exactly):
#
#   mode      DRAFT_CONTIGS (screened)        DECONTAM_CONTIGS (written)
#   ────────  ──────────────────────────────  ──────────────────────────────────
#   illumina  SPAdes contigs_filt.fasta       = FINAL_CONTIGS (decontam is last)
#   contigs   filtered input contigs_filt     = FINAL_CONTIGS (decontam is last)
#   nanopore  reoriented {sample}_fixed       assembly_decontam.fasta → Medaka
#   hybrid    the ILLUMINA SPAdes draft       contigs_sel.fasta → Snippy reference
#                                             + the QC comparator genome
#
# That is why select_contigs writes DECONTAM_CONTIGS and not FINAL_CONTIGS: in
# nanopore and hybrid the delivered genome is produced later by the mode's front
# end (rules/nanopore/, rules/hybrid/), and hard-coding FINAL_CONTIGS here would
# make nanopore circular (select → final → Medaka → select).
#
# Everything referenced here is defined once in 00_common.smk and never
# re-derived: DRAFT_CONTIGS, DECONTAM_CONTIGS, FINAL_CONTIGS, DECONTAM_DIR,
# DECONTAM_BAM, BLOB_PREFIX/JSON/COV, BLOB_TABLE(_PREFIX), CONTIG_LIST,
# CONTIG_DECISIONS, COMPOSITION, BLASTOUT, PLASMID_BLASTOUT, TRIM_R1/TRIM_R2,
# FILT_LONG, BLASTDB, NT_VERSION, SELECT_TAXONOMY_SCRIPT, DECONTAMINATION, LOGS,
# CPUS, capped_cpus, and the capability flags.
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/x.yaml" climbs shared/ → rules/ → workflow/ → workflow/envs/x.yaml.
#
# Resource convention: cpu-bound rules declare Snakemake's built-in
# `threads: capped_cpus(N)` and refer to `{threads}` in the shell. Using the
# BUILT-IN keyword (rather than a custom `resources: cpus`) is what makes
# `--cores N` actually enforce the limit, so a plain `snakemake --cores N` is
# safe on its own and no extra `--resources` flag is needed.


# ────────────── Coverage track (Bowtie2 / minimap2) ────────────
# map_contigs produces the SAME two files in every mode (DECONTAM_BAM plus its
# .bai) but gets there with a different aligner, from a different read type, in a
# different conda environment. The choice is made here at parse time — while
# Snakemake reads the workflow, before any job runs — with a plain if/elif/else
# that defines exactly ONE rule body. MODE is fixed for the whole run, so a reader
# working in nanopore mode sees one 15-line map_contigs and nothing else.
#
# Why not one rule with an input function? The branches differ in their `conda:`
# environment (bowtie.yaml vs minimap.yaml), and Snakemake resolves a rule's conda
# env when it deploys environments — an input function cannot reach it. The other
# alternative, a single rule with a union env carrying bowtie2 AND minimap2 plus a
# branchy shell, would make every mode build a bigger environment and turn a short
# shell into a conditional block. This matches the house pattern already used for
# `if PHAGE_CALLER == "genomad":` in shared/60_plasmid.smk and shared/70_phage.smk.

if HAS_SHORT_READS:

    # Bowtie2 needs its own index of the reference before it can map anything.
    # Module-local constant: both the producer (this rule) and the only consumer
    # (map_contigs, immediately below) live in this file, so the index prefix is
    # not a cross-module contract and does not belong in 00_common. It is built
    # off DECONTAM_DIR, so no stage number is re-derived.
    _BT2_PREFIX = DECONTAM_DIR + "/{sample}_contigs"

    # ── index_contigs — the Bowtie2 index of the draft assembly ──
    # Takes in:  DRAFT_CONTIGS — the length/coverage-filtered SPAdes contigs
    #            written by rule filter_contigs in illumina/20_assembly.smk or
    #            hybrid/20_assembly.smk.
    # Does:      bowtie2-build, which writes six binary index files.
    # Produces:  six temp() .bt2 files; they exist only to let the next rule map.
    # Consumed by: map_contigs.
    #
    # (v1 message: "--- Bowtie2: Build contig db. ---")
    rule index_contigs:
        input:
            contigs = DRAFT_CONTIGS,
        output:
            # All six are temp(): the index is regenerated far more cheaply than
            # it is worth storing, and it is useless once the BAM exists.
            idx1  = temp(_BT2_PREFIX + ".1.bt2"),
            idx2  = temp(_BT2_PREFIX + ".2.bt2"),
            idx3  = temp(_BT2_PREFIX + ".3.bt2"),
            idx4  = temp(_BT2_PREFIX + ".4.bt2"),
            ridx1 = temp(_BT2_PREFIX + ".rev.1.bt2"),
            ridx2 = temp(_BT2_PREFIX + ".rev.2.bt2"),
        params:
            # bowtie2-build takes the shared PREFIX of those six files, not a
            # filename. v1 wrapped this in a pointless temp() — temp() in params
            # does nothing at all — dropped here.
            basename = _BT2_PREFIX,
        conda:
            "../../envs/bowtie.yaml"
        log:
            LOGS + "/index_contigs_{sample}.log"
        shell:
            """
            bowtie2-build \
              -f {input.contigs} \
              {params.basename} > {log} 2>&1
            """

    # ── map_contigs (short-read modes) — reads back onto the draft ──
    # Mapping the sample's own trimmed reads onto its own assembly gives per-contig
    # read depth. A contig from a minor contaminant is usually covered at a very
    # different depth from the isolate's chromosome, and that depth is the second
    # axis — alongside taxonomy — that BlobTools separates organisms on.
    #
    # Takes in: the six Bowtie2 index files (the DAG edge back to index_contigs)
    #           and the fastp-trimmed pairs TRIM_R1/TRIM_R2 written by the mode's
    #           read front end (illumina/10_reads.smk, hybrid/10_reads.smk).
    # Does:     bowtie2 → SAM, converted to BAM, coordinate-sorted, then indexed.
    # Produces: DECONTAM_BAM + .bai, both temp().
    # Consumed by: blob_json (the coverage leg) and map_evaluation (Qualimap, in
    #              shared/20_qc.smk). temp() keeps the BAM alive until both are done.
    #
    # v1→v2: `--write-index` is dropped from samtools sort. It wrote a .csi index
    # alongside the .bai that `samtools index -b` writes on the next line, and
    # BlobTools only reads the .bai. One fewer temp output, no functional change;
    # this unifies on the BacFluxL+ form.
    #
    # (v1 message: "--- Bowtie2: Map reads against contigs. ---")
    rule map_contigs:
        input:
            idx1  = _BT2_PREFIX + ".1.bt2",
            idx2  = _BT2_PREFIX + ".2.bt2",
            idx3  = _BT2_PREFIX + ".3.bt2",
            idx4  = _BT2_PREFIX + ".4.bt2",
            ridx1 = _BT2_PREFIX + ".rev.1.bt2",
            ridx2 = _BT2_PREFIX + ".rev.2.bt2",
            r1 = TRIM_R1,
            r2 = TRIM_R2,
        output:
            bam = temp(DECONTAM_BAM),
            bai = temp(DECONTAM_BAM + ".bai"),
        params:
            db = _BT2_PREFIX,
        conda:
            "../../envs/bowtie.yaml"
            # Uncapped, as in v1: bowtie2 and samtools sort both scale well, and
            # this is one of the few rules where the full budget is worth giving.
        threads: CPUS
        log:
            LOGS + "/map_contigs_{sample}.log"
        shell:
            """
            bowtie2 \
              -x {params.db} \
              -1 {input.r1} \
              -2 {input.r2} \
              -p {threads} \
              -t 2> {log} | \
            samtools view \
              -@ {threads} \
              -hbS - | \
            samtools sort \
              -@ {threads} \
              -o {output.bam} - >> {log} 2>&1

            samtools index \
              -@ {threads} \
              -b {output.bam} >> {log} 2>&1
            """

elif HAS_LONG_READS:

    # ── map_contigs (nanopore) — ONT reads back onto the draft ──
    # Same purpose as the short-read version: a per-contig depth track for
    # BlobTools. minimap2's map-ont preset handles the higher error rate of raw
    # ONT reads, which bowtie2 cannot. No Bowtie2 index rule here — minimap2
    # indexes the reference on the fly.
    #
    # Takes in: FILT_LONG (filtlong-filtered ONT reads, from nanopore/10_reads.smk)
    #           and DRAFT_CONTIGS (the dnaapler-reoriented Flye assembly, from
    #           nanopore/20_assembly.smk).
    # Produces: DECONTAM_BAM + .bai, both temp().
    # Consumed by: blob_json and map_evaluation (Qualimap, shared/20_qc.smk).
    #
    # (v1 message: "--- Minimap2: Map reads against contigs. ---";
    #  v1 rule was also called map_contigs, its Qualimap rule map_qc → D5 renames
    #  that one to map_evaluation everywhere.)
    rule map_contigs:
        input:
            reads = FILT_LONG,
            contigs = DRAFT_CONTIGS,
        output:
            bam = temp(DECONTAM_BAM),
            bai = temp(DECONTAM_BAM + ".bai"),
        conda:
            "../../envs/minimap.yaml"
        threads: CPUS
        log:
            LOGS + "/map_contigs_{sample}.log"
        shell:
            """
            minimap2 \
              -ax map-ont {input.contigs} \
              {input.reads} 2> {log} | \
            samtools view \
              -S \
              -b \
              -u \
              -@ {threads} | \
            samtools sort \
              -o {output.bam} \
              -@ {threads} 2>> {log}

            samtools index \
              {output.bam} \
              -@ {threads} 2>> {log}
            """

else:

    # ── map_contigs (contigs mode) — a deliberate FAKE coverage track ──
    # There are no reads in this mode: the user hands over finished assemblies.
    # `blobtools create -b` still wants a BAM, so the contigs are mapped against
    # THEMSELVES. The resulting depth is near-uniform and carries no information.
    #
    # Preserved from v1 on purpose, and the consequence has to be understood before
    # anyone tries to "use" this BAM: coverage-based separation in BlobTools is
    # MEANINGLESS in contigs mode. Only the taxonomy leg (blastn → BlobTools →
    # selector) is doing real work. That is also why this mode has no Qualimap
    # rule — map_evaluation in shared/20_qc.smk is gated on HAS_READS, and a
    # mapping-quality report on a self-alignment would be a chart of nothing.
    #
    # Takes in: DRAFT_CONTIGS — the input contigs after rule filter_contigs in
    #           contigs/20_assembly.smk: short and low-coverage records are dropped
    #           only if the headers are SPAdes-style (length and coverage are read
    #           out of them); any other style just gets its headers trimmed.
    # Produces: DECONTAM_BAM + .bai, both temp().
    # Consumed by: blob_json only.
    #
    # (v1 message: "--- Minimap2: Map reads against contigs. ---")
    rule map_contigs:
        input:
            contigs = DRAFT_CONTIGS,
        output:
            bam = temp(DECONTAM_BAM),
            bai = temp(DECONTAM_BAM + ".bai"),
        conda:
            "../../envs/minimap.yaml"
        threads: CPUS
        log:
            LOGS + "/map_contigs_{sample}.log"
        shell:
            # No -x preset: contigs against themselves, exactly as v1 FastaFlux.
            """
            minimap2 \
              -a {input.contigs} \
              {input.contigs} 2> {log} | \
            samtools view \
              -S \
              -b \
              -u \
              -@ {threads} | \
            samtools sort \
              -o {output.bam} \
              -@ {threads} 2>> {log}

            samtools index \
              {output.bam} \
              -@ {threads} 2>> {log}
            """


# ─────────────── Taxonomy screen (megablast vs nt) ─────────────
# blast_contigs — megablast every draft contig against the NCBI nucleotide
# database, keeping the top hits WITH their taxids and subject titles. This is the
# taxonomy leg BlobTools turns into a per-contig genus call. megablast rather than
# plain blastn because near-identical matches to known genomes are what an isolate
# assembly is expected to produce, and it is far faster.
#
# Takes in: DRAFT_CONTIGS — the mode's draft assembly, i.e. the same contigs the
#           selector will filter. The screen is structurally pinned to the draft:
#           BlobTools must see every contig it is being asked to judge, including
#           the ones about to be thrown away.
# Does:     one blastn -task megablast per sample against {blast_db}/{nt_version}
#           (nt_version comes from parameters.nt_version, defaulting to core_nt).
# Produces: BLASTOUT — 15 tab-separated columns, subject title (stitle) LAST.
# Consumed by: blob_json in every mode, and — in ILLUMINA AND CONTIGS MODE ONLY —
#              the supplementary "does the nt hit say plasmid?" check in
#              shared/60_plasmid.smk, which greps that last column. The two
#              long-read modes grep blast_final_contigs' output instead; see
#              NEEDS_FINAL_BLAST in 00_common.smk for why.
#
# Do NOT change the -outfmt string. All four v1 modes agree on it, and the plasmid
# check greps the subject title for the word "plasmid". Drop or move stitle and
# that grep silently matches nothing: every Platon call comes back "not verified
# by BLAST search", with no error and no warning.
#
# (v1 message: "--- BLAST: Contigs against NCBI nt db. ---")
rule blast_contigs:
    input:
        contigs = DRAFT_CONTIGS,
    output:
        blast = BLASTOUT,
    params:
        # BLASTDB env var points at the directory (so the taxonomy .dmp files and
        # the taxid mapping are found); -db points at the versioned subfolder.
        dir = BLASTDB,
        db = os.path.join(BLASTDB, NT_VERSION),
    conda:
        "../../envs/blast.yaml"
    threads: capped_cpus(24)
    log:
        LOGS + "/blast_contigs_{sample}.log"
    shell:
        """
        BLASTDB={params.dir} \
        blastn \
          -task megablast \
          -query {input.contigs} \
          -db {params.db} \
          -outfmt '6 qseqid staxids bitscore pident evalue length qlen slen qcovs qcovhsp sskingdoms scomnames sscinames sblastnames stitle' \
          -num_threads {threads} \
          -evalue 1e-5 \
          -max_target_seqs 50 \
          -max_hsps 5 \
          -out {output.blast} > {log} 2>&1
        """


# ────────── Second nt screen of the delivered genome ───────────
# Defined only when NEEDS_FINAL_BLAST is true, which 00_common.smk sets to
# HAS_LONG_READS — so nanopore and hybrid get this rule and illumina and contigs
# never see it (there PLASMID_BLASTOUT simply IS BLASTOUT and no second rule
# enters the DAG).
#
# Both long-read modes need it for the same reason: plasmid_search looks each
# Platon-called contig up in a BLAST table BY CONTIG ID, so the table has to have
# been computed over the contigs Platon actually reported on. In hybrid the screen
# above runs on the ILLUMINA draft while Platon runs on the delivered ONT genome —
# SPAdes names its contigs NODE_1_length_… and Flye names them contig_1, so
# plasmid_search's `grep -m 1 "$contig" <blastout>` would never match a single ID
# and every plasmid would be reported as "not verified by BLAST search". In
# nanopore the screen runs on the pre-Medaka assembly while Platon runs on the
# post-Medaka consensus, and nothing guarantees Medaka preserves contig headers.
#
# v1 BacFluxL+ solved this by running its own blastn inside plasmid_search; v2
# keeps the second BLAST but defines it here, next to the identical command it
# duplicates, and routes it through the PLASMID_BLASTOUT constant.
#
# Cost: the long-read modes therefore BLAST against nt twice per sample (draft +
# final). For hybrid that is exactly what v1 did — not a new cost. For nanopore it
# is one extra blastn, buying immunity from an unverified assumption about Medaka.
if NEEDS_FINAL_BLAST:

    # ── blast_final_contigs — nt screen of the DELIVERED genome ──
    # Takes in: FINAL_CONTIGS — the genome the long-read front end actually
    #           delivers (nanopore: the Medaka consensus, via finalize_contigs in
    #           nanopore/30_polish.smk; hybrid: the ONT+Polypolish genome) — NOT
    #           the draft that blast_contigs screened.
    # Produces: PLASMID_BLASTOUT, same 15-column outfmt as blast_contigs.
    # Consumed by: plasmid_search in shared/60_plasmid.smk, and nothing else —
    #              BlobTools always uses the draft table.
    rule blast_final_contigs:
        input:
            contigs = FINAL_CONTIGS,
        output:
            blast = PLASMID_BLASTOUT,
        params:
            dir = BLASTDB,
            db = os.path.join(BLASTDB, NT_VERSION),
        conda:
            "../../envs/blast.yaml"
        threads: capped_cpus(24)
        log:
            LOGS + "/blast_final_contigs_{sample}.log"
        shell:
            """
            BLASTDB={params.dir} \
            blastn \
              -task megablast \
              -query {input.contigs} \
              -db {params.db} \
              -outfmt '6 qseqid staxids bitscore pident evalue length qlen slen qcovs qcovhsp sskingdoms scomnames sscinames sblastnames stitle' \
              -num_threads {threads} \
              -evalue 1e-5 \
              -max_target_seqs 50 \
              -max_hsps 5 \
              -out {output.blast} > {log} 2>&1
            """


# ────────────── Blob database and per-contig table ─────────────
# blob_json — BlobTools takes the assembly, the read-depth track and the BLAST
# hits and builds the "blobplot" database: per contig, its length, GC, coverage
# and a taxonomic assignment resolved through the NCBI taxonomy dump.
#
# Takes in:
#   contigs     = DRAFT_CONTIGS   (the same contigs both other legs used)
#   bam / bai   = DECONTAM_BAM    (coverage; a self-map in contigs mode)
#   blast       = BLASTOUT        (taxonomy)
#   nodes/names = the NCBI taxonomy dump shipped inside the BLAST database dir.
#                 They are declared as INPUTS, not params, so a missing/incomplete
#                 taxonomy dump is reported up front instead of half-way through.
# Produces: BLOB_JSON and BLOB_COV, both temp() — intermediate binaries that only
#           blob_table reads.
# Consumed by: blob_table.
#
# BLOB_COV's filename is chosen by BlobTools from the BAM's basename; that is why
# 00_common derives it from DECONTAM_BAM instead of typing it out again.
#
# (v1 message: "--- BlobTools: Screen BLAST hits for contaminants. ---")
rule blob_json:
    input:
        contigs = DRAFT_CONTIGS,
        bam = DECONTAM_BAM,
        bai = DECONTAM_BAM + ".bai",
        blast = BLASTOUT,
        nodes = os.path.join(BLASTDB, "nodes.dmp"),
        names = os.path.join(BLASTDB, "names.dmp"),
    output:
        json = temp(BLOB_JSON),
        cov = temp(BLOB_COV),
    params:
        # blobtools create takes an output PREFIX and appends .blobDB.json itself.
        basename = BLOB_PREFIX,
    conda:
        "../../envs/blobtools.yaml"
    log:
        # v1 bug fixed here: BacFlux and BacFluxL+ pointed BOTH this rule and
        # blob_table at logs/blob_table_{sample}.log, and this one truncated it
        # with `>` — so whichever ran second destroyed the other's log. v2 gives
        # each rule its own file.
        LOGS + "/blob_json_{sample}.log"
    shell:
        """
        blobtools create \
          -i {input.contigs} \
          -b {input.bam} \
          -t {input.blast} \
          --nodes {input.nodes} \
          --names {input.names} \
          -o {params.basename} > {log} 2>&1
        """


# blob_table — a contig usually has many BLAST hits pointing at several taxa. The
# "bestsum" tax rule sums bitscores per taxon and keeps the winner, at every rank.
# The result is the flat table the selector reads: one row per contig with its
# length, coverage, GC and its assigned taxonomy from superkingdom down to species.
#
# Takes in: BLOB_JSON from blob_json.
# Does:     blobtools view --taxrule bestsum --rank all --hits.
# Produces: BLOB_TABLE (kept, not temp — it is the evidence behind every
#           keep/discard decision and is worth being able to re-read).
# Consumed by: select_contigs.
#
# Completely mode-independent: same input shape, same command, same output, in all
# four modes.
#
# (v1 message: "--- BlobTools: Collapse taxonomic assignment of BLAST hits
#  according to sum of best scores. ---")
rule blob_table:
    input:
        json = BLOB_JSON,
    output:
        bestscore = BLOB_TABLE,
    params:
        # blobtools view takes an output PREFIX too, and appends
        # ".blob.blobDB.table.txt"; 00_common derives BLOB_TABLE from this prefix.
        basename = BLOB_TABLE_PREFIX,
    conda:
        "../../envs/blobtools.yaml"
    log:
        LOGS + "/blob_table_{sample}.log"
    shell:
        """
        blobtools view \
          --input {input.json} \
          --out {params.basename} \
          --taxrule bestsum \
          --rank all \
          --hits > {log} 2>&1
        """


# ──────────────── Contig selection and audit trail ─────────────
# select_contigs — where contigs are actually kept or dropped. The helper script
# reads the BlobTools table, resolves each contig to a genus, and applies the
# policy the user configured in parameters.decontamination.mode:
#   auto     — keep the genus carried by the most CONTIGS. A count, not base
#              pairs: many short contigs can outvote the genus that holds the
#              genome, and choose_auto_genus() records a case where that sent a
#              4.5 Mb chromosome out as contamination.
#   include  — keep only the listed genera
#   exclude  — drop only the listed genera
#   off      — keep everything (still writes the audit files)
# plus discard_no_hit, which decides what happens to contigs nt could not place.
#
# Genus matching is deliberately fuzzy, and it is worth knowing before reading an
# audit file: the selector treats a curated set of split genera as one target
# (Paenibacillus/Peribacillus/Priestia → Bacillus, Pseudarthrobacter/
# Paenarthrobacter → Arthrobacter, Paraburkholderia → Burkholderia) and also
# strips the prefixes brady/meso/neo/sino/aeri/caldi/geo. BLAST and BlobTools
# routinely scatter one organism across those related names, and without the
# aliases a genuine isolate contig gets dropped as a contaminant. It is a
# safeguard against false removal, not taxonomic reconciliation — see
# GENUS_EQUIVALENCE_ALIASES and GENUS_EQUIVALENCE_PREFIXES in
# scripts/10_decontam/select_contigs_by_taxonomy.py.
#
# In HYBRID mode this rule can cost more than a taxonomy row. Only the Illumina
# reads mapping to the SELECTED contigs become SEL_R1/SEL_R2, and filtlong scores
# ONT reads against those, so a contig dropped here can take its ONT reads with it
# and vanish from the assembly entirely. Plasmids are the usual casualty, because
# their best nt hit is often a different genus from the host. The full worked case
# and the fixes are in the decontamination block of config/config.yaml.
#
# Takes in:
#   bestscore = BLOB_TABLE    (the per-contig taxonomy + coverage table)
#   contigs   = DRAFT_CONTIGS (the sequences themselves)
# Produces (all four paths are the same in every mode except the last):
#   abund     = COMPOSITION       one line per genus, with its share of the DNA
#                                 and its share of the contig count —
#                                 "Bacillus: bases 0.29  contigs 0.30" — over
#                                 EVERY contig in the BlobTools table, kept or
#                                 dropped → read by rule annotation in
#                                 shared/40_annotation.smk as Bakta's --genus
#                                 hint. Read that rule's shell preface before
#                                 trusting the hint: the two labelled figures
#                                 broke its numeric sort.
#   list      = CONTIG_LIST       the kept contig IDs, one per line
#   decisions = CONTIG_DECISIONS  the audit TSV: every contig, its genus, and the
#                                 REASON it was kept or dropped. Required by the
#                                 project convention that every filtering decision
#                                 is auditable (CLAUDE.md).
#   contigs   = DECONTAM_CONTIGS  the kept sequences — the per-mode table in this
#                                 file's header says who consumes them next
#
# conda: NONE — inherited from all four v1 modes. The selector is stdlib-only
# Python and runs in the environment Snakemake was launched from. Adding an env
# would be new behaviour; it is the same deferred decision as cazyme_db_download
# (see docs/README_notes.md item 3).
#
# The seven policy values — mode, discard_no_hit and the five genus/override
# fields — come verbatim from the DECONTAMINATION dict that _decontam_settings()
# resolves once in 00_common.smk. Two of them are the per-sample escape hatch:
# include_genera_by_sample and sample_overrides are FILE PATHS keyed by sample
# name, which is how one awkward isolate deviates without changing the batch
# policy. All 14 flags are passed with {...:q} quoting, so a genus list containing
# spaces or a path with odd characters survives the shell intact.
#
# v1→v2 (additive): v1 had no log:, so the selector's warnings — notably
# "WARNING: N contigs from the BlobTools table were not found in the FASTA" and
# the per-sample kept/dropped counts — went to the console and were lost on a
# large batch. They now go to the log file instead. Trade-off worth knowing: you
# have to open the log to see them.
#
# v1→v2 (removed): FastaFlux used to re-linearise the FASTA into
# contigs_filt_lin.fasta before calling the selector. That was redundant — the
# selector's read_fasta() accumulates sequence chunks and joins them, so it
# already handles wrapped FASTA, and it always writes 80-column-wrapped output
# regardless. The selector's output is byte-identical without the pre-step, so the
# intermediate file simply disappears from the contigs-mode output tree.
rule select_contigs:
    input:
        bestscore = BLOB_TABLE,
        contigs = DRAFT_CONTIGS,
    output:
        abund = COMPOSITION,
        list = CONTIG_LIST,
        decisions = CONTIG_DECISIONS,
        contigs = DECONTAM_CONTIGS,
    params:
        selector = SELECT_TAXONOMY_SCRIPT,
        mode = DECONTAMINATION["mode"],
        include_genera = DECONTAMINATION["include_genera"],
        include_genera_by_sample = DECONTAMINATION["include_genera_by_sample"],
        exclude_genera = DECONTAMINATION["exclude_genera"],
        exclude_genera_file = DECONTAMINATION["exclude_genera_file"],
        sample_overrides = DECONTAMINATION["sample_overrides"],
        discard_no_hit = DECONTAMINATION["discard_no_hit"],
    log:
        LOGS + "/select_contigs_{sample}.log"
    shell:
        """
        python {params.selector:q} \
          --bestscore {input.bestscore:q} \
          --contigs {input.contigs:q} \
          --sample {wildcards.sample:q} \
          --mode {params.mode:q} \
          --include-genera {params.include_genera:q} \
          --include-genera-by-sample {params.include_genera_by_sample:q} \
          --exclude-genera {params.exclude_genera:q} \
          --exclude-genera-file {params.exclude_genera_file:q} \
          --sample-overrides {params.sample_overrides:q} \
          --discard-no-hit {params.discard_no_hit:q} \
          --output-list {output.list:q} \
          --output-fasta {output.contigs:q} \
          --composition {output.abund:q} \
          --decisions {output.decisions:q} > {log} 2>&1
        """


# ══════════════ Contaminant screen of the DELIVERED long-read genome ══════════════
#
# Hybrid screens the Illumina SPAdes draft (the HAS_SHORT_READS branch above). The
# genome it delivers comes from the other side of the mode — Flye, Medaka, dnaapler,
# Polypolish — and passes through none of that.
#
# That was defensible for as long as filtlong scored every long read against the
# decontaminated short reads: a contaminant read was removed before Flye ever saw it,
# and the rule that did it said so. With parameters.hybrid.short_read_guidance now
# defaulting to false, nothing removes contaminant long reads, so the assembly they
# build is screened here instead.
#
# The costly half is already paid for. blast_final_contigs megablasts FINAL_CONTIGS
# against nt so that plasmid_search can read it, and PLASMID_BLASTOUT carries the same
# 15 columns BlobTools wants. Only a coverage track and the BlobTools join are new.
if LONGREAD_SCREEN:

    # ── map_final_contigs — ONT reads onto the delivered genome ──
    # Takes in: FILT_LONG (the filtered ONT reads) and FINAL_CONTIGS.
    # Does:     minimap2 -ax map-ont, exactly as the nanopore screen maps its draft.
    #           The ONT reads are used rather than the Illumina ones on purpose: the
    #           whole reason this screen exists is that the short reads may no longer
    #           be involved in the ONT leg, and a coverage track derived from them
    #           would put the short-read bias back into the decision.
    # Produces: FINAL_BAM (+ .bai), both temp — only BlobTools reads them.
    rule map_final_contigs:
        input:
            reads = FILT_LONG,
            contigs = FINAL_CONTIGS,
        output:
            bam = temp(FINAL_BAM),
            bai = temp(FINAL_BAM + ".bai"),
        conda:
            "../../envs/minimap.yaml"
        threads: CPUS
        log:
            LOGS + "/map_final_contigs_{sample}.log"
        shell:
            """
            minimap2 \
              -ax map-ont {input.contigs} \
              {input.reads} 2> {log} | \
            samtools view -S -b -u -@ {threads} | \
            samtools sort -o {output.bam} -@ {threads} 2>> {log}

            samtools index {output.bam} -@ {threads} 2>> {log}
            """

    # ── blob_json_final / blob_table_final — the same join, on the delivered genome ──
    # Takes in: FINAL_CONTIGS, the coverage track above, and PLASMID_BLASTOUT, which
    #           blast_final_contigs already produced for plasmid_search.
    # Produces: FINAL_BLOB_TABLE, one row per delivered contig with its coverage and
    #           the genus BlobTools settled on.
    rule blob_json_final:
        input:
            contigs = FINAL_CONTIGS,
            bam = FINAL_BAM,
            bai = FINAL_BAM + ".bai",
            blast = PLASMID_BLASTOUT,
            nodes = os.path.join(BLASTDB, "nodes.dmp"),
            names = os.path.join(BLASTDB, "names.dmp"),
        output:
            json = temp(FINAL_BLOB_JSON),
            cov = temp(FINAL_BLOB_COV),
        params:
            basename = FINAL_BLOB_PREFIX,
        conda:
            "../../envs/blobtools.yaml"
        log:
            LOGS + "/blob_json_final_{sample}.log"
        shell:
            """
            blobtools create \
              -i {input.contigs} \
              -b {input.bam} \
              -t {input.blast} \
              --nodes {input.nodes} \
              --names {input.names} \
              -o {params.basename} > {log} 2>&1
            """

    rule blob_table_final:
        input:
            json = FINAL_BLOB_JSON,
        output:
            table = FINAL_BLOB_TABLE,
        params:
            basename = FINAL_BLOB_TABLE_PREFIX,
        conda:
            "../../envs/blobtools.yaml"
        log:
            LOGS + "/blob_table_final_{sample}.log"
        shell:
            """
            blobtools view \
              --input {input.json} \
              --out {params.basename} \
              --taxrule bestsum \
              --rank all \
              --hits > {log} 2>&1
            """

    # ── screen_final_contigs — verdict per delivered contig, discarding nothing ──
    # Takes in: FINAL_BLOB_TABLE and FINAL_CONTIGS.
    # Does:     the same selector the Illumina screen uses, with the same genus
    #           configuration — the organism does not change because the assembler
    #           did — but its own mode, from parameters.decontamination.long_read_mode,
    #           defaulting to "off".
    #
    #           "off" makes the selector return keep/mode_off for every contig, so the
    #           audit is written and nothing is removed. That default is deliberate:
    #           BlobTools identifies outliers within a cloud of contigs, and a
    #           long-read assembly has two to five. A false positive there deletes a
    #           whole replicon — usually the plasmid — rather than trimming a
    #           fragment, and this project has twice seen BLAST bestsum follow
    #           database composition instead of biology on exactly that call.
    # Produces: FINAL_TAXO_DECISIONS, one audited row per contig. The FASTA it writes
    #           is a byproduct; FINAL_CONTIGS remains what every later stage reads.
    rule screen_final_contigs:
        input:
            bestscore = FINAL_BLOB_TABLE,
            contigs = FINAL_CONTIGS,
        output:
            decisions = FINAL_TAXO_DECISIONS,
            abund = DECONTAM_DIR + "/{sample}_final_composition.tsv",
            list = temp(DECONTAM_DIR + "/{sample}_final_keep.list"),
            contigs = temp(DECONTAM_DIR + "/{sample}_final_screened.fasta"),
        params:
            selector = SELECT_TAXONOMY_SCRIPT,
            mode = LONGREAD_SCREEN_MODE,
            include_genera = DECONTAMINATION["include_genera"],
            include_genera_by_sample = DECONTAMINATION["include_genera_by_sample"],
            exclude_genera = DECONTAMINATION["exclude_genera"],
            exclude_genera_file = DECONTAMINATION["exclude_genera_file"],
            sample_overrides = DECONTAMINATION["sample_overrides"],
            discard_no_hit = DECONTAMINATION["discard_no_hit"],
        log:
            LOGS + "/screen_final_contigs_{sample}.log"
        shell:
            """
            python {params.selector:q} \
              --bestscore {input.bestscore:q} \
              --contigs {input.contigs:q} \
              --sample {wildcards.sample:q} \
              --mode {params.mode:q} \
              --include-genera {params.include_genera:q} \
              --include-genera-by-sample {params.include_genera_by_sample:q} \
              --exclude-genera {params.exclude_genera:q} \
              --exclude-genera-file {params.exclude_genera_file:q} \
              --sample-overrides {params.sample_overrides:q} \
              --discard-no-hit {params.discard_no_hit:q} \
              --output-list {output.list:q} \
              --output-fasta {output.contigs:q} \
              --composition {output.abund:q} \
              --decisions {output.decisions:q} > {log} 2>&1
            """
