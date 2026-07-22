# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — Stage 10 decontamination module (rules/shared/10_decontam.smk)
#
# The biology this module does: a bacterial isolate assembly is not guaranteed to
# be one organism. Culture contaminants, index hopping, carry-over from a
# neighbouring library and adapter/vector debris all end up as extra contigs. We
# therefore ask, per contig, "what does this look like taxonomically, and how
# well is it covered by this sample's own reads?" and drop the contigs that do
# not belong to the isolate. Doing this BEFORE annotation matters: contaminant
# contigs inflate CheckM's contamination estimate, pull GTDB-Tk off the right
# lineage, and pollute every downstream annotation and AMR call.
#
# How that is done (BlobTools' classic recipe):
#   1. map this sample's reads back onto its own draft assembly  -> coverage
#   2. megablast every contig against NCBI nt                    -> taxonomy
#   3. BlobTools joins the two into one per-contig table
#   4. a small Python selector applies the user's decontamination policy and
#      writes the kept contigs plus an audit trail of every decision
#
# Data flow (top to bottom):
#
#   DRAFT_CONTIGS ──┬─► index_contigs ─► map_contigs ─► {sample}_map.bam ──┐
#   (from Stage 4)  │   (short reads only)                    │            │
#                   │                                          ▼           │
#                   │                                  map_evaluation      │
#                   │                                  (Qualimap, 20_qc)   │
#                   ├─► blast_contigs ─► {sample}_blastout ────────────────┤
#                   │        │                                             ▼
#                   │        └────────► plasmid_search (60_plasmid) ─► blob_json
#                   │            (non-hybrid modes)                        │
#                   │                                                      ▼
#                   │                                                 blob_table
#                   │                                                      │
#                   └──────────────────────────────────────────────► select_contigs
#                                                                          │
#              ┌───────────────────────────────────────────────────────────┤
#              ▼                    ▼                   ▼                  ▼
#     {sample}_composition   contigs.list   contig_taxonomy_    DECONTAM_CONTIGS
#      -> Bakta --genus                     decisions.tsv        -> see below
#         (40_annotation)                    (audit trail)
#
#   HYBRID ONLY: blast_final_contigs ─► {sample}_final_blastout
#                (FINAL_CONTIGS)          -> plasmid_search (60_plasmid)
#
# WHERE DECONTAMINATION SITS PER MODE (decision D3 — v1 order preserved exactly):
#
#   mode      DRAFT_CONTIGS (screened)        DECONTAM_CONTIGS (written)
#   ────────  ──────────────────────────────  ──────────────────────────────────
#   illumina  SPAdes contigs_filt.fasta       = FINAL_CONTIGS (decontam is last)
#   contigs   filtered input contigs_filt     = FINAL_CONTIGS (decontam is last)
#   nanopore  reoriented {sample}_fixed       assembly_decontam.fasta -> Medaka
#   hybrid    the ILLUMINA SPAdes draft       contigs_sel.fasta -> Snippy ref +
#                                             the QC comparator genome
#
# That is why select_contigs writes DECONTAM_CONTIGS, not FINAL_CONTIGS: in
# nanopore and hybrid the delivered genome is produced LATER by Stage 4, and
# hard-coding FINAL_CONTIGS here would make nanopore circular
# (select -> final -> Medaka -> select).
#
# Everything referenced here is defined once in 00_common.smk and never
# re-derived: DRAFT_CONTIGS, DECONTAM_CONTIGS, FINAL_CONTIGS, DECONTAM_DIR,
# DECONTAM_BAM, BLOB_PREFIX/JSON/COV, BLOB_TABLE(_PREFIX), CONTIG_LIST,
# CONTIG_DECISIONS, COMPOSITION, BLASTOUT, PLASMID_BLASTOUT, TRIM_R1/TRIM_R2,
# FILT_LONG, BLASTDB, NT_VERSION, SELECT_TAXONOMY_SCRIPT, DECONTAMINATION, LOGS,
# CPUS, capped_cpus, and the capability flags.
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


# ── The read-mapping leg: one rule, three possible bodies ────────────────────
# map_contigs produces the SAME two files in every mode (DECONTAM_BAM + its .bai)
# but gets there with a different aligner, from a different read type, in a
# different conda environment. The choice is made HERE, at parse time, with a
# plain if/elif/else that defines exactly ONE rule body. MODE is fixed for the
# whole run, so a reader working in nanopore mode sees exactly one 15-line
# map_contigs and nothing else.
#
# Why not one rule with an input function? Because the branches differ in their
# `conda:` environment (bowtie.yaml vs minimap.yaml), and Snakemake resolves a
# rule's conda env when it deploys environments — an input function cannot reach
# it. The other alternative, a single rule with a union env carrying bowtie2 AND
# minimap2 plus a branchy shell, would make every mode build a bigger environment
# and turn a short shell into a conditional block. This also matches the house
# pattern already used for `if PHAGE_CALLER == "genomad":` in 60/70.

if HAS_SHORT_READS:

    # Bowtie2 needs its own index of the reference before it can map anything.
    # Module-local constant: both the producer (this rule) and the only consumer
    # (map_contigs, immediately below) live in this file, so the index prefix is
    # not a cross-module contract and does not belong in 00_common. It is built
    # off DECONTAM_DIR, so no stage number is re-derived.
    _BT2_PREFIX = DECONTAM_DIR + "/{sample}_contigs"

    # ── Rule: index_contigs — build the Bowtie2 index of the draft assembly ──
    # Takes in:  DRAFT_CONTIGS, the mode's draft assembly (from Stage 4).
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
        priority: 8
        shell:
            """
            bowtie2-build \
              -f {input.contigs} \
              {params.basename} > {log} 2>&1
            """

    # ── Rule: map_contigs (short-read modes) — reads back onto the draft ─────
    # Biology: mapping the sample's own trimmed reads onto its own assembly gives
    # per-contig read depth. A contig from a minor contaminant is usually covered
    # at a very different depth from the isolate's chromosome, which is the second
    # axis (alongside taxonomy) that BlobTools separates organisms on.
    #
    # Takes in: the six Bowtie2 index files (the DAG edge to index_contigs) and
    #           the fastp-trimmed pairs TRIM_R1/TRIM_R2 from the Stage-4 front end.
    # Does:     bowtie2 -> SAM, converted to BAM, coordinate-sorted, then indexed.
    # Produces: DECONTAM_BAM + .bai, both temp().
    # Consumed by: blob_json (the coverage leg) and map_evaluation (Qualimap, in
    #              shared/20_qc.smk). temp() keeps the BAM alive until both are done.
    #
    # v1->v2: `--write-index` is dropped from samtools sort. It wrote a .csi index
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
        priority: 8
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

    # ── Rule: map_contigs (nanopore) — ONT reads back onto the draft ─────────
    # Same purpose as the short-read version: a per-contig depth track for
    # BlobTools. minimap2's map-ont preset handles the higher error rate of raw
    # ONT reads, which bowtie2 cannot.
    #
    # Takes in: FILT_LONG (filtlong-filtered ONT reads) and DRAFT_CONTIGS (the
    #           reoriented Flye assembly), both from the Stage-4 front end.
    # Produces: DECONTAM_BAM + .bai, both temp().
    # Consumed by: blob_json and map_evaluation (Qualimap, shared/20_qc.smk).
    #
    # (v1 message: "--- Minimap2: Map reads against contigs. ---";
    #  v1 rule was also called map_contigs, its Qualimap rule map_qc -> D5 renames
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
        priority: 8
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

    # ── Rule: map_contigs (contigs mode) — a deliberate FAKE coverage track ──
    # There are no reads in this mode: the user hands us finished assemblies. But
    # `blobtools create -b` still wants a BAM, so we map the contigs against
    # THEMSELVES. The resulting depth is near-uniform and carries no information.
    #
    # This is preserved from v1 on purpose, and the consequence must be understood
    # before anyone tries to "use" this BAM: coverage-based separation in
    # BlobTools is MEANINGLESS in contigs mode. Only the taxonomy leg (blastn ->
    # BlobTools -> selector) is doing real work here. That is also why this mode
    # has no Qualimap rule — a mapping-quality report on a self-alignment would be
    # a chart of nothing.
    #
    # Takes in: DRAFT_CONTIGS (the header-fixed / length-filtered input contigs).
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
        priority: 8
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


# ── Rule: blast_contigs — taxonomic identity of every draft contig (megablast) ─
# Biology: megablast each contig against the NCBI nucleotide database and keep the
# top hits WITH their taxids and subject titles. This is the taxonomy leg that
# BlobTools turns into a per-contig genus call. megablast (not blastn) because we
# expect near-identical matches to known genomes, and it is far faster.
#
# Takes in: DRAFT_CONTIGS — the mode's draft assembly, i.e. the same contigs the
#           selector will filter. The screen is structurally pinned to the draft:
#           BlobTools must see every contig it is being asked to judge, including
#           the ones we are about to throw away.
# Does:     one blastn -task megablast per sample against {blast_db}/{nt_version}.
# Produces: BLASTOUT — 15 tab-separated columns, subject title (stitle) LAST.
# Consumed by: blob_json (all modes) and, in every mode EXCEPT hybrid, the
#              supplementary "does the nt hit say plasmid?" check in
#              shared/60_plasmid.smk (which greps that last column).
#
# DO NOT CHANGE THE -outfmt STRING. All four v1 modes agree on it, and the plasmid
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
    priority: 7
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


# ── Hybrid only: the same BLAST screen re-run on the DELIVERED genome ─────────
# In hybrid mode the decontamination screen above runs on the ILLUMINA draft,
# because that is what is being decontaminated, but Platon runs on the ONT genome.
# SPAdes names its contigs NODE_1_length_… and Flye names them contig_1, so
# plasmid_search's `grep -m 1 "$contig" <blastout>` would never match a single
# contig ID and every plasmid would be reported as "not verified by BLAST search".
# v1 BacFluxL+ solved this by running its own blastn inside plasmid_search; v2
# keeps the second BLAST but defines it here, next to the identical command it
# duplicates, and routes it through the PLASMID_BLASTOUT constant.
#
# Cost note: the long-read modes therefore BLAST against nt twice per sample
# (draft + final). For hybrid that is exactly what v1 did — not a new cost. For
# nanopore it is one extra blastn, buying immunity from an unverified assumption
# that Medaka preserves contig headers (see NEEDS_FINAL_BLAST in 00_common).
if NEEDS_FINAL_BLAST:

    # ── Rule: blast_final_contigs — nt screen of the DELIVERED genome ────────
    # Takes in: FINAL_CONTIGS — the genome the Stage-4 long-read front end
    #           actually delivers (hybrid: ONT+Polypolish; nanopore: the Medaka
    #           consensus) — NOT the draft that blast_contigs screened.
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
        priority: 7
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


# ── Rule: blob_json — join coverage + taxonomy into one BlobTools database ────
# Biology: BlobTools takes the assembly, the read-depth track and the BLAST hits
# and builds the "blobplot" database — per contig: length, GC, coverage, and a
# taxonomic assignment resolved through the NCBI taxonomy dump.
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
    priority: 7
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


# ── Rule: blob_table — collapse the hits into one call per contig ─────────────
# Biology: a contig usually has many BLAST hits pointing at several taxa. The
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
# This rule is completely mode-independent: same input shape, same command, same
# output, in all four modes.
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
    priority: 7
    shell:
        """
        blobtools view \
          --input {input.json} \
          --out {params.basename} \
          --taxrule bestsum \
          --rank all \
          --hits > {log} 2>&1
        """


# ── Rule: select_contigs — apply the decontamination policy, and record why ───
# Biology: this is where contigs are actually kept or dropped. The helper script
# reads the BlobTools table, resolves each contig to a genus, and applies the
# policy the user configured:
#   auto     — infer the isolate's genus from the assembly itself and keep it
#   include  — keep only the listed genera
#   exclude  — drop only the listed genera
#   off      — keep everything (still writes the audit files)
# plus discard_no_hit, which decides what happens to contigs nt could not place.
#
# Takes in:
#   bestscore = BLOB_TABLE    (the per-contig taxonomy + coverage table)
#   contigs   = DRAFT_CONTIGS (the sequences themselves)
# Produces (all four are the same in every mode except the last path):
#   abund     = COMPOSITION       "Genus: 0.87" lines -> read by Bakta (40) to
#                                 pick --genus
#   list      = CONTIG_LIST       the kept contig IDs, one per line
#   decisions = CONTIG_DECISIONS  the audit TSV: every contig, its genus, and the
#                                 REASON it was kept or dropped. Required by the
#                                 project convention that every filtering decision
#                                 is auditable (CLAUDE.md).
#   contigs   = DECONTAM_CONTIGS  the kept sequences (see the per-mode table in
#                                 this file's banner for who consumes them next)
#
# conda: NONE — inherited from all four v1 modes. The selector is stdlib-only
# Python and runs in the environment Snakemake was launched from. Adding an env
# would be new behaviour; it is the same deferred decision as cazyme_db_download
# (see docs/README_notes.md item 3).
#
# All 14 selector flags are passed verbatim from the DECONTAMINATION dict resolved
# once in 00_common, with {...:q} quoting so a genus list containing spaces or a
# path with odd characters survives the shell intact.
#
# v1->v2 (additive): v1 had no log:, so the selector's warnings — notably
# "WARNING: N contigs from the BlobTools table were not found in the FASTA" and
# the per-sample kept/dropped counts — went to the console and were lost on a
# large batch. They now go to the log file instead. Trade-off worth knowing: you
# have to open the log to see them.
#
# v1->v2 (removed): FastaFlux used to re-linearise the FASTA into
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
    priority: 6
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
