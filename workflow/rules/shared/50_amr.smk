# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — Stage 05 AMR module (rules/shared/50_amr.smk)
#
# This module screens each finished genome for antimicrobial-resistance (AMR) and
# virulence genes on the ASSEMBLED CONTIGS, using ABRicate against several curated
# databases, then collates the per-database hits into one summary table. Like the
# rest of the shared tail it consumes only FINAL_CONTIGS (D2), so it is identical
# in every mode.
#
# Scope note: BacFlux detects AMR on THREE complementary legs (that is deliberate
# — different methods catch different things). This file holds TWO of them:
#   * the assembly/contig-based ABRicate leg, identical in all four modes; and
#   * the read-based CARD mapping leg (download_amr_db + map_amr_db) at the bottom
#     of this file, gated on HAS_SHORT_READS because it needs raw reads. Mapping
#     READS against CARD is immune to the assembly collapsing repeated genes, so
#     it catches resistance genes the contig leg can miss.
# (The third leg, AMRFinderPlus, belongs to the future mobilome module — see
#  docs/mobilome_module_SPEC.md §4.)
#
# Data flow:
#
#   contigs_final.fasta ─► amr_contigs (once per database) ─► abricate/{sample}/{db}.tsv
#                                                                      │  (all dbs)
#                                              AMR_summary ◄───────────┘
#                                                   └─► abricate/{sample}/AMR_summary.txt
#
#   [short-read modes only]
#   CARD tarball ─► download_amr_db ─► card_db/ ─┐
#                                                 ├─► map_amr_db ─► mapping/{sample}/
#   {sample}_trim_R{1,2}.fastq (from Stage 4) ───┘        {sample}_covstats.tsv
#                                                        {sample}_AMR_legend.tsv
#
# Inherited from 00_common.smk (never re-derived here): FINAL_CONTIGS, DIR_AMR,
# LOGS, DATABASES, DATABASE_PATTERN, and — for the CARD leg — CARD_TARBALL,
# CARD_DB_DIR, TRIM_R1, TRIM_R2, RAM, capped_cpus, HAS_SHORT_READS. The
# 8-database list lives ONCE in 00_common and drives (a) the per-db fan-out in
# rule all, (b) the {db} wildcard constraint, and (c) the summary input — one
# list, one source of truth (D7).
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/abricate.yaml" climbs shared/ -> rules/ -> workflow/ and lands on
# workflow/envs/abricate.yaml (the single shared env copy).
#
# Resource convention: these ABRicate rules are single-threaded in v1 and declare
# no `threads:` — kept as-is (no thread request added), so Snakemake counts them
# as 1 core each.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: amr_contigs — AMR / virulence screen, one run per database (ABRicate) ──
# Biology: ABRicate BLASTs the assembled contigs against a curated resistance /
# virulence gene database and reports every hit passing the identity & coverage
# thresholds. Different databases (CARD, ResFinder, NCBI, VFDB, …) cover
# different gene sets, so the SAME contigs are screened against each one
# separately and the results kept side by side.
#
# Takes in: contigs = FINAL_CONTIGS — the finished, decontaminated assembly.
# Does: one ABRicate run for the database named by the {db} wildcard, applying the
#       EFSA reporting thresholds (>=80% identity, >=70% coverage). These
#       thresholds are kept ONLY on this contig-based leg; the read-based CARD leg
#       (Stage-3) uses its own, stricter mapping criteria, and AMRFinderPlus (see
#       the mobilome spec) would use its curated per-gene cutoffs — applying a
#       blanket 80/70 across all three would fight those.
# Produces: 05.amr/abricate/{sample}/{db}.tsv — one hit table per (sample, db).
# Consumed by: AMR_summary (below), which collates all per-db tables for a sample.
#
# FAN-OUT (unchanged mechanism, relocated source): this rule declares ONE
# templated output and contains NO expand(). The fan-out over the 8 databases
# lives in 00_common._downstream_targets(), which requests
#   expand(DIR_AMR + "/abricate/{sample}/{db}.tsv", sample=SAMPLES, db=DATABASES).
# rule all therefore asks for one {db}.tsv per (sample x database), and Snakemake
# instantiates this rule once per pair. wildcard_constraints bounds {db} to the 8
# legal names (DATABASE_PATTERN); params.db feeds each instance its own database.
# In v1 that expand() sat in the Snakefile's own `rule all`; in v2 it is already
# in 00_common, so nothing extra is needed here.
#
# (v1 message: "--- ABRicate: AMR detection ({wildcards.db}). ---")
rule amr_contigs:
    input:
        contigs = FINAL_CONTIGS,
    output:
        amr_tab = DIR_AMR + "/abricate/{sample}/{db}.tsv",
    wildcard_constraints:
        # Restrict {db} to the 8 known database names so this rule cannot try to
        # match, e.g., the AMR_summary.txt path below.
        db = DATABASE_PATTERN,
    params:
        db = lambda wc: wc.db,
        minid = 80,
        mincov = 70,
    conda:
        "../../envs/abricate.yaml"
    log:
        LOGS + "/amr_{db}_in_{sample}_contigs.log"
    priority: 4
    shell:
        # :q lets Snakemake shell-quote each value safely (paths / db name / log).
        """
        abricate \
          --db {params.db:q} \
          {input.contigs:q} \
          --minid {params.minid} \
          --mincov {params.mincov} \
          --nopath \
          --quiet > {output.amr_tab:q} 2> {log:q}
        """


# ── Rule: AMR_summary — collate the per-database hits into one matrix (ABRicate) ─
# Biology: turns a sample's several per-database hit tables into a single
# presence/absence summary (gene x database) — the human-readable overview of the
# whole AMR / virulence screen for that isolate.
#
# Takes in: all per-database {db}.tsv tables for THIS sample. The double-brace
#           {{sample}} keeps `sample` a wildcard while expand() fills in only the
#           `db` list, so we get exactly this sample's 8 tables, in DATABASES order.
# Does: ABRicate's --summary over that file list, passed POSITIONALLY as {input}.
# Produces: 05.amr/abricate/{sample}/AMR_summary.txt.
# Consumed by: the report / the user (a terminal AMR product).
#
# D7 single-source: v1 hard-coded all 8 database filenames TWICE (as 8 named
# input: keys AND again in the shell). Here input: is driven by the shared
# DATABASES list and the shell passes {input} positionally, so adding or removing
# a database in 00_common.DATABASES now propagates automatically to rule all, the
# {db} wildcard constraint, and this summary. `abricate --summary` accepts the
# files as a plain positional list; their order is the DATABASES order.
#
# DELIBERATE v1→v2 improvement (additive, no behavior change): v1 AMR_summary had
# NO log: and redirected only stdout, so ABRicate's summary-step stderr was lost.
# v2 adds a log and captures stderr, matching every other rule in these two
# modules. The ABRicate command itself is unchanged.
rule AMR_summary:
    input:
        # Unnamed positional list: this sample's per-database tables, DATABASES order.
        # `abricate --summary` column order follows input-file order, so it tracks
        # DATABASES order — do not reorder DATABASES in 00_common expecting the same
        # summary layout.
        expand(DIR_AMR + "/abricate/{{sample}}/{db}.tsv", db=DATABASES),
    output:
        amr_summary = DIR_AMR + "/abricate/{sample}/AMR_summary.txt",
    conda:
        "../../envs/abricate.yaml"
    log:
        LOGS + "/AMR_summary_{sample}.log"
    priority: 3
    shell:
        """
        abricate --summary {input} > {output.amr_summary} 2> {log}
        """


# ── CARD read-mapping leg — short-read modes only ────────────────────────────
# Gated on the CAPABILITY flag HAS_SHORT_READS (illumina + hybrid), not on the
# mode name (D7). This is the SAME flag _downstream_targets() in 00_common uses to
# request these two files, so rule-all and the rule definitions can never
# disagree: either both exist or neither does. nanopore and contigs get no CARD
# leg at all — v1 had none for them, and neither mode's config even declares
# links.card_link. Because envs/bbmap.yaml is referenced only from inside this
# block, those two modes never build it either.
#
# Why a read-based leg at all, when ABRicate already screened the contigs? Because
# assembly COLLAPSES repeats. A resistance gene present in several copies, or one
# sitting on a repeat-rich mobile element that broke the assembly, can be
# under-represented (or missing) in the contigs while still being plainly present
# in the reads. Mapping the trimmed reads straight onto CARD's reference sequences
# side-steps the assembly entirely.
if HAS_SHORT_READS:

    # ── Rule: download_amr_db — fetch the CARD reference sequences ───────────
    # Takes in: nothing — a pure download from the URL in links.card_link.
    # Does:     wget the .tar.bz2, then extract it.
    # Produces: card_tarball (temp) and card_dir (temp DIRECTORY) — both are
    #           deleted once the last map_amr_db job has finished, since nothing
    #           downstream needs the raw database.
    # Consumed by: map_amr_db (every sample waits on this one download).
    #
    # params.link is read from the config INSIDE the rule, deliberately: 00_common
    # §7 leaves card_link and phix_link to their gated consumers, so a nanopore or
    # contigs config that omits card_link still parses cleanly.
    #
    # Structural change from v1: v1 declared the tarball INSIDE the directory()
    # output of the same rule (06.AMR/AMR_db/card.tar.bz2 inside
    # directory("06.AMR/AMR_db")). They are SIBLINGS here — the house rule stated
    # in 40_annotation.smk is that temp scratch never nests inside a directory
    # output. Because of that the extraction directory now has to be created
    # explicitly (v1 got it for free as the tarball's parent).
    #
    # conda: NONE — uses wget and tar from the launch environment, as v1 did. Same
    # deferred decision as cazyme_db_download; see docs/README_notes.md item 3.
    #
    # (v1 message: "--- Download AMR features from CARD repository. ---")
    rule download_amr_db:
        output:
            card_tarball = temp(CARD_TARBALL),
            card_dir = temp(directory(CARD_DB_DIR)),
        params:
            # Resolved and validated once in 00_common (§7), so a missing key is
            # reported by name at parse time rather than as a bare KeyError here.
            link = CARD_LINK,
        log:
            LOGS + "/download_amr.log"
        priority: 9
        shell:
            """
            mkdir -p {output.card_dir}

            wget {params.link} -O {output.card_tarball} > {log} 2>&1
            tar -xjvf {output.card_tarball} -C {output.card_dir} >> {log} 2>&1
            """

    # ── Rule: map_amr_db — map trimmed reads onto CARD (BBMap) ───────────────
    # Biology: align this sample's quality-trimmed Illumina pairs against CARD's
    # protein-homolog-model nucleotide sequences at >=99% identity, and report how
    # much of each reference gene was actually covered by reads. A resistance gene
    # is called present when a large fraction of its length is covered — length
    # coverage, not just "some reads hit it", which is what keeps short conserved
    # domains from producing false positives.
    #
    # Takes in:
    #   r1 / r2  = TRIM_R1 / TRIM_R2, the fastp-trimmed pairs written by the
    #              Stage-4 illumina / hybrid front end. Both short-read modes write
    #              them to the same path, hence one shared constant. They are
    #              declared temp() by their producer, so Snakemake keeps them until
    #              the assembler, map_contigs and this rule are all done.
    #   card_dir = the extracted CARD database from download_amr_db.
    # Does: BBMap at idfilter=0.99, then two post-processing steps —
    #       (1) re-sort covstats by descending Covered_percent, header preserved;
    #       (2) build a human-readable legend: for every feature covered >=70%,
    #           pull its row out of CARD's aro_index.tsv so the output names the
    #           drug class and mechanism rather than an ARO accession.
    # Produces:
    #   covstats     = 05.amr/mapping/{sample}/{sample}_covstats.tsv
    #   amr_legend   = 05.amr/mapping/{sample}/{sample}_AMR_legend.tsv
    #   (plus two temp() intermediates: BBMap's ref/ index and the unsorted
    #    covstats). Those two kept paths are byte-identical to what
    #    _downstream_targets() already requests.
    # Consumed by: the user (terminal AMR products; not fed into MultiQC).
    #
    # PRESERVED FROM v1, ALL KNOWN WARTS, DELIBERATELY NOT FIXED:
    #   * path={output.bbmap_temp} where bbmap_temp already ends in /ref, so BBMap
    #     creates ref/ref/. Odd, harmless, kept.
    #   * the legend's `grep $i {card_dir}/aro_index.tsv` is unquoted, unanchored
    #     and not -F, so an ARO accession that is a substring of another one can
    #     pull in extra rows. Kept as-is.
    #   * NO `set -euo pipefail` in this rule. The legend's `for ... grep ... done`
    #     loop legitimately exits non-zero when nothing clears 70%, and `set -e`
    #     would turn "this isolate has no strongly covered AMR gene" — a perfectly
    #     normal result — into a hard pipeline failure.
    #
    # One deliberate v1->v2 change: -Xmx32g was hard-coded, which fails outright on
    # a machine with less RAM than that. It now asks for min(RAM, 32) GB, the same
    # class of fix as Qualimap's --java-mem-size in shared/20_qc.smk (both taken
    # together, for consistency). Like Qualimap's, this is a GB figure passed
    # through to the JVM, not something the scheduler reserves.
    #
    # (v1 message: "--- Map trimmed reads against CARD db. ---")
    rule map_amr_db:
        input:
            card_dir = CARD_DB_DIR,
            r1 = TRIM_R1,
            r2 = TRIM_R2,
        output:
            bbmap_temp = temp(directory(DIR_AMR + "/mapping/{sample}/ref")),
            covstats_temp = temp(DIR_AMR + "/mapping/{sample}/{sample}_covstats_temp.tsv"),
            covstats = DIR_AMR + "/mapping/{sample}/{sample}_covstats.tsv",
            amr_legend = DIR_AMR + "/mapping/{sample}/{sample}_AMR_legend.tsv",
        params:
            # CARD ships several models; the protein homolog model is the one that
            # holds acquired resistance genes (not the mutation-based models).
            card_target = "nucleotide_fasta_protein_homolog_model.fasta",
            min_id = 0.99,
            max_ram = min(RAM, 32),
        conda:
            "../../envs/bbmap.yaml"
        threads: capped_cpus(24)
        log:
            LOGS + "/map_amr_{sample}.log"
        priority: 5
        shell:
            # Column 5 of BBMap's covstats is Covered_percent; the >=70 threshold
            # and the legend header text below are v1's, kept verbatim. The awk
            # braces are doubled because this awk is written directly in the shell.
            """
            bbmap.sh \
              -in={input.r1} \
              -in2={input.r2} \
              ref={input.card_dir}/{params.card_target} \
              path={output.bbmap_temp} \
              idfilter={params.min_id} \
              idtag \
              -Xmx{params.max_ram}g \
              threads={threads} \
              ambiguous=best \
              secondary=f \
              covstats={output.covstats_temp} > {log} 2>&1

            (head -n 1 {output.covstats_temp} > {output.covstats}) && \
            tail -n +2 {output.covstats_temp} | awk -F'\t' '{{print $5 "\t" $0}}' | sort -t$'\t' -k1,1nr | cut -f2- >> {output.covstats}

            echo "#AMR features with a covered length of at least 70%" > {output.amr_legend}
            (head -n 1 {input.card_dir}/aro_index.tsv >> {output.amr_legend}) && \
            for i in $(tail -n +2 {output.covstats} | awk -F'\t' '$5 >=70' | cut -f1 | awk -F'|' '{{print $5}}'); do \
                grep $i {input.card_dir}/aro_index.tsv; \
            done >> {output.amr_legend}
            """
