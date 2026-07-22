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
# — different methods catch different things). This file holds ONLY the
# assembly/contig-based ABRicate leg, which is the same in all four modes. The
# read-based CARD mapping leg (download_amr_db + map_amr_db) needs raw reads, so
# it exists only in the short-read modes and is added later in Stage 3 (still in
# this file, gated on HAS_SHORT_READS). It is intentionally NOT here yet.
#
# Data flow:
#
#   contigs_final.fasta ─► amr_contigs (once per database) ─► abricate/{sample}/{db}.tsv
#                                                                      │  (all dbs)
#                                              AMR_summary ◄───────────┘
#                                                   └─► abricate/{sample}/AMR_summary.txt
#
# Inherited from 00_common.smk (never re-derived here): FINAL_CONTIGS, DIR_AMR,
# LOGS, DATABASES, DATABASE_PATTERN. The 8-database list lives ONCE in 00_common
# and drives (a) the per-db fan-out in rule all, (b) the {db} wildcard constraint,
# and (c) the summary input — one list, one source of truth (D7).
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/abricate.yaml" climbs shared/ -> rules/ -> workflow/ and lands on
# workflow/envs/abricate.yaml (the single shared env copy).
#
# Resource convention: these ABRicate rules are single-threaded in v1 and declare
# no `resources: cpus` — kept as-is (no thread request added).
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
