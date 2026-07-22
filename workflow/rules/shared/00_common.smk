# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — shared foundation module (rules/shared/00_common.smk)
#
# This is the first file the Snakefile includes and the one place that knows
# about MODE. Everything else — the mode-specific front ends and the shared
# downstream rules — reads the names defined here and never re-derives them.
#
# What this module does, top to bottom:
#   1. read and validate config["mode"] (illumina | nanopore | hybrid | contigs)
#   2. set the working directory and resolve the database / helper paths
#   3. expose resource accessors and mode-capability flags used for gating
#   4. define the unified stage-number layout (D1) and the single canonical
#      final-assembly hand-off file (D2)
#   5. define every de-duplicated helper that used to be copy-pasted into the
#      four separate Snakefiles (locus tags, decontamination settings, boolean
#      coercion, Medaka/Flye resolution, CheckV/dbCAN link parsing)
#   6. discover the samples for the ACTIVE mode only — the glob and the
#      "no input / bad name" exits for the other three modes never run
#   7. build the rule-all target list
#
# House-rule reminder for anyone editing this file: every parse-time side effect
# (a print, a glob of the input directory, a sys.exit guard) that depends on the
# input data MUST sit behind the MODE dispatch, so running one mode never fires
# another mode's discovery or error. Mode-INDEPENDENT setup (database links,
# decontamination policy) may run unconditionally.
# ─────────────────────────────────────────────────────────────────────────────

import glob
import os
import re
import sys
from collections import namedtuple
from pathlib import Path
from urllib.parse import urlparse

# glob_wildcards is a Snakemake helper; import it explicitly so the reader can
# see where it comes from rather than relying on an injected global.
from snakemake.io import glob_wildcards


# ─────────────────────────── 1. Mode dispatch ───────────────────────────────
# The config is REQUIRED on the command line (Snakefile.v2 deliberately declares
# no default `configfile:` — see the note there). Fail with an actionable message
# rather than a bare KeyError if nothing was supplied.
if not config:
    sys.exit(
        "[BacFlux] No configuration supplied. Pass one explicitly, e.g.:\n"
        "  snakemake --sdm conda --cores N "
        "--configfile config/config_v2.yaml"
    )

# The mode is REQUIRED. Bracket access means a config with no "mode" fails
# immediately with a clear KeyError rather than silently guessing a pipeline.
# .lower() normalises case so "Illumina" and "illumina" behave the same.
MODE = config["mode"].lower()
if MODE not in ("illumina", "nanopore", "hybrid", "contigs"):
    sys.exit(
        "config.mode must be one of illumina | nanopore | hybrid | contigs "
        f"(got: {config['mode']!r})"
    )

# One friendly line per mode, printed once at parse time so the log header says
# which pipeline is running. Purely cosmetic.
MODE_TAGLINE = {
    "illumina": "Genomic analysis of bacterial Illumina reads",
    "nanopore": "Genomic analysis of bacterial ONT reads",
    "hybrid":   "Illumina-guided ONT assembly and analysis of bacterial genomes",
    "contigs":  "Downstream analysis of bacterial whole-genome assemblies",
}
print(f"Mode: {MODE} — {MODE_TAGLINE[MODE]}.")


# ─────────────────── 1b. Phage/virus caller selection (D8) ───────────────────
# Which tool calls viruses/prophages. VirSorter2 is the DEFAULT; geNomad is a
# SELECTABLE opt-in (config.phage.caller: genomad). Resolved here next to the MODE
# dispatch (it is mode-independent), following the same validate-or-exit pattern.
#
# LICENSING — why geNomad is opt-in, NOT the default (decided 2026-07-22):
# geNomad is under a Berkeley Lab ACADEMIC / NON-COMMERCIAL-USE-ONLY licence
# ("User must be an accredited academic institution"; commercial use reserved).
# BacFlux is MIT and must not force a non-commercial restriction on downstream
# users, so the DEFAULT pipeline uses only permissively-licensed tools
# (VirSorter2 GPLv2 + Platon GPLv3 + CheckV LBNL-BSD, all commercial-use-OK).
# Selecting geNomad is the user's own informed opt-in — the same treatment BacFlux
# gives licence-encumbered databases.
#
# What the choice controls: geNomad does viruses AND plasmids in one run, so
# opting in turns on BOTH (a) geNomad as the virus caller here, and (b) the
# Platon+geNomad plasmid concordance in 60_plasmid.smk. On the default the plasmid
# stage is Platon-only and geNomad never runs. CheckV runs downstream of whichever
# virus caller ran, either way. `config.get("phage", {}) or {}` guards a phage
# block that is present but null.
PHAGE_CALLER = str((config.get("phage", {}) or {}).get("caller") or "virsorter2").strip().lower()
if PHAGE_CALLER not in ("genomad", "virsorter2"):
    sys.exit(f"config.phage.caller must be genomad | virsorter2 (got: {PHAGE_CALLER!r})")
if PHAGE_CALLER == "genomad":
    print(
        "Phage caller: geNomad (opt-in). NOTE: geNomad is licensed for ACADEMIC / "
        "NON-COMMERCIAL use only (Berkeley Lab). Also enables the Platon+geNomad "
        "plasmid concordance."
    )
else:
    print("Phage caller: virsorter2 (default). Plasmid stage: Platon-only (geNomad off).")


# ─────────────────── 2. Output root, databases, paths ────────────────────────
# MetaFlux-style: resolve the output directory to one absolute root (OUT) and
# build every stage path off it, instead of Snakemake's `workdir:` directive.
#
# Why not workdir (v1 BacFlux's approach, and Stage 1's first draft of this
# file): workdir: changes the process's working directory at PARSE time —
# before any rule runs, even for a dry run — which caused a real problem during
# Stage 1 review: a placeholder output_dir that doesn't exist/isn't writable
# crashed the parse before Snakemake could even build the DAG. It also forces
# every input directory in the config to be absolute too, since a relative one
# would silently resolve under output_dir instead of wherever the user launched
# from. Resolving OUT once here avoids both: the working directory never
# changes, so input paths behave the way a user naturally expects (relative to
# the launch dir, or absolute if given), and only output_dir strictly needs to
# be a real, writable location — checked when a rule first writes there, not at
# parse time for every dry run.
#
# Kept as a plain STRING (not a pathlib.Path object), even though MetaFlux uses
# Path objects: every stage constant below and every rule this project writes
# builds paths with "+" / os.path.join on plain strings, and switching to Path
# objects would mean rewriting that "+" style to Path's "/" everywhere, in this
# file and every rule module still to come. An absolute string root gets the
# actual goal (no workdir side effect, no forced-absolute inputs) without that
# churn.
OUT = str(Path(config["directories"]["output_dir"]).resolve())

# Database directories. Required for every mode (each mode annotates, assigns
# taxonomy, screens for plasmids, etc.), so bracket access is intentional: a
# missing database key should stop the run, not default to something wrong.
BAKTADB  = config["directories"]["bakta_db"]
BLASTDB  = config["directories"]["blast_db"]
DMNDDB   = config["directories"]["eggnog_db"]
GTDBTKDB = config["directories"]["gtdbtk_db"]
PLATONDB = config["directories"]["platon_db"]

# CheckV is the ONE database BacFlux can also fetch for itself, so unlike the five
# above its path is OPTIONAL — hence .get() rather than bracket access.
#
# Why it exists: the official CheckV database lives on portal.nersc.gov, which goes
# down often enough to cost real time (it was unreachable for the whole of
# 2026-07-22, blocking three validation runs). Pointing this key at a copy you
# already hold on disk removes that dependency completely.
#
# When set, it takes PRECEDENCE over links.checkv_link and over CheckV's own
# downloader: nothing is fetched.
#
# BacFlux does NOT hand this path straight to CheckV, for two reasons found the
# hard way on a real shared database:
#
#   1. CheckV needs a DIAMOND index (genome_db/checkv_reps.dmnd) that the official
#      archive does not ship — it is built locally after unpacking. So whatever
#      index sits in a shared database was built by whichever DIAMOND that site
#      happened to have. DIAMOND's database format is versioned, and 2.0.4 cannot
#      run blastp against a format-1 index built by an older build: CheckV dies at
#      "[3/8] Running DIAMOND blastp search... DIAMOND task failed", AFTER the
#      contamination stage has already succeeded, which makes it look like a
#      CheckV bug rather than an index mismatch.
#   2. A shared database is typically NOT writable by the person running the
#      workflow (and must not be rewritten anyway — other people rely on it), so
#      "just rebuild the index in place" is not available.
#
# So rule checkv_db_local (70_phage.smk) builds a local VIEW instead: symlinks to
# the big read-only files, plus a DIAMOND index built by THIS workflow's own
# DIAMOND. Costs ~950 MB and a couple of minutes once per output directory, and
# makes the feature independent of who built the shared index, or with what.
# Your database is only ever read.
#
# Point it at the PARENT directory holding the versioned DB folder, i.e. the same
# shape BacFlux would have created itself:
#     /path/to/checkv/                 <- put THIS in the config
#         checkv-db-v1.5/
#             genome_db/  hmm_db/  README.txt
CHECKVDB = str((config["directories"].get("checkv_db") or "")).strip()

# workflow.basedir is the absolute path of the workflow/ directory (independent
# of workdir). It anchors the helper script and the on-disk rule-module lookup.
WORKFLOW_DIR = workflow.basedir

# Shared post-BlobTools helper: reads contig taxonomy assignments and decides
# which contigs are kept, discarded, or left untouched. Consumed by the shared
# select_contigs rule in every mode.
SELECT_TAXONOMY_SCRIPT = os.path.join(WORKFLOW_DIR, "scripts", "select_contigs_by_taxonomy.py")

# D9 plasmid-concordance helper: joins Platon's plasmid calls and geNomad's plasmid
# calls into one confidence-tiered TSV. Stdlib-only Python; consumed by the
# plasmid_concordance rule in shared/60_plasmid.smk (run in the platon env, which
# already ships a Python — so this adds no dependency).
PLASMID_CONCORDANCE_SCRIPT = os.path.join(WORKFLOW_DIR, "scripts", "plasmid_concordance.py")


# ───────────────────────── 3a. Resource accessors ───────────────────────────
# Single source of truth for compute limits. .get with a default means a config
# that omits a key still parses (nanopore/contigs configs legitimately omit
# ram_gb, since only SPAdes and the JVM tools actually consume large RAM).
_resources = config.get("resources", {}) or {}
CPUS = int(_resources.get("threads", 16))
RAM  = int(_resources.get("ram_gb", 64))


def capped_cpus(n):
    # Cap a rule's thread request at n. Some tools stop scaling (or misbehave)
    # past a certain thread count; rules call capped_cpus(N) instead of spelling
    # out min(CPUS, N) in dozens of places.
    return min(CPUS, n)


# Which NCBI nucleotide database subfolder the contamination BLAST screens
# against. Centralised here (was read inline in blast_contigs) with a safe
# default so a config that omits it still parses.
NT_VERSION = str((config.get("parameters", {}) or {}).get("nt_version") or "core_nt").strip()


# ─────────────────────── 3b. Mode-capability flags ──────────────────────────
# These booleans — not the mode name — drive every downstream gate, so the
# intent reads plainly ("this leg needs short reads") rather than enumerating
# modes at each call site.
IS_HYBRID       = MODE == "hybrid"
HAS_SHORT_READS = MODE in ("illumina", "hybrid")   # phix removal, fastp, SPAdes, CARD read-mapping leg
HAS_LONG_READS  = MODE in ("nanopore", "hybrid")   # filtlong, NanoPlot, Flye, Medaka, dnaapler
HAS_READS       = MODE != "contigs"                # anything populating 01.reads/ at all


# ──────────────── 4. Unified stage layout (D1) + hand-off (D2) ───────────────
# In v1 the shared downstream stages landed at different numbers per mode
# (taxonomy was 03/04/09 depending on how many front-end stages preceded it).
# v2 groups ALL tech-specific front-end work under two fixed parents (01.reads,
# 02.assembly) so every shared stage below has the same number in every mode.
# Every constant is now an ABSOLUTE path (built off OUT, section 2) rather than
# relative-under-workdir. Rules reference these constants and append their own
# "/{sample}/..." — no rule hard-codes a stage number or the output root.
DIR_READS      = OUT + "/01.reads"       # read QC + filtering; sub-dirs illumina/ , ont/  (empty in contigs mode)
DIR_ASSEMBLY   = OUT + "/02.assembly"    # assembly, polishing, reorientation, decontamination, assembly QC
DIR_TAXONOMY   = OUT + "/03.taxonomy"    # GTDB-Tk
DIR_ANNOTATION = OUT + "/04.annotation"  # bakta/ , eggnog/ , antismash/ , dbcan/
DIR_AMR        = OUT + "/05.amr"         # abricate/ , mapping/  (mapping only in short-read modes)
DIR_PLASMIDS   = OUT + "/06.plasmids"    # platon + geNomad concordance
DIR_PHAGES     = OUT + "/07.phages"      # geNomad|virsorter2 caller + checkv
DIR_MOBILOME   = OUT + "/08.mobilome"    # mobilome module, only when config.mobilome.run
DIR_REPORT     = OUT + "/09.report"      # multiqc

# Cross-cutting output locations, sitting alongside the numbered stages rather
# than inside any one of them. Needed as of this file (not before): v1's rules
# wrote log:/benchmark: as bare relative strings (e.g. "logs/foo_{sample}.log"),
# which only landed under output_dir because workdir: put the process there.
# Now that workdir: is gone (see section 2), every rule must build these
# explicitly off an absolute root, or a log would silently land wherever the
# user happened to launch from instead of next to the run's actual output.
LOGS  = OUT + "/logs"
BENCH = OUT + "/benchmarks"

# D2 — the ONE canonical hand-off. Every mode's front end ends by writing this
# exact file, and every rule in 03.* – 08.* consumes ONLY this file, so the
# downstream half never needs to know which front end (or assembler) produced
# it. Built from DIR_ASSEMBLY so the stage number lives in exactly one place.
FINAL_CONTIGS = DIR_ASSEMBLY + "/{sample}/contigs_final.fasta"

# Second cross-stage hand-off, same single-source philosophy as FINAL_CONTIGS:
# the genus-composition table. It is WRITTEN by the decontamination selector
# (select_contigs, in shared/10_decontam.smk) and READ by the Bakta annotation
# rule to pick the isolate's most likely genus. Defining the one canonical path
# here means producer and consumer can never drift onto different strings (the
# failure the Stage-2a review flagged). Lives under 02.assembly/ next to the
# contaminant-screening outputs, per the D1 layout.
#
# CONTENT (verified against workflow/scripts/select_contigs_by_taxonomy.py, not
# inferred): one line per genus, written as
#     Genus: <relative frequency>
# i.e. a COLON followed by a SPACE, then a fraction in 0.00-1.00 with 2 decimals
# — NOT a percentage. Lines are sorted by descending contig count then genus
# name; counts cover ALL contigs in the BlobTools table (not only the kept ones);
# contigs with no taxonomic hit appear under the literal genus name "no-hit".
# The Bakta rule's `sort -t':' -k2 -nr | cut -d':' -f1 | sed -n '1p'` idiom reads
# exactly this shape (numeric sort tolerates the leading space).
COMPOSITION = DIR_ASSEMBLY + "/{sample}/contaminants/{sample}_composition.txt"

# Shared antiSMASH reference-database directory: downloaded once by
# secondary_metabolites_db and read by every per-sample secondary_metabolites_
# analysis run. Single-sourced here (like DBCAN_DB_DIR) so the download rule's
# output and the analysis rule's input reference the exact same path.
ANTISMASH_DB_DIR = DIR_ANNOTATION + "/antismash/databases"

# ── Plasmid (06) + phage (07) stage paths (D8/D9) ────────────────────────────
# When the user OPTS IN to geNomad (PHAGE_CALLER == "genomad"), ONE end-to-end run
# per sample (defined in shared/70_phage.smk) produces BOTH the virus calls
# (→ CheckV) AND the plasmid calls (→ the D9 Platon/geNomad concordance). On the
# default (virsorter2) geNomad does not run at all, and the geNomad constants below
# are simply never referenced by any active rule. 60_plasmid.smk and 70_phage.smk
# consume the files below ONLY through these constants — the same anti-drift,
# single-source rule as FINAL_CONTIGS/COMPOSITION: no rule re-derives a path.
#
# GENOMAD_PREFIX is the basename of FINAL_CONTIGS with the .fasta stripped, i.e.
# the fixed literal "contigs_final". BOTH geNomad and Platon name every output
# file after their input's basename, so this one constant names both tools'
# per-sample outputs (hence the deliberately generic reuse in 60_plasmid.smk).
GENOMAD_PREFIX = os.path.splitext(os.path.basename(FINAL_CONTIGS))[0]   # "contigs_final" — derived so it can't drift if FINAL_CONTIGS is renamed (the {sample} token lives in the directory part, not the basename); geNomad & Platon both name outputs after the input basename

# geNomad — shared virus+plasmid caller (rule genomad_end_to_end, 70_phage.smk).
GENOMAD_DB_DIR          = DIR_PHAGES + "/genomad_db"            # produced by rule genomad_db
GENOMAD_DIR             = DIR_PHAGES + "/genomad/{sample}"      # produced by rule genomad_end_to_end (a DIRECTORY)
GENOMAD_VIRUS_FASTA     = GENOMAD_DIR + "/" + GENOMAD_PREFIX + "_summary/" + GENOMAD_PREFIX + "_virus.fna"
GENOMAD_PLASMID_SUMMARY = GENOMAD_DIR + "/" + GENOMAD_PREFIX + "_summary/" + GENOMAD_PREFIX + "_plasmid_summary.tsv"

# VirSorter2 — selectable alternate virus caller (only when PHAGE_CALLER=="virsorter2").
VS2_DB_DIR = DIR_PHAGES + "/vs2_db"                             # produced by rule virsorter2_db
VS2_DIR    = DIR_PHAGES + "/virsorter/{sample}"                 # produced by rule viral_identification_virsorter2 (a DIRECTORY)

# CheckV — completeness/contamination QC of whichever caller's virus calls.
# ALWAYS this path: "the CheckV database this run uses". Exactly one of two
# mutually exclusive rules in 70_phage.smk produces it —
#   * checkv_db       downloads a database into it (no directories.checkv_db), or
#   * checkv_db_local builds a symlink VIEW of the user's own database into it.
# viral_quality consumes this one name either way and resolves the versioned
# sub-folder at runtime, so nothing downstream knows or cares which it got.
CHECKV_DB_DIR = DIR_PHAGES + "/checkv_db"

# Platon — primary plasmid caller (rule plasmid_search, 60_plasmid.smk). v2 moves
# Platon's output into a platon/ sub-dir (v1 wrote it straight into {sample}/) so
# the concordance TSV can sit as a clean SIBLING — otherwise two rules would write
# inside one rule's directory() output, which Snakemake forbids.
PLATON_DIR          = DIR_PLASMIDS + "/{sample}/platon"                             # a DIRECTORY (rule plasmid_search)
PLATON_TABLE        = PLATON_DIR + "/" + GENOMAD_PREFIX + ".tsv"                    # Platon per-plasmid table (contigs_final.tsv)
PLATON_CHROMOSOME   = PLATON_DIR + "/" + GENOMAD_PREFIX + ".chromosome.fasta"       # Platon chromosome contigs
PLATON_VERIFIED     = PLATON_DIR + "/verified_plasmids.txt"                         # kept v1 BLAST-text check (supplementary)
PLASMID_CONCORDANCE = DIR_PLASMIDS + "/{sample}/{sample}_plasmid_concordance.tsv"   # D9 terminal deliverable (rule plasmid_concordance)

# Cross-stage input contract (same single-source idea as COMPOSITION): the
# contamination-screen BLAST table. WRITTEN by blast_contigs in
# shared/10_decontam.smk, READ by plasmid_search's supplementary "does the nt hit
# say plasmid?" check. CONTRACT: this file must be BLAST outfmt 6 whose LAST
# column is the subject title (stitle), because the check greps that title for
# the word "plasmid" (v1's blast_contigs already emits exactly this outfmt — keep
# it). The screen runs on the DRAFT assembly (DRAFT_CONTIGS, below).
#
# HYBRID CAVEAT — read PLASMID_BLASTOUT below before touching this: in hybrid mode
# the decontamination screen runs on the ILLUMINA draft while Platon runs on the
# delivered ONT genome. The two assemblies have completely different contig names
# (SPAdes NODE_… vs Flye contig_…), so the plasmid check CANNOT use this file in
# hybrid — it uses PLASMID_BLASTOUT instead.
BLASTOUT = DIR_ASSEMBLY + "/{sample}/contaminants/{sample}_blastout"


# ────────── 4b. Decontamination, read and QC hand-offs (Stage 3/4) ───────────
# Everything below is a PATH (or a tiny parse-time derivation of one). No rule
# logic lives here: the decontamination module (shared/10_decontam.smk), the QC
# module (shared/20_qc.smk), taxonomy (shared/30_taxonomy.smk), the CARD leg in
# shared/50_amr.smk and the report (shared/90_report.smk) all reference these
# names, and the per-mode front ends added in Stage 4 must declare the ones
# marked "Stage-4 contract" as the `output:` of whatever rule they choose. Same
# anti-drift pattern as FINAL_CONTIGS and COMPOSITION: the file is named ONCE.

# ── Where the contamination screen keeps its working files ───────────────────
# This is the same directory COMPOSITION and BLASTOUT already live in; naming it
# once keeps every file below from spelling out the stage path again.
DECONTAM_DIR = DIR_ASSEMBLY + "/{sample}/contaminants"

# Read alignment used ONLY to give BlobTools a coverage track (temp; deleted once
# BlobTools and Qualimap are done with it).
#
# LANDMINE: BlobTools names its coverage file after the BAM's BASENAME, so
# BLOB_COV is DERIVED from DECONTAM_BAM rather than typed out again. Renaming the
# BAM would otherwise silently break blob_json's declared `cov` output. Keeping
# the v1 basename ({sample}_map.bam) in every mode is deliberate for the same
# reason.
DECONTAM_BAM = DECONTAM_DIR + "/{sample}_map.bam"

# BlobTools writes <prefix>.blobDB.json (create) and <prefix>.blob.blobDB.table.txt
# (view), so both rules pass a PREFIX and we derive the real filenames from it.
BLOB_PREFIX       = DECONTAM_DIR + "/blob"
BLOB_JSON         = BLOB_PREFIX + ".blobDB.json"                                  # temp
BLOB_COV          = BLOB_PREFIX + "." + os.path.basename(DECONTAM_BAM) + ".cov"   # temp
BLOB_TABLE_PREFIX = DECONTAM_DIR + "/bestscore"
BLOB_TABLE        = BLOB_TABLE_PREFIX + ".blob.blobDB.table.txt"

# Selector outputs other than COMPOSITION. CONTIG_DECISIONS keeps its v1 name (no
# {sample} prefix): it is the per-contig audit trail named in CLAUDE.md and the
# mobilome spec, and other tooling looks for it by that exact name.
CONTIG_LIST      = DECONTAM_DIR + "/contigs.list"
CONTIG_DECISIONS = DECONTAM_DIR + "/contig_taxonomy_decisions.tsv"

# ── The two per-mode assembly hand-offs (D3, Stage-4 contract) ───────────────
# DRAFT_CONTIGS    — what goes INTO the contamination screen (BLAST + BlobTools +
#                    selector). Stage 4's front end must declare this exact string
#                    as an output.
# DECONTAM_CONTIGS — what the selector WRITES. In illumina/contigs decontamination
#                    is the last assembly step, so this IS the canonical
#                    FINAL_CONTIGS. In nanopore/hybrid it is an intermediate and
#                    FINAL_CONTIGS is produced later by Stage 4 (Medaka in
#                    nanopore; the ONT+Polypolish genome in hybrid).
#
# This split is why select_contigs must not hard-code FINAL_CONTIGS on its output
# side: doing so would create a cycle in nanopore (select → final → Medaka →
# select). Verified acyclic in all four modes with the values below.
# NOTE on the short-read paths: contigs_filt.fasta sits directly under
# 02.assembly/{sample}/, NOT inside the spades/ sub-directory. The assembler
# declares directory(SPADES_DIR) as an output, and Snakemake wipes a directory()
# output before re-running its rule — so a filtered file placed inside spades/
# would be destroyed by any re-run of the assembler.
if MODE == "illumina":
    DRAFT_CONTIGS    = DIR_ASSEMBLY + "/{sample}/contigs_filt.fasta"
    DECONTAM_CONTIGS = FINAL_CONTIGS                       # decontam IS the last assembly step
elif MODE == "hybrid":
    DRAFT_CONTIGS    = DIR_ASSEMBLY + "/{sample}/contigs_filt.fasta"          # the ILLUMINA draft
    DECONTAM_CONTIGS = DECONTAM_DIR + "/contigs_sel.fasta"                    # Stage-4 Snippy reference + QC comparator
elif MODE == "nanopore":
    DRAFT_CONTIGS    = DIR_ASSEMBLY + "/{sample}/fix_start/{sample}_fixed.fasta"
    DECONTAM_CONTIGS = DECONTAM_DIR + "/assembly_decontam.fasta"              # -> Stage-4 Medaka
else:  # contigs
    DRAFT_CONTIGS    = DIR_ASSEMBLY + "/{sample}/contigs_filt.fasta"
    DECONTAM_CONTIGS = FINAL_CONTIGS                       # decontam IS the last assembly step

# ── Which BLAST table the plasmid check greps ────────────────────────────────
# plasmid_search looks each Platon-called contig up in a BLAST table BY CONTIG ID
# (`grep -m 1 -F "$contig"`). That only works if the table was computed over the
# SAME contig set Platon reported on. Whether it was depends on whether anything
# between the contamination screen and the delivered genome can change contig
# names — so we gate on that STRUCTURAL fact, not on a mode name:
#
#   illumina / contigs — the screen runs on the draft and the selector's output IS
#       FINAL_CONTIGS (a subset with identical names). Reuse BLASTOUT: no cost.
#   nanopore — the screen runs on the pre-Medaka fix_start assembly, but Platon
#       runs on the post-Medaka consensus. Nothing guarantees Medaka preserves
#       headers.
#   hybrid  — the screen runs on the ILLUMINA draft (SPAdes NODE_… names) while
#       Platon runs on the delivered ONT genome (Flye contig_… names). These can
#       NEVER match.
#
# In the latter two, reusing BLASTOUT would make `grep` match nothing and EVERY
# plasmid would come back "not verified by BLAST search" — silently, with no error.
# So both long-read modes get a second blastn over FINAL_CONTIGS (v1 BacFluxL+ did
# exactly this inside its own plasmid_search; v2 keeps it in 10_decontam, rule
# blast_final_contigs, so the blastn command text lives in one place).
# Consumed by plasmid_search in shared/60_plasmid.smk.
NEEDS_FINAL_BLAST = HAS_LONG_READS
PLASMID_BLASTOUT = (DECONTAM_DIR + "/{sample}_final_blastout") if NEEDS_FINAL_BLAST else BLASTOUT

# ── Read hand-offs (Stage-4 contract) ────────────────────────────────────────
# Gated exactly like the Flye/Medaka block in section 8, so referencing TRIM_R1 in
# nanopore mode raises a clean NameError instead of silently building a path no
# rule will ever produce. Stage 4's front ends must declare these strings as the
# `output:` of their fastp / filtlong / NanoPlot rules.
if HAS_SHORT_READS:
    # fastp-trimmed pairs. Consumed by the assembler (Stage 4), by map_contigs in
    # 10_decontam and by the CARD read-mapping leg in 50_amr. v1 wrote them to one
    # shared path in BOTH short-read modes, so one constant is faithful. Stage 4
    # declares them temp(); Snakemake keeps them until the last consumer is done.
    TRIM_R1 = DIR_READS + "/{sample}/illumina/{sample}_trim_R1.fastq"
    TRIM_R2 = DIR_READS + "/{sample}/illumina/{sample}_trim_R2.fastq"
    # fastp's JSON report — MultiQC input only.
    FASTP_JSON = DIR_READS + "/{sample}/illumina/{sample}_fastp.json"

if HAS_LONG_READS:
    # filtlong-filtered ONT reads: assembled by Flye (Stage 4) and mapped back
    # onto the draft by map_contigs in 10_decontam.
    FILT_LONG = DIR_READS + "/{sample}/ont/{sample}_filt.fastq"
    # NanoPlot read-QC directories, before and after filtering — MultiQC inputs.
    NANOPLOT_RAW_DIR  = DIR_READS + "/{sample}/ont/raw_qc"
    NANOPLOT_FILT_DIR = DIR_READS + "/{sample}/ont/filt_qc"

# ── Front-end working paths (Stage-4 contract) ───────────────────────────────
# Every file each mode's front end creates on the way from raw input to
# DRAFT_CONTIGS / FINAL_CONTIGS. They are named HERE, not inside the front-end
# modules, for two reasons: (a) shared code consumes some of them —
# shared/15_replicons.smk reads FLYE_INFO and DNAAPLER_SUMMARY, and
# shared/40_annotation.smk reads BAKTA_REPLICONS; and (b) the illumina and hybrid
# front ends duplicate the same SPAdes rules, so a single definition means the
# two copies can drift on flags but never on paths.
#
# Gated exactly like the read hand-offs above: referencing FLYE_DIR in illumina
# mode raises a clean NameError instead of silently building a path no rule will
# ever produce.

if HAS_SHORT_READS:
    # fastp's HTML report — the human-readable twin of FASTP_JSON (which is the
    # one MultiQC reads). Written by trim_adapters; nothing downstream consumes it.
    FASTP_HTML = DIR_READS + "/{sample}/illumina/{sample}_fastp.html"

    # PhiX control genome + its Bowtie2 index. No {sample} in these paths: the
    # download and the index happen ONCE per run and are shared by every sample.
    # All of them are declared temp() by the front end.
    PHIX_DIR        = DIR_READS + "/phix"
    PHIX_FASTA      = PHIX_DIR + "/phix.fna.gz"
    PHIX_BT2_PREFIX = PHIX_DIR + "/phix"

    # SPAdes working directory and its raw output.
    #
    # DELIBERATELY NOT derived from DRAFT_CONTIGS (an earlier draft used
    # os.path.dirname(DRAFT_CONTIGS) as an anti-drift trick). That only worked
    # because DRAFT_CONTIGS then lived INSIDE this directory — and that nesting is
    # a data-loss bug: illumina_assembly declares directory(SPADES_DIR) as an
    # output, and Snakemake removes a directory() output before re-running its
    # rule, so any re-run of the assembler would silently delete DRAFT_CONTIGS,
    # the hand-off the whole contamination screen keys on.
    #
    # DRAFT_CONTIGS therefore sits one level up, beside the assembler directory
    # rather than within it (the same place contigs mode already puts it), and
    # these two are spelled out independently.
    SPADES_DIR     = DIR_ASSEMBLY + "/{sample}/spades"
    SPADES_CONTIGS = SPADES_DIR + "/contigs.fasta"

if HAS_LONG_READS:
    # Flye. assembly.fasta and assembly_info.txt are Flye's own names; the ignore
    # list is ours, built by an awk pass in the same rule (contigs Flye did NOT
    # call circular, which dnaapler must not rotate).
    FLYE_DIR         = DIR_ASSEMBLY + "/{sample}/flye"
    FLYE_CONTIGS     = FLYE_DIR + "/assembly.fasta"
    FLYE_INFO        = FLYE_DIR + "/assembly_info.txt"       # -> shared/15_replicons.smk (topology)
    FLYE_IGNORE_LIST = FLYE_DIR + "/ignore_list.txt"

    # dnaapler (replicon reorientation). DNAAPLER_REORIENTED is dnaapler's own
    # output; DNAAPLER_FIXED is that file with headers trimmed to one token and
    # sequences linearised, and in NANOPORE mode it is also DRAFT_CONTIGS (see the
    # assert below). DNAAPLER_SUMMARY is new in v2 as a declared output: v1 wrote
    # it and never used it, but its Gene_Reoriented / Coverage /
    # Identity_Percentage columns are what shared/15_replicons.smk turns into the
    # Bakta replicon type.
    DNAAPLER_DIR        = DIR_ASSEMBLY + "/{sample}/fix_start"
    DNAAPLER_REORIENTED = DNAAPLER_DIR + "/{sample}_reoriented.fasta"
    DNAAPLER_FIXED      = DNAAPLER_DIR + "/{sample}_fixed.fasta"
    DNAAPLER_SUMMARY    = DNAAPLER_DIR + "/{sample}_all_reorientation_summary.tsv"

    # Medaka (ONT consensus polishing). consensus.fasta is Medaka's own name.
    MEDAKA_DIR       = DIR_ASSEMBLY + "/{sample}/medaka"
    MEDAKA_CONSENSUS = MEDAKA_DIR + "/consensus.fasta"

    # The Bakta --replicons table and its audit trail, both written by
    # build_replicons in shared/15_replicons.smk. They live next to the assembly
    # because they describe it.
    BAKTA_REPLICONS       = DIR_ASSEMBLY + "/{sample}/{sample}_replicons.tsv"
    BAKTA_REPLICONS_AUDIT = DIR_ASSEMBLY + "/{sample}/{sample}_replicons_audit.tsv"
    REPLICONS_SCRIPT      = os.path.join(WORKFLOW_DIR, "scripts", "build_bakta_replicons.py")

    # The nanopore screen runs on the reoriented assembly, so those two names must
    # be the SAME file. Asserted rather than assumed: a future edit to either line
    # would otherwise silently split the DAG into two parallel chains.
    if MODE == "nanopore":
        assert DNAAPLER_FIXED == DRAFT_CONTIGS, (
            "nanopore: DNAAPLER_FIXED and DRAFT_CONTIGS must be the same path "
            f"({DNAAPLER_FIXED!r} vs {DRAFT_CONTIGS!r})"
        )

    # The remaining long-read wiring constants — MEDAKA_INPUT, FINALIZE_SOURCE and
    # POLISH_INPUT — are defined at the END of section 8, because two of them
    # depend on USE_MEDAKA, which is only resolved there.

if IS_HYBRID:
    # The decontaminated Illumina pairs: the reads that mapped as proper pairs to
    # the CLEAN Illumina assembly. Kept (not temp) exactly as in v1 — they have
    # standalone value, and deleting them would force SPAdes and the whole screen
    # to re-run just to re-polish with Polypolish.
    SEL_R1 = DIR_READS + "/{sample}/illumina/{sample}_sel_R1.fastq"
    SEL_R2 = DIR_READS + "/{sample}/illumina/{sample}_sel_R2.fastq"
    # Polypolish working directory (its files are declared individually as temp();
    # the directory itself is deliberately NOT a declared output — see the house
    # rule in rules/hybrid/50_polish.smk).
    POLYPOLISH_DIR = DIR_ASSEMBLY + "/{sample}/polypolish"
    # Snippy star comparison of every ONT stage against the Illumina assembly.
    SNPS_DIR     = DIR_ASSEMBLY + "/{sample}/snps"
    SNPS_SUMMARY = SNPS_DIR + "/SNPs_summary.txt"

# Bakta --replicons wiring. Mode-INDEPENDENT (rule annotation in
# shared/40_annotation.smk is one rule shared by all four modes), but it has to be
# defined AFTER BAKTA_REPLICONS above, not next to FINAL_CONTIGS, or the name
# would not exist yet in the long-read modes.
#
# In illumina/contigs mode this is an EMPTY LIST, which Snakemake renders as an
# empty string in the shell — so the annotation rule's `[ -s "{input.replicons}" ]`
# test is simply false and no --replicons flag is added. One static shell text,
# valid in all four modes, no duplicated rule body and no DAG edge where there is
# no producer.
BAKTA_REPLICON_INPUT = [BAKTA_REPLICONS] if HAS_LONG_READS else []

# ── Assembly QC + taxonomy paths (shared/20_qc.smk, shared/30_taxonomy.smk) ───
# Everything genome-QC-ish lives under one 02.assembly/{sample}/eval/ parent.
QC_GENOMES_DIR   = DIR_ASSEMBLY + "/{sample}/eval/genomes"                    # temp staging dir (see QC_GENOMES)
QC_GENOME_TABLE  = DIR_ASSEMBLY + "/{sample}/eval/{sample}_qc_genomes.tsv"    # kept: which genome is which
QUAST_DIR        = DIR_ASSEMBLY + "/{sample}/eval/quast"
CHECKM_DIR       = DIR_ASSEMBLY + "/{sample}/eval/checkm"
# NOTE the {sample}_ prefix — a deliberate rename from v1's bare checkm_stats.tsv.
# The MultiQC staging loop used to recover the sample with `basename $checkm_dir`;
# under the D1 layout that basename is now the literal "checkm", so the sample name
# has to be carried by the FILE name instead.
CHECKM_STATS     = CHECKM_DIR + "/{sample}_checkm_stats.tsv"
CHECKM_LINEAGE   = CHECKM_DIR + "/lineage.ms"                                 # name chosen by CheckM itself
QUALIMAP_DIR     = DIR_ASSEMBLY + "/{sample}/eval/qualimap"
GTDBTK_DIR       = DIR_TAXONOMY + "/{sample}"

# ── Which genomes get QC'd and classified, and what each one is called ───────
# Three of the four modes deliver ONE genome per sample. Hybrid delivers the
# ONT+Polypolish genome but ALSO keeps the decontaminated Illumina draft, and v1
# BacFluxL+ ran CheckM and GTDB-Tk over both so the two could be compared. This
# list is the single source for: the staged FASTA names, the QUAST assembly
# labels, the CheckM bin ids, the GTDB-Tk bin ids, the MultiQC relabel keys and
# the primary/comparator role column. In v1 the "_illumina"/"_ont" suffixes were
# spelled out by hand in three places hundreds of lines apart.
#
#   suffix — appended to the sample name to build the staged FASTA name, which is
#            what CheckM / GTDB-Tk / QUAST then use as the bin id
#   path   — the {sample}-templated assembly to stage
#   role   — "primary" (the delivered genome) or "comparator" (kept for contrast)
#   label  — technology tag used in the MultiQC report ("" = no tag needed)
QcGenome = namedtuple("QcGenome", "suffix path role label")

if IS_HYBRID:
    QC_GENOMES = [
        QcGenome("_illumina", DECONTAM_CONTIGS, "comparator", "Illumina"),
        QcGenome("_ont",      FINAL_CONTIGS,    "primary",    "ONT"),
    ]
else:
    QC_GENOMES = [QcGenome("", FINAL_CONTIGS, "primary", "")]


def qc_genome_fastas(wildcards):
    # Input function for stage_qc_genomes: this sample's 1 (or, in hybrid, 2)
    # assemblies. An input FUNCTION is needed because the number of files varies
    # by mode, which a plain templated string cannot express.
    return [genome.path.format(sample=wildcards.sample) for genome in QC_GENOMES]


def qc_stage_commands(wildcards):
    # One `cp` line per genome: source assembly -> staged FASTA named after its
    # bin id. Generated at parse time so the literal commands show up in the dry
    # run and in the log, instead of a loop over two hidden bash arrays.
    dest = QC_GENOMES_DIR.format(sample=wildcards.sample)
    lines = []
    for genome in QC_GENOMES:
        src = genome.path.format(sample=wildcards.sample)
        lines.append(f"cp {src} {dest}/{wildcards.sample}{genome.suffix}.fasta")
    return "\n".join(lines)


def qc_genome_table_text(wildcards):
    # The complete text of QC_GENOME_TABLE (header + one row per genome), built
    # from the same QC_GENOMES list. The rule drops it into a quoted heredoc, so
    # the tabs and newlines below are the literal file content — no shell quoting
    # or printf format strings to get wrong. Answers, per sample and in writing,
    # "which of these two hybrid rows is the delivered genome?".
    rows = ["bin_id\trole\ttechnology\tsource_path"]
    for genome in QC_GENOMES:
        bin_id = f"{wildcards.sample}{genome.suffix}"
        technology = genome.label or "NA"
        source = genome.path.format(sample=wildcards.sample)
        rows.append(f"{bin_id}\t{genome.role}\t{technology}\t{source}")
    return "\n".join(rows)


def _relabel_awk(kind):
    # awk body that rewrites a CheckM / GTDB-Tk bin id into its MultiQC report
    # label. `kind` is "completeness" or "taxonomy". Simple modes emit one match
    # on the bare sample name; hybrid emits one per technology, using the SAME
    # suffixes stage_qc_genomes copied the FASTAs under — one source of truth for
    # both halves, which is what v1 lacked.
    #
    # The two match rules deliberately do NOT use `next`, and the report rule
    # appends a bare `{ print }` after them, so a row matching neither key still
    # passes through unchanged (v1 behaviour).
    #
    # BRACES: the text below is a params VALUE, and Snakemake formats only the
    # shell TEMPLATE — substituted values are never re-scanned for placeholders.
    # So these awk braces are SINGLE, the opposite of the {{ }} doubling required
    # for awk written directly in a shell: block. Getting this backwards produces
    # silently corrupted output rather than an error.
    lines = []
    for genome in QC_GENOMES:
        label = (kind + " " + genome.label).strip()
        key = f'sample "{genome.suffix}"' if genome.suffix else "sample"
        lines.append(f'$1 == {key} {{ $1 = "{label} | " sample }}')
    return "\n".join(lines)


CHECKM_RELABEL_AWK = _relabel_awk("completeness")
GTDBTK_RELABEL_AWK = _relabel_awk("taxonomy")

# ── CARD read-mapping leg (shared/50_amr.smk, short-read modes only) ─────────
# The tarball is a SIBLING of the extracted-database directory, not a file inside
# it: the house rule (see 40_annotation.smk) is never to nest one declared output
# inside another rule's directory() output. Both are temp() — the database is only
# needed while BBMap runs. card_link itself is deliberately NOT resolved here (see
# section 7); the download rule reads it, so a config without it still parses in
# the two modes that never use it.
CARD_TARBALL = DIR_AMR + "/card.tar.bz2"
CARD_DB_DIR  = DIR_AMR + "/card_db"

# dbCAN reference database directory + sentinel are defined in section 7, once
# the download link has been resolved — the version folder name is derived from
# the link so the two can never disagree (see DBCAN_DB_ID there).


# ───────────────────────── 5a. Shared constants ─────────────────────────────
# ABRicate is run once per database in this list; AMR_summary reads the same
# list. Defined ONCE here (v1 hard-coded it in both rules, in all four repos)
# and both rules are driven by expand(..., db=DATABASES).
DATABASES = ["argannot", "card", "ecoh", "ecoli_vf", "megares", "ncbi", "resfinder", "vfdb"]
DATABASE_PATTERN = "|".join(DATABASES)

# Sample names become filenames, locus tags, wildcards and report labels. These
# characters break paths or make wildcard matching ambiguous, so a sample whose
# name contains any of them is rejected. NOTE: underscore is intentionally
# ALLOWED — it is common in real sample names and inside the {sample} token of
# {sample}_R1 — matching the hybrid v1 behaviour (the other three v1 workflows
# used to reject it; that stricter rule is dropped in v2).
BAD_CHARS = set("*#@%^/! ?&:;|<>")

# Small awk programs shared by several front-end rules. Kept here as raw strings
# so the fragile one-liners are written once and referenced as `awk {NAME:q}`
# (the :q lets Snakemake quote them safely for the shell).
#   FASTA_LIN_CMD   — turn a wrapped multi-line FASTA into one line per record.
#   FASTA_SEL_CMD   — keep SPAdes contigs with coverage >=2.0 and length >=500
#                     (fields 6 and 4 of the "_"-split SPAdes header).
#   FASTA_HEAD_CMD  — trim a FASTA header down to its first whitespace token.
#   IGNORE_LIST_CMD — from a Flye assembly_info.txt, list contigs NOT flagged
#                     circular (column 4 != "Y"); dnaapler skips reorienting them.
FASTA_LIN_CMD   = r"""{if(NR==1) {printf "%s\n", $0} else {if(/^>/) {printf "\n%s\n", $0} else {printf $0}}}"""
FASTA_SEL_CMD   = r"""{if(/^>/ && $6>=2.0 && $4>=500) {printf "%s\n", $0; getline; print}}"""
FASTA_HEAD_CMD  = r"""BEGIN{FS=" "} /^>/{print $1; next} {print}"""
IGNORE_LIST_CMD = r"""BEGIN{FS="[[:space:]]+"} NR>1 && $4!="Y" {print $1}"""


# ───────────────────────────── 5b. Helpers ──────────────────────────────────
# All of these were byte-identical (or trivially different) across the four v1
# Snakefiles. Deduplicated here; behaviour preserved exactly.

def bakta_locus_tag(sample):
    # Bakta needs a short, safe locus-tag prefix. Derive it from the sample name
    # by keeping alphanumerics plus _ . - and replacing anything else with "_",
    # then truncating to 24 characters. Raise if nothing usable survives — this
    # runs when the annotation rule builds params.locus_tag, so the error points
    # at the offending sample. The user's visible sample name is unchanged.
    tag = "".join(char if char.isalnum() or char in "_.-" else "_" for char in str(sample))
    tag = tag[:24]
    if not tag:
        raise ValueError(f"Unable to derive a valid Bakta locus tag from sample name '{sample}'.")
    return tag


def _as_decontam_text(value):
    # The decontamination selector script expects plain text on the command line.
    # YAML may hand us None, a scalar, or a list — normalise all three:
    #   None            -> ""            (nothing set)
    #   list/tuple/set  -> "A;B;C"       (semicolon-joined, blanks dropped)
    #   anything else   -> its stripped string form
    # NOTE: there is deliberately no dict branch. Three decontamination fields
    # (include_genera_by_sample, exclude_genera_file, sample_overrides) are FILE
    # PATHS that the selector opens on disk, so they must arrive as scalar paths.
    # A YAML mapping ({}) would stringify to "{}" and the selector would try to
    # open a file literally named "{}". The config schema keeps those three as
    # null/scalar paths precisely to avoid that.
    if value is None:
        return ""
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(item).strip() for item in value if str(item).strip())
    return str(value).strip()


def _config_bool(value, default=False):
    # Interpret boolean-like YAML values predictably, so "yes"/"true"/"on"/1 all
    # mean True and "no"/"false"/"off"/0 all mean False. Anything unrecognised is
    # a config mistake and raises rather than being silently coerced.
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "on"}:
        return True
    if text in {"false", "0", "no", "off"}:
        return False
    raise ValueError(f"Invalid boolean config value: {value!r}")


def _decontam_settings(default_mode):
    # Collect every decontamination choice into one dict passed verbatim to the
    # selector script, so all samples follow the same filtering policy.
    #
    # Input:  config["parameters"]["decontamination"] (the v2 schema) OR, for
    #         backwards compatibility, an old top-level parameters["genus"].
    # Output: a normalised dict with a validated mode, five text fields, and a
    #         lowercased "true"/"false" string for discard_no_hit (it is passed
    #         as a CLI argument, hence a string not a bool).
    #
    # This block is mode-INDEPENDENT (decontamination policy is the same for all
    # four modes), so it runs unconditionally.
    params = config.get("parameters", {}) or {}
    settings = dict(params.get("decontamination") or {})

    # Legacy back-compat: a pre-v2 config expressed the target as a single
    # "genus" key. Translate it so those configs keep working, even though
    # "genus" is undocumented in the v2 schema.
    if not settings:
        legacy_genus = params.get("genus", None)
        if legacy_genus is not None and str(legacy_genus).strip():
            settings = {"mode": "include", "include_genera": legacy_genus}
        elif "genus" in params:
            settings = {"mode": "auto"}
        else:
            settings = {"mode": default_mode}

    mode = str(settings.get("mode", default_mode) or default_mode).strip().lower()
    if mode not in {"auto", "include", "exclude", "off"}:
        raise ValueError("Invalid decontamination.mode. Expected one of: auto, include, exclude, off.")
    settings["mode"] = mode

    # include_genera / exclude_genera are inline genus lists; the *_by_sample,
    # *_file and sample_overrides fields are file paths — all normalised to text.
    settings["include_genera"] = _as_decontam_text(settings.get("include_genera", ""))
    settings["include_genera_by_sample"] = _as_decontam_text(settings.get("include_genera_by_sample", ""))
    settings["exclude_genera"] = _as_decontam_text(settings.get("exclude_genera", ""))
    settings["exclude_genera_file"] = _as_decontam_text(settings.get("exclude_genera_file", ""))
    settings["sample_overrides"] = _as_decontam_text(settings.get("sample_overrides", ""))
    settings["discard_no_hit"] = str(_config_bool(settings.get("discard_no_hit"), False)).lower()
    return settings


def validate_extensions(label, extensions, allowed_suffixes):
    # Validate one class of discovered input files. Every sample in a batch must
    # share ONE file extension, otherwise the {sample}.{extn} wildcard paths
    # become ambiguous. Lifted from the hybrid v1 validator (the cleanest of the
    # four) and generalised with an allowed_suffixes tuple so the same helper
    # checks FASTQ inputs (fastq/fq/fastq.gz/fq.gz) and the contigs-mode FASTA
    # inputs (fasta/fa/fna).
    #   Input:  a label for messages, the list of extensions glob_wildcards found,
    #           and the tuple of acceptable suffixes.
    #   Output: the single validated extension string (or a clean exit on error).
    if not extensions:
        sys.stderr.write(f"No suitable {label} input files found.\n")
        sys.exit(0)
    if len(set(extensions)) != 1:
        sys.stderr.write(f"More than one {label} file extension detected:\n\t")
        sys.stderr.write("\n\t".join(sorted(set(extensions))))
        sys.stderr.write("\n")
        sys.exit(0)
    extension = extensions[0]
    if not extension.endswith(allowed_suffixes):
        sys.stderr.write(f"{label.capitalize()} file format '{extension}' not recognized.\n")
        sys.exit(0)
    return extension


def medaka_disabled(value):
    # True only when the user explicitly turned Medaka OFF. Crucially, None does
    # NOT count as disabled — a missing/empty medaka_model means "let Medaka
    # infer the model", not "skip Medaka". This is why USE_MEDAKA is defined as
    # `not medaka_disabled(...)` below. Only used in the long-read modes.
    if isinstance(value, bool):
        return value is False
    if value is None:
        return False
    return str(value).strip().lower() in {"false", "0", "no", "off", "none"}


def has_explicit_medaka_model(value):
    # True only when value is a concrete Medaka model name. The aliases for
    # "infer automatically" (auto/true/yes/on/empty) and the aliases for "skip"
    # (false/no/off/none) all return False. Only used in the long-read modes.
    if isinstance(value, bool) or value is None:
        return False
    return str(value).strip().lower() not in {
        "", "true", "1", "yes", "on", "auto", "false", "0", "no", "off", "none"
    }


# ─────────────── 6. Decontamination policy (mode-independent) ────────────────
# Resolve once and print, so the log shows the active filtering policy up front.
DECONTAMINATION = _decontam_settings("auto")
print(
    "Decontamination mode: "
    f"{DECONTAMINATION['mode']} (discard_no_hit={DECONTAMINATION['discard_no_hit']})."
)


# ─────────────── 7. Reference-database download links (parse-time) ───────────
# CheckV and dbCAN links are resolved here because both databases are used by
# every mode; card_link and phix_link are resolved here too, but only checked when
# the mode actually has short reads. All of them fail EARLY with a message naming
# the config key, rather than as a bare KeyError deep inside a rule.
_links = config.get("links") or {}

# CheckV: the link is OPTIONAL. When absent/empty, CheckV downloads its own
# default database and we keep a sensible folder id. When present it must be a
# .tar.gz, and the id is derived from the archive name (strip .tar then .gz).
# (Unified from four divergent v1 forms: the nanopore workflow used to sys.exit
# on a missing link — that user-hostile hard-exit is dropped here.)
CHECKV_LINK = str(_links.get("checkv_link") or "").strip()
CHECKV_DB_ID = "checkv-db-v1.5"

# Three ways to get a CheckV database, in strict order of precedence:
#   1. directories.checkv_db  — a copy you already hold. Nothing is downloaded.
#   2. links.checkv_link      — fetch this .tar.gz, unpack it, build the diamond DB.
#   3. neither                — let CheckV fetch its own default database.
# Case 1 wins outright, and we say so out loud when a link was ALSO set, because a
# silently ignored config key is exactly the kind of thing that wastes an afternoon.
if CHECKVDB:
    # Fail here, at parse time, rather than an hour into a run when viral_quality
    # finally opens the directory. Checking for the actual reference FASTA (not
    # just that the directory exists) also catches the common mistake of pointing
    # at the versioned folder's PARENT's parent, or at a half-unpacked archive.
    _checkv_reps = glob.glob(os.path.join(CHECKVDB, "*", "genome_db", "checkv_reps.faa"))
    if not os.path.isdir(CHECKVDB):
        sys.exit(
            f"[BacFlux] directories.checkv_db points at '{CHECKVDB}', which is not a "
            "directory. Give the PARENT directory that holds the versioned database "
            "folder, e.g. /path/to/checkv/ containing checkv-db-v1.5/."
        )
    if len(_checkv_reps) != 1:
        sys.exit(
            f"[BacFlux] directories.checkv_db is '{CHECKVDB}', but that directory holds "
            f"{len(_checkv_reps)} CheckV database(s) (looking for */genome_db/checkv_reps.faa). "
            "Exactly one is required. Give the PARENT directory that holds a single "
            "versioned database folder, e.g. /path/to/checkv/ containing checkv-db-v1.5/."
        )
    CHECKV_DB_ID = os.path.basename(os.path.dirname(os.path.dirname(_checkv_reps[0])))
    print(f"Using the local CheckV database at '{CHECKVDB}' (db_id='{CHECKV_DB_ID}'). Nothing will be downloaded.")
    if CHECKV_LINK:
        print(
            "  NOTE: links.checkv_link is also set and is being IGNORED — "
            "directories.checkv_db takes precedence."
        )
elif not CHECKV_LINK:
    print(
        "The link to the CheckV database is not specified (or empty). "
        "CheckV will download the database automatically."
    )
else:
    _checkv_name = os.path.basename(urlparse(CHECKV_LINK).path)
    if not _checkv_name.endswith(".tar.gz"):
        raise ValueError(f"Invalid checkv_link: expected a .tar.gz archive, got '{_checkv_name}'.")
    CHECKV_DB_ID = os.path.splitext(os.path.splitext(_checkv_name)[0])[0]
    print(f"Using CheckV database from link: '{CHECKV_LINK}' (db_id='{CHECKV_DB_ID}').")

# dbCAN: the link is REQUIRED (every mode annotates CAZymes) and must be a
# .tar.gz. The matching checksum URL is derived by swapping the suffix. Resolved
# once here (v1 did this inline per rule) so a bad link fails fast at parse time.
DBCAN_LINK = str(_links.get("dbcan_link") or "").strip()
if not DBCAN_LINK:
    raise ValueError("Missing required 'links.dbcan_link' value in the config file.")
if not DBCAN_LINK.endswith(".tar.gz"):
    raise ValueError(
        f"Invalid dbcan_link: expected a .tar.gz archive, got '{os.path.basename(DBCAN_LINK)}'."
    )
DBCAN_SHA_URL = DBCAN_LINK.replace(".tar.gz", ".sha256")

# Name the extracted-DB directory after the archive itself (strip .tar then .gz),
# so the folder version always matches whatever dbcan_link points at. The same
# approach is used for CheckV above. (v1 hard-coded "dbcan_db_v5.1.2", which
# would silently mismatch a link pointing at any other dbCAN version.) The
# sentinel is an empty marker file written after checksum verification so
# Snakemake does not repeat the slow download/verify step.
DBCAN_DB_ID    = os.path.splitext(os.path.splitext(os.path.basename(urlparse(DBCAN_LINK).path))[0])[0]
DBCAN_DB_DIR   = DIR_ANNOTATION + "/dbcan/" + DBCAN_DB_ID
DBCAN_SENTINEL = DBCAN_DB_DIR + "/.verified.sha256"

# CARD: needed only by the read-based AMR leg, which exists only where the mode
# produces short reads (D7). Resolved here so a missing key reports the config key
# by name at parse time instead of surfacing as a bare KeyError from inside the
# download rule. Modes without short reads never look at it.
CARD_LINK = str(_links.get("card_link") or "").strip()
if HAS_SHORT_READS and not CARD_LINK:
    sys.exit(
        f"[BacFlux] mode={MODE} runs the read-based CARD AMR leg, which requires "
        "'links.card_link' to be set in the config."
    )

# PhiX: the small bacteriophage genome Illumina spikes into essentially every
# lane as a sequencing control. Its reads are real sequence from a different
# organism, so the short-read front ends map them out before assembly. Resolved
# here — rather than inline in the download rule as v1 did — because BOTH the
# illumina and hybrid front ends need it and both need the same check, so there is
# one copy of the validation instead of two. Modes without short reads never look
# at it. (This deliberately revises the earlier note in this section that said
# phix_link would stay unresolved here.)
PHIX_LINK = str(_links.get("phix_link") or "").strip()
if HAS_SHORT_READS and not PHIX_LINK:
    sys.exit(
        f"[BacFlux] mode={MODE} removes PhiX spike-in reads before assembly, which "
        "requires 'links.phix_link' to be set in the config."
    )


# ──────────── 8. Long-read assembly / polishing options (gated) ──────────────
# Flye and Medaka only exist in the long-read modes, so their parameter
# resolution — and its several prints — is gated on HAS_LONG_READS. The
# short-read and contigs modes leave the three exposed names as None.
#
# v1→v2 change: these parameters used to be flat (parameters.flye_input_mode,
# parameters.medaka_model). In v2 they are namespaced per mode, so we read
# config["parameters"][MODE] (i.e. parameters.nanopore.* or parameters.hybrid.*).
if HAS_LONG_READS:
    _lr_params = (config.get("parameters", {}) or {}).get(MODE, {}) or {}
    _medaka_model_raw = _lr_params.get("medaka_model")             # None => infer, NOT skip
    _flye_input_mode_raw = _lr_params.get("flye_input_mode", "auto")

    # Decide whether Medaka runs at all. Only an explicit false-like value skips
    # it; None/auto/empty mean "run Medaka and infer the model".
    USE_MEDAKA = not medaka_disabled(_medaka_model_raw)
    MEDAKA_MODEL = None
    FLYE_INPUT_MODE = "--nano-hq"

    # An empty Flye mode is treated as auto (safest default for ONT assemblies).
    if _flye_input_mode_raw is None or str(_flye_input_mode_raw).strip() == "":
        _flye_input_mode_raw = "auto"
    FLYE_MODE_KEY = str(_flye_input_mode_raw).strip().lower()
    FLYE_MODE_MAP = {"nano-raw": "--nano-raw", "nano-hq": "--nano-hq"}
    if FLYE_MODE_KEY not in {"auto", *FLYE_MODE_MAP.keys()}:
        raise ValueError(
            "Invalid 'flye_input_mode' parameter. Expected one of: auto, nano-raw, nano-hq."
        )

    # Announce the Medaka decision and capture an explicit model name if given.
    if not USE_MEDAKA:
        print("Medaka disabled via config. Skipping ONT consensus polishing.")
    elif has_explicit_medaka_model(_medaka_model_raw):
        # str() first: a bare number in YAML (medaka_model: 507) loads as int, and
        # calling .strip() on it directly would crash. has_explicit_medaka_model
        # already stringifies, so a numeric-looking model tag is a valid name here.
        MEDAKA_MODEL = str(_medaka_model_raw).strip()
        print(f"The 'medaka_model' parameter is specified with value: '{MEDAKA_MODEL}'.")
    else:
        print("Medaka model will be inferred automatically from FASTQ headers.")

    # In auto mode only an explicit *fast* Medaka model switches Flye to
    # --nano-raw (fast basecalling => noisier reads); otherwise assume
    # high-quality ONT reads and use --nano-hq.
    if FLYE_MODE_KEY == "auto":
        if MEDAKA_MODEL is not None and "fast" in MEDAKA_MODEL.lower():
            FLYE_INPUT_MODE = "--nano-raw"
            print("Flye input mode set automatically to '--nano-raw' based on the explicit Medaka fast model.")
        else:
            FLYE_INPUT_MODE = "--nano-hq"
            print("Flye input mode set automatically to '--nano-hq'.")
    else:
        FLYE_INPUT_MODE = FLYE_MODE_MAP[FLYE_MODE_KEY]
        print(f"Flye input mode overridden via config with value: '{FLYE_INPUT_MODE}'.")

    # ── The two long-read modes differ ONLY in these wiring constants ─────────
    # Expressing the difference here, once, is what lets the Flye / dnaapler /
    # Medaka rules be byte-identical in rules/nanopore/ and rules/hybrid/. They
    # live in THIS section, not with the other front-end paths in section 4b,
    # because two of them depend on USE_MEDAKA, which is resolved just above.
    #
    # All three are PARSE-TIME constants, so the DAG is fixed before the run
    # starts and `snakemake -n` shows the real chain. v1 used input FUNCTIONS
    # (final_contigs(wc), polishing_input_contigs(wc)) that re-read the config
    # while the DAG was being built.
    if MODE == "nanopore":
        # Medaka polishes AFTER the contamination screen (D3, v1 behaviour).
        MEDAKA_INPUT = DECONTAM_CONTIGS
        # What finalize_contigs copies to FINAL_CONTIGS.
        FINALIZE_SOURCE = MEDAKA_CONSENSUS if USE_MEDAKA else DECONTAM_CONTIGS
    else:  # hybrid
        # The ONT leg is never decontaminated: its reads were already filtered
        # against the decontaminated Illumina reads by filtlong. So Medaka
        # polishes the PRE-screen reoriented assembly.
        MEDAKA_INPUT = DNAAPLER_FIXED
        # What Polypolish corrects, and therefore what becomes FINAL_CONTIGS.
        POLISH_INPUT = MEDAKA_CONSENSUS if USE_MEDAKA else DNAAPLER_FIXED
else:
    # Short-read-only and contigs modes never touch Flye/Medaka.
    USE_MEDAKA = None
    MEDAKA_MODEL = None
    FLYE_INPUT_MODE = None


# ──────────────── 9. Per-mode sample discovery (gated) ───────────────────────
# Only the ACTIVE mode's branch runs: only it touches the filesystem (globs its
# input directory) and only it can exit on "no input" / "bad sample name". The
# other three modes' discovery never fires. This is the gated pattern the
# migration plan requires.
FASTQ_SUFFIXES = ("fastq", "fq", "fastq.gz", "fq.gz")
FASTA_SUFFIXES = ("fasta", "fa", "fna")


def _require_input_dir(config_key):
    # Resolve one input directory from config["input"][config_key], and fail
    # clearly if it is unset or does not exist. Called once per required input
    # directory for the active mode (hybrid needs two).
    inputs = config.get("input", {}) or {}
    value = str(inputs.get(config_key) or "").strip()
    if not value:
        sys.exit(f"[BacFlux] mode={MODE} requires 'input.{config_key}' to be set in the config.")
    if not os.path.isdir(value):
        sys.exit(f"[BacFlux] input.{config_key} is not an existing directory: {value!r}")
    return value


def _check_sample_names(samples):
    # Reject sample names that would break wildcards or paths, and echo each
    # accepted sample so the user can confirm the batch before it runs.
    for sample in sorted(samples):
        if any(char in sample for char in BAD_CHARS):
            sys.stderr.write(f"Sample name '{sample}' contains unsupported characters.\n")
            sys.exit(0)
        print(f"Sample {sample} will be processed.")


def _require_r2_mates(samples, illumina_dir, short_extn):
    # Illumina input is paired; discovery keys off R1, so verify each sample also
    # has its R2 mate before the run rather than failing deep in a rule later.
    # (Deliberate v1->v2 change: BacFlux v1 illumina did NOT check this and would
    # fail later in map_phix; the hybrid v1 workflow did check it. v2 applies the
    # friendlier early check to both short-read modes.)
    for sample in samples:
        r2_path = os.path.join(illumina_dir, f"{sample}_R2.{short_extn}")
        if not os.path.exists(r2_path):
            sys.stderr.write(f"Missing Illumina mate pair for sample '{sample}': {r2_path}\n")
            sys.exit(0)


if MODE == "illumina":
    # Illumina PE: discover samples from R1 files; R2 path is rebuilt per sample.
    ILLUMINA_DIR = _require_input_dir("illumina_dir")
    _samples, _exts = glob_wildcards(os.path.join(ILLUMINA_DIR, "{sample}_R1.{extn}"))
    SHORT_EXTN = validate_extensions("illumina", _exts, FASTQ_SUFFIXES)
    SAMPLES = sorted(set(_samples))
    _check_sample_names(SAMPLES)
    _require_r2_mates(SAMPLES, ILLUMINA_DIR, SHORT_EXTN)
    R1 = "{sample}_R1." + SHORT_EXTN
    R2 = "{sample}_R2." + SHORT_EXTN

elif MODE == "nanopore":
    # ONT long reads: discover from files named {sample}_ont.{extn}.
    NANOPORE_DIR = _require_input_dir("nanopore_dir")
    _samples, _exts = glob_wildcards(os.path.join(NANOPORE_DIR, "{sample}_ont.{extn}"))
    LONG_EXTN = validate_extensions("nanopore", _exts, FASTQ_SUFFIXES)
    SAMPLES = sorted(set(_samples))
    _check_sample_names(SAMPLES)
    ONT = "{sample}_ont." + LONG_EXTN

elif MODE == "hybrid":
    # Hybrid: Illumina and ONT are discovered independently (from two separate
    # directories in v2) and then intersected — a sample must have BOTH.
    ILLUMINA_DIR = _require_input_dir("illumina_dir")
    NANOPORE_DIR = _require_input_dir("nanopore_dir")
    _short_samples, _short_exts = glob_wildcards(os.path.join(ILLUMINA_DIR, "{sample}_R1.{extn}"))
    _long_samples, _long_exts = glob_wildcards(os.path.join(NANOPORE_DIR, "{sample}_ont.{extn}"))
    if not _short_samples:
        sys.exit(f"[BacFlux] no Illumina R1 files found in {ILLUMINA_DIR}.")
    if not _long_samples:
        sys.exit(f"[BacFlux] no ONT files found in {NANOPORE_DIR}.")
    SHORT_EXTN = validate_extensions("illumina", _short_exts, FASTQ_SUFFIXES)
    LONG_EXTN = validate_extensions("nanopore", _long_exts, FASTQ_SUFFIXES)

    # Make a mismatched pairing explicit instead of letting a missing mate become
    # a confusing rule error later.
    _missing_long = sorted(set(_short_samples) - set(_long_samples))
    _missing_short = sorted(set(_long_samples) - set(_short_samples))
    if _missing_long or _missing_short:
        if _missing_long:
            sys.stderr.write("Missing ONT files for the following samples:\n")
            sys.stderr.write("\n".join(f"  - {s}" for s in _missing_long) + "\n")
        if _missing_short:
            sys.stderr.write("Missing Illumina R1 files for the following samples:\n")
            sys.stderr.write("\n".join(f"  - {s}" for s in _missing_short) + "\n")
        sys.exit(0)

    SAMPLES = sorted(set(_short_samples) & set(_long_samples))
    _check_sample_names(SAMPLES)
    _require_r2_mates(SAMPLES, ILLUMINA_DIR, SHORT_EXTN)
    R1 = "{sample}_R1." + SHORT_EXTN
    R2 = "{sample}_R2." + SHORT_EXTN
    ONT = "{sample}_ont." + LONG_EXTN

elif MODE == "contigs":
    # Pre-assembled genomes: discover from {sample}.{extn}. glob_wildcards is
    # GREEDY on {sample}, so the LAST dot splits off the extension: a dotted name
    # like "my.genome.fasta" is read as sample="my.genome", extn="fasta" and is
    # ACCEPTED (the dot stays in the sample name). This matches v1 FastaFlux,
    # which used the same pattern. Dotted sample names flow safely through paths
    # and locus tags; add an explicit check here if you ever need to forbid them.
    CONTIGS_DIR = _require_input_dir("contigs_dir")
    _samples, _exts = glob_wildcards(os.path.join(CONTIGS_DIR, "{sample}.{extn}"))
    CONTIGS_EXTN = validate_extensions("contigs", _exts, FASTA_SUFFIXES)
    SAMPLES = sorted(set(_samples))
    _check_sample_names(SAMPLES)
    CONTIGS = "{sample}." + CONTIGS_EXTN


# ──────────────── 9b. Global {sample} wildcard constraint ────────────────────
# Pin the {sample} wildcard to the exact set of discovered sample names. Without
# this, {sample} defaults to ".+" (matches anything, even "/"), which lets a
# fixed path component collide with a per-sample rule — e.g. the shared
# 04.annotation/antismash/databases directory could be "claimed" by the
# per-sample antismash rule with sample="databases", and likewise the dbCAN DB
# folder. Restricting {sample} to real sample names removes that whole class of
# accidental rule collisions across every module at once, in one place.
# Fallback to ".+" when no samples are discovered (e.g. the bare skeleton at
# parse time) so the directive is always valid; no {sample} job runs in that case.
if SAMPLES:
    _SAMPLE_CONSTRAINT = "|".join(re.escape(s) for s in SAMPLES)
else:
    _SAMPLE_CONSTRAINT = r".+"


wildcard_constraints:
    sample=_SAMPLE_CONSTRAINT


# ─────────────────────────── 10. Rule-all targets ───────────────────────────
# all_targets() is what `rule all` requests. It is the one place 00_common has
# to know downstream output names (acceptable, and matches MetaFlux).
#
# STAGE-1 SAFETY: during the initial migration only this file exists — the rule
# modules that would PRODUCE these targets have not been written yet. Asking
# Snakemake to build a file that no rule produces makes even a dry run fail. So
# until the active mode's front-end rule directory exists on disk, all_targets()
# returns an empty list and `snakemake -n` reports "Nothing to be done" — the
# clean parse the Stage-1 gate checks for. Once the modules land, the full list
# below is returned. (Front-end modules are implemented last in the migration
# plan, so their presence is a good proxy for "the pipeline is fully wired".)

def _rule_modules_present():
    # True once rules/<MODE>/ contains at least one .smk file. glob on a missing
    # directory returns [], so this is False for the bare Stage-1 skeleton.
    return bool(glob.glob(os.path.join(WORKFLOW_DIR, "rules", MODE, "*.smk")))


def _frontend_targets_for(mode):
    # Front-end leaf targets. Every mode guarantees the canonical finished
    # assembly (D2, FINAL_CONTIGS); on top of that, ANYTHING a front end produces
    # that no other rule consumes has to be asked for BY NAME here, or Snakemake
    # will simply never build it.
    #
    # We branch on the capability FLAGS rather than on the `mode` argument (D7),
    # so the intent reads as "this mode has long reads" instead of enumerating
    # mode names. The parameter is kept only so the call site in all_targets()
    # does not change.
    targets = list(expand(FINAL_CONTIGS, sample=SAMPLES))

    if HAS_SHORT_READS:
        # Listed explicitly for clarity even though multiqc already depends on it.
        targets += expand(FASTP_JSON, sample=SAMPLES)

    if HAS_LONG_READS:
        # The two NanoPlot report directories (multiqc pulls these in too).
        targets += expand(NANOPLOT_RAW_DIR, sample=SAMPLES)
        targets += expand(NANOPLOT_FILT_DIR, sample=SAMPLES)
        # ORPHAN without this line: the replicon AUDIT table is terminal — nothing
        # reads it — so it would never be produced. (BAKTA_REPLICONS itself is not
        # an orphan: rule annotation consumes it. FLYE_INFO and DNAAPLER_SUMMARY
        # are not orphans either, now that build_replicons reads them.)
        targets += expand(BAKTA_REPLICONS_AUDIT, sample=SAMPLES)

    if IS_HYBRID:
        # ORPHAN without this line: the Snippy comparison of the ONT stages
        # against the Illumina assembly is a terminal report.
        targets += expand(SNPS_SUMMARY, sample=SAMPLES)

    return targets


def _downstream_targets():
    # The shared 03.* – 09.* deliverables. Identical across modes (that is the
    # whole point of the D1 layout), with two conditional legs. These match the
    # rules that Stages 2–3 add under rules/shared/.
    targets = [
        # 02.assembly/eval — assembly + genome QC leaves. Listed explicitly
        # because v1's `rule all` listed them in all four workflows, and because
        # CheckM is NOT on the path to anything else: GTDB-Tk now reads the staged
        # genomes directly (not CheckM's output), so without this line CheckM
        # would only be reached indirectly, via MultiQC.
        *expand(QC_GENOME_TABLE, sample=SAMPLES),
        *expand(QUAST_DIR,       sample=SAMPLES),
        *expand(CHECKM_STATS,    sample=SAMPLES),
        # 03.taxonomy — GTDB-Tk (one output directory per sample)
        *expand(DIR_TAXONOMY + "/{sample}", sample=SAMPLES),
        # 04.annotation — Bakta, eggNOG, antiSMASH, dbCAN (shared DB + per sample)
        *expand(DIR_ANNOTATION + "/bakta/{sample}", sample=SAMPLES),
        *expand(DIR_ANNOTATION + "/eggnog/{sample}", sample=SAMPLES),
        *expand(DIR_ANNOTATION + "/antismash/{sample}", sample=SAMPLES),
        DBCAN_DB_DIR,
        *expand(DIR_ANNOTATION + "/dbcan/{sample}", sample=SAMPLES),
        # 05.amr — ABRicate: one TSV per database plus the per-sample summary
        *expand(DIR_AMR + "/abricate/{sample}/{db}.tsv", sample=SAMPLES, db=DATABASES),
        *expand(DIR_AMR + "/abricate/{sample}/AMR_summary.txt", sample=SAMPLES),
        # 06.plasmids — terminal plasmid deliverable depends on the phage caller:
        #   geNomad opt-in → the Platon+geNomad concordance TSV (D9; pulls in both
        #                    Platon via PLATON_DIR and geNomad via GENOMAD_DIR).
        #   default        → Platon's verified_plasmids.txt (Platon-only; geNomad
        #                    never runs, so there is no concordance to build).
        # Whichever is the terminal target, plasmid_search always runs.
        *expand(PLASMID_CONCORDANCE if PHAGE_CALLER == "genomad" else PLATON_VERIFIED, sample=SAMPLES),
        # 07.phages — CheckV (unconditional; it pulls in whichever caller ran:
        # geNomad's virus FASTA when opted in, else VirSorter2's).
        *expand(DIR_PHAGES + "/checkv/{sample}", sample=SAMPLES),
        # 09.report — MultiQC
        DIR_REPORT + "/multiqc_report.html",
    ]
    # Mapping QC exists wherever real reads were mapped. Contigs mode maps the
    # contigs against themselves purely to give BlobTools a coverage track, and a
    # Qualimap report on a self-alignment says nothing — so it has none.
    if HAS_READS:
        targets += [*expand(QUALIMAP_DIR, sample=SAMPLES)]
    # CARD read-mapping AMR leg exists only where short reads are produced.
    if HAS_SHORT_READS:
        targets += [
            *expand(DIR_AMR + "/mapping/{sample}/{sample}_covstats.tsv", sample=SAMPLES),
            *expand(DIR_AMR + "/mapping/{sample}/{sample}_AMR_legend.tsv", sample=SAMPLES),
        ]
    # Mobilome module (08.mobilome) is opt-in and default OFF. The module itself
    # is not written yet, so asking for its outputs would abort the DAG with a
    # MissingInputException naming a directory rather than the config key that
    # caused it. Fail with an actionable message instead, and only request the
    # targets once the module actually exists on disk.
    if _config_bool((config.get("mobilome", {}) or {}).get("run"), False):
        if not glob.glob(os.path.join(WORKFLOW_DIR, "rules", "shared", "80_mobilome.smk")):
            sys.exit(
                "[BacFlux] config.mobilome.run is true, but the mobilome module "
                "(workflow/rules/shared/80_mobilome.smk) is not implemented yet. "
                "Set mobilome.run: false."
            )
        targets += [*expand(DIR_MOBILOME + "/{sample}", sample=SAMPLES)]
    return targets


def all_targets():
    if not _rule_modules_present():
        return []
    return _frontend_targets_for(MODE) + _downstream_targets()