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
from pathlib import Path
from urllib.parse import urlparse

# glob_wildcards is a Snakemake helper; import it explicitly so the reader can
# see where it comes from rather than relying on an injected global.
from snakemake.io import glob_wildcards


# ─────────────────────────── 1. Mode dispatch ───────────────────────────────
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
# the "Genus:percentage" composition table. It is WRITTEN by the decontamination
# selector (select_contigs, in the future shared/10_decontam.smk) and READ by the
# Bakta annotation rule to pick the isolate's most likely genus. Defining the one
# canonical path here means producer and consumer can never drift onto different
# strings (the failure the Stage-2a review flagged). Lives under 02.assembly/ next
# to the contaminant-screening outputs, per the D1 layout.
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
CHECKV_DB_DIR = DIR_PHAGES + "/checkv_db"                       # produced by rule checkv_db

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
# contamination-screen BLAST table. WRITTEN by blast_contigs in the future
# shared/10_decontam.smk (Stage 3), READ here by plasmid_search's supplementary
# "does the nt hit say plasmid?" check. CONTRACT for Stage 3: this file must be
# BLAST outfmt 6 whose LAST column is the subject title (stitle), because the
# check greps that title for the word "plasmid" (v1's blast_contigs already emits
# exactly this outfmt — keep it). Cannot be verified until 10_decontam lands.
BLASTOUT = DIR_ASSEMBLY + "/{sample}/contaminants/{sample}_blastout"

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
# every mode. phix_link and card_link are NOT resolved here — they are only
# needed by the short-read front end and are read inside their rules, gated on
# HAS_SHORT_READS.
_links = config.get("links") or {}

# CheckV: the link is OPTIONAL. When absent/empty, CheckV downloads its own
# default database and we keep a sensible folder id. When present it must be a
# .tar.gz, and the id is derived from the archive name (strip .tar then .gz).
# (Unified from four divergent v1 forms: the nanopore workflow used to sys.exit
# on a missing link — that user-hostile hard-exit is dropped here.)
CHECKV_LINK = str(_links.get("checkv_link") or "").strip()
CHECKV_DB_ID = "checkv-db-v1.5"
if not CHECKV_LINK:
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
    # Front-end leaf targets. Every mode guarantees the canonical decontaminated
    # assembly (D2); each front-end module (added in Stage 4) extends this with
    # its own read-QC and assembly-QC leaves, which the report module also pulls
    # in as its inputs.
    return list(expand(FINAL_CONTIGS, sample=SAMPLES))


def _downstream_targets():
    # The shared 03.* – 09.* deliverables. Identical across modes (that is the
    # whole point of the D1 layout), with two conditional legs. These match the
    # rules that Stages 2–3 add under rules/shared/.
    targets = [
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
    # CARD read-mapping AMR leg exists only where short reads are produced.
    if HAS_SHORT_READS:
        targets += [
            *expand(DIR_AMR + "/mapping/{sample}/{sample}_covstats.tsv", sample=SAMPLES),
            *expand(DIR_AMR + "/mapping/{sample}/{sample}_AMR_legend.tsv", sample=SAMPLES),
        ]
    # Mobilome module (08.mobilome) is opt-in and default OFF.
    if _config_bool((config.get("mobilome", {}) or {}).get("run"), False):
        targets += [*expand(DIR_MOBILOME + "/{sample}", sample=SAMPLES)]
    return targets


def all_targets():
    if not _rule_modules_present():
        return []
    return _frontend_targets_for(MODE) + _downstream_targets()