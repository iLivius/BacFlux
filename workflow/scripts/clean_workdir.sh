#!/usr/bin/env bash
# clean_workdir.sh
#
# Remove bulky intermediate folders and downloaded workflow-local databases from
# completed bacterial BioFlux outputs. The script is dry-run by default and
# deletes only when called with --run.
#
# WHEN TO RUN THIS: after you have checked the results and are ready to archive
# or share the output directory. It is deliberately aggressive -- it strips
# intermediate reads, assembler working directories and every downloaded
# database, keeping only results, reports and audit trails. Rerunning the
# workflow on a cleaned directory therefore recomputes a lot (assembly and
# polishing included). That is the intended trade-off, not an oversight.
#
# Raw input reads are never touched: they live in input.illumina_dir /
# input.nanopore_dir and are read in place, never copied into the output tree.
#
# Supported layouts:
#   BacFlux v2 (all four modes: illumina, nanopore, hybrid, contigs -- one
#   directory numbering scheme, so no per-mode branching is needed)
#   plus the retired v1 family kept for any output dirs still around from
#   before the v2 cutover: BacFlux, FastaFlux, BacFluxL, BacFluxL+
#
# Usage:
#   clean_workdir.sh [OUTPUT_DIR]
#   clean_workdir.sh --target OUTPUT_DIR
#   clean_workdir.sh --run [OUTPUT_DIR]
#   clean_workdir.sh --run --target OUTPUT_DIR
#   clean_workdir.sh --run --include-snakemake OUTPUT_DIR
#
# If OUTPUT_DIR is omitted, the current directory is used.
#
# --include-snakemake (opt-in, OFF by default)
#   Also delete a .snakemake/ working directory if one is found INSIDE the
#   output directory. Whether there is one depends entirely on where the
#   workflow was launched from: Snakemake creates .snakemake/ in the CURRENT
#   WORKING DIRECTORY, not in output_dir. Launch from the repo (the documented
#   way: `snakemake --sdm conda --configfile config/config.yaml`) and it lands
#   next to the Snakefile, out of this script's reach. Launch from inside the
#   output folder and it lands there, where this flag can reach it.
#
#   It is opt-in because it is categorically heavier than everything else here.
#   The rest of this script deletes regenerable intermediates; .snakemake/ holds
#   the per-rule CONDA ENVIRONMENTS (tens of GB -- routinely the largest single
#   item in a finished project), the provenance metadata that drives Snakemake's
#   rerun decisions, and the run logs. Deleting it means the next run on that
#   directory rebuilds every environment from scratch, which needs network
#   access and time.
#
#   Use it when archiving or handing over a finished project. Do not use it on
#   a directory you still intend to rerun soon.

set -euo pipefail

DO_RUN=0
TARGET_DIR=""
INCLUDE_SNAKEMAKE=0

usage() {
  # Print the header comment block: everything from the title line down to the
  # last line before `set -euo pipefail`. Derived rather than a fixed line
  # range, so --help cannot silently truncate when the header is edited.
  sed -n '2,/^set -euo pipefail/p' "$0" | sed '$d'
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run)
      DO_RUN=1
      shift
      ;;
    --dry-run)
      DO_RUN=0
      shift
      ;;
    --include-snakemake)
      INCLUDE_SNAKEMAKE=1
      shift
      ;;
    -t|--target|--output-dir)
      [[ $# -ge 2 ]] || die "$1 requires a directory argument."
      TARGET_DIR="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -*)
      die "Unknown option: $1"
      ;;
    *)
      [[ -z "$TARGET_DIR" ]] || die "Only one output directory can be specified."
      TARGET_DIR="$1"
      shift
      ;;
  esac
done

[[ $# -eq 0 ]] || die "Unexpected extra argument(s): $*"

TARGET_DIR="${TARGET_DIR:-.}"
[[ -d "$TARGET_DIR" ]] || die "Target is not a directory: $TARGET_DIR"

TARGET_DIR="$(cd "$TARGET_DIR" && pwd -P)"
[[ "$TARGET_DIR" != "/" ]] || die "Refusing to operate on /"
[[ "$TARGET_DIR" != "$HOME" ]] || die "Refusing to operate on HOME: $TARGET_DIR"

cd "$TARGET_DIR"

has_dirs() {
  local d
  for d in "$@"; do
    [[ -d "$d" ]] || return 1
  done
  return 0
}

detect_workflow() {
  # v2 unified all four modes (illumina/nanopore/hybrid/contigs) onto one
  # directory numbering scheme, so a single marker set covers every mode --
  # unlike v1, where each repo in the family renumbered its own stages (see
  # the elif chain below). 08.mobilome is deliberately NOT part of this check:
  # the module defaults to off, so its absence does not mean "not v2 output".
  if has_dirs 01.reads 02.assembly 04.annotation 09.report; then
    echo "BacFlux-v2"
  elif has_dirs 02.Illumina_assembly 03.post-processing 04.ONT_assembly 10.annotation 13.phages 14.report; then
    echo "BacFluxL+"
  elif has_dirs 01.pre-processing 02.assembly 03.post-processing 04.taxonomy 05.annotation 08.phages 09.report; then
    echo "BacFluxL"
  elif has_dirs 01.pre-processing 02.assembly 03.post-processing 04.taxonomy 05.annotation 08.phages; then
    echo "BacFlux"
  elif has_dirs 01.pre-processing 02.post-processing 03.taxonomy 04.annotation 07.phages 08.report; then
    echo "FastaFlux"
  else
    die "Target does not look like a supported BacFlux/FastaFlux/BacFluxL/BacFluxL+ output directory."
  fi
}

WORKFLOW="$(detect_workflow)"
TARGETS=()

case "$WORKFLOW" in
  BacFlux-v2)
    TARGETS+=(
      # ── Intermediate read files (01.reads) ────────────────────────────────
      # EVERY sequence file under 01.reads is an intermediate: trimmed
      # (fastp), decontaminated (bowtie2-selected pairs) or length/quality
      # filtered (filtlong) reads. The RAW reads are never copied into the
      # output tree -- they are read in place from input.illumina_dir /
      # input.nanopore_dir -- so nothing here is anyone's only copy.
      #
      # Deliberately a blanket sequence-file glob rather than a list of the
      # specific names (*_sel_R1, *_filt, ...): it needs no maintenance when a
      # mode gains a new read-processing step, and it cannot silently miss one.
      # It is also, in practice, the single biggest win in a hybrid run --
      # decontaminated Illumina pairs alone routinely outweigh every other
      # remaining file combined.
      #
      # This DOES retrigger assembly and polishing if the workflow is rerun
      # later, which is accepted: this script is an archive/share step run once
      # the results are known good, not a between-runs tidy.
      #
      # NOT touched here: fastp JSON/HTML and the NanoPlot raw_qc/filt_qc
      # directories in the same folders. They are QC REPORTS (and MultiQC
      # inputs), not sequence data.
      "01.reads/*/*/*.fastq"
      "01.reads/*/*/*.fq"
      "01.reads/*/*/*.fastq.gz"
      "01.reads/*/*/*.fq.gz"

      # PhiX spike-in reference plus its bowtie2 index, downloaded per project
      # (links.phix_link) for the short-read decontamination step. A database,
      # not a result -- re-fetched automatically if the workflow is rerun.
      "01.reads/phix"

      # ── Per-sample assembly scaffolding (02.assembly) ──────────────────────
      # Every mode writes its raw assembler's working directory alongside the
      # final, decontaminated, polished genome (contigs_final.fasta). Once that
      # file exists the working directories below are redundant -- they were
      # only ever inputs to later steps, never something a colleague reads
      # directly.
      "02.assembly/*/spades"                              # illumina/hybrid: SPAdes K-mer dirs, misc, tmp, assembly graphs
      "02.assembly/*/flye"                                 # nanopore/hybrid: Flye 00-assembly .. 40-polishing
      "02.assembly/*/medaka"                                # nanopore/hybrid: long-read polishing BAM + consensus (already copied into contigs_final.fasta)
      "02.assembly/*/snps/*_snps_dir"                       # hybrid: per-stage Snippy BAM/VCF/HTML dirs -- keeps snps/SNPs_summary.txt
      "02.assembly/*/fix_start/*_fixed.fasta"               # hybrid: dnaapler circularisation intermediates --
      "02.assembly/*/fix_start/*_fixed.fasta.fai"           # keeps *_all_reorientation_summary.tsv, the actual
      "02.assembly/*/fix_start/*_fixed.fasta.map-ont.mmi"   # audit record of what was reoriented and why
      "02.assembly/*/fix_start/*_reoriented.fasta"
      "02.assembly/*/fix_start/*_MMseqs2_output.txt"
      "02.assembly/*/fix_start/dnaapler_*.log"
      "02.assembly/*/fix_start/logs"
      # minimap2/samtools indexes built ad hoc against the decontaminated
      # assembly (nanopore/hybrid decontam + read-mapping steps). These are
      # NOT declared Snakemake outputs, so nothing invalidates them between
      # reruns automatically -- they are also, concretely, the exact files
      # that let Medaka silently reuse a stale alignment index after the
      # decontamination step is rerun, polishing the OLD contig set while
      # reporting success (observed in practice, not hypothetical). Routine
      # cleanup removes the trap, not just the disk usage.
      "02.assembly/*/contaminants/*.fai"
      "02.assembly/*/contaminants/*.mmi"

      # CheckM's own scratch (bins/storage/lineage.ms) -- the kept result is
      # *_checkm_stats.tsv.
      "02.assembly/*/eval/checkm/bins"
      "02.assembly/*/eval/checkm/storage"
      "02.assembly/*/eval/checkm/*.ms"

      # GTDB-Tk (03.taxonomy/{sample}) has NO cleanup target here. Unlike v1,
      # where align/identify were separate scratch dirs distinct from wherever
      # the final summary landed, v2's classify_wf (GTDB-Tk >=2.7, per
      # workflow/rules/shared/30_taxonomy.smk's own rule comment) nests the
      # real result INSIDE classify/: classify/gtdbtk.bac120.summary.tsv is
      # the file MultiQC reads, reached from the top level only via a
      # convenience symlink. A "03.taxonomy/*/classify" target here would
      # delete that real result and leave the symlink dangling -- confirmed
      # the hard way, and recoverable only by rerunning taxonomic_assignment
      # for every affected sample. Do not re-add a classify/ target without
      # checking where classify_wf's own output actually lives for the GTDB-Tk
      # version pinned in workflow/envs/gtdbtk.yaml.

      # antiSMASH's bundled reference database, fetched into the output tree
      # per project rather than pointed at a shared path. Reliably the single
      # largest item here (multi-GB).
      "04.annotation/antismash/databases"

      # dbCAN database mirror (symlinks back to directories.dbcan_db in the
      # config -- the real, shared copy lives there, not here).
      "04.annotation/dbcan/dbcan_db_v5.1.2"

      # eggNOG-mapper's temp working directory (emptied by the tool itself,
      # left behind as an empty shell).
      "04.annotation/eggnog/*_tmp"

      # Phage-calling stage (07.phages): shared databases, whichever caller is
      # in use. checkv_db and genomad_db mirror directories.checkv_db /
      # directories.genomad_db; the glob also covers virsorter_db if
      # VirSorter2 is the active caller. The two virsorter patterns below are
      # carried over from the v1 target list and have not been re-verified
      # against a v2 VirSorter2 run.
      "07.phages/*_db"
      "07.phages/virsorter/*/iter-*"                        # unverified against a v2 VirSorter2 run; harmless no-op via nullglob if absent
      "07.phages/virsorter/*/log"

      # Mobilome module (08.mobilome): shared reference-model/database
      # directories, fetched once per project rather than per sample.
      # NOTE -- conjscan_models and icescan_models are third-party,
      # CC BY-NC-SA-licensed model packages fetched from GitHub's API, which
      # has a 60-request/hour UNAUTHENTICATED limit shared by the whole host.
      # Deleting them here means the next BacFlux run on this machine pays
      # that cost again, and a host that has recently run several mobilome
      # jobs can genuinely exhaust the hourly budget. If another mobilome run
      # is expected soon after this cleanup, consider keeping these two and
      # pointing mobilome.icescan.dir at the retained copy instead.
      "08.mobilome/conjscan_models"
      "08.mobilome/icescan_models"
      "08.mobilome/iceberg_db"
      "08.mobilome/tncentral_db"
      # ISOSDB: the IS nucleotide database behind the read-mapping copy-number
      # leg, so it only appears in illumina/hybrid runs (nanopore has no
      # read-mapping step to attach it to). Refetched from its configured URL.
      "08.mobilome/isosdb_db"
      # Per-sample mobilome tool directories (conjscan/, icescan/, isescan/
      # under 08.mobilome/{sample}/) are deliberately NOT listed: like
      # 04.annotation/bakta/{sample} or 04.annotation/dbcan/{sample}, they are
      # per-tool RESULTS, not scratch, and the same restraint v1 already
      # applied to Bakta/dbCAN output applies here too.
    )
    ;;
  BacFlux)
    TARGETS+=(
      "02.assembly/*/K*"
      "02.assembly/*/misc"
      "02.assembly/*/pipeline_state"
      "02.assembly/*/tmp"
      "03.post-processing/completeness_evaluation/*/bins"
      "03.post-processing/completeness_evaluation/*/storage"
      "03.post-processing/completeness_evaluation/*/*.ms"
      "04.taxonomy/*/align"
      "04.taxonomy/*/identify"
      "05.annotation/antismash/databases"
      "05.annotation/dbcan/dbcan_db_v5.1.2"
      "08.phages/*_db"
      "08.phages/checkv/*/tmp"
      "08.phages/virsorter/*/iter-*"
      "08.phages/virsorter/*/log"
    )
    ;;
  FastaFlux)
    TARGETS+=(
      "02.post-processing/completeness_evaluation/*/bins"
      "02.post-processing/completeness_evaluation/*/storage"
      "02.post-processing/completeness_evaluation/*/*.ms"
      "03.taxonomy/*/align"
      "03.taxonomy/*/identify"
      "04.annotation/antismash/databases"
      "04.annotation/dbcan/dbcan_db_v5.1.2"
      "07.phages/*_db"
      "07.phages/checkv/*/tmp"
      "07.phages/virsorter/*/iter-*"
      "07.phages/virsorter/*/log"
    )
    ;;
  BacFluxL)
    TARGETS+=(
      "01.pre-processing/*.fastq"
      "01.pre-processing/*.fq"
      "01.pre-processing/*.fastq.gz"
      "01.pre-processing/*.fq.gz"
      "02.assembly/*/00-assembly"
      "02.assembly/*/10-consensus"
      "02.assembly/*/20-repeat"
      "02.assembly/*/30-contigger"
      "02.assembly/*/40-polishing"
      "03.post-processing/completeness_evaluation/*/bins"
      "03.post-processing/completeness_evaluation/*/storage"
      "03.post-processing/completeness_evaluation/*/*.ms"
      "03.post-processing/consensus/*/*.bam*"
      "03.post-processing/consensus/*/*.bed"
      "03.post-processing/consensus/*/*.hdf"
      "03.post-processing/*/*.fai"
      "03.post-processing/*/*.mmi"
      "03.post-processing/contaminants/*/*.fai"
      "03.post-processing/contaminants/*/*.mmi"
      "04.taxonomy/*/align"
      "04.taxonomy/*/identify"
      "05.annotation/antismash/databases"
      "05.annotation/dbcan/dbcan_db_v5.1.2"
      "08.phages/*_db"
      "08.phages/checkv/*/tmp"
      "08.phages/virsorter/*/iter-*"
      "08.phages/virsorter/*/log"
    )
    ;;
  BacFluxL+)
    TARGETS+=(
      "02.Illumina_assembly/*/K*"
      "02.Illumina_assembly/*/misc"
      "02.Illumina_assembly/*/pipeline_state"
      "02.Illumina_assembly/*/tmp"
      "03.post-processing/*_sel_R1.fastq"
      "03.post-processing/*_sel_R2.fastq"
      "03.post-processing/*_ont_filt.fastq"
      "04.ONT_assembly/*/00-assembly"
      "04.ONT_assembly/*/10-consensus"
      "04.ONT_assembly/*/20-repeat"
      "04.ONT_assembly/*/30-contigger"
      "04.ONT_assembly/*/40-polishing"
      "05.ONT_consensus/*/*.bam*"
      "05.ONT_consensus/*/*.bed"
      "05.ONT_consensus/*/*.hdf"
      "08.SNPs/*/*_snps_dir"
      "09.taxonomy/*/align"
      "09.taxonomy/*/identify"
      "10.annotation/antismash/databases"
      "10.annotation/dbcan/dbcan_db_v5.1.2"
      "11.AMR/AMR_db"
      "11.AMR/AMR_mapping/*/ref"
      "13.phages/*_db"
      "13.phages/checkv/*/tmp"
      "13.phages/virsorter/*/iter-*"
      "13.phages/virsorter/*/log"
    )
    ;;
esac

shopt -s nullglob dotglob

to_delete=()
for pattern in "${TARGETS[@]}"; do
  matches=( $pattern )
  for match in "${matches[@]}"; do
    # The -e test is required, not belt-and-braces: `nullglob` only drops
    # patterns that CONTAIN a wildcard. A literal path with no glob character
    # (e.g. "01.reads/phix", "08.mobilome/iceberg_db" for an optional database
    # that was never enabled) is not subject to pathname expansion at all, so it
    # passes through as itself even when nothing is there. Without this check the
    # dry run lists paths that do not exist and will not be deleted -- which is
    # the wrong way round for the preview people rely on before running --run.
    [[ -e "$match" || -L "$match" ]] || continue
    to_delete+=( "$match" )
  done
done

# Snakemake drops a zero-byte .snakemake_timestamp marker inside every
# directory() output to track whether that directory is up to date. They are
# bookkeeping, not results: meaningless to anyone opening the archive, and they
# survive the patterns above whenever they sit in a directory this script KEEPS
# (Bakta, isescan, platon, ... ). Collected with `find` rather than a glob
# because they occur at several different depths.
#
# Note the `dotglob` above already makes the "*" patterns match hidden entries,
# and anything inside a deleted directory goes with it -- this step exists only
# for the markers left behind in directories that are kept.
#
# Removing them makes Snakemake treat those directory outputs as out of date on
# a later rerun. That is the same accepted trade-off as the rest of this script.
#
# These markers are unrelated to the .snakemake/ WORKING DIRECTORY handled
# separately below -- same prefix, different thing.
#
# Kept in their own array so the listing below can report them as a single
# count. Printing all of them (86 in a routine two-sample run) would bury the
# entries that actually matter in the preview people read before --run.
timestamp_markers=()
while IFS= read -r marker; do
  timestamp_markers+=( "$marker" )
done < <(find . -name ".snakemake_timestamp" -type f 2>/dev/null | sed 's|^\./||')

# ── The .snakemake/ working directory (--include-snakemake only) ─────────────
# Whether one is here at all depends on where the workflow was launched from,
# NOT on any config setting: Snakemake creates .snakemake/ in the current
# working directory. Launching from the repo (the documented way) puts it next
# to the Snakefile, where this script never looks. Launching from inside the
# output folder puts it here.
#
# Only the top level is checked, deliberately. A `find` for .snakemake at any
# depth could match one belonging to a DIFFERENT project that happens to live
# under this tree, and deleting tens of GB of someone else's conda
# environments is not something this script should ever do by inference.
snakemake_workdir=()
if (( INCLUDE_SNAKEMAKE )) && [[ -d ".snakemake" ]]; then
  snakemake_workdir=( ".snakemake" )
fi

echo "== Bacterial BioFlux output cleanup =="
echo "Workflow: $WORKFLOW"
echo "Target: $TARGET_DIR"
if [[ "$DO_RUN" -eq 1 ]]; then
  echo "Mode: RUN (will delete)"
else
  echo "Mode: DRY-RUN (no deletion)"
fi
echo

if (( ${#to_delete[@]} == 0 && ${#timestamp_markers[@]} == 0 && ${#snakemake_workdir[@]} == 0 )); then
  echo "Nothing to clean; no matching files or directories were found."
  # If the only thing here was a .snakemake/ the user did not opt into, say so
  # rather than reporting a clean directory and leaving tens of GB in place.
  if [[ -d ".snakemake" ]]; then
    echo
    echo "NOTE: a .snakemake/ working directory ($(du -sh .snakemake 2>/dev/null | cut -f1)) is present"
    echo "      but was NOT touched. Re-run with --include-snakemake to remove it."
  fi
  exit 0
fi

echo "Targets:"
for path in "${to_delete[@]}"; do
  echo "  - $path"
done
if (( ${#timestamp_markers[@]} )); then
  echo "  - ${#timestamp_markers[@]} x .snakemake_timestamp (zero-byte Snakemake directory markers)"
fi
if (( ${#snakemake_workdir[@]} )); then
  echo "  - .snakemake/ working directory ($(du -sh .snakemake 2>/dev/null | cut -f1)) -- conda envs, metadata and logs"
fi
echo

# Flag the untouched .snakemake/ whenever there IS one and it was not opted
# into, so its size is visible at exactly the moment someone is deciding what
# to archive -- otherwise it is the largest thing left and the easiest to miss.
if (( ! INCLUDE_SNAKEMAKE )) && [[ -d ".snakemake" ]]; then
  echo "NOTE: .snakemake/ ($(du -sh .snakemake 2>/dev/null | cut -f1)) is present and will NOT be removed."
  echo "      It holds this project's conda environments, provenance metadata and logs."
  echo "      Add --include-snakemake to delete it too (rebuilt from scratch on the next run)."
  echo
fi

if [[ "$DO_RUN" -eq 1 ]]; then
  # All three arrays in one call. Overlap is harmless: a marker inside a
  # directory that is itself being deleted simply no longer exists by the time
  # rm reaches it, and -f makes that a no-op rather than an error.
  rm -rf -- "${to_delete[@]}" "${timestamp_markers[@]}" "${snakemake_workdir[@]}"
  echo "Cleanup finished."
else
  echo "Dry-run only. Re-run with --run to delete."
fi
