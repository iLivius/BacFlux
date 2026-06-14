#!/usr/bin/env bash
# clean_workdir.sh
#
# Remove bulky intermediate folders and downloaded workflow-local databases from
# completed bacterial BioFlux outputs. The script is dry-run by default and
# deletes only when called with --run.
#
# Supported layouts:
#   BacFlux, FastaFlux, BacFluxL, BacFluxL+
#
# Usage:
#   clean_workdir.sh [OUTPUT_DIR]
#   clean_workdir.sh --target OUTPUT_DIR
#   clean_workdir.sh --run [OUTPUT_DIR]
#   clean_workdir.sh --run --target OUTPUT_DIR
#
# If OUTPUT_DIR is omitted, the current directory is used.

set -euo pipefail

DO_RUN=0
TARGET_DIR=""

usage() {
  sed -n '1,25p' "$0"
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
  if has_dirs 02.Illumina_assembly 03.post-processing 04.ONT_assembly 10.annotation 13.phages 14.report; then
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
  if (( ${#matches[@]} )); then
    to_delete+=( "${matches[@]}" )
  fi
done

echo "== Bacterial BioFlux output cleanup =="
echo "Workflow: $WORKFLOW"
echo "Target: $TARGET_DIR"
if [[ "$DO_RUN" -eq 1 ]]; then
  echo "Mode: RUN (will delete)"
else
  echo "Mode: DRY-RUN (no deletion)"
fi
echo

if (( ${#to_delete[@]} == 0 )); then
  echo "Nothing to clean; no matching files or directories were found."
  exit 0
fi

echo "Targets:"
for path in "${to_delete[@]}"; do
  echo "  - $path"
done
echo

if [[ "$DO_RUN" -eq 1 ]]; then
  rm -rf -- "${to_delete[@]}"
  echo "Cleanup finished."
else
  echo "Dry-run only. Re-run with --run to delete."
fi
