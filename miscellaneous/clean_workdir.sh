#!/usr/bin/env bash
# clean_workdir.sh
#
# Purpose:
#   Remove bulky temporary folders and downloaded databases that are not required to keep the final results of the workflow.
#   Supports both BacFlux and FastaFlux output directories.
#
# Usage:
#   ./clean_workdir.sh            # dry-run (prints what would be removed)
#   ./clean_workdir.sh --run      # actually delete
#
# Notes:
#   Run from the workflow output directory.

set -euo pipefail

DO_RUN=0
if [[ "${1:-}" == "--run" ]]; then
  DO_RUN=1
elif [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  sed -n '1,80p' "$0"
  exit 0
fi

detect_workflow() {
  if [[ -d "02.assembly" && -d "03.post-processing" && -d "04.taxonomy" ]]; then
    echo "BacFlux"
  elif [[ -d "01.pre-processing" && -d "02.post-processing" && -d "03.taxonomy" ]]; then
    echo "FastaFlux"
  else
    echo "ERROR: current directory does not look like a BacFlux or FastaFlux output directory." >&2
    exit 1
  fi
}

WORKFLOW="$(detect_workflow)"

POSTPROC_DIR=""
TAXONOMY_DIR=""
ANNOTATION_DIR=""
PHAGES_DIR=""
TARGETS=()

if [[ "$WORKFLOW" == "BacFlux" ]]; then
  POSTPROC_DIR="03.post-processing"
  TAXONOMY_DIR="04.taxonomy"
  ANNOTATION_DIR="05.annotation"
  PHAGES_DIR="08.phages"

  TARGETS+=(
    # 02.assembly: SPAdes intermediates and temporary folders
    "02.assembly/*/K*"
    "02.assembly/*/misc"
    "02.assembly/*/pipeline_state"
    "02.assembly/*/tmp"
  )
else
  POSTPROC_DIR="02.post-processing"
  TAXONOMY_DIR="03.taxonomy"
  ANNOTATION_DIR="04.annotation"
  PHAGES_DIR="07.phages"
fi

TARGETS+=(
  # Completeness evaluation intermediates
  "${POSTPROC_DIR}/completeness_evaluation/*/bins"
  "${POSTPROC_DIR}/completeness_evaluation/*/storage"
  "${POSTPROC_DIR}/completeness_evaluation/*/lineage.ms"

  # GTDB-Tk subfolders
  "${TAXONOMY_DIR}/*/align"
  "${TAXONOMY_DIR}/*/identify"

  # Downloaded databases
  "${ANNOTATION_DIR}/dbcan/dbcan_db_v5.1.2"
  "${ANNOTATION_DIR}/antismash/databases"

  # Phage databases and temporary folders
  "${PHAGES_DIR}/*_db"
  "${PHAGES_DIR}/checkv/*/tmp"
  "${PHAGES_DIR}/virsorter/*/iter-*"
)

echo "== Cleanup utility =="
echo "Workflow: $WORKFLOW"
if [[ "$DO_RUN" -eq 1 ]]; then
  echo "Mode: RUN (will delete)"
else
  echo "Mode: DRY-RUN (no deletion)"
fi
echo

# Expand targets safely; ignore non-existing matches
shopt -s nullglob dotglob

to_delete=()
for pat in "${TARGETS[@]}"; do
  # If it looks like a glob, expand it; otherwise, only keep it if it exists.
  if [[ "$pat" == *[\*\?\[]* ]]; then
    matches=( $pat )   # glob expansion happens here
    if (( ${#matches[@]} )); then
      to_delete+=( "${matches[@]}" )
    fi
  else
    if [[ -e "$pat" ]]; then
      to_delete+=( "$pat" )
    fi
  fi
done

if (( ${#to_delete[@]} == 0 )); then
  echo "Nothing to clean (no matches)."
  exit 0
fi

echo "Targets:"
for p in "${to_delete[@]}"; do
  echo "  - $p"
done
echo

if [[ "$DO_RUN" -eq 1 ]]; then
  rm -rf -- "${to_delete[@]}"
  echo "Done."
else
  echo "Dry-run only. Re-run with --run to delete."
fi
