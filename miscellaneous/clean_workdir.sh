#!/usr/bin/env bash
# clean_workdir.sh
#
# Purpose:
#   Remove bulky temporary folders and downloaded databases that are not required to keep the final results of the workflow.
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

need_dir() {
  local d="$1"
  [[ -d "$d" ]] || { echo "ERROR: expected directory not found: $d" >&2; exit 1; }
}

# ensure we are in the right place (security check)
need_dir "02.assembly"

# All targets
TARGETS=(
  # 02.assembly: SPAdes intermediates and temporary folders
  "02.assembly/*/K*"
  "02.assembly/*/misc"
  "02.assembly/*/pipeline_state"
  "02.assembly/*/tmp"

  # 03.post-processing: completeness evaluation intermediates
  "03.post-processing/completeness_evaluation/*/bins"
  "03.post-processing/completeness_evaluation/*/storage"
  "03.post-processing/completeness_evaluation/*/lineage.ms"

  # 04.taxonomy: GTDB-Tk subfolders
  "04.taxonomy/*/align"
  "04.taxonomy/*/identify"

  # 05.annotation: downloaded databases
  "05.annotation/dbcan/dbcan_db_v5.1.2"
  "05.annotation/antismash/databases"

  # 08.phages: databases and temporary folders
  "08.phages/*_db"
  "08.phages/checkv/*/tmp"
  "08.phages/virsorter/*/iter-*"
)

echo "== Cleanup utility =="
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
