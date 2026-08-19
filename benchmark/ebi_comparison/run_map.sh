#!/usr/bin/env bash
#
# Run the EBI mobilome-annotation-pipeline (MAP) on BacFlux benchmark genomes,
# so its ICE/IME calls can be compared against BacFlux's own caller
# (workflow/scripts/mobilome/conjscan_to_ice.py) on identical input.
#
# Usage:  ./run_map.sh <samplesheet.csv> <run_name>
# Example: ./run_map.sh samplesheet_pilot3.csv pilot3
#
# Reads : a MAP samplesheet (sample,assembly,...) pointing at the .fna files in
#         phase7_benchmark/genomes/ - the same FASTAs our own caller was scored on.
# Writes: runs/<run_name>/results/<sample>/prediction/icefinder2lite/<sample>_ices.tsv
#         which is the table the comparison reads.

set -euo pipefail

SAMPLESHEET="${1:?usage: run_map.sh <samplesheet.csv> <run_name>}"
RUN_NAME="${2:?usage: run_map.sh <samplesheet.csv> <run_name>}"

EBI_MAP_ROOT=/media/data/antonielli_dir/BacFlux_v2_validation/ebi_map
source "${EBI_MAP_ROOT}/nf_env.sh"

# Pin Nextflow to 24.10.6.
#
# This is NOT cosmetic. On Nextflow 26.04.6 (the current release, which the conda
# package installs) MAP v5.0.0 dies inside its own ICEfinder2-lite subworkflow:
#
#   ERROR ~ Invalid method invocation `call` with arguments:
#           [[id:pos_test], null, .../pos_test_uniprot_names.tsv] on _closure11 type
#
# That is the `.join(..., remainder: true)` feeding REFINE_BOUNDARIES in
# subworkflows/local/icefinder2lite.nf handing a 3-element tuple to a closure
# written for 6. Newer Nextflow stopped tolerating that arity mismatch. MAP's
# manifest only asks for '!>=24.04.0' and its CI does not pin a version, so this
# is a real forward-compatibility gap in MAP, not a local misconfiguration.
# On 24.10.6 the same commit runs REFINE_BOUNDARIES to completion.
export NXF_VER=24.10.6

RUN_DIR="${EBI_MAP_ROOT}/runs/${RUN_NAME}"
mkdir -p "${RUN_DIR}"

cd "${EBI_MAP_ROOT}/map"

nextflow -log "${RUN_DIR}/nextflow.log" \
    run main.nf \
    -profile singularity \
    --input "${SAMPLESHEET}" \
    --outdir "${RUN_DIR}/results" \
    -work-dir "${RUN_DIR}/work" \
    -c "${EBI_MAP_ROOT}/dbs.config" \
    -c "${EBI_MAP_ROOT}/skip_functional.config" \
    -resume \
    -ansi-log false
