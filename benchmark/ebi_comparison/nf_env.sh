# Shared environment for every EBI mobilome-annotation-pipeline (MAP) run.
#
# Why this file exists: MAP is a Nextflow pipeline that expects a container
# engine. This host has Singularity 2.6.1, which predates the SIF container
# format that all current biocontainers use, so it cannot run them at all.
# Apptainer (the renamed, modern continuation of Singularity) was installed
# into a private conda env instead, and it works "rootless" here because the
# kernel has unprivileged user namespaces enabled. The conda package ships a
# `singularity` -> `apptainer` symlink, which is what lets Nextflow's stock
# `singularity` profile work unmodified.
#
# Source this before any nextflow command:  source nf_env.sh

EBI_MAP_ROOT=/media/data/antonielli_dir/BacFlux_v2_validation/ebi_map

# Put our nextflow + apptainer ahead of the system singularity 2.6.1.
export PATH="${EBI_MAP_ROOT}/envs/nf/bin:$PATH"
export JAVA_HOME="${EBI_MAP_ROOT}/envs/nf"

# Keep every large, regenerable artefact on /media/data (2 TB free), not in
# $HOME or /tmp, which are small here.
export NXF_SINGULARITY_CACHEDIR="${EBI_MAP_ROOT}/cache/singularity"
export APPTAINER_CACHEDIR="${EBI_MAP_ROOT}/cache/apptainer"
export APPTAINER_TMPDIR="${EBI_MAP_ROOT}/cache/apptainer_tmp"
export NXF_HOME="${EBI_MAP_ROOT}/cache/nxf_home"
export NXF_TEMP="${EBI_MAP_ROOT}/cache/nxf_tmp"

mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" \
         "$NXF_HOME" "$NXF_TEMP"
