# Shared environment for every EBI mobilome-annotation-pipeline (MAP) run.
#
# Why this file exists: MAP is a Nextflow pipeline that expects a container
# engine. Singularity 2.x predates the SIF container format that all current
# biocontainers use, so it cannot run them at all. If your system Singularity
# is that old, install Apptainer (the renamed, modern continuation of
# Singularity) into a conda env of its own, together with nextflow; this file
# expects that env at $EBI_MAP_ROOT/envs/nf. Apptainer works "rootless"
# (without root) only if the kernel has unprivileged user namespaces enabled.
# The conda package ships a `singularity` -> `apptainer` symlink, which is what
# lets Nextflow's stock `singularity` profile work unmodified.
#
# Source this before any nextflow command:  source nf_env.sh

# EBI_MAP_ROOT is your MAP work directory: it holds the envs/nf conda env and
# the caches below. Set it before sourcing this file; there is no default.
# This is an explicit test, not ${EBI_MAP_ROOT:?...}, because when this file
# is sourced from an interactive shell, bash prints that error but still runs
# the next lines, which would point PATH, JAVA_HOME and the caches at the top
# of the file system.
if [ -z "${EBI_MAP_ROOT:-}" ]; then
    echo "set EBI_MAP_ROOT to your MAP work directory" >&2
    return 1 2>/dev/null || exit 1   # return if sourced, exit if run as a script
fi

# Put our nextflow + apptainer ahead of any older system singularity.
export PATH="${EBI_MAP_ROOT}/envs/nf/bin:$PATH"
export JAVA_HOME="${EBI_MAP_ROOT}/envs/nf"

# Keep every large, regenerable artefact under EBI_MAP_ROOT, and put that
# folder on a disk with plenty of free space: $HOME and /tmp are often small.
export NXF_SINGULARITY_CACHEDIR="${EBI_MAP_ROOT}/cache/singularity"
export APPTAINER_CACHEDIR="${EBI_MAP_ROOT}/cache/apptainer"
export APPTAINER_TMPDIR="${EBI_MAP_ROOT}/cache/apptainer_tmp"
export NXF_HOME="${EBI_MAP_ROOT}/cache/nxf_home"
export NXF_TEMP="${EBI_MAP_ROOT}/cache/nxf_tmp"

mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" \
         "$NXF_HOME" "$NXF_TEMP"
