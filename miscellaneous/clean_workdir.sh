#!/usr/bin/env bash
# Compatibility wrapper, kept because this is the path older notes, scripts and
# habits point at. The cleanup script itself moved to
# workflow/scripts/clean_workdir.sh in v1.3.0, next to the other workflow
# scripts; this file only forwards to it.
#
# Every argument is passed straight through, so the real script's behaviour is
# unchanged here — including the important part: it is DRY-RUN unless --run is
# given. Read the header of workflow/scripts/clean_workdir.sh (or run this with
# --help) for what it deletes and what it keeps.

set -euo pipefail

# Resolve the sibling repo path from this file's own location rather than the
# caller's working directory, so the wrapper works from anywhere. pwd -P resolves
# symlinks to the real path, so it also works when the repo is reached through one.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
TARGET="${SCRIPT_DIR}/../workflow/scripts/clean_workdir.sh"

# The reminder goes to stderr, not stdout, so piping or capturing this script's
# output still gives only the cleanup listing. `exec` replaces this shell with
# the real script, which keeps its exit status as ours.
echo "NOTE: clean_workdir.sh moved to workflow/scripts/clean_workdir.sh" >&2
exec "$TARGET" "$@"
