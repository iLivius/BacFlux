#!/usr/bin/env bash
# Compatibility wrapper. The cleanup script now lives in workflow/scripts.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
TARGET="${SCRIPT_DIR}/../workflow/scripts/clean_workdir.sh"

echo "NOTE: clean_workdir.sh moved to workflow/scripts/clean_workdir.sh" >&2
exec "$TARGET" "$@"
