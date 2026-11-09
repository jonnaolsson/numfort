#!/usr/bin/env bash

# Run the complete numfort CI matrix in the Fortran CI container.
#
# - Mount repository root and build-ci directory into the container.
# - Execute run-matrix.sh with oneAPI environment initialized.
# - Run up to 8 configurations in parallel by default.
# - Forward any command-line arguments to run-matrix.sh.

# ==============================================================================
# Execution environment
# ==============================================================================
set -euo pipefail
export LC_ALL=C

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../.." && pwd)

mkdir -p "${repo_root}/build-ci"

args=(--jobs 8)
if (( $# > 0 )); then
    args=("${@}")
fi

exec docker run --rm --ulimit core=0 --user "$(id -u):$(id -g)" \
    -v "${repo_root}:/src/numfort:ro,z" \
    -v "${repo_root}/build-ci:/build:z" \
    fortran-ci:ubuntu26.04 with-oneapi \
    bash /src/numfort/ci/ubuntu26.04/run-matrix.sh "${args[@]}"
