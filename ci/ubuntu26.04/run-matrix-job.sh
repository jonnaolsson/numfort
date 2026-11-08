#!/usr/bin/env bash

# Execute one matrix configuration as a CTest worker.
#
# - Called by the matrix scheduler; use run-matrix.sh as the entry point.
# - Receive the run directory, job name, compiler, and build options as
#   arguments.
# - Delegate configuration, builds, installation, and tests to run-tests.sh.
# - Archive previous CTest logs and preserve current logs.
# - Write a separate result row and log for the parent matrix summary.
# - Return the delegated runner's exit status to CTest.

# ==============================================================================
# Execution environment
# ==============================================================================
set -euo pipefail
export LC_ALL=C
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

# ==============================================================================
# Worker arguments
# ==============================================================================
if (( $# != 12 )); then
    printf 'Internal worker: invoke through run-matrix.sh.\n' >&2
    exit 2
fi
run_dir=${1}
name=${2}
compiler=${3}
config=${4}
shared=${5}
tests=${6}
examples=${7}
openmp=${8}
ipo=${9}
blas95=${10}
mkl=${11}
native=${12}

# ==============================================================================
# Job paths and previous logs
# ==============================================================================
log=${run_dir}/${name}.log
toolchain=gnu-${compiler}
if [[ ${compiler} == ifx ]]; then
    toolchain=ifx
fi
job_dir=${BUILD_ROOT:-/build}/${toolchain}/${config}/shared-${shared}-tests-${tests}-examples-${examples}-openmp-${openmp}-ipo-${ipo}-blas95-${blas95}-mkl-${mkl}-native-${native}
test_log=${job_dir}/numfort/Testing/Temporary/LastTest.log
if [[ -f ${test_log} ]]; then
    mv -- "${test_log}" "${run_dir}/${name}-numfort-previous-LastTest.log"
fi

# ==============================================================================
# Configuration execution
# ==============================================================================
code=0
BUILD_SHARED_LIBS=${shared} NUMFORT_BUILD_TESTS=${tests} \
    NUMFORT_BUILD_EXAMPLES=${examples} NUMFORT_USE_OPENMP=${openmp} \
    NUMFORT_USE_LTO=${ipo} NUMFORT_BUILD_BLAS95=${blas95} NUMFORT_USE_MKL=${mkl} \
    NUMFORT_NATIVE_ARCH=${native} \
    bash "${script_dir}/run-tests.sh" "${compiler}" "${config}" > "${log}" 2>&1 || code=$?

# ==============================================================================
# CTest log snapshots
# ==============================================================================
if [[ ${tests} == ON && -f ${test_log} ]]; then
    cp -- "${test_log}" "${run_dir}/${name}-numfort-LastTest.log"
fi

# ==============================================================================
# Result reporting
# ==============================================================================
status=PASS
if (( code != 0 )); then
    status=FAIL
fi
printf '%s\t%s\t%s\t%s\n' "${name}" "${status}" "${code}" "${log}" \
    > "${run_dir}/${name}.tsv"
printf '%s: %s (exit %s)\nLog: %s\n' \
    "${name}" "${status}" "${code}" "${log}"
exit "${code}"
