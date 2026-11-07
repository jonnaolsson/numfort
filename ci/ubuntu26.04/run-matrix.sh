#!/usr/bin/env bash

# Schedule numfort's compiler and build-option matrix with CTest.
#
# - Run inside the Fortran CI image with numfort sources mounted.
# - List configurations with --list, or select jobs by their matrix names.
# - Select all configurations when no job names are given.
# - Limit concurrent configurations with --jobs N (default: 1).
# - Use with-oneapi for selections that include ifx or MKL.
# - Save the manifest, worker logs, and ordered summary under BUILD_ROOT/matrix.
# - Return nonzero for failed or incomplete jobs.

# ==============================================================================
# Execution environment
# ==============================================================================
set -euo pipefail
export LC_ALL=C
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

# ==============================================================================
# Matrix definition
# ==============================================================================
# Deduplicate identical configurations while keeping readable, stable job names.
declare -A seen=()
declare -A jobs=()
order=()
add_job() {
    local name=${1} key row
    shift
    key="${*}"
    if [[ -n ${seen[${key}]:-} ]]; then
        return
    fi
    seen[${key}]=1
    printf -v row '%s\t' "${@}"
    jobs[${name}]=${row%$'\t'}
    order+=("${name}")
}

# Every compiler runs full Debug/static and Release/shared builds.
for major in 11 12 13 14 15 16; do
    add_job "gnu${major}-debug-static" "${major}" Debug OFF ON ON ON OFF OFF OFF
    add_job "gnu${major}-release-shared" "${major}" Release ON ON ON ON OFF OFF OFF
done

# On Ubuntu's default compiler, cross tests/examples/shared/OpenMP in Debug.
for shared in OFF ON; do
    for tests in OFF ON; do
        for examples in OFF ON; do
            for openmp in OFF ON; do
                add_job "gnu15-options-${shared}-${tests}-${examples}-${openmp}" \
                    15 Debug "${shared}" "${tests}" "${examples}" "${openmp}" OFF OFF OFF
            done
        done
    done
done

# Exercise IPO archives and explicit wrappers without tests/examples consumers.
add_job gnu15-release-static-ipo 15 Release OFF ON ON ON ON OFF OFF
add_job gnu15-minimal-blas95 15 Debug OFF OFF OFF OFF OFF ON OFF

# Separate compiler coverage from MKL integration, keeping the option cross
# small.
add_job gnu15-mkl-debug-static 15 Debug OFF ON ON ON OFF OFF ON
add_job gnu15-mkl-release-shared 15 Release ON ON ON ON OFF OFF ON
add_job ifx-debug-static ifx Debug OFF ON ON ON OFF OFF OFF
add_job ifx-release-shared ifx Release ON ON ON ON OFF OFF OFF
add_job ifx-mkl-debug-static ifx Debug OFF ON ON ON OFF OFF ON
add_job ifx-mkl-release-shared ifx Release ON ON ON ON OFF OFF ON

# ==============================================================================
# Job selection
# ==============================================================================
parallel=1
if [[ ${1:-} == --jobs ]]; then
    if [[ ! ${2:-} =~ ^[1-9][0-9]*$ ]]; then
        printf 'Usage: run-matrix.sh [--jobs N] [job ...]\n' >&2
        exit 2
    fi
    parallel=${2}
    shift 2
fi
if [[ ${1:-} == --list ]]; then
    if (( $# != 1 )); then
        printf 'Use --list by itself.\n' >&2
        exit 2
    fi
    printf 'job\tcompiler\tconfiguration\tshared\ttests\texamples\topenmp\tipo\tblas95\tmkl\n'
    for name in "${order[@]}"; do
        printf '%s\t%s\n' "${name}" "${jobs[${name}]}"
    done
    exit
fi
selected=("${order[@]}")
if (( $# > 0 )); then
    selected=("${@}")
fi
# Validate the whole selection before starting any builds.
declare -A selected_names=()
for name in "${selected[@]}"; do
    if [[ -z ${jobs[${name}]:-} ]]; then
        printf 'Unknown matrix job: %s (use --list).\n' "${name}" >&2
        exit 2
    fi
    if [[ -n ${selected_names[${name}]:-} ]]; then
        printf 'Duplicate matrix job: %s\n' "${name}" >&2
        exit 2
    fi
    selected_names[${name}]=1
done

# ==============================================================================
# Run directories and job manifest
# ==============================================================================
mkdir -p "${BUILD_ROOT:-/build}"
BUILD_ROOT=$(cd -- "${BUILD_ROOT:-/build}" && pwd)
export BUILD_ROOT
output=${BUILD_ROOT}/matrix
mkdir -p "${output}"
run_dir=$(mktemp -d "${output}/run.XXXXXX")
run_dir=$(cd -- "${run_dir}" && pwd)
for name in "${selected[@]}"; do
    printf '%s\t%s\n' "${name}" "${jobs[${name}]}"
done > "${run_dir}/jobs.tsv"

# ==============================================================================
# CTest scheduling
# ==============================================================================
cmake -S "${script_dir}/matrix" -B "${run_dir}/scheduler" \
    "-DMATRIX_RUN_DIR=${run_dir}"
result=0
ctest --test-dir "${run_dir}/scheduler" --parallel "${parallel}" \
    --output-on-failure || result=1

# ==============================================================================
# Ordered result summary
# ==============================================================================
# Workers write separate rows; merge them in selection order after CTest exits.
summary=${run_dir}/summary.tsv
printf 'job\tresult\texit_code\tlog\n' > "${summary}"
for name in "${selected[@]}"; do
    if [[ -f ${run_dir}/${name}.tsv ]]; then
        cat -- "${run_dir}/${name}.tsv" >> "${summary}"
    else
        printf '%s\tINCOMPLETE\t-\t%s\n' "${name}" "${run_dir}/${name}.log" >> "${summary}"
        result=1
    fi
done
cat -- "${summary}"
printf 'Matrix summary: %s\n' "${summary}"
exit "${result}"
