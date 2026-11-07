#!/usr/bin/env bash

# Build and validate one numfort configuration inside the Fortran CI image.
#
# - Select a GNU compiler major or ifx, followed by Debug or Release.
# - Read build options from the environment: library type, tests, examples,
#   OpenMP, IPO, BLAS95, and MKL.
# - Use mounted sources; keep builds, installs, and temporary files in
#   BUILD_ROOT.
# - Build a matching fcore installation when numfort tests are enabled.
# - Require the selected BLAS provider, including MKL through with-oneapi.
# - Install numfort, relocate its package, and build and run an installed
#   client.
# - Run enabled fcore and numfort tests sequentially and retain CTest logs.
# - Return build/test failures; run-matrix-job.sh also checks assertion
#   messages.

# ==============================================================================
# Execution environment
# ==============================================================================
set -euo pipefail
export LC_ALL=C
# Container crashes otherwise reach the host's desktop crash reporter.
# Keep CTest failures and logs, but do not collect core dumps for matrix runs.
ulimit -c 0

# ==============================================================================
# Job configuration
# ==============================================================================
if (( $# < 1 || $# > 2 )); then
    printf 'Usage: run-tests.sh <gfortran-major|ifx> [Debug|Release]\n' >&2
    exit 2
fi
compiler=${1}
config=${2:-Debug}
shared=${BUILD_SHARED_LIBS:-OFF}
tests=${NUMFORT_BUILD_TESTS:-ON}
examples=${NUMFORT_BUILD_EXAMPLES:-ON}
openmp=${NUMFORT_USE_OPENMP:-ON}
ipo=${NUMFORT_USE_LTO:-OFF}
blas95=${NUMFORT_BUILD_BLAS95:-OFF}
mkl=${NUMFORT_USE_MKL:-OFF}
if [[ ( ! ${compiler} =~ ^[0-9]+$ && ${compiler} != ifx ) || ! ${config} =~ ^(Debug|Release)$ ]]; then
    printf 'Expected a GNU compiler major or ifx, and Debug/Release.\n' >&2
    exit 2
fi
for value in "${shared}" "${tests}" "${examples}" "${openmp}" "${ipo}" "${blas95}" "${mkl}"; do
    if [[ ! ${value} =~ ^(ON|OFF)$ ]]; then
        printf 'Build options must use ON/OFF.\n' >&2
        exit 2
    fi
done
if [[ ${compiler} == ifx ]]; then
    cc=gcc-15
    fc=ifx
    toolchain=ifx
else
    cc=gcc-${compiler}
    fc=gfortran-${compiler}
    toolchain=gnu-${compiler}
fi
command -v "${cc}" >/dev/null
command -v "${fc}" >/dev/null
blas_options=(-DBLA_VENDOR=Generic)
provider=LAPACK
if [[ ${mkl} == ON ]]; then
    if [[ ! -d ${MKLROOT:-} ]]; then
        printf 'MKL jobs require MKLROOT; invoke the runner through with-oneapi.\n' >&2
        exit 2
    fi
    blas_options=("-DMKL_ROOT=${MKLROOT}" -DMKL_STATIC=OFF)
    provider=MKL
fi

# ==============================================================================
# Source mounts and build directories
# ==============================================================================
# Sources can be read-only; every configuration has its own caches/module files.
numfort_source=${NUMFORT_SOURCE_DIR:-/src/numfort}
fcore_source=${FCORE_SOURCE_DIR:-/src/fcore}
job_dir=${BUILD_ROOT:-/build}/${toolchain}/${config}/shared-${shared}-tests-${tests}-examples-${examples}-openmp-${openmp}-ipo-${ipo}-blas95-${blas95}-mkl-${mkl}
fcore_build=${job_dir}/fcore
fcore_prefix=${job_dir}/fcore-install
numfort_build=${job_dir}/numfort
numfort_prefix=${job_dir}/numfort-install
relocated_prefix=${job_dir}/relocated-install

# Tests use fixed filenames in the system temporary directory.
# Isolate all supported temporary-directory variables across configurations.
mkdir -p "${job_dir}/tmp"
export TMP=${job_dir}/tmp TEMP=${job_dir}/tmp TMPDIR=${job_dir}/tmp

# ==============================================================================
# Compiler and package isolation
# ==============================================================================
# Avoid host flags/search paths and match the C compiler to Fortran explicitly.
unset FFLAGS FCFLAGS CFLAGS CXXFLAGS LDFLAGS CMAKE_PREFIX_PATH fcore_DIR numfort_DIR
common=(-G Ninja "-DCMAKE_C_COMPILER=${cc}" "-DCMAKE_Fortran_COMPILER=${fc}"
    "-DCMAKE_BUILD_TYPE=${config}")
"${fc}" --version
cmake --version

# ==============================================================================
# Build and install fcore
# ==============================================================================
# Only tests require fcore. Other cases verify numfort builds without finding
# it.
# Its module files must be produced by the selected Fortran compiler.
prefix_path=""
if [[ ${tests} == ON ]]; then
    cmake -S "${fcore_source}" -B "${fcore_build}" "${common[@]}" \
        "-DCMAKE_INSTALL_PREFIX=${fcore_prefix}" -DBUILD_SHARED_LIBS=OFF \
        -DFCORE_BUILD_TESTS=ON -DFCORE_BUILD_EXAMPLES=OFF
    cmake --build "${fcore_build}"
    cmake --install "${fcore_build}"
    prefix_path=${fcore_prefix}
fi

# ==============================================================================
# Build and install numfort
# ==============================================================================
# Set every switch, including OFF, and require the selected BLAS provider.
cmake -S "${numfort_source}" -B "${numfort_build}" "${common[@]}" \
    "-DCMAKE_PREFIX_PATH=${prefix_path}" "-DCMAKE_INSTALL_PREFIX=${numfort_prefix}" \
    "-DBUILD_SHARED_LIBS=${shared}" "-DNUMFORT_USE_MKL=${mkl}" "${blas_options[@]}" \
    "-DNUMFORT_BUILD_TESTS=${tests}" "-DNUMFORT_BUILD_EXAMPLES=${examples}" \
    "-DNUMFORT_USE_OPENMP=${openmp}" "-DNUMFORT_USE_LTO=${ipo}" \
    "-DNUMFORT_BUILD_BLAS95=${blas95}"
# The generated package records the resolved provider, including any fallback.
if ! grep -Fqx "if(\"${provider}\" STREQUAL \"MKL\")" "${numfort_build}/numfort-config.cmake"; then
    printf 'Expected BLAS provider %s; inspect the configure log.\n' "${provider}" >&2
    exit 1
fi
cmake --build "${numfort_build}"
# Refresh only generated install trees to exclude stale modules/artifacts.
rm -rf -- "${numfort_prefix}" "${relocated_prefix}"
cmake --install "${numfort_build}"

# ==============================================================================
# Check a relocated installed-package consumer
# ==============================================================================
# Move the package so the original prefix cannot accidentally satisfy discovery.
mv -- "${numfort_prefix}" "${relocated_prefix}"
cmake -S "${numfort_source}/examples/client" -B "${job_dir}/client" "${common[@]}" \
    "-DCMAKE_PREFIX_PATH=${relocated_prefix}" "${blas_options[@]}"
cmake --build "${job_dir}/client"
"${job_dir}/client/client"

# ==============================================================================
# Run tests and retain diagnostics
# ==============================================================================
# Run both suites even if one fails; assertion-only failures may still exit
# zero.
result=0
if [[ ${tests} == ON ]]; then
    ctest --test-dir "${fcore_build}" --parallel 1 --output-on-failure || result=1
    ctest --test-dir "${numfort_build}" --parallel 1 --output-on-failure || result=1
    printf 'CTest logs: %s/Testing/Temporary/LastTest.log\n' "${fcore_build}" "${numfort_build}"
else
    # Disabled tests must leave no registrations, including on repeated runs.
    registrations=$(ctest --test-dir "${numfort_build}" --show-only)
    printf '%s\n' "${registrations}"
    if [[ ${registrations} != *'Total Tests: 0'* ]]; then
        printf 'Unexpected CTest registrations with NUMFORT_BUILD_TESTS=OFF.\n' >&2
        result=1
    fi
fi
exit "${result}"
