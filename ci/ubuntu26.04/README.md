# Fortran CI image: Ubuntu 26.04

- Reusable tools image for Fortran projects; sources and test scripts are mounted.
- One tools image for gfortran **11, 12, 13, 14, 15, and 16**, with matching GCC.
- Intel ifx and CPU MKL development libraries are included in the same image.
- Ubuntu repository packages: CMake, Ninja, Git, reference BLAS/LAPACK, and headers.
- Compiler 15 is Ubuntu's default; jobs select versioned executables explicitly.
- Compiler 16 is the `16-20260322` development snapshot.
- Installed package versions are recorded in `/opt/fortran-ci/packages.tsv`.

The image definition and package inventory live in
[fortran-ci](https://github.com/richardfoltyn/fortran-ci). Build it from a sibling
checkout:

```sh
docker build -t fortran-ci:ubuntu26.04 ../fortran-ci/ubuntu26.04
```

Run a job with initialized numfort submodules:

```sh
mkdir -p build-ci
docker run --rm --ulimit core=0 --user "$(id -u):$(id -g)" \
    -v "$(pwd):/src/numfort:ro,z" \
    -v "$(pwd)/build-ci:/build:z" \
    fortran-ci:ubuntu26.04 bash /src/numfort/ci/ubuntu26.04/run-tests.sh 15 Debug
```

- `:z` supports SELinux hosts and sharing mounts across jobs.
- Core dumps are disabled to avoid host desktop crash-report popups; CTest
  still records crashes as failures.
- Set `-e BUILD_SHARED_LIBS=ON` for shared-library jobs; the default is static.
- Select `ifx` instead of a GNU major version for Intel Fortran jobs; these use
  `gcc-15` for C. numfort and the installed client use the same Fortran compiler.
- For ifx or MKL jobs, insert `with-oneapi` before `bash` in the container command
  to initialize compiler paths, runtime library paths, and `MKLROOT`.
- Set `-e NUMFORT_USE_MKL=ON` to test MKL (default OFF). MKL jobs require MKLROOT
  and fail if CMake falls back to LAPACK. MKL uses LP64, sequential, shared
  libraries; `BUILD_SHARED_LIBS` selects the numfort library type independently.
- Set `-e NUMFORT_NATIVE_ARCH=ON` to build with `-march=native` host optimizations.
- BLAS95 discovery checks compiler compatibility; numfort builds its own wrappers
  when compatible MKL wrappers are unavailable.
- Set `-e CMAKE_BUILD_PARALLEL_LEVEL=N` to change build parallelism (default 4).
- Build output is separated by compiler, configuration, and every build option.
- Each job builds numfort, installs and relocates its package, and runs an
  installed client. Tests/examples/OpenMP follow the selected switches.
- Keep `build-ci` as the job artifact, including `Testing/Temporary/LastTest.log`.
- Executable exit status and CTest results determine test failure.

## Build-option matrix

With the same Docker mounts, replace `run-tests.sh 15 Debug` with:

```sh
with-oneapi bash /src/numfort/ci/ubuntu26.04/run-matrix.sh --list
with-oneapi bash /src/numfort/ci/ubuntu26.04/run-matrix.sh
```

- The matrix has 40 unique jobs: Debug/static and Release/shared for compilers
  11–16; all tests/examples/shared/OpenMP combinations in Debug on compiler 15;
  separate static IPO and explicit BLAS95 jobs; six ifx/MKL jobs; and five
  native host-optimized Release jobs (GNU 15, GNU 16, GNU 15 MKL, ifx, ifx MKL).
- Run selected jobs by passing their names, e.g. `run-matrix.sh gnu11-debug-static`.
  Use this mode for one job per CI matrix worker.
- Each option is also available as an ON/OFF environment variable when running
  `run-tests.sh`: `BUILD_SHARED_LIBS`, `NUMFORT_BUILD_TESTS`,
  `NUMFORT_BUILD_EXAMPLES`, `NUMFORT_USE_OPENMP`, `NUMFORT_USE_LTO`,
  `NUMFORT_BUILD_BLAS95`, `NUMFORT_USE_MKL`, and `NUMFORT_NATIVE_ARCH`.
- Default single jobs enable tests, examples, and OpenMP; shared, IPO, and
  explicit BLAS95 default to OFF (BLAS95 is enabled by CMake when needed).
- The driver continues after failures and returns nonzero if any job exits
  nonzero. Each invocation saves per-job logs, CTest log snapshots, the selected
  job manifest, and `summary.tsv` under `/build/matrix/run.XXXXXX/`.

## Local parallel execution

- Pass `--jobs N` before the job names to run up to N configurations concurrently.
  The default is 1; omitting job names selects the full matrix.
- A standalone CMake project registers the selected configurations with CTest.
  CTest schedules complete build/install/client/test sequences; no additional
  packages or image rebuild are needed.
- Use `--jobs 8` and `CMAKE_BUILD_PARALLEL_LEVEL=4`: at most eight
  configurations, each with up to four compilation workers. The inner
  numfort CTest runs explicitly use one worker, regardless of inherited
  `CTEST_PARALLEL_LEVEL`. OpenMP thread settings remain independent.
- Each configuration retains its own build and temporary directories. The runner
  sets `TMP`, `TEMP`, and `TMPDIR` so tests with fixed filenames remain isolated.
  Duplicate job selections
  are rejected; separate matrix invocations must use different `BUILD_ROOT`
  directories if their selected configurations overlap.
- Summaries preserve selection order. `INCOMPLETE` marks jobs that did not write
  a result, for example after interruption. CTest handles worker scheduling and
  interruption; its diagnostics remain under the invocation's `scheduler/`.

With the same Docker mounts, use this container command for two selected jobs:

```sh
with-oneapi bash /src/numfort/ci/ubuntu26.04/run-matrix.sh --jobs 8 \
    gnu15-mkl-debug-static ifx-mkl-debug-static
```

## ifx and MKL jobs

All six jobs enable tests, examples, and OpenMP, including a relocated installed
client check. Reference jobs explicitly select repository BLAS/LAPACK.

| Compiler | Backend | Debug/static | Release/shared |
| --- | --- | --- | --- |
| GNU 15 | MKL | `gnu15-mkl-debug-static` | `gnu15-mkl-release-shared` |
| ifx | Reference | `ifx-debug-static` | `ifx-release-shared` |
| ifx | MKL | `ifx-mkl-debug-static` | `ifx-mkl-release-shared` |

For example, run the Intel MKL Debug job with the same mounts:

```sh
mkdir -p build-ci
docker run --rm --ulimit core=0 --user "$(id -u):$(id -g)" \
    -v "$(pwd):/src/numfort:ro,z" \
    -v "$(pwd)/build-ci:/build:z" \
    fortran-ci:ubuntu26.04 with-oneapi \
    bash /src/numfort/ci/ubuntu26.04/run-matrix.sh ifx-mkl-debug-static
```
