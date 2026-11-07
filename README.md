# Numfort

_numfort_ implements some of the functionality of 
[NumPy](https://numpy.org/) and [SciPy](https://scipy.org/) in
Fortran, albeit on a much less ambitious scale. The main goal
is to provide routines that make solving quantitative macroeconomic
models in Fortran easier.

## Table of Contents

- [Obtaining the code](#obtaining-the-code)
- [Build instructions](#build-instructions)
- [Usage](#usage)
- [License](#license)
- [Authors](#authors)

***
## Obtaining the code

1. To clone the git repository, run 
    
    ```bash
    git clone https://github.com/richardfoltyn/numfort.git --recursive
    ```
        
1.  _numfort_ contains a git submodule with CMake scripts required to build the 
    library. It needs to be initialized using
    ```bash
    git submodule update --init --recursive
    ```
    This checks out the version recorded by the repository. After a `git pull`,
    run the same command to update the submodule to the recorded revision.

***
## Build instructions

### Linux

- Requires [CMake](https://cmake.org/) 4.1 or later.
- Static and shared builds require BLAS/LAPACK with 32-bit integer interfaces
  (LP64). Link dependencies propagate through `numfort::numfort`.
- Prefers Intel MKL, with fallback to CMake's LAPACK discovery.

Set the GCC major version in your shell. On Fedora this sets the build path;
the default compiler must be the matching version.

Bash:
```bash
GCC_VERSION=16
```

Fish:
```fish
set GCC_VERSION 16
```

The remaining commands work in both shells and build outside the source tree.

#### Debian/Ubuntu ####

On Debian/Ubuntu, the compiler executables have version suffixes:

```sh
cmake -S "$HOME/repos/numfort" -B "$HOME/build/gnu/$GCC_VERSION/numfort" \
    -DCMAKE_Fortran_COMPILER=gfortran-$GCC_VERSION \
    -DCMAKE_C_COMPILER=gcc-$GCC_VERSION \
    -DCMAKE_INSTALL_PREFIX="$HOME/.local"
```

#### Fedora ####

On Fedora, the default GCC compiler does not have a version suffix:

```sh
cmake -S "$HOME/repos/numfort" -B "$HOME/build/gnu/$GCC_VERSION/numfort" \
    -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_INSTALL_PREFIX="$HOME/.local"
```

#### Building and installing ####

```sh
cmake --build "$HOME/build/gnu/$GCC_VERSION/numfort" --parallel 8
cmake --install "$HOME/build/gnu/$GCC_VERSION/numfort"
```

For GCC 16 on Fedora, the CMake package is installed under
`~/.local/lib64/numfort-0.1-gnu-16/cmake/` and modules under
`~/.local/include/numfort-0.1-gnu-16/`.

#### Using a specific compiler ####

You may want to specify an alternative Fortran compiler or compile flags (`FFLAGS`).
For example, to build using Intel's Fortran compiler `ifx` and optimize the
code for the host machine architecture, you could run:

```sh
cmake -S "$HOME/repos/numfort" -B "$HOME/build/intel/numfort" \
    -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_Fortran_FLAGS=-xHost \
    -DCMAKE_INSTALL_PREFIX=/path/to/install
```
The current code is tested with:

- GNU `gfortran` 16.2.1
- Intel `ifx` 2026.2.0

***
### Advanced build instructions

The library uses the dependencies below:

1.  BLAS and LAPACK libraries: _numfort_ works either with Intel MKL
    or some other implementation such as OpenBLAS.
    1.  For Intel MKL, the variable `MKL_ROOT` should point to the
        desired MKL installation directory.
        Optionally, `MKL_FORTRAN95_ROOT` can be specified if the Fortran 95
        wrappers for BLAS are available. Tests and examples use MKL's BLAS95
        wrappers when available; otherwise numfort builds its own wrappers.
        LAPACK95 is not required.
        
        On Fedora, after setting `GCC_VERSION` as above and installing
        MKL95, set the MKL version in your shell:

        Bash:
        ```bash
        MKL_VERSION=2026.1
        ```

        Fish:
        ```fish
        set MKL_VERSION 2026.1
        ```

        Then configure with either shell:

        ```sh
        cmake -S "$HOME/repos/numfort" -B "$HOME/build/gnu/$GCC_VERSION/numfort" \
            -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc \
            -DCMAKE_INSTALL_PREFIX="$HOME/.local" \
            -DCMAKE_PREFIX_PATH="$HOME/.local" \
            -DMKL_ROOT="/opt/intel/oneapi/mkl/$MKL_VERSION" \
            -DMKL_FORTRAN95_ROOT="$HOME/.local/share/mkl/$MKL_VERSION/gnu/$GCC_VERSION" \
            -DNUMFORT_BUILD_TESTS=ON -DNUMFORT_BUILD_EXAMPLES=ON
        ```
                
    2.  If another BLAS/LAPACK implementation should be used, set
        `NUMFORT_USE_MKL=OFF` and _numfort_ will use whichever library is found
        by CMake.
         
        On Fedora, configure with:

        ```sh
        cmake -S "$HOME/repos/numfort" -B "$HOME/build/gnu/$GCC_VERSION/numfort" \
            -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc \
            -DCMAKE_INSTALL_PREFIX="$HOME/.local" \
            -DCMAKE_PREFIX_PATH="$HOME/.local" \
            -DNUMFORT_USE_MKL=OFF \
            -DNUMFORT_BUILD_TESTS=ON -DNUMFORT_BUILD_EXAMPLES=ON
        ```

### Build options

- Standalone builds default to Release and a static library.
- Project options use the `NUMFORT_` prefix:
  - `NUMFORT_BUILD_TESTS=OFF`, `NUMFORT_BUILD_EXAMPLES=OFF`: optional programs;
    tests run with `ctest --test-dir <build-dir>`.
  - `NUMFORT_USE_MKL=ON`, `NUMFORT_USE_OPENMP=ON`: dependency selection;
    OpenMP-specific programs follow the OpenMP option.
  - `NUMFORT_BUILD_BLAS95=OFF`: build/install wrappers; enabled automatically
    when tests/examples need them and MKL BLAS95 is unavailable.
  - `NUMFORT_USE_LTO=OFF`, `NUMFORT_OPT_REPORT=OFF`: IPO and GNU/Intel
    optimization reports (Release/RelWithDebInfo). Static Intel IPO builds pass
    the required linker option to consumers.
  - `NUMFORT_INSTALL_DIR_NAME`: override the version/compiler subdirectory.
- `BUILD_SHARED_LIBS=ON` selects a shared library; tests/examples link that
  library.
- Multi-config generators use `--config` for builds and `-C` for CTest. Build
  configurations sequentially, or use separate build trees.
- Standard `CMAKE_INSTALL_{LIB,INCLUDE,BIN}DIR` overrides are supported. Use
  relative paths for relocation; DLLs go below bindir, libraries/packages below libdir.

***
## Usage

- Enable Fortran, then find an installed package and link its target:

```cmake
find_package(numfort 0.1 CONFIG REQUIRED)
target_link_libraries(<target> PRIVATE numfort::numfort)
```

- Set `CMAKE_PREFIX_PATH` to the install prefix, or `numfort_DIR` to its
  compiler-specific `cmake` directory. Set `MKL_ROOT` if MKL needs a search hint.
- Installed packages require the same Fortran compiler family and major version;
  version requests accept `0.1.x`. Use `EXACT` for an exact version.
- For source-tree integration, use the same target:

```cmake
add_subdirectory(path/to/numfort numfort-build)
target_link_libraries(<target> PRIVATE numfort::numfort)
```

- Module paths and link dependencies propagate automatically; compiler flags,
  output/module directories, and the default configuration stay local to numfort.
- Parents running numfort's optional tests must call `enable_testing()`.

***
## License

This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. See the GNU General Public License for more details.

_numfort_ includes several bundled projects in the `src/external` directory
which were distributed under various different licenses by their original authors.
Any changes and additions made to the code in `src/external` are licensed
under a project's respective original license.
See `LICENSES_bundled.txt` and the subdirectories in `src/external` for details.

## Authors

With the exception of the files in `src/external`, _numfort_ was written
by Richard Foltyn. See the `src/external` directories for information
about the authors of these bundled components.
