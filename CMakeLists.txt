cmake_minimum_required(VERSION 4.1)

project(numfort VERSION 0.1.0 LANGUAGES C Fortran)

if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build configuration" FORCE)
endif()

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(NUMFORT_USE_MKL "Prefer Intel MKL for BLAS/LAPACK" ON)
option(NUMFORT_BUILD_BLAS95 "Build numfort's Fortran 95 BLAS wrappers" OFF)
option(NUMFORT_BUILD_TESTS "Build numfort unit tests" OFF)
option(NUMFORT_BUILD_EXAMPLES "Build numfort examples" OFF)
option(NUMFORT_OPT_REPORT "Generate compiler optimization reports" OFF)
option(NUMFORT_USE_OPENMP "Enable OpenMP in numfort" ON)
option(NUMFORT_USE_LTO "Enable interprocedural optimization in numfort" OFF)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(cmake/CompilerOptions.cmake)

# Resolve link requirements for both static and shared consumers.
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
    "${CMAKE_CURRENT_SOURCE_DIR}/submodules/cmake-modules/Modules")
set(NUMFORT_BLAS_PROVIDER LAPACK)
if(NUMFORT_USE_MKL)
    find_package(NumfortMKL QUIET MODULE)
    if(NumfortMKL_FOUND)
        set(NUMFORT_BLAS_PROVIDER MKL)
    else()
        message(STATUS "Intel MKL not found; using CMake's LAPACK discovery")
    endif()
endif()
if(NUMFORT_BLAS_PROVIDER STREQUAL "LAPACK")
    # The interfaces in src/interfaces use default (32-bit) integers.
    set(BLA_SIZEOF_INTEGER 4)
    find_package(LAPACK REQUIRED)
endif()
if(NUMFORT_USE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()
if(NUMFORT_USE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT numfort_ipo_supported OUTPUT numfort_ipo_error
        LANGUAGES Fortran)
    if(NOT numfort_ipo_supported)
        message(FATAL_ERROR "NUMFORT_USE_LTO requested but unsupported: ${numfort_ipo_error}")
    endif()
endif()

# Only optional programs use BLAS95; prefer MKL's wrappers when available.
if((NUMFORT_BUILD_TESTS OR NUMFORT_BUILD_EXAMPLES) AND NOT NUMFORT_BUILD_BLAS95)
    if(NUMFORT_BLAS_PROVIDER STREQUAL "MKL")
        find_package(NumfortMKL QUIET MODULE OPTIONAL_COMPONENTS BLAS95)
    endif()
    if(NOT TARGET NumfortMKL::BLAS95)
        set(NUMFORT_BUILD_BLAS95 ON)
    endif()
endif()

string(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" NUMFORT_COMPILER_ID)
string(REGEX MATCH "^[0-9]+" NUMFORT_COMPILER_MAJOR "${CMAKE_Fortran_COMPILER_VERSION}")
set(NUMFORT_INSTALL_DIR_NAME
    "numfort-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}-${NUMFORT_COMPILER_ID}-${NUMFORT_COMPILER_MAJOR}"
    CACHE STRING "Version- and compiler-specific installation subdirectory")
set(NUMFORT_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}/${NUMFORT_INSTALL_DIR_NAME}")
set(NUMFORT_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/${NUMFORT_INSTALL_DIR_NAME}")
set(NUMFORT_INSTALL_CMAKEDIR "${NUMFORT_INSTALL_LIBDIR}/cmake")
set(NUMFORT_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
set(NUMFORT_EXTERNAL_MODULE_DIRECTORY "${NUMFORT_MODULE_DIRECTORY}/external")

add_subdirectory(src)
# Generate an empty manifest too, so disabling tests clears old registrations.
enable_testing()
if(NUMFORT_BUILD_TESTS)
    add_subdirectory(tests)
endif()
if(NUMFORT_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

configure_package_config_file(install/config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/numfort-config.cmake"
    INSTALL_DESTINATION "${NUMFORT_INSTALL_CMAKEDIR}")
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/numfort-config-version-base.cmake"
    VERSION "${PROJECT_VERSION}" COMPATIBILITY SameMinorVersion)
configure_file(install/config-version.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/numfort-config-version.cmake" @ONLY)

install(TARGETS numfort EXPORT numfort_targets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}/${NUMFORT_INSTALL_DIR_NAME}"
    LIBRARY DESTINATION "${NUMFORT_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${NUMFORT_INSTALL_LIBDIR}")
install(EXPORT numfort_targets FILE numfortTargets.cmake NAMESPACE numfort::
    DESTINATION "${NUMFORT_INSTALL_CMAKEDIR}" COMPONENT devel)

# Intel's Visual Studio integration appends the configuration automatically.
if(CMAKE_GENERATOR MATCHES "^Visual Studio")
    set(NUMFORT_INSTALL_MODULE_DIRECTORY "${NUMFORT_MODULE_DIRECTORY}/$<CONFIG>")
    set(NUMFORT_INSTALL_EXTERNAL_MODULE_DIRECTORY "${NUMFORT_EXTERNAL_MODULE_DIRECTORY}/$<CONFIG>")
else()
    set(NUMFORT_INSTALL_MODULE_DIRECTORY "${NUMFORT_MODULE_DIRECTORY}")
    set(NUMFORT_INSTALL_EXTERNAL_MODULE_DIRECTORY "${NUMFORT_EXTERNAL_MODULE_DIRECTORY}")
endif()
install(DIRECTORY "${NUMFORT_INSTALL_MODULE_DIRECTORY}/"
    DESTINATION "${NUMFORT_INSTALL_INCLUDEDIR}" COMPONENT devel
    FILES_MATCHING PATTERN "numfort*.mod"
    PATTERN "blas_interfaces.mod" PATTERN "lapack_interfaces.mod"
    PATTERN "external" EXCLUDE)
if(NUMFORT_BUILD_BLAS95)
    install(FILES "${NUMFORT_INSTALL_MODULE_DIRECTORY}/blas95.mod"
        DESTINATION "${NUMFORT_INSTALL_INCLUDEDIR}" COMPONENT devel)
endif()
# Public numfort modules reference modules supplied by the bundled libraries.
install(DIRECTORY "${NUMFORT_INSTALL_EXTERNAL_MODULE_DIRECTORY}/"
    DESTINATION "${NUMFORT_INSTALL_INCLUDEDIR}/external" COMPONENT devel
    FILES_MATCHING PATTERN "*.mod")
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/numfort-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/numfort-config-version.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/numfort-config-version-base.cmake"
    "${CMAKE_CURRENT_SOURCE_DIR}/install/compiler_test.f90"
    DESTINATION "${NUMFORT_INSTALL_CMAKEDIR}" COMPONENT devel)
if(NUMFORT_BLAS_PROVIDER STREQUAL "MKL")
    install(FILES cmake/FindNumfortMKL.cmake
        submodules/cmake-modules/Modules/FindMKL.cmake
        submodules/cmake-modules/Modules/FindTargetArch.cmake
        DESTINATION "${NUMFORT_INSTALL_CMAKEDIR}/modules" COMPONENT devel)
endif()
