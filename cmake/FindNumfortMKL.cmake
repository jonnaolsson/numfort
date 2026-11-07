# Adapt the variable-based finder without exporting machine-specific paths.
# A function confines its compiler/search variables to this dependency lookup.
function(_numfort_find_mkl)
    list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_FUNCTION_LIST_DIR}")
    # Avoid the old architecture probe writing/running in a parent's build root.
    math(EXPR TARGET_ARCH_BITS "${CMAKE_SIZEOF_VOID_P} * 8")
    set(TARGET_ARCH "${TARGET_ARCH_BITS}")
    if(NOT CMAKE_C_COMPILER_LOADED AND NOT CMAKE_CXX_COMPILER_LOADED)
        # FindThreads requires C/CXX. Fortran-only consumers can link the
        # platform thread library directly without enabling another language.
        set(CMAKE_DISABLE_FIND_PACKAGE_Threads TRUE)
        if(UNIX)
            set(CMAKE_THREAD_LIBS_INIT -lpthread)
        else()
            set(CMAKE_THREAD_LIBS_INIT "")
        endif()
    endif()
    find_package(MKL QUIET MODULE COMPONENTS SEQUENTIAL LP64
        OPTIONAL_COMPONENTS ${NumfortMKL_FIND_COMPONENTS})
    set(NumfortMKL_FOUND ${MKL_FOUND} PARENT_SCOPE)
    if(NOT MKL_FOUND)
        return()
    endif()
    if(NOT TARGET NumfortMKL::MKL)
        add_library(NumfortMKL::MKL INTERFACE IMPORTED)
        set(libraries ${MKL_LIBRARIES})
        if(MKL_BLAS95_FOUND)
            list(REMOVE_ITEM libraries "${MKL_BLAS95_LIBRARY}")
        endif()
        set_target_properties(NumfortMKL::MKL PROPERTIES
            INTERFACE_LINK_LIBRARIES "${libraries}")
    endif()
    if(MKL_BLAS95_FOUND AND "BLAS95" IN_LIST NumfortMKL_FIND_COMPONENTS
            AND NOT TARGET NumfortMKL::BLAS95)
        add_library(NumfortMKL::BLAS95 INTERFACE IMPORTED)
        set_target_properties(NumfortMKL::BLAS95 PROPERTIES
            INTERFACE_LINK_LIBRARIES "${MKL_BLAS95_LIBRARY};NumfortMKL::MKL"
            INTERFACE_INCLUDE_DIRECTORIES "${MKL_BLAS95_INCLUDE_DIR}")
    endif()
endfunction()

_numfort_find_mkl()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NumfortMKL REQUIRED_VARS NumfortMKL_FOUND)
