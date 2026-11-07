configure_file(numfort/numfort.f90.in "${CMAKE_CURRENT_BINARY_DIR}/numfort.f90")

add_library(numfort "${CMAKE_CURRENT_BINARY_DIR}/numfort.f90")
add_library(numfort::numfort ALIAS numfort)
set_target_properties(numfort PROPERTIES
    VERSION "${PROJECT_VERSION}"
    Fortran_MODULE_DIRECTORY "${NUMFORT_MODULE_DIRECTORY}")
target_include_directories(numfort
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/numfort/include"
    PUBLIC
        "$<BUILD_INTERFACE:${NUMFORT_MODULE_DIRECTORY}>"
        "$<BUILD_INTERFACE:${NUMFORT_EXTERNAL_MODULE_DIRECTORY}>"
        "$<INSTALL_INTERFACE:${NUMFORT_INSTALL_INCLUDEDIR}>"
        "$<INSTALL_INTERFACE:${NUMFORT_INSTALL_INCLUDEDIR}/external>")
if(CMAKE_GENERATOR MATCHES "^Visual Studio")
    target_include_directories(numfort PUBLIC
        "$<BUILD_INTERFACE:${NUMFORT_MODULE_DIRECTORY}/$<CONFIG>>"
        "$<BUILD_INTERFACE:${NUMFORT_EXTERNAL_MODULE_DIRECTORY}/$<CONFIG>>")
endif()
numfort_configure_compiler_options(numfort modern)
if(NUMFORT_USE_LTO AND CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    get_target_property(numfort_library_type numfort TYPE)
    if(numfort_library_type STREQUAL "STATIC_LIBRARY")
        # Intel's bitcode archives require IPO when linking the final executable.
        if(WIN32)
            target_link_options(numfort INTERFACE
                "$<$<LINK_LANG_AND_ID:Fortran,IntelLLVM>:-Qipo>")
        else()
            target_link_options(numfort INTERFACE
                "$<$<LINK_LANG_AND_ID:Fortran,IntelLLVM>:-ipo>")
        endif()
    endif()
endif()
if(NUMFORT_BLAS_PROVIDER STREQUAL "MKL")
    target_link_libraries(numfort PUBLIC "$<LINK_ONLY:NumfortMKL::MKL>")
else()
    target_link_libraries(numfort PUBLIC "$<LINK_ONLY:LAPACK::LAPACK>")
endif()
if(NUMFORT_USE_OPENMP)
    target_link_libraries(numfort PRIVATE OpenMP::OpenMP_Fortran)
    target_link_libraries(numfort INTERFACE "$<LINK_ONLY:OpenMP::OpenMP_Fortran>")
endif()

add_subdirectory(external)
add_subdirectory(interfaces)
add_subdirectory(wrappers)
add_subdirectory(numfort)
