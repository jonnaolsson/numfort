# Keep legacy compiler options and module files separate from numfort's sources.
function(numfort_add_external name)
    set(target numfort_external_${name})
    add_library(${target} OBJECT ${ARGN})
    get_target_property(library_type numfort TYPE)
    if(library_type STREQUAL "SHARED_LIBRARY")
        set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif()
    set_target_properties(${target} PROPERTIES
        Fortran_MODULE_DIRECTORY "${NUMFORT_EXTERNAL_MODULE_DIRECTORY}")
    target_include_directories(${target} PRIVATE
        "${NUMFORT_EXTERNAL_MODULE_DIRECTORY}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include")
    if(CMAKE_GENERATOR MATCHES "^Visual Studio")
        target_include_directories(${target} PRIVATE
            "${NUMFORT_EXTERNAL_MODULE_DIRECTORY}/$<CONFIG>")
    endif()
    numfort_configure_compiler_options(${target} legacy)
    # Object libraries contribute objects and ordering without install exports.
    target_link_libraries(numfort PRIVATE "$<BUILD_INTERFACE:${target}>")
endfunction()

foreach(component IN ITEMS fitpack lbfgsb_bmnz simplex_csiro minpack newuoa2
        orderpack quadpack CDFLIB90 random sobol slsqp)
    add_subdirectory(${component})
endforeach()
