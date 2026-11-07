function(numfort_add_test name source)
    set(target numfort_test_${name})
    add_executable(${target} "${source}" ${ARGN})
    set_target_properties(${target} PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${target}-modules")
    target_link_libraries(${target} PRIVATE numfort::numfort)
    target_include_directories(${target} PRIVATE
        "${PROJECT_SOURCE_DIR}/src/numfort/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/include")
    numfort_configure_compiler_options(${target} modern)
    if(NUMFORT_USE_OPENMP)
        target_link_libraries(${target} PRIVATE OpenMP::OpenMP_Fortran)
    endif()
    if(TARGET NumfortMKL::BLAS95 AND NOT NUMFORT_BUILD_BLAS95)
        target_link_libraries(${target} PRIVATE NumfortMKL::BLAS95)
    endif()
    add_test(NAME numfort_${name} COMMAND ${target})
endfunction()

foreach(component IN ITEMS arrays common core integrate interpolate io linalg
        optimize polynomial stats)
    add_subdirectory(${component})
endforeach()
if(NUMFORT_BUILD_BLAS95)
    add_subdirectory(wrappers)
endif()
if(NUMFORT_USE_OPENMP)
    add_subdirectory(openmp)
endif()
