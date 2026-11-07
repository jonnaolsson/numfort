# These settings are private to numfort targets and never modify parent flags.
function(numfort_configure_compiler_options target profile)
    target_compile_options(${target} PRIVATE
        "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-frecursive>"
        "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<CONFIG:Debug>>:-fcheck=all;-O0>")
    if(profile STREQUAL "legacy")
        target_compile_options(${target} PRIVATE
            "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=legacy;-Wall;-Wno-uninitialized;-Wno-conversion;-Wno-unused;-Wno-tabs;-Wno-unused-dummy-argument>")
        set(intel_warnings all,nounused,nousage,nodeclarations)
    else()
        target_compile_options(${target} PRIVATE
            "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2018;-fimplicit-none;-Wall;-Wextra;-Wno-compare-reals;-Warray-temporaries;-Wrealloc-lhs;-pedantic;-Wuninitialized;-Wno-do-subscript>")
        set(intel_warnings all)
    endif()
    if(WIN32)
        target_compile_options(${target} PRIVATE
            "$<$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>:/auto;/assume:recursion;/warn:${intel_warnings}>"
            "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>,$<CONFIG:Debug>>:/check:all;/traceback;/Od>")
        if(NOT profile STREQUAL "legacy")
            target_compile_options(${target} PRIVATE
                "$<$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>:/stand:f18>")
        endif()
        if(NUMFORT_OPT_REPORT)
            target_compile_options(${target} PRIVATE
                "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>,$<CONFIG:Release,RelWithDebInfo>>:/Qopt-report:2>")
        endif()
    else()
        target_compile_options(${target} PRIVATE
            "$<$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>:-auto;SHELL:-assume recursion;SHELL:-warn ${intel_warnings}>"
            "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>,$<CONFIG:Debug>>:SHELL:-check all;-traceback;-O0>")
        if(NOT profile STREQUAL "legacy")
            target_compile_options(${target} PRIVATE
                "$<$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>:-std18>")
        endif()
        if(NUMFORT_OPT_REPORT)
            target_compile_options(${target} PRIVATE
                "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>,$<CONFIG:Release,RelWithDebInfo>>:-qopt-report=2>")
        endif()
    endif()
    if(NUMFORT_OPT_REPORT)
        target_compile_options(${target} PRIVATE
            "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<CONFIG:Release,RelWithDebInfo>>:-fopt-info-optimized>")
    endif()
    if(NUMFORT_USE_LTO)
        set_target_properties(${target} PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endfunction()
