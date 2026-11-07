target_sources(numfort PRIVATE
    io.f90
    io_common.F90
    io_fixed.F90
)

if(UNIX)
    target_compile_definitions(numfort PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:_LINUX>")
elseif(WIN32)
    target_compile_definitions(numfort PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:_WIN32>")
endif()
