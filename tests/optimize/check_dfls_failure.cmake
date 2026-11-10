# Check both the child exit status and a line printed after all solver calls.
# A bare Fortran STOP can exit successfully without completing the assertions.
# CTest appends the test executable to the launcher command.
execute_process(COMMAND "${CMAKE_ARGV3}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors
    TIMEOUT 30)
if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "DFLS child failed (${result}):\n${output}${errors}")
endif()
if(NOT output MATCHES "DFLS failure returned to caller")
    message(FATAL_ERROR "DFLS child terminated before returning:\n${output}${errors}")
endif()
message(STATUS "${output}")
