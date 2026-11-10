numfort_add_test(optimize_minimize_dfls minimize_dfls.f90)
numfort_add_test(optimize_minimize_lbfgsb minimize_lbfgsb.f90)
numfort_add_test(optimize_minimize_simplex minimize_simplex.f90)
numfort_add_test(optimize_minimize_simplex_regression minimize_simplex_regression.F90)
set_tests_properties(numfort_optimize_minimize_simplex
    numfort_optimize_minimize_simplex_regression PROPERTIES TIMEOUT 15)
numfort_add_test(optimize_minimize_slsqp minimize_slsqp.f90)
numfort_add_test(optimize_minimize_slsqp_ng minimize_slsqp_ng.f90)
numfort_add_test(optimize_minpack minpack.f90)
numfort_add_test(optimize_root_broyden root_broyden.f90)
numfort_add_test(optimize_root_newton root_newton.f90)
numfort_add_test(optimize_solver_map solver_map.f90)

# Check that numerical difficulties return through the public API.
foreach(kind IN ITEMS 32 64)
    numfort_add_test(optimize_minimize_dfls_failure_${kind} minimize_dfls_failure.F90)
    set(target numfort_test_optimize_minimize_dfls_failure_${kind})
    target_compile_definitions(${target} PRIVATE TEST_REAL_KIND=${kind})
    # Detect a successful but premature STOP in the library.
    set_property(TARGET ${target} PROPERTY TEST_LAUNCHER
        "${CMAKE_COMMAND};-P;${CMAKE_CURRENT_SOURCE_DIR}/check_dfls_failure.cmake")
endforeach()
