module simplex_regression_real32
    !* Simplex termination fixtures in single precision.
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env, only: real32
    use numfort_optimize, optim_result => optim_result_real32
    implicit none
    private
    integer, parameter :: PREC = real32
    integer :: calls, callback_limit
    real(PREC) :: best_value, objective_scale
    real(PREC) :: constant_value = 1.0_PREC
    public :: test_constant, test_partial_flat, test_search_budget
    public :: test_fit_budget, test_quadratic_fit, test_degenerate, test_nonfinite
contains
#include "minimize_simplex_regression.inc"
end module

module simplex_regression_real64
    !* Simplex termination fixtures in double precision.
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env, only: real64
    use numfort_optimize, optim_result => optim_result_real64
    implicit none
    private
    integer, parameter :: PREC = real64
    integer :: calls, callback_limit
    real(PREC) :: best_value, objective_scale
    real(PREC) :: constant_value = 1.0_PREC
    public :: test_constant, test_partial_flat, test_search_budget
    public :: test_fit_budget, test_quadratic_fit, test_degenerate, test_nonfinite
contains
#include "minimize_simplex_regression.inc"
end module

program test_simplex_regression
    use simplex_regression_real32, only: &
        test_constant32 => test_constant, &
        test_partial_flat32 => test_partial_flat, &
        test_search_budget32 => test_search_budget, &
        test_fit_budget32 => test_fit_budget, &
        test_quadratic_fit32 => test_quadratic_fit, &
        test_degenerate32 => test_degenerate, &
        test_nonfinite32 => test_nonfinite
    use simplex_regression_real64, only: &
        test_constant64 => test_constant, &
        test_partial_flat64 => test_partial_flat, &
        test_search_budget64 => test_search_budget, &
        test_fit_budget64 => test_fit_budget, &
        test_quadratic_fit64 => test_quadratic_fit, &
        test_degenerate64 => test_degenerate, &
        test_nonfinite64 => test_nonfinite
    implicit none

    call test_constant32 ()
    call test_partial_flat32 ()
    call test_search_budget32 ()
    call test_fit_budget32 ()
    call test_quadratic_fit32 ()
    call test_degenerate32 ()
    call test_nonfinite32 ()
    call test_constant64 ()
    call test_partial_flat64 ()
    call test_search_budget64 ()
    call test_fit_budget64 ()
    call test_quadratic_fit64 ()
    call test_degenerate64 ()
    call test_nonfinite64 ()
end program
