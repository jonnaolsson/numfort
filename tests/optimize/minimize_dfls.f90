

program test_optimize_minimize_dfls
    !*  Unit tests for Derivative-free least squares (DFLS) minimizer

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_testing
    use numfort_optimize, optim_result => optim_result_real64

    use numfort_stats

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_quadratic ()

end subroutine

subroutine test_quadratic ()

    integer :: n, m
    real (PREC), dimension(:), allocatable :: x, fx_ok
    real (PREC) :: rhobeg, rhoend
    type (optim_result) :: res

    n = 5
    m = 5
    allocate (x(n), fx_ok(m))

    rhobeg = 1.0
    rhoend = 1.0d-8

    call set_seed (1234)
    call random_number (x)

    x(:) = (x - 0.5) * 10.0_PREC

    call minimize_dfls (fobj_quad, x, m, rhobeg, rhoend, res=res)
    fx_ok(:) = 0.0_PREC

    if (.not. (all_close (res%fx, fx_ok, atol = 1.0d-3))) &
        error stop "Smooth quadratic function"

end subroutine


subroutine fobj_quad (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    fx = x ** 2.0
end subroutine

end program
