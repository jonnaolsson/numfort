
program test_optimize_simplex
    !*  Unit tests for Simplex wrapper.

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common_testing
    use numfort_optimize, workspace => workspace_real64, &
        optim_result => optim_result_real64

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_quad ()

end subroutine


subroutine test_quad ()

    type (optim_result) :: res
    integer, parameter :: NDIM = 4
    real (PREC), dimension(NDIM) :: x, x_desired

    x_desired = 0.0

    ! Initial starting point
    call arange (x)

    call minimize_simplex (fcn_quad, x, tol=1.0d-5, res=res)

    if (.not. (res%status == NF_STATUS_OK)) error stop "Exit status"
    if (.not. (all (res%x == x))) error stop "X == RES%X"
    if (.not. (all_close (res%x, x_desired, atol = 1.0d-3))) &
        error stop "Correct minimizer X"
    if (.not. (abs (res%fx(1)) < 1.0d-3)) &
        error stop "Correct minimal function value f(X)"

end subroutine


subroutine fcn_quad (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out) :: fx

    fx = sum(x**2.0_PREC)
end subroutine


end program
