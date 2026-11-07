

program test_numfort_core_libc

    use, intrinsic :: iso_fortran_env

    use numfort_arrays, only: linspace
    use numfort_core
    use numfort_common_testing

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_logaddexp ()

end subroutine

subroutine test_logaddexp ()

    real (PREC), dimension(:), allocatable :: xx, yy, fxx, fxx1
    integer :: n

    ! Test for non-critical values, ie not for large negative values
    n = 100
    allocate (xx(n), yy(n), fxx(n), fxx1(n))

    call linspace (xx, -20.0_PREC, 10.0_PREC)
    call linspace (yy, -10.0_PREC, 15.0_PREC)

    fxx(:) = log(exp(xx) + exp(yy))
    fxx1(:) = logaddexp(xx, yy)

    if (.not. (all_close (fxx1, fxx, atol=1.0e-10_PREC))) &
        error stop "Test LOGADDEXP(x,y) = LOG(EXP(X) + EXP(Y)) for non-critical " // &
            "values"

    fxx(:) = log(exp(yy) + exp(xx))
    fxx1(:) = logaddexp(yy, xx)
    if (.not. (all_close (fxx1, fxx, atol=1.0e-10_PREC))) &
        error stop "Test LOGADDEXP(y,x) = LOG(EXP(Y) + EXP(X)) for non-critical " // &
            "values"

end subroutine

end
