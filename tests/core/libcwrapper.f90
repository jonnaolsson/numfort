

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

    call test_log1p ()

end subroutine

subroutine test_log1p ()

    real (PREC), dimension(:), allocatable :: xx, fxx, fxx1
    integer :: n

    n = 20

    allocate (xx(n), fxx(n), fxx1(n))

    ! 1. Compute for values where LOG1P(x) = LOG(1 + x) as x is large enough
    call linspace (xx, 1.0_PREC, 10.0_PREC)

    fxx(:) = log(1.0_PREC + xx)
    fxx1(:) = log1p (xx)

    if (.not. (all_close (fxx1, fxx, atol = 1.0e-12_PREC))) &
        error stop "Checking LOG1P(x) = LOG(1 + X) for sufficiently large x"

end subroutine

end
