program test_arrays_kron
    !*  Unit tests for Vandermonde matrix creation routines.

    use, intrinsic :: iso_fortran_env
    use numfort_arrays, only: kron
    use numfort_common_testing, only: all_close

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_kron ()

end subroutine

subroutine test_kron ()

    real (PREC), dimension(:,:), allocatable :: m1, m2, mm, mm_ok

    allocate (m1(2,2), m2(2,2))
    m1(:,:) = reshape([1,2,3,4], shape=[2,2])
    m2(:,:) = reshape([1,0,0,1], shape=[2,2])

    allocate (mm(4,4), mm_ok(4,4))
    call kron (m1, m2, mm)

    mm_ok(1:2,1:2) = m2
    mm_ok(3:4,1:2) = m1(2,1) * m2
    mm_ok(1:2,3:4) = m1(1,2) * m2
    mm_ok(3:4,3:4) = m1(2,2) * m2

    if (.not. (all_close (mm, mm_ok, atol = 0.0_PREC, rtol = 0.0_PREC))) &
        error stop "Kronecker product of two 2x2 matrices"

end subroutine

end program
