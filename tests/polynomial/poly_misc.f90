

program test_polynomial_misc
    !*  Unit tests for misc. polynomial (helper) routines from
    !   NUMFORT_POLYNOMIAL_MISC module.

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common
    use numfort_common_workspace, workspace => workspace_real64
    use numfort_common_testing
    use numfort_polynomial

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_polyshift ()

end subroutine

subroutine test_polyshift ()
    !*  Unit tests for POLYSHIFT routine

    integer, parameter :: KMAX = 4
    real (PREC), parameter :: coefs0(0:KMAX) = [ real (PREC) :: &
            -1.23d0, 2.34d0, -5.678d0, 8.234d0, 1.789d0 ]
    real (PREC), dimension(:), allocatable :: coefs1
    real (PREC), dimension(:), allocatable :: x0, x1, y0, y1
    real (PREC) :: xshift
    integer :: k, n
    type (status_t) :: status

    ! Test input checking
    allocate (coefs1(1))
    xshift = 0.0

    call polyshift (coefs0, xshift, coefs1, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Input checks: size(COEFS0) /= size(COEFS1)"

    deallocate (coefs1)

    ! Test with no shift
    xshift = 0.0
    do k = 0, KMAX
        allocate (coefs1(k+1))

        call polyshift (coefs0(0:k), xshift, coefs1, status)

        if (.not. (status == NF_STATUS_OK .and. all_close (coefs1, coefs0(0:k), &
            atol=1.0d-10, rtol=0.0d0))) &
            error stop "POLYSHIFT with zero shift changed coefficients"

        deallocate (coefs1)
    end do

    ! Test with positive shift
    xshift = 1.345d0
    n = 21
    do k = 0, KMAX
        allocate (coefs1(k+1))
        allocate (x0(n), y0(n), x1(n), y1(n))

        call linspace (x0, -2.0d0, 2.0d0)
        x1(:) = x0 - xshift

        call polyshift (coefs0(0:k), xshift, coefs1, status)

        ! Check that polynomials evaluate to the same number
        call polyval (coefs0(0:k), x0, y0)
        call polyval (coefs1, x1, y1)

        if (.not. (status == NF_STATUS_OK .and. all_close (y0, y1, atol=1.0d-10, rtol=0.0d0))) &
            error stop "POLYSHIFT with positive shift changed polynomial values"

        deallocate (coefs1)
        deallocate (x0, x1, y0, y1)
    end do

    ! Test with negative shift
    xshift = -7.234d0
    n = 31
    do k = 0, KMAX
        allocate (coefs1(k+1))
        allocate (x0(n), y0(n), x1(n), y1(n))

        call linspace (x0, -10.0d0, 10.0d0)
        x1(:) = x0 - xshift

        call polyshift (coefs0(0:k), xshift, coefs1, status)

        ! Check that polynomials evaluate to the same number
        call polyval (coefs0(0:k), x0, y0)
        call polyval (coefs1, x1, y1)

        if (.not. (status == NF_STATUS_OK .and. all_close (y0, y1, atol=1.0d-10, rtol=0.0d0))) &
            error stop "POLYSHIFT with negative shift changed polynomial values"

        deallocate (coefs1)
        deallocate (x0, x1, y0, y1)
    end do

end subroutine

end
