

program numfort_test_core_cephes_stats
    !*  Unit tests for selected stats-related functions ported from the
    !   the Cephes math library.

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_arrays, only: linspace
    use numfort_core, only: erfinv, ndtr, ndtri
    use numfort_common_testing, only: all_close

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_erfinv ()
    call test_ndtr ()

end subroutine

subroutine test_erfinv ()
    !*  Test ERFINV inverse function by checking for approximate identity
    !   of ERFINV(ERF(x)) where ERF() is the built-in ERF function.

    real (PREC) :: xi, yi, x1
    real (PREC), dimension(:), allocatable :: xx
    integer :: n, i
    real (PREC) :: tol = 1.0e-10_PREC

    n = 20
    allocate (xx(20))
    call linspace (xx, 0.0_PREC, 3.5_PREC)

    do i = 1, size(xx)
        xi = xx(i)
        yi = erf (xi)
        x1 = erfinv (yi)
        if (.not. (abs (x1 - xx(i)) < tol)) &
            error stop "ERFINV round trip did not recover the original value"
    end do

end subroutine

subroutine test_ndtr ()
    !*  Unit tests for Gaussian CDF and inverse CDF routines NDTR and NDTRI.

    real (PREC) :: xi, yi, x1, y1, diff
    real (PREC), dimension(:), allocatable :: xx, yy
    integer :: n, i

    real (PREC) :: tol = 1.0e-10_PREC

    n = 20
    allocate (xx(20))
    call linspace (xx, -4.0_PREC, 4.0_PREC)

    do i = 1, size(xx)
        xi = xx(i)
        yi = ndtr (xi)
        x1 = ndtri (yi)
        diff = abs(x1-xi)
        if (.not. (diff < tol)) &
            error stop "NDTRI applied to NDTR did not recover the input"
    end do

    allocate (yy(20))
    call linspace (yy, 1.0e-10_PREC, 1.0_PREC-1.0e-10_PREC)

    do i = 1, size(yy)
        yi = yy(i)
        xi = ndtri (yi)
        y1 = ndtr (xi)
        diff = abs(y1-yi)
        if (.not. (diff < tol)) &
            error stop "NDTR applied to NDTRI did not recover the probability"
    end do

    ! Test at some well-known values
    yi = 0.0
    xi = ndtri (yi)
    if (.not. (ieee_class (xi) == IEEE_NEGATIVE_INF)) error stop "NDTRI(0.0) == -inf"

    yi = 1.0
    xi = ndtri (yi)
    if (.not. (ieee_class (xi) == IEEE_POSITIVE_INF)) error stop "NDTRI(1.0) == inf"

    xi = ieee_value (xi, IEEE_NEGATIVE_INF)
    yi = ndtr (xi)
    if (.not. (yi == 0.0_PREC)) error stop "NDTR(-inf) == 0.0"

    xi = ieee_value (xi, IEEE_POSITIVE_INF)
    yi = ndtr (xi)
    if (.not. (yi == 1.0_PREC)) error stop "NDTR(inf) == 1.0"

    xi = 0.0
    yi = ndtr (xi)
    if (.not. (abs (yi - 0.5_PREC) < 1.0e-15_PREC)) error stop "NDTR(0.0) == 0.5"

    yi = 0.5
    xi = ndtri (yi)
    if (.not. (abs (xi) < 1.0e-15_PREC)) error stop "NDTRI(0.5) == 0.0"

end subroutine

end program
