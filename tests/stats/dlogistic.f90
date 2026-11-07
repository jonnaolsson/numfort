program test_numfort_stats_dlogistic

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_stats, only: logistic, dlogistic => dlogistic_real64, &
        cdf, pdf, ppf

    integer, parameter :: PREC = real64

    call test_all ()

contains

subroutine test_all ()

    call test_cdf ()
    call test_pdf ()

end subroutine

subroutine test_cdf ()

    integer :: i
    real (PREC), dimension(3) :: mu, sigma
    real (PREC) :: diff
    real (PREC), dimension(11) :: x, fx, x2, z
    real (PREC), parameter :: tol = 1d-10

    mu = [ real(PREC) :: 0.0, 1.0, -1.0]
    sigma = [ real (PREC) :: 1.0, 0.01, 100.0]

    call linspace (z, - 1.0d1, 1.0d1)

    do i = 1, size(mu)
        x = z * sigma(i) + mu(i)
        fx = cdf (logistic, x, loc=mu(i), scale=sigma(i))
        x2 = ppf (logistic, fx, loc=mu(i), scale=sigma(i))

        diff = maxval(abs(x-x2))
        if (.not. (diff < tol)) error stop "Logistic CDF/inverse-CDF round trip failed"
    end do

end subroutine

subroutine test_pdf ()

    real (PREC) :: fx, x
    type (dlogistic) :: logistic0 = dlogistic(scale=1.0d0)

    ! Test with every small x that should return 0.0
    x = -1.0d5
    fx = pdf (logistic0, x)
    if (.not. (ieee_is_finite (fx) .and. (fx == 0.0))) &
        error stop "PDF at very large negative values"

    ! Test with very large x that should return 0.0
    x = 1.0d5
    fx = pdf (logistic0, x)
    if (.not. (ieee_is_finite (fx) .and. (fx == 0.0))) &
        error stop "PDF at very large positive values"

end subroutine

end program
