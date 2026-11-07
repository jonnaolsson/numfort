program test_numfort_stats_dfrechet

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env

    use numfort_common_testing, only: all_close
    use numfort_arrays
    use numfort_stats, only: dfrechet => dfrechet_real64, &
        cdf, pdf, ppf

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

contains

subroutine test_all ()

    call test_pdf ()
    call test_cdf ()
    call test_ppf ()

end subroutine

subroutine test_pdf ()

    type (dfrechet) :: frechet
    real (PREC) :: x, fx, shp, desired
    logical :: res

    ! Compute to values calculated in Mathematica, taken from Scipy project
    x = 2.5d0
    shp = 0.5d0
    desired = 0.067202904517205343834d0
    fx = pdf (frechet, x, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (fx, desired, rtol=1d-12)
    if (.not. (res)) error stop "PPF(2.500; shape=.500) == Mathematica value"

    x = 3.0d0
    shp = 13.0d0
    desired = 2.7179753508900831584d-6
    fx = pdf (frechet, x, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (fx, desired, rtol=1d-12)
    if (.not. (res)) error stop "PPF(3.000; shape=13.000) == Mathematica value"
end subroutine

subroutine test_cdf ()

    integer, parameter :: NX = 11
    real (PREC), dimension(100) :: m, s, alpha
    real (PREC) :: diff
    real (PREC), dimension(NX) :: z
    real (PREC) :: x, fx
    real (PREC), parameter :: tol=1d-10
    real (PREC) :: desired, shp
    real (PREC), dimension(100) :: power_grid
    real (PREC), dimension(NX) :: power_grid2
    logical :: res

    type (dfrechet) :: frechet

    call linspace (power_grid, 0.0d0, 1.0d0)
    power_grid = power_grid ** 3

    call linspace (power_grid2, 0.0d0, 1.0d0)
    power_grid2 = power_grid2 ** 3

    m = power_grid * 1.0d5
    s = 1.0d-4 + power_grid * 1.0d3
    alpha = 1.0d-4 + power_grid * 1.0d3

    z = 1d-5 + power_grid2 * 1.0d2

    diff = 0

!     do i = 1, size(m)
!         do j = 1, size(s)
!             do k = 1, size(alpha)
!                 x = z * s(j) + m(i)
!                 fx = cdf (frechet, x, loc=m(i), scale=s(j), shape=alpha(k))
!                 x2 = ppf (frechet, fx, loc=m(i), scale=s(j), shape=alpha(k))
!
!                 mask = ieee_is_finite (x2)
!
!                 where (mask)
!                     diff_x = abs(x-x2)
!                 else where
!                     diff_x = 0
!                 end where
!
!                 diff = max(diff, maxval(diff_x))
!             end do
!         end do
!     end do
!

    ! Compare to values calculated in Mathematica, taken from the Scipy project
    x = 2.5d0
    shp = 0.5d0
    fx = cdf(frechet, x, loc=0.0d0, scale=1.0d0, shape=shp)
    desired = 0.53128560913296781152d0
    res = all_close (fx, desired, rtol=1d-12)
    if (.not. (res)) error stop "CDF(2.500; shape=.500) == Mathematica value"

    x = 1d-3
    fx = cdf(frechet, x, loc=0.0d0, scale=1.0d0, shape=shp)
    desired = 1.84672666240969310047d-14
    res = all_close (fx, desired, rtol=1d-12)
    if (.not. (res)) error stop "CDF(.001; shape=.500) == Mathematica value"

    x = 3.0d0
    shp = 13.0d0
    fx = cdf(frechet, x, loc=0.0d0, scale=1.0d0, shape=shp)
    desired = 0.99999937277472231955d0
    res = all_close (fx, desired, rtol=1d-12)
    if (.not. (res)) error stop "CDF(3.000; shape=13.000) == Mathematica value"

end subroutine

subroutine test_ppf ()

    type (dfrechet) :: frechet
    real (PREC) :: q, x, shp, desired
    logical :: res

    ! Compute to values calculated in Mathematica, taken from Scipy project
    q = 0.975d0
    shp = 0.5d0
    desired = 1560.0833306626037620d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    if (.not. (res)) error stop "PPF(.975; shape=.500) == Mathematica value"

    q = 0.5d0
    shp = 0.5d0
    desired = 2.0813689810056077979d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    if (.not. (res)) error stop "PPF(.500; shape=.500) == Mathematica value"

    q = 0.975d0
    shp = 13.0d0
    desired = 1.3268241778612848659d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    if (.not. (res)) error stop "PPF(.975; shape=13.000) == Mathematica value"

    q = 0.5d0
    shp = 13.0d0
    desired = 1.0285944941498834248d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    if (.not. (res)) error stop "PPF(.500; shape=13.000) == Mathematica value"

    q = 0.01d0
    shp = 13.0d0
    desired = 0.88916242414515439302d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    if (.not. (res)) error stop "PPF(.010; shape=13.000) == Mathematica value"

    q = 1d-14
    shp = 13.0d0
    desired =  0.76554999794716080164d0
    x = ppf (frechet, q, loc=0.0d0, scale=1.0d0, shape=shp)
    res = all_close (x, desired, rtol=1d-12)
    if (.not. (res)) error stop "PPF(0.10E-13; shape=13.000) == Mathematica value"
end subroutine

end program
