program test_numfort_stats_dnorm

    use iso_fortran_env

    use numfort_arrays
    use numfort_common
    use numfort_stats, only: dnorm => dnorm_real64, norm, mean, std, cdf, pdf, rvs, &
        set_seed

    implicit none

    integer, parameter :: PREC = real64

    real (PREC), parameter :: tol = 1d-12

    call test_all ()

contains

subroutine test_all ()

    call test_cdf ()
    call test_cdf_rvs ()

end subroutine

subroutine test_cdf ()

    real (real64), dimension(:), allocatable :: x, fx
    real (real64) :: x1, fx1
    real (real64), parameter, dimension(4) :: &
        means = [0.0d0, -21.0d0, 1d-5, 1.234d4], &
        sd = [1.0d0, 1d-4, 1d1, 1.2345d4]

    integer :: i, n

    do i = 1, size(means)

        x1 = means(i)
        fx1 = cdf (norm, x1, loc=means(i), scale=sd(i))
        if (.not. (abs (fx1 - 0.5d0) < tol)) &
            error stop "Normal CDF at the mean differs from 0.5"

        ! Check boundary behavior
        x1 = huge(x1)
        fx1 = cdf (norm, x1, loc=means(i), scale=sd(i))
        if (.not. (abs (fx1 - 1.0d0) < tol)) &
            error stop "Normal CDF at positive infinity differs from one"

        x1 = -huge(x1)
        fx1 = cdf (norm, x1, loc=means(i), scale=sd(i))
        if (.not. (abs (fx1) < tol)) &
            error stop "Normal CDF at negative infinity differs from zero"

        ! Test that CDF is non-decreasing
        n = 100
        allocate (x(n), fx(n))

        call linspace (x, -1d10, 1d10)
        fx = cdf (norm, x, loc=means(i), scale=sd(i))
        if (.not. (all (fx(2:) - fx(1:n - 1) >= 0.0d0))) &
            error stop "Normal CDF is decreasing"

        deallocate (x, fx)
    end do

end subroutine

subroutine test_cdf_rvs ()

    real (PREC), dimension(:), allocatable :: x, fx
    real (PREC), parameter, dimension(4) :: &
        means = [real (PREC) :: 0.0d0, -21.0d0, 1d-5, 1.234d4], &
        sd = [real (PREC) :: 1.0d0, 1d-4, 1d1, 1.2345d4]
    ! Tolerances calibrated to > 6.5 sigma for n = 100,000 draws to prevent false
    ! rejections across compilers (SE(mean) ~= 9.1e-4, SE(var) ~= 2.4e-4)
    real (PREC), parameter :: tol_mean = 6.0e-3_PREC
    real (PREC), parameter :: tol_var = 2.0e-3_PREC
    real (PREC), parameter :: var_expected = 1.0_PREC / 12.0_PREC
    real (PREC) :: m, s, s2
    integer :: i, n
    type (status_t) :: status

    ! Fix seed for reproducible sampling on a given compiler
    call set_seed (1234)

    n = 100000
    allocate (x(n), fx(n))
    do i = 1, size(means)

        ! draw random sample of normally distributed RV
        call rvs (norm, x, loc=means(i), scale=sd(i))
        ! compute CDF(x): is standard uniformly distributed RV
        fx = cdf (norm, x, loc=means(i), scale=sd(i))

        if (.not. (all (fx >= 0.0_PREC) .and. all (fx <= 1.0_PREC))) &
            error stop "Normal CDF of random samples out of [0, 1] range"

        ! compute mean and variance
        call std (fx, s, m, status=status)
        if (status /= NF_STATUS_OK) &
            error stop "Failed to compute moments of CDF-transformed samples"

        s2 = s ** 2
        if (abs (m - 0.5_PREC) >= tol_mean) &
            error stop "Normal CDF of random samples has unexpected mean"

        if (abs (s2 - var_expected) >= tol_var) &
            error stop "Normal CDF of random samples has unexpected variance"

    end do

    deallocate (x, fx)

end subroutine

end program
