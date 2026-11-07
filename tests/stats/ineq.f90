
program test_numfort_stats_ineq
    !*  Module with unit tests for inequality measures in stats module.

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common
    use numfort_common_testing
    use numfort_stats

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_gini ()

end subroutine

subroutine test_gini ()

    real (PREC), dimension(10) :: x, pmf
    real (PREC) :: g, B, g_target
    integer :: i

    ! Degenerate distribution (this will not work if x = )
    x = 0.0
    pmf = 1
    call gini (x(:1), pmf(:1), g, assume_sorted=.true.)
    if (.not. (g == 0.0)) error stop "Degenerate distribution: Gini = 0"

    ! 2-point distribution where one half has nothing
    x(1:2) = [0.0, 1.0]
    pmf = 0.5
    call gini (x(:2), pmf(:2), g, assume_sorted=.true.)
    if (.not. (abs (g - 0.5_PREC) < 1.0e-10)) &
        error stop "Two-point distribution: Gini = 0.5"

    ! 10-point distribution with same logic as before
    x(1:5) = 0.0
    x(6:10) = 0.1
    pmf = 1.0_PREC / size(pmf)
    call gini (x, pmf, g, assume_sorted=.true.)
    if (.not. (abs (g - 0.5_PREC) < 1.0e-10)) &
        error stop "10-point distribution: Gini = 0.5"

    ! Equal distribution among the top N, rest has nothing
    pmf = 1.0_PREC / size(pmf)

    do i = 0, size(x) - 1
        x(1:i) = 0
        x(i+1:) = 1.234
        B = (1.0_PREC - i / real(size(x), PREC)) / 2.0
        g_target = 1.0_PREC - 2.0 * B
        call gini (x, pmf, g, assume_sorted=.true.)
        if (.not. (abs (g - g_target) < 1.0e-10)) &
            error stop "Gini differs from the reference for equal top-income shares"
    end do

end subroutine

end
