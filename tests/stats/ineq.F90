program test_numfort_stats_ineq
    !* Unit tests for discrete Gini coefficients in both real precisions.

    use, intrinsic :: iso_fortran_env, only: real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: ieee_exceptions, only: ieee_invalid, ieee_divide_by_zero, &
        ieee_set_flag, ieee_get_flag

    use numfort_stats, only: gini

    implicit none

#ifdef NUMFORT_TEST_REAL32
    integer, parameter :: PREC = real32
#else
    integer, parameter :: PREC = real64
#endif
    real (PREC), parameter :: TOL = 20.0_PREC * epsilon (1.0_PREC)

    call test_gini_zero (.true.)
    call test_gini_zero (.false.)
    call test_gini_zero ()
    call test_gini_positive_degenerate (.true.)
    call test_gini_positive_degenerate (.false.)
    call test_gini_positive_degenerate ()
    call test_gini_weighted (.true.)
    call test_gini_weighted (.false.)
    call test_gini_weighted ()
    call test_gini_top_shares ()
    call test_gini_undefined (.true.)
    call test_gini_undefined (.false.)
    call test_gini_undefined ()
    call test_gini_unsorted (.false.)
    call test_gini_unsorted ()

    contains

subroutine test_gini_zero (assume_sorted)
    !* - Preserve zero inequality for all-zero occupied support.
    !  - Allow repeated support and empty bins at arbitrary positions.

    logical, intent(in), optional :: assume_sorted
        !* Sorting option, forwarded unchanged when absent.

    real (PREC) :: g
    logical :: flags(2)

    call ieee_set_flag ([ieee_invalid, ieee_divide_by_zero], .false.)

    ! A single zero-valued point mass defines zero inequality.
    call gini ([0.0_PREC], [1.0_PREC], g, assume_sorted)
    if (.not. (g == 0.0_PREC)) error stop "Zero point mass: Gini must be zero"

    ! An empty bin after the occupied zero leaves the probability law unchanged.
    call gini ([0.0_PREC, 1.0_PREC], [1.0_PREC, 0.0_PREC], g, assume_sorted)
    if (.not. (g == 0.0_PREC)) error stop "Trailing empty bin: Gini must be zero"

    ! An empty negative bin before the occupied zero also leaves Gini at zero.
    call gini ([-1.0_PREC, 0.0_PREC], [0.0_PREC, 1.0_PREC], g, assume_sorted)
    if (.not. (g == 0.0_PREC)) error stop "Leading empty bin: Gini must be zero"

    ! Several occupied zero entries still represent the same point mass.
    call gini ([0.0_PREC, 0.0_PREC, 0.0_PREC], &
        [0.25_PREC, 0.5_PREC, 0.25_PREC], g, assume_sorted)
    if (.not. (g == 0.0_PREC)) error stop "Repeated zero support: Gini must be zero"

    ! Empty bins precede, separate, and follow occupied zero entries.
    call gini ([-2.0_PREC, 0.0_PREC, 0.0_PREC, 0.0_PREC, 4.0_PREC], &
        [0.0_PREC, 0.25_PREC, 0.0_PREC, 0.75_PREC, 0.0_PREC], g, assume_sorted)
    if (.not. (g == 0.0_PREC)) error stop "Interspersed empty bins: Gini must be zero"

    call ieee_get_flag ([ieee_invalid, ieee_divide_by_zero], flags)
    if (any (flags)) error stop "Zero point mass raised an IEEE exception"
end subroutine


subroutine test_gini_positive_degenerate (assume_sorted)
    !* - Preserve zero inequality for positive point masses.
    !  - Cover repeated entries and empty bins.

    logical, intent(in), optional :: assume_sorted
        !* Sorting option, forwarded unchanged when absent.

    real (PREC) :: g

    ! A positive singleton has no pairwise differences, so Gini is zero.
    call gini ([3.0_PREC], [1.0_PREC], g, assume_sorted)
    if (.not. (abs (g) < TOL)) error stop "Positive singleton: Gini must be zero"

    ! Unequal weights on identical positive values leave all differences at zero.
    call gini ([3.0_PREC, 3.0_PREC], [0.75_PREC, 0.25_PREC], g, assume_sorted)
    if (.not. (abs (g) < TOL)) error stop "Repeated positive support: Gini must be zero"

    ! Empty bins on either side do not change the positive point mass.
    call gini ([-2.0_PREC, 3.0_PREC, 3.0_PREC, 8.0_PREC], &
        [0.0_PREC, 0.25_PREC, 0.75_PREC, 0.0_PREC], g, assume_sorted)
    if (.not. (abs (g) < TOL)) error stop "Positive point mass with empty bins"
end subroutine


subroutine test_gini_weighted (assume_sorted)
    !* - Check unequal weights and signed support with a positive mean.
    !  - Use independent pairwise differences for expected coefficients.

    logical, intent(in), optional :: assume_sorted
        !* Sorting option, forwarded unchanged when absent.

    real (PREC) :: g, x(10), pmf(10)

    ! Equal masses at zero and one: G = (0.5 * 0.5 * 1) / 0.5 = 0.5.
    call gini ([0.0_PREC, 1.0_PREC], [0.5_PREC, 0.5_PREC], g, assume_sorted)
    if (.not. (abs (g - 0.5_PREC) < TOL)) error stop "Equal two-point Gini"

    ! Splitting each mass into five bins preserves G = 0.5.
    x(:5) = 0.0_PREC
    x(6:) = 0.1_PREC
    pmf = 0.1_PREC
    call gini (x, pmf, g, assume_sorted)
    if (.not. (abs (g - 0.5_PREC) < TOL)) error stop "Equal ten-point Gini"

    ! Unequal masses: G = (0.75 * 0.25 * 2) / 1.5 = 0.25.
    call gini ([1.0_PREC, 3.0_PREC], [0.75_PREC, 0.25_PREC], g, assume_sorted)
    if (.not. (abs (g - 0.25_PREC) < TOL)) error stop "Unequal two-point Gini"

    ! Signed support with mean 2: G = (0.25 * 0.75 * 4) / 2 = 0.375.
    call gini ([-1.0_PREC, 3.0_PREC], [0.25_PREC, 0.75_PREC], g, assume_sorted)
    if (.not. (abs (g - 0.375_PREC) < TOL)) error stop "Signed two-point Gini"

    ! Signed support can exceed one: G = (0.25 * 0.75 * 4) / 0.5 = 1.5.
    call gini ([-2.5_PREC, 1.5_PREC], [0.25_PREC, 0.75_PREC], g, assume_sorted)
    if (.not. (abs (g - 1.5_PREC) < TOL)) error stop "Signed Gini above one"

    ! Sum over unordered pairs: G = (0.125 + 0.75 + 0.5) / 1.5 = 11/12.
    call gini ([-2.0_PREC, 0.0_PREC, 4.0_PREC], &
        [0.25_PREC, 0.25_PREC, 0.5_PREC], g, assume_sorted)
    if (.not. (abs (g - 11.0_PREC / 12.0_PREC) < TOL)) &
        error stop "Signed three-point Gini"

    ! Empty bins preserve the same pairwise numerator and mean, hence G = 11/12.
    call gini ([-5.0_PREC, -2.0_PREC, 0.0_PREC, 1.0_PREC, 4.0_PREC, 8.0_PREC], &
        [0.0_PREC, 0.25_PREC, 0.25_PREC, 0.0_PREC, 0.5_PREC, 0.0_PREC], &
        g, assume_sorted)
    if (.not. (abs (g - 11.0_PREC / 12.0_PREC) < TOL)) &
        error stop "Signed three-point Gini with empty bins"
end subroutine


subroutine test_gini_top_shares ()
    !* - Check equal positive values among the top population shares.

    real (PREC) :: x(10), pmf(10), g, expected
    integer :: i

    pmf = 1.0_PREC / real (size (pmf), PREC)
    do i = 0, size (x) - 1
        ! Case i: i zero bins, with the remaining bins sharing one positive value.
        ! G = q * (1-q) * value / ((1-q) * value) = q, where q = i/10.
        x(:i) = 0.0_PREC
        x(i+1:) = 1.234_PREC
        expected = real (i, PREC) / real (size (x), PREC)
        call gini (x, pmf, g, assume_sorted=.true.)
        if (.not. (abs (g - expected) < TOL)) error stop "Equal top-income shares"
    end do
end subroutine


subroutine test_gini_undefined (assume_sorted)
    !* - Return NaN for zero means with nonzero occupied support.
    !  - Reject empty distributions and mismatched arrays without IEEE exceptions.

    logical, intent(in), optional :: assume_sorted
        !* Sorting option, forwarded unchanged when absent.

    real (PREC) :: g, empty(0)
    logical :: flags(2)

    call ieee_set_flag ([ieee_invalid, ieee_divide_by_zero], .false.)

    ! Opposite values with equal masses cancel the mean but have nonzero dispersion.
    call gini ([-1.0_PREC, 1.0_PREC], [0.5_PREC, 0.5_PREC], g, assume_sorted)
    if (.not. ieee_is_nan (g)) error stop "Signed zero mean: expected NaN"

    ! Unequal masses also cancel: -1 * 0.75 + 3 * 0.25 = 0.
    call gini ([-1.0_PREC, 3.0_PREC], [0.75_PREC, 0.25_PREC], g, assume_sorted)
    if (.not. ieee_is_nan (g)) error stop "Unequal masses with zero mean: expected NaN"

    ! Zero total mass does not define a probability distribution, even at zero.
    call gini ([0.0_PREC, 0.0_PREC], [0.0_PREC, 0.0_PREC], g, assume_sorted)
    if (.not. ieee_is_nan (g)) error stop "Zero total mass: expected NaN"

    ! A singleton with no probability mass also does not define a distribution.
    call gini ([3.0_PREC], [0.0_PREC], g, assume_sorted)
    if (.not. ieee_is_nan (g)) error stop "Empty singleton: expected NaN"

    ! Empty arrays contain no occupied support.
    call gini (empty, empty, g, assume_sorted)
    if (.not. ieee_is_nan (g)) error stop "Empty arrays: expected NaN"

    ! Each support value must have a corresponding probability.
    call gini ([0.0_PREC, 1.0_PREC], [1.0_PREC], g, assume_sorted)
    if (.not. ieee_is_nan (g)) error stop "Mismatched arrays: expected NaN"

    call ieee_get_flag ([ieee_invalid, ieee_divide_by_zero], flags)
    if (any (flags)) error stop "Undefined Gini raised an IEEE exception"
end subroutine


subroutine test_gini_unsorted (assume_sorted)
    !* - Sort support and probabilities together, including empty bins.
    !  - Exercise explicit false and omitted sorting arguments.

    logical, intent(in), optional :: assume_sorted
        !* Sorting option, forwarded unchanged when absent.

    real (PREC) :: g
    logical :: flags(2)

    call ieee_set_flag ([ieee_invalid, ieee_divide_by_zero], .false.)

    ! Reordered unequal masses preserve G = (0.75 * 0.25 * 2) / 1.5 = 0.25.
    call gini ([3.0_PREC, 1.0_PREC], [0.25_PREC, 0.75_PREC], g, assume_sorted)
    if (.not. (abs (g - 0.25_PREC) < TOL)) error stop "Unsorted unequal two-point Gini"

    ! Permuting the signed three-point fixture preserves G = 11/12.
    call gini ([4.0_PREC, -5.0_PREC, 0.0_PREC, 8.0_PREC, -2.0_PREC, 1.0_PREC], &
        [0.5_PREC, 0.0_PREC, 0.25_PREC, 0.0_PREC, 0.25_PREC, 0.0_PREC], &
        g, assume_sorted)
    if (.not. (abs (g - 11.0_PREC / 12.0_PREC) < TOL)) &
        error stop "Unsorted signed Gini with empty bins"

    ! Reordered empty bins leave the occupied support at zero.
    call gini ([4.0_PREC, 0.0_PREC, -2.0_PREC, 0.0_PREC, 0.0_PREC], &
        [0.0_PREC, 0.25_PREC, 0.0_PREC, 0.0_PREC, 0.75_PREC], g, assume_sorted)
    if (.not. (g == 0.0_PREC)) error stop "Unsorted zero point mass"

    ! Reordered empty bins also preserve zero inequality at a positive point mass.
    call gini ([8.0_PREC, 3.0_PREC, -2.0_PREC, 3.0_PREC], &
        [0.0_PREC, 0.25_PREC, 0.0_PREC, 0.75_PREC], g, assume_sorted)
    if (.not. (abs (g) < TOL)) error stop "Unsorted positive point mass"

    ! Reordering a signed zero-mean law leaves its coefficient undefined.
    call gini ([3.0_PREC, -1.0_PREC], [0.25_PREC, 0.75_PREC], g, assume_sorted)
    if (.not. ieee_is_nan (g)) error stop "Unsorted signed zero mean: expected NaN"

    call ieee_get_flag ([ieee_invalid, ieee_divide_by_zero], flags)
    if (any (flags)) error stop "Unsorted Gini raised an IEEE exception"
end subroutine

end program
