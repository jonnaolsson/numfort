program test_numfort_stats_core_quantile
    !* Test discrete quantiles in both supported real precisions.

    use, intrinsic :: iso_fortran_env, only: real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_negative_inf
    use numfort_common_status
    use numfort_stats, only: quantile

    implicit none

#ifdef TEST_REAL32
    integer, parameter :: PREC = real32
#else
    integer, parameter :: PREC = real64
#endif

    call test_quantile_discrete_sort ()
    call test_quantile_discrete_support ()
    call test_quantile_dispatch ()
    call test_quantile_discrete_invalid ()
    call test_quantile_discrete_nonfinite ()

contains


subroutine test_quantile_discrete_sort ()
    !* - Keep weights aligned with unsorted support.
    !  - Preserve the upper choice at an exact interior CDF jump.

    real (PREC) :: x(3), pmf(3), rnk(3), q(3)
    type (status_t) :: status

    ! Issue #30: equal mass at 1 and 3, with zero mass at 2.
    x = [3.0_PREC, 1.0_PREC, 2.0_PREC]
    pmf = [0.5_PREC, 0.5_PREC, 0.0_PREC]
    rnk = [0.25_PREC, 0.75_PREC, 1.0_PREC]
    call quantile (x, pmf, rnk, q, status=status)
    call check_discrete_quantiles (x, pmf, q, &
        [1.0_PREC, 3.0_PREC, 3.0_PREC], status)

    ! The sorted control describes the same distribution without sorting.
    x = [1.0_PREC, 2.0_PREC, 3.0_PREC]
    pmf = [0.5_PREC, 0.0_PREC, 0.5_PREC]
    call quantile (x, pmf, rnk, q, sort=.false., status=status)
    call check_discrete_quantiles (x, pmf, q, &
        [1.0_PREC, 3.0_PREC, 3.0_PREC], status)

    ! At F(1) = 0.5, select 2; rank 1 must also exclude the zero-mass tail.
    pmf = [0.5_PREC, 0.5_PREC, 0.0_PREC]
    rnk = [0.0_PREC, 0.5_PREC, 1.0_PREC]
    call quantile (x, pmf, rnk, q, status=status)
    call check_discrete_quantiles (x, pmf, q, &
        [1.0_PREC, 2.0_PREC, 2.0_PREC], status)

    ! The exact-jump convention is identical with sorting disabled.
    call quantile (x, pmf, rnk, q, sort=.false., status=status)
    call check_discrete_quantiles (x, pmf, q, &
        [1.0_PREC, 2.0_PREC, 2.0_PREC], status)
end subroutine


subroutine test_quantile_discrete_support ()
    !* - Skip zero weights at both ends and in the interior.
    !  - Normalize nonuniform weights and handle singleton support.
    !  - Accept repeated support values and ranks in arbitrary order.

    real (PREC) :: x(7), pmf(7), rnk(7), q(7), q_ok(7)
    type (status_t) :: status

    ! Weights 2 and 6 give P(X=1) = 2/(2+6) = 1/4 and P(X=4) = 3/4.
    x = [0.0_PREC, 1.0_PREC, 2.0_PREC, 3.0_PREC, &
        4.0_PREC, 5.0_PREC, 6.0_PREC]
    pmf = [0.0_PREC, 2.0_PREC, 0.0_PREC, 0.0_PREC, &
        6.0_PREC, 0.0_PREC, 0.0_PREC]
    rnk = [1.0_PREC, 0.25_PREC, 0.0_PREC, 0.125_PREC, &
        0.5_PREC, 0.25_PREC, 1.0_PREC]
    q_ok = [4.0_PREC, 4.0_PREC, 1.0_PREC, 1.0_PREC, &
        4.0_PREC, 4.0_PREC, 4.0_PREC]
    call quantile (x, pmf, rnk, q, sort=.false., status=status)
    call check_discrete_quantiles (x, pmf, q, q_ok, status)

    ! Default sorting must preserve the same masses after this permutation.
    x = [4.0_PREC, 6.0_PREC, 3.0_PREC, 0.0_PREC, &
        5.0_PREC, 1.0_PREC, 2.0_PREC]
    pmf = [6.0_PREC, 0.0_PREC, 0.0_PREC, 0.0_PREC, &
        0.0_PREC, 2.0_PREC, 0.0_PREC]
    call quantile (x, pmf, rnk, q, status=status)
    call check_discrete_quantiles (x, pmf, q, q_ok, status)

    ! A singleton has the same quantile at every rank, including 0 and 1.
    x(1) = 7.0_PREC
    pmf(1) = 2.0_PREC
    q_ok = 7.0_PREC
    call quantile (x(:1), pmf(:1), rnk, q, status=status)
    call check_discrete_quantiles (x(:1), pmf(:1), q, q_ok, status)

    ! A singleton also works with sorting disabled.
    call quantile (x(:1), pmf(:1), rnk, q, sort=.false., status=status)
    call check_discrete_quantiles (x(:1), pmf(:1), q, q_ok, status)

    ! Only one positive weight remains among unsorted zero-mass states.
    x(:3) = [3.0_PREC, 1.0_PREC, 2.0_PREC]
    pmf(:3) = [0.0_PREC, 0.0_PREC, 5.0_PREC]
    q_ok = 2.0_PREC
    call quantile (x(:3), pmf(:3), rnk, q, status=status)
    call check_discrete_quantiles (x(:3), pmf(:3), q, q_ok, status)

    ! Repeated support values with zero weights still give equal mass at 1, 2.
    x(:5) = [1.0_PREC, 1.0_PREC, 2.0_PREC, 2.0_PREC, 3.0_PREC]
    pmf(:5) = [0.0_PREC, 1.0_PREC, 0.0_PREC, 1.0_PREC, 0.0_PREC]
    q_ok = [2.0_PREC, 1.0_PREC, 1.0_PREC, 1.0_PREC, &
        2.0_PREC, 1.0_PREC, 2.0_PREC]
    call quantile (x(:5), pmf(:5), rnk, q, status=status)
    call check_discrete_quantiles (x(:5), pmf(:5), q, q_ok, status)
end subroutine


subroutine test_quantile_dispatch ()
    !* - Equal-length arrays select discrete quantiles and ignore interpolation.
    !  - An extra support value selects the interpolated bin-edge API.
    !  - Scalar ranks use the same dispatch and exact-jump convention.

    real (PREC) :: x(3), pmf(2), rnk(3), q(3), qs
    type (status_t) :: status

    ! Discrete mass at 0 and 2: midpoint must not interpolate between states.
    x = [0.0_PREC, 2.0_PREC, 4.0_PREC]
    pmf = [0.5_PREC, 0.5_PREC]
    rnk = [0.25_PREC, 0.5_PREC, 0.75_PREC]
    call quantile (x(:2), pmf, rnk, q, interp="midpoint", status=status)
    call check_discrete_quantiles (x(:2), pmf, q, &
        [0.0_PREC, 2.0_PREC, 2.0_PREC], status)

    ! Bin midpoints are (0+2)/2 = 1 and (2+4)/2 = 3.
    call quantile (x, pmf, rnk, q, interp="midpoint", status=status)
    if (status /= NF_STATUS_OK) error stop "Bin dispatch status"
    if (any(q /= [1.0_PREC, 3.0_PREC, 3.0_PREC])) &
        error stop "Bin dispatch midpoint values"

    ! Scalar discrete dispatch also ignores midpoint at the exact CDF jump.
    call quantile (x(:2), pmf, 0.5_PREC, qs, interp="midpoint", &
        sort=.false., status=status)
    call check_discrete_quantiles (x(:2), pmf, [qs], [2.0_PREC], status)

    ! Scalar bin dispatch returns the first bin midpoint, (0+2)/2 = 1.
    call quantile (x, pmf, 0.25_PREC, qs, interp="midpoint", status=status)
    if (status /= NF_STATUS_OK) error stop "Scalar bin dispatch status"
    if (qs /= 1.0_PREC) error stop "Scalar bin dispatch midpoint value"
end subroutine


subroutine test_quantile_discrete_invalid ()
    !* - Reject empty, zero-mass, negative-weight, and unordered distributions.
    !  - Validate array sizes and rank bounds; allow an empty rank request.

    real (PREC) :: x(2), pmf(2), rnk(2), q(2)
    type (status_t) :: status

    x = [1.0_PREC, 2.0_PREC]
    pmf = [1.0_PREC, 1.0_PREC]
    rnk = [0.0_PREC, 1.0_PREC]

    ! An empty distribution has no quantiles even with valid requested ranks.
    call quantile (x(:0), pmf(:0), rnk, q, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Empty discrete distribution"

    ! Empty ranks do not make an empty distribution valid.
    call quantile (x(:0), pmf(:0), rnk(:0), q(:0), status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Empty distribution and ranks"

    ! A valid distribution permits an empty rank request.
    call quantile (x, pmf, rnk(:0), q(:0), status=status)
    if (status /= NF_STATUS_OK) error stop "Empty discrete rank request"

    ! All-zero weights cannot be normalized to a probability distribution.
    call quantile (x, [0.0_PREC, 0.0_PREC], rnk, q, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Zero discrete mass"

    ! Negative weights are invalid even when the total weight is positive.
    call quantile (x, [-1.0_PREC, 2.0_PREC], rnk, q, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Negative discrete mass"

    ! sort=false requires nondecreasing support.
    call quantile (x(2:1:-1), pmf, rnk, q, sort=.false., status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Unordered discrete support"

    ! The output must have one element per requested rank.
    call quantile (x, pmf, rnk, q(:1), status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Discrete output size"

    ! A negative rank is outside the probability interval.
    call quantile (x, pmf, [-0.1_PREC, 1.0_PREC], q, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Negative discrete rank"

    ! A rank above one is outside the probability interval.
    call quantile (x, pmf, [0.0_PREC, 1.1_PREC], q, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Discrete rank above one"

    ! Neither discrete support nor bin edges permit fewer values than weights.
    call quantile (x(:1), pmf, rnk, q, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Quantile support size"
end subroutine


subroutine test_quantile_discrete_nonfinite ()
    !* - Reject NaN and infinite support, weights, and ranks.
    !  - Reject a total weight that cannot be represented as a finite value.

    real (PREC) :: x(2), pmf(2), rnk(2), q(2), bad(3)
    type (status_t) :: status
    integer :: i

    x = [1.0_PREC, 2.0_PREC]
    pmf = [1.0_PREC, 1.0_PREC]
    rnk = [0.0_PREC, 1.0_PREC]
    ! Each input is checked with NaN, positive infinity, and negative infinity.
    bad = [ieee_value (0.0_PREC, ieee_quiet_nan), &
        ieee_value (0.0_PREC, ieee_positive_inf), &
        ieee_value (0.0_PREC, ieee_negative_inf)]

    do i = 1, size (bad)
        ! A nonfinite support value does not define a finite discrete quantile.
        call quantile ([bad(i), x(2)], pmf, rnk, q, status=status)
        if (status /= NF_STATUS_INVALID_ARG) error stop "Nonfinite discrete support"

        ! Nonfinite weights cannot define a normalized probability distribution.
        call quantile (x, [bad(i), pmf(2)], rnk, q, status=status)
        if (status /= NF_STATUS_INVALID_ARG) error stop "Nonfinite discrete weight"

        ! Nonfinite ranks do not belong to the closed interval [0,1].
        call quantile (x, pmf, [bad(i), rnk(2)], q, status=status)
        if (status /= NF_STATUS_INVALID_ARG) error stop "Nonfinite discrete rank"
    end do

    ! Two maximal finite weights overflow their sum and must not produce NaNs.
    pmf = huge (1.0_PREC)
    call quantile (x, pmf, rnk, q, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Nonfinite discrete total mass"
end subroutine


subroutine check_discrete_quantiles (x, pmf, q, q_ok, status)
    !* Check exact expected values and membership in positive-mass support.

    real (PREC), intent(in) :: x(:)
        !* Support values supplied to quantile.
    real (PREC), intent(in) :: pmf(:)
        !* Weights associated with the support values.
    real (PREC), intent(in) :: q(:)
        !* Computed quantiles.
    real (PREC), intent(in) :: q_ok(:)
        !* Independently computed expected quantiles.
    type (status_t), intent(in) :: status
        !* Status returned by quantile.

    integer :: i

    if (status /= NF_STATUS_OK) error stop "Discrete quantile status"
    if (any(q /= q_ok)) error stop "Discrete quantile values"
    do i = 1, size (q)
        if (.not. any(x == q(i) .and. pmf > 0.0_PREC)) &
            error stop "Discrete quantile outside positive-mass support"
    end do
end subroutine

end program
