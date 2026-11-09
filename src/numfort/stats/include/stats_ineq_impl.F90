

subroutine gini_impl (x, pmf, g, assume_sorted)
    !* Return the Gini coefficient for a discrete probability distribution.
    !
    ! For nonzero mean, G = sum_i sum_j p_i p_j abs(x_i-x_j) / (2 * mean).
    ! Signed support is allowed; the coefficient need not lie in [0, 1].
    ! All-zero occupied support has G = 0, regardless of empty bins.
    ! A zero mean with nonzero occupied support is undefined and returns a
    ! quiet NaN. Empty arrays, mismatched sizes, and zero total mass also
    ! return a quiet NaN.
    !
    real (PREC), intent(in), dimension(:), contiguous, target :: x
        !* Finite support values, with one entry per probability.
    real (PREC), intent(in), dimension(:), contiguous, target :: pmf
        !* Finite, nonnegative probabilities that sum to one; zero entries are allowed.
    real (PREC), intent(out) :: g
        !* Gini coefficient, or a quiet NaN when undefined.
    logical, intent(in), optional :: assume_sorted
        !* Whether X is already in ascending order. Defaults to false.

    integer, dimension(:), allocatable :: iorder
    real (PREC), dimension(:), allocatable, target :: x_sort, pmf_sort
    real (PREC), dimension(:), contiguous, pointer :: ptr_x, ptr_pmf
    logical :: lassume_sorted

    real (PREC) :: s1, s2
    integer :: n, i
    lassume_sorted = .false.


    if (present(assume_sorted)) lassume_sorted = assume_sorted

    g = ieee_value (0.0_PREC, IEEE_QUIET_NAN)

    n = size (x)
    if (n == 0 .or. n /= size (pmf)) return

    ! --- Ensure sorted inputs ---

    if (lassume_sorted) then
        ! Arrays are already sorted, point to dummy arguments
        ptr_x => x
        ptr_pmf => pmf
    else
        ! Sort arrays
        allocate (iorder(n), x_sort(n), pmf_sort(n))

        ! Get sort order
        call argsort (x, iorder)

        do i = 1, n
            x_sort(i) = x(iorder(i))
            pmf_sort(i) = pmf(iorder(i))
        end do
        
        ptr_x => x_sort
        ptr_pmf => pmf_sort
    end if

    ! --- Compute Gini ---

    g = 0.0_PREC
    s1 = 0.0_PREC
    s2 = 0.0_PREC

    do i = 1, n
        s2 = s1 + ptr_x(i) * ptr_pmf(i)
        g = g + (s1 + s2) * ptr_pmf(i)
        s1 = s2
    end do

    if (s2 == 0.0_PREC) then
        ! Distinguish zero occupied support from cancellation of signed values.
        g = ieee_value (0.0_PREC, IEEE_QUIET_NAN)
        if (any (pmf > 0.0_PREC) &
            .and. all (x == 0.0_PREC .or. pmf == 0.0_PREC)) g = 0.0_PREC
        return
    end if

    g = 1.0_PREC - g / s2

end subroutine
