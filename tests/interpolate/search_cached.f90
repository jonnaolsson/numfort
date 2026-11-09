program test_interp_find_cached
    !* Regression tests for caller-provided interpolation search caches.

    use, intrinsic :: iso_fortran_env, only: real64
    use numfort_interpolate
    implicit none

    integer, parameter :: PREC = real64
    real (PREC), parameter :: TOL = 1.0e-12_PREC
    real (PREC), parameter :: KNOTS(*) = &
        [-2.0_PREC, -1.0_PREC, 1.0_PREC, 5.0_PREC, 11.0_PREC]
    ! Nonmonotone queries, with more queries than knots.
    real (PREC), parameter :: X(*) = [ &
        6.5_PREC, & ! Interior of [5, 11]: lower weight (11 - 6.5) / 6 = 3/4.
        -3.0_PREC, & ! Below the grid: clamp to -2, or extrapolate with weight 2.
        -1.0_PREC, & ! Interior knot: bracket 2, lower weight 1.
        17.0_PREC, & ! Above the grid: clamp to 11, or use (11 - 17) / 6 = -1.
        1.0_PREC, & ! Interior knot: bracket 3, lower weight 1.
        0.0_PREC, & ! Midpoint of [-1, 1]: lower weight 1/2.
        11.0_PREC, & ! Last knot: bracket 4, lower weight 0.
        -2.0_PREC, & ! First knot: bracket 1, lower weight 1.
        6.5_PREC] ! Repeat the first query after visiting both boundaries.
    integer, parameter :: LOWER(*) = [4, 1, 2, 4, 3, 2, 4, 1, 4]
    real (PREC), parameter :: WEIGHT_BOUNDARY(*) = &
        [0.75_PREC, 1.0_PREC, 1.0_PREC, 0.0_PREC, 1.0_PREC, &
        0.5_PREC, 0.0_PREC, 1.0_PREC, 0.75_PREC]
    real (PREC), parameter :: WEIGHT_EXTRAPOLATE(*) = &
        [0.75_PREC, 2.0_PREC, 1.0_PREC, -1.0_PREC, 1.0_PREC, &
        0.5_PREC, 0.0_PREC, 1.0_PREC, 0.75_PREC]

    call test_scalar_boundary ()
    call test_vector_boundary ()
    call test_vector_extrapolate ()
    call test_cache_indices ()
    call test_cache_grid_resize ()

contains

subroutine test_scalar_boundary ()
    !* - Check scalar cached and uncached brackets and lower weights.
    !  - Reuse one cache through nonmonotone queries and repeated points.

    type (search_cache) :: cache
    type (status_t) :: status
    integer :: i, ilb
    real (PREC) :: weight

    do i = 1, size (X)
        ! Each query above has an independently specified bracket and weight.
        call interp_find (X(i), KNOTS, ilb, weight, &
            ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
        call check_scalar (status, ilb, weight, LOWER(i), WEIGHT_BOUNDARY(i))
        call interp_find (X(i), KNOTS, ilb, weight, &
            ext=NF_INTERP_EVAL_BOUNDARY, status=status)
        call check_scalar (status, ilb, weight, LOWER(i), WEIGHT_BOUNDARY(i))
    end do
    if (cache%i /= 4) error stop "Scalar search did not update caller cache"
end subroutine


subroutine test_vector_boundary ()
    !* - Check vector boundary evaluation with long, short, and empty queries.
    !  - Reuse the cache for reversed queries after an empty call.

    type (search_cache) :: cache
    type (status_t) :: status
    integer :: ilb(size (X)), n
    real (PREC) :: weight(size (X)), queries(size (X))

    n = size (X)
    ! More queries than knots, including repeated points and both exteriors.
    call interp_find (X, KNOTS, ilb, weight, &
        ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
    call check_vector (status, ilb, weight, LOWER, WEIGHT_BOUNDARY)
    if (cache%i /= 4) error stop "Vector search did not update caller cache"
    call interp_find (X, KNOTS, ilb, weight, &
        ext=NF_INTERP_EVAL_BOUNDARY, status=status)
    call check_vector (status, ilb, weight, LOWER, WEIGHT_BOUNDARY)

    ! Fewer queries than knots: interior, lower exterior, and exact knot.
    call interp_find (X(:3), KNOTS, ilb(:3), weight(:3), &
        ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
    call check_vector (status, ilb(:3), weight(:3), LOWER(:3), WEIGHT_BOUNDARY(:3))
    call interp_find (X(:3), KNOTS, ilb(:3), weight(:3), &
        ext=NF_INTERP_EVAL_BOUNDARY, status=status)
    call check_vector (status, ilb(:3), weight(:3), LOWER(:3), WEIGHT_BOUNDARY(:3))

    ! Empty queries must succeed with conformable empty output arrays.
    call interp_find (X(:0), KNOTS, ilb(:0), weight(:0), &
        ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
    call check_vector (status, ilb(:0), weight(:0), LOWER(:0), WEIGHT_BOUNDARY(:0))
    call interp_find (X(:0), KNOTS, ilb(:0), weight(:0), &
        ext=NF_INTERP_EVAL_BOUNDARY, status=status)
    call check_vector (status, ilb(:0), weight(:0), LOWER(:0), WEIGHT_BOUNDARY(:0))

    ! Reverse the query order and reuse the cache left by the empty call.
    queries = X(n:1:-1)
    call interp_find (queries, KNOTS, ilb, weight, &
        ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
    call check_vector (status, ilb, weight, LOWER(n:1:-1), WEIGHT_BOUNDARY(n:1:-1))
    call interp_find (queries, KNOTS, ilb, weight, &
        ext=NF_INTERP_EVAL_BOUNDARY, status=status)
    call check_vector (status, ilb, weight, LOWER(n:1:-1), WEIGHT_BOUNDARY(n:1:-1))
end subroutine


subroutine test_vector_extrapolate ()
    !* - Check the separate vector extrapolation path against explicit weights.
    !  - Cover unequal query/grid sizes, empty calls, and reversed queries.

    type (search_cache) :: cache
    type (status_t) :: status
    integer :: ilb(size (X)), n
    real (PREC) :: weight(size (X)), queries(size (X))

    n = size (X)
    ! Extrapolation retains the boundary brackets but permits weights outside [0, 1].
    call interp_find (X, KNOTS, ilb, weight, &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, cache=cache, status=status)
    call check_vector (status, ilb, weight, LOWER, WEIGHT_EXTRAPOLATE)
    if (cache%i /= 4) error stop "Extrapolating search did not update caller cache"
    call interp_find (X, KNOTS, ilb, weight, &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
    call check_vector (status, ilb, weight, LOWER, WEIGHT_EXTRAPOLATE)

    ! Three queries on five knots exercise the short vector path.
    call interp_find (X(:3), KNOTS, ilb(:3), weight(:3), &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, cache=cache, status=status)
    call check_vector (status, ilb(:3), weight(:3), LOWER(:3), &
        WEIGHT_EXTRAPOLATE(:3))
    call interp_find (X(:3), KNOTS, ilb(:3), weight(:3), &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
    call check_vector (status, ilb(:3), weight(:3), LOWER(:3), &
        WEIGHT_EXTRAPOLATE(:3))

    ! Empty queries must also succeed in the specialized extrapolation path.
    call interp_find (X(:0), KNOTS, ilb(:0), weight(:0), &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, cache=cache, status=status)
    call check_vector (status, ilb(:0), weight(:0), LOWER(:0), &
        WEIGHT_EXTRAPOLATE(:0))
    call interp_find (X(:0), KNOTS, ilb(:0), weight(:0), &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
    call check_vector (status, ilb(:0), weight(:0), LOWER(:0), &
        WEIGHT_EXTRAPOLATE(:0))

    ! Reverse the full query sequence after the empty call.
    queries = X(n:1:-1)
    call interp_find (queries, KNOTS, ilb, weight, &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, cache=cache, status=status)
    call check_vector (status, ilb, weight, LOWER(n:1:-1), &
        WEIGHT_EXTRAPOLATE(n:1:-1))
    call interp_find (queries, KNOTS, ilb, weight, &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
    call check_vector (status, ilb, weight, LOWER(n:1:-1), &
        WEIGHT_EXTRAPOLATE(n:1:-1))
end subroutine


subroutine test_cache_indices ()
    !* - Recover from negative, zero, and oversized initial cache indices.
    !  - Check both boundary evaluation and extrapolation with scalar/vector APIs.

    integer, parameter :: INITIAL(*) = [-7, 0, 100]
    integer (NF_ENUM_KIND), parameter :: MODES(*) = &
        [NF_INTERP_EVAL_BOUNDARY, NF_INTERP_EVAL_EXTRAPOLATE]
    type (search_cache) :: scalar_cache, vector_cache
    type (status_t) :: status
    integer :: i, mode, ilb, bounds(1)
    real (PREC) :: weight, weights(1)

    do mode = 1, size (MODES)
        do i = 1, size (INITIAL)
            ! Invalid seeds must still locate 0 in [-1, 1], with lower weight 1/2.
            scalar_cache%i = INITIAL(i)
            vector_cache%i = INITIAL(i)
            call interp_find (0.0_PREC, KNOTS, ilb, weight, &
                ext=MODES(mode), cache=scalar_cache, status=status)
            call check_scalar (status, ilb, weight, 2, 0.5_PREC)
            call interp_find ([0.0_PREC], KNOTS, bounds, weights, &
                ext=MODES(mode), cache=vector_cache, status=status)
            call check_vector (status, bounds, weights, [2], [0.5_PREC])
        end do
        call interp_find (0.0_PREC, KNOTS, ilb, weight, &
            ext=MODES(mode), status=status)
        call check_scalar (status, ilb, weight, 2, 0.5_PREC)
        call interp_find ([0.0_PREC], KNOTS, bounds, weights, &
            ext=MODES(mode), status=status)
        call check_vector (status, bounds, weights, [2], [0.5_PREC])
    end do
end subroutine


subroutine test_cache_grid_resize ()
    !* - Warm scalar/vector caches before shrinking and growing a nonuniform grid.
    !  - Check boundary evaluation and extrapolation after each grid change.

    integer (NF_ENUM_KIND), parameter :: MODES(*) = &
        [NF_INTERP_EVAL_BOUNDARY, NF_INTERP_EVAL_EXTRAPOLATE]
    real (PREC), parameter :: SMALL_KNOTS(*) = [-4.0_PREC, 2.0_PREC, 10.0_PREC]
    type (search_cache) :: scalar_cache, vector_cache
    type (status_t) :: status
    integer :: mode, ilb, bounds(1)
    real (PREC) :: weight, weights(1)

    do mode = 1, size (MODES)
        ! Warm both caches in bracket 4, which does not exist on the smaller grid.
        call interp_find (6.5_PREC, KNOTS, ilb, weight, &
            ext=MODES(mode), cache=scalar_cache, status=status)
        call check_scalar (status, ilb, weight, 4, 0.75_PREC)
        call interp_find ([6.5_PREC], KNOTS, bounds, weights, &
            ext=MODES(mode), cache=vector_cache, status=status)
        call check_vector (status, bounds, weights, [4], [0.75_PREC])
        if (scalar_cache%i /= 4 .or. vector_cache%i /= 4) &
            error stop "Search caches were not warmed before grid shrink"

        ! Shrink to three knots: -1 is halfway between -4 and 2, in bracket 1.
        call interp_find (-1.0_PREC, SMALL_KNOTS, ilb, weight, &
            ext=MODES(mode), cache=scalar_cache, status=status)
        call check_scalar (status, ilb, weight, 1, 0.5_PREC)
        call interp_find ([-1.0_PREC], SMALL_KNOTS, bounds, weights, &
            ext=MODES(mode), cache=vector_cache, status=status)
        call check_vector (status, bounds, weights, [1], [0.5_PREC])
        call interp_find (-1.0_PREC, SMALL_KNOTS, ilb, weight, &
            ext=MODES(mode), status=status)
        call check_scalar (status, ilb, weight, 1, 0.5_PREC)
        call interp_find ([-1.0_PREC], SMALL_KNOTS, bounds, weights, &
            ext=MODES(mode), status=status)
        call check_vector (status, bounds, weights, [1], [0.5_PREC])

        ! Grow back to five knots: 0 again has bracket 2 and lower weight 1/2.
        call interp_find (0.0_PREC, KNOTS, ilb, weight, &
            ext=MODES(mode), cache=scalar_cache, status=status)
        call check_scalar (status, ilb, weight, 2, 0.5_PREC)
        call interp_find ([0.0_PREC], KNOTS, bounds, weights, &
            ext=MODES(mode), cache=vector_cache, status=status)
        call check_vector (status, bounds, weights, [2], [0.5_PREC])
    end do
end subroutine


subroutine check_scalar (status, ilb, weight, expected_ilb, expected_weight)
    !* Check scalar search outputs against independent expected values.

    type (status_t), intent(in) :: status
        !* Search exit status.
    integer, intent(in) :: ilb
        !* Returned lower bracket index.
    real (PREC), intent(in) :: weight
        !* Returned lower weight.
    integer, intent(in) :: expected_ilb
        !* Independently specified bracket index.
    real (PREC), intent(in) :: expected_weight
        !* Hand-computed lower weight.

    if (status /= NF_STATUS_OK) error stop "Scalar search status"
    if (ilb /= expected_ilb) error stop "Scalar search bracket"
    if (.not. (abs (weight - expected_weight) <= TOL)) &
        error stop "Scalar search lower weight"
end subroutine


subroutine check_vector (status, ilb, weight, expected_ilb, expected_weight)
    !* Check vector search outputs against independent expected values.

    type (status_t), intent(in) :: status
        !* Search exit status.
    integer, intent(in) :: ilb(:)
        !* Returned lower bracket indices.
    real (PREC), intent(in) :: weight(:)
        !* Returned lower weights.
    integer, intent(in) :: expected_ilb(:)
        !* Independently specified bracket indices.
    real (PREC), intent(in) :: expected_weight(:)
        !* Hand-computed lower weights.

    if (status /= NF_STATUS_OK) error stop "Vector search status"
    if (any (ilb /= expected_ilb)) error stop "Vector search brackets"
    if (.not. all (abs (weight - expected_weight) <= TOL)) &
        error stop "Vector search lower weights"
end subroutine

end program
