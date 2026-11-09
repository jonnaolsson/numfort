program test_interp_linear_cached
    !* Regression tests for linear interpolation with caller-provided caches.

    use, intrinsic :: iso_fortran_env, only: real64
    use numfort_interpolate
    implicit none

    integer, parameter :: PREC = real64
    real (PREC), parameter :: TOL = 1.0e-12_PREC
    real (PREC), parameter :: KNOTS(*) = &
        [-2.0_PREC, -1.0_PREC, 1.0_PREC, 5.0_PREC, 11.0_PREC]
    ! Tabulated f(x) = 2*x + 3 on an explicit nonuniform grid.
    real (PREC), parameter :: FP(*) = &
        [-1.0_PREC, 1.0_PREC, 5.0_PREC, 13.0_PREC, 25.0_PREC]
    real (PREC), parameter :: X(*) = [ &
        6.5_PREC, & ! Interior of [5, 11]: f(6.5) = 16.
        -3.0_PREC, & ! Lower exterior: f(-2) = -1 clamped, f(-3) = -3 extrapolated.
        -1.0_PREC, & ! Interior knot: f(-1) = 1.
        17.0_PREC, & ! Upper exterior: f(11) = 25 clamped, f(17) = 37 extrapolated.
        1.0_PREC, & ! Interior knot: f(1) = 5.
        0.0_PREC, & ! Midpoint of [-1, 1]: f(0) = 3.
        11.0_PREC, & ! Last knot: f(11) = 25.
        -2.0_PREC, & ! First knot: f(-2) = -1.
        6.5_PREC] ! Repeated interior point: f(6.5) = 16.
    real (PREC), parameter :: FX_BOUNDARY(*) = &
        [16.0_PREC, -1.0_PREC, 1.0_PREC, 25.0_PREC, 5.0_PREC, &
        3.0_PREC, 25.0_PREC, -1.0_PREC, 16.0_PREC]
    real (PREC), parameter :: FX_EXTRAPOLATE(*) = &
        [16.0_PREC, -3.0_PREC, 1.0_PREC, 37.0_PREC, 5.0_PREC, &
        3.0_PREC, 25.0_PREC, -1.0_PREC, 16.0_PREC]

    call test_scalar_boundary ()
    call test_vector_boundary ()
    call test_vector_extrapolate ()
    call test_cache_indices ()
    call test_cache_grid_resize ()

contains

subroutine test_scalar_boundary ()
    !* - Check scalar cached and uncached values against f(x) = 2*x + 3.
    !  - Reuse one cache through nonmonotone queries and boundary clamping.

    type (search_cache) :: cache
    type (status_t) :: status
    integer :: i
    real (PREC) :: fx

    do i = 1, size (X)
        ! Check each query against its explicit affine or clamped value above.
        call interp_linear (X(i), KNOTS, FP, fx, &
            ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
        call check_scalar (status, fx, FX_BOUNDARY(i))
        call interp_linear (X(i), KNOTS, FP, fx, &
            ext=NF_INTERP_EVAL_BOUNDARY, status=status)
        call check_scalar (status, fx, FX_BOUNDARY(i))
    end do
    if (cache%i /= 4) error stop "Scalar interpolation did not update caller cache"
end subroutine


subroutine test_vector_boundary ()
    !* - Check vector boundary evaluation with long, short, and empty queries.
    !  - Reuse the cache for reversed queries after an empty call.

    type (search_cache) :: cache
    type (status_t) :: status
    real (PREC) :: fx(size (X)), queries(size (X))
    integer :: n

    n = size (X)
    ! Nine nonmonotone queries on five knots, including repeated points.
    call interp_linear (X, KNOTS, FP, fx, &
        ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
    call check_vector (status, fx, FX_BOUNDARY)
    if (cache%i /= 4) error stop "Vector interpolation did not update caller cache"
    call interp_linear (X, KNOTS, FP, fx, &
        ext=NF_INTERP_EVAL_BOUNDARY, status=status)
    call check_vector (status, fx, FX_BOUNDARY)

    ! Fewer queries than knots: interior, lower exterior, and exact knot.
    call interp_linear (X(:3), KNOTS, FP, fx(:3), &
        ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
    call check_vector (status, fx(:3), FX_BOUNDARY(:3))
    call interp_linear (X(:3), KNOTS, FP, fx(:3), &
        ext=NF_INTERP_EVAL_BOUNDARY, status=status)
    call check_vector (status, fx(:3), FX_BOUNDARY(:3))

    ! Empty queries must succeed with a conformable empty output array.
    call interp_linear (X(:0), KNOTS, FP, fx(:0), &
        ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
    call check_vector (status, fx(:0), FX_BOUNDARY(:0))
    call interp_linear (X(:0), KNOTS, FP, fx(:0), &
        ext=NF_INTERP_EVAL_BOUNDARY, status=status)
    call check_vector (status, fx(:0), FX_BOUNDARY(:0))

    ! Reverse the full query sequence and reuse the cache after the empty call.
    queries = X(n:1:-1)
    call interp_linear (queries, KNOTS, FP, fx, &
        ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
    call check_vector (status, fx, FX_BOUNDARY(n:1:-1))
    call interp_linear (queries, KNOTS, FP, fx, &
        ext=NF_INTERP_EVAL_BOUNDARY, status=status)
    call check_vector (status, fx, FX_BOUNDARY(n:1:-1))
end subroutine


subroutine test_vector_extrapolate ()
    !* - Check vector extrapolation against the affine function beyond both ends.
    !  - Cover unequal query/grid sizes, empty calls, and reversed queries.

    type (search_cache) :: cache
    type (status_t) :: status
    real (PREC) :: fx(size (X)), queries(size (X))
    integer :: n

    n = size (X)
    ! Affine interpolation and extrapolation both reproduce f(x) = 2*x + 3.
    call interp_linear (X, KNOTS, FP, fx, &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, cache=cache, status=status)
    call check_vector (status, fx, FX_EXTRAPOLATE)
    if (cache%i /= 4) error stop "Extrapolation did not update caller cache"
    call interp_linear (X, KNOTS, FP, fx, &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
    call check_vector (status, fx, FX_EXTRAPOLATE)

    ! Three queries on five knots exercise the short vector path.
    call interp_linear (X(:3), KNOTS, FP, fx(:3), &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, cache=cache, status=status)
    call check_vector (status, fx(:3), FX_EXTRAPOLATE(:3))
    call interp_linear (X(:3), KNOTS, FP, fx(:3), &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
    call check_vector (status, fx(:3), FX_EXTRAPOLATE(:3))

    ! The specialized extrapolation path must also accept empty queries.
    call interp_linear (X(:0), KNOTS, FP, fx(:0), &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, cache=cache, status=status)
    call check_vector (status, fx(:0), FX_EXTRAPOLATE(:0))
    call interp_linear (X(:0), KNOTS, FP, fx(:0), &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
    call check_vector (status, fx(:0), FX_EXTRAPOLATE(:0))

    ! Reverse the query order and continue using the same caller cache.
    queries = X(n:1:-1)
    call interp_linear (queries, KNOTS, FP, fx, &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, cache=cache, status=status)
    call check_vector (status, fx, FX_EXTRAPOLATE(n:1:-1))
    call interp_linear (queries, KNOTS, FP, fx, &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
    call check_vector (status, fx, FX_EXTRAPOLATE(n:1:-1))
end subroutine


subroutine test_cache_indices ()
    !* - Recover from negative, zero, and oversized initial cache indices.
    !  - Check both boundary evaluation and extrapolation with scalar/vector APIs.

    integer, parameter :: INITIAL(*) = [-7, 0, 100]
    integer (NF_ENUM_KIND), parameter :: MODES(*) = &
        [NF_INTERP_EVAL_BOUNDARY, NF_INTERP_EVAL_EXTRAPOLATE]
    type (search_cache) :: scalar_cache, vector_cache
    type (status_t) :: status
    integer :: i, mode
    real (PREC) :: fx, values(1)

    do mode = 1, size (MODES)
        do i = 1, size (INITIAL)
            ! Invalid seeds must still produce the interior value f(0) = 3.
            scalar_cache%i = INITIAL(i)
            vector_cache%i = INITIAL(i)
            call interp_linear (0.0_PREC, KNOTS, FP, fx, &
                ext=MODES(mode), cache=scalar_cache, status=status)
            call check_scalar (status, fx, 3.0_PREC)
            call interp_linear ([0.0_PREC], KNOTS, FP, values, &
                ext=MODES(mode), cache=vector_cache, status=status)
            call check_vector (status, values, [3.0_PREC])
        end do
        call interp_linear (0.0_PREC, KNOTS, FP, fx, &
            ext=MODES(mode), status=status)
        call check_scalar (status, fx, 3.0_PREC)
        call interp_linear ([0.0_PREC], KNOTS, FP, values, &
            ext=MODES(mode), status=status)
        call check_vector (status, values, [3.0_PREC])
    end do
end subroutine


subroutine test_cache_grid_resize ()
    !* - Warm scalar/vector caches before shrinking and growing a nonuniform grid.
    !  - Check boundary evaluation and extrapolation after each grid change.

    integer (NF_ENUM_KIND), parameter :: MODES(*) = &
        [NF_INTERP_EVAL_BOUNDARY, NF_INTERP_EVAL_EXTRAPOLATE]
    real (PREC), parameter :: SMALL_KNOTS(*) = [-4.0_PREC, 2.0_PREC, 10.0_PREC]
    ! The same affine function sampled on the smaller nonuniform grid.
    real (PREC), parameter :: SMALL_FP(*) = [-5.0_PREC, 7.0_PREC, 23.0_PREC]
    type (search_cache) :: scalar_cache, vector_cache
    type (status_t) :: status
    integer :: mode
    real (PREC) :: fx, values(1)

    do mode = 1, size (MODES)
        ! f(6.5) = 16 warms bracket 4, which does not exist on the smaller grid.
        call interp_linear (6.5_PREC, KNOTS, FP, fx, &
            ext=MODES(mode), cache=scalar_cache, status=status)
        call check_scalar (status, fx, 16.0_PREC)
        call interp_linear ([6.5_PREC], KNOTS, FP, values, &
            ext=MODES(mode), cache=vector_cache, status=status)
        call check_vector (status, values, [16.0_PREC])
        if (scalar_cache%i /= 4 .or. vector_cache%i /= 4) &
            error stop "Interpolation caches were not warmed before grid shrink"

        ! Shrink to three knots: -1 is interior to [-4, 2], and f(-1) = 1.
        call interp_linear (-1.0_PREC, SMALL_KNOTS, SMALL_FP, fx, &
            ext=MODES(mode), cache=scalar_cache, status=status)
        call check_scalar (status, fx, 1.0_PREC)
        call interp_linear ([-1.0_PREC], SMALL_KNOTS, SMALL_FP, values, &
            ext=MODES(mode), cache=vector_cache, status=status)
        call check_vector (status, values, [1.0_PREC])
        call interp_linear (-1.0_PREC, SMALL_KNOTS, SMALL_FP, fx, &
            ext=MODES(mode), status=status)
        call check_scalar (status, fx, 1.0_PREC)
        call interp_linear ([-1.0_PREC], SMALL_KNOTS, SMALL_FP, values, &
            ext=MODES(mode), status=status)
        call check_vector (status, values, [1.0_PREC])

        ! Grow back to five knots: 0 is now in [-1, 1], and f(0) = 3.
        call interp_linear (0.0_PREC, KNOTS, FP, fx, &
            ext=MODES(mode), cache=scalar_cache, status=status)
        call check_scalar (status, fx, 3.0_PREC)
        call interp_linear ([0.0_PREC], KNOTS, FP, values, &
            ext=MODES(mode), cache=vector_cache, status=status)
        call check_vector (status, values, [3.0_PREC])
    end do
end subroutine


subroutine check_scalar (status, fx, expected)
    !* Check a scalar interpolated value against the affine reference.

    type (status_t), intent(in) :: status
        !* Interpolation exit status.
    real (PREC), intent(in) :: fx
        !* Returned interpolated value.
    real (PREC), intent(in) :: expected
        !* Hand-computed affine or clamped value.

    if (status /= NF_STATUS_OK) error stop "Scalar interpolation status"
    if (.not. (abs (fx - expected) <= TOL)) error stop "Scalar interpolation value"
end subroutine


subroutine check_vector (status, fx, expected)
    !* Check vector interpolated values against the affine reference.

    type (status_t), intent(in) :: status
        !* Interpolation exit status.
    real (PREC), intent(in) :: fx(:)
        !* Returned interpolated values.
    real (PREC), intent(in) :: expected(:)
        !* Hand-computed affine or clamped values.

    if (status /= NF_STATUS_OK) error stop "Vector interpolation status"
    if (.not. all (abs (fx - expected) <= TOL)) error stop "Vector interpolation values"
end subroutine

end program
