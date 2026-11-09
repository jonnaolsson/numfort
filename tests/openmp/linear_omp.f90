program test_interp_linear

    use, intrinsic :: iso_fortran_env
    use omp_lib, only: omp_get_thread_num, omp_set_dynamic
    use numfort_arrays
    use numfort_interpolate
    implicit none

    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    call test_extrap ()
    call test_cached_find ()
    call test_cached_linear ()

contains

!>  Tests for linear extrapolation when x-coordinate of interpolated value
!   is outside of data points boundaries.
subroutine test_extrap ()

    integer, parameter :: n = int(1e5)
    real (PREC), dimension(10) :: xp, fp
    real (PREC), dimension(:), allocatable :: fx_hat, x, fx
    integer :: i

    ! function f(x) = 2x + 1 on [10, 20]
    call linspace (xp, 1.0d1, 2.0d1)
    call func1 (xp, fp)

    allocate (x(n), fx(n), fx_hat(n))
    ! evaluate on [0,100], ie using extrapolation
    call linspace (x, 0.0d0, 1.0d2)
    fx = 0.0

    !$omp parallel default(shared) private(i)
    !$omp do
    do i = 1, size(x)
        call interp_linear (x(i), xp, fp, fx_hat(i), ext=NF_INTERP_EVAL_EXTRAPOLATE)
    end do
    !$omp end do

    !$omp end parallel

    ! compute true values f(x)
    call func1 (x, fx)

    if (.not. (all (abs (fx - fx_hat) < 1d-12))) &
        error stop "Extrapolation, scalar argument"

end subroutine


subroutine test_cached_find ()
    !* - Check scalar boundary searches with one cache per OpenMP worker.
    !  - Repeatedly search nonmonotone points on distinct nonuniform grids.

    real (PREC), parameter :: BASE_KNOTS(*) = &
        [-2.0_PREC, -1.0_PREC, 1.0_PREC, 5.0_PREC, 11.0_PREC]
    real (PREC), parameter :: BASE_X(*) = [ &
        6.5_PREC, & ! Interior: (11 - 6.5) / (11 - 5) = 3/4 on the lower knot.
        -3.0_PREC, & ! Lower exterior: clamp to the first knot with weight 1.
        -1.0_PREC, & ! Interior knot: bracket 2, lower weight 1.
        17.0_PREC, & ! Upper exterior: clamp to the last knot with weight 0.
        1.0_PREC, & ! Interior knot: bracket 3, lower weight 1.
        0.0_PREC, & ! Midpoint of [-1, 1]: lower weight 1/2.
        11.0_PREC, & ! Last knot: bracket 4, lower weight 0.
        -2.0_PREC, & ! First knot: bracket 1, lower weight 1.
        6.5_PREC] ! Repeated interior point, lower weight 3/4.
    integer, parameter :: LOWER(*) = [4, 1, 2, 4, 3, 2, 4, 1, 4]
    real (PREC), parameter :: WEIGHTS(*) = &
        [0.75_PREC, 1.0_PREC, 1.0_PREC, 0.0_PREC, 1.0_PREC, &
        0.5_PREC, 0.0_PREC, 1.0_PREC, 0.75_PREC]
    type (search_cache) :: cache
    type (status_t) :: status
    real (PREC) :: knots(5), queries(9), shift, weight
    integer :: i, repeat, ilb

    call omp_set_dynamic (.false.)
    !$omp parallel default(none) num_threads(4) &
    !$omp private(cache, status, knots, queries, shift, weight, i, repeat, ilb)
    ! Translation gives each worker its own grid without changing brackets/weights.
    shift = real (omp_get_thread_num (), PREC)
    knots = BASE_KNOTS + shift
    queries = BASE_X + shift
    cache%i = 100
    do repeat = 1, 32
        do i = 1, size (queries)
            call interp_find (queries(i), knots, ilb, weight, &
                ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
            if (status /= NF_STATUS_OK) error stop "OpenMP cached search status"
            if (ilb /= LOWER(i)) error stop "OpenMP cached search bracket"
            if (.not. (abs (weight - WEIGHTS(i)) <= 1.0e-12_PREC)) &
                error stop "OpenMP cached search weight"
            call interp_find (queries(i), knots, ilb, weight, &
                ext=NF_INTERP_EVAL_BOUNDARY, status=status)
            if (status /= NF_STATUS_OK) error stop "OpenMP uncached search status"
            if (ilb /= LOWER(i)) error stop "OpenMP uncached search bracket"
            if (.not. (abs (weight - WEIGHTS(i)) <= 1.0e-12_PREC)) &
                error stop "OpenMP uncached search weight"
        end do
    end do
    !$omp end parallel
end subroutine


subroutine test_cached_linear ()
    !* - Check cached vector interpolation with one cache per OpenMP worker.
    !  - Verify boundary and extrapolated values against f(x) = 2*x + 3.

    real (PREC), parameter :: BASE_KNOTS(*) = &
        [-2.0_PREC, -1.0_PREC, 1.0_PREC, 5.0_PREC, 11.0_PREC]
    real (PREC), parameter :: BASE_X(*) = [ &
        6.5_PREC, & ! Interior: f(6.5) = 16.
        -3.0_PREC, & ! Lower exterior: -1 clamped, -3 extrapolated.
        -1.0_PREC, & ! Interior knot: f(-1) = 1.
        17.0_PREC, & ! Upper exterior: 25 clamped, 37 extrapolated.
        1.0_PREC, & ! Interior knot: f(1) = 5.
        0.0_PREC, & ! Interior midpoint: f(0) = 3.
        11.0_PREC, & ! Last knot: f(11) = 25.
        -2.0_PREC, & ! First knot: f(-2) = -1.
        6.5_PREC] ! Repeat the interior value 16.
    real (PREC), parameter :: BOUNDARY(*) = &
        [16.0_PREC, -1.0_PREC, 1.0_PREC, 25.0_PREC, 5.0_PREC, &
        3.0_PREC, 25.0_PREC, -1.0_PREC, 16.0_PREC]
    real (PREC), parameter :: EXTRAPOLATE(*) = &
        [16.0_PREC, -3.0_PREC, 1.0_PREC, 37.0_PREC, 5.0_PREC, &
        3.0_PREC, 25.0_PREC, -1.0_PREC, 16.0_PREC]
    type (search_cache) :: cache
    type (status_t) :: status
    real (PREC) :: knots(5), fp(5), queries(9), fx(9), expected(9), shift
    integer :: repeat

    call omp_set_dynamic (.false.)
    !$omp parallel default(none) num_threads(4) &
    !$omp private(cache, status, knots, fp, queries, fx, expected, shift, repeat)
    shift = real (omp_get_thread_num (), PREC)
    knots = BASE_KNOTS + shift
    queries = BASE_X + shift
    fp = 2.0_PREC * knots + 3.0_PREC
    cache%i = -7
    do repeat = 1, 32
        ! f(x + shift) = f(x) + 2*shift, including translated boundary values.
        expected = BOUNDARY + 2.0_PREC * shift
        call interp_linear (queries, knots, fp, fx, &
            ext=NF_INTERP_EVAL_BOUNDARY, cache=cache, status=status)
        if (status /= NF_STATUS_OK) error stop "OpenMP cached boundary status"
        if (.not. all (abs (fx - expected) <= 1.0e-12_PREC)) &
            error stop "OpenMP cached boundary values"
        call interp_linear (queries, knots, fp, fx, &
            ext=NF_INTERP_EVAL_BOUNDARY, status=status)
        if (status /= NF_STATUS_OK) error stop "OpenMP uncached boundary status"
        if (.not. all (abs (fx - expected) <= 1.0e-12_PREC)) &
            error stop "OpenMP uncached boundary values"

        ! Reuse each worker's cache for extrapolation, with the same affine shift.
        expected = EXTRAPOLATE + 2.0_PREC * shift
        call interp_linear (queries, knots, fp, fx, &
            ext=NF_INTERP_EVAL_EXTRAPOLATE, cache=cache, status=status)
        if (status /= NF_STATUS_OK) error stop "OpenMP cached extrapolation status"
        if (.not. all (abs (fx - expected) <= 1.0e-12_PREC)) &
            error stop "OpenMP cached extrapolation values"
        call interp_linear (queries, knots, fp, fx, &
            ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
        if (status /= NF_STATUS_OK) error stop "OpenMP uncached extrapolation status"
        if (.not. all (abs (fx - expected) <= 1.0e-12_PREC)) &
            error stop "OpenMP uncached extrapolation values"
    end do
    !$omp end parallel
end subroutine


pure subroutine func1 (x, fx)
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: fx

    fx = 2.0d0 * x + 1.0d0
end subroutine

end program
