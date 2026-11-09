program test_linear_optional
    !* Test optional boundary constants in both interpolation precisions.

    use, intrinsic :: iso_fortran_env, only: TEST_PREC, int32
    use numfort_common_status
    use numfort_interpolate

    implicit none

    integer, parameter :: PREC = TEST_PREC
    real (PREC), parameter :: TOL = 16.0_PREC * epsilon (1.0_PREC)

    ! The affine function f(x) = 2*x + 1 on a nonuniform grid.
    real (PREC), parameter :: XP(3) = [0.0_PREC, 1.0_PREC, 3.0_PREC]
    real (PREC), parameter :: FP(3) = [1.0_PREC, 3.0_PREC, 7.0_PREC]
    ! Unsorted queries cover both exterior sides, endpoints, and an interior point.
    real (PREC), parameter :: X(5) = &
        [4.0_PREC, 0.0_PREC, 0.5_PREC, -1.0_PREC, 3.0_PREC]
    ! Boundary evaluation clamps f(x) to f(0) = 1 and f(3) = 7.
    real (PREC), parameter :: BOUNDARY(5) = &
        [7.0_PREC, 1.0_PREC, 2.0_PREC, 1.0_PREC, 7.0_PREC]
    ! Linear extrapolation of f(x) = 2*x + 1 is exact at every query.
    real (PREC), parameter :: EXTRAPOLATED(5) = &
        [9.0_PREC, 1.0_PREC, 2.0_PREC, -1.0_PREC, 7.0_PREC]
    ! CONST substitutes LEFT and RIGHT only for strictly exterior queries.
    real (PREC), parameter :: LEFT = -11.0_PREC, RIGHT = 13.0_PREC
    real (PREC), parameter :: CONSTANT(5) = &
        [RIGHT, 1.0_PREC, 2.0_PREC, LEFT, 7.0_PREC]

    call test_linear_omitted (NF_INTERP_EVAL_BOUNDARY, BOUNDARY)
    call test_linear_omitted (NF_INTERP_EVAL_EXTRAPOLATE, EXTRAPOLATED)
    call test_eval_omitted (NF_INTERP_EVAL_BOUNDARY, BOUNDARY)
    call test_eval_omitted (NF_INTERP_EVAL_EXTRAPOLATE, EXTRAPOLATED)
    call test_default_omitted ()
    call test_linear_const ()
    call test_eval_const ()
    call test_eval_const_missing ()

contains

subroutine test_linear_omitted (ext, expected)
    !* - Omit both constants in boundary and extrapolation modes.
    !  - Check scalar/vector calls with and without a search cache.

    integer (NF_ENUM_KIND), intent(in) :: ext
        !* Extrapolation mode for this fixture.
    real (PREC), intent(in) :: expected(:)
        !* Independent affine or clamped values at the query points.

    real (PREC) :: fx(size(X))
    type (search_cache) :: cache
    type (status_t) :: status
    integer :: i

    ! Scalar calls without a cache omit LEFT and RIGHT at every query.
    do i = 1, size(X)
        call interp_linear (X(i), XP, FP, fx(i), ext=ext, status=status)
        if (status /= NF_STATUS_OK) error stop "Scalar omission status"
    end do
    if (.not. all (abs (fx - expected) <= TOL)) &
        error stop "Scalar omission values"

    ! Scalar calls reuse a cache across unsorted queries with both constants absent.
    do i = 1, size(X)
        call interp_linear (X(i), XP, FP, fx(i), ext=ext, cache=cache, status=status)
        if (status /= NF_STATUS_OK) error stop "Cached scalar omission status"
    end do
    if (.not. all (abs (fx - expected) <= TOL)) &
        error stop "Cached scalar omission values"

    ! A vector call without a cache omits both constants.
    call interp_linear (X, XP, FP, fx, ext=ext, status=status)
    if (status /= NF_STATUS_OK) error stop "Vector omission status"
    if (.not. all (abs (fx - expected) <= TOL)) &
        error stop "Vector omission values"

    ! A vector call reuses the scalar cache with both constants absent.
    call interp_linear (X, XP, FP, fx, ext=ext, cache=cache, status=status)
    if (status /= NF_STATUS_OK) error stop "Cached vector omission status"
    if (.not. all (abs (fx - expected) <= TOL)) &
        error stop "Cached vector omission values"

end subroutine


subroutine test_eval_omitted (ext, expected)
    !* - Omit both constants when evaluating precomputed brackets and weights.
    !  - Check the scalar and vector evaluator interfaces.

    integer (NF_ENUM_KIND), intent(in) :: ext
        !* Extrapolation mode used for both searching and evaluation.
    real (PREC), intent(in) :: expected(:)
        !* Independent affine or clamped values at the query points.

    integer (int32) :: ilbound(size(X))
    real (PREC) :: weight(size(X)), fx(size(X))
    type (status_t) :: status
    integer :: i

    call interp_find (X, XP, ilbound, weight, ext=ext, status=status)
    if (status /= NF_STATUS_OK) error stop "Omission fixture search status"

    ! Scalar evaluation omits both constants at exterior and interior queries.
    do i = 1, size(X)
        call interp_linear_eval (ilbound(i), weight(i), FP, fx(i), &
            ext=ext, status=status)
        if (status /= NF_STATUS_OK) error stop "Scalar eval omission status"
    end do
    if (.not. all (abs (fx - expected) <= TOL)) &
        error stop "Scalar eval omission values"

    ! Vector evaluation omits both constants for the same precomputed queries.
    call interp_linear_eval (ilbound, weight, FP, fx, ext=ext, status=status)
    if (status /= NF_STATUS_OK) error stop "Vector eval omission status"
    if (.not. all (abs (fx - expected) <= TOL)) &
        error stop "Vector eval omission values"

end subroutine


subroutine test_default_omitted ()
    !* - Default extrapolation accepts absent constants and status.
    !  - Cover interpolation and evaluation, each with scalar/vector input.

    integer (int32) :: ilbound(size(X))
    real (PREC) :: weight(size(X)), fx(size(X))
    integer :: i

    ! Scalar interpolation omits EXT, LEFT, RIGHT, and STATUS.
    do i = 1, size(X)
        call interp_linear (X(i), XP, FP, fx(i))
    end do
    if (.not. all (abs (fx - EXTRAPOLATED) <= TOL)) &
        error stop "Default scalar interpolation"

    ! Vector interpolation omits all optional arguments.
    call interp_linear (X, XP, FP, fx)
    if (.not. all (abs (fx - EXTRAPOLATED) <= TOL)) &
        error stop "Default vector interpolation"

    call interp_find (X, XP, ilbound, weight)

    ! Scalar evaluation omits EXT, LEFT, RIGHT, and STATUS.
    do i = 1, size(X)
        call interp_linear_eval (ilbound(i), weight(i), FP, fx(i))
    end do
    if (.not. all (abs (fx - EXTRAPOLATED) <= TOL)) &
        error stop "Default scalar evaluation"

    ! Vector evaluation also accepts an absent STATUS argument.
    call interp_linear_eval (ilbound, weight, FP, fx)
    if (.not. all (abs (fx - EXTRAPOLATED) <= TOL)) &
        error stop "Default vector evaluation"

end subroutine


subroutine test_linear_const ()
    !* - Preserve explicit constants in CONST mode.
    !  - Check scalar/vector interpolation with and without a cache.

    real (PREC) :: fx(size(X))
    type (search_cache) :: cache
    type (status_t) :: status
    integer :: i

    ! Scalar CONST calls use distinct constants on the left and right.
    do i = 1, size(X)
        call interp_linear (X(i), XP, FP, fx(i), ext=NF_INTERP_EVAL_CONST, &
            left=LEFT, right=RIGHT, status=status)
        if (status /= NF_STATUS_OK) error stop "Scalar CONST status"
    end do
    if (.not. all (abs (fx - CONSTANT) <= TOL)) error stop "Scalar CONST values"

    ! Cached scalar CONST calls preserve endpoint and interior values.
    do i = 1, size(X)
        call interp_linear (X(i), XP, FP, fx(i), ext=NF_INTERP_EVAL_CONST, &
            left=LEFT, right=RIGHT, cache=cache, status=status)
        if (status /= NF_STATUS_OK) error stop "Cached scalar CONST status"
    end do
    if (.not. all (abs (fx - CONSTANT) <= TOL)) &
        error stop "Cached scalar CONST values"

    ! Vector CONST calls use both explicit constants without a cache.
    call interp_linear (X, XP, FP, fx, ext=NF_INTERP_EVAL_CONST, &
        left=LEFT, right=RIGHT, status=status)
    if (status /= NF_STATUS_OK) error stop "Vector CONST status"
    if (.not. all (abs (fx - CONSTANT) <= TOL)) error stop "Vector CONST values"

    ! Cached vector CONST calls reuse the search state from scalar calls.
    call interp_linear (X, XP, FP, fx, ext=NF_INTERP_EVAL_CONST, &
        left=LEFT, right=RIGHT, cache=cache, status=status)
    if (status /= NF_STATUS_OK) error stop "Cached vector CONST status"
    if (.not. all (abs (fx - CONSTANT) <= TOL)) &
        error stop "Cached vector CONST values"

end subroutine


subroutine test_eval_const ()
    !* - Preserve explicit constants when evaluating precomputed brackets.
    !  - Cover both scalar and vector CONST evaluation.

    integer (int32) :: ilbound(size(X))
    real (PREC) :: weight(size(X)), fx(size(X))
    type (status_t) :: status
    integer :: i

    call interp_find (X, XP, ilbound, weight, ext=NF_INTERP_EVAL_CONST, status=status)
    if (status /= NF_STATUS_OK) error stop "CONST fixture search status"

    ! Scalar evaluation uses LEFT and RIGHT only at exterior queries.
    do i = 1, size(X)
        call interp_linear_eval (ilbound(i), weight(i), FP, fx(i), &
            ext=NF_INTERP_EVAL_CONST, left=LEFT, right=RIGHT, status=status)
        if (status /= NF_STATUS_OK) error stop "Scalar eval CONST status"
    end do
    if (.not. all (abs (fx - CONSTANT) <= TOL)) &
        error stop "Scalar eval CONST values"

    ! Vector evaluation uses the same distinct constants.
    call interp_linear_eval (ilbound, weight, FP, fx, ext=NF_INTERP_EVAL_CONST, &
        left=LEFT, right=RIGHT, status=status)
    if (status /= NF_STATUS_OK) error stop "Vector eval CONST status"
    if (.not. all (abs (fx - CONSTANT) <= TOL)) &
        error stop "Vector eval CONST values"

end subroutine


subroutine test_eval_const_missing ()
    !* - Reject absent CONST arguments before evaluating the interpolant.
    !  - Require both constants even for an interior query.

    integer (int32), parameter :: ILBOUND(1) = [1_int32]
    real (PREC), parameter :: WEIGHT(1) = [0.5_PREC]
    real (PREC) :: fx(1)
    type (status_t) :: status

    ! Scalar CONST evaluation rejects omission of both constants.
    call interp_linear_eval (ILBOUND(1), WEIGHT(1), FP, fx(1), &
        ext=NF_INTERP_EVAL_CONST, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Scalar CONST missing both"

    ! Scalar CONST evaluation rejects omission of RIGHT alone.
    call interp_linear_eval (ILBOUND(1), WEIGHT(1), FP, fx(1), &
        ext=NF_INTERP_EVAL_CONST, left=LEFT, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Scalar CONST missing right"

    ! Scalar CONST evaluation rejects omission of LEFT alone.
    call interp_linear_eval (ILBOUND(1), WEIGHT(1), FP, fx(1), &
        ext=NF_INTERP_EVAL_CONST, right=RIGHT, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Scalar CONST missing left"

    ! Vector CONST evaluation rejects omission of both constants.
    call interp_linear_eval (ILBOUND, WEIGHT, FP, fx, &
        ext=NF_INTERP_EVAL_CONST, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Vector CONST missing both"

    ! Vector CONST evaluation rejects omission of RIGHT alone.
    call interp_linear_eval (ILBOUND, WEIGHT, FP, fx, &
        ext=NF_INTERP_EVAL_CONST, left=LEFT, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Vector CONST missing right"

    ! Vector CONST evaluation rejects omission of LEFT alone.
    call interp_linear_eval (ILBOUND, WEIGHT, FP, fx, &
        ext=NF_INTERP_EVAL_CONST, right=RIGHT, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Vector CONST missing left"

end subroutine

end program
