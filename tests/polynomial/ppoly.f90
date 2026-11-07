

program test_ppoly
    !*  Unit tests for function approximation with piecewise polynomials

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common
    use numfort_common_workspace, workspace => workspace_real64
    use numfort_common_testing
    use numfort_polynomial
    use numfort_interpolate

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_bernstein_ppolyfit_input ()
    call test_bernstein_ppolyfit ()

    call test_power_ppolyval_input ()
    call test_power_ppolyval_eval_input ()
    call test_power_ppolyval_cubic ()

end subroutine

subroutine test_bernstein_ppolyfit ()
    !*  Unit tests for fitting and evaluating piecewise polynomials wrt. the
    !   Bernstein basis.

    real (PREC), dimension(:), allocatable :: knots, coefs, coefs1, xx, xp, yp, yy
    real (PREC), dimension(:), allocatable :: fx, fpx
    real (PREC), dimension(:,:), allocatable :: ydat
    integer :: k, ncoefs, nknots, nx, n, i, ilb
    real (PREC) :: xlb, xub, left, right, ylb, yub, wgt
    real (PREC), dimension(:), allocatable :: weight
    integer, dimension(:), allocatable :: ilbound
    type (ppoly_bernstein) :: pp, pp1
    logical :: status_ok, values_ok
    type (status_t) :: status

    n = 11
    xlb = 1.0d-8
    xub = 100
    k = 3
    allocate (xx(n), ydat(2, n))
    allocate (fx(n), fpx(n))

    ! Compute function values and first derivatives
    call linspace (xx, xlb, xub)
    call fcn2 (xx, fx, fpx)
    ydat(1,:) = fx
    ydat(2,:) = fpx

    ncoefs = ppoly_get_ncoefs (pp, n, k)
    nknots = ppoly_get_nknots (pp, n, k)
    allocate (knots(nknots), coefs(ncoefs))

    ! Fit piecewise polynomial to 0th and 1st derivative data
    call ppolyfit (pp, xx, ydat, k, knots, coefs, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "1d API: Cubic: fit piecewise polynomial"

    allocate (yy(n))
    call ppolyval (pp, knots, coefs, xx, yy, status=status)
    values_ok = all_close (yy, fx, atol=1.0e-8_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: Cubic: evaluate at original knots"

    ! 1d API: Test separate find/eval routines
    allocate (ilbound(n), weight(n))

    status = NF_STATUS_UNDEFINED
    call interp_find (xx, knots, ilbound, weight, status=status)
    status_ok = (status == NF_STATUS_OK)

    status = NF_STATUS_UNDEFINED
    yy(:) = 0.0
    call ppolyval_eval (pp, ilbound, weight, coefs, yy, status=status)
    status_ok = status_ok .and. (status == NF_STATUS_OK)
    values_ok = all_close (yy, fx, atol=1.0e-8_PREC)

    if (.not. (status_ok .and. values_ok)) &
        error stop "1d API: Cubic: find/eval at original knots"

    deallocate (yy, ilbound, weight)

    ! Evaluate using scalar interface
    allocate (yy(n))
    status_ok = .true.
    do i = 1, size(xx)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, xx(i), yy(i), status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, fx, atol=1.0e-8_PREC)
    if (.not. (status_ok .and. values_ok)) &
        error stop "0d API: Cubic: evaluate at original knots"

    ! 0d API: Separate find/eval routines
    status_ok = .true.
    do i = 1, size(xx)
        status = NF_STATUS_UNDEFINED
        call interp_find (xx(i), knots, ilb, wgt, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)

        call ppolyval_eval (pp, ilb, wgt, coefs, yy(i), status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do

    values_ok = all_close (yy, fx, atol=1.0e-8_PREC)
    if (.not. (status_ok .and. values_ok)) &
        error stop "0d API: Cubic: find/eval at original knots"

    deallocate (yy)

    ! === Evaluate at different points than knots ===
    nx = 107
    allocate (yp(nx), xp(nx), yy(nx))
    call linspace (xp, xlb, xub)
    ! True function values
    call fcn2 (xp, yp)
    call ppolyval (pp, knots, coefs, xp, yy, status=status)
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: Cubic: Evaluate at points other than knots"

    ! 1d API: separate find/eval routines
    allocate (ilbound(nx), weight(nx))

    status = NF_STATUS_UNDEFINED
    call interp_find (xp, knots, ilbound, weight, status=status)
    status_ok = (status == NF_STATUS_OK)

    status = NF_STATUS_UNDEFINED
    call ppolyval_eval (pp, ilbound, weight, coefs, yy, status=status)
    status_ok = status_ok .and. (status == NF_STATUS_OK)
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)

    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: Cubic: find/eval at points other than knots"

    deallocate (yy, ilbound, weight)

    ! Evaluate using scalar interface
    allocate (yy(nx), source=0.0_PREC)
    status_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, xp(i), yy(i), status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    if (.not. (status_ok .and. values_ok)) &
        error stop "0d API: Cubic: Evaluate at points other than knots"

    ! 0d API: separate find/eval routines
    status_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call interp_find (xp(i), knots, ilb, wgt, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)

        status = NF_STATUS_UNDEFINED
        call ppolyval_eval (pp, ilb, wgt, coefs, yy(i), status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do

    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    if (.not. (status_ok .and. values_ok)) &
        error stop "0d API: Cubic: find/eval at points other than knots"

    deallocate (yp, xp, yy)

    ! === Extrapolation ===
    nx = 17
    allocate (yp(nx), xp(nx), yy(nx))
    call linspace (xp, xlb-1.0, xub+1.0)
    ! Some additional non-interior points
    xp(5) = xlb - 10.0_PREC
    xp(15) = xub + 10.0_PREC
    ! Correct function values
    call fcn2 (xp, yp)
    call ppolyval (pp, knots, coefs, xp, yy, ext=NF_INTERP_EVAL_EXTRAPOLATE, &
        status=status)
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "0d API: Cubic: Evaluate at non-interior points"
    deallocate (yy)

    ! Scalar interface
    allocate (yy(nx), source=0.0_PREC)
    status_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, xp(i), yy(i), &
            ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    if (.not. (status_ok .and. values_ok)) &
        error stop "0d API: Cubic: Evaluate at non-interior points"

    ! === Replace non-interior points with some constant values ===
    left = -10000.0
    right = 10000.0
    where (xp < xlb)
        yp = left
    else where (xp > xub)
        yp = right
    end where
    call ppolyval (pp, knots, coefs, xp, yy, ext=NF_INTERP_EVAL_CONST, &
        left=left, right=right, status=status)
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "0d API: Cubic: Evaluate, assign const to non-interior points"
    deallocate (yy)

    ! Scalar interface
    allocate (yy(nx))
    status_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, xp(i), yy(i), ext=NF_INTERP_EVAL_CONST, &
            left=left, right=right, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    if (.not. (status_ok .and. values_ok)) &
        error stop "0d API: Cubic: Evaluate, assign const to non-interior points"

    ! === Raise error on encountering non-interior points ===
    call ppolyval (pp, knots, coefs, xp, yy, ext=NF_INTERP_EVAL_ERROR, &
        status=status)
    if (.not. (status == NF_STATUS_BOUNDS_ERROR)) &
        error stop "1d API: Cubic: Evaluate, raise error on non-interior points"

    ! Replace non-interior points with boundary values
    call fcn2 (xp, yp)
    call fcn2 (xlb, ylb)
    call fcn2 (xub, yub)
    where (xp < xlb)
        yp = ylb
    else where (xp > xub)
        yp = yub
    end where

    call ppolyval (pp, knots, coefs, xp, yy, ext=NF_INTERP_EVAL_BOUNDARY, &
        status=status)
    values_ok = all_close (yp, yy, atol=1.0e-8_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: Cubic: Evaluate, assign const to non-interior points"

    deallocate (xp, yy, yp)

    ! === Test evaluating first derivative ===
    ncoefs = ppoly_get_ncoefs (pp, k=k-1)
    allocate (coefs1(ncoefs))
    call ppolyder (pp, knots, coefs, pp1, coefs1, m=1, status=status)

    nx = 13
    status = NF_STATUS_UNDEFINED
    allocate (xp(nx), yy(nx), yp(nx))
    call linspace (xp, xlb, xub)
    ! True first derivative
    call fcn2 (xp, fpx=yp)
    call ppolyval (pp1, knots, coefs1, xp, yy, status=status)
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: Cubic: Evaluate first derivative"
    deallocate (yy)

    ! Scalar interface
    allocate (yy(nx), source=0.0_PREC)
    status_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp1, knots, coefs1, xp(i), yy(i), status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    if (.not. (status_ok .and. values_ok)) &
        error stop "0d API: Cubic: Evaluate first derivative"

    ! === Extrapolate first derivative ===
    status = NF_STATUS_UNDEFINED
    call linspace (xp, xlb-10.0, xub+10.0)
    call fcn2 (xp, fpx=yp)
    ! Note: extrapolates by default
    call ppolyval (pp1, knots, coefs1, xp, yy, status=status)
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: Cubic: Evaluate first derivative at non-interior points"
    deallocate (yy)

    ! Scalar interface
    allocate (yy(nx), source=0.0_PREC)
    status_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp1, knots, coefs1, xp(i), yy(i), status=status)
        status_ok = status_ok .and. (status == NF_STATUS_OK)
    end do
    values_ok = all_close (yy, yp, atol=1.0e-8_PREC)
    if (.not. (status_ok .and. values_ok)) &
        error stop "0d API: Cubic: Evaluate first derivaitve at non-interior points"

    deallocate (xp, yp, yy)

end subroutine

subroutine test_bernstein_ppolyfit_input ()
    !*  Unit tests for input checking for Bernstein basis fitting routines

    real (PREC), dimension(:), allocatable :: knots, coefs
    real (PREC), dimension(3) :: x
    real (PREC), dimension(2,3) :: y
    integer :: n, k, ncoefs, nknots
    type (ppoly_bernstein) :: pp
    type (status_t) :: status

    ! Incompatible data point arrays
    k = 3
    x = 0.0
    y = 0.0
    n = size(x)
    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))
    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x(1:2), y, k, knots, coefs, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Incompatible number of data points"

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x, y(1:1,:), k, knots, coefs, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Too few derivatives for requested polynomial degree"

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x, y(:0,:), -1, knots, coefs, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid polynomial degree"

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x, y, k, knots(1:nknots-1), coefs, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Incompatible KNOTS array size"

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x, y, k, knots, coefs(1:ncoefs-1), status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Incompatible COEFS array size"
end subroutine

elemental subroutine fcn1 (x, fx, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), optional :: fpx

    if (present(fx)) fx = log(x)
    if (present(fpx)) fpx = 1.0_PREC / x
end subroutine

elemental subroutine fcn2 (x, fx, fpx)
    real (PREC), intent(in) :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), optional :: fpx

    real (PREC) :: lfx, lfpx
    real (PREC), parameter :: coefs(*) = [-123.34d0, 4.56d0, -7.23d0, 0.234d0]
    real (PREC) :: coefs_d(3)

    call polyval (coefs, x, lfx)
    call polyder (coefs, coefs_d, m=1)
    call polyval (coefs_d, x, lfpx)

    if (present(fx)) fx = lfx
    if (present(fpx)) fpx = lfpx
end subroutine

subroutine test_power_ppolyval_cubic ()
    !*  Unit tests for evaluating cubic polynomials and their derivatives
    !   using a power basis.

    type (ppoly) :: pp, pp_d1
    type (ppoly_bernstein) :: pp_bern, pp_bern_d1

    real (PREC), dimension(:), allocatable :: xx, yy, yy_bern, fx, fpx
    real (PREC), dimension(:,:), allocatable :: ydat
    real (PREC), dimension(:), allocatable :: knots, coefs_bern, coefs
    real (PREC), dimension(:), allocatable :: coefs_bern_d1, coefs_d1
    real (PREC), dimension(:), allocatable :: rwork
    real (PREC) :: lb, ub
    integer :: n, deg, nknots, ncoefs, ncoefs_bern
    type (status_t) :: status
    logical :: is_ok, values_ok

    ! Fit some polynomial using the Bernstein basis
    nknots = 20
    deg = 3
    allocate (xx(nknots), ydat(2,nknots), fx(nknots), fpx(nknots))

    lb = 0.0
    ub = 10.0
    call linspace (xx, 0.0_PREC, 10.0_PREC)
    fx(:) = sin(xx)
    fpx(:) = cos(xx)
    ydat(1,:) = fx
    ydat(2,:) = fpx

    ncoefs = ppoly_get_ncoefs (pp_bern, n=nknots, k=deg)
    allocate (knots(nknots), coefs_bern(ncoefs))

    call ppolyfit (pp_bern, xx, ydat, deg, knots, coefs_bern, status)

    deallocate (xx)

    ncoefs = ppoly_get_ncoefs (pp, n=nknots, k=deg)
    allocate (coefs(ncoefs))

    ! Convert to power-basis
    call ppoly_transform_basis (pp_bern, coefs_bern, pp, coefs, status)

    ! Evaluate at original points
    allocate (yy(nknots))
    call ppolyval (pp, knots, coefs, knots, yy, status=status)

    values_ok = all_close (yy, fx, atol=1.0e-10_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Eval at original points"

    deallocate (yy)

    ! Compare to points evaluated via bernstein basis at points other
    ! than original
    n = 101
    allocate (xx(n), yy(n), yy_bern(n))

    call linspace (xx, lb, ub)
    call ppolyval (pp_bern, knots, coefs_bern, xx, yy_bern, status=status)
    call ppolyval (pp, knots, coefs, xx, yy, status=status)

    is_ok = all_close (yy, yy_bern, atol=1.0e-10_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. is_ok)) &
        error stop "Eval at non-original points, comparing to Bernstein basis"

    deallocate (xx, yy, yy_bern)

    ! Compare 1st derivatives
    ncoefs_bern = ppoly_get_ncoefs (pp_bern_d1, n=nknots, k=deg-1)
    ncoefs = ppoly_get_ncoefs (pp_d1, n=nknots, k=deg-1)
    allocate (coefs_d1(ncoefs), coefs_bern_d1(ncoefs_bern))

    call ppolyder (pp_bern, knots, coefs_bern, pp_bern_d1, coefs_bern_d1, m=1, status=status)
    call ppolyder (pp, knots, coefs, pp_d1, coefs_d1, m=1, status=status)

    ! Compare values to bernstein values
    n = 201
    allocate (xx(n), yy(n), yy_bern(n))
    call linspace (xx, lb, ub)

    call ppolyval (pp_bern_d1, knots, coefs_bern_d1, xx, yy_bern, status=status)
    call ppolyval (pp_d1, knots, coefs_d1, xx, yy, status=status)

    is_ok = all_close (yy, yy_bern, atol=1.0e-10_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. is_ok)) &
        error stop "D1: Comparing to Bernstein value at non-original points"

    deallocate (xx, yy, yy_bern)

    ! Compute first deratives to true values at original knots as these
    ! were used to fit the data.
    allocate (yy(nknots))
    call ppolyval (pp_d1, knots, coefs_d1, knots, yy, status=status)
    ! Create contiguous array to be used as argument
    allocate (rwork(size(ydat, 2)), source=ydat(2,:))
    is_ok = all_close (yy, rwork, atol=1.0e-8_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. is_ok)) &
        error stop "D1: Comparing to true values at original knots"

    deallocate (rwork)

end subroutine

subroutine test_power_ppolyval_input ()
    !*  Unit tests for input checking for Bernstein basis fitting routines

    real (PREC), dimension(:), allocatable :: knots, coefs
    real (PREC), dimension(:), allocatable :: x, y
    integer :: n, k, ncoefs, nknots, nx
    real (PREC) :: xi, yi, wgt
    integer (NF_ENUM_KIND) :: ext
    integer :: ilb
    type (ppoly) :: pp
    type (status_t) :: status
    logical :: status_ok

    xi = 0.0

    ! === "Constant" extrapolation and missing LEFT, RIGHT ===
    k = 1
    nknots = 2
    ncoefs = ppoly_get_ncoefs (pp, k=k, n=nknots)
    nx = 1
    allocate (knots(nknots), coefs(ncoefs))
    allocate (x(nx), y(nx))
    ext = NF_INTERP_EVAl_CONST

    ! --- Missing both LEFT, RIGHT ---

    ! 0d API
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, xi, yi, ext=ext, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "0d API: missing LEFT, RIGHT"

    ! 1d API
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x, y, ext=ext, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "1d API: missing LEFT, RIGHT"

    ! --- Missing LEFT ---
    ! 0d API
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, xi, yi, ext=ext, right=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "0d API: missing LEFT"

    ! 1d API
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x, y, ext=ext, right=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "1d API: missing LEFT"

    ! --- Missing RIGHT ---
    ! 0d API
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, xi, yi, ext=ext, left=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "0d API: missing RIGHT"

    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x, y, ext=ext, left=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "1d API: missing RIGHT"

    deallocate (knots, coefs, x, y)

    ! === Invalid KNOTS array size ===
    k = 3
    nknots = 1
    ncoefs = ppoly_get_ncoefs (pp, k=3, n=nknots)
    allocate (knots(nknots), coefs(ncoefs))
    ! 0d API
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, xi, yi, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "0d API: Invalid KNOTS array size"

    ! 1d API
    nx = 1
    allocate (x(nx), y(nx))
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x, y, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "1d API: Invalid KNOTS array size"

    deallocate (knots, coefs, x, y)

    ! === KNOTS, COEFS input arrays ===

    nknots = 2
    allocate (knots(nknots))

    ! 0d API
    status_ok = .true.
    do k = 0, 3
        ncoefs = ppoly_get_ncoefs (pp, k=k, n=nknots)
        allocate (coefs(ncoefs+1))

        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, xi, yi, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_INVALID_ARG)

        deallocate (coefs)

        allocate (coefs(max(0, ncoefs-1)))
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, xi, yi, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_INVALID_ARG)

        deallocate (coefs)
    end do

    if (.not. (status_ok)) error stop "0d API: Non-conformable KNOTS, COEFS"

    ! 1d API
    nx = 2
    allocate (x(nx), y(nx))
    status_ok = .true.

    do k = 0, 3
        ncoefs = ppoly_get_ncoefs (pp, k=k, n=nknots)
        allocate (coefs(ncoefs+1))

        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, x, y, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_INVALID_ARG)

        deallocate (coefs)

        allocate (coefs(max(0, ncoefs-1)))
        status = NF_STATUS_UNDEFINED
        call ppolyval (pp, knots, coefs, x, y, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_INVALID_ARG)

        deallocate (coefs)
    end do

    if (.not. (status_ok)) error stop "1d API: Non-conformable KNOTS, COEFS"

    deallocate (knots, x, y)

    ! === Non-conformable X, Y arrays ===

    k = 2
    nknots = 10
    ncoefs = ppoly_get_ncoefs (pp, k=k, n=nknots)
    nx = 5
    allocate (knots(nknots), coefs(ncoefs))
    allocate (x(nx), y(nx+1))

    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x, y, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "1d API: Non-conformable X, Y"

end subroutine

subroutine test_power_ppolyval_eval_input ()
    !*  Unit tests for input checking for Bernstein basis fitting routines

    real (PREC), dimension(:), allocatable :: coefs, weight
    real (PREC), dimension(:), allocatable :: y
    integer, dimension(:), allocatable :: ilbound
    integer :: n, k, ncoefs, nknots, nx
    real (PREC) :: yi, wgt
    integer (NF_ENUM_KIND) :: ext
    integer :: ilb
    type (ppoly) :: pp
    type (status_t) :: status
    logical :: status_ok

    ilb = 1
    wgt = 0.0

    ! === "Constant" extrapolation and missing LEFT, RIGHT ===
    k = 3
    nknots = 2
    ncoefs = ppoly_get_ncoefs (pp, k=k, n=nknots)
    nx = 1
    allocate (coefs(ncoefs), ilbound(nx), weight(nx), y(nx))
    ext = NF_INTERP_EVAl_CONST

    ! --- Missing both LEFT, RIGHT ---

    ! 0d API
    status = NF_STATUS_UNDEFINED
    call ppolyval_eval (pp, ilb, wgt, coefs, yi, ext=ext, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "0d API: missing LEFT, RIGHT"

    ! 1d API
    status = NF_STATUS_UNDEFINED
    call ppolyval_eval (pp, ilbound, weight, coefs, y, ext=ext, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "1d API: missing LEFT, RIGHT"

    ! --- Missing LEFT ---
    ! 0d API
    status = NF_STATUS_UNDEFINED
    call ppolyval_eval (pp, ilb, wgt, coefs, yi, ext=ext, right=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "0d API: missing LEFT"

    ! 1d API
    status = NF_STATUS_UNDEFINED
    call ppolyval_eval (pp, ilbound, weight, coefs, y, ext=ext, right=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "1d API: missing LEFT"

    ! --- Missing RIGHT ---
    ! 0d API
    status = NF_STATUS_UNDEFINED
    call ppolyval_eval (pp, ilb, wgt, coefs, yi, ext=ext, left=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "0d API: missing RIGHT"

    ! 1d API
    status = NF_STATUS_UNDEFINED
    call ppolyval_eval (pp, ilbound, weight, coefs, y, ext=ext, left=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "1d API: missing RIGHT"

    deallocate (ilbound, weight, coefs, y)

    ! === Invalid COEFS array size ===

    nknots = 2

    ! 0d API
    status_ok = .true.
    do k = 0, 3
        call ppoly_init (pp, nknots, k)
        ncoefs = ppoly_get_ncoefs (pp, k=k, n=nknots)
        allocate (coefs(ncoefs+1))

        status = NF_STATUS_UNDEFINED
        call ppolyval_eval (pp, ilb, wgt, coefs, yi, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_INVALID_ARG)

        deallocate (coefs)

        allocate (coefs(max(0, ncoefs-1)))
        status = NF_STATUS_UNDEFINED
        call ppolyval_eval (pp, ilb, wgt, coefs, yi, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_INVALID_ARG)

        deallocate (coefs)
    end do

    if (.not. (status_ok)) error stop "0d API: Invalid COEFS array size"

    ! 1d API
    nx = 5
    allocate (ilbound(nx), weight(nx), y(nx))
    status_ok = .true.

    do k = 0, 3
        call ppoly_init (pp, nknots, k)
        ncoefs = ppoly_get_ncoefs (pp, k=k, n=nknots)
        allocate (coefs(ncoefs+1))

        status = NF_STATUS_UNDEFINED
        call ppolyval_eval (pp, ilbound, weight, coefs, y, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_INVALID_ARG)

        deallocate (coefs)

        allocate (coefs(max(0, ncoefs-1)))
        status = NF_STATUS_UNDEFINED
        call ppolyval_eval (pp, ilbound, weight, coefs, y, status=status)
        status_ok = status_ok .and. (status == NF_STATUS_INVALID_ARG)

        deallocate (coefs)
    end do

    if (.not. (status_ok)) error stop "1d API: Invalid COEFS array size"

    deallocate (ilbound, weight, y)

    ! === Non-conformable ILBOUND, WEIGHT arrays ===

    k = 2
    nknots = 10
    call ppoly_init (pp, nknots, k)
    ncoefs = ppoly_get_ncoefs (pp, k=k, n=nknots)
    nx = 5
    allocate (coefs(ncoefs))
    allocate (ilbound(nx), weight(nx+1), y(nx))

    status = NF_STATUS_UNDEFINED
    call ppolyval_eval (pp, ilbound, weight, coefs, y, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "1d API: Non-conformable ILBOUND, WEIGHT"

    deallocate (coefs, ilbound, weight, y)

    ! === Non-conformable ILBOUND, Y arrays ===

    k = 1
    nknots = 4
    ncoefs = ppoly_get_ncoefs (pp, k=k, n=nknots)
    call ppoly_init (pp, nknots, k)
    nx = 1
    allocate (coefs(ncoefs))
    allocate (ilbound(nx), weight(nx), y(nx-1))

    status = NF_STATUS_UNDEFINED
    call ppolyval_eval (pp, ilbound, weight, coefs, y, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "1d API: Non-conformable ILBOUND, Y"

    deallocate (coefs, ilbound, weight, y)

end subroutine

end
