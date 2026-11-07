program test_interp_linear
    !*  Unit tests for linear interpolation in various dimensions.

    use iso_fortran_env
    use numfort_arrays
    use numfort_common_testing
    use numfort_interpolate

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

contains

subroutine test_all ()

    call test_linear_input_checks ()
    call test_exact ()
    call test_extrap ()
    call test_truncate ()
    call test_random ()

    call test_bilinear_input_checks ()
    call test_bilinear_interpolate ()
    call test_bilinear_extrapolate ()

end subroutine

subroutine test_linear_input_checks ()

    real (PREC), dimension(:), allocatable :: x, xp, fp, fx
    real (PREC) :: xi, fxi, left, right
    integer (NF_ENUM_KIND) :: ext
    integer :: np, nx, i
    type (status_t) :: status

    ! === Invalid XP, FP array sizes ===
    nx = 10
    allocate (x(nx), fx(nx))
    xi = 0.0

    do i = 0, 1
        allocate (xp(i), fp(i))

        ! Scalar interface
        status = NF_STATUS_UNDEFINED
        call interp_linear (xi, xp, fp, fxi, status=status)
        if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "0d API: size(XP) < 2"

        ! Array interface
        status = NF_STATUS_UNDEFINED
        call interp_linear (x, xp, fp, fx, status=status)
        if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "1d API: size(XP) < 2"

        deallocate (xp, fp)
    end do

    deallocate (x, fx)

    ! === Non-conformable XP, FP arrays
    nx = 5
    allocate (x(nx), fx(nx))
    allocate (xp(2), fp(3))

    ! Scalar interface
    status = NF_STATUS_UNDEFINED
    call interp_linear (xi, xp, fp, fxi, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "0d API: size(XP) != size(FP)"

    ! Array interface
    status = NF_STATUS_UNDEFINED
    call interp_linear (x, xp, fp, fx, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "1d API: size(XP) != size(FP)"

    deallocate (xp, fp, x, fx)

    ! === Non-conformable X, FX array ===
    np = 10
    allocate (xp(np), fp(np))
    allocate (x(1), fx(2))

    status = NF_STATUS_UNDEFINED
    call interp_linear (x, xp, fp, fx, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "1d API: Non-conformable X, FX"

    deallocate (xp, fp, x, fx)

    ! === Missing LEFT, RIGHT ===
    ext = NF_INTERP_EVAL_CONST
    nx = 10
    np = 9
    allocate (xp(np))
    call arange (xp, 1.0_PREC)
    allocate (fp(np), source=xp)

    allocate (x(nx), fx(nx))
    call random_number (x)

    ! Scalar interface
    xi = 1.0
    status = NF_STATUS_UNDEFINED
    call interp_linear (xi, xp, fp, fxi, ext=ext, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "0d API: Missing LEFT, RIGHT arguments"

    status = NF_STATUS_UNDEFINED
    call interp_linear (xi, xp, fp, fxi, ext=ext, left=1.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "0d API: Missing RIGHT argument"

    status = NF_STATUS_UNDEFINED
    call interp_linear (xi, xp, fp, fxi, ext=ext, right=1.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "1d API: Missing LEFT argument"

    ! Array interface
    status = NF_STATUS_UNDEFINED
    call interp_linear (x, xp, fp, fx, ext=ext, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "0d API: Missing LEFT, RIGHT arguments"

    status = NF_STATUS_UNDEFINED
    call interp_linear (x, xp, fp, fx, ext=ext, left=1.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "0d API: Missing RIGHT argument"

    status = NF_STATUS_UNDEFINED
    call interp_linear (x, xp, fp, fx, ext=ext, right=1.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "1d API: Missing LEFT argument"

end subroutine

subroutine test_extrap ()
    !*  Tests for linear extrapolation when x-coordinate of interpolated value
    !   is outside of data points boundaries.

    real (PREC), dimension(10) :: xp, fp, fx
    integer :: i, j, k

    ! function f(x) = 2x + 1 on [0, 1]
    call linspace (xp, 0.0d0, 1.0d0)
    fp = 2 * xp + 1.0d0

    ! ===== Scalar arguments =====
    ! evaluate outside of domain with extrapolation
    fx = 0.0
    ! limit xp/fp to subrange of values
    j = size(xp) / 2 - 1
    k = size(xp) / 2 + 1

    do i = 1, size(xp)
        call interp_linear (xp(i), xp(j:k), fp(j:k), fx(i), &
            ext=NF_INTERP_EVAL_EXTRAPOLATE)
    end do

    if (.not. (all (abs (fx - fp) < 1d-12))) error stop "Extrapolation, scalar argument"

    ! ===== Array arguments =====
    fx = 0.0
    call interp_linear (xp, xp(j:k), fp(j:k), fx, ext=NF_INTERP_EVAL_EXTRAPOLATE)
    if (.not. (all (abs (fx - fp) < 1d-12))) error stop "Extrapolation, array argument"

end subroutine

subroutine test_exact ()
    !*  Tests for "interpolating" exactly at function data points.

    real (PREC), dimension(10) :: xp, fp, fx, x
    integer :: i

    ! function f(x) = 2x + 1 on [0, 1]
    call linspace (xp, 0.0d0, 1.0d0)
    fp = 2 * xp + 1.0d0

    ! ===== Scalar arguments =====
    ! evaluate at all points in xp
    fx = 0.0
    do i = 1, size(xp)
        call interp_linear (xp(i), xp, fp, fx(i))
    end do
    if (.not. (all (abs (fx - fp) < 1d-12))) &
        error stop "Exact interpolation, scalar argument"

    ! ===== Array arguments =====
    x = xp
    fx = 0.0d0

    call interp_linear (x, xp, fp, fx)
    if (.not. (all (abs (fp - fx) < 1d-12))) &
        error stop "Exact interpolation, array argument"

end subroutine

subroutine test_random ()
    !*  Tests for interpolating at random points within the boundaries defined
    !   but data points.

    real (PREC), dimension(10) :: xp, fp, fx, x, fx1
    integer :: i
    real (PREC) :: x1, xn

    x1 = 0.0d0
    xn = 1.0d0

    ! function f(x) = 2x + 1 on [0, 1]
    call linspace (xp, x1, xn)
    fp = 2 * xp + 1.0d0

    call random_number (x)
    ! correct function values
    x = x * (xn-x1) + x1
    fx1 = 2 * x + 1.0d0

    ! ===== Scalar arguments =====
    ! evaluate at all points in xp
    fx = 0.0
    do i = 1, size(xp)
        call interp_linear (x(i), xp, fp, fx(i))
    end do

    if (.not. (all (abs (fx - fx1) < 1d-12))) error stop "Scalar argument"

    ! ===== Array arguments =====
    fx = 0.0d0

    call interp_linear (x, xp, fp, fx)
    if (.not. (all (abs (fx - fx1) < 1d-12))) error stop "Array argument"

end subroutine

subroutine test_truncate ()
    !*  Tests for interpolation at x-coordinates outside of data point boundaries
    !   when linear extrapolation is disabled.

    real (PREC), dimension(10) :: xp, fp, fx, fx1
    integer :: i, ilb, iub
    real (PREC) :: left, right
    logical :: all_ok
    type (status_t) :: status

    ! function f(x) = 2x + 1 on [0, 1]
    call linspace (xp, 0.0d0, 1.0d0)
    fp = 2 * xp + 1.0d0

    ! ===== Scalar arguments =====

    ! >>> Truncate outside of specified domain, default left/right values
    fx = 0.0
    ! limit xp/fp to subrange of values
    ilb = size(xp) / 2 - 1
    iub = size(xp) / 2 + 1

    all_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call interp_linear (xp(i), xp(ilb:iub), fp(ilb:iub), fx(i), &
            ext=NF_INTERP_EVAL_BOUNDARY, status=status)
        all_ok = all_ok .and. (status == NF_STATUS_OK)
    end do

    ! Compute expected result
    fx1 = fp
    call truncate (xp, xp(ilb), xp(iub), fx1, fp(ilb), fp(iub))

    all_ok = all_ok .and. all_close (fx, fx1, rtol=0.0_PREC, atol=1.0e-12_PREC)
    if (.not. (all_ok)) &
        error stop "Truncation: scalar argument, default left/right values"

    ! === Truncation, exlicitly specify LEFT, RIGHT
    fx = 0.0
    left = -1.0
    right = 10.123_PREC
    all_ok = .true.
    do i = 1, size(xp)
        status = NF_STATUS_UNDEFINED
        call interp_linear (xp(i), xp(ilb:iub), fp(ilb:iub), fx(i), &
            ext=NF_INTERP_EVAL_CONST, left=left, right=right, status=status)
        all_ok = all_ok .and. (status == NF_STATUS_OK)
    end do

    ! compute expected result
    fx1 = fp
    call truncate (xp, xp(ilb), xp(iub), fx1, left, right)

    all_ok = all_ok .and. all_close (fx, fx1, rtol=0.0_PREC, atol=1.0e-12_PREC)
    if (.not. (all_ok)) error stop "Truncation: scalar argument, explicit left value"

    ! ===== Array arguments =====

    ! >>> Truncate outside of specified domain, default left/right values
    fx = 0.0
    ! limit xp/fp to subrange of values
    ilb = size(xp) / 2 - 1
    iub = size(xp) / 2 + 1

    call interp_linear (xp, xp(ilb:iub), fp(ilb:iub), fx, &
        ext=NF_INTERP_EVAL_BOUNDARY, status=status)

    ! Compute expected result
    fx1 = fp
    call truncate (xp, xp(ilb), xp(iub), fx1, fp(ilb), fp(iub))

    all_ok = all_close (fx, fx1, rtol=0.0_PREC, atol=1.0e-12_PREC)
    if (.not. (all_ok .and. (status == NF_STATUS_OK))) &
        error stop "Truncation: array argument, default left/right values"

    ! >>> Truncation, exlicitly specify left value
    fx = 0.0
    left = -1.0
    right = 234.345_PREC

    call interp_linear (xp, xp(ilb:iub), fp(ilb:iub), fx, ext=NF_INTERP_EVAL_CONST, &
        left=left, right=right, status=status)

    ! compute expected result
    fx1 = fp
    call truncate (xp, xp(ilb), xp(iub), fx1, left, right)

    all_ok = all_close (fx, fx1, rtol=0.0_PREC, atol=1.0e-12_PREC)
    if (.not. (all_ok .and. (status == NF_STATUS_OK))) &
        error stop "Truncation: array argument, explicit left value"

end subroutine

subroutine truncate (x, lb, ub, fx, left, right)
    real (PREC), dimension(:) :: x, fx
    real (PREC), intent(in) :: lb, ub, left, right

    where (x < lb)
        fx = left
    else where (x  > ub)
        fx = right
    end where
end subroutine

subroutine test_bilinear_input_checks ()
    !*  Checks input verification for bilinear interpolation routines.

    real (PREC), dimension(:), allocatable :: xp1, xp2, x1, x2, fx1d
    real (PREC) :: fx
    real (PREC), dimension(:,:), allocatable :: fp
    type (status_t) :: status

    ! Test with grid arrays with less than 2 elements
    allocate (xp1(1), xp2(1), fp(1,1))
    status = NF_STATUS_OK
    call interp_bilinear (1.0_PREC, 1.0_PREC, xp1, xp2, fp, fx, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Grid point arrays have fewer than two elements"
    deallocate (xp1, xp2, fp)

    ! Test with incompatible grid arrays
    allocate (xp1(2), xp2(3), fp(2,2))
    status = NF_STATUS_OK
    call interp_bilinear (1.0_PREC, 1.0_PREC, xp1, xp2, fp, fx, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Grid point arrays are incompatible"
    deallocate (xp1, xp2, fp)

    ! Test with eval input arrays of different size
    allocate (xp1(2), xp2(2), fp(2,2))
    allocate (x1(1), x2(2), fx1d(1))
    status = NF_STATUS_OK
    call interp_bilinear (x1, x2, xp1, xp2, fp, fx1d, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Eval input arrays are incompatible"
    deallocate (xp1, xp2, fp, x1, x2, fx1d)

    ! Test with eval output array of different size
    allocate (xp1(2), xp2(2), fp(2,2))
    allocate (x1(2), x2(2), fx1d(1))
    status = NF_STATUS_OK
    call interp_bilinear (x1, x2, xp1, xp2, fp, fx1d, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Eval input arrays are incompatible"
    deallocate (xp1, xp2, fp, x1, x2, fx1d)

end subroutine

subroutine test_bilinear_interpolate ()
    !*  Tests for bilinear interpolation at interior points

    real (PREC), dimension(:), allocatable :: xp1, xp2, x1, x2
    real (PREC), dimension(:,:), allocatable :: fp
    real (PREC), dimension(:), allocatable :: fx, fx_ok
    integer :: n1, n2, n, i, j
    logical :: res_ok
    type (status_t) :: status

    n1 = 10
    n2 = 11
    allocate (xp1(n1), xp2(n2), fp(n1,n2))

    call linspace (xp1, 0.0_PREC, 10.0_PREC)
    call powerspace (xp2, 0.0_PREC, 5.0_PREC, 2.0_PREC)

    do j = 1, size(xp2)
        fp(:,j) = func1_2d (xp1, xp2(j))
    end do

    ! Test with 1d array interface
    n = 15
    allocate (x1(n), x2(n), fx_ok(n), fx(n))
    call random_number (x1)
    call random_number (x2)

    ! Rescale to approx. cover most parts of rectangular grid
    x1(:) = x1 * 10.0
    x2(:) = x2 * 5.0

    ! Correct function values
    fx_ok(:) = func1_2d(x1, x2)

    ! Interpolated values
    call interp_bilinear (x1, x2, xp1, xp2, fp, fx, status=status)
    res_ok = (status == NF_STATUS_OK) .and. all_close (fx, fx_ok, atol=1.0e-10_PREC)
    if (.not. (res_ok)) error stop "Interpolate at random interior points; array input"

    ! Test scalar interface
    fx(:) = 0.0
    fx_ok(1) = func1_2d (x1(1), x2(1))
    call interp_bilinear (x1(1), x2(1), xp1, xp2, fp, fx(1), status=status)
    res_ok = (status == NF_STATUS_OK) .and. all_close (fx(1), fx_ok(1), atol=1.0e-10_PREC)
    if (.not. (res_ok)) error stop "Interpolate at random interior point; scalar input"

    deallocate (x1, x2, fx, fx_ok)

    ! Test evaluating exactly at grid points
    n = n1 * n2
    allocate (x1(n), x2(n), fx(n), fx_ok(n))
    do j = 1, size(xp2)
        i = (j-1) * size(xp1)
        x1(i+1:i+n1) = xp1
        x2(i+1:i+n1) = xp2(j)
    end do

    fx_ok(:) = func1_2d (x1, x2)

    call interp_bilinear (x1, x2, xp1, xp2, fp, fx, status=status)

    res_ok = .true.
    do j = 1, size(xp2)
        i = (j-1) * n1
        res_ok = res_ok .and. all_close (fx(i+1:i+n1), fp(:,j), atol=1.0e-10_PREC)
    end do
    if (.not. (res_ok)) error stop "Interpolate exactly at grid points"

    deallocate (x1, x2, fx, fx_ok)
end subroutine

subroutine test_bilinear_extrapolate ()
    !*  Tests for bilinear interpolation when not all points are interior

    real (PREC), dimension(:), allocatable :: xp1, xp2, x1, x2
    real (PREC), dimension(:,:), allocatable :: fp
    real (PREC), dimension(:), allocatable :: fx, fx_ok
    integer :: n1, n2, n, i, j
    logical :: res_ok
    type (status_t) :: status
    real (PREC) :: fill_value

    n1 = 10
    n2 = 11
    allocate (xp1(n1), xp2(n2), fp(n1,n2))

    call linspace (xp1, -1.0_PREC, 1.0_PREC)
    call powerspace (xp2, -1.0_PREC, 1.0_PREC, 2.0_PREC)

    do j = 1, size(xp2)
        fp(:,j) = func1_2d (xp1, xp2(j))
    end do

    ! Test with 1d array interface
    n = 10
    allocate (x1(n), x2(n), fx_ok(n), fx(n))
    call random_number (x1)
    call random_number (x2)

    ! Rescale to approx. cover most parts of rectangular grid
    x1(:) = x1 * 2.0 - 1.0
    x2(:) = x2 * 2.0 - 1.0

    ! Set some non-interior points (both or either dimension)
    x1(1) = -1.5
    x2(2) = -1.1
    x1(3) = 1.234d0
    x2(4) = 1.456d0
    x1(5) = -1.234d0
    x2(5) = -2.345d0
    x1(6) = 2.3456d0
    x2(6) = 3.4456d0

    fill_value = -1000.0

    ! Correct function values
    fx_ok(:) = func1_2d (x1, x2)

    ! Overwrite non-interior points with FILL_VALUE
    do i = 1,n
        if (x1(i) < xp1(1) .or. x1(i) > xp1(n) .or. x2(i) < xp2(1) .or. x2(i) > xp2(n)) then
            fx_ok(i) = fill_value
        end if
    end do

    ! Interpolated values; set non-internal values to -1000
    call interp_bilinear (x1, x2, xp1, xp2, fp, fx, fill_value=fill_value, status=status)
    res_ok = (status == NF_STATUS_OK) .and. all_close (fx, fx_ok, atol=1.0e-10_PREC)
    if (.not. (res_ok)) error stop "Interpolation at non-int. values using FILL_VALUE"

    ! Test with extrapolation
    fx_ok(:) = func1_2d (x1, x2)
    call interp_bilinear (x1, x2, xp1, xp2, fp, fx, status=status)
    res_ok = (status == NF_STATUS_OK) .and. all_close (fx, fx_ok, atol=1.0e-10_PREC)
    if (.not. (res_ok)) &
        error stop "Interpolation at non-int. values using extrapolation"

    deallocate (x1, x2, fx, fx_ok)

end subroutine

elemental function func1_2d(x1, x2) result(fx)
    !*  Function f: R^2 -> R that can be interpolated exactly by bilinear
    !   interpolation
    real (PREC), intent(in) :: x1
    real (PREC), intent(in) :: x2
    real (PREC) :: fx

    fx = 2.0 * x1 - 3.0 * x2 + 0.5d0 * x1 * x2 + 1.0
end function

end program
