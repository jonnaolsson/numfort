

program test_ppoly_tensor
    !*  Unit tests for function approximation with tensor products of
    !   piecewise polynomials

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common
    use numfort_common_workspace, workspace => workspace_real64
    use numfort_common_testing
    use numfort_polynomial
    use numfort_interpolate
    use numfort_polynomial

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_ppoly2d_fit_input ()
    call test_ppoly2d_val_input ()
    call test_ppoly2d_val_bilinear ()
    call test_ppoly2d_val_bilinear_der ()

end subroutine test_all

subroutine test_ppoly2d_fit_input ()

    real (PREC), allocatable :: x1(:), x2(:), y(:,:)
    real (PREC), dimension(:), allocatable :: knots, coefs
    integer :: nknots, ncoefs, n1, n2, n(2), k
    type (ppoly2d) :: pp
    type (status_t) :: status

    ! Check with arrays X1, X2 that are both too small
    n1 = 0
    n2 = 1
    n = [n1, n2]
    k = 1
    allocate (x1(n1), x2(n2), y(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1, x2, y, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "X1, x2 too small"
    deallocate (x1, x2, y, knots, coefs)

    ! Check with array X1 too small
    n1 = 1
    n2 = 2
    n = [n1, n2]
    k = 1
    allocate (x1(n1), x2(n2), y(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1, x2, y, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "X1 too small"
    deallocate (x1, x2, y, knots, coefs)

    ! Check with array X2 too small
    n1 = 2
    n2 = 1
    n = [n1, n2]
    k = 1
    allocate (x1(n1), x2(n2), y(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1, x2, y, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "X2 too small"
    deallocate (x1, x2, y, knots, coefs)

    ! Check with non-conformable array Y
    n1 = 2
    n2 = 3
    n = [n1, n2]
    k = 1
    allocate (x1(n1), x2(n2), y(n1-1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1, x2, y, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Y not conformable with X1, X2"
    deallocate (x1, x2, y, knots, coefs)

    ! Check with KNOTS array too small
    n1 = 2
    n2 = 3
    n = [n1, n2]
    k = 1
    allocate (x1(n1), x2(n2), y(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots-1), coefs(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1, x2, y, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "KNOTS array too small"
    deallocate (x1, x2, y, knots, coefs)

    ! Check with COEFS array too small
    n1 = 2
    n2 = 3
    n = [n1, n2]
    k = 1
    allocate (x1(n1), x2(n2), y(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs-1))

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1, x2, y, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "COEFS array too small"
    deallocate (x1, x2, y, knots, coefs)

    ! Check with KNOTS, COEFS arrays having non-conformable size
    n1 = 5
    n2 = 3
    n = [n1, n2]
    k = 1
    allocate (x1(n1), x2(n2), y(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs+1))

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1, x2, y, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "KNOTS, COEFS have non-conformable size"
    deallocate (x1, x2, y, knots, coefs)

    ! Invalid polynomial degree
    n1 = 5
    n2 = 3
    n = [n1, n2]
    k = -1
    allocate (x1(n1), x2(n2), y(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1, x2, y, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid polynomial degree"
    deallocate (x1, x2, y, knots, coefs)

end subroutine

subroutine test_ppoly2d_val_input ()

    real (PREC), allocatable :: x1a(:), x2a(:), ya(:), x1in(:), x2in(:), yin(:,:)
    real (PREC), dimension(:), allocatable :: knots, coefs
    real (PREC) :: x1, x2, y
    integer :: nknots, ncoefs, n1, n2, n(2), k, nx, i
    integer (NF_ENUM_KIND) :: ext
    type (ppoly2d) :: pp
    type (status_t) :: status

    n1 = 3
    n2 = 6
    n = [n1, n2]
    k = 1

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (x1in(n1), x2in(n2), yin(n1,n2))
    allocate (knots(nknots), coefs(ncoefs))

    ! Create some arbitrary domain and function values
    call linspace (x1in, 0.0_PREC, 10.0_PREC)
    call linspace (x2in, -1.234_PREC, 2.3456_PREC)

    do i = 1, n2
        call fcn1 (x1in, x2in(i), yin(:,i))
    end do

    ! Fit polynomial, otherwise attributes of PPOLY2D object will not be
    ! properly set.
    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1in, x2in, yin, k, knots, coefs, status)
    deallocate (x1in, x2in, yin)

    ! Non-conformable X1, X2
    nx = 2
    allocate (x1a(nx), x2a(nx+1), ya(nx))
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X1, X2"
    deallocate (x1a, x2a, ya)

    ! Non-conformable X1, Y
    nx = 2
    allocate (x1a(nx), x2a(nx), ya(nx-1))
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X1, Y"
    deallocate (x1a, x2a, ya)

    ! KNOTS has wrong size
    nx = 5
    allocate (x1a(nx), x2a(nx), ya(nx))

    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots(1:nknots-1), coefs, x1a, x2a, ya, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "KNOTS array has wrong size"
    deallocate (x1a, x2a, ya)

    ! COEFS has wrong size
    nx = 5
    allocate (x1a(nx), x2a(nx), ya(nx))
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs(1:ncoefs-1), x1a, x2a, ya, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "COEFS array has wrong size"
    deallocate (x1a, x2a, ya)

    ! EXT argument has unsupported value
    nx = 5
    allocate (x1a(nx), x2a(nx), ya(nx))

    ext = -1
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, ext=ext, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "EXT argument has unsupported value"
    deallocate (x1a, x2a, ya)

    ! EXT = NF_INTEPR_EVAL_CONST without EXT_VAL specified
    nx = 5
    allocate (x1a(nx), x2a(nx), ya(nx))

    ext = NF_INTERP_EVAL_CONST
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, ext=ext, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "EXT=NF_INTERP_EVAL_CONST but EXT_VAL missing"
    deallocate (x1a, x2a, ya)

    ! Scalar interface: KNOTS array has wrong size
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots(1:nknots-1), coefs, x1, x2, y, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Scalar API: KNOTS array has wrong size"

    ! Scalar interface: COEFS array has wrong size
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs(1:ncoefs-1), x1, x2, y, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Scalar API: COEFS array has wrong size"

    ! Scalar interface: EXT argument has unsupported value
    ext = -1
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, ext=ext, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Scalar API: EXT argument has unsupported value"

    ! EXT = NF_INTEPR_EVAL_CONST without EXT_VAL specified
    ext = NF_INTERP_EVAL_CONST
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, ext=ext, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Scalar API: EXT=NF_INTERP_EVAL_CONST but EXT_VAL missing"

    deallocate (knots, coefs)

end subroutine

subroutine test_ppoly2d_val_bilinear ()
    !*  Unit tests for evaluating piecewise bilinear polynomials

    real (PREC), allocatable :: x1a(:), x2a(:), ya(:), x1in(:), x2in(:), yin(:,:)
    real (PREC), allocatable :: ya_ok(:)
    real (PREC), dimension(:), allocatable :: knots, coefs
    integer :: nknots, ncoefs, n1, n2, n(2), k, nx, i
    real (PREC) :: x1bnd(2), x2bnd(2)
    integer (NF_ENUM_KIND) :: ext
    type (ppoly2d) :: pp
    type (status_t) :: status

    ! Test 1: min. number of knots, bilinear function
    n1 = 2
    n2 = 2
    n = [n1, n2]
    k = 1
    x1bnd = [0.0_PREC, 5.0_PREC]
    x2bnd = [-1.0_PREC, 1.234_PREC]
    allocate (x1in(n1), x2in(n2), yin(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    call linspace (x1in, x1bnd(1), x1bnd(2))
    call linspace (x2in, x2bnd(1), x2bnd(2))

    do i = 1, n2
        call fcn1 (x1in, x2in(i), yin(:,i))
    end do

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1in, x2in, yin, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Knots: [2, 2]; fit to bilinear polynomial"
    deallocate (x1in, x2in, yin)

    ! === Evaluate bilinear function at interior points ===
    nx = 11
    allocate (x1a(nx), x2a(nx), ya(nx), ya_ok(nx))
    call linspace (x1a, x1bnd(1)+0.01d0, x1bnd(2)-0.01d0)
    call linspace (x2a, x2bnd(1)+0.01d0, x2bnd(2)-0.01d0)

    ! Evaluate true values
    call fcn1 (x1a, x2a, ya_ok)
    ! Evaluate approx
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (ya, ya_ok, atol=1.0d-10))) &
        error stop "Knots: [2, 2]; evaluate at interior points"

    deallocate (x1a, x2a, ya, ya_ok)

    ! === Evaluate using extrapolation ===
    nx = 15
    allocate (x1a(nx), x2a(nx), ya(nx), ya_ok(nx))
    call linspace (x1a, x1bnd(1)-1.0d0, x1bnd(2) + 1.0d0)
    call linspace (x2a, x2bnd(1)-1.0d0, x2bnd(2) + 1.0d0)

    ! Evaluate true values
    call fcn1 (x1a, x2a, ya_ok)
    ! Evaluate approx
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (ya, ya_ok, atol=1.0d-10))) &
        error stop "Knots: [2, 2]; evaluate at non-interior points"

    deallocate (x1a, x2a, ya, ya_ok)
    deallocate (knots, coefs)

    ! === Grid with more than 2x2 knots ===
    n1 = 11
    n2 = 7
    n = [n1, n2]
    k = 1
    x1bnd = [2.0_PREC, 5.0_PREC]
    x2bnd = [-1.0_PREC, 1.234_PREC]
    allocate (x1in(n1), x2in(n2), yin(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    call linspace (x1in, x1bnd(1), x1bnd(2))
    call linspace (x2in, x2bnd(1), x2bnd(2))

    do i = 1, n2
        call fcn1 (x1in, x2in(i), yin(:,i))
    end do

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1in, x2in, yin, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Knots: [11, 7]; fit bilinear polynomial"
    deallocate (x1in, x2in, yin)

    ! === Evaluate bilinear function at interior points ===
    nx = 11
    allocate (x1a(nx), x2a(nx), ya(nx), ya_ok(nx))
    call linspace (x1a, x1bnd(1)+0.01d0, x1bnd(2)-0.01d0)
    call linspace (x2a, x2bnd(1)+0.01d0, x2bnd(2)-0.01d0)

    ! Evaluate true values
    call fcn1 (x1a, x2a, ya_ok)
    ! Evaluate approx
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (ya, ya_ok, atol=1.0d-10))) &
        error stop "Knots: [11, 7]; evaluate at interior points"

    deallocate (x1a, x2a, ya, ya_ok)

    ! === Evaluate using extrapolation ===
    nx = 15
    allocate (x1a(nx), x2a(nx), ya(nx), ya_ok(nx))
    call linspace (x1a, x1bnd(1)-1.0d0, x1bnd(2) + 1.0d0)
    call linspace (x2a, x2bnd(1)-1.0d0, x2bnd(2) + 1.0d0)

    ! Evaluate true values
    call fcn1 (x1a, x2a, ya_ok)
    ! Evaluate approx
    status = NF_STATUS_UNDEFINED
    ext = NF_INTERP_EVAL_EXTRAPOLATE
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, ext=ext, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (ya, ya_ok, atol=1.0d-10))) &
        error stop "Knots: [11, 7]; evaluate at non-interior points"

    deallocate (x1a, x2a, ya, ya_ok)
    deallocate (knots, coefs)

    ! === Linear test function ===
    n1 = 8
    n2 = 15
    n = [n1, n2]
    k = 1
    x1bnd = [-2.0_PREC, 5.0_PREC]
    x2bnd = [-1.0_PREC, 1.234_PREC]
    allocate (x1in(n1), x2in(n2), yin(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    call linspace (x1in, x1bnd(1), x1bnd(2))
    call linspace (x2in, x2bnd(1), x2bnd(2))

    do i = 1, n2
        call fcn2 (x1in, x2in(i), yin(:,i))
    end do

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1in, x2in, yin, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Knots: [8, 15]; fit to linear polynomial"
    deallocate (x1in, x2in, yin)

    ! === Evaluate using extrapolation ===
    nx = 15
    allocate (x1a(nx), x2a(nx), ya(nx), ya_ok(nx))
    call linspace (x1a, x1bnd(1)-1.0d0, x1bnd(2) + 1.0d0)
    call linspace (x2a, x2bnd(1)-1.0d0, x2bnd(2) + 1.0d0)

    ! Evaluate true values
    call fcn2 (x1a, x2a, ya_ok)
    ! Evaluate approx
    status = NF_STATUS_UNDEFINED
    ext = NF_INTERP_EVAL_EXTRAPOLATE
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, ext=ext, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (ya, ya_ok, atol=1.0d-10))) &
        error stop "Knots: [8, 15]; evaluate at non-interior points"

    deallocate (x1a, x2a, ya, ya_ok)
    deallocate (knots, coefs)

    ! === Scalar interface ===
    n1 = 21
    n2 = 33
    n = [n1, n2]
    k = 1
    x1bnd = [-5.0_PREC, 5.0_PREC]
    x2bnd = [-10.0_PREC, 10.234_PREC]
    allocate (x1in(n1), x2in(n2), yin(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    call linspace (x1in, x1bnd(1), x1bnd(2))
    call linspace (x2in, x2bnd(1), x2bnd(2))

    do i = 1, n2
        call fcn1 (x1in, x2in(i), yin(:,i))
    end do

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1in, x2in, yin, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Knots: [21, 33]; fit to bilinear polynomial"
    deallocate (x1in, x2in, yin)

    ! === Evaluate using extrapolation ===
    nx = 55
    allocate (x1a(nx), x2a(nx), ya(nx), ya_ok(nx))
    call linspace (x1a, x1bnd(1)-10.0d0, x1bnd(2) + 10.0d0)
    call linspace (x2a, x2bnd(1)-10.0d0, x2bnd(2) + 10.0d0)

    ! Evaluate true values
    call fcn1 (x1a, x2a, ya_ok)
    ! Evaluate approx
    ext = NF_INTERP_EVAL_EXTRAPOLATE
    do i = 1, nx
        call ppolyval (pp, knots, coefs, x1a(i), x2a(i), ya(i), ext=ext)
    end do
    if (.not. (all_close (ya, ya_ok, atol=1.0d-10))) &
        error stop "Knots: [21, 33]; evaluate at non-interior points using scalar " &
            // "API"

    deallocate (x1a, x2a, ya, ya_ok)
    deallocate (knots, coefs)

    ! === Test with func. where bilinear approx. is not exact ===
    n1 = 11
    n2 = 7
    n = [n1, n2]
    k = 1
    x1bnd = [-2.0_PREC, 5.0_PREC]
    x2bnd = [1.0_PREC, 5.234_PREC]
    allocate (x1in(n1), x2in(n2), yin(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    call linspace (x1in, x1bnd(1), x1bnd(2))
    call linspace (x2in, x2bnd(1), x2bnd(2))

    do i = 1, n2
        call fcn3 (x1in, x2in(i), yin(:,i))
    end do

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1in, x2in, yin, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Knots: [11, 7]; fit to higher-order polynomial"

    ! === Evaluate at non-interior points ===
    nx = 1
    allocate (x1a(nx), x2a(nx), ya(nx), ya_ok(nx))
    call linspace (x1a, x1bnd(1)-10.0d0, x1bnd(2)+10.0d0)
    call linspace (x2a, x2bnd(1)-10.0d0, x2bnd(2)+10.0d0)

    ! Evaluate using alternative bilinear implementation
    call interp_bilinear (x1a, x2a, x1in, x2in, yin, ya_ok, &
        ext=NF_INTERP_EVAL_EXTRAPOLATE, status=status)
    ! Evaluate approx
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, ya, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (ya, ya_ok, atol=1.0d-10))) &
        error stop "Knots: [11, 7]; eval at non-interior points; compare to " // &
            "INTERP_BILINEAR"

end subroutine

subroutine test_ppoly2d_val_bilinear_der ()
    !*  Unit tests for evaluating partial derivatives of piecewise bilinear
    !   polynomials

    real (PREC), allocatable :: x1in(:), x2in(:), yin(:,:)
    real (PREC), dimension(:), allocatable :: x1a, x2a, dx1, dx2
    real (PREC), dimension(:), allocatable :: dx1_ok, dx2_ok
    real (PREC), dimension(:), allocatable :: knots, coefs
    integer :: nknots, ncoefs, n1, n2, n(2), k, nx, i
    real (PREC) :: x1bnd(2), x2bnd(2)
    type (ppoly2d) :: pp
    type (status_t) :: status

    ! === Fit to bilinear function ===
    n1 = 2
    n2 = 2
    n = [n1, n2]
    k = 1
    x1bnd = [0.0_PREC, 5.0_PREC]
    x2bnd = [-1.0_PREC, 1.234_PREC]
    allocate (x1in(n1), x2in(n2), yin(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    call linspace (x1in, x1bnd(1), x1bnd(2))
    call linspace (x2in, x2bnd(1), x2bnd(2))

    do i = 1, n2
        call fcn1 (x1in, x2in(i), yin(:,i))
    end do

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1in, x2in, yin, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Knots: [2, 2]; fit to bilinear polynomial"
    deallocate (x1in, x2in, yin)

    ! === Evaluate bilinear function at interior points ===
    nx = 11
    allocate (x1a(nx), x2a(nx), dx1(nx), dx2(nx), dx1_ok(nx), dx2_ok(nx))
    call linspace (x1a, x1bnd(1)+0.01d0, x1bnd(2)-0.01d0)
    call linspace (x2a, x2bnd(1)+0.01d0, x2bnd(2)-0.01d0)

    ! Evaluate true values
    call fcn1 (x1a, x2a, dx1=dx1_ok, dx2=dx2_ok)
    ! Evaluate approx
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, dx1, m=1, dim=1, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dx1, dx1_ok, atol=1.0d-10))) &
        error stop "Knots: [2, 2]; eval derivative for dim=1 at interior points"

    call ppolyval (pp, knots, coefs, x1a, x2a, dx2, m=1, dim=2, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dx2, dx2_ok, atol=1.0d-10))) &
        error stop "Knots: [2, 2]; eval derivative for dim=2 at interior points"

    deallocate (x1a, x2a, dx1, dx2, dx1_ok, dx2_ok)

    ! === Evaluate bilinear function at interior points ===
    nx = 16
    allocate (x1a(nx), x2a(nx), dx1(nx), dx2(nx), dx1_ok(nx), dx2_ok(nx))
    call linspace (x1a, x1bnd(1)-10.0d0, x1bnd(2)+10.0d0)
    call linspace (x2a, x2bnd(1)-10.0d0, x2bnd(2)+10.0d0)

    ! Evaluate true values
    call fcn1 (x1a, x2a, dx1=dx1_ok, dx2=dx2_ok)
    ! Evaluate approx
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, dx1, m=1, dim=1, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dx1, dx1_ok, atol=1.0d-10))) &
        error stop "Knots: [2, 2]; eval derivative for dim=1 at non-interior points"

    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, dx2, m=1, dim=2, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dx2, dx2_ok, atol=1.0d-10))) &
        error stop "Knots: [2, 2]; eval derivative for dim=2 at non-interior points"

    deallocate (x1a, x2a, dx1, dx2, dx1_ok, dx2_ok)
    deallocate (knots, coefs)

    ! === Fit to linear function ===
    n1 = 11
    n2 = 7
    n = [n1, n2]
    k = 1
    x1bnd = [-10.0_PREC, 5.0_PREC]
    x2bnd = [-1.0_PREC, 1.234_PREC]
    allocate (x1in(n1), x2in(n2), yin(n1,n2))

    nknots = ppoly_get_nknots (pp, n, k)
    ncoefs = ppoly_get_ncoefs (pp, n, k)

    allocate (knots(nknots), coefs(ncoefs))

    call linspace (x1in, x1bnd(1), x1bnd(2))
    call linspace (x2in, x2bnd(1), x2bnd(2))

    do i = 1, n2
        call fcn2 (x1in, x2in(i), yin(:,i))
    end do

    status = NF_STATUS_UNDEFINED
    call ppolyfit (pp, x1in, x2in, yin, k, knots, coefs, status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Knots: [11, 7]; fit to linear polynomial"
    deallocate (x1in, x2in, yin)

    ! === Evaluate linear function at interior points ===
    nx = 11
    allocate (x1a(nx), x2a(nx), dx1(nx), dx2(nx), dx1_ok(nx), dx2_ok(nx))
    call linspace (x1a, x1bnd(1)+0.01d0, x1bnd(2)-0.01d0)
    call linspace (x2a, x2bnd(1)+0.01d0, x2bnd(2)-0.01d0)

    ! Evaluate true values
    call fcn2 (x1a, x2a, dx1=dx1_ok, dx2=dx2_ok)
    ! Evaluate approx
    status = NF_STATUS_UNDEFINED
    call ppolyval (pp, knots, coefs, x1a, x2a, dx1, m=1, dim=1, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dx1, dx1_ok, atol=1.0d-10))) &
        error stop "Knots: [11, 7]; eval derivative for dim=1 at interior points"

    call ppolyval (pp, knots, coefs, x1a, x2a, dx2, m=1, dim=2, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dx2, dx2_ok, atol=1.0d-10))) &
        error stop "Knots: [11, 7]; eval derivative for dim=2 at interior points"

    deallocate (x1a, x2a, dx1, dx2, dx1_ok, dx2_ok)

end subroutine

elemental subroutine fcn1 (x1, x2, y, dx1, dx2)
    real (PREC), intent(in) :: x1, x2
    real (PREC), intent(out), optional :: y
    real (PREC), intent(out), optional :: dx1
    real (PREC), intent(out), optional :: dx2

    real (PREC), parameter :: COEFS(0:*) = [1.23d0, 4.56d0, -7.234d0, 2.345d0]

    if (present(y)) then
        y = COEFS(0) + COEFS(1)*x1 + COEFS(2)*x2 + COEFS(3)*x1*x2
    end if

    if (present(dx1)) then
        dx1 = COEFS(1) + COEFS(3)*x2
    end if

    if (presenT(dx2)) then
        dx2 = COEFS(2) + COEFS(3)*x1
    end if
end subroutine

elemental subroutine fcn2 (x1, x2, y, dx1, dx2)
    !*  Linear test function
    real (PREC), intent(in) :: x1, x2
    real (PREC), intent(out), optional :: y
    real (PREC), intent(out), optional :: dx1, dx2

    real (PREC), parameter :: COEFS(0:*) = [-.123D0, -4.56d0, 2.345d0]

    if (present(y)) then
        y = COEFS(0) + COEFS(1)*x1 + COEFS(2)*x2
    end if

    if (present(dx1)) dx1 = COEFS(1)
    if (present(dx2)) dx2 = COEFS(2)
end subroutine

elemental subroutine fcn3 (x1, x2, y)
    !*  Function that cannot be exactly approximated by bilinear polynomial
    real (PREC), intent(in) :: x1, x2
    real (PREC), intent(out) :: y

    y = -1.23d0 + 2.56d0*x1 - 3.234d0 * x2 + 2.345 * x1 * x2 &
        - 0.234d0*x1**2.0 + 2.345d0 * x2**2.0
end subroutine

end program
