

program test_interp_pchip
    !*  Unit tests for Piecewise Cubic Hermite Interpolating Polynomials
    !   module.

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_core
    use numfort_common, only: workspace => workspace_real64
    use numfort_common_testing
    use numfort_interpolate

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_input ()
    call test_interp ()
    call test_extrap ()

end subroutine

subroutine test_input ()

    real (PREC), dimension(3) :: x, y, xp
    real (PREC), dimension(:), allocatable :: coef
    integer :: n, ncoef

    type (status_t) :: status

    ! FIT: Check with too small input data arrays
    n = 1
    ncoef = interp_pchip_get_ncoef (n)
    allocate (coef(ncoef))
    call interp_pchip_fit (x(1:n), y(1:n), coef, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "FIT: Insufficient input data"
    deallocate (coef)

    ! FIT: Different-size input arrays
    n = 2
    ncoef = interp_pchip_get_ncoef(n)
    allocate (coef(ncoef))
    status = NF_STATUS_OK
    call interp_pchip_fit (x(1:n), y(1:1), coef, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "FIT: different-size input data arrays"
    deallocate (coef)

    ! FIT: COEF output array too small
    n = 3
    ncoef = interp_pchip_get_ncoef (n-1)
    allocate (coef(ncoef))
    status = NF_STATUS_OK
    call interp_pchip_fit (x(1:n), y(1:n), coef, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "FIT: COEF output array too small"
    deallocate (coef)

    ! EVAL: incompatible input array size
    status = NF_STATUS_OK
    n = 3
    ncoef = interp_pchip_get_ncoef (n)
    allocate (coef(ncoef))
    call interp_pchip_eval (xp, coef, x(1:n), y(1:n-1), status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "EVAL: different sizes of input arrays X, Y"
    deallocate (coef)

    ! EVAL: XP array too small
    n = 3
    ncoef = interp_pchip_get_ncoef (n)
    allocate (coef(ncoef))
    status = NF_STATUS_OK
    call interp_pchip_eval (xp(1:1), coef, x, y, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "EVAL: XP input array too small"
    deallocate (coef)

    ! EVAL: XP and COEF not compatible
    ncoef = interp_pchip_get_ncoef (size(xp)-1)
    allocate (coef(ncoef))
    status = NF_STATUS_OK
    call interp_pchip_eval (xp, coef, x, y, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "EVAL: XP and COEF input arrays incompatible"
    deallocate (coef)

    ! EVAL: invalid ORDER argument
    n = 3
    ncoef = interp_pchip_get_ncoef (n)
    allocate (coef(ncoef))
    status = NF_STATUS_OK
    call interp_pchip_eval (xp, coef, x, y, order=3, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "EVAL: ORDER argument larger than admissible range"

    status = NF_STATUS_OK
    call interp_pchip_eval (xp, coef, x, y, order=-1, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "EVAL: ORDER argument smaller than admissible range"

    ! EVAL: LEFT, RIGHT present if EXT=NF_INTERP_EVAL_CONST
    status = NF_STATUS_OK
    call interp_pchip_eval (xp, coef, x, y, ext=NF_INTERP_EVAL_CONST, &
        right=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "EVAL: missing LEFT when EXT=NF_INTERP_EVAL_CONST"

    status = NF_STATUS_OK
    call interp_pchip_eval (xp, coef, x, y, ext=NF_INTERP_EVAL_CONST, &
        left=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "EVAL: missing RIGHT when EXT=NF_INTERP_EVAL_CONST"

end subroutine

subroutine test_interp ()
    !*  Unit tests for interpolating various functions at interior points.

    type (status_t) :: status
    type (workspace) :: work
    integer :: n, ncoef, np, i, j
    real (PREC), dimension(:), allocatable :: x, y, xp, yp, coef, y_ok, yp_ok
    real (PREC), parameter :: TOL = 1.0d-7

    real (PREC), parameter :: SCIPY_TEST2(*,0:*) = reshape([ real(PREC) :: &
            2.        ,  1.        ,  0.        ,  1.        ,  2., &
           -1.        , -1.        ,  0.        ,  1.        ,  1., &
            0.        ,  0.        ,  5.3333333333d0,  0.        ,  0. &
        ], shape=[5,3])

    real (PREC), parameter :: SCIPY_TEST3(*,0:*) = reshape([ real(PREC) :: &
             0.48775   ,  0.048     ,  1.20815278d0,  4., &
            -1.3475    ,  0.42      ,  2.20986111d0,  3.97222222d0, &
             2.85      ,  1.2       ,  2.39166667d0,  2.47619048d0 &
        ], shape=[4, 3])

    real (PREC), parameter :: SCIPY_TEST4(*,0:*) = reshape([ real(PREC) :: &
            3.25      ,  27.51923077d0, 125.12362637d0, &
            3.25      ,  30.01923077d0,  78.66208791d0, &
            1.5       ,  16.96153846d0,  29.75274725d0 &
        ], shape=[3,3])

    real (PREC), parameter :: SCIPY_TEST5(*,0:*) = reshape([ real(PREC) :: &
        -240.        , -120.        ,  -48.        ,  -12.        ,    0., &
         147.83688817d0, 93.83510204d0, 51.83115207d0, 21.81913043d0, 3.72705882d0, &
         -55.10396055d0,-43.04764988d0,-30.92091382d0,-18.519595d0  ,-2.54117647d0 &
        ], shape=[5,3])

    real (PREC), parameter :: SCIPY_TEST6(*,0:*) = reshape([ real(PREC) :: &
        -1.9807561532d+01,  4.8000005816d-01, &
         2.4447860565d+02,  8.0003569999d-03, &
        -1.4953653148d+03, -4.8480020624d-03 &
        ], shape=[2,3])

    real (PREC), parameter :: SCIPY_TEST7(*,0:*) = reshape([ real(PREC) :: &
        -1.43700000d+00, -1.21650000d+00, -9.96000000d-01, -7.75500000d-01, -5.55000000d-01, &
         2.10000000d-01,  2.10000000d-01,  2.10000000d-01,  2.10000000d-01,  2.10000000d-01, &
         0.0, 0.0, 0.0, 0.0, 0.0 &
        ], shape=[5,3])

    ! Test 1: linear monotinic function
    np = 10
    n = 5
    ncoef = interp_pchip_get_ncoef (np)
    allocate (xp(np), yp(np), yp_ok(np), coef(ncoef), x(n), y(n), y_ok(n))
    call linspace (xp, 1.0_PREC, 10.0_PREC)
    call fcn1 (xp, yp, order=0)
    call linspace (x, 2.0_PREC, 7.0_PREC)
    call interp_pchip_fit (xp, yp, coef, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "FIT: linear function"

    call interp_pchip_eval (xp, coef, xp, yp_ok, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (yp, yp_ok, atol=TOL))) &
        error stop "EVAL: Linear function, interpolated values at XP"

    do i = 0, 2
        call fcn1 (x, y_ok, order=i)
        call interp_pchip_eval (xp, coef, x, y, order=i, status=status)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, y_ok, atol=TOL))) &
            error stop "PCHIP check failed: Linear function"
    end do

    deallocate (xp, yp, yp_ok, coef, x, y, y_ok)

    ! Test 2: f(x) = abs(x)
    np = 9
    n = 5
    ncoef = interp_pchip_get_ncoef (np)
    allocate (xp(np), yp(np), yp_ok(np), coef(ncoef), x(n), y(n), y_ok(n))
    ! Note: XP includes 0.0
    call linspace (xp, -3.0_PREC, 3.0_PREC)
    call fcn2 (xp, yp, order=0)
    call linspace (x, -2.0_PREC, 2.0_PREC)
    call interp_pchip_fit (xp, yp, coef, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "FIT: f(x) = abs(x)"

    call interp_pchip_eval (xp, coef, xp, yp_ok, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (yp, yp_ok, atol=TOL))) &
        error stop "EVAL: p(x) ~ abs(x), check interpolated values at XP"

    do i = 0, 2
        call interp_pchip_eval (xp, coef, x, y, order=i, status=status)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, SCIPY_TEST2(:,i) , &
            atol=TOL))) &
            error stop "PCHIP check failed: p(x) ~ abs(x)"
    end do

    deallocate (xp, yp, yp_ok, coef, x, y, y_ok)

    ! Test 3: f(x) = x**2.0
    np = 13
    n = 4
    ncoef = interp_pchip_get_ncoef (np)
    allocate (xp(np), yp(np), yp_ok(np), coef(ncoef), x(n), y(n), y_ok(n))
    call linspace (xp, -1.0_PREC, 3.0_PREC)
    call fcn3 (xp, yp, order=0)
    call linspace (x, -0.7_PREC, 2.0_PREC)
    call interp_pchip_fit (xp, yp, coef, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "FIT: f(x) = x^2"

    call interp_pchip_eval (xp, coef, xp, yp_ok, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (yp, yp_ok, atol=TOL))) &
        error stop "EVAL: p(x) ~ x^2, check interpolated values at XP"

    do i = 0, 2
        call interp_pchip_eval (xp, coef, x, y, order=i, status=status)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, SCIPY_TEST3(:,i) , &
            atol=TOL))) &
            error stop "PCHIP check failed: p(x) ~ x^2"
    end do

    deallocate (xp, yp, yp_ok, coef, x, y, y_ok)

    ! Test 4: f(x) = x^3
    np = 11
    n = 3
    ncoef = interp_pchip_get_ncoef (np)
    allocate (xp(np), yp(np), yp_ok(np), coef(ncoef), x(n), y(n), y_ok(n))
    call linspace (xp, -10.0_PREC, 10.0_PREC)
    call fcn4 (xp, yp, order=0)
    call linspace (x, 1.0_PREC, 5.0_PREC)
    call interp_pchip_fit (xp, yp, coef, status=status, work=work)
    if (.not. (status == NF_STATUS_OK)) error stop "FIT: f(x) = x^3"

    call interp_pchip_eval (xp, coef, xp, yp_ok, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (yp, yp_ok, atol=TOL))) &
        error stop "EVAL: p(x) ~ x^3, check interpolated values at XP"

    do i = 0, 2
        call interp_pchip_eval (xp, coef, x, y, order=i, status=status)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, SCIPY_TEST4(:,i) , &
            atol=TOL))) &
            error stop "PCHIP check failed: p(x) ~ x^3"
    end do

    deallocate (xp, yp, yp_ok, coef, x, y, y_ok)

    ! Test 5: f(x) = 2*(x^3-x)
    np = 101
    n = 5
    ncoef = interp_pchip_get_ncoef (np)
    allocate (xp(np), yp(np), yp_ok(np), coef(ncoef), x(n), y(n), y_ok(n))
    call linspace (xp, -10.0_PREC, 10.0_PREC)
    call fcn5 (xp, yp, order=0)
    call linspace (x, -5.0_PREC, -1.0_PREC)
    call interp_pchip_fit (xp, yp, coef, status=status, work=work)
    if (.not. (status == NF_STATUS_OK)) error stop "FIT: f(x) = 2(x^3-x)"

    call interp_pchip_eval (xp, coef, xp, yp_ok, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (yp, yp_ok, atol=TOL))) &
        error stop "EVAL: p(x) ~ 2(x^3-x), check interpolated values at XP"

    do i = 0, 2
        call interp_pchip_eval (xp, coef, x, y, order=i, status=status)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, SCIPY_TEST5(:,i) , &
            atol=TOL))) &
            error stop "PCHIP check failed: p(x) ~ 2(x^3-x)"
    end do

    deallocate (xp, yp, yp_ok, coef, x, y, y_ok)

    ! Test 6: CRRA utility
    np = 51
    n = 2
    ncoef = interp_pchip_get_ncoef (np)
    allocate (xp(np), yp(np), yp_ok(np), coef(ncoef), x(n), y(n), y_ok(n))
    call linspace (xp, 0.1_PREC, 10.0_PREC)
    call fcn6 (xp, yp, order=0)
    call linspace (x, 0.2_PREC, 5.0_PREC)
    call interp_pchip_fit (xp, yp, coef, status=status, work=work)
    if (.not. (status == NF_STATUS_OK)) error stop "FIT: CRRA utility"

    call interp_pchip_eval (xp, coef, xp, yp_ok, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (yp, yp_ok, atol=TOL))) &
        error stop "EVAL: CRRA utility, check interpolated values at XP"

    do i = 0, 2
        call interp_pchip_eval (xp, coef, x, y, order=i, status=status)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, SCIPY_TEST6(:,i) , atol=TOL))) &
            error stop "PCHIP check failed: CRRA utility"
    end do

    deallocate (xp, yp, yp_ok, coef, x, y, y_ok)

    ! Test 7: Check that linear function was fit if only 2 data points are
    ! provided.
    np = 2
    n = 5
    ncoef = interp_pchip_get_ncoef (np)
    allocate (xp(np), yp(np), yp_ok(np), coef(ncoef), x(n), y(n), y_ok(n))
    call linspace (xp, 0.5_PREC, 10.0_PREC)
    call fcn6 (xp, yp, order=0)
    call linspace (x, 0.8_PREC, 5.0_PREC)
    call interp_pchip_fit (xp, yp, coef, status=status, work=work)
    if (.not. (status == NF_STATUS_OK)) error stop "FIT: CRRA utility (2 points)"

    call interp_pchip_eval (xp, coef, xp, yp_ok, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (yp, yp_ok, atol=TOL))) &
        error stop "EVAL: CRRA utility (2 points), check interpolated values at XP"

    do i = 0, 2
        call interp_pchip_eval (xp, coef, x, y, order=i, status=status)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, SCIPY_TEST7(:,i) , atol=TOL))) &
            error stop "PCHIP check failed: CRRA utility (2 points)"
    end do

    ! Test 8: Test scalar interface using same setup as in test 7.

    ! Evaluate using scalar interface
    do j = 1, size(yp)
        call interp_pchip_eval (xp, coef, xp(j), yp_ok(j), status=status)
    end do
    if (.not. (status == NF_STATUS_OK .and. all_close (yp, yp_ok, atol=TOL))) &
        error stop "EVAL (scalar): CRRA utility (2 points), check interpolated " // &
            "values at XP"

    do i = 0, 2
        ! Evaluate using scalar interface
        do j = 1, size(y)
            call interp_pchip_eval (xp, coef, x(j), y(j), order=i, status=status)
        end do
        if (.not. (status == NF_STATUS_OK .and. all_close (y, SCIPY_TEST7(:,i) , &
            atol=TOL))) &
            error stop "PCHIP scalar check failed: CRRA utility (2 points)"
    end do

    deallocate (xp, yp, yp_ok, coef, x, y, y_ok)

end subroutine

subroutine test_extrap ()

    type (status_t) :: status
    type (workspace) :: work
    integer :: n, ncoef, np, i
    real (PREC), dimension(:), allocatable :: x, y, xp, yp, coef, y_ok, yp_ok
    real (PREC), parameter :: TOL = 1.0d-7
    real (PREC) :: left, right

    real (PREC), parameter :: SCIPY_TEST5(*,0:*) = reshape([ real(PREC) :: &
             3.21428571,  0.01234568,  3.21428571, &
            -1.78571429,  0.        ,  1.78571429, &
            -2.14285714,  0.        , -2.14285714 &
        ], shape=[3,3])

    real (PREC), parameter :: SCIPY_TEST6(*,0:*) = reshape([ real(PREC) :: &
            -3.2799424026d+01,  0.0,  3.2799424026d+01, &
             1.3643273611d+01, -2.3317027720d+00,  1.3643273611d+01, &
             2.7803696622d+01, 0.0, -2.7803696622d+01 &
        ], shape=[3,3])

    ! Test 1: linear monotinic function, extrapolation
    np = 10
    n = 5
    ncoef = interp_pchip_get_ncoef (np)
    allocate (xp(np), yp(np), yp_ok(np), coef(ncoef), x(n), y(n), y_ok(n))
    call linspace (xp, 1.0_PREC, 10.0_PREC)
    call fcn1 (xp, yp, order=0)
    call linspace (x, 0.0_PREC, 15.0_PREC)
    call interp_pchip_fit (xp, yp, coef, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "FIT: linear function"

    do i = 0, 2
        call fcn1 (x, y_ok, order=i)
        call interp_pchip_eval (xp, coef, x, y, order=i, status=status, &
            ext=NF_INTERP_EVAL_EXTRAPOLATE)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, y_ok, atol=TOL))) &
            error stop "PCHIP check failed: Linear function, extrapolate"
    end do

    ! Test 2: reuse fitted polynomial from above, replace non-interior points
    ! with const
    left = -1.0d2
    right = 1.0d2

    do i = 0, 2
        call fcn1 (x, y_ok, order=i)
        ! replace function values at non-interior points
        where (x > xp(np))
            y_ok = right
        else where (x < xp(1))
            y_ok = left
        end where

        call interp_pchip_eval (xp, coef, x, y, order=i, status=status, &
            ext=NF_INTERP_EVAL_CONST, left=left, right=right)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, y_ok, atol=TOL))) &
            error stop "PCHIP check failed: Linear function, return const"
    end do

    ! Test 3:  reuse fitted polynomial from above, replace non-interior points
    ! with boundary values
    do i = 0, 2
        call fcn1 (x, y_ok, order=i)
        call fcn1 (xp, yp_ok, order=i)
        ! replace function values at non-interior points
        where (x > xp(np))
            y_ok = yp_ok(np)
        else where (x < xp(1))
            y_ok = yp_ok(1)
        end where

        call interp_pchip_eval (xp, coef, x, y, order=i, status=status, &
            ext=NF_INTERP_EVAL_BOUNDARY)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, y_ok, atol=TOL))) &
            error stop "PCHIP check failed: Linear function, use boundary"
    end do

    ! Test 4: reuse fitted polynomial from above, raise error on non-interior
    do i = 0, 2
        call interp_pchip_eval (xp, coef, x, y, order=i, status=status, &
            ext=NF_INTERP_EVAL_ERROR)
        if (.not. (status == NF_STATUS_BOUNDS_ERROR)) &
            error stop "PCHIP check failed: Linear function, raise error"
    end do

    deallocate (xp, yp, yp_ok, coef, x, y, y_ok)

    ! Test 5: quadratic function
    np = 10
    n = 3
    ncoef = interp_pchip_get_ncoef (np)
    allocate (xp(np), yp(np), yp_ok(np), coef(ncoef), x(n), y(n), y_ok(n))
    call linspace (xp, -1.0_PREC, 1.0_PREC)
    call fcn3 (xp, yp, order=0)
    call linspace (x, -2.0_PREC, 2.0_PREC)
    call interp_pchip_fit (xp, yp, coef, work=work, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "FIT: f(x) = x^2"

    do i = 0, 2
        call interp_pchip_eval (xp, coef, x, y, order=i, status=status, &
            ext=NF_INTERP_EVAL_EXTRAPOLATE)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, SCIPY_TEST5(:,i), &
            atol=TOL))) &
            error stop "PCHIP check failed: p(x) ~ x^2, extrapolate"
    end do

    deallocate (xp, yp, yp_ok, coef, x, y, y_ok)

    ! Test 6: cubic function
    np = 10
    n = 3
    ncoef = interp_pchip_get_ncoef (np)
    allocate (xp(np), yp(np), yp_ok(np), coef(ncoef), x(n), y(n), y_ok(n))
    call linspace (xp, -2.0_PREC, 2.0_PREC)
    call fcn5 (xp, yp, order=0)
    call linspace (x, -3.0_PREC, 3.0_PREC)
    call interp_pchip_fit (xp, yp, coef, work=work, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "FIT: f(x) = 2(x^3 - x)"

    do i = 0, 2
        call interp_pchip_eval (xp, coef, x, y, order=i, status=status, &
            ext=NF_INTERP_EVAL_EXTRAPOLATE)
        if (.not. (status == NF_STATUS_OK .and. all_close (y, SCIPY_TEST6(:,i), atol=TOL))) &
            error stop "PCHIP check failed: p(x) ~ 2(x^3 - x), extrapolate"
    end do

    deallocate (xp, yp, yp_ok, coef, x, y, y_ok)

end subroutine

elemental subroutine fcn1 (x, y, order)
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer, intent(in) :: order

    select case (order)
    case (0)
        y = x
    case (1)
        y = 1.0
    case (2)
        y = 0.0
    end select
end subroutine

elemental subroutine fcn2 (x, y, order)
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer, intent(in) :: order

    select case (order)
    case (0)
        y = abs(x)
    case (1)
        y = signum (x)
    case (2)
        y = 0.0
    end select
end subroutine

elemental subroutine fcn3 (x, y, order)
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer, intent(in) :: order

    select case (order)
    case (0)
        y = x ** 2.0
    case (1)
        y = 2.0 * x
    case (2)
        y = 2.0
    end select
end subroutine

elemental subroutine fcn4 (x, y, order)
    !*  Monotone cubic function
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer, intent(in) :: order

    select case (order)
    case (0)
        y = x ** 3.0
    case (1)
        y = 3.0 * x ** 2.0
    case (2)
        y = 6.0
    end select
end subroutine

elemental subroutine fcn5 (x, y, order)
    !*  Non-monotone cubic function
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer, intent(in) :: order

    select case (order)
    case (0)
        y = 2.0 * x ** 3.0 - 2.0 * x
    case (1)
        y = 6.0 * x ** 2.0 - 2.0
    case (2)
        y = 12.0 * x
    end select
end subroutine

elemental subroutine fcn6 (x, y, order)
    !*  Monotone power function
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer, intent(in) :: order

    real (PREC), parameter :: SIGMA = 3.0

    select case (order)
    case (0)
        y = (x ** (1.0-SIGMA) - 1.0) / (1.0 - SIGMA)
    case (1)
        y = x ** (-SIGMA)
    case (2)
        y = x ** (-1.0 - SIGMA) / (-SIGMA)
    end select
end subroutine

end
