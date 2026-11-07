
program test_interpolate_concon

    use, intrinsic :: iso_fortran_env

    use numfort_arrays, only: linspace
    use numfort_interpolate, workspace => workspace_real64

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_cubic ()

end subroutine

subroutine test_cubic ()
    !*   Unit tests for cubic splines

    integer, parameter :: m = 20, ms = 50
    real (PREC), dimension(:), allocatable :: x, fx, fpx, fppx, v
    real (PREC), dimension(:), allocatable :: fx_s, fpx_s, fppx_s
    real (PREC), dimension(:), allocatable :: knots, coefs
    type (workspace) :: work
    real (PREC) :: xlb, xub
    integer, parameter :: k = 3, ik = 2
    real (PREC) :: s, ssr, diff
    type (status_t) :: status
    integer :: nest, n
    integer, parameter :: iopt = 0

    xlb = -2.0
    xub = 2.0

    allocate (x(m), fx(m), v(m))
    call linspace (x, xlb, xub)
    call fcn (x, ik, fx)

    nest = concon_get_nest (m)
    allocate (knots(nest), coefs(nest))

    ! Fit interpolating spline
    s = 1d-10
    v(:) = 1.0
    call concon (x, fx, v, s, knots, coefs, n, iopt=iopt, work=work, ssr=ssr, &
        status=status)

    ! Check spline values and derivatives
    deallocate (x, fx)

    allocate (fx_s(ms), fpx_s(ms), fppx_s(ms))
    allocate (x(ms), fx(ms), fpx(ms), fppx(ms))

    call linspace (x, xlb, xub)

    call fcn (x, ik, fx, fpx, fppx)

    call splev (knots(1:n), coefs(1:n), k, x, fx_s, status=status)
    diff = maxval(abs(fx_s-fx))
    if (.not. (diff < 1e-10)) &
        error stop "Polynomial of degree 2: sup norm (s(x)-f(x)) < 1e-10"

    call splder (knots(1:n), coefs(1:n), k, 1, x, fpx_s, work=work, status=status)
    diff = maxval(abs(fpx-fpx_s))
    if (.not. (diff < 1e-10)) &
        error stop "Polynomial of degree 2: sup norm (s'(x)-f'(x)) < 1e-10"

    call splder (knots(1:n), coefs(1:n), k, 2, x, fppx_s, work=work, status=status)
    diff = maxval(abs(fppx-fppx_s))
    if (.not. (diff < 1e-10)) &
        error stop "Polynomial of degree 2: sup norm (s''(x)-f''(x)) < 1e-10"

    deallocate (x, fx, fpx, fppx)
    deallocate (fx_s, fpx_s, fppx_s)
    deallocate (knots, coefs)

end subroutine

elemental subroutine fcn (x, k, fx, fpx, fppx)
    !*  Polynomial test function (guaranteed to be concave for k=2)
    real (PREC), intent(in) :: x
    integer, intent(in) :: k
        !*  Polynomial degree to use
    real (PREC), intent(out), optional :: fx, fpx, fppx

    integer :: i
    real (PREC), parameter :: coefs(0:*) = &
        [-0.123d0, 1.234d0, -2.345d0, 3.456d0, 3.4583d0, -0.123d0]
        !   Polynomial coefficients in increasing order

    if (present(fx)) then
        fx = 0.0
        do i = 0, k
            fx = fx + coefs(i) * x**i
        end do
    end if

    if (present(fpx)) then
        fpx = 0.0
        do i = 1, k
            fpx = fpx + coefs(i) * i * x**(i-1)
        end do
    end if

    if (present(fppx)) then
        fppx = 0.0
        do i = 2, k
            fppx = fppx + coefs(i) * i * (i-1) * x**(i-2)
        end do
    end if
end subroutine

end program
