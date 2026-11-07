

program test_optimize_minimize_slsqp

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_optimize, workspace => workspace_real64, &
        optim_result => optim_result_real64

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_rosenbrock_scipy ()
    call test_quadratic_scipy ()

end subroutine

subroutine test_rosenbrock_scipy ()

    type (workspace) :: work
    type (optim_result) :: res

    real (PREC), parameter :: tol = 1d-8
    integer, parameter :: n = 2, m = 1
    real (PREC), dimension(n) :: x, lbounds, ubounds
    real (PREC) :: dx, dfx

    ! Results taken from Scipy's wrapper for SLSQP
    ! Problem #1
    real (PREC), parameter :: x1_scipy(2) = &
        [0.786415150971839d0, 0.617698316595411d0]
    real (PREC), parameter :: fx1_scipy = 0.045674808719160d0

    ! Problem #2
    real (PREC), parameter :: x2_scipy(2) = &
        [0.826247374068299d0, 0.682072165838992d0]
    real (PREC), parameter :: fx2_scipy = 0.030227497664672d0

    ! Problem #1: Rosenbrock with inequality constr.
    x = 0.1d0
    lbounds = -1.0d0
    ubounds = 1.0d0
    call minimize_slsqp (fobj, x, lbounds, ubounds, m=1, f_ieqcons=fconstr_ieq1, &
        work=work, res=res, tol=1d-8)

    dx = maxval(abs(x-x1_scipy))
    dfx = abs(res%fx(1)-fx1_scipy)

    if (.not. (dx < tol .and. dfx < tol)) &
        error stop "Problem 1: Rosenbrock with ineq. and box constraints"

    ! Problem #2: Rosenbrock with equality constr.
    x = 0.1d0
    lbounds = -1.0d0
    ubounds = 1.0
    call minimize_slsqp (fobj, x, lbounds, ubounds, m=1, f_eqcons=fconstr_eq2, &
        work=work, res=res, tol=1d-8)

    dx = maxval(abs(x-x2_scipy))
    dfx = abs(res%fx(1)-fx2_scipy)

    if (.not. (dx < tol .and. dfx < tol)) &
        error stop "Problem 2: Rosenbrock with eq. and box constraints"

end subroutine


subroutine test_quadratic_scipy ()

    type (workspace) :: work
    type (optim_result) :: res
    real (PREC), parameter :: tol = 1d-10
    real (PREC), dimension(11) :: x, lbounds, ubounds
    real (PREC) :: dx, dfx
    real (PREC) :: POS_INF, NEG_INF

    ! Scipy results for problem #3
    real (PREC), parameter :: x1_scipy(11) = [ &
         5.773502691896936d-01,  5.773502691896940d-01,  5.773502691896936d-01, &
        -9.765189259803117d-07,  3.082206907816453d+00,  5.000005774156433d-01, &
        -9.765189257569339d-07, -9.765189254581599d-07, -9.765189250856281d-07, &
        -9.765189253218452d-07, -9.765189256346833d-07 ]
    real (PREC), parameter :: fx1_scipy = 10.750000000013397d0

    ! Scipy results for problem #4
    real (PREC), parameter :: x2_scipy(11) = [ &
         5.773502691897339d-01,  5.773502691897370d-01,  5.773502691897344d-01, &
        -1.897200341762655d-07,  3.082206629424660d+00,  5.000022936265232d-01, &
         5.000000000000284d+00, -1.897200340142180d-07, -1.897200341364482d-07, &
        -1.897200340510266d-07, -5.000000000000293d+00 ]
    real (PREC), parameter :: fx2_scipy = 60.750000000107406d0

    ! Problem 3
    x = 2.0d0
    call minimize_slsqp (fobj3, x, m=2, f_ieqcons=fconstr_ieq3, &
        work=work, res=res, tol=1d-8)

    dx = maxval(abs(x-x1_scipy))
    dfx = abs(res%fx(1)-fx1_scipy)

    if (.not. (dx < tol .and. dfx < tol)) &
        error stop "Problem 3: Quadratic obj. with ineq. constraints"

    ! Problem 4: Same as problem 3, but impose additional box constraints
    ! on some of the dimensions of x.
    x = 2.0d0
    x(7) = 7.0d0
    x(11) = -7.0d0

    POS_INF = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    NEG_INF = ieee_value (0.0_PREC, IEEE_NEGATIVE_INF)

    lbounds = NEG_INF
    ubounds = POS_INF
    lbounds(7) = 5.0d0
    ubounds(7) = 10.0d0
    lbounds(11) = -10.0d0
    ubounds(11) = -5.0d0

    call minimize_slsqp (fobj3, x, lbounds=lbounds, ubounds=ubounds, m=2, &
        f_ieqcons=fconstr_ieq3, work=work, res=res, tol=1d-8)

    dx = maxval(abs(x-x2_scipy))
    dfx = abs(res%fx(1)-fx2_scipy)

    if (.not. (dx < tol .and. dfx < tol)) &
        error stop "Problem 4: Quadratic obj. with ineq. and partial box constraints"
end subroutine


subroutine fobj (x, fx, fpx)
    !*  Computes the value and derivative of the Rosenbrock function
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), dimension(:), optional, contiguous :: fpx

    if (present(fx)) then
        ! Compute objective
        fx = 1.0d2 * (x(2)-x(1)**2.0d0)**2.0d0 + (1.0d0-x(1))**2.0d0
    end if

    if (present(fpx)) then
        fpx(1) = - 4.0d2 * (x(2)-x(1)**2.0d0) * x(1) - 2.0d0*(1.0d0-x(1))
        fpx(2) = 2.0d2 * (x(2)-x(1)**2.0d0)
    end if
end subroutine


subroutine fconstr_ieq1 (x, fx, fpx)
    !*  Function evaluating inequality constraints
    !   Constraints needs to be formulated such that C(x) >= 0
    real (real64), intent(in), dimension(:), contiguous :: x
    real (real64), intent(out), dimension(:), contiguous, optional :: fx
    real (real64), intent(out), dimension(:,:), contiguous, optional :: fpx

    if (present(fx)) then
        fx = -x(1)**2.0d0 - x(2) ** 2.0d0 + 1.0d0
    end if

    ! Jacobian of constraint function
    if (present(fpx)) then
        fpx(1,:) = -2.0d0 * x
    end if
end subroutine


subroutine fconstr_eq2 (x, fx, fpx)
    real (real64), intent(in), dimension(:), contiguous :: x
    real (real64), intent(out), dimension(:), contiguous, optional :: fx
    real (real64), intent(out), dimension(:,:), contiguous, optional :: fpx

    if (present(fx)) then
        fx(1) = - x(1)**2.0d0 - x(2)**3.0d0 + 1.0d0
    end if

    if (present(fpx)) then
        fpx(1, 1) = -2.0d0 * x(1)
        fpx(1, 2) = -3.0d0 * x(2) ** 2.0d0
    end if
end subroutine


subroutine fobj3 (x, fx, fpx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), dimension(:), optional, contiguous :: fpx

    if (present(fx)) then
        fx = sum(x**2.0d0)
    end if

    if (present(fpx)) then
        fpx = 2.0d0 * x
    end if
end subroutine


subroutine fconstr_ieq3 (x, fx, fpx)
    real (real64), intent(in), dimension(:), contiguous :: x
    real (real64), intent(out), dimension(:), contiguous, optional :: fx
    real (real64), intent(out), dimension(:,:), contiguous, optional :: fpx

    if (present(fx)) then
        fx(1) = sum(x(1:3) ** 2.0d0) - 1.0d0
        fx(2) = x(5)**2.0d0 + x(6) - 10.0d0
    end if

    if (present(fpx)) then
        fpx(1,:) = 0.0d0
        fpx(1,1:3) = 2.0d0 * x(1:3)

        fpx(2,:) = 0.0d0
        fpx(2,5) = 2.0d0 * x(5)
        fpx(2,6) = 1.0d0
    end if
end subroutine

end program
