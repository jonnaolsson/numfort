

program test_optimize_root_broyden

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_common
    use numfort_common_testing
    use numfort_common_workspace, workspace => workspace_real64

    use numfort_optimize, only: root_broyden, optim_result => optim_result_real64

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_quad ()
    call test_func2 ()

    call test_maxiter ()
    call test_invalid_residual ()

end subroutine

subroutine test_quad ()

    type (optim_result) :: res
    type (workspace) :: work
    real (PREC), dimension(:), allocatable :: x, fx, x0
    integer :: n

    n = 2
    allocate (x(n), fx(n), x0(n))
    x0 = [-2.0d0, 10.0d0]
    x = x0
    call root_broyden (quad_fcn_jac, x, tol=1.0d-8, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Quadratic function; x0=[-2.000, 10.000]"
    deallocate (x, fx, x0)

    ! User-provided WORKSPACE object
    n = 2
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0]
    x = x0
    call root_broyden (quad_fcn_jac, x, tol=1.0d-8, work=work, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Quadratic function; WORK provided; x0=[4.000, 1.000]"
    deallocate (x, fx, x0)

    ! Limit absolute step size
    n = 2
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0]
    x = x0
    call root_broyden (quad_fcn_jac, x, tol=1.0d-8, xstep=0.1d0, work=work, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Quadratic function; WORK, XSTEP provided; x0=[4.000, 1.000]"
    deallocate (x, fx, x0)

    ! Limit relative step size
    n = 2
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, -10.0d0]
    x = x0
    call root_broyden (quad_fcn_jac, x, tol=1.0d-8, rstep=0.1d0, work=work, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Quadratic function; WORK, RSTEP provided; x0=[4.000, -10.000]"
    deallocate (x, fx, x0)

end subroutine


subroutine test_func2 ()

    type (optim_result) :: res
    real (PREC), dimension(:), allocatable :: x, fx, x0
    integer :: n

    ! Separate FCN and JAC routines
    n = 3
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0, 1.0d0]
    x = x0
    call root_broyden (func2_fcn, func2_jac, x, tol=1.0d-8, xstep=0.1d0, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Separate FCN, JAC; XSTEP provided; x0=[4.000, 1.000, 1.000]"
    deallocate (x, fx, x0)

    ! Numerical differentiation
    n = 3
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0, 1.0d0]
    x = x0
    call root_broyden (func2_fcn, x, ndiff=.true., tol=1.0d-8, dstep=1.0d-5, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Numerical diff.; x0=[4.000, 1.000, 1.000]"
    deallocate (x, fx, x0)

end subroutine


subroutine test_maxiter ()
    !*  Test that

    type (optim_result) :: res
    real (PREC), dimension(:), allocatable :: x, fx, x0
    integer :: n
    logical :: lres

    ! Separate FCN and JAC routines
    n = 3
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0, 1.0d0]
    x = x0
    call root_broyden (func2_fcn, func2_jac, x, tol=1.0d-8, xstep=0.1d0, &
        maxiter=1, res=res)
    lres = (.not. res%success) .and. res%status == NF_STATUS_MAX_ITER .and. &
        res%nit == 1
    if (.not. (lres)) error stop "Separate FCN, JAC; MAXITER=1"

    x = x0
    call root_broyden (func2_fcn, func2_jac, x, tol=1.0d-8, xstep=0.1d0, &
        maxfev=1, res=res)
    lres = (.not. res%success) .and. res%status == NF_STATUS_MAX_EVAL .and. &
        res%nit == 1
    if (.not. (lres)) error stop "Separate FCN, JAC; MAXFUN=1"

    deallocate (x, fx, x0)

    ! Numerical differentiation
    n = 3
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0, 1.0d0]
    x = x0
    call root_broyden (func2_fcn, x, ndiff=.true., tol=1.0d-8, dstep=1.0d-5, &
        maxiter=1, res=res)
    lres = (.not. res%success) .and. res%status == NF_STATUS_MAX_ITER .and. &
        res%nit == 1
    if (.not. (lres)) error stop "Numerical diff.; MAXITER=1"

    x = x0
    call root_broyden (func2_fcn, x, ndiff=.true., tol=1.0d-8, dstep=1.0d-5, &
        maxfev=1, res=res)
    lres = (.not. res%success) .and. res%status == NF_STATUS_MAX_EVAL .and. &
        res%nit == 1
    if (.not. (lres)) error stop "Numerical diff.; MAXFUN=1"

    deallocate (x, fx, x0)

end subroutine

subroutine quad_fcn_jac (x, fx, fpx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: fpx

    if (present(fx)) then
        fx(1) = (x(1)-2.0)**2.0
        fx(2) = (x(2)+1.0)**2.0
    end if

    if (present(fpx)) then
        fpx = 0.0
        fpx(1,1) = 2.0 * (x(1) - 2.0)
        fpx(2,2) = 2.0 * (x(2) + 1.0)
    end if
end subroutine


subroutine func2_fcn (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    fx(1) = sin(x(1))
    fx(2) = (x(2) - 2.0) ** 2.0
    fx(3) = cos(x(1)) * x(3)
end subroutine


subroutine func2_jac (x, fpx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:,:), contiguous :: fpx

    fpx = 0.0

    fpx(1,1) = cos(x(1))
    fpx(2,2) = 2.0 * (x(2) - 2.0)
    fpx(3,1) = -sin(x(1)) * x(3)
    fpx(3,3) = cos(x(1))
end subroutine


subroutine test_invalid_residual ()
    !*  Test handling of non-finite (NaN and Inf) line-search residuals,
    !   workspace reuse, and recovery from non-finite trial steps.

    type (optim_result) :: res
    type (workspace) :: work
    real (PREC), dimension(1) :: x, fx_eval
    real (PREC), dimension(1,1) :: jac_inv_init
    logical :: lres

    jac_inv_init(1,1) = 1.0_PREC

    ! Case 1: All line-search trials return NaN with fresh workspace.
    ! At initial point x=0, residual is f(0)=1. Away from x=0, callback returns NaN.
    ! Line search has no finite admissible trial, so the routine must report failure
    ! (NF_STATUS_INVALID_STATE) without reporting false convergence.
    ! Returned point must preserve initial x=0 and residual f(x)=1.
    x = [0.0_PREC]
    call root_broyden (nan_residual_fcn, x, ndiff=.true., &
        jac_inv_init=jac_inv_init, maxiter=3, res=res)
    call nan_residual_fcn (res%x, fx_eval)
    lres = (.not. res%success) .and. (res%status == NF_STATUS_INVALID_STATE) &
        .and. (res%x(1) == 0.0_PREC) .and. (res%fx(1) == fx_eval(1)) &
        .and. (fx_eval(1) == 1.0_PREC)
    if (.not. lres) error stop "All-NaN residuals; fresh workspace"

    ! Case 2: All line-search trials return NaN with reused workspace.
    ! Scratch memory is pre-populated with nonzero values to ensure that
    ! untouched scratch memory is never misinterpreted as converged residual.
    call assert_alloc (work, nrwrk=100)
    work%rwrk = 42.0_PREC
    x = [0.0_PREC]
    call root_broyden (nan_residual_fcn, x, ndiff=.true., &
        jac_inv_init=jac_inv_init, maxiter=3, work=work, res=res)
    call nan_residual_fcn (res%x, fx_eval)
    lres = (.not. res%success) .and. (res%status == NF_STATUS_INVALID_STATE) &
        .and. (res%x(1) == 0.0_PREC) .and. (res%fx(1) == fx_eval(1)) &
        .and. (fx_eval(1) == 1.0_PREC)
    if (.not. lres) error stop "All-NaN residuals; reused workspace"

    ! Case 3: All line-search trials return positive infinity with fresh workspace.
    ! At initial point x=0, residual is f(0)=1. Away from x=0, callback returns +Inf.
    ! Line search must report failure rather than treating scratch space as converged.
    x = [0.0_PREC]
    call root_broyden (inf_residual_fcn, x, ndiff=.true., &
        jac_inv_init=jac_inv_init, maxiter=3, res=res)
    call inf_residual_fcn (res%x, fx_eval)
    lres = (.not. res%success) .and. (res%status == NF_STATUS_INVALID_STATE) &
        .and. (res%x(1) == 0.0_PREC) .and. (res%fx(1) == fx_eval(1)) &
        .and. (fx_eval(1) == 1.0_PREC)
    if (.not. lres) error stop "All-Inf residuals; fresh workspace"

    ! Case 4: All line-search trials return positive infinity with reused workspace.
    work%rwrk = -99.0_PREC
    x = [0.0_PREC]
    call root_broyden (inf_residual_fcn, x, ndiff=.true., &
        jac_inv_init=jac_inv_init, maxiter=3, work=work, res=res)
    call inf_residual_fcn (res%x, fx_eval)
    lres = (.not. res%success) .and. (res%status == NF_STATUS_INVALID_STATE) &
        .and. (res%x(1) == 0.0_PREC) .and. (res%fx(1) == fx_eval(1)) &
        .and. (fx_eval(1) == 1.0_PREC)
    if (.not. lres) error stop "All-Inf residuals; reused workspace"

    ! Case 5: Recovery when an earlier line search trial is invalid (NaN) but a
    ! later trial is valid and finite.
    ! With initial x=0 and initial inverse Jacobian jac_inv_init=[2.0], the full
    ! step gives x=1.0 (|x| > 0.8), producing NaN.
    ! Backtracking evaluates step_size=0.75 giving x=0.75 (|x| <= 0.8), where
    ! residual f(0.75) = 0.75 - 0.5 = 0.25 is finite and has smaller norm than f(0) = -0.5.
    ! Algorithm should recover and converge to closed-form root x = 0.5.
    jac_inv_init(1,1) = 2.0_PREC
    x = [0.0_PREC]
    call root_broyden (recovery_fcn, x, ndiff=.true., &
        jac_inv_init=jac_inv_init, tol=1.0e-8_PREC, res=res)
    call recovery_fcn (res%x, fx_eval)
    lres = res%success .and. (res%status == NF_STATUS_OK) &
        .and. (abs(res%x(1) - 0.5_PREC) < 1.0e-6_PREC) &
        .and. (res%fx(1) == fx_eval(1)) &
        .and. (abs(fx_eval(1)) < 1.0e-6_PREC)
    if (.not. lres) error stop "Recovery when earlier trial is invalid"

end subroutine


subroutine nan_residual_fcn (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    ! Residual is 1.0 at x = 0; quiet NaN everywhere else
    if (abs(x(1)) < 1.0e-12_PREC) then
        fx(1) = 1.0_PREC
    else
        fx(1) = ieee_value (0.0_PREC, ieee_quiet_nan)
    end if
end subroutine


subroutine inf_residual_fcn (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    ! Residual is 1.0 at x = 0; +infinity everywhere else
    if (abs(x(1)) < 1.0e-12_PREC) then
        fx(1) = 1.0_PREC
    else
        fx(1) = ieee_value (0.0_PREC, ieee_positive_inf)
    end if
end subroutine


subroutine recovery_fcn (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    ! Callback has domain restriction |x| <= 0.8; returns NaN outside.
    ! Target closed-form solution is x = 0.5 where f(0.5) = 0.5 - 0.5 = 0.0.
    if (abs(x(1)) > 0.8_PREC) then
        fx(1) = ieee_value (0.0_PREC, ieee_quiet_nan)
    else
        fx(1) = x(1) - 0.5_PREC
    end if
end subroutine


end program
