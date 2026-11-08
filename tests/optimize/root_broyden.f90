

program test_optimize_root_broyden

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: ieee_exceptions, only: ieee_get_flag, ieee_set_flag, &
        ieee_divide_by_zero, ieee_invalid

    use numfort_common
    use numfort_common_testing
    use numfort_common_workspace, workspace => workspace_real64

    use numfort_optimize, only: root_broyden, optim_result => optim_result_real64

    implicit none

    integer, parameter :: PREC = real64
    integer :: callback_count = 0

    call test_all ()

    contains

subroutine test_all ()

    call test_quad ()
    call test_func2 ()

    call test_maxiter ()
    call test_invalid_residual ()
    call test_update_breakdown ()
    call test_initial_root ()
    call test_maxfev ()

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
        res%nit == 0 .and. res%nfev == 1
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
        res%nit == 0 .and. res%nfev == 1
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


subroutine test_update_breakdown ()
    !*  Test detection and handling of zero, non-finite, and singular Broyden
    !   inverse-Jacobian update denominators (Issue #24).
    !
    !   Verifies clean termination with NF_STATUS_INVALID_STATE, consistency of
    !   returned x and fx, retention of uncorrupted inverse Jacobians, and
    !   absence of IEEE divide-by-zero or invalid-operation flags.

    type (optim_result) :: res
    real (PREC), dimension(1) :: x1, fx1_eval
    real (PREC), dimension(1,1) :: jac_inv1
    real (PREC), dimension(2) :: x2, fx2_eval
    real (PREC), dimension(2,2) :: jac_inv2
    logical :: lres, flag_div_zero, flag_invalid

    ! Case 1: Deterministic constant residual (Issue #24 reproduction)
    ! f(x) = 1 everywhere. Initial x=0, initial inverse Jacobian [1.0], xstep=0.2.
    ! Newton step is dx = -1.0, scaled by xstep to dx = -0.2.
    ! Line search accepts x = -0.2 since f(-0.2) = 1.0 has norm 1.0 <= f(0) norm 1.0.
    ! Residual change dfx = f(-0.2) - f(0) = 0.
    ! The Sherman-Morrison denominator denom = dx^T J^{-1} dfx is identically 0.
    ! Expected outcome: NF_STATUS_INVALID_STATE, x = -0.2, fx = 1.0, finite J^{-1},
    ! and clean IEEE exception flags (no divide-by-zero, no invalid).
    call ieee_set_flag (ieee_divide_by_zero, .false.)
    call ieee_set_flag (ieee_invalid, .false.)

    jac_inv1(1,1) = 1.0_PREC
    x1 = [0.0_PREC]
    call root_broyden (const_residual_fcn, x1, ndiff=.true., &
        xstep=0.2_PREC, maxiter=3, jac_inv_init=jac_inv1, res=res)

    call ieee_get_flag (ieee_divide_by_zero, flag_div_zero)
    call ieee_get_flag (ieee_invalid, flag_invalid)
    call const_residual_fcn (res%x, fx1_eval)

    lres = (.not. res%success) .and. (res%status == NF_STATUS_INVALID_STATE) &
        .and. (abs (res%x(1) - (-0.2_PREC)) < 1.0e-12_PREC) &
        .and. (res%fx(1) == fx1_eval(1)) .and. (fx1_eval(1) == 1.0_PREC) &
        .and. all (ieee_is_finite (res%jac_inv)) &
        .and. (.not. flag_div_zero) .and. (.not. flag_invalid)
    if (.not. lres) error stop "Constant residual; singular Broyden denominator"

    ! Case 2: Rank-deficient / orthogonal secant breakdown in 2D
    ! f(x) = [1.0, 1.0 - x(1)]^T. Initial x = [1.0, 0.0]^T, so f(x) = [1.0, 0.0]^T.
    ! Initial inverse Jacobian is identity I_2.
    ! Full search direction is dx = - I_2 f(x) = [-1.0, 0.0]^T.
    ! Along x + step * dx = [1 - step, 0]^T, f = [1.0, step]^T with norm sqrt(1 + step^2) > 1.
    ! Backtracking line search selects least bad step_size = 0.25 (norm ~ 1.03077).
    ! Accepted step dx = [-0.25, 0.0]^T, advancing x to [0.75, 0.0]^T.
    ! Residual becomes f(x) = [1.0, 0.25]^T, so dfx = [0.0, 0.25]^T is non-zero.
    ! However, J^{-1} dfx = [0.0, 0.25]^T is strictly orthogonal to dx = [-0.25, 0.0]^T:
    ! denom = dot_product(dx, J^{-1} dfx) = (-0.25)*0 + 0*(0.25) = 0.
    ! Expected outcome: NF_STATUS_INVALID_STATE, x = [0.75, 0.0]^T, fx = [1.0, 0.25]^T,
    ! uncorrupted finite inverse Jacobian, and clean IEEE exception flags.
    call ieee_set_flag (ieee_divide_by_zero, .false.)
    call ieee_set_flag (ieee_invalid, .false.)

    jac_inv2 = 0.0_PREC
    jac_inv2(1,1) = 1.0_PREC
    jac_inv2(2,2) = 1.0_PREC
    x2 = [1.0_PREC, 0.0_PREC]
    call root_broyden (secant_breakdown_fcn, x2, ndiff=.true., &
        maxiter=3, jac_inv_init=jac_inv2, res=res)

    call ieee_get_flag (ieee_divide_by_zero, flag_div_zero)
    call ieee_get_flag (ieee_invalid, flag_invalid)
    call secant_breakdown_fcn (res%x, fx2_eval)

    lres = (.not. res%success) .and. (res%status == NF_STATUS_INVALID_STATE) &
        .and. (abs (res%x(1) - 0.75_PREC) < 1.0e-12_PREC) &
        .and. (abs (res%x(2) - 0.0_PREC) < 1.0e-12_PREC) &
        .and. (abs (res%fx(1) - fx2_eval(1)) < 1.0e-12_PREC) &
        .and. (abs (res%fx(2) - fx2_eval(2)) < 1.0e-12_PREC) &
        .and. (abs (fx2_eval(1) - 1.0_PREC) < 1.0e-12_PREC) &
        .and. (abs (fx2_eval(2) - 0.25_PREC) < 1.0e-12_PREC) &
        .and. all (ieee_is_finite (res%jac_inv)) &
        .and. (.not. flag_div_zero) .and. (.not. flag_invalid)
    if (.not. lres) error stop "Orthogonal secant breakdown; rank deficient residual change"

end subroutine


subroutine test_initial_root ()
    !*  Test that Broyden does not return arbitrary/stale inverse Jacobians when
    !   exiting at the initial root (nit=0) or on early failures (Issue #25).
    !
    !   Verifies that res%jac_inv remains unallocated unless a valid inverse was
    !   explicitly supplied by caller, counters (nit, nfev) are initialized, and
    !   reused workspace scratch contents are never leaked into the result.

    type (optim_result) :: res
    type (workspace) :: work
    real (PREC), dimension(1) :: x, fx_eval
    real (PREC), dimension(1,1) :: jac_inv_init
    logical :: lres

    ! Case 1: Exact initial root with fresh workspace, no jac_inv_init
    ! f(x) = x - 1 evaluated at x = 1 yields |f(1)| = 0 < tol.
    ! Broyden exits at iteration 0 (fast early convergence).
    ! Inverse Jacobian must remain unallocated rather than advertising scratch memory.
    x = [1.0_PREC]
    call root_broyden (root_fcn, x, ndiff=.true., res=res)
    call root_fcn (res%x, fx_eval)

    lres = res%success .and. (res%status == NF_STATUS_OK) &
        .and. (res%nit == 0) .and. (res%nfev == 1) &
        .and. (abs (res%x(1) - 1.0_PREC) < 1.0e-12_PREC) &
        .and. (res%fx(1) == fx_eval(1)) .and. (abs (fx_eval(1)) < 1.0e-12_PREC) &
        .and. (.not. allocated (res%jac_inv))
    if (.not. lres) error stop "Exact initial root; fresh workspace; no jac_inv"

    ! Case 2: Exact initial root with reused workspace pre-filled with sentinels
    ! Pre-fill workspace storage with 12345.0 to detect any stale-memory leak.
    call assert_alloc (work, nrwrk=100)
    work%rwrk = 12345.0_PREC
    x = [1.0_PREC]
    call root_broyden (root_fcn, x, ndiff=.true., work=work, res=res)

    lres = res%success .and. (res%status == NF_STATUS_OK) &
        .and. (res%nit == 0) .and. (res%nfev == 1) &
        .and. (abs (res%x(1) - 1.0_PREC) < 1.0e-12_PREC) &
        .and. (.not. allocated (res%jac_inv))
    if (.not. lres) error stop "Exact initial root; reused workspace with sentinels"

    ! Case 3: Exact initial root with caller-supplied valid inverse Jacobian
    ! When caller explicitly provides jac_inv_init, it must be preserved and returned
    ! in res%jac_inv even when exiting at the initial root.
    jac_inv_init(1,1) = 2.0_PREC
    x = [1.0_PREC]
    call root_broyden (root_fcn, x, ndiff=.true., jac_inv_init=jac_inv_init, res=res)

    lres = res%success .and. (res%status == NF_STATUS_OK) &
        .and. (res%nit == 0) .and. (res%nfev == 1) &
        .and. (abs (res%x(1) - 1.0_PREC) < 1.0e-12_PREC) &
        .and. (allocated (res%jac_inv)) &
        .and. (abs (res%jac_inv(1,1) - 2.0_PREC) < 1.0e-12_PREC)
    if (.not. lres) error stop "Exact initial root; caller-supplied jac_inv_init"

    ! Case 4: Early failure before inverse construction (invalid initial residual)
    ! f(x) returns NaN at initial point; must exit with NF_STATUS_INVALID_STATE,
    ! nit=0, nfev=1, and unallocated jac_inv.
    x = [0.0_PREC]
    call root_broyden (nan_initial_fcn, x, ndiff=.true., work=work, res=res)

    lres = (.not. res%success) .and. (res%status == NF_STATUS_INVALID_STATE) &
        .and. (res%nit == 0) .and. (res%nfev == 1) &
        .and. (.not. allocated (res%jac_inv))
    if (.not. lres) error stop "Early failure; invalid initial residual"

    ! Case 5: Early failure on invalid argument before workspace allocation
    ! Invalid maxiter must set NF_STATUS_INVALID_ARG with initialized nit=0 and nfev=0.
    x = [1.0_PREC]
    call root_broyden (root_fcn, x, ndiff=.true., maxiter=-1, res=res)

    lres = (.not. res%success) .and. (res%status == NF_STATUS_INVALID_ARG) &
        .and. (res%nit == 0) .and. (res%nfev == 0) &
        .and. (.not. allocated (res%jac_inv))
    if (.not. lres) error stop "Early failure; invalid argument"

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


subroutine const_residual_fcn (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    ! Piecewise-flat residual that remains constant regardless of x
    if (size (x) > 0) fx(1) = 1.0_PREC
end subroutine


subroutine secant_breakdown_fcn (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    ! Residual change along x(1) direction is orthogonal to x(1) step:
    ! f1 is constant 1.0; f2 changes as 1.0 - x(1).
    fx(1) = 1.0_PREC
    fx(2) = 1.0_PREC - x(1)
end subroutine


subroutine root_fcn (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    ! Linear root function f(x) = x - 1
    fx(1) = x(1) - 1.0_PREC
end subroutine


subroutine nan_initial_fcn (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    ! Always returns quiet NaN
    if (size (x) > 0) fx(1) = ieee_value (0.0_PREC, ieee_quiet_nan)
end subroutine


subroutine test_maxfev ()
    !*  Test strict enforcement of maxfev across all stages of Broyden's method:
    !   - Entry budget check (maxfev = 0).
    !   - Early convergence at initial guess (maxfev = 1).
    !   - Numerical differentiation budget check for scalar and multivariate problems.
    !   - Analytical Jacobian differentiation budget check.
    !   - Line search trial exhaustion before and during trials.
    !   - Exact convergence when budget matches required evaluations.

    type (optim_result) :: res
    real (PREC), dimension(:), allocatable :: x
    real (PREC), dimension(1,1) :: jac_inv_init
    logical :: lres

    ! Case 1: Entry budget check with maxfev = 0
    ! Must terminate before the initial evaluation with NF_STATUS_MAX_EVAL.
    allocate (x(1))
    x(1) = 0.0_PREC
    callback_count = 0
    call root_broyden (counted_root_fcn, x, ndiff=.true., maxfev=0, res=res)
    lres = (.not. res%success) .and. (res%status == NF_STATUS_MAX_EVAL) &
        .and. (res%nit == 0) .and. (res%nfev == 0) &
        .and. (callback_count == 0) .and. (res%x(1) == 0.0_PREC) &
        .and. (.not. allocated (res%jac_inv))
    if (.not. lres) error stop "maxfev = 0 budget check at entry"
    deallocate (x)

    ! Case 2: Early convergence at initial guess with maxfev = 1
    ! f(x) = x - 1 evaluated at x = 1 gives f(1) = 0 < tol. Exits at nit = 0.
    allocate (x(1))
    x(1) = 1.0_PREC
    callback_count = 0
    call root_broyden (counted_root_fcn, x, ndiff=.true., maxfev=1, res=res)
    lres = res%success .and. (res%status == NF_STATUS_OK) &
        .and. (res%nit == 0) .and. (res%nfev == 1) &
        .and. (callback_count == 1) &
        .and. (abs (res%x(1) - 1.0_PREC) < 1.0e-12_PREC)
    if (.not. lres) error stop "maxfev = 1 initial root convergence"
    deallocate (x)

    ! Case 3: Numerical differentiation budget check with n = 1 and maxfev = 1 (GH #26)
    ! Initial evaluation takes 1 call. Numerical differentiation requires n = 1 call.
    ! Total required is 1 + 1 = 2 > maxfev. Must terminate before differentiation.
    allocate (x(1))
    x(1) = 0.0_PREC
    callback_count = 0
    call root_broyden (counted_root_fcn, x, ndiff=.true., maxfev=1, res=res)
    lres = (.not. res%success) .and. (res%status == NF_STATUS_MAX_EVAL) &
        .and. (res%nit == 0) .and. (res%nfev == 1) &
        .and. (callback_count == 1) .and. (res%x(1) == 0.0_PREC) &
        .and. (abs (res%fx(1) - (-1.0_PREC)) < 1.0e-12_PREC)
    if (.not. lres) error stop "maxfev = 1 scalar numdiff budget check"
    deallocate (x)

    ! Case 4: Multivariate numerical differentiation budget check (n = 3, maxfev = 2)
    ! Initial evaluation takes 1 call. Numerical differentiation requires n = 3 calls.
    ! Total required is 1 + 3 = 4 > maxfev. Must terminate before differentiation.
    allocate (x(3))
    x(:) = 0.0_PREC
    callback_count = 0
    call root_broyden (counted_multivar_fcn, x, ndiff=.true., maxfev=2, res=res)
    lres = (.not. res%success) .and. (res%status == NF_STATUS_MAX_EVAL) &
        .and. (res%nit == 0) .and. (res%nfev == 1) &
        .and. (callback_count == 1) .and. all (res%x == 0.0_PREC)
    if (.not. lres) error stop "Multivariate numdiff budget check"
    deallocate (x)

    ! Case 5: Analytical Jacobian budget check (n = 2, maxfev = 1)
    ! Initial evaluation takes 1 call. Analytical Jacobian evaluation takes 1 call.
    ! Total required is 1 + 1 = 2 > maxfev. Must terminate before Jacobian evaluation.
    allocate (x(2))
    x(:) = 0.0_PREC
    callback_count = 0
    call root_broyden (counted_multivar_fcn, counted_analytical_jac, x, &
        maxfev=1, res=res)
    lres = (.not. res%success) .and. (res%status == NF_STATUS_MAX_EVAL) &
        .and. (res%nit == 0) .and. (res%nfev == 1) &
        .and. (callback_count == 1) .and. all (res%x == 0.0_PREC)
    if (.not. lres) error stop "Analytical Jacobian budget check"
    deallocate (x)

    ! Case 6: Line search budget exhaustion before any trial step can be evaluated
    ! jac_inv_init bypasses differentiation. Initial residual takes 1 call (nfev = 1).
    ! maxfev = 1 allows 0 line search evaluations; dumb_line_search exits immediately.
    allocate (x(1))
    x(1) = 0.0_PREC
    jac_inv_init(1,1) = 10.0_PREC
    callback_count = 0
    call root_broyden (counted_root_fcn, x, ndiff=.true., &
        jac_inv_init=jac_inv_init, maxfev=1, res=res)
    lres = (.not. res%success) .and. (res%status == NF_STATUS_MAX_EVAL) &
        .and. (res%nit == 0) .and. (res%nfev == 1) &
        .and. (callback_count == 1) .and. (res%x(1) == 0.0_PREC)
    if (.not. lres) error stop "Line search exhaustion before trial evaluation"
    deallocate (x)

    ! Case 7: Line search budget exhaustion after evaluating first trial step
    ! Initial evaluation takes 1 call (f(0) = -1). Step dx = 10 produces f(10) = 9.
    ! Trial 1 evaluation takes 1 call (total nfev = 2 = maxfev). Backtracking halts.
    allocate (x(1))
    x(1) = 0.0_PREC
    jac_inv_init(1,1) = 10.0_PREC
    callback_count = 0
    call root_broyden (counted_root_fcn, x, ndiff=.true., &
        jac_inv_init=jac_inv_init, maxfev=2, res=res)
    lres = (.not. res%success) .and. (res%status == NF_STATUS_MAX_EVAL) &
        .and. (res%nit == 1) .and. (res%nfev == 2) &
        .and. (callback_count == 2)
    if (.not. lres) error stop "Line search exhaustion after trial evaluation"
    deallocate (x)

    ! Case 8: Exact convergence when budget matches required evaluations
    ! Initial evaluation takes 1 call (f(0) = -1). Step dx = 1.0 reaches f(1) = 0.
    ! Trial 1 evaluation takes 1 call (total nfev = 2 = maxfev). Converges at nit = 1.
    allocate (x(1))
    x(1) = 0.0_PREC
    jac_inv_init(1,1) = 1.0_PREC
    callback_count = 0
    call root_broyden (counted_root_fcn, x, ndiff=.true., &
        jac_inv_init=jac_inv_init, maxfev=2, res=res)
    lres = res%success .and. (res%status == NF_STATUS_OK) &
        .and. (res%nit == 1) .and. (res%nfev == 2) &
        .and. (callback_count == 2) &
        .and. (abs (res%x(1) - 1.0_PREC) < 1.0e-12_PREC)
    if (.not. lres) error stop "Exact convergence with tight budget"
    deallocate (x)

end subroutine


subroutine counted_root_fcn (x, fx)
    !*  Scalar linear test function f(x) = x - 1 with independent callback counter
    !   used to verify maxfev enforcement.

    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Evaluation point
    real (PREC), intent(out), dimension(:), contiguous :: fx
        !*  Residual vector f(x) = x - 1

    callback_count = callback_count + 1
    if (size (x) > 0) fx(1) = x(1) - 1.0_PREC
end subroutine


subroutine counted_multivar_fcn (x, fx)
    !*  Multivariate test function f(x) = x - 1 with independent callback counter
    !   used to verify differentiation evaluation budgets.

    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Evaluation point
    real (PREC), intent(out), dimension(:), contiguous :: fx
        !*  Residual vector f(x) = x - 1

    callback_count = callback_count + 1
    fx(:) = x(:) - 1.0_PREC
end subroutine


subroutine counted_analytical_jac (x, jac)
    !*  Analytical Jacobian for multivariate f(x) = x - 1 (identity matrix) with
    !   independent callback counter used to verify analytical Jacobian budgets.

    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Evaluation point
    real (PREC), intent(out), dimension(:,:), contiguous :: jac
        !*  Jacobian matrix (identity matrix)

    integer :: i

    callback_count = callback_count + 1
    jac(:,:) = 0.0_PREC
    do i = 1, min (size (jac, 1), size (jac, 2))
        jac(i,i) = 1.0_PREC
    end do
    if (size (x) < 0) jac(1,1) = 0.0_PREC
end subroutine


end program
