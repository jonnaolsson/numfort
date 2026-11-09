program test_optimize_minimize_dfls
    !* Regression tests for the derivative-free least-squares minimizer.

    use, intrinsic :: iso_fortran_env, only: real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan

    use numfort_common
    use numfort_optimize, optim_result => optim_result_real64

    implicit none

    integer, parameter :: PREC = real64
    integer :: calls, invalid_call

    call test_evaluation_budget ()
    call test_convergence ()
    call test_overdetermined ()
    call test_numerical_failure ()
    call test_invalid_input ()
    call test_real32 ()

contains

subroutine test_evaluation_budget ()
    !* - Enforce tiny budgets using an independent callback counter.
    !  - Return matching coordinates and residuals before and after initialization.

    real (PREC) :: x(2)
    type (optim_result) :: res
    type (workspace_real64) :: work
    integer :: budget

    ! Cases 1-2: Neither of the first two sampled points improves on [3, 4].
    do budget = 1, 2
        x = [3.0_PREC, 4.0_PREC]
        calls = 0
        call minimize_dfls (identity, x, 2, 1.0_PREC, 1.0e-8_PREC, &
            maxfev=budget, work=work, res=res)
        if (res%success .or. res%status /= NF_STATUS_MAX_EVAL) &
            error stop "Tiny budget must report evaluation exhaustion"
        if (calls /= budget .or. res%nfev /= calls) &
            error stop "Tiny budget callback count"
        if (any(x /= [3.0_PREC, 4.0_PREC])) error stop "Tiny budget best point"
        if (any(res%x /= x) .or. any(res%fx /= x)) &
            error stop "Tiny budget residual correspondence"
    end do

    ! Case 3: The second sample improves the first coordinate from -3 to -2.
    x = [-3.0_PREC, 4.0_PREC]
    calls = 0
    call minimize_dfls (identity, x, 2, 1.0_PREC, 1.0e-8_PREC, &
        maxfev=2, work=work, res=res)
    if (res%success .or. res%status /= NF_STATUS_MAX_EVAL) &
        error stop "Partial initialization must report evaluation exhaustion"
    if (calls /= 2 .or. res%nfev /= calls) error stop "Partial initialization count"
    if (any(x /= [-2.0_PREC, 4.0_PREC])) &
        error stop "Partial initialization must retain the best point"
    if (any(res%x /= x) .or. any(res%fx /= x)) &
        error stop "Partial initialization residual correspondence"

    ! Case 4: Exhaust the budget after constructing the five-point model.
    x = [3.0_PREC, 4.0_PREC]
    calls = 0
    call minimize_dfls (identity, x, 2, 1.0_PREC, 1.0e-8_PREC, &
        maxfev=6, work=work, res=res)
    if (res%success .or. res%status /= NF_STATUS_MAX_EVAL) &
        error stop "Iteration budget must report evaluation exhaustion"
    if (calls /= 6 .or. res%nfev /= calls) error stop "Iteration budget count"
    if (any(res%x /= x) .or. any(res%fx /= x)) &
        error stop "Iteration budget residual correspondence"

    ! Case 5: Omitting the result must also respect the evaluation budget.
    x = [3.0_PREC, 4.0_PREC]
    calls = 0
    call minimize_dfls (identity, x, 2, 1.0_PREC, 1.0e-8_PREC, maxfev=2)
    if (calls /= 2) error stop "Budget without result"
end subroutine


subroutine test_convergence ()
    !* - Distinguish convergence at the budget from evaluation exhaustion.
    !  - Check a zero minimum and the existing nonlinear residual problem.

    real (PREC) :: x(2), xquad(5)
    type (optim_result) :: res
    integer :: budget

    ! Cases 1-2: Zero residual is reached on the first or second allowed call.
    do budget = 1, 2
        x = [real(1 - budget, PREC), 0.0_PREC]
        calls = 0
        call minimize_dfls (identity, x, 2, 1.0_PREC, 1.0e-8_PREC, &
            maxfev=budget, res=res)
        if (.not. res%success .or. res%status /= NF_STATUS_OK) &
            error stop "Convergence at exact budget"
        if (calls /= budget .or. res%nfev /= calls) &
            error stop "Convergence at exact budget count"
        if (any(x /= 0.0_PREC) .or. any(res%fx /= 0.0_PREC)) &
            error stop "Exact zero minimum"
        if (any(res%x /= x)) error stop "Exact zero result coordinates"
    end do

    ! Case 3: The sum of squared identity residuals has its minimum at x = 0.
    x = [3.0_PREC, 4.0_PREC]
    calls = 0
    call minimize_dfls (identity, x, 2, 1.0_PREC, 1.0e-8_PREC, &
        maxfev=100, res=res)
    if (.not. res%success .or. res%status /= NF_STATUS_OK) &
        error stop "Identity convergence"
    if (calls > 100 .or. res%nfev /= calls) error stop "Identity convergence count"
    if (any(abs(x) > 1.0e-6_PREC)) error stop "Identity minimum"
    if (any(res%x /= x) .or. any(res%fx /= x)) &
        error stop "Identity residual correspondence"

    ! Case 4: Squared residuals give sum(x**4), whose unique minimum is x = 0.
    xquad = [-2.0_PREC, -1.0_PREC, 0.5_PREC, 1.0_PREC, 2.0_PREC]
    calls = 0
    call minimize_dfls (quadratic, xquad, 5, 1.0_PREC, 1.0e-8_PREC, res=res)
    if (.not. res%success .or. res%status /= NF_STATUS_OK) &
        error stop "Quadratic residual convergence"
    if (res%nfev /= calls .or. calls > 600) error stop "Quadratic residual count"
    if (any(abs(res%fx) > 1.0e-3_PREC)) error stop "Quadratic residual minimum"
    if (any(res%x /= xquad) .or. any(res%fx /= xquad**2)) &
        error stop "Quadratic residual correspondence"
end subroutine


subroutine test_overdetermined ()
    !* - Solve an overdetermined linear problem with a nonzero minimum.
    !  - Check trust-region convergence with ample and exact evaluation budgets.

    real (PREC) :: x(2), expected(3)
    type (optim_result) :: res
    integer :: budget, run

    ! The normal equations give x = [4/3, 7/3], residuals [1/3, 1/3, -1/3],
    ! and a minimum sum of squares of 1/3.
    budget = 200
    do run = 1, 2
        ! Case 1 uses ample budget; case 2 repeats with exactly the observed count.
        x = [3.0_PREC, 4.0_PREC]
        calls = 0
        call minimize_dfls (overdetermined, x, 3, 1.0_PREC, 1.0e-8_PREC, &
            maxfev=budget, res=res)
        if (.not. res%success .or. res%status /= NF_STATUS_OK) &
            error stop "Overdetermined convergence"
        if (calls > budget .or. res%nfev /= calls) &
            error stop "Overdetermined callback count"
        if (run == 2 .and. calls /= budget) error stop "Exact convergence budget"
        if (any(abs(x - [4.0_PREC, 7.0_PREC] / 3.0_PREC) > 1.0e-6_PREC)) &
            error stop "Overdetermined minimizer"
        if (abs(sum(res%fx**2) - 1.0_PREC / 3.0_PREC) > 1.0e-12_PREC) &
            error stop "Overdetermined nonzero minimum"
        expected = [x(1) - 1.0_PREC, x(2) - 2.0_PREC, sum(x) - 4.0_PREC]
        if (any(res%x /= x) .or. any(abs(res%fx - expected) > 1.0e-14_PREC)) &
            error stop "Overdetermined residual correspondence"
        budget = calls
    end do
end subroutine


subroutine test_numerical_failure ()
    !* - Reject nonfinite residuals without certifying convergence.
    !  - Preserve the best valid point when a later evaluation fails.

    real (PREC) :: x(2)
    type (optim_result) :: res

    ! Cases 1-2: A nonfinite first or second evaluation terminates immediately.
    do invalid_call = 1, 2
        x = [3.0_PREC, 4.0_PREC]
        calls = 0
        call minimize_dfls (invalid_residual, x, 2, 1.0_PREC, 1.0e-8_PREC, &
            maxfev=100, res=res)
        if (res%success .or. res%status /= NF_STATUS_INVALID_STATE) &
            error stop "Nonfinite residual must report numerical failure"
        if (calls /= invalid_call .or. res%nfev /= calls) &
            error stop "Numerical failure count"
        if (any(x /= [3.0_PREC, 4.0_PREC]) .or. any(res%x /= x)) &
            error stop "Numerical failure coordinates"
        if (invalid_call == 1) then
            if (all(ieee_is_finite(res%fx))) error stop "Initial invalid residual"
        else
            if (any(res%fx /= x)) error stop "Numerical failure best residual"
        end if
    end do
end subroutine


subroutine test_invalid_input ()
    !* - Reject unsupported dimensions and interpolation counts before callbacks.

    real (PREC) :: x(2)
    type (optim_result) :: res

    ! Case 1: NEWUOA requires at least two variables.
    x = 1.0_PREC
    calls = 0
    call minimize_dfls (identity, x(:1), 1, 1.0_PREC, 1.0e-8_PREC, res=res)
    if (res%success .or. res%status /= NF_STATUS_INVALID_ARG) &
        error stop "Invalid dimension status"
    if (calls /= 0 .or. res%nfev /= 0) error stop "Invalid dimension count"

    ! Case 2: The least-squares variant permits at most 2*n+1 = 5 points.
    call minimize_dfls (identity, x, 2, 1.0_PREC, 1.0e-8_PREC, nint=6, res=res)
    if (res%success .or. res%status /= NF_STATUS_INVALID_ARG) &
        error stop "Invalid interpolation count status"
    if (calls /= 0 .or. res%nfev /= 0) error stop "Invalid interpolation count"
end subroutine


subroutine test_real32 ()
    !* - Check the single-precision wrapper's budget and convergence statuses.

    real (real32) :: x(2)
    type (optim_result_real32) :: res
    integer :: budget

    ! Cases 1-2: Tiny budgets must return the initial point and its residual.
    do budget = 1, 2
        x = [3.0_real32, 4.0_real32]
        calls = 0
        call minimize_dfls (identity_real32, x, 2, 1.0_real32, 1.0e-4_real32, &
            maxfev=budget, res=res)
        if (res%success .or. res%status /= NF_STATUS_MAX_EVAL) &
            error stop "real32 tiny budget status"
        if (calls /= budget .or. res%nfev /= calls) error stop "real32 budget count"
        if (any(x /= [3.0_real32, 4.0_real32])) error stop "real32 best point"
        if (any(res%x /= x) .or. any(res%fx /= x)) &
            error stop "real32 residual correspondence"
    end do

    ! Case 3: With enough evaluations, sum(x**2) reaches its zero minimum.
    x = [3.0_real32, 4.0_real32]
    calls = 0
    call minimize_dfls (identity_real32, x, 2, 1.0_real32, 1.0e-4_real32, &
        maxfev=100, res=res)
    if (.not. res%success .or. res%status /= NF_STATUS_OK) &
        error stop "real32 convergence"
    if (calls > 100 .or. res%nfev /= calls) error stop "real32 convergence count"
    if (any(abs(x) > 1.0e-5_real32)) error stop "real32 minimum"
    if (any(res%x /= x) .or. any(res%fx /= x)) &
        error stop "real32 convergence residual"
end subroutine


subroutine identity (x, fx)
    !* Evaluate identity residuals and count actual callback invocations.
    real (PREC), intent(in), contiguous :: x(:)
        !* Evaluation point.
    real (PREC), intent(out), contiguous :: fx(:)
        !* Residual vector.

    calls = calls + 1
    fx = x
end subroutine


subroutine quadratic (x, fx)
    !* Evaluate squared residuals and count actual callback invocations.
    real (PREC), intent(in), contiguous :: x(:)
        !* Evaluation point.
    real (PREC), intent(out), contiguous :: fx(:)
        !* Residual vector.

    calls = calls + 1
    fx = x**2
end subroutine


subroutine overdetermined (x, fx)
    !* Evaluate three linear residuals in two variables.
    real (PREC), intent(in), contiguous :: x(:)
        !* Evaluation point.
    real (PREC), intent(out), contiguous :: fx(:)
        !* Residual vector.

    calls = calls + 1
    fx = [x(1) - 1.0_PREC, x(2) - 2.0_PREC, sum(x) - 4.0_PREC]
end subroutine


subroutine invalid_residual (x, fx)
    !* Inject a nonfinite residual on the selected evaluation.
    real (PREC), intent(in), contiguous :: x(:)
        !* Evaluation point.
    real (PREC), intent(out), contiguous :: fx(:)
        !* Residual vector.

    calls = calls + 1
    fx = x
    if (calls == invalid_call) fx(1) = ieee_value (0.0_PREC, ieee_quiet_nan)
end subroutine


subroutine identity_real32 (x, fx)
    !* Evaluate single-precision identity residuals and count invocations.
    real (real32), intent(in), contiguous :: x(:)
        !* Evaluation point.
    real (real32), intent(out), contiguous :: fx(:)
        !* Residual vector.

    calls = calls + 1
    fx = x
end subroutine

end program
