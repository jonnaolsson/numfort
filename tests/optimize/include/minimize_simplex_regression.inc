subroutine test_constant ()
    !* - Stop constant objectives with the model's default quadratic fit.
    !  - Preserve the successful result when the fit is explicitly disabled.
    type(optim_result) :: res
    real(PREC) :: x(2)
    integer :: i

    ! Positive, zero, and negative constants have every point as a minimizer.
    do i = 1, 3
        constant_value = real(2 - i, PREC)
        x = 0.0_PREC
        call reset_callback (10000)
        call minimize_simplex (fcn_constant, x, rstep=0.0_PREC, xstep=1.0_PREC, &
            maxfun=10000, res=res)
        if (res%status /= NF_STATUS_OK + NF_STATUS_APPROX) &
            error stop "Constant objective must report unavailable quadratic fit"
        if (index(res%msg, "quadratic fit unavailable") == 0) &
            error stop "Missing quadratic-fit diagnostic"
        if (calls >= 200) error stop "Constant objective used too many callbacks"
        call check_result (x, res, constant_value)
    end do

    ! A small budget expires during expansion before flatness is established.
    x = 0.0_PREC
    call reset_callback (22)
    call minimize_simplex (fcn_constant, x, rstep=0.0_PREC, xstep=1.0_PREC, &
        maxfun=22, res=res)
    if (res%status /= NF_STATUS_MAX_EVAL) error stop "Constant fit budget status"
    if (calls /= 22) error stop "Constant fit callback count"
    call check_result (x, res, constant_value)

    ! Disabling the optional fit still gives an ordinary successful result.
    x = 0.0_PREC
    call reset_callback (10000)
    call minimize_simplex (fcn_constant, x, rstep=0.0_PREC, xstep=1.0_PREC, &
        maxfun=10000, quad=.false., res=res)
    if (res%status /= NF_STATUS_OK) error stop "Constant simplex search failed"
    call check_result (x, res, constant_value)
end subroutine


subroutine test_partial_flat ()
    !* - Stop a quadratic that is independent of the second coordinate.
    type(optim_result) :: res
    real(PREC) :: x(2)

    ! f(x) = x(1)**2 has minimum zero for any x(2) when x(1) = 0.
    x = 0.0_PREC
    call reset_callback (10000)
    call minimize_simplex (fcn_partial, x, rstep=0.0_PREC, xstep=1.0_PREC, &
        maxfun=10000, res=res)
    if (res%status /= NF_STATUS_OK + NF_STATUS_APPROX) &
        error stop "Partially flat objective must report unavailable fit"
    if (calls >= 300) error stop "Partially flat objective used too many callbacks"
    if (abs(x(1)) > 1.0e-3_PREC) error stop "Incorrect partially flat minimizer"
    call check_result (x, res, x(1)**2)
end subroutine


subroutine test_search_budget ()
    !* - Enforce small budgets during initialization and simplex iterations.
    !  - Return the best evaluated point with matching objective and call count.
    type(optim_result) :: res
    real(PREC) :: x(2)
    integer :: budget

    ! Budgets 1 through 8 interrupt initialization, reflection, and contraction.
    do budget = 1, 8
        x = [2.0_PREC, -1.0_PREC]
        objective_scale = 1.0_PREC
        call reset_callback (budget)
        call minimize_simplex (fcn_quadratic, x, rstep=0.0_PREC, xstep=1.0_PREC, &
            maxfun=budget, res=res)
        if (res%status /= NF_STATUS_MAX_EVAL) error stop "Search budget status"
        if (calls /= budget) error stop "Search did not use its callback budget"
        if (res%fx(1) /= best_value) error stop "Search lost best evaluated value"
        call check_result (x, res, sum(x**2))
    end do
end subroutine


subroutine test_fit_budget ()
    !* - Enforce budgets at the centroid, expansion, and fit interpolation points.
    !  - Allow success when the last permitted evaluation completes the fit.
    type(optim_result) :: res
    real(PREC) :: x(2)
    integer :: search_calls, fit_calls, budgets(6), i, budget

    ! A shallow positive definite quadratic requires additional expansion probes.
    objective_scale = 1.0e-20_PREC
    x = 0.0_PREC
    call reset_callback (10000)
    call minimize_simplex (fcn_quadratic, x, rstep=0.0_PREC, xstep=1.0_PREC, &
        maxfun=10000, quad=.false., res=res)
    if (res%status /= NF_STATUS_OK) error stop "Shallow simplex search failed"
    search_calls = calls

    ! Establish the completed-fit count without assuming a compiler-specific path.
    x = 0.0_PREC
    call reset_callback (10000)
    call minimize_simplex (fcn_quadratic, x, rstep=0.0_PREC, xstep=1.0_PREC, &
        maxfun=10000, res=res)
    if (res%status /= NF_STATUS_OK) error stop "Shallow quadratic fit failed"
    fit_calls = calls
    if (fit_calls <= search_calls + 3) error stop "Quadratic expansion was skipped"
    call check_result (x, res, objective_scale * sum(x**2))

    ! Stop before/after the centroid, in expansion, and at the last fit samples.
    budgets = [search_calls - 1, search_calls, search_calls + 1, &
        fit_calls - 2, fit_calls - 1, fit_calls]
    do i = 1, size(budgets)
        budget = budgets(i)
        x = 0.0_PREC
        call reset_callback (budget)
        call minimize_simplex (fcn_quadratic, x, rstep=0.0_PREC, xstep=1.0_PREC, &
            maxfun=budget, res=res)
        if (budget < fit_calls) then
            if (res%status /= NF_STATUS_MAX_EVAL) error stop "Fit budget status"
            if (res%fx(1) /= best_value) error stop "Fit lost best evaluated value"
        else
            if (res%status /= NF_STATUS_OK) error stop "Exact fit budget failed"
        end if
        if (calls /= budget) error stop "Fit did not use its callback budget"
        call check_result (x, res, objective_scale * sum(x**2))
    end do
end subroutine


subroutine test_quadratic_fit ()
    !* - Preserve a complete quadratic fit for an identifiable objective.
    type(optim_result) :: res
    real(PREC) :: x(2)
    integer :: search_calls

    ! The sum of squares has unique minimum zero at the origin.
    objective_scale = 1.0_PREC
    x = [2.0_PREC, -1.0_PREC]
    call reset_callback (10000)
    call minimize_simplex (fcn_quadratic, x, tol=1.0e-7_PREC, &
        rstep=0.0_PREC, xstep=1.0_PREC, maxfun=10000, quad=.false., res=res)
    if (res%status /= NF_STATUS_OK) error stop "Quadratic simplex search failed"
    search_calls = calls

    ! Default quad must perform extra evaluations and finish without approximation.
    x = [2.0_PREC, -1.0_PREC]
    call reset_callback (10000)
    call minimize_simplex (fcn_quadratic, x, tol=1.0e-7_PREC, &
        rstep=0.0_PREC, xstep=1.0_PREC, maxfun=10000, res=res)
    if (res%status /= NF_STATUS_OK) error stop "Identifiable quadratic fit failed"
    if (calls <= search_calls) error stop "Optional quadratic fit was skipped"
    if (any(abs(x) > 1.0e-3_PREC)) error stop "Incorrect quadratic minimizer"
    call check_result (x, res, sum(x**2))
end subroutine


subroutine test_degenerate ()
    !* - Stop when expansion cannot change a representable coordinate.
    !  - Stop before expanded coordinates overflow.
    type(optim_result) :: res
    real(PREC) :: x(2), large_step

    ! The unit step rounds away at 2**(digits+2), collapsing the simplex.
    x = scale(1.0_PREC, digits(1.0_PREC) + 2)
    constant_value = 1.0_PREC
    call reset_callback (10000)
    call minimize_simplex (fcn_constant, x, rstep=0.0_PREC, xstep=1.0_PREC, &
        maxfun=10000, res=res)
    if (res%status /= NF_STATUS_OK + NF_STATUS_APPROX) &
        error stop "Collapsed simplex must report unavailable fit"
    call check_result (x, res, constant_value)

    ! Large finite steps reach the coordinate limit before the flat-probe limit.
    large_step = huge(1.0_PREC) / 64.0_PREC
    x = 0.0_PREC
    call ieee_set_flag (ieee_overflow, .false.)
    call reset_callback (10000)
    call minimize_simplex (fcn_constant, x, rstep=0.0_PREC, xstep=large_step, &
        maxfun=10000, res=res)
    if (res%status /= NF_STATUS_OK + NF_STATUS_APPROX) &
        error stop "Unbounded expansion must report unavailable fit"
    call check_result (x, res, constant_value)
    block
        logical :: overflowed
        call ieee_get_flag (ieee_overflow, overflowed)
        if (overflowed) error stop "Quadratic expansion overflowed"
    end block
end subroutine


subroutine test_nonfinite ()
    !* - Preserve a finite result when a search or fit probe leaves the domain.
    type(optim_result) :: res
    real(PREC) :: x(2)

    ! The initial point is finite, but the next simplex vertex leaves the domain.
    x = [1.5_PREC, 0.0_PREC]
    call reset_callback (10000)
    call minimize_simplex (fcn_bounded, x, rstep=0.0_PREC, xstep=1.0_PREC, &
        maxfun=10000, res=res)
    if (res%status /= NF_STATUS_DOMAIN_ERROR) error stop "Nonfinite search status"
    if (calls /= 2) error stop "Nonfinite search callback count"
    call check_result (x, res, 1.0_PREC)

    ! The objective is constant on a bounded domain and NaN outside it.
    x = 0.0_PREC
    call reset_callback (10000)
    call minimize_simplex (fcn_bounded, x, rstep=0.0_PREC, xstep=1.0_PREC, &
        maxfun=10000, res=res)
    if (res%status /= NF_STATUS_OK + NF_STATUS_APPROX) &
        error stop "Nonfinite fit probe must report unavailable fit"
    if (any(abs(x) > 2.0_PREC)) error stop "Returned point outside finite domain"
    call check_result (x, res, 1.0_PREC)
end subroutine


subroutine reset_callback (limit)
    !* Reset independent callback accounting for one solver invocation.
    integer, intent(in) :: limit
        !* Hard callback limit, independent of the solver's counter.

    calls = 0
    callback_limit = limit
    best_value = huge(1.0_PREC)
end subroutine


subroutine record_callback (x, fx)
    !* Count actual calls and reject nonfinite coordinates before a regression hangs.
    real(PREC), intent(in) :: x(:)
        !* Evaluated coordinates.
    real(PREC), intent(in) :: fx
        !* Observed objective value.

    calls = calls + 1
    if (calls > callback_limit) error stop "Simplex exceeded callback guard"
    if (.not. all(ieee_is_finite(x))) error stop "Nonfinite callback coordinates"
    if (ieee_is_finite(fx)) best_value = min(best_value, fx)
end subroutine


subroutine check_result (x, res, expected)
    !* Check finite, consistent results and independently counted callbacks.
    real(PREC), intent(in) :: x(:)
        !* Returned minimizer.
    type(optim_result), intent(in) :: res
        !* Solver result.
    real(PREC), intent(in) :: expected
        !* Objective computed independently at the returned point.

    if (res%nfev /= calls) error stop "Incorrect simplex evaluation count"
    if (calls > callback_limit) error stop "Callback budget exceeded"
    if (.not. all(ieee_is_finite(x))) error stop "Nonfinite simplex result"
    if (.not. all(res%x == x)) error stop "Inconsistent simplex coordinates"
    if (.not. ieee_is_finite(res%fx(1))) error stop "Nonfinite simplex objective"
    if (abs(res%fx(1) - expected) > &
        8.0_PREC * epsilon(expected) * max(tiny(expected), abs(expected))) &
        error stop "Inconsistent simplex objective"
end subroutine


subroutine fcn_constant (x, fx)
    !* Evaluate a constant objective.
    real(PREC), intent(in), contiguous :: x(:)
        !* Candidate coordinates.
    real(PREC), intent(out) :: fx
        !* Objective value.

    fx = constant_value
    call record_callback (x, fx)
end subroutine


subroutine fcn_partial (x, fx)
    !* Evaluate a quadratic with a flat second coordinate.
    real(PREC), intent(in), contiguous :: x(:)
        !* Candidate coordinates.
    real(PREC), intent(out) :: fx
        !* Objective value.

    fx = x(1)**2
    call record_callback (x, fx)
end subroutine


subroutine fcn_quadratic (x, fx)
    !* Evaluate a scaled positive definite quadratic.
    real(PREC), intent(in), contiguous :: x(:)
        !* Candidate coordinates.
    real(PREC), intent(out) :: fx
        !* Objective value.

    fx = objective_scale * sum(x**2)
    call record_callback (x, fx)
end subroutine


subroutine fcn_bounded (x, fx)
    !* Evaluate a constant objective with a bounded finite domain.
    real(PREC), intent(in), contiguous :: x(:)
        !* Candidate coordinates.
    real(PREC), intent(out) :: fx
        !* Objective value or NaN outside the domain.

    fx = 1.0_PREC
    if (any(abs(x) > 2.0_PREC)) fx = ieee_value(fx, ieee_quiet_nan)
    call record_callback (x, fx)
end subroutine
