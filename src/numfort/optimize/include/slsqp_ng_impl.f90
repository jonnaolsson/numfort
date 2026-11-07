



recursive subroutine slsqp_jac (fobj, x, work, fcon, m, meq, lbounds, ubounds, &
        tol, xtol, maxiter, maxfev, exact_lsearch, maximize, autoscale, iprint, res)
    procedure (fvs_fcn_jac_opt) :: fobj
    real (PREC), intent(inout), dimension(:), contiguous :: x
    type (workspace), intent(inout) :: work
    procedure (fvv_fcn_jac_opt), optional :: fcon
    integer, intent(in), optional :: m
    integer, intent(in), optional :: meq
    real (PREC), intent(in), dimension(:), contiguous, optional :: lbounds, ubounds
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    logical, intent(in), optional :: exact_lsearch
    logical, intent(in), optional :: maximize
        !*  If present and .TRUE., perform maximization instead of minimization.
    logical, intent(in), optional :: autoscale
        !*  If present and .TRUE., automatically rescale problem so that
        !   the objective and the elements of its gradient are all one
        !   in absolute value.
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    type (optim_result), intent(out), optional :: res

    type (fwrapper_vs) :: fobj_wrapper
    type (fwrapper_vv) :: fcon_wrapper
    type (optim_result) :: lres
    integer :: lm, lmeq, lmaxiter, lmaxfev
    real (PREC) :: ltol, lxtol
    logical :: lexact_lsearch, lautoscale, lmaximize
    integer (NF_ENUM_KIND) :: liprint

    lres%status = NF_STATUS_OK
    lres%msg = ''

    call check_input (present(fcon), m, meq, x, lbounds, ubounds, tol, xtol, &
        maxiter, maxfev, exact_lsearch, iprint, status=lres%status, msg=lres%msg)
    if (lres%status /= NF_STATUS_OK) goto 100

    call wrap_procedure (fobj_wrapper, fcn_jac_opt=fobj)
    call wrap_procedure (fcon_wrapper, fcn_jac_opt=fcon)

    ! === Default values ===

    call set_optional_arg (meq, 0, lmeq)
    call set_optional_arg (m, lmeq, lm)
    call set_optional_arg (tol, 1.0e-4_PREC, ltol)
    call set_optional_arg (xtol, 1.0e-6_PREC, lxtol)
    call set_optional_arg (exact_lsearch, .false., lexact_lsearch)
    call set_optional_arg (maximize, .false., lmaximize)
    call set_optional_arg (autoscale, .false., lautoscale)
    call set_optional_arg (maxiter, 100, lmaxiter)
    call set_optional_arg (maxfev, huge(1), lmaxfev)
    call set_optional_arg (iprint, NF_PRINT_NONE, liprint)

    ! === Call implementation routine ===

    call slsqp_impl (fobj_wrapper, x, work, fcon_wrapper, lm, lmeq, &
        lbounds, ubounds, ltol, lxtol, lmaxiter, lmaxfev, lexact_lsearch, &
        lmaximize, lautoscale, liprint, lres)

100 continue

    if (present(res)) res = lres

end subroutine



recursive subroutine slsqp (fobj, x, ndiff, work, fcon, m, meq, lbounds, ubounds, &
        tol, xtol, maxiter, maxfev, exact_lsearch, maximize, autoscale, iprint, xstep, res)
    procedure (fvs_fcn) :: fobj
    real (PREC), intent(inout), dimension(:), contiguous :: x
    logical, intent(in) :: ndiff
    type (workspace), intent(inout) :: work
    procedure (fvv_fcn), optional :: fcon
    integer, intent(in), optional :: m
    integer, intent(in), optional :: meq
    real (PREC), intent(in), dimension(:), contiguous, optional :: lbounds, ubounds
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    logical, intent(in), optional :: exact_lsearch
    logical, intent(in), optional :: maximize
        !*  If present and .TRUE., perform maximization instead of minimization.
    logical, intent(in), optional :: autoscale
        !*  If present and .TRUE., automatically rescale problem so that
        !   the objective and the elements of its gradient are all one
        !   in absolute value.
        integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: xstep
    type (optim_result), intent(out), optional :: res

    type (fwrapper_vs) :: fobj_wrapper
    type (fwrapper_vv) :: fcon_wrapper
    type (optim_result) :: lres
    integer :: lm, lmeq, lmaxiter, lmaxfev
    real (PREC) :: ltol, lxtol
    logical :: lexact_lsearch, lautoscale, lmaximize
    integer (NF_ENUM_KIND) :: liprint

    lres%status = NF_STATUS_OK
    lres%msg = ''

    call check_input (present(fcon), m, meq, x, lbounds, ubounds, tol, xtol, &
        maxiter, maxfev, exact_lsearch, iprint, xstep, lres%status, lres%msg)
    if (lres%status /= NF_STATUS_OK) goto 100

    call wrap_procedure (fobj_wrapper, fcn=fobj, eps=xstep)
    call wrap_procedure (fcon_wrapper, fcn=fcon, eps=xstep)

    ! === Default values ===

    call set_optional_arg (meq, 0, lmeq)
    call set_optional_arg (m, lmeq, lm)
    call set_optional_arg (tol, 1.0e-4_PREC, ltol)
    call set_optional_arg (xtol, 1.0e-6_PREC, lxtol)
    call set_optional_arg (exact_lsearch, .false., lexact_lsearch)
    call set_optional_arg (maximize, .false., lmaximize)
    call set_optional_arg (autoscale, .false., lautoscale)
    call set_optional_arg (maxiter, 100, lmaxiter)
    call set_optional_arg (maxfev, huge(1), lmaxfev)
    call set_optional_arg (iprint, NF_PRINT_NONE, liprint)

   ! === Call implementation routine ===

    call slsqp_impl (fobj_wrapper, x, work, fcon_wrapper, lm, lmeq, &
        lbounds, ubounds, ltol, lxtol, lmaxiter, lmaxfev, lexact_lsearch, &
        lmaximize, lautoscale, liprint, lres)

100 continue

    if (present(res)) res = lres

end subroutine



subroutine slsqp_query (m, meq, n, autoscale, w, iwork, sub_nw, sub_niwork)
    !*  SLSQP_QUERY  returns the minimum workspace requirements for the
    !   SLSQP implementation routine.
    integer, intent(in) :: m, meq
        !*  Total number of constraints M, and number of equality constraints MEQ
    integer, intent(in) :: n
        !*  Dimension of solution vector.
    logical, intent(in), optional :: autoscale
        !*  If true, apply autoscaling to problem. This includes rescaling
        !   the problem to perform maximization.
    integer, intent(out) :: w
        !*  Minimum size for the real working array
    integer, intent(out) :: iwork
        !*  Minimum size for the integer working array
    integer, intent(out), optional :: sub_nw, sub_niwork
        !*  Optional minimum working array sizes for subroutines called
        !   from SLSQP (which is only LSQ).

    integer :: la, nw, niwork

    ! Allow for additional space for augmented LSQ problem which will
    ! append one more column
    call lsq_query (m, meq, n+1, nw, niwork)

    if (present(sub_nw)) sub_nw = nw
    if (present(sub_niwork)) sub_niwork = sub_niwork

    ! Additional arrays needed in SLSQP itself
    la = max(1, m)

    nw = nw + la                ! MU
    nw = nw + (n+1)*n/2 + 1     ! L
    nw = nw + n                 ! X0
    nw = nw + 2*n + la          ! R
    nw = nw + (n+1)             ! S
    nw = nw + (n+1)             ! U
    nw = nw + (n+1)             ! V
    nw = nw + m                 ! C
    nw = nw + (n+1)             ! G
    nw = nw + la*(n+1)          ! A

    if (autoscale) then
        ! local copies of lower and upper bounds
        nw = nw + 2 * n
        ! Add X_SCALE and X_ORIG arrays
        nw = nw + 2 * n
    end if

    w = nw
    iwork = niwork

end subroutine



subroutine check_input (has_fcon, m, meq, x, xl, xu, tol, xtol, maxiter, maxfev, &
        exact_lsearch, iprint, xstep, status, msg)
    logical, intent(in) :: has_fcon
    integer, intent(in), optional :: m, meq
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous, optional :: xl, xu
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    logical, intent(in), optional :: exact_lsearch
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    real (PREC), intent(in), optional :: xstep
    type (status_t), intent(out) :: status
    character (*), intent(out), optional :: msg

    integer, parameter :: IPRINT_VALID(*) = [NF_PRINT_NONE, NF_PRINT_MINIMAL, &
        NF_PRINT_VERBOSE, NF_PRINT_ALL]
    integer :: n

    n = size(x)

    status = NF_STATUS_OK

    call check_nonneg (1, m, 'm', status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_nonneg (1, meq, 'meq', status, msg)
    if (status /= NF_STATUS_OK) goto 100

    if (present(m) .and. present(meq)) then
        call check_cond (m >= meq, 'SLSQP', 'Invalid argument: M >= MEQ required', &
            status, msg)
        if (status /= NF_STATUS_OK) goto 100
    end if

    if (has_fcon) then
        if (.not. present(m)) then
            if (present(msg)) msg = "SLSQP: Argument 'm' missing"
            status = NF_STATUS_INVALID_ARG
            goto 100
        else
            call check_positive (1, m, 'm', status, msg)
            if (status /= NF_STATUS_OK) goto 100
        end if
    end if

    if (present(xl)) then
        call check_cond (size(xl) == n, 'SLSQP', &
            'Non-conformable arrays X, LBOUNDS', status, msg)
        if (status /= NF_STATUS_OK) goto 100
    end if

    if (present(xu)) then
        call check_cond (size(xu) == n, 'SLSQP', &
            'Non-conformable arrays X, UBOUNDS', status, msg)
        if (status /= NF_STATUS_OK) goto 100
    end if

    call check_positive (1.0_PREC, tol, 'tol', status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (1.0_PREC, xtol, 'xtol', status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (1, maxiter, 'maxiter', status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (1, maxfev, 'maxfev', status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (1.0_PREC, xstep, 'xstep', status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_enum (iprint, IPRINT_VALID, 'iprint', status, msg)
    if (status /= NF_STATUS_OK) goto 100

100 continue

end subroutine



recursive subroutine slsqp_impl (fobj, x, work, fcon, m, meq, xl, xu, tol, xtol, &
        maxiter, maxfev, exact_lsearch, maximize, autoscale, iprint, res)
    type (fwrapper_vs), intent(inout) :: fobj
    real (PREC), intent(inout), dimension(:), contiguous :: x
    type (workspace), intent(inout), target :: work
    type (fwrapper_vv), intent(inout) :: fcon
    integer, intent(in) :: m, meq
    real (PREC), intent(in), dimension(:), contiguous, optional, target :: xl, xu
    real (PREC), intent(in) :: tol
    real (PREC), intent(in) :: xtol
    integer, intent(in) :: maxiter
    integer, intent(in) :: maxfev
    logical, intent(in) :: exact_lsearch
    logical, intent(in) :: maximize
    logical, value :: autoscale
    integer (NF_ENUM_KIND), intent(in) :: iprint
    type (optim_result), intent(out) :: res

    type (status_t) :: status

    real (PREC), dimension(:), pointer, contiguous :: mu, c, r, l, x0, s, u, v, g
    real (PREC), dimension(:), pointer, contiguous :: x_scale, x_orig, xsl, xsu
    real (PREC), dimension(:), pointer, contiguous :: lsq_w
    real (PREC), dimension(:,:), pointer, contiguous :: a
    integer :: nrwork, niwork, n, nn1, la, lsq_nrwork
    integer :: i, iter, iline, nreset_BFGS
    real (PREC) :: fx, gs, h1, h2, h3, h4, t0, xnorm, rlxtol, fx0
    real (PREC) :: fx_scale
    logical :: bad_linear

    ! Problem dimensions
    n = size(x)
    la = max(1, m)
    ! Size of data block in L corresponding to BFGS matrix;
    nn1 = (n+1) * n / 2

    ! Initialize values used to perform automatic rescaling
    fx_scale = 1.0
    nullify (x_orig, x_scale, xsl, xsu)

    ! === Workspace arrays ===

    call slsqp_query (m, meq, n, autoscale .or. maximize, nrwork, niwork, lsq_nrwork)

    ! Make sure WORKSPACE attributes are allocated to min. required size
    call assert_alloc (work, nrwrk=nrwork, niwrk=niwork)

    ! Debug: Intialize with NaNs to spot errors
    work%rwrk(:) = ieee_value (0.0_PREC, IEEE_SIGNALING_NAN)

    i = 1
    mu => work%rwrk(i:i+la-1)
    i = i + la
    l => work%rwrk(i:i+nn1)
    i = i + nn1 + 1
    x0 => work%rwrk(i:i+n-1)
    i = i + n
    r => work%rwrk(i:i+2*n+la-1)
    i = i + 2*n + la
    s => work%rwrk(i:i+n)
    i = i + n + 1
    u => work%rwrk(i:i+n)
    i = i + n + 1
    v => work%rwrk(i:i+n)
    i = i + n + 1
    c => work%rwrk(i:i+m-1)
    i = i + m
    g => work%rwrk(i:i+n)
    i = i + n + 1
    a(1:la,1:n+1) => work%rwrk(i:i+la*(n+1)-1)
    i = i + la * (n + 1)
    if (autoscale .or. maximize) then
        x_orig(1:n) => work%rwrk(i:i+n-1)
        i = i + n
        x_scale(1:n) => work%rwrk(i:i+n-1)
        i = i + n
    end if

    if (present(xl)) then
        if (autoscale) then
            xsl(1:n) => work%rwrk(i:i+n-1)
            i = i + n
        else
            xsl(1:n) => xl
        end if
    end if

    if (present(xu)) then
        if (autoscale) then
            xsu => work%rwrk(i:i+n-1)
            i = i + n
        else
            xsu => xu
        end if
    end if

    lsq_w => work%rwrk(i:i+lsq_nrwork-1)

    ! === SLSQP implementation ===

    status = NF_STATUS_OK

    ! --- Initial setup ---
    s = 0.0
    mu = 0.0
    nreset_BFGS = 0
    ! Set ILINE such that fused function value/Jacobian calls will not
    ! be triggered in linesearch the first time it is called.
    iline = MAX_INEXACT_LINESEARCH + 1

    call print_section ('Starting main routine', iprint, NF_PRINT_ALL, prefix='SLSQP')

    ! Reset BFGS matrix
    call reset_BFGS ()

    ! Clip initial guess X to satisfy any boundary constraints
    call clip (xl, xu, x)

    ! Compute initial function values and their Jacobians (using unscaled inputs)
    call eval_and_check (fobj, fcon, m, x, fx, c, g(1:n), a(:,1:n), &
        autoscale=.false., status=status, msg=res%msg)
    if (status /= NF_STATUS_OK) goto 100

    ! --- Apply automatic rescaling ---

    if (autoscale) then
        ! Apply the following scaling to objective and input vector:
        !   fs(xs) = f(H(xs))/abs(f(x0))
        ! where
        !   xs = H(x) = x * df(x0)/dx / abs(f(x0))
        ! and the initial scaled guess is given by
        !   xs0 = x0 * df(x0)/dx / abs(f(x0))
        ! and at x0 we have
        !   dfs(xs0)/dxs = df(x0)/dx dx/dH(x) / abs(f(x0))
        !                = df(x0)/dx / df(x0)/dx * abs(f(x0)) / abs(f(x0))
        !                = 1
        fx_scale = 1.0_PREC / max(1.0e-12_PREC, abs(fx))
        if (maximize) then
            ! In addition to autoscaling also change to maximization problem
            fx_scale = -fx_scale
        end if

        do i = 1, n
            x_scale(i) = max(1.0e-12_PREC, abs(g(i))) * abs(fx_scale)
        end do

        call apply_scaling (autoscale, fx_scale, x_scale, fx, x, g(1:n), a(:,1:n))

        if (present(xl)) then
            do i = 1, n
                xsl(i) = xl(i)
            end do
            call apply_scaling (autoscale, fx_scale, x_scale, x=xsl)
        end if

        if (present(xu)) then
            do i = 1, n
                xsu(i) = xu(i)
            end do
            call apply_scaling (autoscale, fx_scale, x_scale, x=xsu)
        end if
    else if (maximize) then
        ! No autoscaling enable, only maximization requested.

        ! Flip sign of objective
        fx_scale = -1.0_PREC
        ! No scaling applied to solution vector
        x_scale(1:n) = 1.0_PREC

        ! Handle remaining adjustsments using general autoscaling code paths
        autoscale = .true.
        call apply_scaling (autoscale, fx_scale, x_scale, fx, x, g(1:n), a(:,1:n))
    end if

    call print_value ('Initial objective gradient: ', g(1:n), iprint, &
        NF_PRINT_ALL, prefix='SLSQP', fmt=FMT_X)

    call print_value ('Initial constr. Jacobian: ', a(:,1:n), iprint, &
        NF_PRINT_ALL, prefix='SLSQP', fmt=FMT_X)

    ! --- Main iteration loop ---

    main: do iter = 1, maxiter

        if (iprint >= NF_PRINT_VERBOSE) then
            write (ERROR_UNIT, '("SLSQP: Iteration #", i0)') iter
        end if

        ! Prepare bounds passed to underlying solver - do this for
        ! every iteration since the arrays U and V get reused below.
        call set_lsq_bounds (x, xsl, u, ieee_value(0.0_PREC, IEEE_NEGATIVE_INF))
        call set_lsq_bounds (x, xsu, v, ieee_value(0.0_PREC, IEEE_POSITIVE_INF))

        h4 = 1.0

        call lsq (meq, l(1:nn1), g(1:n), a(:,1:n), c, u(1:n), v(1:n), &
            s(1:n), r, lsq_w, work%iwrk, status, res%msg)
        ! Admit only these two return states, since we will handle
        ! the invalid state explicitly below.
        if (.not. (status .in. [NF_STATUS_OK, NF_STATUS_INVALID_STATE])) goto 100

        ! --- Augmented LSQ step ---

        ! Augmented problem for inconsistent linearization:
        ! If it turns out that the original QP problem is inconsistent,
        ! disallow termination with convergence on this iteration,
        ! even if the augmented problem was solved.

        select case (status%code_orig)
        case (STATUS_INCOMPAT_CONSTR)
            ! This flag is applied in the Lawson/Hanson module and propagated
            bad_linear = .true.
        case (STATUS_LSEI_SINGULAR_CONSTR)
            bad_linear = (n == meq)
        case default
            bad_linear = .false.
        end select

        if (bad_linear) then

            call print_msg ("SLSQP: Bad linearization. Solving augmented LS problem", &
                iprint, NF_PRINT_VERBOSE)

            call augmented_lsq (meq, x, xsl, xsu, l, g, a, c, u, v, s, r, lsq_w, &
                work%iwrk, status, res%msg)

            if (status  /= NF_STATUS_OK) then
                res%msg = 'SLSQP: Failed to solve (augmented) LSQ problem'
                call print_msg (trim(res%msg), iprint, NF_PRINT_MINIMAL)
                goto 100
            end if

            h4 = 1.0_PREC - s(n+1)

        else if (status /= NF_STATUS_OK) then
            ! Unrecoverable error
            status = NF_STATUS_INVALID_STATE
            goto 100
        end if

        ! Lagrange multipliers of LSQ problem are stored in vector R
        ! We won't need V until the BFGS update, but in the original code
        ! this is done here, probably because we want G, A and R to be
        ! the ones from the LSQ problem, and G and A will be recomputed
        ! right after line search.
        do i = 1, n
            v(i) = g(i) - dot_product (a(:,i), r(1:la))
        end do

        ! --- Update multipliers for L1 test ---

        gs = dot_product (g(1:n), s(1:n))
        h2 = sum_constr (m, meq, c)
        h1 = abs(gs)
        do i = 1, m
            h3 = abs(r(i))
            mu(i) = max(h3, (mu(i) + h3) / 2.0_PREC)
            h1 = h1 + h3 * abs(c(i))
        end do

        ! --- Check convergence ---

        ! Prevent "successful" exit on iteration with bad linearization
        if (h1 < tol .and. h2 < tol .and. .not. bad_linear) then
            status = NF_STATUS_OK
            goto 100
        end if

        h1 = sum_constr (m, meq, c, mu)
        h3 = gs - h1*h4

        if (h3 >= 0.0_PREC) then
            if (iter == 1 .or. nreset_BFGS >= 5) then
                ! Check RELAXED convergence in case of positive directional
                ! derivative if resetting estimate of the Hessian will not
                ! help, or if we exceeded the max. reset count.
                ! Avoid this step in the first iteration, which already
                ! started with an identity matrix!
                h3 = sum_constr (m, meq, c)
                xnorm = NRM2 (n, s, 1)

                ! Note: Original code checks abs(fx-fx0) < tol, but at this
                ! point these are identical since they differ only during
                ! the line search step!
                rlxtol = 10.0_PREC * tol
                if (xnorm < (xtol * 10.0) .and. h3 < rlxtol .and. .not. bad_linear) then
                    status = NF_STATUS_OK
                    status%code_orig = 0
                else
                    res%msg = 'SLSQP: Positive directional derivative'
                    status = NF_STATUS_INVALID_STATE
                    status%code_orig = 8

                    call print_msg (trim(res%msg), iprint, NF_PRINT_MINIMAL)
                end if

                ! Terminate algorithm in any case
                goto 100
            end if

            ! Positive directional derivative in line search
            call reset_BFGS ()
            status%code_orig = 8
            cycle main
        end if

        ! --- Line search ---

        t0 = fx + sum_constr (m, meq, c, mu)
        ! Store old value of objective function
        fx0 = fx

        if (exact_lsearch) then
            ! TODO: implement this
            continue
        else
            ! Inexact line search
            call lsearch_inexact (fobj, fcon, meq, mu, t0, h3, tol, xtol, &
                iline, x(1:n), fx, c, g(1:n), a(:,1:n), s(1:n), x0, &
                autoscale, x_orig, fx_scale, x_scale, iprint, status, res%msg)
        end if

        if (status /= NF_STATUS_OK .and. fobj%nfev >= maxfev) then
            ! Exceeded max. function evaluations
            status = NF_STATUS_MAX_EVAL
            res%msg = 'Max. number of function evaluations exceeded'
            call print_msg (trim(res%msg), iprint, NF_PRINT_MINIMAL, prefix='SLSQP')

            ! Check if old candidate point should be used
            if (fx0 < fx) then
                fx = fx0
                x(:) = x0
            end if
            goto 100

        else if (NF_STATUS_NOT_CONVERGED .in. status) then
            ! Perform BFGS update, move to next iteration

            call print_msg ('Performing BFGS update', iprint, NF_PRINT_ALL, prefix='SLSQP')

            ! Latest gradient/Jacobian were already computed by
            ! linesearch, just need to do the update.
            call update_BFGS (g(1:n), a(:,1:n), r(1:la), s(1:n), u(1:n), &
                v(1:n), l(1:nn1), status)

            if (status /= NF_STATUS_OK) then
                if (nreset_BFGS <= MAX_RESET_BFGS) then
                    call reset_BFGS ()
                else
                    res%msg = 'Failed to perform BFGS update'
                    call print_msg (trim(res%msg), iprint, NF_PRINT_MINIMAL, &
                        indent=2, prefix='SLSQP')
                    goto 100
                end if
            end if

            ! Move not next iteration
            cycle main
        else
            ! Exit routine with either success or error encountered in LS
            goto 100
        end if

    end do main

    if (iter >= maxiter) then
        ! Max. iteration count exceeded
        status = NF_STATUS_MAX_ITER
        status%code_orig = 8
        res%msg = 'Number of max. iterations exceeded'
        call print_msg (trim(res%msg), iprint, NF_PRINT_MINIMAL, prefix='SLSQP')
    end if


100 continue

    ! Undo internal scaling
    call undo_scaling (autoscale, fx_scale, x_scale, fx, x)

    call result_update (res, x, fx, status=status, nit=iter, nfev=fobj%nfev)
    if (status == NF_STATUS_OK) then
        res%msg = 'Optimization terminated successfully'
        call print_msg (trim(res%msg), iprint, NF_PRINT_ALL, prefix='SLSQP')
    end if

    call print_section ('Exiting main routine', iprint, NF_PRINT_ALL, prefix='SLSQP')

    contains

    subroutine reset_BFGS ()
        !*  Reset BFGS matrix which is stored as lower triangular
        !   LDL' factorization in array L.
        !   Overwrites array contents with identity matrix.
        integer :: i, j

        l = 0.0
        j = 1
        do i = 1, n
            l(j) = 1.0
            j = j + (n+1) - i
        end do

        nreset_BFGS = nreset_BFGS + 1
    end subroutine

end subroutine



subroutine clip (xl, xu, x)
    !*  CLIP given vector to satisfy the element-wise boundary constraints.
    real (PREC), intent(in), dimension(:), contiguous, optional :: xl, xu
        !*  Optional boundary constraints. Any non-finite elements in XL and
        !   XU are ignored.
    real (PREC), intent(inout), dimension(:), contiguous :: x

    if (present(xl)) then
        where (isfinite (xl))
            x = max(x, xl)
        end where
    end if

    if (present(xu)) then
        where (isfinite (xu))
            x = min(x, xu)
        end where
    end if
end subroutine



subroutine set_lsq_bounds (x, src, dst, default)
    !*  SET_LSQ_BOUNDS constructs a vector of lower or upper bounds from
    !   optional user input and given default values in the format
    !   required by the routine LSQ.
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous, optional :: src
        !*  Bounds provided by user. Could be missing or non-finite if
        !   not applicable.
    real (PREC), intent(out), dimension(:), contiguous :: dst
    real (PREC), intent(in) :: default
        !*  Should be either -inf or +inf, depending on whether the lower
        !   or the upper bounds are computed.

    integer :: i

    dst = default
    if (present(src)) then
        do i = 1, min(size(src), size(dst))
            if (isfinite (src(i))) then
                dst(i) = src(i) - x(i)
            end if
        end do
    end if
end subroutine



subroutine apply_scaling (autoscale, fx_scale, x_scale, fx, x, g, a)
    !*  APPLY_SCALING rescale the problem (objective, solution vector,
    !   gradient, Jocobian of constraints) given the scaling data.
    logical, intent(in) :: autoscale
        !*  If true, apply scaling
    real (PREC), intent(in) :: fx_scale
        !*  Scaling factor for objective function such that
        !   fx_scaled = fx * fx_scale
    real (PREC), intent(in), dimension(:), contiguous, optional :: x_scale
        !*  Scaling factors for each dimension of the solution vector
        !   such that x_scaled(i) = x * x_scale(i) for each i.
    real (PREC), intent(inout), optional :: fx
        !*  If present, contains the unscaled objective function on entry
        !   and the scaled objective function on exit.
    real (PREC), intent(inout), dimension(:), contiguous, optional :: x
        !*  If present, contains the unscaled solution vector on entry
        !   and the scaled solution vector on exit.
    real (PREC), intent(inout), dimension(:), contiguous, optional :: g
        !*  If present, contains the unscaled gradient of the objective
        !   on entry and the scaled version on exit.
    real (PREC), intent(inout), dimension(:,:), contiguous, optional :: a
        !*  If present, contains the unscaled Jacobian of constraints
        !   on entry and the scaled version on exit.

    integer :: i

    if (autoscale) then
        if (present(fx)) then
            fx = fx * fx_scale
        end if

        if (present(x) .and. present(x_scale)) then
            x = x * x_scale
        end if

        if (present(g) .and. present(x_scale)) then
            g = g / x_scale * fx_scale
        end if

        if (present(a) .and. present(x_scale)) then
            do i = 1, size(x_scale)
                a(:,i) = a(:,i) / x_scale(i)
            end do
        end if
    end if
end subroutine



subroutine undo_scaling (autoscale, fx_scale, x_scale, fx, x, g)
    !*  UNDO_SCALING converts data of the scaled problem back into the
    !   the unscaled version.
    logical, intent(in) :: autoscale
        !*  If true, undo scaling.
    real (PREC), intent(in) :: fx_scale
        !*  Scaling factor for objective function such that
        !   fx_scaled = fx * fx_scale
    real (PREC), intent(in), dimension(:), contiguous, optional :: x_scale
        !*  Scaling factors for each dimension of the solution vector
        !   such that x_scaled(i) = x * x_scale(i) for each i.
    real (PREC), intent(inout), optional :: fx
        !*  If present, contains the scaled objective on entry and the
        !   unscaled version of exit.
    real (PREC), intent(inout), dimension(:), contiguous, optional :: x
        !*  If present, contains the scaled solution vector on entry and
        !   the unscaled version on exit.
    real (PREC), intent(inout), dimension(:), contiguous, optional :: g
        !*  If present, contains the scaled gradient on entry and the unscaled
        !   version on exit.

    if (autoscale) then
        if (present(fx)) then
            fx = fx / fx_scale
        end if

        if (present(x) .and. present(x_scale)) then
            x = x / x_scale
        end if

        if (present(g) .and. present(x_scale)) then
            g = g * x_scale / fx_scale
        end if
    end if
end subroutine



subroutine augmented_lsq (meq, x, xl, xu, l, g, a, c, u, v, s, r, rwork, &
        iwork, status, msg)
    !*  AUGMENTED_LSQ sets up the array arguments required to performed
    !   an augmented LSQ call.
    integer, intent(in) :: meq
    real (PREC), intent(in), dimension(:), contiguous, optional :: x, xl, xu
    real (PREC), intent(inout), dimension(:), contiguous :: l, g
    real (PREC), intent(inout), dimension(:,:), contiguous :: a
    real (PREC), intent(inout), dimension(:), contiguous :: c, u, v, s, r
    real (PREC), intent(inout), dimension(:), contiguous :: rwork
        !*  Real working array passed to LSQ
    integer, intent(inout), dimension(:), contiguous :: iwork
        !*  Integer working array passed to LSQ
    type (status_t), intent(out) :: status
    character (*), intent(out) :: msg

    integer :: i, m, n, nn1
    integer, parameter :: MAX_LSQ = 5
    logical :: bad_linear

    m = size(c)
    ! G is passed in with size (n+1) to store the augmented problem
    n = size(g) - 1

    nn1 = (n+1)*n / 2

    do i = 1, m
        if (i <= meq) then
            a(i,n+1) = -c(i)
        else
            a(i,n+1) = max(-c(i), 0.0_PREC)
        end if
    end do

    g(n+1) = 0.0
    ! Add element packed L which represents an additional column in the
    ! unpacked L with all zeros except for the last row.
    l(nn1+1) = 100.0_PREC

    ! Update bounds
    call set_lsq_bounds (x, xl, u, ieee_value(0.0_PREC, IEEE_NEGATIVE_INF))
    call set_lsq_bounds (x, xu, v, ieee_value(0.0_PREC, IEEE_POSITIVE_INF))

    u(n+1) = 0.0
    v(n+1) = 1.0

    do i = 1, MAX_LSQ
        ! Attempt 5 tries to solve augmented system
        call lsq (meq, l, g, a, c, u, v, s, r, rwork, iwork, status, msg)
        if (.not. (status .in. [NF_STATUS_OK, NF_STATUS_INVALID_STATE])) return

        select case (status%code_orig)
        case (STATUS_INCOMPAT_CONSTR)
            ! This flag is applied in the Lawson/Hanson module and propagated
            bad_linear = .true.
        case (STATUS_LSEI_SINGULAR_CONSTR)
            bad_linear = (n == meq)
        case default
            bad_linear = .false.
        end select

        ! Augmented problem solved without issues
        if (bad_linear) then
            l(nn1+1) = 10.0_PREC * l(nn1+1)
            cycle
        else if ((status /= NF_STATUS_OK) .or. (i == MAX_LSQ)) then
            ! Unrecoverable error, or we exhausted max LSQ iterations
            status = NF_STATUS_INVALID_STATE
            return
        else
            return
        end if
    end do

end subroutine



function sum_constr (m, meq, c, mult) result(h)
    !*  SUM_CONSTR returns the sum over all constraints, optionally
    !   weighted by their Lagrange multipliers. For equality constraints,
    !   the abs. constraint values are used. For inequality constraints,
    !   non-binding constraints are ignored.
    integer, intent(in) :: m, meq
    real (PREC), intent(in), dimension(:), contiguous :: c
    real (PREC), intent(in), dimension(:), contiguous, optional :: mult
    real (PREC) :: h

    real (PREC) :: val
    integer :: i

    h = 0.0
    do i = 1, m
        if (i <= meq) then
            val = abs(c(i))
        else
            val = max(0.0_PREC, -c(i))
        end if

        if (present(mult)) val = val * mult(i)

        h = h + val
    end do
end function



recursive subroutine eval_and_check (fobj, fcon, m, x, fx, c, g, a, &
        autoscale, x_orig, fx_scale, x_scale, status, msg)
    !*  Evaluate objective function values and/or their derivatives,
    !   and, check results for NaNs
    type (fwrapper_vs) :: fobj
    type (fwrapper_vv) :: fcon
    integer, intent(in) :: m
        !*  Number of equality and inequality constraints, excluding
        !   any boundary constraints.
    real (PREC), intent(in), dimension(:), contiguous, target :: x
        !*  Point at which to evaluate objective and constraints
    real (PREC), intent(out), optional :: fx
        !*  Objective function value
    real (PREC), intent(out), dimension(:), contiguous, optional :: c
        !*  Constraint values
    real (PREC), intent(out), dimension(:), contiguous, optional :: g
        !*  Gradient of objective function
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: a
        !*  Jacobian of constraints function
    logical, intent(in) :: autoscale
        !*  If true, apply autoscaling transformation to inputs and outputs
    real (PREC), intent(inout), dimension(:), contiguous, optional, target :: x_orig
        !*  Temporary array to store the unscaled candidate solution
    real (PREC), intent(in), optional :: fx_scale
        !*  Scaling factor for objective function
    real (PREC), intent(in), dimension(:), contiguous, optional :: x_scale
        !*  Scaling factors for each dimension of the solution vector
    type (status_t), intent(out) :: status
    character (*), intent(out) :: msg

    logical :: do_fcn, do_jac
    real (PREC), dimension(:), pointer, contiguous :: ptr_x

    do_fcn = present(fx) .and. present(c)
    do_jac = present(g) .and. present(a)

    status = NF_STATUS_OK

    if (autoscale .and. present(x_orig)) then
        x_orig = x
        call undo_scaling (autoscale, fx_scale, x_scale, x=x_orig)
        ptr_x => x_orig
    else
        ptr_x => x
    end if

    if (do_fcn .and. do_jac) then
        ! --- Evaluate both function values and Jacobians ---
        ! Objective
        call dispatch_fcn_jac (fobj, ptr_x, fx, g)

        if (.not. isfinite (fx) .or. .not. all(isfinite(g))) then
            msg = MSG_OBJECTIVE_NONFINITE
            status = NF_STATUS_INVALID_STATE
        end if

        ! Constraints
        if (m > 0) then
            call dispatch_fcn_jac (fcon, ptr_x, c, a)
            if (.not. all(isfinite(c)) .or. .not. all(isfinite(a))) then
                msg = MSG_CONSTRAINTS_NONFINITE
                status = NF_STATUS_INVALID_STATE
            end if
        end if
    else if (do_fcn) then
        ! --- Function values only ---
        ! Objective
        call dispatch (fobj, ptr_x, fx)
        if (.not. isfinite(fx)) then
            msg = MSG_OBJECTIVE_NONFINITE
            status = NF_STATUS_INVALID_STATE
        end if

        ! Constraints
        if (m > 0) then
            call dispatch (fcon, ptr_x, c)
            if (.not. all(isfinite(c))) then
                msg = MSG_CONSTRAINTS_NONFINITE
                status = NF_STATUS_INVALID_STATE
            end if
        end if
    else
        ! --- Jacobians only ---
        ! Objective
        call dispatch_jac (fobj, ptr_x, g)
        if (.not. all(isfinite(g))) then
            msg = MSG_OBJECTIVE_NONFINITE
            status = NF_STATUS_INVALID_STATE
        end if

        ! Constraints
        if (m > 0) then
            call dispatch_jac (fcon, ptr_x, a)
            if (.not. all(isfinite(a))) then
                msg = MSG_CONSTRAINTS_NONFINITE
                status = NF_STATUS_INVALID_STATE
            end if
        end if
    end if

    ! Unscaled evaluations may omit FX_SCALE; do not pass it to a required dummy.
    if (autoscale .and. status /= NF_STATUS_INVALID_STATE) then
        call apply_scaling (autoscale, fx_scale, x_scale, fx, g=g, a=a)
    end if

end subroutine



recursive subroutine lsearch_inexact (fobj, fcon, meq, mu, t0, h3, tol, xtol, &
        iline, x, fx, c, g, a, s, x0, autoscale, x_orig, fx_scale, x_scale, &
        iprint, status, msg)
    type (fwrapper_vs) :: fobj
    type (fwrapper_vv) :: fcon
    integer, intent(in) :: meq
    real (PREC), intent(in), dimension(:), contiguous :: mu
    real (PREC), intent(in) :: t0
    real (PREC), value :: h3
    real (PREC), intent(in) :: tol
    real (PREC), intent(in) :: xtol
    integer, intent(inout) :: iline
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(inout) :: fx
    real (PREC), intent(out), dimension(:), contiguous :: c, g
    real (PREC), intent(out), dimension(:,:), contiguous :: a
    real (PREC), intent(out), dimension(:), contiguous :: s
    real (PREC), intent(out), dimension(:), contiguous :: x0
    logical, intent(in) :: autoscale
    real (PREC), intent(inout), dimension(:), contiguous, optional :: x_orig
    real (PREC), intent(in) :: fx_scale
    real (PREC), intent(in), dimension(:), contiguous, optional :: x_scale
    integer (NF_ENUM_KIND), intent(in) :: iprint
    type (status_t), intent(out) :: status
    character (*), intent(out) :: msg

    integer :: n, m
    real (PREC) :: alpha, t, xnorm, h1, fx0, dfx, alpha_total
    integer, parameter :: INDENT = 2

    integer, parameter :: MAXITER = MAX_INEXACT_LINESEARCH + 1
        !*  Increment by +1. See comment at the beginning of the loop.
    logical :: fuse_calls, has_jac

    status = NF_STATUS_OK

    call print_msg ('Starting inexact LINE SEARCH', iprint, NF_PRINT_ALL, prefix='SLSQP')

    m = size(c)
    n = size(x)

    alpha = 1.0
    ! Keep track of total shrinking factor applied over several line search
    ! iterations.
    alpha_total = 1.0

    ! Fuse calls to evaluate function values and gradients if last
    ! linesearch step took only one iteration
    fuse_calls = (iline == 1)
!    fuse_calls = .false.

    ! Store function value and solution vector prior to beginning line
    ! search
    fx0 = fx
    x0 = x

    ! Note: in the original code the line search loop gets evaluated an
    ! ADDITIONAL time since the termination condition is only checked
    ! AFTER the function/constraints evaluation.
    ! Set the constant MAXITER accordingly above!
    do iline = 1, MAXITER

        call print_value ('Iteration #', iline, iprint, NF_PRINT_ALL, &
            prefix='SLSQP', indent=INDENT)

        h3 = alpha * h3
        s = s * alpha
        alpha_total = alpha_total * alpha
        x = x0 + s

        call print_value ('Relative step size: ', alpha_total, iprint, &
            NF_PRINT_ALL, prefix='SLSQP', indent=INDENT+2, fmt='es9.2e2')
        call print_value ('Evaluating at X: ', x, iprint, NF_PRINT_ALL, &
            prefix='SLSQP', indent=INDENT+2, fmt=FMT_X)

        ! Evaluate objective and constraint functions
        if (fuse_calls .or. iline == MAXITER) then
            call eval_and_check (fobj, fcon, m, x, fx, c, g, a, autoscale, &
                x_orig, fx_scale, x_scale, status, msg)
            has_jac = .true.
        else
            call eval_and_check (fobj, fcon, m, x, fx, c, autoscale=autoscale, &
                x_orig=x_orig, fx_scale=fx_scale, x_scale=x_scale, &
                status=status, msg=msg)
            has_jac = .false.
        end if

        call print_value ('Objective: ', fx, iprint, NF_PRINT_ALL, &
            prefix='SLSQP', indent=INDENT+2, fmt=FMT_FX)
        call print_value ('d(Objective): ', fx-fx0, iprint, NF_PRINT_ALL, &
            prefix='SLSQP', indent=INDENT+4, fmt=FMT_DIFF)
        call print_value ('Constraints: ', c, iprint, NF_PRINT_ALL, &
            prefix='SLSQP', indent=INDENT+2, fmt=FMT_X)

        ! Error while evaluating objective or constraints
        if (status /= NF_STATUS_OK) goto 100

        t = fx + sum_constr (m, meq, c, mu)
        h1 = t - t0

        if (h1 <= h3/10.0_PREC .or. iline == MAXITER) then
            ! Check convergence
            h3 = sum_constr (m, meq, c)
            xnorm = NRM2 (n, s, 1)
            dfx = fx - fx0
            if ((abs(dfx) < tol .or. xnorm < xtol) .and. h3 < tol) then
                ! Line search found a satisfactory minimum.
                status = NF_STATUS_OK
                ! Set original SLSQP status code.
                status%code_orig = 0

                ! Check if old candidate point should be used
                if (fx0 < fx) then
                    fx = fx0
                    x = x0

                    call print_msg ('No improvement found', iprint, NF_PRINT_ALL, &
                        prefix='SLSQP', indent=INDENT)
                else
                    call print_msg ('Convergence achieved', iprint, NF_PRINT_ALL, &
                        prefix='SLSQP', indent=INDENT)
                end if

                goto 100
            else
                ! Update Jacobian of objective and constraints such
                ! that on non-converged termination Jacobian is guaranteed
                ! to be updated
                if (.not. has_jac) then
                    call eval_and_check (fobj, fcon, m, x, g=g, a=a, &
                        autoscale=autoscale, x_orig=x_orig, fx_scale=fx_scale, &
                        x_scale=x_scale, status=status, msg=msg)
                    if (status /= NF_STATUS_OK) goto 100
                end if

                call print_value ('Updated objective gradient: ', g(1:n), iprint, &
                    NF_PRINT_ALL, prefix='SLSQP', indent=INDENT+2, fmt=FMT_X)

                call print_value ('Updated constr. Jacobian: ', a(:,1:n), iprint, &
                    NF_PRINT_ALL, prefix='SLSQP', indent=INDENT+2, fmt=FMT_X)

                ! Signal to caller that line search has not found
                ! a satisfactory minimum
                status = NF_STATUS_OK + NF_STATUS_NOT_CONVERGED
                goto 100
            end if
        else
            ! Rescale line search distance
            alpha = max(h3/(2.0_PREC * (h3 - h1)), 1.0e-1_PREC)
            call print_value ('Backtracking #', iline, iprint, NF_PRINT_MINIMAL, &
                prefix='SLSQP', indent=INDENT)
            call print_value ('Rescaling step size by ', alpha, iprint, &
                NF_PRINT_MINIMAL, prefix='SLSQP', indent=INDENT+2, fmt='es9.2e2')
        end if
    end do

100 continue

    call print_msg ('Exiting inexact LINE SEARCH', iprint, NF_PRINT_ALL, prefix='SLSQP')

end subroutine



subroutine update_BFGS (g, a, r, s, u, v, l, status)
    !*  UPDATE_BFGS updates the Cholesky factors of Hessian matrix by modified
    !   BFGS formula.
    real (PREC), intent(in), dimension(:), contiguous :: g
    real (PREC), intent(in), dimension(:,:), contiguous :: a
    real (PREC), intent(in), dimension(:), contiguous :: r, s
    real (PREC), intent(inout), dimension(:), contiguous :: u, v
    real (PREC), intent(inout), dimension(:), contiguous :: l
        !*  On entry, contains current factorization LDL'. On termination,
        !   contains updated LDL' factorization. Stores only rows
        !   i <= j for each column j in packed column-major format.
    type (status_t), intent(out) :: status

    integer :: i, j, k, n
    real (PREC) :: h1, h2, h3, h4

    status = NF_STATUS_OK

    n = size(g)

    do i = 1, n
        u(i) = g(i) - dot_product (a(:,i), r) - v(i)
    end do

    ! L' * s
    k = 0
    do i = 1, n
        h1 = 0.0
        k = k + 1
        do j = i + 1, n
            k = k + 1
            h1 = h1 + l(k)*s(j)
        end do
        v(i) = s(i) + h1
    end do

    ! D * L' * s
    k = 1
    do i = 1, n
        v(i) = l(k) * v(i)
        k = k + (n + 1) - i
    end do

    ! L*D*L' * s
    do i = n, 1, -1
        h1 = 0.0
        k = i
        do j = 1, i-1
            h1 = h1 + l(k)*v(j)
            k = k + n - j
        end do
        v(i) = v(i) + h1
    end do

    h1 = dot_product (s, u)
    h2 = dot_product (s, v)
    h3 = 0.2_PREC * h2
    if (h1 < h3) then
        h4 = (h2 - h3) / (h2 - h1)
        h1 = h3
        do i = 1, n
            u(i) = (1.0_PREC - h4) * v(i) + h4 * u(i)
        end do
    end if

    ! If BFGS update fails just exit, calling code will cycle main loop
    ! anyways.
    if (h1 == 0.0_PREC .or. h2 == 0.0_PREC) then
        status = NF_STATUS_INVALID_STATE
        return
    end if

    call ldl (n, l, u,  1.0_PREC/h1, v)
    call ldl (n, l, v, -1.0_PREC/h2, u)

end subroutine



subroutine lsq_query (m, meq, n, w, iwork, sub_nw, sub_niwork)
    !*  LSQ_QUERY returns the required minimum workspace array sizes for
    !   the routine LSQ, for given problem dimensions.
    integer, intent(in) :: m, meq
        !*  Total number of constraints M, and number of equality constraints MEQ
    integer, intent(in) :: n
        !*  Dimension of solution vector
    integer, intent(out) :: w, iwork
        !*  Minimum sizes for real working array (W) and integer working array
        !   (IWORK).
    integer, intent(out), optional :: sub_nw, sub_niwork
        !*  Optional sizes of working arrays required by routines called
        !   by LSQ.

    integer :: mg, mc

    ! Number of inequality constraints
    mg = m - meq + 2 * n
    mc = meq

    ! LSEI called by LSQ, add to workspace requirement
    call lsei_query (n, n, mg, mc, w, iwork)

    if (present(sub_nw)) sub_nw = w
    if (present(sub_niwork)) sub_niwork = iwork

    ! Additional requirements for LSQ itself
    w = w + (3*n+m) * (n + 1)

end subroutine



subroutine lsq (meq, l, g, a, b, xl, xu, x, y, w, iwork, status, msg)
    integer, intent(in) :: meq
    real (PREC), intent(in), dimension(:), contiguous :: l
    real (PREC), intent(in), dimension(:), contiguous :: g
    real (PREC), intent(in), dimension(:,:), contiguous :: a
    real (PREC), intent(in), dimension(:), contiguous :: b
    real (PREC), intent(in), dimension(:), contiguous, optional :: xl, xu
        !*  Arrays of length N containing lower and upper bounds
    real (PREC), intent(out), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: y
        !*  On exit, contains M + N + N Lagrange multipliers for
        !   M (in)equality containts, N lower and N upper bounds.
    real (PREC), intent(out), dimension(:), contiguous, target :: w
        !*  Worspace array
    integer, intent(out), dimension(:), contiguous :: iwork
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    integer :: m, mineq, mineq_bnd, nw, n, niwork, lsei_nw, nlb, nub
    integer :: ic, id, ig, ih, j
    real (PREC) :: xnorm
    real (PREC), dimension(:), pointer, contiguous :: lsei_w, vecf, vecd, vech
    real (PREC), dimension(:,:), pointer, contiguous :: matE, matC, matG
    type (status_t) :: lstatus
    logical :: cond

    lstatus = NF_STATUS_OK

    ! Length of (potentially augmented) solution vector
    n = size(x)
    ! Number of equality and inequality constraints (excl. boundary constraints)
    m = size(a, 1)
    ! Number of inequality constraints
    mineq = m - meq

    ! === Workspace ===

    call lsq_query (m, meq, n, nw, niwork, lsei_nw)

    ! == Check inputs ===

    ! Check that l has the correct size for the applicable problem
    ! (non-augmented or augmented).
    ! Augmented problem has one additional (diagonal) element, but
    ! N is incremented by one, so L will have size n*(n-1)/2 + 1
    cond = (size(l) == (n+1)*n/2) .or. (size(l) == (n*(n-1)/2 + 1))
    call check_cond (cond, 'LSQ', 'Non-conformable array L', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    cond = (size(w) >= nw)
    call check_cond (cond, 'LSQ', 'Workspace array W too small', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    cond = (size(iwork) >= niwork)
    call check_cond (cond, 'LSQ', 'Workspace array IWORK too small', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! === LSQ implementation ===

    ! Workspace for LSEI assumed to be at very end
    lsei_w(1:lsei_nw) => w(nw-lsei_nw+1:nw)

    ! --- Recover matrix E and vector f ---
    ! Assume that the [N,N] matrix E and f with size N are the
    ! first two elements in W
    nullify (matE, vecf)
    call recover_ef (w(1:n*n), w(n*n+1:n*n+n))

    ! --- Recover matrix C ---
    ! C with dimension [MEQ,N] is stored in upper part of A
    ! Initial index of C on W
    ic = n*n + n + 1
    matC(1:meq,1:n) => w(ic:ic+meq*n-1)
    ! Copy non-contiguous memory in A by column (original code copies by row)
    do j = 1, n
        matC(1:meq,j) = a(1:meq,j)
    end do

    ! --- Recover vector d ---
    id = ic + meq * n
    vecd(1:meq) => w(id:id+meq-1)
    vecd(1:meq) = - b(1:meq)

    ! --- Recover matrix G and vector h ---
    ! G with dimensions of at [M-MEQ,N] is stored in lower part of A.
    ! We augment G by the number of actual (finite) lower and upper bounds.
    ig = id + meq

    ! Determine number of actual lower bounds
    nlb = 0
    if (present(xl)) then
        nlb = count (isfinite (xl))
    end if

    ! Determine number of actual upper bounds
    nub = 0
    if (present(xu)) then
        nub = count (isfinite (xu))
    end if

    mineq_bnd = mineq + nlb + nub

    matG(1:mineq_bnd,1:n) => w(ig:ig+mineq_bnd*n-1)
    do j = 1, n
        matG(1:mineq,j) = a(meq+1:m,j)
    end do

    ! h is stored in lower part of array B
    ! We augment h by appending all actual lower and upper bounds
    ih = ig + mineq_bnd*n
    vech(1:mineq_bnd) => w(ih:ih+mineq_bnd-1)
    vech(1:mineq) = - b(meq+1:m)

    ! --- Append bounds inequalities ---

    ! Append boundary inequality constraints to G and h
    call append_bounds ()

    ! --- Solve constrained least squares problem ---

    if (meq > 0) then
        ! If there are both equality and inequality constraints, call
        ! LSEI -> LSI -> ...
        call lsei (matC, vecd, matE, vecf, matG, vech, x, xnorm, lsei_w, &
            iwork, lstatus, msg)
    else
        ! If there are NO equality constraints, skip LSEI which just
        ! eliminates non-present equality constraints, and call LSI directly.
        call lsi (matE, vecf, matG, vech, x, xnorm, lsei_w, iwork, lstatus, msg)
    end if
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Restore Lagrange multipliers (ignoring those on boundary constraints)
    ! LSEI/LSI store these as the first elements in its working array.
    y(1:m) = lsei_w(1:m)

    ! Set remaining multipliers to NaN (they are not used)
    y(m+1:) = ieee_value (0.0_PREC, IEEE_QUIET_NAN)

    ! Enforce bounds
    call clip (xl, xu, x)

100 continue

    if (present(status)) status = lstatus

    contains

    subroutine recover_ef (e, f)
        real (PREC), intent(out), dimension(:), contiguous, target :: e
        real (PREC), intent(out), dimension(:), contiguous, target :: f

        integer :: nn1, illb, ilub, i, norig
        real (PREC) :: diag, tmp
        logical :: augmented

        nn1 = (n+1)*n/2
        ! Detect if this is the augmented LSQ problem from the size of L
        select case(n)
        case (1)
            augmented = (size(l) == 2)
        case default
            augmented = (size(l) == (n*(n-1)/2 + 1))
        end select

        ! Assumed dimension of square matrix stored in L.
        if (augmented) then
            norig = n - 1
        else
            norig = n
        end if

        matE(1:n,1:n) => e
        vecf(1:n) => f

        matE(:,:) = 0.0_PREC
        vecf(:) = 0.0_PREC

        ! Indices of current column of lower-triangular matrix in L
        illb = 1
        ilub = norig

        do i = 1, norig
            diag = sqrt(l(illb))
            ! Write row which is the transpose of the ith column of L
            matE(i,i+1:norig) = l(illb+1:ilub)
            matE(i,i) = 1.0
            ! Pre-multiplying with diagonal matrix corresponds to
            ! rescaling i-th row by diag(i,i)
            matE(i,:) = matE(i,:) * diag

            ! Update element f(i)
            tmp = g(i) - dot_product (matE(1:i-1,i), vecf(1:i-1))
            vecf(i) =  tmp / matE(i,i)

            ! Update indices pointing to start and end point of
            ! (partial) column on L
            illb = ilub + 1
            ilub = ilub + (norig-i)
        end do

        vecf = - vecf

        ! Process last columns separately
        if (augmented) then
            illb = (norig+1)*norig / 2 + 1
            matE(n,n) = l(illb)
        end if

        ! Set any augmented block to zero
        vecf(norig+1:n) = 0.0

    end subroutine


    subroutine append_bounds ()
        !*  Append finite boundary constraints to inequality constraints.

        integer :: i
        ! Append data for lower bounds
        i = 0
        if (present(xl)) then
            do j = 1, n
                if (.not. isfinite (xl(j))) cycle
                i = i + 1
                matG(mineq+i,:) = 0.0
                matG(mineq+i,j) = 1.0
                vech(mineq+i) = xl(j)
            end do
        end if

        ! Append data for upper bounds (convert to -x >= -xb constraints)
        ! i at this point contains the number of appended lower bounds
        if (present(xu)) then
            do j = 1, n
                if (.not. isfinite (xu(j))) cycle
                i = i + 1
                matG(mineq+i,:) = 0.0
                matG(mineq+i,j) = -1.0
                vech(mineq+i) = -xu(j)
            end do
        end if
    end subroutine

end subroutine




subroutine lsei_query (m, n, mg, mc, w, iwork, sub_nw, sub_niwork)
    !*  LSEI_QUERY returns the minimal workspace array dimensions
    !   required by LSEI.
    integer, intent(in) :: m, n
        !*  Dimensions of matrix A in LS problem Ax = b
    integer, intent(in) :: mg
        !*  First dimension of constraint matrix G in Gx >= h
    integer, intent(in) :: mc
        !*  First dimension of constraint matrix C in Cx >= d
    integer, intent(out) :: w, iwork
        !*  Minimum dimenions of W and IWORK workspace arrays.
    integer, intent(out), optional :: sub_nw, sub_niwork
        !*  Optional maximum of dimensions of workspace arrays required
        !   for subroutines called by LSEI.

    integer :: l, hfti_w, hfti_iwork

    l = max(0, n-mc)
    call lsi_query (m, l, mg, w, iwork)
    ! HFTI will almost surely require a smaller working array, but include
    ! this for completeness
    call hfti_query (m, l, hfti_w, hfti_iwork)

    w = max(hfti_iwork, w)
    iwork = max(iwork, hfti_iwork)

    ! Store max. dimensions required by subroutines, assuming that
    ! these need not be stored at the same time.
    if (present(sub_nw)) sub_nw = w
    if (present(sub_niwork)) sub_niwork = iwork

    ! Additional workspace requirements for LSEI itself
    !   2 * MC + M + (M + MG) * (N - MC)
    w = w + 2 * mc + m + (m + mg) * (n- mc)
end subroutine



subroutine lsei (c, d, e, f, g, h, x, xnorm, w, iwork, status, msg)
    !*  LSEI computes the solution to a least squares problem with
    !   equality and inequality constraints given by
    !       min_{x}  ||E*x - f||
    !           s.t.   C*x = d, G*x >= h
    !
    !   using the approach from section 23.6 of
    !       Lawson, Hanson (1995): Solving least squares problems.
    !   Original implementation in Fortran 77 by Dieter Kraft.
    real (PREC), intent(inout), dimension(:,:), contiguous :: c
        !*  Coefficient matrix C with dimensions [MC,N]. Routine requires
        !   that MC <= N.
    real (PREC), intent(inout), dimension(:), contiguous :: d
        !*  Vector d of size MC.
    real (PREC), intent(inout), dimension(:,:), contiguous :: e
        !*  Coefficient matrix E with dimensions [M,N].
    real (PREC), intent(inout), dimension(:), contiguous :: f
        !*  Vector f of size M.
    real (PREC), intent(inout), dimension(:,:), contiguous :: g
        !*  Coefficient matrix G with dimensions [MG,N].
    real (PREC), intent(inout), dimension(:), contiguous :: h
        !*  Vector h of size MG.
    real (PREC), intent(out), dimension(:), contiguous :: x
        !*  Stores solution vector X on successful termination.
    real (PREC), intent(out) :: xnorm
        !*  Stores Euclidean norm of solution.
    real (PREC), intent(out), dimension(:), contiguous, target :: w
        !*  Workspace array. The first MC + MG elements store the
        !   Lagrange multipliers for the equality and inequality constraints.
    integer, intent(out), dimension(:), contiguous :: iwork
        !*  Integer workspace array.
    type (status_t), intent(out), optional :: status
        !*  Optional status code. Returns NUMFORT status codes as well as
        !   the following numerical codes of the original implementation
        !   in the CODE_ORIG attribute of STATUS_T:
        !       1.  NF_STATUS_OK on successful termination;
        !       2.  NF_STATUS_INVALID_ARG if any of the input arrays has
        !           an invalid dimensions.
        !       3.  NF_STATUS_MAX_ITER if the max. iteration count in NNLS
        !           was exceeded.
        !       4.  NF_STATUS_INVALID_STATE if the inequality constraints
        !           are incompatible.
        !       5.  NF_STATUS_INVALID_STATE if the matrix E is not of full
        !           rank.
        !       6.  NF_STATUS_INVALID_STATE if the matrix C is not of full
        !           rank.
        !       7.  NF_STATUS_INVALID_STATE on rank defects in HFTI.
    character (*), intent(out), optional :: msg
        !*  Optional error message. Not modified on successful termination.

    type (status_t) :: lstatus
    real (PREC), dimension(:), pointer, contiguous :: sub_w, h12_h, ptr_b
    real (PREC), dimension(:), pointer, contiguous :: mult_eq, mult_ineq
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_a, ptr_g
    real (PREC) :: tau, t
    integer :: i
    integer :: n, m, mg, mc, l, nw, niwork, krank
    integer :: sub_nw
    logical :: cond

    lstatus = NF_STATUS_OK
    lstatus%code_orig = STATUS_OK

    m = size(e, 1)
    n = size(e, 2)
    mg = size(g, 1)
    mc = size(c, 1)
    l = n - mc

    ! === Workspace requirements ===

    call lsei_query (m, n, mg, mc, nw, niwork, sub_nw)

    ! === Check inputs ===

    call check_cond (m == size(f), 'LSEI', 'Non-conformable arrays E, F', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (mg == size(h), 'LSEI', 'Non-conformable arrays G, H', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (mc == size(d), 'LSEI', 'Non-conformable arrays C, D', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    cond = (n == size(g, 2) .and. n == size(c, 2))
    call check_cond (cond, 'LSEI', 'Non-conformable arrays E, G, C', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (n == size(x), 'LSEI', 'Non-conformable arrays E, X', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(w) >= nw, 'LSEI', 'Array W too small', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(iwork) >= niwork, 'LSEI', 'Array IWORK too small', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Cannot have more equality constraints than N
    call check_cond (mc <= n, 'LSEI', 'Too many equality constraints', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! === Working arrays ===

    ! Working array layout in the original code, where L = N - MC:
    !   MC                  elements: Langrange multipliers on equality constraints
    !   (L+1)*(MG+2) + 2*MG elements: Workspace for LSI
    !                                   - MC1 points to first element
    !                                   - IW points to last element
    !   MC                  elements: Used to store H argument to H12
    !   M*L                 elements: Used to store modified matrix A
    !                                   - IE points to first element
    !   M                   elements: Stores modified B
    !                                   - IF points to first element
    !   MG*L                elements: Stores modified G
    !                                   - IG points to first element

    ! Vector of Lagrange multipliers for equality constraints, used only
    ! at the very end
    i = 1
    mult_eq(1:mc) => w(i:i+mc-1)
    ! Array used to store Lagrange multipliers for inequality constraints
    ! NOTE THAT THIS OVERLAPS WITH WORKING ARRAY FOR LSI ON PURPOSE,
    ! because LSI returns these in the first MG elements of its working array.
    i = i + mc
    mult_ineq(1:mg) => w(i:i+mg-1)
    ! Working array for LSI or HFTI; only at most one of these will be
    ! called in this routine, so this block is shared between both.
    sub_w(1:sub_nw) => w(i:i+sub_nw-1)
    ! Working array holding H output from H12
    i = i + sub_nw
    h12_h(1:mc) => w(i:i+mc-1)
    ! Transformed matrix A (pointed to by IE index in original code)
    i = i + mc
    ptr_a(1:m,1:l) => w(i:i+m*l-1)
    ! Pointer storing transformed f (pointed to by index IF in original code)
    i = i + m*l
    ptr_b(1:m) => w(i:i+m-1)
    ! Transformed matrix G (pointed to by IG index in original code)
    i = i + m
    ptr_g(1:mg,1:l) => w(i:i+mg*l-1)

    ! === LSEI implementation ===

    ! Triangularize C and apply factors to A and G
    do i = 1, mc
        call h12 (1, i, i+1, c(i,:), h12_h(i), c(i+1:mc,:), trans=.true., &
            status=lstatus, msg=msg)
        if (lstatus /= NF_STATUS_OK) goto 100

        call h12 (2, i, i+1, c(i,:), h12_h(i), e, trans=.true., &
            status=lstatus, msg=msg)
        if (lstatus /= NF_STATUS_OK) goto 100

        call h12 (2, i, i+1, c(i,:), h12_h(i), g, trans=.true., &
            status=lstatus, msg=msg)
        if (lstatus /= NF_STATUS_OK) goto 100
    end do

    ! Solve Cx = d and modify f
    do i = 1, mc
        if (abs(c(i,i)) < epsilon(1.0_PREC)) then
            lstatus = NF_STATUS_INVALID_STATE
            lstatus%code_orig = STATUS_LSEI_SINGULAR_CONSTR
            goto 100
        end if

        x(i) = (d(i) - dot(c(i,1:i-1), x(1:i-1))) / c(i,i)
    end do

    ! Default value if we have no inequality constraints
    ! The original code only erases w(mc+1:mc+mg-mc) which might be an error.
    mult_ineq(:) = 0.0

    if (mc == n) then
        ! X is fully determined by the quality constraint system
        ! Cx = d so there is nothing left to do.

        call solve_orig ()
    else
        ! Either no equality constraints, for Cx=d does not fully determine x
        do i = 1, m
            ! Recall that MC <= N!
            ptr_b(i) = f(i) - dot (e(i,1:mc), x(1:mc))
        end do

        ! Store transformed A
        ptr_a(:,:) = e(:,mc+1:mc+l)

        if (mg == 0) then
            ! No inequality constraints
            tau = sqrt(epsilon(1.0_PREC))
            call hfti (ptr_a, ptr_b, tau, krank, xnorm, sub_w, iwork, &
                status=lstatus, msg=msg)
            if (lstatus /= NF_STATUS_OK) goto 100

            x(mc+1:mc+l) = ptr_b

            if (krank /= l) then
                lstatus = NF_STATUS_INVALID_STATE
                lstatus%code_orig = STATUS_HFTI_RANK_DEFECT
                if (present(msg)) msg = 'LSEI: Matrix rank defect in HFTI'
                goto 100
            end if

            call solve_orig ()
        else
            ! Inequality constraints present.
            ! Modify H and solve the inequality-constrained LS problem.

            ! Implementation note:
            ! LSI/LDP/NNLS store dual vector in the first MG+1 elements of
            ! the working array, which is passed directly from NNLS
            ! through LDP.
            ! These can be read from MULT_INEQ which points at the same
            ! memory block.
            do i = 1, mg
                h(i) = h(i) - dot (g(i,1:mc),x(1:mc))
            end do
            ptr_g(:,:) = g(:,mc+1:mc+l)
            call lsi (ptr_a, ptr_b, ptr_g, h, x(mc+1:n), xnorm, sub_w, &
                iwork, status=lstatus, msg=msg)
            if (lstatus /= NF_STATUS_OK) goto 100

            t = NRM2 (mc, x(1:mc), 1)
            xnorm = sqrt(xnorm**2.0_PREC + t**2.0_PREC)

            if (mc > 0) then
                call solve_orig ()
            end if
        end if
    end if

100 continue

    ! Write back status code
    if (present(status)) status = lstatus

    contains

    subroutine solve_orig ()
        !*  SOLVE_ORIG computes the solution of the original problem as well
        !   as the Lagrange multipliers.
        integer :: i
        real (PREC) :: tmp

        ! TODO: replace with GEMV call
        do i = 1, m
            f(i) = dot (e(i,:), x) - f(i)
        end do

        do i = 1, mc
            tmp = dot_product (e(1:m,i), f(1:m))
            tmp = tmp - dot_product (g(1:mg,i), mult_ineq)
            d(i) = tmp
        end do

        ! Solution vector X
        do i = mc, 1, -1
            call h12 (2, i, i+1, c(i,:), h12_h(i), x, status=lstatus, msg=msg)
            if (lstatus /= NF_STATUS_OK) return
        end do

        ! --- Lagrange multipliers for equality constraints ---

        ! Solve for multipliers, write them to first MC positions of
        ! working array, which is pointed to by MULT_EQ.
        ! Note that we do not use any values in MULT_EQ => W(1:MC)
        ! as inputs!
        !
        ! Recall that we are solving a triangular system backwards.
        ! The last multiplier will be just d(mc)/c(mc,mc),
        ! and recall that MC <= N, so these are valid indices on C.
        do i = mc, 1, -1
            tmp = dot_product (c(i+1:mc,i), mult_eq(i+1:mc))
            mult_eq(i) = (d(i) - tmp) / c(i,i)
        end do

    end subroutine

end subroutine



subroutine lsi_query (m, n, mg, w, iwork)
    integer, intent(in) :: m, n
        !*  Dimensions of matrix A in LS problem Ax = b
    integer, intent(in) :: mg
        !*  First dimension of constraint matrix G in Gx >= h
    integer, intent(out) :: w, iwork

    ! Only workspace requirements are those by LDP()
    call ldp_query (mg, n, w, iwork)

end subroutine



subroutine lsi (e, f, g, h, x, xnorm, w, iwork, status, msg)
    !*  LSI solves the least squares problem with inequality constraints
    !   given by
    !       min_x  ||E*x - f||   s.t. G*x >= h
    !
    !   The algorithm is based on a QR decomposition as described in
    !   section 23.5 of
    !       Lawson, Hanson (1995): Solving least squares problems.
    !
    !   Original implementation in Fortran 77 by Dieter Kraft.
    !
    real (PREC), intent(inout), dimension(:,:), contiguous :: e
        !*  Matrix E with dimension [M,N].
    real (PREC), intent(inout), dimension(:), contiguous :: f
        !*  Vector f of length M.
    real (PREC), intent(inout), dimension(:,:), contiguous :: g
        !*  Matrix G with dimensions [MG,N].
    real (PREC), intent(inout), dimension(:), contiguous :: h
        !*  Vector h of length MG.
    real (PREC), intent(out), dimension(:), contiguous :: x
        !*  Solution vector x of length N.
    real (PREC), intent(out) :: xnorm
        !*  Euclidean norm of the solution.
    real (PREC), intent(out), dimension(:), contiguous :: w
        !*  Working array. On successful termination, the first MG
        !   elements contain the Lagrange multipliers of the inequality
        !   constraints.
    integer, intent(out), dimension(:), contiguous :: iwork
        !*  Integer working array.
    type (status_t), intent(out), optional :: status
        !*  Status code. Returns the NUMFORT status codes,
        !   as well as the following numerical codes that correspond
        !   to those in the original implementation and are stored
        !   in CODE_ORIG of the STATUS_T object.
        !       1.  NF_STATUS_OK on successful termination;
        !       2.  NF_STATUS_INVALID_ARG if any of the input arrays has
        !           an invalid dimensions.
        !       3.  NF_STATUS_MAX_ITER if the max. iteration count in NNLS
        !           was exceeded.
        !       4.  NF_STATUS_INVALID_STATE if the inequality constraints
        !           are incompatible.
        !       5.  NF_STATUS_INVALID_STATE if the matrix E is not of full
        !           rank.
    character (*), intent(out), optional :: msg

    type (status_t) :: lstatus
    integer :: m, n, mg, i, j, k, nw, niwork
    real (PREC) :: d, t

    lstatus = NF_STATUS_OK
    ! Status code used in original SLSQP
    lstatus%code_orig = STATUS_OK

    m = size(e, 1)
    n = size(e, 2)
    mg = size(g, 1)

    ! === Workspace requirements ===

    call ldp_query (mg, n, nw, niwork)

    ! === Input checks ===

    call check_cond (size(f) == m, 'LSI', 'Non-conformable arrays E, F', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(h) == mg, 'LSI', 'Non-conformable arrays G, H', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(g, 2) == n, 'LSI', 'Non-conformable arrays A, G', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(x) == n, 'LSI', 'Non-conformable arrays A, X', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(w) >= nw, 'LSI', 'Working array W too small', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(iwork) >= niwork, 'LSI', 'Working array IWORK too small', lstatus, msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! === Algorithm for LSI ===

    ! QR factors of E, apply to f
    do j = 1, n
        call h12 (1, j, j+1, e(:,j), t, e(:,j+1:n))
        call h12 (2, j, j+1, e(:,j), t, f, status=lstatus, msg=msg)
        if (lstatus /= NF_STATUS_OK) goto 100
    end do

    ! Check that transformed E has full rank
    do j = 1, n
        if (abs(e(j,j)) < epsilon(1.0_PREC)) then
            lstatus = NF_STATUS_INVALID_STATE
            lstatus%code_orig = STATUS_LSI_SINGULAR
            if (present(msg)) msg = 'LSI: Matrix E does not have full rank!'
            goto 100
        end if
    end do

    ! Transform G and h to get least distance problem.
    ! The original code was along the lines of
    !    do j = 1, n
    !        do i = 1, mg
    !            d = dot (g(i,1:j-1), e(1:j-1,j))
    !            g(i,j) = (g(i,j) - d) / e(j,j)
    !        end do
    !    end do
    ! which creates non-contiguous memory access in DOT. Adding an extra
    ! loop speeds things up.
    do j = 1, n
        do k = 1, j - 1
            t = e(k,j)
            do i = 1, mg
                g(i,j) = g(i,j) - g(i,k) * t
            end do
        end do
        g(:,j) = g(:,j) / e(j,j)
    end do

    ! Compute h = - G*f + h
    ! using GEMV with alpha = -1.0, beta = 1.0
!    call GEMV ('N', mg, n, - 1.0_PREC, g, mg, f, 1, 1.0_PREC, h, 1)
    do i = 1, mg
        d = dot (g(i,:), f)
        h(i) = h(i) - d
    end do

    ! Solve least distable problem
    call ldp (g, h, x, xnorm, w, iwork, status=lstatus, msg=msg)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Recover solution of original problem
    x(1:n) = x(1:n) + f(1:n)
    do i = n, 1, -1
        d = dot (e(i,i+1:n),x(i+1:n))
        x(i) = (x(i) - d)/e(i,i)
    end do

    j = min(n+1, m)
    t = NRM2 (m-n, f(j:j+m-n-1), 1)
    xnorm = sqrt(xnorm**2.0 + t**2.0_PREC)

100 continue

    ! Write back status code
    if (present(status)) status = lstatus

end subroutine



function dot (x, y) result(res)
    !*  DOT implements the dot product for assumed-shape arguments
    !   that are not required to be CONTIGUOUS.
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: y
    real (PREC) :: res

    integer :: i, n

    n = size(x)

    res = 0.0_PREC
    do i = 1, n
        res = res + x(i) * y(i)
    end do

end function



subroutine ldl (n, a, z, sigma, w, status, msg)
    integer, intent(in) :: n
    real (PREC), intent(inout), dimension(:), contiguous :: a
    real (PREC), intent(inout), dimension(:), contiguous :: z
    real (PREC), intent(in) :: sigma
    real (PREC), intent(out), dimension(:), contiguous :: w
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    type (status_t) :: lstatus
    integer :: i, j, ij
    real (PREC) :: v, u, t, tp, beta, alpha, delta, gamma

    lstatus = NF_STATUS_OK

    ! === Input checks ===

    call assert_input (size(a) == n*(n+1)/2, 'Non-conformable array A')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input (size(z) == n, 'Non-conformable array Z')
    if (lstatus /= NF_STATUS_OK) goto 100

    call assert_input ((sigma >= 0.0_PREC) .or. (size(w) == n), &
        'Non-conformable array W')
    if (lstatus /= NF_STATUS_OK) goto 100

    ! === LDL implementation ===

    if (sigma == 0.0_PREC) goto 100

    ij = 1
    t = 1.0_PREC / sigma

    if (sigma < 0) then
        ! Prepare negative update
        w(1:n) = z(1:n)

        ! Solve triangular system Lv = z.
        ! Recall that L is lower triangular with unit diagonal elements,
        ! so we have that
        !   v(1) = z(1)
        !   v(2) = z(2) - v(1) * L(2,1)
        !   v(3) = z(3) - v(1) * L(3,1) - v(2) * L(3,2)
        ! etc.
        do i = 1, n
            v = w(i)
            t = t + v*v/a(ij)
            ! Subtract v(i) * L(j+1,i) from z(j+1), etc. for remaining j > i
            do j = i + 1, n
                ij = ij + 1
                w(j) = w(j) - v*a(ij)
            end do

            ij = ij + 1
        end do

        if (t >= 0.0_PREC) then
            t = epsilon(1.0_PREC) / sigma
        end if

        do i = 1, n
            j = n + 1 - i
            ij = ij - i
            u = w(j)
            w(j) = t
            t = t - u*u / a(ij)
        end do
    end if

    ! Update matrix A
    do i = 1, n
        v = z(i)
        delta = v / a(ij)
        if (sigma < 0.0_PREC) then
            tp = w(i)
        else
            ! This covers the case sigma > 0, since sigma = 0 exits immediately
            tp = t + delta*v
        end if

        alpha = tp / t
        a(ij) = alpha * a(ij)

        if (i == n) goto 100

        beta = delta / tp
        if (alpha <= 4.0_PREC) then
            do j = i + 1, n
                ij = ij + 1
                z(j) = z(j) - v*a(ij)
                a(ij) = a(ij) + beta * z(j)
            end do
        else
            gamma = t / tp
            do j = i + 1, n
                ij = ij + 1
                u = a(ij)
                a(ij) = gamma * u + beta*z(j)
                z(j) = z(j) - v*u
            end do
        end if

        ij = ij + 1
        t = tp
    end do


100 continue

    if (present(status)) status = lstatus

    contains

    subroutine assert_input (cond, lmsg)
        logical, intent(in) :: cond
        character (*), intent(in) :: lmsg

        lstatus = NF_STATUS_OK
        if (.not. cond) then
            lstatus = NF_STATUS_INVALID_ARG
            if (present(msg)) then
                msg = 'LDL: Invalid input: ' // lmsg
            end if
        end if
    end subroutine

end subroutine
