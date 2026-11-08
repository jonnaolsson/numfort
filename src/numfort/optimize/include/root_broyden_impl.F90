


subroutine root_broyden_check_input (maxiter, maxfev, &
        tol, xtol, rstep, xstep, status, msg)
    integer, intent(in):: maxiter, maxfev
    real (PREC), intent(in) :: tol, xtol
    real (PREC), intent(in), optional :: rstep, xstep
    type (status_t), intent(inout) :: status
    character (*), intent(inout) :: msg

    status = NF_STATUS_OK

    call check_nonneg (0, maxiter, "maxiter", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_nonneg (0, maxfev, "maxfev", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, tol, "tol", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, xtol, "xtol", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, rstep, "rstep", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, xstep, "xstep", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    return

100 continue
    status = NF_STATUS_INVALID_ARG
end subroutine



recursive subroutine root_broyden (fcn, x, ndiff, tol, xtol, &
        maxiter, maxfev, rstep, xstep, dstep, iprint, work, jac_inv_init, res)
    procedure (fvv_fcn) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    logical, intent(in) :: ndiff
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    real (PREC), intent(in), optional :: dstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
        !*  If present, debug info will be printed according to the value
        !   of IPRINT.
    type (workspace), intent(inout), optional :: work
    real (PREC), intent(in), dimension(:,:), contiguous, optional :: jac_inv_init
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    call wrap_procedure (fwrapper, fcn=fcn, eps=dstep)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, jac_inv_init, res)

end subroutine



recursive subroutine root_broyden_jac (fcn, fjac, x, tol, &
        xtol, maxiter, maxfev, rstep, xstep, iprint, work, jac_inv_init, res)
    procedure (fvv_fcn) :: fcn
    procedure (fvv_jac) :: fjac
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    !*  If present, debug info will be printed according to the value
    !   of IPRINT.
    type (workspace), intent(inout), optional :: work
    real (PREC), intent(in), dimension(:,:), contiguous, optional :: jac_inv_init
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper

    call wrap_procedure (fwrapper, fcn=fcn, jac=fjac)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, jac_inv_init, res)

end subroutine



recursive subroutine root_broyden_fcn_jac_opt (fcn, x, tol, &
        xtol, maxiter, maxfev, rstep, xstep, iprint, work, jac_inv_init, res)
    procedure (fvv_fcn_jac_opt) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    !*  If present, debug info will be printed according to the value
    !   of IPRINT.
    type (workspace), intent(inout), optional :: work
    real (PREC), intent(in), dimension(:,:), contiguous, optional :: jac_inv_init
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac_opt=fcn)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, jac_inv_init, res)

end subroutine



recursive subroutine root_broyden_args (fcn, x, args, ndiff, &
        tol, xtol, maxiter, maxfev, rstep, xstep, dstep, iprint, work, jac_inv_init, res)
    procedure (fvv_fcn_args) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    logical, intent(in) :: ndiff
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    real (PREC), intent(in), optional :: dstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    !*  If present, debug info will be printed according to the value
    !   of IPRINT.
    type (workspace), intent(inout), optional :: work
    real (PREC), intent(in), dimension(:,:), contiguous, optional :: jac_inv_init
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    call wrap_procedure (fwrapper, fcn_args=fcn, args=args, eps=dstep)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, jac_inv_init, res)

end subroutine



recursive subroutine root_broyden_jac_args (fcn, fjac, x, &
        args, tol, xtol, maxiter, maxfev, rstep, xstep, iprint, jac_inv_init, work, res)
    procedure (fvv_fcn_args) :: fcn
    procedure (fvv_jac_args) :: fjac
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    !*  If present, debug info will be printed according to the value
    !   of IPRINT.
    type (workspace), intent(inout), optional :: work
    real (PREC), intent(in), dimension(:,:), contiguous, optional :: jac_inv_init
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper

    call wrap_procedure (fwrapper, fcn_args=fcn, jac_args=fjac, args=args)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, jac_inv_init, res)

end subroutine



recursive subroutine root_broyden_fcn_jac_opt_args (fcn, x, &
        args, tol, xtol, maxiter, maxfev, rstep, xstep, iprint, work, jac_inv_init, res)
    procedure (fvv_fcn_jac_opt_args) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
    !*  Max. number of function evaluations (includes evaluations of
    !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    !*  If present, debug info will be printed according to the value
    !   of IPRINT.
    type (workspace), intent(inout), optional :: work
    real (PREC), intent(in), dimension(:,:), contiguous, optional :: jac_inv_init
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac_opt_args=fcn, args=args)

    call root_broyden_impl (fwrapper, x, tol, xtol, maxiter, maxfev, &
        rstep, xstep, iprint, work, jac_inv_init, res)

end subroutine



recursive subroutine root_broyden_impl (fcn, x, tol, xtol, &
        maxiter, maxfev, rstep, xstep, iprint, work, jac_inv_init, res)

    type (fwrapper_vv), intent(inout) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), optional :: tol
    real (PREC), intent(in), optional :: xtol
    integer, intent(in), optional :: maxiter
    integer, intent(in), optional :: maxfev
        !*  Max. number of function evaluations (includes evaluations of
        !   Jacobian obtained by numerical differentiation, if applicable)
    real (PREC), intent(in), optional :: rstep
        !*  Max. step size in search direction, relative to the current
        !   point (default: unbounded)
    real (PREC), intent(in), optional :: xstep
        !*  Max. absolute step size in search direction (default: unbounded)
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    type (workspace), intent(inout), optional, target :: work
    real (PREC), intent(in), dimension(:,:), contiguous, optional :: jac_inv_init
    type (optim_result), intent(inout), optional, target :: res

    real (PREC) :: ltol, lxtol
    real (PREC) :: dx_scale, denom, denom_scale, nrm, nrmp1, nrm_last, nrm_upd
    integer :: lmaxiter, lmaxfev, req_fev, k, nit, n, i, nrwrk, niwrk, liprint
    integer :: lwork_inv, liwork_inv
        !   Workspace array sizes use for INV routine
    integer, dimension(2) :: shp2d
    real (PREC), dimension(:), pointer, contiguous :: fx, fxlast, dx, dfx
    real (PREC), dimension(:), pointer, contiguous :: vec1, vec2
        !   Temporary arrays used for various BLAS routines
    real (PREC), dimension(:), pointer, contiguous :: x_ls, fx_ls
        !   Working arrays used for line search
    real (PREC), dimension(:,:), pointer, contiguous :: jac, jac_inv
    real (PREC), dimension(:), pointer, contiguous :: rwork_inv
    integer, dimension(:), pointer, contiguous :: iwork_inv
    logical :: has_jac_inv, jac_inv_valid

    type (workspace), pointer :: ptr_work
    type (optim_result), pointer :: ptr_res
    type (status_t) :: status

    ! Arguments for GEMV
    character (1) :: trans
    real (PREC) :: alpha, beta
    integer, parameter :: incx = 1, incy = 1
    integer :: lda, m
    
    status = NF_STATUS_OK
    ! Outer iteration counter (nit) reported in result
    nit = 0
    ! Tracks whether jac_inv holds a valid caller-supplied or computed inverse
    jac_inv_valid = .false.

    nullify (ptr_work, ptr_res)
    nullify (fx, fxlast, dx, dfx)
    nullify (vec1, vec2)
    nullify (x_ls, fx_ls)
    nullify (jac, jac_inv)
    nullify (rwork_inv, iwork_inv)

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    ! Default arguments
    lmaxiter = 100
    ltol = sqrt(epsilon(0.0_PREC))
    lxtol = sqrt(epsilon(0.0_PREC))
    liprint = PRINT_NONE

    ! Overwrite defaults with optionally provided user arguments
    ! Note: RSTEP and XSTEP have no default values, we skip limiting the
    ! step size completly whenever they are not present.
    if (present(maxiter)) lmaxiter = maxiter
    if (present(tol)) ltol = tol
    if (present(xtol)) lxtol = xtol
    if (present(iprint)) liprint = iprint

    ! Max. number of function evaluations given by max. iteration count and
    ! max. number of backtracking steps performed during linesearch.
    ! Add 1 for initial function evaluation, 1 to ensure that root finder
    ! will not exit due to NFEV exceeding MAXFUN when it in fact should
    ! terminate due to exceeding MAXITER.
    lmaxfev = lmaxiter * LINESEARCH_MAX_STEPS + 2
    if (fcn%num_diff) then
        ! If numerical differentiation is performed
        lmaxfev = lmaxfev + size(x)
    end if
    if (present(maxfev)) lmaxfev = maxfev

    ! Validate inputs
    call root_broyden_check_input (lmaxiter, lmaxfev, ltol, lxtol, &
        rstep, xstep, status, ptr_res%msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    if (lmaxfev == 0) then
        ptr_res%msg = "Max. number of function evaluations exceeded"
        status = NF_STATUS_MAX_EVAL
        goto 100
    end if

    n = size(x)

    ! Arguments to GEMV that never change
    m = n
    lda = n

    ! Workspace array size
    nrwrk = 2*n*n + 8*n
    ! Working arrays needed for INV
    call inv_work_query (0.0_PREC, n, lwork_inv, liwork_inv)
    nrwrk = nrwrk + lwork_inv
    niwrk = liwork_inv

    ! Allocate workspace arrays
    call assert_alloc_ptr (work, ptr_work)
    ! Clear any internal state in workspace object, in particular index offsets
    ! (this does not deallocate working arrays)
    call workspace_reset (ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk, niwrk=niwrk)

    call workspace_get_ptr (ptr_work, n, vec1, status)
    if (status /= NF_STATUS_OK) goto 100

    call workspace_get_ptr (ptr_work, n, vec2, status)
    if (status /= NF_STATUS_OK) goto 100

    call workspace_get_ptr (ptr_work, n, fx, status)
    if (status /= NF_STATUS_OK) goto 100

    call workspace_get_ptr (ptr_work, n, fxlast, status)
    if (status /= NF_STATUS_OK) goto 100

    call workspace_get_ptr (ptr_work, n, dx, status)
    if (status /= NF_STATUS_OK) goto 100

    call workspace_get_ptr (ptr_work, n, dfx, status)
    if (status /= NF_STATUS_OK) goto 100

    call workspace_get_ptr (ptr_work, n, fx_ls, status)
    if (status /= NF_STATUS_OK) goto 100

    call workspace_get_ptr (ptr_work, n, x_ls, status)
    if (status /= NF_STATUS_OK) goto 100

    shp2d = n
    call workspace_get_ptr (ptr_work, shp2d, jac, status)
    if (status /= NF_STATUS_OK) goto 100

    ! Working arrays for INV
    call workspace_get_ptr (ptr_work, shp2d, jac_inv, status)
    if (status /= NF_STATUS_OK) goto 100

    call workspace_get_ptr (ptr_work, n, rwork_inv, status)
    if (status /= NF_STATUS_OK) goto 100

    call workspace_get_ptr (ptr_work, n, iwork_inv, status)
    if (status /= NF_STATUS_OK) goto 100


    if (liprint /= PRINT_NONE) then
        print '(tr1, a)', "ROOT_BROYDEN: routine start"
    end if

    has_jac_inv = present (jac_inv_init)
    if (has_jac_inv) then
        ! TODO: check whether initial step size implied by Jacobian is smaller than
        ! xtol, otherwise no point in even trying.
        has_jac_inv = all (ieee_is_finite (jac_inv_init)) &
            .and. any (abs (jac_inv_init) > 1.0e-3_PREC)
    end if

    if (has_jac_inv) then
        jac_inv(:,:) = jac_inv_init
        jac_inv_valid = .true.
        if (iand(liprint, PRINT_JAC) == PRINT_JAC) then
            print '(tr1, a)', "ROOT_BROYDEN: Using caller-provided initial inverse Jacobian"
            do i = 1, n
                print 200, "", jac_inv(i,:)
            end do
        end if
    end if

    call dispatch (fcn, x, fxlast)

    if (.not. all(ieee_is_finite(fxlast))) then
        ptr_res%msg = "Invalid function value encountered"
        status = NF_STATUS_INVALID_STATE
        jac_inv_valid = .false.
        goto 100
    end if

    ! Check whether initial point satisfied convergence criterion
    nrm_last = norm2 (fxlast)

    if (liprint /= PRINT_NONE) then
        print '(tr1, a)', "ROOT_BROYDEN: initial values"
        print 200, "x: ", x
        print 200, "f(x): ", fxlast
        print 201, "2-norm: ", nrm_last
    end if

    if (nrm_last < ltol) then
        ptr_res%msg = "Convergence achieved, func. value smaller than tol"
        status = NF_STATUS_OK
        goto 100
    end if

    if (.not. has_jac_inv) then
        if (fcn%num_diff) then
            req_fev = n
        else
            req_fev = 1
        end if

        if (fcn%nfev + req_fev > lmaxfev) then
            ptr_res%msg = "Max. number of function evaluations exceeded"
            status = NF_STATUS_MAX_EVAL
            goto 100
        end if

        ! Compute Jacobian at initial point. Reuse f(X) evaluated above in
        ! case of numerical differentiation.
        call dispatch_jac (fcn, x, jac, fxlast)

        if (.not. all(ieee_is_finite(jac))) then
            ptr_res%msg = "Invalid value in initial Jacobian encountered"
            status = NF_STATUS_INVALID_STATE
            goto 100
        end if

        if (iand(liprint, PRINT_JAC) == PRINT_JAC) then
            print '(tr1, a)', "ROOT_BROYDEN: initial Jacobian"
            do i = 1, n
                print 200, "", jac(i,:)
            end do
        end if

        ! Initial inverted Jacobian
        call inv (jac, jac_inv, rwork=rwork_inv, iwork=iwork_inv, status=status)
        if (status /= NF_STATUS_OK) then
            ptr_res%msg = "Could not invert initial Jacobian"
            goto 100
        end if

        jac_inv_valid = .true.
    end if

   do k = 1, lmaxiter

        ! 1. Find next candidate point
        alpha = -1.0_PREC
        beta = 0.0_PREC
        trans = 'N'
        ! Compute dx_k = - J^{-1}_{k-1} * f(x_{k-1})
        call GEMV (trans, m, n, alpha, jac_inv, lda, fxlast, incx, beta, dx, incy)

        dx_scale = 1.0_PREC
        if (present(rstep)) then
            ! Perform "dampened" update step, allowing point-wise steps of
            ! at most rstep * abs(x_k).
            ! Standard Broyden would set x_k = x_{k-1} + dx_k
            do i = 1, n
                dx_scale = min(dx_scale, rstep * abs(x(i)) / abs(dx(i)))
            end do
        end if

        if (present(xstep)) then
            ! Check that step size does not violate max. absolute step
            ! size for any dimension.
            do i = 1, n
                dx_scale = min(dx_scale, xstep / abs(dx(i)))
            end do
        end if

        if (iand(liprint, PRINT_STEP) == PRINT_STEP) then
            print '(tr1,a)', "ROOT_BROYDEN: STEP direction"
            print 200, "Unrestricted step direction: ", dx
            print 201, "Rescaling factor: ", dx_scale
        end if

        ! Rescale step size such that change in x is at most EPSX percent
        ! for each element.
        dx(:) = dx_scale * dx

        call dumb_line_search (fcn, x, nrm_last, liprint, x_ls, fx_ls, dx, fx, &
            lmaxfev, status)

        if (status == NF_STATUS_MAX_EVAL) then
            ptr_res%msg = "Max. number of function evaluations exceeded"
            goto 100
        else if (status /= NF_STATUS_OK) then
            ptr_res%msg = "Invalid function value encountered during line search"
            goto 100
        end if

        ! Record completed iteration count once line search finds an admissible step
        nit = k

         ! 2. Check convergence in terms of function values
        nrm_upd = norm2(fx)
        if (nrm_upd < ltol) then
            ptr_res%msg = 'Convergence achieved, func. value smaller than tol'
            status = NF_STATUS_OK
            ! Code below expects the final function value to be stored in FXLAST
            call COPY (n, fx, 1, fxlast, 1)
            x(:) = x + dx

            goto 100
        end if

        ! 3. Check if x_k and x_{k-1} are sufficiently close to each other to
        ! terminate algorithm
        nrm = norm2(dx)
        nrmp1 = norm2(x) + 1.0_PREC
        if (nrm < lxtol * nrmp1) then
            ptr_res%msg = 'Not making sufficient progress, change in x smaller than xtol'
            status = NF_STATUS_NOT_CONVERGED

            ! Determine whether to return (X,FXLAST) or (X+DX,FX) is a
            ! better approximation
            if (nrm_last > nrm_upd) then
                x(:) = x + dx
                call COPY (n, fx, 1, fxlast, 1)
            end if

            goto 100
        end if

        ! 4. Check whether max. number of func. evaluations was exceeded.
        if (fcn%nfev >= lmaxfev) then
            ! Set corresponding exit status.
            ! At this point the last "best" guess for the root is stored in X,
            ! the corresponding function value in FXLAST.
            ptr_res%msg = 'Max. number of function evaluations exceeded'
            status = NF_STATUS_MAX_EVAL

            ! Determine whether to return (X,FXLAST) or (X+DX,FX) is a
            ! better approximation
            if (nrm_last > nrm_upd) then
                x(:) = x + dx
                call COPY (n, fx, 1, fxlast, 1)
            end if

            goto 100
        end if

        ! Update X with the "best" DX identified by line search
        x(:) = x + dx

        ! Compute difference to last objective (fx contains the best linesearch
        ! result)
        ! Use helper routine as otherwise compiler creates temp. arrays
        call diff1d (fx, fxlast, dfx)

        ! Keep FXLAST and NRMLAST synchronized with updated X
        call COPY (n, fx, 1, fxlast, 1)
        nrm_last = nrm_upd

        ! 4. Update inverse of J_k
        ! Compute J^{-1}_{k-1} dfx_k
        alpha = 1.0_PREC
        beta = 0.0_PREC
        trans = 'N'
        call GEMV (trans, m, n, alpha, jac_inv, lda, dfx, incx, beta, vec1, incy)
        ! Denominator in expression for J^{-1}_k, dx_k^T J^{-1}_{k-1} dfx_k
        denom = DOT (n, vec1, incx, dx, incy)

        ! Cauchy-Schwarz upper bound ||dx||_2 * ||vec1||_2 used as relative scale
        ! to detect near-orthogonal or cancelled update directions
        denom_scale = norm2 (dx) * norm2 (vec1)
        if (.not. ieee_is_finite (denom) .or. abs (denom) <= tiny (1.0_PREC) .or. &
            abs (denom) <= epsilon (1.0_PREC) * denom_scale) then
            ptr_res%msg = "Singular or unsafe Broyden update denominator"
            status = NF_STATUS_INVALID_STATE
            goto 100
        end if

        ! Compute transpose of dx_k^T * J^{-1}_{k-1}
        trans = 'T'
        alpha = 1.0_PREC
        beta = 0.0_PREC
        call GEMV (trans, m, n, alpha, jac_inv, lda, dx, incx, beta, vec2, incy)

        ! Compute X argument in call to DGER, ie. dx_k - J^{-1}_{k-1} dfx_k
        ! Set VEC1 = - J^{-1}_{k-1} dfx_k
        alpha = -1.0_PREC
        call SCAL (n, alpha, vec1, incx)
        ! Add dx_k to VEC1
        alpha = 1.0_PREC
        call AXPY (n, alpha, dx, incx, vec1, incy)
        ! Compute update to inverse Jacobian using BLAS rank-1 update routine
        ! which adds alpha * VEC1 * VEC2^T to J^{-1}_{k-1}
        alpha = 1.0_PREC / denom
        call GER (m, n, alpha, vec1, incx, vec2, incy, jac_inv, lda)
    end do

    ptr_res%msg = 'Max. number of iterations exceeded'
    status = NF_STATUS_MAX_ITER
    ! Correct number of iterations
    nit = lmaxiter

100 continue

    if (present(res)) then
        if (jac_inv_valid) then
            call result_update (ptr_res, x, fxlast, jac_inv=jac_inv, &
                status=status, nit=nit, nfev=fcn%nfev)
        else
            call result_update (ptr_res, x, fxlast, &
                status=status, nit=nit, nfev=fcn%nfev)
        end if
    end if

    ! Clean up local WORKSPACE object if none was passed by client code
    call assert_dealloc_ptr (work, ptr_work)

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

    if (liprint /= PRINT_NONE) then
        print '(tr1, a)', "ROOT_BROYDEN: routine complete"
    end if

200 format (tr3, a, *(es12.5e2,:,", "))
201 format (tr3, a, es12.5e2)

end subroutine



recursive subroutine dumb_line_search (fcn, x, nrm, iprint, &
        x_ls, fx_ls, dx, fx, maxfev, status)
    !*  DUMB_LINE_SEARCH uses backtracking to identify the step size in a
    !   given direction that yields a lower Euclidean norm of the objective
    !   function than at a given point X.
    !
    !   If the Euclidean norm is larger at all candidate points in direction
    !   DX, the least bad step size is returned.
    type (fwrapper_vv), intent(inout) :: fcn
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in) :: nrm
        !*  Euclidean norm of the last function value, |FX|_2
    integer (NF_ENUM_KIND), intent(in) :: iprint
    real (PREC), intent(inout), dimension(:), contiguous :: x_ls, fx_ls
        !*  Working arrays to store X, FX evaluated during line search
    real (PREC), intent(inout), dimension(:) :: dx
        !*  On entry, contains the max. step size in a desired direction.
        !   On exit, the value is updated to contain the "optimal" step size.
    real (PREC), intent(inout), dimension(:), contiguous :: fx
        !*  On entry, contains the function value evaluated at X. On exit,
        !   contains the function value at X + DX where DX is the updated
        !   value computed by the routine.
    integer, intent(in) :: maxfev
        !*  Max. number of function evaluations allowed
    type (status_t), intent(out) :: status
        !*  Status code: NF_STATUS_OK if an admissible finite trial was
        !   identified, NF_STATUS_MAX_EVAL if budget was exhausted before any
        !   admissible trial was evaluated, or NF_STATUS_INVALID_STATE if all
        !   evaluated trials were invalid.

    integer :: i, n
    real (PREC) :: step_size, nrm_ls, nrm_ls_best, step_best
    logical :: has_admissible, do_print

    n = size(x)

    do_print = iand(iprint, PRINT_LSEARCH) == PRINT_LSEARCH

    ! Line search
    has_admissible = .false.
    nrm_ls_best = huge(1.0_PREC)
    step_best = 0.0_PREC

    if (do_print) then
        print '(tr1, a)', "ROOT_BROYDEN: LINE SEARCH routine"
    end if

    do i = 1, LINESEARCH_MAX_STEPS
        if (fcn%nfev >= maxfev) exit

        step_size = (LINESEARCH_MAX_STEPS - i + 1.0_PREC) / LINESEARCH_MAX_STEPS
        x_ls(:) = x + step_size * dx

        call dispatch (fcn, x_ls, fx_ls)

        if (do_print) then
            if (i > 1) then
                print '(tr3, a, i0)', "ROOT_BROYDEN: BACKTRACKING step #", i - 1
            end if
            print 200, "x: ", x_ls
            print 200, "f(x): ", fx_ls
        end if

        if (.not. all(ieee_is_finite(fx_ls))) cycle

        nrm_ls = norm2(fx_ls)
        if (.not. ieee_is_finite(nrm_ls)) cycle

        if (do_print) then
            print 201, "2-norm: ", nrm_ls
        end if

        has_admissible = .true.

        if (nrm_ls <= nrm) then
            fx(:) = fx_ls
            dx = dx * step_size
            status = NF_STATUS_OK
            goto 100
        else if (nrm_ls < nrm_ls_best) then
            fx(:) = fx_ls
            step_best = step_size
            nrm_ls_best = nrm_ls
        end if
    end do

    if (has_admissible) then
        ! No step size with decreasing objective found, take the least bad
        dx = step_best * dx
        status = NF_STATUS_OK
    else if (fcn%nfev >= maxfev) then
        status = NF_STATUS_MAX_EVAL
    else
        status = NF_STATUS_INVALID_STATE
    end if

100 continue

    if (do_print) then
        print '(tr1, a)', "ROOT_BROYDEN: LINE SEARCH complete"
    end if

200 format (tr3, a, *(es12.5e2, :, ", "))
201 format (tr3, a, es12.5e2)

end subroutine



subroutine diff1d (x1, x2, dx)
    !*  Helper routine to avoid creating temporary arrays
    real (PREC), intent(in), dimension(:), contiguous :: x1, x2
    real (PREC), intent(out), dimension(:), contiguous :: dx

    dx = x1 - x2

end subroutine
