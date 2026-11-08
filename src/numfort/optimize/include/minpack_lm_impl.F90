

subroutine root_lm (fcn, x, fx, ndiff, ftol, xtol, gtol, maxfev, &
        factor, diag, dstep, drstep, work, res)

    procedure (fvv_fcn) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx
    logical, intent(in) :: ndiff
    real (PREC), intent(in), optional :: ftol
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: gtol
    integer, intent(in), optional :: maxfev
    real (PREC), intent(in), optional :: factor
    real (PREC), intent(in), dimension(:), optional, contiguous :: diag
    real (PREC), intent(in), optional :: dstep
    real (PREC), intent(in), optional :: drstep
        !*  Relative step size for numerical forward differencing. The
        !   step size for the derivative wrt. X(j) is computed as
        !       h = DRSTEP * abs(X(j))
        !   Note: Ignored if DSTEP is present.
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper
    real (PREC) :: reps

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            call result_reset (res)
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    ! MINPACK compatibility: If neither DSTEP nor DRSTEP are present, use
    ! rel. step size of SQRT(machine eps).
    if (.not. present(dstep) .and. .not. present(drstep)) then
        reps = sqrt(epsilon(0.0_PREC))
        call wrap_procedure (fwrapper, fcn=fcn, reps=reps)
    else
        call wrap_procedure (fwrapper, fcn=fcn, eps=dstep, reps=drstep)
    end if
    
    call root_lm_impl (fwrapper, x, fx, ftol, xtol, gtol, maxfev, &
        factor, diag, work, res)

end subroutine



subroutine root_lm_jac (fcn, fjac, x, fx, ftol, xtol, gtol, maxfev, &
        factor, diag, work, res)

    procedure (fvv_fcn) :: fcn
    procedure (fvv_jac) :: fjac
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx
    real (PREC), intent(in), optional :: ftol
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: gtol
    integer, intent(in), optional :: maxfev
    real (PREC), intent(in), optional :: factor
    real (PREC), intent(in), dimension(:), optional, contiguous :: diag
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper

    call wrap_procedure (fwrapper, fcn=fcn, jac=fjac)
    
    call root_lm_impl (fwrapper, x, fx, ftol, xtol, gtol, maxfev, &
        factor, diag, work, res)

end subroutine


subroutine root_lm_fcn_jac (fcn, x, fx, ftol, xtol, gtol, maxfev, &
        factor, diag, work, res)

    procedure (fvv_fcn_jac_opt) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx
    real (PREC), intent(in), optional :: ftol
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: gtol
    integer, intent(in), optional :: maxfev
    real (PREC), intent(in), optional :: factor
    real (PREC), intent(in), dimension(:), optional, contiguous :: diag
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac_opt=fcn)
    
    call root_lm_impl (fwrapper, x, fx, ftol, xtol, gtol, maxfev, &
        factor, diag, work, res)

end subroutine



subroutine root_lm_args (fcn, x, args, fx, ndiff, ftol, xtol, gtol, maxfev, &
        factor, diag, dstep, drstep, work, res)

    procedure (fvv_fcn_args) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:), contiguous :: fx
    logical, intent(in) :: ndiff
    real (PREC), intent(in), optional :: ftol
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: gtol
    integer, intent(in), optional :: maxfev
    real (PREC), intent(in), optional :: factor
    real (PREC), intent(in), dimension(:), optional, contiguous :: diag
    real (PREC), intent(in), optional :: dstep
    real (PREC), intent(in), optional :: drstep
        !*  Relative step size for numerical forward differencing. The
        !   step size for the derivative wrt. X(j) is computed as
        !       h = DRSTEP * abs(X(j))
        !   Note: Ignored if DSTEP is present.
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper
    real (PREC) :: reps

    ! Force NDIFF argument to be TRUE
    if (.not. ndiff) then
        if (present(res)) then
            call result_reset (res)
            res%status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    ! MINPACK compatibility: If neither DSTEP nor DRSTEP are present, use
    ! rel. step size of SQRT(machine eps).
    if (.not. present(dstep) .and. .not. present(drstep)) then
        reps = sqrt(epsilon(0.0_PREC))
        call wrap_procedure (fwrapper, fcn_args=fcn, args=args, reps=reps)
    else
        call wrap_procedure (fwrapper, fcn_args=fcn, args=args, eps=dstep, reps=drstep)
    end if

    call root_lm_impl (fwrapper, x, fx, ftol, xtol, gtol, maxfev, &
        factor, diag, work, res)

end subroutine



subroutine root_lm_jac_args (fcn, fjac, x, args, fx, ftol, xtol, gtol, &
        maxfev, factor, diag, work, res)

    procedure (fvv_fcn_args) :: fcn
    procedure (fvv_jac_args) :: fjac
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:), contiguous :: fx
    real (PREC), intent(in), optional :: ftol
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: gtol
    integer, intent(in), optional :: maxfev
    real (PREC), intent(in), optional :: factor
    real (PREC), intent(in), dimension(:), optional, contiguous :: diag
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper

    call wrap_procedure (fwrapper, fcn_args=fcn, jac_args=fjac, args=args)
    
    call root_lm_impl (fwrapper, x, fx, ftol, xtol, gtol, maxfev, &
        factor, diag, work, res)

end subroutine


subroutine root_lm_fcn_jac_args (fcn, x, args, fx, ftol, xtol, gtol, &
        maxfev, factor, diag, work, res)

    procedure (fvv_fcn_jac_opt_args) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:), contiguous :: fx
    real (PREC), intent(in), optional :: ftol
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: gtol
    integer, intent(in), optional :: maxfev
    real (PREC), intent(in), optional :: factor
    real (PREC), intent(in), dimension(:), optional, contiguous :: diag
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vv) :: fwrapper

    call wrap_procedure (fwrapper, fcn_jac_opt_args=fcn, args=args)
    
    call root_lm_impl (fwrapper, x, fx, ftol, xtol, gtol, maxfev, &
        factor, diag, work, res)

end subroutine





subroutine root_lm_impl (fcn, x, fx, ftol, xtol, gtol, maxfev, &
        factor, diag, work, res)

    type (fwrapper_vv), intent(inout) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx
    real (PREC), intent(in), optional :: ftol
    real (PREC), intent(in), optional :: xtol
    real (PREC), intent(in), optional :: gtol
    integer, intent(in), optional :: maxfev
    real (PREC), intent(in), optional :: factor
    real (PREC), intent(in), dimension(:), optional, target, contiguous :: diag
    type (workspace), intent(inout), target, optional :: work
    type (optim_result), intent(inout), optional :: res

    ! local default values for optional arguments
    real (PREC) :: lxtol, lftol, lgtol, lfactor
    integer :: lmaxfev, lnprint
    integer :: m, n, nrwrk, niwrk, mode, info, nfev, njev
    type (optim_result), pointer :: ptr_res
    type (workspace), pointer :: ptr_work
    ! pointers to various arrays that need to be passed to lmder() that are
    ! segments of memory allocated in workspace
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_fjac
    real (PREC), dimension(:), pointer, contiguous :: ptr_qtf, &
        ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4, ptr_diag
    integer, dimension(:), pointer, contiguous :: ptr_ipvt
    integer, dimension(2) :: shp

    nullify (ptr_work, ptr_res)
    nullify (ptr_fjac, ptr_qtf, ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4, ptr_diag)
    nullify (ptr_ipvt)
    
    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    n = size(x)
    m = size(fx)

    call lm_check_input (m, n, ftol, xtol, gtol, maxfev, factor, &
        status=ptr_res%status, msg=ptr_res%msg)
    if (NF_STATUS_INVALID_ARG .in. ptr_res%status) goto 100

    ! set default values
    lmaxfev = 200 * (n + 1)
    lfactor = 100.0_PREC
    ! use default from scipy
    lxtol = sqrt(epsilon(0.0_PREC))
    lftol = sqrt(epsilon(0.0_PREC))
    lgtol = 0.0_PREC
    lnprint = 0

    mode = 1
    info = 0

    ! override defaults if arguments specified by user
    if (present(xtol)) lxtol = xtol
    if (present(ftol)) lftol = ftol
    if (present(gtol)) lgtol = gtol
    if (present(factor)) lfactor = factor
    if (present(maxfev)) lmaxfev = maxfev

    ! Workspace size requirements:
    !   size(fjac) = m * n
    !   size(wa1) = size(wa2) = size(wa3) = n
    !   size(wa4) = m
    !   size(qtf) = n
    nrwrk = 4*n + m + n*m
    ! Add space to store diag
    if (.not. present(diag)) nrwrk = nrwrk + n
    ! used to store ipvt
    niwrk = n

    call assert_alloc_ptr (work, ptr_work)
    ! Clear any internal state in workspace object, in particular index offsets
    ! (this does not deallocate working arrays)
    call workspace_reset (ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk, niwrk=niwrk)

    ! map array arguments to lmder() onto workspace
    ! First n elements reserved for diag
    shp(1) = m
    shp(2) = n
    call workspace_get_ptr (ptr_work, shp, ptr_fjac)
    call workspace_get_ptr (ptr_work, n, ptr_qtf)
    call workspace_get_ptr (ptr_work, n, ptr_wa1)
    call workspace_get_ptr (ptr_work, n, ptr_wa2)
    call workspace_get_ptr (ptr_work, n, ptr_wa3)
    call workspace_get_ptr (ptr_work, m, ptr_wa4)

    if (present(diag)) then
        ptr_diag => diag
        ! set mode such that use-provided diag is used
        mode = 2
    else
        ! no scaling of diagonal elements; initialize diag to 1.0 even though
        ! this is not needed
        call workspace_get_ptr (ptr_work, n, ptr_diag)
        ptr_diag = 1.0_PREC
    end if

    call workspace_get_ptr (ptr_work, n, ptr_ipvt)

    call minpack_lmder (fcn_wrapper, m, n, x, fx, ptr_fjac, m, &
        lftol, lxtol, lgtol, lmaxfev, &
        ptr_diag, mode, lfactor, lnprint, info, nfev, njev, &
        ptr_ipvt, ptr_qtf, ptr_wa1, ptr_wa2, ptr_wa3, ptr_wa4)

    select case (info)
    case (0)
        ptr_res%status = NF_STATUS_INVALID_ARG
        ptr_res%msg = "Invalid input parameters"
    case (1)
        ptr_res%status = NF_STATUS_OK
        ptr_res%msg = "Convergence in terms of ftol"
    case (2)
        ptr_res%status = NF_STATUS_OK
        ptr_res%msg = "Convergence in terms of xtol"
    case (3)
        ptr_res%status = NF_STATUS_OK
        ptr_res%msg = "Convegence in terms of ftol and xtol"
    case (4)
        ptr_res%status = NF_STATUS_OK
        ptr_res%msg = "Convergence in terms of gtol"
    case (5)
        ptr_res%status = NF_STATUS_MAX_EVAL
        ptr_res%msg = "Exceeded max. number of function evaluations"
    case (6)
        ptr_res%status = NF_STATUS_NOT_CONVERGED
        ptr_res%msg = "ftol too small. No further reduction possible"
    case (7)
        ptr_res%status = NF_STATUS_NOT_CONVERGED
        ptr_res%msg = "xtol is too small. No further improvement in solution possible"
    case (8)
        ptr_res%status = NF_STATUS_NOT_CONVERGED
        ptr_res%msg = "gtol is too small. fvec is orthogonal to columns of Jacobian"
    end select


100 continue

    call result_update (ptr_res, x=x, fx=fx, nfev=nfev)

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

    ! Clean up local WORKSPACE object if none was passed by client code
    call assert_dealloc_ptr (work, ptr_work)

    contains

    ! wrapper function: MINPACK's lmder requires a somewhat inconvenient
    ! function signature, so provide a wrapper around user-supplied function
    subroutine fcn_wrapper (m, n, x, fx, fjac, ldfjac, iflag)
        integer :: m, n, iflag, ldfjac
        real (PREC) :: x(n), fx(m), fjac(ldfjac, n)

        select case (iflag)
        case (1)
            ! Need to evaluate objective function, leave Jacobian unchanged
            call dispatch (fcn, x, fx)
        case (2)
            ! Need to evaluate Jacobian, leave function value unchanged
            call dispatch_jac (fcn, x, fjac)
        end select
    end subroutine
end subroutine


subroutine lm_check_input (m, n, ftol, xtol, gtol, maxfev, factor, eps, &
        status, msg)
    !*  LM_CHECK_INPUT performs input validation for LMDER and LMDIF
    !   routines.

    integer, intent(in) :: m, n
    real (PREC), intent(in), optional :: ftol, xtol, gtol
    integer, intent(in), optional :: maxfev
    real (PREC), intent(in), optional :: factor
    real (PREC), intent(in), optional :: eps
    type (status_t), intent(out) :: status
    character (*), intent(out) :: msg

    status = NF_STATUS_OK

    if (m < n) then
        msg = "Least-squares minimization requires m >= n"
        status = NF_STATUS_INVALID_ARG
        return
    end if

    call check_nonneg (1.0_PREC, ftol, 'ftol', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) return

    call check_nonneg (1.0_PREC, xtol, 'xtol', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) return

    call check_nonneg (1.0_PREC, gtol, 'gtol', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) return

    call check_nonneg (1, maxfev, 'maxfev', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) return

    call check_positive (1.0_PREC, factor, 'factor', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) return

    call check_positive (1.0_PREC, eps, 'eps', status, msg)
    if (NF_STATUS_INVALID_ARG .in. status) return


end subroutine

