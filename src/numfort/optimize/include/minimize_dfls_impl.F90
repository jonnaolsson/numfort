
pure function get_print_level (val) result(res)
    integer (NF_ENUM_KIND), intent(in) :: val
    integer :: res

    select case (val)
    case (NF_PRINT_ALL)
        res = 3
    case (NF_PRINT_VERBOSE)
        res = 2
    case (NF_PRINT_MINIMAL)
        res = 1
    case default
        res = 0
    end select
end function



subroutine dfls_check_input (x, m, rhobeg, rhoend, nint, iprint, maxfev, &
        status, msg)
    real (PREC), intent(in), dimension(:) :: x
    integer, intent(in) :: m
    real (PREC), intent(in) :: rhobeg, rhoend
    integer, intent(in), optional :: nint, iprint, maxfev
    type (status_t), intent(out) :: status
    character (*), intent(out) :: msg

    integer :: n
    integer (NF_ENUM_KIND), parameter :: IPRINT_VALID(*) = &
        [NF_PRINT_NONE, NF_PRINT_MINIMAL, NF_PRINT_VERBOSE, NF_PRINT_ALL]

    msg = ""
    status = NF_STATUS_OK

    n = size(x)

    if (n < 2) then
        msg = "Argument X must contain at least two variables"
        goto 100
    end if

    if (m < 1) then
        msg = 'Argument M has invalid value'
        goto 100
    end if

    if (present(nint)) then
        if (nint < n + 2 .or. nint > 2 * n + 1) then
            msg = "Argument NINT has invalid value"
            goto 100
        end if
    end if

    call check_positive (1, maxfev, "MAXFEV", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (1.0_PREC, rhobeg, "RHOBEG", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (1.0_PREC, rhoend, "RHOEND", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    if (rhobeg < rhoend) then
        msg = "Invalid arguments: RHOEND <= RHOBEG required"
        goto 100
    end if

    call check_enum (iprint, IPRINT_VALID, "IPRINT", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine



subroutine minimize_dfls (fcn, x, m, rhobeg, rhoend, nint, iprint, maxfev, &
        work, res)
    !*  Minime function f:R^n->R^m in the least-square sense using a
    !   modified version of Powell's NEWUOA algorithm.
    procedure (fvv_fcn) :: fcn
        !*  Objective function f:R^n->R^m to be minimized in a least-square
        !   sense.
    real (PREC), intent(inout), dimension(:), contiguous :: x
        !*  Vector containing the initial value on entry and the
        !   minimizer on exit.
    integer, intent(in) :: m
        !*  Dimension of function range R^m
    real (PREC), intent(in) :: rhobeg
        !*  Initial value for the trust region radius. Typically,
        !   RHOBEG should be about one tenth of the greatest
        !   expected change to a variable.
    real (PREC), intent(in) :: rhoend
        !*  Terminal value of the trust region radius, such that
        !   RHOEND <= RHOBEG. RHOEND should indicate the accuracy
        !   that is required in the final values of the variables.
    integer, intent(in), optional :: nint
        !*  Number of interpolation conditions.  Value must be
        !   in the interval [N+2,2*N+1] where
        !   N = SIZE(X)
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
        !*  Option log level for diagnostic messages
    integer, intent(in), optional :: maxfev
        !*  Optional maximum number of function evaluations, including the
        !   residual returned in RES. Exhaustion sets NF_STATUS_MAX_EVAL.
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res
        !*  Result with the evaluated residual at X. Numerical failure sets
        !   NF_STATUS_INVALID_STATE; only convergence sets success to true.

    type (fwrapper_vv) :: fcn_wrapper

    call wrap_procedure (fcn_wrapper, fcn=fcn)

    call minimize_dfls_impl (fcn_wrapper, x, m, rhobeg, rhoend, nint, iprint, &
        maxfev, work, res)

end subroutine


subroutine minimize_dfls_args (fcn, x, args, m, rhobeg, rhoend, nint, iprint, &
        maxfev, work, res)
    !*  Minime function f:R^n->R^m in the least-square sense using a
    !   modified version of Powell's NEWUOA algorithm.
    procedure (fvv_fcn) :: fcn
        !*  Objective function f:R^n->R^m to be minimized in a least-square
        !   sense.
    real (PREC), intent(inout), dimension(:), contiguous :: x
        !*  Vector containing the initial value on entry and the
        !   minimizer on exit.
    class (args_data), intent(inout) :: args
        !*  User-supplied arguments to objective function
    integer, intent(in) :: m
        !*  Dimension of function range R^m
    real (PREC), intent(in) :: rhobeg
        !*  Initial value for the trust region radius. Typically,
        !   RHOBEG should be about one tenth of the greatest
        !   expected change to a variable.
    real (PREC), intent(in) :: rhoend
        !*  Terminal value of the trust region radius, such that
        !   RHOEND <= RHOBEG. RHOEND should indicate the accuracy
        !   that is required in the final values of the variables.
    integer, intent(in), optional :: nint
        !*  Number of interpolation conditions.  Value must be
        !   in the interval [N+2,2*N+1] where
        !   N = SIZE(X)
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
        !*  Option log level for diagnostic messages
    integer, intent(in), optional :: maxfev
        !*  Optional maximum number of function evaluations, including the
        !   residual returned in RES. Exhaustion sets NF_STATUS_MAX_EVAL.
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res
        !*  Result with the evaluated residual at X. Numerical failure sets
        !   NF_STATUS_INVALID_STATE; only convergence sets success to true.

    type (fwrapper_vv) :: fcn_wrapper

    call wrap_procedure (fcn_wrapper, fcn=fcn, args=args)

    call minimize_dfls_impl (fcn_wrapper, x, m, rhobeg, rhoend, nint, iprint, &
        maxfev, work, res)

end subroutine



subroutine minimize_dfls_impl (fcn, x, m, rhobeg, rhoend, nint, iprint, maxfev, &
        work, res)
    type (fwrapper_vv), intent(inout) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
        !*  Vector containing the initial value on entry and the
        !   minimizer on exit.
    integer, intent(in) :: m
    real (PREC), intent(in) :: rhobeg
        !*  Initial value for the trust region radius. Typically,
        !   RHOBEG should be about one tenth of the greatest
        !   expected change to a variable.
    real (PREC), intent(in) :: rhoend
        !*  Terminal value of the trust region radius, such that
        !   RHOEND <= RHOBEG. RHOEND should indicate the accuracy
        !   that is required in the final values of the variables.
    integer, intent(in), optional :: nint
        !*  Number of interpolation conditions.  Value must be
        !   in the interval [N+2,2*N+1] where
        !   N = SIZE(X)
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
        !*  Option log level for diagnostic messages
    integer, intent(in), optional :: maxfev
        !*  Optional maximum number of function evaluations
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    integer :: lmaxfev, lnint, nfev, info
    integer (NF_ENUM_KIND) :: liprint, iprint_impl
    type (workspace), pointer :: ptr_work
    type (optim_result), pointer :: ptr_res

    integer :: nrwrk, n
    real (PREC), dimension(:), pointer, contiguous :: ptr_w
    type (status_t) :: status
    real (PREC), dimension(:), allocatable :: fx

    status = NF_STATUS_OK
    nfev = 0

    nullify (ptr_work, ptr_res, ptr_w)

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    call dfls_check_input (x, m, rhobeg, rhoend, nint, iprint, maxfev, status, &
        ptr_res%msg)
    if (NF_STATUS_INVALID_ARG .in. status) goto 100

    n = size(x)
    lmaxfev = 100 * (n+1)
    liprint = NF_PRINT_NONE
    lnint = 2 * n + 1
    if (present(maxfev)) lmaxfev = maxfev
    if (present(iprint)) liprint = iprint
    if (present(nint)) lnint = nint

    ! Allocate workspace arrays
    nrwrk = (lnint+11)*(lnint+n) + n*(3*n+11)/2 + 3*m*n + m*(n+1)*n/2 &
        + m*lnint + 7*m + n

    call assert_alloc_ptr (work, ptr_work)
    call workspace_reset (ptr_work)

    call assert_alloc (ptr_work, nrwrk=nrwrk)
    call workspace_get_ptr (ptr_work, nrwrk, ptr_w, status)
    if (status /= NF_STATUS_OK) goto 100

    ! Convert IPRINT to value expected by underlying implementation
    iprint_impl = get_print_level (liprint)

    allocate (fx(m))
    call newuoa2 (fcn_wrapper, n, m, lnint, x, rhobeg, rhoend, iprint_impl, &
        lmaxfev, nfev, ptr_w, info=info, fx=fx)

    select case (info)
    case (NEWUOA_CONVERGED)
        status = NF_STATUS_OK
    case (NEWUOA_MAXFUN)
        status = NF_STATUS_MAX_EVAL
        ptr_res%msg = "Maximum number of function evaluations reached"
    case (NEWUOA_INVALID_INPUT)
        status = NF_STATUS_INVALID_ARG
        ptr_res%msg = "Invalid NEWUOA input"
    case default
        status = NF_STATUS_INVALID_STATE
        ptr_res%msg = "NEWUOA encountered a nonfinite value or a failed model step"
    end select
    status%code_orig = info

100 continue

    if (present(res)) then
        call result_update (ptr_res, x, fx, nit=0, nfev=nfev, status=status)
    end if

    ! Clean up local WORKSPACE object if none was passed by client code
    call assert_dealloc_ptr (work, ptr_work)

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

    contains

    subroutine fcn_wrapper (x, fx)
        real (PREC), intent(in), dimension(:), contiguous :: x
        real (PREC), intent(out), dimension(:), contiguous :: fx

        call dispatch (fcn, x, fx)
    end subroutine

end subroutine


