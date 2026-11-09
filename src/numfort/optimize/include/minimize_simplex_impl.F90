


subroutine check_input (tol, maxfun, rstep, xstep, status, msg)
    real (PREC), intent(in) :: tol
    integer, intent(in) :: maxfun
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    type (status_t), intent(inout) :: status
    character (*), intent(inout) :: msg

    status = NF_STATUS_OK

    call check_positive (0, maxfun, "maxfun", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, tol, "tol", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_nonneg (0.0_PREC, xstep, "xstep", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    call check_nonneg (0.0_PREC, rstep, "rstep", status, msg)
    if (status /= NF_STATUS_OK) goto 100

    ! Check that not both xstep and rstep are zero
    if (present(rstep) .and. present(xstep)) then
        if (rstep == 0.0_PREC .and. xstep == 0.0_PREC) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine



pure function map_iprint (i) result(k)
    !*  MAP_IPRINT converts NUMFORT print flags to values expected by
    !   wrapped routine.
    integer (NF_ENUM_KIND), intent(in) :: i
    integer :: k

    select case (i)
    case (NF_PRINT_MINIMAL)
        k = 0
    case (NF_PRINT_VERBOSE)
        k = 1
    case (NF_PRINT_ALL)
        k = 1
    case default
        ! do not print anything by default
        k = -1
    end select
end function



pure subroutine map_ifault (ifault, status, msg)
    !*  MAP_IFAULT maps the status code returned by the underlying implementation
    !   into corresponding NUMFORT status codes.
    integer, intent(in) :: ifault
    !!  Status code as returned by simplex routine
    type (status_t), intent(out) :: status
    !!  On exit, contains the corresponding NF status code
    character (len=*), intent(out), optional :: msg

    status = NF_STATUS_UNKNOWN
    if (ifault == 0) then
        status = NF_STATUS_OK
        if (present(msg)) msg = "Simplex: successful termination"
    else if (ifault == 1) then
        status = NF_STATUS_MAX_EVAL
        if (present(msg)) msg = "Simplex: function evaluation limit reached"
    else if (ifault == 5) then
        status = NF_STATUS_OK + NF_STATUS_APPROX
        if (present(msg)) msg = "Simplex: converged; quadratic fit unavailable"
    else if (ifault == 6) then
        status = NF_STATUS_DOMAIN_ERROR
        if (present(msg)) msg = "Simplex: nonfinite search point or function value"
    else if (ifault == 3 .or. ifault == 4) then
        status = NF_STATUS_INVALID_ARG
        if (present(msg)) msg = "Simplex: invalid input argument(s)"
    end if
end subroutine



subroutine minimize_simplex (fcn, x, tol, maxfun, rstep, &
        xstep, quad, iprint, work, res)

    procedure (fvs_fcn) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxfun
    real (PREC), intent(in), optional :: rstep
        !*  Optional relative initial step size (default: 5% of initial x-value)
        !   The actual step size is computed as (xstep + rstep * x).
    real (PREC), intent(in), optional :: xstep
        !*  Optional absolute initial step size (default: 2.5e-4)
        !   The actual step size is computed as (xstep + rstep * x).
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    logical, intent(in), optional :: quad
        !* Fit a quadratic surface after convergence (default: true).
        !  If the fit is unavailable, return the best finite evaluated point
        !  with NF_STATUS_OK + NF_STATUS_APPROX. An exhausted callback budget
        !  returns NF_STATUS_MAX_EVAL, including during the fit.
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vs) :: fwrapper

    call wrap_procedure (fwrapper, fcn=fcn)

    call simplex_impl (fwrapper, x, tol, maxfun, rstep, xstep, quad, iprint, work, res)
end subroutine


subroutine minimize_simplex_args (fcn, x, args, tol, maxfun, &
        rstep, xstep, quad, iprint, work, res)

    procedure (fvs_fcn_args) :: fcn
    real (PREC), intent(inout), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxfun
    real (PREC), intent(in), optional :: rstep
        !*  Optional relative initial step size (default: 5% of initial x-value)
        !   The actual step size is computed as (xstep + rstep * x).
    real (PREC), intent(in), optional :: xstep
        !*  Optional absolute initial step size (default: 2.5e-4)
        !   The actual step size is computed as (xstep + rstep * x).
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    logical, intent(in), optional :: quad
        !* Fit a quadratic surface after convergence (default: true).
        !  If the fit is unavailable, return the best finite evaluated point
        !  with NF_STATUS_OK + NF_STATUS_APPROX. An exhausted callback budget
        !  returns NF_STATUS_MAX_EVAL, including during the fit.
    type (workspace), intent(inout), optional :: work
    type (optim_result), intent(inout), optional :: res

    type (fwrapper_vs) :: fwrapper

    call wrap_procedure (fwrapper, fcn_args=fcn, args=args)

    call simplex_impl (fwrapper, x, tol, maxfun, rstep, xstep, quad, iprint, work, res)
end subroutine


subroutine simplex_impl (fcn, x, tol, maxfun, rstep, xstep, &
        quad, iprint, work, res)

    type (fwrapper_vs), intent(inout) :: fcn

    real (PREC), intent(inout), dimension(:), contiguous :: x
    real (PREC), intent(in), optional :: tol
    integer, intent(in), optional :: maxfun
    real (PREC), intent(in), optional :: rstep
    real (PREC), intent(in), optional :: xstep
    integer (NF_ENUM_KIND), intent(in), optional :: iprint
    logical, intent(in), optional :: quad
    type (workspace), intent(inout), target, optional :: work
    type (optim_result), intent(inout), target, optional :: res

    integer :: n, lmaxfun, nloop, iquad, liprint, ifault, nrwrk
    logical :: lquad
    real (PREC) :: lrstep, lxstep, simp, fopt, ltol
    real (PREC), dimension(:), pointer, contiguous :: ptr_step, ptr_var

    type (workspace), pointer :: ptr_work
    type (optim_result), pointer :: ptr_res

    nullify (ptr_work, ptr_res)
    nullify (ptr_step, ptr_var)

    call assert_alloc_ptr (res, ptr_res)
    call result_reset (ptr_res)

    n = size(x)
    ltol = 1.0e-4_PREC
    lmaxfun = 100 * n
    lquad = .true.
    iquad = 0
    ! Default relative step size
    lrstep = 0.05_PREC
    ! Default absolute step size
    lxstep = 2.5e-4_PREC

    if (present(tol)) ltol = tol
    if (present(maxfun)) lmaxfun = maxfun
    if (present(quad)) lquad = quad
    if (lquad) iquad = 1
    if (present(rstep)) lrstep = rstep
    if (present(xstep)) lxstep = xstep

    ! disable diagnostic printing
    liprint = map_iprint(NF_PRINT_NONE)
    if (present(iprint)) liprint = map_iprint (iprint)

    call check_input (ltol, lmaxfun, lrstep, lxstep, ptr_res%status, ptr_res%msg)
    if (NF_STATUS_INVALID_ARG .in. ptr_res%status) goto 100

    nrwrk = 2 * n

    ! Allocate workspace arrays
    call assert_alloc_ptr (work, ptr_work)
    ! Clear any internal state in workspace object, in particular index offsets
    ! (this does not deallocate working arrays)
    call workspace_reset (ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk)

    call workspace_get_ptr (ptr_work, n, ptr_step)
    call workspace_get_ptr (ptr_work, n, ptr_var)

    ! set up initial step size
    ptr_step(:) = max(lxstep + lrstep * x, 1.0e-8_PREC)

    nloop = 2 * n
    simp = 1.0e-6_PREC

    call minim (x, ptr_step, n, fopt, lmaxfun, liprint, ltol, nloop, iquad, &
        simp, ptr_var, fcn_wrapper, ifault)

    ! map mimum status into numfort_optimize status
    call map_ifault (ifault, ptr_res%status, ptr_res%msg)

100 continue

    if (present(res)) then
        call result_update (ptr_res, x, fopt, nfev=fcn%nfev)
    end if

    ! Clean up local WORKSPACE object if none was passed by client code
    call assert_dealloc_ptr (work, ptr_work)

    ! Clean up local OPTIM_RESULT object if none was passed by client code
    call assert_dealloc_ptr (res, ptr_res)

    contains

    subroutine fcn_wrapper (x, fx)
        real (PREC), intent(in), dimension(:), contiguous :: x
        real (PREC), intent(out) :: fx

        call dispatch (fcn, x, fx)
    end subroutine

end subroutine
