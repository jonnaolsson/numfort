module numfort_interpolate_fitpack

    use iso_fortran_env
    use numfort_common
    use numfort_common_workspace
    use numfort_interpolate_common, mask_check_input_ext => check_input_ext
    use fitpack_real64, only: fitpack_curfit_real64 => curfit, &
        fitpack_concon_real64 => concon, fitpack_splev_real64 => splev, &
        fitpack_splder_real64 => splder

    implicit none
    private

    public :: NF_STATUS_CURFIT_INTERP
    public :: NF_STATUS_CURFIT_WLS

    integer, parameter :: PREC = real64

    integer (NF_ENUM_KIND), parameter :: NF_STATUS_CURFIT_INTERP = ishft(1, 0)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_CURFIT_WLS = ishft(1, 1)

    integer, parameter :: MIN_SPLINE_DEGREE = 1, MAX_SPLINE_DEGREE = 5, &
        DEFAULT_SPLINE_DEGREE = 3

    interface splev
        module procedure splev_real64, splev_scalar_real64
    end interface

    interface curfit
        module procedure curfit_real64
    end interface

    interface concon
        module procedure concon_real64
    end interface

    interface splder
        module procedure splder_real64, splder_scalar_real64
    end interface

    public :: curfit_get_nest, curfit, splev, splder
    public :: concon_get_nest, concon

contains

! ******************************************************************************
! INPUT CHECKING used by various routines.

!>  Verify that arrays ore of equal length, conditional on being present.
!   If arr2 is not present, array length is considered to be equal.
pure function check_length (arr1, arr2) result(status)
    real (PREC), intent(in), dimension(:) :: arr1
    real (PREC), intent(in), dimension(:), optional :: arr2
    type (status_t) :: status

    status = NF_STATUS_OK
    if (present(arr2)) then
        if (size(arr1) /= size(arr2)) then
            status = NF_STATUS_INVALID_ARG
        end if
    end if

end function

pure subroutine check_input (x, y, w, k, knots, coefs, maxiter, v, sx, bind, status, msg)
    real (PREC), intent(in), dimension(:) :: x, y, knots, coefs
    integer, intent(in), optional :: k
    integer, intent(in), optional :: maxiter
    real (PREC), intent(in), dimension(:), optional :: w
        !!  Vector of weights, defaults to 1.0
    real (PREC), intent(in), dimension(:), optional :: v
        !!  Concavity restrictions, used by concon
    real (PREC), intent(in), dimension(:), optional :: sx
        !!  Used by concon
    logical, intent(in), dimension(:), optional :: bind
        !!  Array containing flag whether constraint is binding: used for concon

    type (status_t), intent(out) :: status
    character (len=*), intent(out), optional :: msg

    status = NF_STATUS_OK

    if (size(x) /= size(y)) then
        if (present(msg)) msg = "Non-conformable arrays x, y"
        goto 100
    end if

    if (check_length (x, w) /= NF_STATUS_OK) then
        if (present(msg)) msg = "Non-conformable arrays x, w"
        goto 100
    end if

    ! check for supported spline degree
    if (present(k)) then
        if (k < MIN_SPLINE_DEGREE .or. k > MAX_SPLINE_DEGREE) then
            if (present(msg)) msg = "Unsupported spline degree k"
            goto 100
        end if
    end if

    if (present(maxiter)) then
        if (maxiter <= 0) then
            if (present(msg)) msg = "Max. number of iterations must be positive integer"
            goto 100
        end if
    end if

    if (size(knots) /= size(coefs)) then
        if (present(msg)) msg = "knots and coefs arrays must be of equal size"
        goto 100
    end if

    if (check_length (x, sx) /= NF_STATUS_OK) then
        if (present(msg)) msg = "x and sx must be of equal size"
        goto 100
    end if

    ! bind is logical array, don't bother to create check_length() that
    ! accepts logical arguments.
    if (present(bind)) then
        if (size(bind) /= size(knots)) then
            if (present(msg)) msg = "knots and bind arrays must be of equal size"
            goto 100
        end if
    end if

    if (check_length (x, v) /= NF_STATUS_OK) then
        if (present(msg)) msg = "x and v must be of equal size"
        goto 100
    end if

    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine

pure subroutine check_eval_input (knots, coefs, k, x, y, status)
    real (PREC), intent(in), dimension(:) :: x, y, knots, coefs
    integer, intent(in) :: k
    type (status_t), intent(out) :: status

    status = NF_STATUS_OK

    if (size(knots) /= size(coefs)) goto 100

    ! check for supported spline degree
    if (k < MIN_SPLINE_DEGREE .or. k > MAX_SPLINE_DEGREE) goto 100

    if (size(x) /= size(y)) goto 100

    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine

pure subroutine check_input_ext (ext, status, msg)
    integer, intent(in) :: ext
    type (status_t), intent(out) :: status
    character (len=*), intent(out), optional :: msg

    integer :: i
    logical :: is_valid
    integer, parameter :: valid(4) = [ &
        NF_INTERP_EVAL_ZERO, NF_INTERP_EVAL_EXTRAPOLATE, &
        NF_INTERP_EVAL_ERROR, NF_INTERP_EVAL_BOUNDARY]

    is_valid = .false.

    do i = 1, size(valid)
        is_valid = is_valid .or. (ext == valid(i))
    end do

    if (.not. is_valid) then
        status = NF_STATUS_INVALID_ARG
        if (present(msg)) msg = "Invalid value extrapolation mode"
    else
        status = NF_STATUS_OK
    end if
end subroutine

! ******************************************************************************
! CURFIT fitting procedures

pure subroutine curfit_real64 (x, y, k, s, knots, coefs, n, &
        iopt, w, xb, xe, work, maxiter, ssr, status, msg)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: y
    integer, intent(in), optional :: k
    real (PREC), intent(in), optional :: s
    real (PREC), intent(inout), dimension(:), contiguous :: knots
    real (PREC), intent(out), dimension(:), contiguous :: coefs
    integer, intent(out) :: n
    integer, intent(in), optional :: iopt
    real (PREC), intent(in), dimension(:), contiguous, optional :: w
    real (PREC), intent(in), optional :: xe
    real (PREC), intent(in), optional :: xb
    type (workspace_real64), intent(inout), optional, target :: work
    integer, intent(in), optional :: maxiter
    real (PREC), intent(out), optional :: ssr
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    integer :: m, liopt, lk, nwrk, nwrk_tot, nest, ier, lmaxiter
    type (status_t) :: lstatus
    real (PREC) :: lxb, lxe, ls, lssr
    type (workspace_real64), pointer :: ptr_work

    nullify(ptr_work)

    call check_input (x, y, w, k, knots, coefs, maxiter=maxiter, status=lstatus, msg=msg)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    ! initialize default values
    lk = DEFAULT_SPLINE_DEGREE
    ! by default to LS on given set of knots
    liopt = -1
    m = size(x)
    lxb = x(1)
    lxe = x(m)
    ls = 0.0_PREC
    lmaxiter = 20
    ! lower bound of recommended interval for s
    if (present(w)) ls = m - sqrt(2*real(m, PREC))

    ! override default values if any of the optional parameters were provided
    if (present(iopt)) liopt = iopt
    if (present(xb)) lxb = xb
    if (present(xe)) lxe = xe
    if (present(k)) lk = k
    if (present(s)) ls = s
    if (present(maxiter)) lmaxiter = maxiter

    call curfit_get_wrk_size (m, lk, nwrk, nest)

    if (size(knots) < nest .or. size(coefs) < nest) then
        if (present(msg)) msg = "knots and coefs arrays must be of length nest"
        goto 100
    end if

    nwrk_tot = nwrk
    ! Add enough space to allocate default weights on workspace array
    if (.not. present(w)) nwrk_tot = nwrk_tot + m

    ! assert that working arrays are sufficiently large
    call assert_alloc_ptr (work, ptr_work)
    call assert_alloc (ptr_work, nrwrk=nwrk_tot, niwrk=nest)

    if (.not. present(w)) then
        ! allocate real working array such that we can append
        ! weights at the end
        ptr_work%rwrk(nwrk+1:nwrk+m) = 1.0_PREC
    end if

    ! Erase values in knots / coefs arrays to eliminate any junk that might
    ! be at those locations.
    if (liopt == 0) then
        knots = 0.0_PREC
        coefs = 0.0_PREC
    end if

    ! call fitpack routine
    if (present(w)) then
        call fitpack_curfit_real64 (liopt, m, x, y, w, lxb, lxe, lk, &
            ls, nest, n, knots, coefs, lssr, &
            ptr_work%rwrk(1:nwrk), nwrk, ptr_work%iwrk, lmaxiter, ier)
    else
        call fitpack_curfit_real64 (liopt, m, x, y, ptr_work%rwrk(nwrk+1:nwrk+m), &
            lxb, lxe, lk, ls, nest, n, knots, coefs, lssr, &
            ptr_work%rwrk(1:nwrk), nwrk, ptr_work%iwrk, lmaxiter, ier)
    end if

    ! Map CURFIT ier code to NF status code
    select case (ier)
    case (0)
        lstatus = NF_STATUS_OK
    case (-1)
        ! retuned spline is an interpolating spline, ie ssr=0
        lstatus = NF_STATUS_OK
        lstatus = lstatus + NF_STATUS_CURFIT_INTERP
    case (-2)
        ! returned spline is an WLS polynomial of degree k
        lstatus = NF_STATUS_OK
        lstatus = lstatus + NF_STATUS_CURFIT_WLS
    case (1)
        ! nest is too small, probably because s is too small.
        ! CURFIT returns an approximation in this case
        lstatus = NF_STATUS_STORAGE_ERROR
        lstatus = lstatus + NF_STATUS_APPROX
    case (2)
        ! Theoretically impossible result encountered, probably because s
        ! is too small
        lstatus = NF_STATUS_INVALID_STATE
        lstatus = lstatus + NF_STATUS_APPROX
    case (3)
        ! Max. number of iterations exceeded, probably because s is too small.
        lstatus = NF_STATUS_MAX_ITER
        lstatus = lstatus + NF_STATUS_APPROX
    case (10)
        lstatus = NF_STATUS_INVALID_ARG
    case default
        lstatus = NF_STATUS_UNKNOWN
    end select

    ! Store original status code
    lstatus%code_orig = ier
    ! Store status message
    if (present(msg)) call curfit_get_msg (ier, msg)
    if (present(ssr)) ssr = lssr

100 continue

    if(present(status)) status = lstatus

    ! Clean up (local) workspace object
    call assert_dealloc_ptr (work, ptr_work)

end subroutine


pure subroutine curfit_get_msg (ier, msg)
    integer, intent(in) :: ier
    character (*), intent(out) :: msg

    select case (ier)
    case (-2)
        msg = "Spline is a weighted LS polynomial of degree k" &
            // " (SSR is upper bound for smoothing factor s)"
    case (-1)
        msg = "Spline is an interpolating spline (SSR=0)."
    case (0)
        msg = "Spline SSR satisfies (SSR-s)/s <= 0.001"
    case (1)
        msg = "Max. storage space exceeded (nest or smoothing parameter too small?)"
    case (2)
        msg = "Theoretically impossible result encountered (s too small?)"
    case (3)
        msg = "Max. number of iterations exceeded"
    case (10)
        msg = "Invalid input argument values"
    case default
        msg = "Unknown error code"
    end select

end subroutine


!-------------------------------------------------------------------------------
! CONCON wrapper routine

pure subroutine concon_real64 (x, y, v, s, &
        knots, coefs, n, iopt, w, maxtr, maxbin, work, ssr, sx, bind, &
        status, msg)

    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous :: y
    real (PREC), intent(in), dimension(:), contiguous :: v
    real (PREC), intent(in), optional :: s
    real (PREC), intent(inout), dimension(:), contiguous :: knots
    real (PREC), intent(out), dimension(:), contiguous :: coefs
    integer, intent(out) :: n
    integer, intent(in), optional :: iopt
    real (PREC), intent(in), dimension(:), contiguous, optional :: w
    integer, intent(in), optional :: maxtr
    integer, intent(in), optional :: maxbin
    type (workspace_real64), intent(inout), optional, target :: work
    real (PREC), intent(out), optional :: ssr
    real (PREC), intent(out), dimension(:), optional, contiguous :: sx
    logical, intent(out), dimension(:), optional, contiguous :: bind
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg

    target :: bind, sx

    integer :: m, liopt, nest, ier, lmaxtr, lmaxbin
    integer :: nrwrk, jrwrk, niwrk, nlwrk, nrwrk_tot
    real (PREC) :: ls, lssr
    type (status_t) :: lstatus
    type (workspace_real64), pointer :: ptr_work
    real (PREC), dimension(:), pointer, contiguous :: ptr_sx, ptr_w, ptr_v
    logical, dimension(:), pointer, contiguous :: ptr_bind

    nullify (ptr_work, ptr_sx, ptr_w, ptr_bind, ptr_v)

    call check_input (x, y, w, knots=knots, coefs=coefs, v=v, bind=bind, sx=sx,&
        status=lstatus, msg=msg)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    ! default values
    liopt = 0
    m = size(x)
    ls = 0.0_PREC
    ! override default values if any of the optional parameters were provided
    ! lower bound of recommended interval for s
    if (present(w)) ls = m - sqrt(2*real(m, PREC))
    if (present(s)) ls = s
    if (present(iopt)) liopt = iopt

    nest = size(knots)
    lmaxbin = nest - 6
    lmaxtr = 200
    if (present(maxbin)) lmaxbin = maxbin
    if (present(maxtr)) lmaxtr = maxtr

    call concon_get_wrk_size (m, lmaxtr, lmaxbin, nest, nrwrk, niwrk)

    ! work size for logical array
    nlwrk = 0
    ! offset for real work array; initially set to end of work array not
    ! including slices used for weights and sx
    nrwrk_tot = nrwrk
    ! elements containing bind array
    if (.not. present(bind)) nlwrk = nlwrk + nest
    ! need to allocate workspace for uniform weights
    if (.not. present(w)) nrwrk_tot = nrwrk_tot + m
    ! need to allocate workspace for sx
    if (.not. present(sx)) nrwrk_tot = nrwrk_tot + m
    ! need to allocate workspace for v as CONCON modifies v in-place
    nrwrk_tot = nrwrk_tot + m

    call assert_alloc_ptr (work, ptr_work)
    call assert_alloc (ptr_work, nrwrk=nrwrk_tot, niwrk=niwrk, nlwrk=nlwrk)

    ! allocate default weights and sx after the actual real working array
    jrwrk = nrwrk + 1
    ! assert that working arrays are sufficiently large
    if (.not. present(w)) then
        ! allocate real working array such that we can append uniform
        ! weights at the end
        ptr_w => ptr_work%rwrk(jrwrk:jrwrk + m - 1)
        ptr_w = 1.0_PREC
        jrwrk = jrwrk + m
    end if

    ! copy v into workspace as it might be modified by CONCON
    ptr_v => ptr_work%rwrk(jrwrk:jrwrk+m-1)
    ptr_v = v
    jrwrk = jrwrk + m
    if (present(sx)) then
        ptr_sx => sx
    else
        ptr_sx => ptr_work%rwrk(jrwrk:jrwrk+m-1)
        jrwrk = jrwrk + m
    end if

    if (present(bind)) then
        ptr_bind => bind
    else
        ptr_bind => ptr_work%lwrk(1:nest)
    end if

    ! call fitpack routine wrapper
    if (present(w)) then
        call fitpack_concon_real64 (liopt, m, x, y, w, ptr_v, ls, nest, lmaxtr, &
            lmaxbin, n, knots, coefs, lssr, ptr_sx, ptr_bind, &
            ptr_work%rwrk(1:nrwrk), nrwrk, ptr_work%iwrk, niwrk, ier)
    else
        call fitpack_concon_real64 (liopt, m, x, y, ptr_w, ptr_v, ls, nest, lmaxtr, &
            lmaxbin, n, knots, coefs, lssr, ptr_sx, ptr_bind, &
            ptr_work%rwrk(1:nrwrk), nrwrk, ptr_work%iwrk, niwrk, ier)
    end if

    if (ier == 0) then
        lstatus = NF_STATUS_OK
    else if (ier == -3) then
        ! this should not happen if nest=m+4, as returned by concon_get_nest()
        lstatus = ior(NF_STATUS_STORAGE_ERROR, NF_STATUS_APPROX)
    else if (ier == -2 .or. ier == -1) then
        lstatus = NF_STATUS_APPROX
    else if (ier == 1 .or. ier == 2) then
        lstatus = NF_STATUS_STORAGE_ERROR
    else if (ier == 3 .or. ier == 4 .or. ier == 5) then
        lstatus = NF_STATUS_OTHER
    else if (ier == 10) then
        lstatus = NF_STATUS_INVALID_ARG
    else
        lstatus = NF_STATUS_UNKNOWN
    end if

    ! Store original status code
    lstatus%code_orig = ier
    ! Store status message if present
    if (present(msg)) call concon_get_msg (ier, msg)
    if (present(ssr)) ssr = lssr

100 continue
    if(present(status)) status = lstatus
    call assert_dealloc_ptr (work, ptr_work)

end subroutine


pure subroutine concon_get_msg (ier, msg)
    integer, intent(in) :: ier
    character (*), intent(out) :: msg

    select case (ier)
    case (-3)
        msg = "Available storage exceeded (nest too small?)"
    case (-2)
        msg = "Maximal number of knots reached (s too small?)"
    case (-1)
        msg = "Adding knots does not reduce SSR (s too small?)"
    case (0)
        msg = "Spline satisfies concavity/convexity constraints"
    case (1)
        msg = "Number of knots with s''(x)=0 exceeds maxbin (maxbin too small?)"
    case (2)
        msg = "Number of records in tree exceeds maxtr"
    case (3)
        msg = "No solution found for quadratic programming problem"
    case (4)
        msg = "Number of knots required to satisfy constraints exceeds nest (nest too small?)"
    case (5)
        msg = "Number of knots required to satisfy constraints exceeds m+4"
    case (10)
        msg = "Invalid input argument values"
    case default
        msg = "Unknown error flag"
    end select

end subroutine


! ******************************************************************************
! CURFIT evaluation

pure subroutine splev_real64 (knots, coefs, k, x, y, ext, status)
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    integer, intent(in), optional :: k
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: y
    integer, intent(in), optional :: ext
    type (status_t), intent(out), optional :: status

    integer :: m, ln, lier, lext, lk
    type (status_t) :: lstatus

    m = size(x)
    ln = size(knots)
    ! assume cubic splines by default
    lk = DEFAULT_SPLINE_DEGREE
    ! by default set function values outside of domain to boundary values
    lext = NF_INTERP_EVAL_EXTRAPOLATE

    if (present(k)) lk = k
    if (present(ext)) lext = ext

    call check_eval_input (knots, coefs, lk, x, y, status=lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    ! Map NUMFORT extrapolation constants to values used by FITPACK
    select case (lext)
    case (NF_INTERP_EVAL_EXTRAPOLATE)
        lext = 0
    case (NF_INTERP_EVAL_ZERO)
        lext = 1
    case (NF_INTERP_EVAL_BOUNDARY)
        lext = 3
    case (NF_INTERP_EVAL_ERROR)
        lext = 2
    case default
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end select

    call fitpack_splev_real64 (knots, ln, coefs, lk, x, y, m, lext, lier)

    ! remap SPLEV status codes to NUMFORT codes
    select case (lier)
    case (0)
        lstatus = NF_STATUS_OK
    case (1)
        lstatus = NF_STATUS_BOUNDS_ERROR
    case (10)
        lstatus = NF_STATUS_INVALID_ARG
    case default
        lstatus = NF_STATUS_UNKNOWN
    end select

100 continue
    if (present(status)) status = lstatus

end subroutine

pure subroutine splev_scalar_real64 (knots, coefs, k, x, y, ext, status)
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    integer, intent(in), optional :: k
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer, intent(in), optional :: ext
    type (status_t), intent(out), optional :: status

    real (PREC) :: lx(1), ly(1)

    lx(1) = x

    call splev_real64 (knots, coefs, k, lx, ly, ext, status)

    ! write back result
    y = ly(1)
end subroutine

! ******************************************************************************
! Evaluation of derivatives

pure subroutine splder_real64 (knots, coefs, k, order, x, y, ext, work, status)
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    integer, intent(in), optional :: k
    integer, intent(in), optional :: order
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: y
    integer, intent(in), optional :: ext
    type (workspace_real64), intent(inout), optional, target :: work
    type (status_t), intent(out), optional :: status

    type (workspace_real64), pointer :: ptr_work

    integer :: lext, lorder, lier, lk, ln, m
    type (status_t) :: lstatus
    
    nullify (ptr_work)

    lk = DEFAULT_SPLINE_DEGREE
    lorder = 1
    lext = NF_INTERP_EVAL_EXTRAPOLATE
    ln = size(knots)
    m = size(x)

    if (present(k)) lk = k
    if (present(ext)) lext = ext
    if (present(order)) lorder = order

    ! Do not bother to check for derivative order since this is performed
    ! by SPLDER.
    call check_eval_input (knots, coefs, lk, x, y, status=lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    select case (lext)
    case (NF_INTERP_EVAL_EXTRAPOLATE)
        lext = 0
    case (NF_INTERP_EVAL_ZERO)
        lext = 1
    case (NF_INTERP_EVAL_ERROR)
        lext = 2
    case default
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end select

    call assert_alloc_ptr (work, ptr_work)
    call assert_alloc (ptr_work, nrwrk=ln)

    call fitpack_splder_real64 (knots, ln, coefs, lk, lorder, x, y, m, &
        lext, ptr_work%rwrk, lier)

    ! remap SPLDER status codes to NUMFORT codes
    select case (lier)
    case (0)
        lstatus = NF_STATUS_OK
    case (1)
        lstatus = NF_STATUS_BOUNDS_ERROR
    case (10)
        lstatus = NF_STATUS_INVALID_ARG
    case default
        lstatus = NF_STATUS_UNKNOWN
    end select

100 continue
    if (present(status)) status = lstatus
    call assert_dealloc_ptr (work, ptr_work)
end subroutine

! splder_scalar_real64 evaluates the derivate of a spline at a single scalar point
! x and stores the result in scalar y.
! Implemented as a wrapper for SPLDER_WRAPPER that creates local size-1 arrays.
pure subroutine splder_scalar_real64 (knots, coefs, k, order, x, y, ext, work, status)
    real (PREC), intent(in), dimension(:), contiguous :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    integer, intent(in), optional :: k
    integer, intent(in), optional :: order
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    integer, intent(in), optional :: ext
    type (workspace_real64), intent(inout), optional, target :: work
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(1) :: lx, ly

    lx(1) = x
    ly(1) = y

    call splder_real64 (knots, coefs, k, order, lx, ly, ext, work, status)

    y = ly(1)
end subroutine

! ******************************************************************************
! Workspace query functions

pure function curfit_get_nest (m, k) result(res)
    integer, intent(in) :: k, m
    integer :: res

    if (k < MIN_SPLINE_DEGREE .or. k > MAX_SPLINE_DEGREE) then
        res = -1
    else
        res = max(m + k + 1, 2*k + 3)
    end if
end function

pure subroutine curfit_get_wrk_size (m, k, nwrk, nest)
    integer, intent(in) :: m, k
    integer, intent(out) :: nwrk, nest

    nest = curfit_get_nest (m, k)
    nwrk = m * (k + 1) + nest * (7 + 3 * k)
end subroutine

pure subroutine concon_get_wrk_size (m, maxtr, maxbin, nest, nrwrk, niwrk)
    integer :: m, maxtr, maxbin, nest, nrwrk, niwrk

    intent (in) :: m, maxtr, maxbin, nest
    intent (out) :: nrwrk, niwrk

    nrwrk = m * 4 + nest * 8 + maxbin * (maxbin + nest + 1)
    niwrk = maxtr * 4 + 2 * (maxbin + 1)
end subroutine

pure function concon_get_nest (m) result(nest)
    integer, intent(in) :: m
    integer :: nest

    ! return maximally required array size
    nest = m + 4
end function

end module
