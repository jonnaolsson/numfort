

pure subroutine linspace (x, xfrom, xto, step, res_n, res_step)
    real (PREC), intent(out), dimension(:) :: x
        !*  Array to store evenly spaced numbers over specified interval.
    real (PREC), intent(in) :: xfrom
        !*  Starting value of sequence.
    real (PREC), intent(in) :: xto
        !*  End value of sequence to be created. If step size is explicitly given
        !   and the array x is not large enough, the last element of X can be
        !   significantly different from XTO.
    real (PREC), intent(in), optional :: step
        !*  Optional step size used to create evenly spaced elements. If not present,
        !   step size is imputed from the boundary values and the size of argument X.
    integer, intent(out), optional :: res_n
        !*  If present, contains the actual number of elements returned. If
        !   the STEP argument is not present, res_n = size(x).
    real (PREC), intent(out), optional :: res_step
        !*  If present, contains the actual step size used. If the STEP argument
        !   is present, then res_step = step.

    real (PREC) :: lstep, xnext, sgn
    integer :: i, nx, n

    n = 0
    nx = size(x)
    lstep = 0.0_PREC

    if (nx == 0) goto 100

    if (.not. present(step)) then
        x(1) = xfrom
        n = 1

        if (nx > 1) then
            ! Infer step size from number boundaries and array size.
            lstep = (xto - xfrom) / (nx - 1.0_PREC)

            do i = 2, nx-1
                x(i) =  xfrom + (i-1) * lstep
            end do

            ! Avoid rounding errors at end point.
            x(nx) = xto
            n = nx
        end if
    else
        sgn = signum (xto-xfrom)
        lstep = step

        xnext = xfrom
        do i = 1, nx
            x(i) = xnext

            xnext = xnext + lstep
            if (xnext*sgn > xto*sgn) exit
        end do

        ! Actual number of elements used need not be identical to array size
        n = min(i, nx)
    end if

100 continue
    if (present(res_n)) res_n = n
    if (present(res_step)) res_step = lstep

end subroutine



pure subroutine powerspace (x, xmin, xmax, pow, sparse_end)
    !*  POWERSPACE returns a sequence of points obtained by taking the power
    !   of a sequence of uniformly spaced points on [0,1] and applying
    !   an affine transformation to match the given start and end point.
    !
    !   More specifically, each element in X is computed as
    !       x(i) = xmin + (xmax-xmin) * u(i) ** pow
    !   where u(i) is the corresponding element on a uniformly-spaced
    !   sequence on [0,1].
    real (PREC), intent(out), dimension(:) :: x
        !*  Array to store "power-spaced" sequence of points
    real (PREC), intent(in) :: xmin
        !*  Starting value
    real (PREC), intent(in) :: xmax
        !*  Endpoint value
    real (PREC), intent(in) :: pow
        !*  Exponent used to create power-spaced sequence
    logical, intent(in), optional :: sparse_end
        !*  If true (default), create a grid that is sparser towards the end.

    integer :: n, i
    real (PREC) :: slope, xi
    logical :: lend

    lend = .true.
    if (present (sparse_end)) lend = sparse_end

    n = size(x)

    ! Silently ignore 0-size arrays
    if (n == 0) return

    ! Singleton grid preserves starting point without power transformation
    if (n == 1) then
        x(1) = xmin
        return
    end if

    call linspace (x, 0.0_PREC, 1.0_PREC)

    x = x ** pow

    if (.not. lend) then
        ! Make grid more sparse at the end
        x = 1.0_PREC - x
        ! Reverse order
        do i = 1, n / 2
            xi = x(i)
            x(i) = x(n-i+1)
            x(n-i+1) = xi
        end do
    end if

    ! Transform to desired range
    slope = xmax - xmin
    x = xmin + slope * x

    ! Explicitly set boundary values to eliminate any rounding issues
    x(1) = xmin
    if (n > 1) x(n) = xmax

end subroutine



pure subroutine logspace (x, logx_min, logx_max, base)
    !*  LOGSPACE returns numbers evenly spaced on a log scale.
    !
    !   To be compatible with Numpy, the start and end points need to be
    !   provided in logs.
    real (PREC), intent(out), dimension(:) :: x
        !*  Array containing samples equally spaced on a log scale
    real (PREC), intent(in) :: logx_min
        !*  Starting value of the sequence in logs, ie. sequence starts
        !   at base ** logx_min.
    real (PREC), intent(in) :: logx_max
        !*  Terminal value of the sequence in logs, ie. sequence ends at
        !   base ** logx_max.
    real (PREC), intent(in), optional :: base
        !*  The base of the log space (default: 10)

    real (PREC) :: lbase
    integer :: i

    lbase = 10.0
    if (present(base)) lbase = base

    call linspace (x, logx_min, logx_max)

    do i = 1, size(x)
        x(i) = lbase ** x(i)
    end do

end subroutine



pure subroutine log_shift_space (x, xmin, xmax, shift, x0, frac_x0, &
        base, status)
    !*  LOG_SHIFT_SPACE returns numbers which, after adding constant, are evenly
    !   spaced in logs.
    real (PREC), intent(out), dimension(:) :: x
        !*  Array to store the resulting numbers
    real (PREC), intent(in) :: xmin
        !*  Starting value (using the non-transformed scale)
    real (PREC), intent(in) :: xmax
        !*  Terminal value (using the non-transformed scale)
    real (PREC), intent(in), optional :: shift
        !*  (Optional) constant added to each point before applying the log
        !   transformation.
    real (PREC), intent(in), optional :: x0
        !*  (Optional) reference value (using the non-transformed scale) which is
        !   used to internally compute additive constant SHIFT (default:
        !   midpoint of interval [XMIN,XMAX]).
    real (PREC), intent(in), optional :: frac_x0
        !*  (Optional) fraction of points in (0,1) that are to be located
        !   on the interval [XMIN,X0]. If present, takes precedence over
        !   the argument SHIFT.
    real (PREC), intent(in), optional :: base
        !*  (Optional) logarithm basis.
    type (status_t), intent(out), optional  :: status
        !*  (Optional) exit code.

    type (status_t) :: lstatus
    real (PREC) :: lx0, lshift
    real (PREC) :: lb, ub, mid, flb, fub, fmid, lxmin, lxmax
    integer :: i

    lstatus = NF_STATUS_OK

    call check_positive (0.0_PREC, base, 'base', status=status)
    if (status /= NF_STATUS_OK) goto 100

    lshift = 0.0
    if (present(shift)) lshift = shift

    ! --- Find shift parameter from (x0, frac_x0) ---

    if (present(frac_x0)) then
        if (frac_x0 <= 0.0_PREC .or. frac_x0 >= 1.0_PREC) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if

        ! Default reference value is interval midpoint
        lx0 = (xmin + xmax) / 2.0_PREC
        if (present(x0)) lx0 = x0

        ! Bisection to find desired shifter
        lb = 1.0e-6_PREC - xmin
        ub = xmax - xmin
        flb = fobj (lx0, lb)
        fub = fobj (lx0, ub)
        ! Use 10 attempts to find a bracket
        do i = 1, 10
            if (flb * fub < 0.0_PREC) exit
            ub = ub * 10.0
            fub = fobj (lx0, ub)
        end do

        if (flb * fub > 0.0_PREC) then
            ! Return error since
            lstatus = NF_STATUS_BOUNDS_ERROR
            goto 100
        end if

        ! Perform actual bisection
        mid = (ub + lb) / 2.0_PREC
        do i = 1, 100
            fmid = fobj (lx0, mid)
            if (abs(fmid) < 1.0e-8_PREC) exit

            if (fmid * flb >= 0.0_PREC) then
                ! f(mid) and f(lb) have the same sign
                lb = mid
            else
                ! If we had an initial bracket, then f(ub) and f(mid) must have
                ! the same sign.
                ub = mid
            end if
        end do

        lshift = mid

    end if

    ! --- Create grid ---

    lxmin = logx (xmin + lshift)
    lxmax = logx (xmax + lshift)

    call linspace (x, lxmin, lxmax)
    if (present(base)) then
        x = base ** x
    else
        x = exp (x)
    end if

    ! Undo shift
    x = x - lshift

    ! Remove any rounding errors in start/end points
    x(1) = xmin
    x(size(x)) = xmax

100 continue

    if (present(status)) status = lstatus

contains

    pure function fobj (x0, shift) result(res)
        real (PREC), intent(in) :: x0, shift
        real (PREC) :: res

        real (PREC) :: dist

        dist = logx (xmin + shift) - logx (xmax + shift)
        res = logx (x0 + shift) - logx (xmin + shift) - frac_x0 * dist
    end function

    elemental function logx (x) result (res)
        real (PREC), intent(in) :: x
        real (PREC) :: res

        res = log(x)
        if (present(base)) then
            res = res / log(base)
        end if
    end function

end subroutine

