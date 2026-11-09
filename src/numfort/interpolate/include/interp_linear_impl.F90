



pure subroutine __APPEND(interp_linear_eval_check_input,__PREC) (ilbound, &
        weight, fx, status)
    !*  INTERP_LINEAR_EVAL_CHECK_INPUT performs input validation for
    !   INTERP_LINEAR_EVAL.
    integer, parameter :: PREC = __PREC

    integer, intent(in), dimension(:) :: ilbound
    real (PREC), intent(in), dimension(:) :: weight
    real (PREC), intent(in), dimension(:) :: fx
        !*  Optional array for storing output, which is not present in
        !   scalar interface.
    type (status_t), intent(out) :: status

    status = NF_STATUS_OK

    if (size(fx) /= size(ilbound)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (size(ilbound) /= size(weight)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

100 continue

end subroutine



pure subroutine __APPEND(interp_linear_check_input,__PREC) (x, xp, fp, fx, status)
    !*  INTERP_LINEAR_CHECK_INPUT performs input validation for INTERP_LINEAR.
    !   This routine will check arguments that are passed to both the
    !   implementations of the FIND and EVAL routines so they don't need
    !   to be checked separately.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:), optional :: x, fx
        !*  Optional interpolation points, since these do not need to be
        !   checked for scalar routines.
    real (PREC), intent(in), dimension(:) :: xp, fp
    type (status_t), intent(out) :: status

    status = NF_STATUS_OK

    if (present(x) .and. present(fx)) then
        if (size(x) /= size(fx)) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    if (size(xp) /= size(fp)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (size(xp) < 2) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

100 continue

end subroutine



pure subroutine __APPEND(interp_linear_eval_impl_scalar,__PREC) (ilbound, &
        weight, fp, fx, ext, left, right)
    !*  Implements linear interpolation for scalar input/return value
    !   and no optional arguments.
    !   Should be called from wrapper routines doing the optional argument handling,
    !   etc.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in) :: ilbound
        !*  Index of the lower bound of the interpolating interval
    real (PREC), intent(in) :: weight
        !*  Weight on the lower bound of the bracketing interval
    real (PREC), intent(in), dimension(:), contiguous :: fp
        !*  y-coordinates of the data points, same length as xp.
    real (PREC), intent(out) :: fx
        !*  Interpolated value
    integer (NF_ENUM_KIND), intent(in) :: ext
        !*  Defines behavious in case function should be evaluated at point x
        !   outside of the range specified by [xp(1), xp(size(xp))].
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'left' for all
        !           x < xp(1) and the value of 'right' for all x > xp(size(xp)).
        !           If 'left' or 'right' are omitted, values fp(1) or fp(size(xp))
        !           are used, respectively.
        !       - NF_INTERP_EVAL_BOUNDARY: Returns fp(1) for all x < xp(1) and
        !           fp(size(xp)) for all x > xp(size(xp)).
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in) :: left
        !*  Value to return if x is smaller than the lower bound of the
        !   function's domain.
    real (PREC), intent(in) :: right
        !*  Value to return if x is larger then the upper bound of the
        !   function'n domain.

    integer (INTSIZE) :: np, ilb
    real (PREC) :: wgt

    if (ext == NF_INTERP_EVAL_CONST) then
        np = size(fp)
        if (ilbound == 0) then
            fx = left
            return
        else if (ilbound == np) then
            fx = right
            return
        end if
    end if

    ! Use linear interpolation / extrapolation according to bracket and
    ! weight identified previously.
    ilb = ilbound
    wgt = weight

    fx = wgt * fp(ilb) + (1.0_PREC - wgt) * fp(ilb+1)

end subroutine



pure subroutine __APPEND(interp_linear_eval_scalar,__PREC) (ilbound, &
        weight, fp, fx, ext, left, right, status)
    !*  Implements linear interpolation for scalar input/return value
    !   and no optional arguments.
    !   Should be called from wrapper routines doing the optional argument handling,
    !   etc.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in) :: ilbound
        !*  Index of the lower bound of the interpolating interval
    real (PREC), intent(in) :: weight
        !*  Weight on the lower bound of the bracketing interval
    real (PREC), intent(in), dimension(:), contiguous :: fp
        !*  y-coordinates of the data points, same length as xp.
    real (PREC), intent(out) :: fx
        !*  Interpolated value
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Defines behavious in case function should be evaluated at point x
        !   outside of the range specified by [xp(1), xp(size(xp))].
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'left' for all
        !           x < xp(1) and the value of 'right' for all x > xp(size(xp)).
        !           If 'left' or 'right' are omitted, values fp(1) or fp(size(xp))
        !           are used, respectively.
        !       - NF_INTERP_EVAL_BOUNDARY: Returns fp(1) for all x < xp(1) and
        !           fp(size(xp)) for all x > xp(size(xp)).
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in), optional :: left
        !*  Value to return if x is smaller than the lower bound of the
        !   function's domain.
    real (PREC), intent(in), optional :: right
        !*  Value to return if x is larger then the upper bound of the
        !   function'n domain.
    type (status_t), intent(out), optional :: status

    integer (NF_ENUM_KIND) :: lext
    type (status_t) :: lstatus
    real (PREC) :: lleft, lright

    lstatus = NF_STATUS_OK

    call check_input_ext (1.0_PREC, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    ! The evaluator requires constants even when the selected mode ignores them.
    lleft = 0.0_PREC
    lright = 0.0_PREC
    if (present(left)) lleft = left
    if (present(right)) lright = right

    call interp_linear_eval_impl (ilbound, weight, fp, fx, lext, lleft, lright)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(interp_linear_eval_impl_1d,__PREC) (ilbound, &
        weight, fp, fx, ext, left, right)
    !*  Implements linear interpolation for scalar input/return value
    !   and no optional arguments.
    !   Should be called from wrapper routines doing the optional argument handling,
    !   etc.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in), dimension(:), contiguous :: ilbound
        !*  Index of the lower bound of the interpolating interval
    real (PREC), intent(in), dimension(:), contiguous :: weight
        !*  Weight on the lower bound of the bracketing interval
    real (PREC), intent(in), dimension(:), contiguous :: fp
        !*  y-coordinates of the data points, same length as xp.
    real (PREC), intent(out), dimension(:), contiguous :: fx
        !*  Interpolated value
    integer (NF_ENUM_KIND), intent(in) :: ext
        !*  Defines behavious in case function should be evaluated at point x
        !   outside of the range specified by [xp(1), xp(size(xp))].
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'left' for all
        !           x < xp(1) and the value of 'right' for all x > xp(size(xp)).
        !           If 'left' or 'right' are omitted, values fp(1) or fp(size(xp))
        !           are used, respectively.
        !       - NF_INTERP_EVAL_BOUNDARY: Returns fp(1) for all x < xp(1) and
        !           fp(size(xp)) for all x > xp(size(xp)).
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in) :: left
        !*  Value to return if x is smaller than the lower bound of the
        !   function's domain.
    real (PREC), intent(in) :: right
        !*  Value to return if x is larger then the upper bound of the
        !   function'n domain.

    integer (INTSIZE) :: np, nx, ilb, i
    real (PREC) :: wgt

    nx = size(fx)

    if (ext == NF_INTERP_EVAL_CONST) then
        ! Deal with constant extrapolation separately
        np = size(fp)

        do i = 1, nx
            ilb = ilbound(i)
            if (ilb == 0) then
                fx(i) = left
            else if (ilb == np) then
                fx(i) = right
            else
                wgt = weight(i)
                fx(i) = wgt * fp(ilb) + (1.0_PREC - wgt) * fp(ilb+1)
            end if
        end do
    else
        ! If we are NOT replacing non-interior interpolant values with some
        ! contants, every other supported extrapolation method
        ! was already applied to bounds/weights as required. Nothing left to
        ! do here.
        do i = 1, nx
            wgt = weight(i)
            ilb = ilbound(i)
            fx(i) = wgt * fp(ilb) + (1.0_PREC - wgt) * fp(ilb+1)
        end do
    end if

end subroutine



pure subroutine __APPEND(interp_linear_eval_1d,__PREC) (ilbound, &
        weight, fp, fx, ext, left, right, status)
    !*  Implements linear interpolation for scalar input/return value
    !   and no optional arguments.
    !   Should be called from wrapper routines doing the optional argument handling,
    !   etc.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), intent(in), dimension(:), contiguous :: ilbound
        !*  Index of the lower bound of the interpolating interval
    real (PREC), intent(in), dimension(:), contiguous :: weight
        !*  Weight on the lower bound of the bracketing interval
    real (PREC), intent(in), dimension(:), contiguous :: fp
        !*  y-coordinates of the data points, same length as xp.
    real (PREC), intent(out), dimension(:), contiguous :: fx
        !*  Interpolated value
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Defines behavious in case function should be evaluated at point x
        !   outside of the range specified by [xp(1), xp(size(xp))].
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'left' for all
        !           x < xp(1) and the value of 'right' for all x > xp(size(xp)).
        !           If 'left' or 'right' are omitted, values fp(1) or fp(size(xp))
        !           are used, respectively.
        !       - NF_INTERP_EVAL_BOUNDARY: Returns fp(1) for all x < xp(1) and
        !           fp(size(xp)) for all x > xp(size(xp)).
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in), optional :: left
        !*  Value to return if x is smaller than the lower bound of the
        !   function's domain.
    real (PREC), intent(in), optional :: right
        !*  Value to return if x is larger then the upper bound of the
        !   function'n domain.
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer (NF_ENUM_KIND) :: lext
    real (PREC) :: lleft, lright

    call interp_linear_eval_check_input (ilbound, weight, fx, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_input_ext (1.0_PREC, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext

    ! The evaluator requires constants even when the selected mode ignores them.
    lleft = 0.0_PREC
    lright = 0.0_PREC
    if (present(left)) lleft = left
    if (present(right)) lright = right

    call interp_linear_eval_impl (ilbound, weight, fp, fx, lext, lleft, lright)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(interp_linear_scalar,__PREC) (x, xp, fp, fx, ext, &
        left, right, cache, status)
    !*  INTERP_LINEAR_SCALAR performs linear interpolation for scalar arguments
    integer, parameter :: INTSIZE = int32
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in) :: x
        !*  The x-coordinate of the interpolated value.
    real (PREC), intent(in), dimension(:), contiguous :: xp
        !*  The x-coordinates of data points in increasing order.
    real (PREC), intent(in), dimension(:), contiguous :: fp
        !*  The y-coordinates of data points, same length as xp.
    real (PREC), intent(out) :: fx
        !*  Interpolated value
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Defines behavious in case function should be evaluated at point x
        !   outside of the range specified by [xp(1), xp(size(xp))].
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'left' for all
        !           x < xp(1) and the value of 'right' for all x > xp(size(xp)).
        !           If 'left' or 'right' are omitted, values fp(1) or fp(size(xp))
        !           are used, respectively.
        !       - NF_INTERP_EVAL_BOUNDARY: Returns fp(1) for all x < xp(1) and
        !           fp(size(xp)) for all x > xp(size(xp)).
        !       - NF_INTERP_EVAL_ERROR: Routine exists with an error code if
        !           a value outside of the function's domain is encountered.
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in), optional :: left
        !*  Value to return if x < xp(1) and extrapolation is disabled (`ext=.false.`).
        !   Default is fp(1).
    real (PREC), intent(in), optional :: right
        !*  Value to return if x > xp(-1) and extrapolation disabled (`ext=.false.`).
        !   Default is fp(-1).
    type (search_cache), intent(inout), optional :: cache
    type (status_t), intent(out), optional :: status

    integer (NF_ENUM_KIND) :: lext
    type (status_t) :: lstatus
    type (search_cache) :: lcache
    integer (INTSIZE) :: ilbound
    real (PREC) :: weight
    real (PREC) :: lleft, lright

    lstatus = NF_STATUS_OK

    call interp_linear_check_input (xp=xp, fp=fp, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_input_ext (1.0_PREC, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext
    if (present(cache)) lcache = cache

    ! The evaluator requires constants even when the selected mode ignores them.
    lleft = 0.0_PREC
    lright = 0.0_PREC
    if (present(left)) lleft = left
    if (present(right)) lright = right

    ! Find index of bracketing interval lower bound and its weight
    call interp_find_impl (x, xp, ilbound, weight, lext, lcache, lstatus)

    if (lstatus /= NF_STATUS_OK) goto 100

    ! Evaluate interpolant and desired location
    call interp_linear_eval_impl (ilbound, weight, fp, fx, lext, lleft, lright)

    if (present(cache)) cache = lcache

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(interp_linear_1d,__PREC) (x, xp, fp, fx, ext, left, &
        right, cache, status)
    !*  INTERP_LINEAR_1D performs linear interpolation for an array of
    !   points.

    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  The x-coordinates of the interpolated values.
    real (PREC), intent(in), dimension(:), contiguous :: xp
        !*  The x-coordinates of data points in increasing order.
    real (PREC), intent(in), dimension(:), contiguous :: fp
        !*  The y-coordinates of data points, same length as xp.
    real (PREC), intent(out), dimension(:), contiguous :: fx
        !*  Array of interpolated values, same shape as x.
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Defines behavious in case function should be evaluated at point x
        !   outside of the range specified by [xp(1), xp(size(xp))].
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'left' for all
        !           x < xp(1) and the value of 'right' for all x > xp(size(xp)).
        !       - NF_INTERP_EVAL_BOUNDARY: Returns fp(1) for all x < xp(1) and
        !           fp(size(xp)) for all x > xp(size(xp)).
        !           If 'left' or 'right' are omitted, values fp(1) or fp(size(xp))
        !           are used, respectively.
        !       - NF_INTERP_EVAL_ERROR: Routine exists with an error code if
        !           a value outside of the function's domain is encountered.
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in), optional :: left
        !*  Value to return if x < xp(1) and extrapolation is disabled (`ext=.false.`).
        !   Default is fp(1)
    real (PREC), intent(in), optional :: right
        !*  Value to return if x > xp(-1) and extrapolation disabled (`ext=.false.`).
        !   Default is fp(-1)
    type (search_cache), intent(inout), optional :: cache
    type (status_t), intent(out), optional :: status

    integer (NF_ENUM_KIND) :: lext
    type (status_t) :: lstatus
    type (search_cache) :: lcache
    integer (INTSIZE), dimension(:), allocatable :: ilbound
    real (PREC), dimension(:), allocatable :: weight
    integer (INTSIZE) :: n
    real (PREC) :: lleft, lright

    lstatus = NF_STATUS_OK

    call interp_linear_check_input (x, xp, fp, fx, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_input_ext (1.0_PREC, ext, left, right, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext
    if (present(cache)) lcache = cache

    ! The evaluator requires constants even when the selected mode ignores them.
    lleft = 0.0_PREC
    lright = 0.0_PREC
    if (present(left)) lleft = left
    if (present(right)) lright = right

    ! Allocate arrays to store location and interpolation weights
    n = size(x)
    allocate (ilbound(n))
    allocate (weight(n))

    ! Find index of bracketing interval lower bound and its weight
    call interp_find_impl (x, xp, ilbound, weight, lext, lcache, lstatus)

    if (lstatus /= NF_STATUS_OK) goto 100

    ! Evaluate interpolant and desired location
    call interp_linear_eval_impl (ilbound, weight, fp, fx, lext, lleft, lright)

    if (present(cache)) cache = lcache

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(interp_bilinear_impl,__PREC) (x1, x2, xp1, xp2, fp, &
        ext, fill_value, fx, wgt)
    !*  INTERP_BILINEAR_IMPL performs bilineare interpolation over a
    !   rectangular grid.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in) :: x1
    real (PREC), intent(in) :: x2
    real (PREC), intent(in), dimension(:), contiguous :: xp1
        !*  Grid points in dimension 1
    real (PREC), intent(in), dimension(:), contiguous :: xp2
        !*  Grid points in dimension 2
    real (PREC), intent(in), dimension(:,:), contiguous :: fp
        !*  Function evaluated at given grid points, ie FP(i,j) = f(xp1(i),xp2(j))
    integer (NF_ENUM_KIND), intent(in) :: ext
        !*  Defines behavious in case function should be evaluated at point (x1,x2)
        !   outside of the range specified by [xpi(1), xpi(size(xpi))] for
        !   i = 1,2.
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'fill_value'
        !           for all xi < xpi(1), i=1,2.
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in) :: fill_value
        !*  Fill value to be used if point (X1,X2) is outside of the domain
        !   defined by XP1, XP2.
    real (PREC), intent(out) :: fx
        !*  Interpolated function value f(x1,x2)
    real (PREC), intent(out), dimension(:) :: wgt
        !*  If present, returns the weights for the lower bound of the bracketing
        !   interval for each dimension.

    integer (INTSIZE) :: ilb1, ilb2
    real (PREC) :: w1, w2, fx1_lb, fx1_ub
    integer :: np1, np2

    np1 = size(xp1)
    np2 = size(xp2)

    select case (ext)
    case (NF_INTERP_EVAL_CONST)
        if (x1 < xp1(1) .or. x1 > xp1(np1) .or. x2 < xp2(1) .or. x2 > xp2(np2)) then
            wgt = 0.0_PREC
            goto 100
        end if
    end select

    ! default: At this point there is either a bracketing interval, or we extrapolate
    ! values outside of domain
    call bsearch (x1, xp1, ilb1)
    call bsearch (x2, xp2, ilb2)

    w1 = (x1 - xp1(ilb1)) / (xp1(ilb1+1) - xp1(ilb1))
    w2 = (x2 - xp2(ilb2)) / (xp2(ilb2+1) - xp2(ilb2))

    ! Interpolate in dimension 1 for upper and lower bracketing value of
    ! dimension 2.
    fx1_lb = (1.0_PREC-w1) * fp(ilb1,ilb2) + w1 * fp(ilb1+1,ilb2)
    fx1_ub = (1.0_PREC-w1) * fp(ilb1,ilb2+1) + w1 * fp(ilb1+1,ilb2+1)

    ! Interpolate in dimension two
    fx = (1.0_PREC-w2) * fx1_lb + w2 * fx1_ub

    ! Return weights on lower bounds
    wgt(1) = 1.0_PREC - w1
    wgt(2) = 1.0_PREC - w2

    return

100 continue
    ! Return constant fill_value if X in any dimension is not "interior".

    fx = fill_value

end subroutine



pure subroutine __APPEND(interp_bilinear_check_input,__PREC) (xp1, xp2, fp, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: xp1, xp2
    real (PREC), intent(in), dimension(:,:) :: fp
    type (status_t), intent(out) :: status

    status = NF_STATUS_OK

    if (size(xp1) /= size(fp, 1) .or. size(xp2) /= size(fp,2)) then
        status = NF_STATUS_INVALID_ARG
        return
    end if

    if (size(xp1) < 2 .or. size(xp2) < 2) then
        status = NF_STATUS_INVALID_ARG
        return
    end if

end subroutine



pure subroutine __APPEND(interp_bilinear_1d,__PREC) (x1, x2, xp1, xp2, fp, fx, &
        ext, fill_value, wgt, status)
    !*  INTERP_BILINEAR_IMPL performs bilineare interpolation over a
    !   rectangular grid.
    integer, parameter :: PREC = __PREC
    integer, parameter :: INTSIZE = int32

    real (PREC), intent(in), dimension(:), contiguous :: x1
        !*  Array of dimension-1 values where function should be interpolated
    real (PREC), intent(in), dimension(:), contiguous :: x2
        !*  Array of dimension-2 values where function should be interpolated
    real (PREC), intent(in), dimension(:), contiguous :: xp1
        !*  Grid points in dimension 1
    real (PREC), intent(in), dimension(:), contiguous :: xp2
        !*  Grid points in dimension 2
    real (PREC), intent(in), dimension(:,:), contiguous :: fp
        !*  Function evaluated at given grid points, ie FP(i,j) = f(xp1(i),xp2(j))
    real (PREC), intent(out), dimension(:), contiguous :: fx
        !*  Array of interpolated values
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Defines behavious in case function should be evaluated at point (x1,x2)
        !   outside of the range specified by [xpi(1), xpi(size(xpi))] for
        !   i = 1,2.
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'fill_value'
        !           for all xi < xpi(1), i=1,2.
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in), optional :: fill_value
        !*  Fill value to be used if point (X1,X2) is outside of the domain
        !   defined by XP1, XP2.
    real (PREC), intent(out), dimension(:,:), optional :: wgt
    type (status_t), intent(out), optional :: status
        !*  Optional exit status.

    type (status_t) :: lstatus
    real (PREC) :: lfill_value
    integer (NF_ENUM_KIND) :: lext
    integer (INTSIZE) :: i
    real (PREC), dimension(2) :: wgt_dummy

    lstatus = NF_STATUS_OK

    if (size(x1) /= size(x2) .or. size(x1) /= size(fx)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! Check remaining inputs using common function
    call interp_bilinear_check_input (xp1, xp2, fp, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(wgt)) then
        if (size(wgt, 1) < 2) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    lfill_value = 0.0_PREC
    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext
    if (present(fill_value)) then
        lfill_value = fill_value
        ! Assume this implies setting non-interior points to constant
        lext = NF_INTERP_EVAL_CONST
    end if

    if (present(wgt)) then
        do i = 1, size(x1)
            call interp_bilinear_impl (x1(i), x2(i), xp1, xp2, fp, lext, &
                lfill_value, fx(i), wgt(:,i))
        end do
    else
        do i = 1, size(x1)
            call interp_bilinear_impl (x1(i), x2(i), xp1, xp2, fp, lext, &
                lfill_value, fx(i), wgt_dummy)
        end do
    end if

100 continue
    if (present(status)) status = lstatus
end subroutine



pure subroutine __APPEND(interp_bilinear_scalar,__PREC) (x1, x2, xp1, xp2, fp, fx, &
        ext, fill_value, wgt, status)
    !*  INTERP_BILINEAR_IMPL performs bilineare interpolation over a
    !   rectangular grid.
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in) :: x1
        !*  Dimension-1 value where function should be interpolated
    real (PREC), intent(in) :: x2
        !*  Dimension-2 value where function should be interpolated
    real (PREC), intent(in), dimension(:), contiguous :: xp1
        !*  Grid points in dimension 1
    real (PREC), intent(in), dimension(:), contiguous :: xp2
        !*  Grid points in dimension 2
    real (PREC), intent(in), dimension(:,:), contiguous :: fp
        !*  Function evaluated at given grid points, ie FP(i,j) = f(xp1(i),xp2(j))
    real (PREC), intent(out) :: fx
        !*  Array of interpolated values
    integer (NF_ENUM_KIND), intent(in), optional :: ext
        !*  Defines behavious in case function should be evaluated at point (x1,x2)
        !   outside of the range specified by [xpi(1), xpi(size(xpi))] for
        !   i = 1,2.
        !   Adminissible constants are
        !       - NF_INTERP_EVAL_CONST: Returns value of argument 'fill_value'
        !           for all xi < xpi(1), i=1,2.
        !       - For all other values, the function value at x is extrapolated.
    real (PREC), intent(in), optional :: fill_value
        !*  Fill value to be used if point (X1,X2) is outside of the domain
        !   defined by XP1, XP2.
    real (PREC), intent(out), dimension(:), optional :: wgt
    type (status_t), intent(out), optional :: status
        !*  Optional exit status.

    type (status_t) :: lstatus
    real (PREC) :: lfill_value
    integer (NF_ENUM_KIND) :: lext
    real (PREC), dimension(2) :: wgt_dummy

    lstatus = NF_STATUS_OK

    ! Check remaining inputs using common function
    call interp_bilinear_check_input (xp1, xp2, fp, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(wgt)) then
    if (size(wgt) < 2) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    end if
    lfill_value = 0.0_PREC
    lext = NF_INTERP_EVAL_EXTRAPOLATE
    if (present(ext)) lext = ext
    if (present(fill_value)) then
        lfill_value = fill_value
        ! Assume this implies setting non-interior points to constant
        lext = NF_INTERP_EVAL_CONST
    end if

    if (present(wgt)) then
        call interp_bilinear_impl (x1, x2, xp1, xp2, fp, lext, lfill_value, &
            fx, wgt)
    else
        call interp_bilinear_impl (x1, x2, xp1, xp2, fp, lext, lfill_value, &
            fx, wgt_dummy)
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine


