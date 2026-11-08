


pure subroutine solver_map_init (self, lb, ub, y0, dy, dx, status)
    !*  MAPPING_INIT initializes the MAPPING object representing the function
    !   y = f(x) such that it maps the X on the real line into one of the
    !   intervals
    !       (1) [-inf, inf]
    !       (2) [y_lb, inf)
    !       (3) [-inf, y_ub] (not yet implemented)
    !       (4) [y_lb, y_ub]
    !
    !   Additionally, a scaling parameter can be computed that ensures
    !   a unit change in x results in a "reasonable" step size in y
    !   around x0 = f^{-1}(y0)
    !
    !   Linear map
    !   ==========
    !   For the mapping of type (1) we apply the transformation
    !       y = x/s
    !   and thus the scale parameter is determined as
    !       s = dx/dy
    !
    !   Exponential map
    !   ===============
    !   For mappings of type (2) we apply the transformation
    !       y = exp(x/s) + y_lb
    !   To pin down the scale parameter s, we use
    !       dy = exp((x0+dx)/s) - exp(x0/s)
    !          = exp(x0/s)*exp(dx/s) - exp(x0/s)
    !   where
    !       x0 = s * log(y0 - y_lb)
    !   and therefore
    !       dy = (y0 - y_lb) * exp(dx/s) - (y0 - y_lb)
    !       dx/s = log(y0 - y_lb + dy) - log(y0 - y_lb)
    !
    !   Logistic map
    !   ============
    !   For bounded mappings of type (4), we apply the transformation
    !       y = y_lb + F(x,s)(y_u - y_lb)
    !   where [y_lb, y_ub] are the bounds to be enforced for y and
    !   F(x,s) is the logistic CDF with the yet-to-be determined scale
    !   coefficient s,
    !       F(x,s) = 1/(1 + exp(-x/s))
    !
    !   We require that at y0
    !       dy = f(x0+dx) - f(x0)
    !   where dy is a "reasonable" step size for y at x0.
    !   Thus
    !       dy  = y_lb + F(x0+dx,s)(y_ub - y_lb) - y_lb - F(x0,s)(y_ub - y_lb)
    !           = (F(x0+dx,s)-F(x0,s))(y_ub - y_lb)
    !
    !   Additionally, the initial point y0 is given by
    !       (y0 - y_lb)/(y_ub - y_lb) = F(x0,s)
    !   so from the definition of the logistic CDF we get
    !       x0 = - s * log(1/z - 1)
    !   where z = (y0 - y_lb)/(y_ub - y_lb)
    !
    !   Combining these expressions, we find that s must satisfy
    !       dy/(y_ub - y_lb) = 1/(1 + (1/z-1)exp(-dx/s)) - z
    !   which implies that
    !       dx/s = - [log(y_ub - y0 - dy) - log(y0 - y_lb + dy)] +
    !               [log(y_ub - y0) - log(y0 - y_lb)]
    type (solver_map), intent(inout) :: self
    real (PREC), intent(in) :: lb
    real (PREC), intent(in), optional :: ub
    real (PREC), intent(in), optional :: y0
    real (PREC), intent(in), optional :: dx
    real (PREC), intent(in), optional :: dy
    type (status_t), intent(out), optional :: status

    real (PREC) :: inf, ninf
    integer :: transform
    real (PREC) :: ylb, yub, a1, a2, s, ldx
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    inf = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    ninf = ieee_value (0.0_PREC, IEEE_NEGATIVE_INF)

    ! Default range boundaries
    ylb = ninf
    yub = inf
    ! Default scale parameter
    s = 1.0_PREC
    ! Default transformation
    transform = TRANSFORM_LINEAR
    ldx = 1.0_PREC

    ylb = lb
    if (present(ub)) yub = ub
    if (present(dx)) ldx = dx

    ! Input validation
    call check_nonzero (0.0_PREC, dy, 'dy', status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, dx, 'dx', status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (ieee_class (ylb) == IEEE_NEGATIVE_INF &
            .and. ieee_class (yub) == IEEE_POSITIVE_INF) then
        ! Maps from R -> R, possibly with rescaling
        transform = TRANSFORM_LINEAR
        if (present(dy)) then
            s = ldx / dy
        end if

    else if (ylb > ninf .and. ieee_class (yub) == IEEE_POSITIVE_INF) then
        ! Maps R -> [lb, inf); use exp-like transformation
        if (present(y0) .and. present(dy)) then
            ! Ensure that y0 in (y_lb, inf)
            if (y0 <= ylb .or. (y0 + dy) <= ylb) then
                lstatus = NF_STATUS_INVALID_ARG
                lstatus = lstatus + NF_STATUS_DOMAIN_ERROR
                goto 100
            end if

            s = ldx / (log(y0 + dy - ylb) - log(y0 - ylb))
        end if

        transform = TRANSFORM_EXP

    else if (ylb > ninf .and. yub < inf) then
        ! Maps R -> [lb,ub]; use logit CDF transformation
        if (present(y0) .and. present(dy)) then
            ! Ensure that y0 in (y_lb, y_ub)
            if (y0 <= ylb .or. yub <= y0 .or. (y0 + dy) <= ylb .or. &
                    (y0 + dy) >= yub) then
                lstatus = NF_STATUS_INVALID_ARG
                lstatus = lstatus + NF_STATUS_DOMAIN_ERROR
                goto 100
            end if

            a1 = log(y0 + dy - ylb) - log(yub - y0 - dy)
            a2 = log(yub - y0) - log(y0 - ylb)

            s = ldx/(a1 + a2)

            if (.not. ieee_is_finite (s)) then
                lstatus = NF_STATUS_INVALID_STATE
                goto 100
            end if
        end if

        transform = TRANSFORM_LOGISTIC
    end if

    self%transform = transform
    self%scale = s
    self%lb = ylb
    self%ub = yub

100 continue
    if (present(status)) status = lstatus

end subroutine


pure subroutine safe_div (num, den, quot, overflows)
    !*  SAFE_DIV computes QUOT = NUM / DEN for finite operands without raising
    !   IEEE divide-by-zero or overflow exception flags.
    real (PREC), intent(in) :: num
        !*  Finite numerator.
    real (PREC), intent(in) :: den
        !*  Finite denominator.
    real (PREC), intent(out) :: quot
        !*  Quotient, or signed infinity when OVERFLOWS is true.
    logical, intent(out) :: overflows
        !*  True if DEN == 0 or |NUM / DEN| exceeds huge(0.0_PREC).

    real (PREC) :: ad, an, dh, qh, sgn

    ad = abs (den)
    an = abs (num)
    sgn = sign (1.0_PREC, num) * sign (1.0_PREC, den)

    if (ad == 0.0_PREC) then
        ! Zero denominator: 0/0 -> quiet NaN, nonzero/0 -> signed Inf.
        if (an == 0.0_PREC) then
            quot = ieee_value (0.0_PREC, IEEE_QUIET_NAN)
        else
            quot = ieee_copy_sign ( &
                ieee_value (0.0_PREC, IEEE_POSITIVE_INF), sgn)
        end if
        overflows = .true.
    else if (ad >= 1.0_PREC) then
        ! |den| >= 1 shrinks or preserves |num|, so NUM / DEN cannot overflow.
        quot = num / den
        overflows = .false.
    else
        ! 0 < |den| < 1: even for the smallest subnormal |den|, dh = |den| * huge
        ! is strictly normal (>= 2^-51 in binary64, >= 2^-22 in binary32).
        dh = ad * huge (0.0_PREC)
        if (an <= 0.5_PREC * dh) then
            ! |num / den| <= 0.5 * huge: direct division cannot overflow.
            quot = num / den
            overflows = .false.
        else if (0.5_PREC * an > 0.75_PREC * dh) then
            ! |num / den| > 1.5 * huge: guaranteed overflow; compare at half/3/4
            ! scale so neither side overflows huge.
            quot = ieee_copy_sign ( &
                ieee_value (0.0_PREC, IEEE_POSITIVE_INF), sgn)
            overflows = .true.
        else
            ! Boundary band 0.5 * dh < |num| <= 1.5 * dh: 0.5 * num is normal and
            ! exact (exponent - 1), so qh = (0.5 * num) / den has |qh| <= 0.75 * huge
            ! with the exact rounded significand of NUM / DEN.
            qh = (0.5_PREC * num) / den
            if (abs (qh) > 0.5_PREC * huge (0.0_PREC)) then
                quot = ieee_copy_sign ( &
                    ieee_value (0.0_PREC, IEEE_POSITIVE_INF), sgn)
                overflows = .true.
            else
                quot = 2.0_PREC * qh
                overflows = .false.
            end if
        end if
    end if

end subroutine


pure subroutine safe_mul (a, b, prod, overflows)
    !*  SAFE_MUL computes PROD = A * B for finite operands without raising the
    !   IEEE overflow exception flag.
    real (PREC), intent(in) :: a
        !*  First finite factor.
    real (PREC), intent(in) :: b
        !*  Second finite factor.
    real (PREC), intent(out) :: prod
        !*  Product, or signed infinity when OVERFLOWS is true.
    logical, intent(out) :: overflows
        !*  True if |A * B| exceeds huge(0.0_PREC).

    real (PREC) :: aa, ab, lim, ph, sgn

    aa = abs (a)
    ab = abs (b)

    ! Guard |a| <= 1 and |b| <= 1 in separate branches so huge / |b| below is
    ! never evaluated when |b| <= 1 (Fortran .and. need not short-circuit).
    if (aa <= 1.0_PREC) then
        prod = a * b
        overflows = .false.
    else if (ab <= 1.0_PREC) then
        prod = a * b
        overflows = .false.
    else
        ! Both |a| > 1 and |b| > 1, so 1 <= lim = huge / |b| < huge is finite.
        sgn = sign (1.0_PREC, a) * sign (1.0_PREC, b)
        lim = huge (0.0_PREC) / ab
        if (aa <= 0.5_PREC * lim) then
            prod = a * b
            overflows = .false.
        else if (0.5_PREC * aa > 0.75_PREC * lim) then
            prod = ieee_copy_sign ( &
                ieee_value (0.0_PREC, IEEE_POSITIVE_INF), sgn)
            overflows = .true.
        else
            ! Boundary band 0.5 * lim < |a| <= 1.5 * lim: 0.5 * a is exact, so
            ! ph = (0.5 * a) * b cannot overflow and tests exact IEEE rounding.
            ph = (0.5_PREC * a) * b
            if (abs (ph) > 0.5_PREC * huge (0.0_PREC)) then
                prod = ieee_copy_sign ( &
                    ieee_value (0.0_PREC, IEEE_POSITIVE_INF), sgn)
                overflows = .true.
            else
                prod = 2.0_PREC * ph
                overflows = .false.
            end if
        end if
    end if

end subroutine


pure subroutine solver_map_eval_scalar (self, x, y, jac)
    !*  SOLVER_MAP_EVAL_SCALAR maps an unconstrained solver coordinate X to the
    !   bounded parameter domain and optionally computes dy/dx.
    type (solver_map), intent(in) :: self
        !*  Initialized mapping object.
    real (PREC), intent(in) :: x
        !*  Unconstrained solver coordinate.
    real (PREC), intent(out) :: y
        !*  Mapped parameter value.
    real (PREC), intent(out), optional :: jac
        !*  Optional derivative dy/dx.

    real (PREC) :: s, lb, ub, e, u, q, span, half_span, w
    logical :: wide_span, ovf

    s = real (self%scale, PREC)
    lb = real (self%lb, PREC)
    ub = real (self%ub, PREC)

    y = 0.0_PREC
    if (present (jac)) jac = 0.0_PREC

    select case (self%transform)
    case (TRANSFORM_LINEAR)
        y = x / s
        if (present (jac)) then
            jac = 1.0_PREC / s
        end if

    case (TRANSFORM_EXP)
        e = exp (x / s)
        y = e + lb
        if (present (jac)) then
            jac = e / s
        end if

    case (TRANSFORM_NEG_EXP)
        continue

    case (TRANSFORM_LOGISTIC)
        if (ieee_is_nan (x)) then
            y = x
            if (present (jac)) jac = x
        else
            ! Evaluate q = expit(-|u|) in [0, 0.5] with u = x / s, and anchor y
            ! to the nearer bound: y = lb + q * (ub - lb) for u <= 0 and
            ! y = ub - q * (ub - lb) for u > 0. This avoids exp overflow for
            ! large |u| and preserves accuracy near ub when |ub| << |lb|.
            u = x / s
            q = expit (-abs (u))

            ! ub - lb can exceed huge only when lb < 0 < ub; check sequentially
            ! so huge + lb never overflows when lb >= 0.
            wide_span = .false.
            if (lb < 0.0_PREC) then
                wide_span = (ub > huge (0.0_PREC) + lb)
            end if

            if (wide_span) then
                ! Factor ub - lb = 2 * half_span with half_span <= huge and
                ! 2 * q <= 1 so neither factor overflows.
                half_span = 0.5_PREC * ub - 0.5_PREC * lb
                if (u <= 0.0_PREC) then
                    y = lb + (2.0_PREC * q) * half_span
                else
                    y = ub - (2.0_PREC * q) * half_span
                end if
            else
                span = ub - lb
                if (u <= 0.0_PREC) then
                    y = lb + q * span
                else
                    y = ub - q * span
                end if
            end if

            if (present (jac)) then
                ! dy/dx = ((ub - lb) / s) * w with w = q * (1 - q) in [0, 0.25].
                ! Prefer dividing span by s first so subnormal w does not lose
                ! precision before division; fall back to safe_div((span * w), s)
                ! when small |s| would overflow span / s. Check |s| >= 1 in its
                ! own branch so |s| * huge is only evaluated for |s| < 1.
                w = q * (1.0_PREC - q)
                if (wide_span) then
                    if (abs (s) >= 2.0_PREC) then
                        jac = (half_span / (0.5_PREC * s)) * w
                    else
                        call safe_div (half_span * (2.0_PREC * w), s, jac, ovf)
                    end if
                else if (abs (s) >= 1.0_PREC) then
                    jac = (span / s) * w
                else if (span <= 0.5_PREC * (abs (s) * huge (0.0_PREC))) then
                    jac = (span / s) * w
                else
                    call safe_div (span * w, s, jac, ovf)
                end if
            end if
        end if
    end select

end subroutine


subroutine solver_map_eval_diag (map, x, y, jac, status)
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: y
    real (PREC), intent(out), dimension(:), optional :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: i, n

    lstatus = NF_STATUS_OK

    n = size(map)

    if (size(y) /= n .or. size(x) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (present(jac)) then
        if (size(jac) /= n) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    ! Forward optional JAC only when requested by the caller.
    if (present(jac)) then
        do i = 1, n
            call solver_map_eval (map(i), x(i), y(i), jac(i))
        end do
    else
        do i = 1, n
            call solver_map_eval (map(i), x(i), y(i))
        end do
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine


subroutine solver_map_eval_matrix (map, x, y, jac, status)
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: y
    real (PREC), intent(out), dimension(:,:) :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: i, n

    lstatus = NF_STATUS_OK

    n = size(map)

    if (size(y) /= n .or. size(x) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (size (jac, 1) /= n .or. size (jac, 2) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    jac = 0.0_PREC
    do i = 1, n
        call solver_map_eval (map(i), x(i), y(i), jac(i,i))
    end do

100 continue
    if (present(status)) status = lstatus

end subroutine


pure subroutine solver_map_eval_inverse_scalar (self, y, x, jac, status)
    !*  SOLVER_MAP_EVAL_INVERSE_SCALAR maps a bounded parameter Y to its
    !   unconstrained solver coordinate X and optionally computes dx/dy.
    type (solver_map), intent(in) :: self
        !*  Initialized mapping object.
    real (PREC), intent(in) :: y
        !*  Parameter value in the target domain.
    real (PREC), intent(out) :: x
        !*  Unconstrained solver coordinate.
    real (PREC), intent(out), optional :: jac
        !*  Optional derivative dx/dy.
    type (status_t), intent(out), optional :: status
        !*  Optional exit status.

    type (status_t) :: lstatus
    real (PREC) :: s, lb, ub, inf, qnan, dl, du, h_dl, h_du, u, t1, t2
    logical :: dl_ovf, du_ovf, mid_band, ovf, ovf1, ovf2

    lstatus = NF_STATUS_OK
    s = real (self%scale, PREC)
    lb = real (self%lb, PREC)
    ub = real (self%ub, PREC)

    inf = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    qnan = ieee_value (0.0_PREC, IEEE_QUIET_NAN)

    x = 0.0_PREC
    if (present (jac)) jac = 0.0_PREC

    select case (self%transform)
    case (TRANSFORM_LINEAR)
        x = y * s
        if (present (jac)) then
            jac = s
        end if

    case (TRANSFORM_EXP)
        ! Check NaN in its own branch before relational `<` so quiet NaNs do
        ! not signal IEEE_INVALID when `.or.` is not short-circuited.
        if (ieee_is_nan (y)) then
            x = qnan
            if (present (jac)) jac = qnan
            lstatus = NF_STATUS_INVALID_ARG + NF_STATUS_DOMAIN_ERROR
            goto 100
        else if (y < lb) then
            x = qnan
            if (present (jac)) jac = qnan
            lstatus = NF_STATUS_INVALID_ARG + NF_STATUS_DOMAIN_ERROR
            goto 100
        end if

        if (y == lb) then
            ! Exact lower bound maps to -sign(s) * Inf; dx/dy = s / (y - lb)
            ! diverges and is only flagged if JAC was explicitly requested.
            x = ieee_copy_sign (inf, -s)
            if (present (jac)) then
                jac = ieee_copy_sign (inf, s)
                lstatus = NF_STATUS_INVALID_STATE
            end if
        else if (.not. ieee_is_finite (y)) then
            ! Upper limit y = +Inf maps to +sign(s) * Inf with dx/dy -> 0.
            x = ieee_copy_sign (inf, s)
            if (present (jac)) then
                jac = ieee_copy_sign (0.0_PREC, s)
            end if
        else
            ! Finite interior y > lb: y - lb can exceed huge only when y > 1
            ! and lb < -1; halve both operands to evaluate log(y - lb) safely.
            dl_ovf = .false.
            if (y > 1.0_PREC .and. lb < -1.0_PREC) then
                dl_ovf = (0.5_PREC * y - 0.5_PREC * lb &
                    > 0.5_PREC * huge (0.0_PREC))
            end if

            if (dl_ovf) then
                h_dl = 0.5_PREC * y - 0.5_PREC * lb
                u = log (h_dl) + log (2.0_PREC)
            else
                dl = y - lb
                u = log (dl)
            end if

            call safe_mul (s, u, x, ovf)
            if (ovf) lstatus = NF_STATUS_INVALID_STATE

            if (present (jac)) then
                if (dl_ovf) then
                    ! Since y - lb > huge, |s / (y - lb)| < |s| / huge <= 1;
                    ! compute via halved numerator/denominator or underflow to 0.
                    if (abs (s) >= 2.0_PREC * tiny (0.0_PREC)) then
                        jac = (0.5_PREC * s) / h_dl
                    else
                        jac = ieee_copy_sign (0.0_PREC, s)
                    end if
                else
                    call safe_div (s, dl, jac, ovf)
                    if (ovf) lstatus = NF_STATUS_INVALID_STATE
                end if
            end if
        end if

    case (TRANSFORM_NEG_EXP)
        continue

    case (TRANSFORM_LOGISTIC)
        ! Check NaN separately so `<` / `>` never signal IEEE_INVALID on NaN.
        if (ieee_is_nan (y)) then
            x = qnan
            if (present (jac)) jac = qnan
            lstatus = NF_STATUS_INVALID_ARG + NF_STATUS_DOMAIN_ERROR
            goto 100
        else if (y < lb .or. y > ub) then
            x = qnan
            if (present (jac)) jac = qnan
            lstatus = NF_STATUS_INVALID_ARG + NF_STATUS_DOMAIN_ERROR
            goto 100
        end if

        ! Handle collapsed intervals and exact bounds before forming gaps.
        ! At y == lb or y == ub, x = -/+sign(s) * Inf is valid when JAC is
        ! omitted, whereas dx/dy diverges (NF_STATUS_INVALID_STATE with JAC).
        if (lb == ub) then
            x = qnan
            if (present (jac)) jac = qnan
            lstatus = NF_STATUS_INVALID_STATE
            goto 100
        else if (y == lb) then
            x = ieee_copy_sign (inf, -s)
            if (present (jac)) then
                jac = ieee_copy_sign (inf, s)
                lstatus = NF_STATUS_INVALID_STATE
            end if
            goto 100
        else if (y == ub) then
            x = ieee_copy_sign (inf, s)
            if (present (jac)) then
                jac = ieee_copy_sign (inf, s)
                lstatus = NF_STATUS_INVALID_STATE
            end if
            goto 100
        end if

        ! Work directly with endpoint gaps dl ~ (y - lb) and du ~ (ub - y)
        ! rather than normalized z = (y - lb) / (ub - lb): normalizing first
        ! rounds z to 1 when |ub| << |lb|, and forming 1 / z or dl / du before
        ! log overflows whenever a gap is smaller than 1 / huge (e.g. subnormal).
        !
        ! Because -huge <= lb < y < ub <= huge, y - lb can exceed huge only when
        ! y > 1 and lb < -1 (where ub - y >= eps(1), so 0.5 * (ub - y) cannot
        ! underflow), and symmetrically ub - y can exceed huge only when ub > 1
        ! and y < -1. Scaling both gaps by 0.5 cancels out in u = log(dl / du).
        dl_ovf = .false.
        if (y > 1.0_PREC .and. lb < -1.0_PREC) then
            dl_ovf = (0.5_PREC * y - 0.5_PREC * lb &
                > 0.5_PREC * huge (0.0_PREC))
        end if

        du_ovf = .false.
        if (ub > 1.0_PREC .and. y < -1.0_PREC) then
            du_ovf = (0.5_PREC * ub - 0.5_PREC * y &
                > 0.5_PREC * huge (0.0_PREC))
        end if

        if (dl_ovf) then
            dl = 0.5_PREC * y - 0.5_PREC * lb
            du = 0.5_PREC * (ub - y)
        else if (du_ovf) then
            dl = 0.5_PREC * (y - lb)
            du = 0.5_PREC * ub - 0.5_PREC * y
        else
            dl = y - lb
            du = ub - y
        end if

        ! Test midpoint band 0.5 * du <= dl <= 2 * du via subtraction of the
        ! smaller gap from the larger: by Sterbenz's lemma this is exact near
        ! the band edges and avoids 0.5 * du underflow and 2 * du overflow.
        if (dl <= du) then
            mid_band = (du - dl <= dl)
        else
            mid_band = (dl - du <= du)
        end if

        if (mid_band) then
            ! Near the midpoint y = (lb + ub) / 2, log(dl) - log(du) suffers
            ! cancellation; use log(dl / du) = 2 * atanh((dl - du) / (dl + du)),
            ! where dl - du is Sterbenz-exact. Halve dl and du if dl + du > huge.
            if (dl > huge (0.0_PREC) - du) then
                h_dl = 0.5_PREC * dl
                h_du = 0.5_PREC * du
                u = 2.0_PREC * atanh ((h_dl - h_du) / (h_dl + h_du))
            else
                u = 2.0_PREC * atanh ((dl - du) / (dl + du))
            end if
        else
            ! Outside [0.5, 2], log(dl) and log(du) are well separated, and
            ! subtracting logs avoids dl / du overflow or underflow in the tails.
            u = log (dl) - log (du)
        end if

        call safe_mul (s, u, x, ovf)
        if (ovf) lstatus = NF_STATUS_INVALID_STATE

        if (present (jac)) then
            ! Evaluate dx/dy = s / (y - lb) + s / (ub - y) as the sum of two
            ! same-sign terms, avoiding (y - lb)^2 underflow near the bounds.
            if (dl_ovf) then
                if (abs (s) >= 2.0_PREC * tiny (0.0_PREC)) then
                    t1 = (0.5_PREC * s) / dl
                else
                    t1 = ieee_copy_sign (0.0_PREC, s)
                end if
                ovf1 = .false.
                call safe_div (s, ub - y, t2, ovf2)
            else if (du_ovf) then
                call safe_div (s, y - lb, t1, ovf1)
                if (abs (s) >= 2.0_PREC * tiny (0.0_PREC)) then
                    t2 = (0.5_PREC * s) / du
                else
                    t2 = ieee_copy_sign (0.0_PREC, s)
                end if
                ovf2 = .false.
            else
                call safe_div (s, dl, t1, ovf1)
                call safe_div (s, du, t2, ovf2)
            end if

            ! Check individual term overflow and sum overflow |t1| + |t2| > huge.
            if (ovf1 .or. ovf2) then
                jac = ieee_copy_sign (inf, s)
                lstatus = NF_STATUS_INVALID_STATE
            else if (abs (t1) > huge (0.0_PREC) - abs (t2)) then
                jac = ieee_copy_sign (inf, s)
                lstatus = NF_STATUS_INVALID_STATE
            else
                jac = t1 + t2
            end if
        end if
    end select

100 continue
    if (present (status)) status = lstatus

end subroutine


subroutine solver_map_eval_inverse_diag (map, y, x, jac, status)
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: y
    real (PREC), intent(out), dimension(:) :: x
    real (PREC), intent(out), dimension(:), optional :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: i, n

    lstatus = NF_STATUS_OK

    n = size(map)

    if (size(y) /= n .or. size(x) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (present(jac)) then
        if (size(jac) /= n) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    ! Forward JAC only when present so no-Jacobian calls do not fail with
    ! NF_STATUS_INVALID_STATE when dx/dy overflows near or at the bounds.
    if (present(jac)) then
        do i = 1, n
            call solver_map_eval_inverse (map(i), y(i), x(i), jac(i), lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100
        end do
    else
        do i = 1, n
            call solver_map_eval_inverse (map(i), y(i), x(i), status=lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100
        end do
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine


subroutine solver_map_eval_inverse_matrix (map, y, x, jac, status)
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: y
    real (PREC), intent(out), dimension(:) :: x
    real (PREC), intent(out), dimension(:,:) :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: i, n

    lstatus = NF_STATUS_OK

    n = size(map)

    if (size(y) /= n .or. size(x) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (size (jac, 1) /= n .or. size (jac, 2) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    jac = 0.0_PREC

    do i = 1, n
        call solver_map_eval_inverse (map(i), y(i), x(i), jac(i,i), lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end do


100 continue
    if (present(status)) status = lstatus

end subroutine


