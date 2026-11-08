


impure elemental function __APPEND(logaddexp,__PREC) (x, y) result(res)
    !*  LOGADDEXP implements the function log(exp(x) + exp(y)) that should
    !   perform better for inputs negative X and Y that are sufficiently
    !   large (in absolute terms) such that exp(X) or exp(Y) evaluates
    !   to zero.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: x
    real (PREC), intent(in) :: y
    real (PREC) :: res

    real (PREC) :: tmp

    real (PREC), parameter :: LOGE2 = log(2.0_PREC)

    if (x == y) then
        ! Handles infinities of the same sign without warnings
        res = x + LOGE2
    else
        tmp = x - y
        if (tmp > 0.0_PREC) then
            res = x + log1p(exp(-tmp))
        else if (tmp <= 0.0_PREC) then
            res = y + log1p(exp(tmp))
        else
            ! NaNs
            res = tmp
        end if
    end if

end function


pure elemental function __APPEND(expit,__PREC) (x) result(res)
    !*  EXPIT evaluates the standard logistic sigmoid function
    !   1 / (1 + exp(-x)) without intermediate exponential overflow.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: x
        !*  Input coordinate on the extended real line.
    real (PREC) :: res
        !*  Logistic function value in [0, 1].

    real (PREC) :: e

    if (ieee_is_nan (x)) then
        res = x
    else if (x < 0.0_PREC) then
        e = exp (x)
        res = e / (1.0_PREC + e)
    else
        e = exp (-x)
        res = 1.0_PREC / (1.0_PREC + e)
    end if

end function


pure elemental function __APPEND(logit,__PREC) (x) result(res)
    !*  LOGIT evaluates the log-odds function log(x / (1 - x)) on [0, 1]
    !   without reciprocal overflow near 0 or cancellation near 0.5.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: x
        !*  Input probability in [0, 1].
    real (PREC) :: res
        !*  Log-odds coordinate in [-Inf, +Inf], or quiet NaN outside [0, 1].

    if (ieee_is_nan (x)) then
        res = ieee_value (0.0_PREC, IEEE_QUIET_NAN)
    else if (x < 0.0_PREC .or. x > 1.0_PREC) then
        res = ieee_value (0.0_PREC, IEEE_QUIET_NAN)
    else if (x == 0.0_PREC) then
        res = ieee_value (0.0_PREC, IEEE_NEGATIVE_INF)
    else if (x == 1.0_PREC) then
        res = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    else if (x >= 0.25_PREC .and. x <= 0.75_PREC) then
        res = 2.0_PREC * atanh (2.0_PREC * x - 1.0_PREC)
    else
        res = log (x / (1.0_PREC - x))
    end if

end function
