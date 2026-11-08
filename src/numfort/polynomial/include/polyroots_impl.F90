

pure subroutine __APPEND(polyroots_check_input,__PREC) (coefs, roots, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: coefs
    real (PREC), intent(in), dimension(:) :: roots
    type (status_t), intent(out) :: status

    integer :: deg

    status = NF_STATUS_INVALID_ARG
    deg = size(coefs) - 1

    if (size(roots) < deg) return

    status = NF_STATUS_OK
end subroutine


pure subroutine __APPEND(polyroots,__PREC) (coefs, roots, n, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: coefs
        !*  Polynomial coefficients in increasing order, ie.
        !   the vector (a_0,a_1,...a_n) that defines the polynomial
        !       p(x) = a_n x^n + ... + a_1 x + a_0
    real (PREC), intent(out), dimension(:) :: roots
        !*  Array of (unique) real roots of a given polynomial.
    integer, intent(out) :: n
        !*  Number of (unique) real roots.
    type (status_t), intent(out), optional :: status
        !*  Status code (optional)

    type (status_t) :: lstatus
    integer :: deg

    call polyroots_check_input (coefs, roots, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    deg = size(coefs) - 1

    select case (deg)
    case (2)
        call polyroots_quad (coefs, roots, n, lstatus)
        goto 100
    case default
        status = NF_STATUS_UNSUPPORTED_OP
    end select


100 continue
    if (present(status)) status = lstatus
end subroutine



pure subroutine __APPEND(polyroots_quad,__PREC) (coefs, roots, n, status)
    !*  POLYROOT_QUAD computes the (unique) real roots of a quadratic
    !   polynomial.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: coefs
    real (PREC), intent(out), dimension(:) :: roots
    integer, intent(out) :: n
    type (status_t), intent(out) :: status

    real (PREC) :: a, b, c, d, x
    real (PREC) :: r1, r2
    real (PREC) :: scale, tol

    status = NF_STATUS_OK

    a = coefs(3)
    b = coefs(2)
    c = coefs(1)

    r1 = 0.0_PREC
    r2 = 0.0_PREC

    ! Compute according to the formula described in
    ! https://www.codeproject.com/Articles/25294/Avoiding-Overflow-Underflow-and-Loss-of-Precision
    if (abs(a) > 0.0_PREC) then
        d = b ** 2.0_PREC - 4.0_PREC * a * c
        ! For roots of multiplicity two, representation roundoff in coefficients
        ! (e.g., c = r^2) or intermediate evaluation via FMA produces a nonzero
        ! residual |d| = O(eps * scale). Use a relative tolerance around zero to
        ! prevent spurious root splitting or classifying real double roots as
        ! complex.
        scale = b ** 2.0_PREC + 4.0_PREC * abs(a * c)
        tol = 4.0_PREC * epsilon(1.0_PREC) * scale
        if (abs(d) <= tol) then
            r1 = - b / (2.0_PREC * a)
            r2 = r1
            n = 1
        else if (d > 0.0_PREC) then
            x = - b - sqrt(d)
            r1 = x / (2.0_PREC * a)
            r2 = 2.0_PREC * c / x
            ! Reorder roots so that they are in increasing order
            if (r1 > r2) then
                x = r1
                r1 = r2
                r2 = x
            end if
            n = 2
        else
            ! All roots are complex
            n = 0
            status = NF_STATUS_INVALID_STATE
        end if
    else if (abs(b) > 0.0_PREC) then
        ! Degenerate quadratic polynomial that is effectively a straight line
        r1 = - c/b
        r2 = r1
        n = 1
    else
        ! Either 0 or infinitily many roots
        n = 0
    end if

    roots(1) = r1
    roots(2) = r2

end subroutine
