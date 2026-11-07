

program test_polynomial_polyval

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common
    use numfort_common_testing
    use numfort_polynomial

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_polyder ()
    call test_polyshift ()

end subroutine

subroutine test_polyder ()

    real (PREC), parameter :: coefs(*) = [ real (PREC) &
            -1.234d0, 0.346d0, 23.345d0, -12.349d0, 2.789d0 ]
    real (PREC), dimension(:), allocatable :: coefs_new, coefs_ok
    integer :: m
    type (status_t) :: status

    ! Invalid order of differentiation
    allocate (coefs_new(1), source=1.0_PREC)
    m = -1
    call polyder (coefs(1:1), coefs_new, m, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Invalid order of differentiation"
    deallocate (coefs_new)

    ! Input array too small
    allocate (coefs_new(1))
    status = NF_STATUS_OK
    m = 0
    call polyder (coefs(1:0), coefs_new, m, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Invalid input coefficient array size"
    deallocate (coefs_new)

    ! Output array too small
    allocate (coefs_new(0))
    status = NF_STATUS_OK
    m = 0
    call polyder (coefs(1:1), coefs_new, m, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Invalid output coefficient array size"
    deallocate (coefs_new)

    ! Contant function
    allocate (coefs_new(1))
    m = 0
    status = NF_STATUS_UNDEFINED
    call polyder (coefs(1:1), coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all (coefs_new == coefs(1:1)))) &
        error stop "Polynomial degree 0, 0-th derivative"

    m = 1
    status = NF_STATUS_UNDEFINED
    call polyder (coefs(1:1), coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all (coefs_new == 0.0_PREC))) &
        error stop "Polynomial degree 0, 1-st derivative"
    deallocate (coefs_new)

    ! Linear function
    allocate (coefs_new(2))
    m = 0
    status = NF_STATUS_UNDEFINED
    call polyder (coefs(1:2), coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all (coefs(1:2) == coefs_new))) &
        error stop "Polynomial degree 1, 0-th derivative"

    m = 1
    status = NF_STATUS_UNDEFINED
    call polyder (coefs(1:2), coefs_new(1:1), m, status)
    if (.not. (status == NF_STATUS_OK .and. coefs_new(1) == coefs(2))) &
        error stop "Polynomial degree 1, 1-st derivative"

    m = 2
    status = NF_STATUS_UNDEFINED
    call polyder (coefs(1:2), coefs_new(1:1), m, status)
    if (.not. (status == NF_STATUS_OK .and. coefs_new(1) == 0.0)) &
        error stop "Polynomial degree 1, 2-nd derivative"
    deallocate (coefs_new)

    ! Quadratic function
    allocate (coefs_new(3), coefs_ok(3))
    m = 0
    status = NF_STATUS_UNDEFINED
    call polyder (coefs(1:3), coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all (coefs(1:3) == coefs_new))) &
        error stop "Polynomial degree 2, 0-th derivative"

    m = 1
    status = NF_STATUS_UNDEFINED
    coefs_ok(:) = [coefs(2), 2.0d0 * coefs(3), 0.0d0]
    call polyder (coefs(1:3), coefs_new(1:2), m, status)
    if (.not. (status == NF_STATUS_OK .and. all_close (coefs_new(1:2), coefs_ok(1:2), atol=1.d-10))) &
        error stop "Polynomial degree 2, 1-st derivative"

    m = 2
    coefs_ok(:) =  [2.0d0 * coefs(3), 0.0d0, 0.0d0]
    status = NF_STATUS_UNDEFINED
    call polyder (coefs(1:3), coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all_close (coefs_new,coefs_ok, atol=1.0d-10))) &
        error stop "Polynomial degree 2, 2-nd derivative"

    m = 3
    coefs_ok(:) =  0.0
    status = NF_STATUS_UNDEFINED
    call polyder (coefs(1:3), coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all_close (coefs_new, coefs_ok, atol=1.0d-10))) &
        error stop "Polynomial degree 2, 3-rd derivative"
    deallocate (coefs_new, coefs_ok)

    ! Cubic function
    allocate (coefs_new(4), coefs_ok(4))
    m = 0
    status = NF_STATUS_UNDEFINED
    call polyder (coefs(1:4), coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all (coefs(1:4) == coefs_new))) &
        error stop "Polynomial degree 3, 0-th derivative"

    m = 2
    status = NF_STATUS_UNDEFINED
    coefs_ok(:) = [2.0 * coefs(3), 3.0 * 2.0 * coefs(4), 0.0d0, 0.0d0]
    call polyder (coefs(1:4), coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all_close (coefs_new, coefs_ok, atol=1.0d-10))) &
        error stop "Polynomial degree 3, 2-nd derivative"

    m = 3
    status = NF_STATUS_UNDEFINED
    coefs_ok(:) = [3.0 * 2.0 * coefs(4), 0.0d0, 0.0d0, 0.0d0]
    call polyder (coefs(1:4), coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all_close (coefs_new, coefs_ok, atol=1.0d-10))) &
        error stop "Polynomial degree 3, 3-rd derivative"

    m = 4
    status = NF_STATUS_UNDEFINED
    coefs_ok(:) = 0.0
    call polyder (coefs(1:4), coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all_close (coefs_new, coefs_ok, atol=1.0d-10))) &
        error stop "Polynomial degree 3, 4-th derivative"
    deallocate (coefs_new, coefs_ok)

    ! Degree 4 polynomial
    allocate (coefs_new(5), coefs_ok(5))
    m = 2
    status = NF_STATUS_UNDEFINED
    coefs_ok(:) = [2.0 * coefs(3), 3.0*2.0*coefs(4), 4.0*3.0*coefs(5), 0.0d0, 0.0d0]
    call polyder (coefs, coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all_close (coefs_new, coefs_ok, atol=1.0d-10))) &
        error stop "Polynomial degree 4, 2-nd derivative"

    m = 3
    status = NF_STATUS_UNDEFINED
    coefs_ok(:) = [3.0*2.0*coefs(4), 4.0*3.0*2.0*coefs(5), 0.0d0, 0.0d0, 0.0d0]
    call polyder (coefs, coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all_close (coefs_new, coefs_ok, atol=1.0d-10))) &
        error stop "Polynomial degree 4, 3-rd derivative"

    m = 4
    status = NF_STATUS_UNDEFINED
    coefs_ok(:) = 0.0
    coefs_ok(1) = 4.0 * 3.0 * 2.0 * coefs(5)
    call polyder (coefs, coefs_new, m, status)
    if (.not. (status == NF_STATUS_OK .and. all_close (coefs_new, coefs_ok, atol=1.0d-10))) &
        error stop "Polynomial degree 4, 4-th derivative"

    deallocate (coefs_new)

end subroutine

subroutine test_polyshift ()
    !*  Unit tests for POLYSHIFT routine

    integer, parameter :: KMAX = 4
    real (PREC), parameter :: coefs0(0:KMAX) = [ real (PREC) :: &
        -1.23d0, 2.34d0, -5.678d0, 8.234d0, 1.789d0 ]
    real (PREC), dimension(:), allocatable :: coefs1
    real (PREC), dimension(:), allocatable :: x0, x1, y0, y1
    real (PREC) :: xshift
    integer :: k, n
    type (status_t) :: status

    ! Test input checking
    allocate (coefs1(1))
    xshift = 0.0

    call polyshift (coefs0, xshift, coefs1, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Input checks: size(COEFS0) /= size(COEFS1)"

    deallocate (coefs1)

    ! Test with no shift
    xshift = 0.0
    do k = 0, KMAX
        allocate (coefs1(k+1))

        call polyshift (coefs0(0:k), xshift, coefs1, status)

        if (.not. (status == NF_STATUS_OK .and. all_close (coefs1, coefs0(0:k), &
            atol=1.0d-10, rtol=0.0d0))) &
            error stop "POLYSHIFT with zero shift changed coefficients"

        deallocate (coefs1)
    end do

    ! Test with positive shift
    xshift = 1.345d0
    n = 21
    do k = 0, KMAX
        allocate (coefs1(k+1))
        allocate (x0(n), y0(n), x1(n), y1(n))

        call linspace (x0, -2.0d0, 2.0d0)
        x1(:) = x0 - xshift

        call polyshift (coefs0(0:k), xshift, coefs1, status)

        ! Check that polynomials evaluate to the same number
        call polyval (coefs0(0:k), x0, y0)
        call polyval (coefs1, x1, y1)

        if (.not. (status == NF_STATUS_OK .and. all_close (y0, y1, atol=1.0d-10, rtol=0.0d0))) &
            error stop "POLYSHIFT with positive shift changed polynomial values"

        deallocate (coefs1)
        deallocate (x0, x1, y0, y1)
    end do

    ! Test with negative shift
    xshift = -7.234d0
    n = 31
    do k = 0, KMAX
        allocate (coefs1(k+1))
        allocate (x0(n), y0(n), x1(n), y1(n))

        call linspace (x0, -10.0d0, 10.0d0)
        x1(:) = x0 - xshift

        call polyshift (coefs0(0:k), xshift, coefs1, status)

        ! Check that polynomials evaluate to the same number
        call polyval (coefs0(0:k), x0, y0)
        call polyval (coefs1, x1, y1)

        if (.not. (status == NF_STATUS_OK .and. all_close (y0, y1, atol=1.0d-10, rtol=0.0d0))) &
            error stop "POLYSHIFT with negative shift changed polynomial values"

        deallocate (coefs1)
        deallocate (x0, x1, y0, y1)
    end do

end subroutine

end program
