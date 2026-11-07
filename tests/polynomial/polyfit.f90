

program test_polynomial_polyfit

    use, intrinsic :: iso_fortran_env
    use numfort_common
    use numfort_common_workspace, workspace => workspace_real64
    use numfort_common_testing
    use numfort_polynomial

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_polyfit ()
    call test_polyfit_deriv ()

end subroutine

subroutine test_polyfit ()

    real (PREC), dimension(:,:), allocatable :: coefs, coefs2, y
    real (PREC), dimension(:), allocatable :: x
    type (workspace) :: work
    type (status_t) :: status
    real (PREC) :: diff
    real (PREC), parameter :: tol=1.0d-10

    integer :: i, deg

    ! Quadratic polynomials
    deg = 2
    allocate (coefs(deg+1,2), coefs2(deg+1,2), y(deg+1,2))
    allocate (x(deg+1))

    coefs(:,1) = [-0.233d0,  1.2354d0, 7.0234d0]
    coefs(:,2) = [ 1.234d0, -7.2389d0, 1.0234d0]

    x(:) = [-1.234d0, 2.234d0, 5.279d0]

    do i = 1, 2
        call polyval (coefs(:,i), x, y(:,i))
    end do

    call polyfit (x, y, deg=deg, coefs=coefs2, work=work, status=status)
    diff = maxval(abs(coefs-coefs2))

    if (.not. (diff < tol .and. status == NF_STATUS_OK)) &
        error stop "Fit polynomials of degree 2"

end subroutine

subroutine test_polyfit_deriv ()

    real (PREC), dimension(:), allocatable :: y, coefs, coefs2, coefs_deriv
    real (PREC) :: x
    type (workspace) :: work
    integer :: deg, n, k
    real (PREC), parameter :: xtol=1.0d-10

    type (status_t) :: status

    ! Quadratic function with all non-zero coefs
    deg = 2
    n = deg + 1
    allocate (y(n), coefs(n), coefs2(n))

    coefs(:) = [-1.234d0, 2.688d0, 3.132d0]
    x = -6.234d0

    allocate (coefs_deriv(deg+1))
    do k = 0, deg
        call polyder (coefs, coefs_deriv(:deg+1-k), k, status)
        call polyval (coefs_deriv(:deg+1-k), x, y(k+1), status)
    end do
    deallocate (coefs_deriv)

    call polyfit_deriv (x, y, coefs2, work, status)
    if (.not. (all_close (coefs, coefs2, atol=xtol))) &
        error stop "Polynomial of degree 2, all non-zero coefs"
    deallocate (y, coefs, coefs2)

    ! Quadratic function with some zero coefs
    deg = 2
    n = deg + 1
    allocate (y(n), coefs(n), coefs2(n))
    coefs(:) = [-123.23d0, 0.0d0, 34.89d0]
    x = 0.1234d0

    allocate (coefs_deriv(deg+1))
    do k = 0, deg
        call polyder (coefs, coefs_deriv(:deg+1-k), k, status)
        call polyval (coefs_deriv(:deg+1-k), x, y(k+1), status)
    end do
    deallocate (coefs_deriv)

    call polyfit_deriv (x, y, coefs2, work, status)
    if (.not. (all_close (coefs, coefs2, atol=xtol))) &
        error stop "Polynomial of degree 2, some zero coefs"
    deallocate (y, coefs, coefs2)

    ! Quadratic function with highest coefficient being zero
    deg = 2
    n = deg + 1
    allocate (y(n), coefs(n), coefs2(n))
    coefs(:) = [-123.23d0, 234.9d0, 0.0d0]
    x = 10.2345d0

    allocate (coefs_deriv(deg+1))
    do k = 0, deg
        call polyder (coefs, coefs_deriv(:deg+1-k), k, status)
        call polyval (coefs_deriv(:deg+1-k), x, y(k+1), status)
    end do
    deallocate (coefs_deriv)

    call polyfit_deriv (x, y, coefs2, work, status)
    if (.not. (all_close (coefs, coefs2, atol=xtol))) &
        error stop "Polynomial of degree 2, some zero coefs"
    deallocate (y, coefs, coefs2)

    ! Cubic polynomial
    deg = 3
    n = deg + 1
    allocate (y(n), coefs(n), coefs2(n))
    coefs(:) = [-123.23d0, 0.0d0, -34.9d0, 3.456d0]
    x = 1.2345d0

    allocate (coefs_deriv(deg+1))
    do k = 0, deg
        call polyder (coefs, coefs_deriv(:deg+1-k), k, status)
        call polyval (coefs_deriv(:deg+1-k), x, y(k+1), status)
    end do
    deallocate (coefs_deriv)

    call polyfit_deriv (x, y, coefs2, work, status)
    if (.not. (all_close (coefs, coefs2, atol=xtol))) &
        error stop "Polynomial of degree 3, some zero coefs"
    deallocate (y, coefs, coefs2)

end subroutine

end program
