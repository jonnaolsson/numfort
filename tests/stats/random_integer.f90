

program test_numfort_stats_random_integer
    !*  Unit tests for PRNG routines that operate on integers

    use iso_fortran_env

    use numfort_stats, only: random_order, random_integer, set_seed

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

contains

subroutine test_all ()

    call test_random_order ()
    call test_random_integer ()
    ! call test_sub2ind ()

end subroutine

subroutine test_random_order ()

    integer, dimension(:), allocatable :: x
    integer :: n, low
    logical correct

    n = 100
    ! call with default low parameter
    allocate (x(n))
    call random_order (x)
    call random_order_verify (x, 1, correct)

    if (.not. (correct)) error stop "argument x(100), default low parameter"
    deallocate (x)

    ! Check with non-default positive argument
    allocate (x(n))
    low = 11
    call random_order (x, low=low)
    call random_order_verify (x, low, correct)
    if (.not. (correct)) error stop "argument x(100), low=11"
    deallocate (x)

    ! Check with non-default negative argument
    allocate (x(n))
    low = -1713
    call random_order (x, low=low)
    call random_order_verify (x, low, correct)
    if (.not. (correct)) error stop "argument x(100), low=-1713"
    deallocate (x)

    ! Check with degenerate size(x) = 1
    n = 1
    allocate (x(n))
    low = 1
    call random_order (x, low=low)
    call random_order_verify (x, low, correct)
    if (.not. (correct)) error stop "argument x(1), low=1"
    deallocate (x)

    ! Check with degenerate size(x) = 0
    n = 0
    allocate (x(n))
    low = 1
    call random_order (x, low=low)
    ! nothing to assert, just check that it does not crash

end subroutine

subroutine random_order_verify (x, low, correct)
    integer, intent(in), dimension(:) :: x
    integer, intent(in) :: low
    logical, intent(out) :: correct

    integer, dimension(:), allocatable :: counter
    integer :: i, j, imin, imax

    allocate (counter(size(x)), source = 0)

    ! check that all integers {low,...,low+n-1} are present exactly once
    do i = 1, size(x)
        j = x(i) - low + 1
        counter(j) = counter(j) + 1
    end do

    imin = minval(x)
    imax = maxval(x)

    correct = all(counter == 1) .and. (imin == low) .and. (imax == low + size(x) - 1)
end subroutine

subroutine test_random_integer ()

    integer, dimension(:), allocatable :: x
    real (PREC), dimension(:), allocatable :: pmf
    integer :: n, lb, ub, nint, i
    logical :: all_ok

    call set_seed (1234)

    n = int(1e6)
    allocate (x(n))

    ! test with degenerate distribution
    call random_integer (x, 0, 0)
    if (.not. (all (x == 0))) error stop "Degenerate draws all equal to 0"

    call random_integer (x, -1, -1)
    if (.not. (all (x == - 1))) error stop "Degenerate draws all equal to -1"

    lb = -1
    ub = 5
    call random_integer (x, lb, ub)

    nint = ub - lb + 1
    allocate (pmf(lb:ub))
    do i = lb, ub
        pmf(i) = count(x == i) / real(n, PREC)
    end do

    all_ok = all(abs(pmf-1.0_PREC/nint) < 1.0d-3) .and. all(x >= lb) &
        .and. all(x <= ub)
    if (.not. (all_ok)) error stop "Uniform distribution on [-1, 5]"

end subroutine

end program
