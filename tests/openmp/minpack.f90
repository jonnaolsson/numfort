

program test_minpack
    !*  Unit tests for MINPACK wrapper routines called from OpenMP
    !   parallel regions.

    use, intrinsic :: iso_fortran_env

    use numfort_common, workspace => workspace_real64
    use numfort_optimize, optim_result => optim_result_real64

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_hybrj ()

end subroutine


subroutine test_hybrj ()

    real (PREC), dimension(:,:), allocatable :: xx, fxx
    real (PREC), dimension(:), allocatable :: x, fx
    type (optim_result) :: res
    integer :: m, n, i
    logical :: values_ok
    real (PREC), parameter :: x0(3) = 50.0

    ! Compute in serial mode

    m = 3
    allocate (x(m), fx(m))
    x(:) = x0

    call root_hybrj (fcn1_jac, x, fx, res=res)

    ! Compute in the same result in parallel
    n = 1000
    allocate (xx(m,n), fxx(m,n))

    do i = 1, n
        xx(:,i) = x0
    end do

    call test_hybrj_dispatch (xx, fxx)

    values_ok = .true.
    do i  = 1, n
        values_ok = values_ok .and. all(abs(xx(:,i)-x) < 1.0e-12_PREC) &
            .and. all(abs(fxx(:,i) - fx) < 1.0e-12_PREC)
    end do

    if (.not. (values_ok)) error stop "Run same problem in parallel"

end subroutine


subroutine test_hybrj_dispatch (xx, fxx)
    real (PREC), intent(out), dimension(:,:), contiguous :: xx, fxx

    integer :: i, n
    type (optim_result) :: res
    type (workspace) :: work

    n = size(xx, 2)

    !$omp parallel default(private) shared(xx, fxx, n)

    !$omp do
    do i = 1, n

        call root_hybrj (fcn1_jac, xx(:,i), fxx(:,i), work=work, res=res)
    end do
    !$omp end do

    !$omp end parallel

end subroutine


pure subroutine fcn1_jac (x, fx, jac, task)
    real (PREC), dimension(:), intent(in), contiguous :: x
    real (PREC), dimension(:), intent(out), contiguous :: fx
    real (PREC), dimension(:,:), intent(out), contiguous :: jac
    integer, intent(inout) :: task

    if (task == 1) then
        fx(1) = x(1)**2 + 2*x(2) - 2*x(3) - 2
        fx(2) = 3*x(1) - x(2)**2 + x(3) - 1
        fx(3) = 5*x(1) - 2*x(2) - x(3)**2 + 12
    else if (task == 2) then
        jac(1, :) = [2*x(1), 2.0d0, -2.0d0]
        jac(2, :) = [3.0d0, -2*x(2), 1.0d0]
        jac(3, :) = [5.0d0, -2.0d0, -2*x(3)]
    end if
end subroutine

end
