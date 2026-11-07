

program test_optimize_root_broyden

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_testing
    use numfort_common_workspace, workspace => workspace_real64

    use numfort_optimize, only: root_broyden, optim_result => optim_result_real64

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_quad ()
    call test_func2 ()

    call test_maxiter ()

end subroutine

subroutine test_quad ()

    type (optim_result) :: res
    type (workspace) :: work
    real (PREC), dimension(:), allocatable :: x, fx, x0
    integer :: n

    n = 2
    allocate (x(n), fx(n), x0(n))
    x0 = [-2.0d0, 10.0d0]
    x = x0
    call root_broyden (quad_fcn_jac, x, tol=1.0d-8, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Quadratic function; x0=[-2.000, 10.000]"
    deallocate (x, fx, x0)

    ! User-provided WORKSPACE object
    n = 2
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0]
    x = x0
    call root_broyden (quad_fcn_jac, x, tol=1.0d-8, work=work, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Quadratic function; WORK provided; x0=[4.000, 1.000]"
    deallocate (x, fx, x0)

    ! Limit absolute step size
    n = 2
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0]
    x = x0
    call root_broyden (quad_fcn_jac, x, tol=1.0d-8, xstep=0.1d0, work=work, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Quadratic function; WORK, XSTEP provided; x0=[4.000, 1.000]"
    deallocate (x, fx, x0)

    ! Limit relative step size
    n = 2
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, -10.0d0]
    x = x0
    call root_broyden (quad_fcn_jac, x, tol=1.0d-8, rstep=0.1d0, work=work, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Quadratic function; WORK, RSTEP provided; x0=[4.000, -10.000]"
    deallocate (x, fx, x0)

end subroutine


subroutine test_func2 ()

    type (optim_result) :: res
    real (PREC), dimension(:), allocatable :: x, fx, x0
    integer :: n

    ! Separate FCN and JAC routines
    n = 3
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0, 1.0d0]
    x = x0
    call root_broyden (func2_fcn, func2_jac, x, tol=1.0d-8, xstep=0.1d0, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Separate FCN, JAC; XSTEP provided; x0=[4.000, 1.000, 1.000]"
    deallocate (x, fx, x0)

    ! Numerical differentiation
    n = 3
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0, 1.0d0]
    x = x0
    call root_broyden (func2_fcn, x, ndiff=.true., tol=1.0d-8, dstep=1.0d-5, res=res)
    fx = res%fx
    if (.not. (maxval (abs (fx)) < 1.0d-7 .and. res%success)) &
        error stop "Numerical diff.; x0=[4.000, 1.000, 1.000]"
    deallocate (x, fx, x0)

end subroutine


subroutine test_maxiter ()
    !*  Test that

    type (optim_result) :: res
    real (PREC), dimension(:), allocatable :: x, fx, x0
    integer :: n
    logical :: lres

    ! Separate FCN and JAC routines
    n = 3
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0, 1.0d0]
    x = x0
    call root_broyden (func2_fcn, func2_jac, x, tol=1.0d-8, xstep=0.1d0, &
        maxiter=1, res=res)
    lres = (.not. res%success) .and. res%status == NF_STATUS_MAX_ITER .and. &
        res%nit == 1
    if (.not. (lres)) error stop "Separate FCN, JAC; MAXITER=1"

    x = x0
    call root_broyden (func2_fcn, func2_jac, x, tol=1.0d-8, xstep=0.1d0, &
        maxfev=1, res=res)
    lres = (.not. res%success) .and. res%status == NF_STATUS_MAX_EVAL .and. &
        res%nit == 1
    if (.not. (lres)) error stop "Separate FCN, JAC; MAXFUN=1"

    deallocate (x, fx, x0)

    ! Numerical differentiation
    n = 3
    allocate (x(n), fx(n), x0(n))
    x0 = [4.0d0, 1.0d0, 1.0d0]
    x = x0
    call root_broyden (func2_fcn, x, ndiff=.true., tol=1.0d-8, dstep=1.0d-5, &
        maxiter=1, res=res)
    lres = (.not. res%success) .and. res%status == NF_STATUS_MAX_ITER .and. &
        res%nit == 1
    if (.not. (lres)) error stop "Numerical diff.; MAXITER=1"

    x = x0
    call root_broyden (func2_fcn, x, ndiff=.true., tol=1.0d-8, dstep=1.0d-5, &
        maxfev=1, res=res)
    lres = (.not. res%success) .and. res%status == NF_STATUS_MAX_EVAL .and. &
        res%nit == 1
    if (.not. (lres)) error stop "Numerical diff.; MAXFUN=1"

    deallocate (x, fx, x0)

end subroutine

subroutine quad_fcn_jac (x, fx, fpx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: fpx

    if (present(fx)) then
        fx(1) = (x(1)-2.0)**2.0
        fx(2) = (x(2)+1.0)**2.0
    end if

    if (present(fpx)) then
        fpx = 0.0
        fpx(1,1) = 2.0 * (x(1) - 2.0)
        fpx(2,2) = 2.0 * (x(2) + 1.0)
    end if
end subroutine


subroutine func2_fcn (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    fx(1) = sin(x(1))
    fx(2) = (x(2) - 2.0) ** 2.0
    fx(3) = cos(x(1)) * x(3)
end subroutine


subroutine func2_jac (x, fpx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:,:), contiguous :: fpx

    fpx = 0.0

    fpx(1,1) = cos(x(1))
    fpx(2,2) = 2.0 * (x(2) - 2.0)
    fpx(3,1) = -sin(x(1)) * x(3)
    fpx(3,3) = cos(x(1))
end subroutine

end program
