

program test_optimize_solver_map
    !*  Unit test for SOLVER_MAP, a helper type that simplifies mapping
    !   unconstrained domains in R into a (possibly one-sidedly) bounded
    !   interval.

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_arrays
    use numfort_common
    use numfort_common_testing
    use numfort_optimize, solver_map => solver_map_real64
    use numfort_stats, only: set_seed

    implicit none

    integer, parameter :: PREC = real64

    real (PREC) :: ninf
    real (PREC) :: inf

    ninf = ieee_value (0.0_PREC, IEEE_NEGATIVE_INF)
    inf = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)

    call test_all ()

    contains

subroutine test_all ()

    call test_linear ()
    call test_exp ()
    call test_logistic ()

end subroutine

subroutine test_linear ()

        ! Note: y0 parameter is irrelevant for linear mappings
    real (PREC), parameter :: dx_all(*) = [1.0d-6, 0.1d0, 1.0d0, 2.34d0]
    real (PREC), parameter :: dy_all(*) = [-1.234d0, -0.123d0, 0.1d0, 1.0d0, 2.234d0]

    type (solver_map), dimension(:), allocatable :: maps
    type (solver_map) :: map
    real (PREC), dimension(:), allocatable :: x, y, x1, y1
    real (PREC), dimension(:), allocatable :: jac_diag, jac_inv_diag
    real (PREC), dimension(:,:), allocatable :: jac, jac_inv, eye, xx
    real (PREC) :: dx, dy
    integer :: i, j, k, nx, ny, n
    type (status_t) :: status

    ! Input checks with invalid arguments
    status = NF_STATUS_UNDEFINED
    call solver_map_init (map, lb=ninf, dy=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "INIT called with invalid dy=0"

    status = NF_STATUS_UNDEFINED
    call solver_map_init (map, lb=ninf, dx=0.0_PREC, dy=1.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "INIT called with invalid dx=0"

    nx = size(dx_all)
    ny = size(dy_all)
    n = nx * ny

    allocate (maps(n))

    do i = 1, size(dx_all)
        do j = 1, size(dy_all)
            k = (i-1) * ny + j

            dx = dx_all(i)
            dy = dy_all(j)

            status = NF_STATUS_UNDEFINED
            call solver_map_init (maps(k), lb=ninf, dy=dy, dx=dx, status=status)
            if (.not. (status == NF_STATUS_OK)) &
                error stop "Linear mapping initialization failed"
        end do
    end do

    ! Test mappings and its inverse
    allocate (x(n), y(n), x1(n), y1(n))

    call set_seed (1234)
    call random_number (x)
    x(:) = x * 10.0

    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps, x, y, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL called at values X"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps, y, x1, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL_INVERSE called at values f(X)"

    if (.not. (status == NF_STATUS_OK .and. all_close (x, x1, atol=1.0e-10_PREC, &
        rtol=0.0_PREC))) &
        error stop "Checking X == f^-1(f(X))"

    status = NF_STATUS_UNDEFINED
    call random_number (y)
    call solver_map_eval_inverse (maps, y, x, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL_INVERSE called at values Y"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps, x, y1, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (y1, y, atol=1.0e-10_PREC, &
        rtol=0.0_PREC))) &
        error stop "Checking Y == f(f^-1(Y))"

    ! Test diagonal Jacobians
    allocate (jac_diag(n), jac_inv_diag(n))
    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps, x, y, jac_diag, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL called with 1-d JAC argument"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps, y, x1, jac_inv_diag, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "EVAL_INVERSE called with 1-d JAC argument"

    if (.not. (all_close (jac_inv_diag, 1.0_PREC / jac_diag))) &
        error stop "Checking inverse of Jacobian == Jacobian of inverse"

    deallocate (jac_diag, jac_inv_diag)

    ! Test Jacobians
    allocate (jac(n,n), jac_inv(n,n), eye(n,n), xx(n,n))
    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps, x, y, jac, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL called with 2-d JAC argument"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps, y, x1, jac_inv, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "EVAL_INVERSE called with 2-d JAC argument"

    xx = matmul (jac, jac_inv)
    call identity (eye)
    if (.not. (all_close (xx, eye, atol=1.0e-10_PREC, rtol=0.0_PREC))) &
        error stop "Checking inverse of Jacobian == Jacobian of inverse"

    deallocate (jac, jac_inv)

end subroutine

subroutine test_exp ()

        ! Note: y0 parameter is irrelevant for linear mappings
    real (PREC), parameter :: dx_all(*) = [1.0d-6, 0.1d0, 1.34d0]
    real (PREC), parameter :: lb_all(*) = [-2.423d0, 0.0d0, 1.234d0]

    type (solver_map), dimension(:), allocatable :: maps
    type (solver_map) :: map
    real (PREC), dimension(:), allocatable :: xx, yy, xx1, yy1, diff_x
    real (PREC), dimension(:), allocatable :: jac_diag, jac_inv_diag
    real (PREC), dimension(:,:), allocatable :: jac, jac_inv, eye, mat
    real (PREC), dimension(:), allocatable :: lb_flat
    logical, dimension(:), allocatable :: mask
    real (PREC) :: dx, dy, lb, y0, x0, y, diff, eps
    integer :: i, j, k, n, nmax
    type (status_t) :: status
    logical :: values_ok

    !  === Input checks ===
    status = NF_STATUS_UNDEFINED
    call solver_map_init (map, lb=0.0_PREC, dy=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "INIT called with invalid dy=0"

    status = NF_STATUS_UNDEFINED
    call solver_map_init (map, lb=0.0_PREC, dx=0.0_PREC, dy=1.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "INIT called with invalid dx=0"

    ! === Verify scaling parameter ===

    ! Check that scaling parameter is computed correctly to yield the
    ! desired dy/dx and some given point.

    nmax = size(dx_all) * size(lb_all) * 20

    allocate (maps(nmax), lb_flat(nmax))

    call set_seed (1234)

    n = 0
    do i = 1, size(lb_all)
        lb = lb_all(i)

        do j = 1, size(dx_all)
            dx = dx_all(j)

            do k = 1, 10

                ! Set some valid points at which to pin down derivative
                call random_number (eps)
                y0 = lb + eps * 2.0

                ! dy needs to satisfy (y0+dy) > lb
                call random_number (eps)
                dy = (lb-y0) + eps * 2.0d0

                status = NF_STATUS_UNDEFINED
                call solver_map_init (map, lb=lb, y0=y0, dy=dy, dx=dx, &
                    status=status)

                n = n + 1
                maps(n) = map
                lb_flat(n) = lb

                ! Check whether dy/dx evaluated at y0 is as expected
                ! First obtain implied x0
                call solver_map_eval_inverse (map, y0, x0)
                ! Evaluate at x0 + dx
                call solver_map_eval (map, x0+dx, y)

                ! Verify that DX moves Y by DY around (X0, Y0)
                diff = (y-y0) - dy
                if (.not. (status == NF_STATUS_OK .and. abs (diff) < 1.0e-8_PREC)) &
                    error stop "Exponential mapping did not produce the requested change"
            end do
        end do
    end do

    ! === Test inverse ===
    allocate (xx(n), yy(n), xx1(n), yy1(n))
    allocate (diff_x(n))
    allocate (mask(n))

    ! Draw some random number in R
    call random_number (xx)
    xx(:) = (xx - 0.5d0) * 10.0d0

    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), xx, yy, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL called at values X"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps(1:n), yy, xx1, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL_INVERSE called at values f(X)"

    ! When checking whether values get mapped back to the same X, ignore
    ! y-values close to lower bound, as there we have less precision.
    mask(:) = (yy - lb_flat(1:n)) > 0.001_PREC * abs(lb_flat(1:n))
    ! Ignore values that get mapped into boundary values +/x Inf
    mask(:) = mask .and. ieee_is_finite (xx)
    mask(:) = mask .and. ieee_is_finite (xx1)
    where (mask)
        diff_x = abs(xx - xx1)
    else where
        diff_x = 0.0
    end where
    values_ok = all(diff_x < 1.0e-10_PREC)

    if (.not. (values_ok)) error stop "Checking X == f^-1(f(X))"

    status = NF_STATUS_UNDEFINED
    call random_number (yy)
    ! Make sure no Y is outside of admissible range
    yy(:) = yy + lb_flat(1:n)
    call solver_map_eval_inverse (maps(1:n), yy, xx, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL_INVERSE called at values Y"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), xx, yy1, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (yy1, yy, atol=1.0e-10_PREC, &
        rtol=0.0_PREC))) &
        error stop "Checking Y == f(f^-1(Y))"

    ! Test diagonal Jacobians
    allocate (jac_diag(n), jac_inv_diag(n))
    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), xx, yy, jac_diag, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL called with 1-d JAC argument"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps(1:n), yy, xx1, jac_inv_diag, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "EVAL_INVERSE called with 1-d JAC argument"

    if (.not. (all_close (jac_inv_diag, 1.0_PREC / jac_diag))) &
        error stop "Checking inverse of Jacobian == Jacobian of inverse"

    deallocate (jac_diag, jac_inv_diag)

    ! Test Jacobians
    allocate (jac(n,n), jac_inv(n,n), eye(n,n), mat(n,n))
    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), xx, yy, jac, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL called with 2-d JAC argument"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps(1:n), yy, xx1, jac_inv, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "EVAL_INVERSE called with 2-d JAC argument"

    mat = matmul (jac, jac_inv)
    call identity (eye)
    if (.not. (all_close (mat, eye, atol=1.0e-10_PREC, rtol=0.0_PREC))) &
        error stop "Checking inverse of Jacobian == Jacobian of inverse"

    deallocate (jac, jac_inv)

end subroutine

subroutine test_logistic ()
    !*  Unit tests for "logistic" mapping f: R -> [lb,ub]

        ! Note: y0 parameter is irrelevant for linear mappings
    real (PREC), parameter :: dx_all(*) = [1.0d-6, 0.1d0, 1.34d0]
    real (PREC), parameter :: lb_all(*) = [-2.423d0, 0.0d0, 1.234d0]
    real (PREC), parameter :: ub_all(*) = [-1.234d0, 1.0d0, 3.450d0]

    type (solver_map), dimension(:), allocatable :: maps
    type (solver_map) :: map
    real (PREC), dimension(:), allocatable :: xx, yy, xx1, yy1, dd, diff_x
    real (PREC), dimension(:), allocatable :: jac_diag, jac_inv_diag
    real (PREC), dimension(:,:), allocatable :: jac, jac_inv, eye, mat
    real (PREC), dimension(:), allocatable :: lb_flat, ub_flat
    logical, dimension(:), allocatable :: mask
    real (PREC) :: dx, dy, lb, ub, y0, x0, y, diff, eps
    integer :: i, j, k, n, nmax, ii
    type (status_t) :: status
    logical :: values_ok

    ! === Input checks ===
    status = NF_STATUS_UNDEFINED
    call solver_map_init (map, lb=0.0_PREC, dy=0.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "INIT called with invalid dy=0"

    status = NF_STATUS_UNDEFINED
    call solver_map_init (map, lb=0.0_PREC, dx=0.0_PREC, dy=1.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "INIT called with invalid dx=0"

    ! === Verify scaling parameter ===

    ! Check that scaling parameter is computed correctly to yield the
    ! desired dy/dx and some given point.
    nmax = size(dx_all) * size(lb_all) * 10

    allocate (maps(nmax))
    allocate (lb_flat(nmax), ub_flat(nmax))

    call set_seed (1234)

    n = 0
    do i = 1, size(lb_all)
        lb = lb_all(i)
        ub = ub_all(i)

        do j = 1, size(dx_all)
            dx = dx_all(j)

            do k = 1, 10

                ii = (i - 1)*size(lb_all) + (j-1)*size(dx_all) + k
                if (modulo(ii, 2) == 0) then
                    ! Position y0 at some random interior point
                    call random_number (eps)
                    y0 = lb + (ub-lb) * eps
                else
                    ! Position y0 exactly at midpoint as this can lead to
                    ! additional complications
                    y0 = (ub + lb) / 2.0_PREC
                end if

                ! dy needs to satisfy: (y0+dy) > lb, (y0+dy) < ub
                ! and thus dy > lb-y0, dy < ub-y0
                call random_number (eps)
                dy = (lb-y0) + (ub-lb) * eps

                status = NF_STATUS_UNDEFINED
                call solver_map_init (map, lb=lb, ub=ub, y0=y0, dy=dy, dx=dx, &
                    status=status)

                n = n + 1
                maps(n) = map
                lb_flat(n) = lb
                ub_flat(n) = ub

                ! Check whether dy/dx evaluated at y0 is as expected
                ! First obtain implied x0
                call solver_map_eval_inverse (map, y0, x0)
                ! Evaluate at x0 + dx
                call solver_map_eval (map, x0+dx, y)

                ! Verify that DX moves Y by DY around (X0, Y0)
                diff = (y-y0) - dy
                if (.not. (status == NF_STATUS_OK .and. abs (diff) < 1.0e-8_PREC)) &
                    error stop "Logistic mapping did not produce the requested change"
            end do
        end do
    end do

    ! === Test inverse ===
    allocate (xx(n), yy(n), xx1(n), yy1(n), dd(n))
    allocate (diff_x(n))
    allocate (mask(n))

    call random_number (xx)
    xx(:) = (xx - 0.50d0) * 10.0

    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), xx, yy, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL called at values X"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps(1:n), yy, xx1, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL_INVERSE called at values f(X)"

    ! When checking whether values get mapped back to the same X, ignore
    ! y-values close to bounds, as there we have less precision.
    mask(:) = (yy - lb_flat(1:n)) > 0.001_PREC * abs(lb_flat(1:n))
    mask(:) = mask .and. (ub_flat(1:n) - yy) > 0.001_PREC * abs(ub_flat(1:n))
    ! Ignore values that get mapped into boundary values +/x Inf
    mask(:) = mask .and. ieee_is_finite (xx)
    mask(:) = mask .and. ieee_is_finite (xx1)
    where (mask)
        diff_x = abs(xx - xx1)
    else where
        diff_x = 0.0
    end where
    values_ok = all(diff_x < 1.0e-10_PREC)

    if (.not. (values_ok)) error stop "Checking X == f^-1(f(X))"

    status = NF_STATUS_UNDEFINED
    call random_number (yy)
    ! Make sure no Y is outside of admissible range
    yy(:) = lb_flat + (ub_flat - lb_flat) * yy
    call solver_map_eval_inverse (maps(1:n), yy, xx, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL_INVERSE called at values Y"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), xx, yy1, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (yy1, yy, atol=1.0e-10_PREC, &
        rtol=0.0_PREC))) &
        error stop "Checking Y == f(f^-1(Y))"

    ! Test diagonal Jacobians
    allocate (jac_diag(n), jac_inv_diag(n))
    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), xx, yy, jac_diag, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL called with 1-d JAC argument"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps(1:n), yy, xx1, jac_inv_diag, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "EVAL_INVERSE called with 1-d JAC argument"

    dd(:) = jac_inv_diag - 1.0_PREC/jac_diag
    if (.not. (all_close (jac_inv_diag, 1.0_PREC / jac_diag, atol=1.0e-8_PREC))) &
        error stop "Checking inverse of Jacobian == Jacobian of inverse"

    deallocate (jac_diag, jac_inv_diag)

    ! Test Jacobians
    allocate (jac(n,n), jac_inv(n,n), eye(n,n), mat(n,n))
    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), xx, yy, jac, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "EVAL called with 2-d JAC argument"

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps(1:n), yy, xx1, jac_inv, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "EVAL_INVERSE called with 2-d JAC argument"

    mat = matmul (jac, jac_inv)
    call identity (eye)
    if (.not. (all_close (mat, eye, atol=1.0e-10_PREC, rtol=0.0_PREC))) &
        error stop "Checking inverse of Jacobian == Jacobian of inverse"

    deallocate (jac, jac_inv)

end subroutine

end
