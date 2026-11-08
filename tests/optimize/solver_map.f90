

program test_optimize_solver_map
    !*  Unit test for SOLVER_MAP, a helper type that simplifies mapping
    !   unconstrained domains in R into a (possibly one-sidedly) bounded
    !   interval.

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: ieee_exceptions, only: ieee_get_flag, ieee_set_flag, &
        ieee_overflow, ieee_divide_by_zero, ieee_invalid

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
    call test_logistic_inverse_near_bounds ()
    call test_logistic_endpoints_and_domain ()
    call test_logistic_forward_extremes ()
    call test_exp_near_bounds ()
    call test_solver_map_real32 ()

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


subroutine test_logistic_inverse_near_bounds ()
    !*  Regression tests for interior near-bound and near-midpoint logistic
    !   inverses (Issue #32):
    !   - Subnormal and tiny normal inputs near lb and ub via vector and scalar APIs.
    !   - Non-unit bounds with negative scale parameter (s < 0).
    !   - Cancellation-free evaluation near the interval midpoint.
    !   - Interior point with representable large Jacobian (1e200) without underflow.
    !   - Verification that IEEE overflow, divide-by-zero, and invalid stay clear.

    type (solver_map) :: map, map_vec(1), map_neg
    real (PREC) :: y_vec(1), x_vec(1), jac_vec(1), x_sc, jac_sc, y_sub, exp_x
    type (status_t) :: status
    logical :: flag_ovf, flag_div0, flag_inv, orig_gradual, restore_uflow

    call solver_map_init (map, lb=0.0_PREC, ub=1.0_PREC, status=status)
    if (status /= NF_STATUS_OK) error stop "Near-bounds map init [0, 1] failed"
    map_vec(1) = map

    if (ieee_support_subnormal (0.0_PREC)) then
        restore_uflow = ieee_support_underflow_control (0.0_PREC)
        if (restore_uflow) then
            call ieee_get_underflow_mode (orig_gradual)
            call ieee_set_underflow_mode (.true.)
        end if

        ! Construct subnormal 1e-310 at runtime after enabling gradual underflow.
        y_sub = 1.0e-300_PREC
        y_sub = y_sub * 1.0e-10_PREC

        ! Case 1 (Issue #32): Vector no-Jacobian inverse at y = [1e-310] on [0, 1];
        ! expected x = log(1e-310) - log(1 - 1e-310) = -713.80137882815416.
        y_vec(1) = y_sub
        status = NF_STATUS_UNDEFINED
        call clear_fp_flags ()
        call solver_map_eval_inverse (map_vec, y_vec, x_vec, status=status)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)
        if (status /= NF_STATUS_OK) &
            error stop "Vector inverse at y=[1e-310]: expected NF_STATUS_OK"
        if (abs (x_vec(1) - (-713.80137882815416_PREC)) > 1.0e-12_PREC) &
            error stop "Vector inverse at y=[1e-310]: coordinate mismatch"
        if (flag_ovf .or. flag_div0 .or. flag_inv) &
            error stop "Vector inverse at y=[1e-310]: raised IEEE exception flag"

        ! Case 2: Scalar no-Jacobian inverse at y = 1e-310 on [0, 1].
        status = NF_STATUS_UNDEFINED
        call clear_fp_flags ()
        call solver_map_eval_inverse (map, y_sub, x_sc, status=status)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)
        if (status /= NF_STATUS_OK) &
            error stop "Scalar inverse at y=1e-310: expected NF_STATUS_OK"
        if (abs (x_sc - (-713.80137882815416_PREC)) > 1.0e-12_PREC) &
            error stop "Scalar inverse at y=1e-310: coordinate mismatch"
        if (flag_ovf .or. flag_div0 .or. flag_inv) &
            error stop "Scalar inverse at y=1e-310: raised IEEE exception flag"

        if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)
    end if

    ! Case 3: Normal near-lower-bound y = [1e-300] on [0, 1];
    ! expected x = -300 * log(10) = -690.77552789821368.
    y_vec(1) = 1.0e-300_PREC
    call clear_fp_flags ()
    call solver_map_eval_inverse (map_vec, y_vec, x_vec, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (status /= NF_STATUS_OK) &
        error stop "Vector inverse at y=[1e-300]: expected NF_STATUS_OK"
    if (abs (x_vec(1) - (-690.77552789821368_PREC)) > 1.0e-12_PREC) &
        error stop "Vector inverse at y=[1e-300]: coordinate mismatch"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Vector inverse at y=[1e-300]: raised IEEE exception flag"

    ! Case 4: Near-upper-bound on [-1, 0] at normal y = [-1e-300] (normalized
    ! z = (y - lb)/(ub - lb) rounds to 1.0, whereas ub - y = 1e-300 is exact);
    ! expected x = log(1) - log(1e-300) = +690.77552789821368.
    call solver_map_init (map_vec(1), lb=-1.0_PREC, ub=0.0_PREC, status=status)
    if (status /= NF_STATUS_OK) error stop "Near-bounds map init [-1, 0] failed"
    y_vec(1) = -1.0e-300_PREC
    call clear_fp_flags ()
    call solver_map_eval_inverse (map_vec, y_vec, x_vec, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (status /= NF_STATUS_OK) &
        error stop "Vector inverse on [-1, 0] at y=-1e-300: expected NF_STATUS_OK"
    if (abs (x_vec(1) - 690.77552789821368_PREC) > 1.0e-12_PREC) &
        error stop "Vector inverse on [-1, 0] at y=-1e-300: coordinate mismatch"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Vector inverse on [-1, 0] at y=-1e-300: raised IEEE flag"

    ! Case 5: Subnormal near-upper-bound on [-1, 0] at y = [-1e-310];
    ! expected x = log(1) - log(1e-310) = +713.80137882815416.
    if (ieee_support_subnormal (0.0_PREC)) then
        restore_uflow = ieee_support_underflow_control (0.0_PREC)
        if (restore_uflow) then
            call ieee_get_underflow_mode (orig_gradual)
            call ieee_set_underflow_mode (.true.)
        end if

        y_sub = -1.0e-300_PREC
        y_sub = y_sub * 1.0e-10_PREC
        y_vec(1) = y_sub

        call clear_fp_flags ()
        call solver_map_eval_inverse (map_vec, y_vec, x_vec, status=status)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)

        if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)

        if (status /= NF_STATUS_OK) &
            error stop "Vector inverse on [-1, 0] at y=-1e-310: status /= OK"
        if (abs (x_vec(1) - 713.80137882815416_PREC) > 1.0e-12_PREC) &
            error stop "Vector inverse on [-1, 0] at y=-1e-310: mismatch"
        if (flag_ovf .or. flag_div0 .or. flag_inv) &
            error stop "Vector inverse on [-1, 0] at y=-1e-310: raised IEEE flag"
    end if

    ! Case 6: Non-unit bounds [-2, 6] with negative scale:
    ! y0 = 2 (midpoint), dy = -2 -> y0 + dy = 0 -> s = 1 / (-log(3)) < 0.
    ! At y = 4: dl = 6, du = 2 -> log(dl/du) = log(3) -> x = s * log(3) = -1.0.
    call solver_map_init (map_neg, lb=-2.0_PREC, ub=6.0_PREC, y0=2.0_PREC, &
        dy=-2.0_PREC, dx=1.0_PREC, status=status)
    if (status /= NF_STATUS_OK) error stop "Negative-scale map init failed"
    call clear_fp_flags ()
    call solver_map_eval_inverse (map_neg, 4.0_PREC, x_sc, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    exp_x = -1.0_PREC
    if (status /= NF_STATUS_OK .or. abs (x_sc - exp_x) > 1.0e-14_PREC) &
        error stop "Negative-scale inverse at y=4 failed"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Negative-scale inverse at y=4 raised IEEE flag"

    ! Case 7: Near-midpoint on [0, 100] (s = 1) at y = fl(50 + 1e-12):
    ! exact 100-digit log((y - 0) / (100 - y)) = 4.0074610296869648e-14.
    call solver_map_init (map_vec(1), lb=0.0_PREC, ub=100.0_PREC, status=status)
    y_vec(1) = 50.0_PREC + 1.0e-12_PREC
    call clear_fp_flags ()
    call solver_map_eval_inverse (map_vec, y_vec, x_vec, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (status /= NF_STATUS_OK) &
        error stop "Midpoint inverse on [0, 100]: status /= OK"
    if (abs (x_vec(1) - 4.0074610296869648e-14_PREC) > 1.0e-28_PREC) &
        error stop "Midpoint inverse on [0, 100]: cancellation loss"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Midpoint inverse on [0, 100]: raised IEEE flag"

    ! Case 8: Interior y = 1e-200 on [0, 1] (s = 1) WITH requested Jacobian
    ! (previous implementation underflowed z**2 = 1e-400 to 0 and divided by 0):
    ! expected x = -200 * log(10) = -460.51701859880916, dx/dy = 1e200 + 1 ~ 1e200.
    map_vec(1) = map
    y_vec(1) = 1.0e-200_PREC
    call clear_fp_flags ()
    call solver_map_eval_inverse (map_vec, y_vec, x_vec, jac_vec, status=status)
    call solver_map_eval_inverse (map, 1.0e-200_PREC, x_sc, jac_sc, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (status /= NF_STATUS_OK) &
        error stop "Inverse with Jacobian at y=1e-200: status /= OK"
    if (abs (x_vec(1) - (-460.51701859880916_PREC)) > 1.0e-12_PREC &
            .or. abs (x_sc - (-460.51701859880916_PREC)) > 1.0e-12_PREC) &
        error stop "Inverse with Jacobian at y=1e-200: coordinate mismatch"
    if (abs (jac_vec(1) - 1.0e200_PREC) / 1.0e200_PREC > 1.0e-15_PREC &
            .or. abs (jac_sc - 1.0e200_PREC) / 1.0e200_PREC > 1.0e-15_PREC) &
        error stop "Inverse with Jacobian at y=1e-200: derivative mismatch"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Inverse with Jacobian at y=1e-200: raised IEEE flag"

end subroutine


subroutine test_logistic_endpoints_and_domain ()
    !*  Endpoint, domain-error, unrepresentable-Jacobian, and shape tests:
    !   - Exact endpoints y = lb and y = ub without jac return +/-Inf and status OK.
    !   - Out-of-domain y < lb, y > ub, and quiet NaN return domain error + NaN.
    !   - Unrepresentable requested jac returns +/-Inf with NF_STATUS_INVALID_STATE.
    !   - Non-square 2D jac matrix returns NF_STATUS_INVALID_ARG.

    type (solver_map) :: map, map_neg, maps2(2), map1(1)
    real (PREC) :: y2(2), x2(2), y1(1), x1(1), jac1(1), jac1_mat(1, 1)
    real (PREC) :: bad_jac(2, 1), x_sc, jac_sc, qnan, y_sub
    type (status_t) :: status
    logical :: flag_ovf, flag_div0, flag_inv, orig_gradual, restore_uflow

    qnan = ieee_value (0.0_PREC, IEEE_QUIET_NAN)

    call solver_map_init (map, lb=0.0_PREC, ub=1.0_PREC, status=status)
    maps2(1) = map
    maps2(2) = map
    map1(1) = map

    ! Case 1: Exact endpoints y = [0, 1] on [0, 1] (s = +1) without Jacobian ->
    ! x = [-Inf, +Inf] with status = NF_STATUS_OK and no IEEE exception flags.
    y2 = [0.0_PREC, 1.0_PREC]
    call clear_fp_flags ()
    call solver_map_eval_inverse (maps2, y2, x2, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (status /= NF_STATUS_OK) &
        error stop "Exact endpoints [0, 1] without jac: status /= OK"
    if (.not. (ieee_class (x2(1)) == IEEE_NEGATIVE_INF &
            .and. ieee_class (x2(2)) == IEEE_POSITIVE_INF)) &
        error stop "Exact endpoints [0, 1] without jac: expected [-Inf, +Inf]"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Exact endpoints [0, 1] without jac: raised IEEE flag"

    ! Case 2: Exact endpoints on [-2, 6] with s < 0 without Jacobian ->
    ! y = -2 maps to +Inf and y = 6 maps to -Inf with NF_STATUS_OK.
    call solver_map_init (map_neg, lb=-2.0_PREC, ub=6.0_PREC, y0=2.0_PREC, &
        dy=-2.0_PREC, dx=1.0_PREC, status=status)
    call solver_map_eval_inverse (map_neg, -2.0_PREC, x2(1), status=status)
    if (status /= NF_STATUS_OK &
            .or. .not. (ieee_class (x2(1)) == IEEE_POSITIVE_INF)) &
        error stop "Exact lb with s < 0: expected +Inf and status OK"
    call solver_map_eval_inverse (map_neg, 6.0_PREC, x2(2), status=status)
    if (status /= NF_STATUS_OK &
            .or. .not. (ieee_class (x2(2)) == IEEE_NEGATIVE_INF)) &
        error stop "Exact ub with s < 0: expected -Inf and status OK"

    ! Case 3: Out-of-domain y < lb (-0.01 < 0) -> NF_STATUS_INVALID_ARG + DOMAIN_ERROR.
    call solver_map_eval_inverse (map, -0.01_PREC, x_sc, jac_sc, status=status)
    if (.not. ((NF_STATUS_INVALID_ARG .in. status) &
            .and. (NF_STATUS_DOMAIN_ERROR .in. status) &
            .and. ieee_is_nan (x_sc) .and. ieee_is_nan (jac_sc))) &
        error stop "Out-of-domain y < lb check failed"

    ! Case 4: Out-of-domain y > ub (1.01 > 1) -> NF_STATUS_INVALID_ARG + DOMAIN_ERROR.
    call solver_map_eval_inverse (map, 1.01_PREC, x_sc, status=status)
    if (.not. ((NF_STATUS_INVALID_ARG .in. status) &
            .and. (NF_STATUS_DOMAIN_ERROR .in. status) &
            .and. ieee_is_nan (x_sc))) &
        error stop "Out-of-domain y > ub check failed"

    ! Case 5: Quiet NaN input -> NF_STATUS_INVALID_ARG + NF_STATUS_DOMAIN_ERROR.
    call clear_fp_flags ()
    call solver_map_eval_inverse (map, qnan, x_sc, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (.not. ((NF_STATUS_INVALID_ARG .in. status) &
            .and. (NF_STATUS_DOMAIN_ERROR .in. status) &
            .and. ieee_is_nan (x_sc))) &
        error stop "Quiet NaN inverse domain check failed"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Quiet NaN inverse raised IEEE exception flag"

    ! Case 6: Unrepresentable requested Jacobian at interior y = 1e-310 on [0, 1]
    ! (s = +1: x = -713.80... is finite, dx/dy = 1e310 > huge -> +Inf with
    ! NF_STATUS_INVALID_STATE across scalar, 1D-jac, and 2D-jac APIs).
    if (ieee_support_subnormal (0.0_PREC)) then
        restore_uflow = ieee_support_underflow_control (0.0_PREC)
        if (restore_uflow) then
            call ieee_get_underflow_mode (orig_gradual)
            call ieee_set_underflow_mode (.true.)
        end if

        y_sub = 1.0e-300_PREC
        y_sub = y_sub * 1.0e-10_PREC
        y1(1) = y_sub

        call clear_fp_flags ()
        call solver_map_eval_inverse (map, y_sub, x_sc, jac_sc, status=status)
        if (status /= NF_STATUS_INVALID_STATE &
                .or. abs (x_sc - (-713.80137882815416_PREC)) > 1.0e-12_PREC &
                .or. .not. (ieee_class (jac_sc) == IEEE_POSITIVE_INF)) then
            if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)
            error stop "Scalar unrepresentable Jacobian at y=1e-310 failed"
        end if

        call solver_map_eval_inverse (map1, y1, x1, jac1, status=status)
        if (status /= NF_STATUS_INVALID_STATE &
                .or. .not. (ieee_class (jac1(1)) == IEEE_POSITIVE_INF)) then
            if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)
            error stop "1D unrepresentable Jacobian at y=[1e-310] failed"
        end if

        call solver_map_eval_inverse (map1, y1, x1, jac1_mat, status=status)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)

        if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)

        if (status /= NF_STATUS_INVALID_STATE &
                .or. .not. (ieee_class (jac1_mat(1, 1)) == IEEE_POSITIVE_INF)) &
            error stop "2D unrepresentable Jacobian at y=[1e-310] failed"
        if (flag_ovf .or. flag_div0 .or. flag_inv) &
            error stop "Unrepresentable Jacobian at y=1e-310 raised IEEE flag"
    end if

    ! Case 7: Requested Jacobian at exact lower endpoint y = lb with s < 0 ->
    ! x = +Inf, jac = -Inf (never NaN), status = NF_STATUS_INVALID_STATE.
    call clear_fp_flags ()
    call solver_map_eval_inverse (map_neg, -2.0_PREC, x_sc, jac_sc, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (status /= NF_STATUS_INVALID_STATE &
            .or. .not. (ieee_class (x_sc) == IEEE_POSITIVE_INF) &
            .or. .not. (ieee_class (jac_sc) == IEEE_NEGATIVE_INF)) &
        error stop "Exact endpoint with s < 0 and requested jac failed"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Exact endpoint with requested jac raised IEEE flag"

    ! Case 8: Non-square 2D jac (shape 2x1 for n=2) rejected in eval & eval_inverse.
    y2 = [0.25_PREC, 0.75_PREC]
    x2 = [0.0_PREC, 1.0_PREC]
    call solver_map_eval (maps2, x2, y2, bad_jac, status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Non-square 2D jac in solver_map_eval not rejected"
    call solver_map_eval_inverse (maps2, y2, x2, bad_jac, status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Non-square 2D jac in solver_map_eval_inverse not rejected"

end subroutine


subroutine test_logistic_forward_extremes ()
    !*  Forward logistic tests at large signed coordinates, infinities, and NaN:
    !   - Large |x| = 800 saturates to exact bounds without IEEE_OVERFLOW.
    !   - Subnormal tails x = -710 on [0, 1] and +710 on [-1, 0] preserve y and jac.
    !   - Infinite coordinates +/-Inf with negative scale s < 0.
    !   - Quiet NaN propagation without IEEE_INVALID.

    type (solver_map) :: map01, map_n10, map_neg, maps2(2), map1(1)
    real (PREC) :: x2(2), y2(2), jac2(2), x1(1), y1(1), jac1(1)
    real (PREC) :: y_sc, jac_sc, qnan
    type (status_t) :: status
    logical :: flag_ovf, flag_div0, flag_inv, orig_gradual, restore_uflow

    qnan = ieee_value (0.0_PREC, IEEE_QUIET_NAN)

    call solver_map_init (map01, lb=0.0_PREC, ub=1.0_PREC, status=status)
    call solver_map_init (map_n10, lb=-1.0_PREC, ub=0.0_PREC, status=status)
    call solver_map_init (map_neg, lb=-2.0_PREC, ub=6.0_PREC, y0=2.0_PREC, &
        dy=-2.0_PREC, dx=1.0_PREC, status=status)

    ! Case 1: Vector no-Jacobian forward eval at x = [-800, +800] on [0, 1] ->
    ! exact endpoints [0.0, 1.0] without IEEE_OVERFLOW.
    maps2(1) = map01
    maps2(2) = map01
    x2 = [-800.0_PREC, 800.0_PREC]
    call clear_fp_flags ()
    call solver_map_eval (maps2, x2, y2, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (status /= NF_STATUS_OK .or. y2(1) /= 0.0_PREC .or. y2(2) /= 1.0_PREC) &
        error stop "Forward eval at [-800, 800] without jac failed"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Forward eval at [-800, 800] without jac raised IEEE flag"

    ! Case 2: Vector 1D-Jacobian forward eval at x = [-800, +800] on [0, 1] ->
    ! y = [0.0, 1.0], jac = [0.0, 0.0] without IEEE_OVERFLOW or IEEE_INVALID.
    call clear_fp_flags ()
    call solver_map_eval (maps2, x2, y2, jac2, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (status /= NF_STATUS_OK .or. jac2(1) /= 0.0_PREC .or. jac2(2) /= 0.0_PREC) &
        error stop "Forward eval at [-800, 800] with jac failed"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Forward eval at [-800, 800] with jac raised IEEE flag"

    if (ieee_support_subnormal (0.0_PREC)) then
        restore_uflow = ieee_support_underflow_control (0.0_PREC)
        if (restore_uflow) then
            call ieee_get_underflow_mode (orig_gradual)
            call ieee_set_underflow_mode (.true.)
        end if

        ! Case 3: Subnormal lower tail x = -710 on [0, 1] (s = 1):
        ! y = exp(-710) = 4.4762862256751298e-309 and dy/dx = exp(-710) > 0
        ! (previous z**2 * exp(710) overflowed exp and underflowed z**2).
        call clear_fp_flags ()
        call solver_map_eval (map01, -710.0_PREC, y_sc, jac_sc)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)
        if (abs (y_sc - 4.4762862256751298e-309_PREC) > 1.0e-321_PREC &
                .or. abs (jac_sc - 4.4762862256751298e-309_PREC) > 1.0e-321_PREC) then
            if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)
            error stop "Forward subnormal lower tail at x=-710 failed"
        end if
        if (flag_ovf .or. flag_div0 .or. flag_inv) then
            if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)
            error stop "Forward subnormal lower tail at x=-710 raised IEEE flag"
        end if

        ! Case 4: Subnormal upper tail x = +710 on [-1, 0] (s = 1):
        ! y = 0 - exp(-710) = -4.4762862256751298e-309 (strictly < 0 = ub) and
        ! dy/dx = 4.4762862256751298e-309.
        map1(1) = map_n10
        x1(1) = 710.0_PREC
        call clear_fp_flags ()
        call solver_map_eval (map1, x1, y1, jac1, status=status)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)

        if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)

        if (status /= NF_STATUS_OK &
                .or. abs (y1(1) - (-4.4762862256751298e-309_PREC)) > 1.0e-321_PREC &
                .or. abs (jac1(1) - 4.4762862256751298e-309_PREC) > 1.0e-321_PREC) &
            error stop "Forward subnormal upper tail at x=+710 on [-1, 0] failed"
        if (flag_ovf .or. flag_div0 .or. flag_inv) &
            error stop "Forward subnormal upper tail at x=+710 raised IEEE flag"
    end if

    ! Case 5: Infinite solver coordinates [-Inf, +Inf] with negative scale s < 0
    ! on [-2, 6] -> y = [6.0, -2.0], jac = [0.0, 0.0], no IEEE exception flags.
    maps2(1) = map_neg
    maps2(2) = map_neg
    x2(1) = ninf
    x2(2) = inf
    call clear_fp_flags ()
    call solver_map_eval (maps2, x2, y2, jac2, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (status /= NF_STATUS_OK .or. y2(1) /= 6.0_PREC .or. y2(2) /= -2.0_PREC &
            .or. jac2(1) /= 0.0_PREC .or. jac2(2) /= 0.0_PREC) &
        error stop "Forward eval at [-Inf, +Inf] with s < 0 failed"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Forward eval at [-Inf, +Inf] with s < 0 raised IEEE flag"

    ! Case 6: Quiet NaN propagates to y = NaN and jac = NaN without IEEE_INVALID.
    call clear_fp_flags ()
    call solver_map_eval (map01, qnan, y_sc, jac_sc)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (.not. (ieee_is_nan (y_sc) .and. ieee_is_nan (jac_sc))) &
        error stop "Forward eval at NaN did not propagate NaN"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "Forward eval at NaN raised IEEE flag"

end subroutine


subroutine test_exp_near_bounds ()
    !*  Near-bound, endpoint, and domain tests for TRANSFORM_EXP:
    !   - Subnormal y = lb + 1e-310 without jac avoids reciprocal overflow.
    !   - Same input with requested jac returns +Inf and NF_STATUS_INVALID_STATE.
    !   - Exact lower bound y = lb maps to -Inf without calling log(0).
    !   - Out-of-domain y < lb and NaN return NF_STATUS_INVALID_ARG + DOMAIN_ERROR.

    type (solver_map) :: map, map1(1)
    real (PREC) :: y1(1), x1(1), jac1(1), x_sc, y_sub, qnan
    type (status_t) :: status
    logical :: flag_ovf, flag_div0, flag_inv, orig_gradual, restore_uflow

    qnan = ieee_value (0.0_PREC, IEEE_QUIET_NAN)

    call solver_map_init (map, lb=0.0_PREC, status=status)
    map1(1) = map

    if (ieee_support_subnormal (0.0_PREC)) then
        restore_uflow = ieee_support_underflow_control (0.0_PREC)
        if (restore_uflow) then
            call ieee_get_underflow_mode (orig_gradual)
            call ieee_set_underflow_mode (.true.)
        end if

        y_sub = 1.0e-300_PREC
        y_sub = y_sub * 1.0e-10_PREC
        y1(1) = y_sub

        ! Case 1: Vector inverse without jac at y = [1e-310] (lb = 0, s = 1) ->
        ! x = log(1e-310) = -713.80137882815416, NF_STATUS_OK, no IEEE_OVERFLOW.
        call clear_fp_flags ()
        call solver_map_eval_inverse (map1, y1, x1, status=status)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)
        if (status /= NF_STATUS_OK &
                .or. abs (x1(1) - (-713.80137882815416_PREC)) > 1.0e-12_PREC) then
            if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)
            error stop "EXP inverse at y=[1e-310] without jac failed"
        end if
        if (flag_ovf .or. flag_div0 .or. flag_inv) then
            if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)
            error stop "EXP inverse at y=[1e-310] without jac raised IEEE flag"
        end if

        ! Case 2: Vector inverse with jac at y = [1e-310] -> jac = +Inf and
        ! status = NF_STATUS_INVALID_STATE without IEEE flags.
        call clear_fp_flags ()
        call solver_map_eval_inverse (map1, y1, x1, jac1, status=status)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)

        if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)

        if (status /= NF_STATUS_INVALID_STATE &
                .or. .not. (ieee_class (jac1(1)) == IEEE_POSITIVE_INF)) &
            error stop "EXP inverse at y=[1e-310] with jac failed"
        if (flag_ovf .or. flag_div0 .or. flag_inv) &
            error stop "EXP inverse at y=[1e-310] with jac raised IEEE flag"
    end if

    ! Case 3: Exact lower endpoint y = lb = 2.5 (s = 1) without jac ->
    ! x = -Inf, status = NF_STATUS_OK, no IEEE divide-by-zero flag.
    call solver_map_init (map, lb=2.5_PREC, status=status)
    call clear_fp_flags ()
    call solver_map_eval_inverse (map, 2.5_PREC, x_sc, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (status /= NF_STATUS_OK &
            .or. .not. (ieee_class (x_sc) == IEEE_NEGATIVE_INF)) &
        error stop "EXP inverse at exact y = lb failed"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "EXP inverse at exact y = lb raised IEEE flag"

    ! Case 4: Out-of-domain y < lb and NaN -> NF_STATUS_INVALID_ARG + DOMAIN_ERROR.
    call solver_map_eval_inverse (map, 2.4_PREC, x_sc, status=status)
    if (.not. ((NF_STATUS_INVALID_ARG .in. status) &
            .and. (NF_STATUS_DOMAIN_ERROR .in. status))) &
        error stop "EXP inverse out-of-domain y < lb failed"
    call solver_map_eval_inverse (map, qnan, x_sc, status=status)
    if (.not. ((NF_STATUS_INVALID_ARG .in. status) &
            .and. (NF_STATUS_DOMAIN_ERROR .in. status))) &
        error stop "EXP inverse quiet NaN domain check failed"

end subroutine


subroutine test_solver_map_real32 ()
    !*  Single-precision (solver_map_real32) boundary and overflow checks:
    !   - Subnormal y = 1e-40_real32 inverse without jac returns finite coordinate.
    !   - Same input with requested jac returns +Inf and NF_STATUS_INVALID_STATE.
    !   - Forward evaluation at x = [-120, +120]_real32 avoids IEEE_OVERFLOW.

    type (solver_map_real32) :: map1(1), maps2(2)
    real (real32) :: y1(1), x1(1), jac1(1), x2(2), y2(2), jac2(2), y_sub
    type (status_t) :: status
    logical :: flag_ovf, flag_div0, flag_inv, orig_gradual, restore_uflow

    call solver_map_init (map1(1), lb=0.0_real32, ub=1.0_real32, status=status)
    if (status /= NF_STATUS_OK) error stop "real32 map init [0, 1] failed"

    if (ieee_support_subnormal (0.0_real32)) then
        restore_uflow = ieee_support_underflow_control (0.0_real32)
        if (restore_uflow) then
            call ieee_get_underflow_mode (orig_gradual)
            call ieee_set_underflow_mode (.true.)
        end if

        y_sub = 1.0e-30_real32
        y_sub = y_sub * 1.0e-10_real32
        y1(1) = y_sub

        ! Case 1: Vector no-Jacobian inverse at y = [1e-40_real32] (1/y overflows
        ! real32 huge ~ 3.4e38); expected log(1e-40) = -92.1034088.
        call clear_fp_flags ()
        call solver_map_eval_inverse (map1, y1, x1, status=status)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)
        if (status /= NF_STATUS_OK &
                .or. abs (x1(1) - (-92.1034088_real32)) > 1.0e-4_real32) then
            if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)
            error stop "real32 inverse at y=[1e-40] without jac failed"
        end if
        if (flag_ovf .or. flag_div0 .or. flag_inv) then
            if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)
            error stop "real32 inverse at y=[1e-40] without jac raised IEEE flag"
        end if

        ! Case 2: Same input with requested Jacobian -> +Inf and INVALID_STATE.
        call clear_fp_flags ()
        call solver_map_eval_inverse (map1, y1, x1, jac1, status=status)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)

        if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)

        if (status /= NF_STATUS_INVALID_STATE &
                .or. .not. (ieee_class (jac1(1)) == IEEE_POSITIVE_INF)) &
            error stop "real32 inverse at y=[1e-40] with jac failed"
        if (flag_ovf .or. flag_div0 .or. flag_inv) &
            error stop "real32 inverse at y=[1e-40] with jac raised IEEE flag"
    end if

    ! Case 3: Forward real32 evaluation at [-120, +120] on [-1, 2] ->
    ! exact [-1.0, 2.0] and zero derivatives without IEEE_OVERFLOW.
    call solver_map_init (maps2(1), lb=-1.0_real32, ub=2.0_real32, status=status)
    maps2(2) = maps2(1)
    x2 = [-120.0_real32, 120.0_real32]
    call clear_fp_flags ()
    call solver_map_eval (maps2, x2, y2, jac2, status=status)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (status /= NF_STATUS_OK .or. y2(1) /= -1.0_real32 .or. y2(2) /= 2.0_real32 &
            .or. jac2(1) /= 0.0_real32 .or. jac2(2) /= 0.0_real32) &
        error stop "real32 forward extremes failed"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "real32 forward extremes raised IEEE flag"

end subroutine


subroutine clear_fp_flags ()
    !*  Clear IEEE overflow, divide-by-zero, and invalid exception flags.

    call ieee_set_flag (ieee_overflow, .false.)
    call ieee_set_flag (ieee_divide_by_zero, .false.)
    call ieee_set_flag (ieee_invalid, .false.)

end subroutine


subroutine get_fp_flags (flag_ovf, flag_div0, flag_inv)
    !*  Query IEEE overflow, divide-by-zero, and invalid exception flags.
    logical, intent(out) :: flag_ovf
        !*  True if IEEE_OVERFLOW is signaling.
    logical, intent(out) :: flag_div0
        !*  True if IEEE_DIVIDE_BY_ZERO is signaling.
    logical, intent(out) :: flag_inv
        !*  True if IEEE_INVALID is signaling.

    call ieee_get_flag (ieee_overflow, flag_ovf)
    call ieee_get_flag (ieee_divide_by_zero, flag_div0)
    call ieee_get_flag (ieee_invalid, flag_inv)

end subroutine


end program
