

program test_stats_lm
    !*  Unit tests for OLS estimation routines in numfort_stats_lm module.

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_negative_inf

    use numfort_common
    use numfort_common_testing
    use numfort_stats, only: set_seed
    use numfort_stats_lm, lm_result => lm_result_real64, lm_config => lm_config_real64
    use numfort_stats_data_helpers
    use numfort_stats_dnorm

    use blas95, only: BLAS95_GEMM => GEMM, BLAS95_GEMV => GEMV

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_OLS_1d_args ()
    call test_OLS_2d_args ()
    call test_OLS_1d ()
    call test_OLS_2d ()
    call test_OLS_fit_equiv ()

    call test_predict_equiv_1d ()
    call test_predict_equiv_2d ()

    call test_WLS ()

    call test_pcr_equiv_1d ()
    call test_pcr_equiv_2d ()

    call test_pcr_ols_equiv_1d ()
    call test_pcr_ols_equiv_2d ()

    call test_pcr_masked ()
    call test_pcr_masked_const ()
    call test_pcr_cv_masked ()
    call test_transform_weights ()
    call test_pcr_weighted_reference ()
    call test_pcr_preprocessing ()
    call test_pcr_invalid_weights ()

end subroutine

subroutine test_OLS_1d_args ()
    !*  Unit tests for OLS argument checking

    real (PREC), dimension(:,:), allocatable :: x
    real (PREC), dimension(:), allocatable :: y
    real (PREC), dimension(:), allocatable :: beta

    integer :: nrhs, nobs, ncoefs
    type (lm_config) :: conf
    type (status_t) :: status

    ! === incompatible X, Y arguments ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs+1), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X, Y"

    deallocate (y, x, beta)

    ! === incompatible X^T, Y ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs, nobs), y(nobs+1), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config(trans_x=.true.)
    call ols (conf, x, y, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X^T, Y"

    deallocate (x, y, beta)

    ! === incompabible X, COEFS ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs-2))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X, COEFS"

    deallocate (x, y, beta)

    ! === incompabible X^T, COEFS ===

    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs+2,nobs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config(trans_x=.true.)
    call ols (conf, x, y, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X^T, COEFS"

    deallocate (x, y, beta)

    ! === incompabible X, COEFS with ADD_INTERCEPT=.FALSE. ===
    nrhs = 3
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config(add_intercept=.false.)
    call ols (conf, x, y, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X, COEFS with ADD_INTERCEPT=.FALSE."

    deallocate (x, y, beta)

    ! === NOBS < NCOEFS ===
    nrhs = 3
    nobs = nrhs
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config(add_intercept=.true.)
    call ols (conf, x, y, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "NOBS < NCOEFS"

    deallocate (x, y, beta)

    ! === Invalid RCOND ===
    nrhs = 3
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config(rcond=-0.1_PREC)
    call ols (conf, x, y, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid RCOND"

    deallocate (x, y, beta)

end subroutine

subroutine test_OLS_2d_args ()
    !*  Unit tests for OLS argument checking

    real (PREC), dimension(:,:), allocatable :: x, y
    real (PREC), dimension(:,:), allocatable :: beta

    integer :: nrhs, nlhs, nobs, ncoefs
    type (lm_config) :: conf
    type (status_t) :: status

    ! === incompatible X, Y arguments ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs+1,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X, Y"

    deallocate (x, y, beta)

    ! === incompatible X^T, Y ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs, nobs), y(nobs+1,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x, y, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X^T, Y"

    deallocate (x, y, beta)

    ! === incompabible X, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(nrhs-1,nlhs))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X, COEFS"

    deallocate (x, y, beta)

    ! === incompabible X^T, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nrhs+2,nobs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x, y, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X^T, COEFS"

    deallocate (x, y, beta)

    ! === incompabible X, COEFS with ADD_INTERCEPT=.FALSE. ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X, COEFS with ADD_INTERCEPT=.FALSE."

    deallocate (x, y, beta)

    ! === incompabible Y, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs-1))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays Y, COEFS"

    deallocate (x, y, beta)

    ! === NOBS < NCOEFS ===
    nrhs = 3
    nlhs = 2
    nobs = nrhs
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    call ols (x=x, y=y, coefs=beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "NOBS < NCOEFS"

    deallocate (x, y, beta)

    ! === Invalid RCOND ===
    nrhs = 3
    nlhs = 2
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = lm_config (rcond=-1.0_PREC)
    call ols (conf, x, y, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid RCOND"

    deallocate (x, y, beta)

end subroutine

subroutine test_OLS_2d ()
    !*  Unit tests for various OLS regression problems using the 2d API

    real (PREC), dimension(:,:), allocatable :: x, y, beta, beta_ok, work
    integer :: nrhs, nlhs, ncoefs, nobs
    integer :: i
    integer :: rank
    real (PREC), parameter :: Rtol=1.0e-4_PREC, ATOL=1.0e-12_PREC
        ! Arguments for ALL_CLOSE()
    real (PREC) :: rsq
    logical :: all_ok, values_ok
    type (lm_config) :: conf
    type (lm_result), allocatable :: lm
    type (status_t) :: status

    call set_seed (123)

    ! === exactly determined system with NOBS = NRHS, no const ===
    nrhs = 10
    nobs = 10
    ncoefs = nrhs
    nlhs = 3

    allocate (x(nobs,nrhs), y(nobs,nlhs), beta(ncoefs,nlhs), beta_ok(ncoefs,nlhs))

    call random_number (x)
    call random_number (beta_ok)

    call BLAS95_GEMM (x, beta_ok, y)

    status = NF_STATUS_UNDEFINED
    allocate (lm)
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, beta, rank=rank, res=lm, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (beta, beta_ok, rtol=RTOL, &
        atol=ATOL))) &
        error stop "Exactly identified system: estimation"
    deallocate (lm)

    ! Test without optional COEFS argument
    status = NF_STATUS_UNDEFINED
    conf = lm_config (add_intercept=.false.)
    allocate (lm)
    call ols (conf, x, y, rank=rank, res=lm, status=status)
    values_ok = all_close (lm%coefs_multi, beta_ok, rtol=RTOL, atol=ATOL)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Exactly identified system: estimation w/o COEFS argument"

    deallocate (x, y, beta, beta_ok)
    deallocate (lm)

    ! === overdetermined system without intercept ===
    nrhs = 5
    nobs = 500
    ncoefs = nrhs
    nlhs = 2

    allocate (x(nobs,nrhs), y(nobs,nlhs), beta(ncoefs,nlhs), beta_ok(ncoefs,nlhs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    call BLAS95_GEMM (x, beta_ok, y)

    status = NF_STATUS_OK
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (beta, beta_ok, rtol=RTOL, &
        atol=ATOL))) &
        error stop "Underdetermined system"

    deallocate (x, y, beta, beta_ok)

    ! === overdetermined system with intercept ===
    nrhs = 5
    nobs = 500
    ncoefs = nrhs + 1
    nlhs = 11

    allocate (x(nobs,nrhs), y(nobs,nlhs), beta(ncoefs,nlhs), beta_ok(ncoefs,nlhs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    do i = 1, nlhs
        y(:,i) = beta_ok(1,i)
    end do
    allocate (work(nrhs,nlhs), source=beta_ok(2:,:))
    call BLAS95_GEMM (x, work, y, beta=1.0_PREC)
    deallocate (work)

    status = NF_STATUS_OK
    conf = lm_config (add_intercept=.true.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    ! Apply less strict check on beta, value of 1.0e-3 does not seem to work
    ! for NLHS=11
    if (.not. (status == NF_STATUS_OK .and. all_close (beta, beta_ok, rtol=&
        1.0e-2_PREC, atol=ATOL))) &
        error stop "Underdetermined system, ADD_INTERCEPT=.TRUE."

    deallocate (x, y, beta, beta_ok)

    ! === overdetermined system with intercept, transposed X ===
    nrhs = 3
    nobs = 500
    ncoefs = nrhs + 1
    nlhs = 7

    allocate (x(nrhs,nobs), y(nobs,nlhs), beta(ncoefs,nlhs), beta_ok(ncoefs,nlhs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    do i = 1, nlhs
        y(:,i) = beta_ok(1,i)
    end do
    allocate (work(nrhs,nlhs), source=beta_ok(2:,:))
    call BLAS95_GEMM (x, work, y, beta=1.0_PREC, transa="T")
    deallocate (work)

    status = NF_STATUS_OK
    conf = lm_config (add_intercept=.true., trans_x=.true.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (beta, beta_ok, rtol=&
        1.0e-2_PREC, atol=ATOL))) &
        error stop "Underdetermined system, ADD_INTERCEPT=.TRUE., TRANS_X=.TRUE."

    deallocate (x, y, beta, beta_ok)

end subroutine

subroutine test_OLS_1d ()
    !*  Unit tests for various OLS regression problems using the 1d API

    real (PREC), dimension(:,:), allocatable :: x
    real (PREC), dimension(:), allocatable :: y, beta, beta_ok, work
    integer :: nrhs, ncoefs, nobs
    integer :: rank
    real (PREC), parameter :: Rtol=1.0e-4_PREC, ATOL=1.0e-12_PREC
        ! Arguments for ALL_CLOSE()
    real (PREC) :: rsq
    type (status_t) :: status
    type (lm_config) :: conf
    type (lm_result), allocatable :: lm

    call set_seed (456)

    ! === exactly determined system with NOBS = NRHS, no const ===

    nrhs = 10
    nobs = 10
    ncoefs = nrhs

    call random_sample (nrhs, nobs, x, y, coefs=beta_ok, add_intercept=.false.)

    status = NF_STATUS_UNDEFINED
    allocate (lm)
    allocate (beta(ncoefs))
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, beta, rank=rank, res=lm, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (beta, beta_ok, rtol=RTOL, &
        atol=ATOL))) &
        error stop "Exactly identified system: estimation"

    ! Test w/o optional COEFS argument, coefs are stored in LM_DATA
    deallocate (lm)
    allocate (lm)
    status = NF_STATUS_UNDEFINED
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, rank=rank, res=lm, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (lm%coefs, beta_ok, rtol=&
        RTOL, atol=ATOL))) &
        error stop "Exactly identified system: estimation w/o COEFS argument"

    ! Compute R^2 which must be 1.0 in this case
    status = NF_STATUS_UNDEFINED
    call post_estim (lm, x, y, rsq=rsq, status=status)
    if (.not. (abs (rsq - 1.0_PREC) < 1.0e-12 .and. status == NF_STATUS_OK)) &
        error stop "Exactly identified system: R^2"

    deallocate (x, y, beta, beta_ok)

    ! === overdetermined system without intercept ===
    nrhs = 5
    nobs = 500
    ncoefs = nrhs

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_ok(ncoefs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    call BLAS95_GEMV (x, beta_ok, y)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (add_intercept=.false.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (beta, beta_ok, rtol=RTOL, &
        atol=ATOL))) &
        error stop "Overdetermined system"

    deallocate (x, y, beta, beta_ok)

    ! === overdetermined system with intercept ===

    nrhs = 5
    nobs = 500
    ncoefs = nrhs + 1

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_ok(ncoefs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    y(:) = beta_ok(1)
    allocate (work(nrhs), source=beta_ok(2:))
    call BLAS95_GEMV (x, work, y, beta=1.0_PREC)
    deallocate (work)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (add_intercept=.true.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (beta, beta_ok, rtol=RTOL, &
        atol=ATOL))) &
        error stop "Overdetermined system, ADD_INTERCEPT=.TRUE."

    deallocate (x, y, beta, beta_ok)

    ! === overdetermined system with intercept, transposed X ===
    nrhs = 2
    nobs = 500
    ncoefs = nrhs + 1

    allocate (x(nrhs,nobs), y(nobs), beta(ncoefs), beta_ok(ncoefs))

    call random_number (beta_ok)

    call rvs (norm, x, loc=0.0_PREC, scale=1.0_PREC)

    ! Create LHS variables
    y(:) = beta_ok(1)
    allocate (work(nrhs), source=beta_ok(2:))
    call BLAS95_GEMV (x, work, y, beta=1.0_PREC, trans="T")
    deallocate (work)

    status = NF_STATUS_OK
    conf = lm_config (add_intercept=.true., trans_x=.true.)
    call ols (conf, x, y, beta, rank=rank, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (beta, beta_ok, rtol=RTOL, &
        atol=ATOL))) &
        error stop "Overdetermined system, ADD_INTERCEPT=.TRUE., TRANS_X=.TRUE."

    deallocate (x, y, beta, beta_ok)

end subroutine

subroutine test_WLS ()
    !*  Unit tests for weightes LS.

    integer :: Nrhs, Nobs, Nlhs, Nadd
    real (PREC), dimension(:,:), allocatable :: X, Y, Xex, Yex
    real (PREC), dimension(:,:), allocatable :: coefs_true, coefs, coefs_w
    real (PREC), dimension(:), allocatable :: intercept_true, intercept, intercept_w
    real (PREC), dimension(:), allocatable :: weights
    type (lm_config) :: conf
    integer :: i, j
    logical :: values_ok
    type (status_t) :: status
    real (PREC), parameter :: Atol=1.0e-10_PREC, Rtol=1.0e-10_PREC

    call set_seed (1234)

    Nobs = 77
    Nlhs = 7
    Nrhs = 20

    call random_sample (nrhs=Nrhs, nobs=Nobs, nlhs=Nlhs, X=X, Y=Y, &
        coefs=coefs_true, intercept=intercept_true, add_intercept=.true., &
        var_error=0.1_PREC, status=status)

    allocate (weights(Nobs), source=1.0_PREC)

    ! --- No weights given vs. uniform weights ---

    allocate (coefs(Nrhs,Nlhs), coefs_w(Nrhs,Nlhs))
    allocate (intercept(Nlhs), intercept_w(Nlhs))

    ! Intercept-only model
    conf%add_intercept = .true.

    status = NF_STATUS_UNDEFINED
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    call ols (conf, X(:,1:0), Y, coefs(1:0,:), intercept, status=status)

    status = NF_STATUS_UNDEFINED
    coefs_w(:,:) = huge(0.0_PREC)
    intercept_w(:) = huge(0.0_PREC)
    call ols (conf, X(:,1:0), Y, coefs_w(1:0,:), intercept_w, weights=weights, &
        status=status)

    values_ok = all_close (intercept_w, intercept, atol=atol, rtol=rtol)

    if (.not. (values_ok .and. status == NF_STATUS_OK)) &
        error stop "Intercept only, uniform weights vs. no weights"

    ! Full model
    status = NF_STATUS_UNDEFINED
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    call ols (conf, X, Y, coefs, intercept, status=status)

    status = NF_STATUS_UNDEFINED
    coefs_w(:,:) = huge(0.0_PREC)
    intercept_w(:) = huge(0.0_PREC)
    call ols (conf, X, Y, coefs_w, intercept_w, weights=weights, &
        status=status)

    values_ok = all_close (intercept_w, intercept, atol=atol, rtol=rtol)

    if (.not. (values_ok .and. status == NF_STATUS_OK)) &
        error stop "Full model, uniform weights vs. no weights"

    ! --- Expanded vs. weighted sample ---

    Nadd = 10
    allocate (Xex(Nobs+Nadd,Nrhs), Yex(Nobs+Nadd,Nlhs))

    weights(:) = 1.0

    ! Expand sample by duplicating observations, increase weight
    Xex(1:Nobs,:) = X
    Yex(1:Nobs,:) = Y
    do i = 1, Nadd
        j = modulo (i * 5, Nobs-1) + 1
        weights(j) = weights(j) + 1.0
        Xex(Nobs+i,:) = X(j,:)
        Yex(Nobs+i,:) = Y(j,:)
    end do

    ! Intercept-only model
    conf%add_intercept = .true.

    status = NF_STATUS_UNDEFINED
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    call ols (conf, Xex(:,1:0), Yex, coefs(1:0,:), intercept, status=status)

    status = NF_STATUS_UNDEFINED
    coefs_w(:,:) = huge(0.0_PREC)
    intercept_w(:) = huge(0.0_PREC)
    call ols (conf, X(:,1:0), Y, coefs_w(1:0,:), intercept_w, weights=weights, &
        status=status)

    values_ok = all_close (intercept_w, intercept, atol=atol, rtol=rtol)

    if (.not. (values_ok .and. status == NF_STATUS_OK)) &
        error stop "Intercept only, weighted vs. expanded sample"

    ! Full model
    status = NF_STATUS_UNDEFINED
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    call ols (conf, Xex, Yex, coefs, intercept, status=status)

    status = NF_STATUS_UNDEFINED
    coefs_w(:,:) = huge(0.0_PREC)
    intercept_w(:) = huge(0.0_PREC)
    call ols (conf, X, Y, coefs_w, intercept_w, weights=weights, &
        status=status)

    values_ok = all_close (intercept_w, intercept, atol=atol, rtol=rtol)

    if (.not. (values_ok .and. status == NF_STATUS_OK)) &
        error stop "Full model, weighted vs. expanded sample"

end subroutine

subroutine test_OLS_fit_equiv ()
    !*  Unit tests for equivalent formulations of the same LS problem,
    !   with OLS() being called with non-transposed and transposed
    !   input and output data.

    type (lm_config) :: conf
    type (lm_result), allocatable :: res
    real (PREC), dimension(:,:), allocatable :: X, X_T, Y, Y_T
    real (PREC), dimension(:,:), allocatable :: coefs, coefs_T, coefs_true
    real (PREC), dimension(:), allocatable ::  intercept, intercept_true
    integer :: Nobs, Nrhs, Nconst, Ncoefs, Nlhs
    type (status_t) :: status
    logical :: values_ok
    real (PREC), parameter :: atol=1.0e-12_PREC, rtol=1.0e-4_PREC
    integer :: jlhs

    call set_seed (1234)

    ! Test data set

    Nobs = 11
    Nrhs = 3
    Nlhs = 5

    allocate (X(Nobs,Nrhs))
    allocate (coefs_true(Nrhs,Nlhs))
    allocate (intercept_true(Nlhs))

    call random_number (X)
    call random_number (coefs_true)
    call random_number (intercept_true)

    coefs_true(:,:) = (coefs_true - 0.5d0) * 5.0d0
    intercept_true(:) = (intercept_true - 0.5d0) * 10.d0

    allocate (Y(Nobs,Nlhs))

    Y(:,:) = spread (intercept_true, dim=1, ncopies=Nobs)

    call BLAS95_GEMM (X, coefs_true, Y, alpha=1.0_PREC, beta=1.0_PREC)

    allocate (X_T(Nrhs,Nobs), Y_T(Nlhs,Nobs))
    X_T(:,:) = transpose(X)
    Y_T(:,:) = transpose(Y)

    ! --- 1d API, no intercept in COEFS ---

    ! trans_x = .false.
    allocate (coefs(Nrhs,1), source=huge(0.0_PREC))
    allocate (intercept(1), source=huge(0.0_PREC))
    allocate (res)

    status = NF_STATUS_UNDEFINED
    call ols (x=X, y=Y(:,1), coefs=coefs(:,1), intercept=intercept(1), res=res, status=status)

    values_ok = all_close (coefs(:,1), coefs_true(:,1), atol=atol, rtol=rtol) &
        .and. all_close (intercept(1), intercept_true(1), atol=atol, rtol=rtol) &
        .and. all_close (res%coefs, coefs_true(:,1), atol=atol, rtol=rtol) &
        .and. all_close (res%intercept, intercept_true(1), atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_X=.FALSE., no intercept in coefs"

    deallocate (res)

    ! trans_x = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x=X_T, y=Y(:,1), coefs=coefs(:,1), intercept=intercept(1), &
        res=res, status=status)

    values_ok = all_close (coefs(:,1), coefs_true(:,1), atol=atol, rtol=rtol) &
        .and. all_close (intercept(1), intercept_true(1), atol=atol, rtol=rtol) &
        .and. all_close (res%coefs, coefs_true(:,1), atol=atol, rtol=rtol) &
        .and. all_close (res%intercept, intercept_true(1), atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_X=.TRUE., no intercept in coefs"

    deallocate (res)
    deallocate (coefs, intercept)

    ! --- 1d API, add intercept to COEFS ---

    ! trans_x = .false.
    Nconst = 1
    Ncoefs = Nrhs + Nconst
    allocate (coefs(Ncoefs,1), source=huge(0.0_PREC))
    allocate (intercept(1), source=huge(0.0_PREC))
    allocate (res)

    jlhs = 2
    status = NF_STATUS_UNDEFINED
    call ols (x=X, y=Y(:,jlhs), coefs=coefs(:,1), intercept=intercept(1), res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,1), coefs_true(:,jlhs), atol=atol, rtol=rtol) &
        .and. all_close (intercept(1), intercept_true(jlhs), atol=atol, rtol=rtol) &
        .and. all_close (res%coefs, coefs_true(:,jlhs), atol=atol, rtol=rtol) &
        .and. all_close (res%intercept, intercept_true(jlhs), atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_X=.FALSE., intercept in coefs"

    deallocate (res)

    ! trans_x = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x=X_T, y=Y(:,jlhs), coefs=coefs(:,1), intercept=intercept(1), &
        res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,1), coefs_true(:,jlhs), atol=atol, rtol=rtol) &
        .and. all_close (coefs(1,1), intercept_true(jlhs), atol=atol, rtol=rtol) &
        .and. all_close (intercept(1), intercept_true(jlhs), atol=atol, rtol=rtol) &
        .and. all_close (res%coefs, coefs_true(:,jlhs), atol=atol, rtol=rtol) &
        .and. all_close (res%intercept, intercept_true(jlhs), atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_X=.TRUE., intercept in coefs"

    deallocate (res)
    deallocate (coefs, intercept)

    ! --- 2d API, no intercept in COEFS ---

    ! trans_x = .false.
    allocate (coefs(Nrhs,Nlhs), source=huge(0.0_PREC))
    allocate (intercept(Nlhs), source=huge(0.0_PREC))
    allocate (res)

    status = NF_STATUS_UNDEFINED
    call ols (x=X, y=Y, coefs=coefs, intercept=intercept, res=res, status=status)

    values_ok = all_close (coefs, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d API: TRANS_X=.FALSE., no intercept in coefs"

    deallocate (res)

    ! trans_x = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x=X_T, y=Y, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_X=.TRUE., no intercept in coefs"

    deallocate (res)

    ! trans_x = .false., trans_y = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.false., trans_y=.true.)
    call ols (conf, x=X, y=Y_T, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_Y=.TRUE., no intercept in coefs"

    deallocate (res)

    ! trans_x = .true., trans_y = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true., trans_y=.true.)
    call ols (conf, x=X_T, y=Y_T, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_X=.TRUE., TRANS_Y=.TRUE., no intercept in coefs"

    deallocate (res)
    deallocate (coefs, intercept)

    ! trans_coefs = .true.
    allocate (coefs_T(Nlhs,Nrhs), source=huge(0.0_PREC))
    allocate (intercept(Nlhs), source=huge(0.0_PREC))
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config(trans_coefs=.true.)
    call ols (conf, x=X, y=Y, coefs=coefs_T, intercept=intercept, res=res, status=status)

    values_ok = all_close (transpose(coefs_T), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d API: TRANS_COEFS=.TRUE., no intercept in coefs"

    deallocate (coefs_T, intercept)
    deallocate (res)

    ! --- 2d API, add intercept to COEFS ---

    ! trans_x = .false.
    Nconst = 1
    Ncoefs = Nrhs + Nconst
    allocate (coefs(Ncoefs,Nlhs), source=huge(0.0_PREC))
    allocate (intercept(Nlhs), source=huge(0.0_PREC))
    allocate (res)

    status = NF_STATUS_UNDEFINED
    call ols (x=X, y=Y, coefs=coefs, intercept=intercept, res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,:), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (coefs(1,:), intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_X=.FALSE., intercept in coefs"

    deallocate (res)

    ! trans_x = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true.)
    call ols (conf, x=X_T, y=Y, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,:), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (coefs(1,:), intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_X=.TRUE., intercept in coefs"

    deallocate (res)

    ! trans_y = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.false., trans_y=.true.)
    call ols (conf, x=X, y=Y_T, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,:), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (coefs(1,:), intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_Y=.TRUE., intercept in coefs"

    deallocate (res)

    ! trans_x = .true., trans_y = .true.
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_x=.true., trans_y=.true.)
    call ols (conf, x=X_T, y=Y_T, coefs=coefs, intercept=intercept, &
        res=res, status=status)

    values_ok = all_close (coefs(1+Nconst:Ncoefs,:), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (coefs(1,:), intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_X=.TRUE., TRANS_Y=.TRUE., intercept in coefs"

    deallocate (res)
    deallocate (coefs, intercept)

    ! trans_coefs = .true.
    Nconst = 1
    Ncoefs = Nrhs + Nconst
    allocate (coefs_T(Nlhs,Ncoefs), source=huge(0.0_PREC))
    allocate (intercept(Nlhs), source=huge(0.0_PREC))
    allocate (res)

    status = NF_STATUS_UNDEFINED
    conf = lm_config (trans_coefs=.true.)
    call ols (conf, x=X, y=Y, coefs=coefs_T, intercept=intercept, res=res, status=status)

    values_ok = all_close (transpose(coefs_T(:,1+Nconst:Ncoefs)), coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (coefs_T(:,1), intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (intercept, intercept_true, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_true, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "1d API: TRANS_COEFS=.TRUE., intercept in coefs"

    deallocate (coefs_T, intercept)
    deallocate (res)

end subroutine

subroutine test_predict_equiv_1d ()
    !*  Test equivalent calls to PREDICT that use transposed variations of
    !   input and output data

    type (lm_config) :: conf
    type (lm_result) :: res
    real (PREC), dimension(:,:), allocatable :: X, Xpred, Xpred_T
    real (PREC), dimension(:), allocatable :: Y, Ypred
    real (PREC), dimension(:), allocatable :: Ypred_true
    real (PREC), dimension(:), allocatable :: coefs
    real (PREC), dimension(:), allocatable :: coefs_intercept
    real (PREC) ::  intercept
    integer :: Nobs, Nrhs, Nconst, Ncoefs, Nlhs, Nobs_pred
    type (status_t) :: status
    logical :: values_ok
    real (PREC), parameter :: atol=1.0e-12_PREC, rtol=1.0e-4_PREC
    integer :: jlhs
    logical :: trans_x, trans_y, trans_coefs

    call set_seed (45678)

    ! --- Training data set ---

    Nobs = 11
    Nrhs = 3
    Nconst = 1
    Ncoefs = Nrhs + Nconst

    allocate (X(Nobs,Nrhs))
    allocate (coefs(Nrhs))

    call random_number (X)
    call random_number (coefs)
    call random_number (intercept)

    coefs(:) = (coefs - 0.5d0) * 5.0d0
    intercept = (intercept - 0.5d0) * 10.d0

    allocate (Y(Nobs))

    Y(:) = intercept

    call BLAS95_GEMV (X, coefs, Y, alpha=1.0_PREC, beta=1.0_PREC)

    ! Coefs + intercept in one array
    allocate (coefs_intercept(Ncoefs))
    coefs_intercept(1) = intercept
    coefs_intercept(1+Nconst:Ncoefs) = coefs

    ! Create OLS result objects so we can test PREDICT(LM_RESULT,...)
    call ols (X=X, Y=Y, res=res, status=status)

    ! --- Test data ---
    Nobs_pred = Nobs + 12

    allocate (Xpred(Nobs_pred,Nrhs), Xpred_T(Nrhs,Nobs_pred))
    allocate (Ypred_true(Nobs_pred))

    call random_number (Xpred)
    Xpred_T(:,:) = transpose (Xpred)

    Ypred_true(:) = intercept
    call BLAS95_GEMV (Xpred, coefs, Ypred_true, alpha=1.0_PREC, beta=1.0_PREC)

    ! Arrays where PREDICT writes its results
    allocate (Ypred(Nobs_pred))

    ! trans_x = .false.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred, coefs, Ypred, intercept, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) error stop "TRANS_X=.FALSE."

    ! trans_x = .false., called via PREDICT(LM_RESULT,...)
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    Ypred(:) = huge(0.0_PREC)
    call predict (res, Xpred, Ypred, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "TRANS_X=.FALSE., using LM_RESULT"

    ! trans_x = .false., extended COEFS, INTERCEPT present
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred, &
        intercept=intercept, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "TRANS_X=.FALSE., extended COEFS, INTERCEPT present"

    ! trans_x = .false., extended COEFS, INTERCEPT not present
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "TRANS_X=.FALSE., extended COEFS"

    ! trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x = .true.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs, Ypred, intercept, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) error stop "TRANS_X=.FALSE."

    ! trans_x = .true., using LM_RESULT
    status = NF_STATUS_UNDEFINED
    trans_x = .true.
    Ypred(:) = huge(0.0_PREC)
    call predict (res, Xpred_T, Ypred, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "TRANS_X=.FALSE., using LM_RESULT"

    ! trans_x = .true., extended COEFS, INTERCEPT present
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred, &
        intercept=intercept, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "TRANS_X=.TRUE., extended COEFS, INTERCEPT present"

    ! trans_x = .true., extended COEFS, INTERCEPT not present
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    Ypred(:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred, trans_x=trans_x, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "TRANS_X=.TRUE., extended COEFS"

end subroutine

subroutine test_predict_equiv_2d ()
    !*  Test equivalent calls to PREDICT that use transposed variations of
    !   input and output data

    type (lm_config) :: conf
    type (lm_result) :: res
    real (PREC), dimension(:,:), allocatable :: X, Xpred, Xpred_T
    real (PREC), dimension(:,:), allocatable :: Y, Ypred, Ypred_T
    real (PREC), dimension(:,:), allocatable :: Ypred_true, Ypred_true_T
    real (PREC), dimension(:,:), allocatable :: coefs, coefs_T
    real (PREC), dimension(:,:), allocatable :: coefs_intercept, coefs_intercept_T
    real (PREC), dimension(:), allocatable ::  intercept
    integer :: Nobs, Nrhs, Nconst, Ncoefs, Nlhs, Nobs_pred
    type (status_t) :: status
    logical :: values_ok
    real (PREC), parameter :: atol=1.0e-12_PREC, rtol=1.0e-4_PREC
    logical :: trans_x, trans_y, trans_coefs

    call set_seed (7890)

    ! --- Training data set ---

    Nobs = 11
    Nrhs = 3
    Nlhs = 5
    Nconst = 1
    Ncoefs = Nrhs + Nconst

    allocate (X(Nobs,Nrhs))
    allocate (coefs(Nrhs,Nlhs))
    allocate (intercept(Nlhs))

    call random_number (X)
    call random_number (coefs)
    call random_number (intercept)

    coefs(:,:) = (coefs - 0.5d0) * 5.0d0
    intercept(:) = (intercept - 0.5d0) * 10.d0

    allocate (Y(Nobs,Nlhs))

    Y(:,:) = spread (intercept, dim=1, ncopies=Nobs)

    call BLAS95_GEMM (X, coefs, Y, alpha=1.0_PREC, beta=1.0_PREC)

    allocate (coefs_T(Nlhs,Nrhs))
    coefs_T(:,:) = transpose(coefs)

    ! Coefs + intercept in one array
    allocate (coefs_intercept(Ncoefs,Nlhs))
    coefs_intercept(1,:) = intercept
    coefs_intercept(1+Nconst:Ncoefs,:) = coefs

    allocate (coefs_intercept_T(Nlhs,Ncoefs))
    coefs_intercept_T(:,:) = transpose (coefs_intercept)

    ! Create OLS result objects so we can test PREDICT(LM_RESULT,...)
    call ols (X=X, Y=Y, res=res, status=status)

    ! --- Test data ---
    Nobs_pred = Nobs + 12

    allocate (Xpred(Nobs_pred,Nrhs), Xpred_T(Nrhs,Nobs_pred))
    allocate (Ypred_true(Nobs_pred,Nlhs), Ypred_true_T(Nlhs,Nobs_pred))

    call random_number (Xpred)
    Xpred_T(:,:) = transpose (Xpred)

    Ypred_true(:,:) = spread (intercept, dim=1, ncopies=Nobs_pred)
    call BLAS95_GEMM (Xpred, coefs, Ypred_true, alpha=1.0_PREC, beta=1.0_PREC)

    Ypred_true_T(:,:) = transpose (Ypred_true)

    ! Arrays where PREDICT writes its results
    allocate (Ypred(Nobs_pred,Nlhs), Ypred_T(Nlhs,Nobs_pred))

    ! --- Tests using LM_RESULT ---

    ! trans_x = .false.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (res, Xpred, Ypred, trans_x=trans_x, trans_y=trans_y, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "LM_RESULT: TRANS_X=.FALSE."

    ! trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x = .true.
    trans_y = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (res, Xpred_T, Ypred, trans_x=trans_x, trans_y=trans_y, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "LM_RESULT: TRANS_X=.TRUE."

    ! trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (res, Xpred, Ypred_T, trans_x=trans_x, trans_y=trans_y, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "LM_RESULT: TRANS_Y=.TRUE."

    ! trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x = .true.
    trans_y = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (res, Xpred_T, Ypred_T, trans_x=trans_x, trans_y=trans_y, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "LM_RESULT: TRANS_X=.TRUE., TRANS_Y=.TRUE."

    ! --- Test using COEFS and INTERCEPT ---

    ! trans_x = .false.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "COEFS: TRANS_X=.FALSE."

    ! trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "COEFS: TRANS_X=.TRUE."

    ! trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "COEFS: TRANS_Y=.TRUE."

    ! trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "COEFS: TRANS_X=.TRUE., TRANS_Y=.TRUE."

    ! trans_coefs = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_T, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "COEFS: TRANS_COEFS=.TRUE."

    ! trans_coefs = .true., trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_T, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "COEFS: TRANS_COEFS=.TRUE., TRANS_X=.TRUE."

    ! trans_coefs = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_T, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "COEFS: TRANS_COEFS=.TRUE., TRANS_Y=.TRUE."

    ! trans_coefs = .true., trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_T, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "COEFS: TRANS_COEFS=.TRUE., TRANS_X=.TRUE., TRANS_Y=.TRUE."

    ! --- Extended COEFS, INTERCEPT present ---

    ! trans_x = .false.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Ext. COEFS, INTERCEPT: TRANS_X=.FALSE."

    ! trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Ext. COEFS, INTERCEPT: TRANS_X=.TRUE."

    ! trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Ext. COEFS, INTERCEPT: TRANS_Y=.TRUE."

    ! trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Ext. COEFS, INTERCEPT: TRANS_X=.TRUE., TRANS_Y=.TRUE."

    ! trans_coefs = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept_T, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Ext. COEFS, INTERCEPT: TRANS_COEFS=.TRUE."

    ! trans_coefs = .true., trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept_T, Ypred, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Ext. COEFS, INTERCEPT: TRANS_COEFS=.TRUE., TRANS_X=.TRUE."

    ! trans_coefs = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept_T, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Ext. COEFS, INTERCEPT: TRANS_COEFS=.TRUE., TRANS_Y=.TRUE."

    ! trans_coefs = .true., trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept_T, Ypred_T, intercept, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Ext. COEFS, INTERCEPT: TRANS_COEFS=.TRUE., TRANS_X=.TRUE., " // &
            "TRANS_Y=.TRUE."

    ! --- Extended COEFS, INTERCEPT not present ---

    ! trans_x = .false.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Extended COEFS: TRANS_X=.FALSE."

    ! trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Extended COEFS: TRANS_X=.TRUE."

    ! trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept, Ypred_T, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Extended COEFS: TRANS_Y=.TRUE."

    ! trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .false.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept, Ypred_T, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Extended COEFS: TRANS_X=.TRUE., TRANS_Y=.TRUE."

    ! trans_coefs = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept_T, Ypred, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Extended COEFS: TRANS_COEFS=.TRUE."

    ! trans_coefs = .true., trans_x = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .false.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept_T, Ypred, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (Ypred, Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Extended COEFS: TRANS_COEFS=.TRUE., TRANS_X=.TRUE."

    ! trans_coefs = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.false.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred, coefs_intercept_T, Ypred_T, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Extended COEFS: TRANS_COEFS=.TRUE., TRANS_Y=.TRUE."

    ! trans_coefs = .true., trans_x = .true., trans_y = .true.
    status = NF_STATUS_UNDEFINED
    trans_x =.true.
    trans_y = .true.
    trans_coefs = .true.
    Ypred(:,:) = huge(0.0_PREC)
    call predict (Xpred_T, coefs_intercept_T, Ypred_T, trans_x=trans_x, &
        trans_y=trans_y, trans_coefs=trans_coefs, status=status)
    values_ok = all_close (transpose(Ypred_T), Ypred_true, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Extended COEFS: TRANS_COEFS=.TRUE., TRANS_X=.TRUE., " // &
            "TRANS_Y=.TRUE."

end subroutine

subroutine test_pcr_ols_equiv_1d ()
    !*  Test models where OLS and PCR yield equivalent results

    type (lm_config) :: conf_pcr, conf_ols
    type (lm_result) :: res_ols, res_pcr
    real (PREC), dimension(:,:), allocatable :: X
    real (PREC), dimension(:), allocatable :: y
    real (PREC), dimension(:), allocatable :: coefs_true, coefs_ols, coefs_pcr
    real (PREC), dimension(:), allocatable :: weights
    real (PREC) ::  intercept_true, intercept_ols, intercept_pcr
    integer :: Nrhs, Nconst, Nobs
    type (status_t) :: status
    logical :: values_ok
    real (PREC), parameter :: atol=1.0e-12_PREC, rtol=1.0e-4_PREC
    integer :: i

    call set_seed (678)

    ! --- Test constant-only model ---

    Nobs = 10

    call random_sample (nrhs=0, nobs=Nobs, X=X, y=y, coefs=coefs_true, &
        intercept=intercept_true, add_intercept=.true., &
        var_error=1.0_PREC, status=status)

    allocate (weights(Nobs))
    call random_number (weights)

    ! Unweighted
    call ols (X=X, Y=Y, res=res_ols, status=status)
    call pcr (X=X, Y=Y, res=res_pcr, status=status)

    values_ok = all_close (res_pcr%intercept, res_ols%intercept, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Constant only, RES present"

    ! Weighted
    call ols (X=X, Y=Y, weights=weights, res=res_ols, status=status)
    call pcr (X=X, Y=Y, weights=weights, res=res_pcr, status=status)

    values_ok = all_close (res_pcr%intercept, res_ols%intercept, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Constant only, weighted, RES present"

    deallocate (X, Y, weights)

    ! --- Test for various number of RHS vars ---

    Nobs = 13

    do i = 1, 10
        Nrhs = i
        call random_sample (nrhs=Nrhs, nobs=Nobs+Nrhs, X=X, y=y, coefs=coefs_true, &
            intercept=intercept_true, add_intercept=.true., &
            var_error=0.1_PREC, status=status)

        allocate (weights(Nobs+Nrhs))
        call random_number (weights)

        ! Force using all available components to get equivalence with OLS
        conf_pcr = lm_config (var_rhs_min=1.0_PREC)

        call ols (X=X, y=y, res=res_ols, status=status)
        status = NF_STATUS_UNDEFINED
        call lm_result_reset (res_pcr)
        call pcr (conf_pcr, X=X, y=y, res=res_pcr, status=status)

        values_ok = all_close (res_pcr%intercept, res_ols%intercept, atol=atol, rtol=rtol) &
            .and. all_close (res_pcr%coefs, res_ols%coefs, atol=atol, rtol=rtol)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS through the result object"

        ! Test with weights
        call lm_result_reset (res_ols)
        call ols (X=X, y=y, weights=weights, res=res_ols, status=status)
        status = NF_STATUS_UNDEFINED
        call lm_result_reset (res_pcr)
        call pcr (conf_pcr, X=X, y=y, weights=weights, res=res_pcr, status=status)

        values_ok = all_close (res_pcr%intercept, res_ols%intercept, atol=atol, rtol=rtol) &
            .and. all_close (res_pcr%coefs, res_ols%coefs, atol=atol, rtol=rtol)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "Weighted PCR disagrees with OLS through the result object"

        ! Check with COEFS, INTERCEPT arguments
        call cond_alloc (coefs_ols, Nrhs)
        call ols (conf_pcr, X=X, y=y, coefs=coefs_ols, intercept=intercept_ols, &
            status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, Nrhs)
        coefs_pcr(:) = huge(0.0_PREC)
        intercept_pcr = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, intercept=intercept_pcr, &
            status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS using coefficients and intercept"

        ! Check with extended COEFS, INTERCEPT arguments
        Nconst = 1
        call cond_alloc (coefs_ols, Nrhs+Nconst)
        call ols (X=X, y=y, coefs=coefs_ols, intercept=intercept_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, Nrhs+Nconst)
        coefs_pcr(:) = huge(0.0_PREC)
        intercept_pcr = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, intercept=intercept_pcr, &
            status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS using extended coefficients and " // &
                "intercept"

        ! Check with extended COEFS, no INTERCEPT argument
        Nconst = 1
        call cond_alloc (coefs_ols, Nrhs+Nconst)
        call ols (X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, Nrhs+Nconst)
        coefs_pcr(:) = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS using extended coefficients"

        ! Check with no ADDED intercept
        call cond_alloc (coefs_ols, Nrhs)
        conf_ols = lm_config (add_intercept=.false.)
        call ols (conf_ols, X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, Nrhs)
        coefs_pcr(:) = huge(0.0_PREC)
        conf_pcr = lm_config (var_rhs_min=1.0_PREC, add_intercept=.false., center=.false.)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (coefs_pcr, coefs_ols, atol=atol, rtol=rtol)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS without an added intercept"

        ! Check with no added intercept, but intercept in X
        Nconst = 1
        call random_sample (nrhs=Nrhs, nobs=Nobs+Nrhs, X=X, y=y, coefs=coefs_true, &
            add_intercept_x=.true., var_error=0.1_PREC, status=status)

        call cond_alloc (coefs_ols, Nrhs + Nconst)
        conf_ols = lm_config (add_intercept=.false.)
        call ols (conf_ols, X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, Nrhs + Nconst)
        coefs_pcr(:) = huge(0.0_PREC)
        conf_pcr = lm_config (var_rhs_min=1.0_PREC, add_intercept=.false., center=.false.)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (coefs_pcr, coefs_ols, atol=atol, rtol=rtol)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS with an intercept column"

        deallocate (weights)

    end do

end subroutine

subroutine test_pcr_ols_equiv_2d ()
    !*  Test models where OLS and PCR yield equivalent results

    type (lm_config) :: conf_pcr, conf_ols
    type (lm_result) :: res_ols, res_pcr
    real (PREC), dimension(:,:), allocatable :: X, Y
    real (PREC), dimension(:,:), allocatable :: coefs_true, coefs_ols, coefs_pcr
    real (PREC), dimension(:), allocatable ::  intercept_true, intercept_ols, intercept_pcr
    real (PREC), dimension(:), allocatable :: weights
    integer :: Nobs, Nrhs, Nlhs, Nconst
    type (status_t) :: status
    logical :: values_ok
    real (PREC), parameter :: atol=1.0e-12_PREC, rtol=1.0e-4_PREC
    integer :: i

    call set_seed (678)

    ! --- Test constant-only model ---

    call random_sample (nrhs=0, nobs=10, Nlhs=3, X=X, y=y, coefs=coefs_true, &
        intercept=intercept_true, add_intercept=.true., &
        var_error=1.0_PREC, status=status)

    call ols (X=X, Y=Y, res=res_ols, status=status)
    status = NF_STATUS_OK
    call pcr (X=X, Y=Y, res=res_pcr, status=status)

    values_ok = all_close (res_pcr%intercept_multi, res_ols%intercept_multi, &
        atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Constant only, RES present"

    allocate (weights(10))
    call random_number (weights)

    call ols (X=X, Y=Y, weights=weights, res=res_ols, status=status)
    status = NF_STATUS_OK
    call pcr (X=X, Y=Y, weights=weights, res=res_pcr, status=status)

    values_ok = all_close (res_pcr%intercept_multi, res_ols%intercept_multi, &
        atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Constant only, weighted, RES present"

    deallocate (X, Y, weights)

    ! --- Test for various number of RHS vars ---

    do i = 1, 10, 3
        Nrhs = i
        Nlhs = modulo (i, 3) + 1
        Nobs = 43 + 2*i
        call random_sample (nrhs=Nrhs, nobs=Nobs, Nlhs=Nlhs, X=X, y=y, &
            coefs=coefs_true, intercept=intercept_true, add_intercept=.true., &
            var_error=0.1_PREC, status=status)

        allocate (weights(Nobs))
        call random_number (weights)

        call cond_alloc (intercept_ols, Nlhs)
        call cond_alloc (intercept_pcr, Nlhs)

        ! Force using all available components to get equivalence with OLS
        conf_pcr = lm_config (var_rhs_min=1.0_PREC)

        call ols (X=X, y=y, res=res_ols, status=status)
        status = NF_STATUS_UNDEFINED
        call lm_result_reset (res_pcr)
        call pcr (conf_pcr, X=X, y=y, res=res_pcr, status=status)

        values_ok = all_close (res_pcr%intercept_multi, res_ols%intercept_multi, &
                atol=atol, rtol=rtol) &
            .and. all_close (res_pcr%coefs_multi, res_ols%coefs_multi, &
                atol=atol, rtol=rtol)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS through the result object"

        ! Test with weights
        call lm_result_reset (res_ols)
        call ols (X=X, y=y, weights=weights, res=res_ols, status=status)
        status = NF_STATUS_UNDEFINED
        call lm_result_reset (res_pcr)
        call pcr (conf_pcr, X=X, y=y, weights=weights, res=res_pcr, status=status)

        values_ok = all_close (res_pcr%intercept_multi, res_ols%intercept_multi, &
                atol=atol, rtol=rtol) &
            .and. all_close (res_pcr%coefs_multi, res_ols%coefs_multi, &
                atol=atol, rtol=rtol)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "Weighted PCR disagrees with OLS through the result object"

        deallocate (weights)

        ! Check with COEFS, INTERCEPT arguments
        call cond_alloc (coefs_ols, [Nrhs,Nlhs])
        call ols (conf_pcr, X=X, y=y, coefs=coefs_ols, intercept=intercept_ols, &
            status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, [Nrhs,Nlhs])
        coefs_pcr(:,:) = huge(0.0_PREC)
        intercept_pcr(:) = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, intercept=intercept_pcr, &
            status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS using coefficients and intercept"

        ! Check with extended COEFS, INTERCEPT arguments
        Nconst = 1
        call cond_alloc (coefs_ols, [Nrhs+Nconst,Nlhs])
        call ols (X=X, y=y, coefs=coefs_ols, intercept=intercept_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, [Nrhs+Nconst,Nlhs])
        coefs_pcr(:,:) = huge(0.0_PREC)
        intercept_pcr(:) = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, intercept=intercept_pcr, &
            status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS using extended coefficients and " // &
                "intercept"

        ! Check with extended COEFS, no INTERCEPT argument
        Nconst = 1
        call cond_alloc (coefs_ols, [Nrhs+Nconst,Nlhs])
        call ols (X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, [Nrhs+Nconst,Nlhs])
        coefs_pcr(:,:) = huge(0.0_PREC)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (intercept_pcr, intercept_ols, atol=atol, rtol=rtol) &
            .and. all_close (coefs_pcr, coefs_ols)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS using extended coefficients"

        ! Check with no ADDED intercept
        call cond_alloc (coefs_ols, [Nrhs,Nlhs])
        conf_ols = lm_config (add_intercept=.false.)
        call ols (conf_ols, X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, [Nrhs,Nlhs])
        coefs_pcr(:,:) = huge(0.0_PREC)
        conf_pcr = lm_config (var_rhs_min=1.0_PREC, add_intercept=.false., center=.false.)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (coefs_pcr, coefs_ols, atol=atol, rtol=rtol)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS without an added intercept"

        ! Check with no added intercept, but intercept in X
        Nconst = 1
        call random_sample (nrhs=Nrhs, nobs=Nobs+Nrhs, nlhs=Nlhs, X=X, y=y, &
            coefs=coefs_true, add_intercept_x=.true., var_error=0.1_PREC, &
            status=status)

        call cond_alloc (coefs_ols, [Nrhs + Nconst,Nlhs])
        conf_ols = lm_config (add_intercept=.false.)
        call ols (conf_ols, X=X, y=y, coefs=coefs_ols, status=status)

        status = NF_STATUS_UNDEFINED
        call cond_alloc (coefs_pcr, [Nrhs + Nconst,Nlhs])
        coefs_pcr(:,:) = huge(0.0_PREC)
        conf_pcr = lm_config (var_rhs_min=1.0_PREC, add_intercept=.false., center=.false.)
        call pcr (conf_pcr, X=X, y=y, coefs=coefs_pcr, status=status)

        values_ok = all_close (coefs_pcr, coefs_ols, atol=atol, rtol=rtol)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "PCR disagrees with OLS with an intercept column"

    end do

end subroutine

subroutine test_pcr_equiv_1d ()
    !*  Test variations of equivalent calls to PCR

    real (PREC), dimension(:,:), allocatable :: X, X_T
    real (PREC), dimension(:), allocatable :: coefs, coefs_ok, coefs_true
    real (PREC), dimension(:), allocatable :: coefs_intercept_ok, coefs_intercept
    real (PREC), dimension(:), allocatable :: y
    real (PREC) :: intercept, intercept_ok, intercept_true
    integer :: Nrhs, Nobs, Ncomp
    real (PREC), parameter :: atol=1.0e-12_PREC, rtol=1.0e-4_PREC
    logical :: values_ok
    type (status_t) :: status
    type (lm_result) :: res_ok, res
    type (lm_config) :: conf

    call set_seed (2345)

    do Nrhs = 1, 5

        Nobs = Nrhs + 4
        call random_sample (nrhs=Nrhs, nobs=Nobs, X=X, y=y, &
            coefs=coefs_true, intercept=intercept_true, add_intercept=.true., &
            var_error=0.1_PREC, status=status)

        call cond_alloc (X_T, [Nrhs, Nobs])
        X_T(:,:) = transpose (X)

        call cond_alloc (coefs_ok, size(coefs_true))
        call cond_alloc (coefs, size(coefs_true))
        call cond_alloc (coefs_intercept_ok, size(coefs_true) + 1)
        call cond_alloc (coefs_intercept, size(coefs_true) + 1)

        do Ncomp = 0, Nrhs
            conf = lm_config (ncomp=Ncomp)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X, y, coefs=coefs_ok, intercept=intercept_ok, status=status)
            call pcr (conf, X, y, coefs=coefs_intercept_ok, status=status)
            call pcr (conf, X, y, res=res_ok, status=status)

            status = NF_STATUS_UNDEFINED
            conf = lm_config (ncomp=Ncomp, trans_x=.true.)
            call pcr (conf, X_T, y, coefs=coefs, intercept=intercept, status=status)

            values_ok = all_close (coefs, coefs_ok, atol=atol, rtol=rtol) .and. &
                all_close (intercept, intercept_ok, atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR coefficients and intercept disagree with " // &
                    "transposed X"

            status = NF_STATUS_UNDEFINED
            call lm_result_reset (res)
            call pcr (conf, X_T, y, res=res, status=status)

            values_ok = all_close (res%coefs, res_ok%coefs, atol=atol, rtol=rtol) .and. &
                all_close (res%intercept, res_ok%intercept, atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR result disagree with transposed X"

            call lm_result_reset (res)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X_T, y, coefs=coefs_intercept, status=status)

            values_ok = all_close (coefs_intercept, coefs_intercept_ok, &
                atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR extended coefficients disagree with transposed X"
        end do
    end do

end subroutine

subroutine test_pcr_equiv_2d ()
    !*  Test variations of equivalent calls to PCR

    real (PREC), dimension(:,:), allocatable :: X, X_T
    real (PREC), dimension(:,:), allocatable :: coefs, coefs_T, coefs_ok, coefs_true
    real (PREC), dimension(:,:), allocatable :: coefs_intercept_ok, &
        coefs_intercept, coefs_intercept_T
    real (PREC), dimension(:,:), allocatable :: Y, Y_T
    real (PREC), dimension(:), allocatable :: intercept, intercept_ok, intercept_true
    integer :: Nrhs, Nobs, Ncomp, Nlhs
    real (PREC), parameter :: atol=1.0e-12_PREC, rtol=1.0e-4_PREC
    logical :: values_ok
    type (status_t) :: status
    type (lm_result) :: res_ok, res
    type (lm_config) :: conf

    call set_seed (2345)

    do Nrhs = 1, 15, 3

        Nlhs = modulo (Nrhs, 5) + 1

        Nobs = Nrhs + 4
        call random_sample (nrhs=Nrhs, nobs=Nobs, nlhs=Nlhs, X=X, y=y, &
            coefs=coefs_true, intercept=intercept_true, add_intercept=.true., &
            var_error=0.1_PREC, status=status)

        call cond_alloc (X_T, [Nrhs, Nobs])
        X_T(:,:) = transpose (X)

        call cond_alloc (Y_T, [Nlhs, Nobs])
        Y_T(:,:) = transpose (Y)

        call cond_alloc (coefs_ok, [Nrhs, Nlhs])
        call cond_alloc (coefs, [Nrhs, Nlhs])
        call cond_alloc (coefs_intercept_ok, [Nrhs + 1, Nlhs])
        call cond_alloc (coefs_intercept, [Nrhs + 1, Nlhs])
        call cond_alloc (coefs_T, [Nlhs, Nrhs])
        call cond_alloc (coefs_intercept_T, [Nlhs, Nrhs + 1])
        call cond_alloc (intercept_ok, Nlhs)
        call cond_alloc (intercept, Nlhs)

        do Ncomp = 0, Nrhs
            conf = lm_config (ncomp=Ncomp)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X, y, coefs=coefs_ok, intercept=intercept_ok, status=status)
            call pcr (conf, X, y, coefs=coefs_intercept_ok, status=status)
            call pcr (conf, X, y, res=res_ok, status=status)

            ! --- Transpose X ---

            status = NF_STATUS_UNDEFINED
            conf = lm_config (ncomp=Ncomp, trans_x=.true.)
            call pcr (conf, X_T, y, coefs=coefs, intercept=intercept, status=status)

            values_ok = all_close (coefs, coefs_ok, atol=atol, rtol=rtol) .and. &
                all_close (intercept, intercept_ok, atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR coefficients and intercept disagree with " // &
                    "transposed X"

            status = NF_STATUS_UNDEFINED
            call lm_result_reset (res)
            call pcr (conf, X_T, y, res=res, status=status)

            values_ok = all_close (res%coefs_multi, res_ok%coefs_multi, &
                    atol=atol, rtol=rtol) .and. &
                all_close (res%intercept_multi, res_ok%intercept_multi, &
                    atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR result disagree with transposed X"

            call lm_result_reset (res)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X_T, y, coefs=coefs_intercept, status=status)

            values_ok = all_close (coefs_intercept, coefs_intercept_ok, &
                atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR extended coefficients disagree with transposed X"

            ! --- Transpose Y ---

            status = NF_STATUS_UNDEFINED
            conf = lm_config (ncomp=Ncomp, trans_y=.true.)
            call pcr (conf, X, Y_T, coefs=coefs, intercept=intercept, status=status)

            values_ok = all_close (coefs, coefs_ok, atol=atol, rtol=rtol) .and. &
                all_close (intercept, intercept_ok, atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR coefficients and intercept disagree with " // &
                    "transposed Y"

            status = NF_STATUS_UNDEFINED
            call lm_result_reset (res)
            call pcr (conf, X, Y_T, res=res, status=status)

            values_ok = all_close (res%coefs_multi, res_ok%coefs_multi, &
                atol=atol, rtol=rtol) .and. &
                all_close (res%intercept_multi, res_ok%intercept_multi, &
                    atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR result disagree with transposed Y"

            call lm_result_reset (res)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X, Y_T, coefs=coefs_intercept, status=status)

            values_ok = all_close (coefs_intercept, coefs_intercept_ok, &
                atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR extended coefficients disagree with transposed Y"

            ! --- Transpose X and Y ---

            status = NF_STATUS_UNDEFINED
            call lm_result_reset (res)
            conf = lm_config (ncomp=Ncomp, trans_x=.true., trans_y=.true.)
            call pcr (conf, X_T, Y_T, res=res, status=status)

            values_ok = all_close (res%coefs_multi, res_ok%coefs_multi, &
                atol=atol, rtol=rtol) .and. &
                all_close (res%intercept_multi, res_ok%intercept_multi, &
                    atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR result disagree with transposed X and Y"

            ! --- Transpose COEFS ---

            status = NF_STATUS_UNDEFINED
            conf = lm_config (ncomp=Ncomp, trans_coefs=.true.)
            call pcr (conf, X, Y, coefs=coefs_T, intercept=intercept, status=status)

            values_ok = all_close (transpose(coefs_T), coefs_ok, &
                    atol=atol, rtol=rtol) .and. &
                all_close (intercept, intercept_ok, atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR coefficients and intercept disagree with " // &
                    "transposed coefficients"

            call lm_result_reset (res)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X, Y, coefs=coefs_intercept_T, status=status)

            values_ok = all_close (transpose(coefs_intercept_T), &
                coefs_intercept_ok, atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR extended coefficients disagree with transposed " // &
                    "coefficients"

            ! --- Transpose X, Y and COEFS ---

            status = NF_STATUS_UNDEFINED
            conf = lm_config (ncomp=Ncomp, trans_x=.true., trans_y=.true., trans_coefs=.true.)
            call pcr (conf, X_T, Y_T, coefs=coefs_T, intercept=intercept, status=status)

            values_ok = all_close (transpose(coefs_T), coefs_ok, &
                atol=atol, rtol=rtol) .and. &
                all_close (intercept, intercept_ok, atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR coefficients and intercept disagree with " // &
                    "transposed X and Y and coefficients"

            call lm_result_reset (res)
            status = NF_STATUS_UNDEFINED
            call pcr (conf, X_T, Y_T, coefs=coefs_intercept_T, status=status)

            values_ok = all_close (transpose(coefs_intercept_T), &
                coefs_intercept_ok, atol=atol, rtol=rtol)

            if (.not. (status == NF_STATUS_OK .and. values_ok)) &
                error stop "PCR extended coefficients disagree with transposed X " // &
                    "and Y and coefficients"

        end do
    end do

end subroutine

subroutine test_pcr_masked_const ()
    !*  TEST_PCR_MASKED_CONST performs unit tests for degenerate intercept-only
    !   PC regressions.

    type (status_t) :: status
    integer :: Nrhs, Nobs, Nlhs

    real (PREC), dimension(:,:), allocatable :: X, Y
    real (PREC), dimension(:,:), allocatable :: coefs_true
    real (PREC), dimension(:), allocatable :: intercept_true
    logical, dimension(:,:), allocatable :: mask

    Nobs = 50
    Nlhs = 10
    Nrhs = 0

    call set_seed (1234)

    call random_sample (nrhs=Nrhs, nobs=Nobs, nlhs=Nlhs, X=X, Y=Y, &
        coefs=coefs_true, intercept=intercept_true, add_intercept=.true., &
        var_error=0.2_PREC, status=status)

    allocate (mask(Nobs,Nlhs), source=.true.)

    ! --- Full mask ---

    call test_pcr_masked_runner (X, Y, mask)

end subroutine

subroutine test_pcr_masked ()

    integer :: Nrhs, Nobs, Nlhs
    real (PREC), dimension(:,:), allocatable :: X, Y
    real (PREC), dimension(:,:), allocatable :: coefs_true
    real (PREC), dimension(:), allocatable :: intercept_true, weights, runif
    logical, dimension(:,:), allocatable :: mask
    integer :: i, ifrom, ito
    type (status_t) :: status, status_mask

    Nobs = 50
    Nlhs = 10
    Nrhs = 20

    call set_seed (1234)

    call random_sample (nrhs=Nrhs, nobs=Nobs, nlhs=Nlhs, X=X, Y=Y, &
        coefs=coefs_true, intercept=intercept_true, add_intercept=.true., &
        var_error=0.1_PREC, status=status)

    allocate (mask(Nobs,Nlhs), source=.true.)
    allocate (weights(Nobs))
    allocate (runif(Nobs))

    call random_number (weights)
    weights(:) = weights / sum(weights)

    ! --- Full mask ---

    ! Unweighted regressions
    call test_pcr_masked_runner (X, Y, mask)
    ! Weighted regression
    call test_pcr_masked_runner (X, Y, mask, weights)

    ! --- Unbalanced obs. count ---

    do i = 1, Nlhs, Nlhs / 3
        call random_number (runif)
        where (runif >= real(i,PREC)/Nlhs)
            mask(:,i) = .true.
        else where
            mask(:,i) = .false.
        end where
    end do

    ! Unweighted regressions
    call test_pcr_masked_runner (X, Y, mask)
    ! Weighted regression
    call test_pcr_masked_runner (X, Y, mask, weights)

    ! --- Zero obs. for some LHS ---

    do i = 1, Nlhs, 2
        mask(:,i) = .false.
    end do

    ! Unweighted regressions
    call test_pcr_masked_runner (X, Y, mask)
    ! Weighted regression
    call test_pcr_masked_runner (X, Y, mask, weights)

end subroutine

subroutine test_pcr_masked_runner (X, Y, mask, weights)
    real (PREC), intent(in), dimension(:,:), contiguous :: X
        !* Regressor observations.
    real (PREC), intent(in), dimension(:,:), contiguous :: Y
        !* Response observations.
    logical, intent(in), dimension(:,:), contiguous :: mask
        !* Mask selecting valid observations for each response.
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
        !* Optional observation weights.

    integer :: Nrhs, Nobs, Nlhs
    real (PREC), dimension(:,:), allocatable :: X_T, Y_T
    real (PREC), dimension(:,:), allocatable :: coefs, coefs_mask, coefs_mask_T
    real (PREC), dimension(:), allocatable :: intercept, intercept_mask
    logical, dimension(:,:), allocatable :: mask_T
    integer :: ncomp_mask, ncomp, n, i, step

    logical :: values_ok
    type (status_t) :: status, status_mask
    type (lm_config) :: conf, conf_mask

    Nobs = size(X, 1)
    Nlhs = size(Y, 2)
    Nrhs = size(X, 2)

    allocate (coefs(Nrhs,Nlhs), coefs_mask(Nrhs,Nlhs), coefs_mask_T(Nlhs,Nrhs))
    allocate (intercept(Nlhs), intercept_mask(Nlhs))
    allocate (X_T(Nrhs,Nobs),Y_T(Nlhs,Nobs))
    allocate (mask_T(Nlhs,Nobs))

    X_T(:,:) = transpose(X)
    Y_T(:,:) = transpose(Y)
    mask_T(:,:) = transpose(mask)

    ! --- Fixed number of components ---

    ! Correctly handle step-size for Nrhs = 0
    step = max(1, Nrhs / 3)

    do n = 0, Nrhs, step
        conf%ncomp = n

        call pcr (conf, X, Y, coefs, ncomp, intercept, weights=weights, &
            status=status)

        status = NF_STATUS_UNDEFINED
        coefs_mask(:,:) = huge(0.0_PREC)
        call pcr (conf, X, Y, mask, coefs_mask, ncomp_mask, intercept_mask, &
            weights=weights, status=status_mask)

        values_ok = all_close_mask (coefs_mask, coefs, mask)

        if (.not. (status_mask == NF_STATUS_OK .and. values_ok)) &
            error stop "Masked PCR disagrees with unmasked PCR with fixed " // &
                "component count"

        ! Compare with quick-and-dirty masked implementation
        coefs(:,:) = huge(0.0_PREC)
        intercept(:) = huge(0.0_PREC)
        call pcr_masked (conf, X, Y, mask, weights, coefs, intercept)

        values_ok = all_close (coefs_mask, coefs, atol=1.0e-10_PREC) .and. &
            all_close (intercept_mask, intercept, atol=1.0e-10_PREC)

        if (.not. (status_mask == NF_STATUS_OK .and. values_ok)) &
            error stop "Masked PCR disagrees with the reference coefficients"
    end do

    ! --- Min. variance ---

    conf%ncomp = -1
    do i = 1, 3
        conf%var_rhs_min = 0.333d0 * i

        call pcr (conf, X, Y, coefs, ncomp, intercept, weights=weights, &
            status=status)
        status = NF_STATUS_UNDEFINED
        coefs_mask(:,:) = huge(0.0_PREC)
        call pcr (conf, X, Y, mask, coefs_mask, ncomp_mask, intercept_mask, &
            weights=weights, status=status_mask)

        values_ok = all_close_mask (coefs_mask, coefs, mask) .and. &
            (ncomp == ncomp_mask)

        if (.not. (status_mask == NF_STATUS_OK .and. values_ok)) &
            error stop "Masked PCR disagrees with unmasked PCR at minimum variance"

        ! Compare with quick-and-dirty masked implementation
        coefs(:,:) = huge(0.0_PREC)
        intercept(:) = huge(0.0_PREC)
        call pcr_masked (conf, X, Y, mask, weights, coefs, intercept)

        values_ok = all_close (coefs_mask, coefs, atol=1.0e-10_PREC) .and. &
            all_close (intercept_mask, intercept, atol=1.0e-10_PREC)

        if (.not. (status_mask == NF_STATUS_OK .and. values_ok)) &
            error stop "Masked PCR disagrees with the reference coefficients"
    end do

    ! --- Transposed X ---

    conf%ncomp = Nrhs / 2
    conf_mask = conf
    conf_mask%trans_x = .true.

    call pcr (conf, X, Y, coefs, ncomp, intercept, weights=weights, &
        status=status)
    status = NF_STATUS_UNDEFINED
    coefs_mask(:,:) = huge(0.0_PREC)
    call pcr (conf_mask, X_T, Y, mask, coefs_mask, ncomp_mask, intercept_mask, &
        weights=weights, status=status_mask)

    values_ok = all_close_mask (coefs_mask, coefs, mask) .and. &
        (ncomp == ncomp_mask)

    if (.not. (status_mask == NF_STATUS_OK .and. values_ok)) &
        error stop "Masked PCR with transposed X disagrees with unmasked PCR"

    ! Compare with quick-and-dirty masked implementation
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    call pcr_masked (conf, X, Y, mask, weights, coefs, intercept)

    values_ok = all_close (coefs_mask, coefs, atol=1.0e-10_PREC) .and. &
        all_close (intercept_mask, intercept, atol=1.0e-10_PREC)

    if (.not. (status_mask == NF_STATUS_OK .and. values_ok)) &
        error stop "Masked PCR with transposed X disagrees with reference " // &
            "coefficients"

    ! --- Transposed Y ---

    conf%ncomp = Nrhs / 2
    conf_mask = conf
    conf_mask%trans_y = .true.

    call pcr (conf, X, Y, coefs, ncomp, intercept, weights=weights, &
        status=status)
    status = NF_STATUS_UNDEFINED
    coefs_mask(:,:) = huge(0.0_PREC)
    call pcr (conf_mask, X, Y_T, mask_T, coefs_mask, ncomp_mask, intercept_mask, &
        weights=weights, status=status_mask)

    values_ok = all_close_mask (coefs_mask, coefs, mask) .and. &
        (ncomp == ncomp_mask)

    if (.not. (status_mask == NF_STATUS_OK .and. values_ok)) &
        error stop "Masked PCR with transposed Y disagrees with unmasked PCR"

    ! Compare with quick-and-dirty masked implementation
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    call pcr_masked (conf, X, Y, mask, weights, coefs, intercept)

    values_ok = all_close (coefs_mask, coefs, atol=1.0e-10_PREC) .and. &
        all_close (intercept_mask, intercept, atol=1.0e-10_PREC)

    if (.not. (status_mask == NF_STATUS_OK .and. values_ok)) &
        error stop "Masked PCR with transposed Y disagrees with reference " // &
            "coefficients"

    ! --- Transposed X, Y ---

    conf%ncomp = Nrhs / 2
    conf_mask = conf
    conf_mask%trans_x = .true.
    conf_mask%trans_y = .true.

    call pcr (conf, X, Y, coefs, ncomp, intercept, weights=weights, &
        status=status)
    status = NF_STATUS_UNDEFINED
    coefs_mask(:,:) = huge(0.0_PREC)
    call pcr (conf_mask, X_T, Y_T, mask_T, coefs_mask, ncomp_mask, intercept_mask, &
        weights=weights, status=status_mask)

    values_ok = all_close_mask (coefs_mask, coefs, mask) .and. &
        (ncomp == ncomp_mask)

    if (.not. (status_mask == NF_STATUS_OK .and. values_ok)) &
        error stop "Masked PCR with transposed X and Y disagrees with unmasked PCR"

    ! Compare with quick-and-dirty masked implementation
    coefs(:,:) = huge(0.0_PREC)
    intercept(:) = huge(0.0_PREC)
    call pcr_masked (conf, X, Y, mask, weights, coefs, intercept)

    values_ok = all_close (coefs_mask, coefs, atol=1.0e-10_PREC) .and. &
        all_close (intercept_mask, intercept, atol=1.0e-10_PREC)

    if (.not. (status_mask == NF_STATUS_OK .and. values_ok)) &
        error stop "Masked PCR with transposed X and Y disagrees with reference " // &
            "coefficients"

    ! --- Transposed COEFS ---

    conf%ncomp = Nrhs / 2
    conf_mask = conf
    conf_mask%trans_coefs = .true.

    call pcr (conf, X, Y, coefs, ncomp, intercept, weights=weights, &
        status=status)
    status = NF_STATUS_UNDEFINED
    coefs_mask_T(:,:) = huge(0.0_PREC)
    call pcr (conf_mask, X, Y, mask, coefs_mask_T, ncomp_mask, intercept_mask, &
        weights=weights, status=status_mask)

    coefs_mask(:,:) = transpose(coefs_mask_T)

    values_ok = all_close_mask (coefs_mask, coefs, mask) .and. &
        (ncomp == ncomp_mask)

    if (.not. (status_mask == NF_STATUS_OK .and. values_ok)) &
        error stop "Masked PCR with transposed coefficients disagrees with " // &
            "unmasked PCR"

end subroutine

subroutine pcr_masked (conf, X, Y, mask, weights, coefs, intercept)
    !*  PCR_MASKED is a quick-and-dirty implementation for masked PCR
    !   which performs PCR separately for each LHS, restricting observations
    !   to those not excluded for the corresponding LHS variable.
    type (lm_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:,:), contiguous :: Y
    logical, intent(in), dimension(:,:), contiguous :: mask
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
    real (PREC), intent(out), dimension(:,:), contiguous :: coefs
    real (PREC), intent(out), dimension(:), contiguous :: intercept

    real (PREC), dimension(:,:), allocatable :: Xpack
    real (PREC), dimension(:), allocatable :: Ypack
    real (PREC), dimension(:), allocatable, target :: Wpack
    real (PREC), dimension(:), pointer, contiguous :: ptr_w
    integer :: j, Nrhs, Nlhs, Nobs, k, ncomp
    type (status_t) :: status

    Nobs = size(X, 1)
    Nrhs = size(X, 2)
    Nlhs = size(Y, 2)

    allocate (Xpack(Nobs,Nrhs), Ypack(Nobs))
    allocate (Wpack(Nobs))

    nullify (ptr_w)

    do j = 1, Nlhs
        k = count (mask(:,j))

        if (k >= Nrhs) then
            call copy (X, Xpack, mask(:,j), dim=1)
            call copy (Y(:,j), Ypack, mask(:,j))
            if (present(weights)) then
                call copy (weights, Wpack, mask(:,j))
                ptr_w => Wpack(1:k)
            end if

            call pcr (conf, Xpack(1:k,:), Ypack(1:k), coefs(:,j), ncomp, &
                intercept(j), weights=ptr_w, status=status)
        else
            coefs(:,j) = 0.0
            intercept(j) = 0.0
        end if
    end do

end subroutine

function all_close_mask (actual, desired, mask, trans) result(res)
    real (PREC), intent(in), dimension(:,:), contiguous :: actual, desired
    logical, intent(in), dimension(:,:), contiguous :: mask
    logical, intent(in), optional :: trans
    logical :: res

    real (PREC), parameter :: Atol=1.0e-10_PREC, Rtol=1.0e-8_PREC
    logical :: ltrans
    integer :: i

    ltrans = .false.
    if (present(trans)) ltrans = trans

    res = .true.

    if (ltrans) then
        do i = 1, size(mask, 1)
            if (all(mask(i,:))) then
                res = res .and. all_close (actual(:,i), desired(:,i), &
                    atol=atol, rtol=rtol)
            end if
        end do
    else
        do i = 1, size(mask, 2)
            if (all(mask(:,i))) then
                res = res .and. all_close (actual(:,i), desired(:,i), &
                    atol=atol, rtol=rtol)
            end if
        end do
    end if
end function

subroutine test_pcr_cv_masked ()

    integer :: Nrhs, Nobs, Nlhs
    real (PREC), dimension(:,:), allocatable :: X, Y
    real (PREC), dimension(:,:), allocatable :: coefs_true
    real (PREC), dimension(:), allocatable :: intercept_true, weights, runif
    logical, dimension(:,:), allocatable :: mask
    integer :: ncomp, ncomp_mask
    real (PREC) :: rmse, rmse_mask
    real (PREC), dimension(:), allocatable :: rmse_lhs, rmse_lhs_mask
    integer, dimension(:), allocatable :: ncomp_lhs, ncomp_lhs_mask
    logical :: values_ok
    type (status_t) :: status, status_mask
    type (lm_config) :: conf

    Nobs = 50
    Nlhs = 10
    Nrhs = 20

    call set_seed (1234)

    call random_sample (nrhs=Nrhs, nobs=Nobs, nlhs=Nlhs, X=X, Y=Y, &
        coefs=coefs_true, intercept=intercept_true, add_intercept=.true., &
        var_error=0.1_PREC, status=status)

    allocate (mask(Nobs,Nlhs), source=.true.)
    allocate (weights(Nobs))
    allocate (runif(Nobs))

    call random_number (weights)
    weights(:) = weights / sum(weights)

    ! --- Full mask, NFOLDS > NLHS ---

    ! Set number of folds above number of LHS
    conf%cv_n = Nlhs + 1

    ! Variant w/o mask, unweighted
    ncomp = -1
    status = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp, rmse=rmse, status=status)

    ! Variant with mask
    ncomp_mask = -1
    status_mask = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_mask, mask=mask, rmse=rmse_mask, status=status_mask)

    values_ok = (ncomp == ncomp_mask) .and. abs(rmse - rmse_mask) < 1.0e-8_PREC

    if (.not. (values_ok .and. status_mask == NF_STATUS_OK)) &
        error stop "Full mask, unweighted, NFOLDS > NLHS"

    ! Variant w/o mask, weighted
    ncomp = -1
    status = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp, weights=weights, rmse=rmse, status=status)

    ! Variant with mask
    ncomp_mask = -1
    status_mask = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_mask, weights=weights, mask=mask, &
        rmse=rmse_mask, status=status_mask)

    values_ok = (ncomp == ncomp_mask) .and. abs(rmse - rmse_mask) < 1.0e-8_PREC

    if (.not. (values_ok .and. status_mask == NF_STATUS_OK)) &
        error stop "Full mask, weighted, NFOLDS > NLHS"

    ! --- Full mask, NLHS > NFOLDS ---

    ! Set number of folds above number of LHS
    conf%cv_n = Nlhs / 2

    ! Variant w/o mask, unweighted
    ncomp = -1
    status = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp, rmse=rmse, status=status)

    ! Variant with mask
    ncomp_mask = -1
    status_mask = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_mask, mask=mask, rmse=rmse_mask, status=status_mask)

    values_ok = (ncomp == ncomp_mask) .and. abs(rmse - rmse_mask) < 1.0e-8_PREC

    if (.not. (values_ok .and. status_mask == NF_STATUS_OK)) &
        error stop "Full mask, unweighted, NLHS > NFOLDS"

    ! Variant w/o mask, weighted
    ncomp = -1
    status = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp, weights=weights, rmse=rmse, status=status)

    ! Variant with mask
    ncomp_mask = -1
    status_mask = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_mask, weights=weights, mask=mask, &
        rmse=rmse_mask, status=status_mask)

    values_ok = (ncomp == ncomp_mask) .and. abs(rmse - rmse_mask) < 1.0e-8_PREC

    if (.not. (values_ok .and. status_mask == NF_STATUS_OK)) &
        error stop "Full mask, weighted, NLHS > NFOLDS"

    ! === CVs with separate result for each LHS ===

    allocate (rmse_lhs(Nlhs), rmse_lhs_mask(Nlhs))
    allocate (ncomp_lhs(Nlhs), ncomp_lhs_mask(Nlhs))

    ! --- Full mask, NFOLDS > NLHS, separate by LHS ---

    ! Set number of folds above number of LHS
    conf%cv_n = Nlhs + 1

    ! Variant w/o mask, unweighted
    ncomp_lhs(:) = -1
    rmse_lhs(:) = huge(0.0_PREC)
    status = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_lhs, rmse=rmse_lhs, status=status)

    ! Variant with mask
    ncomp_lhs_mask(:) = -1
    rmse_lhs_mask(:) = huge(0.0_PREC)
    status_mask = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_lhs_mask, mask=mask, rmse=rmse_lhs_mask, status=status_mask)

    values_ok = all(ncomp_lhs == ncomp_lhs_mask) .and. &
        all_close (rmse_lhs_mask, rmse_lhs, atol=1.0e-10_PREC)

    if (.not. (values_ok .and. status_mask == NF_STATUS_OK)) &
        error stop "By LHS: Full mask, unweighted, NFOLDS > NLHS"

    ! Variant w/o mask, weighted
    ncomp_lhs(:) = -1
    rmse_lhs(:) = huge(0.0_PREC)
    status = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_lhs, weights=weights, rmse=rmse_lhs, status=status)

    ! Variant with mask
    ncomp_lhs_mask(:) = -1
    rmse_lhs_mask(:) = huge(0.0_PREC)
    status_mask = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_lhs_mask, weights=weights, mask=mask, &
        rmse=rmse_lhs_mask, status=status_mask)

    values_ok = all(ncomp_lhs == ncomp_lhs_mask) .and. &
        all_close (rmse_lhs_mask, rmse_lhs, atol=1.0e-10_PREC)

    if (.not. (values_ok .and. status_mask == NF_STATUS_OK)) &
        error stop "By LHS: Full mask, weighted, NFOLDS > NLHS"

    ! --- Full mask, NLHS > NFOLDS ---

    ! Set number of folds above number of LHS
    conf%cv_n = Nlhs / 2

    ! Variant w/o mask, unweighted
    ncomp_lhs(:) = -1
    rmse_lhs(:) = huge(0.0_PREC)
    status = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_lhs, rmse=rmse_lhs, status=status)

    ! Variant with mask
    ncomp_lhs_mask(:) = -1
    rmse_lhs_mask(:) = huge(0.0_PREC)
    status_mask = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_lhs_mask, mask=mask, rmse=rmse_lhs_mask, status=status_mask)

    values_ok = all(ncomp_lhs == ncomp_lhs_mask) .and. &
        all_close (rmse_lhs_mask, rmse_lhs, atol=1.0e-10_PREC)

    if (.not. (values_ok .and. status_mask == NF_STATUS_OK)) &
        error stop "By LHS: Full mask, unweighted, NLHS > NFOLDS"

    ! Variant w/o mask, weighted
    ncomp_lhs(:) = -1
    rmse_lhs(:) = huge(0.0_PREC)
    status = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_lhs, weights=weights, rmse=rmse_lhs, status=status)

    ! Variant with mask
    ncomp_lhs_mask(:) = -1
    rmse_lhs_mask(:) = huge(0.0_PREC)
    status_mask = NF_STATUS_UNDEFINED
    call pcr_cv (conf, X, Y, ncomp_lhs_mask, weights=weights, mask=mask, &
        rmse=rmse_lhs_mask, status=status_mask)

    values_ok = all(ncomp_lhs == ncomp_lhs_mask) .and. &
        all_close (rmse_lhs_mask, rmse_lhs, atol=1.0e-10_PREC)

    if (.not. (values_ok .and. status_mask == NF_STATUS_OK)) &
        error stop "By LHS, Full mask, weighted, NLHS > NFOLDS"
end subroutine


subroutine test_transform_weights ()
    !* Check weighted moments and invalid inputs for both helper overloads.
    real (PREC) :: x(4,2), xout(4,2), xinplace(4,2), weights(4)
    real (PREC) :: means(2), scales(2), expected(4,2), bad_weights(4,5)
    type (status_t) :: status
    integer :: i

    ! Nonuniform weights change the moments; a zero weight is valid.
    ! The second column is affine in the first, keeping both references exact.
    x(:,1) = [0.0_PREC, 1.0_PREC, 2.0_PREC, 6.0_PREC]
    x(:,2) = 3.0_PREC * x(:,1) + 1.0_PREC
    weights = [1.0_PREC, 0.0_PREC, 2.0_PREC, 1.0_PREC]
    ! Weighted population moments divide by sum(weights), without a DOF correction.
    means = [2.5_PREC, 8.5_PREC]
    scales = sqrt(4.75_PREC) * [1.0_PREC, 3.0_PREC]
    do i = 1, 2
        expected(:,i) = (x(:,i) - means(i)) / scales(i)
    end do

    ! Check the transformed data and the moments returned to callers.
    call transform_regressors (x, xout, center=.true., scale=.true., &
        weights=weights, mean_x=means, std_x=scales, status=status)
    if (status /= NF_STATUS_OK) error stop "Weighted transform status"
    if (.not. all_close (xout, expected, atol=1.0e-12_PREC)) &
        error stop "Weighted transform moments"
    if (.not. all_close (means, [2.5_PREC, 8.5_PREC], atol=1.0e-12_PREC)) &
        error stop "Weighted means"
    if (.not. all_close (scales, sqrt(4.75_PREC) * [1.0_PREC, 3.0_PREC], &
        atol=1.0e-12_PREC)) error stop "Weighted standard deviations"

    ! Multiplying every weight by the same constant must not change the transform.
    xinplace = x
    call transform_regressors (xinplace, center=.true., scale=.true., &
        weights=7.0_PREC * weights, status=status)
    if (status /= NF_STATUS_OK) error stop "In-place weighted transform status"
    if (.not. all_close (xinplace, expected, atol=1.0e-12_PREC)) &
        error stop "In-place weighted transform"
    ! Changing input orientation must preserve observation weights and moments.
    call transform_regressors (transpose(x), xout, trans=.true., &
        center=.true., scale=.true., weights=weights, status=status)
    if (status /= NF_STATUS_OK) error stop "Transposed weighted transform status"
    if (.not. all_close (xout, expected, atol=1.0e-12_PREC)) &
        error stop "Transposed weighted transform"

    ! Each column isolates an invalid case. The negative case has a positive sum.
    bad_weights = 1.0_PREC
    bad_weights(2,1) = -1.0_PREC
    bad_weights(2,2) = ieee_value (0.0_PREC, ieee_quiet_nan)
    bad_weights(2,3) = ieee_value (0.0_PREC, ieee_positive_inf)
    bad_weights(2,4) = ieee_value (0.0_PREC, ieee_negative_inf)
    bad_weights(:,5) = 0.0_PREC
    do i = 1, size(bad_weights, 2)
        call transform_regressors (x, xout, center=.true., scale=.true., &
            weights=bad_weights(:,i), status=status)
        if (status /= NF_STATUS_INVALID_ARG) error stop "Invalid transform weights"
        xinplace = x
        call transform_regressors (xinplace, center=.true., scale=.true., &
            weights=bad_weights(:,i), status=status)
        if (status /= NF_STATUS_INVALID_ARG) error stop "Invalid in-place weights"
        ! In-place validation must finish before changing the regressor matrix.
        if (any(xinplace /= x)) error stop "Invalid weights changed input data"
    end do

    ! Reject weight vectors whose length differs from the observation count.
    call transform_regressors (x, xout, weights=weights(1:3), status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Transform weight dimensions"
    xinplace = x
    call transform_regressors (xinplace, weights=weights(1:3), status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "In-place weight dimensions"
end subroutine


subroutine test_pcr_weighted_reference ()
    !* Compare ordinary and partially masked PCR with analytical weighted LS.
    real (PREC) :: x(6,1), y(6,2), xzero(6,0), weights(6), noise(6)
    real (PREC) :: coefs(1,2), intercept(2), b(1), a, mean_x, mean_y
    real (PREC) :: ref_b, ref_a, w(6), rmse
    logical :: mask(6,2)
    type (lm_config) :: conf
    type (lm_result) :: res
    type (status_t) :: status
    integer :: j, imask, ncomp

    x(:,1) = [1, 2, 3, 4, 5, 6]
    weights = [1, 2, 3, 4, 5, 6]
    ! Fixed residuals make the weighted fit differ from the unweighted fit.
    noise = [0.5_PREC, -0.2_PREC, 0.1_PREC, -0.7_PREC, 0.3_PREC, 0.4_PREC]
    y(:,1) = 3.0_PREC + 2.0_PREC * x(:,1) + noise
    y(:,2) = -1.0_PREC + 4.0_PREC * x(:,1) - 2.0_PREC * noise
    ! With one regressor, keeping its sole component must reproduce WLS.
    conf = lm_config(ncomp=1)

    do imask = 0, 1
        mask = .true.
        ! Drop the largest-weight observation for one response; retain the other.
        if (imask == 1) mask(6,1) = .false.
        if (imask == 0) then
            call pcr (conf, x, y, coefs=coefs, intercept=intercept, &
                weights=weights, res=res, status=status)
        else
            call pcr_by_lhs (conf, x, y, mask=mask, weights=weights, &
                coefs=coefs, intercept=intercept, res=res, status=status)
        end if
        if (status /= NF_STATUS_OK) error stop "Reference PCR status"
        do j = 1, 2
            ! The closed-form WLS reference is independent of numfort's OLS.
            ! Zeroing excluded weights applies the response mask to every moment.
            w = merge(weights, 0.0_PREC, mask(:,j))
            mean_x = sum(w * x(:,1)) / sum(w)
            mean_y = sum(w * y(:,j)) / sum(w)
            ref_b = sum(w * (x(:,1) - mean_x) * (y(:,j) - mean_y)) &
                / sum(w * (x(:,1) - mean_x)**2)
            ref_a = mean_y - ref_b * mean_x
            if (abs(coefs(1,j) - ref_b) > 1.0e-12_PREC &
                .or. abs(intercept(j) - ref_a) > 1.0e-12_PREC) &
                error stop "PCR disagrees with analytical WLS"
            ! Result objects must agree with the separate coefficient outputs.
            if (abs(res%coefs_multi(1,j) - ref_b) > 1.0e-12_PREC &
                .or. abs(res%intercept_multi(j) - ref_a) > 1.0e-12_PREC) &
                error stop "PCR result disagrees with analytical WLS"
            if (imask == 0) then
                ! The 1D wrapper must use the same weighted centering as 2D PCR.
                call pcr (conf, x, y(:,j), coefs=b, intercept=a, &
                    weights=weights, status=status)
                if (status /= NF_STATUS_OK) error stop "1D reference PCR status"
                if (abs(b(1) - ref_b) > 1.0e-12_PREC &
                    .or. abs(a - ref_a) > 1.0e-12_PREC) &
                    error stop "1D PCR disagrees with analytical WLS"
            end if
        end do
    end do

    ! With no regressors, the fitted intercept is the weighted response mean.
    ! This checks the branch that skips PCA entirely.
    conf%ncomp = 0
    call pcr (conf, xzero, y, weights=weights, res=res, status=status)
    if (status /= NF_STATUS_OK) error stop "Constant-only reference PCR status"
    if (.not. all_close (res%intercept_multi, matmul(weights, y) / sum(weights), &
        atol=1.0e-12_PREC)) error stop "Constant-only weighted means"
    call pcr (conf, xzero, y(:,1), weights=weights, res=res, status=status)
    if (status /= NF_STATUS_OK) error stop "1D constant-only reference PCR status"
    if (abs(res%intercept - sum(weights * y(:,1)) / sum(weights)) > 1.0e-12_PREC) &
        error stop "1D constant-only weighted mean"

    ! Exact linear data must also predict perfectly across masked CV folds.
    ! Retain the partial mask above to exercise response-specific training samples.
    y(:,1) = 3.0_PREC + 2.0_PREC * x(:,1)
    y(:,2) = y(:,1)
    conf = lm_config(cv_n=3)
    call pcr_cv (conf, x, y, ncomp, weights=weights, mask=mask, &
        rmse=rmse, status=status)
    if (status /= NF_STATUS_OK) error stop "Partially masked reference CV status"
    if (ncomp /= 1 .or. rmse > 1.0e-12_PREC) &
        error stop "Partially masked CV of exact linear data"
end subroutine


subroutine test_pcr_preprocessing ()
    !* Check preprocessing options against separate fits and explicit CV folds.
    integer, parameter :: Nobs = 12, Nrhs = 2, Nlhs = 4
    real (PREC) :: x(Nobs,Nrhs), y(Nobs,Nlhs), xp(Nobs,Nrhs), yp(Nobs)
    real (PREC) :: xt(Nrhs,Nobs), yt(Nlhs,Nobs)
    real (PREC), target :: weights(Nobs), wp(Nobs)
    real (PREC), pointer, contiguous :: ptr_weights(:), ptr_wp(:)
    logical, target :: mask(Nobs,Nlhs)
    logical :: mask_t(Nlhs,Nobs)
    logical, pointer, contiguous :: ptr_mask(:,:)
    logical :: train(Nobs)
    real (PREC) :: coefs(Nrhs,Nlhs), intercept(Nlhs), b(Nrhs), a
    real (PREC) :: coefs_t(Nlhs,Nrhs), intercept_t(Nlhs), residual(Nobs)
    real (PREC) :: ssr(0:Nrhs,Nlhs), total_weight(Nlhs), ref_rmse(Nlhs)
    real (PREC) :: rmse_lhs(Nlhs), rmse, pooled_rmse
    integer :: ncomp, ncomp_lhs(Nlhs), ref_ncomp(Nlhs), pooled_ncomp
    integer :: i, j, k, h, n, ifrom, ito, icenter, iscale, iw, imask, nfolds
    integer :: requested(Nlhs)
    type (lm_config) :: conf, ref_conf, trans_conf
    type (status_t) :: status

    ! Different predictor scales make scaling relevant when only one PC is kept.
    do i = 1, Nobs
        x(i,1) = real(i, PREC)
        x(i,2) = 10.0_PREC * real(mod(i * i, 7), PREC)
    end do
    ! Deterministic nonlinear residuals give CV nonzero prediction errors.
    do j = 1, Nlhs
        y(:,j) = real(j, PREC) * x(:,1) - 0.2_PREC * x(:,2) &
            + sin(x(:,1) * real(j, PREC)) + 3.0_PREC
    end do
    ! GCC PR125998 affects TRANSPOSE expressions passed to CONTIGUOUS,TARGET
    ! dummies. Our GNU 13-16 CI builds supplied the original shape, e.g. (12,2)
    ! instead of (2,12), with incorrect strides; PCR rejected the arguments.
    ! A standalone descriptor probe reproduced this without calling numfort.
    ! Materialize X, Y, and the mask to retain transposed-input test coverage.
    ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=125998
    xt = transpose(x)
    yt = transpose(y)
    ! Exercise zero, truncated, and full component counts across responses.
    requested = [0, 1, 2, 1]

    do iw = 0, 1
        weights = 1.0_PREC
        ! Null pointers omit optional weights; ones still weight the reference SSE.
        nullify (ptr_weights, ptr_wp)
        if (iw == 1) then
            weights = x(:,1)
            weights(2) = 0.0_PREC
            ptr_weights => weights
        end if
        do icenter = 0, 1
            do iscale = 0, 1
                ! The uncentered case fits through the origin; vary scaling separately.
                conf = lm_config(center=icenter == 1, scale=iscale == 1, &
                    add_intercept=icenter == 1)
                do imask = 0, 2
                    ! Compare omitted, all-true, and response-specific masks.
                    mask = .true.
                    nullify (ptr_mask)
                    if (imask > 0) ptr_mask => mask
                    if (imask == 2) then
                        ! Alternate partial and full samples to exercise both paths.
                        mask([2, 9],1) = .false.
                        mask([4, 11],3) = .false.
                    end if
                    call pcr_by_lhs (conf, x, y, mask=ptr_mask, &
                        weights=ptr_weights, ncomp=requested, coefs=coefs, &
                        intercept=intercept, status=status)
                    if (status /= NF_STATUS_OK) error stop "By-LHS options status"
                    ! Explicitly pack each sample and fit it with ordinary PCR.
                    do j = 1, Nlhs
                        n = count(mask(:,j))
                        do i = 1, Nrhs
                            xp(1:n,i) = pack(x(:,i), mask(:,j))
                        end do
                        yp(1:n) = pack(y(:,j), mask(:,j))
                        wp(1:n) = pack(weights, mask(:,j))
                        if (iw == 1) ptr_wp => wp(1:n)
                        ref_conf = conf
                        ref_conf%ncomp = requested(j)
                        call pcr (ref_conf, xp(1:n,:), yp(1:n), coefs=b, &
                            intercept=a, weights=ptr_wp, status=status)
                        if (status /= NF_STATUS_OK) error stop "Packed options status"
                        if (.not. all_close (coefs(:,j), b, atol=1.0e-11_PREC) &
                            .or. abs(intercept(j) - a) > 1.0e-11_PREC) &
                            error stop "By-LHS options disagree with packed PCR"
                    end do

                    ! Transpose both inputs and the coefficient output for the same fit.
                    trans_conf = conf
                    trans_conf%trans_x = .true.
                    trans_conf%trans_y = .true.
                    trans_conf%trans_coefs = .true.
                    ! Use a named mask array for the GCC workaround described above.
                    mask_t = transpose(mask)
                    call pcr_by_lhs (trans_conf, xt, yt, &
                        mask=mask_t, weights=ptr_weights, ncomp=requested, &
                        coefs=coefs_t, intercept=intercept_t, status=status)
                    if (status /= NF_STATUS_OK) error stop "Transposed by-LHS status"
                    if (.not. all_close (transpose(coefs_t), coefs, &
                        atol=1.0e-11_PREC) &
                        .or. .not. all_close (intercept_t, intercept, &
                        atol=1.0e-11_PREC)) error stop "Transposed by-LHS options"

                    ! Four responses exceed three folds, but not six folds:
                    ! exercise parallelization over responses and over folds.
                    do nfolds = 3, 6, 3
                        conf%cv_n = nfolds
                        ssr = 0.0_PREC
                        do j = 1, Nlhs
                            total_weight(j) = sum(weights, mask=mask(:,j))
                            do h = 1, nfolds
                                ! Match the default consecutive folds.
                                ifrom = (h - 1) * (Nobs / nfolds) + 1
                                ito = h * (Nobs / nfolds)
                                ! Exclude held-out and missing data from training.
                                train = mask(:,j)
                                train(ifrom:ito) = .false.
                                n = count(train)
                                do i = 1, Nrhs
                                    xp(1:n,i) = pack(x(:,i), train)
                                end do
                                yp(1:n) = pack(y(:,j), train)
                                wp(1:n) = pack(weights, train)
                                if (iw == 1) ptr_wp => wp(1:n)
                                ! Build the error curve, including zero PCs.
                                do k = 0, Nrhs
                                    ref_conf = conf
                                    ref_conf%ncomp = k
                                    call pcr (ref_conf, xp(1:n,:), yp(1:n), &
                                        coefs=b, intercept=a, weights=ptr_wp, &
                                        status=status)
                                    if (status /= NF_STATUS_OK) &
                                        error stop "Explicit CV fit status"
                                    ! Score only unmasked held-out observations.
                                    residual = matmul(x, b) + a - y(:,j)
                                    ssr(k,j) = ssr(k,j) &
                                        + sum(weights(ifrom:ito) &
                                        * residual(ifrom:ito)**2, &
                                        mask=mask(ifrom:ito,j))
                                end do
                            end do
                        end do
                        ! Convert one-based MINLOC positions to zero-based PC counts.
                        ref_ncomp = minloc(ssr, dim=1) - 1
                        ref_rmse = sqrt(minval(ssr, dim=1) / total_weight)
                        ! Pool raw weighted errors, so larger samples carry more weight.
                        pooled_ncomp = minloc(sum(ssr, dim=2), dim=1) - 1
                        pooled_rmse = sqrt(minval(sum(ssr, dim=2)) &
                            / sum(total_weight))

                        ! Check the selected PC count and its prediction error.
                        call pcr_cv (conf, x, y, ncomp_lhs, weights=ptr_weights, &
                            mask=ptr_mask, rmse=rmse_lhs, status=status)
                        if (status /= NF_STATUS_OK) error stop "Separate CV status"
                        if (any(ncomp_lhs /= ref_ncomp) &
                            .or. .not. all_close (rmse_lhs, ref_rmse, &
                            atol=1.0e-11_PREC)) error stop "Separate CV fold reference"
                        call pcr_cv (conf, x, y, ncomp, weights=ptr_weights, &
                            mask=ptr_mask, rmse=rmse, status=status)
                        if (status /= NF_STATUS_OK) error stop "Pooled CV status"
                        if (ncomp /= pooled_ncomp &
                            .or. abs(rmse - pooled_rmse) > 1.0e-11_PREC) &
                            error stop "Pooled CV fold reference"
                    end do
                end do
            end do
        end do
    end do
end subroutine


subroutine test_pcr_invalid_weights ()
    !* Reject bad weights before square roots, including zero-weight subsets.
    real (PREC) :: x(6,1), y(6,4), weights(6), bad_weights(6,4)
    real (PREC) :: coefs(1,4), intercept(4)
    logical :: mask(6,4)
    type (lm_config) :: conf
    type (lm_result) :: res
    type (status_t) :: status
    integer :: i, j, ncomp, ncomp_lhs(4)

    x(:,1) = [1, 2, 3, 4, 5, 6]
    do j = 1, 4
        y(:,j) = 3.0_PREC + 2.0_PREC * x(:,1)
    end do
    mask = .true.
    conf = lm_config(ncomp=1, cv_n=3)
    ! Valid data isolates weight validation in each public regression entry point.
    ! The negative case retains a positive total; the others test NaN, Inf, and zero.
    bad_weights = 1.0_PREC
    bad_weights(2,1) = -1.0_PREC
    bad_weights(2,2) = ieee_value (0.0_PREC, ieee_quiet_nan)
    bad_weights(2,3) = ieee_value (0.0_PREC, ieee_positive_inf)
    bad_weights(:,4) = 0.0_PREC
    do i = 1, size(bad_weights, 2)
        call pcr (conf, x, y, weights=bad_weights(:,i), res=res, status=status)
        if (status /= NF_STATUS_INVALID_ARG) error stop "Invalid PCR weights"
        call pcr_by_lhs (conf, x, y, weights=bad_weights(:,i), res=res, &
            status=status)
        if (status /= NF_STATUS_INVALID_ARG) error stop "Invalid by-LHS weights"
        call pcr_cv (conf, x, y, ncomp, weights=bad_weights(:,i), status=status)
        if (status /= NF_STATUS_INVALID_ARG) error stop "Invalid CV weights"
        call pcr_cv (conf, x, y, ncomp_lhs, mask=mask, &
            weights=bad_weights(:,i), status=status)
        if (status /= NF_STATUS_INVALID_ARG) error stop "Invalid separate CV weights"
    end do

    ! The total is positive, but the first fold has no positive training weights.
    ! Later folds are valid, so their success must not erase the earlier failure.
    weights = [1.0_PREC, 1.0_PREC, 0.0_PREC, 0.0_PREC, 0.0_PREC, 0.0_PREC]
    ! CV combines status flags from successful and unsuccessful folds.
    call pcr_cv (conf, x, y, ncomp, weights=weights, status=status)
    if (.not. (NF_STATUS_INVALID_ARG .in. status)) &
        error stop "Zero-weight CV training fold"
    call pcr_cv (conf, x, y, ncomp_lhs, weights=weights, mask=mask, status=status)
    if (.not. (NF_STATUS_INVALID_ARG .in. status)) &
        error stop "Zero-weight masked training fold"

    ! Remove all positive weights from one response while keeping the others valid.
    ! Both masked APIs must reject that subset before storing results.
    mask(1:2,1) = .false.
    call pcr_by_lhs (conf, x, y, mask=mask, weights=weights, res=res, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Zero-weight by-LHS subset"
    call pcr (conf, x, y, mask, coefs=coefs, intercept=intercept, &
        weights=weights, status=status)
    if (status /= NF_STATUS_INVALID_ARG) error stop "Zero-weight masked PCR subset"
end subroutine


end
