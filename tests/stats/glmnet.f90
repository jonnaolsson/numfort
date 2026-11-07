

program test_numfort_stats_glmnet

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common
    use numfort_common_testing
    use numfort_stats
    use numfort_stats_lm, lm_result => lm_result_real64, lm_config => lm_config_real64
    use numfort_stats_glmnet, enet_config => enet_config_real64, &
        enet_result => enet_result_real64
    use numfort_stats_data_helpers

    use blas95, only: BLAS_GEMV => GEMV

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_equiv ()
    ! Test special case that triggers different code paths
    call test_equiv_2d (nlhs=1)
    call test_equiv_2d (nlhs=3)
    ! Test with more LHS than RHS variables (routine uses 10 RHS variables)
    call test_equiv_2d (nlhs=13)

    call test_ridge_1d_args ()
    call test_ridge_2d_args ()
    call test_ridge_1d ()

end subroutine

subroutine test_equiv ()
    !*  Unit tests that should all yield same results for different
    !   sets of X and y inputs which are equivalent after applying
    !   the required transformations, e.g.
    !       - Pre-scaling variables vs. letting ENET_* scale them
    !       - Transposed X
    !       - Adding constants to X which are dropped in ENET_*
    !       - Manually passing the grid of ALPHAS which would be created in
    !           ENET_CV

    real (PREC), dimension(:,:), allocatable :: X, X_alt, X_T, Y_multi
    real (PREC), dimension(:), allocatable :: coefs, coefs_true, coefs_default
    real (PREC), dimension(:), allocatable, target :: y, y_alt
    real (PREC), dimension(:), allocatable :: coefs_ols, rwork, alphas
    real (PREC), dimension(:), allocatable :: mean_x, std_x, mean_x_alt, std_x_alt
    integer, dimension(:), allocatable :: iorder, folds_ifrom, folds_size
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Y, ptr_Y_T
    real (PREC), dimension(:,:), allocatable :: coefs_multi
    real (PREC), dimension(:), allocatable :: intercept_multi
    real (PREC) :: rmse_default, rmse
    real (PREC) :: intercept, intercept_true, l1_ratio
    real (PREC) :: intercept_ols, intercept_default
    real (PREC) :: alpha, alpha_default
    real (PREC) :: const, mean_y
    integer :: Nobs, Nvars, Nlhs, Nconst, i, j
    type (status_t) :: status
    type (lm_config) :: conf_lm
    type (enet_config) :: conf
    real (PREC), parameter :: atol=1.0e-12_PREC, rtol=1.0e-4_PREC
    logical :: values_ok
    logical, dimension(:), allocatable :: mask, mask_lhs

    call set_seed (1234)

    ! --- Create data ---

    Nvars = 10
    Nobs = 103

    allocate (X(Nobs, Nvars))
    allocate (mean_x(Nvars), std_x(Nvars))

    call random_number (X)

    call random_number (mean_x)
    call random_number (std_x)

    mean_x(:) = (mean_x - 0.50d0) * 2.0
    std_x(:) = std_x * 5.0d0

    X = X - 0.5d0
    do i = 1, Nvars
        X(:,i) = X(:,i) * std_x(i) + mean_x(i)
    end do

    ! --- True coefs ---

    allocate (coefs_true(Nvars))
    do i = 1, Nvars
        coefs_true(i) = (-1.0_PREC)**i * exp(-i/10.0_PREC)
    end do

    ! Set half to 0
    coefs_true(Nvars/2+1:) = 0.0

    ! True intercept
    intercept_true = 1.234d0

    ! --- Outcome variable ---

    allocate (y(Nobs), source=intercept_true)
    call BLAS_GEMV (X, coefs_true, y, beta=1.0_PREC)

    ! Add some normally-distributed noise
    allocate (rwork(Nobs))

    call rvs (norm, rwork, scale=0.01_PREC)

    y(:) = y + rwork

    deallocate (rwork)

    ! Pointer to degenerate multi-LHS outcome
    ptr_Y(1:Nobs,1:1) => y
    ptr_Y_T(1:1,1:Nobs) => y

    ! --- OLS ---

    allocate (coefs_ols(Nvars))
    conf_lm = lm_config (add_intercept=.false.)
    call ols (conf_lm, X, y, coefs_ols, intercept=intercept_ols, status=status)

    ! --- Default CV ---

    l1_ratio = 0.9

    conf%cv_n = 50
    conf%alpha_eps = 1.0e-5_PREC
    conf%cv_parsimonious = .false.
    status = NF_STATUS_UNDEFINED
    call enet_cv (conf, X, y, alpha_default, l1_ratio, rmse=rmse_default, &
        status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "ENET_CV: default config"

    ! Fit coefs using CV'd alpha
    allocate (coefs_default(Nvars))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X, y, alpha_default, l1_ratio, intercept_default, &
        coefs_default, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "ENET_FIT: default config"

    ! --- Default CV, multi-LHS ---

    ! Use default internal path that is the same as for single-outcome model
    status = NF_STATUS_UNDEFINED
    call enet_cv (conf, X, ptr_Y, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (rmse, rmse_default, atol=atol, rtol=rtol) .and. &
        all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV_MULTI: default config"

    ! Fit coefs using CV'd alpha
    allocate (coefs_multi(Nvars,1), source=huge(0.0_PREC))
    allocate (intercept_multi(1), source=huge(0.0_PREC))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X, ptr_Y, alpha, l1_ratio, intercept_multi, &
        coefs_multi, status=status)
    values_ok = all_close(intercept_multi(1), intercept_default, atol=atol, rtol=rtol) .and. &
        all_close (coefs_multi(:,1), coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT_MULTI: default config"

    ! Force multi-outcome code path even for single-outcome model
    conf%force_multi = .true.
    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X, ptr_Y, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (rmse, rmse_default, atol=atol, rtol=rtol) .and. &
        all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV_MULTI: default config, multi-outcome code path"

    ! Fit coefs using CV'd alpha
    status = NF_STATUS_UNDEFINED
    intercept_multi(:) = huge(0.0_PREC)
    coefs_multi(:,:) = huge(0.0_PREC)
    call enet_fit (conf, X, ptr_Y, alpha, l1_ratio, intercept_multi, &
        coefs_multi, status=status)
    values_ok = all_close(intercept_multi(1), intercept_default, atol=atol, rtol=rtol) .and. &
        all_close (coefs_multi(:,1), coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT_MULTI: default config, multi-outcome code path"

    deallocate (coefs_multi, intercept_multi)

    ! --- Transposed X ---

    conf%trans_x = .true.

    allocate (X_T(Nvars, Nobs))
    X_T(:,:) = transpose (X)

    status = NF_STATUS_UNDEFINED
    call enet_cv (conf, X_T, y, alpha, l1_ratio, rmse=rmse, status=status)

    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol) .and. &
        all_close (rmse, rmse_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: transposed X"

    ! Fit coefs using CV'd alpha
    allocate (coefs(Nvars))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_T, y, alpha, l1_ratio, intercept, coefs, status=status)

    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol) .and. &
        all_close (rmse, rmse_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed X"

    ! Multi-outcomes, transposed X and Y
    ! Force multi-outcome code path even for single-outcome model
    conf%force_multi = .true.
    conf%trans_y = .true.
    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X_T, ptr_Y_T, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (rmse, rmse_default, atol=atol, rtol=rtol) .and. &
        all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV_MULTI: transposed X and Y, multi-outcome code path"

    ! Fit coefs using CV'd alpha
    allocate (coefs_multi(Nvars,1), source=huge(0.0_PREC))
    allocate (intercept_multi(1), source=huge(0.0_PREC))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_T, ptr_Y_T, alpha, l1_ratio, intercept_multi, &
        coefs_multi, status=status)
    values_ok = all_close(intercept_multi(1), intercept_default, atol=atol, rtol=rtol) .and. &
        all_close (coefs_multi(:,1), coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT_MULTI: transposed X and Y, multi-outcome code path"

    deallocate (coefs, X_T)
    deallocate (coefs_multi, intercept_multi)
    conf%trans_x = .false.
    conf%trans_y = .false.

    ! --- Pre-specified ALPHA grid ---

    allocate (alphas(conf%alpha_n))

    allocate (folds_ifrom(conf%cv_n), folds_size(conf%cv_n))

    ! Create CV chunk indices
    call split_uniform (Nobs, conf%cv_n, folds_ifrom, folds_size, status)

    ! Create the same grid of ALPHAS used for CV
    call create_alpha_grid_cv (conf, X, y, l1_ratio, folds_ifrom, folds_size, &
        alphas, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "Create alpha grid for CV"

    ! Call with pre-specified ALPHAS
    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X, y, alpha, l1_ratio, alphas=alphas, rmse=rmse, &
        status=status)

    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol) .and. &
        all_close (rmse, rmse_default, atol=atol, rtol=rtol)

    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: Use pre-specified ALPHAS"

    ! Call with randomly-ordered alphas
    allocate (iorder(conf%alpha_n))
    call random_order (iorder)

    allocate (rwork(conf%alpha_n))
    do i = 1, conf%alpha_n
        rwork(i) = alphas(iorder(i))
    end do

    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X, y, alpha, l1_ratio, alphas=rwork, rmse=rmse, &
        status=status)

    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol) .and. &
        all_close (rmse, rmse_default, atol=atol, rtol=rtol)

    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: Use pre-specified, randomly ordered ALPHAS"

    deallocate (rwork)

    ! --- X with constants ---

    Nconst = 2
    allocate (mask(Nvars + Nconst), source=.true.)
    mask(4) = .false.
    mask(9) = .false.
    allocate (X_alt(Nobs, Nvars + Nconst))

    j = 0
    do i = 1, size(mask)
        if (mask(i)) then
            j = j + 1
            X_alt(:,i) = X(:,j)
        else
            call random_number (const)
            X_alt(:,i) = const
        end if
    end do

    status = NF_STATUS_UNDEFINED
    call enet_cv (conf, X_alt, y, alpha, l1_ratio, rmse=rmse, status=status)

    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol) .and. &
        all_close (rmse, rmse_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: add. const. in X"

    allocate (coefs(Nvars + Nconst))

    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_alt, y, alpha, l1_ratio, intercept, coefs, status=status)

    values_ok = all(pack(coefs, .not. mask) == 0.0_PREC) .and. &
        all_close (pack(coefs, mask), coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X"

    ! Multi-outcomes
    ! Force multi-outcome code path even for single-outcome model
    conf%force_multi = .true.
    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X_alt, ptr_Y, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (rmse, rmse_default, atol=atol, rtol=rtol) .and. &
        all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV_MULTI: add. const. in X, multi-outcome code path"

    ! Fit coefs using CV'd alpha
    allocate (coefs_multi(Nvars+Nconst,1), source=huge(0.0_PREC))
    allocate (intercept_multi(1), source=huge(0.0_PREC))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_alt, ptr_Y, alpha, l1_ratio, intercept_multi, &
        coefs_multi, status=status)
    values_ok = all_close(intercept_multi(1), intercept_default, atol=atol, rtol=rtol) .and. &
        all_close (pack(coefs_multi(:,1), mask), coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT_MULTI: add. const. in X, multi-outcome code path"

    deallocate (coefs_multi, intercept_multi)

    ! --- X with constants, transposed ---
    conf%trans_x = .true.

    allocate (X_T(Nvars + Nconst, Nobs))
    X_T(:,:) = transpose (X_alt)

    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X_T, y, alpha, l1_ratio, rmse=rmse, status=status)

    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol) .and. &
        all_close (rmse, rmse_default, atol=atol, rtol=rtol)

    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: transposed X, add. const."

    ! Fit coefs using CV'd alpha
    status = NF_STATUS_UNDEFINED
    coefs(:) = huge(0.0_PREC)
    intercept = huge(0.0_PREC)
    call enet_fit (conf, X_T, y, alpha, l1_ratio, intercept, coefs, status=status)

    values_ok = all(pack(coefs, .not. mask) == 0.0_PREC) .and. &
        all_close (pack(coefs, mask), coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, intercept_default, atol=atol, rtol=rtol)

    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed X, add. const."

    ! Multi-outcomes
    ! Force multi-outcome code path even for single-outcome model
    conf%force_multi = .true.
    conf%trans_y = .true.
    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X_T, ptr_Y_T, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (rmse, rmse_default, atol=atol, rtol=rtol) .and. &
        all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV_MULTI: X^T with constants, Y^T, multi-outcome code path"

    ! Fit coefs using CV'd alpha
    allocate (coefs_multi(Nvars+Nconst,1), source=huge(0.0_PREC))
    allocate (intercept_multi(1), source=huge(0.0_PREC))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_T, ptr_Y_T, alpha, l1_ratio, intercept_multi, &
        coefs_multi, status=status)
    values_ok = all_close(intercept_multi(1), intercept_default, atol=atol, rtol=rtol) .and. &
        all_close (pack(coefs_multi(:,1), mask), coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT_MULTI: X^T with constants, Y^T, multi-outcome code path"

    deallocate (coefs_multi, intercept_multi)

    conf%trans_x = .false.
    conf%trans_y = .false.

    ! --- X, y with constants and pre-centered variables ---

    ! For the remaining calls we need to pass the pre-computed ALPHAS,
    ! since otherwise the max. grid point will differ as X and y are not
    ! standardized within CV chunks!

    allocate (y_alt(Nobs), source=y)

    allocate (mean_x_alt(Nvars + Nconst), std_x_alt(Nvars + Nconst))
    call standardize (X_alt, dim=1, mean_x=mean_x_alt, std_x=std_x_alt, &
        center=.true., scale=.false., status=status)
    call standardize (y_alt, mean_x=mean_y, center=.true., scale=.false., status=status)

    ! Variables are pre-centered, no need to center in ENET routines
    conf%center = .false.

    alpha = 0.0
    rmse = 0.0
    status = NF_STATUS_UNDEFINED
    call enet_cv (conf, X_alt, y_alt, alpha, l1_ratio, rmse=rmse, &
        alphas=alphas, status=status)
    ! It's not entirely clear that these should be identical since the
    ! RMSEs across chunks will differ, and hence in principle a different
    ! point on the ALPHAS grid could be optiomal.
    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: add. const. in X, pre-centered variables"

    status = NF_STATUS_UNDEFINED
    intercept = huge(0.0_PREC)
    coefs(:) = huge(0.0_PREC)
    call enet_fit (conf, X_alt, y_alt, alpha, l1_ratio, intercept, coefs, status=status)

    values_ok = all(pack(coefs, .not. mask) == 0.0_PREC) .and. &
        all_close (pack(coefs, mask), coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, 0.0_PREC, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X, pre-centered variables"

    ! Multi-outcomes

    ptr_Y(1:Nobs,1:1) => y_alt
    ! Force multi-outcome code path even for single-outcome model
    conf%force_multi = .true.
    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X_alt, ptr_Y, alpha, l1_ratio, rmse=rmse, &
        alphas=alphas, status=status)
    values_ok = all_close (rmse, rmse_default, atol=atol, rtol=rtol) .and. &
        all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV_MULTI: X with constants, pre-centered, multi-outcome " // &
            "code path"

    ! Fit coefs using CV'd alpha
    allocate (coefs_multi(Nvars+Nconst,1), source=huge(0.0_PREC))
    allocate (intercept_multi(1), source=huge(0.0_PREC))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_alt, ptr_Y, alpha, l1_ratio, intercept_multi, &
        coefs_multi, status=status)
    values_ok = all_close(intercept_multi(1), 0.0_PREC, atol=atol, rtol=rtol) .and. &
        all_close (pack(coefs_multi(:,1), mask), coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT_MULTI: X^T with constants, Y^T, multi-outcome code path"

    deallocate (coefs_multi, intercept_multi)

    ! --- X pre-centered and pre-scaled, y pre-centered ---

    ! X_alt is already center from before
    call standardize (X_alt, dim=1, std_x=std_x_alt, center=.false., &
        scale=.true., status=status)

    ! Eliminate all previous constants since these will now show up as NaN
    do i = 1, size(X_alt, 2)
        if (.not. mask(i)) then
            X_alt(:,i) = 0.0
        end if
    end do

    ! Variables are pre-centered, no need to center in ENET routines
    conf%center = .false.
    conf%scale = .false.

    alpha = 0.0
    rmse = 0.0
    status = NF_STATUS_UNDEFINED
    call enet_cv (conf, X_alt, y_alt, alpha, l1_ratio, rmse=rmse, &
        alphas=alphas, status=status)
    ! It's not entirely clear that these should be identical since the
    ! RMSEs across chunks will differ, and hence in principle a different
    ! point on the ALPHAS grid could be optiomal.
    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: add. const. in X, pre-centered, pre-scaled variables"

    status = NF_STATUS_UNDEFINED
    intercept = huge(0.0_PREC)
    coefs(:) = huge(0.0_PREC)
    call enet_fit (conf, X_alt, y_alt, alpha, l1_ratio, intercept, coefs, status=status)

    ! Undo scaling of coefficients
    where (std_x_alt > 0.0)
        coefs = coefs / std_x_alt
    end where

    values_ok = all(pack(coefs, .not. mask) == 0.0_PREC) .and. &
        all_close (pack(coefs, mask), coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, 0.0_PREC, atol=atol, rtol=rtol)

    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X, pre-centered, pre-scaled variables"

    ! Multi-outcomes

    ! Force multi-outcome code path even for single-outcome model
    conf%force_multi = .true.
    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X_alt, ptr_Y, alpha, l1_ratio, rmse=rmse, &
        alphas=alphas, status=status)
    values_ok = all_close (rmse, rmse_default, atol=atol, rtol=rtol) .and. &
        all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV_MULTI: X with constants, pre-centered, multi-outcome " // &
            "code path"

    ! Fit coefs using CV'd alpha
    allocate (coefs_multi(Nvars+Nconst,1), source=huge(0.0_PREC))
    allocate (intercept_multi(1), source=huge(0.0_PREC))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_alt, ptr_Y, alpha, l1_ratio, intercept_multi, &
        coefs_multi, status=status)

    ! Undo scaling of coefficients
    do i = 1, size(coefs_multi, 1)
        if (std_x_alt(i) > 0.0) then
            coefs_multi(i,1) = coefs_multi(i,1) / std_x_alt(i)
        end if
    end do

    values_ok = all_close(intercept_multi(1), 0.0_PREC, atol=atol, rtol=rtol) .and. &
        all_close (pack(coefs_multi(:,1), mask), coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT_MULTI: X^T with constants, Y^T, multi-outcome code path"

    deallocate (coefs_multi, intercept_multi)

    ! --- Multiple LHS variables, some constant ---

    Nlhs = 3
    allocate (Y_multi(Nlhs,Nobs))

    Y_multi(1,:) = 0.56789_PREC
    Y_multi(2,:) = y
    Y_multi(3,:) = 0.0

    conf%trans_x = .true.
    conf%trans_y = .true.
    conf%drop_const = .true.
    conf%center = .true.
    conf%scale = .true.

    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X_T, Y_multi, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (rmse, rmse_default, atol=atol, rtol=rtol) .and. &
        all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV_MULTI: X^T with constants, Y^T with constants"

    allocate (coefs_multi(Nvars+Nconst,Nlhs), source=huge(0.0_PREC))
    allocate (intercept_multi(Nlhs), source=huge(0.0_PREC))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_T, Y_multi, alpha, l1_ratio, intercept_multi, &
        coefs_multi, status=status)
    values_ok = all_close (intercept_multi(2), intercept_default, atol=atol, rtol=rtol) &
        .and. all_close (intercept_multi(1), 0.56789_PREC, atol=atol, rtol=rtol) &
        .and. all_close (intercept_multi(3), 0.0_PREC, atol=atol, rtol=rtol) &
        .and. all_close (pack(coefs_multi(:,2), mask), coefs_default, atol=atol, rtol=rtol) &
        .and. all(coefs_multi(:,1) == 0.0_PREC) &
        .and. all(coefs_multi(:,3) == 0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT_MULTI: X^T with constants, Y^T with constants"

end subroutine

subroutine test_equiv_2d (nlhs)
    !*  Unit tests that should all yield same results for different
    !   sets of X and multiple Y inputs which are equivalent after applying
    !   the required transformations, e.g.
    !       - Pre-scaling variables vs. letting ENET_* scale them
    !       - Transposed X
    !       - Adding constants to X which are dropped in ENET_*
    !       - Manually passing the grid of ALPHAS which would be created in
    !           ENET_CV
    integer, intent(in) :: Nlhs

    real (PREC), dimension(:,:), allocatable :: X, X_alt, X_T
    real (PREC), dimension(:,:), allocatable :: coefs, coefs_T, coefs_true, coefs_default
    real (PREC), dimension(:,:), allocatable :: coefs_all, coefs_const
    real (PREC), dimension(:), allocatable :: intercept, intercept_default, intercept_true
    real (PREC), dimension(:,:), allocatable :: Y, Y_T, Y_alt
    real (PREC), dimension(:), allocatable :: rwork, alphas
    real (PREC), dimension(:), allocatable :: mean_x, std_x, mean_x_alt, scale_x_alt, mean_y
    integer, dimension(:), allocatable :: iorder, folds_ifrom, folds_size
    integer, dimension(:), allocatable :: irhs
    real (PREC) :: rmse_default, rmse
    real (PREC) :: l1_ratio
    real (PREC) :: alpha, alpha_default
    real (PREC) :: const
    integer :: Nobs, Nrhs, Nconst, i, j
    type (status_t) :: status
    type (lm_config) :: conf_lm
    type (enet_result) :: res
    type (enet_config) :: conf
    real (PREC), parameter :: atol=1.0e-12_PREC, rtol=1.0e-4_PREC
    logical :: values_ok
    logical, dimension(:), allocatable :: mask, mask_lhs
    integer, parameter :: CV_N = 50
    real (PREC), parameter :: ALPHA_EPS = 1.0e-5_PREC

    call set_seed (1234)

    ! --- Create data ---

    Nrhs = 10
    Nobs = 103

    call random_sample (nobs=Nobs, nrhs=Nrhs, nlhs=Nlhs, X=X, Y=Y, &
        coefs=coefs_true, add_intercept=.true., intercept=intercept_true, &
        status=status)

    allocate (Y_T(Nlhs,Nobs))
    Y_T(:,:) = transpose (Y)

    ! --- Default CV ---

    l1_ratio = 0.9

    conf = enet_config (cv_n=CV_N, alpha_eps=ALPHA_EPS)
    status = NF_STATUS_UNDEFINED
    call enet_cv (conf, X, Y, alpha_default, l1_ratio, rmse=rmse_default, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "ENET_CV: default config"

    ! Fit coefs using CV'd alpha
    status = NF_STATUS_UNDEFINED
    allocate (coefs_default(Nrhs,Nlhs), intercept_default(Nlhs))
    call enet_fit (conf, X, Y, alpha_default, l1_ratio, intercept_default, &
        coefs_default, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "ENET_FIT: default config"

    ! --- Pre-specified ALPHA grid ---

    allocate (alphas(conf%alpha_n))
    allocate (folds_ifrom(conf%cv_n), folds_size(conf%cv_n))

    ! Create CV chunk indices
    call split_uniform (Nobs, conf%cv_n, folds_ifrom, folds_size, status)

    ! Create the same grid of ALPHAS used for CV
    conf = enet_config (cv_n=CV_N, alpha_eps=ALPHA_EPS)
    call create_alpha_grid_cv (conf, X, Y, l1_ratio, folds_ifrom, folds_size, &
        alphas, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "Create alpha grid for CV"

    ! --- Transposed X ---

    allocate (X_T(Nrhs, Nobs))
    X_T(:,:) = transpose (X)

    status = NF_STATUS_UNDEFINED
    conf = enet_config (cv_n=CV_N, alpha_eps=ALPHA_EPS, trans_x=.true.)
    call enet_cv (conf, X_T, Y, alpha, l1_ratio, rmse=rmse, status=status)

    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol) .and. &
        all_close (rmse, rmse_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: transposed X"

    ! Fit coefs using CV'd alpha
    allocate (coefs(Nrhs,Nlhs), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_T, Y, alpha, l1_ratio, intercept, coefs, status=status)

    values_ok = all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed X"
    deallocate (coefs, intercept)

    ! Fit using result object
    status = NF_STATUS_UNDEFINED
    ! Reset result object
    res = enet_result (conf=conf)
    call enet_fit (conf, X_T, Y, alpha, l1_ratio, res=res, status=status)
    values_ok = all_close (res%coefs_multi, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (res%intercept_multi, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed X, RES present"

    ! --- Transposed Y ---

    conf = enet_config (cv_n=CV_N, alpha_eps=ALPHA_EPS, trans_y=.true.)
    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X, Y_T, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (rmse, rmse_default, atol=atol, rtol=rtol) .and. &
        all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: transposed Y"

    ! Fit coefs using CV'd alpha
    allocate (coefs(Nrhs,Nlhs), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X, Y_T, alpha, l1_ratio, intercept, coefs, status=status)
    values_ok = all_close(intercept, intercept_default, atol=atol, rtol=rtol) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed Y"
    deallocate (coefs, intercept)

    ! Fit using result object
    status = NF_STATUS_UNDEFINED
    res = enet_result (conf=conf)
    call enet_fit (conf, X, Y_T, alpha, l1_ratio, res=res, status=status)
    values_ok = all_close (res%coefs_multi, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (res%intercept_multi, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed Y, RES present"

    ! --- Transposed X and Y ---

    conf = enet_config (cv_n=CV_N, alpha_eps=ALPHA_EPS, trans_x=.true., trans_y=.true.)
    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    call enet_cv (conf, X_T, Y_T, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (rmse, rmse_default, atol=atol, rtol=rtol) .and. &
        all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: transposed X and Y"

    ! Fit coefs using CV'd alpha
    allocate (coefs(Nrhs,Nlhs), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_T, Y_T, alpha, l1_ratio, intercept, coefs, status=status)
    values_ok = all_close(intercept, intercept_default, atol=atol, rtol=rtol) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed X and Y"
    deallocate (coefs, intercept)

    ! Fit using result object
    status = NF_STATUS_UNDEFINED
    res = enet_result (conf=conf)
    call enet_fit (conf, X_T, Y_T, alpha, l1_ratio, res=res, status=status)
    values_ok = all_close (res%coefs_multi, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (res%intercept_multi, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed X and Y, RES present"

    ! ---- Transposed COEFS ----

    allocate (coefs_T(Nlhs,Nrhs), coefs(Nrhs,Nlhs), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    conf = enet_config (trans_coefs=.true.)
    call enet_fit (conf, X, Y, alpha_default, l1_ratio, intercept, coefs_T, status=status)
    coefs(:,:) = transpose (coefs_T)
    values_ok = all_close(intercept, intercept_default, atol=atol, rtol=rtol) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed COEFS"

    ! Transposed COEFS, RES present
    coefs_T(:,:) = huge(0.0_PREC)
    res = enet_result (conf=conf)
    call enet_fit (conf, X, Y, alpha_default, l1_ratio, intercept, coefs_T, res=res, status=status)
    coefs(:,:) = transpose (coefs_T)
    values_ok = all_close(intercept, intercept_default, atol=atol, rtol=rtol) &
        .and. all_close (coefs, coefs_default, atol=atol, rtol=rtol) &
        .and. all_close (res%intercept_multi, intercept, atol=atol, rtol=rtol) &
        .and. all_close (res%coefs_multi, coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed COEFS; RES, COEFS and INTERCEPT present"

    deallocate (coefs_T, coefs, intercept)

    ! --- Transposed COEFS, transposed X ---

    allocate (coefs_T(Nlhs,Nrhs), coefs(Nrhs,Nlhs), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    conf = enet_config (trans_x=.true., trans_coefs=.true.)
    call enet_fit (conf, X_T, Y, alpha, l1_ratio, intercept, coefs_T, status=status)
    coefs(:,:) = transpose (coefs_T)
    values_ok = all_close(intercept, intercept_default, atol=atol, rtol=rtol) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed COEFS, transposed X"
    deallocate (coefs_T, coefs, intercept)

    ! --- Transposed COEFS, transposed Y ---

    allocate (coefs_T(Nlhs,Nrhs), coefs(Nrhs,Nlhs), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    conf = enet_config (trans_y=.true., trans_coefs=.true.)
    call enet_fit (conf, X, Y_T, alpha, l1_ratio, intercept, coefs_T, status=status)
    coefs(:,:) = transpose (coefs_T)
    values_ok = all_close(intercept, intercept_default, atol=atol, rtol=rtol) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed COEFS, transposed Y"
    deallocate (coefs_T, coefs, intercept)

    ! --- Transposed COEFS, transposed X and Y---

    allocate (coefs_T(Nlhs,Nrhs), coefs(Nrhs,Nlhs), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    conf = enet_config (trans_x=.true., trans_y=.true., trans_coefs=.true.)
    call enet_fit (conf, X_T, Y_T, alpha, l1_ratio, intercept, coefs_T, status=status)
    coefs(:,:) = transpose (coefs_T)
    values_ok = all_close(intercept, intercept_default, atol=atol, rtol=rtol) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed COEFS, transposed X and Y"
    deallocate (coefs_T, coefs, intercept)

    ! --- Create X with constants ---

    Nconst = 2
    allocate (mask(Nrhs + Nconst), source=.true.)
    mask(4) = .false.
    mask(9) = .false.
    allocate (X_alt(Nobs, Nrhs + Nconst))

    j = 0
    do i = 1, size(mask)
        if (mask(i)) then
            j = j + 1
            X_alt(:,i) = X(:,j)
        else
            call random_number (const)
            X_alt(:,i) = const
        end if
    end do

    if (allocated(X_T)) deallocate (X_T)
    allocate (X_T(Nrhs + Nconst, Nobs))
    X_T(:,:) = transpose (X_alt)

    ! --- X with constants ---

    status = NF_STATUS_UNDEFINED
    conf = enet_config (cv_n=CV_N, alpha_eps=ALPHA_EPS)
    call enet_cv (conf, X_alt, Y, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol) .and. &
        all_close (rmse, rmse_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: add. const. in X"

    allocate (coefs_all(Nrhs+Nconst, Nlhs), intercept(Nlhs))

    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_alt, Y, alpha, l1_ratio, intercept, coefs_all, status=status)
    ! Reduce to coefs of non-constant variables
    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all(coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X"

    deallocate (coefs, coefs_all, coefs_const, intercept)

    ! --- X with constants, transposed ---

    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    conf = enet_config (cv_n=CV_N, alpha_eps=ALPHA_EPS, trans_x=.true.)
    call enet_cv (conf, X_T, Y, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol) .and. &
        all_close (rmse, rmse_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: transposed X, add. const."

    allocate (coefs_all(Nrhs+Nconst,Nlhs), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_T, Y, alpha, l1_ratio, intercept, coefs_all, status=status)

    ! Copy coefs of non-constant variables
    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed X, add. const."
    deallocate (coefs, coefs_all, coefs_const, intercept)

    ! Call with RESULT object
    status = NF_STATUS_UNDEFINED
    res = enet_result (conf=conf)
    call enet_fit (conf, X_T, Y, alpha, l1_ratio, res=res, status=status)

    ! Copy coefs of non-constant variables
    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (res%coefs_multi, coefs, mask, dim=1)
    call copy (res%coefs_multi, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (res%intercept_multi, intercept_default, atol=atol, rtol=rtol)

    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed X, add. const., RES present"

    deallocate (coefs, coefs_const)

    ! --- X with constants, transposed Y ---

    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    conf = enet_config (cv_n=CV_N, alpha_eps=ALPHA_EPS, trans_y=.true.)
    call enet_cv (conf, X_alt, Y_T, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol) .and. &
        all_close (rmse, rmse_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: transposed Y, add. const. in X"

    allocate (coefs_all(Nrhs+Nconst,Nlhs), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_alt, Y_T, alpha_default, l1_ratio, intercept, coefs_all, status=status)

    ! Copy coefs of non-constant variables
    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed Y, add. const. in X"
    deallocate (coefs, coefs_all, coefs_const, intercept)

    ! --- X with constants, transposed X and Y ---

    status = NF_STATUS_UNDEFINED
    alpha = huge(0.0_PREC)
    rmse = huge(0.0_PREC)
    conf = enet_config (cv_n=CV_N, alpha_eps=ALPHA_EPS, trans_x=.true., trans_y=.true.)
    call enet_cv (conf, X_T, Y_T, alpha, l1_ratio, rmse=rmse, status=status)
    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol) .and. &
        all_close (rmse, rmse_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: transposed Y, add. const. in X"

    allocate (coefs_all(Nrhs+Nconst,Nlhs), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    call enet_fit (conf, X_T, Y_T, alpha_default, l1_ratio, intercept, coefs_all, status=status)
    ! Copy coefs of non-constant variables
    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed X and Y; add. const. in X"
    deallocate (coefs, coefs_all, coefs_const, intercept)

    ! --- X with constants, transposed COEFS ---

    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    conf = enet_config (trans_coefs=.true.)
    call enet_fit (conf, X_alt, Y, alpha_default, l1_ratio, intercept, coefs_T, status=status)

    ! Copy coefs of non-constant variables
    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    coefs_all(:,:) = transpose (coefs_T)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed COEFS, add. const. in X"
    deallocate (coefs, coefs_all, coefs_T, coefs_const, intercept)

    ! --- X with constants, transposed COEFS and X ---

    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    conf = enet_config (trans_x=.true., trans_coefs=.true.)
    call enet_fit (conf, X_T, Y, alpha_default, l1_ratio, intercept, coefs_T, status=status)

    ! Copy coefs of non-constant variables
    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    coefs_all(:,:) = transpose (coefs_T)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed COEFS and X, add. const. in X"
    deallocate (coefs, coefs_all, coefs_T, coefs_const, intercept)

    ! --- X with constants, transposed COEFS and Y ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    conf = enet_config (trans_y=.true., trans_coefs=.true.)
    call enet_fit (conf, X_alt, Y_T, alpha_default, l1_ratio, intercept, coefs_T, status=status)

    ! Copy coefs of non-constant variables
    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    coefs_all(:,:) = transpose (coefs_T)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed COEFS and Y, add. const. in X"
    deallocate (coefs, coefs_all, coefs_T, coefs_const, intercept)

    ! --- X with constants, transposed COEFS, X and Y ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    status = NF_STATUS_UNDEFINED
    conf = enet_config (trans_x=.true., trans_y=.true., trans_coefs=.true.)
    call enet_fit (conf, X_T, Y_T, alpha_default, l1_ratio, intercept, coefs_T, status=status)

    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    coefs_all(:,:) = transpose (coefs_T)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all_close (intercept, intercept_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: transposed COEFS, X and Y; add. const. in X"
    deallocate (coefs, coefs_all, coefs_T, coefs_const, intercept)

    ! --- Pre-center variables ---

    ! For the remaining calls we need to pass the pre-computed ALPHAS,
    ! since otherwise the max. grid point will differ as X and Y are not
    ! standardized within CV chunks!

    allocate (Y_alt(Nobs,Nlhs), source=Y)

    allocate (mean_x_alt(Nrhs + Nconst), scale_x_alt(Nrhs + Nconst))
    allocate (mean_y(Nlhs))
    call standardize (X_alt, dim=1, mean_x=mean_x_alt, scale_x=scale_x_alt, &
        center=.true., scale=.false., status=status)
    call standardize (Y_alt, mean_x=mean_y, center=.true., scale=.false., &
        dim=1, status=status)

    X_T(:,:) = transpose(X_alt)
    Y_T(:,:) = transpose(Y_alt)

    ! --- X, Y with constants and pre-centered variables ---

    alpha = 0.0
    rmse = 0.0
    status = NF_STATUS_UNDEFINED
    ! Variables are pre-centered, no need to center in ENET routines
    conf = enet_config (cv_n=CV_N, alpha_eps=ALPHA_EPS, center=.false.)
    call enet_cv (conf, X_alt, Y_alt, alpha, l1_ratio, rmse=rmse, &
        alphas=alphas, status=status)
    ! It's not entirely clear that these should be identical since the
    ! RMSEs across chunks will differ, and hence in principle a different
    ! point on the ALPHAS grid could be optiomal.
    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: add. const. in X, pre-centered variables"

    status = NF_STATUS_UNDEFINED
    allocate (coefs_all(Nrhs+Nconst,Nlhs), intercept(Nlhs))
    call enet_fit (conf, X_alt, Y_alt, alpha, l1_ratio, intercept, coefs_all, status=status)

    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X, pre-centered variables"
    deallocate (coefs_all, coefs, coefs_const, intercept)

    ! Test with RESULT present
    status = NF_STATUS_UNDEFINED
    conf = enet_config (center=.false.)
    call enet_fit (conf, X_alt, Y_alt, alpha, l1_ratio, res=res, status=status)

    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (res%coefs_multi, coefs, mask, dim=1)
    call copy (res%coefs_multi, coefs_const, .not. mask, dim=1)
    values_ok = all(coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(res%intercept_multi) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X, pre-centered variables, RES present"
    deallocate (coefs, coefs_const)

    ! --- X^T, Y with constants and pre-centered variables ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_all(Nrhs+Nconst,Nlhs), intercept(Nlhs))
    conf = enet_config (center=.false., trans_x=.true.)
    call enet_fit (conf, X_T, Y_alt, alpha, l1_ratio, intercept, coefs_all, status=status)

    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X^T, pre-centered variables"
    deallocate (coefs_all, coefs, coefs_const, intercept)

    ! --- X, Y^T with constants and pre-centered variables ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_all(Nrhs+Nconst,Nlhs), intercept(Nlhs))
    conf = enet_config (center=.false., trans_y=.true.)
    call enet_fit (conf, X_alt, Y_T, alpha, l1_ratio, intercept, coefs_all, status=status)

    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X, Y^T, pre-centered variables"
    deallocate (coefs_all, coefs, coefs_const, intercept)

    ! --- X^T, Y^T with constants and pre-centered variables --

    status = NF_STATUS_UNDEFINED
    allocate (coefs_all(Nrhs+Nconst,Nlhs), intercept(Nlhs))
    conf = enet_config (center=.false., trans_x=.true., trans_y=.true.)
    call enet_fit (conf, X_T, Y_T, alpha, l1_ratio, intercept, coefs_all, status=status)

    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X^T, Y^T, pre-centered variables"
    deallocate (coefs_all, coefs, coefs_const, intercept)

    ! --- COEFS^T, X, Y with constants and pre-centered variables ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    conf = enet_config (trans_coefs=.true., center=.false.)
    call enet_fit (conf, X_alt, Y_alt, alpha, l1_ratio, intercept, coefs_T, status=status)

    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    coefs_all(:,:) = transpose (coefs_T)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: COEFS^T, add. const. in X, pre-centered variables"
    deallocate (coefs_all, coefs, coefs_T, coefs_const, intercept)

    ! --- COEFS^T, X^T, Y with constants and pre-centered variables ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    conf = enet_config (trans_coefs=.true., center=.false., trans_x=.true.)
    call enet_fit (conf, X_T, Y_alt, alpha, l1_ratio, intercept, coefs_T, status=status)

    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    coefs_all(:,:) = transpose (coefs_T)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: COEFS^T, add. const. in X^T, pre-centered variables"
    deallocate (coefs_all, coefs, coefs_T, coefs_const, intercept)

    ! --- X, Y^T with constants and pre-centered variables ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    conf = enet_config (trans_coefs=.true., center=.false., trans_y=.true.)
    call enet_fit (conf, X_alt, Y_T, alpha, l1_ratio, intercept, coefs_T, status=status)

    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    coefs_all(:,:) = transpose (coefs_T)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: COEFS^T, add. const. in X, Y^T, pre-centered variables"
    deallocate (coefs_all, coefs, coefs_T, coefs_const, intercept)

    ! --- COEFS^T, X^T, Y^T with constants and pre-centered variables --

    status = NF_STATUS_UNDEFINED
    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    conf = enet_config (trans_coefs=.true.,center=.false., trans_x=.true., trans_y=.true.)
    call enet_fit (conf, X_T, Y_T, alpha, l1_ratio, intercept, coefs_T, status=status)

    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    coefs_all(:,:) = transpose (coefs_T)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)
    values_ok = all (coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X, Y^T, pre-centered variables"
    deallocate (coefs_all, coefs, coefs_T, coefs_const, intercept)

    ! --- Additionally pre-scale X ---

    ! X_alt is already center from before
    call standardize (X_alt, dim=1, scale_x=scale_x_alt, center=.false., &
        scale=.true., skip_const=.true., status=status)
    X_T(:,:) = transpose (X_alt)

    ! --- X pre-centered and pre-scaled, Y pre-centered ---

    alpha = 0.0
    rmse = 0.0
    status = NF_STATUS_UNDEFINED
    conf = enet_config (cv_n=CV_N, alpha_eps=ALPHA_EPS, center=.false., scale=.false.)
    call enet_cv (conf, X_alt, Y_alt, alpha, l1_ratio, rmse=rmse, &
        alphas=alphas, status=status)
    ! It's not entirely clear that these should be identical since the
    ! RMSEs across chunks will differ, and hence in principle a different
    ! point on the ALPHAS grid could be optiomal.
    values_ok = all_close (alpha, alpha_default, atol=atol, rtol=rtol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_CV: add. const. in X, pre-centered, pre-scaled variables"

    status = NF_STATUS_UNDEFINED
    allocate (coefs_all(Nrhs+Nconst,Nlhs), intercept(Nlhs))
    call enet_fit (conf, X_alt, Y_alt, alpha, l1_ratio, intercept, coefs_all, status=status)
    ! Undo scaling of coefficients
    coefs_all(:,:) = coefs_all / spread (scale_x_alt, dim=2, ncopies=Nlhs)
    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)

    values_ok = all(coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X, pre-centered, pre-scaled variables"
    deallocate (coefs_all, coefs, coefs_const, intercept)

    ! Present RES
    status = NF_STATUS_UNDEFINED
    res = enet_result (conf=conf)
    call enet_fit (conf, X_alt, Y_alt, alpha, l1_ratio, res=res, status=status)
    ! Undo scaling of coefficients
    allocate (coefs_all, source=res%coefs_multi)
    coefs_all(:,:) = coefs_all / spread (scale_x_alt, dim=2, ncopies=Nlhs)
    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)

    values_ok = all(coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(res%intercept_multi) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X, pre-centered, pre-scaled variables"
    deallocate (coefs_all, coefs, coefs_const)

    ! --- X^T pre-centered and pre-scaled, Y pre-centered ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_all(Nrhs+Nconst,Nlhs), intercept(Nlhs))
    conf = enet_config (trans_x=.true., center=.false., scale=.false.)
    call enet_fit (conf, X_T, Y_alt, alpha, l1_ratio, intercept, coefs_all, status=status)
    ! Undo scaling of coefficients
    coefs_all(:,:) = coefs_all / spread (scale_x_alt, dim=2, ncopies=Nlhs)
    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)

    values_ok = all(coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X^T, pre-centered, pre-scaled variables"
    deallocate (coefs_all, coefs, coefs_const, intercept)

    ! --- X pre-centered and pre-scaled, Y^T pre-centered ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_all(Nrhs+Nconst,Nlhs), intercept(Nlhs))
    conf = enet_config (trans_y=.true., center=.false., scale=.false.)
    call enet_fit (conf, X_alt, Y_T, alpha, l1_ratio, intercept, coefs_all, status=status)
    ! Undo scaling of coefficients
    coefs_all(:,:) = coefs_all / spread (scale_x_alt, dim=2, ncopies=Nlhs)
    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)

    values_ok = all(coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X, Y^T, pre-centered, pre-scaled " // &
            "variables"
    deallocate (coefs_all, coefs, coefs_const, intercept)

    ! --- X^T pre-centered and pre-scaled, Y^T pre-centered ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_all(Nrhs+Nconst,Nlhs), intercept(Nlhs))
    conf = enet_config (trans_x=.true., trans_y=.true., center=.false., scale=.false.)
    call enet_fit (conf, X_T, Y_T, alpha, l1_ratio, intercept, coefs_all, status=status)
    ! Undo scaling of coefficients
    coefs_all(:,:) = coefs_all / spread (scale_x_alt, dim=2, ncopies=Nlhs)
    allocate (coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)

    values_ok = all(coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: add. const. in X^T, Y^T, pre-centered, pre-scaled " // &
            "variables"
    deallocate (coefs_all, coefs, coefs_const, intercept)

    ! --- COEFS^T, X pre-centered and pre-scaled, Y pre-centered ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    conf = enet_config (trans_coefs=.true., center=.false., scale=.false.)
    call enet_fit (conf, X_alt, Y_alt, alpha, l1_ratio, intercept, coefs_T, status=status)
    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    ! Undo scaling of coefficients
    coefs_all(:,:) = transpose(coefs_T) / spread (scale_x_alt, dim=2, ncopies=Nlhs)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)

    values_ok = all(coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: COEFS^T, add. const. in X, pre-centered, pre-scaled " &
            // "variables"
    deallocate (coefs_all, coefs, coefs_T, coefs_const, intercept)

    ! --- COEFS^T, X^T pre-centered and pre-scaled, Y pre-centered ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    conf = enet_config (trans_coefs=.true., trans_x=.true., center=.false., scale=.false.)
    call enet_fit (conf, X_T, Y_alt, alpha, l1_ratio, intercept, coefs_T, status=status)
    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    ! Undo scaling of coefficients
    coefs_all(:,:) = transpose(coefs_T) / spread (scale_x_alt, dim=2, ncopies=Nlhs)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)

    values_ok = all(coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: COEFS^T, add. const. in X^T, pre-centered, " // &
            "pre-scaled variables"
    deallocate (coefs_all, coefs, coefs_T, coefs_const, intercept)

    ! --- X pre-centered and pre-scaled, Y^T pre-centered ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    conf = enet_config (trans_coefs=.true., trans_y=.true., center=.false., scale=.false.)
    call enet_fit (conf, X_alt, Y_T, alpha, l1_ratio, intercept, coefs_T, status=status)
    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    ! Undo scaling of coefficients
    coefs_all(:,:) = transpose (coefs_T) / spread (scale_x_alt, dim=2, ncopies=Nlhs)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)

    values_ok = all(coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: COEFS^T, add. const. in X, Y^T, pre-centered, " // &
            "pre-scaled variables"
    deallocate (coefs_all, coefs, coefs_T, coefs_const, intercept)

    ! --- COEFS^T, X^T pre-centered and pre-scaled, Y^T pre-centered ---

    status = NF_STATUS_UNDEFINED
    allocate (coefs_T(Nlhs,Nrhs+Nconst), intercept(Nlhs))
    conf = enet_config (trans_coefs=.true., trans_x=.true., trans_y=.true., center=.false., scale=.false.)
    call enet_fit (conf, X_T, Y_T, alpha, l1_ratio, intercept, coefs_T, status=status)
    allocate (coefs_all(Nrhs+Nconst,Nlhs), coefs(Nrhs,Nlhs), coefs_const(Nconst,Nlhs))
    ! Undo scaling of coefficients
    coefs_all(:,:) = transpose (coefs_T) / spread (scale_x_alt, dim=2, ncopies=Nlhs)
    call copy (coefs_all, coefs, mask, dim=1)
    call copy (coefs_all, coefs_const, .not. mask, dim=1)

    values_ok = all(coefs_const == 0.0_PREC) .and. &
        all_close (coefs, coefs_default, atol=atol, rtol=rtol) .and. &
        all(abs(intercept) < atol)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "ENET_FIT: COEFS^T, add. const. in X^T, Y^T, pre-centered, " // &
            "pre-scaled variables"
    deallocate (coefs_all, coefs, coefs_T, coefs_const, intercept)

end subroutine

subroutine test_ridge_1d_args ()
    !*  Unit tests for rigde regression argument checking

    real (PREC), dimension(:,:), allocatable :: x
    real (PREC), dimension(:), allocatable :: y
    real (PREC), dimension(:), allocatable :: beta

    real (PREC) :: alpha
    integer :: nrhs, nobs, ncoefs
    type (status_t) :: status
    type (enet_config) :: conf

    alpha = 0.0

    ! === incompatible X, Y arguments ===
    nrhs = 3
    nobs = 5
    ! Note: constant added by default
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs+1), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    call ridge (conf, x, y, alpha, beta, status=status)
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
    conf = enet_config ()
    conf%trans_x = .true.
    call ridge (conf, x, y, alpha, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X^T, Y"

    deallocate (x, y, beta)

    ! === incompabible X, COEFS ===
    nrhs = 3
    nobs = 5
    ncoefs = nrhs
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs-1))

    status = NF_STATUS_UNDEFINED
    conf = enet_config()
    call ridge (conf, x, y, alpha, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X, COEFS"

    deallocate (x, y, beta)

    ! === incompabible X^T, COEFS ===
    nrhs = 3
    nobs = 5
    ncoefs = nrhs
    allocate (x(nrhs+1,nobs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = enet_config ()
    conf%trans_x = .true.
    call ridge (conf, x, y, alpha, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X^T, COEFS"

    deallocate (x, y, beta)

    ! === Invalid RCOND ===
    nrhs = 3
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    conf = enet_config ()
    call ridge (conf, x, y, alpha, beta, rcond=-0.1_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid RCOND"

    deallocate (x, y, beta)

    ! === Invalid alpha ===

    nrhs = 3
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs), beta(ncoefs))

    status = NF_STATUS_UNDEFINED
    alpha = -1.0_PREC
    conf = enet_config ()
    call ridge (conf, x, y, alpha, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid LAMBDA"

    deallocate (x, y, beta)

end subroutine

subroutine test_ridge_2d_args ()
    !*  Unit tests for ridge regression argument checking

    real (PREC), dimension(:,:), allocatable :: x, y
    real (PREC), dimension(:,:), allocatable :: beta

    integer :: nrhs, nlhs, nobs, ncoefs
    real (PREC) :: alpha
    type (status_t) :: status
    type (enet_config) :: conf

    alpha = 0.0

    ! === incompatible X, Y arguments ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs
    allocate (x(nobs, nrhs), y(nobs+1,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = enet_config ()
    call ridge (conf, x, y, alpha, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X, Y"

    deallocate (x, y, beta)

    ! === incompatible X^T, Y ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs
    allocate (x(nrhs, nobs), y(nobs+1,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = enet_config (trans_x=.true.)
    call ridge (conf, x, y, alpha, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X^T, Y"

    deallocate (x, y, beta)

    ! === incompabible X, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs-1,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = enet_config ()
    call ridge (conf, x, y, alpha, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X, COEFS"

    deallocate (x, y, beta)

    ! === incompabible X^T, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs
    allocate (x(nrhs+1,nobs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = enet_config(trans_x=.true.)
    call ridge (conf, x, y, alpha, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays X^T, COEFS"

    deallocate (x, y, beta)

    ! === incompabible Y, COEFS ===
    nrhs = 3
    nlhs = 1
    nobs = 5
    ncoefs = nrhs
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs-1))

    status = NF_STATUS_UNDEFINED
    conf = enet_config ()
    call ridge (conf, x, y, alpha, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable arrays Y, COEFS"

    deallocate (x, y, beta)

    ! === Invalid RCOND ===
    nrhs = 3
    nlhs = 2
    nobs = 5
    ncoefs = nrhs + 1
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    conf = enet_config ()
    call ridge (conf, x, y, alpha, beta, rcond=-0.1_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid RCOND"

    deallocate (x, y, beta)

    ! === Invalid LAMBDA ===
    nrhs = 3
    nlhs = 2
    nobs = 5
    ncoefs = nrhs
    allocate (x(nobs, nrhs), y(nobs,nlhs), beta(ncoefs,nlhs))

    status = NF_STATUS_UNDEFINED
    alpha = -1.0
    conf = enet_config ()
    call ridge (conf, x, y, alpha, beta, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid LAMBDA"

    deallocate (x, y, beta)

end subroutine

subroutine test_ridge_1d ()
    !*  Unit tests for various OLS regression problems using the 1d API

    real (PREC), dimension(:,:), allocatable :: x
    real (PREC), dimension(:), allocatable :: y, beta, beta_true, work
    real (PREC), dimension(:), allocatable :: beta_ols
    integer :: nrhs, ncoefs, nobs
    integer :: rank, rank_ols
    real (PREC) :: intercept, intercept_ols, intercept_true
    real (PREC), parameter :: Rtol=1.0e-3_PREC, ATOL=0.0_PREC
    logical :: values_ok, trans_x
    ! Arguments for ALL_CLOSE()
    real (PREC) :: rsq, alpha, ssq, ssq_ols
    type (status_t) :: status
    type (lm_config) :: conf_lm
    type (lm_result), allocatable :: lm_ols
    type (enet_result), allocatable :: glm
    type (enet_config) :: conf

    call set_seed (456)

    ! === exactly determined system with NOBS = NRHS, no intercept ===

    nrhs = 10
    nobs = 10
    ncoefs = nrhs

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_true(ncoefs))
    allocate (beta_ols(ncoefs))

    call random_sample (nobs, nrhs, X, y, beta_true, add_intercept=.false., &
        var_error=0.0_PREC, status=status)

    allocate (glm, lm_ols)

    conf_lm = lm_config (add_intercept=.false.)
    call ols (conf_lm, x, y, beta_ols, rank=rank_ols, res=lm_ols)

    ! Test with alpha = 0.0
    status = NF_STATUS_UNDEFINED
    alpha = 0.0
    conf = enet_config (center=.false., scale=.false.)
    call ridge (conf, x, y, alpha, beta, rank=rank, res=glm, status=status)

    values_ok = all_close (beta, beta_true, rtol=RTOL, atol=ATOL)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Estim: Exactly identified system, LAMBDA=0"

    ! Test with alpha > 0.0
    status = NF_STATUS_UNDEFINED
    alpha = 0.1
    beta(:) = 0.0
    call ridge (conf, x, y, alpha, beta, rank=rank, res=glm, status=status)

    ssq = sum(beta**2.0)
    ssq_ols = sum(beta_ols**2.0)
    if (.not. (status == NF_STATUS_OK .and. ssq <= ssq_ols)) &
        error stop "Estim: Exactly identified system, LAMBDA>0"

    deallocate (x, y, beta, beta_true, beta_ols)

    ! === overdetermined system without intercept ===

    nrhs = 5
    nobs = 500
    ncoefs = nrhs

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_true(ncoefs))
    allocate (beta_ols(ncoefs))

    call random_sample (nobs, nrhs, X, y, beta_true, add_intercept=.false., &
        var_error=1.0e-3_PREC, status=status)

    conf_lm = lm_config (add_intercept=.false.)
    call ols (conf_lm, x, y, beta_ols, rank=rank_ols)

    ! Lambda = 0
    status = NF_STATUS_UNDEFINED
    alpha = 0.0
    conf = enet_config (center=.false., scale=.false.)
    call ridge (conf, x, y, alpha, beta, rank=rank, status=status)

    values_ok = all_close (beta, beta_ols, rtol=RTOL, atol=ATOL)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Estim: Overdet. system with small residuals, LAMBDA=0"

    ! Lambda > 0
    status = NF_STATUS_UNDEFINED
    alpha = 0.1d0
    call ridge (conf, x, y, alpha, beta, rank=rank, status=status)

    ssq = sum(beta**2.0)
    ssq_ols = sum(beta_ols**2.0)
    if (.not. (status == NF_STATUS_OK .and. ssq <= ssq_ols)) &
        error stop "Estim: Overdet. system with small residuals, LAMBDA>0"

    deallocate (x, y, beta, beta_true, beta_ols)

    ! === overdetermined system with intercept ===

    nrhs = 5
    nobs = 500
    ncoefs = nrhs

    allocate (x(nobs,nrhs), y(nobs), beta(ncoefs), beta_true(ncoefs))
    allocate (beta_ols(ncoefs))

    call random_sample (nobs, nrhs, X, y, beta_true, add_intercept=.true., &
        intercept=intercept_true, var_error=1.0e-3_PREC, status=status)

    conf_lm = lm_config (add_intercept=.true.)
    call ols (conf_lm, x, y, beta_ols, intercept=intercept_ols, rank=rank_ols)

    ! Lambda = 0.0
    status = NF_STATUS_UNDEFINED
    alpha = 0.0
    conf = enet_config (center=.true., scale=.false.)
    call ridge (conf, x, y, alpha, beta, intercept=intercept, rank=rank, status=status)

    values_ok = all_close (beta, beta_ols, rtol=RTOL, atol=ATOL) .and. &
        all_close (intercept, intercept_ols, atol=ATOL, rtol=RTOL)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Estim: Overdet. system with small residuals, LAMBDA=0, " // &
            "CENTER=.TRUE."

    ! Lambda > 0.0
    status = NF_STATUS_UNDEFINED
    alpha = 0.5
    conf = enet_config (center=.true., scale=.true.)
    call ridge (conf, x, y, alpha, beta, intercept=intercept, rank=rank, status=status)

    ssq = sum(beta**2.0)
    ssq_ols = sum(beta_ols**2.0)
    if (.not. (status == NF_STATUS_OK .and. ssq <= ssq_ols)) &
        error stop "Estim: Overdet. system with small residuals, LAMBDA>0, " // &
            "CENTER=.TRUE."

    deallocate (x, y, beta, beta_true, beta_ols)

    ! === overdetermined system with intercept, transposed X ===

    nrhs = 2
    nobs = 500
    ncoefs = nrhs

    allocate (x(nrhs,nobs), y(nobs), beta(ncoefs), beta_true(ncoefs))
    allocate (beta_ols(ncoefs))

    call random_sample (nobs, nrhs, X, y, beta_true, add_intercept=.true., &
        intercept=intercept_true, trans_x=.true., var_error=1.0e-2_PREC, &
        status=status)

   trans_x = .true.

    conf_lm = lm_config (add_intercept=.true., trans_x=trans_x)
    call ols (conf_lm, x, y, beta_ols, intercept=intercept_ols, rank=rank_ols)

    ! Lambda = 0.0
    status = NF_STATUS_UNDEFINED
    alpha = 0.0
    conf = enet_config (center=.true., trans_x=trans_x)
    call ridge (conf, x, y, alpha, beta, rank=rank, intercept=intercept, status=status)

    values_ok = all_close (beta, beta_ols, rtol=RTOL, atol=ATOL) .and. &
        all_close (intercept, intercept_ols, atol=ATOL, rtol=RTOL)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Estim: Overdet. system with small resid, LAMBDA=0, " // &
            "CENTER=.TRUE., TRANS_X=.TRUE."

    ! Lambda > 0.0
    status = NF_STATUS_UNDEFINED
    alpha = 0.1
    call ridge (conf, x, y, alpha, beta, rank=rank, status=status)

    ssq = sum(beta**2.0)
    ssq_ols = sum(beta_ols**2.0)
    if (.not. (status == NF_STATUS_OK .and. ssq <= ssq_ols)) &
        error stop "Estim: Overdet. system with small resid, LAMBDA>0, " // &
            "CENTER=.TRUE., TRANS_X=.TRUE."

    deallocate (x, y, beta, beta_true, beta_ols)

end subroutine

end
