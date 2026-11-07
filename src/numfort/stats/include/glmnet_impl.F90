



pure subroutine check_result (res, status, X, Y_multi, trans_x, trans_y, msg)
    type (enet_result), intent(in) :: res
    type (status_t), intent(out) :: status
    real (PREC), intent(in), dimension(:,:), optional :: X
    real (PREC), intent(in), dimension(:,:), optional :: Y_multi
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
    character (*), intent(out), optional :: msg

    integer :: Nobs, Nrhs, Nlhs
    logical :: ltrans_x, ltrans_y
    character (*), parameter :: NAME = 'GLMNET'

    call set_optional_arg (trans_x, res%conf%trans_x, ltrans_x)
    call set_optional_arg (trans_y, res%conf%trans_y, ltrans_y)

    call check_cond (res%l1_ratio >= 0.0_PREC .and. res%l1_ratio <= 1.0_PREC, &
        NAME, 'Invalid L1_RATIO attribute', status, msg)
    if (status /= NF_STATUS_OK) return

    call check_cond (res%alpha >= 0.0_PREC, NAME, 'Invalid ALPHA attribute', &
        status, msg)
    if (status /= NF_STATUS_OK) return

    call check_cond (res%nobs > 0, NAME, 'Uninitialized NOBS attribute', &
        status, msg)
    if (status /= NF_STATUS_OK) return

    if (present(X)) then
        ! Check that X and the model for which RES was estimated have the
        ! same number of variabled (but allow for different number of
        ! obs., e.g. for prediction)
        call get_regressor_dims (X, ltrans_x, nvars=Nrhs, nobs=Nobs)
        call check_cond (Nrhs == res%nrhs, NAME, 'X: Non-conformable array', &
            status, msg)
    end if

    if (present(Y_multi)) then
        ! Check that Y_MULTI and the model for which RES was estimated have
        ! the same number of LHS variables.
        call get_regressor_dims (Y_multi, ltrans_y, nvars=Nlhs)
        call check_cond (Nlhs == res%nlhs, NAME, 'Y_MULTI: Non-conformable array', &
            status, msg)
    end if

    if (res%nlhs == 1) then
        call check_cond (allocated(res%coefs), NAME, 'COEFS attribute not allocated', &
            status, msg)
        if (status /= NF_STATUS_OK) return

        call check_cond (size(res%coefs) == res%nrhs, NAME, &
            'COEFS attribute has wrong size', status, msg)
        if (status /= NF_STATUS_OK) return
    else
        call check_cond (allocated(res%coefs_multi), NAME, &
            'COEFS_MULTI attribute not allocated', status, msg)
        if (status /= NF_STATUS_OK) return

        call check_cond (size(res%coefs_multi,1) == res%nrhs &
            .and. size(res%coefs_multi,2) == res%nlhs, NAME, &
            'COEFS_MULTI attribute has wrong size', status, msg)
        if (status /= NF_STATUS_OK) return

        call check_cond (allocated(res%intercept_multi), NAME, &
            'INTERCEPT_MULTI attribute not allocated', status, msg)
        if (status /= NF_STATUS_OK) return

        call check_cond (size(res%intercept_multi) == res%nlhs, NAME, &
            'INTERCEPT_MULTI attribute has wrong size', status, msg)
        if (status /= NF_STATUS_OK) return
    end if

    status = NF_STATUS_OK
    if (present(msg)) msg = ''

end subroutine



pure subroutine check_conf (conf, status, msg)
    type (enet_config), intent(in) :: conf
    type (status_t), intent(out) :: status
    character (*), intent(out), optional :: msg

    character (*), parameter :: NAME = 'GLMNET'

    call check_cond (conf%tol_const > 0.0, NAME, 'CONST_TOL: Invalid value', status, msg)
    if (status /= NF_STATUS_OK) return

    call check_cond (conf%alpha_n > 0, NAME, 'ALPHA_N: Invalid value', status, msg)
    if (status /= NF_STATUS_OK) return

    call check_cond (conf%alpha_eps > 0.0, NAME, 'ALPHA_EPS: Invalid value', status, msg)
    if (status /= NF_STATUS_OK) return

    call check_cond (conf%cv_n > 0, NAME, 'CV_N: Invalid value', status, msg)
    if (status /= NF_STATUS_OK) return

    call check_cond (conf%maxiter > 0, NAME, 'MAXITER: Invalid value', status, msg)
    if (status /= NF_STATUS_OK) return

    call check_cond (conf%tol > 0.0, NAME, 'TOL: Invalid value', status, msg)
    if (status /= NF_STATUS_OK) return

    status = NF_STATUS_OK
    if (present(msg)) msg = ''

end subroutine


pure function conf_transform_data (conf, nlhs) result (res)
    type (enet_config), intent(in) :: conf
    integer, intent(in), optional :: nlhs
    logical :: res

    res = conf%trans_x .or. conf%center .or. conf%scale .or. conf%drop_const

    if (present(nlhs)) then
        res = res .or. conf%trans_y
    end if

end function



subroutine check_dims (X, y, status, trans_x, Xout, yout, XX, Xy)
    real (PREC), intent(in), dimension(:,:) :: X
    real (PREC), intent(in), dimension(:) :: y
    type (status_t), intent(out) :: status
    logical, intent(in), optional :: trans_x
    real (PREC), intent(in), dimension(:,:), optional :: Xout
    real (PREC), intent(in), dimension(:), optional :: yout
    real (PREC), intent(in), dimension(:,:), optional :: XX
    real (PREC), intent(in), dimension(:), optional :: Xy

    integer :: n, k

    status = NF_STATUS_OK

    call get_regressor_dims (X, trans_x, k, n)

    if (size(y) /= n) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (present(Xout)) then
        if (size(Xout, 1) /= n .or. size(Xout, 2) /= k) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    if (present(yout)) then
        if (size(yout) /= n) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    if (present(XX)) then
        if (size(XX, 1) /= k .or. size(XX, 2) /= k) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    if (present(Xy)) then
        if (size(Xy) /= k) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

100 continue

end subroutine



subroutine check_dims_multi (X, y, status, trans_x, trans_y, Xout, yout)
    real (PREC), intent(in), dimension(:,:) :: X
    real (PREC), intent(in), dimension(:,:) :: y
    type (status_t), intent(out) :: status
    logical, intent(in), optional :: trans_x, trans_y
    real (PREC), intent(in), dimension(:,:), optional :: Xout
    real (PREC), intent(in), dimension(:,:), optional :: yout

    integer :: Nobs, Nobs_y, Nobs_x, Npred, Nlhs

    status = NF_STATUS_OK

    call get_regressor_dims (X, trans_x, Npred, Nobs_x)
    call get_regressor_dims (Y, trans_y, Nlhs, Nobs_y)

    if (Nobs_y /= Nobs_x) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    Nobs = Nobs_x

    if (present(Xout)) then
        if (size(Xout, 1) /= Nobs .or. size(Xout, 2) /= Npred) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    if (present(Yout)) then
        if (size(Yout,1) /= Nobs .or. size(Yout,2) /= Nlhs) then
            status = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

100 continue

end subroutine



subroutine compute_moments (X, y, XX, Xy, status)
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous :: y
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: XX
    real (PREC), intent(out), dimension(:), contiguous, optional :: Xy
    type (status_t), intent(out), optional :: status

    integer :: n, k
    type (status_t) ::lstatus

    lstatus = NF_STATUS_OK

    n = size(X, 1)
    k = size(X, 2)

    if (present(XX)) then
        call gram (X, XX, trans='T', status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(Xy)) then
        if (size(Xy) /= k) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
        call BLAS_GEMV ('T', n, k, 1.0_PREC, X, n, y, 1, 0.0_PREC, Xy, 1)
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine get_alpha_grid_max (Xy, nobs, l1_ratio, alpha_max, status)
    !*  GET_ALPHA_GRID_MAX returns the value of the regularization parameter
    !   alpha such that all coefficients beta are estimated to be zero.
    !
    !   Consider the the minimization problem given by
    !       min_{b} 1/(2*N) * (y-Xb)'(y-Xb) + alpha * L1 * |b|_1
    !               + alpha * (1-L1) * 1/2 * |b|_2^2
    !   where |.|_1 and |.|_2 are the L1 and L2 norms, respectively, and
    !   the data y and X are assumed to be be centered. We additionally
    !   require that L1 > 0.
    !
    !   Assuming that all elements in b are non-zero, the first-order condition
    !   is given by
    !       0 = - X'y/N + X'Xb/N + alpha * L1 * sign(b) + alpha * (1-L1) * b
    !   Consider the case when all elements in b but one are zero, and
    !   the j-th element is b = eps > 0.
    !   The above condition is then
    !       0 = -X_j'y/N + (X_j)'(X_j) b_j/N + alpha * L1 + alpha + (1-L1) * b
    !   We want to choose alpha such that the RHS is positive for any
    !   b_j > 0 such that it is optimal to decrease b_j towards 0.
    !   A sufficient condition is that
    !       RHS > -X_j'y/N + alpha * L1 > 0
    !   which implies that we need to impose
    !       alpha > X_j'y / (L1 * N)
    !
    !   Conversely, assume that b_j < 0: the first-order condition is then
    !        0 = - X_j'y/N + (X_j)'(X_j)b_j/N - alpha * L1 + alpha * (1-L1) * b_j
    !   We want to enforce that the RHS < 0, such that increasing b_j
    !   towards zero lowers the objective. A sufficient condition is thus
    !       RHS < -X_j'y/N - alpha * L1 < 0
    !   and consequently we need to impose
    !       alpha > - X_j'y / (L1 * N)
    !   Combining these two sufficient conditions, we see that
    !       alpha > max_j {|X_j'y| / (L1 * N)}
    !   is sufficient to shrink all b_j to zero.
    real (PREC), intent(in), dimension(:), contiguous, target :: Xy
        !*  Pre-computed moment X'y (possibly for a sub-sample, such as
        !   the training sample).
    integer, intent(in) :: nobs
        !*  Number of observations used to compute X'y (e.g. number of
        !   observations in training sample).
    real (PREC), intent(in) :: l1_ratio
        !*  Elastic net mixing parameter, which is required to satisfy
        !   0 < L1_RATIO <= 1. L1_RATIO = 0 is not supported as then no grid
        !   upper bound can be computed.
    real (PREC), intent(out) :: alpha_max
        !*  Alpha required to shrink all coefficients to zero.
    type (status_t), intent(out), optional :: status
        !*  Exit code

    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Xy

    ptr_Xy(1:size(Xy),1:1) => Xy

    ! Call multi-output version
    call get_alpha_grid_max_multi (ptr_Xy, nobs, l1_ratio, alpha_max, status)

end subroutine



subroutine get_alpha_grid_max_multi (Xy, nobs, l1_ratio, alpha_max, status)
    !*  GET_ALPHA_GRID_MAX returns the value of the regularization parameter
    !   alpha such that all coefficients beta are estimated to be zero.
    !
    !   Implementation for multiple outcome variables ("tasks"). See
    !   Single-outcome variant for an exposition of the algorithm.
    real (PREC), intent(in), dimension(:,:), contiguous :: Xy
    !*  Pre-computed moment X'Y (possibly for a sub-sample, such as
    !   the training sample).
    integer, intent(in) :: nobs
    !*  Number of observations used to compute X'y (e.g. number of
    !   observations in training sample).
    real (PREC), intent(in) :: l1_ratio
    !*  Elastic net mixing parameter, which is required to satisfy
    !   0 < L1_RATIO <= 1. L1_RATIO = 0 is not supported as then no grid
    !   upper bound can be computed.
    real (PREC), intent(out) :: alpha_max
    !*  Alpha required to shrink all coefficients to zero.
    type (status_t), intent(out), optional :: status
    !*  Exit code

    real (PREC) :: Xy_max
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'GET_ALPHA_GRID_MAX'

    lstatus = NF_STATUS_OK
    ! Prevent annoying uninitialized variable compiler warnings
    alpha_max = 0.0_PREC

    ! --- Input checks ---

    call check_cond (l1_ratio > 0.0_PREC .and. l1_ratio <= 1.0_PREC, &
        NAME, 'L1_RATIO: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Implementation ---

    Xy_max = maxval(sum(Xy**2.0_PREC, dim=2))

    alpha_max = sqrt(Xy_max) / (nobs * l1_ratio)

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine create_alpha_grid (conf, X, y, l1_ratio, grid, Xy, status)
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:), contiguous, target :: y
    real (PREC), intent(in) :: l1_ratio
    real (PREC), intent(out), dimension(:), contiguous :: grid
    real (PREC), intent(in), dimension(:), contiguous, optional, target :: Xy
    type (status_t), intent(out), optional :: status

    integer :: i, nobs, k, Nvars, Nalpha
    real (PREC) :: alpha_max, tmp, resolution
    type (status_t) :: lstatus
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Xp
    real (PREC), dimension(:), pointer, contiguous :: ptr_Xy, ptr_yp

    lstatus = NF_STATUS_OK

    call check_conf (conf, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call get_regressor_dims (X, conf%trans_x, nobs, k)

    ! --- Compute X'y (optional) ---

    if (.not. present(Xy)) then
        if (conf_transform_data (conf)) then
            allocate (ptr_Xp(nobs,k))
            call transform_regressors (X, ptr_Xp, trans=conf%trans_x, &
                center=conf%center, scale=conf%scale, nkeep=Nvars, status=lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100

            allocate (ptr_yp(nobs), source=y)
            call standardize (ptr_yp, center=conf%center, scale=.false.)
        else
            Nvars = k
            ptr_Xp => X
            ptr_yp => y
        end if

        allocate (ptr_Xy(Nvars))
        call compute_moments (ptr_Xp(:,1:Nvars), ptr_yp, Xy=ptr_Xy, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        ptr_Xy => Xy
    end if

    ! --- Find grid upper bound ---

    call get_alpha_grid_max (ptr_Xy, nobs, l1_ratio, alpha_max)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Create grid ---

    ! Replicate same grid for tiny alpha_max as in scikit-learn
    resolution = 10.0_PREC ** (-precision(0.0_PREC))
    if (alpha_max < resolution) then
        grid = resolution
        goto 100
    end if

    call logspace (grid, log10(alpha_max*conf%alpha_eps), log10(alpha_max))
    ! Invert sort order such that smallest alpha is last
    Nalpha = size(grid)
    do i = 1, Nalpha / 2
        tmp = grid(i)
        grid(i) = grid(Nalpha-i+1)
        grid(Nalpha-i+1) = tmp
    end do

100 continue

    call assert_dealloc_ptr (X, ptr_Xp)
    call assert_dealloc_ptr (y, ptr_yp)
    call assert_dealloc_ptr (Xy, ptr_Xy)

    if (present(status)) status = lstatus

end subroutine



subroutine create_alpha_grid_cv (conf, X, y, l1_ratio, folds_ifrom, folds_size, &
        grid, status)
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous, target :: y
    real (PREC), intent(in) :: l1_ratio
    integer, intent(in), dimension(:), contiguous :: folds_ifrom
    integer, intent(in), dimension(:), contiguous :: folds_size
    real (PREC), intent(out), dimension(:), contiguous :: grid
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer, contiguous :: ptr_y

    ptr_y(1:size(y),1:1) => y

    call create_alpha_grid_cv_multi (conf, X, ptr_y, l1_ratio, folds_ifrom, &
        folds_size, grid, status)

end subroutine



subroutine create_alpha_grid_cv_multi (conf, X, Y, l1_ratio, folds_ifrom, &
        folds_size, grid, status)
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:,:), contiguous :: Y
    real (PREC), intent(in) :: l1_ratio
    integer, intent(in), dimension(:), contiguous :: folds_ifrom
    integer, intent(in), dimension(:), contiguous :: folds_size
    real (PREC), intent(out), dimension(:), contiguous :: grid
    type (status_t), intent(out), optional :: status

    integer :: i, ifrom, n, Nobs, Nlhs, k, nchunks, Nvars, Nalpha, Ntrain
    integer :: shp(2)
    real (PREC) :: alpha_max, tmp, resolution
    real (PREC), dimension(:,:), allocatable :: X_train, Y_train, Xy
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    ! --- Input checks ---

    call check_conf (conf, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_dims (X, Y, trans_x=conf%trans_x, trans_y=conf%trans_y, &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (size(folds_ifrom) /= size(folds_size)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! --- Iteratore over chunks ---

    alpha_max = 0.0_PREC
    tmp = 0.0_PREC

    nchunks = size(folds_ifrom)

    call get_regressor_dims (X, conf%trans_x, k, Nobs)
    call get_regressor_dims (Y, conf%trans_y, nvars=Nlhs)

    do i = 1, nchunks
        n = folds_size(i)
        ifrom = folds_ifrom(i)
        Ntrain = Nobs - n

        ! Split into training and test samples
        call extract_block_alloc (X, ifrom, n, trans=conf%trans_x, x_rest=X_train)
        call extract_block_alloc (Y, ifrom, n, trans=conf%trans_y, x_rest=Y_train)

        call transform_regressors (X_train, center=.true., scale=.true., &
            drop_const=.true., tol_const=conf%tol_const, drop_na=.true., &
            nkeep=Nvars, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        call standardize (Y_train, center=.true., scale=.false., dim=1)

        shp(1) = Nvars
        shp(2) = Nlhs
        call cond_alloc (Xy, shp)
        ! Compute X'Y
        call BLAS_GEMM ('T', 'N', Nvars, Nlhs, Ntrain, 1.0_PREC, &
            X_train(:,1:Nvars), Ntrain, Y_train, Ntrain, 0.0_PREC, Xy, Nvars)

        ! Find alpha grid upper bound for this chunk
        call get_alpha_grid_max (Xy, Ntrain, l1_ratio, tmp)

        alpha_max = max(alpha_max, tmp)

    end do

    if (alpha_max == 0.0_PREC) then
        lstatus = NF_STATUS_INVALID_STATE
        goto 100
    end if

    ! --- Create grid ---

    ! Replicate same grid for tiny alpha_max as in scikit-learn
    resolution = 10.0_PREC ** (-precision(0.0_PREC))
    if (alpha_max < resolution) then
        grid = resolution
        goto 100
    end if

    call logspace (grid, log10(alpha_max*conf%alpha_eps), log10(alpha_max))
    ! Invert sort order such that smallest alpha is last
    Nalpha = size(grid)
    do i = 1, Nalpha / 2
        tmp = grid(i)
        grid(i) = grid(Nalpha-i+1)
        grid(Nalpha-i+1) = tmp
    end do

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine enet_path (conf, X, y, l1_ratio, alphas, coefs, &
        coefs_init, XX, Xy, dual_gaps, status)
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:), contiguous, target :: y
    real (PREC), intent(in) :: l1_ratio
    real (PREC), intent(in), dimension(:), contiguous :: alphas
    real (PREC), intent(out), dimension(:,:), contiguous :: coefs
    logical, intent(in), optional :: coefs_init
        !*  If true, use the values in COEFS present on entry as initial guesses
        !   (default: initialize COEFS with zeros)
    real (PREC), intent(in), dimension(:,:), contiguous, optional, target :: XX
    real (PREC), intent(in), dimension(:), contiguous, optional, target :: Xy
    real (PREC), intent(out), dimension(:), contiguous, optional :: dual_gaps
    type (status_t), intent(out), optional :: status

    logical :: lcoefs_init
    real (PREC) :: l1_reg, l2_reg, alpha, gap
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_XX
    real (PREC), dimension(:), pointer, contiguous :: ptr_Xy
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Xp
    real (PREC), dimension(:), pointer, contiguous :: ptr_yp
    type (status_t) :: lstatus, status_all
    integer :: Nobs, k, Nalphas, i, Nvars
    character (*), parameter :: NAME = 'ENET_PATH'

    lstatus = NF_STATUS_OK

    nullify (ptr_Xp, ptr_yp)
    nullify (ptr_XX, ptr_Xy)

    call get_regressor_dims (X, conf%trans_x, k, Nobs)
    Nalphas = size(alphas)

    ! --- Input checks ---

    call set_optional_arg (coefs_init, .false., lcoefs_init)

    call check_conf (conf, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (all(alphas > 0.0_PREC), NAME, &
        'ALPHAS: Invalid value, ALPHA > 0 required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (l1_ratio >= 0.0_PREC .and. l1_ratio <= 1.0_PREC, NAME, &
        'L1_RATIO: Invalid value, 0 <= L1_RATIO <= 1 required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_dims (X, y, lstatus, conf%trans_x, XX=XX, Xy=Xy)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(coefs,1) == k .and. size(coefs,2) == Nalphas, &
        NAME, 'Non-conformable COEFS argument', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(dual_gaps)) then
        call check_cond (size(dual_gaps) == Nalphas, NAME, &
            'Non-conformable DUAL_GAPS argument', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Pre-compute X'X and X'y ---

    if (.not. present(XX) .or. .not. present(Xy)) then
        if (conf_transform_data (conf)) then
            allocate (ptr_Xp(Nobs,k))

            call transform_regressors (X, ptr_Xp, trans=conf%trans_x, &
                center=conf%center, scale=conf%scale, drop_const=conf%drop_const, &
                tol_const=conf%tol_const, drop_na=.true., nkeep=Nvars, &
                status=lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100

            allocate (ptr_yp(Nobs), source=y)
            call standardize (ptr_yp, center=conf%center, scale=.false.)
        else
            Nvars = k
            ptr_Xp => X
            ptr_yp => y
        end if

        allocate (ptr_XX(Nvars,Nvars), ptr_Xy(Nvars))
        call compute_moments (ptr_Xp(:,1:Nvars), ptr_yp, ptr_XX, ptr_Xy, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        Nvars = k
        ptr_XX => XX
        ptr_Xy => Xy
    end if


    ! --- Estimate coefficients ---

    if (.not. lcoefs_init) then
        ! No initial values given in COEFS, set estimate for first alpha
        ! to zero.
        coefs(:,1) = 0.0
    end if

    status_all = NF_STATUS_OK
    do i = 1, Nalphas
        alpha = alphas(i)
        l1_reg = alpha * l1_ratio * Nobs
        l2_reg = alpha * (1.0_PREC - l1_ratio) * Nobs

        if (.not. lcoefs_init .and. i > 1) then
            ! No initial values given in COEFS, use estimate from previous
            ! iteration
            coefs(:,i) = coefs(:,i-1)
        end if

        call coord_descent_gram (ptr_XX, ptr_Xy, y, l1_reg, l2_reg, coefs(:,i), &
            lstatus, conf%maxiter, conf%tol, conf%positive, gap)
        if (.not. (NF_STATUS_OK .in. lstatus)) goto 100

        status_all = status_all + lstatus

        if (present(dual_gaps)) dual_gaps(i) = gap
    end do

    lstatus = status_all

100 continue

    call assert_dealloc_ptr (X, ptr_Xp)
    call assert_dealloc_ptr (y, ptr_yp)

    call assert_dealloc_ptr (XX, ptr_XX)
    call assert_dealloc_ptr (Xy, ptr_Xy)

    if (present(status)) status = lstatus

end subroutine



subroutine enet_path_multi (conf, X, Y, l1_ratio, alphas, coefs, &
        coefs_init, dual_gaps, status)
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
        !*  Array of NPRED predictors containing NOBS observations,
        !   shaped [NOBS,NPRED]. If CONF%TRANS_X=.TRUE., X stores
        !   the data in transposed form with shape [NPRED,NOBS].
    real (PREC), intent(in), dimension(:,:), contiguous, target :: Y
        !*  Array of NLHS outcomes (responses) containing NOBS observations,
        !   shaped [NOBS,NLHS]. If CONF%TRANS_Y=.TRUE., Y stores
        !   the data in transposed form with shape [NLHS,NOBS].
    real (PREC), intent(in) :: l1_ratio
        !*  Mixture weight on Lasso component.
    real (PREC), intent(in), dimension(:), contiguous :: alphas
        !*  Vector of alphas (regularization parameters) of length NALPHAS
        !   which represent the path along which elastic net should be fitted.
    real (PREC), intent(out), dimension(:,:,:), contiguous :: coefs
        !*  Array storing the fitted coefficients, shaped
        !   [NPRED,NLHS,NALPHAS], or [NLHS,NPRED,NALPHAS] if
        !   CONF%TRANS_COEFS = .TRUE.
    logical, intent(in), optional :: coefs_init
        !*  If true, use the values in COEFS present on entry as initial guesses
        !   (default: initialize COEFS with zeros)
    real (PREC), intent(out), dimension(:), contiguous, optional :: dual_gaps
    type (status_t), intent(out), optional :: status

    logical :: lcoefs_init
    real (PREC) :: l1_reg, l2_reg, alpha, gap
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Xp
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Yp
    type (status_t) :: lstatus, status_all
    real (PREC), dimension(:,:), allocatable :: coefs_T
    integer, dimension(:), allocatable :: ilhs
    integer :: Nobs, k, Nalphas, i, j, iv, Nrhs, Nlhs, Nlhs_in, dim_coefs
    character (*), parameter :: NAME = 'ENET_PATH_MULTI'

    lstatus = NF_STATUS_OK

    nullify (ptr_Xp, ptr_Yp)

    call get_regressor_dims (X, conf%trans_x, k, Nobs)
    call get_regressor_dims (Y, conf%trans_y, nvars=Nlhs)
    Nalphas = size(alphas)

    ! --- Input checks ---

    call set_optional_arg (coefs_init, .false., lcoefs_init)

    call check_conf (conf, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (all(alphas > 0.0_PREC), NAME, &
        'ALPHAS: Invalid value, ALPHA > 0 required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (l1_ratio >= 0.0_PREC .and. l1_ratio <= 1.0_PREC, NAME, &
        'L1_RATIO: Invalid value, 0 <= L1_RATIO <= 1 required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_dims (X, Y, lstatus, trans_x=conf%trans_x, trans_y=conf%trans_y)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Check dimensions of COEFS array
    ! DIM_COEFS is the dimension along which coefficients are aligned
    dim_coefs = 1
    if (conf%trans_coefs) dim_coefs = 2

    call check_cond (size(coefs,dim_coefs) == k .and. &
        size(coefs,3-dim_coefs) == Nlhs .and. size(coefs,3) == Nalphas, &
        NAME, 'Non-conformable COEFS argument', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(dual_gaps)) then
        call check_cond (size(dual_gaps) == Nalphas, NAME, &
            'Non-conformable DUAL_GAPS argument', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Prepare data ---

    if (conf_transform_data (conf, Nlhs)) then
        allocate (ptr_Xp(Nobs,k))

        call transform_regressors (X, ptr_Xp, trans=conf%trans_x, &
            center=conf%center, scale=conf%scale, drop_const=conf%drop_const, &
            tol_const=conf%tol_const, drop_na=.true., nkeep=Nrhs, &
            status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        allocate (ptr_yp(Nobs,Nlhs), source=y)
        allocate (ilhs(Nlhs))

        call transform_regressors (Y, ptr_Yp, trans=conf%trans_y, &
            center=conf%center, scale=.false., drop_const=conf%drop_const, &
            tol_const=conf%tol_const, drop_na=.false., ikeep=ilhs, &
            nkeep=Nlhs_in, status=status)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        Nrhs = k
        Nlhs_in = Nlhs
        allocate (ilhs(Nlhs))
        call arange (ilhs, 1)
        ptr_Xp => X
        ptr_yp => y
    end if

    ! --- Estimate coefficients ---

    allocate (coefs_T(Nlhs_in,Nrhs))

    if (lcoefs_init) then
        do j = 1, Nlhs_in
            iv = ilhs(j)
            if (conf%trans_coefs) then
                coefs_T(j,:) = coefs(iv,:,1)
            else
                coefs_T(j,:) = coefs(:,iv,1)
            end if
        end do
    else
        ! No initial values given in COEFS, use zeros instead.
        ! Any user-provided initial values will be copied over in the loop below.
        coefs_T(:,:) = 0.0
    end if

    status_all = NF_STATUS_OK
    do i = 1, Nalphas
        alpha = alphas(i)
        l1_reg = alpha * l1_ratio * Nobs
        l2_reg = alpha * (1.0_PREC - l1_ratio) * Nobs

        call coord_descent_multi_task (ptr_Xp, ptr_Yp, l1_reg, l2_reg, &
            coefs_T, lstatus, conf%maxiter, conf%tol, gap)
        if (.not. (NF_STATUS_OK .in. lstatus)) goto 100

        if (Nlhs_in == Nlhs) then
            ! Quick code path if all LHS variables are included
            if (conf%trans_coefs) then
                coefs(:,:,i) = coefs_T
            else
                coefs(:,:,i) = transpose (coefs_T)
            end if
        else
            coefs(:,:,i) = 0.0
            do j = 1, Nlhs_in
                iv = ilhs(j)
                if (conf%trans_coefs) then
                    coefs(iv,:,i) = coefs_T(j,:)
                else
                    coefs(:,iv,i) = coefs_T(j,:)
                end if
            end do
        end if

        status_all = status_all + lstatus

        if (present(dual_gaps)) dual_gaps(i) = gap
    end do

    lstatus = status_all

100 continue

    call assert_dealloc_ptr (X, ptr_Xp)
    call assert_dealloc_ptr (Y, ptr_Yp)

    if (present(status)) status = lstatus

end subroutine



subroutine enet_path_mse (conf, X, y, l1_ratio, alphas, test_ifrom, test_n, &
        path_mse, status)
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous :: y
    real (PREC), intent(in) :: l1_ratio
    real (PREC), intent(in), dimension(:), contiguous :: alphas
    integer, intent(in) :: test_ifrom
    integer, intent(in) :: test_n
    real (PREC), intent(out), dimension(:), contiguous :: path_mse
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    type (enet_config) :: lconf
    integer, dimension(:), allocatable :: ivars
    real (PREC), dimension(:,:), allocatable :: X_train, X_test
    real (PREC), dimension(:), allocatable :: y_train, y_test
    real (PREC), dimension(:,:), allocatable :: coefs
    real (PREC), dimension(:,:), allocatable :: XX, resid
    real (PREC), dimension(:), allocatable :: Xy, intercepts
    real (PREC), dimension(:), allocatable :: mean_x, std_x
    real (PREC) :: mean_y, ssr
    integer :: Nalphas, Nobs, k, i, j, Nvars, Ntrain

    lstatus = NF_STATUS_OK

    call get_regressor_dims (X, conf%trans_x, k, Nobs)
    Nalphas = size(alphas)

    ! --- Input checks ---

    call check_conf (conf, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_dims (X, y, lstatus, conf%trans_x)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(alphas) == size(path_mse), 'GLMNET', &
        'Non-conformable arrays ALPHAS, PATH_MSE', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Split sample ---

    Ntrain = Nobs - test_n

    ! Create train and test subsamples in contiguous arrays.
    ! This already applies X^T as needed!
    call extract_block_alloc (X, test_ifrom, test_n, conf%trans_x, x_test, x_train)
    call extract_block_alloc (y, test_ifrom, test_n, y_test, y_train)

    ! --- Preprocess data ---

    ! SPLIT_SAMPLE takes care of transposing X
    lconf = conf
    lconf%trans_x = .false.

    allocate (mean_x(k), std_x(k), ivars(k))

    if (conf_transform_data (lconf)) then
        ! Transform regressors in training set
        call transform_regressors (X_train, center=conf%center, scale=conf%scale, &
            drop_const=conf%drop_const, tol_const=conf%tol_const, &
            drop_na=.true., mean_x=mean_x, std_x=std_x, ikeep=ivars, &
            nkeep=Nvars, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        ! Transform outcome in training set
        call standardize (y_train, center=conf%center, scale=.false., &
            mean_x=mean_y)
    else
        Nvars = k
        call arange (ivars, 1)
        call std_impl (X_train, s=std_x, m=mean_x, dim=1, dof=1)
        call mean_impl (y_train, m=mean_y)
    end if

    ! --- Pre-compute X'X, X'y ---

    allocate (XX(Nvars,Nvars), Xy(Nvars))

    call compute_moments (X_train(:,1:Nvars), y_train, XX, Xy, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Estimate coefs ---

    ! These transformations have all been taken care of above
    lconf%trans_x = .false.
    lconf%scale = .false.
    lconf%center = .false.
    lconf%drop_const = .false.

    allocate (coefs(Nvars, Nalphas))

    ! Estimate coefs along path of alphas
    call enet_path (lconf, X_train(:,1:Nvars), y_train, l1_ratio, alphas, &
        coefs, XX=XX, Xy=Xy, status=lstatus)
    if (.not. (NF_STATUS_OK .in. lstatus)) then
        goto 100
    else
        ! Do not propagate non-convergence status downwards
        lstatus = NF_STATUS_OK
    end if

    ! Rescale coefs to undo standardisation
    if (conf%scale) then
        do i = 1, Nalphas
            do j = 1, Nvars
                coefs(j,i) = coefs(j,i) / std_x(j)
            end do
        end do
    end if

    ! Compute intercepts for each alpha:
    !   const = mean(y) - mean(X) * coefs
    allocate (intercepts(Nalphas), source=mean_y)
    call BLAS_GEMV ('T', Nvars, Nalphas, -1.0_PREC, coefs, Nvars, &
        mean_x(1:Nvars), 1, 1.0_PREC, intercepts, 1)

    ! --- Compute residuals ---

    ! Drop irrelevant variables from X_test; relevant variables
    ! will be copied to the first NVARS columns.
    call replay_transform (X_test, ikeep=ivars(1:Nvars))

    ! resid = y - yhat = y - intercept - Xb
    allocate (resid(test_n, Nalphas))

    ! Compute predicted values -Xb
    call BLAS_GEMM ('N', 'N', test_n, Nalphas, Nvars, -1.0_PREC, X_test(:,1:Nvars), &
        test_n, coefs, Nvars, 0.0_PREC, resid, test_n)

    ! Complete residuals by adding y - intercept
    do i = 1, Nalphas
        resid(:,i) = resid(:,i) + y_test - intercepts(i)
    end do

    ! --- MSEs for each alpha ---

    do i = 1, Nalphas
        ! Sum of squared residuals
        ssr = BLAS_NRM2 (test_n, resid(:,i), 1)
        ! Undo the square root applied by NRM2 and compute average
        path_mse(i) = ssr ** 2.0_PREC / test_n
    end do

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine enet_path_mse_multi (conf, X, Y, l1_ratio, alphas, test_ifrom, &
        test_n, path_mse, status)
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:,:), contiguous :: Y
    real (PREC), intent(in) :: l1_ratio
    real (PREC), intent(in), dimension(:), contiguous :: alphas
    integer, intent(in) :: test_ifrom
    integer, intent(in) :: test_n
    real (PREC), intent(out), dimension(:), contiguous :: path_mse
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    type (enet_config) :: lconf
    integer, dimension(:), allocatable :: irhs, ilhs
    real (PREC), dimension(:,:), allocatable :: X_train, X_test
    real (PREC), dimension(:,:), allocatable :: Y_train, Y_test
    real (PREC), dimension(:,:,:), allocatable :: coefs
    real (PREC), dimension(:,:), allocatable, target :: resid
    real (PREC), dimension(:,:), allocatable :: intercepts
    real (PREC), dimension(:), allocatable :: mean_x, std_x, mean_y
    real (PREC), dimension(:), pointer, contiguous :: ptr_resid
    real (PREC) :: ssr
    integer :: Nalphas, Nobs, Nrhs, Nrhs_in, Ntrain, Nlhs, Nlhs_in
    integer :: i, j

    lstatus = NF_STATUS_OK

    call get_regressor_dims (X, conf%trans_x, Nrhs, Nobs)
    call get_regressor_dims (Y, conf%trans_y, nvars=Nlhs)
    Nalphas = size(alphas)

    ! --- Input checks ---

    call check_conf (conf, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_dims (X, Y, lstatus, trans_x=conf%trans_x, trans_y=conf%trans_y)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(alphas) == size(path_mse), 'GLMNET', &
        'Non-conformable arrays ALPHAS, PATH_MSE', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Split sample ---

    Ntrain = Nobs - test_n

    ! Create train and test subsamples in contiguous arrays.
    ! This already applies X^T as needed!
    call extract_block_alloc (X, test_ifrom, test_n, conf%trans_x, x_test, x_train)
    call extract_block_alloc (Y, test_ifrom, test_n, conf%trans_y, Y_test, Y_train)

    ! --- Preprocess data ---

    ! SPLIT_SAMPLE takes care of transposing X
    lconf = conf
    lconf%trans_x = .false.
    lconf%trans_y = .false.

    allocate (mean_x(Nrhs), std_x(Nrhs), mean_y(Nlhs))
    allocate (irhs(Nrhs), ilhs(Nlhs))

    if (conf_transform_data (lconf, Nlhs)) then
        ! Transform regressors in training set
        call transform_regressors (X_train, center=conf%center, scale=conf%scale, &
            drop_const=conf%drop_const, tol_const=conf%tol_const, &
            drop_na=.true., mean_x=mean_x, std_x=std_x, ikeep=irhs, &
            nkeep=Nrhs_in, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        ! Transform outcomes in training set; transform has already been
        ! taken care of by EXTRACT_BLOCK_ALLOC
        call transform_regressors (Y_train, center=conf%center, scale=.false., &
            drop_const=conf%drop_const, tol_const=conf%tol_const, &
            drop_na=.false., mean_x=mean_y, ikeep=ilhs, nkeep=Nlhs_in, &
            status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        Nrhs_in = Nrhs
        call arange (irhs, 1)
        Nlhs_in = Nlhs
        call arange (ilhs, 1)
        call std_impl (X_train, s=std_x, m=mean_x, dim=1, dof=1)
        call mean_impl (Y_train, m=mean_y, dim=1)
    end if

    ! --- Estimate coefs ---

    ! These transformations have all been taken care of above
    lconf%trans_x = .false.
    lconf%trans_y = .false.
    lconf%scale = .false.
    lconf%center = .false.
    lconf%drop_const = .false.
    ! We allocate a non-transposed COEFS array in this routine
    lconf%trans_coefs = .false.

    allocate (coefs(Nrhs_in, Nlhs_in, Nalphas))

    ! Estimate coefs along path of alphas
    call enet_path (lconf, X_train(:,1:Nrhs_in), Y_train(:,1:Nlhs_in), &
        l1_ratio, alphas, coefs, status=lstatus)
    if (.not. (NF_STATUS_OK .in. lstatus)) then
        goto 100
    else
        ! Do not propagate non-convergence status downwards
        lstatus = NF_STATUS_OK
    end if

    ! Rescale coefs to undo standardisation
    if (conf%scale) then
        do i = 1, Nalphas
            do j = 1, Nlhs_in
                coefs(:,j,i) = coefs(:,j,i) / std_x(1:Nrhs_in)
            end do
        end do
    end if

    ! Compute intercepts for each alpha indexed by i:
    !   const = mean(Y) - mean(X) * coefs
    ! computed as
    !   const = mean(Y) - coefs(:,:,i)' * mean(x)
    allocate (intercepts(Nlhs_in,Nalphas))
    intercepts(:,:) = spread(mean_y(1:Nlhs_in), dim=2, ncopies=Nalphas)

    do i = 1, Nalphas
        call BLAS_GEMV ('T', Nrhs_in, Nlhs_in, -1.0_PREC, coefs(:,:,i), Nrhs_in, &
            mean_x(1:Nrhs_in), 1, 1.0_PREC, intercepts(:,i), 1)
    end do

    ! --- Compute residuals ---

    ! Drop irrelevant variables from X_test; relevant variables
    ! will be copied to the first NVARS columns.
    call replay_transform (X_test, ikeep=irhs(1:Nrhs_in))
    call replay_transform (Y_test, ikeep=ilhs(1:Nlhs_in))

    allocate (resid(test_n, Nlhs_in))

    ! 1d-pointer used to compute MSE using NRM2
    ptr_resid(1:test_n * Nlhs_in) => resid

    do i = 1, Nalphas
        ! Compute predicted values on X_test, -XB, which is shaped [Ntest,Nlhs_in]
        call BLAS_GEMM ('N', 'N', test_n, Nlhs_in, Nrhs_in, -1.0_PREC, X_test, &
            test_n, coefs(:,:,i), Nrhs_in, 0.0_PREC, resid, test_n)

        ! Compute residuals shaped [Ntest,Nlhs_in]
        !   resid = y - yhat = y - intercept - XB
        resid(:,:) = resid(:,:) + Y_test(:,1:Nlhs_in)
        ! Subtract intercept
        do j = 1, Nlhs_in
            resid(:,j) = resid(:,j) - intercepts(j,i)
        end do

        ! Sum of squared residuals
        ssr = BLAS_NRM2 (test_n * Nlhs_in, ptr_resid, 1)
        ! Undo the square root applied by NRM2 and compute average
        path_mse(i) = ssr ** 2.0_PREC / test_n / Nlhs_in
    end do

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine enet_cv (conf, X, y, alpha, l1_ratio, alphas, rmse, status)
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous :: y
    real (PREC), intent(out) :: alpha
    real (PREC), intent(in), optional :: l1_ratio
    real (PREC), intent(in), dimension(:), contiguous, optional :: alphas
    real (PREC), intent(out), optional :: rmse
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    real (PREC) :: ll1_ratio, ubound
    integer :: Nalphas, Nfolds, Ntest, Nobs, k, i, ifrom, imin
    real (PREC), dimension(:), allocatable :: lalphas
    integer, dimension(:), allocatable :: folds_ifrom, folds_size
    real (PREC), dimension(:,:), allocatable :: path_mse
    real (PREC), dimension(:), allocatable :: mean_mse, se_mse, resid
    integer, dimension(:), allocatable :: iwork
    real (PREC) :: min_mse
    integer (NF_ENUM_KIND) :: status_code
    type (enet_config) :: lconf
    character (*), parameter :: NAME = 'ENET_CV'

    lstatus = NF_STATUS_OK

    ! Initialize to something to avoid uninitialized variable compiler warnings
    alpha = huge(0.0_PREC)

    ! --- Input checks ---

    call set_optional_arg (l1_ratio, DEFAULT_L1_RATIO, ll1_ratio)

    call check_conf (conf, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (l1_ratio >= 0.0_PREC .and. l1_ratio <= 1.0_PREC, NAME, &
        'L1_RATIO: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (l1_ratio == 0.0_PREC) then
        ! Require user-provided grid of ALPHAS we cannot determine
        ! grid upper bound for model without Lasso component.
        if (.not. present(alphas)) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    call check_dims (X, y, trans_x=conf%trans_x, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call get_regressor_dims (X, conf%trans_x, k, Nobs)

    ! Ensure that there are not enough obs. for # of CV chunks
    call check_cond (conf%cv_n <= Nobs, NAME, &
        'CONF: Too few obs. for requested number of CV folds', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- k-fold CV chunks ---

    Nfolds = conf%cv_n
    allocate (folds_ifrom(Nfolds), folds_size(Nfolds))

    call split_uniform (Nobs, Nfolds, folds_ifrom, folds_size, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Alpha grid ---

    if (.not. present (alphas)) then
        Nalphas = conf%alpha_n
        allocate (lalphas(Nalphas))
        call create_alpha_grid_cv (conf, X, y, ll1_ratio, folds_ifrom, &
            folds_size, lalphas, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        Nalphas = size(alphas)
        allocate (lalphas(Nalphas))
        ! Make sure that user-provided ALPHAS are in descending order
        allocate (iwork(Nalphas))
        call argsort (alphas, iwork)
        do i = 1, Nalphas
            lalphas(i) = alphas(iwork(Nalphas-i+1))
        end do
        deallocate (iwork)
    end if

    ! --- Perform cross-validation ---

    ! Local copy of config file
    lconf = conf
    ! We always want to standardize during cross-validation. Even if the
    ! *whole* sample is standardized, individual CV'd chunks will not be
    ! in general.
    lconf%center = .true.
    lconf%scale = .true.

    allocate (path_mse(Nalphas, Nfolds))

    ! Workaround: OpenMP in gfortran 7.x cannot use user-defined
    ! operators, keep track using integer instead.
    status_code = NF_STATUS_OK

    !$omp parallel default(none) &
    !$omp shared(Nfolds,folds_ifrom,folds_size,lconf,X,y,ll1_ratio,lalphas) &
    !$omp shared(path_mse) &
    !$omp private(i,ifrom,Ntest,lstatus) &
    !$omp reduction(ior: status_code)

    !$omp do schedule(auto)
    do i = 1, Nfolds
        ifrom = folds_ifrom(i)
        Ntest = folds_size(i)
        call enet_path_mse (lconf, X, y, ll1_ratio, lalphas, ifrom, Ntest, &
            path_mse(:,i), lstatus)
        ! Convert to integer to which IOR reduction can be applied
        status_code = lstatus
    end do
    !$omp end do
    !$omp end parallel

    ! Convert back from integer to STATUS_T
    lstatus = status_code

    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Select optimal alpha ---

    allocate (mean_mse(Nalphas), se_mse(Nalphas))
    allocate (resid(Nfolds))

    ! MSE for each alpha, averaged across cross-validation samples
    mean_mse(:) = sum(path_mse, dim=2) / Nfolds

    do i = 1, Nalphas
        resid(:) =  path_mse(i,:) - mean_mse(i)
        ! SE of mean: 1/sqrt(Nobs) * std(resid) = 1/sqrt(Nobs) * sqrt(sum(resid**2.0)/(Nobs-1))
        se_mse(i) = sqrt(sum(resid**2.0_PREC)) / sqrt(real(Nfolds,PREC)) &
            / sqrt(Nfolds-1.0_PREC)
    end do

    min_mse = huge(0.0_PREC)
    imin = 0

    do i = 1, Nalphas
        if (mean_mse(i) < min_mse) then
            imin = i
            min_mse = mean_mse(i)
        end if
    end do

    ! Pick most parsimonious model within one standard deviation of
    ! the minimum.
    if (conf%cv_parsimonious) then
        ubound = min_mse + se_mse(imin)
        do i = 1, imin - 1
            if (mean_mse(i) <= ubound) then
                imin = i
                exit
            end if
        end do

    end if

    alpha = lalphas(imin)
    if (present(rmse)) rmse = sqrt(min_mse)

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine enet_cv_multi (conf, X, Y, alpha, l1_ratio, alphas, rmse, status)
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:,:), contiguous, target :: Y
    real (PREC), intent(out) :: alpha
    real (PREC), intent(in), optional :: l1_ratio
    real (PREC), intent(in), dimension(:), contiguous, optional :: alphas
    real (PREC), intent(out), optional :: rmse
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    real (PREC) :: ll1_ratio, ubound
    integer :: Nalphas, Nfolds, Ntest, Nobs, Nlhs, k, i, ifrom, imin
    real (PREC), dimension(:), allocatable :: lalphas
    integer, dimension(:), allocatable :: folds_ifrom, folds_size
    real (PREC), dimension(:,:), allocatable :: path_mse
    real (PREC), dimension(:), allocatable :: mean_mse, se_mse, resid
    integer, dimension(:), allocatable :: iwork
    real (PREC), dimension(:), pointer, contiguous :: ptr_Y
    real (PREC) :: min_mse, std_resid
    integer (NF_ENUM_KIND) :: status_code
    type (enet_config) :: lconf
    character (*), parameter :: NAME = 'ENET_CV'

    lstatus = NF_STATUS_OK

    nullify (ptr_Y)

    ! Initialize to something to avoid uninitialized variable compiler warnings
    alpha = huge(0.0_PREC)

    ! --- Input checks ---

    call set_optional_arg (l1_ratio, DEFAULT_L1_RATIO, ll1_ratio)

    call check_conf (conf, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (l1_ratio >= 0.0_PREC .and. l1_ratio <= 1.0_PREC, NAME, &
        'L1_RATIO: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (l1_ratio == 0.0_PREC) then
        ! Require user-provided grid of ALPHAS we cannot determine
        ! grid upper bound for model without Lasso component.
        if (.not. present(alphas)) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    call check_dims (X, Y, trans_x=conf%trans_x, trans_y=conf%trans_y, &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call get_regressor_dims (X, conf%trans_x, k, Nobs)
    call get_regressor_dims (Y, conf%trans_y, nvars=Nlhs)

    ! Ensure that there are not enough obs. for # of CV chunks
    call check_cond (conf%cv_n <= Nobs, NAME, &
        'CONF: Too few obs. for requested number of CV folds', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    Nfolds = conf%cv_n
    call check_cond (Nfolds >= 1 .or. Nfolds <= Nobs, NAME, &
        'CONF: Invalid number of folds', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- k-fold CV chunks ---

    allocate (folds_ifrom(Nfolds), folds_size(Nfolds))

    call split_uniform (Nobs, Nfolds, folds_ifrom, folds_size, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Alpha grid ---

    if (.not. present (alphas)) then
        Nalphas = conf%alpha_n
        allocate (lalphas(Nalphas))
        call create_alpha_grid_cv (conf, X, Y, ll1_ratio, folds_ifrom, &
            folds_size, lalphas, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        Nalphas = size(alphas)
        allocate (lalphas(Nalphas))
        ! Make sure that user-provided ALPHAS are in descending order
        allocate (iwork(Nalphas))
        call argsort (alphas, iwork)
        do i = 1, Nalphas
            lalphas(i) = alphas(iwork(Nalphas-i+1))
        end do
        deallocate (iwork)
    end if

    ! --- Perform cross-validation ---

    ! Local copy of config file
    lconf = conf
    ! We always want to standardize during cross-validation. Even if the
    ! *whole* sample is standardized, individual CV'd chunks will not be
    ! in general.
    lconf%center = .true.
    lconf%scale = .true.

    allocate (path_mse(Nalphas, Nfolds))

    ! Workaround: OpenMP in gfortran 7.x cannot use user-defined
    ! operators, keep track using integer instead.
    status_code = NF_STATUS_OK

    !$omp parallel default(none) &
    !$omp shared(Nfolds,folds_ifrom,folds_size,lconf,X,Y,ll1_ratio,lalphas) &
    !$omp shared(path_mse,Nlhs) &
    !$omp private(i,ifrom,Ntest,lstatus,ptr_Y) &
    !$omp reduction(ior: status_code)

    !$omp do schedule(auto)
    do i = 1, Nfolds
        ifrom = folds_ifrom(i)
        Ntest = folds_size(i)
        if (Nlhs == 1 .and. .not. lconf%force_multi) then
            ! If there is only one LHS variable, skip the more complex
            ! multi-task CV and channel execution.
            ptr_Y(1:size(Y)) => Y
            call enet_path_mse (lconf, X, ptr_Y, ll1_ratio, lalphas, ifrom, &
                Ntest, path_mse(:,i), lstatus)
        else
            call enet_path_mse (lconf, X, Y, ll1_ratio, lalphas, ifrom, &
                Ntest, path_mse(:,i), lstatus)
        end if
        ! Convert to integer to which IOR reduction can be applied
        status_code = lstatus
    end do
    !$omp end do
    !$omp end parallel

    ! Convert back from integer to STATUS_T
    lstatus = status_code

    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Select optimal alpha ---

    allocate (mean_mse(Nalphas), se_mse(Nalphas))
    allocate (resid(Nfolds))

    ! MSE for each alpha, averaged across cross-validation samples
    mean_mse(:) = sum(path_mse, dim=2) / Nfolds

    do i = 1, Nalphas
        resid(:) =  path_mse(i,:) - mean_mse(i)
        ! SE of mean: 1/sqrt(Nobs) * std(resid) = 1/sqrt(Nobs) * sqrt(sum(resid**2.0)/(Nobs-1))
        std_resid = sqrt(sum(resid**2.0)) / sqrt(Nfolds-1.0_PREC)
        se_mse(i) = std_resid / sqrt(Nfolds-1.0_PREC)
    end do

    min_mse = huge(0.0_PREC)
    imin = 0

    do i = 1, Nalphas
        if (mean_mse(i) < min_mse) then
            imin = i
            min_mse = mean_mse(i)
        end if
    end do

    ! Pick most parsimonious model within one standard deviation of
    ! the minimum.
    if (conf%cv_parsimonious) then
        ubound = min_mse + se_mse(imin)
        do i = 1, imin - 1
            if (mean_mse(i) <= ubound) then
                imin = i
                exit
            end if
        end do
    end if

    alpha = lalphas(imin)
    if (present(rmse)) rmse = sqrt(min_mse)

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine enet_fit (conf, X, y, alpha, l1_ratio, intercept, coefs, res, status)
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:), contiguous, target :: y
    real (PREC), intent(in) :: alpha
    real (PREC), intent(in) :: l1_ratio
    real (PREC), intent(out), optional :: intercept
    real (PREC), intent(out), dimension(:), contiguous, optional :: coefs
    type (enet_result), intent(out), optional :: res
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    type (enet_config) :: lconf
    real (PREC) :: ll1_ratio
    integer :: Nobs, k, Nvars
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Xp
    real (PREC), dimension(:), pointer, contiguous :: ptr_yp
    real (PREC), dimension(:,:), allocatable :: XX
    real (PREC), dimension(:), allocatable :: Xy, mean_x, scale_x
    integer, dimension(:), allocatable :: ivars
    real (PREC), dimension(:,:), allocatable :: lcoefs
    real (PREC), dimension(1) :: alphas
    real (PREC) :: mean_y, lintercept
    character (*), parameter :: NAME = 'ENET_FIT'

    lstatus = NF_STATUS_OK

    nullify (ptr_Xp, ptr_Yp)

    call get_regressor_dims (X, conf%trans_x, k, Nobs)

    ! --- Input checks ---

    call set_optional_arg (l1_ratio, DEFAULT_L1_RATIO, ll1_ratio)

    call check_conf (conf, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (alpha > 0.0_PREC, NAME, &
        'ALPHA: Invalid value, ALPHA > 0 required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (l1_ratio >= 0.0_PREC .and. l1_ratio <= 1.0_PREC, NAME, &
        'L1_RATIO: Invalid value, 0 <= L1_RATIO <= 1 required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_dims (X, y, trans_x=conf%trans_x, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(coefs)) then
        call check_cond (size(coefs) >= k, NAME, 'COEFS: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Pre-process input data ---

    allocate (mean_x(k), scale_x(k), ivars(k))

    if (conf_transform_data (conf)) then
        allocate (ptr_Xp(Nobs, k))
        call transform_regressors (X, ptr_Xp, trans=conf%trans_x, &
            center=conf%center, scale=conf%scale, drop_const=conf%drop_const, &
            tol_const=conf%tol_const, drop_na=.true., mean_x=mean_x, &
            scale_x=scale_x, ikeep=ivars, nkeep=Nvars, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        allocate (ptr_yp(Nobs), source=y)
        call standardize (ptr_yp, center=conf%center, scale=.false., mean_x=mean_y)
    else
        ptr_Xp => X
        ptr_yp => y
        Nvars = k
        call arange (ivars, 1)
        call mean_impl (y, m=mean_y)
    end if

    allocate (XX(Nvars, Nvars), Xy(Nvars))

    call compute_moments (ptr_Xp(:,1:Nvars), ptr_yp, XX, Xy, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Unset all transformation-related settings
    lconf = conf
    lconf%trans_x = .false.
    lconf%center = .false.
    lconf%scale = .false.
    lconf%drop_const = .false.

    ! --- Fit elastic net ---

    alphas(1) = alpha
    allocate (lcoefs(Nvars,1))

    ! Transpose and standardization are done, so disable these in ENET_PATH
    call enet_path (lconf, ptr_Xp(:,1:Nvars), ptr_yp, ll1_ratio, alphas, lcoefs, &
        XX=XX, Xy=Xy, status=lstatus)

    if (.not. (NF_STATUS_OK .in. lstatus)) goto 100

    ! Undo scaling
    if (conf%scale) then
        lcoefs(:,1) = lcoefs(:,1) / scale_x(1:Nvars)
    end if

    ! Recover intercept
    lintercept = mean_y - dot_product (mean_x(1:Nvars), lcoefs(:,1))

    if (present(intercept)) intercept = lintercept

    if (present(coefs)) then
        ! Unpack coefficients for non-constant variables
        call unpack_coefs (ivars(1:Nvars), lcoefs(:,1), coefs)
    endif

    if (present(res)) then
        ! Copy result data into ENET_RESULT object
        res = enet_result (conf=conf, nrhs=k, Nobs=Nobs, l1_ratio=l1_ratio, &
            alpha=alpha, intercept=lintercept)
        call copy_alloc (ivars(1:Nvars), res%irhs)
        call cond_alloc (res%coefs, k)
        call unpack_coefs (ivars(1:Nvars), lcoefs(:,1), res%coefs)
    end if

100 continue

    call assert_dealloc_ptr (X, ptr_Xp)
    call assert_dealloc_ptr (y, ptr_yp)

    if (present(status)) status = lstatus

contains

    subroutine unpack_coefs (ivars, src, dst)
        integer, intent(in), dimension(:), contiguous :: ivars
        real (PREC), intent(in), dimension(:), contiguous :: src
        real (PREC), intent(out), dimension(:), contiguous :: dst

        integer :: i, ik
        dst = 0.0
        do i = 1, size(ivars)
            ik = ivars(i)
            dst(ik) = src(i)
        end do
    end subroutine

end subroutine



subroutine enet_fit_multi (conf, X, Y, alpha, l1_ratio, intercept, coefs, &
        coefs_init, res, status)
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:,:), contiguous, target :: Y
    real (PREC), intent(in) :: alpha
    real (PREC), intent(in) :: l1_ratio
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: intercept
    real (PREC), intent(inout), dimension(:,:), contiguous, optional :: coefs
    logical, intent(in), optional :: coefs_init
    type (enet_result), intent(out), optional :: res
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    type (enet_config) :: lconf
    real (PREC) :: ll1_ratio
    integer :: Nobs, Nrhs, Nrhs_in, Nlhs, Nlhs_in, dim_coefs
    integer :: i
    integer :: shp(2)
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Xp
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Yp
    real (PREC), dimension(:), allocatable :: mean_x, scale_x, mean_y
    real (PREC), dimension(:,:,:), allocatable, target :: coefs_T
    real (PREC), dimension(:), allocatable :: lintercept
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_coefs
    integer, dimension(:), allocatable :: irhs, ilhs
    real (PREC), dimension(1) :: alphas
    logical :: lcoefs_init
    character (*), parameter :: NAME = 'ENET_FIT'

    lstatus = NF_STATUS_OK

    nullify (ptr_Xp, ptr_Yp)

    call get_regressor_dims (X, conf%trans_x, Nrhs, Nobs)
    call get_regressor_dims (Y, conf%trans_y, nvars=Nlhs)

    ! --- Input checks ---

    call set_optional_arg (l1_ratio, DEFAULT_L1_RATIO, ll1_ratio)
    call set_optional_arg (coefs_init, .false., lcoefs_init)

    call check_conf (conf, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (alpha > 0.0_PREC, NAME, &
        'ALPHA: Invalid value, ALPHA > 0 required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (l1_ratio >= 0.0_PREC .and. l1_ratio <= 1.0_PREC, NAME, &
        'L1_RATIO: Invalid value, 0 <= L1_RATIO <= 1 required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_dims (X, Y, trans_x=conf%trans_x, trans_y=conf%trans_y, &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Check dimensions of COEFS array
    dim_coefs = 1
    if (conf%trans_coefs) dim_coefs = 2

    if (present(coefs)) then
        call check_cond (size(coefs,dim_coefs) == Nrhs .and. &
            size(coefs,3-dim_coefs) == Nlhs, NAME, &
            'COEFS: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(intercept)) then
        call check_cond (size(intercept) >= Nlhs, NAME, &
            'INTERCEPT: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Pre-process input data ---

    allocate (mean_x(Nrhs), scale_x(Nrhs), mean_y(Nlhs))
    allocate (irhs(Nrhs), ilhs(Nlhs))

    if (conf_transform_data (conf, Nlhs)) then
        allocate (ptr_Xp(Nobs, Nrhs))
        call transform_regressors (X, ptr_Xp, trans=conf%trans_x, &
            center=conf%center, scale=conf%scale, drop_const=conf%drop_const, &
            tol_const=conf%tol_const, drop_na=.true., mean_x=mean_x, &
            scale_x=scale_x, ikeep=irhs, nkeep=Nrhs_in, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        allocate (ptr_Yp(Nobs,Nlhs))
        call transform_regressors (Y, ptr_Yp, trans=conf%trans_y, &
            center=conf%center, scale=.false., drop_const=conf%drop_const, &
            tol_const=conf%tol_const, drop_na=.false., mean_x=mean_y, &
            ikeep=ilhs, nkeep=Nlhs_in, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        ptr_Xp => X
        ptr_Yp => Y
        Nrhs_in = Nrhs
        call arange (irhs, 1)
        Nlhs_in = Nlhs
        call arange (ilhs, 1)
        call mean_impl (Y, m=mean_y, dim=1)
    end if

    ! Unset all transformation-related settings
    lconf = conf
    lconf%trans_x = .false.
    lconf%trans_y = .false.
    lconf%center = .false.
    lconf%scale = .false.
    lconf%drop_const = .false.

    ! --- No LHS variables ---

    ! Gracefully handle degenerate case when there are no (non-constant)
    ! LHS variables

    if (Nlhs_in == 0) then
        ! Allocate zero-size arrays to not break unpacking code below
        allocate (coefs_T(0,Nrhs_in,1))
        allocate (lintercept(0))
        goto 10
    end if

    ! --- Fit elastic net ---

    ! We fit the elastic net on the subset of non-constant LHS and RHS
    ! variables. The coefficients (other than the intercepts) of constant
    ! predictors or constant response variables will be set to zero.
    ! The intercepts of constant response variables will be set to
    ! the corresponding constant values.

    ! We need a contiguous array in case N*hs_in < N*hs, as then user-provided
    ! COEFS(1:Nlhs_in,1:Nrhs_in) will not be contiguous.
    allocate (coefs_T(Nlhs_in,Nrhs_in,1))

    ! --- Populate with initial guesses ---

    if (lcoefs_init) then
        if (Nrhs == Nrhs_in .and. Nlhs == Nlhs_in) then
            if (conf%trans_coefs) then
                coefs_T(:,:,1) = coefs
            else
                coefs_T(:,:,1) = transpose (coefs)
            end if
        else
            if (conf%trans_coefs) then
                call pack_indexed (coefs, ilhs(1:Nlhs_in), irhs(1:Nrhs_in), &
                    coefs_T(:,:,1), trans=.false.)
            else
                call pack_indexed (coefs, irhs(1:Nrhs_in), ilhs(1:Nlhs_in), &
                    coefs_T(:,:,1), trans=.true.)
            end if
        end if
    end if

    ! --- Fit elastic net ---

    alphas(1) = alpha

    if (Nlhs_in == 1 .and. .not. lconf%force_multi) then
        ! Since the last two dimensions of LCOEFS both have size 1, we
        ! here interchange the LHS and ALPHA dimensions, since this
        ! makes no difference.
        ptr_coefs(1:size(coefs_T),1:1) => coefs_T
        lconf%trans_coefs = .false.
        call enet_path (lconf, ptr_Xp(:,1:Nrhs_in), ptr_Yp(:,1), ll1_ratio, &
            alphas, ptr_coefs, coefs_init=coefs_init, status=lstatus)
    else
        lconf%trans_coefs = .true.
        call enet_path (lconf, ptr_Xp(:,1:Nrhs_in), ptr_Yp(:,1:Nlhs_in), &
            ll1_ratio, alphas, coefs_T, coefs_init=coefs_init, status=lstatus)
    end if

    if (.not. (NF_STATUS_OK .in. lstatus)) goto 100

    ! --- Undo scaling ---

    if (conf%scale) then
        do i = 1, Nrhs_in
            coefs_T(:,i,1) = coefs_T(:,i,1) / scale_x(i)
        end do
    end if

    ! --- Recover intercept ---

    ! Only compute intercept if it needs to be returned to caller
    if (present(intercept) .or. present(res)) then
        allocate (lintercept(Nlhs_in))
        ! For multi-outcome Y we have
        !   intercept = mean(Y) - B' mean(X)
        ! where B is shaped [Nrhs_in,Nlhs]
        call BLAS_COPY (Nlhs_in, mean_y, 1, lintercept, 1)
        call BLAS_GEMV ('N', Nlhs_in, Nrhs_in, -1.0_PREC, coefs_T(:,:,1), Nlhs_in, &
            mean_x(1:Nrhs_in), 1, 1.0_PREC, lintercept, 1)
    end if

    ! Jump to this point of no LHS variables are included in the model
10  continue

    ! --- Unpack intercept ---

    if (present(intercept)) then
        call unpack_intercept (lintercept, intercept)
    end if

    ! --- Unpack coefficients ---

    if (present(coefs)) then
        ! Unpack coefficients for non-constant LHS and RHS variables
        if (conf%trans_coefs) then
            call unpack_indexed (coefs_T(:,:,1), ilhs(1:Nlhs_in), irhs(1:Nrhs_in), &
                coefs, fill=0.0_PREC, trans=.false.)
        else
            call unpack_indexed (coefs_T(:,:,1), irhs(1:Nrhs_in), ilhs(1:Nlhs_in), &
                coefs, fill=0.0_PREC, trans=.true.)
        end if
    endif

    if (present(res)) then
        ! Copy result data into ENET_RESULT object
        res = enet_result (conf=conf, nrhs=Nrhs, nlhs=Nlhs, Nobs=Nobs, &
            l1_ratio=l1_ratio, alpha=alpha)
        call copy_alloc (irhs(1:Nrhs_in), res%irhs)

        call cond_alloc (res%intercept_multi, Nlhs)
        call unpack_intercept (lintercept, res%intercept_multi)

        shp(1) = Nrhs
        shp(2) = Nlhs
        call cond_alloc (res%coefs_multi, shp)
        call unpack_indexed (coefs_T(:,:,1), irhs(1:Nrhs_in), ilhs(1:Nlhs_in), &
            res%coefs_multi, fill=0.0_PREC, trans=.true.)
    end if

100 continue

    call assert_dealloc_ptr (X, ptr_Xp)
    call assert_dealloc_ptr (Y, ptr_Yp)

    if (present(status)) status = lstatus

contains

    subroutine unpack_intercept (src, dst)
        !*  UNPACK_INTERCEPT expands the intercepts computed on the set
        !   of non-constant response variables to the set of all response
        !   variables. This is done by filling in the blacks with the
        !   corresponding constant values.
        real (PREC), intent(in), dimension(:), contiguous :: src
        real (PREC), intent(out), dimension(:), contiguous :: dst

        dst = 0.0

        ! Populate with first obs. of Y. For constant responses, this
        ! will also be their mean and hence their intercept.
        if (conf%trans_y) then
            dst(1:Nlhs) = Y(1:Nlhs,1)
        else
            dst(1:Nlhs) = Y(1,1:Nlhs)
        end if

        call unpack_indexed (src(1:Nlhs_in), ilhs(1:Nlhs_in), dst(1:Nlhs))
    end subroutine

end subroutine



subroutine enet_predict (res, X, y_pred, status)
    type (enet_result), intent(in), target :: res
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(out), dimension(:), contiguous :: y_pred
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer, contiguous :: ptr_X
    real (PREC), dimension(:), pointer, contiguous :: ptr_coefs
    integer :: Nobs, k, Nvars, i, iv
    logical :: drop_vars
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'ENET_PREDICT'

    lstatus = NF_STATUS_OK

    ! --- Input checks ---

    call get_regressor_dims (X, res%conf%trans_x, k, Nobs)
    Nvars = size(res%irhs)

    call check_result (res, X=X, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(y_pred) == Nobs, NAME, 'Y_PRED: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Determine whether any variables in X need to be dropped from the model
    drop_vars = .not. has_all_vars (X, res%irhs, trans=res%conf%trans_x)
    if (drop_vars .or. res%conf%trans_x) then
        allocate (ptr_X(Nobs,Nvars))
        call replay_transform (X, ptr_X, add_intercept=.false., &
            trans=res%conf%trans_x, ikeep=res%irhs, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        ! Eliminate irrelevant coefficients.
        ! Note that we don't want to apply the unpacked coef array to the
        ! unpacked regressor matrix X, as X could contain Infs and NaNs
        ! which were dropped when fitting. The corresponding variable would
        ! then have received a zero coefficient, but NaN * 0 or Inf * 0
        ! will corrupt the prediction.
        allocate (ptr_coefs(Nvars))
        do i = 1, Nvars
            iv = res%irhs(i)
            ptr_coefs(iv) = res%coefs(i)
        end do
    else
        ! No transformations need to be applied, use input data directly
        ptr_X => X
        ptr_coefs => res%coefs
    end if

    ! Compute predicted values
    call BLAS_GEMV ('N', Nobs, Nvars, 1.0_PREC, ptr_X, Nobs, ptr_coefs, 1, &
        0.0_PREC, y_pred, 1)

    ! Add intercept
    y_pred = y_pred + res%intercept

100 continue

    call assert_dealloc_ptr (X, ptr_X)
    call assert_dealloc_ptr (res%coefs, ptr_coefs)

    if (present(status)) status = lstatus

end subroutine



subroutine enet_predict_multi (res, X, Y_pred, trans_x, trans_y, status)
    type (enet_result), intent(in), target :: res
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(out), dimension(:,:), contiguous :: Y_pred
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
    type (status_t), intent(out), optional :: status

    integer :: Nobs, k, Nrhs
    logical :: ltrans_x, ltrans_y
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    ! --- Input checks ---

    call set_optional_arg (trans_x, res%conf%trans_x, ltrans_x)
    call set_optional_arg (trans_y, res%conf%trans_y, ltrans_y)

    call get_regressor_dims (X, res%conf%trans_x, k, Nobs)
    Nrhs = size(res%irhs)

    call check_result (res, X=X, Y_multi=Y_pred, trans_x=ltrans_x, &
        trans_y=ltrans_y, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (size(res%irhs) == k) then
        call predict (X, res%coefs_multi, Y_pred, res%intercept_multi, &
            trans_x=ltrans_x, trans_y=ltrans_y, trans_coefs=.false., &
            status=lstatus)
    else
        call predict (X, res%coefs_multi, Y_pred, res%intercept_multi, &
            irhs=res%irhs, trans_x=ltrans_x, trans_y=ltrans_y, &
            trans_coefs=.false., status=lstatus)
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine enet_post_estim_impl (res, X, y, rsq, status)
    type (enet_result), intent(in) :: res
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous :: y
    real (PREC), intent(out), optional :: rsq
    type (status_t), intent(out), optional :: status

    integer :: k, Nobs, Nvars
    real (PREC) :: std_y, var_y, ssr
    type (status_t) :: lstatus
    real (PREC), dimension(:), allocatable :: y_pred, resid

    lstatus = NF_STATUS_OK

    ! --- Input checks ---

    call check_dims (X, y, trans_x=res%conf%trans_x, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call get_regressor_dims (X, res%conf%trans_x, k, Nobs)
    Nvars = size(res%irhs)

    call check_result (res, X=X, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Compute R-squared ---

    if (present(rsq)) then
        call compute_predicted (lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        call std (y, s=std_y)
        var_y = std_y ** 2.0

        ! Compute residuals
        call compute_resid ()

        ! Compute SSR / N = 1/N sum(resid_i**2.0)
        ssr = BLAS_NRM2 (Nobs, resid, 1) ** 2.0 / Nobs

        rsq = 1.0_PREC - ssr/var_y
    end if

100 continue

    if (present(status)) status = lstatus

contains

    subroutine compute_predicted (status)
        type (status_t), intent(out) :: status
        if (.not. allocated (y_pred)) then
            allocate (y_pred(Nobs))
            call enet_predict (res, X, y_pred, status=status)
        end if
    end subroutine

    subroutine compute_resid ()
        if (.not. allocated (resid)) then
            allocate (resid(Nobs))
            resid(:) = y - y_pred
        end if
    end subroutine
end subroutine



subroutine coord_descent_gram (xx, xy, y, alpha, beta, w, status, &
        maxiter, tol, positive, gap, niter, gtol)
    !*  COORD_DESCENT_GRAM implements the coordinate descent algorithm
    !   for the elastic net regression which uses the Gram matrix
    !   G = X'X as an input.
    !
    !   The routine computes the solution to
    !       min_w 1/2 w'X'Xw - y'Xw + alpha ||w||_1 + 1/2 beta ||w||_2
    !   where ||.||_1 and ||.||_2 are the l1 and l2 norms, prespectively.
    !
    !   The code is a direct port of the Cython implementation from
    !   scikit-learn.
    real (PREC), intent(in), dimension(:,:), contiguous :: xx
        !*  Array containing the matrix X'X, shaped [Npred,Npred]
        !   where Npred is the number of predictors in the model.
    real (PREC), intent(in), dimension(:), contiguous :: xy
        !*  Array containing the vector X'y of length Npred.
    real (PREC), intent(in), dimension(:), contiguous :: y
        !*  Array containing the outcome variable y of length Nobs.
    real (PREC), intent(in) :: alpha
        !*  L1 regularization parameter.
    real (PREC), intent(in) :: beta
        !*  L2 regularization parameter
    real (PREC), intent(out), dimension(:), contiguous :: w
        !*  Solution vector of length Npred.
    type (status_t), intent(out) :: status
        !*  Exit code
    integer, intent(in), optional :: maxiter
        !*  Optional maximum number of iterations.
    real (PREC), intent(in), optional :: tol
        !*  Optional termination tolerance
    logical, intent(in), optional :: positive
    real (PREC), intent(out), optional :: gap
        !*  If present, contains the duality gap on successful termination.
    integer, intent(out), optional :: niter
        !*  If present, contains the actual number of iterations on exit.
    real (PREC), intent(out), optional :: gtol
        !*  If present, contains the actual termination tolerance for the
        !   duality gap.

    logical :: lpositive
    integer :: lmaxiter, iter, ik
    real (PREC) :: ltol, lgap, dw_tol, w_max, dw_max, y_norm2, w_ik, dw_ik, tmp, sgn
    real (PREC) :: qw, dual_norm_XtA, R_norm2, w_norm2, w_norm1, A_norm2, const
    integer :: Nobs, k
    real (PREC), dimension(:), allocatable :: H, XtA

    status = NF_STATUS_OK
    iter = 0

    call set_optional_arg (tol, DEFAULT_COORD_DESCENT_TOL, ltol)
    call set_optional_arg (maxiter, DEFAULT_COORD_DESCENT_MAXITER, lmaxiter)
    call set_optional_arg (positive, .false., lpositive)

    ! Number of samples
    Nobs = size(y)
    ! Number of features (regressors)
    k = size(xx, 1)

    ! Compute y'y as |y|_2^2
    y_norm2 = BLAS_NRM2 (Nobs, y, 1) ** 2.0_PREC

    lgap = ltol + 1.0
    dw_tol = ltol
    ltol = ltol * y_norm2

    allocate (H(k))
    allocate (XtA(k), source=0.0_PREC)

    ! Compute H = X'X * w
    call BLAS_GEMV ('N', k, k, 1.0_PREC, xx, k, w, 1, 0.0_PREC, H, 1)

    do iter = 1, lmaxiter
        w_max = 0.0
        dw_max = 0.0

        ! --- Loop over coordinates ---
        do ik = 1, k
            if (xx(ik, ik) == 0.0_PREC) cycle

            ! Store previous vablue
            w_ik = w(ik)

            if (w_ik /= 0.0_PREC) then
                ! Compute H = H - w_ik * X'X(ik,:) = H - w_ik * X'X(:,ik)
                ! since X'X is symmetric.
                call BLAS_AXPY (k, -w_ik, xx(:,ik), 1, H, 1)
            end if

            tmp = xy(ik) - H(ik)
            
            if (lpositive .and. tmp < 0.0_PREC) then
                w(ik) = 0.0
            else
                sgn = signum (tmp)
                w(ik) = sgn * max(abs(tmp)-alpha, 0.0_PREC) / (xx(ik,ik) + beta)
            end if

            if (w(ik) /= 0.0_PREC) then
                ! Update H = X'X w
                ! H = H + w(ik) * X'X(ik,:) = H + w(ik) * X'X(:,ik)
                call BLAS_AXPY (k, w(ik), xx(:,ik), 1, H, 1)
            end if

            ! Perform maximum absolute coefficient update
            dw_ik = abs(w(ik) - w_ik)
            if (dw_ik > dw_max) then
                dw_max = dw_ik
            end if

            if (abs(w_ik) > w_max) then 
                w_max = abs(w_ik)
            end if
        end do

        if ((w_max == 0.0_PREC) .or. (dw_max < dw_tol*w_max) &
                .or. (iter == lmaxiter)) then
            
            qw = dot_product (w, xy)

            do ik = 1, k
                XtA(ik) = xy(ik) - H(ik) - beta * w(ik)
            end do

            if (lpositive) then
                dual_norm_XtA = maxval(XtA)
            else
                ik = BLAS_IAMAX (k, XtA, 1)
                dual_norm_XtA = abs(XtA(ik))
            end if

            ! tmp = sum(w * H)
            tmp = sum(w * H)
            R_norm2 = y_norm2 + tmp - 2.0_PREC * qw

            w_norm2 = BLAS_NRM2 (k, w, 1) ** 2.0_PREC

            if (dual_norm_XtA > alpha) then
                const = alpha / dual_norm_XtA
                A_norm2 = R_norm2 * const**2.0_PREC
                lgap = 0.5_PREC * (R_norm2 + A_norm2)
            else
                const = 1.0
                lgap = R_norm2
            end if

            ! Call to ASUM is equivalent to L1 norm
            w_norm1 = BLAS_ASUM (k, w, 1)
            lgap = lgap + (alpha * w_norm1 - const * y_norm2 + const * qw & 
                + 0.5_PREC * beta * (1.0_PREC + const**2.0_PREC) * w_norm2)

            ! Return if the desired tolerance was achieved
            if (lgap < ltol) exit
        end if

    end do

    if (iter > lmaxiter) then
        status = NF_STATUS_OK
        status = status + NF_STATUS_NOT_CONVERGED
    end if

    if (present(gap)) gap = lgap
    if (present(niter)) niter = min(iter, lmaxiter)
    if (present(gtol)) gtol = ltol

end subroutine



subroutine coord_descent_multi_task (X, Y, alpha, beta, W, status, maxiter, &
        tol, gap, niter, gtol)
    !*  COORD_DESCENT_MULTI_TASK implements the coordinate descent algorithm
    !   for multile outcome variables ("task"), use for elastic net
    !   multi-task regression.
    !
    !   The routine computes the minimum
    !       min_W 1/2 ||Y-XW'||_F^2 + alpha ||W'||_21 + 1/2 * beta * ||W||_F^2
    !   where ||.||_F is the Frobenius norm, and ||.||_21 is the l1-norm
    !   of the column-wise l2-norms of a given matrix.
    !
    !   The code is a direct port of the Cython implementation from
    !   scikit-learn.
    real (PREC), intent(in), dimension(:,:), contiguous :: X
        !*  Array of predictors of shape [Nobs,Npred], ie. each column of X
        !   contains one predictor variable.
    real (PREC), intent(in), dimension(:,:), contiguous, target :: Y
        !*  Array of outcomies of shape [Nobs,Ntask], ie. each column of Y
        !   contains one outcome ("task").
    real (PREC), intent(in) :: alpha
        !*  L1 regularization parameter
    real (PREC), intent(in) :: beta
        !*  L2 regularization parameter
    real (PREC), intent(out), dimension(:,:), contiguous, target :: W
        !*  Coefficient array of shape [Ntask,Npred]
    type (status_t), intent(out) :: status
        !*  Exit code.
    integer, intent(in), optional :: maxiter
        !*  Optional maximum number of iteration
    real (PREC), intent(in), optional :: tol
        !*  Optional termination tolerance.
    real (PREC), intent(out), optional :: gap
        !*  If present, contains the duality gap on successful exit.
    integer, intent(out), optional :: niter
        !*  If present, contains the actual number of iterations on exit.
    real (PREC), intent(out), optional :: gtol
        !*  If present, contains the actual termination tolerance for the
        !   duality gap.

    integer :: lmaxiter, Nobs, Npred, Ntask
    real (PREC) :: ltol, w_ij
    real (PREC), dimension(:,:), allocatable :: XtA
    real (PREC), dimension(:,:), allocatable, target :: Resid
    real (PREC), dimension(:), allocatable :: norm_cols_X, tmp, w_i
    real (PREC) :: XtA_dim1_norm, dual_norm_XtA
    real (PREC) :: dw_max, w_max, dw_i, nn, w_i_abs_max, lgap, dw_tol, const
    real (PREC) :: R_norm, W_norm, A_norm, ry_sum, l21_norm, Y_norm_F, factr
    real (PREC), dimension(:), pointer, contiguous :: ptr_Y, ptr_Resid, ptr_W
    integer :: i, j, iter

    status = NF_STATUS_OK
    iter = 0

    ! --- Input processing ---

    call set_optional_arg (tol, DEFAULT_COORD_DESCENT_TOL, ltol)
    call set_optional_arg (maxiter, DEFAULT_COORD_DESCENT_MAXITER, lmaxiter)

    ! Number of observations (samples)
    Nobs = size(X, 1)
    ! Number of predictors ("features", independent regressors)
    Npred = size(X, 2)
    ! Number of outcome variables
    Ntask = size(Y, 2)

    ! Note: We flip the dimensions of XtA compared to sklearn implementation
    allocate (XtA(Ntask, Npred))

    allocate (norm_cols_X(Npred))
    allocate (tmp(Ntask))
    allocate (w_i(Ntask))

    do i = 1, Npred
        norm_cols_X(i) = BLAS_NRM2 (Nobs, X(:,i), 1) ** 2.0_PREC
    end do

    ! Initial value of residuals Y - XW'
    allocate (Resid(Nobs, Ntask), source=Y)
    call BLAS_GEMM ('N', 'T', Nobs, Ntask, Npred, -1.0_PREC, X, Nobs, W, Ntask, &
        1.0_PREC, Resid, Nobs)

    lgap = ltol + 1.0
    dw_tol = ltol

    ! tol = tol * ||Y||_F ** 2
    ! We can compute the Frobenius norm of Y simply as the l2 norm
    ! of a flattened vector Y
    ptr_Y(1:Ntask*Nobs) => Y
    Y_norm_F = BLAS_NRM2 (Ntask*Nobs, ptr_Y, 1) ** 2.0_PREC
    ltol = ltol * Y_norm_F

    ! 1d pointers to arrays so we can compute the Frobenius norm directly
    ! as L2 norm on vectors
    ptr_W(1:Ntask*Npred) => W
    ptr_Resid(1:Ntask*Nobs) => Resid

    do iter = 1, lmaxiter
        w_max = 0.0
        dw_max = 0.0

        do i = 1, Npred

            if (norm_cols_X(i) == 0.0_PREC) cycle

            ! Store previous value
            w_i(:) = W(:,i)

            ! Perform rank-1 update of residual
            !  R = R + X(:,i) * W(:,i)'
            ! The corresponding call to BLAS GER is
            !   call BLAS_GER (Nobs, Ntask, 1.0_PREC, X(:,i), 1, w_i, 1, Resid, Nobs)
            ! We do this using AXPY instead of GER to avoid multiplying
            ! columns of X with zeros.
            do j = 1, Ntask
                w_ij = w_i(j)
                if (w_ij /= 0.0_PREC) then
                    call BLAS_AXPY (Nobs, w_ij, X(:,i), 1, Resid(:,j), 1)
                end if
            end do

            ! Compute tmp' = X(:,i)' R
            ! We instead compute the transpose tmp = R'X(:,i)
            call BLAS_GEMV ('T', Nobs, Ntask, 1.0_PREC, Resid, Nobs, X(:,i), &
                1, 0.0_PREC, tmp, 1)

            nn = BLAS_NRM2 (Ntask, tmp, 1)

            ! W(:,i) = tmp * max(1-alpha/nn, 0) / (norm_cols_X(i) + beta)
            factr = max(1.0_PREC-alpha/nn, 0.0_PREC) / (norm_cols_X(i) + beta)
            W(:,i) = tmp * factr

            ! If factor = 0, there is no point in executing the rank-1
            ! update as W(:,i) will be all zeros.
            if (factr > 0.0_PREC) then
                ! Perform rank-1 update of residual
                !   R = R - X(:,i) * W(:,i)'
                ! The corresponding call to BLAS GER is
                !   call BLAS_GER (Nobs, Ntask, -1.0_PREC, X(:,i), 1, W(:,i), 1, Resid, Nobs)
                ! We do this using AXPY instead of GER to avoid multiplying
                ! columns of X with zeros.
                do j = 1, Ntask
                    w_ij = W(j,i)
                    if (w_ij /= 0.0_PREC) then
                        call BLAS_AXPY (Nobs, -w_ij, X(:,i), 1, Resid(:,j), 1)
                    end if
                end do
            end if

            ! Update maximum absolute coefficient
            dw_i = maxval(abs(W(:,i) - w_i))

            if (dw_i > dw_max) then
                dw_max = dw_i
            end if

            w_i_abs_max = maxval(abs(W(:,i)))
            if (w_i_abs_max > w_max) then
                w_max = w_i_abs_max
            end if
        end do

        if ((w_max == 0.0_PREC) .or. (dw_max < dw_tol*w_max) &
                .or. (iter == lmaxiter)) then
            ! The biggest coordinate update of this iteration was smaller than
            ! the tolerance: check the duality lgap as ultimate stopping criterion

            ! XtA = R'X - beta * W
            ! sklearn code computes the transpose, XtA' = X'R - beta * W'
            XtA(:,:) = W
            call BLAS_GEMM ('T', 'N', Ntask, Npred, Nobs, 1.0_PREC, Resid, Nobs, &
                X, Nobs, -beta, XtA, Ntask)

            ! Compute dual_norm_XtA = maxval(sqrt(sum(XtA**2, dim=1)))
            ! It's actually possible to use the above statement directly,
            ! but we prefer to use NRM2 on each column.
            dual_norm_XtA = 0.0
            do i = 1, Npred
                XtA_dim1_norm = BLAS_NRM2 (Ntask, XtA(:,i), 1)
                if (XtA_dim1_norm > dual_norm_XtA) then
                    dual_norm_XtA = XtA_dim1_norm
                end if
            end do

            R_norm = BLAS_NRM2 (Nobs * Ntask, ptr_Resid, 1)
            W_norm = BLAS_NRM2 (Ntask * Npred, ptr_W, 1)

            if (dual_norm_XtA > alpha) then
                const = alpha / dual_norm_XtA
                A_norm = R_norm * const
                lgap = 0.5 * (R_norm**2.0_PREC + A_norm**2.0_PREC)
            else
                const = 1.0
                lgap = R_norm**2.0
            end if

            ! ry_sum = sum(R * Y)
            ry_sum = BLAS_DOT (Nobs * Ntask, ptr_Resid, 1, ptr_Y, 1)

            ! l21_norm = sum(sqrt(sum(W**2, dim=1))
            l21_norm = 0.0
            do i = 1, Npred
                l21_norm = l21_norm + BLAS_NRM2 (Ntask, W(:,i), 1)
            end do

            lgap = lgap + alpha * l21_norm - const * ry_sum &
                + 0.5_PREC * beta * (1.0 + const**2.0_PREC) * W_norm**2.0_PREC

            if (lgap < ltol) exit
        end if

    end do

    if (iter > lmaxiter) then
        status = NF_STATUS_OK
        status = status + NF_STATUS_NOT_CONVERGED
    end if

    if (present(gap)) gap = lgap
    if (present(niter)) niter = min(iter, lmaxiter)
    if (present(gtol)) gtol = ltol

end subroutine



subroutine ridge_multi (conf, X, Y, alpha, coefs, rcond, intercept, rank, res, &
        status)
    !*  RIDGE_2D computes the ridge regression for given independent
    !   data X and (potentially multiple) dependent variables Y.
    !
    !   Note that a regression constant is added by default.
    !
    !   The problem is solved using SVD as implemented in LAPACK's GESDD
    !   routine. The optional argument RCOND can be used to control the
    !   effective rank of the matrix X'X.
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
        !*  Array of RHS (predictor) variables
    real (PREC), intent(in), dimension(:,:), contiguous, target :: Y
        !*  Array of LHS (response) variables (separate regression is performed
        !   for each LHS variables using the same set of RHS variables)
    real (PREC), intent(in) :: alpha
        !*  Ridge parameter controlling the penalty for large coefficients.
    real (PREC), intent(out), dimension(:,:), contiguous, optional, target :: coefs
        !*  Array of estimated coefficients.
    real (PREC), intent(in), optional :: rcond
        !*  Optional argument to control the effective rank of the regressor
        !   matrix.
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: intercept
    integer, intent(out), optional :: rank
        !*  Contains effective rank of regressor matrix
    type (enet_result), optional :: res
        !*  Optional result objects for linear models. Note that a separate
        !   object is returned for each LHS variable.
    type (status_t), intent(out), optional :: status
        !*  Optional exit code

    real (PREC) :: lrcond, var_rhs
    integer :: Nobs, ncoefs, nconst, i, Nrhs, Nrhs_in, Nlhs
    type (status_t) :: lstatus
    real (PREC), dimension(:,:), allocatable :: Xp
    real (PREC), dimension(:,:), allocatable :: lcoefs
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Yp
    real (PREC), dimension(:), allocatable :: mean_x, scale_x, mean_y
    real (PREC), dimension(:), allocatable :: lintercept
    integer, dimension(:), allocatable :: irhs
    integer :: shp(2)
    character (*), parameter :: NAME = 'RIDGE'

    ! GESDD arguments
    real (PREC), dimension(:), allocatable :: sval, work
    real (PREC), dimension(:,:), allocatable :: vt
    real (PREC), dimension(1) :: qwork
    integer, dimension(:), allocatable :: iwork
    character (1), parameter :: jobz = 'O'
    real (PREC), dimension(0,0) :: u
    integer :: lrank
    integer :: lwork, info, m, n, lda, ldvt, ldu, mn

    ! GEMM arguments
    real (PREC), dimension(:,:), allocatable :: mat_Uty

    lstatus = NF_STATUS_OK

    nullify (ptr_Yp)

    ! --- Input processing ---

    call set_optional_arg (rcond, 0.0_PREC, lrcond)

    call check_conf (conf, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (alpha >= 0.0, NAME, 'ALPHA: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (lrcond >= 0.0, NAME, 'RCOND: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_dims (X, Y, trans_x=conf%trans_x, trans_y=conf%trans_y, &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call get_regressor_dims (X, conf%trans_x, Nrhs, Nobs)
    call get_regressor_dims (Y, conf%trans_y, nvars=Nlhs)

    Nconst = 0
    if (conf%add_intercept) Nconst = 1
    Ncoefs = Nconst + Nrhs

    if (present(coefs)) then
        call check_cond (size(coefs,1) == Ncoefs .and. size(coefs,2) >= Nlhs, &
            NAME, 'COEFS: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(intercept)) then
        call check_cond (size(intercept) >= Nlhs, NAME, &
            'INTERCEPT: Array too small', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Prepare input data ---

    allocate (irhs(Nrhs))
    allocate (mean_x(Nrhs), scale_x(Nrhs), mean_y(Nlhs))

    if (conf_transform_data (conf, Nlhs)) then
        allocate (Xp(Nobs,Nrhs))
        call transform_regressors (X, Xp, center=conf%center, scale=conf%scale, &
            trans=conf%trans_x, drop_const=conf%drop_const, &
            tol_const=conf%tol_const, drop_na=.true., mean_x=mean_x, &
            scale_x=scale_x, ikeep=irhs, nkeep=Nrhs_in, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        allocate (ptr_Yp(Nobs,Nlhs))
        call transform_regressors (Y, ptr_Yp, center=conf%center, scale=.false., &
            trans=conf%trans_y, drop_const=.false., drop_na=.false., &
            mean_x=mean_y, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        Nrhs_in = Nrhs
        call arange (irhs, 1)
        ! X will be overwritten by GESDD in place, so we need to create a
        ! copy in any case
        allocate (Xp, source=X)
        ptr_Yp => Y

        call mean_impl (Y, dim=1, m=mean_y)
        call mean_impl (X, dim=1, m=mean_x)
    end if

    ! --- Check for degenerate problem ---

    if (Nrhs_in == 0 .or. Nlhs == 0) then
        ! We can skip the entire SVD part and immediately move to
        ! computing the intercepts
        allocate (lcoefs(Nrhs_in, Nlhs))
        goto 50
    end if

    ! --- SVD decomposition via GESDD ---

    ! Set up GESDD call
    m = Nobs
    n = Nrhs_in
    mn = max(1, min(m, n))
    lda = Nobs
    ldvt = n
    ldu = m

    allocate (iwork(8*mn))
    allocate (vt(n, n))
    allocate (sval(mn))

    ! workspace query
    lwork = -1
    call LAPACK_GESDD (jobz, m, n, Xp, lda, sval, u, ldu, vt, ldvt, qwork, &
        lwork, iwork, info)

    ! Recover minimal work space size
    if (info /= 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if
    lwork = int(qwork(1))

    ! perform actual SVD
    allocate (work(lwork))

    call LAPACK_GESDD (jobz, m, n, Xp, lda, sval, u, ldu, vt, ldvt, work, &
        lwork, iwork, info)

    if (info /= 0) then
        if (info < 0) then
            ! This should not happen since we checked arguments before
            lstatus = NF_STATUS_INVALID_ARG
        else
            lstatus = NF_STATUS_NOT_CONVERGED
        end if
        goto 100
    end if

    ! --- Determine effective rank of X ---

    ! We have rank <= ncoefs
    lrank = count(sval > lrcond)

    ! --- Compute U'y ---

    ! At this point we have: first NRHS_IN columns of U stored in matrix X

    allocate (mat_Uty(lrank, Nlhs))

    call BLAS_GEMM (transa='T', transb='N', m=lrank, n=Nlhs, k=Nobs, &
        alpha=1.0_PREC, a=Xp, lda=Nobs, b=ptr_Yp, ldb=Nobs, beta=0.0_PREC, &
        c=mat_Uty, ldc=lrank)

    ! --- Rescale by inverse of diagonal matrix ---

    ! Compute (diag(s^2) + alpha I)^{-1} diag(s)
    ! and then (diag(s^2) + alpha I)^{-1} diag(s) U'y

    do i = 1, lrank
        sval(i) = sval(i) / (sval(i)**2.0_PREC + alpha)
    end do

    ! Apply scaling factors by columns
    do i = 1, Nlhs
        mat_Uty(:,i) = mat_Uty(:,i) * sval(1:lrank)
    end do

    ! --- Last step: premulitply with (V')^(-1) ---

    ! Due to properties of V we have
    ! V'V = I => V' = V^(-1) => (V')^(-1) = V
    ! Note that GESDD returns V' which is stored in the matrix VT

    allocate (lcoefs(Nrhs_in, Nlhs))

    call BLAS_GEMM (transa='T', transb='N', m=Nrhs_in, n=Nlhs, k=lrank, &
        alpha=1.0_PREC, a=vt, lda=ncoefs, b=mat_Uty, ldb=lrank, beta=0.0_PREC, &
        c=lcoefs, ldc=size(lcoefs,1))

    ! --- Undo rescaling ---

    ! 1. Rescale regression coefficients
    if (conf%scale) then
        do i = 1, Nlhs
            lcoefs(:,i) = lcoefs(:,i) / scale_x(1:Nrhs_in)
        end do
    end if

50 continue

    ! --- Recover intercept ---

    ! Recover intercept only if there is a way to return it to the caller
    if (conf%add_intercept .or. present(intercept) .or. present(res)) then
        allocate (lintercept(Nlhs), source=mean_y(1:Nlhs))

        ! Intercept is given by
        !   intercept = mean(Y) - coefs' * mean(X)
        call BLAS_GEMV ('T', Nrhs_in, Nlhs, -1.0_PREC, lcoefs, Nrhs_in, &
            mean_x(1:Nrhs_in), 1, 1.0_PREC, lintercept, 1)
    end if

    ! --- Process return args ---

    if (present(coefs)) then
        if (conf%add_intercept) then
            coefs(1,:) = lintercept
        end if
        ! Unpack remaining coefficients, fill gaps with zeros.
        call unpack_indexed (lcoefs, irhs(1:Nrhs_in), 1, coefs(Nconst+1:,:), &
            fill=0.0_PREC)
    end if

    if (present(intercept)) then
        intercept = lintercept
    end if

    if (present(res)) then
        ! Update ENET_DATA result objects for OLS model

        ! Fraction of RHS variance used, analogous to PCR regression
        ! (only applicable if RHS matrix does not have full rank)
        var_rhs = 0.0
        if (allocated(sval)) then
            var_rhs = sum(sval(1:lrank) ** 2.0_PREC) / sum(sval(1:Nrhs_in)**2.0_PREC)
        end if

        res = enet_result (conf=conf, nrhs=Nrhs, nlhs=Nlhs, Nobs=Nobs, &
            l1_ratio=0.0_PREC, alpha=alpha, var_rhs=var_rhs)

        call copy_alloc (irhs(1:Nrhs_in), res%irhs)

        call copy_alloc (lintercept, res%intercept_multi)

        shp(1) = Nrhs
        shp(2) = Nlhs
        call cond_alloc (res%coefs_multi, shp)
        call unpack_indexed (lcoefs, irhs(1:Nrhs_in), 1, res%coefs_multi, &
            fill=0.0_PREC)
    end if

    ! Copy over optional output arguments
    if (present(rank)) rank = lrank

100 continue

    call assert_dealloc_ptr (Y, ptr_Yp)

    if (present(status)) status = lstatus

end subroutine



subroutine ridge_single (conf, X, y, alpha, coefs, rcond, intercept, rank, &
        res, status)
    !*  RIDGE_SINGLE computes the ridge regression for given independent
    !   data X and a single response variable Y.
    !
    !   The problem is solved using SVD as implemented in LAPACK's GESDD
    !   routine. The optional argument RCOND can be used to control the
    !   effective rank of the matrix X'X.
    type (enet_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
        !*  Array of RHS variables
    real (PREC), intent(in), dimension(:), contiguous, target :: y
        !*  Array of LHS variables (separate regression is performed for each
        !   LHS variables using the same set of RHS variables)
    real (PREC), intent(in) :: alpha
        !*  Ridge parameter controlling the penalty for large coefficients.
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: coefs
        !*  Array of estimated coefficients.
    real (PREC), intent(in), optional :: rcond
        !*  Optional argument to control the effective rank of the regressor
        !   matrix.
    real (PREC), intent(out), optional :: intercept
    integer, intent(out), optional :: rank
        !*  Contains effective rank of regressor matrix
    type (enet_result), optional :: res
        !*  Optional result objects for linear models.
    type (status_t), intent(out), optional :: status
        !*  Optional exit code

    real (PREC), dimension(:,:), contiguous, pointer :: ptr_coefs, ptr_y
    real (PREC), dimension(:), contiguous, pointer :: ptr_intercept

    integer :: nobs
    type (status_t) :: lstatus
    type (enet_config) :: lconf

    lstatus = NF_STATUS_OK

    nullify (ptr_coefs, ptr_y, ptr_intercept)

    lconf = conf
    lconf%trans_y = .false.

    nobs = size(y)

    ptr_y(1:nobs,1:1) => y

    if (present(coefs)) then
        ptr_coefs(1:size(coefs),1:1) => coefs
    end if

    if (present(intercept)) then
        allocate (ptr_intercept(1))
    end if

    call ridge (lconf, X, ptr_y, alpha, ptr_coefs, rcond, ptr_intercept, &
        rank, res, lstatus)

    if (lstatus == NF_STATUS_OK) then
        if (present(res)) then
            if (allocated(res%coefs_multi)) then
                call copy_alloc (res%coefs_multi(:,1), res%coefs)
            end if

            res%intercept = 0.0
            if (allocated(res%intercept_multi)) then
                res%intercept = res%intercept_multi(1)
            end if
        end if

        if (present(intercept)) then
            intercept = ptr_intercept(1)
        end if
    end if

    if (associated (ptr_intercept)) deallocate (ptr_intercept)

    if (present(status)) status = lstatus

end subroutine



