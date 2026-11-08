

pure subroutine transform_regr (X, xout, add_intercept, trans, &
        center, scale, drop_const, tol_const, drop_na, skip_const, mean_x, std_x, &
        shift_x, scale_x, ikeep, nkeep, idrop, ndrop, weights, status)
    !*  TRANSFORM_REGRESSORS processes the regressor matrix to obtain
    !   the data expected by estimation routines.
    real (PREC), intent(in), dimension(:,:), contiguous :: X
        !*  User-provided regressor matrix
    real (PREC), intent(out), dimension(:,:), contiguous :: xout
        !*  Processed regressor matrix
    logical, intent(in), optional :: add_intercept
        !*  Add intercept to regressor matrix
    logical, intent(in), optional :: trans
        !*  Transpose regressor matrix
    logical, intent(in), optional :: center
        !*  Center variables in X.
    logical, intent(in), optional :: scale
        !*  Scale variables in X to have unit variance.
    logical, intent(in), optional :: drop_const
        !*  Drop any constant variables in X.
    real (PREC), intent(in), optional :: tol_const
        !*  Tolerance threshold determining whether a variable is considered
        !   a constant whenever Var(x) < TOL_CONST.
    logical, intent(in), optional :: drop_na
        !*  Drop any variables with non-finite (NaN, +/- infinity) values.
    logical, intent(in), optional :: skip_const
        !*  Skip rescaling of constant "variables" since this would create
        !   NaNs or Infs. Defaults to true if DROP_CONST=.FALSE.
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: mean_x
        !*  Means of variables in X (excluding dropped variables)
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: std_x
        !*  Std. dev. of variables in X (excluding dropped variables)
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: shift_x
        !*  Values by which variables in X where shifted.
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: scale_x
        !*  Values by which variables in X where scaled.
    integer, intent(out), dimension(:), optional, contiguous, target :: ikeep
        !*  Indices of variables in X which were NOT dropped.
    integer, intent(out), optional :: nkeep
        !*  Number of variables from X which were NOT dropped.
    integer, intent(out), dimension(:), optional :: idrop
        !*  Optional array of variable indices that were dropped from the model.
    integer, intent(out), optional :: ndrop
        !*  Number of variables that were dropped from the model.
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
        !*  Finite, nonnegative moment weights with a positive total.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: nconst, nobs, k, ncoefs, i, ik, dim, lnkeep
    logical :: ladd_intercept, ltrans, ldrop_const, lcenter, lscale, ldrop_na
    logical :: need_std, need_mean_only, lskip_const
    real (PREC) :: ltol_const
    real (PREC), dimension(:), allocatable :: std_x_all, mean_x_all, w_norm
    real (PREC), dimension(:), contiguous, pointer :: ptr_std_x, ptr_mean_x
    real (PREC), dimension(:), allocatable :: lstd_x
    logical, dimension(:), allocatable :: mask_vars
    integer, dimension(:), pointer, contiguous :: ptr_ikeep
    character (*), parameter :: NAME = 'TRANSFORM_REGRESSORS'
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    nullify (ptr_ikeep, ptr_mean_x, ptr_std_x)

    ! --- Input checks ---

    call set_optional_arg (trans, .false., ltrans)
    call set_optional_arg (center, .false., lcenter)
    call set_optional_arg (scale, .false., lscale)
    call set_optional_arg (add_intercept, .false., ladd_intercept)
    call set_optional_arg (drop_const, .false., ldrop_const)
    call set_optional_arg (drop_na, .true., ldrop_na)
    call set_optional_arg (tol_const, 1.0e-14_PREC, ltol_const)
    call set_optional_arg (skip_const, .not. ldrop_const, lskip_const)

    ! --- Determine required moments ---

    call get_regressor_dims (X, ltrans, k, nobs)

    if (present(weights)) then
        call check_cond (size(weights) == nobs, NAME, &
            "WEIGHTS: Non-conformable array", lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
        ! Check finiteness before comparisons or weighted moment calculations.
        call check_cond (all(ieee_is_finite(weights)), NAME, &
            "WEIGHTS: Values must be finite", lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
        call check_cond (all(weights >= 0.0_PREC), NAME, &
            "WEIGHTS: Values must be nonnegative", lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
        call check_cond (sum(weights) > 0.0_PREC, NAME, &
            "WEIGHTS: Sum must be positive", lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! Determine which, if any, moments are needed, and precompute those
    need_std = ldrop_const .or. lscale .or. present(std_x)
    need_mean_only = (lcenter .or. present(mean_x)) .and. .not. need_std

    dim = 1
    if (ltrans) dim = 2

    if (need_std) then
        ! Compute std. dev., get mean for free.
        allocate (std_x_all(k), mean_x_all(k))
        if (present(weights)) then
            allocate (w_norm(nobs))
            w_norm(:) = weights / sum(weights)
            if (dim == 1) then
                do i = 1, k
                    mean_x_all(i) = sum(w_norm * X(:, i))
                    std_x_all(i) = &
                        sqrt(sum(w_norm * (X(:, i) - mean_x_all(i))**2.0_PREC))
                end do
            else
                do i = 1, k
                    mean_x_all(i) = sum(w_norm * X(i, :))
                    std_x_all(i) = &
                        sqrt(sum(w_norm * (X(i, :) - mean_x_all(i))**2.0_PREC))
                end do
            end if
        else
            call std_impl (X, s=std_x_all, m=mean_x_all, dim=dim, dof=1)
        end if
    else if (need_mean_only) then
        allocate (mean_x_all(k))
        if (present(weights)) then
            allocate (w_norm(nobs))
            w_norm(:) = weights / sum(weights)
            if (dim == 1) then
                do i = 1, k
                    mean_x_all(i) = sum(w_norm * X(:, i))
                end do
            else
                do i = 1, k
                    mean_x_all(i) = sum(w_norm * X(i, :))
                end do
            end if
        else
            call mean_impl (X, m=mean_x_all, dim=dim)
        end if
    end if

    ! --- Filter variables ---

    call assert_alloc_ptr (ikeep, k, ptr_ikeep)

    call filter_vars (X, ptr_ikeep, lnkeep, trans, ldrop_const, ltol_const, &
        ldrop_na, std_x_all, idrop, ndrop, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Set mask for variables to be included in model to .TRUE.
    allocate (mask_vars(k), source=.false.)
    do i = 1, lnkeep
        mask_vars(ptr_ikeep(i)) = .true.
    end do

    ! --- Intercept ---

    nconst = 0
    if (ladd_intercept) nconst = 1

    ncoefs = lnkeep + nconst

    ! --- Check remaining args ---

    ! For these checks we need to know how many variables in X are included
    ! in the model

    call check_cond (size(xout, 1) == nobs .and. size(xout, 2) >= ncoefs, &
        NAME, 'XOUT: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(mean_x)) then
        call check_cond (size(mean_x) >= lnkeep, NAME, 'MEAN_X: Array too small', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(std_x)) then
        call check_cond (size(std_x) >= lnkeep, NAME, 'STD_X: Array too small', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(shift_x)) then
        call check_cond (size(shift_x) >= lnkeep, NAME, 'SHIFT_X: Array too small', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(scale_x)) then
        call check_cond (size(scale_x) >= lnkeep, NAME, 'SCALE_X: Array too small', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Copy & transpose data ---

    if (ltrans) then
        do i = 1, lnkeep
            ik = ptr_ikeep(i)
            xout(:,nconst+i) = X(ik,:)
        end do
    else
        do i = 1, lnkeep
            ik = ptr_ikeep(i)
            xout(:, nconst+i) = X(:,ik)
        end do
    end if

    if (ladd_intercept) xout(:, 1) = 1.0_PREC

    ! --- Standardize ---

    if (lscale .or. present(std_x)) then
        call assert_alloc_ptr (std_x, lnkeep, ptr_std_x)
        ptr_std_x(:) = 0.0
        ptr_std_x(1:lnkeep) = pack (std_x_all, mask_vars)
    end if

    if (lcenter .or. present(mean_x)) then
        call assert_alloc_ptr (mean_x, lnkeep, ptr_mean_x)
        ptr_mean_x(:) = 0.0
        ptr_mean_x(1:lnkeep) = pack (mean_x_all, mask_vars)
    end if

    if (lskip_const .and. lscale) then
        ! We need to present NaNs or Infs that would arise from division by
        ! zero, but we want to return the correct std. devs. to user
        ! if STD_X is present, so we create a second array containing
        ! adjusted STD_X with 1.0 for constant variables.
        allocate (lstd_x(lnkeep), source=ptr_std_x(1:lnkeep))
        where (lstd_x == 0.0_PREC)
            lstd_x = 1.0_PREC
        end where

        ! Use helper routine to center and scale X, if applicable
        call standardize_impl (xout(:,1+nconst:ncoefs), ptr_mean_x, lstd_x, &
            dim=1, center=lcenter, scale=lscale, shift_x=shift_x, scale_x=scale_x)
    else
        ! Use helper routine to center and scale X, if applicable
        call standardize_impl (xout(:,1+nconst:ncoefs), ptr_mean_x, ptr_std_x, &
            dim=1, center=lcenter, scale=lscale, shift_x=shift_x, scale_x=scale_x)
    end if


    if (present(nkeep)) nkeep = lnkeep

100 continue

    call assert_dealloc_ptr (std_x, ptr_std_x)
    call assert_dealloc_ptr (mean_x, ptr_mean_x)
    call assert_dealloc_ptr (ikeep, ptr_ikeep)

    if (present(status)) status = lstatus

end subroutine



pure subroutine transform_regr_in_place (X, center, scale, drop_const, &
        tol_const, drop_na, skip_const, mean_x, std_x, shift_x, scale_x, ikeep, &
        nkeep, idrop, ndrop, weights, status)
    !*  TRANSFORM_REGRESSORS processes the regressor matrix to obtain
    !   the data expected by estimation routines.
    real (PREC), intent(inout), dimension(:,:), contiguous :: X
        !*  User-provided regressor matrix
    logical, intent(in), optional :: center
        !*  Center variables in X.
    logical, intent(in), optional :: scale
        !*  Scale variables in X to have unit variance.
    logical, intent(in), optional :: drop_const
        !*  Drop any constant variables in X.
    real (PREC), intent(in), optional :: tol_const
        !*  Tolerance threshold determining whether a variable is considered
        !   a constant whenever Var(x) < TOL_CONST.
    logical, intent(in), optional :: drop_na
        !*  Drop any variables with non-finite (NaN, +/- infinity) values.
    logical, intent(in), optional :: skip_const
        !*  Skip rescaling of constant "variables" since this would create
        !   NaNs or Infs. Defaults to true if DROP_CONST=.FALSE.
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: mean_x
        !*  Means of variables in X (excluding dropped variables)
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: std_x
        !*  Std. dev. of variables in X (excluding dropped variables)
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: shift_x
        !*  Values by which variables in X where shifted.
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: scale_x
        !*  Values by which variables in X where scaled.
    integer, intent(out), dimension(:), optional, contiguous, target :: ikeep
        !*  Indices of variables in X which were NOT dropped.
    integer, intent(out), optional :: nkeep
        !*  Number of variables from X which were NOT dropped.
    integer, intent(out), dimension(:), optional :: idrop
        !*  Optional array of variable indices that were dropped from the model.
    integer, intent(out), optional :: ndrop
        !*  Number of variables that were dropped from the model.
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
        !*  Finite, nonnegative moment weights with a positive total.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: Nobs, k, i, ik, lnkeep
    logical :: ldrop_const, lcenter, lscale, ldrop_na, lskip_const
    real (PREC) :: ltol_const
    real (PREC), dimension(:), allocatable :: std_x_all, mean_x_all, w_norm
    real (PREC), dimension(:), pointer, contiguous :: ptr_std_x, ptr_mean_x
    logical, dimension(:), allocatable :: mask_vars
    integer, dimension(:), pointer, contiguous :: ptr_ikeep
    character (*), parameter :: NAME = 'TRANSFORM_REGRESSORS'
    type (status_t) :: lstatus
    logical :: need_std, need_mean_only
    logical, parameter :: trans = .false.
    real (PREC), dimension(:), allocatable :: lstd_x

    lstatus = NF_STATUS_OK

    nullify (ptr_std_x, ptr_mean_x, ptr_ikeep)

    ! --- Input checks ---

    call set_optional_arg (center, .false., lcenter)
    call set_optional_arg (scale, .false., lscale)
    call set_optional_arg (drop_const, .false., ldrop_const)
    call set_optional_arg (drop_na, .true., ldrop_na)
    call set_optional_arg (tol_const, 1.0e-14_PREC, ltol_const)
    call set_optional_arg (skip_const, .not. ldrop_const, lskip_const)

    ! --- Determine required moments ---

    Nobs = size(X, 1)
    k = size(X, 2)

    if (present(weights)) then
        call check_cond (size(weights) == Nobs, NAME, &
            "WEIGHTS: Non-conformable array", lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
        ! Check finiteness before comparisons or weighted moment calculations.
        call check_cond (all(ieee_is_finite(weights)), NAME, &
            "WEIGHTS: Values must be finite", lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
        call check_cond (all(weights >= 0.0_PREC), NAME, &
            "WEIGHTS: Values must be nonnegative", lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
        call check_cond (sum(weights) > 0.0_PREC, NAME, &
            "WEIGHTS: Sum must be positive", lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! Determine which, if any, moments are needed, and precompute those
    need_std = ldrop_const .or. lscale .or. present(std_x)
    need_mean_only = (lcenter .or. present(mean_x)) .and. .not. need_std

    if (need_std) then
        ! Compute std. dev., get mean for free.
        allocate (std_x_all(k), mean_x_all(k))
        if (present(weights)) then
            allocate (w_norm(Nobs))
            w_norm(:) = weights / sum(weights)
            do i = 1, k
                mean_x_all(i) = sum(w_norm * X(:, i))
                std_x_all(i) = sqrt(sum(w_norm * (X(:, i) - mean_x_all(i))**2.0_PREC))
            end do
        else
            call std_impl (X, s=std_x_all, m=mean_x_all, dim=1, dof=1)
        end if
    else if (need_mean_only) then
        allocate (mean_x_all(k))
        if (present(weights)) then
            allocate (w_norm(Nobs))
            w_norm(:) = weights / sum(weights)
            do i = 1, k
                mean_x_all(i) = sum(w_norm * X(:, i))
            end do
        else
            call mean_impl (X, m=mean_x_all, dim=1)
        end if
    end if

    ! --- Filter variables ---

    call assert_alloc_ptr (ikeep, k, ptr_ikeep)

    call filter_vars (X, ptr_ikeep, lnkeep, trans, ldrop_const, ltol_const, &
        ldrop_na, std_x_all, idrop, ndrop, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Set mask for variables to be included in model to .TRUE.
    allocate (mask_vars(k), source=.false.)
    do i = 1, lnkeep
        mask_vars(ptr_ikeep(i)) = .true.
    end do

    ! --- Check remaining args ---

    ! For these checks we need to know how many variables in X are included
    ! in model

    if (present(mean_x)) then
        call check_cond (size(mean_x) >= lnkeep, NAME, &
            'MEAN_X: Array size too small', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(std_x)) then
        call check_cond (size(std_x) >= lnkeep, NAME, &
            'STD_X: Array size too small', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(shift_x)) then
        call check_cond (size(shift_x) >= lnkeep, NAME, &
            'SHIFT_X: Array size too small', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(scale_x)) then
        call check_cond (size(scale_x) >= lnkeep, NAME, &
            'SCALE_X: Array size too small', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Move data to contiguous block ---

    ! Skip this step if no variables need to be dropped
    if (k /= lnkeep) then
        do i = 1, lnkeep
            ik = ptr_ikeep(i)
            if (i /= ik) then
                call copy_1d (X(:,ik), X(:,i))
            end if
        end do
    end if

    ! --- Standardize ---

    if (lscale .or. present(std_x)) then
        call assert_alloc_ptr (std_x, lnkeep, ptr_std_x)
        ptr_std_x(:) = 0.0
        ptr_std_x(1:lnkeep) = pack (std_x_all, mask_vars)
    end if

    if (lcenter .or. present(mean_x)) then
        call assert_alloc_ptr (mean_x, lnkeep, ptr_mean_x)
        ptr_mean_x(:) = 0.0
        ptr_mean_x(1:lnkeep) = pack (mean_x_all, mask_vars)
    end if

    if (lskip_const .and. lscale) then
        ! We need to present NaNs or Infs that would arise from division by
        ! zero, but we want to return the correct std. devs. to user
        ! if STD_X is present, so we create a second array containing
        ! adjusted STD_X with 1.0 for constant variables.
        allocate (lstd_x(lnkeep), source=ptr_std_x(1:lnkeep))
        where (lstd_x == 0.0_PREC)
            lstd_x = 1.0_PREC
        end where
        ! Use helper routine to center and scale X, if applicable
        call standardize_impl (X(:,1:lnkeep), ptr_mean_x, lstd_x, dim=1, &
            center=lcenter, scale=lscale, shift_x=shift_x, scale_x=scale_x)
    else
        ! Use helper routine to center and scale X, if applicable
        call standardize_impl (X(:,1:lnkeep), ptr_mean_x, ptr_std_x, dim=1, &
            center=lcenter, scale=lscale, shift_x=shift_x, scale_x=scale_x)
    end if

    if (present(nkeep)) nkeep = lnkeep

100 continue

    call assert_dealloc_ptr (std_x, ptr_std_x)
    call assert_dealloc_ptr (mean_x, ptr_mean_x)
    call assert_dealloc_ptr (ikeep, ptr_ikeep)

    if (present(status)) status = lstatus

end subroutine




pure subroutine filter_vars (X, ikeep, nkeep, trans, drop_const, tol_const, &
        drop_na, std_x, idrop, ndrop, status)
    !*  FILTER_VARS compiles lists of variable indices that should be
    !   included and/or excluded in the regressor matrix.
    real (PREC), intent(in), dimension(:,:), contiguous :: X
        !*  User-provided regressor matrix
    integer, intent(out), dimension(:), contiguous :: ikeep
        !*  Indices of variables in X which were NOT dropped.
    integer, intent(out) :: nkeep
        !*  Number of variables from X which were NOT dropped.
    logical, intent(in), optional :: trans
        !*  If presenet and true, transpose regressor matrix
    logical, intent(in), optional :: drop_const
        !*  Drop any constant variables in X.
    real (PREC), intent(in), optional :: tol_const
        !*  Tolerance threshold determining whether a variable is considered
        !   a constant whenever Var(x) < TOL_CONST.
    logical, intent(in), optional :: drop_na
        !*  Drop any variables with non-finite (NaN, +/- infinity) values.
    real (PREC), intent(in), dimension(:), contiguous, optional :: std_x
        !*  Std. dev. of all variables in X
    integer, intent(out), dimension(:), optional :: idrop
        !*  Optional list of variable indices to be dropped.
    integer, intent(out), optional :: ndrop
        !*  Optional number of variables to be dropped.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: Nobs, k, i, ik, n, dim
    logical :: ltrans, ldrop_const, ldrop_na
    logical :: is_finite
    real (PREC) :: ltol_const
    logical, dimension(:), allocatable :: mask_vars
    character (*), parameter :: NAME = 'FILTER_VARS'
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    ! --- Input checks ---

    call set_optional_arg (trans, .false., ltrans)
    call set_optional_arg (drop_const, .false., ldrop_const)
    call set_optional_arg (drop_na, .true., ldrop_na)
    call set_optional_arg (tol_const, 1.0e-14_PREC, ltol_const)

    ! --- Quick exit ---

    call get_regressor_dims (X, ltrans, k, Nobs)

    if (.not. (ldrop_const .or. ldrop_na)) then
        nkeep = k
        call check_cond (size(ikeep) >= nkeep, NAME, 'IKEEP: Array too small', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        ikeep = 0
        call arange (ikeep(1:nkeep), 1)

        goto 50
    end if

    ! --- Implementation ---

    dim = 1
    if (ltrans) dim = 2

    ! Default: all variables are included
    allocate (mask_vars(k), source=.true.)

    if (ldrop_na) then
        do ik = 1, k
            if (ltrans) then
                is_finite = all(ieee_is_finite (X(ik, :)))
            else
                is_finite = all(ieee_is_finite (X(:,k)))
            end if

            mask_vars(ik) = is_finite
        end do
    end if

    if (ldrop_const) then
        mask_vars(:) = mask_vars .and. (std_x(1:k) >= ltol_const)
    end if

    ! Number of admissible variables in X
    nkeep = count (mask_vars)

    ! Check IVARS array size: here we require only that it is sufficiently
    ! large to hold indices of all non-constant variables, but not
    ! all variables in X.
    call check_cond (size(ikeep) >= nkeep, NAME, 'IKEEP: Array too small', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ikeep = 0

    if (nkeep < k) then
        nkeep = 0
        do ik = 1, k
            if (mask_vars (ik)) then
                nkeep = nkeep + 1
                ikeep(nkeep) = ik
            end if
        end do
    else
        call arange (ikeep(1:k), 1)
    end if

50 continue

    ! Process optional variables to be dropped

    if (present(ndrop)) ndrop = k - nkeep
    if (present(idrop)) then
        ! Compile indices of variables which were dropped from model
        n = 0
        do i = 1, k
            if (.not. mask_vars(i)) then
                n = n + 1
                ! Check for enough space, and truncate list of dropped indices
                ! if there is not enough.
                if (size(idrop) > n) exit
                idrop(n) = i
            end if
        end do
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine replay_transform (X, Xout, add_intercept, trans, shift_x, &
        scale_x, ikeep, status)
    !*  REPLAY_TRANSFORM (re)applies a transformation to array X, defined
    !   by the indices of variables to keep, the amount by which variables
    !   in X should be shifted and scaled, etc.
    real (PREC), intent(in), dimension(:,:), contiguous :: X
        !*  Array to which transformyyation should be applied
    real (PREC), intent(out), dimension(:,:), contiguous :: Xout
        !*  Output array generating by transforming X
    logical, intent(in), optional :: add_intercept
        !*  If true, insert an additional column of ones in the first position.
        !   (default: .FALSE.)
    logical, intent(in), optional :: trans
        !*  Assume that data in X is provided in transposed form, ie.
        !   variables are aligned along the first dimension (default: .FALSE.)
    real (PREC), intent(in), dimension(:), contiguous, optional :: shift_x
        !*  Array containing the amounts by which each variable in X
        !   should be shifted (default: do not shift variables)
    real (PREC), intent(in), dimension(:), contiguous, optional :: scale_x
        !*  Array containing the amounts by which each variable in X
        !   should be rescaled (default: do not scale variables)
    integer, intent(in), dimension(:), contiguous, target, optional :: ikeep
        !*  Array of indices of variables in X that should be included
        !   in output array (default: include all variables).
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    logical :: ltrans, ladd_intercept
    integer :: Nvars, Ncoefs, Nconst, Nobs
    integer :: imin, imax, i, ik, k
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'REPLAY_TRANSFORM'

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (trans, .false., ltrans)
    call set_optional_arg (add_intercept, .false., ladd_intercept)

    call get_regressor_dims (X, ltrans, k, Nobs)

    Nconst = 0
    Ncoefs = k
    Nvars = k

    if (ladd_intercept) then
        Nconst = 1
        Ncoefs = k + 1
    end if

    if (present(ikeep)) then
        ! Check that we have only valid indices in IVARS
        imin = minval(ikeep)
        imax = maxval(ikeep)
        call check_cond (imin >= 1 .and. imin <= k, NAME, &
            'IVARS: Invalid value encountered', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        Nvars = size(ikeep)
        Ncoefs = Nvars + Nconst
    end if

    ! Check output array dimensions
    call check_cond (size(Xout,1) == Nobs .and. size(Xout,2) >= Ncoefs, NAME, &
        'XOUT: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(shift_x)) then
        call check_cond (size(shift_x) == Nvars, NAME, &
            'SHIFT_X: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(scale_x)) then
        call check_cond (size(scale_x) == Nvars, NAME, &
            'SCALE_X: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Copy & transpose data ---

    if (ltrans) then
        do i = 1, Nvars
            if (present(ikeep)) then
                ik = ikeep(i)
            else
                ik = i
            end if
            Xout(:,Nconst+i) = X(ik,:)
        end do
    else
        do i = 1, Nvars
            if (present(ikeep)) then
                ik = ikeep(i)
            else
                ik = i
            end if
            Xout(:, Nconst+i) = X(:,ik)
        end do
    end if

    ! Standardize, if applicable.
    ! Pass CENTER=.TRUE. and SCALE=.TRUE., but these operations will
    ! only be performed if SHIFT_X and SCALE_X are present.
    call standardize_impl (Xout, dim=1, mean_x=shift_x, std_x=scale_x, &
        center=.true., scale=.true.)

    if (ladd_intercept) Xout(:,1) = 1.0_PREC


100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine replay_transform_in_place (X, shift_x, scale_x, ikeep, status)
    !*  REPLAY_TRANSFORM_IN_PLACE (re)applies a transformation to array X, defined
    !   by the indices of variables to keep, and the amount by which variables
    !   X should be shifted and scaled.
    real (PREC), intent(inout), dimension(:,:), contiguous :: X
        !*  Array to which transformyyation should be applied
    real (PREC), intent(in), dimension(:), contiguous, optional :: shift_x
        !*  Array containing the amounts by which each variable in X
        !   should be shifted (default: do not shift variables)
    real (PREC), intent(in), dimension(:), contiguous, optional :: scale_x
        !*  Array containing the amounts by which each variable in X
        !   should be rescaled (default: do not scale variables)
    integer, intent(in), dimension(:), contiguous, target, optional :: ikeep
        !*  Array of indices of variables in X that should be included
        !   in output array (default: include all variables).
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: nkeep, Nobs
    integer :: imin, imax, i, ik, k
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'REPLAY_TRANSFORM'

    lstatus = NF_STATUS_OK

    Nobs = size(X, 1)
    k = size(X, 2)

    nkeep = k

    if (present(ikeep)) then
        ! Check that we have only valid indices in IVARS
        imin = minval(ikeep)
        imax = maxval(ikeep)
        call check_cond (imin >= 1 .and. imin <= k, NAME, &
            'IVARS: Invalid value encountered', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        nkeep = size(ikeep)
    end if

    if (present(shift_x)) then
        call check_cond (size(shift_x) == nkeep, NAME, &
            'SHIFT_X: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(scale_x)) then
        call check_cond (size(scale_x) == nkeep, NAME, &
            'SCALE_X: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Copy variables to contiguous block ---

    if (present(ikeep)) then
        do i = 1, nkeep
            ik = ikeep(i)
            if (i /= ik) then
                call copy_1d (X(:,ik), X(:,i))
            end if
        end do
    end if

    ! --- Standardize ---

    ! Pass CENTER=.TRUE. and SCALE=.TRUE., since these operations will
    ! only be performed if SHIFT_X and SCALE_X are present.
    call standardize_impl (X, dim=1, mean_x=shift_x, std_x=scale_x, &
        center=.true., scale=.true.)

100 continue

    if (present(status)) status = lstatus

end subroutine



function has_all_vars (X, ivars, trans) result (res)
    !*  HAS_ALL_VARS returns .TRUE. if the indices in IVARS include all the
    !   variables in X, and .FALSE. otherwise.
    real (PREC), intent(in), dimension(:,:) :: X
        !*  Array containing multiple variables.
    integer, intent(in), dimension(:), contiguous, optional :: ivars
        !*  Array of variable indices.
    logical, intent(in), optional :: trans
        !*  If present and true, assume that X is provided in transposed form
        !   (variables in first dimension).
    logical :: res

    integer :: Nuniq, imin, imax, Nobs, Nvars

    if (.not. present(ivars)) then
        ! If no indices are given, assume that all variables are included
        ! in the model
        res = .true.
        return
    end if

    call get_regressor_dims (X, trans, Nvars, Nobs)

    res = (size(ivars) == Nvars)
    if (res) then
        ! Check that all indices are present in IVARS
        call unique (ivars, Nuniq)
        imin = minval(ivars)
        imax = maxval(ivars)
        res = (Nuniq == Nvars .and. imin == 1 .and. imax == Nvars)
    end if

end function



pure subroutine get_regressor_dims (X, trans, nvars, nobs)
    !*  GET_REGRESSOR_DIMS returns the dimensions of the matrix of regressors,
    !   depending on whether the matrix is provided in transposed form or not.
    real (PREC), intent(in), dimension(:,:) :: X
    logical, intent(in), optional :: trans
        !*  If true, regressor matrix is provided in transposed form
        !   (variables along dimension 1, observations along dimension 2)
    integer, intent(out), optional :: nvars
        !*  Number of variables
    integer, intent(out), optional :: nobs
        !*  Number of observations

    logical :: ltrans
    integer :: lnobs, lnvars, dim

    ltrans = .false.
    if (present(trans)) ltrans = trans

    ! Dimension along with observations are aligned
    dim = 1
    if (ltrans) dim = 2

    lnobs = size(X, dim)
    lnvars = size(X, 3-dim)

    if (present(nvars)) nvars = lnvars
    if (present(nobs)) nobs = lnobs

end subroutine



subroutine extract_block_2d (x, ifrom, n, trans, x_block, x_rest, status)
    !*  EXTRACT_BLOCK splits a contiguous block from a given array X,
    !   and optionally returns this block and/or the remainder of the array.
    real (PREC), intent(in), dimension(:,:), contiguous :: x
        !*  Array to be split
    integer, intent(in) :: ifrom
        !*  Initial index of the block to be extracted. The index is assumed
        !   to correspond to the first dimension, unless TRANS=.TRUE.
        !   is specified.
    integer, intent(in) :: n
        !*  Size of the contiguous block to be extracted.
    logical, intent(in), optional :: trans
        !*  If present and true, store data in transposed form in output
        !   arrays (default: false)
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: x_block
        !*  Array containing the contiguous block extracted from X.
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: x_rest
        !*  Array containing the remained of X, with X_BLOCK removed.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: i, Nvars, Nobs, Nrest
    logical :: ltrans
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'EXTRACT_BLOCK'

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (trans, .false., ltrans)

    call get_regressor_dims (x, ltrans, Nvars, Nobs)
    Nrest = Nobs - n

    if (present(x_block)) then
        call check_cond (size(x_block,1) == n .and. size(x_block,2) == Nvars, &
            NAME, 'X_BLOCK: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(x_rest)) then
        call check_cond (size(x_rest,1) == Nrest .and. size(x_rest,2) == Nvars, &
            NAME, 'X_REST: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Copy data ---

    if (present(x_block)) then
        if (ltrans) then
            do i = 1, Nvars
                x_block(1:n,i) = x(i,ifrom:ifrom+n-1)
            end do
        else
            do i = 1, Nvars
                x_block(1:n,i) = x(ifrom:ifrom+n-1,i)
            end do
        end if
    end if

    if (present(x_rest)) then
        if (ltrans) then
            do i = 1, Nvars
                x_rest(1:ifrom-1,i) = x(i,1:ifrom-1)
                x_rest(ifrom:Nrest,i) = x(i,ifrom+n:Nobs)
            end do
        else
            do i = 1, Nvars
                x_rest(1:ifrom-1,i) =  x(1:ifrom-1,i)
                x_rest(ifrom:Nrest,i) = x(ifrom+n:Nobs,i)
            end do
        end if
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine extract_block_alloc_2d (x, ifrom, n, trans, x_block, x_rest, status)
    !*  EXTRACT_BLOCK_ALLOC splits a contiguous block from a given array X,
    !   and optionally returns this block and/or the remainder of the array.
    !
    !   The arrays X_BLOCK and X_REST are (re)allocated to store the
    !   data if necessary.
    real (PREC), intent(in), dimension(:,:), contiguous :: x
        !*  Array to be split
    integer, intent(in) :: ifrom
        !*  Initial index of the block to be extracted. The index is assumed
        !   to correspond to the first dimension, unless TRANS=.TRUE.
        !   is specified.
    integer, intent(in) :: n
        !*  Size of the contiguous block to be extracted.
    logical, intent(in), optional :: trans
        !*  If present and true, store data in transposed form in output
        !   arrays (default: false)
    real (PREC), intent(inout), dimension(:,:), allocatable, optional :: x_block
        !*  Array containing the contiguous block extracted from X.
        !   If the input array does not have the exact size required to
        !   hold the data, it is (re)allocated.
    real (PREC), intent(inout), dimension(:,:), allocatable, optional :: x_rest
        !*  Array containing the contiguous block extracted from X.
        !   If the input array does not have the exact size required to
        !   hold the data, it is (re)allocated.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: Nvars, Nobs, Nrest
    integer :: shp(2)

    call get_regressor_dims (x, trans, Nvars, Nobs)
    Nrest = Nobs - n

    if (present(x_block)) then
        shp(1) = n
        shp(2) = Nvars
        call cond_alloc (x_block, shp)
    end if

    if (present(x_rest)) then
        shp(1) = Nrest
        shp(2) = Nvars
        call cond_alloc (x_rest, shp)
    end if

    ! ifort-2019 seems to have troubles handling optional, allocatable
    ! dummy arguments when passing them to routines. Split out calls
    ! depending on which arguments are actually present.
    if (present(x_block) .and. present(x_rest)) then
        call extract_block (x, ifrom, n, trans, x_block, x_rest, status)
    else if (present(x_block)) then
        call extract_block (x, ifrom, n, trans, x_block=x_block, status=status)
    else if (present(x_rest)) then
        call extract_block (x, ifrom, n, trans, x_rest=x_rest, status=status)
    end if

end subroutine



subroutine extract_block_1d (x, ifrom, n, x_block, x_rest, status)
    !*  EXTRACT_BLOCK splits a contiguous block from a given array X,
    !   and optionally returns this block and/or the remainder of the array.
    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Array to be split
    integer, intent(in) :: ifrom
        !*  Initial index of the block to be extracted.
    integer, intent(in) :: n
        !*  Size of the contiguous block to be extracted.
    real (PREC), intent(out), dimension(:), contiguous, optional :: x_block
        !*  Array containing the contiguous block extracted from X.
    real (PREC), intent(out), dimension(:), contiguous, optional :: x_rest
        !*  Array containing the remained of X, with X_BLOCK removed.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: Nobs, Nrest
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'EXTRACT_BLOCK'

    lstatus = NF_STATUS_OK

    Nobs = size(x)
    Nrest = Nobs - n

    if (present(x_block)) then
        call check_cond (size(x_block) == n, NAME, &
            'X_BLOCK: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(x_rest)) then
        call check_cond (size(x_rest) == Nrest, NAME, &
            'X_REST: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Copy data ---

    if (present(x_block)) then
        x_block(1:n) = x(ifrom:ifrom+n-1)
    end if

    if (present(x_rest)) then
        x_rest(1:ifrom-1) =  x(1:ifrom-1)
        x_rest(ifrom:Nrest) = x(ifrom+n:Nobs)
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine extract_block_alloc_1d (x, ifrom, n, x_block, x_rest, status)
    !*  EXTRACT_BLOCK_ALLOC splits a contiguous block from a given array X,
    !   and optionally returns this block and/or the remainder of the array.
    !
    !   The arrays X_BLOCK and X_REST are (re)allocated to store the
    !   data if necessary.
    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Array to be split
    integer, intent(in) :: ifrom
        !*  Initial index of the block to be extracted.
    integer, intent(in) :: n
        !*  Size of the contiguous block to be extracted.
    real (PREC), intent(inout), dimension(:), allocatable, optional :: x_block
        !*  Array containing the contiguous block extracted from X.
        !   If the input array does not have the exact size required to
        !   hold the data, it is (re)allocated.
    real (PREC), intent(inout), dimension(:), allocatable, optional :: x_rest
        !*  Array containing the contiguous block extracted from X.
        !   If the input array does not have the exact size required to
        !   hold the data, it is (re)allocated.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: Nobs, Nrest

    Nobs = size(x)
    Nrest = Nobs - n

    if (present(x_block)) then
        call cond_alloc (x_block, n)
    end if

    if (present(x_rest)) then
        call cond_alloc (x_rest, Nrest)
    end if

    ! ifort-2019 seems to have troubles handling optional, allocatable
    ! dummy arguments when passing them to routines. Split out calls
    ! depending on which arguments are actually present.
    if (present(x_block) .and. present(x_rest)) then
        call extract_block (x, ifrom, n, x_block, x_rest, status)
    else if (present(x_block)) then
        call extract_block (x, ifrom, n, x_block=x_block, status=status)
    else if (present(x_rest)) then
        call extract_block (x, ifrom, n, x_rest=x_rest, status=status)
    end if

end subroutine



pure subroutine copy_1d (src, dst)
    !*  Quick and dirty workaround to BLAS's COPY which is not PURE.
    real (PREC), intent(in), dimension(:), contiguous :: src
    real (PREC), intent(out), dimension(:), contiguous :: dst

    dst = src
end subroutine



subroutine random_sample_2d (Nobs, Nrhs, Nlhs, X, y, coefs, &
        add_intercept, add_intercept_coefs, add_intercept_x, trans_x, trans_y, &
        trans_coefs, var_error, intercept, status)
    !*  RANDOM_SAMPLE creates a random sample and coefficients
    !   with multiple outcome variables.
    !   Used for testing statistical routines.
    integer, intent(in) :: Nobs
        !*  Number of observations
    integer, value :: Nrhs
        !*  Number of RHS variables (predictors)
    integer, intent(in) :: Nlhs
        !*  Number of LHS variables (outcomes)
    real (PREC), intent(inout), dimension(:,:), allocatable :: X
        !*  Matrix of regressors, by default shaped [Nobs,Nrhs]
    real (PREC), intent(inout), dimension(:,:), allocatable :: Y
        !*  Matrix of outcome variables, by default shaped [Nobs,Nlhs]
    real (PREC), intent(inout), dimension(:,:), allocatable :: coefs
    logical, intent(in), optional :: add_intercept
    logical, intent(in), optional :: add_intercept_coefs
        !*  If present and true, add intercept as the first row or column
        !   of the COEFS array.
    logical, intent(in), optional :: add_intercept_x
        !*  If present and true, prepend a row / column of ones to X
    logical, intent(in), optional :: trans_x
        !*  If present and true, return regressor matrix in transposed form,
        !   ie. shaped [Nrhs,Nobs].
    logical, intent(in), optional :: trans_y
        !*  If present and true, return matrix of outcomes in transposed form,
        !   ie. shaped [Nlhs,Nobs]
    logical, intent(in), optional :: trans_coefs
        !*  If present and true, return coefficients matrix in transposed form,
        !   ie. shaped [Nlhs,Ncoefs]
    real (PREC), intent(in), optional :: var_error
        !*  If present, a uniformly distributed error with given variance
        !   will be added to the outcome variables.
    real (PREC), intent(inout), dimension(:), allocatable, optional :: intercept
        !*  Model intercept. Will be non-zero only if ADD_INTERCEPT=.TRUE.
        !   add ADD_INTERCEPT_X=.FALSE.
    type (status_t), intent(out), optional :: status

    logical :: ladd_intercept, ladd_intercept_coefs, ltrans_x, ltrans_y, ltrans_coefs
    logical :: ladd_intercept_x
    real (PREC), dimension(:,:), allocatable :: lX, lY, lcoefs
    real (PREC), dimension(:), allocatable :: lintercept
    real (PREC) :: lvar_error
    integer :: Nconst_coefs, i
    integer :: shp(2)
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'CREATE_RANDOM_SAMPLE'

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (add_intercept, .true., ladd_intercept)
    call set_optional_arg (add_intercept_coefs, .false., ladd_intercept_coefs)
    call set_optional_arg (add_intercept_x, .false., ladd_intercept_x)
    call set_optional_arg (trans_x, .false., ltrans_x)
    call set_optional_arg (trans_y, .false., ltrans_y)
    call set_optional_arg (trans_coefs, .false., ltrans_coefs)
    call set_optional_arg (var_error, 0.0_PREC, lvar_error)

    call check_cond (Nobs >= 0, NAME, 'NOBS: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (Nrhs >= 0, NAME, 'NVARS: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (Nlhs >= 0, NAME, 'NLHS: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (lvar_error >= 0.0_PREC, NAME, 'VAR_ERROR: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    Nconst_coefs = 0
    if (ladd_intercept_x) then
        ladd_intercept = .false.
        ladd_intercept_coefs = .false.
        Nrhs = Nrhs + 1
    else if (ladd_intercept_coefs) then
        Nconst_coefs = 1
    end if

    ! --- Create random sample ---

    allocate (lX(Nobs,Nrhs))
    allocate (lY(Nobs,Nlhs), source=0.0_PREC)
    allocate (lcoefs(Nrhs,Nlhs))
    allocate (lintercept(Nlhs), source=0.0_PREC)

    call random_number (lX)
    lx(:,:) = (lx - 0.5_PREC) * 2.0_PREC
    if (ladd_intercept_x) then
        lx(:,1) = 1.0_PREC
    end if

    if (lvar_error > 0.0_PREC) then
        call random_number (lY)
        ! Variance of standard uniform is 1/12, so rescale accordingly
        ! to get desired variance
        lY(:,:) = (lY - 0.5_PREC) * sqrt(12.0_PREC * lvar_error)
    end if

    call random_number (lcoefs)
    lcoefs(:,:) = (lcoefs - 0.5_PREC) * 5.0_PREC

    if (ladd_intercept) then
        call random_number (lintercept)
        lintercept(:) = (lintercept - 0.5_PREC) * 10.0_PREC
        do i = 1, Nlhs
            lY(:,i) = lY(:,i) + lintercept(i)
        end do
    end if

    if (Nrhs > 0) then
        call BLAS_GEMM ('N', 'N', Nobs, Nlhs, Nrhs, 1.0_PREC, lX, Nobs, &
            lcoefs, Nrhs, 1.0_PREC, lY, Nobs)
    end if

    ! --- Store output ---

    if (ltrans_x) then
        shp(1) = Nrhs
        shp(2) = Nobs
        call cond_alloc (X, shp)
        X(:,:) = transpose (lx)
    else
        shp(1) = Nobs
        shp(2) = Nrhs
        call cond_alloc (X, shp)
        X(:,:) = lX
    end if

    if (ltrans_y) then
        shp(1) = Nlhs
        shp(2) = Nobs
        call cond_alloc (Y, shp)
        Y(:,:) = transpose (lY)
    else
        shp(1) = Nobs
        shp(2) = Nlhs
        call cond_alloc (Y, shp)
        Y(:,:) = lY
    end if

    if (ltrans_coefs) then
        shp(1) = Nlhs
        shp(2) = Nrhs + Nconst_coefs
        call cond_alloc (coefs, shp)
        do i = 1, Nlhs
            coefs(1:Nlhs,Nconst_coefs+i) = lcoefs(1:Nrhs,i)
        end do
        if (ladd_intercept_coefs) then
            coefs(1:Nlhs,1) = lintercept
        end if
    else
        shp(1) = Nrhs + Nconst_coefs
        shp(2) = Nlhs
        call cond_alloc (coefs, shp)
        do i = 1, Nlhs
            coefs(1+Nconst_coefs:,i) = lcoefs(1:Nrhs,i)
        end do
        if (ladd_intercept_coefs) then
            coefs(1,1:Nlhs) = lintercept
        end if
    end if

    if (present(intercept)) then
        call cond_alloc (intercept, Nlhs)
        intercept(:) = lintercept
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine random_sample_1d (Nobs, Nrhs, X, y, coefs, &
        add_intercept, add_intercept_coefs, add_intercept_x, trans_x, &
        var_error, intercept, status)
    !*  RANDOM_SAMPLE creates a random sample and coefficients.
    !   Used for testing statistical routines.
    integer, intent(in) :: Nobs
        !*  Number of observations
    integer, value :: Nrhs
        !*  Number of RHS variables (predictors)
    real (PREC), intent(inout), dimension(:,:), allocatable :: X
        !*  Matrix of regressors, by default shaped [Nobs,Nrhs]
    real (PREC), intent(inout), dimension(:), allocatable :: y
        !*  Vector of length Nobs containing the outcome variable
    real (PREC), intent(inout), dimension(:), allocatable :: coefs
        !*  Vector of coefficients
    logical, intent(in), optional :: add_intercept
        !*  If true (default), add an intercept to the model. This will be
        !   stored in INTERCEPT.
    logical, intent(in), optional :: add_intercept_coefs
        !*  If present and true, add intercept as the first element
        !   of the COEFS array.
    logical, intent(in), optional :: add_intercept_x
        !*  If present and true, prepend a row / column of ones to X
    logical, intent(in), optional :: trans_x
        !*  If present and true, return regressor matrix in transposed form,
        !   ie. shaped [Nrhs,Nobs].
    real (PREC), intent(in), optional :: var_error
        !*  If present, a uniformly distributed error with given variance
        !   will be added to the outcome variables.
    real (PREC), intent(out), optional :: intercept
        !*  Model intercept. Will be non-zero only if ADD_INTERCEPT=.TRUE.
        !   add ADD_INTERCEPT_X=.FALSE.
    type (status_t), intent(out), optional :: status

    logical :: ladd_intercept, ladd_intercept_coefs, ladd_intercept_x, ltrans_x
    real (PREC), dimension(:,:), allocatable :: lX
    real (PREC), dimension(:), allocatable :: ly, lcoefs
    real (PREC) :: lvar_error
    real (PREC) :: lintercept
    integer :: Nconst_coefs
    integer :: shp(2)
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'CREATE_RANDOM_SAMPLE'

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (add_intercept, .true., ladd_intercept)
    call set_optional_arg (add_intercept_coefs, .false., ladd_intercept_coefs)
    call set_optional_arg (add_intercept_x, .false., ladd_intercept_x)
    call set_optional_arg (trans_x, .false., ltrans_x)
    call set_optional_arg (var_error, 0.0_PREC, lvar_error)

    call check_cond (Nobs >= 0, NAME, 'NOBS: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (Nrhs >= 0, NAME, 'NVARS: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (lvar_error >= 0.0_PREC, NAME, 'VAR_ERROR: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    Nconst_coefs = 0
    if (ladd_intercept_x) then
        ladd_intercept_coefs = .false.
        ladd_intercept = .false.
        Nrhs = Nrhs + 1
    else if (ladd_intercept_coefs) then
        Nconst_coefs = 1
    end if

    ! --- Create random sample ---

    allocate (lX(Nobs,Nrhs))
    allocate (lY(Nobs), source=0.0_PREC)
    allocate (lcoefs(Nrhs))

    call random_number (lX)
    lx(:,:) = (lx - 0.5_PREC) * 2.0_PREC
    if (ladd_intercept_x) then
        lx(:,1) = 1.0_PREC
    end if

    if (lvar_error > 0.0_PREC) then
        call random_number (lY)
        ! Variance of standard uniform is 1/12, so rescale accordingly
        ! to get desired variance
        lY(:) = (lY - 0.5_PREC) * sqrt(12.0_PREC * lvar_error)
    end if

    call random_number (lcoefs)
    lcoefs(:) = (lcoefs - 0.5_PREC) * 5.0_PREC

    lintercept = 0.0
    if (ladd_intercept) then
        call random_number (lintercept)
        lintercept = (lintercept - 0.5_PREC) * 10.0_PREC
        ly(:) = ly + lintercept
    end if

    if (Nrhs > 0) then
        call BLAS_GEMV ('N', Nobs, Nrhs, 1.0_PREC, lX, Nobs, &
            lcoefs, 1, 1.0_PREC, lY, 1)
    end if

    ! --- Store output ---

    if (ltrans_x) then
        shp(1) = Nrhs
        shp(2) = Nobs
        call cond_alloc (X, shp)
        X(:,:) = transpose (lx)
    else
        shp(1) = Nobs
        shp(2) = Nrhs
        call cond_alloc (X, shp)
        X(:,:) = lX
    end if

    call cond_alloc (Y, Nobs)
    Y(:) = lY

    call cond_alloc (coefs, Nrhs + Nconst_coefs)
    coefs(1+Nconst_coefs:) = lcoefs

    if (ladd_intercept_coefs) then
        coefs(1) = lintercept
    end if

    if (present(intercept)) then
        intercept = lintercept
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine
