

subroutine lm_result_update (self, conf, model, coefs, nobs, nrhs, nlhs, &
        ncomp, var_rhs, rank_rhs)

    type (lm_result), intent(inout) :: self
    integer, intent(in), optional :: model
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:), optional :: coefs
    integer, intent(in), optional :: nobs
    integer, intent(in), optional :: nrhs
    integer, intent(in), optional :: nlhs
    integer, intent(in), optional :: ncomp
    real (PREC), intent(in), optional :: var_rhs
    integer, intent(in), optional :: rank_rhs

    if (present(coefs)) then
        call copy_alloc (coefs, self%coefs)
    end if

    if (present(model)) self%model = model
    if (present(conf)) self%conf = conf
    if (present(nobs)) self%nobs = nobs
    if (present(nrhs)) self%nrhs = nrhs
    if (present(nlhs)) self%nlhs = nlhs
    if (present(ncomp)) self%ncomp = ncomp
    if (present(var_rhs)) self%var_rhs = var_rhs
    if (present(rank_rhs)) self%rank_rhs = rank_rhs
end subroutine



subroutine lm_result_reset (self)
    type (lm_result), intent(out) :: self

    ! INTENT(OUT) makes sure that all allocatable attribute arrays are
    ! deallocated.

    ! Overwrite with default LM_CONFIG, whatever that may be
    type (lm_config) :: conf
    self%conf = conf

    self%intercept = 0.0

    self%model = 0

    self%nrhs = 0
    self%nlhs = 0
    self%nobs = 0
    self%ncomp = 0

    self%rank_rhs = 0
    self%var_rhs = 0.0

end subroutine



pure subroutine get_dims_2d (X, Y, coefs, intercept, add_intercept, trans_x, &
        trans_y, trans_coefs, weights, mask, nobs, nrhs, nlhs, ncoefs, nconst_add, &
        nconst_coefs, status)
    !*  GET_DIMS returns the dimensions corresponding to various arrays
    !   used in linear models.
    real (PREC), intent(in), dimension(:,:) :: X
    real (PREC), intent(in), dimension(:,:) :: Y
    real (PREC), intent(in), dimension(:,:), optional :: coefs
    real (PREC), intent(in), dimension(:), optional :: intercept
    logical, intent(in), optional :: add_intercept
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
    logical, intent(in), optional :: trans_coefs
    real (PREC), intent(in), dimension(:), optional :: weights
    logical, intent(in), dimension(:,:), optional :: mask
    integer, intent(out), optional :: nobs, nrhs, nlhs
    integer, intent(out), optional :: ncoefs
        !*  Number of coefficients stored in COEFS array.
    integer, intent(out), optional :: nconst_add
        !*  Number of constants ADDED to the model (irrespective of what is
        !   in the user-provided regressor matrix). Either 0 or 1, the latter
        !   being the case when ADD_INTERCEPT=.TRUE.
    integer, intent(out), optional :: nconst_coefs
        !*  Number of constants added to COEFS arrays. This value satisfies
        !   0 <= Nconst_coefs <= Nconst <= 1.
        !   Non-zero value can only be returned if COEFS is present, otherwise
        !   it is impossible to determine from the inputs.
    type (status_t), intent(out), optional :: status
        !*  Exit code. Returns NF_STATUS_INVALID_ARG if any of the dimensions
        !   are non-conformable with the given array sizes.

    logical :: ltrans_x, ltrans_y, ltrans_coefs
    logical :: dim_ok
    integer :: lnrhs, lnlhs, lncoefs, lnobs, lnobs_lhs, dim
    integer :: lnconst_add, lnconst_coefs
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'GET_DIMS'
    type (lm_config) :: conf

    lstatus = NF_STATUS_OK

    ! Use CONF values as default values
    call set_optional_arg (trans_x, conf%trans_x, ltrans_x)
    call set_optional_arg (trans_y, conf%trans_y, ltrans_y)
    call set_optional_arg (trans_coefs, conf%trans_coefs, ltrans_coefs)

    lnobs = 0
    lnrhs = 0
    lnlhs = 0
    lncoefs = 0
    lnconst_add = 0
    lnconst_coefs = 0

    if (present(add_intercept)) then
        if (add_intercept) lnconst_add = 1
    end if

    if (.not. ltrans_x) then
        lnobs = size(x, 1)
        lnrhs = size(x, 2)
    else
        lnobs = size(x, 2)
        lnrhs = size(x, 1)
    end if

    if (ltrans_y) then
        lnlhs = size(Y, 1)
        lnobs_lhs = size(Y, 2)
    else
        lnlhs = size(Y, 2)
        lnobs_lhs = size(Y, 1)
    end if

    if (present(coefs)) then
        ! Dimensions of COEFS array
        dim = 1
        if (ltrans_coefs) dim = 2

        lncoefs = size(coefs,dim)
        lnconst_coefs = lncoefs - lnrhs

        if (lnconst_coefs < 0 .or. lnconst_coefs > 1) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if

        if (present(add_intercept)) then
            if (.not. add_intercept) then
                ! Intercept should not be added, so don't allow COEFS
                ! to have one additional row/column
                dim_ok = (lncoefs == lnrhs) .and. (size(coefs,3-dim) == lnlhs)
            else
                ! Intercept to be added to model, but user still has the option
                ! of having it stored in COEFS or not, so both are acceptable
                dim_ok = ((lncoefs == lnrhs) .or. (lncoefs == (lnrhs + 1))) &
                    .and. (size(coefs,3-dim) == lnlhs)
            end if
        else
            ! ADD_INTERCEPT not present, so we don't know whether COEFS
            ! has one row/column more than NRHS.
            ! This is relevant when called from the PREDICT routines.
            dim_ok = (lncoefs == (lnrhs + 1) .or. lncoefs == lnrhs) &
                .and. size(coefs,3-dim) == lnlhs
        end if

        call check_cond (dim_ok, NAME, 'COEFS: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    call check_cond (lnobs == lnobs_lhs, NAME, 'Y: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(intercept)) then
        call check_cond (size(intercept) >= lnlhs, NAME, &
            'INTERCEPT: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(weights)) then
        call check_cond (size(weights) == lnobs, NAME, &
            'WEIGHTS: Non-conformable arrays', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
        ! Reject invalid weights before callers normalize or take square roots.
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

    if (present(mask)) then
        call check_cond (size(mask,1) == size(Y,1) .and. size(mask,2) == size(Y,2), &
            NAME, 'MASK: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

100 continue

    if (present(nobs)) nobs = lnobs
    if (present(nrhs)) nrhs = lnrhs
    if (present(nlhs)) nlhs = lnlhs
    if (present(ncoefs)) ncoefs = lncoefs
    if (present(nconst_add)) nconst_add = lnconst_add
    if (present(nconst_coefs)) nconst_coefs = lnconst_coefs

    if (present(status)) status = lstatus

end subroutine



subroutine ols_2d (conf, X, Y, coefs, intercept, weights, rank, res, status)
    !*  OLS_2D computes the ordinary least-squares problem for given independent
    !   data X and (potentially multiple) dependent variables Y.
    !
    !   The LS problem is solved using SVD as implemented in LAPACK's GELSD
    !   routine. The optional arguments RCOND and RANK are passed directly
    !   to/from GELSD.
    !
    !   Intercepts
    !   ----------
    !   The routine supports the following scenarios:
    !       1.  CONF%ADD_INTERCEPT = .FALSE. and INTERCEPT not present.
    !           No intercept should be added or returned: the model is
    !           estimated with an intercept if and only if the user
    !           manually adds a constant to X.
    !       2.  CONF%ADD_INTERCEPT = .TRUE. and INTERCEPT not present:
    !           The routine adds an intercept to the model.
    !       3.  CONF%ADD_INTERCEPT = .FALSE. and INTERCEPT is present:
    !           The routine adds an intercept to the model and its estimates
    !           are stored in INTERCEPT.
    !       4.  ADD_INTERCEPT = .TRUE. and INTERCEPT is present:
    !           The routine adds an intercept to the model and its estimates
    !           are stored INTERCEPT.
    !   Consequently, in cases (2) - (4), the regressor matrix is augmented
    !   to be [1, X], and the RANK output argument will be the rank of
    !   this matrix.
    !
    !   Additionally, if the array COEFS has shape [NRHS+1,NLHS], the
    !   estimated intercepts will ADDITIONALLY be stored in the first
    !   row of COEFS.
    !
    !   Note that the routine creates copies for input arrays X and Y as these
    !   will be overwritten by GELSD.
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
        !*  Array of RHS variables
    real (PREC), intent(in), dimension(:,:), contiguous :: Y
        !*  Array of LHS variables (separate regression is performed for each
        !   LHS variables using the same set of RHS variables)
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: coefs
        !*  Array of estimated coefficients. Argument is optional in case
        !   coefficients should only be stored in RES argument.
    real (PREC), intent(out), dimension(:), contiguous, optional :: intercept
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
    integer, intent(out), optional :: rank
        !*  Contains effective rank of regressor matrix
    type (lm_result), intent(inout), optional :: res
        !*  Optional result object for linear models.
    type (status_t), intent(out), optional :: status
        !*  Optional exit code

    real (PREC) :: var_rhs
    integer :: nobs, Nrhs, Ncoefs, nconst_add, nlhs, i
    type (status_t) :: lstatus
    type (lm_config) :: lconf
    real (PREC), dimension(:), allocatable :: weights_sqrt
    character (*), parameter :: NAME = 'OLS'

    ! GELSD arguments
    real (PREC), dimension(:), allocatable :: work, sval
    integer, dimension(:), allocatable :: iwork
    real (PREC), dimension(:,:), allocatable :: Xp, Yp
    integer :: lrank
    integer :: lwork, info, m, n, lda, ldb, mn, liwork

    lstatus = NF_STATUS_OK

    ! --- Input processsing ---

    if (present(conf)) lconf = conf

    call check_cond (lconf%rcond >= 0.0_PREC, NAME, 'CONF: Invalid value for RCOND', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call get_dims (X, Y, coefs, intercept, lconf%add_intercept, lconf%trans_x, &
        lconf%trans_y, trans_coefs=lconf%trans_coefs, weights=weights, nobs=nobs, &
        nrhs=Nrhs, nlhs=Nlhs, ncoefs=Ncoefs, nconst_add=Nconst_add, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (Nobs >= (Nrhs + Nconst_add), NAME, 'Too few observations', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(intercept)) then
        call check_cond (size(intercept) >= nlhs, NAME, &
            'INTERCEPT: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

   ! --- Pre-process data ---

    ! Allocate (possibly transposed) array to store X variables; will be
    ! overwritten by GELSD, so we have to allocate in any case.
    ! Add constant as requested.
    allocate (Xp(nobs,Nrhs+Nconst_add))

    call transform_regressors (X, Xp, trans=lconf%trans_x, drop_na=.false., &
        add_intercept=lconf%add_intercept, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (lconf%trans_y) then
        allocate (Yp(nobs,Nlhs))
        call transform_regressors (Y, Yp, trans=lconf%trans_y, drop_na=.false., &
            status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        ! Just create a copy since Y would otherwise be overwritten in place
        allocate (Yp(nobs,nlhs), source=Y)
    end if

    ! --- Apply weights ---

    if (present(weights)) then
        allocate (weights_sqrt(Nobs), source=weights)
        weights_sqrt(:) = weights_sqrt / sum(weights_sqrt)
        weights_sqrt(:) = sqrt(weights_sqrt)

        do i = 1, Nrhs + Nconst_add
            Xp(:,i) = Xp(:,i) * weights_sqrt
        end do

        do i = 1, Nlhs
            Yp(:,i) = Yp(:,i) * weights_sqrt
        end do
    end if

    ! --- GELSD ---

    m = nobs
    n = Nrhs + Nconst_add
    mn = max(1, min(m, n))
    lda = nobs
    ldb = nobs
    ! GELSD solves Ax=b, we solve Y = Xb, so our nlhs = nrhs in GELSD notation
    allocate (sval(mn))

    ! 1. Workspace query
    lwork = -1
    allocate (work(1))
    allocate (iwork(1))

    call LAPACK_GELSD (m, n, nlhs, Xp, lda, Yp, ldb, sval, lconf%rcond, lrank, &
        work, lwork, iwork, info)

    if (info < 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! 2. Allocate optimal / minimal work arrays
    lwork = int(work(1))
    liwork = iwork(1)

    deallocate (work, iwork)
    allocate (work(lwork), iwork(liwork))

    ! 3. Actual call to perform LS
    call LAPACK_GELSD (m, n, nlhs, Xp, lda, Yp, ldb, sval, lconf%rcond, lrank, &
        work, lwork, iwork, info)

    ! Check whether algorithm for computing SVD failed to converge
    if (info > 0) then
        lstatus = NF_STATUS_NOT_CONVERGED
        goto 100
    else if (info < 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! --- Store results ---

    ! Copy coefficients
    if (present(coefs)) then
        if (lconf%trans_coefs) then
            if (Ncoefs == Nrhs) then
                ! No intercept should be added to COEFS
                do i = 1, Nrhs
                    coefs(1:Nlhs,i) = Yp(Nconst_add+i,1:Nlhs)
                end do
            else
                ! Include added intercept in COEFS
                do i = 1, Ncoefs
                    coefs(1:Nlhs,i) = Yp(i,1:Nlhs)
                end do
            end if
        else
            if (Ncoefs == Nrhs) then
                ! No intercept should be added to COEFS
                do i = 1, nlhs
                    coefs(1:Nrhs,i) = Yp(1+Nconst_add:Nrhs+Nconst_add,i)
                end do
            else
                ! Include added intercept in COEFS
                do i = 1, nlhs
                    coefs(1:Ncoefs,i) = Yp(1:Ncoefs,i)
                end do
            end if
        end if
    end if

    if (present(intercept)) then
        ! Costant was added to X, intercept is in the first row of Yp
        intercept = 0.0
        if (lconf%add_intercept) then
            intercept(1:Nlhs) = Yp(1,:)
        end if
    end if

    ! Copy over optional output arguments
    if (present(rank)) rank = lrank

    if (present(res)) then
        ! Update LM_RESULT object for OLS model

        ! Fraction of RHS variance used, analogous to PCR regression
        ! (only applicable if RHS matrix does not have full rank)
        var_rhs = sum(sval(1:lrank) ** 2.0_PREC) / sum(sval**2.0_PREC)

        call lm_result_update (res, conf=lconf, model=NF_STATS_LM_OLS, &
            nobs=nobs, nrhs=Nrhs, nlhs=Nlhs, var_rhs=var_rhs, rank_rhs=lrank)

        ! Never store intercept in first row in RES%COEFS_MULTI array
        call copy_alloc (Yp(1+Nconst_add:Nrhs+Nconst_add,1:Nlhs), res%coefs_multi)

        if (lconf%add_intercept) then
            call copy_alloc (Yp(1,1:Nlhs), res%intercept_multi)
        else
            ! We cannot leave the array unallocated or initialized, make
            ! sure that we set intercepts to zero.
            call cond_alloc (res%intercept_multi, Nlhs)
            res%intercept_multi(:) = 0.0
        end if
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine



subroutine ols_1d (conf, X, y, coefs, intercept, weights, rank, res, status)
    !*  OLS_1D provides a convenient wrapper for OLS_2D for one-dimensional
    !   input data (ie for regressions with a single dependent variable).
    !
    !   See the documentation for OLS_2D for details.
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous, target :: y
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: coefs
    real (PREC), intent(out), optional :: intercept
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
    integer, intent(out), optional :: rank
    type (lm_result), intent(inout), optional :: res
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), contiguous, pointer :: ptr_y
    real (PREC), dimension(:,:), contiguous, pointer :: ptr_coefs
    real (PREC), dimension(:), allocatable :: lintercept
    type (status_t) :: lstatus
    type (lm_config) :: lconf

    nullify (ptr_y, ptr_coefs)

    ptr_y(1:size(y),1:1) => y

    if (present(conf)) then
        ! Enforce no transposing of data in 1d arrays
        lconf = conf
        lconf%trans_y = .false.
        lconf%trans_coefs = .false.
    end if

    if (present(coefs)) then
        ptr_coefs(1:size(coefs),1:1) => coefs
    end if

    if (present(intercept)) then
        allocate (lintercept(1))
    end if

    call ols (conf, X, ptr_y, ptr_coefs, lintercept, weights, rank, res, lstatus)

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
    end if

    if (present(intercept)) then
        intercept = lintercept(1)
    end if

    if (present(status)) status = lstatus

end subroutine



subroutine pca (X, trans_x, center, scale, scores, sval, loadings, shift_x, &
        scale_x, propvar, status)

    real (PREC), intent(in), dimension(:,:), contiguous :: X
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: center
    logical, intent(in), optional :: scale
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: scores
    real (PREC), intent(out), dimension(:), contiguous, optional :: sval
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: loadings
    real (PREC), intent(out), dimension(:), contiguous, optional :: shift_x
    real (PREC), intent(out), dimension(:), contiguous, optional :: scale_x
    real (PREC), intent(out), dimension(:), contiguous, optional :: propvar
    type (status_t), intent(out), optional :: status

    logical :: lscale, lcenter, ltrans_x
    type (status_t) :: lstatus
    real (PREC), dimension(:,:), allocatable :: Xp
    real (PREC), dimension(:), allocatable :: work
    integer :: Nvars, Nobs, i, Ncomp
    character (*), parameter :: NAME = 'PCA'

    ! Argument for GESDD
    real (PREC), dimension(:), allocatable :: ls
    real (PREC), dimension(:,:), allocatable :: vt
    real (PREC), dimension(1) :: qwork
    integer, dimension(:), allocatable :: iwork
    character (1), parameter :: jobz = 'O'
    real (PREC), dimension(0,0) :: u
    integer :: lda, ldu, ldvt, lwork, info, m, n, mn

    lstatus = NF_STATUS_OK

    ! --- Input processing ---
    call set_optional_arg (center, .true., lcenter)
    call set_optional_arg (scale, .true., lscale)
    call set_optional_arg (trans_x, .false., ltrans_x)

    call get_regressor_dims (X, ltrans_x, Nvars=Nvars, Nobs=Nobs)

    ! --- Quick termination ---

    call check_cond (Nobs >= Nvars, NAME, 'X: too few observations', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Check remaining inputs ---

    if (present(scores)) then
        ! Do not check second dimension, we will simply fill in as many scores
        ! as there are columns
        call check_cond (size(scores, 1) == Nobs, NAME, &
            'SCORES: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(loadings)) then
        ! Don't check second dimension, we will simply fill in as many
        ! principal components as there are columns
        call check_cond (size(loadings,1) == Nvars, NAME, &
            'LOADINGS: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Pre-process data ---

    ! 1. Standardize data in place; hence we first need to create a copy,
    ! accounting for the possibility that X needs to be transformed.
    ! GESDD overwrites its inputs, so we need a copy of X in any case!

    allocate (Xp(Nobs,Nvars))

    call transform_regressors (X, Xp, trans=ltrans_x, center=lcenter, &
        scale=lscale, shift_x=shift_x, scale_x=scale_x, drop_na=.false., &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- GESDD call ---


    ! Notes on calling GESDD: we are using JOBZ='O' so that only the first
    ! NVARS columns of the matrix U in the SVD X=USV' will be computed and
    ! stored in the array X_TMP.
    ! The input checks ensure that NOBS >= NVARS, so the matrix V' will always
    ! be written into the array VT and never into X_TMP.
    ! For JOBZ='O' and m >= n, the array U will never be referenced, so we
    ! can just pass a size-0 dummy array.
    m = Nobs
    n = Nvars
    lda = Nobs
    ldvt = n
    ldu = m
    mn = max(1, min(m, n))

    ! These arrays need to be allocated irrespective of whether GESDD is
    ! called or not
    allocate (vt(n, n))
    allocate (ls(mn))

    ! Skip the remaining SVD decomposition of no variables are present
    if (Nvars == 0) goto 50

    allocate (iwork(8*mn))

    ! workspace query
    lwork = -1
    call LAPACK_GESDD (jobz, m, n, Xp, lda, ls, u, ldu, vt, ldvt, qwork, &
        lwork, iwork, info)

    ! Recover minimal work space size
    if (info /= 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if
    lwork = int(qwork(1))

    ! perform actual SVD
    allocate (work(lwork))

    call LAPACK_GESDD (jobz, m, n, Xp, lda, ls, u, ldu, vt, ldvt, work, &
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

    ! --- Store results ---

    ! Jump to this point if SVD decomposition is skipped
50  continue

    if (present(scores)) then
        ! Compute scores as S(:,i) = U(:,i) * sval(i) = X * V(:,i)
        ncomp = min(Nvars, size(scores, 2))
        do i = 1, ncomp
            scores(1:Nobs,i) = Xp(1:Nobs,i) * ls(i)
        end do

        ! Set any excess columns to zero
        scores(:,ncomp+1:) = 0.0
    end if

    if (present(sval)) then
        sval = 0.0
        ncomp = min(Nvars, size(sval))
        sval(1:ncomp) = ls(1:ncomp)
    end if

    if (present(loadings)) then
        ncomp = min(Nvars, size(loadings, 2))
        ! V^T returned by GESDD has shape [NCOMP,NVARS] where in that routine
        ! NCOMP = NVARS. If we want to restrict to NCOMP < NVARS, we have
        ! to use the first 1:NCOMP rows.
        do i = 1, ncomp
            loadings(1:Nvars,i) = vt(i,1:Nvars)
        end do
        ! Set any excess columns to zero
        loadings(:,ncomp+1:) = 0.0
    end if

    ! compute vector of variance share captured by each PC.
    ! Variance share is given by relative size of each eigenvalue (ie. squared
    ! singular value) relative to sum of all EVs.
    if (present(propvar)) then
        propvar(:) = 0.0
        if (allocated(ls)) then
            ncomp = min(Nvars, size(propvar), size(ls))
            propvar(1:ncomp) = ls(1:ncomp) ** 2 / sum(ls(1:Nvars) ** 2)
        end if
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine



pure subroutine pcr_get_dims (Y, coefs, scores, sval, loadings, mean_x, &
        scale_x, add_intercept, trans_y, trans_coefs, nobs, nrhs, &
        ncomp, nlhs, ncoefs, nconst_add, nconst_coefs, status)
    !*  PCR_GET_DIMS returns the number of various dimensions of the PCR
    !   problem, depending on input data arrays and whether an intercept
    !   shoud be added, etc.
    real (PREC), intent(in), dimension(:,:) :: Y
    real (PREC), intent(in), dimension(:,:) :: coefs
    real (PREC), intent(in), dimension(:,:) :: scores
    real (PREC), intent(in), dimension(:) :: sval
    real (PREC), intent(in), dimension(:,:) :: loadings
    real (PREC), intent(in), dimension(:), optional :: mean_x, scale_x
    logical, intent(in), optional :: add_intercept
    logical, intent(in) :: trans_y
    logical, intent(in) :: trans_coefs
    integer, intent(out), optional :: nobs, nrhs, nlhs, ncoefs, nconst_add, ncomp
    integer, intent(out), optional :: nconst_coefs
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: lnobs, lnrhs, lnlhs, lncoefs, lnconst_add, lncomp, dim
    integer :: lnconst_coefs
    logical :: dim_ok
    character (*), parameter :: NAME = 'PCR_GET_DIMS'

    lstatus = NF_STATUS_INVALID_ARG

    lnobs = 0
    lnrhs = 0
    lnlhs = 0
    lncoefs = 0
    lnconst_add = 0
    lnconst_coefs = 0
    lncomp = 0

    if (present(add_intercept)) then
        if (add_intercept) lnconst_add = 1
    end if

    dim = 1
    if (trans_y) dim = 2
    lnobs = size(Y, dim)
    lnlhs = size(Y, 3-dim)

    lnrhs = size(loadings, 1)
    lncomp = size(scores, 2)

    if (size(scores,1) /= lnobs) goto 100
    if (size(loadings,2) /= lncomp) goto 100

    ! --- Check COEFS arrays ---

    dim = 1
    if (trans_coefs) dim = 2
    lncoefs = size(coefs,dim)
    lnconst_coefs = lncoefs - lnrhs

    if (lnconst_coefs < 0 .or. lnconst_coefs > 1) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (present(add_intercept)) then
        if (.not. add_intercept) then
            ! Intercept should not be added, so don't allow COEFS
            ! to have one additional row/column
            dim_ok = (lncoefs == lnrhs) .and. (size(coefs,3-dim) == lnlhs)
        else
            ! Intercept to be added to model, but user still has the option
            ! of having it stored in COEFS or not, so both are acceptable
            dim_ok = ((lncoefs == lnrhs) .or. (lncoefs == (lnrhs + 1))) &
                .and. (size(coefs,3-dim) == lnlhs)
        end if
    else
        ! ADD_INTERCEPT not present, so we don't know whether COEFS
        ! has one row/column more than NRHS.
        ! This is relevant when called from the PREDICT routines.
        dim_ok = (lncoefs == (lnrhs + 1) .or. lncoefs == lnrhs) &
            .and. size(coefs,3-dim) == lnlhs
    end if

    call check_cond (dim_ok, NAME, 'COEFS: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(mean_x)) then
        if (size(mean_x) /= lnrhs) goto 100
    end if

    if (present(scale_x)) then
        if (size(scale_x) /= lnrhs) goto 100
    end if

    lstatus = NF_STATUS_OK

100 continue

    if (present(nobs)) nobs = lnobs
    if (present(nrhs)) nrhs = lnrhs
    if (present(nlhs)) nlhs = lnlhs
    if (present(ncoefs)) ncoefs = lncoefs
    if (present(nconst_add)) nconst_add = lnconst_add
    if (present(nconst_coefs)) nconst_coefs = lnconst_coefs
    if (present(ncomp)) ncomp = lncomp

    if (present(status)) status = lstatus

end subroutine




subroutine pcr_2d (Y, scores, sval, loadings, coefs, shift_x, scale_x, &
        add_intercept, trans_y, trans_coefs, center, intercept, weights_sqrt, status)
    real (PREC), intent(in), dimension(:,:), contiguous, target :: Y
    real (PREC), intent(in), dimension(:,:), contiguous :: scores
    real (PREC), intent(in), dimension(:), contiguous :: sval
    real (PREC), intent(in), dimension(:,:), contiguous :: loadings
    real (PREC), intent(out), dimension(:,:), contiguous :: coefs
    real (PREC), intent(in), dimension(:), contiguous, optional, target :: shift_x
    real (PREC), intent(in), dimension(:), contiguous, optional, target :: scale_x
    logical, intent(in), optional :: add_intercept
    logical, intent(in), optional :: trans_y
    logical, intent(in), optional :: trans_coefs
    logical, intent(in), optional :: center
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: intercept
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights_sqrt
    type (status_t), intent(out), optional :: status

    logical :: lcenter, ltrans_y, ltrans_coefs
    integer :: Nrhs, Nobs, ncomp, Nlhs, Ncoefs, Nconst_coefs, dim
    real (PREC), dimension(:), pointer, contiguous :: ptr_intercept
    integer :: i
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'PCR'

    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Yp
    real (PREC), dimension(:,:), allocatable :: PCy, lcoefs
    real (PREC), dimension(:), allocatable :: shift_y, lweights
    ! Arguments to GEMM
    character (1) :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    real (PREC) :: alpha, beta
    ! Arguments to GEMV

    nullify (ptr_intercept, ptr_Yp)

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (center, .true., lcenter)
    call set_optional_arg (trans_y, .false., ltrans_y)
    call set_optional_arg (trans_coefs, .false., ltrans_coefs)

    call pcr_get_dims (Y, coefs, scores, sval, loadings, shift_x, scale_x, &
        add_intercept, ltrans_y, ltrans_coefs, nobs=Nobs, nrhs=Nrhs, &
        ncomp=Ncomp, nlhs=Nlhs, ncoefs=Ncoefs, nconst_coefs=Nconst_coefs, &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(weights_sqrt)) then
        call check_cond (size(weights_sqrt) == Nobs, NAME, &
            'WEIGHTS: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Pre-process data ---

    if (ltrans_y .or. lcenter .or. present(weights_sqrt)) then
        allocate (ptr_Yp(Nobs,Nlhs))
        allocate (shift_y(Nlhs))
        if (present(weights_sqrt)) then
            allocate (lweights(Nobs))
            lweights(:) = weights_sqrt**2.0_PREC
            call transform_regressors (Y, ptr_Yp, center=lcenter, scale=.false., &
                trans=ltrans_y, shift_x=shift_y, weights=lweights, status=lstatus)
        else
            call transform_regressors (Y, ptr_Yp, center=lcenter, scale=.false., &
                trans=ltrans_y, shift_x=shift_y, status=lstatus)
        end if
        if (lstatus /= NF_STATUS_OK) goto 100

        if (present(weights_sqrt)) then
            do i = 1, Nlhs
                ptr_Yp(:,i) = ptr_Yp(:,i) * weights_sqrt
            end do
        end if
    else
        ptr_Yp => Y
    end if

    ! If no components are present, we can still run a constant-only regression
    ! if requested by user. Otherwise exit doing nothing.
    if (ncomp == 0) then
        allocate (lcoefs(Nrhs,Nlhs), source=0.0_PREC)
        goto 50
    end if

    ! --- Compute slope coefficients ---

    ! Compute regression coefficients of regression Y on PC
    ! 1. Compute PC'y
    m = ncomp
    n = Nlhs
    k = Nobs
    lda = Nobs
    ldb = Nobs
    ldc = ncomp
    transa = 'T'
    transb = 'N'
    alpha = 1.0_PREC
    beta = 0.0_PREC
    allocate (PCy(ncomp, Nlhs))
    call BLAS_GEMM (transa, transb, m, n, k, alpha, scores, lda, ptr_Yp, ldb, beta, PCy, ldc)

    ! 2. Regression coefficients are given by (diag(sval**2))^{-1} PC'y,
    ! can be computed directly without solving equation system.
    do i = 1, Nlhs
        PCy(1:ncomp, i) = PCy(1:ncomp, i) * (sval(1:ncomp) ** (-2.0_PREC))
    end do

    ! PCR coefficients derived from regressing on k PC are given by
    ! beta_pcr = V_k * beta_ols
    m = Nrhs
    n = Nlhs
    k = ncomp
    lda = Nrhs
    ldb = ncomp
    ldc = Nrhs
    transa = 'N'
    transb = 'N'
    alpha = 1.0_PREC
    beta = 0.0_PREC
    allocate (lcoefs(Nrhs, Nlhs))
    call BLAS_GEMM (transa, transb, m, n, k, alpha, loadings, lda, PCy, ldb, beta, lcoefs, ldc)

    ! Scale back betas if explanatory variables were normalized
    if (present(scale_x)) then
        do i = 1, Nlhs
            lcoefs(1:Nrhs, i) = lcoefs(1:Nrhs, i) / scale_x(1:Nrhs)
        end do
    end if

50 continue

    ! --- Recover intercept ---

    dim = 1
    if (ltrans_coefs) dim = 2

    if (present(intercept) .or. size(coefs,1) > Nrhs) then

        call assert_alloc_ptr (intercept, Nlhs, ptr_intercept)

        ptr_intercept(:) = 0.0

        if (allocated(shift_y)) then
            ! If LHS variables were centered in this routine, initialize
            ! intercept using their mean (ie the value by which they were
            ! shifted).
            ! SHIFT_Y will be zeros if Y is pre-centered by TRANS_Y=.TRUE.
            ptr_intercept(1:Nlhs) = shift_y
        end if

        if (present(shift_x) .and. ncomp > 0) then
            ! Compute intercept as a = mean(y) - mean(x)'b
            ! where mean(x)'b is stored in temporary array work.
            ! Manually create temporary array as coefs(1:1,:) is not contiguous
            ! in memory and ifort issues runtime warning in debug mode every time.

            call BLAS_GEMV ('T', Nrhs, Nlhs, -1.0_PREC, lcoefs, Nrhs, &
                shift_x, 1, 1.0_PREC, ptr_intercept, 1)
        end if

        if (Nconst_coefs == 1) then
            if (ltrans_coefs) then
                coefs(1:Nlhs,1) = ptr_intercept
            else
                coefs(1,1:Nlhs) = ptr_intercept
            end if
        end if
    end if

    ! --- Store output ---

    if (ltrans_coefs) then
        do i = 1, Nrhs
            coefs(1:Nlhs,Nconst_coefs+i) = lcoefs(i,1:Nlhs)
        end do
    else
        do i = 1, Nlhs
            coefs(1+Nconst_coefs:Ncoefs,i) = lcoefs(:,i)
        end do
    end if

100 continue

    call assert_dealloc_ptr (Y, ptr_Yp)
    call assert_dealloc_ptr (intercept, ptr_intercept)

    if (present(status)) status = lstatus

end subroutine



subroutine pcr_1d (y, scores, sval, loadings, coefs, shift_x, scale_x, &
        add_intercept, center, intercept, weights_sqrt, status)

    real (PREC), intent(in), dimension(:), contiguous, target :: y
    real (PREC), intent(in), dimension(:,:), contiguous :: scores
    real (PREC), intent(in), dimensioN(:), contiguous :: sval
    real (PREC), intent(in), dimension(:,:), contiguous :: loadings
    real (PREC), intent(out), dimension(:), contiguous, target :: coefs
    real (PREC), intent(in), dimension(:), contiguous, optional :: shift_x
    real (PREC), intent(in), dimension(:), contiguous, optional :: scale_x
    logical, intent(in), optional :: add_intercept
    logical, intent(in), optional :: center
    real (PREC), intent(out), optional :: intercept
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights_sqrt
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), contiguous, pointer :: ptr_y, ptr_coefs
    real (PREC), dimension(:), allocatable :: lintercept
    integer :: ncoefs, nobs
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    nobs = size(y)
    ncoefs = size(coefs)

    ptr_y(1:nobs,1:1) => y
    ptr_coefs(1:ncoefs,1:1) => coefs

    if (present(intercept)) then
        allocate (lintercept(1))
    end if

    call pcr (ptr_y, scores, sval, loadings, ptr_coefs, shift_x, scale_x, &
        add_intercept, center=center, intercept=lintercept, &
        weights_sqrt=weights_sqrt, status=lstatus)

    if (lstatus == NF_STATUS_OK) then
        if (present(intercept)) then
            intercept = lintercept(1)
        end if
    end if

    if (present(status)) status = lstatus

end subroutine




subroutine pcr_pca_2d (conf, X, Y, coefs, ncomp, intercept, weights, res, status)
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:,:), contiguous :: Y
    real (PREC), intent(out), dimension(:,:), contiguous, optional, target :: coefs
    integer, intent(out), optional :: ncomp
        !   If present, contains the actual number of PCs used on exit.
    real (PREC), intent(out), dimension(:), optional, contiguous, target :: intercept
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
    type (lm_result), intent(inout), optional, target :: res
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    type (lm_config) :: lconf
    real (PREC), dimension(:,:), allocatable :: lcoefs
    real (PREC), dimension(:), allocatable :: lintercept
    real (PREC), dimension(:,:), allocatable :: scores, loadings
    real (PREC), dimension(:), allocatable :: sval, shift_x, scale_x, propvar
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Xp
    real (PREC), dimension(:), allocatable :: weights_sqrt
    integer :: Nrhs, Nconst_add, nobs, ncoefs, Nlhs, lncomp, Nconst_coefs
    logical :: trans_coefs
    real (PREC) :: var_rhs
    integer :: k, i
    character (*), parameter :: NAME = 'PCR_PCA'

    lstatus = NF_STATUS_OK

    nullify (ptr_Xp)

    if (present(conf)) lconf = conf

    ! --- Input checks ---

    call get_dims (X, Y, coefs, intercept, lconf%add_intercept, lconf%trans_x, &
        lconf%trans_y, lconf%trans_coefs, weights=weights, nobs=nobs, nrhs=Nrhs, &
        nlhs=nlhs, ncoefs=ncoefs, nconst_add=Nconst_add, nconst_coefs=Nconst_coefs, &
        status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (Nobs >= (Nrhs + Nconst_add), NAME, 'Too few observations', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (lconf%var_rhs_min <= 1.0_PREC, NAME, &
        'CONF: Invalid value for VAR_RHS_MIN', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (lconf%ncomp_min >= 0 .and. lconf%ncomp_min <= Nrhs, &
        NAME, 'CONF: Invalid value for NCOMP_MIN', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! NCOMP will be ignored if negative
    call check_cond (lconf%ncomp <= Nrhs, NAME, &
        'CONF: Invalid value for NCOMP', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (lconf%ncomp /= PCR_NCOMP_UNDEFINED) then
        call check_cond (lconf%ncomp >= 0 .and. lconf%ncomp <= Nrhs, NAME, &
            'CONF: Invalid value for NCOMP', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(weights)) then
        call check_cond (size(weights) == Nobs, NAME, &
            'WEIGHTS: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Run PCA ---

    ! Allocate working arrays
    ! Preliminary number of components: use all available, restrict below.
    lncomp = Nrhs
    allocate (scores(nobs, lncomp))
    allocate (loadings(Nrhs, lncomp))
    allocate (sval(lncomp))
    allocate (shift_x(Nrhs), scale_x(Nrhs))
    allocate (propvar(Nrhs), source=0.0_PREC)

    if (present(weights)) then
        allocate (weights_sqrt(Nobs))
        weights_sqrt(:) = weights / sum(weights)
        weights_sqrt(:) = sqrt(weights_sqrt)
    end if

    ! Skip PCA if no RHS are present
    if (Nrhs == 0) goto 50

    ! --- Transform X ---

    if (lconf%trans_x .or. lconf%center .or. lconf%scale .or. present(weights)) then
        allocate (ptr_Xp(Nobs,Nrhs))
        if (present(weights)) then
            call transform_regressors (X, ptr_Xp, trans=lconf%trans_x, &
                center=lconf%center, scale=lconf%scale, shift_x=shift_x, &
                scale_x=scale_x, drop_na=.false., weights=weights, status=lstatus)
        else
            call transform_regressors (X, ptr_Xp, trans=lconf%trans_x, &
                center=lconf%center, scale=lconf%scale, shift_x=shift_x, &
                scale_x=scale_x, drop_na=.false., status=lstatus)
        end if
        if (lstatus /= NF_STATUS_OK) goto 100

        if (present(weights)) then
            do i = 1, Nrhs
                ptr_Xp(:,i) = ptr_Xp(:,i) * weights_sqrt
            end do
        end if
    else
        ptr_Xp => X
        scale_x(:) = 1.0
        shift_x(:) = 0.0
    end if

    ! Perform principal component analysis
    call pca (ptr_Xp, center=.false., scale=.false., trans_x=.false., &
        scores=scores, sval=sval, loadings=loadings, propvar=propvar, status=lstatus)
    if (.not. (NF_STATUS_OK .in. lstatus)) goto 100

    ! Select number of principal components to use:
    ! 1. Select by min. variance, if applicable;
    if (lconf%var_rhs_min >= 0.0_PREC) then
        var_rhs = 0.0
        do k = 1, Nrhs
            var_rhs = var_rhs + propvar(k)
            if ((var_rhs - lconf%var_rhs_min) >= 0.0_PREC) exit
        end do
        lncomp = min(k, Nrhs)
    end if

    ! 2. User-imposed fixed number of components
    if (lconf%ncomp >= 0 .and. lconf%ncomp <= Nrhs) then
        lncomp = lconf%ncomp
    end if

    ! 3. Apply bounds
    if (lconf%ncomp_min > 0) then
        lncomp = max(lconf%ncomp_min, lncomp)
    end if

    if (lconf%ncomp_max >= 0) then
        lncomp = min(lconf%ncomp_max, lncomp)
    end if

    ! This should be satisfied, but still make sure it is.
    lncomp = min(lncomp, Nrhs)

    ! Jump to this point of PCA is skipped because of no RHS variables.
50  continue

    allocate (lcoefs(Nrhs,Nlhs))
    allocate (lintercept(Nlhs))
    ! Save user-provided TRANS_COEFS
    trans_coefs = lconf%trans_coefs
    lconf%trans_coefs = .false.

    ! Run principal component regression
    call pcr (Y, scores(:,1:lncomp), sval(1:lncomp), loadings(:,1:lncomp), &
        lcoefs, shift_x, scale_x, lconf%add_intercept, lconf%trans_y,  &
        lconf%trans_coefs, lconf%center, lintercept, weights_sqrt, lstatus)

    ! Restore user-provided TRANS_COEFS
    lconf%trans_coefs = trans_coefs

    if (.not. (NF_STATUS_OK .in. lstatus)) goto 100

    ! --- Store result object ---

    if (present(intercept)) then
        intercept = 0.0_PREC
        if (lconf%add_intercept) then
            intercept(1:Nlhs) = lintercept(1:Nlhs)
        end if
    end if

    if (present(coefs)) then
        if (lconf%trans_coefs) then
            do i = 1, Nlhs
                coefs(i,1+Nconst_coefs:Ncoefs) = lcoefs(1:Nrhs,i)
            end do
            if (Nconst_coefs == 1) then
                coefs(1:Nlhs,1) = lintercept
            end if
        else
            do i = 1, Nlhs
                coefs(1+Nconst_coefs:Ncoefs,i) = lcoefs(1:Nrhs,i)
            end do
            if (Nconst_coefs == 1) then
                coefs(1,1:Nlhs) = lintercept
            end if
        end if
    end if

    if (present(res)) then
        var_rhs = sum(propvar(1:lncomp))
        call lm_result_update (res, lconf, model=NF_STATS_LM_PCR, &
            nobs=nobs, nrhs=Nrhs, nlhs=Nlhs, ncomp=lncomp, var_rhs=var_rhs)

        if (allocated(res%coefs_multi)) then
            call copy_alloc (lcoefs, res%coefs_multi)
        else
            call move_alloc (lcoefs, res%coefs_multi)
        end if

        if (lconf%add_intercept) then
            if (allocated(res%intercept_multi)) then
                call copy_alloc (lintercept, res%intercept_multi)
            else
                call move_alloc (lintercept, res%intercept_multi)
            end if
        else
            call cond_alloc (res%intercept_multi, Nlhs)
            res%intercept_multi(:) = 0.0_PREC
        end if
    end if

    if (present(ncomp)) ncomp = lncomp

100 continue

    call assert_dealloc_ptr (X, ptr_Xp)

    if (present(status)) status = lstatus

end subroutine



subroutine pcr_pca_1d (conf, X, y, coefs, ncomp, intercept, weights, res, status)
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous, target :: y
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: coefs
    integer, intent(inout), optional :: ncomp
    real (PREC), intent(out), optional :: intercept
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
    type (lm_result), intent(inout), optional :: res
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer, contiguous :: ptr_y
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_coefs
    real (PREC), dimension(:), allocatable :: lintercept

    type (status_t) :: lstatus
    type (lm_config) :: lconf

    lstatus = NF_STATUS_OK

    nullify (ptr_coefs, ptr_y)

    if (present(conf)) lconf = conf
    lconf%trans_y = .false.
    lconf%trans_coefs = .false.

    ptr_y(1:size(y),1:1) => y

    if (present(coefs)) then
        ptr_coefs(1:size(coefs),1:1) => coefs
    end if

    if (present(intercept)) then
        allocate (lintercept(1))
    end if

    call pcr (lconf, X, ptr_y, ptr_coefs, ncomp, lintercept, weights, res, lstatus)

    if (lstatus == NF_STATUS_OK) then
        if (present(res)) then
            if (allocated(res%coefs_multi)) then
                call copy_alloc (res%coefs_multi(:,1), res%coefs)
            end if

            if (allocated(res%intercept_multi)) then
                res%intercept = res%intercept_multi(1)
            end if
        end if

        if (present(intercept)) then
            intercept = lintercept(1)
        end if
    end if

    if (present(status)) status = lstatus

end subroutine



subroutine pcr_pca_masked_2d (conf, X, Y, mask, coefs, ncomp, intercept, &
        weights, res, status)
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:,:), contiguous, target :: Y
    logical, intent(in), dimension(:,:), contiguous, target :: mask
        !*  Array identifying observations to be included in estimation. Only
        !   observations where MASK is .TRUE. will be included. Observations
        !   are permitted to vary across LHS variables.
    real (PREC), intent(out), dimension(:,:), contiguous, optional, target :: coefs
    integer, intent(out), optional :: ncomp
        !   If present, contains the actual number of PCs used on exit.
    real (PREC), intent(out), dimension(:), optional, contiguous, target :: intercept
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
    type (lm_result), intent(inout), optional, target :: res
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    type (lm_config) :: lconf, lconf_i
    type (lm_result) :: lres
    real (PREC), dimension(:), allocatable :: lcoefs
    real (PREC) :: lintercept
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Xp, ptr_Yp
    real (PREC), dimension(:,:), allocatable :: Xpack
    real (PREC), dimension(:), allocatable :: Ypack, Wpack
    logical, dimension(:,:), pointer, contiguous :: ptr_mask
    integer :: Nrhs, Nconst_add, nobs, ncoefs, Nlhs, Nconst_coefs
    integer :: lncomp, ncomp_i, Nobs_total, Nobs_i
    integer :: shp(2)
    real (PREC) :: Nobs_avg, var_rhs, var_rhs_avg
    integer :: i
    character (*), parameter :: NAME = 'PCR_PCA_MASKED'

    lstatus = NF_STATUS_OK

    nullify (ptr_Xp, ptr_Yp, ptr_mask)

    if (present(conf)) lconf = conf

    ! --- Input checks ---

    call get_dims (X, Y, coefs, intercept, lconf%add_intercept, lconf%trans_x, &
        lconf%trans_y, lconf%trans_coefs, weights=weights, mask=mask, &
        nobs=nobs, nrhs=Nrhs, nlhs=nlhs, ncoefs=ncoefs, nconst_add=Nconst_add, &
        nconst_coefs=Nconst_coefs, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (Nobs >= (Nrhs + Nconst_add), NAME, 'Too few observations', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (lconf%var_rhs_min <= 1.0_PREC, NAME, &
        'CONF: Invalid value for VAR_RHS_MIN', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (lconf%ncomp_min >= 0 .and. lconf%ncomp_min <= Nrhs, &
        NAME, 'CONF: Invalid value for NCOMP_MIN', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! NCOMP will be ignored if negative
    call check_cond (lconf%ncomp <= Nrhs, NAME, &
        'CONF: Invalid value for NCOMP', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (lconf%ncomp /= PCR_NCOMP_UNDEFINED) then
        call check_cond (lconf%ncomp >= 0 .and. lconf%ncomp <= Nrhs, NAME, &
            'CONF: Invalid value for NCOMP', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Transpose data into required shape ---

    if (.not. lconf%trans_x) then
        allocate (ptr_Xp(Nrhs,Nobs))
        call transpose_no_copy (X, ptr_Xp)
    else if (Nrhs == 0) then
        ! Work around problem that ASSOCIATED() returns false for
        ! zero-sized arrays
        allocate (ptr_Xp(0,Nobs))
    else
        ptr_Xp => X
    end if

    if (lconf%trans_y) then
        allocate (ptr_Yp(Nobs,Nlhs), ptr_mask(Nobs,Nlhs))
        call transpose_no_copy (Y, ptr_Yp)
        call transpose_no_copy (mask, ptr_mask)
    else
        ptr_Yp => Y
        ptr_mask => mask
    end if


    allocate (Xpack(Nrhs,Nobs))
    allocate (Ypack(Nobs))

    if (present(weights)) then
        allocate (Wpack(Nobs))
    end if

    ! Force transpose of X in actual PCA/PCR routine
    lconf%trans_x = .true.
    lconf%trans_y = .false.

    allocate (lcoefs(Ncoefs))
    lncomp = 0

    if (present(res)) then
        shp(1) = Nrhs
        shp(2) = Nlhs
        call cond_alloc (res%coefs_multi, shp)
        call cond_alloc (res%intercept_multi, Nlhs)

        res%coefs_multi(:,:) = 0.0
        res%intercept_multi(:) = 0.0
    end if

    Nobs_total = count (mask)
    var_rhs_avg = 0.0
    Nobs_avg = real(Nobs_total, PREC) / Nlhs

    do i = 1, Nlhs
        Nobs_i = count (ptr_mask(:,i))

        if (Nobs_i >= Ncoefs) then
            lconf_i = lconf
            lconf_i%ncomp = min(lconf%ncomp, max(0, Nobs_i - Nconst_coefs))

            call copy (ptr_Xp, Xpack(:,1:Nobs_i), ptr_mask(:,i), dim=2, status=lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100

            Ypack(1:Nobs_i) = pack (ptr_Yp(:,i), ptr_mask(:,i))
            if (present(weights)) then
                Wpack(1:Nobs_i) = pack (weights, ptr_mask(:,i))
                ! The nested fit validates and normalizes the subset weights.
            end if

            if (present(weights)) then
                call pcr (lconf, Xpack(:,1:Nobs_i), Ypack(1:Nobs_i), lcoefs, &
                    ncomp_i, lintercept, Wpack(1:Nobs_i), res=lres, status=lstatus)
            else
                call pcr (lconf, Xpack(:,1:Nobs_i), Ypack(1:Nobs_i), lcoefs, &
                    ncomp_i, lintercept, res=lres, status=lstatus)
            end if
            if (lstatus /= NF_STATUS_OK) goto 100

            var_rhs = lres%var_rhs
        else
            lcoefs(:) = 0.0
            lintercept = 0.0
            var_rhs = 0.0
            ncomp_i = 0
        end if

        lncomp = max(lncomp, ncomp_i)
        var_rhs_avg = var_rhs_avg + (Nobs_i / real(Nobs_total, PREC)) * var_rhs

        if (present(intercept)) then
            intercept(i) = lintercept
        end if

        if (present(coefs)) then
            if (lconf%trans_coefs) then
                coefs(i,1:Ncoefs) = lcoefs(1:Ncoefs)
            else
                coefs(1:Ncoefs,i) = lcoefs(1:Ncoefs)
            end if
        end if

    end do

    if (present(ncomp)) ncomp = lncomp

    if (present(res)) then
        call lm_result_update (res, lconf, model=NF_STATS_LM_PCR, &
            nobs=int(Nobs_avg), nrhs=Nrhs, ncomp=lncomp, var_rhs=var_rhs_avg)
    end if

100 continue

    call assert_dealloc_ptr (X, ptr_Xp)
    call assert_dealloc_ptr (Y, ptr_Yp)
    call assert_dealloc_ptr (mask, ptr_mask)

    if (present(status)) status = lstatus

end subroutine



subroutine pcr_pca_by_lhs (conf, X, Y, mask, weights, ncomp, coefs, intercept, &
        ncomp_used, res, status)
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:,:), contiguous, target :: Y
    logical, intent(in), dimension(:,:), contiguous, optional, target :: mask
        !*  Array identifying observations to be included in estimation. Only
        !   observations where MASK is .TRUE. will be included. Observations
        !   are permitted to vary across LHS variables.
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
    integer, intent(in), dimension(:), contiguous, optional :: ncomp
        !*  If present, contains the number of princial components to be used
        !   for each LHS variable.
    real (PREC), intent(out), dimension(:,:), contiguous, optional, target :: coefs
    real (PREC), intent(out), dimension(:), optional, contiguous, target :: intercept
    integer, intent(out), dimension(:), optional, contiguous :: ncomp_used
    type (lm_result), intent(inout), optional, target :: res
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    type (lm_config) :: lconf
    real (PREC), dimension(:,:), allocatable :: lcoefs
    real (PREC), dimension(:), allocatable :: lintercept
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Xp, ptr_Yp
    real (PREC), dimension(:,:), allocatable :: X_lhs_T, Xw_all, xp_lhs
    real (PREC), dimension(:), allocatable :: Y_lhs, W_lhs
    real (PREC), dimension(:), allocatable :: W_all
    real (PREC), dimension(:,:), allocatable :: scores, loadings
    real (PREC), dimension(:), allocatable :: sval, shift_x, scale_x
    real (PREC), dimension(:,:), allocatable :: scores_all, loadings_all
    real (PREC), dimension(:), allocatable :: sval_all, shift_x_all, scale_x_all
    logical, dimension(:,:), pointer, contiguous :: ptr_mask
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_scores, ptr_loadings
    real (PREC), dimension(:), pointer, contiguous :: ptr_sval, ptr_shift_x, ptr_scale_x
    real (PREC), dimension(:), pointer, contiguous :: ptr_W_lhs, ptr_Y_lhs
    integer :: Nrhs, Nconst_add, nobs, ncoefs, Nlhs, Nconst_coefs
    integer :: Nobs_lhs
    integer :: shp(2)
    real (PREC) :: Nobs_avg, Ncomp_avg, var_rhs_avg, var_rhs
    real (PREC) :: weights_lhs, weights_total
    logical :: has_all_obs, has_weights
    logical, dimension(:), allocatable :: lwork
    integer :: i, j, k
    integer (NF_ENUM_KIND) :: status_code
    character (*), parameter :: NAME = 'PCR_PCA_MASKED'

    target :: scores_all, loadings_all, sval_all, shift_x_all, scale_x_all
    target :: scores, loadings, sval, shift_x, scale_x
    target :: Y_lhs, W_lhs, W_all, Xw_all

    lstatus = NF_STATUS_OK

    nullify (ptr_Xp, ptr_Yp, ptr_mask)

    if (present(conf)) lconf = conf

    has_weights = present(weights)

    ! --- Input checks ---

    call get_dims (X, Y, coefs, intercept, lconf%add_intercept, lconf%trans_x, &
        lconf%trans_y, lconf%trans_coefs, weights=weights, mask=mask, &
        nobs=nobs, nrhs=Nrhs, nlhs=nlhs, ncoefs=ncoefs, nconst_add=Nconst_add, &
        nconst_coefs=Nconst_coefs, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (Nobs >= (Nrhs + Nconst_add), NAME, 'Too few observations', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (lconf%var_rhs_min <= 1.0_PREC, NAME, &
        'CONF: Invalid value for VAR_RHS_MIN', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (lconf%ncomp_min >= 0 .and. lconf%ncomp_min <= Nrhs, &
        NAME, 'CONF: Invalid value for NCOMP_MIN', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(ncomp)) then
        call check_cond (size(ncomp) == Nlhs, NAME, 'NCOMP: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        call check_cond (all((ncomp >= 0 .and. ncomp <= Nrhs) &
            .or. ncomp == PCR_NCOMP_UNDEFINED), NAME, &
            'NCOMP: Invalid value encountered', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        if (lconf%ncomp /= PCR_NCOMP_UNDEFINED) then
            call check_cond (lconf%ncomp >= 0 .and. lconf%ncomp <= Nrhs, NAME, &
                'CONF: Invalid value for NCOMP', lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100
        end if
    end if

    if (present(ncomp_used)) then
        call check_cond (size(ncomp_used) == Nlhs, NAME, &
            'NCOMP_USED: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Precompute PCA if all obs. are to be used ---

    has_all_obs = .not. present(mask)
    if (present(mask)) then
        allocate (lwork(Nlhs))
        if (lconf%trans_y) then
            lwork(:) = all(mask, dim=2)
        else
            lwork(:) = all(mask, dim=1)
        end if
        has_all_obs = any (lwork)
        deallocate (lwork)
    end if

    if (has_all_obs) then
        allocate (scores_all(Nobs,Nrhs))
        allocate (loadings_all(Nrhs,Nrhs))
        allocate (sval_all(Nrhs))
        allocate (shift_x_all(Nrhs), scale_x_all(Nrhs))

        allocate (Xw_all(Nobs,Nrhs))
        if (has_weights) then
            weights_total = sum(weights)
            allocate (W_all, source=weights)
            W_all(:) = W_all / weights_total
            W_all(:) = sqrt(W_all)

            call transform_regressors (X, Xw_all, trans=lconf%trans_x, &
                center=lconf%center, scale=lconf%scale, shift_x=shift_x_all, &
                scale_x=scale_x_all, drop_na=.false., weights=weights, status=lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100

            do i = 1, Nrhs
                Xw_all(:,i) = Xw_all(:,i) * W_all
            end do
        else
            call transform_regressors (X, Xw_all, trans=lconf%trans_x, &
                center=lconf%center, scale=lconf%scale, shift_x=shift_x_all, &
                scale_x=scale_x_all, drop_na=.false., status=lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100
        end if

        call pca (Xw_all, center=.false., scale=.false., &
            trans_x=.false., scores=scores_all, sval=sval_all, &
            loadings=loadings_all, status=lstatus)

        deallocate (Xw_all)

        if (.not. (NF_STATUS_OK .in. lstatus)) then
            status_code = lstatus
            goto 100
        end if
    end if

    ! --- Transpose data into required shape ---

    if (.not. lconf%trans_x) then
        allocate (ptr_Xp(Nrhs,Nobs))
        call transpose_no_copy (X, ptr_Xp)
    else if (Nrhs == 0) then
        ! Work around problem that ASSOCIATED() returns false for
        ! zero-sized arrays
        allocate (ptr_Xp(0,Nobs))
    else
        ptr_Xp => X
    end if

    if (lconf%trans_y) then
        allocate (ptr_Yp(Nobs,Nlhs), ptr_mask(Nobs,Nlhs))
        call transpose_no_copy (Y, ptr_Yp)
        if (present(mask)) then
            call transpose_no_copy (mask, ptr_mask)
        end if
    else
        ptr_Yp => Y
        if (present(mask)) then
            ptr_mask => mask
        end if
    end if

    ! --- Process each LHS variable ---

    Nobs_avg = 0.0
    Ncomp_avg = 0.0
    var_rhs_avg = 0.0
    weights_total = 0.0

    if (present(ncomp_used)) then
        ! Default value for LHS variables that will be skipped due to low
        ! obs. count
        ncomp_used = PCR_NCOMP_UNDEFINED
    end if

    allocate (lcoefs(Nrhs,Nlhs), lintercept(Nlhs))

    allocate (X_lhs_T(Nrhs,Nobs))
    allocate (Y_lhs(Nobs))

    if (present(weights)) then
        allocate (W_lhs(Nobs))
    end if

    allocate (loadings(Nrhs,Nrhs))
    allocate (shift_x(Nrhs), scale_x(Nrhs))
    allocate (sval(Nrhs))

    nullify (ptr_Y_lhs, ptr_W_lhs)
    nullify (ptr_scores, ptr_loadings, ptr_sval, ptr_shift_x, ptr_scale_x)

    ! Initialize here to avoid gfortran warnings
    status_code = NF_STATUS_OK

    loop_lhs: do j = 1, Nlhs

        ! Skip remaining tasks if we previously encountered an error
        if (status_code /= NF_STATUS_OK) cycle

        if (present(mask)) then
            Nobs_lhs = count (ptr_mask(:,j))
        else
            Nobs_lhs = Nobs
        end if

        if (Nobs_lhs <= Ncoefs) then
            lcoefs(:,j) = 0.0
            lintercept(j) = 0.0
            cycle loop_lhs
        end if

        ! --- Create contiguous input data ---

        weights_lhs = Nobs_lhs

        if (Nobs_lhs /= Nobs) then

            ! Create contiguous subset of RHS variables
            call copy (ptr_Xp, X_lhs_T(:,1:Nobs_lhs), ptr_mask(:,j), dim=2, status=lstatus)
            if (.not. (NF_STATUS_OK .in. lstatus)) then
                status_code = lstatus
                cycle loop_lhs
            end if

            ! Create contiguous array of LHS obs.
            call copy (ptr_Yp(:,j), Y_lhs(1:Nobs_lhs), ptr_mask(:,j))
            ptr_Y_lhs => Y_lhs(1:Nobs_lhs)

            ! Keep raw subset weights until weighted centering and scaling are done.
            nullify (ptr_W_lhs)
            if (has_weights) then
                call copy (weights, W_lhs(1:Nobs_lhs), ptr_mask(:,j))
                weights_lhs = sum(W_lhs(1:Nobs_lhs))
                ptr_W_lhs => W_lhs(1:Nobs_lhs)
            end if

            call cond_alloc (xp_lhs, [Nobs_lhs, Nrhs])
            call transform_regressors (X_lhs_T(:,1:Nobs_lhs), xp_lhs, &
                trans=.true., center=lconf%center, scale=lconf%scale, &
                shift_x=shift_x, scale_x=scale_x, drop_na=.false., &
                weights=ptr_W_lhs, status=lstatus)
            if (lstatus /= NF_STATUS_OK) then
                status_code = lstatus
                cycle loop_lhs
            end if

            if (has_weights) then
                ! Apply weights after computing moments in the original units.
                W_lhs(1:Nobs_lhs) = sqrt(W_lhs(1:Nobs_lhs) / weights_lhs)
                do i = 1, Nrhs
                    xp_lhs(:,i) = xp_lhs(:,i) * W_lhs(1:Nobs_lhs)
                end do
            end if

            ! --- Run PCA ---

            shp(1) = Nobs_lhs
            shp(2) = Nrhs
            call cond_alloc (scores, shp(1:2))

            ! PCA must not center or scale the already weighted matrix again.
            call pca (xp_lhs, center=.false., scale=.false., &
                trans_x=.false., scores=scores, sval=sval, &
                loadings=loadings, status=lstatus)
            if (.not. (NF_STATUS_OK .in. lstatus)) then
                status_code = lstatus
                cycle loop_lhs
            end if

            ptr_scores => scores
            ptr_loadings => loadings
            ptr_sval => sval
            ptr_shift_x => shift_x
            ptr_scale_x => scale_x
        else
            ! Use precomputed PCA data derived from all RHS observations
            ptr_Y_lhs => ptr_Yp(:,j)
            nullify (ptr_W_lhs)
            if (has_weights) ptr_W_lhs => W_all

            ptr_scores => scores_all
            ptr_loadings => loadings_all
            ptr_sval => sval_all
            ptr_shift_x => shift_x_all
            ptr_scale_x => scale_x_all

            if (present(weights)) then
                weights_lhs = sum(weights)
            end if
        end if

        ! --- Select number of PCs ---

        k = 0

        ! 1. Select by min variance, if applicable
        if (lconf%var_rhs_min >= 0.0_PREC) then
            var_rhs = 0.0
            do k = 1, Nrhs
                var_rhs = var_rhs + ptr_sval(k)**2.0_PREC / sum(ptr_sval**2.0_PREC)
                if ((var_rhs - lconf%var_rhs_min) >= 0.0_PREC) exit
            end do
            k = min(k, Nrhs)
        end if

        ! 2. User-imposed overall number of components
        if (lconf%ncomp >= 0 .and. lconf%ncomp <= Nrhs) then
            k = lconf%ncomp
        end if

        ! 3. Run PCR on using given PCA data
        if (present(ncomp)) then
            if (ncomp(j) /= PCR_NCOMP_UNDEFINED) k = ncomp(j)
        end if

        ! 4. User-imposed lower bound
        if (lconf%ncomp_min > 0) then
            k = max(lconf%ncomp_min, k)
        end if

        ! 5. User-imposed upper bound
        if (lconf%ncomp_max >= 0) then
            k = min(lconf%ncomp_max, k)
        end if

        ! Center Y consistently with the X preprocessing above.

        call pcr (ptr_Y_lhs, ptr_scores(:,1:k), ptr_sval(1:k), &
            ptr_loadings(:,1:k), lcoefs(:,j), shift_x=ptr_shift_x, &
            scale_x=ptr_scale_x, add_intercept=lconf%add_intercept, &
            center=lconf%center, &
            intercept=lintercept(j), weights_sqrt=ptr_W_lhs, status=lstatus)
        if (lstatus /= NF_STATUS_OK) then
            status_code = lstatus
            cycle loop_lhs
        end if

        ! Fraction of RHS variance captured by components used
        var_rhs = sum(ptr_sval(1:k)**2.0_PREC) / sum(ptr_sval**2.0_PREC)

        weights_total = weights_total + weights_lhs
        Nobs_avg = Nobs_avg + Nobs_lhs * weights_lhs
        Ncomp_avg = Ncomp_avg + k * weights_lhs
        var_rhs_avg = var_rhs_avg + var_rhs * weights_lhs

        if (present(ncomp_used)) then
            ncomp_used(j) = k
        end if

    end do loop_lhs

    ! Do not publish results or average weights after a failed response fit.
    lstatus = status_code
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(intercept)) then
        intercept(1:Nlhs) = lintercept(1:Nlhs)
    end if

    if (present(coefs)) then
        do j = 1, Nlhs
            if (lconf%trans_coefs) then
                coefs(j,1+Nconst_coefs:Ncoefs) = lcoefs(1:Nrhs,j)
                if (Nconst_coefs == 1) then
                    coefs(j,1) = lintercept(j)
                end if
            else
                coefs(1+Nconst_coefs:Ncoefs,j) = lcoefs(1:Nrhs,j)
                if (Nconst_coefs == 1) then
                    coefs(1,j) = lintercept(j)
                end if
            end if
        end do
    end if

    if (present(res)) then
        shp(1) = Nrhs
        shp(2) = Nlhs
        call cond_alloc (res%coefs_multi, shp)
        call cond_alloc (res%intercept_multi, Nlhs)

        var_rhs_avg = var_rhs_avg / weights_total
        Ncomp_avg = Ncomp_avg / weights_total
        Nobs_avg = Nobs_avg / weights_total

        call copy_alloc (lcoefs, res%coefs_multi)
        call copy_alloc (lintercept, res%intercept_multi)

        call lm_result_update (res, lconf, model=NF_STATS_LM_PCR, &
            nobs=int(Nobs_avg), nrhs=Nrhs, ncomp=int(Ncomp_avg), &
            var_rhs=var_rhs_avg)
    end if

100 continue

    call assert_dealloc_ptr (X, ptr_Xp)
    call assert_dealloc_ptr (Y, ptr_Yp)
    call assert_dealloc_ptr (mask, ptr_mask)

    if (present(status)) status = lstatus

end subroutine



subroutine pcr_cv_pooled_2d (conf, X, Y, ncomp, weights, mask, rmse, status)
    !*  PCR_CV implements cross-validation for principal component regression,
    !   selecting the number of principal components which minimize the
    !   average mean squared error computed on test sub-samples.
    type (lm_config), intent(in), optional :: conf
        !*  Configuration object
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
        !*  Predictor variables (regressors)
    real (PREC), intent(in), dimension(:,:), contiguous, target :: Y
        !*  Response variables (outcomes)
    integer, intent(out) :: ncomp
        !   Contains the cross-validated number of PCs used on exit.
    logical, intent(in), dimension(:,:), contiguous, optional, target :: mask
        !*  If present, only include observations where MASK is .TRUE. when
        !   estimating projection coefficients and evaluating goodness of fit.
    real (PREC), intent(in), dimension(:), contiguous, optional, target :: weights
    real (PREC), intent(out), optional :: rmse
        !*  Root MSE of the MSE-minimizing model
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: Nrhs, Nobs, Nlhs, Nconst_add, Nfolds
    integer, dimension(:), allocatable :: folds_ifrom, folds_size
    real (PREC), dimension(:,:,:), allocatable :: mse_ncomp_lhs, mse_weights_lhs
    real (PREC), dimension(:), allocatable :: mean_mse, se_mse, resid
    real (PREC), dimension(:,:), allocatable :: mse_ncomp, mse_weights
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_X, ptr_Y
    real (PREC), dimension(:), pointer, contiguous :: ptr_weights
    logical, dimension(:,:), pointer, contiguous :: ptr_mask
    integer, dimension(:), allocatable :: iorder
    integer :: i, imin, dim
    real (PREC) :: min_mse, ubound, std_resid, sum_w
    type (status_t) :: lstatus
    type (lm_config) :: lconf
    character (*), parameter :: NAME = 'PCR_CV'

    nullify (ptr_X, ptr_Y, ptr_weights, ptr_mask)

    ! --- Input processing ---

    lconf = conf

    call get_dims (X, Y, mask=mask, add_intercept=lconf%add_intercept, &
        trans_x=lconf%trans_x, trans_y=lconf%trans_y, nobs=nobs, nrhs=Nrhs, &
        nlhs=nlhs, nconst_add=Nconst_add, weights=weights, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(weights)) then
        call check_cond (size(weights) == Nobs, NAME, &
            'WEIGHTS: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Randomly reorder X and Y ---

    if (lconf%cv_shuffle_obs) then
        allocate (iorder(Nobs))
        call random_order (iorder)

        allocate (ptr_X(Nobs,Nrhs))
        allocate (ptr_Y(Nobs,Nlhs))

        dim = 1
        if (lconf%trans_x) dim = 2
        call pack_indexed (X, iorder, ptr_X, dim, lconf%trans_x, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        dim = 1
        if (lconf%trans_y) dim = 2
        call pack_indexed (Y, iorder, ptr_Y, dim, lconf%trans_y, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        if (present(weights)) then
            allocate (ptr_weights(Nobs))
            call pack_indexed (weights, iorder, ptr_weights)
        end if

        if (present(mask)) then
            dim = 1
            if (lconf%trans_y) dim = 2
            allocate (ptr_mask(Nobs,Nlhs))
            call pack_indexed (mask, iorder, ptr_mask, dim, lconf%trans_y, lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100
        end if

        lconf%trans_x = .false.
        lconf%trans_y = .false.
    else
        ptr_X => X
        ptr_Y => Y
        if (present(weights)) then
            ptr_weights => weights
        end if
        if (present(mask)) then
            ptr_mask => mask
        end if
    end if

    ! --- k-fold CV chunks ---

    Nfolds = lconf%cv_n
    allocate (folds_ifrom(Nfolds), folds_size(Nfolds))

    call split_uniform (Nobs, Nfolds, folds_ifrom, folds_size, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Perform cross-validation ---

    ! We have the following parallelization scenarios:
    !   1.  No masked obs: PCA is computed once per fold, so we parallelize
    !       over folds.
    !   2.  Masked obs: PCA needs to be done for each LHS variance, so if
    !       2a. NLHS > NFOLDS: parallelization should be over LHS
    !       2b. NFOLDS > NLHS: parallelization should be over folds

    ! Hold MSEs for number of components from 0-Nrhs, for each LHS and each fold
    allocate (mse_ncomp_lhs(0:Nrhs,Nlhs,Nfolds))
    allocate (mse_weights_lhs(0:Nrhs,Nlhs,Nfolds))

    if (.not. present(mask)) then
        call pcr_cv_dispatch_folds_omp (lconf, ptr_X, ptr_Y, folds_ifrom, &
            folds_size, mse_ncomp_lhs, mse_weights_lhs, lstatus, ptr_weights)
    else if (Nlhs > Nfolds) then
        ! Case 2a: parallelize over LHS
        call pcr_cv_dispatch_folds (lconf, ptr_X, ptr_Y, folds_ifrom, &
            folds_size, mse_ncomp_lhs, mse_weights_lhs, lstatus, ptr_weights, ptr_mask)
    else
        ! Case 2b, parallelize over folds
        call pcr_cv_dispatch_folds_omp (lconf, ptr_X, ptr_Y, folds_ifrom, &
            folds_size, mse_ncomp_lhs, mse_weights_lhs, lstatus, ptr_weights, ptr_mask)
    end if

    if (lstatus /= NF_STATUS_OK) goto 100

    ! Sum across LHS to create MSE and sum of weights for each #PC/fold
    allocate (mse_ncomp(0:Nrhs,Nfolds), mse_weights(0:Nrhs,Nfolds))
    do i = 1, Nfolds
        mse_ncomp(:,i) = sum(mse_ncomp_lhs(:,:,i), dim=2)
        mse_weights(:,i) = sum(mse_weights_lhs(:,:,i), dim=2)
    end do

    ! Compute weighted average MSE for each # of PCs:
    ! Each fold h returns, for each number of PCs k,
    !   MSE(k,h) = sum_{j} sum_{i in I_{jh}} (w_ij * resid_ijk^2)
    ! where j in 1,...,NLHS and i in 1,...,Nobs(j,h) is the observation index
    ! for LHS variable j and fold h.
    ! w_ij is the weight attached to the i-th residual of the j-th LHS variable.
    !
    ! For each fold h and each PC k, the sum of weights across all LHS and
    ! observations is returned in MSE_WEIGHTS(k,h)
    !
    ! Thus the average MSE across *all* folds is for #PC = k is given by
    !   sum_h MSE(k,h) / sum_h MSE_WEIGHTS(k,h)
    ! Note that the sum of weights should actually be identical for all k,
    ! since we have that Nobs >= Nrhs >= #PC.

    allocate (mean_mse(0:Nrhs), se_mse(0:Nrhs))
    allocate (resid(Nfolds))

    mean_mse(:) = sum(mse_ncomp, dim=2) / sum(mse_weights, dim=2)

    do i = 0, Nrhs
        sum_w = sum(mse_weights(i,:))
        ! "Residual" MSE for each fold around mean MSE for #PC = k
        resid(:) = mse_ncomp(i,:)/mse_weights(i,:) - mean_mse(i)
        ! Weight by relative fold size
        resid(:) = resid * sqrt(mse_weights(i,:)/sum_w)
        std_resid = sqrt(sum(resid**2.0_PREC))
        se_mse(i) = std_resid / sqrt(real(Nfolds,PREC))
    end do

    min_mse = huge(0.0_PREC)
    imin = 0

    do i = 0, Nrhs
        if (mean_mse(i) < min_mse) then
            imin = i
            min_mse = mean_mse(i)
        end if
    end do

    ! Pick most parsimonious model within one standard deviation of
    ! the minimum.
    if (lconf%cv_parsimonious) then
        ubound = min_mse + se_mse(imin) * max(lconf%cv_se_mult, 0.0_PREC)
        do i = 1, imin-1
            if (mean_mse(i) < ubound) then
                imin = i
            end if
        end do
    end if

    ncomp = imin
    if (present(rmse)) rmse = sqrt(min_mse)

100 continue

    call assert_dealloc_ptr (X, ptr_X)
    call assert_dealloc_ptr (Y, ptr_Y)
    call assert_dealloc_ptr (weights, ptr_weights)
    call assert_dealloc_ptr (mask, ptr_mask)

    if (present(status)) status = lstatus

end subroutine



subroutine pcr_cv_separate_2d (conf, X, Y, ncomp, weights, mask, rmse, status)
    !*  PCR_CV implements cross-validation for principal component regression,
    !   selecting the number of principal components which minimize the
    !   average mean squared error computed on test sub-samples.
    type (lm_config), intent(in), optional :: conf
        !*  Configuration object
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
        !*  Predictor variables (regressors)
    real (PREC), intent(in), dimension(:,:), contiguous, target :: Y
        !*  Response variables (outcomes)
    integer, intent(out), dimension(:), contiguous :: ncomp
        !   Contains the cross-validated number of PCs for each LHS variable.
    logical, intent(in), dimension(:,:), contiguous, optional, target :: mask
        !*  If present, only include observations where MASK is .TRUE. when
        !   estimating projection coefficients and evaluating goodness of fit.
    real (PREC), intent(in), dimension(:), contiguous, optional, target :: weights
    real (PREC), intent(out), dimension(:), contiguous, optional :: rmse
        !*  Root MSE of the MSE-minimizing model
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    integer :: Nrhs, Nobs, Nlhs, Nconst_add, Nfolds
    integer, dimension(:), allocatable :: folds_ifrom, folds_size
    real (PREC), dimension(:,:,:), allocatable :: mse_ncomp_lhs, mse_weights_lhs
    real (PREC), dimension(:), allocatable :: mean_mse, se_mse, resid
    real (PREC), dimension(:,:), allocatable :: mse_ncomp, mse_weights
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_X, ptr_Y
    real (PREC), dimension(:), pointer, contiguous :: ptr_weights
    logical, dimension(:,:), pointer, contiguous :: ptr_mask
    integer, dimension(:), allocatable :: iorder
    integer :: i, j, imin, dim
    real (PREC) :: min_mse, ubound, std_resid, sum_w, Neff
    type (status_t) :: lstatus
    type (lm_config) :: lconf
    character (*), parameter :: NAME = 'PCR_CV'

    nullify (ptr_X, ptr_Y, ptr_weights, ptr_mask)

    ! --- Input processing ---

    lconf = conf

    call get_dims (X, Y, mask=mask, add_intercept=lconf%add_intercept, &
        trans_x=lconf%trans_x, trans_y=lconf%trans_y, nobs=nobs, nrhs=Nrhs, &
        nlhs=nlhs, nconst_add=Nconst_add, weights=weights, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(ncomp) == Nlhs, NAME, 'NCOMP: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(weights)) then
        call check_cond (size(weights) == Nobs, NAME, &
            'WEIGHTS: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Randomly reorder X and Y ---

    if (lconf%cv_shuffle_obs) then
        allocate (iorder(Nobs))
        call random_order (iorder)

        allocate (ptr_X(Nobs,Nrhs))
        allocate (ptr_Y(Nobs,Nlhs))

        dim = 1
        if (lconf%trans_x) dim = 2
        call pack_indexed (X, iorder, ptr_X, dim, lconf%trans_x, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        dim = 1
        if (lconf%trans_y) dim = 2
        call pack_indexed (Y, iorder, ptr_Y, dim, lconf%trans_y, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        if (present(weights)) then
            allocate (ptr_weights(Nobs))
            call pack_indexed (weights, iorder, ptr_weights)
        end if

        if (present(mask)) then
            dim = 1
            if (lconf%trans_y) dim = 2
            allocate (ptr_mask(Nobs,Nlhs))
            call pack_indexed (mask, iorder, ptr_mask, dim, lconf%trans_y, lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100
        end if

        lconf%trans_x = .false.
        lconf%trans_y = .false.
    else
        ptr_X => X
        ptr_Y => Y
        if (present(weights)) then
            ptr_weights => weights
        end if
        if (present(mask)) then
            ptr_mask => mask
        end if
    end if

    ! --- k-fold CV chunks ---

    Nfolds = lconf%cv_n
    allocate (folds_ifrom(Nfolds), folds_size(Nfolds))

    call split_uniform (Nobs, Nfolds, folds_ifrom, folds_size, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Perform cross-validation ---

    ! We have the following parallelization scenarios:
    !   1.  No masked obs: PCA is computed once per fold, so we parallelize
    !       over folds.
    !   2.  Masked obs: PCA needs to be done for each LHS variance, so if
    !       2a. NLHS > NFOLDS: parallelization should be over LHS
    !       2b. NFOLDS > NLHS: parallelization should be over folds

    ! Hold MSEs for number of components from 0-Nrhs, for each LHS and each fold
    allocate (mse_ncomp_lhs(0:Nrhs,Nlhs,Nfolds))
    allocate (mse_weights_lhs(0:Nrhs,Nlhs,Nfolds))

    if (.not. present(mask)) then
        call pcr_cv_dispatch_folds_omp (lconf, ptr_X, ptr_Y, folds_ifrom, &
            folds_size, mse_ncomp_lhs, mse_weights_lhs, lstatus, ptr_weights)
    else if (Nlhs > Nfolds) then
        ! Case 2a: parallelize over LHS
        call pcr_cv_dispatch_folds (lconf, ptr_X, ptr_Y, folds_ifrom, &
            folds_size, mse_ncomp_lhs, mse_weights_lhs, lstatus, ptr_weights, ptr_mask)
    else
        ! Case 2b, parallelize over folds
        call pcr_cv_dispatch_folds_omp (lconf, ptr_X, ptr_Y, folds_ifrom, &
            folds_size, mse_ncomp_lhs, mse_weights_lhs, lstatus, ptr_weights, ptr_mask)
    end if

    if (lstatus /= NF_STATUS_OK) goto 100

    ! Compute weighted average MSE for each # of PCs:
    ! Each fold h returns, for each number of PCs k,
    !   MSE(k,h) = sum_{j} sum_{i in I_{jh}} (w_ij * resid_ijk^2)
    ! where j in 1,...,NLHS and i in 1,...,Nobs(j,h) is the observation index
    ! for LHS variable j and fold h.
    ! w_ij is the weight attached to the i-th residual of the j-th LHS variable.
    !
    ! For each fold h and each PC k, the sum of weights across all LHS and
    ! observations is returned in MSE_WEIGHTS(k,h)
    !
    ! Thus the average MSE across *all* folds is for #PC = k is given by
    !   sum_h MSE(k,h) / sum_h MSE_WEIGHTS(k,h)
    ! Note that the sum of weights should actually be identical for all k,
    ! since we have that Nobs >= Nrhs >= #PC.

    ! Sum across LHS to create MSE and sum of weights for each #PC/fold
    allocate (mse_ncomp(0:Nrhs,Nfolds), mse_weights(0:Nrhs,Nfolds))

    allocate (mean_mse(0:Nrhs), se_mse(0:Nrhs))
    allocate (resid(Nfolds))

    do j = 1, Nlhs
        mse_ncomp(:,:) = mse_ncomp_lhs(:,j,:)
        mse_weights(:,:) = mse_weights_lhs(:,j,:)

        do i = 0, Nrhs
            sum_w = sum(mse_weights(i,:))

            if (sum_w == 0.0_PREC) then
                mean_mse(i) = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
                se_mse(i) = 0.0_PREC
                cycle
            end if

            mean_mse(i) = sum(mse_ncomp(i,:)) / sum_w

            ! "Residual" MSE for each fold around mean MSE for #PC = k
            where (mse_weights(i,:) > 0)
                resid = mse_ncomp(i,:)/mse_weights(i,:) - mean_mse(i)
            else where
                resid = 0.0
            end where
            ! Weight by relative fold size
            resid(:) = resid * sqrt(mse_weights(i,:)/sum_w)
            std_resid = sqrt(sum(resid**2.0_PREC))

            ! Effective number of folds
            Neff = sum(mse_weights(i,:))**2.0 / sum(mse_weights(i,:)**2.0)

            se_mse(i) = std_resid / sqrt(real(Neff,PREC))
        end do

        if (.not. any(ieee_is_finite(mean_mse))) then
            ! Quick exit if there we don't have a single valid MSE
            ncomp(j) = PCR_NCOMP_UNDEFINED
            if (present(rmse)) rmse(j) = huge(0.0_PREC)
            cycle
        end if

        min_mse = huge(0.0_PREC)
        imin = 0

        do i = 0, Nrhs
            if (mean_mse(i) < min_mse) then
                imin = i
                min_mse = mean_mse(i)
            end if
        end do

        ! Pick most parsimonious model within one standard deviation of
        ! the minimum.
        if (lconf%cv_parsimonious) then
            ubound = min_mse + se_mse(imin) * max(lconf%cv_se_mult, 0.0_PREC)
            do i = 0, imin - 1
                if (mean_mse(i) < ubound) then
                    imin = i
                    exit
                end if
            end do
        end if

        ncomp(j) = imin
        if (present(rmse)) rmse(j) = sqrt(min_mse)
    end do

100 continue

    call assert_dealloc_ptr (X, ptr_X)
    call assert_dealloc_ptr (Y, ptr_Y)
    call assert_dealloc_ptr (weights, ptr_weights)
    call assert_dealloc_ptr (mask, ptr_mask)

    if (present(status)) status = lstatus

end subroutine



subroutine pcr_cv_1d (conf, X, Y, ncomp, rmse, weights, status)
    !*  PCR_CV implements cross-validation for principal component regression,
    !   selecting the number of principal components which minimize the
    !   average mean squared error computed on test sub-samples.
    !
    !   Routine implements a wrapper which accepts a single outcome variable.
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous, target :: Y
    integer, intent(out), optional :: ncomp
        !*  Contains the cross-validated number of PCs used on exit.
    real (PREC), intent(out), optional :: rmse
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer, contiguous :: ptr_Y
    type (lm_config) :: lconf

    lconf = conf
    lconf%trans_y = .false.

    ptr_Y(1:size(Y),1:1) => Y

    call pcr_cv (lconf, X, ptr_Y, ncomp, weights=weights, rmse=rmse, status=status)

end subroutine



subroutine pcr_cv_dispatch_folds_omp (conf, X, Y, folds_ifrom, folds_size, &
        mse_ncomp, mse_weights, status, weights, mask)
    type (lm_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X, Y
    integer, intent(in), dimension(:), contiguous :: folds_ifrom, folds_size
    real (PREC), intent(out), dimension(:,:,:), contiguous :: mse_ncomp
    real (PREC), intent(out), dimension(:,:,:), contiguous :: mse_weights
    type (status_t), intent(out) :: status
    real (PREC), intent(in), dimension(:), contiguous, optional, target :: weights
    logical, intent(in), dimension(:,:), contiguous, optional :: mask

    integer :: status_code, fold_status, i, ifrom, Ntest, Nfolds
    logical :: has_mask
    type (status_t) :: lstatus
    real (PREC), dimension(:), contiguous, pointer :: ptr_weights

    nullify (ptr_weights)
    if (present(weights)) ptr_weights => weights

    has_mask = present(mask)
    Nfolds = size(folds_ifrom)

    ! Workaround: OpenMP in gfortran 7.x cannot use user-defined
    ! operators, keep track using integer instead.
    status_code = NF_STATUS_OK

    !$omp parallel default(none) &
    !$omp shared(Nfolds,folds_ifrom,folds_size,conf,X,Y,mse_ncomp) &
    !$omp shared(ptr_weights,mask,has_mask,mse_weights) &
    !$omp private(i,ifrom,Ntest,lstatus,fold_status) &
    !$omp reduction(ior: status_code)

    ! Initalize again to avoid gfortran warnings
    status_code = NF_STATUS_OK

    !$omp do schedule (auto)
    do i = 1, Nfolds
        ifrom = folds_ifrom(i)
        Ntest = folds_size(i)
        if (has_mask) then
            call pcr_cv_mse_mask (conf, X, Y, mask, ifrom, Ntest, mse_ncomp(:,:,i), &
                mse_weights(:,:,i), weights=ptr_weights, status=lstatus)
        else
            call pcr_cv_mse (conf, X, Y, ifrom, Ntest, mse_ncomp(:,:,i), &
                mse_weights(:,:,i), weights=ptr_weights, status=lstatus)
        end if

        ! Retain failures from earlier folds handled by this thread.
        fold_status = lstatus
        status_code = ior(status_code, fold_status)
    end do
    !$omp end do
    !$omp end parallel

    ! Convert back from integer to STATUS_T
    status = status_code

end subroutine



subroutine pcr_cv_dispatch_folds (conf, X, Y, folds_ifrom, folds_size, &
        mse_ncomp, mse_weights, status, weights, mask)
    type (lm_config), intent(in) :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X, Y
    integer, intent(in), dimension(:), contiguous :: folds_ifrom, folds_size
    real (PREC), intent(out), dimension(:,:,:), contiguous :: mse_ncomp
    real (PREC), intent(out), dimension(:,:,:), contiguous :: mse_weights
    type (status_t), intent(out) :: status
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
    logical, intent(in), dimension(:,:), contiguous, optional :: mask

    integer :: i, ifrom, Ntest, Nfolds
    logical :: has_mask
    type (status_t) :: lstatus

    has_mask = present(mask)
    Nfolds = size(folds_ifrom)

    status = NF_STATUS_OK
    do i = 1, Nfolds
        ifrom = folds_ifrom(i)
        Ntest = folds_size(i)
        if (has_mask) then
            call pcr_cv_mse_mask (conf, X, Y, mask, ifrom, Ntest, mse_ncomp(:,:,i), &
                mse_weights(:,:,i), weights=weights, status=lstatus)
        else
            call pcr_cv_mse (conf, X, Y, ifrom, Ntest, mse_ncomp(:,:,i), &
                mse_weights(:,:,i), weights=weights, status=lstatus)
        end if

        ! Convert to integer to which IOR reduction can be applied
        status = status + lstatus
    end do

end subroutine



subroutine pcr_cv_mse (conf, X, Y, test_ifrom, Ntest, mse_ncomp, mse_weights, &
        weights, status)
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:,:), contiguous :: Y
    integer, intent(in) :: test_ifrom
    integer, intent(in) :: Ntest
    real (PREC), intent(out), dimension(:,:), contiguous :: mse_ncomp
    real (PREC), intent(out), dimension(:,:), contiguous :: mse_weights
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
    type (status_t), intent(out), optional :: status


    integer :: Ntrain, Nobs, Nrhs, Nlhs, Nconst_add, Ncomp
    integer :: ncomp_min, ncomp_max, i, j
    real (PREC) :: ssr, mse, sum_W_test
    real (PREC), dimension(:,:), allocatable :: X_test, X_train, Y_test, Y_train
    real (PREC), dimension(:,:), allocatable :: Xp
    real (PREC), dimension(:), allocatable :: weights_train, weights_test
    real (PREC), dimension(:), allocatable :: shift_x, scale_x, sval
    real (PREC), dimension(:,:), allocatable :: scores, loadings
    real (PREC), dimension(:,:), allocatable :: coefs
    real (PREC), dimension(:), allocatable :: intercept
    real (PREC), dimension(:,:), allocatable :: resid
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call get_dims (X, Y, add_intercept=conf%add_intercept, trans_x=conf%trans_x, &
        trans_y=conf%trans_y, nobs=nobs, nrhs=Nrhs, nlhs=nlhs, &
        nconst_add=Nconst_add, weights=weights, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Split sample ---

    Ntrain = Nobs - Ntest

    ! Create train and test subsamples in contiguous arrays.
    ! This already applies X^T and Y^T as needed!
    call extract_block_alloc (X, test_ifrom, Ntest, conf%trans_x, x_test, x_train)
    call extract_block_alloc (Y, test_ifrom, Ntest, conf%trans_y, Y_test, Y_train)

    allocate (weights_test(Ntest), source=1.0_PREC)
    sum_W_test = Ntest

    if (present(weights)) then
        call extract_block_alloc (weights, test_ifrom, Ntest, &
            x_block=weights_test, x_rest=weights_train)

        ! Sum of weights in test sample (for single LHS variable, single PC)
        sum_W_test = sum(weights_test)
    end if

   ! PC range to process
    ncomp_min = 0
    if (size(mse_ncomp) <= Nrhs)  ncomp_min = 1
    ncomp_max = min(Nrhs, size(mse_ncomp))

    ! --- Pre-compute effective weights ---

    ! Compute square root once since that's what will be needed below.
    ! Do NOT normalize weights, this will be done in caller routine, taking
    ! into account potentially different test sample sizes.
    weights_test(:) = sqrt(weights_test)

    ! --- Process regressors ---

    allocate (Xp(Ntrain,Nrhs))
    allocate (shift_x(Nrhs), scale_x(Nrhs))
    ! Fit the same preprocessing as the final model, using training data only.
    call transform_regressors (X_train, Xp, trans=.false., &
        center=conf%center, scale=conf%scale, shift_x=shift_x, &
        scale_x=scale_x, drop_na=.false., weights=weights_train, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(weights)) then
        ! Keep raw weights until their moments and positive total are checked.
        weights_train(:) = sqrt(weights_train / sum(weights_train))
        do i = 1, Nrhs
            Xp(:,i) = Xp(:,i) * weights_train
        end do
    end if

    ! --- PCA on training sample ---

    allocate (scores(Ntrain, ncomp_max))
    allocate (loadings(Nrhs, ncomp_max))
    allocate (sval(ncomp_max))

    ! Perform principal component analysis
    call pca (Xp, center=.false., scale=.false., trans_x=.false., &
        scores=scores, sval=sval, loadings=loadings, status=lstatus)
    if (.not. (NF_STATUS_OK .in. lstatus)) goto 100

    ! --- Run PCR for each number of PCs ---

    allocate (coefs(Nrhs,Nlhs))
    allocate (intercept(Nlhs))
    allocate (resid(Ntest,Nlhs))

    do Ncomp = ncomp_min, ncomp_max
        i = Ncomp - ncomp_min + 1

        ! Use the same centering for Y as for the training regressors.
        call pcr (Y_train, scores(:,1:Ncomp), sval(1:Ncomp), loadings(:,1:Ncomp), &
            coefs, shift_x=shift_x, scale_x=scale_x, &
            add_intercept=conf%add_intercept, center=conf%center, intercept=intercept, &
            weights_sqrt=weights_train, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        ! Compute predicted values on test sample
        call predict (X_test, coefs, resid, intercept, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        ! Compute (negative) residuals in place
        resid(:,:) = resid - Y_test
        do j = 1, Nlhs
            resid(:,j) = resid(:,j) * weights_test

            ssr = BLAS_NRM2 (Ntest, resid(:,j), 1)
            mse = ssr ** 2.0_PREC

            mse_ncomp(i,j) = mse
            mse_weights(i,j) = sum_W_test
        end do
    end do

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine pcr_cv_mse_mask (conf, X, Y, mask, test_ifrom, Ntest, mse_ncomp, &
        mse_weights, weights, status)
    type (lm_config), intent(in), optional :: conf
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:,:), contiguous :: Y
    logical, intent(in), dimension(:,:), contiguous :: mask
    integer, intent(in) :: test_ifrom
    integer, intent(in) :: Ntest
    real (PREC), intent(out), dimension(:,:), contiguous :: mse_ncomp
    real (PREC), intent(out), dimension(:,:), contiguous :: mse_weights
    real (PREC), intent(in), dimension(:), contiguous, optional :: weights
    type (status_t), intent(out), optional :: status


    integer :: Ntrain, Nobs, Nrhs, Nlhs, Nconst_add, Ncomp, Ntest_lhs, Ntrain_lhs
    integer :: ncomp_min, ncomp_max, ncomp_max_lhs, i, j
    real (PREC) :: ssr, mse_lhs, intercept, weights_lhs, s
    real (PREC), dimension(:,:), allocatable :: X_test, X_train, Y_test, Y_train
    real (PREC), dimension(:,:), allocatable :: X_train_T, X_test_T, X_lhs_T
    real (PREC), dimension(:,:), allocatable :: xp_lhs
    real (PREC), dimension(:), allocatable :: W_test_lhs, W_train_lhs
    real (PREC), dimension(:), allocatable :: Y_train_pack, Y_test_pack
    real (PREC), dimension(:), allocatable :: weights_train, weights_test
    logical, dimension(:,:), allocatable :: mask_train, mask_test
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_X_test_lhs_T
    real (PREC), dimension(:), pointer, contiguous :: ptr_Y_train_lhs, ptr_W_train_lhs
    real (PREC), dimension(:), pointer, contiguous :: ptr_Y_test_lhs
    real (PREC), dimension(:), allocatable :: W_train_all
    real (PREC), dimension(:), allocatable :: coefs
    real (PREC), dimension(:), allocatable :: resid
    real (PREC), dimension(:,:), allocatable :: scores, loadings
    real (PREC), dimension(:), allocatable :: sval, scale_x, shift_x
    real (PREC), dimension(:,:), allocatable :: scores_all, loadings_all
    real (PREC), dimension(:), allocatable :: sval_all, scale_x_all, shift_x_all
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_scores, ptr_loadings
    real (PREC), dimension(:), pointer, contiguous :: ptr_sval, ptr_shift_x, ptr_scale_x
    logical, dimension(:), allocatable :: lwork
    logical :: has_all_obs
    integer :: shp(2)
    integer (NF_ENUM_KIND) :: status_code
    type (status_t) :: lstatus
    type (lm_config) :: lconf
    logical :: has_weights

    target :: Y_train_pack, Y_test_pack, W_train_lhs, Y_train, Y_test
    target :: weights_train, X_lhs_T, X_test_T
    target :: scores, loadings, sval, shift_x, scale_x
    target :: scores_all, loadings_all, sval_all, shift_x_all, scale_x_all
    target :: W_train_all

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    lconf = conf

    has_weights = present(weights)

    call get_dims (X, Y, add_intercept=lconf%add_intercept, trans_x=lconf%trans_x, &
        trans_y=lconf%trans_y, nobs=nobs, nrhs=Nrhs, nlhs=nlhs, &
        nconst_add=Nconst_add, weights=weights, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! PC range to process
    ncomp_max = min(Nrhs, size(mse_ncomp, 1))
    ncomp_min = max(0, ncomp_max - size(mse_ncomp, 1) + 1)

    ! --- Split sample ---

    Ntrain = Nobs - Ntest

    ! Create train and test subsamples in contiguous arrays.
    ! This already applies X^T and Y^T as needed!
    call extract_block_alloc (X, test_ifrom, Ntest, conf%trans_x, x_test, x_train)
    call extract_block_alloc (Y, test_ifrom, Ntest, conf%trans_y, Y_test, Y_train)

    call extract_block_alloc (mask, test_ifrom, Ntest, conf%trans_y, &
        mask_test, mask_train)

    if (present(weights)) then
        call extract_block_alloc (weights, test_ifrom, Ntest, &
            x_block=weights_test, x_rest=weights_train)
    else
        ! Default uniform weights. Use this to determine relative sample sizes
        ! even for unweighted PCR.
        allocate (weights_test(Ntest), source=1.0_PREC)
    end if

    allocate (X_train_T(Nrhs,Ntrain))
    allocate (X_test_T(Nrhs,Ntest))

    X_train_T(:,:) = transpose (X_train)
    X_test_T(:,:) = transpose (X_test)

    ! Detect if there is any LHS variable which includes all obs
    allocate (lwork(Nlhs))
    lwork(:) = all(mask_train, dim=1)
    has_all_obs = any (lwork)
    deallocate (lwork)

    if (has_all_obs) then
        allocate (scores_all(Ntrain,Ncomp_max))
        allocate (loadings_all(Nrhs,Ncomp_max))
        allocate (sval_all(Ncomp_max))
        allocate (shift_x_all(Nrhs), scale_x_all(Nrhs))

        allocate (X_lhs_T(Ntrain,Nrhs))
        ! Validate raw training weights before normalizing them.
        call transform_regressors (X_train, X_lhs_T, trans=.false., &
            center=conf%center, scale=conf%scale, shift_x=shift_x_all, &
            scale_x=scale_x_all, drop_na=.false., weights=weights_train, &
            status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        if (has_weights) then
            s = sum(weights_train)
            allocate (W_train_all, source=weights_train)
            W_train_all(:) = sqrt(W_train_all / s)
            do i = 1, Nrhs
                X_lhs_T(:,i) = X_lhs_T(:,i) * W_train_all
            end do
        end if

        call pca (X_lhs_T, center=.false., scale=.false., &
            trans_x=.false., scores=scores_all, sval=sval_all, &
            loadings=loadings_all, status=lstatus)

        deallocate (X_lhs_T)

        if (.not. (NF_STATUS_OK .in. lstatus)) then
            status_code = lstatus
            goto 100
        end if
    end if

    ! Initialize here to avoid gfortran warnings
    status_code = NF_STATUS_OK

    mse_ncomp = 0.0
    mse_weights = 0.0

    !$omp parallel default(none) &
    !$omp shared(conf,X_train_T,X_test_T,mask_train,mask_test,weights_train) &
    !$omp shared (weights_test,Y_test,Y_train,has_weights,Nrhs,Ntest,Ntrain) &
    !$omp shared (Nlhs,ncomp_min,ncomp_max) &
    !$omp shared (scores_all,loadings_all,sval_all,shift_x_all,scale_x_all) &
    !$omp shared (W_train_all,mse_weights,mse_ncomp) &
    !$omp private (coefs,resid,X_lhs_T,Y_train_pack,Y_test_pack,W_train_lhs) &
    !$omp private (xp_lhs) &
    !$omp private (W_test_lhs,ptr_X_test_lhs_T,ptr_Y_train_lhs,ptr_Y_test_lhs) &
    !$omp private (ptr_W_train_lhs,shift_x,scale_x,scores,loadings,sval) &
    !$omp private (i,j,Ntrain_lhs,Ntest_lhs,ncomp_max_lhs,lstatus,s,shp) &
    !$omp private (intercept,weights_lhs,Ncomp,ssr,mse_lhs) &
    !$omp private (ptr_scores,ptr_loadings,ptr_sval,ptr_shift_x,ptr_scale_x) &
    !$omp reduction(ior: status_code)
    allocate (coefs(Nrhs))
    allocate (resid(Ntest))

    allocate (X_lhs_T(Nrhs,max(Ntrain,Ntest)))
    allocate (Y_train_pack(Ntrain), Y_test_pack(Ntest))
    allocate (W_train_lhs(Ntrain), W_test_lhs(Ntest))

    nullify (ptr_X_test_lhs_T)
    nullify (ptr_Y_train_lhs, ptr_Y_test_lhs)
    nullify (ptr_W_train_lhs)

    allocate (shift_x(Nrhs), scale_x(Nrhs))
    allocate (sval(ncomp_max))

    nullify (ptr_scores, ptr_loadings, ptr_sval, ptr_shift_x, ptr_scale_x)

    ! Initialize here to avoid gfortran warnings
    status_code = NF_STATUS_OK

    !$omp do
    loop_lhs: do j = 1, Nlhs

        ! Skip remaining tasks if we previously encountered an error
        if (status_code /= NF_STATUS_OK) cycle

        Ntrain_lhs = count (mask_train(:,j))
        Ntest_lhs = count (mask_test(:,j))

        if (Ntrain_lhs < Nrhs .or. Ntest_lhs == 0) cycle

        ncomp_max_lhs = min(ncomp_max, Ntrain_lhs)

        ! --- Contiguous training data ---

        if (Ntrain_lhs /= Ntrain) then

            call copy (X_train_T, X_lhs_T(:,1:Ntrain_lhs), mask_train(:,j), &
                dim=2, status=lstatus)
            if (.not. (NF_STATUS_OK .in. lstatus)) then
                status_code = lstatus
                cycle loop_lhs
            end if

            call copy (Y_train(:,j), Y_train_pack(1:Ntrain_lhs), mask_train(:,j))
            ptr_Y_train_lhs => Y_train_pack(1:Ntrain_lhs)

            ! Each response can have a different training sample and weighted mean.
            nullify (ptr_W_train_lhs)
            if (has_weights) then
                call copy (weights_train, W_train_lhs(1:Ntrain_lhs), mask_train(:,j))
                s = sum(W_train_lhs(1:Ntrain_lhs))
                ptr_W_train_lhs => W_train_lhs(1:Ntrain_lhs)
            end if

            call cond_alloc (xp_lhs, [Ntrain_lhs, Nrhs])
            call transform_regressors (X_lhs_T(:,1:Ntrain_lhs), xp_lhs, &
                trans=.true., center=conf%center, scale=conf%scale, &
                shift_x=shift_x, scale_x=scale_x, drop_na=.false., &
                weights=ptr_W_train_lhs, status=lstatus)
            if (lstatus /= NF_STATUS_OK) then
                status_code = lstatus
                cycle loop_lhs
            end if

            if (has_weights) then
                ! Weight the transformed data, retaining moments in original units.
                W_train_lhs(1:Ntrain_lhs) = sqrt(W_train_lhs(1:Ntrain_lhs) / s)
                do i = 1, Nrhs
                    xp_lhs(:,i) = xp_lhs(:,i) * W_train_lhs(1:Ntrain_lhs)
                end do
            end if

            ! --- Run PCA ---

            shp(1) = Ntrain_lhs
            shp(2) = ncomp_max_lhs
            call cond_alloc (scores, shp(1:2))

            shp(1) = Nrhs
            shp(2) = ncomp_max_lhs
            call cond_alloc (loadings, shp(1:2))

            ! Do not recenter or rescale after applying least-squares weights.
            call pca (xp_lhs, center=.false., scale=.false., &
                trans_x=.false., scores=scores, sval=sval(1:ncomp_max_lhs), &
                loadings=loadings, status=lstatus)
            if (.not. (NF_STATUS_OK .in. lstatus)) then
                status_code = lstatus
                cycle loop_lhs
            end if

            ptr_scores => scores
            ptr_loadings => loadings
            ptr_sval => sval(1:ncomp_max_lhs)
            ptr_shift_x => shift_x
            ptr_scale_x => scale_x
        else
            ! Use precomputed PCA data derived from all RHS observations
            ptr_Y_train_lhs => Y_train(:,j)
            nullify (ptr_W_train_lhs)
            if (has_weights) ptr_W_train_lhs => W_train_all

            ptr_scores => scores_all
            ptr_loadings => loadings_all
            ptr_sval => sval_all
            ptr_shift_x => shift_x_all
            ptr_scale_x => scale_x_all
        end if

        ! --- Prepare test sample ---

        if (Ntest_lhs /= Ntest) then
            call copy (X_test_T, X_lhs_T(:,1:Ntest_lhs), mask_test(:,j), &
                dim=2, status=lstatus)
            if (.not. (NF_STATUS_OK .in. lstatus)) then
                status_code = lstatus
                cycle loop_lhs
            end if
            ptr_X_test_lhs_T => X_lhs_T(:,1:Ntest_lhs)

            call copy (Y_test(:,j), Y_test_pack(1:Ntest_lhs), mask_test(:,j))
            ptr_Y_test_lhs => Y_test_pack(1:Ntest_lhs)
            if (has_weights) then
                W_test_lhs(1:Ntest_lhs) = pack (weights_test, mask_test(:,j))
            else
                ! Use uniform weights to "automatically" keep track of
                ! potentially different test sample sizes.
                ! Do NOT normalize weights here.
                W_test_lhs(1:Ntest_lhs) = 1.0_PREC
            end if
        else
            ptr_X_test_lhs_T => X_test_T
            ptr_Y_test_lhs => Y_test(:,j)
            if (has_weights) then
                W_test_lhs(:) = weights_test
            else
                ! Use uniform weights to "automatically" keep track of
                ! potentially different test sample sizes
                ! Do NOT normalize weights here.
                W_test_lhs(:) = 1.0_PREC
            end if
        end if

        ! Store sqrt. of test weights for later use
        weights_lhs = sum(W_test_lhs(1:Ntest_lhs))
        ! Do NOT normalie test weights, total MSE will be divided by sum
        ! of total weights.
        W_test_lhs(1:Ntest_lhs) = sqrt(W_test_lhs(1:Ntest_lhs))

        ! --- Run PCR for each number of PCs ---

        do Ncomp = ncomp_min, ncomp_max_lhs
            i = Ncomp - ncomp_min + 1

            ! Match Y centering to the preprocessing of this training sample.
            call pcr (ptr_Y_train_lhs, ptr_scores(:,1:Ncomp), ptr_sval(1:Ncomp), &
                ptr_loadings(:,1:Ncomp), coefs, shift_x=ptr_shift_x, &
                scale_x=ptr_scale_x, add_intercept=conf%add_intercept, &
                center=conf%center, intercept=intercept, &
                weights_sqrt=ptr_W_train_lhs, status=lstatus)
            if (lstatus /= NF_STATUS_OK) then
                status_code = lstatus
                cycle loop_lhs
            end if

            ! Compute predicted values on test sample
            call predict (ptr_X_test_lhs_T, coefs, resid(1:Ntest_lhs), &
                intercept, trans_x=.true., status=lstatus)
            if (lstatus /= NF_STATUS_OK) then
                status_code = lstatus
                cycle loop_lhs
            end if

            ! Compute residuals in place
            !   resid = (yhat - y)
            resid(1:Ntest_lhs) = resid(1:Ntest_lhs) - ptr_Y_test_lhs
            ! Squaring these residuals gives the weighted sum of squared errors.
            resid(1:Ntest_lhs) = resid(1:Ntest_lhs) * W_test_lhs(1:Ntest_lhs)

            ! Compute Euclidean norm
            !   |resid|_2 = sqrt(sum(wresid**2.0))
            ssr = BLAS_NRM2 (Ntest_lhs, resid(1:Ntest_lhs), 1)
            ! Undo square root applied by NRM2, so now MSE_LHS stores
            !   MSE_LHS = sum(wresid**2.0) = sum(weights) * MSE
            mse_lhs = ssr ** 2.0_PREC

            mse_ncomp(i,j) = mse_lhs
            ! Store sum of weights for given LHS test sample
            mse_weights(i,j) = weights_lhs
        end do
    end do loop_lhs
    !$omp end do

    deallocate (coefs,resid)
    deallocate (X_lhs_T,Y_train_pack,Y_test_pack,W_train_lhs,W_test_lhs)
    deallocate (shift_x,scale_x,sval)

    if (allocated(scores)) deallocate (scores)
    if (allocated(loadings)) deallocate (loadings)
    if (allocated(xp_lhs)) deallocate (xp_lhs)

    !$omp end parallel

    lstatus = status_code

100 continue

    if (present(status)) status = lstatus

end subroutine




subroutine lm_post_estim (model, X, y, trans_x, rsq, status)
    !*  LM_POST_ESTIM computes common post-estimation statistics such as
    !   R^2 for a given linear model. The model must have been estimated
    !   using before invoking this routine.
    type (lm_result), intent(in) :: model
        !*  Estimated model
    real (PREC), intent(in), dimension(:,:), contiguous :: X
        !*  Array of explanatory (RHS) variables used to estimate model.
        !   Array is assumed to have shape [NOBS, NRHS], unless
        !   TRANS_RHS is present and has value .TRUE.
    real (PREC), intent(in), dimension(:), contiguous :: y
        !*  Vector of observations of the dependent variable used to estimate
        !   model.
    logical, intent(in), optional :: trans_x
    real (PREC), intent(out), optional :: rsq
        !*  If present, contains the model's R^2 on exit.
    type (status_t), intent(out), optional :: status
        !*  Exit code

    logical :: ltrans_x
    real (PREC), dimension(:), allocatable :: resid, y_pred
    character (*), parameter :: NAME = 'POST_ESTIM'

    type (status_t) :: lstatus
    real (PREC) :: var_u, var_y, norm_u
    integer :: Nobs

    lstatus = NF_STATUS_OK

    call set_optional_arg (trans_x, model%conf%trans_x, ltrans_x)

    ! Check whether model has been estimated
    if (.not. allocated(model%coefs)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    call get_regressor_dims (X, trans=ltrans_x, nobs=Nobs)

    call check_cond (size(y) == Nobs, NAME, 'Y: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(rsq)) then

        call compute_resid (lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        ! Compute variance of residuals
        norm_u = BLAS_NRM2 (model%nobs, resid, 1)
        var_u = norm_u**2.0_PREC / model%nobs

        ! Compute variable of LHS variable Y
        call std (y, s=var_y, dof=0)
        var_y = var_y ** 2.0_PREC

        rsq = 0.0_PREC
        if (var_y > 0.0_PREC) then
            rsq = 1.0_PREC - var_u / var_y
        end if
    end if


100 continue

    if (allocated(y_pred)) deallocate (y_pred)
    if (allocated(resid)) deallocate (resid)

    if (present(status)) status  = lstatus

contains

    subroutine compute_predicted (status)
        type (status_t), intent(out) :: status
        if (.not. allocated (y_pred)) then
            allocate (y_pred(Nobs))
            call predict (X, model%coefs, y_pred, model%intercept, &
                trans_x=ltrans_x, trans_y=.false., trans_coefs=.false., &
                status=status)
        end if
    end subroutine

    subroutine compute_resid (status)
        type (status_t), intent(out) :: status

        call compute_predicted (status)
        if (status /= NF_STATUS_OK) return

        if (.not. allocated (resid)) then
            allocate (resid, source=y)
            resid(:) = resid - y_pred
        end if
    end subroutine

end subroutine



subroutine lm_predict_2d (model, X, Y_pred, trans_x, trans_y, status)
    type (lm_result), intent(in), target :: model
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(out), dimension(:,:), contiguous :: Y_pred
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
    type (status_t), intent(out), optional :: status

    logical :: ltrans_x, ltrans_y
    integer :: Nlhs
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'PREDICT'

    lstatus = NF_STATUS_OK

    call set_optional_arg (trans_x, model%conf%trans_x, ltrans_x)
    call set_optional_arg (trans_y, model%conf%trans_y, ltrans_y)

    call check_cond (allocated(model%coefs_multi), NAME, &
        'MODEL: COEFS_MULTI not allocated', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (allocated(model%intercept_multi), NAME, &
        'MODEL: INTERCEPT_MULTI not allocated', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call get_dims (X, Y_pred, model%coefs_multi, trans_x=ltrans_x, &
        trans_y=ltrans_y, nlhs=Nlhs, status=lstatus)

    call check_cond (model%nlhs == Nlhs, NAME, 'Y_PRED: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call predict (X, model%coefs_multi, Y_pred, intercept=model%intercept_multi, &
        trans_x=ltrans_x, trans_y=ltrans_y, trans_coefs=.false., status=lstatus)

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine lm_predict_1d (model, X, y_pred, trans_x, trans_y, status)
    type (lm_result), intent(in), target :: model
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(out), dimension(:), contiguous :: y_pred
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
        !*  Ignored. Present only for API-compatibility with 2d-variant
        !   of this routine.
    type (status_t), intent(out), optional :: status

    logical :: ltrans_x
    type (status_t) :: lstatus
    character (*), parameter :: NAME = 'PREDICT'
    real (PREC), dimension(:), pointer, contiguous :: ptr_coefs

    lstatus = NF_STATUS_OK

    call set_optional_arg (trans_x, model%conf%trans_x, ltrans_x)

    call check_cond (model%nlhs == 1, NAME, 'Y_PRED: Non-conformable array', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (allocated(model%coefs)) then
        ptr_coefs => model%coefs
    else if (allocated(model%coefs_multi)) then
        ptr_coefs => model%coefs_multi(:,1)
    else
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! Prevent intercept being added twice, so pass INTERCEPT argument
    ! only if intercept is not included in COEFS already.
    call predict (X, model%coefs, y_pred, intercept=model%intercept, &
        trans_x=ltrans_x, status=lstatus)

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine predict_2d (X, coefs, Y_pred, intercept, irhs, trans_x, &
        trans_y, trans_coefs, status)
    !*  PREDICT computes predicted values for a linear model, ie. it evaluates
    !   some variant of
    !       Y_pred = op(X) * op(B)
    !   or
    !       Y_pred^T = op(B) * op(X)
    !   where all three arrays can be optionally provided in transposed form.
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:,:), contiguous, target :: coefs
    !*  Coefficient matrix shaped [NCOEFS,NLHS]
    real (PREC), intent(out), dimension(:,:), contiguous :: Y_pred
    real (PREC), intent(in), dimension(:), contiguous, optional :: intercept
    integer, intent(in), dimension(:), contiguous, optional :: irhs
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
    logical, intent(in), optional :: trans_coefs
    type (status_t), intent(out), optional :: status

    logical :: ltrans_y

    call set_optional_arg (trans_y, .false., ltrans_y)

    if (ltrans_y) then
        call predict_trans_2d (X, coefs, Y_pred, intercept, irhs, trans_x, &
            trans_coefs, status)
    else
        call predict_notrans_2d (X, coefs, Y_pred, intercept, irhs, trans_x, &
            trans_coefs, status)
    end if

end subroutine



subroutine predict_notrans_2d (X, coefs, Y_pred, intercept, irhs, trans_x, &
        trans_coefs, status)
    !*  PREDICT_NOTRANS computes predicted values for linear models when
    !   the outcome variables are not transposed, ie. it returns
    !       Y = op(X) * op(B)
    !   where Y has shape [Nobs,Nlhs].
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:,:), contiguous, target :: coefs
        !*  Coefficient matrix shaped [NCOEFS,NLHS]
    real (PREC), intent(out), dimension(:,:), contiguous :: Y_pred
    real (PREC), intent(in), dimension(:), contiguous, optional :: intercept
    integer, intent(in), dimension(:), contiguous, optional :: irhs
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_coefs
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    logical :: ltrans_x, ltrans_coefs
    integer :: Nrhs, Nlhs, Nobs, Ncoefs, Nconst_coefs, Nrhs_in
    integer :: i, j, jrhs
    character (1) :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    real (PREC) :: beta, b0, b_ij
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_b
    real (PREC), dimension(:,:), allocatable :: X_in
    logical :: deallocate_coefs

    character (*), parameter :: NAME = 'PREDICT'

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (trans_x, .false., ltrans_x)
    call set_optional_arg (trans_coefs, .false., ltrans_coefs)

    call get_dims (X, Y_pred, coefs, trans_x=ltrans_x, &
        trans_y=.false., trans_coefs=ltrans_coefs, nrhs=Nrhs, nlhs=Nlhs, &
        nobs=Nobs, ncoefs=Ncoefs, nconst_coefs=Nconst_coefs, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(intercept)) then
        call check_cond (size(intercept) == Nlhs, NAME, &
            'INTERCEPT: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Add intercept ---

    if (present(intercept)) then
        do j = 1, Nlhs
            Y_pred(:,j) = intercept(j)
        end do
        ! Add X'B without intercept in GEMM, GEMV call
        beta = 1.0
    else if (Nconst_coefs == 1) then
        ! Initialize Y with intercepts from first row in COEFS
        do j = 1, Nlhs
            if (ltrans_coefs) then
                b0 = coefs(j,1)
            else
                b0 = coefs(1,j)
            end if
            Y_pred(:,j) = b0
        end do

        ! Add X'B without intercept in GEMM calls
        beta = 1.0
    else
        beta = 0.0
    end if

    ! --- Add XB in some form ---

    if (present(irhs)) then
        Nrhs_in = size(irhs)
    else
        Nrhs_in = Nrhs
    end if

    ! Exit immediately for degenerate matrix X
    if (Nrhs_in == 0) goto 100

    ! We need to account for the following variants:
    !   1.  X and COEFS not transposed: Y += X*B(i:,:)
    !   2.  X transposed, COEFS not transposed: Y += X^T * B(i:,:)
    !   3.  X not transposed, COEFS transposed: Y += X * B(:,i:)^T
    !   4.  X transposed, COEFS transposed: Y += X^T * B(:,i:)^T
    !
    !   where the initial row i depends on whether COEFS includes
    !   an intercept that is not in X.

    if (present(irhs)) then
        ! Include only selected RHS indices.

        ! If there is no intercept, Y_PREC needs to be initialized to zero.
        if (.not. (Nconst_coefs == 1 .or. present(intercept))) then
            Y_pred(:,:) = 0.0
        end if

        if (.not. ltrans_x) then
            ! X not tranposed: simple case with mostly memory-contiguous access,
            ! so that we can just loop through columns of Y and X.
            do i = 1, Nlhs
                do j = 1, Nrhs_in
                    jrhs = irhs(j)
                    if (ltrans_coefs) then
                        b_ij = coefs(i,Nconst_coefs+jrhs)
                    else
                        b_ij = coefs(Nconst_coefs+jrhs,i)
                    end if

                    if (b_ij == 0.0_PREC) cycle

                    call BLAS_AXPY (Nobs, b_ij, X(:,jrhs), 1, Y_pred(:,i), 1)
                end do
            end do
        else
            ! We need to compute
            !   Y = X(irhs,:)^T * op(B)(Nconst+irhs,:)
            ! which cannot be accomplished in a contiguous fashion in any
            ! reasonable way. Instead create a contiguous versions of X.
            ! and compute
            !   Y = X_in * op(B)(Nconst+irhs,:)

            ! Allocate packed X and COEFS
            allocate (X_in(Nobs, Nrhs_in))

            call replay_transform (X, X_in, trans=.true., ikeep=irhs)

            do i = 1, Nlhs
                do j = 1, Nrhs_in
                    jrhs = irhs(j)
                    if (ltrans_coefs) then
                        b_ij = coefs(i,Nconst_coefs+jrhs)
                    else
                        b_ij = coefs(Nconst_coefs+jrhs,i)
                    end if

                    if (b_ij == 0.0_PREC) cycle

                    call BLAS_AXPY (Nobs, b_ij, X_in(:,j), 1, Y_pred(:,i), 1)
                end do
            end do
        end if
    else
        ! Include all RHS variables.


        ! Y is not transposed, we compute some variant of
        ! Y = op(X) * op(B)(1+Nconst:,:)
        ldc = Nobs
        m = Nobs
        n = Nlhs
        k = Nrhs

        if (ltrans_x) then
            ! Case 2 + 4
            ! X provided in transposed form, needs to be transposed back to
            ! [Nobs,Nrhs]
            transa = 'T'
            lda = Nrhs
        else
            ! X provided in non-transposed form, remains as is
            ! Case 1 + 3
            transa = 'N'
            lda = Nobs
        end if

        ! Shift initial relevant row or column depending on whether intecept is
        ! (additionally) stored in COEFS.
        deallocate_coefs = .false.
        if (ltrans_coefs) then
            ! Case  3 + 4
            ! COEFS transposed by caller, needs to be transposed back
            ! to [NCOEFS,NLHS]
            transb = 'T'
            ldb = Nlhs
            ptr_b => coefs(:,1+Nconst_coefs:Ncoefs)
        else
            ! Case 1 + 2
            ! COEFS provided in non-transposed form, remains as is
            transb = 'N'
            ldb = Nrhs
            if (Nconst_coefs == 1) then
                deallocate_coefs = .true.
                allocate (ptr_b(Nrhs,Nlhs))
                do i = 1, Nlhs
                    call BLAS_COPY (Nrhs, coefs(2:Ncoefs,i), 1, ptr_b(:,i), 1)
                end do
            else
                ptr_b => coefs
            end if
        end if

        call BLAS_GEMM (transa, transb, m, n, k, 1.0_PREC, X, lda, &
            ptr_b, ldb, beta, Y_pred, ldc)

        if (deallocate_coefs) then
            deallocate (ptr_b)
        end if
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine predict_trans_2d (X, coefs, Y_pred, intercept, irhs, trans_x, &
        trans_coefs, status)
    !*  PREDICT_TRANS implements prediction for linear models when
    !   the outcome array is transposed, ie.
    !       Y^T = op(B) * op(X)
    real (PREC), intent(in), dimension(:,:), contiguous, target :: X
    real (PREC), intent(in), dimension(:,:), contiguous, target :: coefs
    !*  Coefficient matrix shaped [NCOEFS,NLHS]
    real (PREC), intent(out), dimension(:,:), contiguous :: Y_pred
    real (PREC), intent(in), dimension(:), contiguous, optional :: intercept
    integer, intent(in), dimension(:), contiguous, optional :: irhs
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_coefs
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    logical :: ltrans_x, ltrans_coefs
    integer :: Nrhs, Nlhs, Nobs, Ncoefs, Nconst_coefs, Nrhs_in
    integer :: i, j, jrhs
    character (1) :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    real (PREC) :: beta
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_b
    real (PREC), dimension(:,:), allocatable :: X_in_T, coefs_in_T
    character (*), parameter :: NAME = 'PREDICT'
    logical :: deallocate_coefs

    lstatus = NF_STATUS_OK

    ! --- Input processing ---

    call set_optional_arg (trans_x, .false., ltrans_x)
    call set_optional_arg (trans_coefs, .false., ltrans_coefs)

    call get_dims (X, Y_pred, coefs, trans_x=ltrans_x, trans_y=.true., &
        trans_coefs=ltrans_coefs, nrhs=Nrhs, nlhs=Nlhs, nobs=Nobs, &
        ncoefs=Ncoefs, nconst_coefs=Nconst_coefs, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(intercept)) then
        call check_cond (size(intercept) == Nlhs, NAME, &
            'INTERCEPT: Non-conformable array', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Add intercept ---

    if (present(intercept)) then
        Y_pred(:,:) = spread(intercept, dim=2, ncopies=Nobs)
        ! Add X'B without intercept in GEMM, GEMV call
        beta = 1.0
    else if (Nconst_coefs == 1) then
        ! Initialize Y with intercepts from first row in COEFS
        if (ltrans_coefs) then
            Y_pred(:,:) = spread(coefs(:,1), dim=2, ncopies=Nobs)
        else
            Y_pred(:,:) = spread(coefs(1,:), dim=2, ncopies=Nobs)
        end if

        ! Add X'B without intercept in GEMM calls
        beta = 1.0
    else
        beta = 0.0
    end if

    ! --- Add XB in some form ---

    if (present(irhs)) then
        Nrhs_in = size(irhs)
    else
        Nrhs_in = Nrhs
    end if

    ! Exit immediately for degenerate matrix X
    if (Nrhs_in == 0) goto 100

    ! We need to account for the following variants:
    !   1.  X and COEFS not transposed: Y^T += B(i:,:)^T * X^T
    !   2.  X transposed, COEFS not transposed: Y^T += B(i:,:)^T * X
    !   3.  X not transposed, COEFS transposed: Y^T += B(:,i:) * X^T
    !   4.  X transposed, COEFS transposed: Y^T += B(:,i:) * X
    !
    !   where the initial row i depends on whether COEFS includes
    !   an intercept that is not in X.

    if (present(irhs)) then
        ! Include only selected RHS indices.
        ! This is impossible to solve in a way that uses contiguous
        ! memory access, so we first create contiguous packed versions
        ! of X and COEFS.
        !
        ! We need to compute Y^T += B^T * X^T
        ! so we need to create transposed arrays of B and X, if these
        ! are not already on entry.

        ! Allocate packed X and COEFS
        allocate (X_in_T(Nrhs_in, Nobs))
        allocate (coefs_in_T(Nlhs,Nrhs_in))

        if (ltrans_x) then
            do i = 1, Nobs
                do j = 1, Nrhs_in
                    jrhs = irhs(j)
                    X_in_T(j,i) = X(jrhs,i)
                end do
            end do
        else
            do i = 1, Nobs
                do j = 1, Nrhs_in
                    jrhs = irhs(j)
                    X_in_T(j,i) = X(i,jrhs)
                end do
            end do
        end if

        if (ltrans_coefs) then
            ! Transposed COEFS with shape [NLHS,NCOEFS]
            do j = 1, Nrhs_in
                jrhs = irhs(j)
                coefs_in_T(:,j) = coefs(1:Nlhs,Nconst_coefs+jrhs)
            end do
        else
            ! Non-transposed COEFS with shape [NCOEFS,NLHS]
            do j = 1, Nrhs_in
                jrhs = irhs(j)
                coefs_in_T(:,j) = coefs(Nconst_coefs+jrhs,1:Nlhs)
            end do
        end if

        call BLAS_GEMM ('N', 'N', Nlhs, Nobs, Nrhs_in, 1.0_PREC, coefs_in_T, &
            Nlhs, X_in_T, Nrhs_in, beta, Y_pred, Nlhs)
    else
        ! Include all RHS variables.

        ldc = Nlhs
        m = Nlhs
        n = Nobs
        k = Nrhs

        if (ltrans_x) then
            ! Case 2 + 4
            ! X provided in transposed form, remains as is
            transb = 'N'
            ldb = Nrhs
        else
            ! Case 1 + 3
            ! X provided in non-transposed form, needs to be transposed to
            ! shape [Nrhs,Nobs]
            transb = 'T'
            ldb = Nobs
        end if

        deallocate_coefs = .false.
        if (ltrans_coefs) then
            ! Case 3 + 4
            ! Already transposed by caller, remains as is
            transa = 'N'
            lda = Nlhs
            ptr_b => coefs(:,1+Nconst_coefs:Ncoefs)
        else
            ! Case 1 + 2
            ! COEFS provided in non-transposed form, needs to be transposed
            ! to shape [Nlhs,Ncoefs]
            transa = 'T'
            lda = Nrhs
            if (Nconst_coefs == 1) then
                deallocate_coefs = .true.
                allocate (ptr_b(Nrhs,Nlhs))
                do i = 1, Nlhs
                    call BLAS_COPY (Nrhs, coefs(2:Ncoefs,i), 1, ptr_b(:,i), 1)
                end do
            else
                ptr_b => coefs
            end if
        end if

        call BLAS_GEMM (transa, transb, m, n, k, 1.0_PREC, ptr_b, lda, &
            X, ldb, beta, Y_pred, ldc)

        if (deallocate_coefs) then
            deallocate (ptr_b)
        end if
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine predict_1d (X, coefs, y_pred, intercept, irhs, trans_x, trans_y, &
        trans_coefs, status)
    real (PREC), intent(in), dimension(:,:), contiguous :: X
    real (PREC), intent(in), dimension(:), contiguous, target :: coefs
        !*  Coefficient vector of length NCOEFS, where NCOEFS = NRHS + 1
        !   if ADD_INTERCEPT=.TRUE., and NCOEFS = NRHS otherwise.
    real (PREC), intent(out), dimension(:), contiguous, target :: y_pred
    real (PREC), intent(in), optional :: intercept
    integer, intent(in), dimension(:), contiguous, optional :: irhs
    logical, intent(in), optional :: trans_x
    logical, intent(in), optional :: trans_y
        !*  Ignored. Present only for API compatibility with 2d variant of this
        !   routine.
    logical, intent(in), optional :: trans_coefs
        !*  Ignored. Present only for API compatibility with 2d variant of this
        !   routine.
    type (status_t), intent(out), optional :: status
        !*  Exit code.

    logical :: ltrans_x, ltrans_y, ltrans_coefs
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_y
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_coefs
    real (PREC), dimension(:), allocatable :: lintercept

    call set_optional_arg (trans_x, .false., ltrans_x)

    nullify (ptr_coefs, ptr_y)

    if (ltrans_x) then
        ! The optimal path to take is to compute
        !   Y^T = B^T * X^T
        ! at least if IRHS is not present, since that is only one GEMM call
        ! and no additional memory allocations.
        ptr_coefs(1:1,1:size(coefs)) => coefs
        ptr_y(1:1,1:size(y_pred)) => y_pred
        ltrans_y = .true.
        ltrans_coefs = .true.
    else
        ! The optimal path to take is
        !   Y = X * B
        ! since that causes no additional allocations and only uses contiguous
        ! memory access.
        ptr_coefs(1:size(coefs),1:1) => coefs
        ptr_y(1:size(y_pred),1:1) => y_pred
        ltrans_y = .false.
        ltrans_coefs = .false.
    end if

    if (present(intercept)) then
        allocate (lintercept(1), source=intercept)
    end if

    call predict (X, ptr_coefs, ptr_y, intercept=lintercept, irhs=irhs, &
        trans_x=ltrans_x, trans_y=ltrans_y, trans_coefs=ltrans_coefs, status=status)

end subroutine
