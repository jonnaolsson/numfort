

subroutine __APPEND(dist_set_params,__PREC) (obj, mean, cov, status)
    !*  DIST_SET_PARAMS copies the given distribution parametes
    !   (arrays in the case of the MV normal distribution) into the
    !   corresponding object attributes
    integer, parameter :: PREC = __PREC
    type (__APPEND(dmvnorm,__PREC)), intent(inout) :: obj
    real (PREC), intent(in), dimension(:) :: mean
        !*  Vector of means
    real (PREC), intent(in), dimension(:,:) :: cov
        !*  Variance-covariance matrix
    type (status_t), intent(out), optional :: status
        !*  Optional exit code. Returns NF_STATUS_OK input arrays are
        !   non-empty and have conformable shapes.

    type (status_t) :: lstatus
    integer :: nd, shp(2)

    lstatus = NF_STATUS_OK

    nd = size(mean)
    if (nd < 1) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (size(cov,1) /= nd .or. size(cov,2) /= nd) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    call cond_alloc (obj%mean, nd)
    obj%mean(:) = mean

    shp = nd
    call cond_alloc (obj%cov, shp)
    obj%cov(:,:) = cov

100 continue
    if (present(status)) status = lstatus

end subroutine



subroutine __APPEND(rvs_check_input,__PREC) (obj, x, mean, cov, dim, tol, status)
    !*  RVS_CHECK_INPUT performs input validation for arguments to
    !   RVS routine for MV-normal distribution.
    integer, parameter :: PREC = __PREC
    type (__APPEND(dmvnorm,__PREC)), intent(in) :: obj
    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(in), dimension(:), contiguous, target, optional :: mean
    real (PREC), intent(in), dimension(:,:), contiguous, target, optional :: cov
    integer, intent(in), optional :: dim
    real (PREC), intent(in), optional :: tol
    type (status_t), intent(out) :: status

    integer :: ldim
    real (PREC), dimension(:), pointer, contiguous :: ptr_mean
    real (PREC), dimension(:,:), pointer, contiguous :: ptr_cov

    call dist_get_params (obj, mean, cov, ptr_mean, ptr_cov)

    if (.not. associated(ptr_mean) .or. .not. associated(ptr_cov)) goto 100

    if (size(ptr_mean) < 1) goto 100
    if (size(ptr_mean) /= size(ptr_cov,1)) goto 100
    if (size(ptr_cov,1) /= size(ptr_cov,2)) goto 100

    ldim = 1
    if (present(dim)) ldim = dim
    if (ldim /= 1 .and. ldim /= 2) goto 100
    if (size(x, 3-ldim) /= size(ptr_mean)) goto 100

    call check_positive (0.0_PREC, tol, 'tol', status)
    if (status /= NF_STATUS_OK) goto 100

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG
end subroutine



subroutine __APPEND(dist_get_params,__PREC) (obj, mean, cov, ptr_mean, ptr_cov)
    !*  DIST_GET_PARAMS returns pointers to the MV-normal parameter
    !   arrays that should be used.
    !
    !   If MEAN or COV arguments are provided by the caller, it is assumed
    !   that these should override the attributes of the OBJ distribution
    !   object, and hence pointers the the former will be returned.
    integer, parameter :: PREC = __PREC
    type (__APPEND(dmvnorm,__PREC)), intent(in), target :: obj
    real (PREC), intent(in), dimension(:), optional, target, contiguous :: mean
    real (PREC), intent(in), dimension(:,:), optional, target, contiguous :: cov
    real (PREC), intent(inout), dimension(:), pointer, contiguous :: ptr_mean
    real (PREC), intent(inout), dimension(:,:), pointer, contiguous :: ptr_cov

    nullify (ptr_mean, ptr_cov)

    if (present(mean)) then
        ptr_mean => mean
    else if (allocated(obj%mean)) then
        ptr_mean => obj%mean
    end if

    if (present(cov)) then
        ptr_cov => cov
    else if (allocated(obj%cov)) then
        ptr_cov => obj%cov
    end if

end subroutine



subroutine __APPEND(dmvnorm_rvs_1d,__PREC) (obj, x, mean, cov, &
        check_cov, tol, status)
    !*  DMVNORM_RVS_1D returns one (vector-valued) draw from a MV-normal
    !   distribution.
    integer, parameter :: PREC = __PREC
    type (__APPEND(dmvnorm,__PREC)), intent(in) :: obj
    real (PREC), intent(out), dimension(:), target, contiguous :: x
    real (PREC), intent(in), dimension(:), contiguous, optional :: mean
    real (PREC), intent(in), dimension(:,:), contiguous, optional :: cov
    logical, intent(in), optional :: check_cov
    real (PREC), intent(in), optional :: tol
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:,:), pointer, contiguous :: x2d
    integer, parameter :: dim = 2

    x2d(1:size(x),1:1) => x
    call rvs (obj, x2d, mean, cov, dim, check_cov, tol, status)

end subroutine



subroutine __APPEND(dmvnorm_rvs_2d,__PREC) (obj, x, mean, cov, dim, &
        check_cov, tol, status)
    !*  DMVNORM_RVS_2D draw samples from a multivariate normal distribution.
    !
    !   The distribution can either be parametrized by the MEAN and COV
    !   attributes of the OBJ distribution object, or by explicitly passing
    !   MEAN and/or COV arguments which override the corresponding attributes
    !   in OBJ.
    integer, parameter :: PREC = __PREC
    type (__APPEND(dmvnorm,__PREC)), intent(in) :: obj
        !*  DMVNORM distribution object
    real (PREC), intent(out), dimension(:,:), contiguous :: x
        !*  Array containing the desired sample draw on exit
    real (PREC), intent(in), dimension(:), contiguous, optional :: mean
        !*  Optional vector of means. If present, overrides values in
        !   OBJ%MEAN.
    real (PREC), intent(in), dimension(:,:), contiguous, optional :: cov
        !*  Optional variance-covariance matrix. If present, overrides values
        !   in OBJ%COV.
    integer, intent(in), optional :: dim
        !*  Optional DIM argument, controlling along which dimension
        !   observations are stacked.
        !   For DIM=1 (default), each row in X corresponds to one
        !   (multivariate) draw, and thus each column represents one
        !   (univariate) random variable.
        !   The interpretation of DIM is analogous to its usage in
        !   various routines computing sample moments, such as MEAN, STD and COV.
    logical, intent(in), optional :: check_cov
        !*  Optional. If true (the default), the given variance-covariance matrix
        !   is verified to be postive semi-definite using the given tolerance
        !   TOL. If the VCV matrix is not p.d.s., the routine exists immediately.
    real (PREC), intent(in), optional :: tol
        !*  Optional tolerance level when checking positive semi-definiteness
        !   of variance-covariance matrix. COV is assumed to be p.s.d. if
        !       |V*S*V' - COV| < TOL
        !   where COV = U*S*V' is the SVD factorization of COV.
    type (status_t), intent(out), optional :: status
        !*  Optional exit code.

    type (status_t) :: lstatus
    logical :: lcheck_cov
    integer :: ldim
    real (PREC), dimension(:), contiguous, pointer :: ptr_mean
    real (PREC), dimension(:,:), contiguous, pointer :: ptr_cov

    integer :: i, j, nvars, nobs, shp(2)
    real (PREC), dimension(:,:), allocatable :: rwork2d
    real (PREC) :: diff, ltol

    ! Arguments used for call to GESDD
    real (PREC), dimension(:), allocatable :: s
    real (PREC), dimension(:,:), allocatable :: vt
    real (PREC), dimension(1,1) :: u
    real (PREC), dimension(:), allocatable :: rwork
    integer, dimension(:), allocatable :: iwork
    integer :: info, lwork

    lstatus = NF_STATUS_OK
    call rvs_check_input (obj, x, mean, cov, dim, tol, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ldim = 1
    lcheck_cov = .true.
    ltol = 1.0e-8_PREC
    if (present(dim)) ldim = dim
    if (present(check_cov)) lcheck_cov = check_cov
    if (present(tol)) ltol = tol

    call dist_get_params (obj, mean, cov, ptr_mean, ptr_cov)

    ! dimension of mv-normally distributed random vector
    shp = shape(x)
    nvars = size(ptr_mean)
    nobs = shp(ldim)

    allocate (s(nvars), vt(nvars,nvars))
    allocate (rwork(1))
    allocate (iwork(8*nvars))

    ! Create copy of VCV matrix as it will be overwritten by GESDD
    allocate (rwork2d(nvars,nvars), source=ptr_cov)

    ! Workspace query for SVD decomposition
    call LAPACK_GESDD (jobz='O', m=nvars, n=nvars, a=rwork2d, lda=nvars, s=s, &
        u=u, ldu=1, vt=vt, ldvt=nvars, work=rwork, lwork=-1, iwork=iwork, &
        info=info)

    lwork = int(rwork(1))
    deallocate (rwork)
    allocate (rwork(lwork))

    ! Compute SVD decomposition (only S, V)
    call LAPACK_GESDD (jobz='O', m=nvars, n=nvars, a=rwork2d, lda=nvars, s=s, &
        u=u, ldu=1, vt=vt, ldvt=nvars, work=rwork, lwork=lwork, iwork=iwork, &
        info=info)

    if (info /= 0) then
        if (info < 0 ) then
            ! This should not happend since we checked input arguments above
            lstatus = NF_STATUS_INVALID_ARG
        else
            lstatus = NF_STATUS_NOT_CONVERGED
        end if
        goto 100
    end if

    ! Convert V' to V * diag(sqrt(s))
    ! We keep calling it VT, though...
    rwork2d(:,:) = transpose(vt)
    do i = 1, nvars
        vt(:,i) = rwork2d(:,i) * sqrt(s(i))
    end do

    ! Check whether COV was a positive semi-definite matrix, if applicable
    if (lcheck_cov) then
        call BLAS_GEMM (transa='N', transb='T', m=nvars, n=nvars, k=nvars, &
            alpha=1.0_PREC, a=vt, lda=nvars, b=vt, ldb=nvars, &
            beta=0.0_PREC, c=rwork2d, ldc=nvars)

        diff = maxval(abs(ptr_cov - rwork2d))
        if (diff > ltol) then
            lstatus = NF_STATUS_INVALID_STATE
            goto 100
        end if
    end if

    deallocate (rwork2d)

    ! Draw independent standard-normal realizations
    allocate (rwork2d(size(x,1), size(x,2)))
    do j = 1, size(x, 2)
        do i = 1, size(x, 1)
            rwork2d(i,j) = __APPEND(random_normal,__PREC) ()
        end do
    end do

    ! Transform to potentially correlated MV-normal variables
    if (ldim == 1) then
        ! Populate X with mean values first
        do i = 1, nvars
            x(:,i) = ptr_mean(i)
        end do

        ! Each column contains one realized variable, compute
        !   X = Z * V'
        ! where Z is the array of iid standard-normal variables.
        call BLAS_GEMM (transa="N", transb="T", m=nobs, n=nvars, k=nvars, &
            alpha=1.0_PREC, a=rwork2d, lda=nobs, b=vt, ldb=nvars, &
            beta=1.0_PREC, c=x, ldc=nobs)
    else
        ! Populate X with mean values first
        do i = 1, nobs
            x(:,i) = ptr_mean
        end do
        ! Each row contains one realized variable, compute
        !   X = V * Z
        ! instead.
        call BLAS_GEMM (transa='N', transb='N', m=nvars, n=nobs, k=nvars, &
            alpha=1.0_PREC, a=vt, lda=nvars, b=rwork2d, ldb=nvars, &
            beta=1.0_PREC, c=x, ldc=nvars)
    end if

    ! Add back mean for each variable

100 continue

    if (present(status)) status = lstatus

end subroutine
