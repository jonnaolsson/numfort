

program test_numfort_stats_dmvnorm

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_stats, dmvnorm => dmvnorm_real64
    use numfort_common
    use numfort_common_testing

    use blas_interfaces, only: BLAS_GEMM => GEMM

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_rvs ()

end subroutine

subroutine test_rvs ()
    !*  Unit tests for MV normal random number generator

    real (PREC), dimension(:,:), allocatable :: xarr, vcv, vcv_hat, corr_x
    real (PREC), dimension(:), allocatable :: mean_x, xvec, std_x, mean_hat
    integer :: nobs, nvars, i, j, dim
    type (dmvnorm) :: mvnorm
    type (status_t) :: status
    logical :: vcv_ok, mean_ok

    call set_seed (5678)

    ! === Input checks (use 1d interface) ===
    nvars = 1
    nobs = 1
    allocate (xvec(nvars), mean_x(nvars), vcv(nvars, nvars))
    mean_x(:) = 0.0
    vcv(:,:) = 0.0
    status = NF_STATUS_UNDEFINED
    ! Call without any mean_x or VCV parameters
    call rvs (mvnorm, xvec, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Calling RVS without MEAN, COV parameters"

    ! Call without any VCV parameter
    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xvec, mean=mean_x, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Calling RVS without MEAN parameter"

    ! Call without any mean_x parameter
    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xvec, cov=vcv, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Calling RVS without COV parameter"

    ! Call with non-conformable dimensions
    deallocate (xvec)
    allocate (xvec(nvars+1))
    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xvec, mean=mean_x, cov=vcv, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Calling RVS with non-conformable X parameter"

    ! Call with non-conformable COV argument
    deallocate (xvec, vcv)
    allocate (xvec(nvars), vcv(nvars,nvars+1))
    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xvec, mean=mean_x, cov=vcv, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Calling RVS with non-conformable VCV parameter"

    deallocate (vcv)
    allocate (vcv(nvars,nvars), source=0.0_PREC)
    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xvec, mean=mean_x, cov=vcv, tol=-1.0_PREC, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Calling RVS with invalid TOL argument"

    deallocate (xvec, mean_x, vcv)

    ! === input checks for 2d-interface ===
    ! Additionally check DIM argument which is not available for 1d interface
    nvars = 3
    nobs = 10
    allocate (xarr(nvars,nobs), mean_x(nvars), vcv(nvars,nvars))
    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xarr, mean=mean_x, cov=vcv, dim=0, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Calling RVS with invalid DIM=0 argument"

    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xarr, mean=mean_x, cov=vcv, dim=3, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Calling RVS with invalid DIM=3 argument"

    ! Check with DIM that does not match the shape of XARR
    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xarr, mean=mean_x, cov=vcv, dim=1, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Calling RVS with DIM argument that does not match X shape"

    deallocate (xarr)
    allocate (xarr(nobs,nvars))
    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xarr, mean=mean_x, cov=vcv, dim=2, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Calling RVS with DIM argument that does not match X shape"

    deallocate (xarr, mean_x, vcv)

    ! === 1d interface ===
    ! Check MV-normal with 1 variable
    nvars = 1
    nobs = 1
    allocate (xvec(nvars), mean_x(nvars), vcv(nvars,nvars))
    mean_x(:) = 0.0
    vcv(1,1) = 2.0
    call dist_set_params (mvnorm, mean_x, vcv)
    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xvec, status=status)
    if (.not. (status == NF_STATUS_oK)) &
        error stop "Calling RVS with 1d-interface for 1 variable"
    deallocate (xvec, mean_x, vcv)

    ! Check call with n-dimensional vector
    nvars = 3
    nobs = 1
    allocate (xvec(nvars), mean_x(nvars), vcv(nvars,nvars), std_x(nvars))
    call random_number (mean_x)
    ! Create independent variables
    call random_number (std_x)
    call diag (std_x, vcv)
    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xvec, mean_x, vcv, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Calling RVS with 1d-interface for 3-dimensional random variable"
    deallocate (xvec, mean_x, vcv, std_x)

    ! === Calls to 2d interface ===
    ! Test with independent RVs
    nvars = 5
    nobs = int(1e6)
    allocate (xarr(nobs,nvars), mean_x(nvars), vcv(nvars,nvars), std_x(nvars))
    call random_number (std_x)
    call random_number (mean_x)
    mean_x(:) = 2.0 * (mean_x - 0.5)
    std_x(:) = std_x * 3.0
    call diag (std_x, vcv)
    status = NF_STATUS_UNDEFINED
    call rvs (mvnorm, xarr, mean_x, vcv, status=status)

    ! Compute sample mean
    allocate (mean_hat(nvars), vcv_hat(nvars,nvars))
    call mean (xarr, mean_hat, dim=1)
    call cov (xarr, vcv_hat, dim=1)

    vcv_ok = all_close (vcv_hat, vcv, rtol=0.0_PREC, atol=1.0e-2_PREC)
    mean_ok = all_close (mean_hat, mean_x, rtol=0.0_PREC, atol=1.0e-2_PREC)

    if (.not. (status == NF_STATUS_OK .and. mean_ok .and. vcv_ok)) &
        error stop "Calling RVS with independent VCV"

    deallocate (xarr, mean_x, std_x, vcv)
    deallocate (mean_hat, vcv_hat)

    ! === Test with partially dependent RVs ===
    nobs = int(1e6)
    nvars = 3
    dim = 2
    allocate (xarr(nvars,nobs), mean_x(nvars), vcv(nvars,nvars), std_x(nvars))
    allocate (corr_x(nvars,nvars), source=0.0_PREC)
    call random_number (mean_x)
    call random_number (std_x)
    std_x(:) = std_x * 2.0
    mean_x(:) = 2.0*(mean_x - 0.5)
    do i = 1, nvars
        corr_x(i,i) = 1.0
    end do
    corr_x(1,2) = 0.5
    corr_x(2,1) = 0.5

    ! Compute implied VCV
    do j = 1, nvars
        do i = 1, nvars
            vcv(i,j) = std_x(i)*std_x(j)*corr_x(i,j)
        end do
    end do

    status = NF_STATUS_UNDEFINED
    call dist_set_params (mvnorm, mean_x, vcv)
    call rvs (mvnorm, xarr, dim=dim, status=status)

    ! Compute sample moments
    allocate (mean_hat(nvars), vcv_hat(nvars,nvars))
    call mean (xarr, mean_hat, dim=dim)
    call cov (xarr, vcv_hat, dim=dim)

    vcv_ok = all_close (vcv_hat, vcv, rtol=0.0_PREC, atol=1.0e-2_PREC)
    mean_ok = all_close (mean_hat, mean_x, rtol=0.0_PREC, atol=1.0e-2_PREC)

    if (.not. (status == NF_STATUS_OK .and. mean_ok .and. vcv_ok)) &
        error stop "Calling RVS with partially interdependent random variables"

    deallocate (xarr, vcv, mean_x, std_x, corr_x)
    deallocate (mean_hat, vcv_hat)

    ! === Test with fully interdependent RVs ===
    nobs = int(1e6)
    nvars = 4
    dim = 1
    allocate (xarr(nobs,nvars), mean_x(nvars), vcv(nvars,nvars), std_x(nvars))
    allocate (corr_x(nvars,nvars))
    allocate (mean_hat(nvars), vcv_hat(nvars,nvars))
    call random_number (mean_x)
    mean_x(:) = 2.0*(mean_x - 0.5)

    ! Abuse VCV_HAT as workspace to construct VCV matrix as XX'
    call random_number (vcv_hat)
    vcv_hat(:,:) = 2.0 * (vcv_hat - 0.5)
    call BLAS_GEMM (transa="N", transb="T", m=nvars, n=nvars, k=nvars, &
        alpha=1.0_PREC, a=vcv_hat, lda=nvars, b=vcv_hat, ldb=nvars, &
        beta=0.0_PREC, c=vcv, ldc=nvars)
    ! Diagnostics: compute implied corr. matrix
    do j = 1, nvars
        do i = 1, nvars
            corr_x(i,j) = vcv(i,j) / sqrt(vcv(i,i)) / sqrt(vcv(j,j))
        end do
    end do

    status = NF_STATUS_UNDEFINED
    call dist_set_params (mvnorm, mean_x, vcv)
    call rvs (mvnorm, xarr, dim=dim, status=status)

    ! Compute sample moments
    call mean (xarr, mean_hat, dim=dim)
    call cov (xarr, vcv_hat, dim=dim)

    vcv_ok = all_close (vcv_hat, vcv, rtol=0.0_PREC, atol=1.0e-2_PREC)
    mean_ok = all_close (mean_hat, mean_x, rtol=0.0_PREC, atol=1.0e-2_PREC)

    if (.not. (status == NF_STATUS_OK .and. mean_ok .and. vcv_ok)) &
        error stop "Calling RVS with fully interdependent random variables"

end subroutine

end program
