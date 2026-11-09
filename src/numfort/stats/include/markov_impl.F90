


pure function is_trans_matrix (mat, tol, transposed) result(res)
    !* Check transition probabilities and row or column sums.
    !  Matrix shape is not validated.
    real (PREC), intent(in), dimension(:,:) :: mat
        !*  Transition matrix to check
    real (PREC), intent(in), optional :: tol
        !* Absolute tolerance on row or column sums (default: 1e-12).
        !  The absolute deviation from one must be strictly less than TOL.
    logical, intent(in), optional :: transposed
        !*  If true, assume that transition matrix is transposed, ie. columns
        !   sum to 1.0 instead of rows.
    logical :: res
        !* True if probabilities are nonnegative and sums satisfy the tolerance.

    real (PREC) :: ltol
    integer :: dim
    logical :: ltransposed
    real (PREC), dimension(:), allocatable :: rwork1d

    ltol = 1e-12_PREC
    ltransposed = .false.

    if (present(tol)) ltol = tol
    if (present(transposed)) ltransposed = transposed

    ! Dimension across which elements should sum to 1.0
    dim = 2
    if (ltransposed) dim = 1

    ! Check that all elements are non-negative and that rows (or columns if
    ! transposed=.true.) approximately sum to 1.0
    ! Explicitly allocate array to prevent temp. array compiler warnings.
    allocate (rwork1d(size(mat,3-dim)))
    rwork1d(:) = sum(mat, dim=dim)
    res = all(mat >= 0.0_PREC) .and. all(abs(rwork1d - 1.0_PREC) < ltol)

    deallocate (rwork1d)

end function



subroutine truncate_trans_matrix (mat, tol, transposed, status)
    !* Truncate small transition probabilities and normalize rows or columns.
    !  If any conditional distribution would have no positive mass left,
    !  leave the entire matrix unchanged, including when STATUS is absent.

    real (PREC), intent(inout), dimension(:,:) :: mat
        !* Markov transition matrix with elements MAT(i,j) = Prob[x'=j | x=i].
    real (PREC), intent(in) :: tol
        !* Probabilities strictly smaller than TOL are truncated to zero.
    logical, intent(in), optional :: transposed
        !* If true, MAT(i,j) = Prob[x'=i | x=j] and columns are normalized
        !  (default: false).
    type (status_t), intent(out), optional :: status
        !* NF_STATUS_OK on success; NF_STATUS_INVALID_ARG if truncation would
        !  leave any row (or column if TRANSPOSED) without positive mass.

    integer :: n, dim, i
    logical :: ltransposed
    real (PREC), dimension(:), allocatable :: rwork

    if (present (status)) status = NF_STATUS_OK

    ltransposed = .false.
    if (present (transposed)) ltransposed = transposed

    dim = 2
    if (ltransposed) then
        dim = 1
    end if

    ! Check retained mass before modifying any entries.
    n = size (mat, 3 - dim)
    allocate (rwork(n))
    rwork(:) = sum (mat, dim=dim, mask=mat >= tol)
    if (any (rwork <= 0.0_PREC)) then
        if (present (status)) status = NF_STATUS_INVALID_ARG
        return
    end if

    where (mat < tol)
        mat = 0.0_PREC
    end where

    ! Normalize such that conditional trans. probabilities sum to 1.0
    if (ltransposed) then
        do i = 1, size (mat, 2)
            mat(:,i) = mat(:,i) / rwork(i)
        end do
    else
        do i = 1, size (mat, 2)
            mat(:,i) = mat(:,i) / rwork
        end do
    end if

end subroutine



subroutine markov_approx_input_checks (rho, sigma, states, tm, status)
    !*  MARKOV_APPROX_INPUT_CHECKS validates inputs for the ROUWENHORST
    !   and TAUCHEN routines.
    real (PREC), intent(in) :: rho
    real (PREC), intent(in) :: sigma
    real (PREC), intent(in), dimension(:) :: states
    real (PREC), intent(in), dimension(:,:) :: tm
    type (status_t), intent(out) :: status

    integer :: n

    status = NF_STATUS_OK
    n = size(states)

    ! Input checks
    if ((abs(rho) >= 1.0_PREC) .and. (n > 1)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if ((sigma < 0.0_PREC) .and. (n > 1)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (n /= size(tm,1) .or. n /= size(tm,2)) then
        status = NF_STATUS_INVALID_ARG
        goto 100
    end if

100 continue

end subroutine



subroutine rouwenhorst (rho, sigma, states, tm, sigma_cond, status)
    !*  ROUWENHORST returns the discretized approximation of an AR(1) process
    !   using the Rouwenhorst method.
    !   Matching empty arrays are a successful no-op. A singleton grid has
    !   state zero and transition probability one; rho and sigma are ignored
    !   for both cases. With more than one state, abs(rho) < 1 and sigma >= 0
    !   are required. Zero sigma sets every state to zero and preserves the
    !   rho-dependent transition matrix.
    real (PREC), intent(in) :: rho
        !*  Autocorrelation coefficient
    real (PREC), intent(in) :: sigma
        !*  Conditional or unconditional standard deviation
        !   (depending on SIGMA_COND)
    real (PREC), intent(out), dimension(:) :: states
        !*  Contains discretized state vector on successful exit.
    real (PREC), intent(out), dimension(:,:) :: tm
        !*  Contains transition matrix on successful exit.
    logical, intent(in), optional :: sigma_cond
        !*  If present and true, argument SIGMA is interpreted as the
        !   conditional standard deviation (default: true).
    type (status_t), intent(out), optional :: status
        !*  Exit status (optional)

    type (status_t) :: lstatus
    integer :: n, i, j, m
    real (PREC), dimension(:,:), allocatable :: work
    real (PREC) :: eps, p, sigma_z
    logical :: lsigma_cond

    lstatus = NF_STATUS_OK

    n = size(states)

    ! Check inputs
    call markov_approx_input_checks (rho, sigma, states, tm, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Permit calling routine with size-0 arrays, do nothing in that case
    if (n == 0) goto 100

    lsigma_cond = .true.
    if (present(sigma_cond)) lsigma_cond = sigma_cond

    ! Treat degenerate case with one state separately
    if (n == 1) then
        tm(1,1) = 1.0_PREC
        states(1) = 0.0_PREC
        goto 100
    end if

    if (lsigma_cond) then
        sigma_z = sigma / sqrt(1.0_PREC - rho ** 2.0_PREC)
    else
        sigma_z = sigma
    end if

    ! Construct discretized state space
    eps = sqrt(n-1.0_PREC) * sigma_z
    call linspace (states, -eps, eps)

    ! At this point the state space has at least N=2
    allocate (work(n+1,n+1), source=0.0_PREC)

    p = (1.0_PREC+rho)/2.0_PREC

    ! Initial 2-by-2 transition matrix
    tm(1,1) = p
    tm(1,2) = 1.0_PREC - p
    tm(2,1) = 1.0_PREC - p
    tm(2,2) = p

    do m = 2, n-1
        call rouwenhorst_pad_matrix (tm(1:m,1:m), work)

        do j = 1, m+1
            do i = 1, m+1
                tm(i,j) = p * work(1+i,1+j) &
                    + (1.0_PREC - p) * work(1+i,j) &
                    + (1.0_PREC - p) * work(i,j+1) &
                    + p * work(i,j)
            end do
        end do

        do j = 1, m+1
            tm(2:m,j) = tm(2:m,j) / 2.0_PREC
        end do
    end do

    ! normalize row sums to 1
    do i = 1, n
        tm(i,:) = tm(i,:) / sum(tm(i,:))
    end do

100 continue
    if (present(status)) status = lstatus

    if (allocated(work)) deallocate (work)

end subroutine


subroutine rouwenhorst_pad_matrix (mat, res)
    real (PREC), intent(in), dimension(:,:) :: mat
    real (PREC), intent(out), dimension(:,:) :: res

    integer :: n, i

    n = size(mat, 1)

    res(1:n+2,1:n+2) = 0.0_PREC
    do i = 1, n
        res(2:n+1,i+1) = mat(1:n,i)
    end do
end subroutine



subroutine tauchen (rho, sigma, states, tm, m, sigma_cond, status)
    !*  TAUCHEN returns the discretized approximation of an AR(1) process
    !   using the Tauche (1986) method.
    real (PREC), intent(in) :: rho
        !*  Autocorrelation coefficient
    real (PREC), intent(in) :: sigma
        !*  Conditional (default) or unconditional standard deviation
        !   (interpretation depends on the argument SIGMA_COND)
    real (PREC), intent(out), dimension(:) :: states
        !*  Contains discretized state space on successful exit.
    real (PREC), intent(out), dimension(:,:) :: tm
        !*  Contains transition matrix on successful exit.
    real (PREC), intent(in), optional :: m
        !*  If present, determines the range of the discretized state space
        !   as (sigma * M), where SIGMA is the unconditional standard deviation
        !   (default: 3).
    logical, intent(in), optional :: sigma_cond
        !*  If present and true, interpret SIGMA argument as the conditional
        !   standard deviation (default: true).
    type (status_t), intent(out), optional :: status
        !*  Exit status (optional)

    logical :: lsigma_cond
    real (PREC) :: lm
    real (PREC) :: sigma_z, sigma_e
    real (PREC) :: zjk, w2, zj1, zjn
    integer :: i, k, n
    type (status_t) :: lstatus

    type (dnorm), parameter :: norm = dnorm (loc=0.0_PREC,scale=1.0_PREC)

    lstatus = NF_STATUS_OK

    ! Input checks
    call markov_approx_input_checks (rho, sigma, states, tm, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    n =  size(states)

    if (n == 0) then
        ! Permit calling routine with size-0 arrays; do nothing in that case
        goto 100
    else if (n == 1) then
        ! Degenerate process with single state and certain transition
        states(1) = 0.0
        tm(1,1) = 1.0
        goto 100
    end if

    lsigma_cond = .true.
    lm = 3

    if (present(sigma_cond)) lsigma_cond = sigma_cond
    if (present(m)) lm = m

    if (lsigma_cond) then
        sigma_z = sigma / sqrt(1.0_PREC - rho ** 2.0_PREC)
        sigma_e = sigma
    else
        sigma_z = sigma
        sigma_e = sigma_z * sqrt(1.0_PREC - rho ** 2.0_PREC)
    end if

    call linspace(states, -sigma_z * lm, sigma_z * lm)
    w2 = (states(2) - states(1)) / 2.0_PREC

    do i = 1, n
        do k = 2, n-1
            zjk = states(k) - rho * states(i)

            tm(i, k) = cdf (norm, (zjk + w2) / sigma_e) &
                - cdf (norm, (zjk - w2) / sigma_e)
        end do

        zj1 = states(1) - rho * states(i)
        zjn = states(n) - rho * states(i)

        tm(i, 1) = cdf (norm, (zj1 + w2) / sigma_e)
        tm(i, n) = 1 - cdf (norm, (zjn - w2) / sigma_e)
    end do

    ! normalize row sums to 1
    do i = 1, n
        tm(i,:) = tm(i,:) / sum(tm(i,:))
    end do

100 continue
    if (present(status)) status = lstatus

end subroutine



subroutine ergodic_dist (tm, edist, inverse, maxiter, &
        transposed, tol, initial, work, status)
    !*  ERGODIC_DIST returns the ergodic distribution implied by a given
    !   Markov transition matrix using one of two methods (see argument INVERSE).
    real (PREC), intent(in), dimension(:,:), contiguous :: tm
        !*  Nonempty, square Markov transition matrix.
    real (PREC), intent(inout), dimension(:), contiguous :: edist
        !*  Contains the ergodic distribution on successful exit. Unchanged
        !   for invalid dimensions or a singular inverse system.
    logical, intent(in), optional :: inverse
        !*  If present and true, compute ergodic distribution using the
        !   "inversion" method (default: true). Otherwise, iterate on the
        !   transition matrix until element-wise convergence.
    integer, intent(in), optional :: maxiter
        !*  Maximum number of iterations if INVERSE is .FALSE.
    logical, intent(in), optional :: transposed
        !*  If present and true, assume that TM is a transposed transition
        !   matrix, ie TM(i,j) = Prob(x'=i|x=j). (default: false)
    real (PREC), intent(in), optional :: tol
        !*  Convergence tolerance, only used for iterative method.
        !   (default: 1e-12)
    logical, intent(in), optional :: initial
        !*  If present and true, use the values in EDIST as the initial
        !   guess for the ergodic distribution (default: false).
        !   Only used for iterative method.
    type (workspace), intent(inout), optional, target :: work
        !*  Optional caller-owned workspace, retained on success or failure.
    type (status_t), intent(out), optional :: status
        !*  Exit code. Invalid dimensions return NF_STATUS_INVALID_ARG.
        !   A singular inverse system, including nonunique stationary
        !   distributions, returns NF_STATUS_INVALID_STATE.

    ! Working arrays
    real (PREC), dimension(:,:), pointer, contiguous :: tm_inv, tm_T
    real (PREC), dimension(:), pointer, contiguous :: mu1, mu2, diff_mu
    real (PREC), dimension(:), pointer, contiguous :: rwork_inv
    integer, dimension(:), pointer, contiguous :: iwork_inv

    type (workspace), pointer :: ptr_work

    integer :: nrwrk, niwrk, nrwrk_inv, niwrk_inv
    integer, dimension(2) :: shp2d

    integer, parameter :: NITER = 100
        !*  Number of iterations to perform before checking for convergence
        !   when iterative method is used.
    integer :: i, j, lmaxiter, n, NBATCH
    logical :: linverse, ltransposed, linitial
    real (PREC) :: ltol
    type (status_t) :: lstatus
    real (PREC) :: mass
    real (PREC), parameter :: zero_trunc = 1.0d-15

    ! Arguments for BLAS routines
    character (1), parameter :: trans = 'N'
    integer :: m, lda
    integer, parameter :: incx = 1, incy = 1
    real (PREC), parameter :: beta = 0.0_PREC
    real (PREC) :: alpha

    lstatus = NF_STATUS_OK

    nullify (ptr_work)
    nullify (tm_T, tm_inv)
    nullify (mu1, mu2, diff_mu)

    n = size(tm, 1)

    if (n < 1 .or. size(tm, 2) /= n .or. size(edist) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    linitial = .false.
    linverse = .true.
    lmaxiter = 10000
    ltransposed = .false.
    ltol = 1.0e-12_PREC
    if (present(inverse)) linverse = inverse
    if (present(maxiter)) lmaxiter = maxiter
    if (present(transposed)) ltransposed = transposed
    if (present(tol)) ltol = tol
    if (present(initial)) linitial = initial

    ! Allocate workspace arrays
    call assert_alloc_ptr (work, ptr_work)
    ! Clear any internal state in workspace object, in particular index offsets
    ! (this does not deallocate working arrays)
    call workspace_reset (ptr_work)

    if (linverse) then
        ! Query workspace requirements for INV
        call inv_work_query (tm, nrwrk_inv, niwrk_inv, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        nrwrk = 2*n*n + nrwrk_inv
        niwrk = niwrk_inv
    else
        nrwrk = n*n + 3*n
        niwrk = n
    end if

    call assert_alloc (ptr_work, nrwrk=nrwrk, niwrk=niwrk)

    shp2d = n
    call workspace_get_ptr (ptr_work, shp2d, tm_T, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (ltransposed) then
        tm_T(:,:) = tm
    else
        tm_T(:,:) = transpose(tm)
    end if

    ! compute ergodic distribution using "inverse" method
    if (linverse) then

        do i = 1, n
            tm_T(i,i) = tm_T(i,i) - 1.0_PREC
        end do

        tm_T(n, :) = 1.0_PREC

        call workspace_get_ptr (ptr_work, shp2d, tm_inv, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        call workspace_get_ptr (ptr_work, nrwrk_inv, rwork_inv, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        call workspace_get_ptr (ptr_work, niwrk_inv, iwork_inv, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        call inv (tm_T, tm_inv, rwork=rwork_inv, iwork=iwork_inv, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
        
        edist = tm_inv(:, n)

        nullify (tm_inv, rwork_inv, iwork_inv, tm_T)
    else
        ! Check for convergence every NBATCH iterations
        NBATCH = int(ceiling(real(lmaxiter,real64)/NITER))

        ! Set GEMV arguments that do not change during iteration
        lda = n
        m = n

        call workspace_get_ptr (ptr_work, n, mu1, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        call workspace_get_ptr (ptr_work, n, mu2, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        call workspace_get_ptr (ptr_work, n, diff_mu, lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100

        ! compute ergodic distribution by iteration
        if (linitial) then
            ! Use EDIST as initial distribution if requested by caller.
            mu1(1:n) = edist
        else
            ! No initial distribution provided, initialize with uniform
            ! distribution over states.
            mu1(1:n) = 1.0_PREC / n
        end if

        do i = 1, NBATCH
            do j = 1, NITER
                alpha = 1.0_PREC
                call BLAS_GEMV (trans, m, n, alpha, tm_T, lda, mu1, incx, &
                    beta, mu2, incy)

                call swap (mu1, mu2)
            end do

            ! Check for convergence
            call BLAS_COPY (n, mu1, incx, diff_mu, incy)
            ! obtain DIFF_MU = MU1 - MU2
            alpha = -1.0_PREC
            call BLAS_AXPY (n, alpha, mu2, incx, diff_mu, incy)

            ! check whether convergence in mu was achieved
            if (all(abs(diff_mu) < ltol)) exit
        end do

        edist = mu1

        nullify (mu1, mu2, diff_mu, tm_T)

        if (i > NBATCH) then
            lstatus = NF_STATUS_NOT_CONVERGED
            goto 100
        end if
    end if

    ! normalize to 1
    do i = 1, n
        if (edist(i) < 0.0_PREC) then
            if (abs(edist(i)) < zero_trunc) then
                edist(i) = 0.0_PREC
            else
                lstatus = NF_STATUS_INVALID_STATE
                goto 100
            end if
        end if
    end do

    mass = sum(edist)
    edist = edist / mass

100 continue

    if (present(status)) status = lstatus

    ! Clean up local WORKSPACE object if none was passed by client code
    call assert_dealloc_ptr (work, ptr_work)

end subroutine



subroutine moments (states, tm, mean, acorr, sigma, sigma_e, &
        edist, inverse, transposed, status)
    !*  MOMENTS computes selected moments implied by a Markov
    !   process with given state space and transition matrix.
    real (PREC), intent(in), dimension(:), contiguous :: states
        !*  Markov chain state space
    real (PREC), intent(in), dimension(:,:), contiguous :: tm
        !*  Transition matrix with element (i,j) being the transition
        !   probability Prob[x_t=s_j|x_{t-1}=s_i]
    real (PREC), intent(out), optional :: mean
        !*  If present, contains unconditional mean on successful exit.
    real (PREC), intent(out), optional :: acorr
        !*  On exit, contains implied autocorrelation coefficient
    real (PREC), intent(out), optional :: sigma
        !*  On exit, contains implied unconditional std. dev.
    real (PREC), intent(out), optional :: sigma_e
        !*  On exit, contains implied std. dev. of AR(1) error
        !   term (ie. conditional std. dev.)
    real (PREC), intent(in), dimension(:), contiguous, optional, target :: edist
        !*  If present, contains the ergodic distribution implied by transition
        !   matrix. Will be computed if missing.
    logical, intent(in), optional :: inverse
        !*  Use inverse of transition matrix when computing ergodic distribution
        !   (ignored if EDIST argument is present).
    logical, intent(in), optional :: transposed
        !*  If present and true, assume that transition matrix TM is provided
        !   in transposed format, ie each element (i,j) represents the
        !   probabolity Prob[x'=i | x=j].
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(:), contiguous, pointer :: ptr_pmf
    real (PREC), dimension(:), allocatable :: x
    real (PREC) :: mean_x, var_x, covar_x, acorr_x
    integer :: i, j, n
    type (status_t) :: lstatus
    logical :: ltransposed

    lstatus = NF_STATUS_OK
    n = size(states)

    ltransposed = .false.
    if (present(transposed)) ltransposed = transposed

    if (.not. present(edist)) then
        allocate (ptr_pmf(n))
        call markov_ergodic_dist (tm, ptr_pmf, inverse=inverse, &
            transposed=transposed, status=lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    else
        ptr_pmf => edist
    end if

    ! Unconditional mean
    mean_x = dot_product (ptr_pmf, states)

    ! Demeaned states
    allocate (x(n), source=states)
    x(:) = x - mean_x

    ! Unconditional variance
    var_x = 0.0
    do i = 1, n
        var_x = var_x + ptr_pmf(i) * x(i) ** 2
    end do

    ! Compute autocovariane
    covar_x = 0.0
    if (ltransposed) then
        do j = 1, n
            do i = 1, n
                covar_x = covar_x + x(i) * x(j) * tm(j,i) * ptr_pmf(i)
            end do
        end do
    else
        do j = 1, n
            do i = 1, n
                covar_x = covar_x + x(i) * x(j) * tm(i,j) * ptr_pmf(i)
            end do
        end do
    end if

    if (var_x > 0.0_PREC) then
        acorr_x = covar_x / var_x
    else
        ! If var_x == 0 we have a degenerate distribution with only one
        ! state, hence there is perfect correlation.
        acorr_x = 1.0_PREC
    end if

    if (present(mean)) mean = mean_x
    if (present(sigma_e)) sigma_e = sqrt((1.0_PREC-acorr_x**2.0_PREC) * var_x)
    if (present(sigma)) sigma = sqrt(var_x)
    if (present(acorr)) acorr = acorr_x

100 continue

    if (present(status)) status = lstatus

    call assert_dealloc_ptr (edist, ptr_pmf)
    if (allocated(x)) deallocate(x)

end subroutine



subroutine pmf_to_histogram (nobs, pmf, hist, status)
    !*  PMF_TO_HISTOGRAM converts the PMF over a set of discrete states into an
    !   "optimal" histograms for the given sample size so that sum(hist) = nobs.
    integer, intent(in) :: nobs
        !*  Number of observations
    real (PREC), intent(in), dimension(:) :: pmf
        !*  PMF over states
    integer, intent(out), dimension(:) :: hist
        !*  Histogram of observations for each state such that sum(hist) = nobs.
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    type (status_t) :: lstatus
    logical, dimension(:), allocatable :: mask
    integer :: i

    lstatus = NF_STATUS_OK

    ! --- Input argument checking ---

    if (.not. all(pmf >= 0.0_PREC) .or. abs(sum(pmf) - 1.0_PREC) > 1.0e-10_PREC) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (nobs <= 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! --- Create histogram ---

    hist(:) = int(nobs * pmf)

    allocate (mask(size(pmf)))
    mask(:) = (pmf > 0.0_PREC)

    do while (sum(hist) > nobs)
        i = maxloc (hist / nobs - pmf, dim=1, mask=mask)
        hist(i) = hist(i) - 1
    end do

    do while (sum(hist) < nobs)
        i = minloc (hist / nobs - pmf, dim=1, mask=mask)
        hist(i) = hist(i) + 1
    end do

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine trans_histogram (nobs, transm, trans_hist, status)
    !*  TRANS_HISTOGRAM returns the optiomal transition "histogram" for a given
    !   sample size so that sum(trans_hist,1) = sum(trans_hist,2), i.e., after
    !   each round of transitions the PMF over states remains unchanged.
    integer, intent(in) :: nobs
        !*  Number of observations
    real (PREC), intent(in), dimension(:,:), contiguous :: transm
        !*  Transition matrix of Markov process
    integer, intent(out), dimension(:,:), contiguous :: trans_hist
        !*  Transition "histogram" for given sample size
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    type (status_t) :: lstatus
    real (PREC), dimension(:), allocatable :: edist
    integer, dimension(:), allocatable :: edist_hist, hist_to
    integer :: nstates, i, ihi, ilo, ibest
    real (PREC) :: eps, eps_min

    lstatus = NF_STATUS_OK

    ! --- Input argument checking ---

    if (.not. all(transm >= 0.0_PREC) .or. &
            all(abs(sum(transm, dim=2) - 1.0_PREC) > 1.0e-10_PREC)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (nobs <= 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    nstates = size (transm, 1)
    allocate (edist(nstates))

    ! Ergodic distribution implied by transition matrix
    call ergodic_dist (transm, edist, inverse=.TRUE., status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Histogram for ergodic distribution
    allocate (edist_hist(nstates))
    call pmf_to_histogram (nobs, edist, edist_hist, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Create approximate histogram for each transition matrix row
    do i = 1, nstates
        call pmf_to_histogram (edist_hist(i), transm(i, :), trans_hist(i, :), lstatus)
    end do

    allocate (hist_to(nstates))

    hist_to(:) = sum (trans_hist, dim=1)

    do while (any (hist_to /= edist_hist))
        ! Column index where sum is too high
        ihi = minloc (edist_hist - hist_to, dim=1)
        ! Column index where sum is too low
        ilo = maxloc (edist_hist - hist_to, dim=1)

        eps_min = ieee_value (1.0_PREC, IEEE_POSITIVE_INF)
        ibest = -1

        ! Find adminissible row to flip one element between columns
        do i = 1, nstates
            ! Do not add to elements which are zero in original matrix
            if (transm(i,ilo) == 0.0_PREC) cycle
            ! Do not subtract from elements that are already zero (can this happen?)
            if (trans_hist(i, ihi) == 0) cycle
            ! Compute deviation from original transition matrix
            eps = abs((trans_hist(i,ihi) - 1) / edist_hist(i) - transm(i,ihi)) &
                + abs((trans_hist(i,ilo) + 1) / edist_hist(i) - transm(i,ilo))

            if (eps < eps_min) then
                ibest = i
                eps_min = eps
            end if
        end do

        if (ibest == -1) then
            lstatus = NF_STATUS_INVALID_STATE
            goto 100
        end if

        trans_hist(ibest,ihi) = trans_hist(ibest,ihi) - 1
        trans_hist(ibest,ilo) = trans_hist(ibest,ilo) + 1

        hist_to(:) = sum (trans_hist, dim=1)
    end do

100 continue

    if (present(status)) status = lstatus

end subroutine
