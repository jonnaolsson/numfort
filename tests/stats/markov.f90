
program test_stats_markov

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: ieee_exceptions, only: ieee_get_flag, ieee_set_flag, &
        ieee_divide_by_zero, ieee_invalid

    use numfort_common
    use numfort_stats
    use numfort_common_testing
    use numfort_common_workspace, workspace => workspace_real64

    implicit none

    integer, parameter :: PREC = real64

    call test_truncate_trans_matrix ()
    call test_truncate_trans_matrix_no_mass ()
    call test_tauchen ()
    call test_rouwenhorst_iid ()
    call test_rouwenhorst_python ()

    call test_histogram ()

    call test_simulate ()
    call test_simulate_advanced ()

    contains


subroutine test_truncate_trans_matrix ()
    !* - Check normal truncation, no-op tolerance, and threshold equality.
    !  - Check singleton and empty matrices.
    !  - Preserve retained probability ratios and normalize both orientations.
    !  - Avoid invalid operations and division by zero.

    real (PREC), parameter :: ATOL = 32.0_PREC * epsilon (1.0_PREC)
    real (PREC), parameter :: CONTROL(2,2) = reshape ([ &
        0.9_PREC, 0.1_PREC, 0.2_PREC, 0.8_PREC], [2,2], order=[2,1])
    real (PREC), parameter :: DESIRED(2,2) = reshape ([ &
        1.0_PREC, 0.0_PREC, 0.2_PREC, 0.8_PREC], [2,2], order=[2,1])
    real (PREC), parameter :: ASYMMETRIC(3,3) = reshape ([ &
        0.6_PREC,   0.3_PREC,   0.1_PREC, &
        0.125_PREC, 0.375_PREC, 0.5_PREC, &
        0.25_PREC,  0.0_PREC,   0.75_PREC], [3,3], order=[2,1])
    ! Retained odds are 2:1 and 3:4, so probabilities are 2/3, 1/3 and 3/7, 4/7.
    real (PREC), parameter :: RATIOS(3,3) = reshape ([ &
        2.0_PREC / 3.0_PREC, 1.0_PREC / 3.0_PREC, 0.0_PREC, &
        0.0_PREC, 3.0_PREC / 7.0_PREC, 4.0_PREC / 7.0_PREC, &
        0.25_PREC, 0.0_PREC, 0.75_PREC], [3,3], order=[2,1])
    real (PREC) :: mat(2,2), expected(2,2), mat3(3,3), singleton(1,1), empty(0,0)
    type (status_t) :: status
    logical :: flag_invalid, flag_div_zero

    call ieee_set_flag (ieee_invalid, .false.)
    call ieee_set_flag (ieee_divide_by_zero, .false.)

    ! The independent control removes only 0.1 and normalizes its row to [1,0].
    mat = CONTROL
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (mat, 0.15_PREC, status=status)
    if (status /= NF_STATUS_OK) error stop "Truncation rejected the normal control"
    if (.not. all_close (mat, DESIRED, rtol=0.0_PREC, atol=ATOL)) &
        error stop "Truncation differs from the normal control"
    if (.not. all (abs (sum (mat, dim=2) - 1.0_PREC) <= ATOL)) &
        error stop "Truncation did not normalize rows"

    ! Transposed input must give the transpose of the independent control.
    mat = transpose (CONTROL)
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (mat, 0.15_PREC, .true., status)
    if (status /= NF_STATUS_OK) error stop "Truncation rejected transposed control"
    if (.not. all_close (mat, transpose (DESIRED), rtol=0.0_PREC, atol=ATOL)) &
        error stop "Truncation differs from the transposed control"
    if (.not. all (abs (sum (mat, dim=1) - 1.0_PREC) <= ATOL)) &
        error stop "Truncation did not normalize columns"

    ! Zero tolerance preserves an already normalized matrix, including zeros.
    expected = reshape ([1.0_PREC, 0.0_PREC, 0.25_PREC, 0.75_PREC], &
        [2,2], order=[2,1])
    mat = expected
    call truncate_trans_matrix (mat, 0.0_PREC)
    if (any (mat /= expected)) error stop "Zero tolerance changed the matrix"

    ! Keep probabilities equal to the threshold; remove only the smaller 0.125.
    mat = reshape ([0.75_PREC, 0.25_PREC, 0.125_PREC, 0.875_PREC], &
        [2,2], order=[2,1])
    expected = reshape ([0.75_PREC, 0.25_PREC, 0.0_PREC, 1.0_PREC], &
        [2,2], order=[2,1])
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (mat, 0.25_PREC, status=status)
    if (status /= NF_STATUS_OK) error stop "Truncation rejected threshold equality"
    if (any (mat /= expected)) error stop "Truncation removed threshold equality"

    ! A singleton retains its unit probability at the exact threshold of one.
    singleton = 1.0_PREC
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (singleton, 1.0_PREC, status=status)
    if (status /= NF_STATUS_OK) error stop "Truncation rejected a singleton"
    if (any (singleton /= 1.0_PREC)) error stop "Truncation changed a singleton"

    ! Rows with different retained masses preserve their original odds.
    mat3 = ASYMMETRIC
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (mat3, 0.2_PREC, transposed=.false., status=status)
    if (status /= NF_STATUS_OK) error stop "Truncation rejected asymmetric rows"
    if (.not. all (ieee_is_finite (mat3))) &
        error stop "Truncation returned nonfinite row probabilities"
    if (.not. all_close (mat3, RATIOS, rtol=0.0_PREC, atol=ATOL)) &
        error stop "Truncation changed retained row probability ratios"
    if (.not. all (abs (sum (mat3, dim=2) - 1.0_PREC) <= ATOL)) &
        error stop "Truncation did not normalize asymmetric rows"

    ! Transposed asymmetric input preserves the same odds in each column.
    mat3 = transpose (ASYMMETRIC)
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (mat3, 0.2_PREC, transposed=.true., status=status)
    if (status /= NF_STATUS_OK) error stop "Truncation rejected asymmetric columns"
    if (.not. all (ieee_is_finite (mat3))) &
        error stop "Truncation returned nonfinite column probabilities"
    if (.not. all_close (mat3, transpose (RATIOS), rtol=0.0_PREC, atol=ATOL)) &
        error stop "Truncation changed retained column probability ratios"
    if (.not. all (abs (sum (mat3, dim=1) - 1.0_PREC) <= ATOL)) &
        error stop "Truncation did not normalize asymmetric columns"

    ! Empty matrices contain no conditional distributions to truncate.
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (empty, 0.6_PREC, status=status)
    if (status /= NF_STATUS_OK) error stop "Truncation rejected an empty matrix"

    call ieee_get_flag (ieee_invalid, flag_invalid)
    call ieee_get_flag (ieee_divide_by_zero, flag_div_zero)
    if (flag_invalid .or. flag_div_zero) &
        error stop "Successful truncation raised an IEEE exception"

end subroutine


subroutine test_truncate_trans_matrix_no_mass ()
    !* - Reject truncation that leaves any conditional distribution without mass.
    !  - Preserve the entire input, including when status is omitted.
    !  - Avoid NaNs, invalid operations, and division by zero.

    real (PREC), parameter :: MIXED(2,2) = reshape ([ &
        0.9_PREC, 0.1_PREC, 0.5_PREC, 0.5_PREC], [2,2], order=[2,1])
    real (PREC) :: mat(2,2), singleton(1,1)
    type (status_t) :: status
    logical :: flag_invalid, flag_div_zero

    call ieee_set_flag (ieee_invalid, .false.)
    call ieee_set_flag (ieee_divide_by_zero, .false.)

    ! Issue #43: truncation at 0.6 would remove every entry of the uniform matrix.
    mat = 0.5_PREC
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (mat, 0.6_PREC, status=status)
    call ieee_get_flag (ieee_invalid, flag_invalid)
    call ieee_get_flag (ieee_divide_by_zero, flag_div_zero)
    if (flag_invalid .or. flag_div_zero) &
        error stop "All-removed truncation raised an IEEE exception"
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Truncation accepted removal of every probability"
    if (.not. all (ieee_is_finite (mat))) &
        error stop "All-removed truncation returned nonfinite probabilities"
    if (any (mat /= 0.5_PREC)) error stop "All-removed truncation changed the matrix"

    ! Existing callers without status receive the same unchanged matrix.
    mat = 0.5_PREC
    call truncate_trans_matrix (mat, 0.6_PREC)
    if (any (mat /= 0.5_PREC)) &
        error stop "All-removed truncation without status changed the matrix"

    ! Only the second row loses all mass; neither row may be modified.
    mat = MIXED
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (mat, 0.6_PREC, status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Truncation accepted removal of one row"
    if (any (mat /= MIXED)) error stop "Rejected truncation partially changed rows"

    ! Only the second column loses all mass in the transposed fixture.
    mat = transpose (MIXED)
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (mat, 0.6_PREC, transposed=.true., status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Truncation accepted removal of one column"
    if (any (mat /= transpose (MIXED))) &
        error stop "Rejected truncation partially changed columns"

    ! A threshold above one would remove the singleton's only transition.
    singleton = 1.0_PREC
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (singleton, 1.1_PREC, status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Truncation accepted removal of the singleton transition"
    if (any (singleton /= 1.0_PREC)) &
        error stop "Rejected truncation changed a singleton"

    ! Retained zero entries provide no positive mass even with zero tolerance.
    mat = 0.0_PREC
    status = NF_STATUS_UNDEFINED
    call truncate_trans_matrix (mat, 0.0_PREC, status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Truncation accepted a zero-mass matrix"
    if (any (mat /= 0.0_PREC)) error stop "Truncation changed a zero-mass matrix"

    call ieee_get_flag (ieee_invalid, flag_invalid)
    call ieee_get_flag (ieee_divide_by_zero, flag_div_zero)
    if (flag_invalid .or. flag_div_zero) &
        error stop "Rejected truncation raised an IEEE exception"

end subroutine


subroutine test_tauchen ()
    !*  TEST_TAUCHEN compares discretized moments of Markov chain generated
    !   using the TAUCHEN method to those tabulated in Tauchen (1986).

    real (PREC) :: rho, sigma
    real (PREC) :: rho_m, sigma_m
    real (PREC), dimension(:), allocatable :: states
    real (PREC), dimension(:,:), allocatable :: tm
    integer :: n
    type (status_t) :: status

    ! Replicate results reported in Tauchen (1986) paper
    rho = 0.10d0
    sigma = 0.101
    n = 9

    allocate (states(n), tm(n,n))
    call tauchen (rho, sigma, states, tm, sigma_cond=.false., status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    if (.not. (all_close (rho_m, 0.1d0, atol=1.0d-3) .and. all_close (sigma_m, &
        0.103d0, atol=1.0d-3))) &
        error stop "Tauchen moments differ from the reference"
    deallocate (states, tm)

    rho = 0.80d0
    sigma = 0.167
    n = 9

    allocate (states(n), tm(n,n))
    call tauchen (rho, sigma, states, tm, sigma_cond=.false., status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    if (.not. (all_close (rho_m, 0.798d0, atol=1.0d-3) .and. all_close (sigma_m, &
        0.176d0, atol=1.0d-3))) &
        error stop "Tauchen moments differ from the reference"
    deallocate (states, tm)

    rho = 0.90d0
    sigma = 0.229
    n = 9

    allocate (states(n), tm(n,n))
    call tauchen (rho, sigma, states, tm, sigma_cond=.false., status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    if (.not. (all_close (rho_m, 0.898d0, atol=1.0d-3) .and. all_close (sigma_m, &
        0.253d0, atol=1.0d-3))) &
        error stop "Tauchen moments differ from the reference"
    deallocate (states, tm)

    rho = 0.90d0
    sigma = 0.229
    n = 5

    allocate (states(n), tm(n,n))
    call tauchen (rho, sigma, states, tm, sigma_cond=.false., status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    if (.not. (all_close (rho_m, 0.932d0, atol=1.0d-3) .and. all_close (sigma_m, &
        0.291d0, atol=1.0d-3))) &
        error stop "Tauchen moments differ from the reference"
    deallocate (states, tm)

    ! Test some degenerate processes
    rho = 0.0d0
    sigma = sqrt(0.0522d0)
    n = 5

    allocate (states(n), tm(n,n))
    call tauchen (rho, sigma, states, tm, status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    if (.not. (all_close (rho_m, 0.0d0, atol=1.0d-8))) &
        error stop "Tauchen returned nonzero autocorrelation for independent states"
    deallocate (states, tm)

end subroutine

subroutine test_rouwenhorst_iid ()

    real (PREC) :: rho, sigma
    real (PREC) :: rho_m, sigma_m
    real (PREC), dimension(:), allocatable :: states
    real (PREC), dimension(:,:), allocatable :: tm
    integer :: n
    type (status_t) :: status

    rho = 0.0d0
    sigma = sqrt(0.0522d0)
    n = 5
    allocate (states(n), tm(n,n))
    call rouwenhorst (rho, sigma, states, tm, status=status)
    call markov_moments (states, tm, acorr=rho_m, sigma=sigma_m, status=status)
    if (.not. (all_close (rho_m, 0.0d0, atol=1.0d-8))) &
        error stop "Rouwenhorst returned nonzero autocorrelation for independent " // &
            "states"
    deallocate (states, tm)

end subroutine

subroutine test_rouwenhorst_python ()

    real (PREC) :: rho, sigma
    real (PREC), dimension(:), allocatable :: states, edist, edist2
    real (PREC), dimension(:,:), allocatable :: tm, tm_T
    type (workspace) :: work
    integer :: n
    type (status_t) :: status

    integer, parameter :: N2 = 9
    real (PREC), parameter :: TM2_desired(*,*) = reshape([ &
        6.63420431d-01,   2.79334918d-01,   5.14564323d-02,   5.41646656d-03, &
        3.56346484d-04,   1.50040625d-05,   3.94843750d-07,   5.93750000d-09, &
        3.90625000d-11,   3.49168648d-02,   6.76284539d-01,   2.46449229d-01, &
        3.87704975d-02,   3.39466914d-03,   1.78469375d-04,   5.63171875d-06, &
        9.87500000d-08,   7.42187500d-10,   1.83772973d-03,   7.04140653d-02, &
        6.85549548d-01,   2.12408226d-01,   2.77697840d-02,   1.94249469d-03, &
        7.65292188d-05,   1.60906250d-06,   1.41015625d-08,   9.67226172d-05, &
        5.53864250d-03,   1.06204113d-01,   6.91139238d-01,   1.77494044d-01, &
        1.85302288d-02,   9.71247344d-04,   2.54956250d-05,   2.67929688d-07, &
        5.09066406d-06,   3.87962188d-04,   1.11079136d-02,   1.41995235d-01, &
        6.93007596d-01,   1.41995235d-01,   1.11079136d-02,   3.87962188d-04, &
        5.09066406d-06,   2.67929688d-07,   2.54956250d-05,   9.71247344d-04, &
        1.85302288d-02,   1.77494044d-01,   6.91139238d-01,   1.06204113d-01, &
        5.53864250d-03,   9.67226172d-05,   1.41015625d-08,   1.60906250d-06, &
        7.65292188d-05,   1.94249469d-03,   2.77697840d-02,   2.12408226d-01, &
        6.85549548d-01,   7.04140653d-02,   1.83772973d-03,   7.42187500d-10, &
        9.87500000d-08,   5.63171875d-06,   1.78469375d-04,   3.39466914d-03, &
        3.87704975d-02,   2.46449229d-01,   6.76284539d-01,   3.49168648d-02, &
        3.90625000d-11,   5.93750000d-09,   3.94843750d-07,   1.50040625d-05, &
        3.56346484d-04,   5.41646656d-03,   5.14564323d-02,   2.79334918d-01, &
        6.63420431d-01], &
        shape=[N2,N2], order=[2, 1])

    real (PREC), parameter :: states2_desired(*) = [ &
        -0.64888568d0, -0.48666426d0, -0.32444284d0, -0.16222142d0, 0.0d0,  &
         0.16222142d0,  0.32444284d0,  0.48666426d0,  0.64888568d0]

    real (PREC), parameter :: EDIST2_DESIRED(*) = [ &
        0.00390625d0, 0.03125d0, 0.109375d0, 0.21875d0, 0.2734375d0, &
        0.21875d0, 0.109375d0, 0.03125d0, 0.00390625d0 ]

    ! Compare to results from Python implementation for N=9, rho=0.9, sigma=0.1
    rho = 0.9d0
    sigma = 0.1d0
    n = N2

    allocate (states(n), tm(n,n), tm_T(n,n), edist(n), edist2(n))

    call rouwenhorst (rho, sigma, states, tm, status=status)

    if (.not. (status == NF_STATUS_OK)) &
        error stop "Rouwenhorst returned an error status"

    if (.not. (all_close (tm, TM2_desired, rtol=1.0d-7))) &
        error stop "Transition matrix differs from the Python reference"

    if (.not. (all_close (states, states2_desired, atol=1.0d-7))) &
        error stop "State grid differs from the Python reference"

    call markov_ergodic_dist (tm, edist, inverse=.true., work=work, status=status)

    if (.not. (status == NF_STATUS_OK)) &
        error stop "Ergodic distribution solver returned an error status"

    if (.not. (all_close (edist, EDIST2_DESIRED, atol=1.0d-4))) &
        error stop "Ergodic distribution differs from the Python reference"

    tm_T(:,:) = transpose(tm)
    status = NF_STATUS_UNDEFINED
    call markov_ergodic_dist (tm_T, edist2, inverse=.true., transposed=.true., &
        work=work, status=status)

    if (.not. (status == NF_STATUS_OK)) &
        error stop "Ergodic distribution solver failed with transposed input"

    if (.not. (all_close (edist2, edist, atol=1.0d-10, rtol=0.0d0))) &
        error stop "Ergodic distribution differs for transposed input"

    edist2(:) = 0.0
    status = NF_STATUS_UNDEFINED
    call markov_ergodic_dist (tm, edist2, inverse=.false., work=work, status=status)

    if (.not. (status == NF_STATUS_OK)) &
        error stop "Direct ergodic distribution solver returned an error status"

    if (.not. (all_close (edist, edist, atol=1.0d-10, rtol=0.0d0))) &
        error stop "Direct ergodic distribution comparison failed"

    edist2(:) = 0.0
    status = NF_STATUS_UNDEFINED
    call markov_ergodic_dist (tm_T, edist2, inverse=.false., transposed=.true., &
        work=work, status=status)

    if (.not. (status == NF_STATUS_OK)) &
        error stop "Direct ergodic solver failed with transposed input"

    if (.not. (all_close (edist, edist, atol=1.0d-10, rtol=0.0d0))) &
        error stop "Transposed direct ergodic distribution comparison failed"

    deallocate (states, tm, tm_T, edist, edist2)

end subroutine

subroutine test_histogram ()

    type (status_t) :: status

    real (PREC), dimension(:,:), allocatable :: transm
    real (PREC), dimension(:), allocatable :: states
    integer, dimension(:,:), allocatable :: trans_hist
    real (PREC) :: rho, sigma
    integer :: nobs, nstates

    rho = 0.9
    sigma = 0.2
    nobs = 10000
    nstates = 7

    allocate (transm(nstates, nstates))
    allocate (states(nstates))

    ! Create Rouwenhorst approximation
    call rouwenhorst (rho, sigma, states, transm, status=status)

    ! Create discretized histogram
    allocate (trans_hist(nstates,nstates))

    call trans_histogram (nobs, transm, trans_hist, status=status)
    if (.not. (status == NF_STATUS_OK .and. all (sum (trans_hist, dim=1) == sum &
        (trans_hist, dim=2)) .and. sum (trans_hist) == nobs)) &
        error stop "Transition histogram failed the count conservation check"

end subroutine

subroutine test_simulate ()

    real (PREC), dimension(:,:), allocatable :: transm, transm_T
    real (PREC), dimension(:), allocatable :: states
    integer, dimension(:), allocatable :: sim, sim_T
    integer :: n, Nobs
    integer :: in, ir, is
    real (PREC) :: rho, sigma
    real (PREC), parameter :: RHO_ALL(*) = [ real (PREC) :: 0.0, 0.2, 0.5, 0.9]
    real (PREC), parameter :: SIGMA_ALL(*) = [ real (PREC) :: 0.01, 0.1, 0.2]
    integer, parameter :: N_ALL(*) = [3, 5, 7]
    logical :: transposed, values_ok
    type (status_t) :: status, status1, status2

    ! === Input checks ===

    ! 1. Invalid transition matrix
    Nobs = 100
    allocate (sim(Nobs))
    n = 2
    allocate (transm(n,n), source=0.0_PREC)
    allocate (transm_T(n,n), source=0.0_PREC)

    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, 1, sim, transposed=.false., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Markov simulation invalid-matrix check failed"

    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, 1, sim, transposed=.true., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Markov simulation invalid-matrix check failed with transposed " &
            // "input"

    transm(:,:) = 1.0
    transm_T(:,:) = 1.0
    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, 1, sim, transposed=.false., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Markov simulation invalid-matrix check failed"

    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, 1, sim, transposed=.true., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Markov simulation invalid-matrix check failed with transposed " &
            // "input"

    ! 2. Invalid initial value
    transm(:,:) = 1.0_PREC / n
    transm_T(:,:) = 1.0_PREC / n
    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, 0, sim, transposed=.false., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Markov simulation invalid-state check failed"

    status = NF_STATUS_UNDEFINED
    call markov_simulate (transm, n+1, sim, transposed=.false., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Markov simulation invalid-state check failed"

    deallocate (transm, transm_T, sim)

    ! === Compare transposed vs. non-transposed TM input ===

    Nobs = 10000
    allocate (sim(Nobs), sim_T(Nobs))

    do in = 1, size(N_ALL)

        n = N_ALL(in)

        allocate (transm(n,n), transm_T(n,n), states(n))

        do ir = 1, size(RHO_ALL)
            do is = 1, size(SIGMA_ALL)
                rho = RHO_ALL(ir)
                sigma = SIGMA_ALL(is)

                call rouwenhorst (rho, sigma, states, transm)
                transm_T(:,:) = transpose(transm)

                call set_seed (1234)
                status1 = NF_STATUS_UNDEFINED
                sim(:) = 0
                call markov_simulate (transm, 1, sim, transposed=.false., &
                    status=status1)

                call set_seed (1234)
                status2 = NF_STATUS_UNDEFINED
                sim_T(:) = 0
                call markov_simulate (transm_T, 1, sim_T, transposed=.true., &
                    status=status2)

                values_ok = all(sim == sim_T)

                if (.not. (values_ok .and. status1 == NF_STATUS_OK .and. status2 == &
                    NF_STATUS_OK)) &
                    error stop "Markov simulations differ for transposed and " // &
                        "non-transposed inputs"

            end do
        end do

        deallocate (transm, transm_T, states)
    end do

end subroutine

subroutine test_simulate_advanced ()

    real (PREC), dimension(:,:), allocatable :: transm, transm_T
    real (PREC), dimension(:,:), allocatable :: pmf_trans, pmf_sample
    real (PREC), dimension(:), allocatable :: states, edist
    integer, dimension(:,:), allocatable :: ntrans
    integer, dimension(:), allocatable :: sim, sim_T
    integer :: n, Nobs
    integer :: in, ir, is, i
    real (PREC) :: rho, sigma
    real (PREC), parameter :: RHO_ALL(*) = [ real (PREC) :: 0.0, 0.2d0, 0.5d0, 0.9d0]
    real (PREC), parameter :: SIGMA_ALL(*) = [ real (PREC) :: 0.01d0, 0.1d0, 0.2d0]
    integer, parameter :: N_ALL(*) = [3, 5]
    real (PREC), parameter :: rtol=5.0e-2_PREC, atol=1.0e-3_PREC
    logical :: transposed, values_ok
    type (status_t) :: status, status1, status2

    ! === Input checks ===

    ! 1. Invalid transition matrix
    Nobs = 100
    allocate (sim(Nobs))
    n = 2
    allocate (transm(n,n), source=0.0_PREC)
    allocate (transm_T(n,n), source=0.0_PREC)

    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, 1, sim, transposed=.false., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Advanced Markov simulation invalid-matrix check failed"

    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, 1, sim, transposed=.true., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Advanced Markov simulation invalid-matrix check failed with " // &
            "transposed input"

    transm(:,:) = 1.0
    transm_T(:,:) = 1.0
    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, 1, sim, transposed=.false., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Advanced Markov simulation invalid-matrix check failed"

    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, 1, sim, transposed=.true., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Advanced Markov simulation invalid-matrix check failed with " // &
            "transposed input"

    ! 2. Invalid initial value
    transm(:,:) = 1.0_PREC / n
    transm_T(:,:) = 1.0_PREC / n
    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, 0, sim, transposed=.false., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Advanced Markov simulation invalid-state check failed"

    status = NF_STATUS_UNDEFINED
    call markov_simulate_advanced (transm, n+1, sim, transposed=.false., status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Advanced Markov simulation invalid-state check failed"

    deallocate (transm, transm_T, sim)

    ! === Compare transposed vs. non-transposed TM input ===

    Nobs = 10000
    allocate (sim(Nobs), sim_T(Nobs))

    do in = 1, size(N_ALL)

        n = N_ALL(in)

        allocate (transm(n,n), transm_T(n,n), states(n), edist(n))
        allocate (pmf_trans(n,n), pmf_sample(n,n))
        allocate (ntrans(n,n))

        do ir = 1, size(RHO_ALL)
            do is = 1, size(SIGMA_ALL)
                rho = RHO_ALL(ir)
                sigma = SIGMA_ALL(is)

                call rouwenhorst (rho, sigma, states, transm)
                transm_T(:,:) = transpose(transm)

                call set_seed (1234)
                status1 = NF_STATUS_UNDEFINED
                sim(:) = 0
                call markov_simulate_advanced (transm, 1, sim, rtol=rtol, &
                    atol=atol, transposed=.false., status=status1)

                call set_seed (1234)
                status2 = NF_STATUS_UNDEFINED
                sim_T(:) = 0
                call markov_simulate_advanced (transm_T, 1, sim_T, rtol=rtol, &
                    atol=atol, transposed=.true., status=status2)

                values_ok = all(sim == sim_T)

                if (.not. (values_ok .and. status1 == NF_STATUS_OK .and. status2 == &
                    NF_STATUS_OK)) &
                    error stop "Advanced Markov simulations differ across " // &
                        "transpose modes"

                ! Check that generated sequence is close to actual
                ! transition matrix
                call count_trans (sim, ntrans)

                ! Theoretical transition PMF
                call markov_ergodic_dist (transm, edist, inverse=.true.)
                do i = 1, n
                    pmf_trans(i,:) = transm(i,:) * edist(i)
                end do

                ! Transition PMF for simulated sample
                do i = 1, n
                    pmf_sample(i,:) = ntrans(i,:) / real(sum(ntrans), PREC)
                end do

                ! This is not exactly the terminal condition that simulation
                ! routine computes, so we need to less restrictive tolerance
                values_ok = all_close (pmf_sample, pmf_trans, &
                    atol=atol + 1.0e-3_PREC, rtol=rtol)

                if (.not. (values_ok)) &
                    error stop "Advanced Markov simulation has unexpected " // &
                        "transition frequencies"

            end do
        end do

        deallocate (transm, transm_T, ntrans, states, edist)
        deallocate (pmf_trans, pmf_sample)
    end do

end subroutine

pure subroutine count_trans (x, num)
    !*  COUNT_TRANS returns the number of observed transitions for each
    !   transition type.
    integer, intent(in), dimension(:), contiguous :: x
    integer, intent(out), dimension(:,:), contiguous :: num
    !*  Matrix containing the number of transitions for each possible
    !   combination of origin and destrination values.

    integer :: i, n
    integer :: ifrom,  ito

    n = size(x)
    num = 0

    do i = 1, n - 1
        ifrom = x(i)
        ito = x(i+1)
        num(ifrom, ito) = num(ifrom,ito) + 1
    end do

end subroutine

end program
