program test_stats_markov_ergodic_dist

    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    use numfort_common
    use numfort_common_workspace, only: workspace => workspace_real64
    use numfort_stats, only: markov_ergodic_dist

    implicit none

    integer, parameter :: PREC = real64
    real (PREC), parameter :: TOL = 64.0_PREC * epsilon (1.0_PREC)
    real (PREC), parameter :: SENTINEL = -7.0_PREC

    ! Asymmetric chain with transition probabilities a = 1/4 and b = 1/2.
    real (PREC), parameter :: TM_TWO(2,2) = reshape ([ &
        0.75_PREC, 0.25_PREC, 0.5_PREC, 0.5_PREC], [2, 2], order=[2, 1])
    ! Balance gives pi = [b, a] / (a + b) = [2/3, 1/3].
    real (PREC), parameter :: PI_TWO(2) = [2.0_PREC, 1.0_PREC] / 3.0_PREC

    call test_ergodic_dist_solutions ()
    call test_ergodic_dist_singular ()
    call test_ergodic_dist_invalid ()
    call test_ergodic_dist_workspace ()

    contains

subroutine test_ergodic_dist_solutions ()
    !* - Check independent stationary distributions with internal workspace.
    !  - Cover default options, transposed input, and a periodic chain.

    real (PREC) :: edist(2), singleton(1), tm_single(1,1), tm_cycle(2,2)
    real (PREC) :: tm_transposed(2,2)
    type (status_t) :: status

    ! Default inverse method and orientation, with no caller workspace.
    edist = SENTINEL
    call markov_ergodic_dist (TM_TWO, edist, status=status)
    call check_distribution (TM_TWO, edist, PI_TWO, status)

    ! Explicit inverse method with the asymmetric matrix stored by columns.
    tm_transposed = transpose (TM_TWO)
    edist = SENTINEL
    call markov_ergodic_dist (tm_transposed, edist, inverse=.true., &
        transposed=.true., status=status)
    call check_distribution (TM_TWO, edist, PI_TWO, status)

    ! A singleton places all probability mass on its only state.
    tm_single = 1.0_PREC
    singleton = SENTINEL
    call markov_ergodic_dist (tm_single, singleton, status=status)
    call check_distribution (tm_single, singleton, [1.0_PREC], status)

    ! Alternation has equal stationary weights even though iteration can cycle.
    tm_cycle = reshape ([0.0_PREC, 1.0_PREC, 1.0_PREC, 0.0_PREC], [2, 2])
    edist = SENTINEL
    call markov_ergodic_dist (tm_cycle, edist, inverse=.true., status=status)
    call check_distribution (tm_cycle, edist, [0.5_PREC, 0.5_PREC], status)
end subroutine


subroutine test_ergodic_dist_singular ()
    !* - Reject nonunique stationary distributions in the inverse method.
    !  - Preserve output on failure and allow subsequent internal allocations.

    real (PREC) :: tm_identity(2,2), tm_reducible(3,3), edist(2), edist3(3)
    type (status_t) :: status

    ! Both identity states are absorbing, so every distribution is stationary.
    tm_identity = reshape ([1.0_PREC, 0.0_PREC, 0.0_PREC, 1.0_PREC], [2, 2])
    edist = SENTINEL
    call markov_ergodic_dist (tm_identity, edist, status=status)
    if (status /= NF_STATUS_INVALID_STATE) &
        error stop "Identity inverse system did not report a singular matrix"
    if (any (edist /= SENTINEL)) &
        error stop "Singular inverse system changed the output"

    ! Two absorbing states and one transient state also give a singular system.
    tm_reducible = reshape ([ &
        1.0_PREC, 0.0_PREC, 0.0_PREC, &
        0.0_PREC, 1.0_PREC, 0.0_PREC, &
        0.25_PREC, 0.75_PREC, 0.0_PREC], [3, 3], order=[2, 1])
    edist3 = SENTINEL
    call markov_ergodic_dist (tm_reducible, edist3, status=status)
    if (status /= NF_STATUS_INVALID_STATE) &
        error stop "Reducible inverse system did not report a singular matrix"
    if (any (edist3 /= SENTINEL)) &
        error stop "Reducible inverse system changed the output"

    ! A valid call after singular failures must replace the failure status.
    call markov_ergodic_dist (TM_TWO, edist, status=status)
    call check_distribution (TM_TWO, edist, PI_TWO, status)
end subroutine


subroutine test_ergodic_dist_invalid ()
    !* - Reject output lengths and matrix dimensions before using LAPACK.
    !  - Preserve output on each argument error without caller workspace.

    real (PREC) :: short_output(1), long_output(3), edist(2), edist3(3)
    real (PREC) :: wide(2,3), tall(3,2), empty(0,0), empty_output(0)
    real (PREC) :: no_columns(2,0)

    ! The output is shorter than the number of states.
    call check_invalid_argument (TM_TWO, short_output)

    ! The output is longer than the number of states.
    call check_invalid_argument (TM_TWO, long_output)

    ! A matrix with extra columns must fail even when its rows sum to one.
    wide = 1.0_PREC / 3.0_PREC
    call check_invalid_argument (wide, edist)

    ! A matrix with missing columns must fail before the transpose assignment.
    tall = 0.5_PREC
    call check_invalid_argument (tall, edist3)

    ! Transposed storage does not make a wide matrix valid.
    call check_invalid_argument (wide, edist, transposed=.true.)

    ! Transposed storage does not make a tall matrix valid.
    call check_invalid_argument (tall, edist3, transposed=.true.)

    ! An empty chain has no normalized stationary distribution.
    call check_invalid_argument (empty, empty_output)

    ! Nonempty rows with no columns must also be rejected.
    call check_invalid_argument (no_columns, edist)
end subroutine


subroutine test_ergodic_dist_workspace ()
    !* - Reuse caller workspace while growing and shrinking the state space.
    !  - Retain caller allocations after invalid and singular calls.
    !  - Recover from failures without stale workspace offsets or status.

    type (workspace) :: work
    type (status_t) :: status
    real (PREC) :: tm_single(1,1), singleton(1), tm_iid(3,3), edist3(3)
    real (PREC) :: tm_identity(2,2), edist(2), wide(2,3)
    integer :: nrwrk, niwrk, i
    ! Identical rows imply pi = [1/2, 1/3, 1/6] after one transition.
    real (PREC), parameter :: PI_IID(3) = [3.0_PREC, 2.0_PREC, 1.0_PREC] / 6.0_PREC

    ! Invalid output with fresh workspace must not allocate working arrays.
    singleton = SENTINEL
    call markov_ergodic_dist (TM_TWO, singleton, work=work, status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Invalid output length with caller workspace was accepted"
    if (any (singleton /= SENTINEL)) &
        error stop "Invalid output length changed the output"
    if (allocated (work%rwrk) .or. allocated (work%iwrk)) &
        error stop "Invalid output length allocated caller workspace"

    ! Start with one state and unrelated caller-owned workspace components.
    allocate (work%lwrk(2), source=.true.)
    work%cwrk = "caller data"
    tm_single = 1.0_PREC
    call markov_ergodic_dist (tm_single, singleton, work=work, status=status)
    call check_distribution (tm_single, singleton, [1.0_PREC], status)
    if (.not. allocated (work%rwrk) .or. .not. allocated (work%iwrk)) &
        error stop "Successful call did not retain caller workspace"

    ! Grow the same workspace to fit three independent states.
    do i = 1, 3
        tm_iid(i,:) = PI_IID
    end do
    edist3 = SENTINEL
    call markov_ergodic_dist (tm_iid, edist3, work=work, status=status)
    call check_distribution (tm_iid, edist3, PI_IID, status)
    nrwrk = size (work%rwrk)
    niwrk = size (work%iwrk)

    ! Shrink the problem while retaining the larger workspace allocations.
    edist = SENTINEL
    call markov_ergodic_dist (TM_TWO, edist, work=work, status=status)
    call check_distribution (TM_TWO, edist, PI_TWO, status)

    ! A singular inverse system must leave caller storage allocated and usable.
    tm_identity = reshape ([1.0_PREC, 0.0_PREC, 0.0_PREC, 1.0_PREC], [2, 2])
    edist = SENTINEL
    call markov_ergodic_dist (tm_identity, edist, work=work, status=status)
    if (status /= NF_STATUS_INVALID_STATE) &
        error stop "Singular system with caller workspace was accepted"
    if (any (edist /= SENTINEL)) &
        error stop "Singular system with caller workspace changed the output"
    if (.not. allocated (work%rwrk) .or. .not. allocated (work%iwrk)) &
        error stop "Singular system deallocated caller workspace"
    if (size (work%rwrk) /= nrwrk .or. size (work%iwrk) /= niwrk) &
        error stop "Smaller problems resized caller workspace"

    ! Invalid shape must leave previously allocated workspace contents intact.
    wide = 1.0_PREC / 3.0_PREC
    work%rwrk = SENTINEL
    work%iwrk = -7
    call markov_ergodic_dist (wide, edist, work=work, status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Invalid matrix shape with caller workspace was accepted"
    if (any (edist /= SENTINEL)) &
        error stop "Invalid matrix shape changed the output"
    if (.not. allocated (work%rwrk) .or. .not. allocated (work%iwrk)) &
        error stop "Invalid matrix shape deallocated caller workspace"
    if (any (work%rwrk /= SENTINEL) .or. any (work%iwrk /= -7)) &
        error stop "Invalid matrix shape changed caller workspace contents"

    ! Reuse the workspace after both failure paths, with a smaller valid chain.
    singleton = SENTINEL
    call markov_ergodic_dist (tm_single, singleton, work=work, status=status)
    call check_distribution (tm_single, singleton, [1.0_PREC], status)
    if (.not. allocated (work%lwrk) .or. .not. allocated (work%cwrk)) &
        error stop "Ergodic solver deallocated unrelated caller storage"
    if (.not. all (work%lwrk) .or. work%cwrk /= "caller data") &
        error stop "Ergodic solver changed unrelated caller storage"
end subroutine


subroutine check_distribution (tm, edist, expected, status)
    !* Check the independent reference and stationary probability identities.

    real (PREC), intent(in) :: tm(:,:)
        !* Transition probabilities stored by rows.
    real (PREC), intent(in) :: edist(:)
        !* Computed stationary distribution.
    real (PREC), intent(in) :: expected(:)
        !* Independently derived stationary probabilities.
    type (status_t), intent(in) :: status
        !* Solver exit status.

    if (status /= NF_STATUS_OK) &
        error stop "Inverse ergodic solver returned an error status"
    if (.not. all (ieee_is_finite (edist))) &
        error stop "Inverse ergodic solver returned nonfinite probabilities"
    if (any (abs (edist - expected) > TOL)) &
        error stop "Stationary distribution differs from independent solution"
    if (any (edist < 0.0_PREC)) &
        error stop "Stationary distribution contains negative probabilities"
    if (abs (sum (edist) - 1.0_PREC) > TOL) &
        error stop "Stationary distribution is not normalized"
    ! Stationarity requires pi * P = pi for a row-stochastic transition matrix.
    if (any (abs (matmul (edist, tm) - edist) > TOL)) &
        error stop "Stationary distribution does not satisfy pi * P = pi"
end subroutine


subroutine check_invalid_argument (tm, edist, transposed)
    !* Check argument rejection and unchanged output with internal workspace.

    real (PREC), intent(in), contiguous :: tm(:,:)
        !* Transition matrix passed to the solver.
    real (PREC), intent(out), contiguous :: edist(:)
        !* Output initialized to a sentinel before the invalid call.
    logical, intent(in), optional :: transposed
        !* Whether the input uses transposed storage.

    type (status_t) :: status

    edist = SENTINEL
    call markov_ergodic_dist (tm, edist, transposed=transposed, status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Invalid ergodic solver argument did not report INVALID_ARG"
    if (any (edist /= SENTINEL)) &
        error stop "Invalid ergodic solver argument changed the output"
end subroutine

end program
