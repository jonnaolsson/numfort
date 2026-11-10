program test_minimize_dfls_failure
    !* Check least-squares results and process survival through the public API.

    use, intrinsic :: iso_fortran_env, only: real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use numfort_common
#if TEST_REAL_KIND == 32
    use numfort_optimize, only: minimize_dfls, optim_result => optim_result_real32
#else
    use numfort_optimize, only: minimize_dfls, optim_result => optim_result_real64
#endif

    implicit none

#if TEST_REAL_KIND == 32
    integer, parameter :: PREC = real32
    real (PREC), parameter :: LARGE_WEIGHT = 2.0_PREC**8
    real (PREC), parameter :: SMALL_WEIGHT = 2.0_PREC**(-4)
    real (PREC), parameter :: RHOEND = 1.0e-4_PREC
#else
    integer, parameter :: PREC = real64
    real (PREC), parameter :: LARGE_WEIGHT = 2.0_PREC**16
    real (PREC), parameter :: SMALL_WEIGHT = 2.0_PREC**(-12)
    real (PREC), parameter :: RHOEND = 1.0e-8_PREC
#endif
    integer, parameter :: MAXFEV = 200
    real (PREC), parameter :: XMIN(3) = [1.0_PREC, 1.0_PREC, 0.0_PREC]
    integer :: calls
    real (PREC) :: weights(4)
    real (PREC) :: points(3,MAXFEV), residuals(4,MAXFEV), norms(MAXFEV)

    call test_well_conditioned ()
    call test_ill_conditioned ()
    print "(a)", "DFLS failure returned to caller"

contains

subroutine test_well_conditioned ()
    !* Check the analytical minimum with equal residual weights.
    !  - Four affine residuals in three variables, starting at zero.

    real (PREC) :: x(3)
    type (optim_result) :: res

    ! The first two residuals vanish at x1=x2=1-x3. The remaining squared
    ! residuals sum to (x3-1)**2+(x3+1)**2=2*x3**2+2, minimized at x3=0.
    weights = 1.0_PREC
    x = 0.0_PREC
    calls = 0
    call minimize_dfls (affine_residual, x, 4, 0.25_PREC, RHOEND, &
        maxfev=MAXFEV, res=res)
    if (.not. res%success .or. res%status /= NF_STATUS_OK) &
        error stop "Well-conditioned affine convergence"
    if (any(abs(x - XMIN) > sqrt(epsilon(1.0_PREC)))) &
        error stop "Well-conditioned affine minimizer"
    if (abs(sum(res%fx**2) - 2.0_PREC) > sqrt(epsilon(1.0_PREC))) &
        error stop "Well-conditioned affine minimum"
    if (res%nfev /= calls) error stop "Well-conditioned callback count"
end subroutine


subroutine test_ill_conditioned ()
    !* Check that poorly scaled residuals return a usable result to the caller.
    !  - Unequal weights with an optimization result.
    !  - The same objective with the optional result omitted.

    real (PREC) :: x(3), minimum
    type (optim_result) :: res
    integer :: run, i, sample

    ! Positive weights preserve the unique minimizer [1,1,0]. Their ratio
    ! is 2**12 in real32 and 2**28 in real64, making the problem ill conditioned.
    weights = [LARGE_WEIGHT, LARGE_WEIGHT, SMALL_WEIGHT, SMALL_WEIGHT]
    ! At the minimum, only the residuals -SMALL_WEIGHT and SMALL_WEIGHT remain.
    minimum = 2.0_PREC * SMALL_WEIGHT**2
    do run = 1, 2
        x = 0.0_PREC
        calls = 0
        if (run == 1) then
            ! Case 1: A numerical failure must return a consistent saved result.
            res%success = .true.
            call minimize_dfls (affine_residual, x, 4, 0.25_PREC, RHOEND, &
                maxfev=MAXFEV, res=res)
        else
            ! Case 2: Omitting RES must also return to this caller.
            call minimize_dfls (affine_residual, x, 4, 0.25_PREC, RHOEND, &
                maxfev=MAXFEV)
        end if
        if (calls < 1 .or. calls > MAXFEV) error stop "Ill-conditioned budget"
        if (.not. all(ieee_is_finite(x))) error stop "Ill-conditioned coordinates"

        ! Compare with actual callback observations, without assuming a step
        ! sequence, an iteration count, or a particular choice among tied minima.
        sample = 0
        do i = 1, calls
            if (all(x == points(:,i))) sample = i
        end do
        if (sample == 0) error stop "Returned point was never evaluated"
        if (norms(sample) > minval(norms(:calls))) &
            error stop "Returned point must retain the best evaluation"

        if (run == 1) then
            if (res%nfev /= calls) error stop "Ill-conditioned callback count"
            if (any(res%x /= x) .or. any(res%fx /= residuals(:,sample))) &
                error stop "Ill-conditioned result correspondence"
            if (res%success) then
                ! Convergence is acceptable if the analytical minimum is reached
                ! to relative objective accuracy 1e-3; failure is not required.
                if (res%status /= NF_STATUS_OK) error stop "Convergence status"
                if (abs(sum(res%fx**2) - minimum) > 1.0e-3_PREC * minimum) &
                    error stop "Ill-conditioned false convergence"
            else
                if (res%status /= NF_STATUS_INVALID_STATE) &
                    error stop "Ill-conditioned numerical failure status"
                if (len_trim(res%msg) == 0) error stop "Missing failure diagnostic"
            end if
        end if
    end do
end subroutine


subroutine affine_residual (x, fx)
    !* Evaluate the weighted affine objective and record every callback.

    real (PREC), intent(in), contiguous :: x(:)
        !* Evaluation point.
    real (PREC), intent(out), contiguous :: fx(:)
        !* Weighted residual vector.

    calls = calls + 1
    if (calls > MAXFEV) error stop "DFLS exceeded the callback budget"
    fx = weights * [x(1) + x(3) - 1.0_PREC, x(2) + x(3) - 1.0_PREC, &
        x(3) - 1.0_PREC, x(3) + 1.0_PREC]
    points(:,calls) = x
    residuals(:,calls) = fx
    norms(calls) = sum(fx**2)
end subroutine

end program
