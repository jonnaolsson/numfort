program test_numfort_core_logexp
    !*  Unit tests for special functions in NUMFORT_CORE_LOGEXP.

    use, intrinsic :: iso_fortran_env, only: real32, real64
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: ieee_exceptions, only: ieee_get_flag, ieee_set_flag, &
        ieee_overflow, ieee_divide_by_zero, ieee_invalid

    use numfort_arrays, only: linspace
    use numfort_core, only: logaddexp, expit, logit
    use numfort_common_testing, only: all_close

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains


subroutine test_all ()
    !*  Run all LOGEXP test families.

    call test_logaddexp ()
    call test_expit ()
    call test_logit ()

end subroutine


subroutine test_logaddexp ()
    !*  Unit tests for LOGADDEXP:
    !   - Verify identity with log(exp(x) + exp(y)) on moderate grids.
    !   - Verify symmetry logaddexp(x, y) == logaddexp(y, x).

    real (PREC), dimension(:), allocatable :: xx, yy, fxx, fxx1
    integer :: n

    ! Test for non-critical values, ie not for large negative values
    n = 100
    allocate (xx(n), yy(n), fxx(n), fxx1(n))

    call linspace (xx, -20.0_PREC, 10.0_PREC)
    call linspace (yy, -10.0_PREC, 15.0_PREC)

    fxx(:) = log (exp (xx) + exp (yy))
    fxx1(:) = logaddexp (xx, yy)

    if (.not. (all_close (fxx1, fxx, atol=1.0e-10_PREC))) &
        error stop "Test LOGADDEXP(x,y) = LOG(EXP(X) + EXP(Y)) for non-critical values"

    fxx(:) = log (exp (yy) + exp (xx))
    fxx1(:) = logaddexp (yy, xx)
    if (.not. (all_close (fxx1, fxx, atol=1.0e-10_PREC))) &
        error stop "Test LOGADDEXP(y,x) = LOG(EXP(Y) + EXP(X)) for non-critical values"

end subroutine


subroutine test_expit ()
    !*  Unit tests for EXPIT:
    !   - Interior closed-form values, symmetry, SciPy references, and real32 inputs.
    !   - Large signed coordinates, subnormal tails, infinities, NaNs, and IEEE flags.

    call test_expit_interior ()
    call test_expit_extremes ()

end subroutine


subroutine test_expit_interior ()
    !*  Interior and moderate-argument tests for EXPIT:
    !   - Exact midpoint expit(0) = 1/2.
    !   - Closed-form points x = +/-log(3) yielding 3/4 and 1/4.
    !   - Reflection identity expit(x) + expit(-x) = 1.
    !   - Hard-coded SciPy reference values in real64 and real32.

    real (PREC), parameter :: X_REF(4) = [-2.0_PREC, -0.5_PREC, 0.5_PREC, 2.0_PREC]
    ! SciPy reference values: scipy.special.expit([-2.0, -0.5, 0.5, 2.0])
    real (PREC), parameter :: Y_REF(4) = [ &
        1.1920292202211755e-01_PREC, &
        3.7754066879814541e-01_PREC, &
        6.2245933120185459e-01_PREC, &
        8.8079707797788231e-01_PREC &
    ]
    real (real32), parameter :: X32_REF(3) = [-1.5_real32, 0.0_real32, 1.5_real32]
    ! SciPy reference values: scipy.special.expit(float32([-1.5, 0.0, 1.5]))
    real (real32), parameter :: Y32_REF(3) = [ &
        1.82425529e-01_real32, 5.0e-01_real32, 8.17574441e-01_real32 &
    ]

    real (PREC) :: y, ln3
    real (PREC) :: y_vec(4), sym_sum(4)
    real (real32) :: y32_vec(3)

    ! Case 1: At x = 0, 1 / (1 + exp(0)) = 1 / 2 = 0.5.
    y = expit (0.0_PREC)
    if (y /= 0.5_PREC) error stop "EXPIT(0.0) /= 0.5"

    ! Case 2: At x = log(3), 1 / (1 + 1/3) = 3 / 4 = 0.75.
    ln3 = log (3.0_PREC)
    y = expit (ln3)
    if (abs (y - 0.75_PREC) > 1.0e-15_PREC) error stop "EXPIT(log(3)) /= 0.75"

    ! Case 3: At x = -log(3), 1 / (1 + 3) = 1 / 4 = 0.25.
    y = expit (-ln3)
    if (abs (y - 0.25_PREC) > 1.0e-15_PREC) error stop "EXPIT(-log(3)) /= 0.25"

    ! Case 4: Elemental evaluation against SciPy real64 reference vector.
    y_vec = expit (X_REF)
    if (.not. all_close (y_vec, Y_REF, atol=1.0e-15_PREC, rtol=1.0e-15_PREC)) &
        error stop "EXPIT mismatch against SciPy real64 references"

    ! Case 5: Reflection symmetry expit(x) + expit(-x) = 1 across moderate and large x.
    sym_sum = expit ([-30.0_PREC, -3.5_PREC, 0.25_PREC, 50.0_PREC]) &
        + expit ([30.0_PREC, 3.5_PREC, -0.25_PREC, -50.0_PREC])
    if (.not. all_close (sym_sum, [1.0_PREC, 1.0_PREC, 1.0_PREC, 1.0_PREC], &
            atol=1.0e-15_PREC, rtol=0.0_PREC)) &
        error stop "EXPIT reflection symmetry failed"

    ! Case 6: Elemental evaluation against SciPy real32 reference vector.
    y32_vec = expit (X32_REF)
    if (.not. all_close (y32_vec, Y32_REF, atol=1.0e-6_real32, rtol=1.0e-6_real32)) &
        error stop "EXPIT mismatch against SciPy real32 references"

end subroutine


subroutine test_expit_extremes ()
    !*  Extreme-coordinate, subnormal, non-finite, and IEEE-flag tests for EXPIT:
    !   - Large negative and positive arguments without IEEE_OVERFLOW.
    !   - Subnormal tail value at x = -710 in real64 when subnormals are supported.
    !   - Exact saturation at +/-Infinity and quiet NaN propagation.
    !   - Single-precision extreme inputs (+/-120_real32) without IEEE_OVERFLOW.

    real (PREC) :: y, inf, ninf, qnan
    real (real32) :: y32(2)
    logical :: flag_ovf, flag_div0, flag_inv, orig_gradual, restore_uflow

    inf = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    ninf = ieee_value (0.0_PREC, IEEE_NEGATIVE_INF)
    qnan = ieee_value (0.0_PREC, IEEE_QUIET_NAN)

    ! Case 1: x = -800 in real64 underflows smoothly to 0.0 without IEEE_OVERFLOW
    ! (naive exp(-(-800)) = exp(800) overflows binary64 huge ~ exp(709.78)).
    call clear_fp_flags ()
    y = expit (-800.0_PREC)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (y /= 0.0_PREC) error stop "EXPIT(-800) /= 0.0"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "EXPIT(-800) raised unexpected IEEE exception flag"

    ! Case 2: x = +800 in real64 saturates to 1.0 without IEEE exception flags.
    call clear_fp_flags ()
    y = expit (800.0_PREC)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (y /= 1.0_PREC) error stop "EXPIT(800) /= 1.0"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "EXPIT(800) raised unexpected IEEE exception flag"

    ! Case 3: Subnormal tail x = -710: 1 / (1 + exp(710)) ~ exp(-710) =
    ! 4.4762862256751298e-309 (subnormal > 0, no IEEE_OVERFLOW).
    if (ieee_support_subnormal (0.0_PREC)) then
        restore_uflow = ieee_support_underflow_control (0.0_PREC)
        if (restore_uflow) then
            call ieee_get_underflow_mode (orig_gradual)
            call ieee_set_underflow_mode (.true.)
        end if

        call clear_fp_flags ()
        y = expit (-710.0_PREC)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)

        if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)

        if (abs (y - 4.4762862256751298e-309_PREC) > 1.0e-321_PREC) &
            error stop "EXPIT(-710) subnormal tail mismatch"
        if (flag_ovf .or. flag_div0 .or. flag_inv) &
            error stop "EXPIT(-710) raised unexpected IEEE exception flag"
    end if

    ! Case 4: x = -Inf maps to 0.0 without IEEE exception flags.
    call clear_fp_flags ()
    y = expit (ninf)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (y /= 0.0_PREC) error stop "EXPIT(-Inf) /= 0.0"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "EXPIT(-Inf) raised unexpected IEEE exception flag"

    ! Case 5: x = +Inf maps to 1.0 without IEEE exception flags.
    call clear_fp_flags ()
    y = expit (inf)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (y /= 1.0_PREC) error stop "EXPIT(+Inf) /= 1.0"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "EXPIT(+Inf) raised unexpected IEEE exception flag"

    ! Case 6: Quiet NaN propagates to NaN without IEEE_INVALID.
    call clear_fp_flags ()
    y = expit (qnan)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (.not. ieee_is_nan (y)) error stop "EXPIT(NaN) did not return NaN"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "EXPIT(NaN) raised unexpected IEEE exception flag"

    ! Case 7: Single-precision extremes [-120, +120] (exp(120) overflows real32
    ! huge ~ exp(88.72), while exp(-120) underflows below 2^-149 ~ exp(-103.28)).
    call clear_fp_flags ()
    y32(1) = expit (-120.0_real32)
    y32(2) = expit (120.0_real32)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (y32(1) /= 0.0_real32 .or. y32(2) /= 1.0_real32) &
        error stop "EXPIT real32 extremes mismatch"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "EXPIT real32 extremes raised unexpected IEEE exception flag"

end subroutine


subroutine test_logit ()
    !*  Unit tests for LOGIT:
    !   - Interior closed-form, near-midpoint, near-0/near-1, and subnormal fixtures.
    !   - Exact endpoints (0, 1) and out-of-domain / NaN inputs.

    call test_logit_interior ()
    call test_logit_boundaries ()

end subroutine


subroutine test_logit_interior ()
    !*  Interior and tail accuracy tests for LOGIT:
    !   - Exact midpoint logit(0.5) = 0 and quartiles 0.25, 0.75 -> -/+log(3).
    !   - Cancellation-sensitive offset x = 0.5 + 1e-12 vs SciPy reference.
    !   - Tiny normal x = 1e-300 and subnormal x = 1e-310 without IEEE_OVERFLOW.
    !   - Near-unity x = 1 - 2^(-52) and single-precision subnormal x = 1e-40_real32.

    real (PREC) :: x_in, y, ln3
    real (real32) :: x32_in, y32
    logical :: flag_ovf, flag_div0, flag_inv, orig_gradual, restore_uflow

    ln3 = log (3.0_PREC)

    ! Case 1: At x = 0.5, log(0.5 / 0.5) = 0.0.
    y = logit (0.5_PREC)
    if (y /= 0.0_PREC) error stop "LOGIT(0.5) /= 0.0"

    ! Case 2: At x = 0.25, log(0.25 / 0.75) = -log(3).
    y = logit (0.25_PREC)
    if (abs (y + ln3) > 1.0e-15_PREC) error stop "LOGIT(0.25) /= -log(3)"

    ! Case 3: At x = 0.75, log(0.75 / 0.25) = +log(3).
    y = logit (0.75_PREC)
    if (abs (y - ln3) > 1.0e-15_PREC) error stop "LOGIT(0.75) /= +log(3)"

    ! Case 4: Near-midpoint x = fl(0.5 + 1e-12) = 0.50000000000099998;
    ! SciPy reference scipy.special.logit(0.5 + 1e-12) = 3.9999115131195140e-12.
    call clear_fp_flags ()
    y = logit (0.5_PREC + 1.0e-12_PREC)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (abs (y - 3.9999115131195140e-12_PREC) > 1.0e-27_PREC) &
        error stop "LOGIT(0.5 + 1e-12) midpoint cancellation mismatch"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "LOGIT(0.5 + 1e-12) raised unexpected IEEE exception flag"

    ! Case 5: Tiny normal x = 1e-300 -> log(1e-300) = -690.77552789821368.
    call clear_fp_flags ()
    y = logit (1.0e-300_PREC)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (abs (y - (-690.77552789821368_PREC)) > 1.0e-12_PREC) &
        error stop "LOGIT(1e-300) mismatch"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "LOGIT(1e-300) raised unexpected IEEE exception flag"

    ! Case 6: Subnormal x = 1e-310 where 1/x overflows binary64;
    ! SciPy reference scipy.special.logit(1e-310) = -713.80137882815416.
    if (ieee_support_subnormal (0.0_PREC)) then
        restore_uflow = ieee_support_underflow_control (0.0_PREC)
        if (restore_uflow) then
            call ieee_get_underflow_mode (orig_gradual)
            call ieee_set_underflow_mode (.true.)
        end if

        ! Construct 1e-310 at runtime after enabling gradual underflow.
        x_in = 1.0e-300_PREC
        x_in = x_in * 1.0e-10_PREC

        call clear_fp_flags ()
        y = logit (x_in)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)

        if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)

        if (abs (y - (-713.80137882815416_PREC)) > 1.0e-12_PREC) &
            error stop "LOGIT(1e-310) subnormal mismatch"
        if (flag_ovf .or. flag_div0 .or. flag_inv) &
            error stop "LOGIT(1e-310) raised unexpected IEEE exception flag"
    end if

    ! Case 7: Near-one x = 1 - 2^(-52): log((1 - 2^(-52)) / 2^(-52)) =
    ! 36.043653389117154 (SciPy reference).
    call clear_fp_flags ()
    y = logit (1.0_PREC - scale (1.0_PREC, -52))
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (abs (y - 36.043653389117154_PREC) > 1.0e-13_PREC) &
        error stop "LOGIT(1 - 2^-52) mismatch"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "LOGIT(1 - 2^-52) raised unexpected IEEE exception flag"

    ! Case 8: Single-precision subnormal x = 1e-40_real32 (1/x overflows real32);
    ! SciPy reference scipy.special.logit(float32(1e-40)) = -92.1034088.
    if (ieee_support_subnormal (0.0_real32)) then
        restore_uflow = ieee_support_underflow_control (0.0_real32)
        if (restore_uflow) then
            call ieee_get_underflow_mode (orig_gradual)
            call ieee_set_underflow_mode (.true.)
        end if

        x32_in = 1.0e-30_real32
        x32_in = x32_in * 1.0e-10_real32

        call clear_fp_flags ()
        y32 = logit (x32_in)
        call get_fp_flags (flag_ovf, flag_div0, flag_inv)

        if (restore_uflow) call ieee_set_underflow_mode (orig_gradual)

        if (abs (y32 - (-92.1034088_real32)) > 1.0e-4_real32) &
            error stop "LOGIT(1e-40_real32) mismatch"
        if (flag_ovf .or. flag_div0 .or. flag_inv) &
            error stop "LOGIT(1e-40_real32) raised unexpected IEEE exception flag"
    end if

end subroutine


subroutine test_logit_boundaries ()
    !*  Endpoint and domain-error tests for LOGIT:
    !   - Exact endpoints x = 0 -> -Inf and x = 1 -> +Inf.
    !   - Out-of-domain inputs x < 0, x > 1, and quiet NaN -> quiet NaN.
    !   - Single-precision vector endpoint and out-of-domain checks.

    real (PREC) :: y, qnan
    real (real32) :: y32(4)
    logical :: flag_ovf, flag_div0, flag_inv

    qnan = ieee_value (0.0_PREC, IEEE_QUIET_NAN)

    ! Case 1: Exact lower endpoint x = 0.0 -> -Inf without divide-by-zero flag.
    call clear_fp_flags ()
    y = logit (0.0_PREC)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (.not. (ieee_class (y) == IEEE_NEGATIVE_INF)) &
        error stop "LOGIT(0.0) /= -Inf"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "LOGIT(0.0) raised unexpected IEEE exception flag"

    ! Case 2: Exact upper endpoint x = 1.0 -> +Inf without divide-by-zero flag.
    call clear_fp_flags ()
    y = logit (1.0_PREC)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (.not. (ieee_class (y) == IEEE_POSITIVE_INF)) &
        error stop "LOGIT(1.0) /= +Inf"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "LOGIT(1.0) raised unexpected IEEE exception flag"

    ! Case 3: Out-of-domain x = -0.1 < 0 -> quiet NaN without IEEE_INVALID.
    call clear_fp_flags ()
    y = logit (-0.1_PREC)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (.not. ieee_is_nan (y)) error stop "LOGIT(-0.1) did not return NaN"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "LOGIT(-0.1) raised unexpected IEEE exception flag"

    ! Case 4: Out-of-domain x = 1.1 > 1 -> quiet NaN without IEEE_INVALID.
    call clear_fp_flags ()
    y = logit (1.1_PREC)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (.not. ieee_is_nan (y)) error stop "LOGIT(1.1) did not return NaN"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "LOGIT(1.1) raised unexpected IEEE exception flag"

    ! Case 5: Quiet NaN input -> quiet NaN without IEEE_INVALID.
    call clear_fp_flags ()
    y = logit (qnan)
    call get_fp_flags (flag_ovf, flag_div0, flag_inv)
    if (.not. ieee_is_nan (y)) error stop "LOGIT(NaN) did not return NaN"
    if (flag_ovf .or. flag_div0 .or. flag_inv) &
        error stop "LOGIT(NaN) raised unexpected IEEE exception flag"

    ! Case 6: Single-precision elemental boundary and out-of-domain vector.
    y32 = logit ([0.0_real32, 1.0_real32, -0.5_real32, 1.5_real32])
    if (.not. (ieee_class (y32(1)) == IEEE_NEGATIVE_INF)) &
        error stop "LOGIT(0.0_real32) /= -Inf"
    if (.not. (ieee_class (y32(2)) == IEEE_POSITIVE_INF)) &
        error stop "LOGIT(1.0_real32) /= +Inf"
    if (.not. (ieee_is_nan (y32(3)) .and. ieee_is_nan (y32(4)))) &
        error stop "LOGIT real32 out-of-domain did not return NaN"

end subroutine


subroutine clear_fp_flags ()
    !*  Clear IEEE overflow, divide-by-zero, and invalid exception flags.

    call ieee_set_flag (ieee_overflow, .false.)
    call ieee_set_flag (ieee_divide_by_zero, .false.)
    call ieee_set_flag (ieee_invalid, .false.)

end subroutine


subroutine get_fp_flags (flag_ovf, flag_div0, flag_inv)
    !*  Query IEEE overflow, divide-by-zero, and invalid exception flags.
    logical, intent(out) :: flag_ovf
        !*  True if IEEE_OVERFLOW is signaling.
    logical, intent(out) :: flag_div0
        !*  True if IEEE_DIVIDE_BY_ZERO is signaling.
    logical, intent(out) :: flag_inv
        !*  True if IEEE_INVALID is signaling.

    call ieee_get_flag (ieee_overflow, flag_ovf)
    call ieee_get_flag (ieee_divide_by_zero, flag_div0)
    call ieee_get_flag (ieee_invalid, flag_inv)

end subroutine


end program
