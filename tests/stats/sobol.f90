program test_numfort_stats_sobol
    !*  Unit tests for Sobol quasi-random number generator routines.

    use, intrinsic :: iso_fortran_env

    use numfort_common_testing, only: all_close
    use numfort_stats, only: sobol_init, sobol_next, sobol_reset, sobol_state

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

contains

    subroutine test_all ()
        !*  Run all Sobol test families.

        call test_sobol_initial_zero ()
        call test_sobol_reference_pts_0_7 ()
        call test_sobol_reference_pts_8_15 ()
        call test_sobol_reference_real32 ()
        call test_sobol_step ()
        call test_sobol_reinit ()
        call test_sobol_independent_streams ()
        call test_sobol_scalar_vector ()

    end subroutine


    subroutine test_sobol_initial_zero ()
        !*  Test that the initial point returned after initialization is zero.

        type (sobol_state) :: state
        real (PREC) :: x1_scalar
        real (PREC), dimension(1) :: x1_vec
        real (PREC), dimension(2) :: x2
        real (PREC), dimension(5) :: x5

        ! 1D scalar interface returns 0.0 at point 0
        call sobol_init (state, ndim=1)
        call sobol_next (state, x1_scalar)
        if (x1_scalar /= 0.0_PREC) error stop "initial zero: 1D scalar must be 0"

        ! 1D vector interface returns 0.0 at point 0
        call sobol_init (state, ndim=1)
        call sobol_next (state, x1_vec)
        if (x1_vec(1) /= 0.0_PREC) error stop "initial zero: 1D vector must be 0"

        ! 2D vector interface returns all zeros at point 0
        call sobol_init (state, ndim=2)
        call sobol_next (state, x2)
        if (any(x2 /= 0.0_PREC)) error stop "initial zero: 2D vector must be 0"

        ! 5D vector interface returns all zeros at point 0
        call sobol_init (state, ndim=5)
        call sobol_next (state, x5)
        if (any(x5 /= 0.0_PREC)) error stop "initial zero: 5D vector must be 0"

    end subroutine


    subroutine test_sobol_reference_pts_0_7 ()
        !*  Test points 0 to 7 of 5D Sobol sequence against independent Joe-Kuo
        !*  reference values (scipy.stats.qmc.Sobol).

        type (sobol_state) :: state
        real (PREC), dimension(5) :: x

        call sobol_init (state, ndim=5)

        ! Point 0: [0, 0, 0, 0, 0]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.0_PREC, 0.0_PREC, 0.0_PREC, &
            0.0_PREC, 0.0_PREC])) then
            error stop "reference 5D: point 0 mismatch"
        end if

        ! Point 1: [1/2, 1/2, 1/2, 1/2, 1/2]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.5_PREC, 0.5_PREC, 0.5_PREC, &
            0.5_PREC, 0.5_PREC])) then
            error stop "reference 5D: point 1 mismatch"
        end if

        ! Point 2: [3/4, 1/4, 1/4, 1/4, 3/4]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.75_PREC, 0.25_PREC, 0.25_PREC, &
            0.25_PREC, 0.75_PREC])) then
            error stop "reference 5D: point 2 mismatch"
        end if

        ! Point 3: [1/4, 3/4, 3/4, 3/4, 1/4]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.25_PREC, 0.75_PREC, 0.75_PREC, &
            0.75_PREC, 0.25_PREC])) then
            error stop "reference 5D: point 3 mismatch"
        end if

        ! Point 4: [3/8, 3/8, 5/8, 7/8, 3/8]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.375_PREC, 0.375_PREC, 0.625_PREC, &
            0.875_PREC, 0.375_PREC])) then
            error stop "reference 5D: point 4 mismatch"
        end if

        ! Point 5: [7/8, 7/8, 1/8, 3/8, 7/8]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.875_PREC, 0.875_PREC, 0.125_PREC, &
            0.375_PREC, 0.875_PREC])) then
            error stop "reference 5D: point 5 mismatch"
        end if

        ! Point 6: [5/8, 1/8, 7/8, 5/8, 5/8]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.625_PREC, 0.125_PREC, 0.875_PREC, &
            0.625_PREC, 0.625_PREC])) then
            error stop "reference 5D: point 6 mismatch"
        end if

        ! Point 7: [1/8, 5/8, 3/8, 1/8, 1/8]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.125_PREC, 0.625_PREC, 0.375_PREC, &
            0.125_PREC, 0.125_PREC])) then
            error stop "reference 5D: point 7 mismatch"
        end if

    end subroutine


    subroutine test_sobol_reference_pts_8_15 ()
        !*  Test points 8 to 15 of 5D Sobol sequence against independent Joe-Kuo
        !*  reference values. Specifically tests dimension 5 at point 8 where
        !*  degree-3 primitive polynomial coefficient a_5 = 2 is required.

        type (sobol_state) :: state
        real (PREC), dimension(5) :: x
        integer :: i

        call sobol_init (state, ndim=5)

        ! Discard points 0 to 7
        do i = 0, 7
            call sobol_next (state, x)
        end do

        ! Point 8: [3/16, 5/16, 15/16, 7/16, 9/16]
        ! Dimension 5 has s_5 = 3, a_5 = 2; direction number m_4 = 11, giving
        ! v_4 = 11/16. In Gray code recurrence, x_8(5) = x_7(5) ^ v_4 =
        ! (2/16) ^ (11/16) = 9/16 = 0.5625 (previously buggy code yielded 0.6875).
        call sobol_next (state, x)
        if (.not. all_close (x, [0.1875_PREC, 0.3125_PREC, 0.9375_PREC, &
            0.4375_PREC, 0.5625_PREC])) then
            error stop "reference 5D: point 8 mismatch (degree-3 coeff extraction)"
        end if

        ! Point 9: [11/16, 13/16, 7/16, 15/16, 1/16]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.6875_PREC, 0.8125_PREC, 0.4375_PREC, &
            0.9375_PREC, 0.0625_PREC])) then
            error stop "reference 5D: point 9 mismatch"
        end if

        ! Point 10: [15/16, 1/16, 11/16, 3/16, 5/16]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.9375_PREC, 0.0625_PREC, 0.6875_PREC, &
            0.1875_PREC, 0.3125_PREC])) then
            error stop "reference 5D: point 10 mismatch"
        end if

        ! Point 11: [7/16, 9/16, 3/16, 11/16, 13/16]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.4375_PREC, 0.5625_PREC, 0.1875_PREC, &
            0.6875_PREC, 0.8125_PREC])) then
            error stop "reference 5D: point 11 mismatch"
        end if

        ! Point 12: [5/16, 3/16, 5/16, 9/16, 15/16]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.3125_PREC, 0.1875_PREC, 0.3125_PREC, &
            0.5625_PREC, 0.9375_PREC])) then
            error stop "reference 5D: point 12 mismatch"
        end if

        ! Point 13: [13/16, 11/16, 13/16, 1/16, 7/16]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.8125_PREC, 0.6875_PREC, 0.8125_PREC, &
            0.0625_PREC, 0.4375_PREC])) then
            error stop "reference 5D: point 13 mismatch"
        end if

        ! Point 14: [9/16, 7/16, 1/16, 13/16, 3/16]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.5625_PREC, 0.4375_PREC, 0.0625_PREC, &
            0.8125_PREC, 0.1875_PREC])) then
            error stop "reference 5D: point 14 mismatch"
        end if

        ! Point 15: [1/16, 15/16, 9/16, 5/16, 11/16]
        call sobol_next (state, x)
        if (.not. all_close (x, [0.0625_PREC, 0.9375_PREC, 0.5625_PREC, &
            0.3125_PREC, 0.6875_PREC])) then
            error stop "reference 5D: point 15 mismatch"
        end if

    end subroutine


    subroutine test_sobol_reference_real32 ()
        !*  Test that real32 single-precision interface matches the exact dyadic
        !*  rational values for the first 16 points in 5 dimensions.

        type (sobol_state) :: state
        real (real32), dimension(5) :: x
        integer :: i

        call sobol_init (state, ndim=5)

        ! Advance past points 0 to 7
        do i = 0, 7
            call sobol_next (state, x)
        end do

        ! Point 8 in real32: [3/16, 5/16, 15/16, 7/16, 9/16]
        call sobol_next (state, x)
        if (abs(x(5) - 0.5625_real32) > 1.0e-6_real32) then
            error stop "reference real32: point 8 dim 5 mismatch"
        end if

    end subroutine


    subroutine test_sobol_step ()
        !*  Test that sequential calls progress within bounds [0, 1).

        type (sobol_state) :: state
        real (PREC), dimension(3) :: x
        integer :: i

        call sobol_init (state, ndim=3)

        do i = 0, 31
            call sobol_next (state, x)
            if (any(x < 0.0_PREC .or. x >= 1.0_PREC)) then
                error stop "step: points must lie in [0, 1)"
            end if
        end do

    end subroutine


    subroutine test_sobol_reinit ()
        !*  Test that reinitialization reproduces the identical sequence.

        type (sobol_state) :: state
        real (PREC), dimension(5, 16) :: seq1, seq2, seq3
        integer :: i

        ! Run 1: initial run
        call sobol_init (state, ndim=5)
        do i = 1, 16
            call sobol_next (state, seq1(:, i))
        end do

        ! Run 2: reinitialize directly with sobol_init
        call sobol_init (state, ndim=5)
        do i = 1, 16
            call sobol_next (state, seq2(:, i))
        end do
        if (.not. all_close (seq1, seq2)) then
            error stop "reinit: sobol_init must reset sequence"
        end if

        ! Run 3: explicit sobol_reset then sobol_init
        call sobol_reset (state)
        call sobol_init (state, ndim=5)
        do i = 1, 16
            call sobol_next (state, seq3(:, i))
        end do
        if (.not. all_close (seq1, seq3)) then
            error stop "reinit: sobol_reset then sobol_init must reset sequence"
        end if

    end subroutine


    subroutine test_sobol_independent_streams ()
        !*  Test that two sobol_state objects advance independently.

        type (sobol_state) :: s1, s2
        real (PREC), dimension(4) :: x1, x2, ref8
        integer :: i

        call sobol_init (s1, ndim=4)
        call sobol_init (s2, ndim=4)

        ! Advance s1 by 9 points (indices 0..8)
        do i = 0, 8
            call sobol_next (s1, x1)
        end do
        ref8 = x1

        ! Advance s2 by only 3 points (indices 0..2)
        do i = 0, 2
            call sobol_next (s2, x2)
        end do

        ! Advance s1 one more point (index 9)
        call sobol_next (s1, x1)

        ! Advance s2 to index 8
        do i = 3, 8
            call sobol_next (s2, x2)
        end do

        ! s2 at index 8 must match s1 at index 8
        if (.not. all_close (x2, ref8)) then
            error stop "independent streams: streams must not interfere"
        end if

    end subroutine


    subroutine test_sobol_scalar_vector ()
        !*  Test that scalar and vector interfaces yield identical 1D sequences.

        type (sobol_state) :: s_scalar, s_vector
        real (PREC) :: xs
        real (PREC), dimension(1) :: xv
        integer :: i

        call sobol_init (s_scalar, ndim=1)
        call sobol_init (s_vector, ndim=1)

        do i = 0, 15
            call sobol_next (s_scalar, xs)
            call sobol_next (s_vector, xv)
            if (xs /= xv(1)) then
                error stop "scalar_vector: scalar and vector 1D outputs must match"
            end if
        end do

    end subroutine

end program
