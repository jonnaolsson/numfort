program test_arrays_grid
    !*  Unit tests for array grid constructors (linspace and powerspace).

    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_exceptions, only: ieee_get_flag, ieee_set_flag, &
        ieee_divide_by_zero, ieee_invalid
    use numfort_arrays, only: linspace, powerspace
    use numfort_common_testing, only: all_close

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

contains


subroutine test_all ()
    !*  Run all grid constructor test suites.

    call test_linspace ()
    call test_powerspace ()

end subroutine


subroutine test_linspace ()
    !*  Unit tests for LINSPACE.
    !   - Empty, singleton, and edge cases.
    !   - General multi-point, exact endpoints, and descending grids.

    call test_linspace_edges ()
    call test_linspace_general ()

end subroutine


subroutine test_linspace_edges ()
    !*  Edge cases for LINSPACE:
    !   - Empty array (size 0).
    !   - Singleton array (size 1) with IEEE flag and res_step verification.
    !   - Singleton array with optional explicit step.
    !   - Singleton array with equal endpoints.
    !   - Two-point interval.

    real (PREC), dimension(:), allocatable :: x
    real (PREC) :: res_step
    integer :: res_n
    logical :: flag_div_zero, flag_invalid

    ! Case 1: Empty array (size 0) returns n=0 and finite step 0.0.
    allocate (x(0))
    call linspace (x, 2.0_PREC, 4.0_PREC, res_n=res_n, res_step=res_step)
    if (res_n /= 0) error stop "Empty linspace: expected res_n = 0"
    if (res_step /= 0.0_PREC) error stop "Empty linspace: expected res_step = 0.0"
    deallocate (x)

    ! Case 2: Singleton array (size 1) without step argument (Issue #29).
    ! Start point xfrom = 2.0 must be preserved.
    ! Step convention must be finite (0.0), res_n = 1.
    ! IEEE_DIVIDE_BY_ZERO and IEEE_INVALID must NOT be raised.
    allocate (x(1))
    call ieee_set_flag (ieee_divide_by_zero, .false.)
    call ieee_set_flag (ieee_invalid, .false.)
    call linspace (x, 2.0_PREC, 4.0_PREC, res_n=res_n, res_step=res_step)
    call ieee_get_flag (ieee_divide_by_zero, flag_div_zero)
    call ieee_get_flag (ieee_invalid, flag_invalid)
    if (flag_div_zero) &
        error stop "Singleton linspace: IEEE divide-by-zero raised"
    if (flag_invalid) &
        error stop "Singleton linspace: IEEE invalid raised"
    if (x(1) /= 2.0_PREC) &
        error stop "Singleton linspace: x(1) /= xfrom"
    if (res_n /= 1) &
        error stop "Singleton linspace: expected res_n = 1"
    if (res_step /= 0.0_PREC) &
        error stop "Singleton linspace: expected res_step = 0.0"

    ! Case 3: Singleton array with explicit step argument.
    ! Step convention returns the supplied step, starting value xfrom = 2.0.
    call linspace (x, 2.0_PREC, 4.0_PREC, step=0.5_PREC, res_n=res_n, &
        res_step=res_step)
    if (x(1) /= 2.0_PREC) &
        error stop "Singleton linspace explicit step: x(1) /= xfrom"
    if (res_n /= 1) &
        error stop "Singleton linspace explicit step: res_n /= 1"
    if (res_step /= 0.5_PREC) &
        error stop "Singleton linspace explicit step: res_step /= step"

    ! Case 4: Singleton array with equal boundaries (xfrom == xto).
    call linspace (x, 5.0_PREC, 5.0_PREC, res_n=res_n, res_step=res_step)
    if (x(1) /= 5.0_PREC) &
        error stop "Singleton linspace equal bounds: x(1) /= 5.0"
    if (res_n /= 1) &
        error stop "Singleton linspace equal bounds: res_n /= 1"
    if (res_step /= 0.0_PREC) &
        error stop "Singleton linspace equal bounds: res_step /= 0.0"
    deallocate (x)

    ! Case 5: Two-point grid on [2.0, 4.0].
    ! Step is (4.0 - 2.0) / (2 - 1) = 2.0. Expected array is [2.0, 4.0].
    allocate (x(2))
    call linspace (x, 2.0_PREC, 4.0_PREC, res_n=res_n, res_step=res_step)
    if (x(1) /= 2.0_PREC .or. x(2) /= 4.0_PREC) &
        error stop "Two-point linspace: incorrect endpoints"
    if (res_n /= 2) error stop "Two-point linspace: res_n /= 2"
    if (res_step /= 2.0_PREC) error stop "Two-point linspace: res_step /= 2.0"
    deallocate (x)

end subroutine


subroutine test_linspace_general ()
    !*  General multi-point cases for LINSPACE:
    !   - Uniform spacing across multiple points.
    !   - Exact boundary reproduction.
    !   - Descending interval (xfrom > xto).
    !   - Degenerate interval (xfrom == xto).
    !   - Explicit step with multi-point output.

    real (PREC), dimension(:), allocatable :: x
    real (PREC) :: res_step
    integer :: res_n

    ! Case 1: Multiple points (N=5) on [0.0, 1.0].
    ! Expected closed-form values: x_i = (i-1) / 4 -> [0.0, 0.25, 0.5, 0.75, 1.0].
    allocate (x(5))
    call linspace (x, 0.0_PREC, 1.0_PREC, res_n=res_n, res_step=res_step)
    if (res_n /= 5) error stop "Multi-point linspace: res_n /= 5"
    if (.not. all_close (x, [0.0_PREC, 0.25_PREC, 0.5_PREC, 0.75_PREC, 1.0_PREC], &
            atol=1.0e-14_PREC)) &
        error stop "Multi-point linspace: incorrect values"
    if (abs (res_step - 0.25_PREC) > 1.0e-14_PREC) &
        error stop "Multi-point linspace: res_step /= 0.25"
    deallocate (x)

    ! Case 2: Exact endpoint verification on non-integer bounds [-3.7, 8.3].
    ! First element must exactly equal xfrom and last element must equal xto.
    allocate (x(11))
    call linspace (x, -3.7_PREC, 8.3_PREC)
    if (x(1) /= -3.7_PREC) error stop "Linspace exact endpoints: x(1) /= xfrom"
    if (x(11) /= 8.3_PREC) error stop "Linspace exact endpoints: x(N) /= xto"
    deallocate (x)

    ! Case 3: Descending bounds [10.0, 1.0] with N=4.
    ! Step is (1.0 - 10.0) / (4 - 1) = -3.0. Expected array is [10.0, 7.0, 4.0, 1.0].
    allocate (x(4))
    call linspace (x, 10.0_PREC, 1.0_PREC, res_n=res_n, res_step=res_step)
    if (res_n /= 4) error stop "Descending linspace: res_n /= 4"
    if (.not. all_close (x, [10.0_PREC, 7.0_PREC, 4.0_PREC, 1.0_PREC], &
            atol=1.0e-14_PREC)) &
        error stop "Descending linspace: incorrect values"
    if (abs (res_step - (-3.0_PREC)) > 1.0e-14_PREC) &
        error stop "Descending linspace: res_step /= -3.0"
    deallocate (x)

    ! Case 4: Equal bounds [5.0, 5.0] with N=4.
    ! Step is (5.0 - 5.0) / 3 = 0.0. All elements equal 5.0.
    allocate (x(4))
    call linspace (x, 5.0_PREC, 5.0_PREC, res_step=res_step)
    if (.not. all (x == 5.0_PREC)) error stop "Equal bounds linspace: x /= 5.0"
    if (res_step /= 0.0_PREC) error stop "Equal bounds linspace: res_step /= 0.0"
    deallocate (x)

    ! Case 5: Explicit step argument (step = 0.5) on [0.0, 2.0] with N=5.
    ! Generates elements by advancing: 0.0, 0.5, 1.0, 1.5, 2.0.
    allocate (x(5))
    call linspace (x, 0.0_PREC, 2.0_PREC, step=0.5_PREC, res_n=res_n, &
        res_step=res_step)
    if (res_n /= 5) error stop "Explicit step linspace: res_n /= 5"
    if (.not. all_close (x, [0.0_PREC, 0.5_PREC, 1.0_PREC, 1.5_PREC, 2.0_PREC], &
            atol=1.0e-14_PREC)) &
        error stop "Explicit step linspace: incorrect values"
    if (res_step /= 0.5_PREC) error stop "Explicit step linspace: res_step /= 0.5"
    deallocate (x)

end subroutine


subroutine test_powerspace ()
    !*  Unit tests for POWERSPACE.
    !   - Empty, singleton, and edge cases.
    !   - Exponent fixtures (pow=1, pow=2), orientations, and descending bounds.

    call test_powerspace_edges ()
    call test_powerspace_powers ()

end subroutine


subroutine test_powerspace_edges ()
    !*  Edge cases for POWERSPACE:
    !   - Empty array (size 0).
    !   - Singleton array (size 1) with sparse_end=.true. and clean IEEE flags.
    !   - Singleton array (size 1) with sparse_end=.false. and clean IEEE flags.
    !   - Singleton array with equal bounds.
    !   - Two-point interval with both orientations.
    !   - Equal bounds across multiple points.

    real (PREC), dimension(:), allocatable :: x
    logical :: flag_div_zero, flag_invalid

    ! Case 1: Empty array (size 0) returns silently without error.
    allocate (x(0))
    call powerspace (x, 2.0_PREC, 4.0_PREC, 2.0_PREC)
    deallocate (x)

    ! Case 2: Singleton array (size 1) with default sparse_end=.true. (Issue #29).
    ! Starting value xmin = 2.0 must be preserved.
    ! IEEE_DIVIDE_BY_ZERO and IEEE_INVALID must NOT be raised.
    allocate (x(1))
    call ieee_set_flag (ieee_divide_by_zero, .false.)
    call ieee_set_flag (ieee_invalid, .false.)
    call powerspace (x, 2.0_PREC, 4.0_PREC, 2.0_PREC, sparse_end=.true.)
    call ieee_get_flag (ieee_divide_by_zero, flag_div_zero)
    call ieee_get_flag (ieee_invalid, flag_invalid)
    if (flag_div_zero) &
        error stop "Singleton powerspace (sparse_end=T): IEEE divide-by-zero raised"
    if (flag_invalid) &
        error stop "Singleton powerspace (sparse_end=T): IEEE invalid raised"
    if (x(1) /= 2.0_PREC) &
        error stop "Singleton powerspace (sparse_end=T): x(1) /= xmin"

    ! Case 3: Singleton array (size 1) with sparse_end=.false. (Issue #29).
    ! Starting value xmin = 2.0 must be preserved.
    ! IEEE_DIVIDE_BY_ZERO and IEEE_INVALID must NOT be raised.
    call ieee_set_flag (ieee_divide_by_zero, .false.)
    call ieee_set_flag (ieee_invalid, .false.)
    call powerspace (x, 2.0_PREC, 4.0_PREC, 2.0_PREC, sparse_end=.false.)
    call ieee_get_flag (ieee_divide_by_zero, flag_div_zero)
    call ieee_get_flag (ieee_invalid, flag_invalid)
    if (flag_div_zero) &
        error stop "Singleton powerspace (sparse_end=F): IEEE divide-by-zero raised"
    if (flag_invalid) &
        error stop "Singleton powerspace (sparse_end=F): IEEE invalid raised"
    if (x(1) /= 2.0_PREC) &
        error stop "Singleton powerspace (sparse_end=F): x(1) /= xmin"

    ! Case 4: Singleton array with equal boundaries (xmin == xmax).
    call powerspace (x, 5.0_PREC, 5.0_PREC, 2.0_PREC)
    if (x(1) /= 5.0_PREC) &
        error stop "Singleton powerspace equal bounds: x(1) /= 5.0"
    deallocate (x)

    ! Case 5: Two points on [2.0, 4.0] with sparse_end=.true. -> [2.0, 4.0].
    allocate (x(2))
    call powerspace (x, 2.0_PREC, 4.0_PREC, 2.0_PREC, sparse_end=.true.)
    if (x(1) /= 2.0_PREC .or. x(2) /= 4.0_PREC) &
        error stop "Two-point powerspace (sparse_end=T): incorrect endpoints"

    ! Case 6: Two points on [2.0, 4.0] with sparse_end=.false. -> [2.0, 4.0].
    call powerspace (x, 2.0_PREC, 4.0_PREC, 2.0_PREC, sparse_end=.false.)
    if (x(1) /= 2.0_PREC .or. x(2) /= 4.0_PREC) &
        error stop "Two-point powerspace (sparse_end=F): incorrect endpoints"
    deallocate (x)

    ! Case 7: Equal bounds [3.0, 3.0] with N=4 points.
    allocate (x(4))
    call powerspace (x, 3.0_PREC, 3.0_PREC, 2.0_PREC)
    if (.not. all (x == 3.0_PREC)) &
        error stop "Equal bounds powerspace: x /= 3.0"
    deallocate (x)

end subroutine


subroutine test_powerspace_powers ()
    !*  Exponent fixtures, grid orientations, and bounds for POWERSPACE:
    !   - Exponent 1 (reduces to linear grid).
    !   - Hand-computable exponent 2 with sparse_end=.true.
    !   - Hand-computable exponent 2 with sparse_end=.false.
    !   - Descending interval (xmin > xmax).
    !   - Exact endpoints for large multi-point grid.

    real (PREC), dimension(:), allocatable :: x, x_lin

    ! Case 1: Exponent 1.0 on [1.0, 5.0] with N=5.
    ! u^1 = u, so powerspace must match linspace identically.
    allocate (x(5), x_lin(5))
    call powerspace (x, 1.0_PREC, 5.0_PREC, 1.0_PREC)
    call linspace (x_lin, 1.0_PREC, 5.0_PREC)
    if (.not. all_close (x, x_lin, atol=1.0e-14_PREC)) &
        error stop "Powerspace exponent 1: mismatch with linspace"
    deallocate (x, x_lin)

    ! Case 2: Hand-computable exponent 2 on [1.0, 5.0] with N=3 and sparse_end=.true.
    ! u = [0.0, 0.5, 1.0]
    ! u^2 = [0.0, 0.25, 1.0]
    ! x = xmin + (xmax - xmin) * u^2 = 1.0 + 4.0 * [0.0, 0.25, 1.0] = [1.0, 2.0, 5.0].
    ! Points are denser near xmin (first spacing = 1.0, second spacing = 3.0).
    allocate (x(3))
    call powerspace (x, 1.0_PREC, 5.0_PREC, 2.0_PREC, sparse_end=.true.)
    if (.not. all_close (x, [1.0_PREC, 2.0_PREC, 5.0_PREC], atol=1.0e-14_PREC)) &
        error stop "Hand-computable powerspace (sparse_end=T): incorrect values"

    ! Case 3: Hand-computable exponent 2 on [1.0, 5.0] with N=3 and sparse_end=.false.
    ! u = [0.0, 0.5, 1.0]
    ! u^2 = [0.0, 0.25, 1.0]
    ! 1 - u^2 = [1.0, 0.75, 0.0], reversed v = [0.0, 0.75, 1.0]
    ! x = xmin + (xmax - xmin) * v = 1.0 + 4.0 * [0.0, 0.75, 1.0] = [1.0, 4.0, 5.0].
    ! Points are denser near xmax (first spacing = 3.0, second spacing = 1.0).
    call powerspace (x, 1.0_PREC, 5.0_PREC, 2.0_PREC, sparse_end=.false.)
    if (.not. all_close (x, [1.0_PREC, 4.0_PREC, 5.0_PREC], atol=1.0e-14_PREC)) &
        error stop "Hand-computable powerspace (sparse_end=F): incorrect values"

    ! Case 4: Descending bounds [5.0, 1.0] with N=3, pow=2.0, sparse_end=.true.
    ! slope = 1.0 - 5.0 = -4.0.
    ! x = 5.0 + (-4.0) * [0.0, 0.25, 1.0] = [5.0, 4.0, 1.0].
    call powerspace (x, 5.0_PREC, 1.0_PREC, 2.0_PREC, sparse_end=.true.)
    if (.not. all_close (x, [5.0_PREC, 4.0_PREC, 1.0_PREC], atol=1.0e-14_PREC)) &
        error stop "Descending powerspace: incorrect values"
    deallocate (x)

    ! Case 5: Exact endpoints for multi-point grid (N=17) on [-2.5, 7.5].
    allocate (x(17))
    call powerspace (x, -2.5_PREC, 7.5_PREC, 2.5_PREC, sparse_end=.true.)
    if (x(1) /= -2.5_PREC) &
        error stop "Powerspace exact endpoints (sparse_end=T): x(1) /= xmin"
    if (x(17) /= 7.5_PREC) &
        error stop "Powerspace exact endpoints (sparse_end=T): x(N) /= xmax"

    call powerspace (x, -2.5_PREC, 7.5_PREC, 2.5_PREC, sparse_end=.false.)
    if (x(1) /= -2.5_PREC) &
        error stop "Powerspace exact endpoints (sparse_end=F): x(1) /= xmin"
    if (x(17) /= 7.5_PREC) &
        error stop "Powerspace exact endpoints (sparse_end=F): x(N) /= xmax"
    deallocate (x)

end subroutine


end program
