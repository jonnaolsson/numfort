

program test_numfort_linalg_core

    use, intrinsic :: iso_fortran_env

    use numfort_linalg
    use numfort_common
    use numfort_common_testing
    use numfort_stats, only: set_seed

    implicit none

    integer, parameter :: PREC = real64

    call test_gram ()
    call test_inv ()
    call test_inv_invalid ()
    call test_inv_singular ()
    call test_inv_work_query ()

    contains

subroutine test_gram ()

    real (PREC), dimension(:,:), allocatable :: a, g, a_T, g_ok
    type (status_t) :: status
    character (1) :: trans
    integer :: m, n
    logical :: values_ok

    call set_seed (1234)

    ! --- Input checks ---

    ! Invalid TRANS

    m = 10
    n = 3
    allocate (a(m,n), g(n,n))

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans="e", status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid TRANS argument"

    ! Try with valid values
    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans="t", status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "TRANS='t'"

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans="t", status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "TRANS='T'"

    deallocate (g)
    allocate (g(m,m))

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans="n", status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "TRANS='n'"

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans="N", status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "TRANS='N'"

    deallocate (a, g)

    ! Non-conformable A, G
    m = 8
    n = 5
    allocate (a(m,n), g(m+1,m+1))

    ! g should be of share (n,n)
    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans="T", status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Non-conformable G, trans=T"

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans="N", status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Non-conformable G, trans=N"

    deallocate (a, g)

    ! --- Valid invocations ---

    ! Call with trans='T'
    m = 6
    n = 11
    trans = "T"
    allocate (a(m,n), g(n,n), a_T(n,m), g_ok(n,n))

    call random_number (a)
    a_T(:,:) = transpose (a)
    g_ok(:,:) = matmul (a_T, a)

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans=trans, status=status)
    values_ok = all_close (g, g_ok, atol=1.0e-12_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) error stop "Call with trans='T'"

    deallocate (a, a_T, g, g_ok)

    ! Call with trans='N'
    m = 7
    n = 11
    trans = "N"
    allocate (a(m,n), g(m,m), a_T(n,m), g_ok(m,m))

    call random_number (a)
    a_T(:,:) = transpose (a)
    g_ok(:,:) = matmul (a, a_T)

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans=trans, status=status)
    values_ok = all_close (g, g_ok, atol=1.0e-12_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) error stop "Call with trans='N'"

    deallocate (a, a_T, g, g_ok)

    ! --- Degenerate arguments ---

    ! SHAPE(a) = [1,1]
    ! Call with trans = 'T'
    m = 1
    n = 1
    trans = "T"
    allocate (a(m,n), g(n,n), a_T(n,m), g_ok(n,n))

    call random_number (a)
    a_T(:,:) = transpose (a)
    g_ok(:,:) = matmul (a_T, a)

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans=trans, status=status)
    values_ok = all_close (g, g_ok, atol=1.0e-12_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Call with trans='T', shape(a) = [1,1]"

    deallocate (a, a_T, g, g_ok)

    ! Call with trans = 'N'
    m = 1
    n = 1
    trans = "N"
    allocate (a(m,n), g(m,m), a_T(n,m), g_ok(m,m))

    call random_number (a)
    a_T(:,:) = transpose (a)
    g_ok(:,:) = matmul (a, a_T)

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans=trans, status=status)
    values_ok = all_close (g, g_ok, atol=1.0e-12_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Call with trans='N', shape(a) = [1,1]"

    deallocate (a, a_T, g, g_ok)

    ! SHAPE(a) = [1,k]
    ! Call with trans = 'T'
    m = 1
    n = 11
    trans = "T"
    allocate (a(m,n), g(n,n), a_T(n,m), g_ok(n,n))

    call random_number (a)
    a_T(:,:) = transpose (a)
    g_ok(:,:) = matmul (a_T, a)

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans=trans, status=status)
    values_ok = all_close (g, g_ok, atol=1.0e-12_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Call with trans='T', shape(a) = [1,k]"

    deallocate (a, a_T, g, g_ok)

    ! Call with trans = 'N'
    m = 1
    n = 7
    trans = "N"
    allocate (a(m,n), g(m,m), a_T(n,m), g_ok(m,m))

    call random_number (a)
    a_T(:,:) = transpose (a)
    g_ok(:,:) = matmul (a, a_T)

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans=trans, status=status)
    values_ok = all_close (g, g_ok, atol=1.0e-12_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Call with trans='N', shape(a) = [1,k]"

    deallocate (a, a_T, g, g_ok)

    ! SHAPE(a) = [k,1]
    ! Call with trans = 'T'
    m = 8
    n = 1
    trans = "T"
    allocate (a(m,n), g(n,n), a_T(n,m), g_ok(n,n))

    call random_number (a)
    a_T(:,:) = transpose (a)
    g_ok(:,:) = matmul (a_T, a)

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans=trans, status=status)
    values_ok = all_close (g, g_ok, atol=1.0e-12_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Call with trans='T', shape(a) = [k,1]"

    deallocate (a, a_T, g, g_ok)

    ! Call with trans = 'N'
    m = 2
    n = 1
    trans = "N"
    allocate (a(m,n), g(m,m), a_T(n,m), g_ok(m,m))

    call random_number (a)
    a_T(:,:) = transpose (a)
    g_ok(:,:) = matmul (a, a_T)

    status = NF_STATUS_UNDEFINED
    call gram (a, g, trans=trans, status=status)
    values_ok = all_close (g, g_ok, atol=1.0e-12_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Call with trans='N', shape(a) = [k,1]"

    deallocate (a, a_T, g, g_ok)

end subroutine


subroutine test_inv ()
    !* - Check scalar, diagonal, and nonsymmetric inverses independently.
    !  - Preserve inputs and verify products with the inverse are identity.
    !  - Reuse explicit arrays with minimum workspace sizes.

    real (PREC) :: a1(1,1), a1_before(1,1), inverse1(1,1), expected1(1,1)
    real (PREC) :: a2(2,2), a2_before(2,2), inverse2(2,2), expected2(2,2)
    real (PREC) :: a3(3,3), a3_before(3,3), inverse3(3,3), expected3(3,3)
    real (PREC) :: rwork(3)
    integer :: iwork(3)
    type (status_t) :: status

    ! Scalar inverse: [4]**(-1) = [1/4], with internally allocated work.
    a1 = 4.0_PREC
    a1_before = a1
    expected1 = 0.25_PREC
    status = NF_STATUS_UNDEFINED
    call inv (a1, inverse1, status=status)
    call check_inverse (a1, a1_before, inverse1, expected1, status)

    ! The scalar also needs only one real and one integer workspace entry.
    status = NF_STATUS_UNDEFINED
    call inv (a1, inverse1, rwork=rwork(:1), iwork=iwork(:1), status=status)
    call check_inverse (a1, a1_before, inverse1, expected1, status)

    ! A diagonal inverse has reciprocal diagonal entries and zero off-diagonals.
    a3 = 0.0_PREC
    a3(1,1) = 2.0_PREC
    a3(2,2) = -4.0_PREC
    a3(3,3) = 0.5_PREC
    a3_before = a3
    expected3 = 0.0_PREC
    expected3(1,1) = 0.5_PREC
    expected3(2,2) = -0.25_PREC
    expected3(3,3) = 2.0_PREC
    status = NF_STATUS_UNDEFINED
    call inv (a3, inverse3, status=status)
    call check_inverse (a3, a3_before, inverse3, expected3, status)

    ! Reuse the explicit arrays for a larger matrix, again with n entries.
    status = NF_STATUS_UNDEFINED
    call inv (a3, inverse3, rwork=rwork, iwork=iwork, status=status)
    call check_inverse (a3, a3_before, inverse3, expected3, status)

    ! The nonsymmetric matrix [[1,2],[3,4]] requires row pivoting.
    a2 = reshape ([1.0_PREC, 3.0_PREC, 2.0_PREC, 4.0_PREC], [2,2])
    a2_before = a2
    ! For [[a,b],[c,d]], the inverse is [[d,-b],[-c,a]]/(a*d-b*c).
    ! Here det(A) = 1*4 - 2*3 = -2.
    expected2 = reshape ([-2.0_PREC, 1.5_PREC, 1.0_PREC, -0.5_PREC], [2,2])
    status = NF_STATUS_UNDEFINED
    call inv (a2, inverse2, status=status)
    call check_inverse (a2, a2_before, inverse2, expected2, status)

    ! Reuse the same arrays for a smaller matrix, with exactly n entries.
    status = NF_STATUS_UNDEFINED
    call inv (a2, inverse2, rwork=rwork(:2), iwork=iwork(:2), status=status)
    call check_inverse (a2, a2_before, inverse2, expected2, status)

end subroutine


subroutine test_inv_invalid ()
    !* - Reject short real or integer workspaces and incompatible shapes.
    !  - Preserve the input on every invalid invocation.

    real (PREC) :: a(2,2), a_before(2,2), ainv(2,2), rwork(2)
    real (PREC) :: rectangular(2,3), wrong_rows(3,2), wrong_columns(2,3)
    integer :: iwork(2)
    type (status_t) :: status

    a = reshape ([1.0_PREC, 3.0_PREC, 2.0_PREC, 4.0_PREC], [2,2])
    a_before = a

    ! The real workspace is one entry short, with a valid integer workspace.
    status = NF_STATUS_UNDEFINED
    call inv (a, ainv, rwork=rwork(:1), iwork=iwork, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG) .or. any (a /= a_before)) &
        error stop "inv: short real workspace"

    ! The integer workspace is one entry short, with a valid real workspace.
    status = NF_STATUS_UNDEFINED
    call inv (a, ainv, rwork=rwork, iwork=iwork(:1), status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG) .or. any (a /= a_before)) &
        error stop "inv: short integer workspace"

    ! An empty real workspace must be rejected before calling LAPACK.
    status = NF_STATUS_UNDEFINED
    call inv (a, ainv, rwork=rwork(:0), iwork=iwork, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG) .or. any (a /= a_before)) &
        error stop "inv: empty real workspace"

    ! An empty integer workspace must be rejected before calling LAPACK.
    status = NF_STATUS_UNDEFINED
    call inv (a, ainv, rwork=rwork, iwork=iwork(:0), status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG) .or. any (a /= a_before)) &
        error stop "inv: empty integer workspace"

    ! A nonsquare input is invalid even when the output shape matches.
    rectangular = 1.0_PREC
    status = NF_STATUS_UNDEFINED
    call inv (rectangular, wrong_columns, rwork=rwork, iwork=iwork, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG) &
        .or. any (rectangular /= 1.0_PREC)) &
        error stop "inv: nonsquare input"

    ! The output has the wrong number of rows for a square input.
    status = NF_STATUS_UNDEFINED
    call inv (a, wrong_rows, rwork=rwork, iwork=iwork, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG) .or. any (a /= a_before)) &
        error stop "inv: mismatched output rows"

    ! The output has the wrong number of columns for a square input.
    status = NF_STATUS_UNDEFINED
    call inv (a, wrong_columns, rwork=rwork, iwork=iwork, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG) .or. any (a /= a_before)) &
        error stop "inv: mismatched output columns"

end subroutine


subroutine test_inv_singular ()
    !* - Report exactly singular matrices as invalid state.
    !  - Preserve singular inputs with internal and explicit workspaces.

    real (PREC) :: a(2,2), a_before(2,2), ainv(2,2), rwork(2)
    integer :: iwork(2)
    type (status_t) :: status

    ! The second row is twice the first, so det(A) = 1*4 - 2*2 = 0.
    a = reshape ([1.0_PREC, 2.0_PREC, 2.0_PREC, 4.0_PREC], [2,2])
    a_before = a
    status = NF_STATUS_UNDEFINED
    call inv (a, ainv, status=status)
    if (.not. (status == NF_STATUS_INVALID_STATE) .or. any (a /= a_before)) &
        error stop "inv: singular matrix with internal workspace"

    ! The same singular input must fail with minimum explicit workspaces.
    status = NF_STATUS_UNDEFINED
    call inv (a, ainv, rwork=rwork, iwork=iwork, status=status)
    if (.not. (status == NF_STATUS_INVALID_STATE) .or. any (a /= a_before)) &
        error stop "inv: singular matrix with explicit workspace"

end subroutine


subroutine test_inv_work_query ()
    !* - Exercise both real64 workspace query overloads at several sizes.
    !  - Invert matrices using exactly the returned workspace lengths.

    integer, parameter :: SIZES(3) = [1, 3, 2]
    real (PREC), allocatable :: a(:,:), a_before(:,:), ainv(:,:), expected(:,:)
    real (PREC), allocatable :: rwork(:)
    integer, allocatable :: iwork(:)
    integer :: k, n, i, query_kind, lwork, liwork
    type (status_t) :: status

    ! Scalar, larger diagonal, then smaller diagonal fixtures.
    do k = 1, size (SIZES)
        n = SIZES(k)
        allocate (a(n,n), a_before(n,n), ainv(n,n), expected(n,n))
        a = 0.0_PREC
        expected = 0.0_PREC
        do i = 1, n
            a(i,i) = real (i + 1, PREC)
            ! The inverse of diag(2,...,n+1) has entries 1/(i+1).
            expected(i,i) = 1.0_PREC / real (i + 1, PREC)
        end do
        a_before(:,:) = a

        do query_kind = 1, 2
            status = NF_STATUS_UNDEFINED
            lwork = -1
            liwork = -1
            if (query_kind == 1) then
                ! Query by matrix, checking that the query preserves its input.
                call inv_work_query (a, lwork, liwork, status=status)
            else
                ! Query by dimension with a real64 kind selector.
                call inv_work_query (0.0_PREC, n, lwork, liwork, status=status)
            end if
            if (.not. (status == NF_STATUS_OK)) error stop "inv_work_query: status"
            if (any (a /= a_before)) error stop "inv_work_query: input changed"
            if (lwork < n .or. liwork < n) error stop "inv_work_query: short workspace"

            ! Use the returned sizes without assuming backend optimal lengths.
            allocate (rwork(lwork), iwork(liwork))
            status = NF_STATUS_UNDEFINED
            call inv (a, ainv, rwork=rwork, iwork=iwork, status=status)
            call check_inverse (a, a_before, ainv, expected, status)
            deallocate (rwork, iwork)
        end do
        deallocate (a, a_before, ainv, expected)
    end do

end subroutine


subroutine check_inverse (a, a_before, ainv, expected, status)
    !* Check a successful inversion against an independent reference.

    real (PREC), intent(in) :: a(:,:)
        !* Input matrix after inversion.
    real (PREC), intent(in) :: a_before(:,:)
        !* Copy of the input before inversion.
    real (PREC), intent(in) :: ainv(:,:)
        !* Computed inverse.
    real (PREC), intent(in) :: expected(:,:)
        !* Inverse computed from a closed-form expression.
    type (status_t), intent(in) :: status
        !* Inversion exit status.

    real (PREC), parameter :: TOL = 1.0e-12_PREC
    real (PREC) :: identity(size (a, 1),size (a, 1))
    real (PREC) :: product(size (a, 1),size (a, 1))
    integer :: i

    if (.not. (status == NF_STATUS_OK)) error stop "inv: expected success"
    if (any (a /= a_before)) error stop "inv: input changed"
    if (.not. all_close (ainv, expected, atol=TOL, rtol=0.0_PREC)) &
        error stop "inv: incorrect inverse"

    ! A*A**(-1) = I, whose diagonal entries are one and others are zero.
    identity = 0.0_PREC
    do i = 1, size (a, 1)
        identity(i,i) = 1.0_PREC
    end do
    product = matmul (a, ainv)
    if (.not. all_close (product, identity, atol=TOL, rtol=0.0_PREC)) &
        error stop "inv: product is not identity"

end subroutine

end
