

program test_numfort_linalg_core

    use, intrinsic :: iso_fortran_env

    use numfort_linalg
    use numfort_common
    use numfort_common_testing
    use numfort_stats, only: set_seed

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_gram ()

end subroutine

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

end
