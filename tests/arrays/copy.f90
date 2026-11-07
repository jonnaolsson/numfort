

program test_numfort_arrays_copy

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common_status

    use numfort_common_testing

    use numfort_stats

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_pack_indexed ()

end subroutine

subroutine test_pack_indexed ()
    !*   Unit tests for COPY_INDEXED.

    real (PREC), dimension(:,:), allocatable :: src, dst, dst_ok
    integer, dimension(:), allocatable :: indices
    integer :: m, n, k
    logical :: values_ok
    type (status_t) :: status

    call set_seed (123)

    ! --- 2d-API: size-1 inputs ---

    m = 1
    n = 1
    allocate (src(m,n), dst(m,n))
    allocate (indices(m), source=1)

    call random_number (src)

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, status=status)
    values_ok = all (src == dst)
    if (.not. (values_ok .and. status == NF_STATUS_OK)) &
        error stop "2d-API: size-1 input array, default args"

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=1, trans=.true., status=status)
    values_ok = all (src == dst)
    if (.not. (values_ok .and. status == NF_STATUS_OK)) &
        error stop "2d-API: size-1 input array, TRANS=T"

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=2, trans=.false., status=status)
    values_ok = all (src == dst)
    if (.not. (values_ok .and. status == NF_STATUS_OK)) &
        error stop "2d-API: size-1 input array, DIM=2, TRANS=F"

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=2, trans=.true., status=status)
    values_ok = all (src == dst)
    if (.not. (values_ok .and. status == NF_STATUS_OK)) &
        error stop "2d-API: size-1 input array, DIM=2, TRANS=T"

    deallocate (indices, dst, src)

    ! --- 2d-API: 1d-inputs, arange'd indices ---

    m = 10
    n = 1
    allocate (src(m,n))
    call random_number (src)

    ! column vector SRC, trans = .false.
    allocate (dst(m,n), indices(m))
    call arange (indices)

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=1, trans=.false., status=status)
    values_ok = all (src == dst)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: column vector SRC, DIM=1, TRANS=F"

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, [1], dst, dim=2, trans=.false., status=status)
    values_ok = all (src == dst)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: column vector SRC, DIM=2, TRANS=F"

    deallocate (dst)

    ! column vector SRC, trans = .true.
    allocate (dst(n,m))
    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=1, trans=.true., status=status)
    values_ok = all (src == transpose(dst))
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: column vector SRC, DIM=1, TRANS=T"

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, [1], dst, dim=2, trans=.true., status=status)
    values_ok = all (src == transpose(dst))
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: column vector SRC, DIM=2, TRANS=T"

    deallocate (src, dst, indices)

    ! row vector SRC
    m = 1
    n = 3
    allocate (src(m,n))
    call random_number (src)

    ! row vector SRC, trans = .false.
    allocate (dst(m,n), indices(n))
    call arange (indices)

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, [1], dst, dim=1, trans=.false., status=status)
    values_ok = all (src == dst)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: row vector SRC, DIM=1, TRANS=F"

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=2, trans=.false., status=status)
    values_ok = all (src == dst)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: row vector SRC, DIM=2, TRANS=F"

    deallocate (dst)

    ! row vector SRC, trans = .true.
    allocate (dst(n,m))
    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, [1], dst, dim=1, trans=.true., status=status)
    values_ok = all (src == transpose(dst))
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: row vector SRC, DIM=1, TRANS=T"

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=2, trans=.true., status=status)
    values_ok = all (src == transpose(dst))
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: row vector SRC, DIM=2, TRANS=T"

    deallocate (src, dst, indices)

    ! --- 2d-API: 2d-inputs, arange'd indices ---

    m = 10
    n = 3
    allocate (src(m,n))
    call random_number (src)

    ! dim = 1, trans = .false.
    allocate (dst(m,n), indices(m))
    call arange (indices)

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=1, trans=.false., status=status)
    values_ok = all (src == dst)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: DIM=1, TRANS=F"

    deallocate (indices)

    ! dim = 2, trans = .false.
    allocate (indices(n))
    call arange (indices)
    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src,indices, dst, dim=2, trans=.false., status=status)
    values_ok = all (src == dst)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: DIM=2, TRANS=F"

    deallocate (dst, indices)

    ! dim = 1, trans = .true.
    allocate (dst(n,m), indices(m))
    call arange (indices)
    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=1, trans=.true., status=status)
    values_ok = all (src == transpose(dst))
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: DIM=1, TRANS=T"

    deallocate (indices)

    ! dim = 2, trans = .true.
    allocate (indices(n))
    call arange (indices)
    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=2, trans=.true., status=status)
    values_ok = all (src == transpose(dst))
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: DIM=2, TRANS=T"

    deallocate (src, dst, indices)

    ! --- 2d-API: 2d-inputs, reordered indices ---

    m = 13
    n = 2

    allocate (src(m,n))
    call random_number (src)

    ! dim = 1, trans = .false.
    allocate (indices(m))
    call random_order (indices)
    allocate (dst(m,n), dst_ok(m,n))

    dst_ok(:,:) = src(indices,:)

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=1, trans=.false., status=status)
    values_ok = all (dst == dst_ok)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: Reordered indices, DIM=1, TRANS=F"

    deallocate (indices)

    ! dim = 2, trans = .false.
    allocate (indices(n))
    call random_order (indices)
    dst_ok(:,:) = src(:,indices)

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=2, trans=.false., status=status)
    values_ok = all (dst == dst_ok)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: Reordered indices, DIM=2, TRANS=F"

    deallocate (indices, dst, dst_ok)

    ! dim = 1, trans = .true.
    allocate (indices(m), dst(n,m), dst_ok(n,m))
    call random_order (indices)

    dst_ok(:,:) = transpose(src(indices,:))

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=1, trans=.true., status=status)
    values_ok = all (dst == dst_ok)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: Reordered indices, DIM=1, TRANS=T"

    deallocate (indices)

    ! dim = 2, trans = .true.
    allocate (indices(n))
    call random_order (indices)

    dst_ok(:,:) = transpose(src(:,indices))

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=2, trans=.true., status=status)
    values_ok = all (dst == dst_ok)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: Reordered indices, DIM=2, TRANS=T"

    deallocate (src, indices, dst, dst_ok)

    ! --- 2d-API: 2d-inputs, subset of reordered indices ---

    m = 7
    n = 16
    k = 4

    allocate (src(m,n))
    call random_number (src)

    ! dim = 1, trans = .false.
    allocate (indices(k))
    call random_order (indices)
    allocate (dst(k,n), dst_ok(k,n))

    dst_ok(:,:) = src(indices(1:k),:)

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=1, trans=.false., status=status)
    values_ok = all (dst == dst_ok)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: Subset of reordered indices, DIM=1, TRANS=F"

    deallocate (indices, dst, dst_ok)

    ! dim = 2, trans = .false.
    allocate (indices(k), dst(m,k), dst_ok(m,k))
    call random_order (indices)
    dst_ok(:,:) = src(:,indices(1:k))

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=2, trans=.false., status=status)
    values_ok = all (dst == dst_ok)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: Subset of reordered indices, DIM=2, TRANS=F"

    deallocate (indices, dst, dst_ok)

    ! dim = 1, trans = .true.
    allocate (indices(k), dst(n,k), dst_ok(n,k))
    call random_order (indices)

    dst_ok(:,:) = transpose(src(indices(1:k),:))

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=1, trans=.true., status=status)
    values_ok = all (dst == dst_ok)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: Subset of reordered indices, DIM=1, TRANS=T"

    deallocate (indices, dst, dst_ok)

    ! dim = 2, trans = .true.
    allocate (indices(k), dst(k,m), dst_ok(k,m))
    call random_order (indices)

    dst_ok(:,:) = transpose(src(:,indices(1:k)))

    status = NF_STATUS_UNDEFINED
    dst(:,:) = huge(0.0_PREC)
    call pack_indexed (src, indices, dst, dim=2, trans=.true., status=status)
    values_ok = all (dst == dst_ok)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "2d-API: Subset of reordered indices, DIM=2, TRANS=T"

    deallocate (src, indices, dst, dst_ok)

end subroutine

end
