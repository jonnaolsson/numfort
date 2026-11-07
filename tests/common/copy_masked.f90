

program test_numfort_common_copy_masked
    ! Unit tests for masked copy routines

    use, intrinsic :: iso_fortran_env

    use numfort_common_copy_masked
    use numfort_common_status
    use numfort_stats, only: set_seed

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_copy_2d_1dim ()

end subroutine

subroutine test_copy_2d_1dim ()

    real (PREC), dimension(:,:), allocatable :: src, dst
    logical, dimension(:), allocatable :: mask
    integer :: m, n, i, j
    logical :: values_ok
    type (status_t) :: status

    call set_seed (123)

    ! === Input checks ===

    ! 1. Non-conformable SRC, MASK along dim = 1
    m = 3
    n = 4
    allocate (src(m,n), dst(m,n), mask(m-1))
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=1, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable SRC, MASK along dim=1"
    deallocate (src, dst, mask)

    ! 2. Non-conformable SRC, MASK along dim=2
    m = 1
    n = 2
    allocate (src(m,n), dst(m,n), mask(n+1))
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=2, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable SRC, MASK along dim=2"
    deallocate (src, dst, mask)

    ! 3. Non-conformable SRC, DST along non-masked dim = 2
    m = 11
    n = 100
    allocate (src(m,n), dst(m,n-1), mask(m))
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=1, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable SRC, DST along non-masked dimension; DIM=1"
    deallocate (src, dst, mask)

    ! 4. Non-conformable SRC, DST along non-masked dim = 1
    m = 2
    n = 1
    allocate (src(m,n), dst(m+1,n), mask(n))
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=2, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Non-conformable SRC, DST along non-masked dimension; DIM=2"
    deallocate (src, dst, mask)

    ! 5. Invalid dim < 1
    m = 1
    n = 1
    allocate (src(m,n), dst(m,n), mask(m))
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=0, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid DIM < 1"
    deallocate (src, dst, mask)

    ! 6. Invalid dim > 2
    m = 12
    n = 1
    allocate (src(m,n), dst(m,n), mask(n))
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=3, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid DIM > 2"
    deallocate (src, dst, mask)

    ! === Degenerate input args ===

    ! 1. 0-sized rows, dim = 1
    m = 0
    n = 10
    allocate (src(m,n), dst(m,n), mask(m))
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=1, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "0-sized row SRC, DST; DIM=1"
    deallocate (src, dst, mask)

    ! 2. 0-sized rows, dim = 2
    m = 0
    n = 1
    allocate (src(m,n), dst(m,n), mask(n))
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=2, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "0-sized row SRC, DST; DIM=2"
    deallocate (src, dst, mask)

    ! 3. 0-sized cols, dim = 1
    m = 102
    n = 0
    allocate (src(m,n), dst(m,n), mask(m))
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=1, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "0-sized col SRC, DST; DIM=1"
    deallocate (src, dst, mask)

    ! 4. 0-sized cols, dim = 2
    m = 1
    n = 0
    allocate (src(m,n), dst(m,n), mask(n))
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=2, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "0-sized col SRC, DST; DIM=2"
    deallocate (src, dst, mask)

    ! 5. 0-sized rows and cols
    m = 0
    n = 0
    allocate (src(m,n), dst(m,n), mask(m))
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "0-sized row/col SRC, DST"
    deallocate (src, dst, mask)

    ! === No masked elements ===

    ! 1. dim = 1
    m = 3
    n = 5
    allocate (src(m,n), dst(m,n), mask(m))
    call random_number (src)
    mask(:) = .true.
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=1, status=status)
    values_ok = all(src == dst)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "No masked elements, DIM=1"
    deallocate (src, dst, mask)

    ! 2. dim = 2
    m = 1
    n = 11
    allocate (src(m,n), dst(m,n), mask(n))
    call random_number (src)
    mask(:) = .true.
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=2, status=status)
    values_ok = all(src == dst)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "No masked elements, DIM=2"
    deallocate (src, dst, mask)

    ! === All elements masked ===

    ! 1. dim = 1
    m = 5
    n = 123
    allocate (src(m,n), dst(m,n), mask(m))
    mask(:) = .false.
    src(:,:) = 1.0
    dst(:,:) = -1.0
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=1, status=status)
    ! Check that values in DST are untouched
    values_ok = all (dst == -1.0)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "All elements masked, DIM=1"
    deallocate (src, dst, mask)

    ! 2. dim = 2
    m = 3
    n = 2
    allocate (src(m,n), dst(m,n), mask(n))
    mask(:) = .false.
    src(:,:) = 1.0
    dst(:,:) = -1.0
    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=2, status=status)
    ! Check that values in DST are untouched
    values_ok = all (dst == -1.0)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "All elements masked, DIM=2"
    deallocate (src, dst, mask)

    ! === One dimension of size 1 ===

    ! 1. 1 column, dim = 1
    m = 10
    n = 1
    allocate (src(m,n), dst(m,n), mask(m))
    call random_number (src)
    dst(:,:) = -1.0
    mask(:) = .false.
    mask([1,4,8]) = .true.

    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=1, status=status)
    values_ok = all(src([1,4,8],1) == dst(1:3,1)) .and. all(dst(4:,1) == -1.0)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Row vector SRC, DST; DIM = 1"

    deallocate (src, dst, mask)

    ! (2) 1 column, dim = 2, MASK = .true.
    m = 4
    n = 1
    allocate (src(m,n), dst(m,n), mask(n))
    dst(:,:) = -1.0
    call random_number (src)
    mask(:) = .true.

    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=2, status=status)
    values_ok = all(src == dst)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Row vector SRC, DST; MASK = .TRUE., DIM=2"
    deallocate (src, dst, mask)

    ! (3) 1 column, dim = 2, MASK = .false.
    m = 11
    n = 1
    allocate (src(m,n), dst(m,n), mask(n))
    dst(:,:) = -1.0
    call random_number (src)
    mask(:) = .false.

    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=2, status=status)
    values_ok = all (dst == -1.0)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Row vector SRC, DST; MASK = .FALSE.; DIM=2"
    deallocate (src, dst, mask)

    ! (4) 1 row, DIM = 2
    m = 1
    n = 11
    allocate (src(m,n), dst(m,n), mask(n))
    call random_number (src)
    dst(:,:) = -1.0
    mask(:) = .false.
    mask([2,4,7,10,11]) = .true.

    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=2, status=status)
    values_ok = all(src(1,[2,4,7,10,11]) == dst(1,1:5)) .and. all (dst(1,6:) == -1.0)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Column vector SRC, DST; DIM = 1"

    deallocate (src, dst, mask)

    ! (5) 1 row, DIM = 1, MASK = .true.
    m = 1
    n = 5
    allocate (src(m,n), dst(m,n), mask(m))
    call random_number (src)
    dst(:,:) = -1.0
    mask(:) = .true.

    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=1, status=status)
    values_ok = all(src == dst)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Column vector SRC, DST; MASK = .TRUE.; DIM = 1"
    deallocate (src, dst, mask)

    ! (6) 1 row, DIM = 1, MASK = .false.
    m = 1
    n = 2
    allocate (src(m,n), dst(m,n), mask(m))
    call random_number (src)
    dst(:,:) = -1.0
    mask(:) = .false.

    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=1, status=status)
    values_ok = all(dst == -1.0)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Column vector SRC, DST; MASK = .FALSE.; DIM=1"
    deallocate (src, dst, mask)

    ! === "Regular" use-cases ===

    ! (1) Mask along DIM = 1, with DST of exact size
    m = 8
    n = 3
    allocate (src(m,n), dst(4,n), mask(m))
    call random_number (src)
    dst(:,:) = -1.0
    mask(:) = .false.
    mask([2,3,4,7]) = .true.

    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=1, status=status)
    values_ok = .true.
    j = 0
    do i = 1, m
        if (.not. mask(i)) cycle
        j = j + 1
        values_ok = all(dst(j,:) == src(i,:))
    end do
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Exact-size DST, DIM=1"
    deallocate (src, dst, mask)

    ! (2) Mask along DIM = 2, with DST of exact size
    m = 5
    n = 11
    allocate (src(m,n), dst(m,n-4), mask(n))
    call random_number (src)
    dst(:,:) = -1.0
    mask(:) = .false.
    mask([4,5,6,10]) = .true.

    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=2, status=status)
    values_ok = .true.
    j = 0
    do i = 1, n
        if (.not. mask(i)) cycle
        j = j + 1
        values_ok = all(dst(:,j) == src(:,i))
    end do
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Exact-size DST, DIM=2"

    deallocate (src, dst, mask)

    ! === Over-sized DST array ===

    ! DST array is permitted to be larger along masked dimension
    ! 1. dim = 1
    m = 5
    n = 3
    allocate (src(m,n), dst(m+1,n), mask(m))
    call random_number (src)
    dst(:,:) = -1.0
    mask(:) = .false.
    mask([1,3,4]) = .true.

    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=1, status=status)

    values_ok = .true.
    j = 0
    do i = 1, m
        if (.not. mask(i)) cycle
        j = j + 1
        values_ok = values_ok .and. all(dst(j,:) == src(i,:))
    end do

    values_ok = values_ok .and. all (dst(m+1:,:) == -1.0)

    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Over-sized DST array, DIM=1"

    deallocate (src, dst, mask)

    ! 2. dim = 2
    m = 7
    n = 11
    allocate (src(m,n), dst(m,n+5), mask(n))
    call random_number (src)
    dst(:,:) = -1.0
    mask(:) = .false.
    mask([4,6,7,11]) = .true.

    status = NF_STATUS_UNDEFINED
    call copy (src, dst, mask, dim=2, status=status)

    values_ok = .true.
    j = 0
    do i = 1, n
        if (.not. mask(i)) cycle
        j = j + 1
        values_ok = values_ok .and. all(dst(:,j) == src(:,i))
    end do

    values_ok = values_ok .and. all(dst(:,n+1:) == -1.0)

    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Over-sized DST array, DIM=2"

    deallocate (src, dst, mask)

end subroutine

end program
