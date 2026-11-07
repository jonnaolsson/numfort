

program test_common_cond_alloc
    !*  Module contains unit test routines for the array INSERT routine.

    use, intrinsic :: iso_fortran_env
    use numfort_common, only: shape_equal, cond_alloc

    use numfort_stats, only: set_seed

    implicit none

    integer, parameter :: PREC = real32

    call test_all ()

    contains

subroutine test_all ()

    call test_1d ()
    call test_1d_array ()
    call test_2d ()
    call test_3d ()
    call test_4d ()
    call test_5d ()
    call test_6d ()
    call test_7d ()

end subroutine

subroutine test_1d ()
    !*  Unit tests for in-place insertion into given argument array.

    integer, parameter :: NDIM = 1
    integer, parameter :: n = 2
    integer, dimension(NDIM) :: shp
    integer, dimension(:), allocatable :: arr
    integer :: stat, source
    integer :: j

    ! Test with scalar interface
    do j = 0, n

        shp = j

        stat = 1000
        call cond_alloc (arr, shp(1), stat=stat)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an unallocated array"

        ! Test with allocated array of correct shape
        stat = 1000
        call cond_alloc (arr, shp(1), stat=stat)
        if (.not. (all (shape (arr) == shp) .and. stat == - 1)) &
            error stop "COND_ALLOC check failed for a correctly shaped array"

        ! Test with allocate array of wrong shape
        shp = shp + 1
        stat = 1000
        call cond_alloc (arr, shp(1), stat=stat)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an array with the wrong shape"

        ! Test with source argument
        deallocate (arr)
        source = 10
        stat = 1000
        call cond_alloc (arr, shp(1), source=source, stat=stat)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == &
            source))) &
            error stop "COND_ALLOC check failed for an unallocated array and scalar source"

        ! Test with source argument, allocated ARR
        stat = 1000
        call cond_alloc (arr, shp(1), source=source, stat=stat)
        if (.not. (all (shape (arr) == shp) .and. stat == - 1)) &
            error stop "COND_ALLOC check failed for a correctly shaped array and " // &
                "scalar source"

        ! Test with source argument, allocated ARR of wrong size
        shp = shp + 1
        stat = 1000
        source = 123
        call cond_alloc (arr, shp(1), source=source, stat=stat)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == &
            source))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and scalar source"
        deallocate (arr)
    end do

end subroutine

subroutine test_1d_array ()

    integer, parameter :: NDIM = 1
    integer, dimension(NDIM) :: shp
    integer, dimension(:), allocatable :: shp2
    real (PREC), dimension(NDIM) :: rwork
    real (PREC), dimension(:), allocatable :: arr, src, xx
    real (PREC) :: source
    integer :: stat, i

    call set_seed (1234)

    do i = 1, 50

        call random_number (rwork)
        shp = int(ceiling(rwork * 10.0))

        ! 1. Unallocated array
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an unallocated array"
        arr(:) = 0.0

        ! 2. Allocated array of correct size
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == - 1)) &
            error stop "COND_ALLOC check failed for a correctly shaped array"
        arr(:) = 0.0

        ! 3. Allocated array of wrong shape
        shp = shp + 1
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an array with the wrong shape"
        arr(:) = 0.0

        ! 4. Unallocated array, scalar source argument
        deallocate (arr)
        source = 123.4
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == &
            source))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "scalar source"
        arr(:) = 0.0

        ! 5. Allocated array of wrong shape, scalar source argument
        shp = shp + 1
        source = 234.5
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == &
            source))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and scalar source"
        arr(:) = 0.0

        ! 6. Unallocated array, array source argument
        deallocate (arr)
        allocate (src(shp(1)), source=0.12345_PREC)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == &
            src))) &
            error stop "COND_ALLOC check failed for an unallocated array and array source"
        arr(:) = 0.0

        ! 7. Allocated array of wrong shape, array source argument
        shp = shp + 10
        call cond_alloc (arr, shp)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == &
            src))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and array source"

        ! 8. Allocated array of correct size, array source argument
        arr(:) = 0.0
        src(:) = 123.45
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == - 1 .and. all (arr &
            == 0.0))) &
            error stop "COND_ALLOC check failed for a correctly shaped array and " // &
                "array source"

        deallocate (arr, src)
    end do

    ! Test with invalid shape argument
    shp = [10]
    do i = 0, 5
        if (i == NDIM) cycle
        allocate (shp2(i), source=1)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. ( .not. allocated (arr) .and. stat == - 2)) &
            error stop "Unallocated ARR, invalid SHP=[10]"

        allocate (arr(1), source=0.0_PREC)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == [1]) .and. stat == - 2)) &
            error stop "Allocated ARR, invalid SHP=[10]"

        deallocate (arr, shp2)
    end do

end subroutine

subroutine test_2d ()

    integer, parameter :: NDIM = 2
    integer, dimension(NDIM) :: shp
    integer, dimension(:), allocatable :: shp2
    real (PREC), dimension(NDIM) :: rwork
    real (PREC), dimension(:,:), allocatable :: arr, src, xx
    real (PREC) :: source
    integer :: stat, i

    call set_seed (1234)

    do i = 1, 50

        call random_number (rwork)
        shp = int(ceiling(rwork * 10.0))

        ! 1. Unallocated array
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an unallocated array"
        arr(:,:) = 0.0

        ! 2. Allocated array of correct size
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == - 1)) &
            error stop "COND_ALLOC check failed for a correctly shaped array"
        arr(:,:) = 0.0

        ! 3. Allocated array of wrong shape
        shp = shp + 1
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an array with the wrong shape"
        arr(:,:) = 0.0

        ! 4. Unallocated array, scalar source argument
        deallocate (arr)
        source = 123.4
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == &
            source))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "scalar source"
        arr(:,:) = 0.0

        ! 5. Allocated array of wrong shape, scalar source argument
        shp = shp + 1
        source = 234.5
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == &
            source))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and scalar source"
        arr(:,:) = 0.0

        ! 6. Unallocated array, array source argument
        deallocate (arr)
        allocate (src(shp(1), shp(2)), source=0.12345_PREC)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == &
            src))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "array source"
        arr(:,:) = 0.0

        ! 7. Allocated array of wrong shape, array source argument
        shp = shp + 10
        call cond_alloc (arr, shp)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == &
            src))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and array source"

        ! 8. Allocated array of correct size, array source argument
        arr(:,:) = 0.0
        src(:,:) = 123.45
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == - 1 .and. all (arr &
            == 0.0))) &
            error stop "COND_ALLOC check failed for a correctly shaped array and array source"

        deallocate (arr, src)
    end do

    ! Test with invalid shape argument
    shp = [10,2]
    do i = 0, 5
        if (i == NDIM) cycle
        allocate (shp2(i), source=1)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. ( .not. allocated (arr) .and. stat == - 2)) &
            error stop "Unallocated ARR, invalid SHP=[10, 2]"

        allocate (arr(1,2), source=0.0_PREC)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == [1,2]) .and. stat == - 2)) &
            error stop "Allocated ARR, invalid SHP=[10, 2]"

        deallocate (arr, shp2)
    end do

end subroutine

subroutine test_3d ()

    integer, parameter :: NDIM = 3
    integer, dimension(NDIM) :: shp
    integer, dimension(:), allocatable :: shp2
    real (PREC), dimension(NDIM) :: rwork
    real (PREC), dimension(:,:,:), allocatable :: arr, src, xx
    real (PREC) :: source
    integer :: stat, i

    call set_seed (1234)

    do i = 1, 50

        call random_number (rwork)
        shp = int(ceiling(rwork * 10.0))

        ! 1. Unallocated array
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an unallocated array"
        arr(:,:,:) = 0.0

        ! 2. Allocated array of correct size
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == - 1)) &
            error stop "COND_ALLOC check failed for a correctly shaped array"
        arr(:,:,:) = 0.0

        ! 3. Allocated array of wrong shape
        shp = shp + 1
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an array with the wrong shape"
        arr(:,:,:) = 0.0

        ! 4. Unallocated array, scalar source argument
        deallocate (arr)
        source = 123.4
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == &
            source))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "scalar source"
        arr(:,:,:) = 0.0

        ! 5. Allocated array of wrong shape, scalar source argument
        shp = shp + 1
        source = 234.5
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == source))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and scalar source"
        arr(:,:,:) = 0.0

        ! 6. Unallocated array, array source argument
        deallocate (arr)
        allocate (src(shp(1), shp(2), shp(3)), source=0.12345_PREC)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == src))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "array source"
        arr(:,:,:) = 0.0

        ! 7. Allocated array of wrong shape, array source argument
        shp = shp + 10
        call cond_alloc (arr, shp)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == src))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and array source"

        ! 8. Allocated array of correct size, array source argument
        arr(:,:,:) = 0.0
        src(:,:,:) = 123.45
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == - 1 .and. all (arr == 0.0))) &
            error stop "COND_ALLOC check failed for a correctly shaped array and " // &
                "array source"

        deallocate (arr, src)
    end do

    ! Test with invalid shape argument
    shp = [10,2,3]
    do i = 0, 5
        if (i == NDIM) cycle
        allocate (shp2(i), source=1)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. ( .not. allocated (arr) .and. stat == - 2)) &
            error stop "Unallocated ARR, invalid SHP=[10, 2, 3]"

        allocate (arr(1,2,3), source=0.0_PREC)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == [1,2,3]) .and. stat == - 2)) &
            error stop "Allocated ARR, invalid SHP=[10, 2, 3]"

        deallocate (arr, shp2)
    end do

end subroutine

subroutine test_4d ()

    integer, parameter :: NDIM = 4
    integer, dimension(NDIM) :: shp
    integer, dimension(:), allocatable :: shp2
    real (PREC), dimension(NDIM) :: rwork
    real (PREC), dimension(:,:,:,:), allocatable :: arr, src, xx
    real (PREC) :: source
    integer :: stat, i

    call set_seed (1234)

    do i = 1, 50

        call random_number (rwork)
        shp = int(ceiling(rwork * 5.0))

        ! 1. Unallocated array
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an unallocated array"
        arr(:,:,:,:) = 0.0

        ! 2. Allocated array of correct size
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == - 1)) &
            error stop "COND_ALLOC check failed for a correctly shaped array"
        arr(:,:,:,:) = 0.0

        ! 3. Allocated array of wrong shape
        shp = shp + 1
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an array with the wrong shape"
        arr(:,:,:,:) = 0.0

        ! 4. Unallocated array, scalar source argument
        deallocate (arr)
        source = 123.4
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == &
            source))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "scalar source"
        arr(:,:,:,:) = 0.0

        ! 5. Allocated array of wrong shape, scalar source argument
        shp = shp + 1
        source = 234.5
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == &
            source))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and scalar source"
        arr(:,:,:,:) = 0.0

        ! 6. Unallocated array, array source argument
        deallocate (arr)
        allocate (src(shp(1),shp(2),shp(3),shp(4)), source=0.12345_PREC)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == src))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "array source"
        arr(:,:,:,:) = 0.0

        ! 7. Allocated array of wrong shape, array source argument
        shp = shp + 10
        call cond_alloc (arr, shp)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == src))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and array source"

        ! 8. Allocated array of correct size, array source argument
        arr(:,:,:,:) = 0.0
        src(:,:,:,:) = 123.45
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == - 1 .and. all (arr == 0.0))) &
            error stop "COND_ALLOC check failed for a correctly shaped array and " // &
                "array source"

        deallocate (arr, src)
    end do

    ! Test with invalid shape argument
    shp = [10,2,3,7]
    do i = 0, 5
        if (i == NDIM) cycle
        allocate (shp2(i), source=1)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. ( .not. allocated (arr) .and. stat == - 2)) &
            error stop "Unallocated ARR, invalid SHP=[10, 2, 3, 7]"

        allocate (arr(1,2,3,7), source=0.0_PREC)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == [1,2,3,7]) .and. stat == - 2)) &
            error stop "Allocated ARR, invalid SHP=[10, 2, 3, 7]"

        deallocate (arr, shp2)
    end do

end subroutine

subroutine test_5d ()

    integer, parameter :: NDIM = 5
    integer, dimension(NDIM) :: shp
    integer, dimension(:), allocatable :: shp2
    real (PREC), dimension(NDIM) :: rwork
    real (PREC), dimension(:,:,:,:,:), allocatable :: arr, src, xx
    real (PREC) :: source
    integer :: stat, i

    call set_seed (1234)

    do i = 1, 50

        call random_number (rwork)
        shp = int(ceiling(rwork * 3.0))

        ! 1. Unallocated array
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an unallocated array"
        arr(:,:,:,:,:) = 0.0

        ! 2. Allocated array of correct size
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == - 1)) &
            error stop "COND_ALLOC check failed for a correctly shaped array"
        arr(:,:,:,:,:) = 0.0

        ! 3. Allocated array of wrong shape
        shp = shp + 1
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an array with the wrong shape"
        arr(:,:,:,:,:) = 0.0

        ! 4. Unallocated array, scalar source argument
        deallocate (arr)
        source = 123.4
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == source))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "scalar source"
        arr(:,:,:,:,:) = 0.0

        ! 5. Allocated array of wrong shape, scalar source argument
        shp = shp + 1
        source = 234.5
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == &
            source))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and scalar source"
        arr(:,:,:,:,:) = 0.0

        ! 6. Unallocated array, array source argument
        deallocate (arr)
        allocate (src(shp(1),shp(2),shp(3),shp(4),shp(5)), source=0.12345_PREC)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == src))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "array source"
        arr(:,:,:,:,:) = 0.0

        ! 7. Allocated array of wrong shape, array source argument
        shp = shp + 10
        call cond_alloc (arr, shp)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == src))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and array source"

        ! 8. Allocated array of correct size, array source argument
        arr(:,:,:,:,:) = 0.0
        src(:,:,:,:,:) = 123.45
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == - 1 .and. all (arr &
            == 0.0))) &
            error stop "COND_ALLOC check failed for a correctly shaped array and " // &
                "array source"

        deallocate (arr, src)
    end do

    ! Test with invalid shape argument
    shp = [10,2,3,7,11]
    do i = 0, 7
        if (i == NDIM) cycle
        allocate (shp2(i), source=1)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. ( .not. allocated (arr) .and. stat == - 2)) &
            error stop "Unallocated ARR, invalid SHP=[10, 2, 3, 7, 11]"

        allocate (arr(1,2,3,7,12), source=0.0_PREC)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == [1,2,3,7,12]) .and. stat == - 2)) &
            error stop "Allocated ARR, invalid SHP=[10, 2, 3, 7, 11]"

        deallocate (arr, shp2)
    end do

end subroutine

subroutine test_6d ()

    integer, parameter :: NDIM = 6
    integer, dimension(NDIM) :: shp
    integer, dimension(:), allocatable :: shp2
    real (PREC), dimension(NDIM) :: rwork
    real (PREC), dimension(:,:,:,:,:,:), allocatable :: arr, src, xx
    real (PREC) :: source
    integer :: stat, i

    call set_seed (1234)

    do i = 1, 50

        call random_number (rwork)
        shp = int(ceiling(rwork * 3.0))

        ! 1. Unallocated array
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an unallocated array"
        arr(:,:,:,:,:,:) = 0.0

        ! 2. Allocated array of correct size
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == - 1)) &
            error stop "COND_ALLOC check failed for a correctly shaped array"
        arr(:,:,:,:,:,:) = 0.0

        ! 3. Allocated array of wrong shape
        shp = shp + 1
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an array with the wrong shape"
        arr(:,:,:,:,:,:) = 0.0

        ! 4. Unallocated array, scalar source argument
        deallocate (arr)
        source = 123.4
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == source))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "scalar source"
        arr(:,:,:,:,:,:) = 0.0

        ! 5. Allocated array of wrong shape, scalar source argument
        shp = shp + 1
        source = 234.5
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == source))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and scalar source"
        arr(:,:,:,:,:,:) = 0.0

        ! 6. Unallocated array, array source argument
        deallocate (arr)
        allocate (src(shp(1),shp(2),shp(3),shp(4),shp(5),shp(6)), source=0.12345_PREC)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == src))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "array source"
        arr(:,:,:,:,:,:) = 0.0

        ! 7. Allocated array of wrong shape, array source argument
        shp = shp + 10
        call cond_alloc (arr, shp)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == src))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and array source"

        ! 8. Allocated array of correct size, array source argument
        arr(:,:,:,:,:,:) = 0.0
        src(:,:,:,:,:,:) = 123.45
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == - 1 .and. all (arr == 0.0))) &
            error stop "COND_ALLOC check failed for a correctly shaped array and " // &
                "array source"

        deallocate (arr, src)
    end do

    ! Test with invalid shape argument
    shp = [10,2,3,7,1,3]
    do i = 0, 7
        if (i == NDIM) cycle
        allocate (shp2(i), source=1)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. ( .not. allocated (arr) .and. stat == - 2)) &
            error stop "Unallocated ARR, invalid SHP=[10, 2, 3, 7, 1, 3]"

        allocate (arr(1,2,3,7,2,1), source=0.0_PREC)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == [1,2,3,7,2,1]) .and. stat == - 2)) &
            error stop "Allocated ARR, invalid SHP=[10, 2, 3, 7, 1, 3]"

        deallocate (arr, shp2)
    end do

end subroutine

subroutine test_7d ()

    integer, parameter :: NDIM = 7
    integer, dimension(NDIM) :: shp
    integer, dimension(:), allocatable :: shp2
    real (PREC), dimension(NDIM) :: rwork
    real (PREC), dimension(:,:,:,:,:,:,:), allocatable :: arr, src, xx
    real (PREC) :: source
    integer :: stat, i

    call set_seed (1234)

    do i = 1, 50

        call random_number (rwork)
        shp = int(ceiling(rwork * 3.0))

        ! 1. Unallocated array
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an unallocated array"
        arr(:,:,:,:,:,:,:) = 0.0

        ! 2. Allocated array of correct size
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == - 1)) &
            error stop "COND_ALLOC check failed for a correctly shaped array"
        arr(:,:,:,:,:,:,:) = 0.0

        ! 3. Allocated array of wrong shape
        shp = shp + 1
        stat = 1000
        call cond_alloc (arr, shp, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0)) &
            error stop "COND_ALLOC check failed for an array with the wrong shape"
        arr(:,:,:,:,:,:,:) = 0.0

        ! 4. Unallocated array, scalar source argument
        deallocate (arr)
        source = 123.4
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == source))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "scalar source"
        arr(:,:,:,:,:,:,:) = 0.0

        ! 5. Allocated array of wrong shape, scalar source argument
        shp = shp + 1
        source = 234.5
        stat = 1000
        call cond_alloc (arr, shp, source=source, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == shp) .and. stat == 0 .and. all (arr == source))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and scalar source"
        arr(:,:,:,:,:,:,:) = 0.0

        ! 6. Unallocated array, array source argument
        deallocate (arr)
        allocate (src(shp(1),shp(2),shp(3),shp(4),shp(5),shp(6),shp(7)), source=0.12345_PREC)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == src))) &
            error stop "COND_ALLOC check failed for an unallocated array and " // &
                "array source"
        arr(:,:,:,:,:,:,:) = 0.0

        ! 7. Allocated array of wrong shape, array source argument
        shp = shp + 10
        call cond_alloc (arr, shp)
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == 0 .and. all (arr == src))) &
            error stop "COND_ALLOC check failed for an array with the wrong shape " &
                // "and array source"

        ! 8. Allocated array of correct size, array source argument
        arr(:,:,:,:,:,:,:) = 0.0
        src(:,:,:,:,:,:,:) = 123.45
        stat = 1000
        call cond_alloc (arr, source=src, stat=stat)
        call cond_alloc (xx, source=src)
        if (.not. (all (shape (arr) == shape (src)) .and. stat == - 1 .and. all (arr == 0.0))) &
            error stop "COND_ALLOC check failed for a correctly shaped array and " // &
                "array source"

        deallocate (arr, src)
    end do

    ! Test with invalid shape argument
    shp = [10,2,3,7,1,3,1]
    do i = 0, 7
        if (i == NDIM) cycle
        allocate (shp2(i), source=1)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. ( .not. allocated (arr) .and. stat == - 2)) &
            error stop "Unallocated ARR, invalid SHP=[10, 2, 3, 7, 1, 3, 1]"

        allocate (arr(1,2,3,7,2,1,2), source=0.0_PREC)
        stat = 1000
        call cond_alloc (arr, shp2, stat=stat)
        call cond_alloc (xx, shp)
        if (.not. (all (shape (arr) == [1,2,3,7,2,1,2]) .and. stat == - 2)) &
            error stop "Allocated ARR, invalid SHP=[10, 2, 3, 7, 1, 3, 1]"

        deallocate (arr, shp2)
    end do

end subroutine

end program
