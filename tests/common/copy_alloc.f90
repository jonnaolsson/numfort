program numfort_common_copy_alloc
    !*  Unit tests for COPY_ALLOC routine

    use, intrinsic :: iso_fortran_env
    use numfort_common, only: shape_equal, copy_alloc

    implicit none

    integer, parameter :: PREC = real64

    ! Distinct literal entries expose misplaced axes and slices in every rank.
    real (PREC), parameter :: VALUES(12) = [ &
        2.0_PREC, -3.0_PREC, 5.0_PREC, 7.0_PREC, -11.0_PREC, 13.0_PREC, &
        17.0_PREC, -19.0_PREC, 23.0_PREC, 29.0_PREC, -31.0_PREC, 37.0_PREC]
    real (PREC), parameter :: REPLACEMENTS(12) = [ &
        -41.0_PREC, 43.0_PREC, 47.0_PREC, -53.0_PREC, 59.0_PREC, 61.0_PREC, &
        -67.0_PREC, 71.0_PREC, 73.0_PREC, -79.0_PREC, 83.0_PREC, 89.0_PREC]

    call test_missing_src ()
    call test_unalloc_dst ()
    call test_wrong_dst_shape ()
    call test_copy ()
    call test_copy_5d ()
    call test_empty_5d ()
    call test_missing_src_5d ()
    call test_copy_6d ()
    call test_empty_6d ()
    call test_missing_src_6d ()
    call test_copy_7d ()
    call test_empty_7d ()
    call test_missing_src_7d ()

    contains


subroutine test_missing_src ()

    real (PREC), allocatable :: src1(:), src2(:,:), src3(:,:,:), src4(:,:,:,:)
    real (PREC), allocatable :: dst1(:), dst2(:,:), dst3(:,:,:), dst4(:,:,:,:)

    ! 1. Unallocated SRC, unallocated DST argument
    call copy_alloc (src1, dst1)
    if (.not. ( .not. allocated (dst1))) error stop "1d: unallocated SRC/DST"

    call copy_alloc (src2, dst2)
    if (.not. ( .not. allocated (dst2))) error stop "2d: unallocated SRC/DST"

    call copy_alloc (src3, dst3)
    if (.not. ( .not. allocated (dst3))) error stop "3d: unallocated SRC/DST"

    call copy_alloc (src4, dst4)
    if (.not. ( .not. allocated (dst4))) error stop "4d: unallocated SRC/DST"

    ! 2. Missing SRC argument
    call copy_alloc (dst=dst1)
    if (.not. ( .not. allocated (dst1))) error stop "1d: missing SRC, unallocated DST"

    call copy_alloc (dst=dst2)
    if (.not. ( .not. allocated (dst2))) error stop "2d: missing SRC, unallocated DST"

    call copy_alloc (dst=dst3)
    if (.not. ( .not. allocated (dst3))) error stop "3d: missing SRC, unallocated DST"

    call copy_alloc (dst=dst4)
    if (.not. ( .not. allocated (dst4))) error stop "4d: missing SRC, unallocated DST"

end subroutine

subroutine test_unalloc_dst ()

    real (PREC), allocatable :: src1(:), src2(:,:), src3(:,:,:), src4(:,:,:,:)
    real (PREC), allocatable :: dst1(:), dst2(:,:), dst3(:,:,:), dst4(:,:,:,:)

    allocate (src1(10), source=1.0_PREC)
    call copy_alloc (src1, dst1)
    if (.not. (shape_equal (src1, dst1) .and. all (src1 == dst1))) &
        error stop "1d: allocated SRC, unalloc. DST"

    allocate (src2(10,10), source=1.0_PREC)
    call copy_alloc (src2, dst2)
    if (.not. (shape_equal (src2, dst2) .and. all (src2 == dst2))) &
        error stop "2d: allocated SRC, unalloc. DST"

    allocate (src3(10,10,10), source=1.0_PREC)
    call copy_alloc (src3, dst3)
    if (.not. (shape_equal (src3, dst3) .and. all (src3 == dst3))) &
        error stop "3d: allocated SRC, unalloc. DST"

    allocate (src4(10,10,10,10), source=1.0_PREC)
    call copy_alloc (src4, dst4)
    if (.not. (shape_equal (src4, dst4) .and. all (src4 == dst4))) &
        error stop "4d: allocated SRC, unalloc. DST"

end subroutine

subroutine test_wrong_dst_shape ()

    real (PREC), allocatable :: src1(:), src2(:,:), src3(:,:,:), src4(:,:,:,:)
    real (PREC), allocatable :: dst1(:), dst2(:,:), dst3(:,:,:), dst4(:,:,:,:)

    allocate (src1(10), source=1.0_PREC)
    allocate (dst1(5))
    call copy_alloc (src1, dst1)
    if (.not. (shape_equal (src1, dst1) .and. all (src1 == dst1))) &
        error stop "1d: allocated SRC, alloc. DST, shape(SRC) != shape(DST)"

    allocate (src2(10,10), source=1.0_PREC)
    allocate (dst2(10,5))
    call copy_alloc (src2, dst2)
    if (.not. (shape_equal (src2, dst2) .and. all (src2 == dst2))) &
        error stop "2d: allocated SRC, alloc. DST, shape(SRC) != shape(DST)"

    allocate (src3(10,10,10), source=1.0_PREC)
    allocate (dst3(10,10,1))
    call copy_alloc (src3, dst3)
    if (.not. (shape_equal (src3, dst3) .and. all (src3 == dst3))) &
        error stop "3d: allocated SRC, alloc. DST, shape(SRC) != shape(DST)"

    allocate (src4(10,10,10,10), source=1.0_PREC)
    allocate (dst4(10,10,10,20))
    call copy_alloc (src4, dst4)
    if (.not. (shape_equal (src4, dst4) .and. all (src4 == dst4))) &
        error stop "4d: allocated SRC, alloc. DST, shape(SRC) != shape(DST)"

end subroutine

subroutine test_copy ()

    real (PREC), allocatable :: src1(:), src2(:,:), src3(:,:,:), src4(:,:,:,:)
    real (PREC), allocatable :: dst1(:), dst2(:,:), dst3(:,:,:), dst4(:,:,:,:)

    allocate (src1(10), source=1.0_PREC)
    allocate (dst1(10), source=0.0_PREC)
    call copy_alloc (src1, dst1)
    if (.not. (shape_equal (src1, dst1) .and. all (src1 == dst1))) &
        error stop "1d: allocated SRC, allocated DST, shape(SRC) == shape(DST)"

    allocate (src2(10,10), source=1.0_PREC)
    allocate (dst2(10,10), source=0.0_PREC)
    call copy_alloc (src2, dst2)
    if (.not. (shape_equal (src2, dst2) .and. all (src2 == dst2))) &
        error stop "2d: allocated SRC, allocated DST, shape(SRC) == shape(DST)"

    allocate (src3(10,10,10), source=1.0_PREC)
    allocate (dst3(10,10,10), source=0.0_PREC)
    call copy_alloc (src3, dst3)
    if (.not. (shape_equal (src3, dst3) .and. all (src3 == dst3))) &
        error stop "3d: allocated SRC, allocated DST, shape(SRC) == shape(DST)"

    allocate (src4(10,10,10,10), source=1.0_PREC)
    allocate (dst4(10,10,10,10), source=0.0_PREC)
    call copy_alloc (src4, dst4)
    if (.not. (shape_equal (src4, dst4) .and. all (src4 == dst4))) &
        error stop "4d: allocated SRC, allocated DST, shape(SRC) == shape(DST)"

end subroutine


subroutine test_copy_5d ()
    !* - Copy rank-5 real64 arrays with unequal extents and singleton axes.
    !  - Check replacement, interior/trailing resizing, and independent storage.

    ! Reshape lists the fixed entries in Fortran array element order.
    real (PREC), parameter :: INITIAL(2,1,3,1,2) = &
        reshape (VALUES, [2,1,3,1,2])
    real (PREC), parameter :: REPLACEMENT(2,1,3,1,2) = &
        reshape (REPLACEMENTS, [2,1,3,1,2])
    real (PREC), parameter :: INTERIOR(2,1,3,2,2) = &
        reshape ([VALUES, REPLACEMENTS], [2,1,3,2,2])
    real (PREC), parameter :: TRAILING(2,1,3,2,1) = &
        reshape (REPLACEMENTS, [2,1,3,2,1])
    real (PREC), allocatable :: src(:,:,:,:,:), dst(:,:,:,:,:)

    ! An unallocated destination receives every entry of the source.
    allocate (src, source=INITIAL)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "5d: initial allocation"
    if (any (shape (dst) /= [2,1,3,1,2])) error stop "5d: initial shape"
    if (any (dst /= INITIAL)) error stop "5d: initial values"
    if (any (src /= INITIAL)) error stop "5d: source changed by copy"

    ! Editing the source after copying leaves every destination entry intact.
    src(:,:,:,:,:) = -99.0_PREC
    if (any (dst /= INITIAL)) error stop "5d: source edit changed destination"

    ! A destination with the same shape has all its old values replaced.
    src(:,:,:,:,:) = REPLACEMENT
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "5d: replacement allocation"
    if (any (shape (dst) /= [2,1,3,1,2])) error stop "5d: replacement shape"
    if (any (dst /= REPLACEMENT)) error stop "5d: replacement values"

    ! Grow the penultimate extent while retaining the leading/trailing extents.
    deallocate (src)
    allocate (src, source=INTERIOR)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "5d: interior resize allocation"
    if (any (shape (dst) /= [2,1,3,2,2])) error stop "5d: interior resize shape"
    if (any (dst /= INTERIOR)) error stop "5d: interior resize values"

    ! Shrink only the trailing extent and verify all remaining entries.
    deallocate (src)
    allocate (src, source=TRAILING)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "5d: trailing resize allocation"
    if (any (shape (dst) /= [2,1,3,2,1])) error stop "5d: trailing resize shape"
    if (any (dst /= TRAILING)) error stop "5d: trailing resize values"
end subroutine


subroutine test_empty_5d ()
    !* - Keep rank-5 destinations allocated for zero-size sources.
    !  - Check interior/trailing zero extents and copying back to nonempty arrays.

    real (PREC), parameter :: EXPECTED(2,1,3,1,2) = &
        reshape (VALUES, [2,1,3,1,2])
    real (PREC), allocatable :: src(:,:,:,:,:), dst(:,:,:,:,:)

    ! An interior zero extent replaces a populated destination with an empty one.
    allocate (src(2,1,3,0,2))
    allocate (dst, source=EXPECTED)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "5d: empty interior allocation"
    if (any (shape (dst) /= [2,1,3,0,2])) error stop "5d: empty interior shape"
    if (size (dst) /= 0) error stop "5d: empty interior size"

    ! A trailing zero extent produces an allocated empty destination from scratch.
    deallocate (src, dst)
    allocate (src(2,1,3,1,0))
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "5d: empty trailing allocation"
    if (any (shape (dst) /= [2,1,3,1,0])) error stop "5d: empty trailing shape"
    if (size (dst) /= 0) error stop "5d: empty trailing size"

    ! Copying a nonempty source resizes the previously empty destination.
    deallocate (src)
    allocate (src, source=EXPECTED)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "5d: empty to nonempty allocation"
    if (any (shape (dst) /= [2,1,3,1,2])) error stop "5d: empty to nonempty shape"
    if (any (dst /= EXPECTED)) error stop "5d: empty to nonempty values"
end subroutine


subroutine test_missing_src_5d ()
    !* - Treat unallocated and omitted rank-5 sources as absent.
    !  - Leave an unallocated destination alone or deallocate a populated one.

    real (PREC), allocatable :: src(:,:,:,:,:), dst(:,:,:,:,:)

    ! An unallocated source leaves an unallocated destination unallocated.
    call copy_alloc (src, dst)
    if (allocated (dst)) error stop "5d: unallocated source and destination"

    ! An unallocated source deallocates a populated destination.
    allocate (dst(2,1,3,1,2), source=1.0_PREC)
    call copy_alloc (src, dst)
    if (allocated (dst)) error stop "5d: unallocated source, allocated destination"

    ! An omitted source leaves an unallocated destination unallocated.
    call copy_alloc (dst=dst)
    if (allocated (dst)) error stop "5d: missing source, unallocated destination"

    ! An omitted source deallocates a populated destination.
    allocate (dst(2,1,3,1,2), source=1.0_PREC)
    call copy_alloc (dst=dst)
    if (allocated (dst)) error stop "5d: missing source, allocated destination"
end subroutine


subroutine test_copy_6d ()
    !* - Copy rank-6 real64 arrays with unequal extents and singleton axes.
    !  - Check replacement, interior/trailing resizing, and independent storage.

    ! Reshape lists the fixed entries in Fortran array element order.
    real (PREC), parameter :: INITIAL(2,1,1,3,1,2) = &
        reshape (VALUES, [2,1,1,3,1,2])
    real (PREC), parameter :: REPLACEMENT(2,1,1,3,1,2) = &
        reshape (REPLACEMENTS, [2,1,1,3,1,2])
    real (PREC), parameter :: INTERIOR(2,1,1,3,2,2) = &
        reshape ([VALUES, REPLACEMENTS], [2,1,1,3,2,2])
    real (PREC), parameter :: TRAILING(2,1,1,3,2,1) = &
        reshape (REPLACEMENTS, [2,1,1,3,2,1])
    real (PREC), allocatable :: src(:,:,:,:,:,:), dst(:,:,:,:,:,:)

    ! An unallocated destination receives every entry of the source.
    allocate (src, source=INITIAL)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "6d: initial allocation"
    if (any (shape (dst) /= [2,1,1,3,1,2])) error stop "6d: initial shape"
    if (any (dst /= INITIAL)) error stop "6d: initial values"
    if (any (src /= INITIAL)) error stop "6d: source changed by copy"

    ! Editing the source after copying leaves every destination entry intact.
    src(:,:,:,:,:,:) = -99.0_PREC
    if (any (dst /= INITIAL)) error stop "6d: source edit changed destination"

    ! A destination with the same shape has all its old values replaced.
    src(:,:,:,:,:,:) = REPLACEMENT
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "6d: replacement allocation"
    if (any (shape (dst) /= [2,1,1,3,1,2])) error stop "6d: replacement shape"
    if (any (dst /= REPLACEMENT)) error stop "6d: replacement values"

    ! Grow the penultimate extent while retaining the leading/trailing extents.
    deallocate (src)
    allocate (src, source=INTERIOR)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "6d: interior resize allocation"
    if (any (shape (dst) /= [2,1,1,3,2,2])) error stop "6d: interior resize shape"
    if (any (dst /= INTERIOR)) error stop "6d: interior resize values"

    ! Shrink only the trailing extent and verify all remaining entries.
    deallocate (src)
    allocate (src, source=TRAILING)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "6d: trailing resize allocation"
    if (any (shape (dst) /= [2,1,1,3,2,1])) error stop "6d: trailing resize shape"
    if (any (dst /= TRAILING)) error stop "6d: trailing resize values"
end subroutine


subroutine test_empty_6d ()
    !* - Keep rank-6 destinations allocated for zero-size sources.
    !  - Check interior/trailing zero extents and copying back to nonempty arrays.

    real (PREC), parameter :: EXPECTED(2,1,1,3,1,2) = &
        reshape (VALUES, [2,1,1,3,1,2])
    real (PREC), allocatable :: src(:,:,:,:,:,:), dst(:,:,:,:,:,:)

    ! An interior zero extent replaces a populated destination with an empty one.
    allocate (src(2,1,1,3,0,2))
    allocate (dst, source=EXPECTED)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "6d: empty interior allocation"
    if (any (shape (dst) /= [2,1,1,3,0,2])) error stop "6d: empty interior shape"
    if (size (dst) /= 0) error stop "6d: empty interior size"

    ! A trailing zero extent produces an allocated empty destination from scratch.
    deallocate (src, dst)
    allocate (src(2,1,1,3,1,0))
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "6d: empty trailing allocation"
    if (any (shape (dst) /= [2,1,1,3,1,0])) error stop "6d: empty trailing shape"
    if (size (dst) /= 0) error stop "6d: empty trailing size"

    ! Copying a nonempty source resizes the previously empty destination.
    deallocate (src)
    allocate (src, source=EXPECTED)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "6d: empty to nonempty allocation"
    if (any (shape (dst) /= [2,1,1,3,1,2])) error stop "6d: empty to nonempty shape"
    if (any (dst /= EXPECTED)) error stop "6d: empty to nonempty values"
end subroutine


subroutine test_missing_src_6d ()
    !* - Treat unallocated and omitted rank-6 sources as absent.
    !  - Leave an unallocated destination alone or deallocate a populated one.

    real (PREC), allocatable :: src(:,:,:,:,:,:), dst(:,:,:,:,:,:)

    ! An unallocated source leaves an unallocated destination unallocated.
    call copy_alloc (src, dst)
    if (allocated (dst)) error stop "6d: unallocated source and destination"

    ! An unallocated source deallocates a populated destination.
    allocate (dst(2,1,1,3,1,2), source=1.0_PREC)
    call copy_alloc (src, dst)
    if (allocated (dst)) error stop "6d: unallocated source, allocated destination"

    ! An omitted source leaves an unallocated destination unallocated.
    call copy_alloc (dst=dst)
    if (allocated (dst)) error stop "6d: missing source, unallocated destination"

    ! An omitted source deallocates a populated destination.
    allocate (dst(2,1,1,3,1,2), source=1.0_PREC)
    call copy_alloc (dst=dst)
    if (allocated (dst)) error stop "6d: missing source, allocated destination"
end subroutine


subroutine test_copy_7d ()
    !* - Copy rank-7 real64 arrays with unequal extents and singleton axes.
    !  - Check replacement, interior/trailing resizing, and independent storage.

    ! Reshape lists the fixed entries in Fortran array element order.
    real (PREC), parameter :: INITIAL(2,1,1,1,3,1,2) = &
        reshape (VALUES, [2,1,1,1,3,1,2])
    real (PREC), parameter :: REPLACEMENT(2,1,1,1,3,1,2) = &
        reshape (REPLACEMENTS, [2,1,1,1,3,1,2])
    real (PREC), parameter :: INTERIOR(2,1,1,1,3,2,2) = &
        reshape ([VALUES, REPLACEMENTS], [2,1,1,1,3,2,2])
    real (PREC), parameter :: TRAILING(2,1,1,1,3,2,1) = &
        reshape (REPLACEMENTS, [2,1,1,1,3,2,1])
    real (PREC), allocatable :: src(:,:,:,:,:,:,:), dst(:,:,:,:,:,:,:)

    ! An unallocated destination receives every entry of the source.
    allocate (src, source=INITIAL)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "7d: initial allocation"
    if (any (shape (dst) /= [2,1,1,1,3,1,2])) error stop "7d: initial shape"
    if (any (dst /= INITIAL)) error stop "7d: initial values"
    if (any (src /= INITIAL)) error stop "7d: source changed by copy"

    ! Editing the source after copying leaves every destination entry intact.
    src(:,:,:,:,:,:,:) = -99.0_PREC
    if (any (dst /= INITIAL)) error stop "7d: source edit changed destination"

    ! A destination with the same shape has all its old values replaced.
    src(:,:,:,:,:,:,:) = REPLACEMENT
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "7d: replacement allocation"
    if (any (shape (dst) /= [2,1,1,1,3,1,2])) error stop "7d: replacement shape"
    if (any (dst /= REPLACEMENT)) error stop "7d: replacement values"

    ! Grow the penultimate extent while retaining the leading/trailing extents.
    deallocate (src)
    allocate (src, source=INTERIOR)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "7d: interior resize allocation"
    if (any (shape (dst) /= [2,1,1,1,3,2,2])) error stop "7d: interior resize shape"
    if (any (dst /= INTERIOR)) error stop "7d: interior resize values"

    ! Shrink only the trailing extent and verify all remaining entries.
    deallocate (src)
    allocate (src, source=TRAILING)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "7d: trailing resize allocation"
    if (any (shape (dst) /= [2,1,1,1,3,2,1])) error stop "7d: trailing resize shape"
    if (any (dst /= TRAILING)) error stop "7d: trailing resize values"
end subroutine


subroutine test_empty_7d ()
    !* - Keep rank-7 destinations allocated for zero-size sources.
    !  - Check interior/trailing zero extents and copying back to nonempty arrays.

    real (PREC), parameter :: EXPECTED(2,1,1,1,3,1,2) = &
        reshape (VALUES, [2,1,1,1,3,1,2])
    real (PREC), allocatable :: src(:,:,:,:,:,:,:), dst(:,:,:,:,:,:,:)

    ! An interior zero extent replaces a populated destination with an empty one.
    allocate (src(2,1,1,1,3,0,2))
    allocate (dst, source=EXPECTED)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "7d: empty interior allocation"
    if (any (shape (dst) /= [2,1,1,1,3,0,2])) error stop "7d: empty interior shape"
    if (size (dst) /= 0) error stop "7d: empty interior size"

    ! A trailing zero extent produces an allocated empty destination from scratch.
    deallocate (src, dst)
    allocate (src(2,1,1,1,3,1,0))
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "7d: empty trailing allocation"
    if (any (shape (dst) /= [2,1,1,1,3,1,0])) error stop "7d: empty trailing shape"
    if (size (dst) /= 0) error stop "7d: empty trailing size"

    ! Copying a nonempty source resizes the previously empty destination.
    deallocate (src)
    allocate (src, source=EXPECTED)
    call copy_alloc (src, dst)
    if (.not. allocated (dst)) error stop "7d: empty to nonempty allocation"
    if (any (shape (dst) /= [2,1,1,1,3,1,2])) error stop "7d: empty to nonempty shape"
    if (any (dst /= EXPECTED)) error stop "7d: empty to nonempty values"
end subroutine


subroutine test_missing_src_7d ()
    !* - Treat unallocated and omitted rank-7 sources as absent.
    !  - Leave an unallocated destination alone or deallocate a populated one.

    real (PREC), allocatable :: src(:,:,:,:,:,:,:), dst(:,:,:,:,:,:,:)

    ! An unallocated source leaves an unallocated destination unallocated.
    call copy_alloc (src, dst)
    if (allocated (dst)) error stop "7d: unallocated source and destination"

    ! An unallocated source deallocates a populated destination.
    allocate (dst(2,1,1,1,3,1,2), source=1.0_PREC)
    call copy_alloc (src, dst)
    if (allocated (dst)) error stop "7d: unallocated source, allocated destination"

    ! An omitted source leaves an unallocated destination unallocated.
    call copy_alloc (dst=dst)
    if (allocated (dst)) error stop "7d: missing source, unallocated destination"

    ! An omitted source deallocates a populated destination.
    allocate (dst(2,1,1,1,3,1,2), source=1.0_PREC)
    call copy_alloc (dst=dst)
    if (allocated (dst)) error stop "7d: missing source, allocated destination"
end subroutine

end program
