program numfort_common_copy_alloc
    !*  Unit tests for COPY_ALLOC routine

    use, intrinsic :: iso_fortran_env
    use numfort_common, only: shape_equal, copy_alloc

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_missing_src ()
    call test_unalloc_dst ()
    call test_wrong_dst_shape ()
    call test_copy ()

end subroutine

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

end program
