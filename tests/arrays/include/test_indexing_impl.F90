subroutine __APPEND(test_indexing,__INTSIZE) ()
    !*  Run all indexing unit tests for __INTSIZE.

    call __APPEND(test_ind2sub,__INTSIZE) ()
    call __APPEND(test_sub2ind,__INTSIZE) ()
    call __APPEND(test_indexing_1d,__INTSIZE) ()
    call __APPEND(test_shape2sub,__INTSIZE) ()
    call __APPEND(test_roundtrip,__INTSIZE) ()

end subroutine


subroutine __APPEND(test_ind2sub,__INTSIZE) ()
    !*  Direct unit tests for ind2sub with __INTSIZE.

    integer, parameter :: INTSIZE = __INTSIZE

    integer (INTSIZE), dimension(:), allocatable :: lin, shp
    integer (INTSIZE), dimension(:, :), allocatable :: sub, sub_ref
    integer (INTSIZE) :: i, n, idx, k1, k2, k3

    ! 1D array with rank-2 sub output
    n = 100_INTSIZE
    allocate (lin(n), sub(n, 1))
    lin = [(i, i=1, n)]
    sub = -1_INTSIZE

    call ind2sub ([n], lin, sub)
    if (.not. (all (lin == sub(:, 1)))) &
        error stop "ind2sub failed for 1D input"
    deallocate (lin, sub)

    ! 2D array [7, 8]
    shp = [7_INTSIZE, 8_INTSIZE]
    n = product (shp)
    allocate (lin(n), sub(n, 2), sub_ref(n, 2))
    lin = [(i, i=1, n)]

    call ind2sub (shp, lin, sub)

    sub_ref(:, 1) = modulo (lin - 1_INTSIZE, shp(1)) + 1_INTSIZE
    sub_ref(:, 2) = (lin - 1_INTSIZE) / shp(1) + 1_INTSIZE

    if (.not. (all (sub == sub_ref))) &
        error stop "ind2sub failed for shape [7, 8]"
    deallocate (lin, sub, sub_ref, shp)

    ! 3D array [3, 4, 5]
    shp = [3_INTSIZE, 4_INTSIZE, 5_INTSIZE]
    n = product (shp)
    allocate (lin(n), sub(n, 3), sub_ref(n, 3))
    lin = [(i, i=1, n)]

    idx = 1_INTSIZE
    do k3 = 1_INTSIZE, shp(3)
        do k2 = 1_INTSIZE, shp(2)
            do k1 = 1_INTSIZE, shp(1)
                sub_ref(idx, 1) = k1
                sub_ref(idx, 2) = k2
                sub_ref(idx, 3) = k3
                idx = idx + 1_INTSIZE
            end do
        end do
    end do

    call ind2sub (shp, lin, sub)
    if (.not. (all (sub == sub_ref))) &
        error stop "ind2sub failed for shape [3, 4, 5]"
    deallocate (lin, sub, sub_ref, shp)

    ! Singleton dimension: [1]
    shp = [1_INTSIZE]
    allocate (lin(1), sub(1, 1))
    lin = [1_INTSIZE]
    call ind2sub (shp, lin, sub)
    if (.not. (sub(1, 1) == 1_INTSIZE)) &
        error stop "ind2sub failed for shape [1]"
    deallocate (lin, sub, shp)

    ! Singleton dimension: [1, 1]
    shp = [1_INTSIZE, 1_INTSIZE]
    allocate (lin(1), sub(1, 2))
    lin = [1_INTSIZE]
    call ind2sub (shp, lin, sub)
    if (.not. (all (sub(1, :) == [1_INTSIZE, 1_INTSIZE]))) &
        error stop "ind2sub failed for shape [1, 1]"
    deallocate (lin, sub, shp)

    ! Singleton dimension: [1, 5]
    shp = [1_INTSIZE, 5_INTSIZE]
    n = product (shp)
    allocate (lin(n), sub(n, 2), sub_ref(n, 2))
    lin = [(i, i=1, n)]
    sub_ref(:, 1) = 1_INTSIZE
    sub_ref(:, 2) = lin
    call ind2sub (shp, lin, sub)
    if (.not. (all (sub == sub_ref))) &
        error stop "ind2sub failed for shape [1, 5]"
    deallocate (lin, sub, sub_ref, shp)

    ! Singleton dimension: [5, 1]
    shp = [5_INTSIZE, 1_INTSIZE]
    n = product (shp)
    allocate (lin(n), sub(n, 2), sub_ref(n, 2))
    lin = [(i, i=1, n)]
    sub_ref(:, 1) = lin
    sub_ref(:, 2) = 1_INTSIZE
    call ind2sub (shp, lin, sub)
    if (.not. (all (sub == sub_ref))) &
        error stop "ind2sub failed for shape [5, 1]"
    deallocate (lin, sub, sub_ref, shp)

    ! Singleton interior dimension: [3, 1, 4]
    shp = [3_INTSIZE, 1_INTSIZE, 4_INTSIZE]
    n = product (shp)
    allocate (lin(n), sub(n, 3), sub_ref(n, 3))
    lin = [(i, i=1, n)]
    idx = 1_INTSIZE
    do k3 = 1_INTSIZE, shp(3)
        do k2 = 1_INTSIZE, shp(2)
            do k1 = 1_INTSIZE, shp(1)
                sub_ref(idx, 1) = k1
                sub_ref(idx, 2) = k2
                sub_ref(idx, 3) = k3
                idx = idx + 1_INTSIZE
            end do
        end do
    end do
    call ind2sub (shp, lin, sub)
    if (.not. (all (sub == sub_ref))) &
        error stop "ind2sub failed for shape [3, 1, 4]"
    deallocate (lin, sub, sub_ref, shp)

end subroutine


subroutine __APPEND(test_sub2ind,__INTSIZE) ()
    !*  Direct unit tests for sub2ind with __INTSIZE.

    integer, parameter :: INTSIZE = __INTSIZE

    integer (INTSIZE), dimension(:), allocatable :: lin, lin_ref, shp
    integer (INTSIZE), dimension(:, :), allocatable :: sub, sub_ref
    integer (INTSIZE) :: i, n, idx, k1, k2, k3
    integer (INTSIZE), dimension(4, 2) :: sub_pts
    integer (INTSIZE), dimension(4) :: lin_pts

    ! 1D array with rank-2 sub input
    n = 100_INTSIZE
    allocate (lin(n), lin_ref(n), sub(n, 1))
    lin_ref = [(i, i=1, n)]
    sub(:, 1) = lin_ref
    call sub2ind ([n], sub, lin)
    if (.not. (all (lin == lin_ref))) &
        error stop "sub2ind failed for 1D input"
    deallocate (lin, lin_ref, sub)

    ! 2D array [7, 8]
    shp = [7_INTSIZE, 8_INTSIZE]
    n = product (shp)
    allocate (lin(n), lin_ref(n), sub_ref(n, 2))
    lin_ref = [(i, i=1, n)]
    sub_ref(:, 1) = modulo (lin_ref - 1_INTSIZE, shp(1)) + 1_INTSIZE
    sub_ref(:, 2) = (lin_ref - 1_INTSIZE) / shp(1) + 1_INTSIZE

    call sub2ind (shp, sub_ref, lin)
    if (.not. (all (lin == lin_ref))) &
        error stop "sub2ind failed for shape [7, 8]"
    deallocate (lin, lin_ref, sub_ref, shp)

    ! 3D array [3, 4, 5]
    shp = [3_INTSIZE, 4_INTSIZE, 5_INTSIZE]
    n = product (shp)
    allocate (lin(n), lin_ref(n), sub_ref(n, 3))
    lin_ref = [(i, i=1, n)]
    idx = 1_INTSIZE
    do k3 = 1_INTSIZE, shp(3)
        do k2 = 1_INTSIZE, shp(2)
            do k1 = 1_INTSIZE, shp(1)
                sub_ref(idx, 1) = k1
                sub_ref(idx, 2) = k2
                sub_ref(idx, 3) = k3
                idx = idx + 1_INTSIZE
            end do
        end do
    end do
    call sub2ind (shp, sub_ref, lin)
    if (.not. (all (lin == lin_ref))) &
        error stop "sub2ind failed for shape [3, 4, 5]"
    deallocate (lin, lin_ref, sub_ref, shp)

    ! Singleton shapes: [1], [1, 1], [1, 5], [5, 1], [3, 1, 4]
    shp = [1_INTSIZE]
    allocate (lin(1), sub(1, 1))
    sub(1, 1) = 1_INTSIZE
    call sub2ind (shp, sub, lin)
    if (.not. (lin(1) == 1_INTSIZE)) &
        error stop "sub2ind failed for shape [1]"
    deallocate (lin, sub, shp)

    shp = [1_INTSIZE, 1_INTSIZE]
    allocate (lin(1), sub(1, 2))
    sub(1, :) = [1_INTSIZE, 1_INTSIZE]
    call sub2ind (shp, sub, lin)
    if (.not. (lin(1) == 1_INTSIZE)) &
        error stop "sub2ind failed for shape [1, 1]"
    deallocate (lin, sub, shp)

    shp = [1_INTSIZE, 5_INTSIZE]
    n = product (shp)
    allocate (lin(n), lin_ref(n), sub_ref(n, 2))
    lin_ref = [(i, i=1, n)]
    sub_ref(:, 1) = 1_INTSIZE
    sub_ref(:, 2) = lin_ref
    call sub2ind (shp, sub_ref, lin)
    if (.not. (all (lin == lin_ref))) &
        error stop "sub2ind failed for shape [1, 5]"
    deallocate (lin, lin_ref, sub_ref, shp)

    shp = [5_INTSIZE, 1_INTSIZE]
    n = product (shp)
    allocate (lin(n), lin_ref(n), sub_ref(n, 2))
    lin_ref = [(i, i=1, n)]
    sub_ref(:, 1) = lin_ref
    sub_ref(:, 2) = 1_INTSIZE
    call sub2ind (shp, sub_ref, lin)
    if (.not. (all (lin == lin_ref))) &
        error stop "sub2ind failed for shape [5, 1]"
    deallocate (lin, lin_ref, sub_ref, shp)

    shp = [3_INTSIZE, 1_INTSIZE, 4_INTSIZE]
    n = product (shp)
    allocate (lin(n), lin_ref(n), sub_ref(n, 3))
    lin_ref = [(i, i=1, n)]
    idx = 1_INTSIZE
    do k3 = 1_INTSIZE, shp(3)
        do k2 = 1_INTSIZE, shp(2)
            do k1 = 1_INTSIZE, shp(1)
                sub_ref(idx, 1) = k1
                sub_ref(idx, 2) = k2
                sub_ref(idx, 3) = k3
                idx = idx + 1_INTSIZE
            end do
        end do
    end do
    call sub2ind (shp, sub_ref, lin)
    if (.not. (all (lin == lin_ref))) &
        error stop "sub2ind failed for shape [3, 1, 4]"
    deallocate (lin, lin_ref, sub_ref, shp)

    ! Explicit coordinate points for shape [5, 10]
    shp = [5_INTSIZE, 10_INTSIZE]
    sub_pts(1, :) = [1_INTSIZE, 1_INTSIZE]
    sub_pts(2, :) = [5_INTSIZE, 1_INTSIZE]
    sub_pts(3, :) = [1_INTSIZE, 2_INTSIZE]
    sub_pts(4, :) = [5_INTSIZE, 10_INTSIZE]
    call sub2ind (shp, sub_pts, lin_pts)
    if (.not. (all (lin_pts == [1_INTSIZE, 5_INTSIZE, 6_INTSIZE, 50_INTSIZE]))) &
        error stop "sub2ind failed for explicit coordinate points"
    deallocate (shp)

end subroutine


subroutine __APPEND(test_indexing_1d,__INTSIZE) ()
    !*  Unit tests for 1D overloads of ind2sub and sub2ind with __INTSIZE.

    integer, parameter :: INTSIZE = __INTSIZE

    integer (INTSIZE), dimension(:), allocatable :: lin, lin_1d, sub_1d
    integer (INTSIZE) :: i, n

    n = 50_INTSIZE
    allocate (lin_1d(n), sub_1d(n), lin(n))
    lin_1d = [(i, i=1, n)]

    call ind2sub ([n], lin_1d, sub_1d)
    if (.not. (all (sub_1d == lin_1d))) &
        error stop "1D overload ind2sub failed"

    call sub2ind ([n], sub_1d, lin)
    if (.not. (all (lin == lin_1d))) &
        error stop "1D overload sub2ind failed"
    deallocate (lin_1d, sub_1d, lin)

end subroutine


subroutine __APPEND(test_shape2sub,__INTSIZE) ()
    !*  Direct unit tests for shape2sub with __INTSIZE.

    integer, parameter :: INTSIZE = __INTSIZE

    integer (INTSIZE), dimension(:), allocatable :: lin_ref, shp
    integer (INTSIZE), dimension(:, :), allocatable :: sub, sub_ref
    integer (INTSIZE) :: i, n, idx, k1, k2, k3

    ! 2D [7, 8]
    shp = [7_INTSIZE, 8_INTSIZE]
    n = product (shp)
    allocate (sub(n, 2), sub_ref(n, 2), lin_ref(n))
    lin_ref = [(i, i=1, n)]
    sub_ref(:, 1) = modulo (lin_ref - 1_INTSIZE, shp(1)) + 1_INTSIZE
    sub_ref(:, 2) = (lin_ref - 1_INTSIZE) / shp(1) + 1_INTSIZE
    deallocate (lin_ref)

    call shape2sub (shp, sub)
    if (.not. (all (sub == sub_ref))) &
        error stop "shape2sub failed for shape [7, 8]"
    deallocate (sub, sub_ref, shp)

    ! 3D with singleton: [3, 1, 4]
    shp = [3_INTSIZE, 1_INTSIZE, 4_INTSIZE]
    n = product (shp)
    allocate (sub(n, 3), sub_ref(n, 3))
    idx = 1_INTSIZE
    do k3 = 1_INTSIZE, shp(3)
        do k2 = 1_INTSIZE, shp(2)
            do k1 = 1_INTSIZE, shp(1)
                sub_ref(idx, 1) = k1
                sub_ref(idx, 2) = k2
                sub_ref(idx, 3) = k3
                idx = idx + 1_INTSIZE
            end do
        end do
    end do
    call shape2sub (shp, sub)
    if (.not. (all (sub == sub_ref))) &
        error stop "shape2sub failed for shape [3, 1, 4]"
    deallocate (sub, sub_ref, shp)

    ! Singleton [1]
    shp = [1_INTSIZE]
    allocate (sub(1, 1))
    call shape2sub (shp, sub)
    if (.not. (sub(1, 1) == 1_INTSIZE)) &
        error stop "shape2sub failed for shape [1]"
    deallocate (sub, shp)

end subroutine


subroutine __APPEND(test_roundtrip,__INTSIZE) ()
    !*  Round trip fixed point checks for __INTSIZE.

    integer, parameter :: INTSIZE = __INTSIZE

    call __APPEND(ind2sub_fp,__INTSIZE) ([11_INTSIZE])
    call __APPEND(ind2sub_fp,__INTSIZE) ([7_INTSIZE, 8_INTSIZE])
    call __APPEND(ind2sub_fp,__INTSIZE) ([3_INTSIZE, 5_INTSIZE, 15_INTSIZE])
    call __APPEND(ind2sub_fp,__INTSIZE) ([3_INTSIZE, 5_INTSIZE, 15_INTSIZE, 7_INTSIZE])
    call __APPEND(ind2sub_fp,__INTSIZE) &
        ([3_INTSIZE, 5_INTSIZE, 15_INTSIZE, 1_INTSIZE, 5_INTSIZE])

end subroutine


subroutine __APPEND(ind2sub_fp,__INTSIZE) (shp)
    !*  Check fixed point x = sub2ind(shp, ind2sub(shp, x)).

    integer, parameter :: INTSIZE = __INTSIZE

    integer (INTSIZE), intent(in), dimension(:) :: shp
        !* Shape used for the indexing round trip.

    integer (INTSIZE), dimension(:), allocatable :: lin, lin2
    integer (INTSIZE), dimension(:, :), allocatable :: sub
    integer (INTSIZE) :: i, n

    n = product (shp)
    allocate (lin(n), lin2(n), sub(n, size (shp)))
    lin = [(i, i=1, n)]

    call ind2sub (shp, lin, sub)
    call sub2ind (shp, sub, lin2)

    if (.not. (all (lin == lin2))) &
        error stop "IND2SUB/SUB2IND round trip failed"

    deallocate (lin, lin2, sub)

end subroutine
