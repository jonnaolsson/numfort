program test_nf_arrays_setops
    !*   Unit tests for array set routines.

    use, intrinsic :: iso_fortran_env
    use numfort_arrays

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

contains

subroutine test_all ()

    call test_setdiff ()
    call test_intersect ()
    call test_union ()
    call test_in ()

end subroutine

subroutine test_setdiff ()

    real (PREC), dimension(:), allocatable :: set_a, set_b, diff
    integer, dimension(:), allocatable :: idx
    integer :: n

    ! 1. Degenerate inputs A, B
    allocate (set_a(0), set_b(0), diff(0), idx(0))
    call setdiff (set_a, set_b, n, diff, idx)
    if (.not. (n == 0)) error stop "Zero-size input sets A, B"
    deallocate (set_a, set_b, diff, idx)

    ! 2. Degenerate input B, some set A
    allocate (set_a(1), set_b(0), diff(1), idx(1))
    set_a = 1.0
    call setdiff (set_a, set_b, n, diff, idx)
    if (.not. (n == 1 .and. all (set_a == diff) .and. all (idx == 1))) &
        error stop "Singleton A, zero-size B"
    deallocate (set_a, set_b, diff, idx)

    ! Degenerate set B, non-sorted set A with duplicates
    allocate (set_a(5), set_b(0), diff(5), idx(5))
    set_a(:) = [real (PREC) :: 1, 1, 5, 4, 7]
    idx = 0
    diff = 0
    call setdiff (set_a, set_b, n, diff, idx)
    ! unclear whether the underlying routine in ORDERPACK always
    ! returns index of last occurence of duplicate entries
    if (.not. (n == 4 .and. all (diff(1:4) == [1,4,5,7]) .and. all (idx(2:4) == &
        [4,3,5]) .and. (idx(1) == 1 .or. idx(1) == 2))) &
        error stop "Zero-size B, unsorted A with duplicates"
    deallocate (set_a, set_b, diff, idx)

    ! 3. Degenerate input A, some set B
    allocate (set_a(0), set_b(1), diff(1), idx(1))
    set_b = 1.0
    call setdiff (set_a, set_b, n, diff, idx)
    if (.not. (n == 0)) error stop "Zero-size A, singleton B"
    deallocate (set_a, set_b, diff, idx)

    ! 4. Identical non-degenerate set A, B
    allocate (set_a(5), set_b(5), diff(5), idx(5))
    call linspace (set_a, 1.0_PREC, 5.0_PREC)
    set_b(:) = set_a
    call setdiff (set_a, set_b, n, diff, idx)
    if (.not. (n == 0)) error stop "Non-degenerate identical A, B"
    deallocate (set_a, set_b, diff, idx)

    ! 5. B is strict superset of A
    allocate (set_a(3), set_b(5), diff(5), idx(5))
    call linspace (set_b, 1.0_PREC, 5.0_PREC)
    set_a(:) = set_b(2:4)
    call setdiff (set_a, set_b, n, diff, idx)
    if (.not. (n == 0)) error stop "Non-degenerate A, B; A strict subset of B"
    deallocate (set_a, set_b, diff, idx)

    ! 6. Some example that yields non-empty A\B
    allocate (set_a(5), set_b(3), diff(5), idx(5))
    diff = 0
    idx = 0
    call linspace (set_a, 1.0_PREC, 5.0_PREC)
    ! Create non-sorted input array
    set_a([1, 5]) = set_a([5, 1])
    set_b(:) = set_a(2:4)
    call setdiff (set_a, set_b, n, diff, idx)
    if (.not. (n == 2 .and. all (diff(1:2) == [1,5]) .and. all (idx(1:2) == [5,1]))) &
        error stop "Non-generate A, B; B strict subset of A"

    ! 7. Test with assume_unique = .true. using the same input arrays
    diff = 0
    idx = 0
    n = 0
    call setdiff (set_a, set_b, n, diff, idx, assume_unique=.true.)
    if (.not. (n == 2 .and. all (diff(1:2) == [1,5]) .and. all (idx(1:2) == [5,1]))) &
        error stop "Non-degenerate A, B; B strict subset of A; assumed unique"
    deallocate (set_a, set_b, diff, idx)

    ! 8. Test with B strict subset of A, B contains initial sorted elements
    ! of A.
    allocate (set_a(5), set_b(3), diff(2), idx(2))
    diff = 0
    idx = 0
    call linspace (set_a, 1.0_PREC, 5.0_PREC)
    set_b(1:3) = set_a(1:3)
    call setdiff (set_a, set_b, n, diff, idx, assume_unique=.true.)
    if (.not. (n == 2 .and. all (diff(1:2) == [4,5]) .and. all (idx(1:2) == [4,5]))) &
        error stop "B strict subset of A; B contains initial elems of A"

    ! Test with B containing last elements of A
    diff = 0
    idx = 0
    n = 0
    set_b(1:3) = set_a(3:5)
    call setdiff (set_a, set_b, n, diff, idx, assume_unique=.true.)
    if (.not. (n == 2 .and. all (diff(1:2) == [1,2]) .and. all (idx(1:2) == [1,2]))) &
        error stop "B strict subset of A; B contains last elemts of A"
    deallocate (diff, idx, set_a, set_b)

    ! Test with A, B disjoint
    allocate (set_a(3), set_b(5), idx(3), diff(3))
    diff = 0
    idx = 0
    n = 0
    set_b = [1,3,5,7,9]
    set_a = [0,2,4]
    call setdiff (set_a, set_b, n, diff, idx)
    if (.not. (n == size (set_a) .and. all (diff == set_a) .and. all (idx == &
        [1,2,3]))) &
        error stop "A, B disjoint"
    deallocate (set_a, set_b, diff, idx)

end subroutine

subroutine test_intersect ()

    real (PREC), dimension(:), allocatable :: set_a, set_b, res
    integer :: n, i

    ! 1. Both operands are empty sets
    allocate (set_a(0), set_b(0), res(0))
    n = -1
    call intersect (set_a, set_b, res, n)
    if (.not. (n == 0)) error stop "A, B empty sets"
    deallocate (set_a, set_b, res)

    ! 2. First operand empty, second operand non-empty
    allocate (set_a(0), set_b(10), res(10))
    set_b = [(i, i = 1, 10)]
    n = -1
    call intersect (set_a, set_b, res, n)
    if (.not. (n == 0)) error stop "A empty, B non-empty"
    deallocate (set_a, set_b, res)

    ! 3. Second operand empty, first non-empty
    allocate (set_a(10), set_b(0), res(0))
    set_a = [(i, i = 1, 10)]
    n = -1
    call intersect (set_a, set_b, res, n)
    if (.not. (n == 0)) error stop "A non-empty, B empty"
    deallocate (set_a, set_b, res)

    ! 4. Nonempty operands, empty intersection
    allocate (set_a(5), set_b(5), res(0))
    set_a = [(i, i=1,5)]
    set_b = [(i, i=6,10)]
    n = -1
    call intersect (set_a, set_b, res, n)
    if (.not. (n == 0)) error stop "A, B non-empty, empty intersection"
    deallocate (set_a, set_b, res)

    ! 5. A strict subset of B
    allocate (set_a(5), set_b(10), res(5))
    set_a = [(i, i = 1, 5)]
    set_b = [(i, i = 1, 10)]
    n = -1
    call intersect (set_a, set_b, res, n)
    if (.not. (n == 5 .and. all (res(1:n) == set_a))) error stop "A, B nonempty, C == A"
    deallocate (set_a, set_b, res)

    ! 6. B strict subset of A
    allocate (set_a(10), set_b(5), res(5))
    set_a = [(i, i = 1, 10)]
    set_b = [(i, i = 6, 10)]
    n = -1
    call intersect (set_a, set_b, res, n)
    if (.not. (n == 5 .and. all (res(1:n) == set_b))) error stop "A, B nonempty, C == B"
    deallocate (set_a, set_b, res)

    ! 7. Test some random non-empty intersection
    allocate (set_a(5), set_b(5), res(3))
    set_a = [(i, i = 1, 5)]
    set_b = [2, 3, 5, 7, 10]
    n = -1
    call intersect (set_a, set_b, res, n)
    if (.not. (n == 3 .and. all (res(1:n) == [2,3,5]))) &
        error stop "Non-empty A, B; non-empty intersection"

    ! Test the same if assuming that elements are unique
    n = -1
    call intersect (set_a, set_b, res, n)
    if (.not. (n == 3 .and. all (res(1:n) == [2,3,5]))) &
        error stop "Non-empty A, B; non-empty intersection; assumed unique"
    deallocate (set_a, set_b, res)

    ! 8. Test with arrays containing duplicates
    allocate (set_a(10), set_b(10), res(10))
    set_a = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5]
    set_b = [2, 2, 3, 4, 4, 6, 6, 7, 8, 9]
    n = -1
    call intersect (set_a, set_b, res, n)
    if (.not. (n == 3 .and. all (res(1:n) == [2,3,4]))) &
        error stop "A, B containing duplicates, non-empty intersection"
    deallocate (set_a, set_b, res)

end subroutine

subroutine test_union ()

    real (PREC), dimension(:), allocatable :: set_a, set_b, res
    integer :: n, i

    ! 1. Both operands are empty sets
    allocate (set_a(0), set_b(0), res(0))
    n = -1
    call union (set_a, set_b, res, n)
    if (.not. (n == 0)) error stop "A, B empty sets"
    deallocate (set_a, set_b, res)

    ! 2. First operand empty, second operand non-empty
    allocate (set_a(0), set_b(10), res(10))
    set_b = [(i, i = 1, 10)]
    n = -1
    call union (set_a, set_b, res, n)
    if (.not. (n == size (set_b) .and. all (res(1:n) == set_b))) &
        error stop "A empty, B non-empty"
    deallocate (set_a, set_b, res)

    ! 3. Second operand empty, first non-empty
    allocate (set_a(10), set_b(0), res(10))
    set_a = [(i, i = 1, 10)]
    n = -1
    call union (set_a, set_b, res, n)
    if (.not. (n == size (set_a) .and. all (res(1:n) == set_a))) &
        error stop "A non-empty, B empty"
    deallocate (set_a, set_b, res)

    ! 4. Nonempty operands
    allocate (set_a(5), set_b(5), res(10))
    set_a = [(i, i=1,5)]
    set_b = [(i, i=6,10)]
    n = -1
    call union (set_a, set_b, res, n)
    if (.not. (n == 10 .and. all (res(1:n) == [(i,i = 1,10)]))) &
        error stop "A, B non-empty"
    deallocate (set_a, set_b, res)

    ! 5. A strict subset of B
    allocate (set_a(5), set_b(10), res(10))
    set_a = [(i, i = 1, 5)]
    set_b = [(i, i = 1, 10)]
    n = -1
    call union (set_a, set_b, res, n)
    if (.not. (n == size (set_b) .and. all (res(1:n) == set_b))) &
        error stop "A, B nonempty, C == B"
    deallocate (set_a, set_b, res)

    ! 6. B strict subset of A
    allocate (set_a(10), set_b(5), res(10))
    set_a = [(i, i = 1, 10)]
    set_b = [(i, i = 6, 10)]
    n = -1
    call union (set_a, set_b, res, n)
    if (.not. (n == size (set_a) .and. all (res(1:n) == set_a))) &
        error stop "A, B nonempty, C == A"
    deallocate (set_a, set_b, res)

    ! 7. Test some random non-empty unionion, res array not large enough
    ! no hold union
    allocate (set_a(5), set_b(5), res(3))
    set_a = [(i, i = 1, 5)]
    set_b = [2, 3, 5, 7, 10]
    n = -1
    call union (set_a, set_b, res, n)
    if (.not. (n == 7 .and. all (res == [1,2,3]))) &
        error stop "Non-empty A, B; result array too small to hold union"

    ! Test the same if assuming that elements are unique
    n = -1
    res = 0
    call union (set_a, set_b, res, n)
    if (.not. (n == 7 .and. all (res == [1,2,3]))) &
        error stop "Non-empty A, B; result array too small to hold union; assumed unique"
    deallocate (set_a, set_b, res)

    ! 8. Test with arrays containing duplicates
    allocate (set_a(10), set_b(10), res(10))
    set_a = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5]
    set_b = [2, 2, 3, 4, 4, 6, 6, 7, 8, 9]
    n = -1
    call union (set_a, set_b, res, n)
    if (.not. (n == 9 .and. all (res(1:n) == [1,2,3,4,5,6,7,8,9]))) &
        error stop "A, B containing duplicates"
    deallocate (set_a, set_b, res)

end subroutine

subroutine test_in ()

    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), dimension(:), allocatable :: haystack
    integer (INTSIZE) :: needle
    logical :: res

    allocate (haystack(0))
    needle = 1
    res = needle .in. haystack
    if (res) error stop "size-zero HAYSTACK"
    deallocate (haystack)

    allocate (haystack(1), source=0_INTSIZE)
    needle = 1
    res = needle .in. haystack
    if (res) error stop "size-one HAYSTACK, NEEDLE not in HAYSTACK"

    needle = 0
    res = needle .in. haystack
    if (.not. (res)) error stop "size-one HAYSTACK, NEEDLE in HAYSTACK"
    deallocate (haystack)

    allocate (haystack(5))
    call arange (haystack, 1_INTSIZE)
    needle = 0
    res = needle .in. haystack
    if (res) error stop "sorted HAYSTACK, NEEDLE not in HAYSTACK"

    needle = 3
    res = needle .in. haystack
    if (.not. (res)) error stop "sorted HAYSTACK, NEEDLE in HAYSTACK"
    deallocate (haystack)

end subroutine

end program
