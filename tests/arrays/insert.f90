
program test_arrays_insert
    !*  Module contains unit test routines for the array INSERT routine.

    use, intrinsic :: iso_fortran_env
    use numfort_arrays, only: insert

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_inplace ()
    call test_copy ()

end subroutine

subroutine test_inplace ()
    !*  Unit tests for in-place insertion into given argument array.

    integer, parameter :: na = 5, nv = 3
    integer :: arr(na), val(nv), arr0(na)
    integer :: i

    arr = [(i, i=1,size(arr))]
    ! store copy of original array
    arr0 = arr
    val = [(10+i, i=1,size(val))]

    ! Invalid index
    call insert (arr, 0, val)
    if (.not. (all (arr == arr0))) error stop "INSERT with idx=0"

    call insert (arr, size(arr)+2, val)
    if (.not. (all (arr == arr0))) error stop "INSERT with idx >= size(arr)+1"

    ! Insert empty VAL array
    call insert (arr, 1, val(:0))
    if (.not. (all (arr == arr0))) error stop "INSERT with empty VAL argument"

    ! Insert VAL at beginning
    call insert (arr, 1, val)
    if (.not. (all (arr(1:nv) == val) .and. all (arr(4:) == arr0(:2)))) &
        error stop "Insert non-empty VAL at the beginning"

    ! Insert VAL at end
    arr = arr0
    call insert (arr, na-nv+1, val)
    if (.not. (all (arr(:na - nv) == arr0(:na - nv)) .and. all (arr(na - nv + 1:) == val))) &
        error stop "Insert non-empty VAL at the end"

    ! Insert VAL in the interior
    arr = arr0
    call insert (arr, 2, val)
    if (.not. (all (arr(2:2 + nv - 1) == val) .and. arr(1) == arr0(1) .and. arr(na) == arr0(2))) &
        error stop "Insert non-empty VAL in the interior of ARR"

    ! Insert VAL such that it won't fit completely
    arr = arr0
    call insert (arr, 4, val)
    if (.not. (all (arr(:3) == arr0(:3)) .and. all (arr(4:) == val(:2)))) &
        error stop "Insert VAL at end that won't fit completely"

    arr = arr0
    call insert (arr(1:2), 1, val)
    if (.not. (all (arr(:2) == val(:2)) .and. all (arr(3:) == arr0(3:)))) &
        error stop "Insert VAL at the beginning that won't fit completely"

end subroutine

subroutine test_copy ()
    !*  Unit tests for inserting elements into array and storing the result
    !   in a new array.

    integer, parameter :: na = 5, nv = 3, nr = 10
    integer :: arr(na), val(nv), res(nr)
    integer :: i

    arr = [(i, i=1,size(arr))]
    ! store copy of original array
    val = [(10+i, i=1,size(val))]

    ! Invalid index
    res = -1
    call insert (arr, 0, val, res)
    if (.not. (all (res == - 1))) error stop "INSERT with idx=0"

    call insert (arr, size(arr)+2, val, res)
    if (.not. (all (res == - 1))) error stop "INSERT with idx >= size(arr)+1"

    ! Insert empty VAL array; this should just create copy of ARR
    res = -1
    call insert (arr, 1, val(:0), res)
    if (.not. (all (arr == res(1:na)) .and. all (res(na + 1:) == - 1))) &
        error stop "INSERT with empty VAL argument"

    ! Insert empty VAL array, but pass RES that is not large enough to hold
    ! copy of ARR
    call insert (arr, 1, val(:0), res(1:na-1))
    if (.not. (all (arr(1:na - 1) == res(1:na - 1)))) &
        error stop "INSERT with empty VAL argument; size(RES) < size(ARR)"

    ! Insert VAL at beginning; RES large enough to hold final array
    res = -1
    call insert (arr, 1, val, res)
    if (.not. (all (res(1:nv) == val) .and. all (res(4:nv + na) == arr) .and. all &
        (res(nv + na + 1:) == - 1))) &
        error stop "Insert VAL at the beginning; RES large enough"

    ! Insert VAL at beginning; RES NOT large enough to hold final array
    call insert (arr, 1, val, res(1:na-1))
    if (.not. (all (res(1:nv) == val) .and. all (res(nv + 1:na - 1) == arr(1)))) &
        error stop "Insert VAL at the beginning; RES is NOT large enough"

    ! Insert VAL at end; RES is large enough to hold final array
    res = -1
    call insert (arr, na-nv+1, val, res)
    if (.not. (all (arr(:na - nv) == res(:na - nv)) .and. all (res(na - nv + 1:na) == &
        val) .and. all (res(na + 1:na + nv) == arr(na - nv + 1:na)) .and. all (res(na &
        + nv + 1:) == - 1))) &
        error stop "Insert VAL at the end; RES is sufficiently large"

    ! Insert VAL at end; RES is NOT large enough to hold final array
    call insert (arr, na-nv+1, val, res(1:na+nv-1))
    if (.not. (all (arr(:na - nv) == res(:na - nv)) .and. all (res(na - nv + 1:na) == &
        val) .and. all (res(na + 1:na + nv - 1) == arr(na - nv + 1:na - 1)))) &
        error stop "Insert VAL at the end; RES is NOT sufficiently large"

    ! Insert VAL in the interior
    res = -1
    call insert (arr, 2, val, res)
    if (.not. (all (res(2:2 + nv - 1) == val) .and. res(1) == arr(1) .and. all &
        (res(nv + 2:na + nv) == arr(2:na)) .and. all (res(na + nv + 1:) == - 1))) &
        error stop "Insert VAL in the interior of ARR; RES is sufficiently large"

    ! Insert VAL such that it won't fit completely
    call insert (arr, 4, val, res(1:na))
    if (.not. (all (res(:3) == arr(:3)) .and. all (res(4:na) == val(:2)))) &
        error stop "Insert VAL at end that won't fit completely"

    ! Insert VAL at beginning such that no element of ARR will fit
    call insert (arr(1:2), 1, val, res(1:nv-1))
    if (.not. (all (res(:2) == val(:2)))) &
        error stop "Insert VAL at the beginning that won't fit completely"

    ! Call with non-defineable arguments (constants)
    call insert ([1, 2, 3], 1, [10, 11], res)
    if (.not. (all (res(:2) == [10,11]) .and. all (res(3:5) == [1,2,3]))) &
        error stop "Call with constant ARR, VAL arguments"

end subroutine

end program
