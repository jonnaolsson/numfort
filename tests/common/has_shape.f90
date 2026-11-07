

program test_common_has_shape
    !*  Module contains unit test routines for the array INSERT routine.

    use, intrinsic :: iso_fortran_env
    use numfort_common, only: has_shape

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

    contains

subroutine test_all ()

    call test_1d ()
    call test_2d ()

end subroutine

subroutine test_1d ()
    !*  Unit tests for in-place insertion into given argument array.

    integer, parameter :: n = 2
    integer, dimension(n) :: arr
    logical :: res, res0
    integer :: i, j

    do i = 0, n+1
        do j = 0, n
            ! Test with scalar interface
            res = has_shape (arr(:j), i)
            res0 = size(arr(:j)) == i
            if (.not. (res .eqv. res0)) &
                error stop "HAS_SHAPE disagrees with a scalar shape comparison"

            ! Test with empty arrays, array interface
            res = has_shape (arr(:j), [i])
            res0 = all(shape(arr(:j)) == [i])
            if (.not. (res .eqv. res0)) &
                error stop "HAS_SHAPE disagrees with an array shape comparison"
        end do
    end do

end subroutine

subroutine test_2d ()
    !*  Unit tests for in-place insertion into given argument array.

    integer, parameter :: n = 2
    integer, dimension(n,n) :: arr
    integer, dimension(:), allocatable :: shp
    logical :: res, res0
    integer :: i, j, k, l

    do i = 0, n+1
        do j = 0, n+1
            do k = 0, n
                do l = 0, n
                    res = has_shape (arr(:k,:l), [i,j])
                    res0 = all(shape(arr(:k,:l)) == [i,j])
                    if (.not. (res .eqv. res0)) &
                        error stop "HAS_SHAPE disagrees with the two-dimensional shape"
                end do
            end do
        end do
    end do

    ! Check with incompatible shape parameter, should return .false.
    ! instead of crashing the application
    do i = 0, 5
        if (i == size(shape(arr))) cycle

        allocate (shp(i), source=0)
        res = has_shape (arr, shp)
        if (res) error stop "HAS_SHAPE accepted an incompatible shape vector"

        deallocate (shp)
    end do

end subroutine

end program
