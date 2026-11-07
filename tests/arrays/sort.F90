

program test_numfort_arrays_sort

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common_status

    use numfort_common_testing

    use numfort_stats

    implicit none

    call test_all ()

    contains

subroutine test_all ()

    call test_argsort_real32_int32 ()
    call test_argsort_real32_int64 ()
    call test_argsort_real64_int32 ()
    call test_argsort_real64_int64 ()

    call test_argsort_int32 ()
    call test_argsort_int64 ()

end subroutine

subroutine test_argsort_real32_int32 ()
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    real (PREC), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine

subroutine test_argsort_real32_int64 ()
    integer, parameter :: PREC = real32
    integer, parameter :: INTSIZE = int32

    real (PREC), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine

subroutine test_argsort_real64_int32 ()
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    real (PREC), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine

subroutine test_argsort_real64_int64 ()
    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int64

    real (PREC), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine

subroutine test_argsort_int32 ()
    integer, parameter :: INTSIZE = int32

    integer (INTSIZE), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine

subroutine test_argsort_int64 ()
    integer, parameter :: INTSIZE = int64

    integer (INTSIZE), dimension(:), allocatable :: arr, arr_sorted
    real (real64), dimension(:), allocatable :: rarr
    integer (INTSIZE), dimension(:), allocatable :: iorder
    integer :: i, k, n
    integer, parameter :: SIZES(*) = [0, 1, 5, 10, 100, 1000, 1000]
    integer (NF_ENUM_KIND) :: status
    logical :: values_ok

    call set_seed (123)

#include "test_argsort_impl.F90"
end subroutine

end
