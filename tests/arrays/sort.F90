

program test_numfort_arrays_sort
    !* Test argsort:
    !  - Check sorted values, index permutations, payloads, and invalid sizes.

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common_status

    use numfort_stats

    implicit none

    integer(int32), parameter :: MERGE_SIZES(*) = [7, 8, 9, 15, 16, 17, 31, 32, 33]

    call test_argsort_real32_int32 ()
    call test_argsort_real32_int64 ()
    call test_argsort_real64_int32 ()
    call test_argsort_real64_int64 ()

    call test_argsort_int32 ()
    call test_argsort_int64 ()

    call test_argsort_small ()
    call test_argsort_three ()
    call test_argsort_four ()
    call test_argsort_five ()
    call test_argsort_sorted_boundaries ()
    call test_argsort_reversed_boundaries ()
    call test_argsort_tied_boundaries ()
    call test_argsort_equal_boundaries ()
    call test_argsort_payload ()
    call test_argsort_invalid_size ()

    contains


subroutine test_argsort_real32_int32 ()
    !* Test real32 keys with int32 indices:
    !  - Check seeded random inputs, including empty and singleton arrays.

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
    !* Test real32 keys with int64 indices:
    !  - Check seeded random inputs, including empty and singleton arrays.

    integer, parameter :: PREC = real32
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


subroutine test_argsort_real64_int32 ()
    !* Test real64 keys with int32 indices:
    !  - Check seeded random inputs, including empty and singleton arrays.

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
    !* Test real64 keys with int64 indices:
    !  - Check seeded random inputs, including empty and singleton arrays.

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
    !* Test int32 keys and indices:
    !  - Check seeded random inputs, including empty and singleton arrays.

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
    !* Test int64 keys and indices:
    !  - Check seeded random inputs, including empty and singleton arrays.

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


subroutine test_argsort_small ()
    !* Test arrays of length zero through two:
    !  - Cover empty inputs, singletons, and ordered or reversed pairs.
    !  - Preserve both indices of an equal pair.

    ! Empty keys and indices form an empty permutation.
    call check_argsort_fixture ([integer(int32) ::], [integer(int32) ::])
    ! A singleton has only one possible index.
    call check_argsort_fixture ([-2], [-2])
    ! An ordered pair needs no exchange.
    call check_argsort_fixture ([-1, 2], [-1, 2])
    ! A reversed pair exchanges both indices.
    call check_argsort_fixture ([2, -1], [-1, 2])
    ! Both indices of an equal pair must survive.
    call check_argsort_fixture ([2, 2], [2, 2])
end subroutine


subroutine test_argsort_three ()
    !* Test three-element merges:
    !  - Cover all six permutations of distinct keys.
    !  - Preserve repeated keys without requiring stable tie ordering.

    ! A sorted triple leaves the trailing singleton after the pair.
    call check_argsort_fixture ([1, 2, 3], [1, 2, 3])
    ! The trailing singleton belongs inside the ordered pair.
    call check_argsort_fixture ([1, 3, 2], [1, 2, 3])
    ! The trailing singleton belongs before the ordered pair.
    call check_argsort_fixture ([2, 3, 1], [1, 2, 3])
    ! An initial pair exchange leaves the singleton at the end.
    call check_argsort_fixture ([2, 1, 3], [1, 2, 3])
    ! An initial pair exchange places the singleton inside the pair.
    call check_argsort_fixture ([3, 1, 2], [1, 2, 3])
    ! A reversed triple also requires the initial pair exchange.
    call check_argsort_fixture ([3, 2, 1], [1, 2, 3])
    ! Equal keys span the initial pair and trailing singleton.
    call check_argsort_fixture ([2, 1, 1], [1, 1, 2])
    ! An equal initial pair follows the smaller trailing singleton.
    call check_argsort_fixture ([2, 2, 1], [1, 2, 2])
    ! Three equal keys require three distinct indices.
    call check_argsort_fixture ([2, 2, 2], [2, 2, 2])
end subroutine


subroutine test_argsort_four ()
    !* Test four-element merges:
    !  - Cover all six interleavings of two ordered pairs.
    !  - Preserve ties within and between pairs.

    ! Both elements of the first pair precede the second pair.
    call check_argsort_fixture ([1, 2, 3, 4], [1, 2, 3, 4])
    ! The ordered pairs alternate, starting with the first pair.
    call check_argsort_fixture ([1, 3, 2, 4], [1, 2, 3, 4])
    ! The second pair lies inside the first pair.
    call check_argsort_fixture ([1, 4, 2, 3], [1, 2, 3, 4])
    ! The first pair lies inside the second pair.
    call check_argsort_fixture ([2, 3, 1, 4], [1, 2, 3, 4])
    ! The ordered pairs alternate, starting with the second pair.
    call check_argsort_fixture ([2, 4, 1, 3], [1, 2, 3, 4])
    ! Both elements of the second pair precede the first pair.
    call check_argsort_fixture ([3, 4, 1, 2], [1, 2, 3, 4])
    ! A reversed array exchanges both initial pairs.
    call check_argsort_fixture ([4, 3, 2, 1], [1, 2, 3, 4])
    ! Repeated keys occur in both initial pairs.
    call check_argsort_fixture ([2, 1, 2, 1], [1, 1, 2, 2])
    ! Each initial pair contains equal keys.
    call check_argsort_fixture ([2, 2, 1, 1], [1, 1, 2, 2])
    ! Four equal keys still require four distinct indices.
    call check_argsort_fixture ([2, 2, 2, 2], [2, 2, 2, 2])
end subroutine


subroutine test_argsort_five ()
    !* Test the first merge beyond four elements:
    !  - Cover sorted, reversed, mixed, repeated, and all-equal keys.

    ! Sorted keys leave the fifth element at the end.
    call check_argsort_fixture ([-2, -1, 0, 1, 2], [-2, -1, 0, 1, 2])
    ! Reversed keys move the fifth element to the front.
    call check_argsort_fixture ([2, 1, 0, -1, -2], [-2, -1, 0, 1, 2])
    ! The fifth element belongs inside the first four sorted keys.
    call check_argsort_fixture ([2, -2, 1, -1, 0], [-2, -1, 0, 1, 2])
    ! Ties cross the four-element merge boundary.
    call check_argsort_fixture ([2, -1, 2, -1, 2], [-1, -1, 2, 2, 2])
    ! Five equal keys require five distinct indices.
    call check_argsort_fixture ([2, 2, 2, 2, 2], [2, 2, 2, 2, 2])
end subroutine


subroutine test_argsort_sorted_boundaries ()
    !* Test sorted arrays at merge boundaries:
    !  - Cover lengths immediately below, at, and above 8, 16, and 32.

    integer(int32) :: i, k, n
    integer(int32) :: keys(maxval (MERGE_SIZES))

    do k = 1, size (MERGE_SIZES)
        n = MERGE_SIZES(k)
        ! Each sorted fixture contains the integers 1 through n.
        keys(1:n) = [(i, i = 1, n)]
        call check_argsort_fixture (keys(1:n), keys(1:n))
    end do
end subroutine


subroutine test_argsort_reversed_boundaries ()
    !* Test reversed arrays at merge boundaries:
    !  - Cover lengths immediately below, at, and above 8, 16, and 32.

    integer(int32) :: i, k, n
    integer(int32) :: keys(maxval (MERGE_SIZES)), expected(maxval (MERGE_SIZES))

    do k = 1, size (MERGE_SIZES)
        n = MERGE_SIZES(k)
        ! Each reversed fixture has the same values as the sequence 1 through n.
        keys(1:n) = [(i, i = n, 1, -1)]
        expected(1:n) = [(i, i = 1, n)]
        call check_argsort_fixture (keys(1:n), expected(1:n))
    end do
end subroutine


subroutine test_argsort_tied_boundaries ()
    !* Test repeated keys at merge boundaries:
    !  - Cross several merge runs with two distinct keys.
    !  - Cover lengths immediately below, at, and above 8, 16, and 32.

    integer(int32) :: i, k, n
    integer(int32) :: keys(maxval (MERGE_SIZES)), expected(maxval (MERGE_SIZES))

    do k = 1, size (MERGE_SIZES)
        n = MERGE_SIZES(k)
        ! Each fixture alternates 1, -1, starting with 1.
        keys(1:n) = [(2 * modulo (i, 2) - 1, i = 1, n)]
        ! There are floor(n/2) negative keys, followed by the remaining ones.
        expected(1:n / 2) = -1
        expected(n / 2 + 1:n) = 1
        call check_argsort_fixture (keys(1:n), expected(1:n))
    end do
end subroutine


subroutine test_argsort_equal_boundaries ()
    !* Test all-equal arrays at merge boundaries:
    !  - Require complete permutations even when sorted values are identical.
    !  - Cover lengths immediately below, at, and above 8, 16, and 32.

    integer(int32) :: k, n
    integer(int32) :: keys(maxval (MERGE_SIZES))

    do k = 1, size (MERGE_SIZES)
        n = MERGE_SIZES(k)
        ! Every key is -2, so sortedness alone cannot detect lost observations.
        keys(1:n) = -2
        call check_argsort_fixture (keys(1:n), keys(1:n))
    end do
end subroutine


subroutine test_argsort_payload ()
    !* Test unequal weights attached to tied keys:
    !  - Preserve each support/weight pair for real64 and int32 keys.
    !  - Check independent group totals without requiring stable tie ordering.

    integer(int32), parameter :: KEYS(*) = [2, -1, 2, 0, -1, 2, 0, -1, 2]
    real(real64), parameter :: SUPPORT(*) = real (KEYS, real64) / 4.0_real64
    real(real64), parameter :: WEIGHTS(*) = &
        [1.0_real64, 2.0_real64, 4.0_real64, 8.0_real64, 16.0_real64, &
        32.0_real64, 64.0_real64, 128.0_real64, 256.0_real64]
    integer(int32) :: iorder(size (KEYS)), i, variant
    real(real64) :: sorted_support(size (KEYS)), sorted_weights(size (KEYS))
    integer(NF_ENUM_KIND) :: status

    do variant = 1, 2
        if (variant == 1) then
            ! Real support has three tied values and distinct positive weights.
            call argsort (SUPPORT, iorder, status=status)
        else
            ! Integer destination keys carry the same unequal weights.
            call argsort (KEYS, iorder, status=status)
        end if
        if (status /= NF_STATUS_OK) error stop "Weighted argsort: invalid status"
        call check_index_permutation (iorder)
        sorted_support = SUPPORT(iorder)
        sorted_weights = WEIGHTS(iorder)
        if (any (sorted_support(1:8) > sorted_support(2:9))) &
            error stop "Weighted argsort: unsorted support"
        do i = 1, size (KEYS)
            if (count (sorted_support == SUPPORT(i) &
                    .and. sorted_weights == WEIGHTS(i)) /= 1) &
                error stop "Weighted argsort: lost or duplicated key/weight pair"
        end do

        ! Sorted group totals are 2+16+128, 8+64, and 1+4+32+256.
        if (sum (sorted_weights(1:3)) /= 146.0_real64 &
                .or. sum (sorted_weights(4:5)) /= 72.0_real64 &
                .or. sum (sorted_weights(6:9)) /= 293.0_real64) &
            error stop "Weighted argsort: incorrect group weights"
    end do
end subroutine


subroutine test_argsort_invalid_size ()
    !* Test mismatched output lengths:
    !  - Reject shorter and longer index arrays for real64 and int32 keys.
    !  - Reject mismatches involving empty inputs or outputs.

    integer(int32), parameter :: KEYS(*) = [2, 1, 2]
    real(real64), parameter :: SUPPORT(*) = [0.5_real64, 0.25_real64, 0.5_real64]
    integer(int32) :: iorder(4)
    integer(NF_ENUM_KIND) :: status

    ! Real keys with a shorter output array.
    call argsort (SUPPORT, iorder(1:2), status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Real argsort accepted short output"
    ! Real keys with a longer output array.
    call argsort (SUPPORT, iorder, status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Real argsort accepted long output"
    ! Empty real input with nonempty output.
    call argsort (SUPPORT(1:0), iorder(1:1), status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Real argsort accepted output for empty input"
    ! Nonempty real input with empty output.
    call argsort (SUPPORT(1:1), iorder(1:0), status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Real argsort accepted empty output"

    ! Integer keys with a shorter output array.
    call argsort (KEYS, iorder(1:2), status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Integer argsort accepted short output"
    ! Integer keys with a longer output array.
    call argsort (KEYS, iorder, status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Integer argsort accepted long output"
    ! Empty integer input with nonempty output.
    call argsort (KEYS(1:0), iorder(1:1), status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Integer argsort accepted output for empty input"
    ! Nonempty integer input with empty output.
    call argsort (KEYS(1:1), iorder(1:0), status=status)
    if (status /= NF_STATUS_INVALID_ARG) &
        error stop "Integer argsort accepted empty output"
end subroutine


subroutine check_argsort_fixture (keys, expected)
    !* Check one fixture for int32 keys and fractional real64 keys.
    !  Equal keys may appear in any order; argsort does not document stability.

    integer(int32), intent(in) :: keys(:)
        !* Input keys for both overloads.
    integer(int32), intent(in) :: expected(:)
        !* Independently specified sorted keys.

    integer(int32) :: iorder(size (keys))
    real(real64) :: support(size (keys)), expected_support(size (keys))
    integer(NF_ENUM_KIND) :: status

    call argsort (keys, iorder, status=status)
    if (status /= NF_STATUS_OK) error stop "Integer argsort: invalid status"
    call check_index_permutation (iorder)
    if (any (keys(iorder) /= expected)) &
        error stop "Integer argsort: incorrect sorted keys"

    ! Positive scaling preserves order; quarters are exactly representable.
    support = real (keys, real64) / 4.0_real64
    expected_support = real (expected, real64) / 4.0_real64
    call argsort (support, iorder, status=status)
    if (status /= NF_STATUS_OK) error stop "Real argsort: invalid status"
    call check_index_permutation (iorder)
    if (any (support(iorder) /= expected_support)) &
        error stop "Real argsort: incorrect sorted keys"
end subroutine


subroutine check_index_permutation (iorder)
    !* Check that indices contain every integer in 1:n exactly once.

    integer(int32), intent(in) :: iorder(:)
        !* Returned indices to validate before indexing any payload.

    integer(int32) :: i

    if (any (iorder < 1) .or. any (iorder > size (iorder))) &
        error stop "Argsort: index outside 1:n"
    do i = 1, size (iorder)
        if (count (iorder == i) /= 1) &
            error stop "Argsort: indices are not a permutation"
    end do
end subroutine


end program
