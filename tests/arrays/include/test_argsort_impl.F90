

! Seeded random fixtures cover empty, singleton, small, and larger arrays.
do k = 1, size (SIZES)
    n = SIZES(k)

    allocate (rarr(n), arr(n), arr_sorted(n))
    allocate (iorder(n), source=0_INTSIZE)

    call random_number (rarr)

    arr(:) = rarr * 1000.0

    ! ARGSORT returns rank array where the value at position i contains the index
    ! the element which would be at position i in a sorted array.
    call argsort (arr, iorder, status=status)

    if (status /= NF_STATUS_OK) error stop "ARGSORT returned an invalid status"
    if (any (iorder < 1) .or. any (iorder > n)) &
        error stop "ARGSORT returned an out-of-range index"
    do i = 1, n
        if (count (iorder == i) /= 1) &
            error stop "ARGSORT indices are not a permutation"
    end do

    do i = 1, n
        arr_sorted(i) = arr(iorder(i))
    end do

    values_ok = all (arr_sorted(1:n-1) <= arr_sorted(2:n))

    if (.not. values_ok) error stop "ARGSORT returned unsorted values"

    deallocate (rarr, arr, arr_sorted, iorder)
end do
