

subroutine __APPEND(test_fixed,__PREC) ()
    integer, parameter :: PREC = __PREC

    real (PREC), dimension(:), allocatable :: rwork1d
    real (PREC), dimension(:), allocatable :: dat1d_out, dat1d_in
    integer, dimension(:), allocatable :: idat1d_out, idat1d_in
    real (PREC), dimension(:,:), allocatable :: dat2d_out, dat2d_in
    real (PREC), dimension(:,:,:), allocatable :: dat3d_out, dat3d_in
    integer :: n, m, k
    character (100) :: iomsg
    character (1024) :: tmpdir, path
    character (*), parameter :: filename = "test_numfort_io_fixed.txt"
    character (20) :: fmt
    type (status_t) :: status
    logical :: values_ok

    if (PREC == real32) then
    else if (PREC == real64) then
    end if

    call get_temp_directory (tmpdir, status)
    path = trim(tmpdir) // "/" // filename

    call set_seed (1234)

    ! === Test with REAL 1d arrays ===
    m = 50
    allocate (dat1d_in(m), dat1d_out(m))

    call random_number (dat1d_out)

    fmt = "(*(5(f8.4),:,/))"
    call write_fixed (path, fmt, dat1d_out, msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Writing 1d REAL data in flat form, reshaping via FMT"

    call read_fixed (path, fmt, dat1d_in, msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dat1d_out, dat1d_in, atol = &
        1.0e-3_PREC))) &
        error stop "Reading 1d REAL data in flat form, reshaping via FMT"

    deallocate (dat1d_in, dat1d_out)

    ! === Test with INTEGER 1d arrays ===
    m = 50
    allocate (idat1d_in(m), idat1d_out(m), rwork1d(m))

    call random_number (rwork1d)
    idat1d_out(:) = int(100.0d0 * rwork1d)

    fmt = "(*(5(i3),:,/))"
    call write_fixed (path, fmt, idat1d_out, msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Writing 1d INTEGER data in flat form, reshaping via FMT"

    call read_fixed (path, fmt, idat1d_in, msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK .and. all (idat1d_out == idat1d_in))) &
        error stop "Reading 1d INTEGER data in flat form, reshaping via FMT"

    deallocate (idat1d_out, idat1d_in, rwork1d)

    ! === Test with 2d arrays ===
    n = 10
    m = 11
    allocate (dat2d_out(m,n), dat2d_in(m,n))

    call random_number (dat2d_out)

    ! Write in transposed form
    call write_fixed (path, "(*(f10.6))", dat2d_out, transform="transpose", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "Writing file in transposed order"

    call read_fixed (path, "(*(f10.6))", dat2d_in, transform="transpose", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dat2d_out, dat2d_in, atol = &
        1.0e-5_PREC))) &
        error stop "Reading written file in transposed order"

    ! Write in non-transposed order such that one line in the file corresponds
    ! to one row in the array
    call write_fixed (path, "(*(f10.6))", dat2d_out, transform="none", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Writing file in non-transposed order"

    call read_fixed (path, "(*(f10.6))", dat2d_in, transform="none", &
        msg=iomsg, status=status)
    values_ok= all_close (dat2d_in, dat2d_out, atol=1.0e-4_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Reading written file in non-transposed order"

    ! Write in flatten (column-major) order
    call write_fixed (path, "(*(f10.6))", dat2d_out, transform="format", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "Writing file in flatten form"

    call read_fixed (path, "(*(f10.6))", dat2d_in, transform="format", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dat2d_out, dat2d_in, atol = &
        1.0e-5_PREC))) &
        error stop "Reading written file in flatten form"

    ! Test with reshaping via format specifier
    fmt = "(*(2(f10.8),:,/))"
    call write_fixed (path, fmt, dat2d_out, transform="format", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Writing file in flat form, reshaping via FMT"

    call read_fixed (path, fmt, dat2d_in, transform="format", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dat2d_out, dat2d_in, atol = &
        1.0e-5_PREC))) &
        error stop "Reading written file in flat form, reshaping via FMT"

    ! === Test with 3d array ===
    m = 5
    n = 4
    k = 2
    allocate (dat3d_out(m,n,k), dat3d_in(m,n,k))

    call random_number (dat3d_out)

    ! --- Full format string ---
    fmt = "(*(4(f15.8),:,/))"
    call write_fixed (path, fmt, dat3d_out, transform="format", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Writing 3d data in flat form, reshaping via FMT"

    call read_fixed (path, fmt, dat3d_in, transform="format", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dat3d_out, dat3d_in, atol = &
        1.0e-7_PREC))) &
        error stop "Reading 3d data in flat form, reshaping via FMT"

    ! --- Element-specific format string ---
    fmt = "f15.9"
    call write_fixed (path, fmt, dat3d_out, transform="format", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Writing 3d data in flat form, element-only FMT argument"

    fmt = "(5(f15.9))"
    call read_fixed (path, fmt, dat3d_in, transform="format", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dat3d_out, dat3d_in, atol = &
        1.0e-7_PREC))) &
        error stop "Reading 3d data in flat form, element-only FMT argument"

    ! --- Default formatting ---
    call write_fixed (path, dat=dat3d_out, transform="format", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "Writing 3d data in flat form, missing FMT argument"

    ! Use default formatting that is chosen in WRITE_FIXED
    fmt = "(5(g16.8))"
    call read_fixed (path, fmt, dat3d_in, transform="format", &
        msg=iomsg, status=status)
    if (.not. (status == NF_STATUS_OK .and. all_close (dat3d_out, dat3d_in, atol = &
        1.0e-7_PREC))) &
        error stop "Reading 3d data in flat form, element-only FMT argument"

end subroutine

subroutine __APPEND(test_fixed_alloc,__PREC) ()
    integer, parameter :: PREC = __PREC

    real (PREC), dimension(:,:), allocatable :: dat, dat_in
    integer :: nrow, ncol, i
    type (status_t) :: status
    logical :: values_ok
    character (20) :: fmt
    character (1024) :: tmpdir, path
    character (*), parameter :: filename = "test_numfort_io_fixed_alloc.txt"

    if (PREC == real32) then
    else if (PREC == real64) then
    end if

    call set_seed (2345)

    call get_temp_directory (tmpdir, status)
    path = trim(tmpdir) // "/" // filename

    ! === Test with 1 row / 1 col ===
    ncol = 1
    nrow = 1
    allocate (dat(nrow, ncol))

    call random_number (dat)
    status = NF_STATUS_UNDEFINED
    fmt = "(f12.6)"
    call write_fixed (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "Writing [1,1] array"

    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat_in, status=status)
    values_ok = all_close (dat, dat_in, atol=1.0e-6_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) error stop "Reading [1,1] array"

    ! Read with already allocated array
    status = NF_STATUS_UNDEFINED
    dat_in(:,:) = 0.0
    call read_fixed_alloc (path, fmt, dat_in, status=status)
    values_ok = all_close (dat, dat_in, atol=1.0e-6_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Reading [1,1] array, allocated DAT argument"

    deallocate (dat, dat_in)

    ! === Test with 1 row, multiple columns ===
    nrow = 1
    ncol = 7
    allocate (dat(nrow, ncol))

    call random_number (dat)
    dat(:,:) = dat * 1.0e4_PREC
    status = NF_STATUS_UNDEFINED
    fmt = "(f16.8)"
    call write_fixed (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "Writing [1,n] array"

    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat_in, status=status)
    values_ok = all_close (dat, dat_in, atol=1.0e-6_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) error stop "Reading [1,n] array"

    ! Read with already allocated DAT array
    status = NF_STATUS_UNDEFINED
    dat_in(:,:) = 0.0
    call read_fixed_alloc (path, fmt, dat_in, status=status)
    values_ok = all_close (dat, dat_in, atol=1.0e-6_PREC, rtol=0.0_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Reading [1,n] array, allocated DAT argument"

    deallocate (dat, dat_in)

    ! === Test with n rows, 1 column ===
    nrow = 1234
    ncol = 1
    allocate (dat(nrow, ncol))

    call random_number (dat)
    dat(:,:) = dat * 1.0e3_PREC
    status = NF_STATUS_UNDEFINED
    fmt = "(es13.6e2)"
    call write_fixed (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_OK)) error stop "Writing [n,1] array"

    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat_in, status=status)
    values_ok = all_close (dat, dat_in, atol=0.0_PREC, rtol=1.0e-4_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) error stop "Reading [n,1] array"

    ! Read with already allocated array
    status = NF_STATUS_UNDEFINED
    dat_in(:,:) = 0.0
    call read_fixed_alloc (path, fmt, dat_in, status=status)
    values_ok = all_close (dat, dat_in, atol=0.0_PREC, rtol=1.0e-4_PREC)
    if (.not. (status == NF_STATUS_OK .and. values_ok)) &
        error stop "Reading [n,1] array, allocated DAT argument"
    deallocate (dat, dat_in)

    ! === Test with n-by-n arrays of various sizes ===
    nrow = 6
    ncol = 3
    do i = 1, 7
        nrow = nrow * 2
        ncol = ncol * 2
        allocate (dat(nrow, ncol))

        call random_number (dat)
        dat(:,:) = dat * 1.0e2_PREC
        fmt = "(f12.6)"
        status = NF_STATUS_UNDEFINED
        call write_fixed (path, fmt, dat, status=status)
        if (.not. (status == NF_STATUS_OK)) &
            error stop "WRITE_FIXED failed the two-dimensional write check"

        status = NF_STATUS_UNDEFINED
        call read_fixed_alloc (path, fmt, dat_in, status=status)
        values_ok = all_close (dat_in, dat, atol=1.0e-3_PREC, rtol=0.0_PREC)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "READ_FIXED_ALLOC failed the two-dimensional read check"

        ! Read using already allocated array
        status = NF_STATUS_UNDEFINED
        dat_in(:,:) = 0.0
        call read_fixed_alloc (path, fmt, dat_in, status=status)
        values_ok = all_close (dat_in, dat, atol=1.0e-3_PREC, rtol=0.0_PREC)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "READ_FIXED_ALLOC failed the preallocated array read check"

        deallocate (dat_in)

        ! Read in as transpose
        status = NF_STATUS_UNDEFINED
        call read_fixed_alloc (path, fmt, dat_in, transform="transpose", status=status)
        values_ok = all_close (transpose(dat_in), dat, atol=1.0e-3_PREC, rtol=0.0_PREC)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "READ_FIXED_ALLOC failed the transposed array read check"

        ! Read in as transpose, array allocated with wrong shape
        status = NF_STATUS_UNDEFINED
        deallocate (dat_in)
        allocate (dat_in(ncol+1, nrow+1), source=0.0_PREC)
        call read_fixed_alloc (path, fmt, dat_in, transform="transpose", status=status)
        values_ok = all_close (transpose(dat_in), dat, atol=1.0e-3_PREC, rtol=0.0_PREC)
        if (.not. (status == NF_STATUS_OK .and. values_ok)) &
            error stop "READ_FIXED_ALLOC failed the resized transposed array read " &
                // "check"

        deallocate (dat, dat_in)
    end do

end subroutine
