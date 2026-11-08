

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


subroutine __APPEND(test_fixed_alloc_malformed,__PREC) ()
    !*  Test error handling in READ_FIXED_ALLOC for malformed first records:
    !   - Line length not divisible by field width.
    !   - Empty newline record and zero-byte empty file (EOF on initial probe).
    !   - Absence of MSG optional argument.
    !   - Retention of preallocated array state on failure.

    integer, parameter :: PREC = __PREC

    real (PREC), dimension(:,:), allocatable :: dat
    type (status_t) :: status
    character (100) :: msg
    character (20) :: fmt
    character (1024) :: tmpdir, path
    character (*), parameter :: filename = "test_numfort_io_fixed_alloc_malformed.txt"
    integer :: uid, iostat

    call get_temp_directory (tmpdir, status)
    path = trim(tmpdir) // "/" // filename

    ! Hand-written fixture with a first line of exactly 5 characters: '12345'.
    ! The field width for '(f4.1)' is 4.
    ! Because 5 mod 4 = 1 /= 0, the record length is not divisible by the field width.
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open temporary fixture file"
    write (uid, "(a)") "12345"
    close (unit=uid)

    fmt = "(f4.1)"

    ! Test case 1: Malformed line width with unallocated DAT and MSG present
    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Status not INVALID_STATE for malformed width with MSG"
    if (trim(msg) /= "Failed to determine record length") &
        error stop "Diagnostic message mismatch for malformed width"
    if (allocated(dat)) &
        error stop "Output array unexpectedly allocated on failure"

    ! Test case 2: Malformed line width with unallocated DAT and MSG omitted
    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Status not INVALID_STATE for malformed width without MSG"
    if (allocated(dat)) &
        error stop "Output array unexpectedly allocated on failure without MSG"

    ! Test case 3: Malformed line width with preallocated DAT and MSG present
    allocate (dat(2, 2))
    dat(:,:) = -999.0_PREC
    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Status not INVALID_STATE with preallocated DAT and MSG"
    if (trim(msg) /= "Failed to determine record length") &
        error stop "Diagnostic message mismatch with preallocated DAT"
    if (.not. (allocated(dat) .and. size(dat, 1) == 2 .and. size(dat, 2) == 2)) &
        error stop "Preallocated array shape altered on failure"
    if (any(dat /= -999.0_PREC)) &
        error stop "Preallocated array contents altered on failure"

    ! Test case 4: Malformed line width with preallocated DAT and MSG omitted
    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Status not INVALID_STATE with preallocated DAT without MSG"
    if (any(dat /= -999.0_PREC)) &
        error stop "Preallocated array contents altered without MSG"

    deallocate (dat)

    ! Hand-written empty record fixture (0 characters).
    ! Field width is 4; line length is 0, so nrec = 0 <= 0 cannot determine record length.
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open temporary fixture file for empty line"
    write (uid, "(a)") ""
    close (unit=uid)

    ! Test case 5: Empty line record with MSG present
    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Status not INVALID_STATE for empty record"
    if (trim(msg) /= "Failed to determine record length") &
        error stop "Diagnostic message mismatch for empty record"
    if (allocated(dat)) &
        error stop "Output array unexpectedly allocated for empty record"

    ! Test case 6: Empty line record without MSG
    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Status not INVALID_STATE for empty record without MSG"

    ! Zero-byte empty file fixture (0 bytes, immediate IOSTAT_END on first probe read).
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open temporary fixture for zero-byte file"
    close (unit=uid)

    ! Test case 7: Zero-byte empty file with MSG present
    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Status not INVALID_STATE for zero-byte file"
    if (trim (msg) /= "Failed to determine record length") &
        error stop "Diagnostic message mismatch for zero-byte file"
    if (allocated (dat)) &
        error stop "Output array unexpectedly allocated for zero-byte file"

    ! Clean up temporary file
    open (newunit=uid, file=path, status="old", action="read", iostat=iostat)
    if (iostat == 0) close (unit=uid, status="delete")

end subroutine


subroutine __APPEND(test_fixed_alloc_long_records,__PREC) ()
    !*  Test READ_FIXED_ALLOC with long records and chunk boundaries:
    !   - Records at the 16 KiB probe buffer boundary (16384 characters).
    !   - Records above the boundary (16400 characters, asserting final fields).
    !   - Records below the boundary (16380 characters).
    !   - Multi-row records exceeding the boundary.
    !   - Malformed first-record lengths exceeding the boundary (16401 characters).
    !   - Chunk-boundary crossing (105 rows exceeding CHUNK_NCOL = 100).
    !   - Multi-chunk probe accumulation (> 32 KiB, 32784 characters).

    integer, parameter :: PREC = __PREC

    real (PREC), dimension(:,:), allocatable :: dat
    type (status_t) :: status
    character (100) :: msg
    character (20) :: fmt
    character (4) :: tail_val
    character (1024) :: tmpdir, path
    character (*), parameter :: filename = "test_numfort_io_fixed_alloc_long.txt"
    integer :: uid, iostat, k

    call get_temp_directory (tmpdir, status)
    path = trim (tmpdir) // "/" // filename

    fmt = "(f4.1)"

    ! Test case 1: Long record exactly at the 16 KiB boundary (16384 characters).
    ! 4096 fields * 4 characters per field = 16384 characters = LINE_BUFFER_LEN.
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for boundary test"
    write (uid, "(a)") repeat (" 1.0", 4095) // " 2.5"
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "READ_FIXED_ALLOC failed at boundary length"
    if (.not. (allocated (dat) .and. size (dat, 1) == 4096 &
            .and. size (dat, 2) == 1)) &
        error stop "Shape mismatch for boundary record length"
    if (.not. (abs (dat(1, 1) - 1.0_PREC) < 1.0e-5_PREC &
            .and. abs (dat(4096, 1) - 2.5_PREC) < 1.0e-5_PREC)) &
        error stop "Value mismatch for boundary record length"
    deallocate (dat)

    ! Test case 2: Long record above boundary (16400 characters, 4100 fields).
    ! 4100 fields * 4 characters = 16400 characters (exceeds 16384 by 16).
    ! Fields 4097-4100 are unique to ensure truncation beyond 16 KiB is detected.
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for above-boundary test"
    write (uid, "(a)") repeat (" 1.0", 4096) // " 2.1 3.2 4.3 5.4"
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "READ_FIXED_ALLOC failed above boundary length"
    if (.not. (allocated (dat) .and. size (dat, 1) == 4100 &
            .and. size (dat, 2) == 1)) &
        error stop "Shape mismatch above boundary length"
    if (.not. (abs (dat(4097, 1) - 2.1_PREC) < 1.0e-5_PREC &
            .and. abs (dat(4098, 1) - 3.2_PREC) < 1.0e-5_PREC &
            .and. abs (dat(4099, 1) - 4.3_PREC) < 1.0e-5_PREC &
            .and. abs (dat(4100, 1) - 5.4_PREC) < 1.0e-5_PREC)) &
        error stop "Final field values truncated or corrupted above boundary"
    deallocate (dat)

    ! Test case 3: Long record below boundary (16380 characters, 4095 fields).
    ! 4095 fields * 4 characters = 16380 characters (4 characters below boundary).
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for below-boundary test"
    write (uid, "(a)") repeat (" 1.0", 4094) // " 9.9"
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "READ_FIXED_ALLOC failed below boundary length"
    if (.not. (allocated (dat) .and. size (dat, 1) == 4095 &
            .and. size (dat, 2) == 1)) &
        error stop "Shape mismatch below boundary length"
    if (abs (dat(4095, 1) - 9.9_PREC) >= 1.0e-5_PREC) &
        error stop "Value mismatch below boundary length"
    deallocate (dat)

    ! Test case 4: Multi-row records above boundary (3 rows * 4100 fields).
    ! Each row is 16400 characters; final field of row k has unique value k * 1.1.
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for multi-row test"
    do k = 1, 3
        write (tail_val, "(f4.1)") real (k, PREC) * 1.1_PREC
        write (uid, "(a)") repeat (" 0.0", 4099) // tail_val
    end do
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "READ_FIXED_ALLOC failed multi-row long records"
    if (.not. (allocated (dat) .and. size (dat, 1) == 4100 &
            .and. size (dat, 2) == 3)) &
        error stop "Shape mismatch multi-row long records"
    if (.not. (abs (dat(4100, 1) - 1.1_PREC) < 1.0e-5_PREC &
            .and. abs (dat(4100, 2) - 2.2_PREC) < 1.0e-5_PREC &
            .and. abs (dat(4100, 3) - 3.3_PREC) < 1.0e-5_PREC)) &
        error stop "Value mismatch in multi-row long records final fields"
    deallocate (dat)

    ! Test case 5: Malformed long record above boundary (16401 characters).
    ! 4100 fields * 4 characters + 1 extra character = 16401 characters.
    ! Since 16401 mod 4 = 1 /= 0, line length is not divisible by field width.
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for malformed long record"
    write (uid, "(a)") repeat (" 1.0", 4100) // "7"
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Status not INVALID_STATE for malformed record above boundary"
    if (trim (msg) /= "Failed to determine record length") &
        error stop "Diagnostic message mismatch for malformed record above boundary"
    if (allocated (dat)) &
        error stop "Array allocated on malformed record above boundary"

    ! Test case 6: Chunk-boundary crossing (105 rows, exceeding CHUNK_NCOL = 100).
    ! Tests linking across chunk 1 (rows 1-100) and chunk 2 (rows 101-105).
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for chunk boundary test"
    do k = 1, 105
        ! For k = 100, mod(100, 10) * 0.5 = 0.0
        ! For k = 101, mod(101, 10) * 0.5 = 0.5
        ! For k = 105, mod(105, 10) * 0.5 = 2.5
        write (uid, "(f4.1,f4.1)") real (mod (k, 10), PREC) * 0.5_PREC, 0.9_PREC
    end do
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "READ_FIXED_ALLOC failed chunk boundary test"
    if (.not. (allocated (dat) .and. size (dat, 1) == 2 &
            .and. size (dat, 2) == 105)) &
        error stop "Shape mismatch chunk boundary test"
    if (.not. (abs (dat(1, 100) - 0.0_PREC) < 1.0e-5_PREC &
            .and. abs (dat(1, 101) - 0.5_PREC) < 1.0e-5_PREC &
            .and. abs (dat(1, 105) - 2.5_PREC) < 1.0e-5_PREC)) &
        error stop "Value mismatch across chunk boundary"
    deallocate (dat)

    ! Test case 7: Multi-chunk probe accumulation (> 32 KiB, 32784 characters).
    ! 8196 fields * 4 characters = 32784 = 2 * 16384 + 16 (3 probe iterations).
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for multi-chunk probe test"
    write (uid, "(a)") repeat (" 1.0", 8192) // " 6.1 7.2 8.3 9.4"
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_OK)) &
        error stop "READ_FIXED_ALLOC failed multi-chunk probe length"
    if (.not. (allocated (dat) .and. size (dat, 1) == 8196 &
            .and. size (dat, 2) == 1)) &
        error stop "Shape mismatch for multi-chunk probe length"
    if (.not. (abs (dat(8193, 1) - 6.1_PREC) < 1.0e-5_PREC &
            .and. abs (dat(8194, 1) - 7.2_PREC) < 1.0e-5_PREC &
            .and. abs (dat(8195, 1) - 8.3_PREC) < 1.0e-5_PREC &
            .and. abs (dat(8196, 1) - 9.4_PREC) < 1.0e-5_PREC)) &
        error stop "Final field values corrupted on multi-chunk probe length"
    deallocate (dat)

    ! Clean up temporary file
    open (newunit=uid, file=path, status="old", action="read", iostat=iostat)
    if (iostat == 0) close (unit=uid, status="delete")

end subroutine


subroutine __APPEND(test_fixed_alloc_row_widths,__PREC) ()
    !*  Test READ_FIXED_ALLOC width enforcement on subsequent records:
    !   - Short second row with fewer fields (rejects silent zero-padding).
    !   - Empty second row.
    !   - Partial field width on second row with omitted MSG argument.
    !   - Overlong second row with extra fields (rejects silent truncation).
    !   - Overlong second row with trailing character and preallocated DAT.
    !   - Truncated final row without trailing newline at EOF.
    !   - Width mismatch past CHUNK_NCOL = 100 chunk boundary.
    !   - Short and overlong subsequent rows exceeding the 16 KiB buffer.

    integer, parameter :: PREC = __PREC

    real (PREC), dimension(:,:), allocatable :: dat
    type (status_t) :: status
    character (100) :: msg
    character (20) :: fmt
    character (1024) :: tmpdir, path
    character (*), parameter :: filename = "test_numfort_io_fixed_alloc_widths.txt"
    integer :: uid, iostat, k

    call get_temp_directory (tmpdir, status)
    path = trim (tmpdir) // "/" // filename

    fmt = "(f4.1)"

    ! Test case 1: Short second row (row 1 has 3 fields = 12 chars; row 2 has 2 = 8).
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for short second row"
    write (uid, "(a)") " 1.0 2.0 3.0"
    write (uid, "(a)") " 4.0 5.0"
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Short second row did not return INVALID_STATE"
    if (trim (msg) /= "Inconsistent record length") &
        error stop "Diagnostic mismatch for short second row"
    if (allocated (dat)) &
        error stop "Output array allocated on short second row"

    ! Test case 2: Empty second row (row 1 has 8 chars; row 2 has 0 chars).
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for empty second row"
    write (uid, "(a)") " 1.0 2.0"
    write (uid, "(a)") ""
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Empty second row did not return INVALID_STATE"
    if (trim (msg) /= "Inconsistent record length") &
        error stop "Diagnostic mismatch for empty second row"
    if (allocated (dat)) &
        error stop "Output array allocated on empty second row"

    ! Test case 3: Partial field on second row (7 chars vs 8) with MSG omitted.
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for partial second row"
    write (uid, "(a)") " 1.0 2.0"
    write (uid, "(a)") " 3.0 4."
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    call read_fixed_alloc (path, fmt, dat, status=status)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Partial second row without MSG did not return INVALID_STATE"
    if (allocated (dat)) &
        error stop "Output array allocated on partial second row"

    ! Test case 4: Extra field on second row (row 1 has 8 chars; row 2 has 12 chars).
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for overlong second row"
    write (uid, "(a)") " 1.0 2.0"
    write (uid, "(a)") " 3.0 4.0 5.0"
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Extra-field second row did not return INVALID_STATE"
    if (trim (msg) /= "Inconsistent record length") &
        error stop "Diagnostic mismatch for extra-field second row"
    if (allocated (dat)) &
        error stop "Output array allocated on extra-field second row"

    ! Test case 5: Extra trailing character on row 2 (9 chars vs 8) with preallocated
    ! DAT; verifies preallocated array shape and values remain untouched on failure.
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for trailing-char row"
    write (uid, "(a)") " 1.0 2.0"
    write (uid, "(a)") " 3.0 4.07"
    close (unit=uid)

    allocate (dat(2, 2))
    dat(:,:) = -999.0_PREC
    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Trailing-char second row did not return INVALID_STATE"
    if (trim (msg) /= "Inconsistent record length") &
        error stop "Diagnostic mismatch for trailing-char second row"
    if (.not. (allocated (dat) .and. size (dat, 1) == 2 &
            .and. size (dat, 2) == 2)) &
        error stop "Preallocated array shape altered on row-width failure"
    if (any (dat /= -999.0_PREC)) &
        error stop "Preallocated array contents altered on row-width failure"
    deallocate (dat)

    ! Test case 6: Short final row without trailing newline at EOF (4 chars vs 8).
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for unterminated short row"
    write (uid, "(a)") " 1.0 2.0"
    write (uid, "(a)", advance="no") " 3.0"
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Unterminated short final row did not return INVALID_STATE"
    if (trim (msg) /= "Inconsistent record length") &
        error stop "Diagnostic mismatch for unterminated short final row"
    if (allocated (dat)) &
        error stop "Output array allocated on unterminated short final row"

    ! Test case 7: Width mismatch on row 102 (after CHUNK_NCOL = 100 chunk boundary).
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for second-chunk mismatch"
    do k = 1, 101
        write (uid, "(a)") " 1.0 2.0"
    end do
    write (uid, "(a)") " 3.0"
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Second-chunk short row did not return INVALID_STATE"
    if (trim (msg) /= "Inconsistent record length") &
        error stop "Diagnostic mismatch on second-chunk short row"
    if (allocated (dat)) &
        error stop "Output array allocated on second-chunk short row"

    ! Test case 8: Short second row above 16 KiB boundary (16400 vs 16396 chars).
    ! Row 1 has 4100 * 4 = 16400 characters; row 2 has 4099 * 4 = 16396 characters.
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for long short-row test"
    write (uid, "(a)") repeat (" 1.0", 4100)
    write (uid, "(a)") repeat (" 1.0", 4099)
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Long short second row did not return INVALID_STATE"
    if (trim (msg) /= "Inconsistent record length") &
        error stop "Diagnostic mismatch for long short second row"
    if (allocated (dat)) &
        error stop "Output array allocated on long short second row"

    ! Test case 9: Overlong second row above 16 KiB boundary (16400 vs 16404 chars).
    ! Row 1 has 4100 * 4 = 16400 characters; row 2 has 4101 * 4 = 16404 characters.
    open (newunit=uid, file=path, status="replace", action="write", iostat=iostat)
    if (iostat /= 0) error stop "Failed to open fixture for long overlong-row test"
    write (uid, "(a)") repeat (" 1.0", 4100)
    write (uid, "(a)") repeat (" 1.0", 4101)
    close (unit=uid)

    status = NF_STATUS_UNDEFINED
    msg = ""
    call read_fixed_alloc (path, fmt, dat, status=status, msg=msg)
    if (.not. (status == NF_STATUS_INVALID_STATE)) &
        error stop "Long overlong second row did not return INVALID_STATE"
    if (trim (msg) /= "Inconsistent record length") &
        error stop "Diagnostic mismatch for long overlong second row"
    if (allocated (dat)) &
        error stop "Output array allocated on long overlong second row"

    ! Clean up temporary file
    open (newunit=uid, file=path, status="old", action="read", iostat=iostat)
    if (iostat == 0) close (unit=uid, status="delete")

end subroutine
