
    ! Implementation of READ_FIXED_ALLOC.
    ! Infers record length from the first record and requires all subsequent
    ! records to have the same length.

    character (*), intent(in) :: path
    character (*), intent(in) :: fmt
    character (*), intent(in), optional :: transform
    type (status_t), intent(out), optional :: status
    character (*), intent(out), optional :: msg


    integer, parameter :: LINE_BUFFER_LEN = 2**14
    integer :: linelen, fwidth, nrec, nlines, i, n, ioffset
    character (:), allocatable :: buf
    integer :: shp(2)
    integer, parameter :: CHUNK_NCOL = 100
    character (:), allocatable :: lfmt, fmt_field

    integer (NF_ENUM_KIND) :: itransform
    integer :: iostat, uid
    character (100) :: lmsg
    type (status_t) :: lstatus

    nullify (ptr_first, ptr_curr, ptr_next)

    ! Default exit status
    lstatus = NF_STATUS_IO_ERROR
    lmsg = ''

    ! ----- Parse and check TRANSFORM argument -----
    itransform = parse_transform (transform)
    if (itransform == 0) then
        lstatus = NF_STATUS_INVALID_ARG
        lmsg = "Invalid value for argument 'transform'"
        goto 100
    end if

    if (itransform /= NF_IO_TRANSFORM_NONE .and. &
            itransform /= NF_IO_TRANSFORM_TRANSPOSE) then
        lstatus = NF_STATUS_UNSUPPORTED_OP
        lmsg = "Unsupported value for argument 'transform'"
        goto 100
    end if

    ! ----- Determine field width implied by format string -----
    fwidth = format_get_field_width (fmt, dummy)
    if (fwidth < 0) then
        lstatus = NF_STATUS_INVALID_ARG
        lmsg = 'Invalid format specifier'
        goto 100
    end if

    open (newunit=uid, file=path, status='old', action='read', &
        access='sequential', iostat=iostat, iomsg=lmsg, err=100)

    ! ---- Determine line length -----
    ! Read through the end of the first record in chunks of LINE_BUFFER_LEN
    ! to determine the complete record length without truncating long records.
    allocate (character (LINE_BUFFER_LEN) :: buf)
    linelen = 0
    do
        read (unit=uid, fmt="(a)", iomsg=lmsg, iostat=iostat, size=n, &
            advance="no", err=50) buf
        if (is_iostat_end (iostat)) exit
        linelen = linelen + n
        if (is_iostat_eor (iostat)) exit
    end do

    deallocate (buf)

    nrec = linelen / fwidth

    if (nrec <= 0 .or. nrec * fwidth /= linelen) then
        lstatus = NF_STATUS_INVALID_STATE
        lmsg = "Failed to determine record length"
        goto 50
    end if

    rewind (unit=uid, iostat=iostat, iomsg=lmsg, err=50)

    ! ----- Read data into chunks -----
    ! We still don't know the number of lines and hance the ultimate array
    ! size, so we allocate chunks of size (NREC,CHUNK_NCOL) and use a
    ! linked list to store all the data read in a dynamic fashion.

    allocate (ptr_first)
    ptr_curr => ptr_first

    allocate (ptr_curr%dat(nrec,CHUNK_NCOL))

    nlines = 0

    ! Add some additional space that will be used below
    n = len_trim (fmt)
    allocate (character (n+20) :: lfmt)
    allocate (character (n) :: fmt_field)
    call format_strip_parenthesis (fmt, fmt_field)
    ! Do not add any EOR character at the end, this breaks ifort
    write (lfmt, '("(*(", a, "))")') trim(fmt_field)

    ! Allocate a buffer of length linelen + 1 so that any record of exact
    ! length linelen ends with IOSTAT_EOR and n == linelen, shorter records
    ! yield n < linelen, and longer records fill linelen + 1 characters
    ! without reaching end-of-record.
    allocate (character (linelen + 1) :: buf)

    do while (.true.)
        do i = 1, CHUNK_NCOL
            read (unit=uid, fmt="(a)", iomsg=lmsg, iostat=iostat, size=n, &
                advance="no", err=10) buf
            if (is_iostat_end (iostat)) goto 5
            if (.not. is_iostat_eor (iostat) .or. n /= linelen) then
                lstatus = NF_STATUS_INVALID_STATE
                lmsg = "Inconsistent record length"
                goto 10
            end if
            read (buf(1:linelen), fmt=lfmt, iomsg=lmsg, iostat=iostat, err=10) &
                ptr_curr%dat(:,i)
            nlines = nlines + 1
        end do

        ! Haven't reached the end of the file yet, need to allocate another
        ! data chunk
        allocate (ptr_curr%ptr_next)
        ptr_curr => ptr_curr%ptr_next
        allocate (ptr_curr%dat(nrec, CHUNK_NCOL))
    end do

5   continue

    if (nlines <= 0) then
        lstatus = NF_STATUS_INVALID_STATE
        lmsg = "Failed to determine record length"
        goto 10
    end if

    ! Processing reached EOF, no error encountered
    ! Copy linked list of chunks into output array, transform as requested

    select case (itransform)
    case (NF_IO_TRANSFORM_NONE)
        ! Create output data such that one column in the data array corresponds
        ! to one row in the data file. This effectively orders the data in
        ! the same way as it is ordered in the data file.
        shp(1) = nrec
        shp(2) = nlines
        call cond_alloc (dat, shp)

        ptr_curr => ptr_first
        ioffset = 0
        do while (associated(ptr_curr))
            n = min(CHUNK_NCOL, nlines - ioffset)
            do i = 1, n
                dat(:,ioffset+i) = ptr_curr%dat(:,i)
            end do
            ioffset = ioffset + n
            ptr_curr => ptr_curr%ptr_next
        end do

    case (NF_IO_TRANSFORM_TRANSPOSE)
        ! Create output data such that one row (line) in the data file
        ! corresponds to one row in the input array. Note that for Fortran this
        ! effectively transposes the data (data is left as-is in C order)
        shp(1) = nlines
        shp(2) = nrec
        call cond_alloc (dat, shp)

        ptr_curr => ptr_first
        ioffset = 0
        do while (associated(ptr_curr))
            n = min(CHUNK_NCOL, nlines - ioffset)
            do i = 1, n
                dat(ioffset+i,:) = ptr_curr%dat(:,i)
            end do
            ioffset = ioffset + n
            ptr_curr => ptr_curr%ptr_next
        end do

    end select

    lstatus = NF_STATUS_OK
    ! Clear local diagnostic so transient runtime EOF messages from the
    ! read loop above are not copied into msg on success.
    lmsg = ""

10  continue

    ! Clean up dynamically allocated memory
    ptr_curr => ptr_first

    do while (associated(ptr_curr))
        if (allocated(ptr_curr%dat)) deallocate (ptr_curr%dat)
        ptr_next => ptr_curr%ptr_next
        ! Deallocate container object itself
        deallocate (ptr_curr)
        ptr_curr => ptr_next
    end do

    nullify (ptr_curr, ptr_first, ptr_next)


50  continue

    ! Close file
    close (unit=uid)

100 continue

    if (allocated(buf)) deallocate (buf)
    if (allocated(fmt_field)) deallocate (fmt_field)
    if (allocated(lfmt)) deallocate (lfmt)

    if (present(status)) status = lstatus
    if (present(msg)) msg = lmsg
