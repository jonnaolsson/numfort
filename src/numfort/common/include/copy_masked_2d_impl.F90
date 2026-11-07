
! Implementation for COPY_MASKED_2D

real (PREC), intent(in), dimension(:,:), contiguous :: src
    !*  Input array
real (PREC), intent(inout), dimension(:,:), contiguous :: dst
    !*  Output array. Dimension along which subset of elements to be selectively
    !   copied must be at least as large the corresponding in dimension in
    !   SRC, while the remaining dimension must have identical size to SRC.
logical, intent(in), dimension(:), contiguous :: mask
    !*  Logical array identifying elements to be copied. Must have the
    !   same length as the dimension in SRC along which elements should
    !   be selectively copied.
integer, intent(in), optional :: dim
    !*  Dimension along which to operate (default: 1)
type (status_t), intent(out), optional :: status
    !*  Optional exit code.

type (status_t) :: lstatus
integer :: ldim, m, n, Nvalid, shp(2), i, j
integer, dimension(:), allocatable :: idims

lstatus = NF_STATUS_OK

! === Input validation ===
! Guess DIM argument if none provided
if (.not. present(dim)) then
    if (size(src, 1) == size(mask)) then
        ldim = 1
    else if (size(src, 2) == size(mask)) then
        ldim = 2
    else
        ! Default: use first dimension; check whether this conforms with
        ! shape of SRC below.
        ldim = 1
    end if
else
    ldim = dim
end if

if (ldim < 1 .or. ldim > 2) then
    lstatus = NF_STATUS_INVALID_ARG
    goto 100
end if

if (size(src, ldim) /= size(mask)) then
    lstatus = NF_STATUS_INVALID_ARG
    goto 100
end if

Nvalid = count(mask)

shp = shape(src)
shp(ldim) = Nvalid

! Check that destination array is large enough to hold the result.
! Require that the size of the dimension which is NOT to be filtered
! is the same between SRC and DST.
if (ldim == 1) then
    if (size(src, 2) /= size(dst, 2) .or. size(dst,1) < Nvalid) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if
else
    if (size(src, 1) /= size(dst, 1) .or. size(dst, 2) < Nvalid) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if
end if

if (Nvalid == 0 .or. size (src) == 0) goto 100

if (Nvalid == size(mask) .and. shape_equal (src, dst)) then
    ! Shortcut: if everything is to be included and arrays are of equal size,
    ! copy over data directly.
    dst(:,:) = src
    goto 100
end if

m = size(src, 1)
n = size(src, 2)

if (ldim == 1) then
    ! Compute array containing the valid indices along dimension 1,
    ! we will use this index array to copy all valid values in each
    ! row at once.
    allocate (idims(NVALID))
    j = 0
    do i = 1, m
        if (.not. mask(i)) cycle
        j = j + 1
        idims(j) = i
    end do

    do j = 1, n
        dst(1:Nvalid,j) = src(idims, j)
    end do
else
    ! Copy by column, skipping those that are marked as invalid
    i = 0
    do j = 1, n
        if (.not. mask(j)) cycle
        i = i + 1
        dst(:,i) = src(:,j)
    end do
end if


100 continue

if (allocated(idims)) deallocate (idims)

if (present(status)) status = lstatus
