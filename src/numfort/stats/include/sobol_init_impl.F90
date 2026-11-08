

pure subroutine __APPEND(sobol_init_check_input,__INTSIZE) (ndim, s, a, m, max_len, status)

    integer, parameter :: INTSIZE = __INTSIZE

    integer (INTSIZE), intent(in) :: ndim
        !*  Number of dimension of each sequence element
    integer (INTSIZE), intent(in), dimension(:), optional :: s
        !*  Array of primitive polynomial degrees, one for each dimension
    integer (INTSIZE), intent(in), dimension(:), optional :: a
        !*  Array of polynomial coefficients stored (either 0 or 1) stored
        !   as integers, one for each dimension.
    integer (INTSIZE), intent(in), dimension(:,:), optional :: m
        !*  Array of direction numbers. Column j should contains 1:s(j) direction
        !   numbers for dimension j.
    integer (INTSIZE), intent(in), optional :: max_len
        !*  Maximum sequence length (max. 2 ** 63 points, which is also the
        !   default value)
    type (status_t), intent(out) :: status

    logical :: user_data
    integer :: lndim, nm

    integer (int64), parameter :: MAX_SUPPORTED_LEN = &
        ishft(1_int64, SOBOL_MAX_BITS - 1_int64)

    status = NF_STATUS_INVALID_ARG

    user_data = present(s) .and. present(a) .and. present(m)

    if (user_data) then
        lndim = size(s)
        if (lndim /= size(a) .or. lndim > size(m, 2)) goto 100
        nm = maxval(s)
        if (nm > size(m, 1)) goto 100

        if (lndim /= ndim) goto 100
    else
        ! Do not permit partial specification of user-provided direction numbers
        if (present(s) .or. present(a) .or. present(m)) goto 100
        lndim = ndim
    end if

    if (lndim < 1) goto 100

    if (present(max_len)) then
        if (max_len < 1 .or. max_len > MAX_SUPPORTED_LEN) goto 100
    end if

    status = NF_STATUS_OK

100 continue

end subroutine


pure subroutine __APPEND(sobol_init_array,__INTSIZE) (self, ndim, s, a, m, max_len, status)
    !*  SOBOL_INIT initializes the state used to sample from a Sobol sequence
    !   generator using a recursive algorithm.

    integer, parameter :: INTSIZE = __INTSIZE

    type (sobol_state), intent(in out) :: self
        !*  Object to store Sobol sequence state
    integer (INTSIZE), intent(in) :: ndim
        !*  Number of dimension of each sequence element
    integer (INTSIZE), intent(in), dimension(:), optional :: s
        !*  Array of primitive polynomial degrees, one for each dimension
    integer (INTSIZE), intent(in), dimension(:), optional :: a
        !*  Array of polynomial coefficients stored (either 0 or 1) stored
        !   as integers, one for each dimension.
    integer (INTSIZE), intent(in), dimension(:,:), optional :: m
        !*  Array of direction numbers. Column j should contains 1:s(j) direction
        !   numbers for dimension j.
    integer (INTSIZE), intent(in), optional :: max_len
        !*  Maximum sequence length (max. 2 ** 63 points, which is also the
        !   default value)
    type (status_t), intent(out), optional :: status
        !*  Status object. If present, contains status information on exit.

    type (status_t) :: lstatus
    integer (INTSIZE) :: lndim

    integer :: max_bits
    integer (int64), dimension(SOBOL_MAX_BITS) :: m_j
        !   Temporary array to store (m_1,m_2,...) for dimension j
    integer (int64) :: s_j, a_j, mk
    integer :: i, j, k, offset
    logical :: user_data

    call sobol_init_check_input (ndim, s, a, m, max_len, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call sobol_reset (self)

    user_data = present(s) .and. present(a) .and. present(m)
    lndim = ndim
    self%ndim = lndim

    ! Maximum number of bits needed, depending on max. lenth
    if (present(max_len)) then
        !
        max_bits = ceiling(log(real(max_len, real64)) / log(2.0d0))
        if (max_bits >= SOBOL_MAX_BITS) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    else
        max_bits = SOBOL_MAX_BITS
    end if

    ! Store dimensions in rows since then evaluating Sobol sequences
    ! in recursive fashion uses column-wise memory access patterns.
    allocate (self%v(lndim, max_bits))

    ! compute m_k for each dimension for sequences of max. length given
    ! by 2 ** max_bits.
    ! Dimension 1: Let all m_k = 1
    do i = 1, max_bits
        self%v(1,i) = ishft(1_int64, int(SOBOL_MAX_BITS - i, int64))
    end do

    do j = 2, lndim
        if (user_data) then
            s_j = s(j)
            m_j(1:s_j) = m(1:s_j,j)
            a_j = a(j)
        else
            s_j = SOBOL_DEFAULT_S(j)
            offset = sum(SOBOL_DEFAULT_S(:j-1))
            m_j(1:s_j) = SOBOL_DEFAULT_M(offset+1:offset+s_j)
            a_j = SOBOL_DEFAULT_A(j)
        end if

        ! scale up by 2 ** 64
        do i = 1, s_j
            m_j(i) = ishft(m_j(i), int(SOBOL_MAX_BITS - i, int64))
        end do

        ! Compute "direction numbers" m_{k,j} for k > s_j
        do k = s_j + 1, max_bits
            mk = ieor(m_j(k-s_j), ishft(m_j(k-s_j), - s_j))
            do i = 1, s_j - 1
                ! Negative shift extracts bit (s_j - i - 1) of a_j; a positive shift
                ! in Fortran shifts left, which drops higher-degree coefficients.
                mk = ieor(mk, m_j(k-i) * iand(ishft(a_j, -(s_j - i - 1)), 1_int64))
            end do
            m_j(k) = mk
        end do

        self%v(j,:) = m_j
    end do

100 continue
    if (present(status)) status = lstatus

end subroutine

pure subroutine __APPEND(sobol_init_scalar,__INTSIZE) (self, ndim, s, a, m, max_len, status)

    integer, parameter :: INTSIZE = __INTSIZE

    type (sobol_state), intent(in out) :: self
    integer (INTSIZE), intent(in), optional :: ndim
    integer (INTSIZE), intent(in) :: s
    integer (INTSIZE), intent(in) :: a
    integer (INTSIZE), intent(in), dimension(:) :: m
    integer (INTSIZE), intent(in), optional :: max_len
    type (status_t), intent(out), optional :: status

    integer (INTSIZE), dimension(1) :: s1, a1
    integer (INTSIZE), dimension(:,:), allocatable :: m2
    integer (INTSIZE), parameter :: lndim = 1

    allocate (m2(size(m), 1))
    s1(1) = s
    a1(1) = a
    call sobol_init (self, lndim, s1, a1, m2, max_len, status)

end subroutine
