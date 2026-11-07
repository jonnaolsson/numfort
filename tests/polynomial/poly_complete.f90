

program test_polynomial_poly_complete
    !*  Unit tests for complete polynomial routines.

    use, intrinsic :: iso_fortran_env

    use numfort_arrays
    use numfort_common
    use numfort_common_testing
    use numfort_core
    use numfort_polynomial

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

contains

subroutine test_all ()

    call test_exponents ()
    call test_basis ()
    call test_basis_jac ()

end subroutine

subroutine test_exponents ()
    !*  Unit tests for creating a list of all permutations of exponents
    !   that occur in a complete polynomial of given degree in a given number
    !   of variables.

    integer :: kmax, k, n, j, ifrom, ito, nterms
    logical :: ok
    integer, dimension(:,:), allocatable :: exponents, exponents_ok
    type (status_t) :: status

    ! Test with invalid arguments
    allocate (exponents(0,0))

    status = NF_STATUS_OK
    call polyexponents_complete (-1, 0, exponents, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid dimension argument"

    status = NF_STATUS_OK
    call polyexponents_complete (0, -1, exponents, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid polynomial degree"

    status = NF_STATUS_OK
    call polyexponents_complete (1, 0, exponents, status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "Invalid EXPONENT array size"

    deallocate (exponents)

    ! Test with 1 dimension
    n = 1
    kmax = 5
    allocate (exponents(n,kmax+1))

    do k = 0, kmax
        status = NF_STATUS_OK
        call polyexponents_complete (n, k, exponents, status)
        ok = .true.
        do j = 1, k+1
            ok = ok .and. (exponents(1,j) == (j-1))
        end do
        if (.not. (ok .and. status == NF_STATUS_OK)) &
            error stop "Complete polynomial exponents differ from the reference"
    end do

    deallocate (exponents)

    ! Test with 2 dimensions
    n = 2
    kmax = 5
    nterms = poly_complete_get_nterms (n, kmax)
    allocate (exponents(n,nterms), exponents_ok(n,nterms))

    do k = 0, kmax

        status = NF_STATUS_OK
        exponents = -1
        exponents_ok = -2

        ifrom = 1
        ito = poly_complete_get_nterms (n, k)

        call polyexponents_complete (n, k, exponents(:,ifrom:ito), status)

        ! Use alternative method to obtain (correct) exponents
        call get_exponents_alt (n, k, exponents_ok(:,ifrom:ito))

        ok = exponents_identical (exponents(:,ifrom:ito), exponents_ok(:,ifrom:ito))

        if (.not. (ok .and. status == NF_STATUS_OK)) &
            error stop "Complete polynomial exponents differ from the reference"
    end do

    deallocate (exponents, exponents_ok)

    ! Test with 3 dimensions
    n = 3
    kmax = 5
    nterms = poly_complete_get_nterms (n, kmax)
    allocate (exponents(n,nterms), exponents_ok(n,nterms))

    do k = 0, kmax
        status = NF_STATUS_OK
        exponents = -1
        exponents_ok = -2

        ifrom = 1
        ito = poly_complete_get_nterms (n, k)

        call polyexponents_complete (n, k, exponents(:,ifrom:ito), status)

        ! Use alternative method to obtain (correct) exponents
        call get_exponents_alt (n, k, exponents_ok(:,ifrom:ito))

        ok = exponents_identical (exponents(:,ifrom:ito), exponents_ok(:,ifrom:ito))
        if (.not. (ok .and. status == NF_STATUS_OK)) &
            error stop "Complete polynomial exponents differ from the reference"
    end do

    deallocate (exponents, exponents_ok)

    ! Test with 4 dimensions
    n = 4
    kmax = 5
    nterms = poly_complete_get_nterms (n, kmax)
    allocate (exponents(n,nterms), exponents_ok(n,nterms))

    do k = 0, kmax
        status = NF_STATUS_OK
        exponents = -1
        exponents_ok = -2

        ifrom = 1
        ito = poly_complete_get_nterms (n, k)

        call polyexponents_complete (n, k, exponents(:,ifrom:ito), status)

        ! Use alternative method to obtain (correct) exponents
        call get_exponents_alt (n, k, exponents_ok(:,ifrom:ito))

        ok = exponents_identical (exponents(:,ifrom:ito), exponents_ok(:,ifrom:ito))
        if (.not. (ok .and. status == NF_STATUS_OK)) &
            error stop "Complete polynomial exponents differ from the reference"
    end do

    deallocate (exponents, exponents_ok)

    ! Test with 5 dimensions
    n = 5
    kmax = 5
    nterms = poly_complete_get_nterms (n, kmax)
    allocate (exponents(n,nterms), exponents_ok(n,nterms))

    do k = 0, kmax
        status = NF_STATUS_OK
        exponents = -1
        exponents_ok = -2

        ifrom = 1
        ito = poly_complete_get_nterms (n, k)

        call polyexponents_complete (n, k, exponents(:,ifrom:ito), status)

        ! Use alternative method to obtain (correct) exponents
        call get_exponents_alt (n, k, exponents_ok(:,ifrom:ito))

        ok = exponents_identical (exponents(:,ifrom:ito), exponents_ok(:,ifrom:ito))
        if (.not. (ok .and. status == NF_STATUS_OK)) &
            error stop "Complete polynomial exponents differ from the reference"
    end do

    deallocate (exponents, exponents_ok)

end subroutine

subroutine get_exponents_alt (n, k, exponents)
    !*  Alternative algorithm to find all exponent permutations, which
    !   however consumes a lot of memory and fails for larger (n, k) as
    !   it would exhaust all computer memory.
    integer, intent(in) :: n, k
    integer, intent(out), dimension(:,:) :: exponents

    integer :: m, xp, i, j
    integer, dimension(:), allocatable :: shp, work
    integer, dimension(:,:), allocatable :: exponents_all

    ! Total number of permutations of [0,...,k] exponents in N variables
    m = (k + 1) ** n

    allocate (work(m))
    call arange (work)

    allocate (shp(n), source=k+1)
    allocate (exponents_all(m, n))

    call ind2sub (shp, work, exponents_all)
    exponents_all = exponents_all - 1

    deallocate (work, shp)
    allocate (work(n))

    i = 1
    do j = 1, m
        work(:) = exponents_all(j, :)
        xp = sum(work)
        if (xp <= k) then
            ! Note: output is expected to be in transposed format
            exponents(:,i) = work
            i = i + 1
        end if
    end do

    deallocate (work)
    deallocate (exponents_all)

end subroutine

function exponents_identical (xp1, xp2) result(res)
    !*  Verify whether two sets of exponent permutations are identical.
    !
    !   The routine sorts both arrays in increasing order of both
    !   exponent and dimension.

    integer, intent(in), dimension(:,:) :: xp1, xp2
    logical :: res

    integer (int64), dimension(:), allocatable :: li1, li2, ii1, ii2
    integer :: nterms, i

    nterms = size(xp1, 2)

    ! We need to sort exponents into the same order to be able to compare them
    allocate (li1(nterms), li2(nterms))

    call create_linear_index (xp1, li1)
    call create_linear_index (xp2, li2)

    allocate (ii1(nterms), ii2(nterms))

    call argsort (li1, ii1)
    call argsort (li2, ii2)

    res = .true.
    do i = 1, nterms
        res = res .and. all(xp1(:,ii1(i)) == xp2(:,ii2(i)))
    end do

    deallocate (li1, li2)
    deallocate (ii1, ii2)

end function

subroutine create_linear_index (xp, idx)
    !*  Map each permutation of exponents to a unique linear index to
    !   establish a sort order for any list of exponent permutations.
    integer, intent(in), dimension(:,:) :: xp
    integer (int64), intent(out), dimension(:) :: idx

    integer :: n, i, nterms, j
    integer (int64) :: v

    n = size(xp, 1)
    nterms = size(xp, 2)

    do j = 1, nterms
        v = 0
        do i = 1, n
            v = v + xp(i, j) * (10_int64 ** (n-i))
        end do
        idx(j) = v
    end do

end subroutine

subroutine test_basis ()
    !* Unit tests for complete polynomial basis construction routines

    real (PREC), dimension(:,:), allocatable :: x, basis, basis_ok
    integer, parameter :: KMAX = 5
    integer :: nterms, k, n, nx, i, j
    type (status_t) :: status
    logical :: ok

    ! Test input checks

    allocate (x(1,1), basis(1,1))
    status = NF_STATUS_OK
    call polybasis_complete (x, -1, basis, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid polynomial degree"
    deallocate (x, basis)

    allocate (x(1,2), basis(1,1))
    status = NF_STATUS_OK
    call polybasis_complete (x, 0, basis, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "BASIS input has too few columns"
    deallocate (x, basis)

    allocate (x(1,1), basis(1,1))
    status = NF_STATUS_OK
    call polybasis_complete (x, 1, basis, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "BASIS input has too few rows"
    deallocate (x, basis)

    ! Test with 1d
    n = 1
    nx = 10
    allocate (x(n, nx))
    call random_number (x)

    do k = 0, KMAX
        nterms = poly_complete_get_nterms (k, n)
        allocate (basis(nterms, nx), basis_ok(nterms, nx))

        call polybasis_complete (x, k, basis, status=status)

        do i = 1, k+1
            do j = 1, nx
                basis_ok(i, j) = x(1,j) ** (i-1)
            end do
        end do

        ok = all_close (basis, basis_ok, atol=1.0e-10_PREC, rtol=0.0_PREC)
        if (.not. (ok .and. status == NF_STATUS_OK)) &
            error stop "Complete polynomial basis differs from the reference"

        deallocate (basis, basis_ok)
    end do

    deallocate (x)

    ! Test with 2d
    n = 2
    nx = 10
    allocate (x(n, nx))
    call random_number (x)

    do k = 0, 3
        nterms = poly_complete_get_nterms (k, n)
        allocate (basis(nterms, nx), basis_ok(nterms, nx))

        call polybasis_complete (x, k, basis, status=status)

        basis_ok(1,:) = 1.0
        if (k >= 1) then
            basis_ok(2,:) = x(1,:)
            basis_ok(3,:) = x(2,:)
        end if

        if (k >= 2) then
            basis_ok(4,:) = x(1,:) ** 2.0
            basis_ok(5,:) = x(1,:) * x(2,:)
            basis_ok(6,:) = x(2,:) ** 2.0
        end if

        if (k >= 3) then
            basis_ok(7,:) = x(1,:) ** 3.0
            basis_ok(8,:) = x(1,:) ** 2.0 * x(2,:)
            basis_ok(9,:) = x(1,:) * x(2,:) ** 2.0
            basis_ok(10,:) = x(2,:) ** 3.0
        end if

        ok = all_close (basis, basis_ok, atol=1.0e-10_PREC, rtol=0.0_PREC)
        if (.not. (ok .and. status == NF_STATUS_OK)) &
            error stop "Complete polynomial basis differs from the reference"

        deallocate (basis, basis_ok)
    end do

    deallocate (x)

    ! Test with 3d
    n = 3
    nx = 10
    allocate (x(n, nx))
    call random_number (x)

    do k = 0, 3
        nterms = poly_complete_get_nterms (k, n)
        allocate (basis(nterms, nx), basis_ok(nterms, nx))

        call polybasis_complete (x, k, basis, status=status)

        basis_ok(1,:) = 1.0
        if (k >= 1) then
            basis_ok(2,:) = x(1,:)
            basis_ok(3,:) = x(2,:)
            basis_ok(4,:) = x(3,:)
        end if

        if (k >= 2) then
            basis_ok(5,:) = x(1,:) ** 2.0
            basis_ok(6,:) = x(1,:) * x(2,:)
            basis_ok(7,:) = x(1,:) * x(3,:)
            basis_ok(8,:) = x(2,:) ** 2.0
            basis_ok(9,:) = x(2,:) * x(3,:)
            basis_ok(10,:) = x(3,:) ** 2.0
        end if

        if (k >= 3) then
            basis_ok(11,:) = x(1,:) ** 3.0
            basis_ok(12,:) = x(1,:) ** 2.0 * x(2,:)
            basis_ok(13,:) = x(1,:) ** 2.0 * x(3,:)
            basis_ok(14,:) = x(1,:) * x(2,:) ** 2.0
            basis_ok(15,:) = x(1,:) * x(2,:) * x(3,:)
            basis_ok(16,:) = x(1,:) * x(3,:) ** 2.0
            basis_ok(17,:) = x(2,:) ** 3.0
            basis_ok(18,:) = x(2,:) ** 2.0 * x(3,:)
            basis_ok(19,:) = x(2,:) * x(3,:) ** 2.0
            basis_ok(20,:) = x(3,:) ** 3.0
        end if

        ok = all_close (basis, basis_ok, atol=1.0e-10_PREC, rtol=0.0_PREC)
        if (.not. (ok .and. status == NF_STATUS_OK)) &
            error stop "Complete polynomial basis differs from the reference"

        deallocate (basis, basis_ok)
    end do

    deallocate (x)

end subroutine

subroutine test_basis_jac ()
    !* Unit tests for complete polynomial basis construction routines

    real (PREC), dimension(:,:), allocatable :: jac, jac_ok
    real (PREC), dimension(:), allocatable :: x
    integer :: nterms, k, n, i
    type (status_t) :: status
    logical :: ok

    ! Test input checks

    allocate (x(1), jac(1,1))
    status = NF_STATUS_OK
    call polybasis_jac_complete (x, -1, jac, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "Invalid polynomial degree"
    deallocate (x, jac)

    allocate (x(2), jac(1,1))
    status = NF_STATUS_OK
    call polybasis_jac_complete (x, 0, jac, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) &
        error stop "JAC input has too few columns"
    deallocate (x, jac)

    allocate (x(1), jac(1,1))
    status = NF_STATUS_OK
    call polybasis_jac_complete (x, 1, jac, status=status)
    if (.not. (status == NF_STATUS_INVALID_ARG)) error stop "JAC input has too few rows"
    deallocate (x, jac)

    ! Test with 1d
    n = 1
    allocate (x(n))
    call random_number (x)

    do k = 0, 5
        nterms = poly_complete_get_nterms (k, n)
        allocate (jac(nterms, n), jac_ok(nterms, n))

        call polybasis_jac_complete (x, k, jac, status=status)

        jac_ok(1,:) = 0.0
        do i = 1, k
            jac_ok(i+1, n) = i * x(1) ** (i-1)
        end do

        ok = all_close (jac, jac_ok, atol=1.0e-10_PREC, rtol=0.0_PREC)
        if (.not. (ok .and. status == NF_STATUS_OK)) &
            error stop "Complete polynomial basis derivatives differ from the " // &
                "reference"

        deallocate (jac, jac_ok)
    end do

    deallocate (x)

    ! Test with 2d
    n = 2
    allocate (x(n))
    call random_number (x)

    do k = 0, 3
        nterms = poly_complete_get_nterms (k, n)
        allocate (jac(nterms, n), jac_ok(nterms, n))

        call polybasis_jac_complete (x, k, jac, status=status)

        jac_ok(:,:) = 0.0
        if (k >= 1) then
            jac_ok(2,1) = 1.0
            jac_ok(3,2) = 1.0
        end if

        if (k >= 2) then
            jac_ok(4,1) = 2.0 * x(1)
            jac_ok(5,1) = x(2)
            jac_ok(5,2) = x(1)
            jac_ok(6,2) = 2.0 * x(2)
        end if

        if (k >= 3) then
            jac_ok(7,1) = 3.0 * x(1) ** 2.0
            jac_ok(8,1) = 2.0 * x(1) * x(2)
            jac_ok(8,2) = x(1) ** 2.0
            jac_ok(9,1) = x(2) ** 2.0
            jac_ok(9,2) = x(1) * 2.0 * x(2)
            jac_ok(10,2) = 3.0 * x(2) ** 2.0
        end if

        ok = all_close (jac, jac_ok, atol=1.0e-10_PREC, rtol=0.0_PREC)
        if (.not. (ok .and. status == NF_STATUS_OK)) &
            error stop "Complete polynomial basis derivatives differ from the reference"

        deallocate (jac, jac_ok)
    end do

    deallocate (x)

    ! Test with 3d
    n = 3
    allocate (x(n))
    call random_number (x)

    do k = 0, 3
        nterms = poly_complete_get_nterms (k, n)
        allocate (jac(nterms, n), jac_ok(nterms, n))

        call polybasis_jac_complete (x, k, jac, status=status)

        jac_ok(:,:) = 0.0
        if (k >= 1) then
            jac_ok(2,1) = 1.0
            jac_ok(3,2) = 1.0
            jac_ok(4,3) = 1.0
        end if

        if (k >= 2) then
            jac_ok(5,1) = 2.0 * x(1)
            jac_ok(6,1) = x(2)
            jac_ok(6,2) = x(1)
            jac_ok(7,1) = x(3)
            jac_ok(7,3) = x(1)
            jac_ok(8,2) = 2.0 * x(2)
            jac_ok(9,2) = x(3)
            jac_ok(9,3) = x(2)
            jac_ok(10,3) = 2.0 * x(3)
        end if

        if (k >= 3) then
            jac_ok(11,1) = 3.0 * x(1) ** 2.0
            jac_ok(12,1) = 2.0 * x(1) * x(2)
            jac_ok(12,2) = x(1) ** 2.0
            jac_ok(13,1) = 2.0 * x(1) * x(3)
            jac_ok(13,3) = x(1) ** 2.0
            jac_ok(14,1) = x(2) ** 2.0
            jac_ok(14,2) = x(1) * 2.0 * x(2)
            jac_ok(15,1) = x(2) * x(3)
            jac_ok(15,2) = x(1) * x(3)
            jac_ok(15,3) = x(1) * x(2)
            jac_ok(16,1) = x(3) ** 2.0
            jac_ok(16,3) = x(1) * 2.0 * x(3)
            jac_ok(17,2) = 3.0 * x(2) ** 2.0
            jac_ok(18,2) = 2.0 * x(2) * x(3)
            jac_ok(18,3) = x(2) ** 2.0
            jac_ok(19,2) = x(3) ** 2.0
            jac_ok(19,3) = x(2) * 2.0 * x(3)
            jac_ok(20,3) = 3.0 * x(3) ** 2.0
        end if

        ok = all_close (jac, jac_ok, atol=1.0e-10_PREC, rtol=0.0_PREC)
        if (.not. (ok .and. status == NF_STATUS_OK)) &
            error stop "Complete polynomial basis derivatives differ from the reference"

        deallocate (jac, jac_ok)
    end do

    deallocate (x)

end subroutine

end
