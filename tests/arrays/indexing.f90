program test_numfort_indexing

    use iso_fortran_env
    use numfort_arrays
    implicit none

    integer, parameter :: INTSIZE = int32

    call test_all ()

contains

subroutine test_all ()

    call test_ind2sub ()
    ! call test_sub2ind ()

end subroutine

subroutine test_ind2sub ()

    integer (INTSIZE), dimension(:), allocatable :: lin, shp
    integer (INTSIZE), dimension(:, :), allocatable :: sub, sub2
    integer :: i, n

    ! test conversion of 1-d arrays
    n = 100
    allocate (lin(n), sub(n, 1))
    lin = [(i, i=1,n)]
    sub = -1

    call ind2sub ([n], lin, sub)
    if (.not. (all (lin == sub(:, 1)))) error stop "Conversion of 1d linear indices"
    deallocate (lin, sub)

    ! check some simple transformations
    allocate (shp(2))
    shp = [7, 8]
    n = product(shp)
    allocate (lin(n), sub(n, size(shp)), sub2(n, size(shp)))
    lin = [(i, i=1,n)]

    call ind2sub (shp, lin, sub)

    ! manually compute subindices
    sub2(:, 1) = modulo (lin - 1, shp(1)) + 1
    sub2(:, 2) = ceiling (lin / real (shp(1)))

    if (.not. (all (sub(1,:) == [1, 1]))) error stop "in2dsub, shp = [7, 8]"
    deallocate (lin, sub, sub2, shp)

    ! check fixed point of ind2sub / sub2ind
    call ind2sub_fp ([11])
    call ind2sub_fp ([7, 8])
    call ind2sub_fp ([3, 5, 15])
    call ind2sub_fp ([3, 5, 15, 7])
    call ind2sub_fp ([3, 5, 15, 1, 5])
end subroutine

subroutine ind2sub_fp (shp)
    !*  Check fixed point x = ind2sub(sub2ind(x))
    integer, intent(in), dimension(:) :: shp
        !* Shape used for the indexing round trip.

    integer (INTSIZE), dimension(:), allocatable :: lin, lin2
    integer (INTSIZE), dimension(:, :), allocatable :: sub
    integer :: i, n

    n = product(shp)
    allocate (lin(n), lin2(n), sub(n, size(shp)))
    lin = [(i, i=1,n)]

    call ind2sub (shp, lin, sub)
    call sub2ind (shp, sub, lin2)

    if (.not. (all (lin == lin2))) error stop "IND2SUB/SUB2IND round trip failed"
end subroutine

end program
