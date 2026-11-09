program test_arrays_create
    !* Unit tests for array creation routines.

    use numfort_arrays, only: arange

    implicit none

    call test_arange ()

    contains

subroutine test_arange ()
    !* - Check default-integer ranges with omitted and explicit start/step values.
    !  - Preserve surrounding sentinels for the calibration index slice.

    integer, parameter :: INTSIZE = kind (0)
    integer (INTSIZE) :: idx_left(0:6), values(4), nleft

    ! An empty calibration slice leaves every surrounding sentinel unchanged.
    idx_left = -99
    nleft = 0
    call arange (idx_left(1:nleft), 2)
    if (any (idx_left /= [-99, -99, -99, -99, -99, -99, -99])) &
        error stop "arange: empty slice changed sentinels"

    ! A singleton calibration slice contains only its starting index, 2.
    idx_left = -99
    nleft = 1
    call arange (idx_left(1:nleft), 2)
    if (any (idx_left /= [-99, 2, -99, -99, -99, -99, -99])) &
        error stop "arange: singleton slice or sentinels"

    ! Four remaining indices start at 2; entries outside the slice stay intact.
    idx_left = -99
    nleft = 4
    call arange (idx_left(1:nleft), 2)
    if (any (idx_left /= [-99, 2, 3, 4, 5, -99, -99])) &
        error stop "arange: calibration slice or sentinels"

    ! Omitting both optional arguments starts at 1 with unit step.
    values = -99
    call arange (values)
    if (any (values /= [1, 2, 3, 4])) &
        error stop "arange: default start and step"

    ! A positive nonunit step of 3 gives four values starting at 2.
    values = -99
    call arange (values, 2, step=3)
    if (any (values /= [2, 5, 8, 11])) &
        error stop "arange: positive nonunit step"

    ! A negative step of -3 gives a descending range crossing zero.
    values = -99
    call arange (values, 2, step=-3)
    if (any (values /= [2, -1, -4, -7])) &
        error stop "arange: negative step"

    ! A zero step repeats the starting value under the current contract.
    values = -99
    call arange (values, 2, step=0)
    if (any (values /= [2, 2, 2, 2])) &
        error stop "arange: zero step"

    ! Supplying only the step retains the default starting value, 1.
    values = -99
    call arange (values, step=2)
    if (any (values /= [1, 3, 5, 7])) &
        error stop "arange: default start with explicit step"
end subroutine

end program
