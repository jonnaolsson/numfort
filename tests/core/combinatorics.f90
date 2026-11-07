
program test_core_combinatorics

    use, intrinsic :: iso_fortran_env
    use numfort_arrays
    use numfort_common
    use numfort_common_testing
    use numfort_core

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()

contains

subroutine test_all ()

    call test_factorial ()
    call test_comb ()
    call test_poch ()

end subroutine

subroutine test_factorial ()

    integer, parameter :: N = 10
    integer, parameter, dimension(N) :: fact_true = &
        [ 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880]
    integer, dimension(N) :: x, fx

    ! Test factorial function for first 10 values
    call arange (x, 0)
    fx = factorial (x)

    if (.not. (all (fx == fact_true))) &
        error stop "Factorial function for values [0,...,9]"

end subroutine

subroutine test_comb ()

    integer, PARAMETER :: NMAX = 10
    logical :: ok
    integer (int64) :: res, res_true, n, k
        ! Note: need to use in64 otherwise the results for higher N won't fit
        ! into default INTEGER!

    do n = 0, NMAX
        ok = .true.
        do k = 0, n
            res = comb (n, k, repetition=.true.)
            res_true = factorial (n+k-1) / factorial(k) / factorial (n-1)
            ok = ok .and. (res == res_true)
        end do

        if (.not. (ok)) error stop "Binomial coefficients differ from the reference"
    end do

    do n = 0, NMAX
        ok = .true.
        do k = 0, n
            res = comb (n, k, repetition=.false.)
            res_true = factorial (n) / factorial(k) / factorial (n-k)
            ok = ok .and. (res == res_true)
        end do

        if (.not. (ok)) error stop "Binomial coefficients differ from the reference"
    end do

end subroutine

subroutine test_poch ()
    !*  Unit tests for Pochhammer symbol

    integer :: n, k, res, res_ok

    ! Test with empty product
    res_ok = 1
    res = poch (0, 0)
    if (.not. (res == res_ok)) error stop "POCH(0,0) == 1"

    do n = 1, 5
        do k = 0, 5
            res = poch (n, k)
            res_ok = factorial (k) * comb (n+k-1, k)
            if (.not. (res == res_ok)) &
                error stop "POCH failed the factorial/binomial identity check"
        end do
    end do

end subroutine

end
