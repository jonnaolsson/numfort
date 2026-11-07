program compiler_test

    use iso_fortran_env, only: real64
    use numfort, only: NUMFORT_VERSION
    use numfort_arrays, only: linspace
    use numfort_linalg, only: inv
    use numfort_stats, only: set_seed
    implicit none

    real(real64) :: a(2, 2), ainv(2, 2), x(2)

    ! Check public modules, bundled modules, and the BLAS/LAPACK link interface.
    if (len(NUMFORT_VERSION) == 0) stop 1
    call set_seed(123)
    call linspace(x, 1.0_real64, 2.0_real64)
    a = 0.0_real64
    a(1, 1) = x(1)
    a(2, 2) = x(2)
    call inv(a, ainv)
    if (any(ainv /= ainv)) stop 1

end program
