

subroutine __APPEND(test_gemv,__PREC) ()
    !* TEST_GEMV implements unit tests for the F95 wrapper for GEMV.
    use blas_interfaces, only: BLAS_GEMV => GEMV

    integer, parameter  :: PREC = __PREC

    integer :: m, n, lda
    integer, parameter :: mrange(*) = [1, 2, 10]
    integer, parameter :: nrange(*) = [1, 2, 10]
    integer, parameter :: incx = 1, incy = 1
    real (PREC), dimension(:,:), allocatable :: a
    real (PREC), dimension(:), allocatable :: x, y, y95
    real (PREC) :: alpha, beta
    character (1) :: trans
    character (1) :: trans_range(2) = ["N", "T"]

    integer :: i, j, k

    call set_seed (123)

    do k = 1, size(trans_range)
        trans = trans_range(k)
        do i = 1, size(mrange)
            m = mrange(i)
            lda = mrange(i)

            do j = 1, size(nrange)
                n = nrange(j)

                allocate (a(m, n))
                if (trans == "N") then
                    allocate (x(n))
                    allocate (y(m), y95(m))
                else
                    allocate (x(m))
                    allocate (y(n), y95(n))
                end if

                call random_number (a)
                call random_number (x)

                ! === Call with default ALPHA, BETA ===
                alpha = 1.0
                beta = 0.0

                ! Call F77 interface
                call BLAS_GEMV (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)

                ! Call F95 interface (with default alpha, beta) and default
                ! TRANS (if applicable)
                if (trans == "N") then
                    call GEMV (a, x, y95)
                else
                    call GEMV (a, x, y95, trans=trans)
                end if

                if (.not. (all (y == y95))) &
                    error stop "GEMV with default arguments differs from BLAS"

                call GEMV (a, x, y95, alpha=alpha, beta=beta, trans=trans)

                if (.not. (all (y == y95))) &
                    error stop "GEMV with explicit default scaling differs from BLAS"

                ! === Call with non-default ALPHA, BETA ===
                y95(:) = y
                call random_number (alpha)
                call random_number (beta)
                ! Call F77 interface
                call BLAS_GEMV (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
                call GEMV (a, x, y95, alpha=alpha, beta=beta, trans=trans)

                if (.not. (all (y == y95))) &
                    error stop "GEMV with scaling differs from BLAS"

                deallocate (a, x, y, y95)

            end do
        end do
    end do

end subroutine
