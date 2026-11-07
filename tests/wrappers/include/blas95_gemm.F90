

subroutine __APPEND(test_gemm,__PREC) ()
    !* TEST_GEMM implements unit tests for the F95 wrapper for GEMM.
    use blas_interfaces, only: BLAS_GEMM => GEMM

    integer, parameter  :: PREC = __PREC

    integer :: m, n, k, lda, ldb, ldc
    integer, parameter :: mrange(*) = [1, 2, 10]
    integer, parameter :: nrange(*) = [1, 2, 10]
    integer, parameter :: krange(*) = [1, 2, 10]
    real (PREC), dimension(:,:), allocatable :: a, b, c, c95
    real (PREC) :: alpha, beta
    character (1) :: transa, transb
    character (1) :: trans_range(2) = ["N", "T"]

    integer :: i, j, ik, ita, itb

    call set_seed (123)

    do ita = 1, size(trans_range)
        transa = trans_range(ita)
        do itb = 1, size(trans_range)
            transb = trans_range(itb)

            do i = 1, size(mrange)
                m = mrange(i)

                do j = 1, size(nrange)
                    n = nrange(j)

                    do ik = 1, size(krange)
                        k = krange(ik)

                        if (transa == "N") then
                            allocate (a(m, k))
                            lda = m
                        else
                            allocate (a(k, m))
                            lda = k
                        end if

                        if (transb == "N") then
                            allocate (b(k, n))
                            ldb = k
                        else
                            allocate (b(n, k))
                            ldb = n
                        end if

                        allocate (c(m, n), c95(m, n))
                        ldc = m

                        call random_number (a)
                        call random_number (b)
                        call random_number (c)

                        c95(:,:) = c

                        ! === Call with default ALPHA, BETA ===
                        alpha = 1.0
                        beta = 0.0

                        ! Call F77 interface
                        call BLAS_GEMM (transa, transb, m, n, k, alpha, a, lda, &
                            b, ldb, beta, c, ldc)

                        ! Call F95 interface (with default alpha, beta) and default
                        ! TRANS{A,B} (if applicable)
                        if (transa == "N" .and. transb == "N") then
                            call GEMM (a, b, c95)
                        else
                            call GEMM (a, b, c95, transa=transa, transb=transb)
                        end if

                        if (.not. (all (c == c95))) &
                            error stop "GEMM with default arguments differs from BLAS"

                        call GEMM (a, b, c95, alpha=alpha, beta=beta, &
                            transa=transa, transb=transb)

                        if (.not. (all (c == c95))) &
                            error stop "GEMM with explicit default scaling " // &
                                "differs from BLAS"

                        ! === Call with non-default ALPHA, BETA ===
                        c95(:,:) = c
                        call random_number (alpha)
                        call random_number (beta)
                        ! Call F77 interface
                        call BLAS_GEMM (transa, transb, m, n, k, alpha, a, lda, &
                            b, ldb, beta, c, ldc)
                        call GEMM (a, b, c95, alpha=alpha, beta=beta, &
                            transa=transa, transb=transb)

                        if (.not. (all (c == c95))) &
                            error stop "GEMM with scaling differs from BLAS"

                        deallocate (a, b, c, c95)

                    end do
                end do
            end do
        end do
    end do

end subroutine
