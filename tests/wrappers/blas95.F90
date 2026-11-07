

#include <numfort.h>

program test_wrappers_blas95

    use, intrinsic :: iso_fortran_env

    use blas95
    use numfort_common
    use numfort_common_testing
    use numfort_stats, only: set_seed

    implicit none

    call test_all ()

    contains

subroutine test_all ()

    call test_gemv_real32 ()
    call test_gemv_real64 ()

    call test_gemm_real32 ()
    call test_gemm_real64 ()

end subroutine

#include <numfort_real32.h>
#include "blas95_gemv.F90"
#include "blas95_gemm.F90"

#include <numfort_real64.h>
#include "blas95_gemv.F90"
#include "blas95_gemm.F90"

end program
