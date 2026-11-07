#include <numfort.h>

program test_numfort_indexing
    !*  Unit tests for array indexing routines (ind2sub, sub2ind, shape2sub).

    use, intrinsic :: iso_fortran_env
    use numfort_arrays

    implicit none

    call test_all ()

    contains

subroutine test_all ()
    !*  Run all indexing unit tests for 32-bit and 64-bit integers.

    call test_indexing_int32 ()
    call test_indexing_int64 ()

end subroutine


#include <numfort_int32.h>
#include "test_indexing_impl.F90"

#include <numfort_int64.h>
#include "test_indexing_impl.F90"

end program
