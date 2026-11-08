
#include <numfort.h>

program test_io_fixed
    !*  Unit tests for routines performing fixed-format text I/O.

    use, intrinsic :: iso_fortran_env

    use numfort_stats
    use numfort_common
    use numfort_common_testing
    use numfort_io

    implicit none

    call test_fixed_real32 ()
    call test_fixed_real64 ()

    call test_fixed_alloc_real32 ()
    call test_fixed_alloc_real64 ()

    call test_fixed_alloc_malformed_real32 ()
    call test_fixed_alloc_malformed_real64 ()

    contains

#include <numfort_real32.h>
#include "io_fixed_impl.F90"

#include <numfort_real64.h>
#include "io_fixed_impl.F90"

end program

