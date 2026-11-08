

#include <numfort.h>


module numfort_core_logexp
    !*  Module implements various special functions related to log and
    !   exponentials.

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, &
        IEEE_QUIET_NAN, IEEE_POSITIVE_INF, IEEE_NEGATIVE_INF

    use numfort_core_libc, only: log1p

    implicit none
    private


    public :: logaddexp
    public :: expit
    public :: logit

#include <numfort_real32.h>
#include "logexp_spec.F90"

#include <numfort_real64.h>
#include "logexp_spec.F90"

    contains


#include <numfort_real32.h>
#include "logexp_impl.F90"


#include <numfort_real64.h>
#include "logexp_impl.F90"


end module
