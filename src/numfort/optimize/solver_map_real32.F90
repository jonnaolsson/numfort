


module numfort_optimize_solver_map_real32

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_common_status
    use numfort_common_input_checks
    use numfort_core_logexp, only: expit
    use numfort_optimize_solver_map_common

    implicit none

    private

    public :: solver_map
    public :: solver_map_init
    public :: solver_map_eval
    public :: solver_map_eval_inverse

    integer, parameter :: PREC = real32

    type :: solver_map
        private
        integer :: transform = TRANSFORM_LINEAR
        real (PREC) :: lb
        real (PREC) :: ub
        real (PREC) :: scale = 1.0_PREC
    end type


    interface solver_map_init
        procedure solver_map_init
    end interface

    interface solver_map_eval
        procedure solver_map_eval_scalar, &
            solver_map_eval_diag, &
            solver_map_eval_matrix
    end interface

    interface solver_map_eval_inverse
        procedure solver_map_eval_inverse_scalar, &
            solver_map_eval_inverse_diag, &
            solver_map_eval_inverse_matrix
    end interface

    contains


#include "solver_map_impl.F90"

end module
