

real (PREC), parameter :: DEFAULT_L1_RATIO = 0.5_PREC

real (PREC), parameter :: DEFAULT_COORD_DESCENT_TOL = 1.0e-4_PREC
integer, parameter :: DEFAULT_COORD_DESCENT_MAXITER = 1000
logical, parameter :: DEFAULT_POSITIVE = .false.

type :: enet_config
    logical :: trans_x = .false.
        !*  Transpose regressor matrix
    logical :: trans_y = .false.
        !*  Transpose outcome matrix (multi-outcome model only)
    logical :: trans_coefs = .false.
        !*  Transpose coefficient matrix (multi-outcome model only)
        !   The transposed coefficient matrix is assumed to have shape
        !   [Nlhs,Ncoefs].
    logical :: center = .true.
        !*  Center regressors and outcome variables
    logical :: scale = .true.
        !*  Scale regressors (features) to have unit variance
    logical :: add_intercept = .false.
        !*  If true, prepend intercept as the first element (or row)
        !   in the coefficient array.

    logical :: drop_const = .true.
        !*  Drop any constant regressors (features)
    real (PREC) :: tol_const = 1.0e-14_PREC
        !*  Tolerance for variance of a variable in X to be considered
        !   as a constant and dropped from the model.

    integer :: cv_n = 20
        !*  Number of chunks to use for k-fold cross-validation
    logical :: cv_parsimonious = .false.
        !*  If true, pick the most parsimonious model within one standard
        !   deviation of the minimal MSE.

    integer :: alpha_n = 100
        !*  Number of grid points used for alpha regularization parameter
        !   for cross-validation
    real (PREC) :: alpha_eps = 1.0e-4_PREC
        !*  Scaling parameter to determine lower bound of alpha grid as a
        !   fraction of the upper bound.

    integer :: maxiter = 1000
        !*  Max. number of iterations for coordinate descent algorithm
    real (PREC) :: tol = 1.0e-4_PREC
        !*  Termination tolerance for coordinate descent algorithm
    logical :: positive = .false.

    logical :: force_multi = .false.
        !*  Force code path for multiple outcomes even if outcome data
        !   contains only one LHS variable (internal use only!)
end type


type :: enet_result
    type (enet_config) :: conf

    real (PREC), dimension(:), allocatable :: coefs
        !*  Coefficient array of estimated model
    real (PREC) :: intercept = 0.0
        !*  Intercept of estimated model

    real (PREC), dimension(:,:), allocatable :: coefs_multi
    real (PREC), dimension(:), allocatable :: intercept_multi

    integer, dimension(:), allocatable :: irhs
        !*  Indices of variables that were actually included in the model
        !   (ie. excludes variables that were dropped because they were
        !   constant, non-finite, etc.)

    real (PREC) :: l1_ratio = -1.0_PREC
    real (PREC) :: alpha = -1.0_PREC

    integer :: nobs = 0
        !*  Number of observations
    integer :: nrhs = 0
        !*  Number of explanatory (RHS) variables
    integer :: nlhs = 1
        !*  Number of outcome (LHS) varialbes

    real (PREC) :: var_rhs = 0.0
        !*  Contains frac. of RHS variance included in the model, if
        !   applicable (used for Ridge regression)
end type



interface enet_cv
    procedure enet_cv, enet_cv_multi
end interface

interface enet_path
    procedure enet_path, enet_path_multi
end interface

interface enet_path_mse
    procedure enet_path_mse, enet_path_mse_multi
end interface

interface enet_fit
    procedure enet_fit, enet_fit_multi
end interface

interface predict
    procedure enet_predict
end interface

interface enet_predict
    procedure enet_predict, enet_predict_multi
end interface

interface post_estim
    procedure enet_post_estim_impl
end interface

interface enet_post_estim
    procedure enet_post_estim_impl
end interface

interface create_alpha_grid
    procedure create_alpha_grid
end interface

interface create_alpha_grid_cv
    procedure create_alpha_grid_cv, create_alpha_grid_cv_multi
end interface

interface check_dims
    procedure check_dims, check_dims_multi
end interface

interface get_alpha_grid_max
    procedure get_alpha_grid_max, get_alpha_grid_max_multi
end interface

interface ridge
    procedure ridge_single, ridge_multi
end interface
