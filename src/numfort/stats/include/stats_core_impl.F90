


pure subroutine get_data_dims (x, dim, nvars, nobs, status)
    !*  GET_DATA_DIMS sets the number of variables and number of observations
    !   of for a given 2-dimensional array.
    real (PREC), intent(in), dimension(:,:) :: x
        !*  Data array
    integer, intent(in) :: dim
        !*  Dimension along with reduction should be performed
    integer, intent(out) :: nvars
        !*  Number of variables
    integer, intent(out) :: nobs
        !*  Number of observations
    type (status_t), intent(out), optional :: status
        !*  Exit code (either NF_STATUS_OK or NF_STATUS_INVALID_ARG)

    nobs = 0
    nvars = 0

    if (present(status)) then
        status = NF_STATUS_OK
        if (dim /= 1 .and. dim /= 2) then
            status = NF_STATUS_INVALID_ARG
            return
        end if
    end if

    nvars = size(x, 3-dim)
    nobs = size(x, dim)
end subroutine



pure subroutine mean_1d (x, m, dim, status)
    !*  MEAN computes the mean of values in array X.
    real (PREC), intent(in), dimension(:), contiguous, target :: x
        !*  Array containing numbers for which mean should be computed.
    real (PREC), intent(out) :: m
        !*  Computed mean
    integer, intent(in), optional :: dim
        !*  Ignored for 1d-arrays, present to support same API as for
        !   higher-rank arrays
    type (status_t), intent(out), optional :: status
        !*  Exit code

    type (status_t) :: lstatus

    ! initialize to default values
    lstatus = NF_STATUS_OK

    ! --- Input checks ---

    call check_cond (size(x) > 0, 'MEAN', 'X: Non-zero number of observations required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Implementation ---

    call mean_impl_1d (x, m, dim=1)

100 continue
    if (present(status)) status = lstatus
end subroutine



pure subroutine mean_2d (x, m, dim, status)
    !*  MEAN computes the mean of array elements along a dimension.
    real (PREC), intent(in), dimension(:,:), contiguous :: x
        !*  Input data. Variables can either be stored column-wise (default)
        !   or row-wise. (see DIM argument)
    real (PREC), intent(out), dimension(:), contiguous :: m
        !*  Array of means for each variable in X
    integer, intent(in), optional :: dim
        !*  If present, specifies the dimension along which to compute means
        !   (default: dim=1 implies that variables are stored column-wise)
    type (status_t), intent(out), optional :: status
        !*  Optional exit status

    type (status_t) :: lstatus
    ! iv indexes variables, iobs individual observations
    integer :: nvars, nobs, ldim

    ! default values
    lstatus = NF_STATUS_OK

    ! initialize to something so gfortran does not complain
    nvars = -1
    nobs = -1

    ! --- Input checks ---

    call set_optional_arg (dim, 1, ldim)

    call get_data_dims (x, ldim, nvars, nobs, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(m) == nvars, 'MEAN', 'M: Non-conformable array size', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (nobs > 0, 'MEAN', 'X: Non-zero number of observations required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Implementation ---

    call mean_impl (x, m, ldim)

100 continue
    if (present(status)) status = lstatus
end subroutine



pure subroutine mean_impl_1d (x, m, dim)
    !*  MEAN computes the mean of values in array X.
    !
    !   Implementation routine, does NOT perform any input checks.
    real (PREC), intent(in), dimension(:), contiguous, target :: x
        !*  Array containing numbers for which mean should be computed.
    real (PREC), intent(out) :: m
        !*  Computed mean
    integer, intent(in), optional :: dim
        !*  Ignored for 1d-arrays, present to support same API as for higher-rank
        !   arrays

    integer :: iobs, nobs
    real (PREC) :: m_old, m_new

    nobs = size(x)

    ! --- Quick termination ---

    if (nobs == 0) return

    ! --- Implementation ---

    ! Implementation note: implement this directly here without calling
    ! MEAN_IMPL_2D, since in a pure procedure we cannot create a 2d pointer
    ! to point to X to construct a 2d argument.
    m_new = x(1)
    m_old = m_new

    do iobs = 2, nobs
        m_new = m_old + (x(iobs) - m_old)  / iobs
        m_old = m_new
    end do

    m = m_new

end subroutine



pure subroutine mean_impl_2d (x, m, dim)
    !*  MEAN computes the mean of array elements along a dimension.
    !
    !   Implementation routine, does NOT perform any input checks.
    real (PREC), intent(in), dimension(:,:), contiguous :: x
        !*  Input data. Variables can either be stored column-wise
        !   or row-wise. (see DIM argument)
    real (PREC), intent(out), dimension(:), contiguous, target :: m
        !*  Array of means, needs to be at least as large as the number of
        !   variables.
    integer, intent(in) :: dim
        !*  Dimension along which to compute means

    real (PREC), dimension(:), allocatable :: m_old

    ! iv indexes variables, iobs individual observations
    integer :: iv, iobs, nvars, nobs

    ! initialize to something so gfortran does not complain
    nvars = size(x, 3-dim)
    nobs = size(x, dim)

    ! --- Quick termination ---

    if (nobs == 0 .or. nvars == 0) return

    ! --- Implementation ---

    ! compute mean and std. dev. using Welford's method; see
    ! http://www.johndcook.com/blog/standard_deviation/

    if (dim == 1) then
        ! If variables are stored in columns, use 1d mean routine
        do iv = 1, nvars
            call mean_impl (x(:,iv), m(iv))
        end do
    else
        ! Compute running mean and std if variables are stored in rows
        allocate (m_old(nvars))
        m = x(:, 1)
        m_old(:) = m

        ! Compute running means for all variables simultaneuosly
        do iobs = 2, nobs
            m = m_old + (x(:, iobs) - m_old) / iobs
            m_old(:) = m
        end do
    end if

end subroutine



pure subroutine std_1d (x, s, m, dof, dim, status)
    !*  Compute standard deviation (and optionally mean) for data given in
    !   1d array.
    real (PREC), intent(in), dimension(:), contiguous :: x
        !   Array containing numbers for which standard deviation should be
        !   computed
    real (PREC), intent(out) :: s
        !   Contains computed standard deviation on exit.
    real (PREC), intent(out), optional :: m
        !   If present, contains computed mean on exit.
    integer, intent(in), optional :: dof
        !*  Degrees of freedom correction to use (default: 1). Computes
        !   standard deviation as sum of squared deviations divided by (N-dof).
    integer, intent(in), optional :: dim
        !*  Ignored for 1d-arrays, present to support same API as for higher-rank
        !   arrays
    type (status_t), intent(out), optional :: status
        !*  Ignored for 1d-arrays, present to support same API as for higher-rank
        !   arrays

    integer :: ldof
    type (status_t) :: lstatus

    ! initiliaze default values
    lstatus = NF_STATUS_OK

    ! --- Input checks ---

    call set_optional_arg (dof, 1, ldof)

    call check_cond (ldof == 0 .or. ldof == 1, 'STD', 'DOF: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(x) > 0, 'STD', 'X: Non-zero number of observations required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Implementation ---

    call std_impl (x, s, m, ldof)

100 continue
    if (present(status)) status = lstatus
end subroutine



pure subroutine std_2d (x, s, m, dof, dim, status)
    !*  Compute standard deviation (and optionally) mean of data in 2d-array
    !   along the desired dimension.
    real (PREC), intent(in), dimension(:,:), contiguous :: x
        !*  2d-array containing data for which standard deviation should be
        !   computed.
    real (PREC), intent(out), dimension(:), contiguous :: s
        !*  On exit, contains standard deviations along desired dimension.
        !   Array size must be sufficient to hold computed values.
    real (PREC), intent(out), dimension(:), contiguous, optional :: m
        !*  If present, contains vector of means computed along desired dimension.
        !   Array size must be sufficient to hold computed values.
    integer, intent(in), optional :: dof
        !*  Degrees of freedom correction to use (default: 1). Computes
        !   standard deviation as sum of squared deviations divided by (N-dof).
    integer, intent(in), optional :: dim
        !*  If present, specifies dimension along which std/mean should be computed.
    type (status_t), intent(out), optional :: status
        !*  If present, returns NF_STATUS_OK on exit if operation was successful,
        !   and status > 0 if an error was encountered

    type (status_t) :: lstatus
    ! iv indexes variables, iobs individual observations
    integer :: nvars, nobs, ldof, ldim

    ! set default values
    lstatus = NF_STATUS_OK

    ! --- Input checks ---

    call set_optional_arg (dof, 1, ldof)
    call set_optional_arg (dim, 1, ldim)

    ! Initialize to something so gfortran does not complain
    nvars = -1
    nobs = -1
    call get_data_dims (x, ldim, nvars, nobs, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (ldof == 0 .or. ldof == 1, 'STD', 'DOF: Invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (nobs > 0, 'STD', 'X: Non-zero number of observations required', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (size(s) == nvars, 'STD', 'S: Non-conformable array size', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(m)) then
        call check_cond (size(m) == nvars, 'STD', 'M: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Implementation ---

    call std_impl (x, s, m, ldof, ldim)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine std_impl_1d (x, s, m, dof, dim)
    !*  Compute standard deviation (and optionally mean) for data given in
    !   1d array.
    !
    !   Implementation routine, does NOT perform any input checks.
    real (PREC), intent(in), dimension(:), contiguous :: x
        !   Array containing numbers for which standard deviation should be
        !   computed
    real (PREC), intent(out) :: s
        !   Contains computed standard deviation on exit.
    real (PREC), intent(out), optional :: m
        !   If present, contains computed mean on exit.
    integer, intent(in), optional :: dof
        !*  Degrees of freedom correction to use (default: 1). Computes
        !   standard deviation as sum of squared deviations divided by (N-dof).
    integer, intent(in), optional :: dim
        !*  Ignored for 1d-arrays, present to support same API as for higher-rank
        !   arrays

    real (PREC) :: m_old, m_new, s_old, s_new, x_i
    integer :: iobs, nobs

    nobs = size(x)

    ! --- Quick termination ---

    if (nobs == 0) return

    ! --- Implementation ----

    ! compute mean and std. dev. using Welford's method; see
    ! http://www.johndcook.com/blog/standard_deviation/

    ! Note on implementation: We cannot create a 2d pointer to X in order
    ! to call the 2d implementation in a PURE procedure, so we
    ! implement the 1d case directly here.
    m_new = x(1)
    m_old = m_new
    s_new = 0.0_PREC
    s_old = 0.0_PREC

    do iobs = 2, nobs
        x_i = x(iobs)
        m_new = m_old + (x_i - m_old) / iobs
        s_new = s_old + (x_i - m_old) * (x_i - m_new)

        m_old = m_new
        s_old = s_new
    end do

    s = s_new

    if (nobs > 1) then
        s = sqrt(s_new / (nobs - dof))
    end if

    if (present(m)) m = m_new

end subroutine



pure subroutine std_impl_2d (x, s, m, dof, dim)
    !*  STD_IMPL computes the standard deviations (and optionally) means of
    !   data in 2d-array along the desired dimension.
    real (PREC), intent(in), dimension(:,:), contiguous :: x
        !*  2d-array containing data for which standard deviation should be
        !   computed.
    real (PREC), intent(out), dimension(:), contiguous :: s
        !*  On exit, contains standard deviations along desired dimension.
        !   Array size must be sufficient to hold computed values.
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: m
        !*  If present, contains vector of means computed along desired dimension.
        !   Array size must be sufficient to hold computed values.
    integer, intent(in) :: dof
        !*  Degrees of freedom correction to use (default: 1). Computes
        !   standard deviation as sum of squared deviations divided by (N-dof).
    integer, intent(in) :: dim
        !*  Dimension along which std/mean should be computed.

    real (PREC), dimension(:), pointer, contiguous :: ptr_m
    real (PREC), dimension(:), allocatable :: m_old, s_old
    integer :: iv, iobs, nvars, nobs

    nullify (ptr_m)

    ! Initialize to something so gfortran does not complain
    nvars = size(x, 3-dim)
    nobs = size(x, dim)

    ! --- Quick termination ---

    if (nobs == 0 .or. nvars == 0) return

    ! --- Implementation ---

    if (present(m)) then
        ptr_m(1:nvars) => m
    else
        allocate (ptr_m(nvars))
    end if

    ! compute mean and std. dev. using Welford's method; see
    ! http://www.johndcook.com/blog/standard_deviation/

    if (dim == 1) then
        ! If data is organized in columns we can use 1d routine without
        ! hurting cache locality
        do iv = 1, nvars
            call std_impl (x(:, iv), s(iv), ptr_m(iv), dof)
        end do
    else if (dim == 2) then
        ! Compute running mean and std if variables are stored in rows
        ! In this case we don't use 1d method as this would result in
        ! inefficient memory access patterns
        ptr_m(:) = x(:, 1)
        s = 0.0_PREC

        allocate (m_old(nvars), source=x(:, 1))
        allocate (s_old(nvars), source=0.0_PREC)

        ! Compute means / std for all variables simultaneously
        do iobs = 2, nobs
            ptr_m(:) = m_old + (x(:, iobs) - m_old) / iobs
            s(:) = s_old + (x(:, iobs) - m_old) * (x(:, iobs) - ptr_m)

            m_old(:) = ptr_m
            s_old(:) = s
        end do

        if (nobs > 1) then
            s = sqrt(s / (nobs - dof))
        end if
    end if

    ! Clean up potentially allocated array
    call assert_dealloc_ptr (m, ptr_m)

end subroutine



pure subroutine normalize_2d (x, m, s, dim, center, scale, &
        dof, status)
    !*  NORMALIZE centers and optionally scales input data
    !   such that the resulting array has zero mean and unity standard deviation
    !   along the specified dimension.
    real (PREC), intent(inout), dimension(:,:), contiguous :: x
        !*  Data that should be normalized. Multiple 'variables' are supported
        !   and can be organized either in rows or columns (goverened by
        !   the dim argument)
    real (PREC), intent(out), dimension(:), contiguous, optional :: m, s
        !*  If present, contain the per-variable means and standard deviations
        !   on exit. Arrays need to be at least as large as the number of variables.
    integer, intent(in), optional :: dim
        !!  If present, specifies the dimension along which to normalize (default: dim=1)
    logical, intent(in), optional :: center
        !!  If present and true, center input data around its mean.
    logical, intent(in), optional :: scale
        !!  If present and true, divide each variable by its standard deviation
    integer, intent(in), optional :: dof
        !*  If present, specifies the degrees-of-freedom correction applied
        !   when computing the standard deviation. Default: 1
    type (status_t), intent(out), optional :: status
        !*  If present, returns 0 on exit if normalization was successful,
        !   and status > 0 if an error was encountered

    real (PREC), dimension(:), allocatable :: lmean, lstd
    integer :: ldim, ldof
    type (status_t) :: lstatus
    ! iv indexes variables, iobs individual observations
    integer :: iv, iobs, nvars, nobs
    logical :: lscale, lcenter
    integer, dimension(2) :: shp

    lcenter = .true.
    lscale = .true.
    ldof = 1
    if (present(center)) lcenter = center
    if (present(scale)) lscale = scale
    if (present(dof)) ldof = dof

    ! Initialize to supress gfortran warnings
    nvars = -1
    nobs = -1

    ! initialize dim, nvars and nobs from input data
    shp = shape(x)
    call mean_std_init (shp, dim, ldim, nvars, nobs, status=lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    ! Delegate all remaining input checking to STD or MEAN routines

    allocate (lmean(nvars), source=0.0_PREC)
    allocate (lstd(nvars), source=1.0_PREC)

    ! Terminate immediately if there's nothing to do
    if (.not. (lscale .or. lcenter)) goto 50

    ! Compute mean (and std) as required.
    ! If only centering is requested (no scaling), compute only
    ! mean and set std = 1.0
    if (lscale) then
        call std (x, s=lstd, m=lmean, dim=ldim, dof=ldof, status=lstatus)
        if (.not. (NF_STATUS_OK .in. lstatus)) goto 100
    else if (lcenter) then
        call mean (x, m=lmean, dim=ldim, status=lstatus)
        if (.not. (NF_STATUS_OK .in. lstatus)) goto 100
    end if

    ! Reset to zero mean correction if case lscale = .true. but lcenter = .false.
    if (.not. lcenter) then
        lmean = 0.0_PREC
    end if

    ! Apply normalization in-place
    if (ldim == 1) then
        ! Compute running mean and std if variables are stored in columns
        do iv = 1, nvars
            ! Apply normalization in-place
            do iobs = 1, nobs
                x(iobs, iv) = (x(iobs, iv) - lmean(iv)) / lstd(iv)
            end do
        end do
    else if (ldim == 2) then
        do iobs = 1, nobs
            do iv = 1, nvars
                x(iv, iobs) = (x(iv, iobs) - lmean(iv)) / lstd(iv)
            end do
        end do
    end if

50  continue

    ! Return mean and std. only if 1) output arrays are present; and 2)
    ! centering / scaling was actually performed.
    if (present(m) .and. lcenter) m = lmean
    if (present(s) .and. lscale) s = lstd

100 continue
    if (present(status)) status = lstatus

end subroutine



pure subroutine standardize_1d (x, dim, center, scale, skip_const, dof, mean_x, &
        std_x, shift_x, scale_x, status)
    !*  STANDARDIZE centers and scales input data such that the resulting
    !   array has zero mean and unit standard deviation
    real (PREC), intent(inout), dimension(:), contiguous :: x
        !*  Data that should be standardized.
    integer, intent(in), optional :: dim
        !*  Ignored for 1d array, only present for API compatibility.
    logical, intent(in), optional :: center
        !*  If present and true, center input data around its mean
        !   (default: true)
    logical, intent(in), optional :: scale
        !*  If true, divide each variable by its standard deviation
        !   (default: true)
    logical, intent(in), optional :: skip_const
        !*  If true, do not scale constant "variables" in X to avoid
        !   division by zero (default: true)
    integer, intent(in), optional :: dof
        !*  Degrees of freedom correction for computing std. deviation
        !   (default: 1)
    real (PREC), intent(out), optional :: mean_x
        !*  If present, stores the mean of X.
    real (PREC), intent(out), optional :: std_x
        !*  If present, stores the std. dev. of X.
    real (PREC), intent(out), optional :: shift_x
        !*  If present, stores the quantity by which the elements in X
        !   were shifted (for center=.true., SHIFT_X contains the mean of X,
        !   and otherwise is zero)
    real (PREC), intent(out), optional :: scale_x
        !*  If present, stores the quantity by which the elements in X
        !   were scaled (for scale=.true., SCALE_X contains the std. dev. of X,
        !   and otherwise is one).
    type (status_t), intent(out), optional :: status
        !*  Exit code

    integer :: ldof
    type (status_t) :: lstatus
    ! iv indexes variables, iobs individual observations
    integer :: nobs
    logical :: lscale, lcenter, lskip_const, need_std, need_mean_only
    real (PREC) :: lmean_x, lstd_x
    character (*), parameter :: NAME = 'STANDARDIZE'

    lstatus = NF_STATUS_OK
    nobs = size(x)

    ! --- Quick termination ---

    if (nobs == 0) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    ! --- Input checks ---

    call set_optional_arg (center, .true., lcenter)
    call set_optional_arg (scale, .true., lscale)
    call set_optional_arg (skip_const, .false., lskip_const)
    call set_optional_arg (dof, 1, ldof)

    call check_cond (ldof == 0 .or. ldof == 1, NAME, 'DOF: invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Implementation ---

    lmean_x = 0.0
    lstd_x = 1.0

    need_std = lscale .or. present(std_x)
    need_mean_only = (lcenter .or. present(mean_x)) .and. .not. need_std

    if (need_std) then
        call std_impl (x, dof=ldof, s=lstd_x, m=lmean_x)
    else if (need_mean_only) then
        call mean_impl (x, m=lmean_x)
    end if

    if (lskip_const .and. lscale) then
        ! We use modified STD_X only for the impl. routine, the
        ! user will get the correct STD_X if the argument is present,
        ! the SCALE_X value will be correctly set to 1.0 for constant
        ! variable.
        call standardize_impl (x, lmean_x, std_x=1.0_PREC, dim=1, center=lcenter, &
            scale=lscale, shift_x=shift_x, scale_x=scale_x)
    else
        call standardize_impl (x, lmean_x, lstd_x, dim=1, center=lcenter, &
            scale=lscale, shift_x=shift_x, scale_x=scale_x)
    end if


    if (present(mean_x)) mean_x = lmean_x
    if (present(std_x)) std_x = lstd_x

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine standardize_2d (x, dim, center, scale, skip_const, dof, mean_x, &
        std_x, shift_x, scale_x, status)
    !*  STANDARDIZE centers and scales input data such that the resulting
    !   array has zero mean and unit standard deviation along the specified
    !   dimension.
    real (PREC), intent(inout), dimension(:,:), contiguous :: x
        !*  Data that should be standardized. Multiple 'variables' are supported
        !   and can be organized either in rows or columns (goverened by
        !   the dim argument)
    integer, intent(in), optional :: dim
        !*  Dimension along which to normalize (default: dim=1)
    logical, intent(in), optional :: center
        !*  If present and true, center input data around its mean
        !   (default: true)
    logical, intent(in), optional :: scale
        !*  If true, divide each variable by its standard deviation
        !   (default: true)
    logical, intent(in), optional :: skip_const
        !*  If true, do not scale constant "variables" in X to avoid
        !   division by zero (default: true)
    integer, intent(in), optional :: dof
        !*  Degrees of freedom correction for computing std. deviation
        !   (default: 1)
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: mean_x
        !*  Array of means of variables in X over dimension DIM.
    real (PREC), intent(out), dimension(:), contiguous, optional, target :: std_x
        !*  Array of standard deviations of variables in X over dimension DIM.
    real (PREC), intent(out), dimension(:), contiguous, optional :: shift_x
        !*  If present, stores the quantities by which each variable in
        !   X was shifted (for center=.true., SHIFT_X contains the means of X,
        !   otherwise zeros).
    real (PREC), intent(out), dimension(:), contiguous, optional :: scale_x
        !*  If present, stores the quantities by which each variable in
        !   X was scaled (for scale=.true., SCALE_X contains the std. dev. of X,
        !   otherwise ones).
    type (status_t), intent(out), optional :: status
        !*  Exit code

    integer :: ldim, ldof
    type (status_t) :: lstatus
    ! iv indexes variables, iobs individual observations
    integer :: nvars, nobs
    logical :: lscale, lcenter, lskip_const, need_std, need_mean_only
    real (PREC), dimension(:), pointer, contiguous :: ptr_mean_x, ptr_std_x
    real (PREC), dimension(:), allocatable :: lstd_x
    character (*), parameter :: NAME = 'STANDARDIZE'

    lstatus = NF_STATUS_OK
    nullify (ptr_mean_x, ptr_std_x)

    ! --- Input checks ---

    call set_optional_arg (dim, 1, ldim)

    ! Get data dimensions; this also checks validity of DIM argument
    call get_data_dims (x, dim, nvars, nobs, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_cond (nobs > 0, NAME, 'X: Too few observations', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! --- Quick termination ---

    ! STATUS should remain at its value NF_STATUS_OK!
    if (nvars == 0) goto 100

    ! --- Remaining input checks

    call set_optional_arg (center, .true., lcenter)
    call set_optional_arg (scale, .true., lscale)
    call set_optional_arg (skip_const, .false., lskip_const)
    call set_optional_arg (dof, 1, ldof)

    call check_cond (ldof == 0 .or. ldof == 1, NAME, 'DOF: invalid value', lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(mean_x)) then
        call check_cond (size(mean_x) == nvars, NAME, &
            'MEAN_X: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(std_x)) then
        call check_cond (size(std_x) == nvars, NAME, &
            'STD_X: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(shift_x)) then
        call check_cond (size(shift_x) == nvars, NAME, &
            'SHIFT_X: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    if (present(scale_x)) then
        call check_cond (size(scale_x) == nvars, NAME, &
            'SCALE_X: Non-conformable array size', lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end if

    ! --- Compute moments ---

    need_std = lscale .or. present(std_x)
    need_mean_only = (lcenter .or. present(mean_x)) .and. .not. need_std

    if (need_std) then
        call assert_alloc_ptr (mean_x, nvars, ptr_mean_x)
        call assert_alloc_ptr (std_x, nvars, ptr_std_x)
        call std_impl (x, dim=ldim, dof=ldof, s=ptr_std_x(1:nvars), m=ptr_mean_x(1:nvars))
    else if (need_mean_only) then
        call assert_alloc_ptr (mean_x, nvars, ptr_mean_x)
        call mean_impl (x, dim=ldim, m=ptr_mean_x(1:nvars))
    end if

    ! --- Implementation ---

    if (lskip_const .and. lscale) then
        ! Adjust STD_X so that constant "variables" in X will not be
        ! modified
        allocate (lstd_x(size(ptr_std_x)), source=ptr_std_x)
        where (lstd_x == 0.0_PREC)
            lstd_x = 1.0
        end where

        ! We use modified STD_X array only for the impl. routine, the
        ! user will get the correct STD_X if the argument is present,
        ! the SCALE_X values will be correctly set to 1.0 for constant
        ! variables.
        call standardize_impl (x, ptr_mean_x, lstd_x, ldim, lcenter, lscale, &
            shift_x, scale_x)
    else
        call standardize_impl (x, ptr_mean_x, ptr_std_x, ldim, lcenter, lscale, &
            shift_x, scale_x)
    end if

    call assert_dealloc_ptr (mean_x, ptr_mean_x)
    call assert_dealloc_ptr (std_x, ptr_std_x)

100 continue

    if (present(status)) status = lstatus

end subroutine



pure subroutine standardize_impl_1d (x, mean_x, std_x, dim, center, scale, &
        shift_x, scale_x)
    !*  STANDARDIZE_IMPL centers and scales input data
    !   such that the resulting array has zero mean and unit standard deviation
    !   along the specified dimension.
    !
    !   Implementation routine, does NOT perform any input checks.
    real (PREC), intent(inout), dimension(:), contiguous :: x
        !*  Data that should be standardized. Multiple 'variables' are supported
        !   and can be organized either in rows or columns (goverened by
        !   the dim argument)
    real (PREC), intent(in), optional :: mean_x
        !*  Array of means of variables in X over dimension DIM.
    real (PREC), intent(in), optional :: std_x
        !*  Array of standard deviations of variables in X over dimension DIM.
    integer, intent(in) :: dim
        !*  Ignored for 1d arrays, present only for API compatibility
    logical, intent(in) :: center
        !*  If true, center input data around its mean.
    logical, intent(in) :: scale
        !*  If true, divide each variable by its standard deviation
    real (PREC), intent(out), optional :: shift_x
        !*  If present, stores the quantity by which elements in
        !   X were shifted (for center=.true., SHIFT_X contains the mean of X,
        !   and otherwise is zero).
    real (PREC), intent(out), optional :: scale_x
        !*  If present, stores the quantity by which elements in
        !   X were scaled (for scale=.true., SCALE_X contains the std. dev. of X,
        !   and otherwise is ones).

    integer :: nobs

    nobs = size(x)

    ! --- Quick termination ---

    if (nobs == 0) return

    ! --- Implementation ---

    if (center .and. present(mean_x)) then
        x = x - mean_x
        if (present(shift_x)) then
            shift_x = mean_x
        end if
    else
        if (present(shift_x)) then
            shift_x = 0.0_PREC
        end if
    end if

    if (scale .and. present(std_x)) then
        x = x / std_x
        if (present(scale_x)) then
            scale_x = std_x
        end if
    else
        if (present(scale_x)) then
            scale_x = 1.0_PREC
        end if
    end if

end subroutine



pure subroutine standardize_impl_2d (x, mean_x, std_x, dim, center, scale, &
        shift_x, scale_x)
    !*  STANDARDIZE_IMPL centers and scales input data
    !   such that the resulting array has zero mean and unit standard deviation
    !   along the specified dimension.
    !
    !   Implementation routine, does NOT perform any input checks.
    real (PREC), intent(inout), dimension(:,:), contiguous :: x
        !*  Data that should be standardized. Multiple 'variables' are supported
        !   and can be organized either in rows or columns (goverened by
        !   the dim argument)
    real (PREC), intent(in), dimension(:), contiguous, optional :: mean_x
        !*  Array of means of variables in X over dimension DIM.
    real (PREC), intent(in), dimension(:), contiguous, optional :: std_x
        !*  Array of standard deviations of variables in X over dimension DIM.
    integer, intent(in) :: dim
        !*  Dimension along which to normalize.
    logical, intent(in) :: center
        !*  If true, center input data around its mean.
    logical, intent(in) :: scale
        !*  If true, divide each variable by its standard deviation
    real (PREC), intent(out), dimension(:), contiguous, optional :: shift_x
        !*  If present, stores the quantities by which each variable in
        !   X was shifted (for center=.true., SHIFT_X contains the means of X,
        !   otherwise zeros).
    real (PREC), intent(out), dimension(:), contiguous, optional :: scale_x
        !*  If present, stores the quantities by which each variable in
        !   X was scaled (for scale=.true., SCALE_X contains the std. dev. of X,
        !   otherwise ones).

    ! iv indexes variables, iobs individual observations
    integer :: iv, iobs, nvars, nobs

    ! Initialize to supress gfortran warnings
    nvars = size(x, 3-dim)
    nobs = size(x, dim)

    ! --- Quick termination ---

    if (nobs == 0 .or. nvars == 0) return

    ! --- Implementation ---

    ! We don't check whether these arrays are possibily too large, so
    ! erase any "excess" data that will not be overwritten below.
    ! These are also the default values if no centering or scaling is performed.
    if (present(shift_x)) then
        shift_x = 0.0_PREC
    end if

    if (present(scale_x)) then
        scale_x = 1.0_PREC
    end if

    if (center .and. present(mean_x)) then
        if (dim == 1) then
            do iv = 1, nvars
                x(:,iv) = x(:,iv) - mean_x(iv)
            end do
        else
            do iobs = 1, nobs
                x(:,iobs) = x(:,iobs) - mean_x
            end do
        end if

        if (present(shift_x)) then
            shift_x(1:nvars) = mean_x(1:nvars)
        end if
    end if

    if (scale .and. present(std_x)) then
        if (dim == 1) then
            do iv = 1, nvars
                x(:,iv) = x(:,iv) / std_x(iv)
            end do
        else
            do iobs = 1, nobs
                x(:,iobs) = x(:,iobs) / std_x
            end do
        end if

        if (present(scale_x)) then
            scale_x(1:nvars) = std_x(1:nvars)
        end if
    end if

end subroutine


! ------------------------------------------------------------------------------
! QUANTILE ROUTINES

pure subroutine quantile_bins_check_input (x, pmf, q, pctl, &
        interp, status)
    !*  QUANTILE_BINS_CHECK_INPUT performs input validation for the
    !   QUANTILE_BINS routine.
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: pmf
    real (PREC), intent(in), dimension(:) :: q
    real (PREC), intent(in), dimension(:) :: pctl
    character (*), intent(in), optional :: interp
    type (status_t), intent(out) :: status

    integer (NF_ENUM_KIND) :: imode
    status = NF_STATUS_INVALID_ARG

    if (present(interp)) then
        ! The following routine returns 0 if INTERP does not match any of the
        ! acceptable values
        imode = quantile_interp_to_enum (interp)
        if (imode == 0) goto 100
    end if

    if ((size(x) - 1) /= size(pmf)) goto 100
    if (size(q) /= size(pctl)) goto 100

    ! Enforce percentiles to be in valid range
    if (any(q < 0.0_PREC)) goto 100
    if (any(q > 1.0_PREC)) goto 100

    status = NF_STATUS_OK

100 continue

end subroutine



pure subroutine quantile_discrete_check_input (x, pmf, rnk, q, status)
    !* Validate discrete support, weights, ranks, and output size.

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    real (PREC), intent(in), dimension(:) :: x
        !* Finite support values.
    real (PREC), intent(in), dimension(:) :: pmf
        !* Finite, nonnegative weights with at least one positive value.
    real (PREC), intent(in), dimension(:) :: rnk
        !* Finite ranks in the closed interval [0,1].
    real (PREC), intent(in), dimension(:) :: q
        !* Output array whose size must match the ranks.
    type (status_t), intent(out) :: status
        !* Validation status.

    status = NF_STATUS_INVALID_ARG

    if (size(x) == 0) goto 100
    if (size(x) /= size(pmf)) goto 100
    if (size(rnk) /= size(q)) goto 100

    if (.not. all (ieee_is_finite (x))) goto 100
    if (.not. all (ieee_is_finite (pmf))) goto 100
    if (any(pmf < 0.0_PREC)) goto 100
    if (.not. any(pmf > 0.0_PREC)) goto 100

    if (.not. all (ieee_is_finite (rnk))) goto 100
    if (any(rnk < 0.0_PREC)) goto 100
    if (any(rnk > 1.0_PREC)) goto 100

    status = NF_STATUS_OK

100 continue

end subroutine



pure subroutine quantile_bins (x, pmf, rnk, q, interp, status)
    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Array of bin edges (or bin midpoints) in increasing order
    real (PREC), intent(in), dimension(:), contiguous :: pmf
        !*  PMF associated with bins defined by X
    real (PREC), intent(in), dimension(:), contiguous :: rnk
        !*  Array of quantiles "ranks" to compute which must be between
        !   0 and 1 inclusive.
    real (PREC), intent(out), dimension(:), contiguous :: q
        !*  Output array storing (interpolated) quantiles corresponding to RNK
    character (*), intent(in), optional :: interp
        !*  Interpolation method to use when desired percentile is between
        !   two data points with indices i and i+1. Valid values are
        !      1.  'linear': linearly interpolate between values at i, i+1
        !      2.  'lower': x(i)
        !      3.  'higher': x(i+1)
        !      4.  'nearest': value at index which is nearest
        !      5.  'midpoint': (x(i)+x(i+1))/2
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    type (status_t) :: lstatus
    integer (NF_ENUM_KIND) :: imode
    integer :: n, nq, i, imax, ncdf, ilb
    real (PREC), dimension(:), allocatable :: cdf
    real (PREC) :: wgt, ri
    type (search_cache) :: cache

    lstatus = NF_STATUS_OK

    call quantile_bins_check_input (x, pmf, rnk, q, interp, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Convert to numeric interpolation enum
    imode = NF_STATS_QUANTILE_LINEAR
    if (present(interp)) then
        imode = quantile_interp_to_enum (interp)
    end if

    n = size(pmf)
    nq = size(rnk)

    ncdf = n + 1
    allocate (cdf(ncdf))
    cdf(1) = 0.0

    do i = 1, n
        cdf(i+1) = cdf(i) + pmf(i)
    end do
    cdf(:) = cdf / cdf(ncdf)

    ! Find last bin with non-zero mass
    do imax = n, 2, -1
        if (pmf(imax) > 0.0_PREC) exit
    end do

    ! CDF has one additional element, so increment max. index
    imax = imax + 1

    select case (imode)
    case (NF_STATS_QUANTILE_LINEAR)

        call interp_linear (rnk, cdf(1:imax), x(1:imax), q, &
            ext=NF_INTERP_EVAL_BOUNDARY)

    case (NF_STATS_QUANTILE_NEAREST)

        do i = 1, nq
            call interp_find (rnk(i), cdf(1:imax), ilb, wgt, cache=cache)
            if (wgt >= 0.50_PREC) then
                q(i) = x(ilb)
            else
                q(i) = x(ilb+1)
            end if
        end do

    case (NF_STATS_QUANTILE_LOWER)

        do i = 1, nq
            ri = rnk(i)
            if (ri == 1.0_PREC) then
                ! Take the max., ignore "lower" argument
                q(i) = x(imax)
            else
                call bsearch_cached (ri, cdf(1:imax), ilb, cache)
                q(i) = x(ilb)
            end if
        end do

    case (NF_STATS_QUANTILE_HIGHER)

        do i = 1, nq
            ri = rnk(i)
            call bsearch_cached (ri, cdf(1:imax), ilb, cache)
            if (ri > 0.0_PREC) then
                q(i) = x(ilb+1)
            else
                ! For 0.0 rank (ie. the min), take the lower bound
                ! Note: we need to call BSEARCH even in this case as the
                ! PMF could contain zero values in the beginning
                q(i) = x(ilb)
            end if
        end do

    case (NF_STATS_QUANTILE_MIDPOINT)

        wgt = 0.5_PREC
        do i = 1, nq
            ri = rnk(i)
            call bsearch_cached (ri, cdf(1:imax), ilb, cache)
            q(i) = wgt * x(ilb) + (1.0_PREC-wgt) * x(ilb+1)
        end do
    end select

100 continue
    if (allocated(cdf)) deallocate (cdf)

    if (present(status)) status = lstatus

end subroutine



subroutine quantile_discrete (x, pmf, rnk, q, sort, status)
    !* Compute discrete quantiles using positive-mass support.
    !  For 0 <= r < 1, select the smallest support value with F(x) > r.
    !  At r = 1, select the largest support value with positive mass.

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    real (PREC), intent(in), dimension(:), contiguous :: x
        !* Nonempty array of finite support values.
    real (PREC), intent(in), dimension(:), contiguous :: pmf
        !* Finite, nonnegative weights associated with X. Their total must
        !  be positive and finite; weights are normalized internally.
    real (PREC), intent(in), dimension(:), contiguous :: rnk
        !*  Array of quantiles "ranks" to compute which must be between
        !   0 and 1 inclusive.
    real (PREC), intent(out), dimension(:), contiguous :: q
        !* Quantiles corresponding to RNK, each in the positive-mass support.
    logical, intent(in), optional :: sort
        !* Sort X internally by default. If false, X must be nondecreasing.
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    type (status_t) :: lstatus
    integer :: n, nq, i, j, imax, ncdf, ilb
    real (PREC), dimension(:), allocatable :: cdf
    real (PREC) :: ri
    type (search_cache) :: cache
    logical :: lsort
    integer, dimension(:), allocatable :: idx_sort

    lstatus = NF_STATUS_OK
    lsort = .true.

    if (present(sort)) lsort = sort

    call quantile_discrete_check_input (x, pmf, rnk, q, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    n = size(pmf)
    nq = size(rnk)

    if (.not. lsort) then
        if (any(x(2:n) < x(1:n-1))) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    ncdf = n + 1
    allocate (cdf(ncdf))
    cdf(1) = 0.0_PREC
    imax = 0

    ! Track the last positive weight in the same order used to build the CDF.
    if (lsort) then
        allocate (idx_sort(n))
        call argsort (x, idx_sort)

        do i = 1, n
            j = idx_sort(i)
            cdf(i+1) = cdf(i) + pmf(j)
            if (pmf(j) > 0.0_PREC) imax = i
        end do
    else
        do i = 1, n
            cdf(i+1) = cdf(i) + pmf(i)
            if (pmf(i) > 0.0_PREC) imax = i
        end do
    end if

    if (.not. ieee_is_finite (cdf(ncdf))) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if
    cdf(:) = cdf / cdf(ncdf)

    do i = 1, nq
        ri = rnk(i)
        call bsearch_cached (ri, cdf(1:imax+1), ilb, cache)
        if (lsort) then
            j = idx_sort(ilb)
        else
            j = ilb
        end if
        q(i) = x(j)
    end do

100 continue
    if (allocated(cdf)) deallocate (cdf)
    if (allocated(idx_sort)) deallocate (idx_sort)

    if (present(status)) status = lstatus

end subroutine



subroutine quantile_dispatch (x, pmf, rnk, q, interp, sort, status)
    !* Compute quantiles from discrete support or bin edges.
    !  Equal-length X and PMF arrays select discrete quantiles; INTERP is
    !  ignored. For 0 <= r < 1, return the smallest positive-mass support
    !  value with F(x) > r, so an exact interior CDF jump selects the next
    !  positive-mass value. Ranks 0 and 1 select the minimum and maximum
    !  positive-mass support values, respectively.
    !  With one more element in X than PMF, use bin-edge interpolation.

    real (PREC), intent(in), dimension(:), contiguous :: x
        !* Discrete support values, or bin edges in increasing order.
    real (PREC), intent(in), dimension(:), contiguous :: pmf
        !* Weights associated with X. Discrete distributions require nonempty,
        !  finite support and finite, nonnegative weights with a positive,
        !  finite total. Weights are normalized internally.
    real (PREC), intent(in), dimension(:), contiguous :: rnk
        !*  Percentiles to compute which must be between 0.0 and 1.0 inclusive.
    real (PREC), intent(out), dimension(:), contiguous :: q
        !*  (Interpolated) quantiles corresponding to RNK
    character (*), intent(in), optional :: interp
        !*  Bin interpolation method; ignored for discrete support.
        !   Interpolation to use when desired percentile is between
        !   two data points with indices i and i+1. Valid values are
        !      1.  "linear": linearly interpolate between values at i, i+1
        !      2.  "lower": x(i)
        !      3.  "higher": x(i+1)
        !      4.  "nearest": value at index which is nearest
        !      5.  "midpoint": (x(i)+x(i+1))/2
    logical, intent(in), optional :: sort
        !* Sort discrete support internally by default. If false, discrete
        !  support must be nondecreasing. Ignored for bin edges.
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    if (size(x) == size(pmf)) then
        call quantile_discrete (x, pmf, rnk, q, sort, status)
    else
        call quantile_bins (x, pmf, rnk, q, interp, status)
    end if

end subroutine



subroutine quantile_dispatch_scalar (x, pmf, rnk, q, interp, sort, status)
    !* Scalar version of QUANTILE_DISPATCH with the same input requirements,
    !  discrete CDF jump convention, and bin-edge interpolation methods.

    real (PREC), intent(in), dimension(:), contiguous :: x
        !* Discrete support values, or bin edges in increasing order.
    real (PREC), intent(in), dimension(:), contiguous :: pmf
        !* Weights associated with the support values or bins defined by X.
    real (PREC), intent(in) :: rnk
        !*  Percentile to compute which must be between 0.0 and 1.0 inclusive.
    real (PREC), intent(out) :: q
        !*  (Interpolated) quantile corresponding to RNK
    character (*), intent(in), optional :: interp
        !*  Bin interpolation method; ignored for discrete support.
        !   Interpolation to use when desired percentile is between
        !   two data points with indices i and i+1. Valid values are
        !      1.  "linear": linearly interpolate between values at i, i+1
        !      2.  "lower": x(i)
        !      3.  "higher": x(i+1)
        !      4.  "nearest": value at index which is nearest
        !      5.  "midpoint": (x(i)+x(i+1))/2
    logical, intent(in), optional :: sort
        !* Sort discrete support internally by default. If false, discrete
        !  support must be nondecreasing. Ignored for bin edges.
    type (status_t), intent(out), optional :: status
        !*  Optional status code

    real (PREC) :: q1(1), rnk1(1)

    rnk1(1) = rnk
    call quantile (x, pmf, rnk1, q1, interp, sort, status)
    q = q1(1)

end subroutine



subroutine cov_check_input (x, cov, dim, dof, status)
    real (PREC), intent(in), dimension(:,:) :: x
    real (PREC), intent(in), dimension(:,:) :: cov
    integer, intent(in), optional :: dim
    integer, intent(in), optional :: dof
    type (status_t), intent(out) :: status

    integer :: shp(2), nvars, nobs, ldim

    ldim = 1
    if (present(dim)) ldim = dim

    if (ldim < 1 .or. ldim > 2) goto 100
    if (present(dof)) then
        if (dof < 0 .or. dof  > 1) goto 100
    end if

    shp = shape(x)
    nvars = shp(3-ldim)
    nobs = shp(ldim)

    if (size(cov,1) /= nvars) goto 100
    if (size(cov,1) /= size(cov,2)) goto 100

    status = NF_STATUS_OK
    return

100 continue

    status = NF_STATUS_INVALID_ARG

end subroutine



subroutine cov (x, vcv, dim, dof, status)
    !*  COV estimates the variance-covariance matrix for a given set of
    !   variables.
    real (PREC), intent(in), dimension(:,:), contiguous :: x
        !*  Array of input data. If DIM=1 (the default), it is assumed that
        !   each column contains a variable, while for DIM=2 each column
        !   contains one (multivariate) observation.
    real (PREC), intent(out), dimension(:,:), contiguous :: vcv
        !*  Estimated variance-covariance matrix
    integer, intent(in), optional :: dim
        !*  Dimension along which covariances should be computed (default: dim=1
        !   which implies that each column corresponds to one variable)
    integer, intent(in), optional :: dof
        !*  Degrees-of-freedom correction to apply (default: dof=1)
    type (status_t), intent(out), optional :: status
        !*  Optional exit status

    integer :: ldim, ldof, nvars, nobs
    type (status_t) :: lstatus
    real (PREC) :: alpha, beta
    integer :: shp(2)
    real (PREC), dimension(:), allocatable :: m

    ldim = 1
    ldof = 1
    if (present(dim)) ldim = dim
    if (present(dof)) ldof = dof

    lstatus = NF_STATUS_OK

    call cov_check_input (x, vcv, dim, dof, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    shp = shape(x)
    nvars = shp(3-ldim)
    nobs = shp(ldim)

    vcv = 0.0
    ! Covariance is zero in degenerate cases
    if (nobs < 2) goto 100

    ! Compute vector of means
    allocate (m(nvars))
    m(:) = sum(x, dim=ldim)
    m(:) = m / nobs

    ! Compute outer product m x m' = E[X] * E[X']
    call BLAS_GER (m=nvars, n=nvars, alpha=1.0_PREC, x=m, incx=1, y=m, &
        incy=1, a=vcv, lda=nvars)

    ! Note on DOF correction: Assume we are computing the sample
    ! Covariance as
    !   Cov(x,y) = 1/(n-1) * sum((x_i-mu_x)(y_i-mu_y))
    ! Then this equals
    !   Cov(x,y) = 1/(n-1) * sum(x_i * y_i) - n/(n-1) * mu_x * mu_y
    ! thus we have to rescale the product of means by n/(n-1).
    alpha = 1.0_PREC / (nobs - ldof)
    ! Scaling factor on product of means, as computed above
    beta = - real(nobs, PREC) / (nobs - ldof)

    ! Compute covariance E[XX'] - E[X]E[X']
    if (ldim == 1) then
        ! X has shape [NOBS, NVARS], need to compute X'X
        call BLAS_GEMM (transa='T', transb='N', m=nvars, n=nvars, k=nobs, &
            alpha=alpha, a=x, lda=nobs, b=x, ldb=nobs, beta=beta, &
            c=vcv, ldc=nvars)
    else
        ! X has shape [NVARS, NOBS], need to compute XX'
        call BLAS_GEMM (transa='N', transb='T', m=nvars, n=nvars, k=nobs, &
            alpha=alpha, a=x, lda=nvars, b=x, ldb=nvars, beta=beta, &
            c=vcv, ldc=nvars)
    end if

100 continue

    if (present(status)) status = lstatus

end subroutine



subroutine corrcoef (x, corr, dim, dof, status)
    !*  CORRCOEF estimates the Pearson's correlation coefficient for a given
    !   set of vaiables.
    real (PREC), intent(in), dimension(:,:), contiguous :: x
        !*  Array of input data. If DIM=1 (the default), it is assumed that
        !   each column contains a variable, while for DIM=2 each column
        !   contains one (multivariate) observation.
    real (PREC), intent(out), dimension(:,:), contiguous :: corr
        !*  Estimated variance-covariance matrix
    integer, intent(in), optional :: dim
        !*  Dimension along which covariances should be computed (default: dim=1
        !   which implies that each column corresponds to one variable)
    integer, intent(in), optional :: dof
        !*  Degrees-of-freedom correction to apply (default: dof=1)
    type (status_t), intent(out), optional :: status
        !*  Optional exit status

    type (status_t) :: lstatus
    real (PREC) :: s
    integer :: i, ldim, shp(2), nvars

    ! Compute VCV matrix (input checks are also performed by COV routine
    call cov (x, corr, dim, dof, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ldim = 1
    if (present(dim)) ldim = dim

    shp = shape(x)
    nvars = shp(3-ldim)

    do i = 1, nvars
        s = sqrt(corr(i,i))
        corr(:,i) = corr(:,i) / s
        corr(i,:) = corr(i,:) / s
        ! Make sure diagonal elements are identically 1.0
        corr(i,i) = 1.0
    end do

100 continue
    if (present(status)) status = lstatus

end subroutine
