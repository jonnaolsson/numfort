
      real (PREC), parameter :: ZERO    = 0.0_PREC
      real (PREC), parameter :: TENTH   = 0.1_PREC
      real (PREC), parameter :: QUART   = 0.25_PREC
      real (PREC), parameter :: HALF    = 0.5_PREC
      real (PREC), parameter :: ONE     = 1.0_PREC
      real (PREC), parameter :: TWO     = 2.0_PREC
      real (PREC), parameter :: TWOPI   = 8.0_PREC * ATAN(1.0_PREC)

C     Termination reasons returned by the numfort interface.
      INTEGER, PARAMETER, PUBLIC :: NEWUOA_CONVERGED = 0
      INTEGER, PARAMETER, PUBLIC :: NEWUOA_MAXFUN = 1
      INTEGER, PARAMETER, PUBLIC :: NEWUOA_NUMERICAL_ERROR = 2
      INTEGER, PARAMETER, PUBLIC :: NEWUOA_INVALID_INPUT = 3

