C  Important Notice:
C   These algorithms are modifications and based on the Software Newuoa, authored by M. J. D. Powell,
C   to minimize sum of squares by taking advantage of the problem structure.
C    min_{x \in R^n }  F(x) := Sum_{i=1}^{mv}  v_err_i(x)^2
C   where v_err(x) : R^n \to R^{mv} is a vector function.
C   This subroutine seeks the least value of sum of the squres of the components of v_err(x)
C   by combing trust region method and Levenberg-Marquardt method 
C
C   References:
C
C   1.  M. J. D. Powell, The NEWUOA software for unconstrained optimization without derivatives,
C       DAMTP 2004/ NA 05
C   2.  H. Zhang, A. R. CONN, AND K. SCHEINBERG, A DERIVATIVE-FREE ALGORITHM FOR THE LEAST-SQUARE
C       MINIMIZATION, technical report, 2009
C
C      -----------------------------------------------------------------
C      | This program is free software; you can redistribute it and/or  |
C      |modify it under the terms of the GNU General Public License as  |
C      |published by the Free Software Foundation; either version 2 of  |
C      |the License, or (at your option) any later version.             |
C      |This program is distributed in the hope that it will be useful, |
C      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
C      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
C      |GNU General Public License for more details.                    |
C      |                                                                |
C      |You should have received a copy of the GNU General Public       |
C      |License along with this program; if not, write to the Free      |
C      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
C      |MA  02110-1301  USA                                             |
C      -----------------------------------------------------------------|

      SUBROUTINE NEWUOA_H(FCN,N,MV,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,
     &  NF,W,INFO,FX)
      IMPLICIT NONE
      PROCEDURE (fobj_if) :: FCN
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: MV
      INTEGER, INTENT(IN) :: NPT
      REAL (PREC), INTENT(INOUT), DIMENSION(:), CONTIGUOUS :: X
      REAL (PREC), INTENT(IN) :: RHOBEG, RHOEND
      INTEGER, INTENT(IN) :: IPRINT
      INTEGER, INTENT(IN) :: MAXFUN
      INTEGER, INTENT(OUT) :: NF
c       Contains the number of function evaluations on exit
      REAL (PREC), INTENT(OUT), DIMENSION(:), TARGET, CONTIGUOUS :: W
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
C       NEWUOA_CONVERGED, NEWUOA_MAXFUN, NEWUOA_NUMERICAL_ERROR, or
C       NEWUOA_INVALID_INPUT. NF counts all calls, including failed ones.
      REAL (PREC), INTENT(OUT), OPTIONAL :: FX(:)
C       Residual at the returned X, without an additional evaluation.
C       Undefined for invalid input; may be nonfinite if the first call fails.
C
C     N must be set to the number of variables and must be at least two.
C     mv must be set to the lengh of the vector function v_err(x):  R^n \to R^{mv}. 
C
C     NPT is the number of interpolation conditions. Its value must be in the
C     interval [N+2,2N+1].  Recommended: NPT = 2*N+1
C
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C     will be changed to the values that give the least calculated F(x) = Sum_{i=1}^{mv} v_err_i(x)^2.
C
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C     region radius, so both must be positive with RHOEND<=RHOBEG. Typically
C     RHOBEG should be about one tenth of the greatest expected change to a
C     variable, and RHOEND should indicate the accuracy that is required in
C     the final values of the variables. Default: RHOBEG = 1.0, RHOEND = 10^{-8} 
C
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C     amount of printing. Specifically, there is no output if IPRINT=0 and
C     there is output only at the return if IPRINT=1. Otherwise, each new
C     value of RHO is printed, with the best vector of variables so far and
C     the corresponding value of the objective function. Further, each new
C     value of F with its variables are output if IPRINT=3.
C
C     MAXFUN must be set to an upper bound on the number of calls of subroutine
C     dfovec(n, mv, x, v_err) which provides the values of the vector function v_err(x).
C     Here: n, mv, x \in R^n are input, v_err \in R^{mv} are output.
C     Default:  MAXFUN= 100(n+1), i.e 100 (simplex) gradients for reasonable accuracy.
C               MAXFUN= infinity, to let the algorithm explore the lowest function value  
C                       as much as it could.
C
C     The array W will be used for working space. Its length must be at least
C     (NPT+11)*(NPT+N)+N*(3*N+11)/2 + 3*MV*N + MV*(N+1)*N/2 + MV*NPT + 7*MV + N
C
C     SUBROUTINE dfovec(n, mv, x, v_err) must be provided by the user. It must provide
C     the values of the vector function v_err(x) : R^n to R^{mv} at the variables X(1),X(2),...,X(N).
C

      INTEGER :: NP, NPTM, NDIM, IERR
      INTEGER :: IXB, IXO, IXN, IXP, IGQ, IHQ, IPQ
      INTEGER :: IBMAT, IZMAT, ID, IVL, IW
      INTEGER :: IGQV, IHQV, IPQV, IWV, IVERR, IVBEG, IVTEMP, IDIFFV
      INTEGER :: IVBASE, IVOPT, IVQUAD, IHD1, IGQVOPT

      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: XBASE
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: XOPT
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: XNEW

      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: GQ, HQ, PQ
      REAL (PREC), DIMENSION(:,:), POINTER, CONTIGUOUS :: XPT
      REAL (PREC), DIMENSION(:,:), POINTER, CONTIGUOUS :: BMAT
      REAL (PREC), DIMENSION(:,:), POINTER, CONTIGUOUS :: ZMAT

      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: D
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: VLAG
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: V_ERR
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: V_BEG
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: V_TEMP
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: DIFFV
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: V_OPT
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: V_VQUAD
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: V_BASE
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: HD1
      REAL (PREC), DIMENSION(:,:), POINTER, CONTIGUOUS :: GQV
      REAL (PREC), DIMENSION(:,:), POINTER, CONTIGUOUS :: HQV
      REAL (PREC), DIMENSION(:,:), POINTER, CONTIGUOUS :: PQV
      REAL (PREC), DIMENSION(:,:), POINTER, CONTIGUOUS :: WV
      REAL (PREC), DIMENSION(:,:), POINTER, CONTIGUOUS :: GQV_OPT

      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: W_REST

      NF=0
      IERR=NEWUOA_INVALID_INPUT
      IF (PRESENT(INFO)) INFO=IERR
      NP=N+1
      NPTM=NPT-NP
      IF (NPT .LT. N+2 .OR. NPT .GT. 2*N+1) THEN
          PRINT 10
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
     1      '[N+2, 2N+1]')
          GO TO 20
      END IF
      NDIM=NPT+N

      IXB=1
      IXO=IXB+N
      IXN=IXO+N
      IXP=IXN+N
      IGQ=IXP+N*NPT
      IHQ=IGQ+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ID=IZMAT+NPT*NPTM
      IVL=ID+N

      IGQV=IVL+NDIM
      IHQV=IGQV + MV*N
      IPQV=IHQV + MV*(N+1)*N/2
      IWV=IPQV + MV*NPT
      IVERR=IWV + MV*N
      IVBEG=IVERR + MV
      IVTEMP=IVBEG + MV
      IDIFFV=IVTEMP + MV
      IVOPT=IDIFFV + MV
      IVQUAD=IVOPT + MV
      IVBASE=IVQUAD + MV
      IHD1=IVBASE + MV
      IGQVOPT=IHD1 + N
      IW=IGQVOPT + MV*N

      XBASE => W(IXB:IXB+N-1)
      XOPT => W(IXO:IXO+N-1)
      XNEW => W(IXN:IXN+N-1)
      XPT(1:NPT,1:N) => W(IXP:IXP+N*NPT-1)
      GQ => W(IGQ:IGQ+N-1)
      HQ => W(IHQ:IHQ+N*(N+1)/2-1)
      PQ => W(IPQ:IPQ+NPT-1)
      BMAT(1:NDIM,1:N) => W(IBMAT:IBMAT+(NDIM*N)-1)
      ZMAT(1:NPT,1:NPTM) => W(IZMAT:IZMAT+(NPT*NPTM)-1)
      D => W(ID:ID+N-1)
      VLAG => W(IVL:IVL+NDIM-1)

      GQV(1:MV,1:N) => W(IGQV:IGQV+(MV*N)-1)
      HQV(1:MV,1:(N+1)*N/2) => W(IHQV:IHQV+MV*(N+1)*N/2-1)
      PQV(1:MV,1:NPT) => W(IPQV:IPQV+(MV*NPT)-1)
      WV(1:MV,1:N) => W(IWV:IWV+(MV*N)-1)
      GQV_OPT(1:MV,1:N) => W(IGQVOPT:IGQVOPT+(MV*N)-1)
      V_ERR => W(IVERR:IVERR+MV-1)
      V_BEG => W(IVBEG:IVBEG+MV-1)
      V_TEMP => W(IVTEMP:IVTEMP+MV-1)
      DIFFV => W(IDIFFV:IDIFFV+MV-1)
      V_OPT => W(IVOPT:IVOPT+MV-1)
      V_VQUAD => W(IVQUAD:IVQUAD+MV-1)
      V_BASE => W(IVBASE:IVBASE+MV-1)
      HD1 => W(IHD1:IHD1+N-1)
      W_REST => W(IW:)

C
C     The above settings provide a partition of W for subroutine NEWUOB_H.
C
      CALL NEWUOB_H (FCN,N,MV,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,
     &  XBASE,XOPT,XNEW,XPT,GQ,HQ,PQ,BMAT,ZMAT,
     &  NDIM,D,VLAG,GQV,HQV,PQV,WV,V_ERR,V_BEG,V_TEMP,DIFFV,
     &  V_OPT,V_VQUAD,V_BASE,HD1,GQV_OPT,W_REST,NF,IERR)
      IF (PRESENT(INFO)) INFO=IERR
      IF (PRESENT(FX)) FX=V_ERR
   20 RETURN
      END

