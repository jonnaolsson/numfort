C
C   Important Notice:
C   This TRSAPP_H are modifications and based on the subroutine TRSAPP in the software NEWUOA, authored by M. J. D. Powell.
C
      SUBROUTINE TRSAPP_H (N,MV,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,STEP,
     1  W,CRVMIN,GQV,HQV,PQV,XBASE,vquad,
     1  GQV_opt,v_opt,v_base,XOPTSQ,model_update,opt_update,INFO)
      IMPLICIT REAL (PREC) (A-H,O-Z)
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: MV
      INTEGER, INTENT(IN) :: NPT
      REAL (PREC), INTENT(IN), DIMENSION(:), CONTIGUOUS :: XOPT
      REAL (PREC), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: XPT
      REAL (PREC), INTENT(INOUT), DIMENSION(:), CONTIGUOUS :: GQ
      REAL (PREC), INTENT(INOUT), DIMENSION(:), CONTIGUOUS :: HQ
      REAL (PREC), INTENT(INOUT), DIMENSION(:), CONTIGUOUS :: PQ
      REAL (PREC), INTENT(IN) :: DELTA
      REAL (PREC), INTENT(OUT), DIMENSION(:), CONTIGUOUS :: STEP
      REAL (PREC), INTENT(OUT), DIMENSION(:), CONTIGUOUS, TARGET :: W
      REAL (PREC), INTENT(INOUT) :: CRVMIN
      REAL (PREC), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: GQV
      REAL (PREC), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: HQV
      REAL (PREC), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: PQV
      REAL (PREC), INTENT(INOUT), DIMENSION(:), CONTIGUOUS :: XBASE
      REAL (PREC), INTENT(INOUT) :: VQUAD
      REAL (PREC), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: GQV_OPT
      REAL (PREC), INTENT(INOUT), DIMENSION(:), CONTIGUOUS :: V_OPT
      REAL (PREC), INTENT(INOUT), DIMENSION(:), CONTIGUOUS :: V_BASE
      REAL (PREC), INTENT(IN) :: XOPTSQ
      LOGICAL, INTENT(INOUT) :: model_update, opt_update
      INTEGER, INTENT(OUT) :: INFO
C       Report a failed trust-region step without terminating the caller.

      logical zero_res, debug
      REAL (PREC) v_gtemp(MV), gbeg(N), gtemp(N)

      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: D
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: G
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: HS
      REAL (PREC), DIMENSION(:), POINTER, CONTIGUOUS :: HD

      D => W(1:N)
      G => W(N+1:N+N)
      HD => W(2*N+1:2*N+N)
      HS => W(3*N+1:3*N+N)

      debug = .false.

      if ((.not.model_update).and.(.not.opt_update)) go to 8
       model_update = .false.
       opt_update = .false.


      if (SQRT(XOPTSQ).gt.0.25_PREC*delta) then
c
c Use the gradient at xopt to formulate J^t J
c
       do m1=1,mv
         do i=1,n
           GQV_opt(m1,i) = GQV(m1,i)
         enddo
         do k=1, npt
           temp = zero
           do j=1,n
             temp = temp+XPT(k,j)*XOPT(j)
           enddo
           temp = temp*PQV(m1,k)
           do i=1,n
             GQV_opt(m1,i)=GQV_opt(m1,i)+temp*XPT(k,i)
           enddo
         enddo
         IH=0
         do j=1,n
           do i=1,j
             IH=IH+1
             if (i .lt. j) GQV_opt(m1,j)=GQV_opt(m1,j)
     &                                   +HQV(m1,IH)*XOPT(i)
                GQV_opt(m1,i)=GQV_opt(m1,i)+HQV(m1,IH)*XOPT(j)
           enddo
         enddo  
       enddo

       call f_grad(mv,v_opt,v_gtemp)
       gnorm2 = zero
       do i=1,n
         GQ(i) = zero
         do m1=1,mv
            GQ(i) = GQ(i) + v_gtemp(m1)*GQV_opt(m1,i)
         enddo
         gnorm2 = gnorm2 + GQ(i)**2
       enddo
c
c Calculate the explicite Hessian.
c
       f_opt = zero
       call f_value(mv,v_opt,f_opt)
       if (gnorm2.ge.1.0_PREC.or.f_opt.le.SQRT(gnorm2)) then
         zero_res = .true.
       else
         zero_res = .false.
       endif

       IH=0
       do j=1,n
         do i=1,j
           IH=IH+1
           if (zero_res) then
             t1 = zero
             do m1=1,mv  
               t1 = t1+GQV_opt(m1,i)*GQV_opt(m1,j)
             enddo
             HQ(IH) = 2.0_PREC*t1
           else
             t1 = zero
             do m1=1,mv
               t2 = zero
               do k=1,npt
                 t2 = t2 + XPT(k,i)*PQV(m1,k)*XPT(k,j)
               enddo
               t2 = t2 + HQV(m1,IH)
               t1 = t1+(GQV_opt(m1,i)*GQV_opt(m1,j)+v_opt(m1)*t2) 
             enddo
             HQ(IH) = 2.0_PREC*t1
           endif
         enddo
       enddo

      else
c
c Use the gradient at xbase to formulate J^t J
c
       call f_grad(mv,v_base,v_gtemp)
       gnorm2 = zero
       do i=1,n
         GQ(i) = zero
         do m1=1,mv
            GQ(i) = GQ(i) + v_gtemp(m1)*GQV(m1,i)
         enddo
         gnorm2 = gnorm2 + GQ(i)**2
       enddo
c
c Calculate the explicite Hessian.
c
       f_base = zero
       call f_value(mv,v_base,f_base)
       if (gnorm2.ge.1.0_PREC.or.f_base.le.SQRT(gnorm2)) then
         zero_res = .true.
       else
         zero_res = .false.
       endif

       IH=0
       do j=1,n
         do i=1,j
           IH=IH+1
           if (zero_res) then
             t1 = zero
             do m1=1,mv           
               t1 = t1+GQV(m1,i)*GQV(m1,j)
             enddo
             HQ(IH) = 2.0_PREC*t1
           else           
             t1 = zero
             do m1=1,mv
               t2 = zero
               do k=1,npt
                 t2 = t2 + XPT(k,i)*PQV(m1,k)*XPT(k,j)
               enddo
               t2 = t2 + HQV(m1,IH)
               t1 = t1+(GQV(m1,i)*GQV(m1,j)+v_base(m1)*t2) 
             enddo
             HQ(IH) = 2.0_PREC*t1
           endif
         enddo
       enddo
c calculte the gradient at xopt
       IH=0
       do j=1,n
         do i=1,j
           IH=IH+1
           if (i .lt. j) GQ(j)=GQ(j)+HQ(IH)*XOPT(i)
              GQ(i)=GQ(i)+HQ(IH)*XOPT(j)
         enddo
       enddo  

      endif

    8 DELSQ=DELTA*DELTA
      ITERC=0
      ITERMAX=N
      ITERSW=ITERMAX
      if (debug) then
        t = zero
        do i=1,n
          t = t + xopt(i)**2
        enddo
        print *, " ||xopt||=",SQRT(t)
      endif
      gnorm2 = zero
      DO 10 I=1,N
      gnorm2 = gnorm2 + GQ(i)**2
   10 D(I)=zero 
      gnorm2 = SQRT(gnorm2)
      if (debug) print *, " gnorm2=", gnorm2 
      GOTO 170
C
C     Prepare for the first line search.
C
   20 continue
      QRED=ZERO
      DD=ZERO
      DO 30 I=1,N
      STEP(I)=ZERO
      HS(I)=ZERO
      G(I) = GQ(I)
      D(I)=-G(I)
      gbeg(i) = G(i)
   30 DD=DD+D(I)**2
      CRVMIN=ZERO
      IF (DD .EQ. ZERO) GOTO 160
      DS=ZERO
      SS=ZERO
      GG=DD
      GGBEG=GG
      if (debug) print *, " GGBEG=", GGBEG
C
C     Calculate the step to the trust region boundary and the product HD.
C
   40 ITERC=ITERC+1
      TEMP=DELSQ-SS
      BSTEP=TEMP/(DS+SQRT(DS*DS+DD*TEMP))
c      BSTEP=(-DS+SQRT(DS*DS+DD*TEMP))/DD
      if (debug) print *, " BSTEP=", BSTEP
      GOTO 170
   50 DHD=ZERO
      DO 60 J=1,N
   60 DHD=DHD+D(J)*HD(J)
C
C     Update CRVMIN and set the step-length ALPHA.
C
      ALPHA=BSTEP
      if (debug) then
         print *, " ITERC=",ITERC
         print *, " DHD/DD=", DHD/DD
      endif
      IF (DHD .GT. ZERO) THEN
          TEMP=DHD/DD
          IF (ITERC .EQ. 1) CRVMIN=TEMP
          CRVMIN=MIN(CRVMIN,TEMP)
          ALPHA=MIN(ALPHA,GG/DHD)
      END IF
      QADD=ALPHA*(GG-HALF*ALPHA*DHD)
      QRED=QRED+QADD
C
C     Update STEP and HS.
C
      GGSAV=GG
      GG=ZERO
      DO 70 I=1,N
      STEP(I)=STEP(I)+ALPHA*D(I)
      HS(I)=HS(I)+ALPHA*HD(I)
   70 GG=GG+(G(I)+HS(I))**2
      if (debug) print *, " GG=",GG

      IF (GG .LE. MIN(1.0e-4_PREC*GGBEG,1.0e-16_PREC)) GOTO 160
      IF (GG .LE. 1.0e-14_PREC*gnorm2) GOTO 160

      IF (ITERC .EQ. ITERMAX) GOTO 160
C
C     Begin another conjugate direction iteration if required.
C
      IF (ALPHA .LT. BSTEP) THEN
          IF (QADD .LE. 1.0e-6_PREC*QRED) GOTO 160
          TEMP=GG/GGSAV
          DD=ZERO
          DS=ZERO
          SS=ZERO
          DO 80 I=1,N
          D(I)=TEMP*D(I)-G(I)-HS(I)
          DD=DD+D(I)**2
          DS=DS+D(I)*STEP(I)
   80     SS=SS+STEP(I)**2
          IF (SS .LT. DELSQ) GOTO 40
      END IF
      CRVMIN=ZERO
      ITERSW=ITERC
C
C     Test whether an alternative iteration is required.
C
   90 IF (GG .LE. 1.0e-4_PREC*GGBEG) GOTO 160
      if (debug) print *, "curve search performed"
      SG=ZERO
      SHS=ZERO
      DO 100 I=1,N
      SG=SG+STEP(I)*G(I)
  100 SHS=SHS+STEP(I)*HS(I)
      SGK=SG+SHS
      ANGTEST=SGK/SQRT(GG*DELSQ)
      IF (ANGTEST .LE. -0.99_PREC) GOTO 160
C
C     Begin the alternative iteration by calculating D and HD and some
C     scalar products.
C
      ITERC=ITERC+1
      TEMP=SQRT(DELSQ*GG-SGK*SGK)
      TEMPA=DELSQ/TEMP
      TEMPB=SGK/TEMP
      DO 110 I=1,N
  110 D(I)=TEMPA*(G(I)+HS(I))-TEMPB*STEP(I)
      GOTO 170
  120 DG=ZERO
      DHD=ZERO
      DHS=ZERO
      DO 130 I=1,N
      DG=DG+D(I)*G(I)
      DHD=DHD+HD(I)*D(I)
  130 DHS=DHS+HD(I)*STEP(I)
C
C     Seek the value of the angle that minimizes Q.
C
      CF=HALF*(SHS-DHD)
      QBEG=SG+CF
      QSAV=QBEG
      QMIN=QBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/REAL(IU+1, PREC)
      DO 140 I=1,IU
      ANGLE=REAL(I, PREC)*TEMP
      CTH=COS(ANGLE)
      STH=SIN(ANGLE)
      QNEW=(SG+CF*CTH)*CTH+(DG+DHS*CTH)*STH
      IF (QNEW .LT. QMIN) THEN
          QMIN=QNEW
          ISAVE=I
          TEMPA=QSAV
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=QNEW
      END IF
  140 QSAV=QNEW
      IF (ISAVE .EQ. ZERO) TEMPA=QNEW
      IF (ISAVE .EQ. IU) TEMPB=QBEG
      ANGLE=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-QMIN
          TEMPB=TEMPB-QMIN
          ANGLE=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(REAL(ISAVE, PREC)+ANGLE)
C
C     Calculate the new STEP and HS. Then test for convergence.
C
      CTH=COS(ANGLE)
      STH=SIN(ANGLE)
      REDUC=QBEG-(SG+CF*CTH)*CTH-(DG+DHS*CTH)*STH
      GG=ZERO
      DO 150 I=1,N
      STEP(I)=CTH*STEP(I)+STH*D(I)
      HS(I)=CTH*HS(I)+STH*HD(I)
  150 GG=GG+(G(I)+HS(I))**2
      QRED=QRED+REDUC
      RATIO=REDUC/QRED
      IF (ITERC .LT. ITERMAX .AND. RATIO .GT. 0.01_PREC) GOTO 90
  160 continue
      do 161 i=1,n
         HD(i) = zero
  161 continue
      IH=0
      DO 162 J=1,N
      DO 162 I=1,J
      IH=IH+1
      IF (I .LT. J) HD(J)=HD(J)+HQ(IH)*step(I)
  162 HD(I)=HD(I)+HQ(IH)*step(J)
c      vquad = zero
c      do 163 i=1,n
c  163 vquad = vquad + step(i)*(GQ(i)+HALF*HD(i))
c     &        + XOPT(i)*HD(i)
      vquad= zero
      do 163 i=1,n
  163 vquad = vquad + step(i)*(gbeg(i)+HALF*HD(i))
      INFO=NEWUOA_CONVERGED
      IF (VQUAD .GT. ZERO) THEN
          T=SUM(STEP(1:N)**2)
          IF (SQRT(T) .GE. HALF*DELTA) INFO=NEWUOA_NUMERICAL_ERROR
      END IF
      RETURN
C
C     The following instructions act as a subroutine for setting the vector
C     HD to the vector D multiplied by the second derivative matrix of Q.
C     They are called from three different places, which are distinguished
C     by the value of ITERC.
C
  170 continue

      DO 315 I=1,N
  315 HD(I) = ZERO
      IH=0
      DO 320 J=1,N
      DO 320 I=1,J
      IH=IH+1
      IF (I .LT. J) HD(J)=HD(J)+HQ(IH)*D(I)
  320 HD(I)=HD(I)+HQ(IH)*D(J) 

      IF (ITERC .EQ. 0) GOTO 20
      IF (ITERC .LE. ITERSW) GOTO 50
      GOTO 120
      END

      subroutine f_grad (mv,v_base,v_gtemp)
      integer mv
      REAL (PREC), INTENT(IN), DIMENSION(:)  :: v_base
      REAL (PREC), INTENT(OUT), DIMENSION(:)  :: v_gtemp
      integer m1

      do 10 m1=1,mv
         v_gtemp(m1)=2.0_PREC*v_base(m1)
   10 continue  

      return
      end

