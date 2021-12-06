c
c THIS CODE IS SUPERIOR TO PDF AND PDFLL.  NO ERROR AT ZETA=(4.,0.1)
c  disp(z) is the plasma z function
c
      complex function disp(z)
      complex z
      real im
      x=real(z)
      y=aimag(z)
      call zzdisp(x,y,re,im)
      DISP=cmplx(re,im)
      return
      end function disp
c
c***************************************************************************
c
c
c  disp1(z) calculates first derivative of disp(z)
c
      complex function disp1(z)
      complex z,dp,disp
      complex zz
      zz=z
      dp=disp(z)
      disp1=-2.0-2.0*zz*dp
      return
      end function
c
c***************************************************************************
c
c
c  calculates plasma z function
c
      subroutine zzdisp(x,y,zzr,zzi)
      x1=abs(x)
      y1=abs(y)
      call wzdisp(x1,y1,wzr1,wzi1)
      if(x.ge.0..and.y.ge.0.) go to 1
      if(x.le.0..and.y.ge.0.) go to 2
      a=2*x1*y1
      b=-(x1*x1-y1*y1)
      abr=2*exp(b)*cos(a)
      abi=-2*exp(b)*sin(a)
      wzr1=abr-wzr1
      wzi1=abi-wzi1
      if(x.le.0..and.y.le.0.) go to 1
    2 wzi1=-wzi1
    1 zzr=-1.7724538509055*wzi1
      zzi=1.7724538509055*wzr1
      return
      end
c
c***************************************************************************
c
cc   needed for z function
cc
      subroutine wzdisp(x,y,re,im)
      equivalence(epsh,epsl,epsy)
      real im,lambda
      integer capn
      logical b
      epsh=1.e-12
      if(y.lt.4.29 .and. x.lt.5.33) go to 10
c
c  (x,y) belongs to q1-r
c
      h=0.
      capn=0
      nu=8
      go to 20
c
c  (x,y) belongs to r
c
   10 s=(1.-y/4.29)*sqrt(1.-x*x/28.41)
      h=1.6*s
      h2=2.*h
      capn=6.+23*s+.5
      nu=9.+21*s+.5
      lambda=h2**capn
   20 b=(h.eq.0. .or. lambda.lt.epsl)
c
c  statement (lambda.lt.epsl) covers the underflow case
c  when h(.gt.0) is very small.
c
      rr=0.
      ri=0.
      sr=0.
      si=0.
      nup=nu+1
      do 100 i=1,nup
      n=nup-i
      np1=n+1
      tr=y+h+np1*rr
      ti=x-np1*ri
      c=.5/(tr*tr+ti*ti)
      rr=c*tr
      ri=c*ti
      if(.not.(h.gt.0..and.n.le.capn)) go to 100
      tr=lambda+sr
      sr=rr*tr-ri*si
      si=ri*tr+rr*si
      lambda=lambda/h2
  100 continue
      cc=1.12837916709551
      if(y.lt.epsy) go to 120
      if(b) go to 110
      re=sr*cc
      go to 130
  110 re=rr*cc
      go to 130
  120 re=exp(-x*x)
  130 if(b) go to 140
      im=si*cc
      go to 150
  140 im=ri*cc
  150 return
      end
c
c***************************************************************************
c
      SUBROUTINE ZFUN(Z,FU)
C
C     ROUTINE WHICH EVALUATES THE PLASMA DISPERSION FUNCTION.  USES
C     NUMERICAL INTEGRATION (ABSOLUTE VALUE OF THE COMPLEX ARGUMENT, Z,
C     LESS THAN 5) OR ASYMPTOTIC EXPANSION.
C
      COMPLEX Z,FU,TEMP1,TEMP2,Z2,TPIIOD
      DIMENSION C(21),D(21),E(21),F(21),W(21)
      DATA DELTA/5.0E-1/, YI/-1.0E+0/, N/21/, ITEST/0/
      DATA ZERO/0.0E+0/, OSQPI/5.6418958355E-1/, TPI/6.28318530718E+0/
C
      IF(ITEST .EQ. 1) GO TO 1
      ITEST=1
C***     DEFINE WEIGHTS AND STORE CONSTANTS USED FOR INTEGRATION.
C***     WEIGHTS ARE DERIVED FROM A 3 POINT INTEGRATION SCHEME.
      NM3=N-3
      CONOI=DELTA*OSQPI
      W(1)=CONOI*3.75E-1
      W(2)=CONOI*1.1666666667E+0
      W(3)=CONOI*9.5833333333E-1
      DO 20 I=1,3
        NMIP1=N-I+1
   20 W(NMIP1)=W(I)
      DO 30 I=4,NM3
   30 W(I)=CONOI
      NO2=N/2
      NO2P1=NO2+1
      X=ZERO
      Y=YI
      E(NO2P1)=X
      F(NO2P1)=Y
      TEMP1=CMPLX(X,Y)
      TEMP2=CEXP(-TEMP1*TEMP1)
      C(NO2P1)=REAL(TEMP2)*W(NO2P1)
      D(NO2P1)=AIMAG(TEMP2)*W(NO2P1)
      DO 200 I=1,NO2
      X=DELTA*FLOAT(I)
      NPI=NO2P1+I
      NMI=NO2P1-I
      TEMP1=CMPLX(X,Y)
      TEMP2=CEXP(-TEMP1*TEMP1)
      C(NPI)=REAL(TEMP2)*W(NPI)
      C(NMI)=REAL(TEMP2)*W(NMI)
      D(NPI)=AIMAG(TEMP2)*W(NPI)
      D(NMI)=AIMAG(TEMP2)*W(NMI)*(-1.0E+0)
      E(NMI)=-X
      E(NPI)=X
      F(NPI)=Y
      F(NMI)=Y
  200 CONTINUE
      TPIOD=TPI/DELTA
      TPIIOD=CMPLX(ZERO,TPIOD)
C***     BEGIN CALCULATIONS.
    1 G=REAL(Z)
      YY=AIMAG(Z)
      H=ABS(YY)
      ZA=Z*CONJG(Z)
      IF(ZA .GE. 2.5E+1) GO TO 5
      Z2=Z*Z
C***     NUMERICAL INTEGRATION.
C***     F=1/SQRT(PI)*SUM OF...W(I)*EXP(-X(I)**2)/(X(I)-Z)...I=1,N.
C***     INTEGRATION IS ALONG A LINE X(I) IN THE COMPLEX PLANE, WHERE
C***     THE IMAGINARY PART OF X(I)=YI AND THE DIFFERENCE BETWEEN
C***     SUCCESSIVE REAL PARTS OF X(I)=DELTA.  LIMITS OF INTEGRATION
C***     ARE FROM -DELTA*N/2 TO DELTA*N/2.
C***     COMPUTE THE INTEGRAL BY TAKING THE SUM FROM 1 TO N OF THE
C***     CONSTANTS DIVIDED BY X(I)-Z.  USES REAL ARITHMETIC.
      ZR=0.0E+0
      ZI=0.0E+0
      DO 7 I=1,N
      A=E(I)-G
      B=F(I)-H
      DEN=A*A+B*B
      ODEN=1.0E+0/DEN
      ZR=ZR+(A*C(I)+B*D(I))*ODEN
      ZI=ZI+(A*D(I)-B*C(I))*ODEN
    7 CONTINUE
C***     ADD THE CORRECTION TERM.
      FU=CMPLX(ZR,ZI)+(0.0E+0,-3.5449077018E+0)*CEXP(-Z2-TPIOD*
     *   (H-YI)+TPIIOD*G)
      IF(YY .GE. ZERO) GO TO 6
C***     IMAGINARY PART OF ARGUMENT IS NEGATIVE.
      FU=CONJG(FU)+(0.0E+0,3.5449077018E+0)*CEXP(-Z2)
      GO TO 6
C***     MAGNITUDE OF ARGUMENT IS GREATER THAN 5, USE
C***     ASYMPTOTIC EXPANSION.
    5 CALL AEXPAN(Z,G,H,YY,FU)
    6 RETURN
      END
      SUBROUTINE AEXPAN(Z,G,H,YY,FU)
C
C     ROUTINE WHICH COMPUTES THE PLASMA DISPERSION FUNCTION USING
C     ASYMPTOTIC EXPANSION.  IF THE IMAGINARY PART OF THE ARGUMENT, YY,
C     IS EQUAL TO ZERO REAL ARITHMETIC IS USED.
C
      COMPLEX Z,FU,A,Z2,OZ2
      DATA ZERO/0.0E+0/, N/8/
C
      IF(YY .EQ. ZERO) GO TO 1
C***     COMPLEX ARITHMETIC.
      Z2=Z*Z
      OZ2=1.0E+0/Z2
      FU=-1.0E+0/Z
      A=FU
      EN=5.0E-1
      DO 10 I=1,N
      A=EN*A*OZ2
      FU=FU+A
   10 EN=EN+1.0E+0
      IF(YY .GT. ZERO) GO TO 30
      IF(H .GT. SQRT(G*G+1.72E+2)) GO TO 20
      FU=FU+(0.0E+0,3.5449077018E+0)*CEXP(-Z2)
      GO TO 30
C***     ERROR STOP TO AVOID OVERFLOW.
   20 WRITE (51,100) Z
  100 FORMAT(//1X,'*** ERROR STOP IN FRDCNT ROUTINE, ARGUMENT IS',
     * ' TOO SMALL,  ARG =',1P2E14.6)
      STOP
C***     REAL ARITHMETIC.
    1 X2=G*G
      OX2=1.0E+0/X2
      F=-1.0E+0/G
      B=F
      EN=5.0E-1
      DO 11 I=1,N
      B=EN*B*OX2
      F=F+B
   11 EN=EN+1.0E+0
      C=1.7724538509E+0*EXP(-X2)
      FU=CMPLX(F,C)
   30 RETURN
      END
