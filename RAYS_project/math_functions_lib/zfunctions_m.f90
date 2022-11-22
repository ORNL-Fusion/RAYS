module zfunctions_m

!   calculates plasma function defined by Eq.(8-82), Stix.
! N.B. This assumes that kz is real and non-zero

    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals

    interface zfun0
        module procedure zfun0_S, zfun0_D
    end interface

    interface zfun
        module procedure zfun_S, zfun_D
    end interface

!    private
!    public zfun, zfun0
     contains

!***********************************************************************
    
 complex function zfun0_S(z, kz)
! Calculates plasma function defined by Eq.(8-82), Stix.
! N.B. This assumes that kz is real and non-zero

    implicit none

    complex :: z
    real :: kz

!   See Eq.(8-82).
    if ( kz > 0. ) then
          zfun0_S = zfun_S(z)
    else if ( kz < 0. ) then
          zfun0_S = -zfun_S(-z)
    else
       write(0, *) 'ZFUN0: Error, kz must be real and non-zero, kz= ', kz
       stop 1
    end if
    return
 end function zfun0_S
!***********************************************************************
                                                                        
!  zfun(z) is the plasma dispersion z function                          
!                                                                       
      complex function zfun_S(z) 
      complex, intent(in) :: z 
      real :: x, y, re, im 
      x=real(z) 
      y=aimag(z) 
      call zzdisp(x,y,re,im) 
      zfun_S=cmplx(re,im) 
      return 
      END function zfun_S                                     
!***********************************************************************
!                                                                       
!                                                                       
!  disp1(z) calculates first derivative of disp(z)                      
!                                                                       
      complex function disp1(z) 
      complex, intent(in) :: z 
!      complex :: zz,dp,zfun 
      complex :: zz,dp 
      zz=z
      dp=zfun_S(z) 
      disp1=-2.0-2.0*zz*dp 
      return 
      END                                           
!                                                                       
!***********************************************************************
                                                                     
!                                                                       
!  calculates plasma z function                                         
!                                                                       
      subroutine zzdisp(x,y,zzr,zzi)
      real, intent(in) :: x, y
      real, intent(out) :: zzr,zzi
      real :: x1, y1, a, b, abr, abi, wzr1, wzi1
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
      END                                           
!                                                                       
!***********************************************************************
!                                                                       
!c   needed for z function                                              
!c                                                                      
      subroutine wzdisp(x,y,re,im) 
      real, intent(in) :: x, y
      real, intent(out) :: re, im
      real :: epsh,epsl,epsy
      equivalence(epsh,epsl,epsy) 
      real h, h2, lambda
      real :: r1, s, sr, si, tr, ti, c, cc, rr, ri
      integer capn, nu, i, n, nup, np1
      logical b 
      epsh=1.e-12 
      if(y.lt.4.29 .and. x.lt.5.33) go to 10 
!                                                                       
!  (x,y) belongs to q1-r                                                
!                                                                       
      h=0. 
      capn=0 
      nu=8 
      go to 20 
!                                                                       
!  (x,y) belongs to r                                                   
!                                                                       
   10 s=(1.-y/4.29)*sqrt(1.-x*x/28.41) 
      h=1.6*s 
      h2=2.*h 
      capn=6.+23*s+.5 
      nu=9.+21*s+.5 
      lambda=h2**capn 
   20 b=(h.eq.0. .or. lambda.lt.epsl) 
!                                                                       
!  statement (lambda.lt.epsl) covers the underflow case                 
!  when h(.gt.0) is very small.                                         
!                                                                       
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
      END                                           
!***********************************************************************
     
 complex(kind = rkind) function zfun0_D(z, kz)
! Calculates plasma function defined by Eq.(8-82), Stix.
! N.B. This assumes that kz is real and non-zero

    implicit none

    complex(kind = rkind) :: z
    real(kind = rkind) ::kz

!    complex, external :: zfun

!   See Eq.(8-82).
    if ( kz > 0. ) then
          zfun0_D = zfun_D(z)
    else if ( kz < 0. ) then
          zfun0_D = -zfun_D(-z)
    else
       write(0, *) 'ZFUN0: Error, kz must be real and non-zero, kz= ', kz
       stop 1
    end if
    return
 end function zfun0_D
!***********************************************************************

!  zfun(z) is the plasma dispersion z function                          
!                                                                       
      complex(kind = rkind) function zfun_D(z) 
      complex(kind = rkind), intent(in) :: z 
      real(kind = rkind) ::x, y, re, im 
      x=real(z) 
      y=aimag(z) 
      call zzdisp_D(x,y,re,im) 
      zfun_D=cmplx(re,im) 
      return 
      END function zfun_D                                  
!                                                                       
!***********************************************************************
!                                                                       
!                                                                       
!  disp1(z) calculates first derivative of disp(z)                      
!                                                                       
      complex(kind = rkind) function disp1_D(z) 
      complex(kind = rkind), intent(in) :: z 
      complex(kind = rkind) :: zz,dp 
      zz=z
      dp=zfun_D(z) 
      disp1_D=-2.0-2.0*zz*dp 
      return 
      END                                           
!                                                                       
!***********************************************************************
!                                                                       
!                                                                       
!  calculates plasma z function                                         
!                                                                       
      subroutine zzdisp_D(x,y,zzr,zzi)
      real(kind = rkind), intent(in) :: x, y
      real(kind = rkind), intent(out) :: zzr,zzi
      real(kind = rkind) ::x1, y1, a, b, abr, abi, wzr1, wzi1
      x1=abs(x) 
      y1=abs(y) 
      call wzdisp_D(x1,y1,wzr1,wzi1) 
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
      END                                           
!                                                                       
!***********************************************************************
!                                                                       
!c   needed for z function                                              
!c                                                                      
      subroutine wzdisp_D(x,y,re,im) 
      real(kind = rkind), intent(in) :: x, y
      real(kind = rkind), intent(out) :: re, im
      real(kind = rkind) ::epsh,epsl,epsy
      equivalence(epsh,epsl,epsy) 
      real h, h2, lambda
      real(kind = rkind) ::r1, s, sr, si, tr, ti, c, cc, rr, ri
      integer capn, nu, i, n, nup, np1
      logical b 
      epsh=1.e-12 
      if(y.lt.4.29 .and. x.lt.5.33) go to 10 
!                                                                       
!  (x,y) belongs to q1-r                                                
!                                                                       
      h=0. 
      capn=0 
      nu=8 
      go to 20 
!                                                                       
!  (x,y) belongs to r                                                   
!                                                                       
   10 s=(1.-y/4.29)*sqrt(1.-x*x/28.41) 
      h=1.6*s 
      h2=2.*h 
      capn=6.+23*s+.5 
      nu=9.+21*s+.5 
      lambda=h2**capn 
   20 b=(h.eq.0. .or. lambda.lt.epsl) 
!                                                                       
!  statement (lambda.lt.epsl) covers the underflow case                 
!  when h(.gt.0) is very small.                                         
!                                                                       
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
      END                                           
    
 end module zfunctions_m

