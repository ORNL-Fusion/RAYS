module zfunctions_m

!   calculates plasma function defined by Eq.(8-82), Stix.
! N.B. This assumes that kz is real and non-zero

    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind), parameter :: zero = 0.0_rkind, one = 1.0_rkind, two = 2.0_rkind
    real(KIND=rkind), parameter :: pi = atan(zero, -one), sqrt_pi = sqrt(pi)

    interface zfun0
        module procedure zfun0_S, zfun0_D
    end interface

    interface zfun
        module procedure zfun_S, zfun_D
    end interface

    interface zfun_real_arg
        module procedure zfun_real_arg_spline_D
    end interface

    interface zfun0_real_arg
        module procedure zfun0_real_arg_D
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
    1 zzr=-sqrt_pi*wzi1
      zzi=sqrt_pi*wzr1
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
      x=real(z,kind=rkind)
      y=aimag(z)
      call zzdisp_D(x,y,re,im)
      zfun_D=cmplx(re,im, rkind)
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
      disp1_D=-two-two*zz*dp
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
      a=2._rkind*x1*y1
      b=-(x1*x1-y1*y1)
      abr=two*exp(b)*cos(a)
      abi=-two*exp(b)*sin(a)
      wzr1=abr-wzr1
      wzi1=abi-wzi1
      if(x.le.0..and.y.le.0.) go to 1
    2 wzi1=-wzi1
    1 zzr=-sqrt_pi*wzi1
      zzi=sqrt_pi*wzr1
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
      real(kind = rkind) h, h2, lambda
      real(kind = rkind) ::r1, s, sr, si, tr, ti, c, cc, rr, ri
      integer capn, nu, i, n, nup, np1
      logical b
      epsh=1.e-12_rkind
      if(y.lt.4.29 .and. x.lt.5.33) go to 10
!
!  (x,y) belongs to q1-r
!
      h=zero
      capn=0
      nu=8
      go to 20
!
!  (x,y) belongs to r
!
   10 s=(1.-y/4.29)*sqrt(1.-x*x/28.41)
      h=1.6_rkind*s
      h2=two*h
      capn=6._rkind+23._rkind*s+.5_rkind
      nu=9._rkind+21._rkind*s+.5_rkind
      lambda=h2**capn
   20 b=(h.eq.0. .or. lambda.lt.epsl)
!
!  statement (lambda.lt.epsl) covers the underflow case
!  when h(.gt.0) is very small.
!
      rr=zero
      ri=zero
      sr=zero
      si=zero
      nup=nu+1
      do 100 i=1,nup
      n=nup-i
      np1=n+1
      tr=y+h+np1*rr
      ti=x-np1*ri
      c=.5_rkind/(tr*tr+ti*ti)
      rr=c*tr
      ri=c*ti
      if(.not.(h.gt.0..and.n.le.capn)) go to 100
      tr=lambda+sr
      sr=rr*tr-ri*si
      si=ri*tr+rr*si
      lambda=lambda/h2
  100 continue
      cc=1.12837916709551_rkind
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

 complex(kind = rkind) function zfun0_real_arg_D(z, kz)
! Calculates plasma function defined by Eq.(8-82), Stix.
! N.B. This assumes that z and kz are real and kz is non-zero

    implicit none

    real(kind = rkind) :: z
    real(kind = rkind) ::kz

!    complex, external :: zfun_real_arg_D

!   See Eq.(8-82).
    if ( kz > 0. ) then
          zfun0_real_arg_D = zfun_real_arg_spline_D(z)
    else if ( kz < 0. ) then
          zfun0_real_arg_D = -zfun_real_arg_spline_D(-z)
    else
       write(0, *) 'ZFUN0_real_arg_D: Error, kz must be real and non-zero, kz= ', kz
       stop 1
    end if
    return
 end function zfun0_real_arg_D

 !***********************************************************************

 complex(kind = rkind) function zfun_real_arg_spline_D(z)

! N.B. !!  This thing is not thread safe.  Maybe work on it more later.

! Calculates plasma dispersion function for real argument.
! For |x| <= spline_range, the real part of Z(x) is calculated by spline interpolation
! of Zfun(x) on grid -> x_grid.
! For |x| > spline_range, the real part of Z(x) is calculated from the asymptotic expansion.
! The imaginary part of Z(x) is analytical = sqrt(pi)*i*exp(-x**2)
! Grid parameters are set below: nx, x_grid_min, x_grid_max.
! N.B. Calling with x outside [x_grid_min, x_grid_max] will produce an error.
!
! On the first call an initialization is performed that calculates the spline coefficients.
! Spline coefficients for both real and imaginary parts are calculated, but only the real
! part is presently splined.  The coding for evaluating the coefficients using splines
! for the imaginary part is is there if someone wants to use it but commented out.
!
! The spline routines can also return the 1st and 2nd derivatives (controlled by iselect
! below).  Right now I don't use these, but maybe later.


	implicit none
	real(kind = rkind) :: z

	logical, save :: initialized = .false.
    integer :: i

! Stuff for Cubic splines
! Grid definition
    integer, parameter :: nx = 2001
    real(KIND=rkind), parameter ::  spline_range = 10.0_rkind
    real(KIND=rkind), parameter ::  x_grid_min = -spline_range, x_grid_max = spline_range
    real(KIND=rkind), save ::  x_grid(nx)

! Stuff for asymptotic expansion
! Order of the expansion is specified by N_asymp = number of terms in expansion.
! N_asymp = 6 (i.e. 1/x**11) gives accuracy ~10**-11. Coefficients are provided up to
! N_asymp = 8.
    integer, parameter :: N_asymp = 6
    real(KIND=rkind), parameter ::  A_coef(8)= (/1._rkind, 1._rkind/2._rkind, &
     & 3._rkind/4._rkind, 15._rkind/8._rkind, 105._rkind/16._rkind, 945._rkind/32._rkind, &
     & 10395._rkind/64._rkind, 135135._rkind/128._rkind/)
    real(KIND=rkind) :: z_inv


! bcspeval arguments
	integer ibcxmin, ibcxmax, ibcthmin, ibcthmax, ilinx, ilnx,ier
    real(KIND=rkind) :: bcxmin(nx),bcxmax(nx)      ! (nx) if used
    real(KIND=rkind) :: bcthmin(nx),bcthmax(nx)  ! (inx) if used
    real(KIND=rkind) ::  fsplRe(4, nx)=0., fsplIm(4, nx)=0.
    real(KIND=rkind) ::  fvalRe(3)=0., fvalIm(3)=0.
    integer :: iselect(3) = (/ 1, 0, 0 /) ! Can change to /1,1,1/ to get derivatives
    real(KIND=rkind) ::  wk(nx)

    complex(KIND=rkind)  :: fC

! If this is first call initialize splines
	init: if (.not. initialized) then

! Generate spline grid and evaluate function
		do i = 1, nx
			x_grid(i) = x_grid_min + (i-1)*(x_grid_max - x_grid_min)/(nx-1)
			fC = zfun(cmplx(x_grid(i), 0.0_rkind, kind=rkind))
			fsplRe(1,i) = real(fC, rkind)
			fsplIm(1,i) = aimag(fC)

		end do

! Get spline coefficients
		ibcxmin = 0
		bcxmin = zero
		ibcxmax = 0
		bcxmax = zero
		ibcthmin = 0
		bcthmin = zero
		ibcthmax = 0
		bcthmax =zero
		ilinx = 0

		call cspline(x_grid,nx,fsplRe,ibcxmin,bcxmin,ibcxmax,bcxmax,wk,nx,ilinx,ier)
		call cspline(x_grid,nx,fsplIm,ibcxmin,bcxmin,ibcxmax,bcxmax,wk,nx,ilinx,ier)

		initialized = .true.
! 		write(*,*) 'initialized = ',initialized
! 		stop
	end if init

	if (abs(z) <= spline_range) then ! spline real part
		call cspeval(z,iselect,fvalRe,x_grid,nx,ilinx,fsplRe,ier)
	else ! asymptotic expansion of real part
		z_inv = one/z
		fvalRe(1) = zero
		do i = 1, N_asymp
			fvalRe(1) = fvalRe(1) - z_inv**(2*i - 1)* A_coef(i)
		end do
	end if

	fvalIm(1) = sqrt_pi*exp(-z**2)
	if (ier .ne. 0) write (*,*) 'cspeval: ier = ', ier

	zfun_real_arg_spline_D = cmplx(fvalRe(1), fvalIm(1), rkind)
	return

 end function zfun_real_arg_spline_D

 end module zfunctions_m

