PROGRAM test_bisect
   USE bisect_m, only : solve_bisection

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
!    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

    INTEGER, PARAMETER :: npoints = 101
    REAL :: x, x_low, x_high, y, eps, delta
    REAL :: a, b, data(2)
    REAL, EXTERNAL :: fun, funX
!    REAL :: fun
    REAL(kind = rkind) :: xD, xD_low, xD_high, yD, epsD, deltaD
    REAL(kind = rkind) :: aD, bD, dataD(2)
    REAL(kind = rkind), EXTERNAL :: funD, funDX
    integer :: ierr


! test single
	x_low = 0.
	x_high = 6.
	eps = 1.e-6
	y = 15.

    write(*,*) ' '
    write(*,*) 'single precision: y = ', y
    CALL solve_bisection(fun, x,  x_low, x_high, y, eps, ierr)
    delta = x - sqrt(y)
    WRITE (*, *) ' x  = ', x, '  delta = ', delta, '  ierr = ', ierr

! test real (i.e. double)
	xD_low = 0.0_rkind
	xD_high = 6.0_rkind
	epsD = 1.e-10_rkind

	yD = 15.0_rkind

    write(*,*) ' '
    write(*,*) 'real precision: y = ', yD
    CALL solve_bisection(funD, xD,  xD_low, xD_high, yD, epsD, ierr)
    deltaD = xD - sqrt(yD)
    WRITE (*, *) ' x  = ', xD, '  delta = ', deltaD, '  ierr = ', ierr

	y = 15.9

    write(*,*) ' '
    write(*,*) 'single precision: y = ', y
    CALL solve_bisection(fun, x,  x_low, x_high, y, eps, ierr)
    delta = x - sqrt(y)
    WRITE (*, *) ' x  = ', x, '  delta = ', delta, '  ierr = ', ierr

	yD = 15.9_rkind

    write(*,*) ' '
    write(*,*) 'real precision: y = ', yD
    CALL solve_bisection(funD, xD,  xD_low, xD_high, yD, epsD, ierr)
    deltaD = xD - sqrt(yD)
    WRITE (*, *) ' x  = ', xD, '  delta = ', deltaD, '  ierr = ', ierr

	y = 0.01

    write(*,*) ' '
    write(*,*) 'single precision: y = ', y
    CALL solve_bisection(fun, x,  x_low, x_high, y, eps, ierr)
!    delta = x - y**2 ! for f = f = x**2
    delta = x - sqrt(y) ! for f = sqrt(x)
    WRITE (*, *) ' x  = ', x, '  delta = ', delta, '  ierr = ', ierr

	yD = 0.01_rkind

    write(*,*) ' '
    write(*,*) 'real precision: y = ', yD
    CALL solve_bisection(funD, xD,  xD_low, xD_high, yD, epsD, ierr)
!     deltaD = xD - yD**2 ! for f = x**2
    deltaD = xD - sqrt(yD)! for f = sqrt(x)
    WRITE (*, *) ' x  = ', xD, '  delta = ', deltaD, '  ierr = ', ierr

!****************************************************************************
    write(*,*) ' '
    write(*,*) 'Test x0 > x1'

	yD = 0.01_rkind

    write(*,*) ' '
    write(*,*) 'real precision: y = ', yD
    CALL solve_bisection(funD, xD, xD_high, xD_low, yD, epsD, ierr)
!     deltaD = xD - yD**2 ! for f = x**2
    deltaD = xD - sqrt(yD)! for f = sqrt(x)
    WRITE (*, *) ' x  = ', xD, '  delta = ', deltaD, '  ierr = ', ierr


!****************************************************************************
! Test extra data
!****************************************************************************

    write(*,*) ' '
    write(*,*) 'Test extra data'

! test single extra data
	x_low = 0.
	x_high = 6.
	eps = 1.e-6
	y = 3.
	a =1.0
	b = 1.0
	data = [a,b]

    write(*,*) ' '
    write(*,*) 'single precision, extra data: y = ', y
    CALL solve_bisection(funX, x,  x_low, x_high, y, eps, ierr, data)
    delta = x - ((y-b)/a)**2
    WRITE (*, *) ' x  = ', x, '  delta = ', delta, '  ierr = ', ierr

!****************************************************************************
! test real extra data
	xD_low = 0._rkind
	xD_high = 6._rkind
	epsD = 1.e-10_rkind
	yD = 3._rkind
	aD =1.0_rkind
	bD = 1.0_rkind
	dataD = [aD,bD]

    write(*,*) ' '
    write(*,*) 'single precision, extra data: y = ', yD
    CALL solve_bisection(funDX, xD,  xD_low, xD_high, yD, epsD, ierr, dataD)
    deltaD = xD - ((yD-bD)/aD)**2
    WRITE (*, *) ' x  = ', xD, '  delta = ', deltaD, '  ierr = ', ierr


!****************************************************************************
! contains
!****************************************************************************

END PROGRAM test_bisect


   function fun(x)
	IMPLICIT NONE
	REAL :: fun, x
	fun = x**2
!	fun = sqrt(x)
	return
  end function fun

  function funD(x)
	IMPLICIT NONE
    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals

    real(KIND=rkind) funD, x
	funD = x**2
!	funD = sqrt(x)
	return
  end function funD

   function funX(x, data)
	IMPLICIT NONE
	REAL :: funX, x, data(2), a, b

	a = data(1)
	b = data(2)
! 	fun = x**2
	funX = a*sqrt(x) + b
	return
  end function funX

  function funDX(x, data)
	IMPLICIT NONE
    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals

    real(KIND=rkind) funDX, x, data(2), a, b

	a = data(1)
	b = data(2)

! 	funD = x**2
	funDX = a*sqrt(x) + b
	return
  end function funDX
