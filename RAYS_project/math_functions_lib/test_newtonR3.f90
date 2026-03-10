PROGRAM test_newtonR3
!   USE newtonR3_m, only : newtonR3_single, newtonR3_real
   USE newtonR3_m, only : solve_newtonR3

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
!    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

    REAL :: x(3), x_start(3)
    REAL :: step_max
    REAL :: y
    REAL :: eps, delta
    REAL :: a, b, Xdata(2)

    REAL(kind = rkind) :: xD(3), xD_start(3)
    REAL(kind = rkind) :: step_maxD
    REAL(kind = rkind) :: yD
    REAL(kind = rkind) :: epsD, deltaD
    REAL(kind = rkind) :: aD, bB, XdataD(2)

    integer :: i_iter, nsig, max_iter = 100

! test single
	x_start = (/0.1, 0.2, 0.3/)
	step_max = 0.1
	eps = 1.e-6
	nsig = 8
	y = 1.

    write(*,*) ' '
    write(*,*) 'single precision: y = ', y
    CALL solve_newtonR3(fun, gradfun, x,  x_start, step_max, y, eps, nsig, i_iter,&
            & max_iter_in=100)
    WRITE (*, *) ' x  = ', x, '  fun(x) = ',fun(x), '  i_iter = ', i_iter

! test real (i.e. double)
	xD_start = (/0.1_rkind, 0.2_rkind, 0.3_rkind/)
	step_maxD = 0.1_rkind
	epsD = 1.0e-6_rkind
	nsig = 8
	yD = 1.0_rkind

    write(*,*) ' '
    write(*,*) 'real precision: y = ', yD
    CALL  solve_newtonR3(funD, gradfunD, xD,  xD_start, step_maxD, yD, epsD, nsig, i_iter,&
            & max_iter_in=100)
    WRITE (*, *) ' xD  = ', xD, '  funD(xD) = ',funD(xD), '  i_iter = ', i_iter



!****************************************************************************
! Test extra data
!****************************************************************************

!     write(*,*) ' '
!     write(*,*) 'Test extra data'
!
! ! test single extra data
! 	x_start = (/0.1, 0.2, 0.3/)
! 	step_max = 0.1
! 	y = 1.
! 	eps = 1.e-6
! 	nsig = 8
! 	Xdata = [1.0, 0.]
!
!     write(*,*) ' '
!     write(*,*) 'single precision: y = ', y
!
!     CALL solve_newtonR3(funX, gradfunX, x,  x_start, step_max, y, eps, nsig, i_iter,&
!             & Xdata)
! !            & Xdata,  max_iter_in=100)
!     WRITE (*, *) ' x  = ', x, '  fun(x) = ',fun(x), '  i_iter = ', i_iter
!
! ! test real (i.e. double) extra data
! 	xD_start = (/0.1_rkind, 0.2_rkind, 0.3_rkind/)
! 	step_maxD = 0.1_rkind
! 	yD = 1.0_rkind
! 	epsD = 1.0e-6_rkind
! 	nsig = 8
! 	XdataD = [1.0_rkind, 0.0_rkind]
!
!     write(*,*) ' '
!     write(*,*) 'real precision: y = ', yD
!     CALL  solve_newtonR3(funDX, gradfunDX, xD,  xD_start, step_maxD, yD, epsD, nsig, i_iter,&
!             & XdataD)
! !            & XdataD, max_iter_in=100)
!     WRITE (*, *) ' xD  = ', xD, '  funD(xD) = ',funD(xD), '  i_iter = ', i_iter
!

!****************************************************************************
contains
!****************************************************************************

   function fun(x)
	IMPLICIT NONE
	real, intent(in) :: x(3)
	real :: fun
	fun = dot_product(x,x) + x(1)
	return
   end function fun

   function gradfun(x)
	IMPLICIT NONE
	real,intent(in) :: x(3)
	REAL :: gradfun(3)
	gradfun(:) = 2.0*x(:) + (/1.,0.,0./)
	return
   end function gradfun

  function funD(x)
	IMPLICIT NONE
    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind), intent(in) :: x(3)
    real(KIND=rkind) ::	funD
    funD = dot_product(x,x)
	return
   end function funD


   function gradfunD(x)
	IMPLICIT NONE
    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind), intent(in) :: x(3)
    real(KIND=rkind) :: gradfunD(3)
	gradfunD = 2.0_rkind*x
	return
   end function gradfunD


   function funX(x, extra_data)
	IMPLICIT NONE
	real, intent(in) :: x(3)
	real, intent(in) :: extra_data(2)
	real :: funX

	real :: a, b

	a = extra_data(1)
	b = extra_data(2)
	funX = a*dot_product(x,x) + b
	return
  end function funX

  function gradfunX(x, extra_data)
	IMPLICIT NONE
	real, intent(in) :: x(3)
	real, intent(in) :: extra_data(2)
	real :: gradfunX(3)

	real :: a, b

	a = extra_data(1)
	b = extra_data(2)
	gradfunX = 2.0_rkind*a*x
	return
  end function gradfunX


   function funDX(x, extra_data)
	IMPLICIT NONE
    real(KIND=rkind), intent(in) :: x(3)
    real(KIND=rkind), intent(in) :: extra_data(2)
    real(KIND=rkind) :: funDX

    real(KIND=rkind) :: a, b

	a = extra_data(1)
	b = extra_data(2)
	funDX = a*dot_product(x,x) + b
	return
  end function funDX

  function gradfunDX(x, extra_data)
	IMPLICIT NONE
    real(KIND=rkind), intent(in) :: x(3)
    real(KIND=rkind), intent(in) :: extra_data(2)
    real(KIND=rkind) :: gradfunDX(3)

    real(KIND=rkind) :: a, b

	a = extra_data(1)
	b = extra_data(2)
	gradfunDX = 2.0_rkind*a*x
	return
  end function gradfunDX


END PROGRAM test_newtonR3
