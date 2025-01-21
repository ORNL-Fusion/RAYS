PROGRAM fun_invert
   USE mono_funct_inversion_m, only : list_invert

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
!    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

!    TYPE(cube_spline_function_1D) :: sp_fun

    INTEGER, PARAMETER :: npoints = 5, n_test = 4

!     REAL, allocatable :: x(:)
!     REAL :: x_low, x_high, y, eps, delta
!     REAL, EXTERNAL :: fun

    REAL(kind = rkind) :: xD(npoints), yD(npoints)
    REAL(kind = rkind) :: xD_low, xD_high, epsD, deltaD
    REAL(kind = rkind), EXTERNAL :: fun_real

    REAL(kind = rkind) :: yD_low, yD_high
    REAL(kind = rkind) :: x_calc(npoints), y_calc(npoints)
    REAL(kind = rkind) :: x_test(n_test), y_test(n_test)
    integer :: i

 	CHARACTER (len = 80), parameter :: f_name = 'f(x)'

! test single
! 	x_low = 0.
! 	x_high = 6.
! 	eps = 1.e-6
! 	y = 15.
!
!     write(*,*) ' '
!     write(*,*) 'single precision: y = ', y
!     CALL solve_bisection(fun, x,  x_low, x_high, y, eps, ierr)
!     delta = x - sqrt(y)
!     WRITE (*, *) ' x  = ', x, '  delta = ', delta, '  ierr = ', ierr

	xD_low = 0.0_rkind
	xD_high = 6.0_rkind
	epsD = 1.e-10_rkind

	yD_low = fun_real(xD_low)
	yD_high = fun_real(xD_high)

    write(*,*) ' '
    write(*,*) 'Test list_invert_explicit'

	do i = 1, npoints
		yD(i) = yD_low + (i-1)*(yD_high-yD_low)/(npoints-1)
	end do

	call list_invert(fun_real, npoints, yD,  xD_low, xD_high, epsD, xD)

	do i = 1, npoints
		y_calc(i) = fun_real(xD(i))
		x_calc(i) = sqrt(yD(i))
	end do

    WRITE (*, *) '     x  = ', xD
    WRITE (*, *) 'x_calc  = ', x_calc
    WRITE (*, *) '     y  = ', yD
    WRITE (*, *) 'y_calc  = ', y_calc

    write(*,*) ' '
    write(*,*) 'Test list_invert_uniform'
	call list_invert(fun_real, npoints,  xD_low, xD_high, epsD, xD)

	do i = 1, npoints
		yD(i) = yD_low + (i-1)*(yD_high-yD_low)/(npoints-1)
		y_calc(i) = fun_real(xD(i))
		x_calc(i) = sqrt(y_calc(i))
	end do

    WRITE (*, *) '     x  = ', xD
    WRITE (*, *) 'x_calc  = ', x_calc
    WRITE (*, *) '     y  = ', yD
    WRITE (*, *) 'y_calc  = ', y_calc

!     write(*,*) ' '
!     write(*,*) 'Test spline_invert_uniform'
! 	call spline_invert_uniform_init(sp_fun, fun_real, f_name, npoints, xD_low, xD_high, epsD)
!
! 	do i = 1, n_test
! 		yD(i) = fun_real(xD(i))
! 		y_test(i) = yD_low + (i-1)*(yD_high-yD_low)/(n_test-1)
! 		call sp_fun%eval_1D_f(y_test(i), x_test(i))
! 		x_calc(i) = sqrt(y_test(i))
! 		y_calc(i) = fun_real(x_calc(i))
! 	end do
!
!     WRITE (*, *) 'x_test = ', x_test
!     WRITE (*, *) 'x_calc  = ', x_calc
!     WRITE (*, *) 'y_test  = ', y_test
!     WRITE (*, *) 'y_calc  = ', y_calc

END PROGRAM fun_invert

!   function fun(x)
! 	IMPLICIT NONE
! 	REAL :: fun, x
! ! 	fun = x**2
! 	fun = sqrt(x)
! 	return
!   end function fun

  function fun_real(x)
	IMPLICIT NONE
    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind) fun_real, x

	fun_real = x**2
! 	fun_real = sqrt(x)
	return
  end function fun_real
