  module oneD_spline_function_m

! Wrapper code for cspline.f90 and cspeval.f90 (themselves extracted from PPPL pspline).
! Provides a derived type, spline_function_1D, which contains the spline coefficients and
! setup data needed for cubic spline interpolation plus type bound module procedures for
! initialization and evaluation.  The is much more flexibility in the pspline routines than
! is kept here, particularly with respect to the grid, which is expected to be uniformly
! spaced, and boundary conditions, here taken to be not-a-knot.

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

	implicit none

    type spline_function_1D

		integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
		real(KIND=rkind), parameter :: zero = 0.0_rkind, one = 1.0_rkind

		integer :: nx ! Number of points to be splined
		real(KIND=rkind), allocatable ::  x(:) ! Grid
        real(KIND=rkind), allocatable ::  fspl(:,:)
	    character (len = 80) :: name

	contains
		procedure :: init, eval

    end type spline_function_1D

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

    subroutine init(this, nx_in, x_in, y_in, name_in )

		implicit none

		class (spline_function_1D) ::  this

		integer :: nx_in
		real(KIND=rkind) ::  x(nx_in), y(nx_in)

		character (len = 80) :: name_in

		integer :: i, j
		real(KIND=rkind) ::  f, fx, fy, fxx, fxy, fyy

		integer ibcxmin,ibcxmax,ilinx, ilinth,ier
		real(KIND=rkind) :: bcxmin(1),bcxmax(1)      ! Would be nx if used
		real(KIND=rkind) :: wk(1) ! Not used except for periodic BC
		real(KIND=rkind) :: fspl(4,nx_in)

		nx = nx_in
		name = name_in
		x = x_in
		y = y_in


	! cspline fixed args -> BC = not-a-knot
		ibcxmin = 0
		bcxmin = 0.
		ibcxmax = 0
		bcxmax = 0.
		ilinx = 0

		fspl1D = zero
		fspl1D(1,:) = y

		call cspline(x,nx,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,wk,1,ilinx,ier)
		write (*,*) 'ilinx = ', ilinx

    end subroutine init(this, nx_in, x_in, y_in, name_in )

end  module oneD_spline_function_m