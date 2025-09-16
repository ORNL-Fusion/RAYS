  module oneD_spline_function_m

! Wrapper code for cspline.f90 and cspeval.f90 (themselves extracted from PPPL pspline).
! Provides a derived type, spline_function_1D, which contains the spline coefficients and
! setup data needed for cubic spline interpolation, plus type bound procedures for
! initialization and evaluation.  There is much more flexibility in the pspline routines
! than is kept here, particularly with respect to the grid, which is expected to be uniformly
! spaced, and boundary conditions, here taken to be not-a-knot.
!
! Exports one derived type -> spline_function_1D, and 4 subroutines:
! spline_1D_init -> Initializes spline coefficients for 'this' instance
! eval_f(this, x, f) -> evaluates splined function f
! eval_fp(this, x, f, fp) -> evaluates splined function f plus 1st derivative
! eval_fpp(this, x, f, fp, fpp) -> evaluates splined function f plus 1st ands 2nd derivatives
! Could bundle this into a generic eval procedure but I think it is clearer to expose the
! arguments.

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

	implicit none

	integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
	real(KIND=rkind), parameter :: zero = 0.0_rkind, one = 1.0_rkind

    type spline_function_1D

		integer :: nx ! Number of points to be splined
		real(KIND=rkind), allocatable ::  x_grid(:) ! Grid
		real(KIND=rkind), allocatable ::  y(:) ! Values on grid
        real(KIND=rkind), allocatable ::  fspl(:,:)
	    character (len = 80) :: name ! Function name, useful for arrays of this type

	contains
		procedure :: spline_1D_init
		procedure :: eval_f, eval_fp, eval_fpp

    end type spline_function_1D

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

    subroutine spline_1D_init(this, nx_in, x_in, y_in, name_in)

		implicit none

		class (spline_function_1D) ::  this

		integer :: nx_in
		real(KIND=rkind) ::  x_in(nx_in), y_in(nx_in)

		character (len = 80) :: name_in

		integer ibcxmin,ibcxmax,ilinx, ilinth,ier
		real(KIND=rkind) :: bcxmin(1),bcxmax(1)      ! Not used, would be nx if used
		real(KIND=rkind) :: wk(1) ! Not used except for periodic BC
		real(KIND=rkind) :: fspl(4,nx_in)


		this%nx = nx_in
		this%name = name_in
		this%x_grid = x_in
		this%y = y_in

	! cspline fixed args -> BC = not-a-knot
		ibcxmin = 0
		bcxmin = 0.
		ibcxmax = 0
		bcxmax = 0.
		ilinx = 0

		allocate(this%fspl(4,nx_in), source = zero)
		this%fspl(1,:) = y_in

		call cspline(this%x_grid,this%nx,this%fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,wk,1,ilinx,ier)
		if (ilinx == 2) then
			write (*,*) 'spline_function_1D: init, grid not evenly spaced'
			stop
		end if
		if (ier /= 0) then
			write (*,*) 'spline_function_1D: init, error return ier = ', ier
			stop
		end if

		return

    end subroutine spline_1D_init

!****************************************************************************

    subroutine eval_f(this, x, f)

		implicit none

		class (spline_function_1D), intent(in) ::  this

		real(KIND=rkind), intent(in) ::  x
		real(KIND=rkind), intent(out) ::  f

		integer :: ier
		integer :: ilinx = 1
        integer :: iselect(3) = (/1,0,0/)  !output selector
		real(KIND=rkind) ::  f_eval(3)

        call cspeval(x, iselect, f_eval, this%x_grid, this%nx, ilinx, this%fspl, ier)
		f = f_eval(1)

		return
    end subroutine eval_f

!****************************************************************************

    module subroutine eval_fp(this, x, f, fp)

		implicit none

		class (spline_function_1D), intent(in) ::  this

		real(KIND=rkind), intent(in) ::  x
		real(KIND=rkind), intent(out) ::  f, fp

		integer :: ier
		integer :: ilinx = 1
        integer :: iselect(3) = (/1,1,0/)   ! output selector
		real(KIND=rkind) ::  f_eval(3)

        call cspeval(x, iselect, f_eval, this%x_grid, this%nx, ilinx, this%fspl, ier)
		f = f_eval(1)
		fp = f_eval(2)

		return
    end subroutine eval_fp

!****************************************************************************


    module subroutine eval_fpp(this, x, f, fp, fpp)

		implicit none

		class (spline_function_1D), intent(in) ::  this

		real(KIND=rkind), intent(in) ::  x
		real(KIND=rkind), intent(out) ::  f, fp, fpp

		integer :: ier
		integer :: ilinx = 1
        integer :: iselect(3) = (/1,1,1/)   ! output selector
		real(KIND=rkind) ::  f_eval(3)

        call cspeval(x, iselect, f_eval, this%x_grid, this%nx, ilinx, this%fspl, ier)
		f = f_eval(1)
		fp = f_eval(2)
		fpp = f_eval(3)

		return

    end subroutine eval_fpp

end  module oneD_spline_function_m