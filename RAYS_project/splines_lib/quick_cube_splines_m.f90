  module quick_cube_splines_m

! Wrapper code for cspline.f90 and cspeval.f90 (themselves extracted from PPPL pspline).
! Provides two derived types, spline_function_1D, which contains the spline coefficients and
! setup data needed for cubic spline interpolation, plus type bound procedures for
! initialization and evaluation.  There is much more flexibility in the pspline routines
! than is kept here, particularly with respect to the grid, which is expected to be uniformly
! spaced, and boundary conditions, here taken to be not-a-knot.
!
! Exports derived type -> spline_function_1D, and 4 subroutines:
! cube_spline_1D_init -> Initializes spline coefficients for 'this' instance
! eval_1D_f(this, x, f) -> evaluates splined function f
! eval_1D_fp(this, x, f, fp) -> evaluates splined function f plus 1st derivative
! eval_1D_fpp(this, x, f, fp, fpp) -> evaluates splined function f plus 1st and 2nd derivatives
!
! Exports derived type -> spline_function_2D, and 4 subroutines:
! cube_spline_2D_init -> Initializes spline coefficients for 'this' instance
! eval_2D_f(this, x, f) -> evaluates splined function f
! eval_2D_fp(this, x, f, fx, fy) -> evaluates splined function f plus 1st derivative
! eval_2D_fpp(this, x, y, f, fx, fy, fxx, fxy, fyy) -> evaluates splined function f plus
! 1st and 2nd derivatives
	implicit none

	integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
	real(KIND=rkind), parameter :: zero = 0.0_rkind, one = 1.0_rkind

    type cube_spline_function_1D
		integer :: nx ! Number of points to be splined
		real(KIND=rkind), allocatable ::  x_grid(:) ! Grid
		real(KIND=rkind), allocatable ::  f_val(:) ! Values on grid
        real(KIND=rkind), allocatable ::  fspl(:,:) ! Spline coeffieients
	    character (len = 80) :: name ! Function name, useful for arrays of this type
	contains
		procedure :: cube_spline_1D_init
		procedure :: eval_1D_f, eval_1D_fp, eval_1D_fpp
    end type cube_spline_function_1D

    type cube_spline_function_2D
		integer :: nx, ny ! Number of points to be splined
		real(KIND=rkind), allocatable ::  x_grid(:), y_grid(:) ! Grid
		real(KIND=rkind), allocatable ::  f_val(:,:) ! Values on grid
        real(KIND=rkind), allocatable ::  fspl(:,:, :, :)
	    character (len = 80) :: name ! Function name, useful for arrays of this type
	contains
		procedure :: cube_spline_2D_init
		procedure :: eval_2D_f, eval_2D_fp, eval_2D_fpp
    end type cube_spline_function_2D

!****************************************************************************

contains

!****************************************************************************
! 1D Stuff

    subroutine cube_spline_1D_init(this, nx_in, x_in, f_in, name_in)

		implicit none

		class (cube_spline_function_1D) ::  this

		integer :: nx_in
		real(KIND=rkind) ::  x_in(nx_in), f_in(nx_in)

		character (len = 80) :: name_in

		integer ibcxmin,ibcxmax,ilinx, ilinth,ier
		real(KIND=rkind) :: bcxmin(1),bcxmax(1)      ! Not used, would be nx if used
		real(KIND=rkind) :: wk(1) ! Not used except for periodic BC
		real(KIND=rkind) :: fspl(4,nx_in)

		this%nx = nx_in
		this%x_grid = x_in
		this%f_val = f_in
		this%name = name_in

	! cspline fixed args -> BC = not-a-knot
		ibcxmin = 0
		bcxmin = 0.
		ibcxmax = 0
		bcxmax = 0.
		ilinx = 0

		allocate(this%fspl(4,nx_in), source = zero)
		this%fspl(1,:) = f_in

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

    end subroutine cube_spline_1D_init

!****************************************************************************

    subroutine eval_1D_f(this, x, f)

		implicit none

		class (cube_spline_function_1D), intent(in) ::  this

		real(KIND=rkind), intent(in) ::  x
		real(KIND=rkind), intent(out) ::  f

		integer :: ier
		integer :: ilinx = 1
        integer :: iselect(3) = (/1,0,0/)  !output selector
		real(KIND=rkind) ::  f_eval(3)

        call cspeval(x, iselect, f_eval, this%x_grid, this%nx, ilinx, this%fspl, ier)
		f = f_eval(1)

		return
    end subroutine eval_1D_f

!****************************************************************************

    module subroutine eval_1D_fp(this, x, f, fp)

		implicit none

		class (cube_spline_function_1D), intent(in) ::  this

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
    end subroutine eval_1D_fp

!****************************************************************************


    module subroutine eval_1D_fpp(this, x, f, fp, fpp)

		implicit none

		class (cube_spline_function_1D), intent(in) ::  this

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

    end subroutine eval_1D_fpp

!****************************************************************************
! 2D Stuff
!****************************************************************************

    subroutine cube_spline_2D_init(this, nx_in, x_in, ny_in, y_in, f_in, name_in)

		implicit none

		class (cube_spline_function_2D) ::  this

		integer :: nx_in, ny_in
		real(KIND=rkind) ::  x_in(nx_in),y_in(ny_in), f_in(nx_in, ny_in)

		character (len = 80) :: name_in


		! bcspline arguments
		integer :: ibcxmin, ibcxmax, ibcymin, ibcymax, ilinx, iliny, nwk, ier
		real(KIND=rkind), allocatable :: fspl(:,:,:,:), wk(:)
		real(KIND=rkind):: bcxmin(1),bcxmax(1)  ! Not used, would be nx if used
		real(KIND=rkind) :: bcymin(1),bcymax(1)  ! Not used, would be ny if used
		real(KIND=rkind) :: fval(6)

		integer :: i, j

		this%nx = nx_in
		this%x_grid = x_in
		this%ny = ny_in
		this%y_grid = y_in
		this%f_val = f_in
		this%name = name_in

		allocate (this%fspl(4,4,nx_in,ny_in))
	    nwk = 4*nx_in*ny_in +5*max(nx_in,ny_in)
		allocate (wk(nwk))

		! cspline fixed args -> BC = not-a-knot
		ibcxmin = 0
		bcxmin = 0.
		ibcxmax = 0
		bcxmax = 0.
		ibcymin = 0
		bcymin = 0.
		ibcymax = 0
		bcymax = 0.
		ilinx = 1
		iliny = 1

	! Set f(1,1,i,j) = f_in
		do i = 1, nx_in; do j = 1, ny_in
			this%fspl(1,1,i,j) = f_in(i,j)
		end do; end do

		call bcspline(x_in,nx_in,y_in,ny_in,this%fspl,nx_in, &
			 ibcxmin,bcxmin,ibcxmax,bcxmax, ibcymin,bcymin,ibcymax,bcymax, &
			 wk,nwk,ilinx,iliny,ier)

 		if (ilinx /= 1 .or. ilinx /= 1) then
			write (*,*) 'cube_spline_2D_init: init, grid not evenly spaced'
			stop
		end if
		if (ier /= 0) then
			write (*,*) 'cube_spline_2D_init: init, error return ier = ', ier
			stop
		end if

		return

    end subroutine cube_spline_2D_init

!****************************************************************************

    subroutine eval_2D_f(this, x, y, f)

		implicit none

		class (cube_spline_function_2D), intent(in) ::  this

		real(KIND=rkind), intent(in) ::  x, y
		real(KIND=rkind), intent(out) ::  f

		integer :: ilinx = 1, iliny = 1, ier
        integer :: iselect(6) = (/1, 0, 0, 0, 0, 0/)  !output selector
		real(KIND=rkind) ::  fval(6)

	! evaluate spline fits
		call bcspeval(x,y,iselect,fval,this%x_grid, this%nx,this%y_grid, this%ny,&
			 & ilinx,iliny,this%fspl,this%nx,ier)
		if (ier .ne. 0) write (*,*) 'bcspeval: ier = ', ier

		f = fval(1)

		return
    end subroutine eval_2D_f

!****************************************************************************

    module subroutine eval_2D_fp(this, x, y, f, fx, fy)

		implicit none

		class (cube_spline_function_2D), intent(in) ::  this

		real(KIND=rkind), intent(in) ::  x, y
		real(KIND=rkind), intent(out) ::  f, fx, fy

		integer :: ilinx = 1, iliny = 1, ier
        integer :: iselect(6) = (/1, 1, 1, 0, 0, 0/)  !output selector
		real(KIND=rkind) ::  fval(6)

	! evaluate spline fits
		call bcspeval(x,y,iselect,fval,this%x_grid, this%nx,this%y_grid, this%ny,&
			 & ilinx,iliny,this%fspl,this%nx,ier)
		if (ier .ne. 0) write (*,*) 'bcspeval: ier = ', ier

		f = fval(1)
		fx = fval(2)
		fy = fval(3)

		return
    end subroutine eval_2D_fp

!****************************************************************************


    module subroutine eval_2D_fpp(this, x, y, f, fx, fy, fxx, fxy, fyy)

		implicit none

		class (cube_spline_function_2D), intent(in) ::  this

		real(KIND=rkind), intent(in) ::  x, y
		real(KIND=rkind), intent(out) ::  f, fx, fy, fxx, fxy, fyy

		integer :: ilinx = 1, iliny = 1, ier
        integer :: iselect(6) = (/1, 1, 1, 1, 1, 1/)  !output selector
		real(KIND=rkind) ::  fval(6)

	! evaluate spline fits
		call bcspeval(x,y,iselect,fval,this%x_grid, this%nx,this%y_grid, this%ny,&
			 & ilinx,iliny,this%fspl,this%nx,ier)
		if (ier .ne. 0) write (*,*) 'bcspeval: ier = ', ier

		f = fval(1)
		fx = fval(2)
		fy = fval(3)
		fxx = fval(4)
		fyy = fval(5)
		fxy = fval(6)

		return

    end subroutine eval_2D_fpp

end  module quick_cube_splines_m