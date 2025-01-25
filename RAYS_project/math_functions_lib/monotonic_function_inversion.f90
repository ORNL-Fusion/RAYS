 MODULE mono_funct_inversion_m

! generic procedures: to invert a monotonic function y(x) -> x(y).
! Exports two functions:
!
! list_invert: For an external function f(x), accepts a strictly increasing list
! of y values, yi, returns a list of x values, xi, such that f(xi) = yi.
!
! spline_invert: For an external function f(x), takes a y list, uses list_invert
! above to produce an inverse x list, then returns a cube_spline_function_1D derived
! type object, obtained from module quick_cube_splines_m, that approximates x(y)
!
! More on the exported routines:
!
! list_invert is the generic interface for two module procedures: list_invert_explicit
! and list_invert_uniform:
!
! list_invert_explicit() takes an explicit list of y values to solve for x within
! the search domain [x_min, x_max]. Values in the y list must be strictly monotonic
! and must lie within the range [f(x_min), f(x_max)]
!
! list_invert_uniform constructs a uniform y grid from [y(xmin), y(x_max)] and inverts
! on that grid

! f = objective function y = f(x)
! x = solution
! x_min = lower search bound
! x_max = upper search bound
! y = equation rhs
! eps = error tolerance: ABS(f(x)-y) <= eps

    use bisect_m, only : solve_bisection
!     use quick_cube_splines_m, only : cube_spline_function_1D

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

    interface list_invert
        module procedure list_invert_explicit, list_invert_uniform
    end interface

    CONTAINS
 !*********************************************************************************

  SUBROUTINE list_invert_explicit(f, npoints, y_in,  x_min, x_max, eps, x)

	IMPLICIT NONE

	REAL(kind = rkind),	EXTERNAL :: f
	INTEGER, INTENT(IN) :: npoints
	REAL(kind = rkind), INTENT(in) :: y_in(npoints)
	REAL(kind = rkind), INTENT(IN) :: x_min
	REAL(kind = rkind), INTENT(IN) :: x_max
	REAL(kind = rkind), INTENT(IN) :: eps
	REAL(kind = rkind), INTENT(out) :: x(npoints)

! Local variables
	REAL(kind = rkind) :: y(npoints)
	REAL(kind = rkind) :: y0, y1, y2, x0, x1, x2
	INTEGER :: i, ierr
	LOGICAL :: increasing

	x0 = x_min
	x1 = x_max
	y0 = f(x0)
	y1 = f(x1)

	increasing = .true.
	if (y0 > y1) increasing = .false.

	DO i = 1, npoints
		call solve_bisection(f, x(i), x0, x1, y_in(i), eps, ierr)

		if (ierr < 1) then
			write(*,*) 'list_invert_explicit: bisection failed, y(i) = ', y(i), &
			&  'ierr = ',ierr
			stop
		end if

		if (increasing .eqv. .true.) then
			x0 = x(i)
		else
			x1 = x(i)
		end if
	END DO

  RETURN
  END SUBROUTINE list_invert_explicit


 !*********************************************************************************

  SUBROUTINE list_invert_uniform(f, npoints,  x_min, x_max, eps, x)

	IMPLICIT NONE

	REAL(kind = rkind),	EXTERNAL :: f
	INTEGER, INTENT(IN) :: npoints
	REAL(kind = rkind), INTENT(IN) :: x_min
	REAL(kind = rkind), INTENT(IN) :: x_max
	REAL(kind = rkind), INTENT(IN) :: eps
	REAL(kind = rkind), INTENT(out) :: x(npoints)

! Local variables
	REAL(kind = rkind) :: y(npoints)
	REAL(kind = rkind) :: y0, y1, y2, x0, x1, x2
	INTEGER :: i, ierr
	LOGICAL :: increasing

	x0 = x_min
	x1 = x_max
	y0 = f(x0)
	y1 = f(x1)

	increasing = .true.
	if (y0 > y1) increasing = .false.

	DO i = 1, npoints
		y(i) = y0 + (i-1)*(y1-y0)/(npoints-1)
		call solve_bisection(f, x(i), x0, x1, y(i), eps, ierr)

		if (ierr < 1) then
			write(*,*) 'list_invert_uniform: bisection failed, y(i) = ', y(i), &
			&  'ierr = ',ierr
			stop
		end if

		if (increasing .eqv. .true.) then
			x0 = x(i)
		else
			x1 = x(i)
		end if
	END DO

  RETURN
  END SUBROUTINE list_invert_uniform


 !*********************************************************************************

  SUBROUTINE spline_invert_uniform_init(sp_fun, f, f_name, npoints, x_min, x_max, eps)

    use quick_cube_splines_m, only : cube_spline_function_1D

	IMPLICIT NONE

    TYPE(cube_spline_function_1D), INTENT(out) :: sp_fun
	REAL(kind = rkind),	EXTERNAL :: f
 	CHARACTER (len = 80), INTENT(IN) :: f_name
	INTEGER, INTENT(IN) :: npoints
	REAL(kind = rkind), INTENT(IN) :: x_min
	REAL(kind = rkind), INTENT(IN) :: x_max
	REAL(kind = rkind), INTENT(IN) :: eps

	REAL(kind = rkind) :: y_low, y_high
	REAL(kind = rkind) :: y_grid(npoints), x(npoints)
	integer :: i

	call list_invert_uniform(f, npoints,  x_min, x_max, eps, x)

! Generate uniform grid of y values and initialize spline coefficients for sp_fun
	y_low = f(x_min)
	y_high = f(x_max)
	do i = 1, npoints
		y_grid(i) = y_low + (i-1)*(y_high-y_low)/(npoints-1)
	end do

    call sp_fun%cube_spline_1D_init(npoints, y_grid, x, f_name)

  RETURN
  END SUBROUTINE spline_invert_uniform_init

  END MODULE mono_funct_inversion_m