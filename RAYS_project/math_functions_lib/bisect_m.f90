 MODULE bisect_m

! generic procedure: solve_bisection() to sole an equation by the bisection method.
! This is slightly different in that it solves f(x)=y rather than f(x)=0
! Usual caveats for non-linear solvers apply.  May fail if inverse function is
! multi-valued in search domain.  There is no attempt to find multiple solutions.
!
! f = objective function y = f(x)
! x = solution
! x0_in = lower search bound
! x1_in = upper search bound
! y = equation rhs
! eps = error tolerance: ABS(f(x)-y) <= eps
! ierr return:
! ierr = -1 => (y(x1)-y(x0))(x1-x0) > 0 => no solution or multiple solutions in [x0,x1]
! ierr = 0 => no convergence after max_iter iterations
! ierr >= 1 => sucessful solution, ierr = number of iterations
!
! There are four routines under the solve_bisection interface. bisect_single, bisect_real
! which accept the arguments listed above. The calling sequence for these is:
!   call solve_bisection(f, x,  x0_in, x1_in, y, eps, ierr)
! The routines bisect_singleX and bisect_realX accept an additional argument,
! extra_data(:), which is an assumed shape vector.  The calling sequence for these is:
!   call solve_bisection(f, x,  x0_in, x1_in, y, eps, ierr, extra_data)
! The actual length of extra_data(:) must be declared in the code calling solve_bisection
! and also declared in f().
!
!______________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision
	INTEGER, PARAMETER :: max_iter = 1000

    interface solve_bisection
        module procedure bisect_single, bisect_real, bisect_singleX, bisect_realX
    end interface

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

  SUBROUTINE bisect_single(f, x,  x0_in, x1_in, y, eps, ierr)

	IMPLICIT NONE

	REAL(kind = skind),	EXTERNAL :: f
	REAL(kind = skind), INTENT(IN) :: x0_in
	REAL(kind = skind), INTENT(IN) :: x1_in
	REAL(kind = skind), INTENT(IN) :: y
	REAL(kind = skind), INTENT(IN) :: eps
	REAL(kind = skind), INTENT(out) :: x
	INTEGER, INTENT(out) :: ierr

! Local variables
	REAL(kind = skind) :: y0, y1, y2, x0, x1, x2
	INTEGER :: i

	x0 = x0_in
	x1 = x1_in

	DO i = 1, max_iter
		ierr = i
		y0 = f(x0) - y
		y1 = f(x1) - y

		if (abs(y0) <= eps) then
			x = x0
			exit
		end if

		if (abs(y1) <= eps) then
			x = x1
			exit
		end if

		if (y0*y1 > 0.) then
			ierr = -1
			exit
		end if

		x2 = (x0 + x1)/2.0
		y2 = f(x2) - y

		if(y0*y2<0)then
		   x1=x2
		else
		   x0=x2
		endif
	END DO

	ierr = i
	if (i >= max_iter)ierr = 0 ! No convergence in max_iter iterations

  RETURN
  END SUBROUTINE bisect_single

 !*********************************************************************************

  SUBROUTINE bisect_real(f, x,  x0_in, x1_in, y, eps, ierr)

	IMPLICIT NONE

	REAL(kind = rkind),	EXTERNAL :: f
	REAL(kind = rkind), INTENT(IN) :: x0_in
	REAL(kind = rkind), INTENT(IN) :: x1_in
	REAL(kind = rkind), INTENT(IN) :: y
	REAL(kind = rkind), INTENT(IN) :: eps
	REAL(kind = rkind), INTENT(out) :: x
	INTEGER, INTENT(out) :: ierr

	REAL(kind = rkind) :: y0, y1, y2, x0, x1, x2
	INTEGER :: i

	ierr = 0
	x0 = x0_in
	x1 = x1_in

	DO i = 1, max_iter
		ierr = i
		y0 = f(x0) - y
		y1 = f(x1) - y

		if (abs(y0) <= eps) then
			x = x0
			exit
		end if

		if (abs(y1) <= eps) then
			x = x1
			exit
		end if

		if (y0*y1 > 0.0_rkind) then
			ierr = -1
			exit
		end if

		x2 = (x0 + x1)/2.0_rkind
		y2 = f(x2) - y

		if(y0*y2<0)then
		   x1=x2
		else
		   x0=x2
		endif
	END DO

	ierr = i
	if (i >= max_iter)ierr = 0 ! No convergence in max_iter iterations

  RETURN
  END SUBROUTINE bisect_real

 !*********************************************************************************

  SUBROUTINE bisect_singleX(f, x,  x0_in, x1_in, y, eps, ierr, extra_data)

	IMPLICIT NONE

	REAL(kind = skind),	EXTERNAL :: f
	REAL(kind = skind), INTENT(IN) :: x0_in
	REAL(kind = skind), INTENT(IN) :: x1_in
	REAL(kind = skind), INTENT(IN) :: y
	REAL(kind = skind), INTENT(IN) :: eps
	REAL(kind = skind), INTENT(IN) :: extra_data(:)
	REAL(kind = skind), INTENT(out) :: x
	INTEGER, INTENT(out) :: ierr

! Local variables
	REAL(kind = skind) :: y0, y1, y2, x0, x1, x2
	INTEGER :: i

	x0 = x0_in
	x1 = x1_in

	DO i = 1, max_iter
		ierr = i
		y0 = f(x0, extra_data(:)) - y
		y1 = f(x1, extra_data(:)) - y

		if (abs(y0) <= eps) then
			x = x0
			exit
		end if

		if (abs(y1) <= eps) then
			x = x1
			exit
		end if

		if (y0*y1 > 0.) then
			ierr = -1
			exit
		end if

		x2 = (x0 + x1)/2.0
		y2 = f(x2, extra_data(:)) - y

		if(y0*y2<0)then
		   x1=x2
		else
		   x0=x2
		endif
	END DO

	ierr = i
	if (i >= max_iter)ierr = 0 ! No convergence in max_iter iterations

  RETURN
  END SUBROUTINE bisect_singleX

 !*********************************************************************************

  SUBROUTINE bisect_realX(f, x,  x0_in, x1_in, y, eps, ierr, extra_data)

	IMPLICIT NONE

	REAL(kind = rkind),	EXTERNAL :: f
	REAL(kind = rkind), INTENT(IN) :: x0_in
	REAL(kind = rkind), INTENT(IN) :: x1_in
	REAL(kind = rkind), INTENT(IN) :: y
	REAL(kind = rkind), INTENT(IN) :: eps
	REAL(kind = rkind), INTENT(IN) :: extra_data(:)
	REAL(kind = rkind), INTENT(out) :: x
	INTEGER, INTENT(out) :: ierr

	REAL(kind = rkind) :: y0, y1, y2, x0, x1, x2
	INTEGER :: i

	ierr = 0
	x0 = x0_in
	x1 = x1_in

	DO i = 1, max_iter
		ierr = i
		y0 = f(x0, extra_data(:)) - y
		y1 = f(x1, extra_data(:)) - y

		if (abs(y0) <= eps) then
			x = x0
			exit
		end if

		if (abs(y1) <= eps) then
			x = x1
			exit
		end if

		if (y0*y1 > 0.0_rkind) then
			ierr = -1
			exit
		end if

		x2 = (x0 + x1)/2.0_rkind
		y2 = f(x2, extra_data(:)) - y

		if(y0*y2<0)then
		   x1=x2
		else
		   x0=x2
		endif
	END DO

	ierr = i
	if (i >= max_iter)ierr = 0 ! No convergence in max_iter iterations

  RETURN
  END SUBROUTINE bisect_realX

  END MODULE bisect_m