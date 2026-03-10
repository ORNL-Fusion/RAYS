 MODULE newtonR3_m

! generic procedure: solve_newtonR3() to solve a scalar equation on R3 by
! Newtons method, f(x)=y. User must supply functions f(x) and gradf(x)
! Usual caveats for non-linear solvers apply. There is no attempt to find multiple
! solutions.
!
! f = objective function y = f(x)
! gradf = gradient of objective function f(x)
! x = solution, dimension = 3
! x0 starting point, dimension = 3
! step_max = upper limit on norm2 of length of x step to take
! y = equation rhs, constant
! eps = error tolerance: ABS(f(x)-y) <= eps
! nsig = Stopping criteria if successive iterations agree to nsig significant figures
! i_iter return:
! i_iter = -2 => error return, no convergence after max_iter iterations
! i_iter = -1 => error return, |gradf| = 0
! i_iter >= 0 => sucessful solution, i_iter = number of iterations
! max_iter_in, optional => can be used to override the default maximum number of iterations.
!
! There are four routines under the solve_newtonR3 interface.
! Routines newtonR3_single and newtonR3_real
! accept the arguments listed above. The calling sequence for these is:
!   call solve_newtonR3(f, x, x0, step_max, y, eps, i_iter, max_iter_in[keyword, optional])
!
! Routines newtonR3_singleX and newtonR3_realX
! accept an additional argument,extra_data(:), which is an assumed shape vector.
! The calling sequence for these is:
!   call solve_newtonR3(f, x,  x0, step_max, y, eps, i_iter, &
!      & extra_data, max_iter_in[keyword, optional])

! N.B. The actual length of extra_data(:) must be declared in the code calling
! solve_newtonR3 and also declared in f().

!______________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________
!
! (1/14/2026) I have not got the extra data routines to work in test_newtonR3.f90 to work
! yet.  It complains that there is no specific interface that matches the calls.  Figure
! it out later.
!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision
	INTEGER :: max_iter = 1000

    interface solve_newtonR3
        module procedure newtonR3_single, newtonR3_real, newtonR3_singleX, newtonR3_realX
    end interface

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

  SUBROUTINE newtonR3_single(f, gradf, x,  x_start, step_max, y, eps, nsig, i_iter,&
            & max_iter_in)

	IMPLICIT NONE

	interface
		function f(x)
			real, intent(in) :: x(3)
			real :: f
		end function
	end interface

	interface
		function gradf(x)
			real, intent(in) :: x(3)
			real :: gradf(3)
		end function
	end interface

! 	REAL(kind = skind),	EXTERNAL :: f
! 	REAL(kind = skind), EXTERNAL :: gradf
	REAL(kind = skind), INTENT(out) :: x(3)
	REAL(kind = skind), INTENT(IN) :: x_start(3)
	REAL(kind = skind), INTENT(IN) :: step_max
	REAL(kind = skind), INTENT(IN) :: y
	REAL(kind = skind), INTENT(IN) :: eps
	INTEGER, INTENT(IN) :: nsig
	INTEGER, INTENT(out) :: i_iter
	INTEGER, OPTIONAL, INTENT(in) :: max_iter_in

! Local variables
	REAL(kind = skind) :: x0(3), delta_x(3), fp(3)
	REAL(kind = skind) :: y0, y1, delta
	INTEGER :: i

	if (present(max_iter_in)) max_iter = max_iter_in

	i_iter = 0
	x0 = x_start
	y0 = f(x0)
	if (abs(y0-y) <= eps) then
		x = x0
		return
	end if

	DO i = 1, max_iter
!  write(*,*) 'newtonR3_single: x0 = ', x0
		i_iter = i
		y0 = f(x0)
		fp = gradf(x0)

		! Check that gradf is non-zero, otherwise error return
		if (norm2(fp) == 0.) then
			i_iter = -1
			x = x0
			return
		end if

		delta = -(f(x0) - y)/dot_product(fp,fp)
		delta_x = fp*delta

		! check if |delta_x| < 10**-nsig
		if (norm2(delta_x) < 10.**(-nsig)) then
			x = x0
			return
		end if

		! Iterate
		x0 = x0 + delta_x
		y1 = f(x0) - y

		! check if |f(x0) - y| < eps
		if (abs(y1) <= eps) then
			x = x0
			return
		end if

	END DO

	if (i >= max_iter)i_iter = -2 ! No convergence in max_iter iterations

  RETURN
  END SUBROUTINE newtonR3_single

 !*********************************************************************************

  SUBROUTINE newtonR3_real(f, gradf, x,  x_start, step_max, y, eps, nsig, i_iter,&
           & max_iter_in)

	IMPLICIT NONE

	interface
		function f(x)
  			integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
			real(kind = rkind), intent(in) :: x(3)
			real(kind = rkind) :: f
		end function
	end interface

	interface
		function gradf(x)
		    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
			real(kind = rkind), intent(in) :: x(3)
			real(kind = rkind) :: gradf(3)
		end function
	end interface

	REAL(kind = rkind), INTENT(out) :: x(3)
	REAL(kind = rkind), INTENT(IN) :: x_start(3)
	REAL(kind = rkind), INTENT(IN) :: step_max
	REAL(kind = rkind), INTENT(IN) :: y
	REAL(kind = rkind), INTENT(IN) :: eps
	INTEGER, INTENT(IN) :: nsig
	INTEGER, INTENT(out) :: i_iter
	INTEGER, OPTIONAL, INTENT(in) :: max_iter_in

! Local variables
	REAL(kind = rkind) :: x0(3), delta_x(3), fp(3)
	REAL(kind = rkind) :: y0, y1, delta
	INTEGER :: i

	if (present(max_iter_in)) max_iter = max_iter_in

	i_iter = 0
	x0 = x_start
	y0 = f(x0)
	if (abs(y0-y) <= eps) then
		x = x0
		return
	end if

	DO i = 1, max_iter
		i_iter = i
		y0 = f(x0)
		fp = gradf(x0)

		! Check that gradf is non-zero, otherwise error return
		if (norm2(fp) == 0.) then
			i_iter = -1
			x = x0
			return
		end if

		delta = -(f(x0) - y)/dot_product(fp,fp)
		delta_x = fp*delta

		! check if |delta_x| < 10**-nsig
		if (norm2(delta_x) < 10.**(-nsig)) then
			x = x0
			return
		end if

		! Iterate
		x0 = x0 + delta_x
		y1 = f(x0) - y

		! check if |f(x0) - y| < eps
		if (abs(y1) <= eps) then
			x = x0
			exit
		end if

	END DO

	if (i >= max_iter)i_iter = -2 ! No convergence in max_iter iterations

  RETURN
  END SUBROUTINE newtonR3_real

!_________________________________________________________________________________________

  SUBROUTINE newtonR3_singleX(f, gradf, x,  x_start, step_max, y, eps, nsig, i_iter,&
         & extra_data, max_iter_in)

	IMPLICIT NONE

	interface
		function f(x, extra_data)
			real, intent(in) :: x(3)
			real, intent(in) :: extra_data(:)
			real :: f
		end function
	end interface

	interface
		function gradf(x, extra_data)
			real, intent(in) :: x(3)
			real, intent(in) :: extra_data(:)
			real :: gradf(3)
		end function
	end interface

	REAL(kind = skind), INTENT(out) :: x(3)
	REAL(kind = skind), INTENT(IN) :: x_start(3)
	REAL(kind = skind), INTENT(IN) :: step_max
	REAL(kind = skind), INTENT(IN) :: y
	REAL(kind = skind), INTENT(IN) :: eps
	INTEGER, INTENT(IN) :: nsig
	INTEGER, INTENT(out) :: i_iter
	REAL(kind = skind), INTENT(IN) :: extra_data(:)
	INTEGER, OPTIONAL, INTENT(in) :: max_iter_in

! Local variables
	REAL(kind = skind) :: x0(3), delta_x(3), fp(3)
	REAL(kind = skind) :: y0, y1, delta
	INTEGER :: i

	if (present(max_iter_in)) max_iter = max_iter_in

	i_iter = 0
	x0 = x_start
	y0 = f(x0, extra_data(:))
	if (abs(y0-y) <= eps) then
		x = x0
		return
	end if

	DO i = 1, max_iter
!  write(*,*) 'newtonR3_singleX: x0 = ', x0
		i_iter = i
		y0 = f(x0, extra_data(:))
		fp = gradf(x0, extra_data(:))

		! Check that gradf is non-zero, otherwise error return
		if (norm2(fp) == 0.) then
			i_iter = -1
			x = x0
			return
		end if

		delta = -(f(x0, extra_data(:)) - y)/dot_product(fp,fp)
		delta_x = fp*delta

		! check if |delta_x| < 10**-nsig
		if (norm2(delta_x) < 10.**(-nsig)) then
			x = x0
			return
		end if

		! Iterate
		x0 = x0 + delta_x
		y1 = f(x0, extra_data(:)) - y

		! check if |f(x0) - y| < eps
		if (abs(y1) <= eps) then
			x = x0
			return
		end if

	END DO

	if (i >= max_iter)i_iter = -2 ! No convergence in max_iter iterations

  RETURN
  END SUBROUTINE newtonR3_singleX

 !*********************************************************************************

  SUBROUTINE newtonR3_realX(f, gradf, x,  x_start, step_max, y, eps, nsig, i_iter,&
         & extra_data, max_iter_in)

	IMPLICIT NONE

	interface
		function f(x, extra_data)
  			integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
			real(kind = rkind), intent(in) :: x(3)
			real(kind = rkind), intent(in) :: extra_data(:)
			real(kind = rkind) :: f
		end function
	end interface

	interface
		function gradf(x, extra_data)
		    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
			real(kind = rkind), intent(in) :: x(3)
			real(kind = rkind), intent(in) :: extra_data(:)
			real(kind = rkind) :: gradf(3)
		end function
	end interface

	REAL(kind = rkind), INTENT(out) :: x(3)
	REAL(kind = rkind), INTENT(IN) :: x_start(3)
	REAL(kind = rkind), INTENT(IN) :: step_max
	REAL(kind = rkind), INTENT(IN) :: y
	REAL(kind = rkind), INTENT(IN) :: eps
	INTEGER, INTENT(IN) :: nsig
	INTEGER, INTENT(out) :: i_iter
	REAL(kind = rkind), INTENT(IN) :: extra_data(:)
	INTEGER, OPTIONAL, INTENT(in) :: max_iter_in

! Local variables
	REAL(kind = rkind) :: x0(3), delta_x(3), fp(3)
	REAL(kind = rkind) :: y0, y1, delta
	INTEGER :: i

	if (present(max_iter_in)) max_iter = max_iter_in

	i_iter = 0
	x0 = x_start
	y0 = f(x0, extra_data(:))
	if (abs(y0-y) <= eps) then
		x = x0
		return
	end if

	DO i = 1, max_iter
		i_iter = i
		y0 = f(x0, extra_data(:))
		fp = gradf(x0, extra_data(:))

		! Check that gradf is non-zero, otherwise error return
		if (norm2(fp) == 0.) then
			i_iter = -1
			x = x0
			return
		end if

		delta = -(f(x0, extra_data(:)) - y)/dot_product(fp,fp)
		delta_x = fp*delta

		! check if |delta_x| < 10**-nsig
		if (norm2(delta_x) < 10.**(-nsig)) then
			x = x0
			return
		end if

		! Iterate
		x0 = x0 + delta_x
		y1 = f(x0, extra_data(:)) - y

		! check if |f(x0) - y| < eps
		if (abs(y1) <= eps) then
			x = x0
			exit
		end if

	END DO

	if (i >= max_iter)i_iter = -2 ! No convergence in max_iter iterations

  RETURN
  END SUBROUTINE newtonR3_realX

  END MODULE newtonR3_m