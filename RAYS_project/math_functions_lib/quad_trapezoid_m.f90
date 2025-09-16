 MODULE trapezoid_quad_m

! generic procedure: trapezoid_quad, calculates integral from a vector of y
! values, y(:), using simple trapezoidal approximation

! generic procedure: trapezoid_quad_cumulative, calculates integral from a vector of y
! values, y(:), using simple trapezoidal approximation and also returns the cumulative
! integral at each grid point in a vector, cum(:) of length size(y).

! Input x is either a fixed increment, h, or a vector of x(:) values (not necessarily
! evenly spaced). x,y vector value can be real or complex, single precision or kind = rkind.

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

    interface trapezoid_quad
        module procedure h_yvec_real, h_yvec_D, h_yvec_cmplx, h_yvec_cmplxD,&
                       & xvec_yvec_real, xvec_yvec_realD, xvec_yvec_cmplx, xvec_yvec_cmplxD
    end interface

    interface trapezoid_quad_cumulative
        module procedure h_yvec_cum_real, h_yvec_cum_D,&
                       & xvec_yvec_cum_real, xvec_yvec_cum_D
    end interface

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

    SUBROUTINE h_yvec_real(h, y, integral)

        IMPLICIT NONE
        REAL(kind = skind), INTENT(IN) :: h, y(:)
        REAL(kind = skind), INTENT(out) :: integral

        integral = h*(SUM(y(2:SIZE(y))) + y(1)/2. - y(SIZE(y))/2.)

    RETURN
    END SUBROUTINE h_yvec_real
 !*********************************************************************************

    SUBROUTINE h_yvec_D(h, y, integral)

        IMPLICIT NONE
        REAL(kind = rkind), INTENT(IN) :: h, y(:)
        REAL(kind = rkind), INTENT(out) :: integral

        integral = h*(SUM(y(2:SIZE(y))) + y(1)/2. - y(SIZE(y))/2.)

    RETURN
    END SUBROUTINE h_yvec_D
 !*********************************************************************************


    SUBROUTINE h_yvec_cmplx(h, y, integral)

        IMPLICIT NONE
        REAL(kind = skind), INTENT(IN) :: h
        COMPLEX(kind = skind), INTENT(IN) :: y(:)
        COMPLEX(kind = skind), INTENT(out) :: integral

        integral = h*(SUM(y(2:SIZE(y))) + y(1)/2. - y(SIZE(y))/2.)

    RETURN
    END SUBROUTINE h_yvec_cmplx
 !*********************************************************************************

    SUBROUTINE h_yvec_cmplxD(h, y, integral)

        IMPLICIT NONE
        REAL(kind = rkind), INTENT(IN) :: h
        COMPLEX(kind = rkind), INTENT(IN) :: y(:)
        COMPLEX(kind = rkind), INTENT(out) :: integral

        integral = h*(SUM(y(2:SIZE(y))) + y(1)/2. - y(SIZE(y))/2.)

    RETURN
    END SUBROUTINE h_yvec_cmplxD
 !*********************************************************************************

    SUBROUTINE xvec_yvec_real(x, y, integral)

        IMPLICIT NONE
        REAL(kind = skind), INTENT(IN) :: x(:), y(:)
        REAL(kind = skind), INTENT(out) :: integral
        integer :: i

        if (size(x) /= size(y)) then
            write(*,*) 'quad_trapezoid: size(x) != size(y)'
            stop
        end if

        integral = 0.
        do i = 1, size(x)-1
        integral = integral + (x(i+1)-x(i))*(y(i+1)+y(i))
        end do
        integral = integral/2.

    RETURN
    END SUBROUTINE xvec_yvec_real
 !*********************************************************************************

    SUBROUTINE xvec_yvec_realD(x, y, integral)

        IMPLICIT NONE
        REAL(kind = rkind), INTENT(IN) :: x(:), y(:)
        REAL(kind = rkind), INTENT(out) :: integral
        integer :: i

        if (size(x) /= size(y)) then
            write(*,*) 'quad_trapezoid: size(x) != size(y)'
            stop
        end if

        integral = 0.
        do i = 1, size(x)-1
        integral = integral + (x(i+1)-x(i))*(y(i+1)+y(i))
        end do
        integral = integral/2.

    RETURN
    END SUBROUTINE xvec_yvec_realD
 !*********************************************************************************

    SUBROUTINE xvec_yvec_cmplx(x, y, integral)

        IMPLICIT NONE
        REAL(kind = skind), INTENT(IN) :: x(:)
        complex(kind = skind), INTENT(IN) :: y(:)
        complex(kind = skind), INTENT(out) :: integral
        integer :: i

        if (size(x) /= size(y)) then
            write(*,*) 'quad_trapezoid: size(x) != size(y)'
            stop
        end if

        integral = 0.
        do i = 1, size(x)-1
        integral = integral + (x(i+1)-x(i))*(y(i+1)+y(i))
        end do
        integral = integral/2.

    RETURN
    END SUBROUTINE xvec_yvec_cmplx
 !*********************************************************************************

    SUBROUTINE xvec_yvec_cmplxD(x, y, integral)

        IMPLICIT NONE
        REAL(kind = rkind), INTENT(IN) :: x(:)
        complex(kind = rkind), INTENT(IN) :: y(:)
        complex(kind = rkind), INTENT(out) :: integral
        integer :: i

        if (size(x) /= size(y)) then
            write(*,*) 'quad_trapezoid: size(x) != size(y)'
            stop
        end if

        integral = 0.
        do i = 1, size(x)-1
        integral = integral + (x(i+1)-x(i))*(y(i+1)+y(i))
        end do
        integral = integral/2.

    RETURN
    END SUBROUTINE xvec_yvec_cmplxD
 !*********************************************************************************

    SUBROUTINE h_yvec_cum_real(h, y, integral, cum)

    IMPLICIT NONE
    REAL, INTENT(IN) :: h, y(:)
    REAL, INTENT(out) :: integral
    REAL, INTENT(out) :: cum(SIZE(y))
    INTEGER  :: i

    cum(1) = 0.
    DO i = 2, SIZE(y)
        cum(i) = cum(i-1) + (y(i-1) + y(i))/2.
    END DO

    cum = h*cum
    integral = cum(SIZE(y))

    RETURN
    END SUBROUTINE h_yvec_cum_real
 !*********************************************************************************

    SUBROUTINE h_yvec_cum_D(h, y, integral, cum)

    IMPLICIT NONE
    REAL(kind = rkind), INTENT(IN) :: h, y(:)
    REAL(kind = rkind), INTENT(out) :: integral
    REAL(kind = rkind), INTENT(out) :: cum(SIZE(y))
    INTEGER  :: i

    cum(1) = 0.
    DO i = 2, SIZE(y)
        cum(i) = cum(i-1) + (y(i-1) + y(i))/2.
    END DO

    cum = h*cum
    integral = cum(SIZE(y))

    RETURN
    END SUBROUTINE h_yvec_cum_D
!*********************************************************************************

    SUBROUTINE xvec_yvec_cum_real(x, y, integral, cum)

        IMPLICIT NONE
        REAL(kind = skind), INTENT(IN) :: x(:), y(:)
        REAL(kind = skind), INTENT(out) :: integral
		REAL(kind = skind), INTENT(out) :: cum(SIZE(y))
        integer :: i

        if (size(x) /= size(y)) then
            write(*,*) 'quad_trapezoid: size(x) != size(y)'
            stop
        end if

		cum(1) = 0.
		DO i = 2, SIZE(y)
			cum(i) = cum(i-1) + (x(i)-x(i-1))*(y(i-1) + y(i))/2.
		END DO

	    integral = cum(SIZE(y))

    RETURN
    END SUBROUTINE xvec_yvec_cum_real
!*********************************************************************************

    SUBROUTINE xvec_yvec_cum_D(x, y, integral, cum)

        IMPLICIT NONE
        REAL(kind = rkind), INTENT(IN) :: x(:), y(:)
        REAL(kind = rkind), INTENT(out) :: integral
		REAL(kind = rkind), INTENT(out) :: cum(SIZE(y))
        integer :: i

        if (size(x) /= size(y)) then
            write(*,*) 'quad_trapezoid: size(x) != size(y)'
            stop
        end if

		cum(1) = 0.
		DO i = 2, SIZE(y)
			cum(i) = cum(i-1) + (x(i)-x(i-1))*(y(i-1) + y(i))/2.
		END DO

	    integral = cum(SIZE(y))

    RETURN
    END SUBROUTINE xvec_yvec_cum_D

 END MODULE trapezoid_quad_m