PROGRAM test_trapezoid
    USE trapezoid_quad_m, only : trapezoid_quad, trapezoid_quad_cumulative

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
!    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

    INTEGER, PARAMETER :: npoints = 101
    REAL :: h, x, xvec(npoints), y(npoints), integral
    REAL(kind = rkind) :: hD, xD, xDvec(npoints),yD(npoints), integralD
    COMPLEX :: yC(npoints), integralC
    COMPLEX(kind = rkind) :: yDC(npoints), integralDC

    REAL :: cumulative(npoints)
    REAL(kind = rkind) :: cumulativeD(npoints)

    real :: pi = 4.0*ATAN(1.0)
    REAL(kind = rkind) :: piD = 4.D0*ATAN(1.D0)
    integer :: i

    h = 1.0/(npoints-1)
    hD = 1.0_rkind/(npoints-1)

    DO i = 1, npoints
        x = h*(i-1)
        xvec(i) = x
        y(i) = x**3
        yC(i) = complex(cos(2.*pi*x), sin(2.*pi*x))

        xD = hD*(i-1)
        xDvec(i) = xD
        yD(i) = xD**3
        yDC(i) = complex(cos(2.*piD*xD), sin(2.*piD*xD))
    END DO

! fixed increment h, hD
    write(*,*) ' '
    write(*,*) 'fixed increment: h = ', h
    write(*,*) 'y = ', y

    write(*,*) ' '
    CALL trapezoid_quad(h, y, integral)
    WRITE (*, *) ' single precision integral = ', integral

    CALL trapezoid_quad(hD, yD, integralD)
    WRITE (*, *) ' double precision integral = ', integralD

    CALL trapezoid_quad(h, yC, integralC)
    WRITE (*, *) ' cmplx single precision integral = ', integralC

    CALL trapezoid_quad(hD, yDC, integralDC)
    WRITE (*, *) ' cmplx double precision integral = ', integralDC

    CALL trapezoid_quad_cumulative(h, y, integral, cumulative)
    WRITE (*, *) ' single precision integral = ', integral
    WRITE (*, *) ' cumulative integral = ', cumulative

    CALL trapezoid_quad_cumulative(hD, yD, integralD, cumulativeD)
    WRITE (*, *) ' double precision integral = ', integralD
    WRITE (*, *) ' cumulative integral = ', cumulativeD


! x vector
    write(*,*) ' '
    write(*,*) 'x vector'
    write(*,*) 'xvec = ', xvec

    write(*,*) ' '
    CALL trapezoid_quad(xvec, y, integral)
    WRITE (*, *) ' xvec single precision integral = ', integral

    CALL trapezoid_quad(xDvec, yD, integralD)
    WRITE (*, *) ' xvec double precision integral = ', integralD

    CALL trapezoid_quad(xvec, yC, integralC)
    WRITE (*, *) ' xvec cmplx single precision integral = ', integralC

    CALL trapezoid_quad(xDvec, yDC, integralDC)
    WRITE (*, *) ' xvec cmplx double precision integral = ', integralDC

    CALL trapezoid_quad_cumulative(xvec, y, integral, cumulative)
    WRITE (*, *) ' xvec single precision integral = ', integral
    WRITE (*, *) ' cumulative integral = ', cumulative

    CALL trapezoid_quad_cumulative(xDvec, yD, integralD, cumulativeD)
    WRITE (*, *) ' xvec double precision integral = ', integralD
    WRITE (*, *) ' cumulative integral = ', cumulativeD


!     write(*,*) ' '
!     WRITE (*, *) 'integral = ', integral
!     WRITE (*, *) 'cumulative integral = ', y_cumulative

END PROGRAM test_trapezoid
