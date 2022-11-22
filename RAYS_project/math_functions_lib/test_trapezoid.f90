PROGRAM test_trapezoid
    USE trapezoid_quad_m, only : trapezoid_quad

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
!    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

    INTEGER, PARAMETER :: npoints = 101
    REAL :: h, x, xvec(npoints), y(npoints), integral
    REAL(kind = rkind) :: hD, xD, xDvec(npoints),yD(npoints), integralD
    COMPLEX :: yC(npoints), integralC
    COMPLEX(kind = rkind) :: yDC(npoints), integralDC

!     REAL :: y_cumulative(npoints)
!     REAL(kind = rkind) :: y_cumulativeD(npoints)
!     COMPLEX :: yC_cumulative(npoints)
!     COMPLEX(kind = rkind) :: y_cumulativeDC(npoints)
    
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


! x vector   
    write(*,*) ' '
    write(*,*) 'x vector'
    write(*,*) 'xvec = ', xvec

    write(*,*) ' '
    CALL trapezoid_quad(xvec, y, integral)    
    WRITE (*, *) ' single precision integral = ', integral
    
    CALL trapezoid_quad(xDvec, yD, integralD)    
    WRITE (*, *) ' double precision integral = ', integralD
    
    CALL trapezoid_quad(xvec, yC, integralC)    
    WRITE (*, *) ' cmplx single precision integral = ', integralC
    
    CALL trapezoid_quad(xDvec, yDC, integralDC)    
    WRITE (*, *) ' cmplx double precision integral = ', integralDC



!     write(*,*) ' '
!     WRITE (*, *) 'integral = ', integral
!     WRITE (*, *) 'cumulative integral = ', y_cumulative
    
END PROGRAM test_trapezoid
    