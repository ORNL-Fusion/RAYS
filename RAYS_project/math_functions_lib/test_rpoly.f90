PROGRAM test_rpoly
!USE Solve_Real_Poly
USE rpoly_m, only : rpoly
IMPLICIT NONE

 integer, parameter :: dp = selected_real_kind(15,307) ! kind parameter for reals


REAL (dp)  :: p(50), zr(50), zi(50)
INTEGER    :: degree, i, istat
LOGICAL    :: fail = .false.

WRITE(*, 5000)
degree = 10
p(1) = 1._dp
p(2) = -55._dp
p(3) = 1320._dp
p(4) = -18150._dp
p(5) = 157773._dp
p(6) = -902055._dp
p(7) = 3416930._dp
p(8) = -8409500._dp
p(9) = 12753576._dp
p(10) = -10628640._dp
p(11) = 3628800._dp

!CALL rpoly(p, degree, zr, zi, fail)
CALL rpoly(p, degree, zr, zi, istat)

IF (fail) THEN
  WRITE(*, *) ' ** Failure by RPOLY **'
ELSE
  WRITE(*, '(a/ (2g23.15))') ' Real part           Imaginary part',  &
                             (zr(i), zi(i), i=1,degree)
END IF

! This test provided by Larry Wigton

WRITE(*, *)
WRITE(*, *) "Now try case where 1 is an obvious root"
degree = 5
p(1) = 8.D0
p(2) = -8.D0
p(3) = 16.D0
p(4) = -16.D0
p(5) = 8.D0
p(6) = -8.D0

!CALL rpoly(p, degree, zr, zi, fail)
CALL rpoly(p, degree, zr, zi, istat)

IF (fail) THEN
  WRITE(*, *) ' ** Failure by RPOLY **'
ELSE
  WRITE(*, *) ' Real part           Imaginary part'
  WRITE(*, '(2g23.15)') (zr(i), zi(i), i=1,degree)
END IF


WRITE(*, *)
WRITE(*, *) "4th order"
degree = 4
p(1) = 1.D0
p(2) = 6.D0
p(3) = -41.D0
p(4) = -6.D0
p(5) = 40.D0

!CALL rpoly(p, degree, zr, zi, fail)
CALL rpoly(p, degree, zr, zi, istat)

IF (fail) THEN
  WRITE(*, *) ' ** Failure by RPOLY **'
ELSE
  WRITE(*, *) ' Real part           Imaginary part'
  WRITE(*, '(2g23.15)') (zr(i), zi(i), i=1,degree)
END IF

STOP

5000 FORMAT (' EXAMPLE 1. POLYNOMIAL WITH ZEROS 1,2,...,10.')

END PROGRAM test_rpoly