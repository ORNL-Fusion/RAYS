 MODULE complete_elliptic_int_m

! Calculates complete elliptic integrals using the argument convention E(m) and K(m).  See
! Abramowitz and Stegan Section 17.3. These routines are extracted a more complete F90
! code by John Burkhardt called "Elliptic_integral" which in turn was adap[ed from F77
! by Bille Carlson, ACM algorithm 577.  The more complete code includes incomplete elliptic
! integrals as well as the "k" argument convention.
!
! N.B. I (DBB) have changed the name of Burkhardt's function "elliptic_fm" to
! elliptic_Km which is more common usage i.e. agrees with Abramoitz and Stegun and
! with Mathematica.
!
! In order to put these functions in a module is was necessary to comment out the
! declarations for functions "rd" and "rf" in subroutines elliptic_em(m) and
! elliptic_fm(m).  The declarations shadowed "rd" and "rf" as being module functions.

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision
    integer, parameter :: rk = rkind

    real(KIND=rkind),parameter :: lolim = 3.D-78
    real(KIND=rkind),parameter :: uplim = 1.D+75

    CONTAINS
 !*********************************************************************************

function elliptic_Em(m)
!*****************************************************************************80
!
!! ELLIPTIC_EM evaluates the complete elliptic integral E(M).
!
!  Discussion:
!
!    The value is computed using Carlson elliptic integrals:
!
!      E(m) = RF ( 0, 1-m, 1 ) - 1/3 m RD ( 0, 1-m, 1 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) M, the argument.
!
!    Output, real ( kind = rk ) ELLIPTIC_EM, the function value.
!
  implicit none

!  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) elliptic_Em
  real ( kind = rk ) errtol
  integer ierr
  real ( kind = rk ) m
!  real ( kind = rk ) rd
!  real ( kind = rk ) rf
  real ( kind = rk ) value
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  x = 0.0D+00
  y = 1.0D+00 - m
  z = 1.0D+00
  errtol = 1.0D-03

  value = rf ( x, y, z, errtol, ierr ) &
    - m * rd ( x, y, z, errtol, ierr ) / 3.0D+00

  elliptic_Em = value

  return
end function elliptic_Em

function elliptic_Km( m )
!*****************************************************************************80
!
!! ELLIPTIC_FM evaluates the complete elliptic integral F(M).
!
!  Discussion:
!
!    The value is computed using Carlson elliptic integrals:
!
!      F(m) = RF ( 0, 1-m, 1 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 May 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) M, the argument.
!
!    Output, real ( kind = rk ) ELLIPTIC_FM, the function value.
!
! N.B.  "elliptic_fm" renamed to "elliptic_Km" (DBB)

  implicit none

!  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) elliptic_Km
  real ( kind = rk ) errtol
  integer ierr
  real ( kind = rk ) m
!  real ( kind = rk ) rf
  real ( kind = rk ) value
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  x = 0.0D+00
  y = 1.0D+00 - m
  z = 1.0D+00
  errtol = 1.0D-03

  value = rf( x, y, z, errtol, ierr )

  elliptic_Km = value

  return
end function elliptic_Km

function rf( x, y, z, errtol, ierr )
!*****************************************************************************
!
!! RF computes an incomplete elliptic integral of the first kind, RF(X,Y,Z).
!
!  Discussion:
!
!    This function computes the incomplete elliptic integral of the first kind.
!
!    RF(X,Y,Z) = Integral ( 0 <= T < oo )
!
!                    -1/2     -1/2     -1/2
!          (1/2)(T+X)    (T+Y)    (T+Z)    DT,
!
!    where X, Y, and Z are nonnegative and at most one of them is zero.
!
!    If X or Y or Z is zero, the integral is complete.
!
!    The duplication theorem is iterated until the variables are
!    nearly equal, and the function is then expanded in Taylor
!    series to fifth order.
!
!    Check by addition theorem:
!
!      RF(X,X+Z,X+W) + RF(Y,Y+Z,Y+W) = RF(0,Z,W),
!      where X, Y, Z, W are positive and X * Y = Z * W.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Bille Carlson,
!    Computing Elliptic Integrals by Duplication,
!    Numerische Mathematik,
!    Volume 33, 1979, pages 1-16.
!
!    Bille Carlson, Elaine Notis,
!    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, pages 398-403, September 1981.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, Z, the arguments in the integral.
!
!    Input, real ( kind = rk ) ERRTOL, the error tolerance.
!    Relative error due to truncation is less than
!      ERRTOL ^ 6 / (4 * (1 - ERRTOL)).
!    Sample choices:
!      ERRTOL   Relative truncation error less than
!      1.D-3    3.D-19
!      3.D-3    2.D-16
!      1.D-2    3.D-13
!      3.D-2    2.D-10
!      1.D-1    3.D-7
!
!    Output, integer IERR, the error flag.
!    0, no error occurred.
!    1, abnormal termination.
!
  implicit none

!  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) c3
  real ( kind = rk ) e2
  real ( kind = rk ) e3
  real ( kind = rk ) epslon
  real ( kind = rk ) errtol
  integer ierr
  real ( kind = rk ) lamda
!  real ( kind = rk ) lolim
  real ( kind = rk ) mu
  real ( kind = rk ) rf
  real ( kind = rk ) s
!  real ( kind = rk ) uplim
  real ( kind = rk ) x
  real ( kind = rk ) xn
  real ( kind = rk ) xndev
  real ( kind = rk ) xnroot
  real ( kind = rk ) y
  real ( kind = rk ) yn
  real ( kind = rk ) yndev
  real ( kind = rk ) ynroot
  real ( kind = rk ) z
  real ( kind = rk ) zn
  real ( kind = rk ) zndev
  real ( kind = rk ) znroot
!
!  LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!  LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
!  UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.
!
!   save lolim
!   save uplim
!
!   data lolim /3.D-78/
!   data uplim /1.D+75/

  if ( &
    x < 0.0D+00 .or. &
    y < 0.0D+00 .or. &
    z < 0.0D+00 .or. &
    x + y < lolim .or. &
    x + z < lolim .or. &
    y + z < lolim .or. &
    uplim <= x .or. &
    uplim <= y .or. &
    uplim <= z ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RF - Error!'
    write ( *, '(a)' ) '  Invalid input arguments.'
    write ( *, '(a,d23.16)' ) '  X = ', x
    write ( *, '(a,d23.16)' ) '  Y = ', y
    write ( *, '(a,d23.16)' ) '  Z = ', z
    write ( *, '(a)' ) ''
    ierr = 1
    rf = 0.0D+00
    return
  end if

  ierr = 0
  xn = x
  yn = y
  zn = z

  do

    mu = ( xn + yn + zn ) / 3.0d0
    xndev = 2.0d0 - ( mu + xn ) / mu
    yndev = 2.0d0 - ( mu + yn ) / mu
    zndev = 2.0d0 - ( mu + zn ) / mu
    epslon = max ( abs ( xndev ), abs ( yndev ), abs ( zndev ) )

    if ( epslon < errtol ) then
      c1 = 1.0d0 / 24.0d0
      c2 = 3.0d0 / 44.0d0
      c3 = 1.0d0 / 14.0d0
      e2 = xndev * yndev - zndev * zndev
      e3 = xndev * yndev * zndev
      s = 1.0d0 + ( c1 * e2 - 0.1d0 - c2 * e3 ) * e2 + c3 * e3
      rf = s / sqrt ( mu )
      return
    end if

    xnroot = sqrt ( xn )
    ynroot = sqrt ( yn )
    znroot = sqrt ( zn )
    lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot
    xn = ( xn + lamda ) * 0.25d0
    yn = ( yn + lamda ) * 0.25d0
    zn = ( zn + lamda ) * 0.25d0

  end do
  return
end function rf

function rd( x, y, z, errtol, ierr )
!*****************************************************************************80
!
!! RD computes an incomplete elliptic integral of the second kind, RD(X,Y,Z).
!
!  Discussion:
!
!    This function computes an incomplete elliptic integral of the second kind.
!
!    RD(X,Y,Z) = Integral ( 0 <= T < oo )
!
!                    -1/2     -1/2     -3/2
!          (3/2)(T+X)    (T+Y)    (T+Z)    DT,
!
!    where X and Y are nonnegative, X + Y is positive, and Z is positive.
!
!    If X or Y is zero, the integral is complete.
!
!    The duplication theorem is iterated until the variables are
!    nearly equal, and the function is then expanded in Taylor
!    series to fifth order.
!
!    Check:
!
!      RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y) = 3 / sqrt ( X * Y * Z ),
!      where X, Y, and Z are positive.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Bille Carlson,
!    Computing Elliptic Integrals by Duplication,
!    Numerische Mathematik,
!    Volume 33, 1979, pages 1-16.
!
!    Bille Carlson, Elaine Notis,
!    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, pages 398-403, September 1981.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, Z, the arguments in the integral.
!
!    Input, real ( kind = rk ) ERRTOL, the error tolerance.
!    The relative error due to truncation is less than
!      3 * ERRTOL ^ 6 / (1-ERRTOL) ^ 3/2.
!    Sample choices:
!      ERRTOL   Relative truncation error less than
!      1.D-3    4.D-18
!      3.D-3    3.D-15
!      1.D-2    4.D-12
!      3.D-2    3.D-9
!      1.D-1    4.D-6
!
!    Output, integer IERR, the error flag.
!    0, no error occurred.
!    1, abnormal termination.
!
  implicit none

!  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) c3
  real ( kind = rk ) c4
  real ( kind = rk ) ea
  real ( kind = rk ) eb
  real ( kind = rk ) ec
  real ( kind = rk ) ed
  real ( kind = rk ) ef
  real ( kind = rk ) epslon
  real ( kind = rk ) errtol
  integer ierr
  real ( kind = rk ) lamda
!   real ( kind = rk ) lolim
  real ( kind = rk ) mu
  real ( kind = rk ) power4
  real ( kind = rk ) rd
  real ( kind = rk ) sigma
  real ( kind = rk ) s1
  real ( kind = rk ) s2
!   real ( kind = rk ) uplim
  real ( kind = rk ) x
  real ( kind = rk ) xn
  real ( kind = rk ) xndev
  real ( kind = rk ) xnroot
  real ( kind = rk ) y
  real ( kind = rk ) yn
  real ( kind = rk ) yndev
  real ( kind = rk ) ynroot
  real ( kind = rk ) z
  real ( kind = rk ) zn
  real ( kind = rk ) zndev
  real ( kind = rk ) znroot
!
!  LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!  LOLIM IS NOT LESS THAN 2 / (MACHINE MAXIMUM) ^ (2/3).
!  UPLIM IS NOT GREATER THAN (0.1 * ERRTOL / MACHINE
!  MINIMUM) ^ (2/3), WHERE ERRTOL IS DESCRIBED BELOW.
!  IN THE FOLLOWING TABLE IT IS ASSUMED THAT ERRTOL WILL
!  NEVER BE CHOSEN SMALLER THAN 1.D-5.
!
!   save lolim
!   save uplim
!
!   data lolim /6.D-51/
!   data uplim /1.D+48/

  if ( &
    x < 0.0D+00 .or. &
    y < 0.0D+00 .or. &
    x + y < lolim .or. &
    z < lolim .or. &
    uplim < x .or. &
    uplim < y .or. &
    uplim < z ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RD - Error!'
    write ( *, '(a)' ) '  Invalid input arguments.'
    write ( *, '(a,d23.16)' ) '  X = ', x
    write ( *, '(a,d23.16)' ) '  Y = ', y
    write ( *, '(a,d23.16)' ) '  Z = ', z
    write ( *, '(a)' ) ''
    ierr = 1
    rd = 0.0D+00
    return
  end if

  ierr = 0
  xn = x
  yn = y
  zn = z
  sigma = 0.0d0
  power4 = 1.0d0

  do

    mu = ( xn + yn + 3.0d0 * zn ) * 0.2d0
    xndev = ( mu - xn ) / mu
    yndev = ( mu - yn ) / mu
    zndev = ( mu - zn ) / mu
    epslon = max ( abs ( xndev ), abs ( yndev ), abs ( zndev ) )

    if ( epslon < errtol ) then
      c1 = 3.0d0 / 14.0d0
      c2 = 1.0d0 / 6.0d0
      c3 = 9.0d0 / 22.0d0
      c4 = 3.0d0 / 26.0d0
      ea = xndev * yndev
      eb = zndev * zndev
      ec = ea - eb
      ed = ea - 6.0d0 * eb
      ef = ed + ec + ec
      s1 = ed * ( - c1 + 0.25d0 * c3 * ed - 1.5d0 * c4 * zndev * ef )
      s2 = zndev * ( c2 * ef + zndev * ( - c3 * ec + zndev * c4 * ea ) )
      rd = 3.0d0 * sigma + power4 * ( 1.0d0 + s1 + s2 ) / ( mu * sqrt ( mu ) )

      return
    end if

    xnroot = sqrt ( xn )
    ynroot = sqrt ( yn )
    znroot = sqrt ( zn )
    lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot
    sigma = sigma + power4 / ( znroot * ( zn + lamda ) )
    power4 = power4 * 0.25d0
    xn = ( xn + lamda ) * 0.25d0
    yn = ( yn + lamda ) * 0.25d0
    zn = ( zn + lamda ) * 0.25d0

  end do
  return

end function rd


 END MODULE complete_elliptic_int_m