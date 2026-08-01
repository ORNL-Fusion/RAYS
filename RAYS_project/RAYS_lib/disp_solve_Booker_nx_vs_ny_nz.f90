!****************************************************************************

  subroutine solve_Booker_nx_vs_theta_ny_nz(eq, theta, ny, nz, nx)

! Booker quartic is a cold plasma model that only includes electron dynamics.  It allows
! for a B to have a component in the x direction which makes it quartic instead of bi-
! quadratic. nx is complex vector containing all 4 roots
!
! Here (x,y,z) is a local coordinate system:
!   x = The distinguished direction, nx, solved for. e.g. the direction of grad(n)
!   y = direction transverse x, perpendicular to both B and x
!   z = direction = x cross y <--> i.e. B direction if B is perpendicular to grad(n).
!       i.e. B lies in the x,z plane.

! N.B. nx output is complex(KIND=rkind).  Calling program must account for that.

! The roots are classified according to
! whether the sign in front of the discriminant is + or -.  Roots are also classified as
! to which |nx| is larger -> slow mode, or smaller -> fast mode.  This returns both choices:
! n1sq(1) -> plus
! n1sq(2) -> minus
! n1sq(3) -> fast
! n1sq(4) -> slow

! N.B. refractive index convention

    use constants_m, only : rkind, zero, one, two
    use diagnostics_m, only : message
    use equilibrium_m, only : eq_point
    use rpoly_m, only : rpoly

    implicit none

!      Equilibrium data for the spatial point in the plasma
    type(eq_point), intent(in) :: eq

    real(KIND=rkind), intent(in) :: theta, ny, nz
    complex(KIND=rkind), intent(out) :: nx(4)

    real(KIND=rkind) :: X, Y
    real(KIND=rkind) :: C4, C3, C2, C1, C0 ! Booker coefficients in descending powers of nx
    real(KIND=rkind) :: coeffs(5) ! Booker coefficient vector
    integer, parameter :: degree = 4 ! order of equation
    real(KIND=rkind) :: nx_re(4), nx_im(4)
    integer :: i, istat

    X = eq%alpha(0)
    Y = abs(eq%gamma(0))

!  write(*,*) 'solve_Booker: theta = ', theta, '  ny = ', ny, '  nz = ', nz
!  write(*,*) 'solve_Booker: X = ', X, '  Y = ', Y

    C4 = (1 - X/(1 - Y**2))*cos(theta)**2 + (1 - X)*sin(theta)**2

    c3 = (-2*nz*X*Y**2*cos(theta)*sin(theta))/(-1 + Y**2)

    C2 = (ny**2*(-1 + X + Y**2) + (-1 + X)*(-1 + X + Y**2)*cos(theta)**4 + &
          &(-(ny**2*(-1 + X)*(-1 + Y**2)) + (-1 + X)*(-1 + X + Y**2) + &
          &   nz**2*(-2 + 2*X + 2*Y**2 - X*Y**2))*sin(theta)**2 + &
          & (-1 + X)*(-1 + X + Y**2)*sin(theta)**4 + &
          & cos(theta)**2*(1 - 2*nz**2 - 2*X + 2*nz**2*X + X**2 - Y**2 + &
          &   2*nz**2*Y**2 - nz**2*X*Y**2 + ny**2*(-1 + X + Y**2) + &
          &  2*(-1 + X)*(-1 + X + Y**2)*sin(theta)**2))/(-1 + Y**2)

    C1 = (2*nz*(-1 + ny**2 + nz**2)*X*Y**2*cos(theta)*sin(theta))/(1 - Y**2)

    C0 = ((-1 + X)*(1 - 2*X + X**2 - Y**2 + nz**2*(-1 + X + Y**2))*cos(theta)**4 + &
         & (ny**2 + nz**2 + (-1 + X)*sin(theta)**2)* &
         &  (ny**2*(-1 + X + Y**2) + &
         &   (1 - 2*X + X**2 - Y**2 + nz**2*(-1 + X + Y**2))*sin(theta)**2) + &
         & cos(theta)**2*(nz**2*(-1 + X)*(-1 + X + Y**2 - nz**2*(-1 + Y**2)) + &
         &    ny**2*(2 + 2*X**2 - 2*Y**2 + X*(-4 + Y**2) - &
         &       nz**2*(-1 + X)*(-1 + Y**2)) + &
         &    2*(-1 + X)*(1 - 2*X + X**2 - Y**2 + nz**2*(-1 + X + Y**2))*sin(theta)**2 &
             ))/(-1 + Y**2)

    coeffs = (/C4, C3, C2, C1, C0/)
    call rpoly(coeffs, degree, nx_re, nx_im, istat)
    do i = 1,4
        nx(i) = cmplx(nx_re(i), nx_im(i),rkind)
! write(*,*) 'nx(i) = ', nx(i)
    end do

    return
 end subroutine solve_Booker_nx_vs_theta_ny_nz
!****************************************************************************

complex(KIND=rkind) function disp_fun_Booker(eq, theta, nx, ny, nz)
! Calculates the dispersion function of the Booker quartic versus nx, ny, nz where the
! local coordinate system is as described above.
! N.B. nx is complex, as is the returned value disp_fun_Booker.  But ny, nz are real.

    use constants_m, only : rkind, one
    use diagnostics_m, only : message
    use equilibrium_m, only : eq_point

    implicit none

!  Derived type containing equilibrium data for a spatial point in the plasma
   type(eq_point), intent(in) :: eq

   complex(KIND=rkind), intent(in) :: nx
   real(KIND=rkind), intent(in) :: theta, ny, nz

    real(KIND=rkind) :: X, Y
    real(KIND=rkind) :: C4, C3, C2, C1, C0 ! Booker coefficients in descending powers of nx

    X = eq%alpha(0)
    Y = abs(eq%gamma(0))

    C4 = (1 - X/(1 - Y**2))*cos(theta)**2 + (1 - X)*sin(theta)**2

    c3 = (-2*nz*X*Y**2*cos(theta)*sin(theta))/(-1 + Y**2)

    C2 = (ny**2*(-1 + X + Y**2) + (-1 + X)*(-1 + X + Y**2)*cos(theta)**4 + &
          &(-(ny**2*(-1 + X)*(-1 + Y**2)) + (-1 + X)*(-1 + X + Y**2) + &
          &   nz**2*(-2 + 2*X + 2*Y**2 - X*Y**2))*sin(theta)**2 + &
          & (-1 + X)*(-1 + X + Y**2)*sin(theta)**4 + &
          & cos(theta)**2*(1 - 2*nz**2 - 2*X + 2*nz**2*X + X**2 - Y**2 + &
          &   2*nz**2*Y**2 - nz**2*X*Y**2 + ny**2*(-1 + X + Y**2) + &
          &  2*(-1 + X)*(-1 + X + Y**2)*sin(theta)**2))/(-1 + Y**2)

    C1 = (2*nz*(-1 + ny**2 + nz**2)*X*Y**2*cos(theta)*sin(theta))/(1 - Y**2)

    C0 = ((-1 + X)*(1 - 2*X + X**2 - Y**2 + nz**2*(-1 + X + Y**2))*cos(theta)**4 + &
         & (ny**2 + nz**2 + (-1 + X)*sin(theta)**2)* &
         &  (ny**2*(-1 + X + Y**2) + &
         &   (1 - 2*X + X**2 - Y**2 + nz**2*(-1 + X + Y**2))*sin(theta)**2) + &
         & cos(theta)**2*(nz**2*(-1 + X)*(-1 + X + Y**2 - nz**2*(-1 + Y**2)) + &
         &    ny**2*(2 + 2*X**2 - 2*Y**2 + X*(-4 + Y**2) - &
         &       nz**2*(-1 + X)*(-1 + Y**2)) + &
         &    2*(-1 + X)*(1 - 2*X + X**2 - Y**2 + nz**2*(-1 + X + Y**2))*sin(theta)**2 &
             ))/(-1 + Y**2)

    disp_fun_Booker = C0 + nx*C1 + nx**2*C2 + nx**3*C3 + nx**4*C4

       return
 end function disp_fun_Booker

