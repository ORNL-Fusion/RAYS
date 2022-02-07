subroutine rk4_RAYS(eqn_ray, nv, v, s, sout, ray_stop)

! Modified for use in RAYS DBB 2/4/2022

!*****************************************************************************80
!
!! rk4 approximates an ODE using a Runge-Kutta fourth order method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    function handle dydt: a function that evaluates the right hand side.
!
!    real (  kind = rkind ) tspan(2): contains the initial and final times.
!
!    real (  kind = rkind ) y0(m): the initial condition.
!
!    integer n: the number of steps to take.
!
!    integer m: the number of variables.
!
!  Output:
!
!    real (  kind = rkind ) t(n+1), y(n+1,m): the times and solution values.
!

  use constants_m, only : rkind
  use ode_m, only : ode_stop

  implicit none

  integer, intent(in) :: nv
  real (  kind = rkind ), intent(inout) :: v(nv)

  real (  kind = rkind ), intent(inout) :: s, sout

  type(ode_stop), intent(out)  :: ray_stop
  
  external eqn_ray 

  real (  kind = rkind )::  ds  
  real (  kind = rkind ) :: f1(nv)
  real (  kind = rkind ) :: f2(nv)
  real (  kind = rkind ) :: f3(nv)
  real (  kind = rkind ) :: f4(nv)
  
  ds = sout-s

    call eqn_ray ( s, v(:), f1(:), ray_stop )
! write(*,*) 'f1 = ', f1
    if (ray_stop%stop_ode .eqv. .true.) return
    call eqn_ray ( s + ds/2.0, v(:)+ ds*f1(:)/2.0, f2, ray_stop )
! write(*,*) 'f2 = ', f2
    call eqn_ray ( s + ds/2.0, v(:)+ ds*f2(:)/2.0, f3, ray_stop )
! write(*,*) 'f3 = ', f3
    if (ray_stop%stop_ode .eqv. .true.) return
    call eqn_ray ( s+ds, v(:)+ ds*f3(:), f4, ray_stop )
! write(*,*) 'f4 = ', f4
    if (ray_stop%stop_ode .eqv. .true.) return

    sout = s+ds
    v = v + ds * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0

! write(*,*) ' v = ', v

  return
end
