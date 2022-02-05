subroutine rk4_RAYS ( dydt, tspan, y0, n, m, t, y, ray_stop )

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

  real (  kind = rkind ), intent(inout) :: t(n+1)
  real (  kind = rkind ), intent(in) :: tspan(2)
  real (  kind = rkind ) y(n+1,m)
  real (  kind = rkind ) y0(m)

  integer, intent(in) :: m
  integer, intent(in) :: n

  type(ode_stop), intent(out)  :: ray_stop
  
  external dydt ! N.B. here dydt is passed in from RK4_ode_m as eqn_ray 

  real (  kind = rkind )::  dt
  
  real (  kind = rkind ) :: f1(m)
  real (  kind = rkind ) :: f2(m)
  real (  kind = rkind ) :: f3(m)
  real (  kind = rkind ) :: f4(m)
  integer :: i
  
  dt = ( tspan(2) - tspan(1) ) / real ( n,  kind = rkind )

  t(1) = tspan(1)
  y(1,1:m) = y0(1:m)

  do i = 1, n

    call dydt ( t(i),            y(i,1:m),                 f1, ray_stop )
    call dydt ( t(i) + dt / 2.0, y(i,1:m) + dt * f1 / 2.0, f2, ray_stop )
    call dydt ( t(i) + dt / 2.0, y(i,1:m) + dt * f2 / 2.0, f3, ray_stop )
    call dydt ( t(i) + dt,       y(i,1:m) + dt * f3,       f4, ray_stop )

    t(i+1) = t(i) + dt
    y(i+1,1:m) = y(i,1:m) + dt * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0

  end do

  return
end
