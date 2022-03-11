 submodule (ode_m) RK4_ode_m
!   contains parameters specific to RK4_ode ODE solver.

    use constants_m, only : rkind

!   Access type definition ode_stop should be available by host association from ode_m.
!   But it confuses make not to use ode_m since it thinks RK4_ode_m doesn't depend on it.
!   Therefore:
    use ode_m, only : ode_stop
    
    implicit none


! namelist /RK4_ode_list/  N.B. RK4 has no parameters so there is no namelist to read

!********************************************************************

contains

!********************************************************************

  module subroutine initialize_RK4_ode

    use constants_m, only : input_unit
    use diagnostics_m, only : message_unit, text_message
    
    implicit none
    
! N.B. RK4 has no parameters so there is no namelist to read    
! Read and write input namelist
!     open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
!     read(input_unit, RK4_ode_list)
!     close(unit=input_unit)
!     write(message_unit, RK4_ode_list)

    call text_message('RK4: RK4 has no initialization')
    return
  end subroutine initialize_RK4_ode

!********************************************************************

  module subroutine ray_init_RK4_ode
  ! Resets RK4_ode solver to appropriate state for start of ray.
  ! N.B. RK4 has no parameters so there is nothing to set. 
    implicit none
    return
  end subroutine ray_init_RK4_ode

!********************************************************************
 
   module subroutine RK4_ode(eqn_ray, nv, v, s, sout, ray_stop)
  ! Simple RK4_ode solver

        use diagnostics_m, only : message_unit, message, verbosity
        use constants_m, only : rkind
        use ode_m, only : ode_stop

        implicit none

      ! Arguments of ODE
        external eqn_ray
        integer, intent(in) :: nv
        real(KIND=rkind), intent(inout) :: v(nv)
        real(KIND=rkind), intent(inout) :: s, sout
        type(ode_stop), intent(out)  :: ray_stop
  
        real (  kind = rkind )::  ds  
        real (  kind = rkind ) :: f1(nv)
        real (  kind = rkind ) :: f2(nv)
        real (  kind = rkind ) :: f3(nv)
        real (  kind = rkind ) :: f4(nv)
  
        ds = sout-s

        call eqn_ray ( s, v(:), f1(:), ray_stop )
        if (ray_stop%stop_ode .eqv. .true.) return
        call eqn_ray ( s + ds/2.0, v(:)+ ds*f1(:)/2.0, f2, ray_stop )
        if (ray_stop%stop_ode .eqv. .true.) return
        call eqn_ray ( s + ds/2.0, v(:)+ ds*f2(:)/2.0, f3, ray_stop )
        if (ray_stop%stop_ode .eqv. .true.) return
        call eqn_ray ( s+ds, v(:)+ ds*f3(:), f4, ray_stop )
        if (ray_stop%stop_ode .eqv. .true.) return

        v = v + ds * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0
        s = sout
        return
  end subroutine RK4_ode

 end submodule RK4_ode_m

