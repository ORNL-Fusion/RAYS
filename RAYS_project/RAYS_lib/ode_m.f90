 module ode_m
!   contains parameters for initial value ODE solver.

    use constants_m, only : rkind
    
    implicit none

!   Switch to select ODE solvers for ray tracing.  Presently supported solvers are:
!   ode_solver= SG_ode: subroutine ODE developed by L. F. Shampine and M. K. Gordon.
!   ode_solver = RK4_ode: Simple Runge-Kutta 4th order integrator.
    character (len = 15) :: ode_solver_name  

!   Name of routine used to calculate RHS derivatives for ray equations. Presently
!   supported routines are:
!   ray_deriv_name = cold
!   ray_deriv_name = numerical
    character(len=60) :: ray_deriv_name
    
!   nv: No. of ODEs to be solved.
    integer :: nv 
    
! Maximum length of ray
     real(KIND=rkind) :: s_max

!   ODE step size.
    real(KIND=rkind) :: ds

!   Maximum no. of steps allowed.
    integer :: nstep_max
   
! Derived type containing information on terminating ode integration
    type ode_stop    
        logical :: stop_ode = .false.
        character(len=60) :: ode_stop_flag = ''
    end type ode_stop

! Interfaces for Shampine and Gordon ODE 
    interface initialize_SG_ode
        module subroutine initialize_SG_ode(read_input)
          logical, intent(in) :: read_input
        end subroutine initialize_SG_ode
    end interface initialize_SG_ode

    interface ray_init_SG_ode
         module subroutine ray_init_SG_ode
         end subroutine ray_init_SG_ode
    end interface ray_init_SG_ode

    interface SG_ode
      module subroutine SG_ode(eqn_ray, nv, v, s, sout, ray_stop)
        use diagnostics_m, only : message_unit, message, verbosity
        use constants_m, only : rkind
      ! Arguments of ODE
        external eqn_ray
        integer, intent(in) :: nv
        real(KIND=rkind), intent(inout) :: v(nv)
        real(KIND=rkind), intent(inout) :: s, sout
        type(ode_stop), intent(out)  :: ray_stop
      end subroutine SG_ode
    end interface SG_ode

! Interfaces for RK4 ODE 
    interface initialize_RK4_ode
        module subroutine initialize_RK4_ode
        end subroutine initialize_RK4_ode
    end interface initialize_RK4_ode

    interface ray_init_RK4_ode
         module subroutine ray_init_RK4_ode
         end subroutine ray_init_RK4_ode
    end interface ray_init_RK4_ode

    interface RK4_ode
      module subroutine RK4_ode(eqn_ray, nv, v, s, sout, ray_stop)
        use diagnostics_m, only : message_unit, message, verbosity
        use constants_m, only : rkind
      ! Arguments of ODE
        external eqn_ray
        integer, intent(in) :: nv
        real(KIND=rkind), intent(inout) :: v(nv)
        real(KIND=rkind), intent(inout) :: s, sout
        type(ode_stop), intent(out)  :: ray_stop
      end subroutine RK4_ode
    end interface RK4_ode

! Namelist   
    namelist /ode_list/ ode_solver_name, ray_deriv_name, nstep_max, s_max, ds

!********************************************************************

 contains

!********************************************************************

  subroutine initialize_ode_solver(read_input)

    use constants_m, only : input_unit
    use damping_m, only : damping_model, multi_spec_damping
    use diagnostics_m, only : message_unit, message, text_message, integrate_eq_gradients
    use species_m, only : nspec
    
    implicit none
    logical, intent(in) :: read_input

    if (read_input .eqv. .true.) then    
    ! Read and write input namelist
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, ode_list)
        close(unit=input_unit)
        write(message_unit, ode_list)
    end if

!   Select ode solver.
    solver: select case (trim(ode_solver_name))

       case ('SG_ODE')
          call initialize_SG_ode(read_input)

       case ('RK4_ODE')
          call initialize_RK4_ode

       case default
          write(0,*) 'read_ode_namelists, invalid ode solver = ', trim(ode_solver_name)
          call text_message ('read_ode_namelists: invalid ode solver', trim(ode_solver_name),0)
          stop 1

    end select solver


!   Calculate number of of ODEs to be solved, nv.

    nv = 7  ! Default, integrate x(s) and k(s) and s
    
    if (trim(damping_model) == 'no_damp') then
        continue
    else
        nv = nv +1  ! Also integrate total damping.
    endif
    
    ! Also integrate damping for each species (electrons + nspec ion species)   
    if (multi_spec_damping .eqv. .true.) nv = nv + 1 + nspec

    ! Also integrate grad(B) tensor, grad(ne/ne0) and grad(te/te0) along ray for 
    ! diagnostic purposes
    if (integrate_eq_gradients .eqv. .true. ) nv = nv + 5
    
    call message ('initialize_ode_solver: ODE vector length nv', nv)

    return
  end subroutine initialize_ode_solver
    
!********************************************************************

  subroutine ray_init_ode_solver
! Generic wrapper for various ODE solvers. Resets state of ode solver to whatever is
! required at the beginning of the ray.

!Presently supported solvers are:
! ode_solver = SG_ode: subroutine ODE developed by by L. F. Shampine and M. K. Gordon.
! ode_solver = RK4_ode: Simple Runge-Kutta 4th order integrator.  The initialization
!              routine ray_init_RK4_ode does nothing since RK4 has no state to be reset.

    use diagnostics_m, only : message, text_message

      
    implicit none

!   Select ode solver.
    solver: select case (trim(ode_solver_name))

       case ('SG_ODE')
          call ray_init_SG_ode

       case ('RK4_ODE')
          call ray_init_RK4_ode

       case default  ! By this point there has to be a solver
           stop 2

    end select solver
 
  return
  end subroutine ray_init_ode_solver
    
!********************************************************************

  subroutine ode_solver(eqn_ray, nv, v, s, sout, ray_stop)
! Generic wrapper for various ODE solvers. Takes one step integrating from s to sout.

!Presently supported solvers are:
! ode_solver = SG_ode: subroutine ODE developed by by L. F. Shampine and M. K. Gordon.

    use diagnostics_m, only : message, text_message
     
    implicit none

    integer, intent(in) :: nv
    real(KIND=rkind), intent(inout) :: v(nv)
    real(KIND=rkind), intent(inout) :: s
    real(KIND=rkind), intent(inout) :: sout
    type(ode_stop), intent(out)  :: ray_stop
    
    external eqn_ray
    
!   Select ode solver.
    solver: select case (trim(ode_solver_name))

       case ('SG_ODE')
          call SG_ode(eqn_ray, nv, v, s, sout, ray_stop)

       case ('RK4_ODE')
          call RK4_ode(eqn_ray, nv, v, s, sout, ray_stop)

       case default
          write(0,*) 'ode_solver, invalid ode solver = ', trim(ode_solver_name)
          call text_message ('ode_solver: invalid ode solver', trim(ode_solver_name),0)
          stop 2

    end select solver
 
  return
  end subroutine ode_solver
    
 
 end module ode_m

