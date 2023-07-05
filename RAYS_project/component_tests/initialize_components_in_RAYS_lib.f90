 subroutine initialize_components_in_RAYS_lib(read_input)

! N.B.  This routine only initializes components that are used in the tests, others are
!       commented out. To add tests that need these components uncomment whatever is needed. 

! Calls the initialize routines for all other components, and opens some output files.
! The argument read_input determines whether or not the initialization routines read the
! input data files.  Normally on startup this will be true, but when used as an internal
! component of a host code, the host code may modify initial values by directly using
! the modules and not want data reread from files.

    use constants_m, only : initialize_constants_m
    use diagnostics_m, only : initialize_diagnostics, date_v, message_unit, message,&
        & text_message, run_description, run_label, ray_list_unit, output_unit
    use equilibrium_m, only : equilib_model, initialize_equilibrium_m
    use rf_m, only : initialize_rf_m
    use damping_m, only : initialize_damping_m
    use species_m, only : initialize_species_m
!    use ode_m, only : initialize_ode_solver_m
!    use ray_init_m, only : initialize_ray_init_m
!    use ray_results_m, only : initialize_ray_results_m

    implicit none
    logical, intent(in) :: read_input

!************* read input data and set up for messages and diagnostic output  **************

    write(*,*) 'Initializing components in RAYS_lib'

!   Read data to set up diagnostic output
    call initialize_diagnostics(read_input)

!   Find date and time
    call date_and_time (values=date_v)

    call text_message('Initializing components in RAYS_lib')
    write (message_unit,fmt="(i2,'-',i2,'-',i4,'   ',i2,':',i2,':',i2,'.',i3)") &
     & date_v(2), date_v(3), date_v(1), date_v(5), date_v(6), date_v(7), date_v(8)
    call message()
    call text_message(trim(run_description))
    call text_message(trim(run_label))
    call message()
   
    write(*,*) trim(run_description)
    write(*,*) trim(run_label)

! ****** Initialize the modules, they read input data from module namelists   ***********

    call initialize_constants_m
    call message()

    call initialize_species_m(read_input)
    call message()
 
   call initialize_rf_m(read_input)
    call message()
    
    call initialize_damping_m(read_input)
    call message()
    call initialize_equilibrium_m(read_input)
    call message()

!     call initialize_ray_init_m(read_input)
!     call message()
!     
!     call initialize_ode_solver_m(read_input)
!     call message()
! 
!     call initialize_ray_results_m(read_input)    
!     call message()

!*************** Open output files ******************************

!   Open a formatted file to receive number of rays and number of steps per ray
    open(unit=ray_list_unit, file='ray_list.'//trim(run_label),action='write', &
                & status='replace', form='formatted')


!   Open a file for formatted ray output.
    open(unit=output_unit, file='ray_out.'//trim(run_label),action='write', &
                & status='replace', form='formatted')

! Nota Bene: For now don't write binary output files since nothing uses them.  In long term
!            should replace with netCDF or something.
!   Open a binary file for post-processing.
!     open(unit=94, file='rays.bin',                    &
!        & action='write', status='replace', form='unformatted')
    
!   Open a binary file to receive number of rays and number of steps per ray
!     open(unit=95, file='ray_list.bin',                  &
!        & action='write', status='replace', form='unformatted')
!    write(95) date_v
    
    return 
 end subroutine initialize_components_in_RAYS_lib
