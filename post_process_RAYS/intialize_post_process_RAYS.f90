 subroutine initialize_post_process_RAYS
!   reads input file and initializes variables.  Called from main program RAYS.

    use constants_m, only : input_unit, output_unit, ray_list_unit, initialize_constants_m, eps0
    use diagnostics_m, only : initialize_diagnostics, date_v, message_unit, message,&
        & text_message, run_description, run_label
    use equilibrium_m, only : equilib_model, initialize_equilibrium, equilibrium
    use ode_m, only : initialize_ode_solver
    use ray_init_m, only : initialize_ray_init
    use rf_m, only : omgrf, initialize_rf, frf, k0
    use damping_m, only : initialize_damping
    use species_m, only : initialize_species_m
    use post_processing_m, only : initialize_post_processing_m
    implicit none

    integer :: is 
!************* read input data and set up for messages and diagnostic output  **************


    write(*,*) 'Starting post_processing_RAYS'

!   Read data to set up diagnostic output
    call initialize_diagnostics

    call message()
    call text_message('Starting post_processing_RAYS')
    call text_message(trim(run_description))
    call text_message(trim(run_label))
    call message()
   
    write(*,*) trim(run_description)
    write(*,*) trim(run_label)


!*************** Open Input files containing ray data ******************************

!   Open a formatted file containing number of rays and number of steps per ray
    open(unit=ray_list_unit, file='ray_list.'//trim(run_label),action='read', &
                & status='old', form='formatted')

!   Open file containing ray data. File generated in rays as ray.out. Here open for read
    open(unit=output_unit, file='ray_out.'//trim(run_label),action='read', &
                & status='old', form='formatted')

! Binary writes. Use formatted instead
!   Open a binary file containing ray data
!     open(unit=94, file='rays.bin',                    &
!        & action='read', status='old', form='unformatted')
    
!   Open a binary file describing number of rays and number of steps in each ray
!     open(unit=95, file='ray_list.bin',                  &
!        & action='read', status='old', form='unformatted')
    

! ****** Initialize the modules, they read input data from module namelists   ***********

    call initialize_constants_m
    call message()

    call initialize_species_m
    call message()
    
    call initialize_rf
    call message()
    
    call initialize_damping
    call message()

    call initialize_equilibrium
    call message()

    call initialize_ray_init
    call message()

    call initialize_ode_solver
    call message()
    
    call initialize_post_processing_m
    call message()
    

    return 
 end subroutine initialize_post_process_RAYS
