 subroutine initialize_post_process_RAYS(read_input)
!   Reads input file and initializes variables.  Called from main program post_process_RAYS
!   or from another host program.
! N.B. If called from main program post_process_RAYS it is necessary to initialize all of
!      the rays_lib components.  If called from a host program wich has already done the
!      component initializations only the initialize_post_processing_m routine should be
!      called.  Therefore this routine should be called from the host program with 
!      initialize_post_processing_m = .false.

    use constants_m, only : initialize_constants_m
    use diagnostics_m, only : initialize_diagnostics, date_v, message_unit, message, &
                & ray_list_unit, output_unit, text_message, run_description, run_label
    use equilibrium_m, only : equilib_model, initialize_equilibrium_m
    use ode_m, only : initialize_ode_solver_m
    use ray_init_m, only : initialize_ray_init_m
    use rf_m, only : initialize_rf_m
    use damping_m, only : initialize_damping_m
    use species_m, only : initialize_species_m
    use ray_results_m, only : initialize_ray_results_m
    use post_processing_m, only : initialize_post_processing_m

    implicit none
    logical, intent(in) :: read_input

    integer :: is 
!************* read input data and set up for messages and diagnostic output  **************


    write(*,*) 'Starting post_processing_RAYS'

!   Read data to set up diagnostic output
    call initialize_diagnostics(read_input)

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

    call initialize_species_m(read_input)
    call message()
    
    call initialize_rf_m(read_input)
    call message()
    
    call initialize_damping_m(read_input)
    call message()

    call initialize_equilibrium_m(read_input)
    call message()

    call initialize_ray_init_m(read_input)
    call message()
    
    call initialize_ode_solver_m(read_input)
    call message()

    call initialize_ray_results_m(read_input)    
    call message()
    
    call initialize_post_processing_m(.true.)
    call message()
    

    return 
 end subroutine initialize_post_process_RAYS
