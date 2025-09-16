 subroutine initialize(read_input)

! Calls the initialize routines for all other components, and opens some output files.
! The argument read_input determines whether or not the initialization routines read the
! input data files.  Normally on startup this will be true, but when used as an internal
! component of a host code, the host code may modify initial values by directly using
! the modules and not want data reread from files.
!
! N.B. This routine also initializes openMP

! Working notes:

    use constants_m, only : initialize_constants_m
    use diagnostics_m, only : initialize_diagnostics, date_v, message_unit, &
        & ray_list_unit, output_unit, verbosity, write_formatted_ray_files,  &
        & messages_to_stdout, message, text_message, run_description, run_label
    use equilibrium_m, only : equilib_model, initialize_equilibrium_m
    use ode_m, only : initialize_ode_solver_m
    use ray_init_m, only : initialize_ray_init_m
    use rf_m, only : initialize_rf_m
    use damping_m, only : initialize_damping_m
    use species_m, only : initialize_species_m
    use ray_results_m, only : initialize_ray_results_m
!$  use openmp_m, only : initialize_openmp_m

    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder

!************* read input data and set up for messages and diagnostic output  **************

!   Read data to set up diagnostic output
    call initialize_diagnostics(read_input)

    call text_message('initializing RAYS', 0)
    if (verbosity > 0) then
		write (message_unit,fmt="(i2,'-',i2,'-',i4,'   ',i2,':',i2,':',i2,'.',i3)") &
		 & date_v(2), date_v(3), date_v(1), date_v(5), date_v(6), date_v(7), date_v(8)
		if (messages_to_stdout) then
			write (*,fmt="(i2,'-',i2,'-',i4,'   ',i2,':',i2,':',i2,'.',i3)") &
			 & date_v(2), date_v(3), date_v(1), date_v(5), date_v(6), date_v(7), date_v(8)
		end if
    end if
    call text_message('run_description = ',trim(run_description), 1)
    call text_message('run_label = ', trim(run_label), 1)
    call message(1)

! ****** Initialize the modules, they read input data from module namelists   ***********

    call initialize_constants_m
    call message(0)

    call initialize_species_m(read_input)
    call message(0)

    call initialize_rf_m(read_input)
    call message(0)

    call initialize_damping_m(read_input)
    call message(0)

    call initialize_equilibrium_m(read_input)
    call message(0)

    call initialize_ray_init_m(read_input)
    call message(0)

    call initialize_ode_solver_m(read_input)
    call message(0)

    call initialize_ray_results_m(read_input)
    call message(0)

!$  call initialize_openmp_m(read_input)
!$  call message(0)

!*************** Open output files ******************************

	if (write_formatted_ray_files) then
	!   Open a formatted file to receive number of rays and number of steps per ray
	!   File written in ray_tracing()
		ray_list_unit = get_unit_number()
		open(unit=ray_list_unit, file='ray_list.'//trim(run_label),action='write', &
					& status='replace', form='formatted')


	!   Open a file for formatted ray output. File written in check_save()
		output_unit = get_unit_number()
		open(unit=output_unit, file='ray_out.'//trim(run_label),action='write', &
					& status='replace', form='formatted')
	end if

    return
 end subroutine initialize
