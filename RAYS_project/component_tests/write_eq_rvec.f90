 program write_eq_rvec
! Stand alone code to write contents of an equilibrium derived type eq_point evaluated
! at position rvec = (x,y,z).
! The equiilibrium is specified by the data contained in a standard 'rays.in' file

   use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message, verbosity, &
                            & output_unit
    use equilibrium_m, only : equilibrium, eq_point, write_eq_point

    implicit none
    logical :: read_input = .true.
	integer :: n_args
	integer :: input_unit, get_unit_number ! External, free unit finder

    real(KIND=rkind) :: x, y, z, rvec(3)

    type(eq_point) :: eq

	write (*,*) 'Enter x, y, z'
    read (*,*) x, y, z
    rvec = (/x, y, z/)
	write (*,*) 'rvec = ', rvec

!   Read input files and initialize variables needed from RAYS_lib
    call initialize_components_needed_for_eq(read_input)

!   Calculate the plasma equilibrium.
    call equilibrium(rvec, eq)

	call write_eq_point(eq)

 contains

!*****************************************************************************************
	  subroutine initialize_components_needed_for_eq(read_input)

! N.B.  This routine only initializes components that are used in to evaulate equilibrium
! Calls the initialize routines for the needed components, and opens output file = "eq_rvec.<run_label>"

		use constants_m, only : initialize_constants_m
		use diagnostics_m, only : initialize_diagnostics, date_v, message_unit, message,&
			& text_message, run_description, run_label, output_unit
		use equilibrium_m, only : equilib_model, initialize_equilibrium_m
		use species_m, only : initialize_species_m

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

		call initialize_equilibrium_m(read_input)
		call message()


!*************** Open output file ******************************

!   Open a file for formatted ray output. File written in check_save()
		output_unit = get_unit_number()
		open(unit=output_unit, file='eq_rvec.'//trim(run_label),action='write', &
					& status='replace', form='formatted')

		return
	 end subroutine initialize_components_needed_for_eq

	 end program write_eq_rvec
