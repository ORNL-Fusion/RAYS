 module post_processing_m
! Generic wrapper for post processing by individual processors, mostly differentiated by
! The specific geometry. It does the initializations for the specific processor.
!
! This also reads the data from a given RAYS run.  For now this is done either of 2 ways.
! 1) Read an ASCII file containing data from the ray_results_m module as written by routine
!    read_results_instance_LD().  The module data is thread safe and the file is only
!    written at the end of the RAYS RUN. The data read here goes back into the
!    ray_results_m module
! 2) Read an netCDF file containing data from the ray_results_m module as written by
!    type bound procedure %read_results_instance_NC().  This is called on an instance of
!    type run_results.  The data from this instance is then loaded into the ray_results_m
!    module using type bound procedure %to_module()
! Which of the input methods to use is selected in the post_process_rays.in file from
! variable ray_data_input_mode = (LD, or NC), input filenames are constructed from
! namelist variable run_label which should match ray_results_m variable RAYS_run_label
!
! There is a third, legacy, method, which is present for the present but may be eliminated.
! 3) Read an ASCII (or binary) pair of files as specified by unit numbers: ray_list_unit and
!    output_unit,  These files are written incrementally as the rays are traced and
!    therefore are still available if the code crashes.  The data read is put into this
!    module variables: npoints, s_vec, and v_vec.  For now the older processors use this
!    data. ray_data_input_mode = ASCII
!
! It requires an input file 'post_process_rays.in' telling it which specific processor to
! use and which ray data input method to use. It also contains a namelist group with data
! for whichever processor is being used.
!
! It also reads the "rays.in" file (done in subroutine initialize_diagnostics) to get
! some metadata from the RAYS run.

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use constants_m, only : rkind

    implicit none

	character(len=256) :: error_message

! Calculated below from data in input files

    integer :: npoints_max
    integer, allocatable :: npoints(:)
    real(KIND=rkind), allocatable :: s_vec(:,:), v_vec(:,:,:)

!_________________________________________________________________________________________
! Namelist data for /post_process_list/
!_________________________________________________________________________________________

! Switch to select specific post processor
    character(len=80) :: processor = ''

! Selector for ray data input mode -> LD, NC, or ASCII (obsolete but still here)
    character(len=80) :: ray_data_input_mode = ''

    namelist /post_process_list/ processor, ray_data_input_mode

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

 subroutine initialize_post_processing_m(read_input)

!   Reads input file and initializes variables.  Called from main program post_process_RAYS
!   or from another host program.
! N.B. If called from main program post_process_RAYS it is necessary to initialize all of
!      the rays_lib components.  If called from a host program wich has already done the
!      component initializations only the initialize_post_processing_m routine should be
!      called.  Therefore this routine should be called from the host program with
!      initialize_post_processing_m = .false.

    use constants_m, only : initialize_constants_m
    use diagnostics_m, only : initialize_diagnostics, date_v, message_unit, message, &
                & messages_to_stdout, ray_list_unit, output_unit, text_message, verbosity, &
                & run_description, run_label
    use equilibrium_m, only : equilib_model, initialize_equilibrium_m
    use ode_m, only : initialize_ode_solver_m
    use ray_init_m, only : initialize_ray_init_m
    use rf_m, only : initialize_rf_m
    use damping_m, only : initialize_damping_m
    use species_m, only : initialize_species_m
    use ray_results_m, only : initialize_ray_results_m

    use slab_processor_m, only : initialize_slab_processor
    use solovev_processor_m, only : initialize_solovev_processor
    use axisym_toroid_processor_m, only : initialize_axisym_toroid_processor
	use mirror_processor_m, only : initialize_mirror_processor
    use ray_init_m, only : nray_ray_init => nray

    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder

    call text_message('initialize_post_process_rays', 1)

! ****** Initialize the RAYS_lib modules, they read input data from module namelists   ****
    call initialize_diagnostics(read_input) ! This generates new date_v, overwrite with
                                            ! original date_v from RAYS run later
    call message(1)

    call initialize_constants_m
    call message(1)

    call initialize_species_m(read_input)
    call message(1)

    call initialize_rf_m(read_input)
    call message(1)

    call initialize_damping_m(read_input)
    call message(1)

    call initialize_equilibrium_m(read_input)
    call message(1)

    call initialize_ray_init_m(read_input)
    call message(1)

    call initialize_ode_solver_m(read_input)
    call message(1)

    call initialize_ray_results_m(read_input)
    call message(1)

! Read post process namelist file
    if (read_input .eqv. .true.) then
	! Read and write input namelist
  		input_unit = get_unit_number()
        open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
        read(input_unit, post_process_list)
        close(unit=input_unit)
	end if
    if (verbosity > 0) then
		write(message_unit, post_process_list)
		if (messages_to_stdout) write(*, post_process_list)
    end if

    select case (trim(processor))

       case ('slab')
          call initialize_slab_processor(read_input)

       case ('solovev')
          call initialize_solovev_processor(read_input)

       case ('axisym_toroid')
          call initialize_axisym_toroid_processor(read_input)

       case ('multiple_mirror')
          call initialize_mirror_processor(read_input)

       case default
          error_message = 'post_process_rays: unimplemented ray_data_input_mode ='//&
                    & trim(ray_data_input_mode)
          write(*,*) trim(error_message)
          call text_message(trim(error_message))
          stop 1

       end select

!****** Read in all ray data *********************************

    select case (trim(ray_data_input_mode))

       case ('LD')
          call read_results_data_LD

       case ('NC')
           call read_results_data_NC

       case ('ASCII')
          call read_ray_data

       case default
          error_message = 'post_process_rays: unimplemented ray_data_input_mode ='//&
                    & trim(ray_data_input_mode)
          write(*,*) trim(error_message)
          call text_message(trim(error_message))
          stop 1

       end select

    return
 end subroutine initialize_post_processing_m

!*************************************************************************

  subroutine post_process_m

    use diagnostics_m, only : message_unit, message, text_message, verbosity
    use slab_processor_m, only : slab_processor
    use solovev_processor_m, only : solovev_processor
    use axisym_toroid_processor_m, only : axisym_toroid_processor
!	use solovev_magnetics_m, only : inner_bound, outer_bound, rmaj, kappa
    use axisym_toroid_eq_m, only : r_axis, z_axis, &
                          & box_rmin, box_rmax, box_zmin, box_zmax, &
                          & inner_bound, outer_bound, upper_bound, lower_bound
	use mirror_processor_m, only : mirror_processor

    implicit none

    select case (trim(processor))

       case ('slab')
          call text_message('calling slab_processor', 1)
          call slab_processor

       case ('solovev')
          call text_message('calling solovev_processor', 1)
          call solovev_processor

       case ('axisym_toroid')
          call text_message('calling axisym_toroid_processor', 1)
          call axisym_toroid_processor

       case ('multiple_mirror')
          call text_message('calling mirror_processor', 1)
          call mirror_processor

       end select

    return
 end subroutine post_process_m

! Read data from RAYS: two routines: read_results_data_LD and read_results_data_LD


!*************************************************************************

! Read the contents of the ray_results_m module as written by write_results_LD(), which is
! formatted ASCII files written with default list directed output format in finalize_run()

 subroutine read_results_data_LD

    use diagnostics_m, only : message_unit, message, text_message, verbosity, run_label
    use ray_results_m, only : read_results_LD, RAYS_run_label, date_vector

    implicit none

 !  File name for input
    character(len=80) :: in_filename

    in_filename = 'run_results.'//trim(run_label)
	call read_results_LD(in_filename)

	call message(1)
	call text_message('read_results_data: RAYS_run_label = ', RAYS_run_label, 1)
	call message('read_results_data: date_vector', date_vector, 3, 1)
	call message(1)

 end subroutine read_results_data_LD

!*************************************************************************

! Read the contents of the ray_results_m module as written by write_results_NC(), which is
! a netCDF file

 subroutine read_results_data_NC

    use diagnostics_m, only : message_unit, message, text_message, verbosity, run_label
    use ray_results_m, only : run_results, RAYS_run_label, date_vector

    implicit none

 !  File name for input
    character(len=80) :: in_filename

    type(run_results) :: res

    in_filename = 'run_results.'//trim(run_label)//'.nc'
 	call res%read_results_instance_NC(in_filename) ! Get instance of run_results -> res

    call res%to_module ! Load data in res to module

	call message(1)
	call text_message('read_results_data_NC: RAYS_run_label = ', RAYS_run_label, 1)
	call message('read_results_data: date_vector', date_vector, 3, 1)
	call message(1)

 end subroutine read_results_data_NC

! *************************************************************************
! N.B This routine is not used any more.  It is superseded by read_results_data_LD and
! read_results_data_NC.  I haven't deleted it yet just in case.  You know!
!
!  Read number of rays and number of points points on each ray from file: ray_list.<run label>
!  Then read the ray data from file: ray_out.<run_label>.  These are formatted ASCII files
!  written with default list directed output format in ray_tracing() and check_save()
!  respectively.
!

 subroutine read_ray_data

    use diagnostics_m, only : message_unit, message, text_message, verbosity, &
                            & output_unit, ray_list_unit, run_label

    implicit none

    integer :: nray, nv, iray, ipoint
    real(KIND=rkind) :: s
    real(KIND=rkind), allocatable :: v(:)
    real(KIND=rkind), allocatable :: end_residuals(:)
    character(len = 20), allocatable :: ray_stop(:)

!	integer :: input_unit, output_unit, get_unit_number ! External, free unit finder

!*************** Open Input files containing ray data ******************************

!   Open a formatted file containing number of rays and number of steps per ray
    open(unit=ray_list_unit, file='ray_list.'//trim(run_label),action='read', &
                & status='old', form='formatted')

!   Open file containing ray data. File generated in rays as ray.out. Here open for read
    open(unit=output_unit, file='ray_out.'//trim(run_label),action='read', &
                & status='old', form='formatted')

    read(ray_list_unit, *) nray
    call message('nray = ', nray, 1)
    allocate(npoints(nray))
    read(ray_list_unit, *) nv
    allocate(v(nv))

    read(ray_list_unit, *) npoints
    npoints_max = maxval(npoints)
    allocate(end_residuals(nray))
    read(ray_list_unit, *) end_residuals
    allocate(ray_stop(nray))
    read(ray_list_unit, *) ray_stop

    allocate(s_vec(nray, npoints_max))
    allocate(v_vec(nray, npoints_max, nv ))

    ray_loop : do iray = 1, nray
        do ipoint = 1, npoints(iray)

            read (output_unit, *) s, v
            s_vec(iray, ipoint) = s
            v_vec(iray, ipoint, :) = v
        end do
    end do ray_loop

! Write data to an ASCII file for diagnostic
    if (verbosity >= 3) then
       call message(1)
       call text_message('Ray data')
       do iray = 1, nray
            call message('ray number = ', iray)
            do ipoint = 1, npoints(iray)
                write(message_unit, *) s_vec(iray, ipoint), v_vec(iray, ipoint, 1:nv)
            end do
            call text_message('stopped ',ray_stop(iray))
        end do
    end if

  end subroutine read_ray_data

!*************************************************************************

 subroutine finalize_post_processing_m
!   Finish up, do post processing if any

    use diagnostics_m, only : message_unit, message, text_message, ray_list_unit, output_unit

    implicit none

    close(ray_list_unit)
    close(output_unit)

    call message(1)
    call text_message('Post processing RAYS finished', 1)
    call message(1)
    close(message_unit)

! Copy messages file to log.RAYS so it won't get clobbered by post processing
    call system('mv messages log.post_process_RAYS')

    return
 end subroutine finalize_post_processing_m

!*************************************************************************

 end module post_processing_m

