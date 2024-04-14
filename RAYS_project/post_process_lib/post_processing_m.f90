! Generic wrapper for post processing by individual processors, mostly differentiated by
! The specific geometry. It does the initializations for the specific processor.
!
! This also reads the data from a given RAYS run.  For now this is done both of 2 ways.
! 1) Read an ASCII (or binary) pair of files as specified by unit numbers: ray_list_unit and
!    output_unit,  These files are written incrementally as the rays are traced and
!    therefore are still available if the code crashes.  The data read is put into this
!    module variables: npoints, s_vec, and v_vec.  For now the older processors use this
!    data.
! 2) Read an ASCII file containing data from the ray_results_m module as written by routine
!    read_results_instance_LD().  The module data is thread safe and the file is only
!    written at the end of the RAYS RUN. The data read here goes back into the
!    ray_results_m module
! 3) Read an netCDF file containing data from the ray_results_m module as written by
!    type bound procedure %read_results_instance_NC().  This is called on an instance of
!    type run_results.  The data from this instance is then loaded into the ray_results_m
!    module using type bound procedure %to_module()
! Which of the input methods to use is selected in the post_process_rays.in file from
! variable ray_data_input_mode = (ASCII, LD, or NC), input filenames are constructed from
! namelist variable run_label which should match ray_results_m variable RAYS_run_label

 module post_processing_m

    use constants_m, only : rkind

    implicit none

	character(len=256) :: error_message

! Calculated below from data in input files
    integer :: npoints_max
    integer, allocatable :: npoints(:)
    real(KIND=rkind), allocatable :: s_vec(:,:), v_vec(:,:,:)

! Switch to select specific post processor
    character(len=80) :: processor = ''

! Selector for ray data input mode
    character(len=80) :: ray_data_input_mode = ''

    namelist /post_process_list/ processor, ray_data_input_mode

 contains

 subroutine initialize_post_processing_m(read_input)

    use diagnostics_m, only : message_unit, message, text_message, verbosity, &
                            & output_unit, ray_list_unit
    use slab_processor_m, only : initialize_slab_processor
    use solovev_processor_m, only : initialize_solovev_processor
    use axisym_toroid_processor_m, only : initialize_axisym_toroid_processor
    use ray_init_m, only : nray_ray_init => nray

    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder

    call text_message('initialize_post_process_rays', 3)

    if (read_input .eqv. .true.) then
	! Read and write input namelist
  		input_unit = get_unit_number()
        open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
        read(input_unit, post_process_list)
        close(unit=input_unit)
        write(message_unit, post_process_list)
	end if

    select case (trim(processor))

       case ('slab')
          call initialize_slab_processor(read_input)

       case ('solovev')
          call initialize_solovev_processor(read_input)

       case ('axisym_toroid')
          call initialize_axisym_toroid_processor(read_input)

       case default
          error_message = 'post_process_rays: unimplemented ray_data_input_mode ='//&
                    & trim(ray_data_input_mode)
          write(*,*) trim(error_message)
          call text_message(trim(error_message))
          stop 1

       end select

!****** Read in all ray data *********************************

    select case (trim(ray_data_input_mode))

       case ('ASCII')
          call read_ray_data

       case ('LD')
          call read_results_data_LD

       case ('NC')
           call read_results_data_NC

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

  subroutine post_process

    use diagnostics_m, only : message_unit, message, text_message, verbosity
    use slab_processor_m, only : slab_processor
    use solovev_processor_m, only : solovev_processor
    use axisym_toroid_processor_m, only : axisym_toroid_processor
!	use solovev_magnetics_m, only : inner_bound, outer_bound, rmaj, kappa
    use axisym_toroid_eq_m, only : r_axis, z_axis, &
                          & box_rmin, box_rmax, box_zmin, box_zmax, &
                          & inner_bound, outer_bound, upper_bound, lower_bound

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

       end select

    return
 end subroutine post_process

!*************************************************************************

!  Read number of rays and number of points points on each ray from file: ray_list.<run label>
!  Then read the ray data from file: ray_out.<run_label>.  These are formatted ASCII files
!  written with default list directed output format in ray_tracing() and check_save()
!  respectively.
!
!  Note: The units for ray_list.<run label> and ray_out.<run label>  are set, and files are
!  opened, in subroutine initialize().  The actual filenames, including .<run_label> are
!  not needed in this subroutine.

 subroutine read_ray_data

    use diagnostics_m, only : message_unit, message, text_message, verbosity, &
                            & output_unit, ray_list_unit

    implicit none

    integer :: nray, nv, iray, ipoint
    real(KIND=rkind) :: s
    real(KIND=rkind), allocatable :: v(:)
    real(KIND=rkind), allocatable :: end_residuals(:)
    character(len = 20), allocatable :: ray_stop(:)

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

!*************************************************************************

 end module post_processing_m

