! Generic wrapper for post processing by individual processors, mostly differentiated by
! The specific geometry. It does the initializations for the specific processor.
! This also reads the data from a given RAYS run.  For now this is done both of 2 ways.
! 1) Read an ASCII (or binary) pair of files as specified by unit numbers: ray_list_unit and
!    output_unit,  These files are written incrementally as the rays are traced and 
!    therefore are still available if the code crashes.  The data read is put into this
!    module variables: npoints, s_vec, and v_vec.  For now the older processors use this
!    data.
! 2) Read an ASCII file containing data from the ray_results_m module.  This file is
!    intended to be thread safe and is only written at the end of the RAYS RUN. The data
!    read goes back into the ray_results_m module
! Eventually I expect the ASCII files will be replaced by NETCDF or the like.  But for now
! I want to avoid external libraries.
!
! Working notes:
! 2/28.2022 (DBB)
! Turn ray_list_unit and output_unit into a subroutine.

 module post_processing_m

    use constants_m, only : rkind

    implicit none
    
! Calculated below from data in input files
    integer :: npoints_max
    integer, allocatable :: npoints(:)
    real(KIND=rkind), allocatable :: s_vec(:,:), v_vec(:,:,:)

! Switch to select specific post processor
    character(len=80) :: processor = ''

! Switches to control reading of ray data
    logical :: read_ray_data_file, read_ray_results_file
    
    namelist /post_process_list/ processor, read_ray_data_file, read_ray_results_file

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
          write(*,*) 'post_process_rays: unimplemented post_processor =', trim(processor)
          call text_message('post_process_rays: unimplemented post_processor', trim(processor))
          stop 1

       end select

!****** Read in all ray data *********************************    

	if (read_ray_data_file) call read_ray_data

	if (read_ray_data_file) call read_results_data
   
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

      case default
          call text_message('post_process_rays: unimplemented post_processor =', trim(processor))
          stop 1

       end select

    return
 end subroutine post_process

!*************************************************************************     

!  Read number of rays and number of points points on each ray from file: ray_list.<run label>
!  Then read the ray data from file: ray_out.<run_label>.
!  Note: The units for ray_list.<run label> (= unit 95) and ray_list.<run label> (= unit 94)
!  are opened in subroutine initialize().  The actual filenames, including .<run_label> are
!  not used here.

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

!  Read number of rays and number of points points on each ray from file: ray_list.<run label>
!  Then read the ray data from file: ray_out.<run_label>.
!  Note: The units for ray_list.<run label> (= unit 95) and ray_list.<run label> (= unit 94)
!  are opened in subroutine initialize().  The actual filenames, including .<run_label> are
!  not used here.

 subroutine read_results_data

    use diagnostics_m, only : message_unit, message, text_message, verbosity, run_label
    use ray_results_m, only : read_results_LD, RAYS_run_label, date_vector

    implicit none
    
 !  File name for input
    character(len=80) :: in_filename

    in_filename = 'run_results.'//trim(run_label)
	call read_results_LD(in_filename)

	call message(1)	
	call text_message('read_results_data: RAYS_run_label = ', RAYS_run_label, 1)
	call message('read_results_data: date_vector = ', date_vector, 3, 1)
	call message(1)
			
 end subroutine read_results_data
 
!*************************************************************************     

 end module post_processing_m

