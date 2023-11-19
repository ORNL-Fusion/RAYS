 module ray_results_m

! A module that contains results from ray tracing for storage in memory or for writing once
! at the end of the run.
!
! N.B. ray_vec(:,:,:) is allocated ray_vec(nv, nstep_max+1, nray) to accomodate the
! initial point (step zero) plus nstep_max steps beyond that.

    use constants_m, only : rkind

    implicit none

    logical :: write_results_list_directed = .true.

! Time and date vector, get from diagnostics_m
    integer :: date_vector(8)

!  Run label (N.B. should be legal in file name, get from diagnostics_m)
    character(len=60) :: RAYS_run_label = ''

    integer :: number_of_rays, max_number_of_steps, dim_v_vector

! ray data
    real(KIND=rkind), allocatable :: ray_vec(:,:,:) ! nv, nstep_max+1, nray
    real(KIND=rkind), allocatable :: residual(:,:)  ! nstep_max, nray
    integer, allocatable :: npoints(:)              ! nray

! Summary data
    real(KIND=rkind), allocatable :: ray_trace_time(:)
    real(KIND=rkind), allocatable :: end_residuals(:)
    real(KIND=rkind), allocatable :: max_residuals(:)
    real(KIND=rkind), allocatable :: end_ray_parameter(:)
    real(KIND=rkind), allocatable :: start_ray_vec(:,:)
    real(KIND=rkind), allocatable :: end_ray_vec(:,:)
    character(len=60), allocatable :: ray_stop_flag(:)

    real(KIND=rkind)  :: total_trace_time

! Derived type containing same data as above so can have multiple run results in memory
    type run_results

	! Time and date vector, get from diagnostics_m
		integer :: date_vector(8)

	!  Run label (N.B. should be legal in file name, get from diagnostics_m)
		character(len=60) :: RAYS_run_label = ''

		integer :: number_of_rays, max_number_of_steps, dim_v_vector

	! ray data
		real(KIND=rkind), allocatable :: ray_vec(:,:,:) ! nv, nstep_max+1, nray
		real(KIND=rkind), allocatable :: residual(:,:)  ! nstep_max, nray
		integer, allocatable :: npoints(:)              ! nray

	! Summary data
		real(KIND=rkind), allocatable :: ray_trace_time(:)
		real(KIND=rkind), allocatable :: end_residuals(:)
		real(KIND=rkind), allocatable :: max_residuals(:)
		real(KIND=rkind), allocatable :: end_ray_parameter(:)
		real(KIND=rkind), allocatable :: start_ray_vec(:,:)
		real(KIND=rkind), allocatable :: end_ray_vec(:,:)
		character(len=60), allocatable :: ray_stop_flag(:)

		real(KIND=rkind)  :: total_trace_time

    end type run_results

    namelist /ray_results_list/ write_results_list_directed

!****************************************************************************

contains

!****************************************************************************

    subroutine initialize_ray_results_m(read_input)

        use diagnostics_m, only : message_unit, message, text_message, run_label, date_v, &
                                & messages_to_stdout, verbosity
        use ray_init_m, only : nray  ! Number of rays initialized
        use ode_m, only : nv, nstep_max ! dimension of ray vector, max number of steps allowed

        implicit none
        logical, intent(in) :: read_input
 		integer :: input_unit, get_unit_number ! External, free unit finder

		call message(1)
		call text_message('Initializing ray_results_m ', 1)

        if (read_input .eqv. .true.) then
        ! Read and write input namelist
  		  	input_unit = get_unit_number()
            open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
            read(input_unit, ray_results_list)
            close(unit=input_unit)

			allocate (ray_vec(nv, nstep_max+1, nray))
			allocate (residual(nstep_max, nray))
			allocate (npoints(nray))
			allocate (ray_trace_time(nray))
			allocate (end_ray_parameter(nray))
			allocate (end_residuals(nray))
			allocate (max_residuals(nray))
			allocate (start_ray_vec(nv, nray))
			allocate (end_ray_vec(nv, nray))
			allocate (ray_stop_flag(nray))
        end if

! Write input namelist
		if (verbosity >= 0) then
			write(message_unit, ray_results_list)
			if (messages_to_stdout) write(*, ray_results_list)
		end if

        date_vector = date_v
        RAYS_run_label = run_label
        number_of_rays = nray
        max_number_of_steps = nstep_max
        dim_v_vector = nv

        ray_vec = 0.
        residual = 0.
        npoints = 0
        ray_trace_time = 0.
        end_ray_parameter = 0.
        end_residuals = 0.
        max_residuals = 0.
        start_ray_vec = 0.
        end_ray_vec = 0.
        ray_stop_flag = ''


    return
    end subroutine initialize_ray_results_m

!****************************************************************************

    subroutine write_results_LD

    use diagnostics_m, only : run_label

    implicit none

    integer :: results_star_unit, get_unit_number

 !  File name for  output
    character(len=80) :: out_filename

    ! Open fortran ascii file for results output
    results_star_unit = get_unit_number
    out_filename = 'run_results.'//trim(run_label)
    open(unit=results_star_unit, file=trim(out_filename), &
       & action='write', status='replace', form='formatted')

    write (results_star_unit,*) 'RAYS_run_label'
    write (results_star_unit,*) RAYS_run_label
    write (results_star_unit,*) 'date_vector'
    write (results_star_unit,*) date_vector
    write (results_star_unit,*) 'number_of_rays'
    write (results_star_unit,*) number_of_rays
    write (results_star_unit,*) 'max_number_of_steps'
    write (results_star_unit,*) max_number_of_steps
    write (results_star_unit,*) 'dim_v_vector'
    write (results_star_unit,*) dim_v_vector
    write (results_star_unit,*) 'npoints'
    write (results_star_unit,*) npoints

    write (results_star_unit,*) 'total_trace_time'
    write (results_star_unit,*) total_trace_time
    write (results_star_unit,*) 'ray_trace_time'
    write (results_star_unit,*) ray_trace_time
    write (results_star_unit,*) 'end_ray_parameter'
    write (results_star_unit,*) end_ray_parameter
    write (results_star_unit,*) 'end_residuals'
    write (results_star_unit,*) end_residuals
    write (results_star_unit,*) 'max_residuals'
    write (results_star_unit,*) max_residuals
    write (results_star_unit,*) 'ray_stop_flag'
    write (results_star_unit,*) ray_stop_flag
    write (results_star_unit,*) 'start_ray_vec'
    write (results_star_unit,*) start_ray_vec
    write (results_star_unit,*) 'end_ray_vec'
    write (results_star_unit,*) end_ray_vec
    write (results_star_unit,*) 'residual'
    write (results_star_unit,*) residual
    write (results_star_unit,*) 'ray_vec'
    write (results_star_unit,*) ray_vec

    close(unit=results_star_unit)

    end subroutine write_results_LD

!****************************************************************************

    subroutine read_results_LD(in_filename)

	use diagnostics_m, only : run_label, date_v
	use ray_init_m, only : nray  ! Number of rays initialized
	use ode_m, only : nv, nstep_max ! dimension of ray vector, max number of steps allowed

    implicit none

 !  File name for input
    character(len=80), intent(in) :: in_filename
    character(len=80) :: var_name

    integer :: results_star_unit, get_unit_number
    ! Open fortran ascii file for results output

    results_star_unit = get_unit_number()
    open(unit=results_star_unit, file=trim(in_filename), &
       & action='read', status='old', form='formatted')

! Read scalars and date_v which is fixed length -> integer :: date_v(8)

    read (results_star_unit,*) var_name
    if (var_name .ne. 'RAYS_run_label') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) RAYS_run_label

    read (results_star_unit,*) var_name
    if (var_name .ne. 'date_vector') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) date_vector

    read (results_star_unit,*) var_name
    if (var_name .ne. 'number_of_rays') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) number_of_rays

    read (results_star_unit,*) var_name
    if (var_name .ne. 'max_number_of_steps') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) max_number_of_steps

    read (results_star_unit,*) var_name
    if (var_name .ne. 'dim_v_vector') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) dim_v_vector

! If the arrays are already allocated deallocate them and reallocate with dimensions just read

	if (allocated(ray_vec)) call deallocate_ray_results_m
	allocate (ray_vec(nv, nstep_max+1, nray))
	allocate (residual(nstep_max, nray))
	allocate (npoints(nray))
	allocate (ray_trace_time(nray))
	allocate (end_ray_parameter(nray))
	allocate (end_residuals(nray))
	allocate (max_residuals(nray))
	allocate (start_ray_vec(nv, nray))
	allocate (end_ray_vec(nv, nray))
	allocate (ray_stop_flag(nray))


! Read arrays

    read (results_star_unit,*) var_name
    if (var_name .ne. 'npoints') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) npoints

    read (results_star_unit,*) var_name
    if (var_name .ne. 'total_trace_time') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) total_trace_time

    read (results_star_unit,*) var_name
    if (var_name .ne. 'ray_trace_time') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) ray_trace_time

    read (results_star_unit,*) var_name
    if (var_name .ne. 'end_ray_parameter') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) end_ray_parameter

    read (results_star_unit,*) var_name
    if (var_name .ne. 'end_residuals') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) end_residuals

    read (results_star_unit,*) var_name
    if (var_name .ne. 'max_residuals') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) max_residuals

    read (results_star_unit,*) var_name
    if (var_name .ne. 'ray_stop_flag') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) ray_stop_flag

    read (results_star_unit,*) var_name
    if (var_name .ne. 'start_ray_vec') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) start_ray_vec

    read (results_star_unit,*) var_name
    if (var_name .ne. 'end_ray_vec') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) end_ray_vec

    read (results_star_unit,*) var_name
    if (var_name .ne. 'residual') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) residual

    read (results_star_unit,*) var_name
    if (var_name .ne. 'ray_vec') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) ray_vec

    close(unit=results_star_unit)

    end subroutine read_results_LD

!********************************************************************

    subroutine deallocate_ray_results_m
        deallocate (ray_vec)
        deallocate (residual)
        deallocate (npoints)
        deallocate (ray_trace_time)
        deallocate (end_ray_parameter)
        deallocate (end_residuals)
        deallocate (max_residuals)
        deallocate (start_ray_vec)
        deallocate (end_ray_vec)
        deallocate (ray_stop_flag)

        return
    end subroutine deallocate_ray_results_m

!****************************************************************************

    subroutine results_type_from_module_data_m(this_run)

    implicit none
    type (run_results), intent(out) :: this_run

	allocate (this_run%ray_vec(dim_v_vector, max_number_of_steps+1, number_of_rays))
	allocate (this_run%residual(max_number_of_steps, number_of_rays))
	allocate (this_run%npoints(number_of_rays))
	allocate (this_run%ray_trace_time(number_of_rays))
	allocate (this_run%end_ray_parameter(number_of_rays))
	allocate (this_run%end_residuals(number_of_rays))
	allocate (this_run%max_residuals(number_of_rays))
	allocate (this_run%start_ray_vec(dim_v_vector, number_of_rays))
	allocate (this_run%end_ray_vec(dim_v_vector, number_of_rays))
	allocate (this_run%ray_stop_flag(number_of_rays))

    this_run%RAYS_run_label = RAYS_run_label
    this_run%date_vector = date_vector
    this_run%number_of_rays = number_of_rays
    this_run%max_number_of_steps = max_number_of_steps
    this_run%dim_v_vector = dim_v_vector
    this_run%npoints = npoints

    this_run%total_trace_time = total_trace_time
    this_run%ray_trace_time = ray_trace_time
    this_run%end_ray_parameter = end_ray_parameter
    this_run%end_residuals = end_residuals
    this_run%max_residuals = max_residuals
    this_run%start_ray_vec = start_ray_vec
    this_run%end_ray_vec = end_ray_vec
    this_run%residual = residual
    this_run%ray_vec = ray_vec

    end subroutine results_type_from_module_data_m



 end module ray_results_m

!********************************************************************
