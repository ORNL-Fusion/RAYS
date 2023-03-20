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

    real(KIND=rkind)  :: run_trace_time

    namelist /ray_results_list/ write_results_list_directed

!****************************************************************************

contains

!****************************************************************************

    subroutine initialize_ray_results_m(read_input)

        use diagnostics_m, only : message_unit, run_label, date_v
        use ray_init_m, only : nray  ! Number of rays initialized
        use ode_m, only : nv, nstep_max ! dimension of ray vector, max number of steps allowed
 
        implicit none
        logical, intent(in) :: read_input
 		integer :: input_unit, get_unit_number ! External, free unit finder   
        
        write(*,*) 'initialize_ray_results'
     
        if (read_input .eqv. .true.) then    
        ! Read and write input namelist
  		  	input_unit = get_unit_number()
            open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
            read(input_unit, ray_results_list)
            close(unit=input_unit)
            write(message_unit, ray_results_list)

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

    write (results_star_unit,*) 'run_trace_time'
    write (results_star_unit,*) run_trace_time
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
    if (var_name .ne. 'run_trace_time') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) run_trace_time

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
     
 end module ray_results_m
