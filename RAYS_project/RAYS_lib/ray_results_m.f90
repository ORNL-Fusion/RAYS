 module ray_results_m

! A module that contains results from ray tracing for storage in memory or for writing
! at the end of the run.  The is a derived type 'run_results' that can accept the results
! module data so we can have multiple instances (maybe useful for analysis). And there are
! routines to write the data to files in either netCDF or list-directed ascii format.
!
! N.B. ray_vec(:,:,:) is allocated ray_vec(nv, max_number_of_points, nray) where
! max_number_of_points = nstep_max+1 to accomodate the initial point (step zero)
! plus nstep_max steps beyond that.  Save for residual(max_number_of_points, nray)
! and presumably the residual at step 0 is 0.
!
! Note also that nstep_max may be much bigger than the actual number of steps taken by
! any of the rays.  We typically set nstep_max large so we won't have to worry too much
! about how big it is.  However when we write out the results to file there is no point in
! keeping the superflous points.  So on output we just keep up to the largest value of
! npoints(:).  Something to be alert to is that when data is read back in and put into the
! module, max_number_of_points will equal 'actual_max_npoints' = the largest value of
! npoints(:).

    use constants_m, only : rkind

    implicit none

    logical :: write_results_list_directed = .false.
    logical :: write_results_netCDF = .false.

! Time and date vector, get from diagnostics_m
    integer :: date_vector(8)

!  Run label (N.B. should be legal in a file name, get from diagnostics_m)
    character(len=60) :: RAYS_run_label = ''

    integer :: number_of_rays, max_number_of_points, dim_v_vector

! ray data
    real(KIND=rkind), allocatable :: ray_vec(:,:,:) ! nv, max_number_of_points=nstep_max+1, nray
    real(KIND=rkind), allocatable :: residual(:,:)  ! nstep_max, nray
    integer, allocatable :: npoints(:)              ! nray

! Summary data
    real(KIND=rkind), allocatable :: initial_ray_power(:) ! nray
    real(KIND=rkind), allocatable :: ray_trace_time(:)    ! nray
    real(KIND=rkind), allocatable :: end_residuals(:)     ! nray
    real(KIND=rkind), allocatable :: max_residuals(:)     ! nray
    real(KIND=rkind), allocatable :: end_ray_parameter(:) ! nray
    real(KIND=rkind), allocatable :: start_ray_vec(:,:)   ! nv, nray
    real(KIND=rkind), allocatable :: end_ray_vec(:,:)     ! nv, nray
    character(len=60), allocatable :: ray_stop_flag(:)    ! nray

    real(KIND=rkind)  :: total_trace_time

! Derived type containing same data as above so can have multiple run results in memory
    type run_results

	! Time and date vector, get from diagnostics_m
		integer :: date_vector(8)

	!  Run label (N.B. should be legal in file name, get from diagnostics_m)
		character(len=60) :: RAYS_run_label = ''

		integer :: number_of_rays, max_number_of_points, dim_v_vector

	! ray data
		real(KIND=rkind), allocatable :: ray_vec(:,:,:) ! nv, nstep_max+1, nray
		real(KIND=rkind), allocatable :: residual(:,:)  ! nstep_max+1, nray
		integer, allocatable :: npoints(:)              ! nray

	! Summary data
		real(KIND=rkind), allocatable :: initial_ray_power(:) ! nray
		real(KIND=rkind), allocatable :: ray_trace_time(:)    ! nray
		real(KIND=rkind), allocatable :: end_residuals(:)     ! nray
		real(KIND=rkind), allocatable :: max_residuals(:)     ! nray
		real(KIND=rkind), allocatable :: end_ray_parameter(:) ! nray
		real(KIND=rkind), allocatable :: start_ray_vec(:,:)   ! nv, nray
		real(KIND=rkind), allocatable :: end_ray_vec(:,:)     ! nv, nray
		character(len=60), allocatable :: ray_stop_flag(:)    ! nray

		real(KIND=rkind)  :: total_trace_time

	contains
		procedure :: from_module, to_module, read_results_instance_NC

    end type run_results

    namelist /ray_results_list/ write_results_list_directed, write_results_netCDF

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

        dim_v_vector = nv
        number_of_rays = nray
        max_number_of_points = nstep_max+1

        if (read_input .eqv. .true.) then
        ! Read and write input namelist
  		  	input_unit = get_unit_number()
            open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
            read(input_unit, ray_results_list)
            close(unit=input_unit)

			allocate (ray_vec(nv, max_number_of_points, nray))
			allocate (residual(max_number_of_points, nray))
			allocate (npoints(nray))
			allocate (initial_ray_power(nray))
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

        ray_vec = 0.
        residual = 0.
        npoints = 0
        initial_ray_power = 0.
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

  subroutine write_results_NC
! writes a netCDF file containing all module data

    use diagnostics_m, only : run_label
    use netcdf

    implicit none

    integer :: actual_max_npoints

! netCDF Declarations
    integer :: ncid
! Declarations: dimensions
    integer, parameter :: n_dims = 5
    integer :: d8, d60
    integer :: number_of_rays_id, max_number_of_points_id, dim_v_vector_id, d8_id, d60_id

! Declarations: variables
    integer, parameter :: n_vars =  12
    integer :: date_vector_id, ray_vec_id, residual_id, npoints_id, initial_ray_power_id,&
             & ray_trace_time_id, end_residuals_id, max_residuals_id, end_ray_parameter_id,&
             & start_ray_vec_id, end_ray_vec_id, ray_stop_flag_id, total_trace_time_id

 !  File name for  output
    character(len=80) :: out_filename

    out_filename = 'run_results.'//trim(run_label)//'.nc'

!   Open NC file
    call check( nf90_create(trim(out_filename), nf90_clobber, ncid) )

    actual_max_npoints = maxval(npoints) ! see discussion above

!   Define NC dimensions
    call check( nf90_def_dim(ncid, 'number_of_rays', number_of_rays, number_of_rays_id))
    call check( nf90_def_dim(ncid, 'max_number_of_points', actual_max_npoints, max_number_of_points_id))
    call check( nf90_def_dim(ncid, 'dim_v_vector', dim_v_vector, dim_v_vector_id))
    call check( nf90_def_dim(ncid, 'd8', 8, d8_id))
    call check( nf90_def_dim(ncid, 'd60', 60, d60_id))

! Define NC variables
    call check( nf90_def_var(ncid, 'date_vector', NF90_INT, [d8_id], date_vector_id))
    call check( nf90_def_var(ncid, 'ray_vec', NF90_DOUBLE, [dim_v_vector_id,max_number_of_points_id,number_of_rays_id], ray_vec_id))
    call check( nf90_def_var(ncid, 'residual', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], residual_id))
    call check( nf90_def_var(ncid, 'npoints', NF90_INT, [number_of_rays_id], npoints_id))
    call check( nf90_def_var(ncid, 'initial_ray_power', NF90_FLOAT, [number_of_rays_id], initial_ray_power_id))
    call check( nf90_def_var(ncid, 'ray_trace_time', NF90_FLOAT, [number_of_rays_id], ray_trace_time_id))
    call check( nf90_def_var(ncid, 'end_residuals', NF90_FLOAT, [number_of_rays_id], end_residuals_id))
    call check( nf90_def_var(ncid, 'max_residuals', NF90_FLOAT, [number_of_rays_id], max_residuals_id))
    call check( nf90_def_var(ncid, 'end_ray_parameter', NF90_FLOAT, [number_of_rays_id], end_ray_parameter_id))
    call check( nf90_def_var(ncid, 'start_ray_vec', NF90_FLOAT, [dim_v_vector_id,number_of_rays_id], start_ray_vec_id))
    call check( nf90_def_var(ncid, 'end_ray_vec', NF90_FLOAT, [dim_v_vector_id,number_of_rays_id], end_ray_vec_id))
    call check( nf90_def_var(ncid, 'ray_stop_flag', NF90_CHAR, [d60_id,number_of_rays_id], ray_stop_flag_id))
    call check( nf90_def_var(ncid, 'total_trace_time', NF90_FLOAT, total_trace_time_id))

! Put global attributes
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'RAYS_run_label', RAYS_run_label))

call check( nf90_enddef(ncid))

! Put NC variables
    call check( nf90_put_var(ncid, date_vector_id, date_vector))
    call check( nf90_put_var(ncid, ray_vec_id, ray_vec(:,1:actual_max_npoints,:)))
    call check( nf90_put_var(ncid, residual_id, residual(1:actual_max_npoints,:)))
    call check( nf90_put_var(ncid, npoints_id, npoints))
    call check( nf90_put_var(ncid, initial_ray_power_id, initial_ray_power))
    call check( nf90_put_var(ncid, ray_trace_time_id, ray_trace_time))
    call check( nf90_put_var(ncid, end_residuals_id, end_residuals))
    call check( nf90_put_var(ncid, max_residuals_id, max_residuals))
    call check( nf90_put_var(ncid, end_ray_parameter_id, end_ray_parameter))
    call check( nf90_put_var(ncid, start_ray_vec_id, start_ray_vec))
    call check( nf90_put_var(ncid, end_ray_vec_id, end_ray_vec))
    call check( nf90_put_var(ncid, ray_stop_flag_id, ray_stop_flag))
    call check( nf90_put_var(ncid, total_trace_time_id, total_trace_time))

!   Close the NC file
    call check( nf90_close(ncid) )

  end subroutine write_results_NC

!****************************************************************************

  subroutine read_results_instance_NC(this, in_filename)
! Reads a netCDF file containing all module data into a run_results instsance 'this'

    use diagnostics_m, only : run_label
    use netcdf

    implicit none
    class (run_results) ::  this

! netCDF Declarations
    integer :: ncid
    integer, parameter :: n_dims = 5
    integer :: number_of_rays_id, max_number_of_points_id, dim_v_vector_id, d8_id, d60_id
    integer :: d8, d60 ! Need these declarations, they are not in the module data

! Declarations: variables
    integer, parameter :: n_vars =  13
    integer :: date_vector_id, ray_vec_id, residual_id, npoints_id, initial_ray_power_id,&
             & ray_trace_time_id, end_residuals_id, max_residuals_id, end_ray_parameter_id,&
             & start_ray_vec_id, end_ray_vec_id, ray_stop_flag_id, total_trace_time_id

 !  File name for  output
    character(len=80) :: in_filename

! If this instance is already allocated deallocate its arrays
  if (allocated(this%ray_vec)) deallocate(this%ray_vec)
  if (allocated(this%residual)) deallocate(this%residual)
  if (allocated(this%npoints)) deallocate(this%npoints)
  if (allocated(this%initial_ray_power)) deallocate(this%initial_ray_power)
  if (allocated(this%ray_trace_time)) deallocate(this%ray_trace_time)
  if (allocated(this%end_residuals)) deallocate(this%end_residuals)
  if (allocated(this%max_residuals)) deallocate(this%max_residuals)
  if (allocated(this%end_ray_parameter)) deallocate(this%end_ray_parameter)
  if (allocated(this%start_ray_vec)) deallocate(this%start_ray_vec)
  if (allocated(this%end_ray_vec)) deallocate(this%end_ray_vec)
  if (allocated(this%ray_stop_flag)) deallocate(this%ray_stop_flag)

!   Open NC file
    call check( nf90_open(in_filename, nf90_nowrite, ncid) )

!   Inquire NC dimensions
    call check( nf90_inq_dimid(ncid, 'number_of_rays', number_of_rays_id))
    call check(nf90_inquire_dimension(ncid, number_of_rays_id, len=this%number_of_rays))
    call check( nf90_inq_dimid(ncid, 'max_number_of_points', max_number_of_points_id))
    call check(nf90_inquire_dimension(ncid, max_number_of_points_id, len=this%max_number_of_points))
    call check( nf90_inq_dimid(ncid, 'dim_v_vector', dim_v_vector_id))
    call check(nf90_inquire_dimension(ncid, dim_v_vector_id, len=this%dim_v_vector))

!   Allocate arrays
    allocate(this%ray_vec(this%dim_v_vector,this%max_number_of_points,this%number_of_rays))
    allocate(this%residual(this%max_number_of_points,this%number_of_rays))
    allocate(this%npoints(this%number_of_rays))
    allocate(this%initial_ray_power(this%number_of_rays))
    allocate(this%ray_trace_time(this%number_of_rays))
    allocate(this%end_residuals(this%number_of_rays))
    allocate(this%max_residuals(this%number_of_rays))
    allocate(this%end_ray_parameter(this%number_of_rays))
    allocate(this%start_ray_vec(this%dim_v_vector,this%number_of_rays))
    allocate(this%end_ray_vec(this%dim_v_vector,this%number_of_rays))
    allocate(this%ray_stop_flag(this%number_of_rays))

!   Get varids and variable values
    call check(nf90_inq_varid(ncid, 'date_vector', date_vector_id))
    call check( nf90_get_var(ncid, date_vector_id, this%date_vector))
    call check(nf90_inq_varid(ncid, 'ray_vec', ray_vec_id))
    call check( nf90_get_var(ncid, ray_vec_id, this%ray_vec))
    call check(nf90_inq_varid(ncid, 'residual', residual_id))
    call check( nf90_get_var(ncid, residual_id, this%residual))
    call check(nf90_inq_varid(ncid, 'npoints', npoints_id))
    call check( nf90_get_var(ncid, npoints_id, this%npoints))
    call check(nf90_inq_varid(ncid, 'npoints', npoints_id))
    call check(nf90_inq_varid(ncid, 'initial_ray_power', initial_ray_power_id))
    call check( nf90_get_var(ncid, initial_ray_power_id, this%initial_ray_power))
    call check(nf90_inq_varid(ncid, 'ray_trace_time', ray_trace_time_id))
    call check( nf90_get_var(ncid, ray_trace_time_id, this%ray_trace_time))
    call check(nf90_inq_varid(ncid, 'end_residuals', end_residuals_id))
    call check( nf90_get_var(ncid, end_residuals_id, this%end_residuals))
    call check(nf90_inq_varid(ncid, 'max_residuals', max_residuals_id))
    call check( nf90_get_var(ncid, max_residuals_id, this%max_residuals))
    call check(nf90_inq_varid(ncid, 'end_ray_parameter', end_ray_parameter_id))
    call check( nf90_get_var(ncid, end_ray_parameter_id, this%end_ray_parameter))
    call check(nf90_inq_varid(ncid, 'start_ray_vec', start_ray_vec_id))
    call check( nf90_get_var(ncid, start_ray_vec_id, this%start_ray_vec))
    call check(nf90_inq_varid(ncid, 'end_ray_vec', end_ray_vec_id))
    call check( nf90_get_var(ncid, end_ray_vec_id, this%end_ray_vec))
    call check(nf90_inq_varid(ncid, 'ray_stop_flag', ray_stop_flag_id))
    call check( nf90_get_var(ncid, ray_stop_flag_id, this%ray_stop_flag))
    call check(nf90_inq_varid(ncid, 'total_trace_time', total_trace_time_id))
    call check( nf90_get_var(ncid, total_trace_time_id, this%total_trace_time))

! Get global attributes
    call check( nf90_get_att(ncid, NF90_GLOBAL, 'RAYS_run_label', this%RAYS_run_label))

!   Close the NC file
    call check( nf90_close(ncid) )

  end subroutine read_results_instance_NC

!****************************************************************************

  subroutine check(status)
    use netcdf
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check

!****************************************************************************

    subroutine write_results_LD

    use diagnostics_m, only : run_label

    implicit none

    integer :: results_star_unit, get_unit_number

 !  File name for  output
    character(len=80) :: out_filename

    ! Open fortran ascii file for results output
    results_star_unit = get_unit_number()
    out_filename = 'run_results.'//trim(run_label)
    open(unit=results_star_unit, file=trim(out_filename), &
       & action='write', status='replace', form='formatted')

    write (results_star_unit,*) 'RAYS_run_label'
    write (results_star_unit,*) RAYS_run_label
    write (results_star_unit,*) 'date_vector'
    write (results_star_unit,*) date_vector
    write (results_star_unit,*) 'number_of_rays'
    write (results_star_unit,*) number_of_rays
    write (results_star_unit,*) 'max_number_of_points'
    write (results_star_unit,*) max_number_of_points
    write (results_star_unit,*) 'dim_v_vector'
    write (results_star_unit,*) dim_v_vector
    write (results_star_unit,*) 'npoints'
    write (results_star_unit,*) npoints

    write (results_star_unit,*) 'total_trace_time'
    write (results_star_unit,*) total_trace_time
    write (results_star_unit,*) 'initial_ray_power'
    write (results_star_unit,*) initial_ray_power
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
    if (var_name .ne. 'max_number_of_points') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) max_number_of_points

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
	allocate (initial_ray_power(nray))
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
    if (var_name .ne. 'initial_ray_power') then
    	write(*,*) 'read_results_LD: inconsistent variable name = ', var_name
    	stop
    end if
    read (results_star_unit,*) initial_ray_power

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

    subroutine allocate_ray_results_m
		allocate (ray_vec(dim_v_vector, max_number_of_points, number_of_rays))
		allocate (residual(max_number_of_points, number_of_rays))
		allocate (npoints(number_of_rays))
		allocate (initial_ray_power(number_of_rays))
		allocate (ray_trace_time(number_of_rays))
		allocate (end_ray_parameter(number_of_rays))
		allocate (end_residuals(number_of_rays))
		allocate (max_residuals(number_of_rays))
		allocate (start_ray_vec(dim_v_vector, number_of_rays))
		allocate (end_ray_vec(dim_v_vector, number_of_rays))
		allocate (ray_stop_flag(number_of_rays))

        return
    end subroutine allocate_ray_results_m

!********************************************************************

    subroutine deallocate_ray_results_m
        deallocate (ray_vec)
        deallocate (residual)
        deallocate (npoints)
        deallocate (initial_ray_power)
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

    subroutine from_module(this)

    implicit none

    class (run_results) ::  this

! 	allocate ( this%ray_vec(dim_v_vector, max_number_of_points, number_of_rays))
! 	allocate ( this%residual(max_number_of_points, number_of_rays))
! 	allocate ( this%npoints(number_of_rays))
! 	allocate ( this%initial_ray_power(number_of_rays))
! 	allocate ( this%ray_trace_time(number_of_rays))
! 	allocate ( this%end_ray_parameter(number_of_rays))
! 	allocate ( this%end_residuals(number_of_rays))
! 	allocate ( this%max_residuals(number_of_rays))
! 	allocate ( this%start_ray_vec(dim_v_vector, number_of_rays))
! 	allocate ( this%end_ray_vec(dim_v_vector, number_of_rays))
! 	allocate ( this%ray_stop_flag(number_of_rays))

     this%RAYS_run_label = RAYS_run_label
     this%date_vector = date_vector
     this%number_of_rays = number_of_rays
     this%max_number_of_points = max_number_of_points
     this%dim_v_vector = dim_v_vector
     this%npoints = npoints

     this%total_trace_time = total_trace_time
     this%initial_ray_power = initial_ray_power
     this%ray_trace_time = ray_trace_time
     this%end_ray_parameter = end_ray_parameter
     this%end_residuals = end_residuals
     this%max_residuals = max_residuals
     this%start_ray_vec = start_ray_vec
     this%end_ray_vec = end_ray_vec
     this%residual = residual
     this%ray_vec = ray_vec
     this%ray_stop_flag = ray_stop_flag

    end subroutine from_module

!****************************************************************************

    subroutine to_module(this)

    implicit none

    class (run_results) ::  this

! 	allocate ( this%ray_vec(dim_v_vector, max_number_of_points, number_of_rays))
! 	allocate ( this%residual(max_number_of_points, number_of_rays))
! 	allocate ( this%npoints(number_of_rays))
! 	allocate ( this%initial_ray_power(number_of_rays))
! 	allocate ( this%ray_trace_time(number_of_rays))
! 	allocate ( this%end_ray_parameter(number_of_rays))
! 	allocate ( this%end_residuals(number_of_rays))
! 	allocate ( this%max_residuals(number_of_rays))
! 	allocate ( this%start_ray_vec(dim_v_vector, number_of_rays))
! 	allocate ( this%end_ray_vec(dim_v_vector, number_of_rays))
! 	allocate ( this%ray_stop_flag(number_of_rays))

    call deallocate_ray_results_m

     RAYS_run_label = this%RAYS_run_label
     date_vector = this%date_vector
     number_of_rays = this%number_of_rays
     max_number_of_points = this%max_number_of_points
     dim_v_vector = this%dim_v_vector
     npoints = this%npoints

     total_trace_time = this%total_trace_time
     initial_ray_power = this%initial_ray_power
     ray_trace_time = this%ray_trace_time
     end_ray_parameter = this%end_ray_parameter
     end_residuals = this%end_residuals
     max_residuals = this%max_residuals
     start_ray_vec = this%start_ray_vec
     end_ray_vec = this%end_ray_vec
     residual = this%residual
     ray_vec = this%ray_vec
     ray_stop_flag = this%ray_stop_flag

    end subroutine to_module

 end module ray_results_m
