module scanner_m

! Modifies data in modules to implement a scan of ray runs with varying input data.
! scan_parameter selects which parameter is to be varied.  Presently only supported are:
!   ode_m -> 'ds'
! scan_algorithm selects the algorithm by which the parameter is varied.  Presently only
! supported are: 'fixed_increment'
!
! Reads standard rays.in namelist file with added /scanner_list/.  Serves as template for
! all data that is not varied during scan.
!
! N.B.  For now we only aggregate data for ray # 1 per run.  A parameter scan varying only
! simulation parameters might reasonably only consist of one ray per run.


    use constants_m, only : rkind
    
    implicit none

    character(len=60) :: scan_parameter
    character(len=60) :: scan_algorithm
    character(len=60) :: scan_id
    integer :: n_runs
    real(KIND=rkind), allocatable :: p_values(:)
    character(len=60), allocatable :: file_name_suffix(:)
    integer :: scan_date_v(8)
     
   
! data for fixed increment scan algorithm
    real(KIND=rkind) :: p_start, p_incr

! Scan summary data
    real(KIND=rkind), allocatable :: trace_time_run(:)
    real(KIND=rkind), allocatable :: end_ray_param_run(:)
    real(KIND=rkind), allocatable :: end_resid_run(:) 
    real(KIND=rkind), allocatable :: max_resid_run(:) 
    real(KIND=rkind), allocatable :: end_ray_vec_run(:,:)
    character(len=60), allocatable :: ray_stop_flag_run(:)

    real(KIND=rkind) :: scan_trace_time

 namelist /scanner_list/ &
     & scan_id, scan_parameter, scan_algorithm, n_runs, p_start, p_incr
     
!********************************************************************

contains

!********************************************************************

  subroutine initialize_scanner_m(read_input)

    use constants_m, only : input_unit
    use diagnostics_m, only : message_unit, verbosity, run_label
	use ode_m, only : nv ! dimension of ray vector
   
    implicit none
    logical, intent(in) :: read_input
    integer :: i_run
    character (len=4) :: chr_iter_number

        
    write(*,*) 'initialize_scanner_m'

    if (read_input .eqv. .true.) then    
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, scanner_list)
        close(unit=input_unit)
        write(message_unit, scanner_list)

		allocate( p_values(n_runs) )
		allocate( file_name_suffix(n_runs) )

		allocate (trace_time_run(n_runs))
		allocate (end_ray_param_run(n_runs))
		allocate (end_resid_run(n_runs))
		allocate (max_resid_run(n_runs))
		allocate (end_ray_vec_run(nv, n_runs))
		allocate (ray_stop_flag_run(n_runs))

    end if
    write(*,*) 'n_runs = ', n_runs

        p_values = 0.
        file_name_suffix = ''

        trace_time_run = 0.
        end_ray_param_run = 0.
        end_resid_run = 0.
        max_resid_run = 0.
        end_ray_vec_run = 0.
        ray_stop_flag_run = ''        
        scan_trace_time = 0.
        
!   Calculate parameter values and 
    algorithm: select case (trim(scan_algorithm))

       case ('fixed_increment')

		do i_run = 1, n_runs	
		   p_values(i_run) = p_start +  i_run * p_incr
		   write (chr_iter_number, '(I4)') i_run
    write(*,*) 'chr_iter_number = ', adjustl(trim(chr_iter_number))
		   file_name_suffix(i_run) = 'run_'//adjustl(trim(chr_iter_number))
		end do

	   case default
		  write(0,*) 'initialize_scanner_m: unknown scan algorithm = ', scan_algorithm
		  stop 1
    
    end select algorithm

  write(*,*) 'p_values = ', p_values
  write(*,*) ' '
  write(*,*) 'file_name_suffix = ', file_name_suffix

    return
  end subroutine initialize_scanner_m

!********************************************************************

 subroutine update_scan_parameter(i_run)

    use diagnostics_m, only : message_unit, verbosity, run_label
    use ode_m, only : ds
    
    implicit none
    
    integer, intent(in) :: i_run

    param: select case (trim(scan_parameter))

       case ('ds')
       	ds = p_values(i_run)
       	run_label = trim(file_name_suffix(i_run))

	   case default
		  write(0,*) 'initialize_scanner_m: unknown scan parameter = ', scan_parameter
		  stop 1
    
    end select param

 end subroutine update_scan_parameter

!********************************************************************

 subroutine aggregate_run_data(i_run)
! N.B.  For now we only aggregate data for ray # 1 per run.  See above.

    use diagnostics_m, only : message_unit, verbosity, run_label
    use ray_results_m, only : end_residuals, max_residuals, end_ray_parameter, end_ray_vec,&
                            & run_trace_time, ray_trace_time, ray_stop_flag



    implicit none
    
    integer, intent(in) :: i_run

    trace_time_run(i_run) = ray_trace_time(1)
	end_ray_param_run(i_run) = end_ray_parameter(1)
	end_resid_run(i_run) = end_residuals(1)
	max_resid_run(i_run) = max_residuals(1)
	end_ray_vec_run(:, i_run) = end_ray_vec(:, 1)
	ray_stop_flag_run(i_run) = ray_stop_flag(1)     
   
    scan_trace_time = scan_trace_time + run_trace_time
    
 end subroutine aggregate_run_data

!********************************************************************

    subroutine write_scan_summary

    use diagnostics_m, only : message_unit, run_label   
 	use ode_m, only : nv ! dimension of ray vector
    use ray_results_m, only : end_residuals, max_residuals, end_ray_parameter, end_ray_vec,&
                            & ray_trace_time
        
    implicit none
    
    integer :: scan_star_unit
    
 !  File name for  output
    character(len=80) :: out_filename
   
    ! Open fortran ascii file for results output
    scan_star_unit = 59 
    out_filename = 'scan_summary.'//trim(scan_id)
    open(unit=scan_star_unit, file=trim(out_filename), &
       & action='write', status='replace', form='formatted')     

    write (scan_star_unit,*) 'scan_id'
    write (scan_star_unit,*) scan_id
    write (scan_star_unit,*) 'scan_parameter'
    write (scan_star_unit,*) scan_parameter
    write (scan_star_unit,*) 'scan_date_v'
    write (scan_star_unit,*) scan_date_v
    write (scan_star_unit,*) 'p_values'
    write (scan_star_unit,*) p_values
    write (scan_star_unit,*) 'dim_v_vector'
    write (scan_star_unit,*) nv

    write (scan_star_unit,*) 'trace_time_run'
    write (scan_star_unit,*) trace_time_run
    write (scan_star_unit,*) 'end_ray_param_run'
    write (scan_star_unit,*) end_ray_param_run
    write (scan_star_unit,*) 'end_resid_run'
    write (scan_star_unit,*) end_resid_run
    write (scan_star_unit,*) 'max_resid_run'
    write (scan_star_unit,*) max_resid_run
    write (scan_star_unit,*) 'ray_stop_flag_run'
    write (scan_star_unit,*) ray_stop_flag_run
    write (scan_star_unit,*) 'end_ray_vec_run'
    write (scan_star_unit,*) end_ray_vec_run

    close(unit=scan_star_unit)

    end subroutine write_scan_summary

!********************************************************************

    
!********************************************************************

    subroutine deallocate_scanner_m
		if (allocated(p_values)) then
			deallocate( p_values )
			deallocate( file_name_suffix )
			deallocate (trace_time_run)
			deallocate (end_ray_param_run)
			deallocate (end_resid_run)
			deallocate (max_resid_run)
			deallocate (end_ray_vec_run)
			deallocate (ray_stop_flag_run)
		end if
		return
    end subroutine deallocate_scanner_m

end module scanner_m
