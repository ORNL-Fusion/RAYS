module scanner_m

! Modifies data in modules to implement a scan of ray runs with varying input data.
! scan_parameter selects which parameter is to be varied.  Presently only supported are:
!   ode_m -> 'ds'
!   openmp_m -> num_threads
! scan_algorithm selects the algorithm by which the parameter is varied.  Presently only
! supported are:
!   ode_m -> 'ds' -> 'fixed_increment', 'pwr_of_2', 'integer_divide', 'algorithm_1'
!   openmp_m -> num_threads -> 'increment_num_threads', 'double_num_threads'
!
! Reads standard rays.in namelist file with added /scanner_list/.  Serves as template for
! all data that is not varied during scan.
!
! N.B.  For now we only aggregate data for ray # 1 per run.  A parameter scan varying only
! simulation parameters might reasonably only consist of one ray per run.
! Except that trace_time_run is for the whole run.


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

! data for pwr_of_2 scan algorithm and integer_divide
    real(KIND=rkind) :: p_max
    integer :: max_divide
    real(KIND=rkind) :: delta = 1.e-14

! data for algorithm_1
    real(KIND=rkind) :: S_max
    integer :: n_max, k_factor

! Scan summary data
    real(KIND=rkind), allocatable :: trace_time_run(:)
    real(KIND=rkind), allocatable :: end_ray_param_run(:)
    real(KIND=rkind), allocatable :: end_resid_run(:)
    real(KIND=rkind), allocatable :: max_resid_run(:)
    real(KIND=rkind), allocatable :: start_ray_vec_run(:,:)
    real(KIND=rkind), allocatable :: end_ray_vec_run(:,:)
    character(len=60), allocatable :: ray_stop_flag_run(:)

    real(KIND=rkind) :: scan_trace_time

 namelist /scanner_list/ &
     & scan_id, scan_parameter, scan_algorithm, n_runs, p_start, p_incr, p_max, max_divide,&
     & S_max, n_max, k_factor

!********************************************************************

contains

!********************************************************************

  subroutine initialize_scanner_m(read_input)

    use diagnostics_m, only : message_unit, verbosity, run_label
	use ode_m, only : nv ! dimension of ray vector

    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder

    integer :: i_run
    character (len=4) :: chr_iter_number


    write(*,*) 'initialize_scanner_m'

    if (read_input .eqv. .true.) then
   		input_unit = get_unit_number()
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
		allocate (start_ray_vec_run(nv, n_runs))
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
        start_ray_vec_run = 0.
        end_ray_vec_run = 0.
        ray_stop_flag_run = ''
        scan_trace_time = 0.

!   Calculate parameter values and
    algorithm: select case (trim(scan_algorithm))

       case ('fixed_increment')
		do i_run = 1, n_runs
		   p_values(i_run) = p_start +  (i_run - 1) * p_incr
		   write (chr_iter_number, '(I4)') i_run
		   file_name_suffix(i_run) = 'run_'//adjustl(trim(chr_iter_number))
		end do

       case ('pwr_of_2') ! p_max is always (an integral number of p_values) - delta
		do i_run = 1, n_runs
		   p_values(i_run) = p_max/2.**(n_runs-i_run) - delta
		   write (chr_iter_number, '(I4)') i_run
		   file_name_suffix(i_run) = 'run_'//adjustl(trim(chr_iter_number))
		end do

       case ('integer_divide') ! p_max is always (an integral number of p_values) - delta
		do i_run = 1, n_runs
		   p_values(i_run) = p_max/(max_divide - i_run +1)
		   write (chr_iter_number, '(I4)') i_run
		   file_name_suffix(i_run) = 'run_'//adjustl(trim(chr_iter_number))
		end do

       case ('algorithm_1') ! p_max is always (an integral number of p_values) - delta
		do i_run = 1, n_runs
		   p_values(i_run) = S_max/(n_max + k_factor*(n_runs-i_run)) - delta
		   write (chr_iter_number, '(I4)') i_run
		   file_name_suffix(i_run) = 'run_'//adjustl(trim(chr_iter_number))
		end do

       case ('increment_num_threads')
		do i_run = 1, n_runs
		   p_values(i_run) = i_run
		   write (chr_iter_number, '(I4)') i_run
		   file_name_suffix(i_run) = 'run_'//adjustl(trim(chr_iter_number))
		end do

       case ('double_num_threads')
		do i_run = 1, n_runs
		   p_values(i_run) = 2.**(i_run-1)
		   write (chr_iter_number, '(I4)') i_run
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
    use openmp_m, only : num_threads

    implicit none

    integer, intent(in) :: i_run
    integer :: stat

    param: select case (trim(scan_parameter))

       case ('ds')
       	ds = p_values(i_run)
       	run_label = trim(file_name_suffix(i_run))

       case ('double_num_threads', 'increment_num_threads')
       	num_threads = nint(p_values(i_run))
       	run_label = trim(file_name_suffix(i_run))
!
!        case ('increment_num_threads')
!        	num_threads = nint(p_values(i_run))
!        	run_label = trim(file_name_suffix(i_run))

	   case default
		  write(0,*) 'initialize_scanner_m: unknown scan parameter = ', scan_parameter
		  stop 1

    end select param

 end subroutine update_scan_parameter

!********************************************************************

 subroutine aggregate_run_data(i_run)
! N.B.  For now we only aggregate data for ray # 1 per run.  See above.

    use diagnostics_m, only : message_unit, verbosity, run_label
    use ray_results_m, only : end_residuals, max_residuals, end_ray_parameter, start_ray_vec,&
                            & end_ray_vec, total_trace_time, ray_trace_time, ray_stop_flag



    implicit none

    integer, intent(in) :: i_run

    trace_time_run(i_run) = total_trace_time
	end_ray_param_run(i_run) = end_ray_parameter(1)
	end_resid_run(i_run) = end_residuals(1)
	max_resid_run(i_run) = max_residuals(1)
	start_ray_vec_run(:, i_run) = start_ray_vec(:, 1)
	end_ray_vec_run(:, i_run) = end_ray_vec(:, 1)
	ray_stop_flag_run(i_run) = ray_stop_flag(1)

    scan_trace_time = scan_trace_time + total_trace_time

 end subroutine aggregate_run_data

!********************************************************************

    subroutine write_scan_summary

    use diagnostics_m, only : message_unit, run_label
 	use ode_m, only : nv ! dimension of ray vector
    use ray_results_m, only : end_residuals, max_residuals, end_ray_parameter, end_ray_vec,&
                            & ray_trace_time

    implicit none

    integer :: scan_star_unit, get_unit_number ! External, free unit finder

 !  File name for  output
    character(len=80) :: out_filename

    ! Open fortran ascii file for results output
    scan_star_unit = get_unit_number()
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
    write (scan_star_unit,*) 'start_ray_vec_run'
    write (scan_star_unit,*) start_ray_vec_run
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
			deallocate (start_ray_vec_run)
			deallocate (end_ray_vec_run)
			deallocate (ray_stop_flag_run)
		end if
		return
    end subroutine deallocate_scanner_m

end module scanner_m
