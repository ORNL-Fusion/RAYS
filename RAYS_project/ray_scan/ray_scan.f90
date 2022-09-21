 program ray_scan
! A code that links to RAYS_lib and performs multiple ray tracing runs
! as determined by the module ray_scanner_m.  
!
! The scan loop is as follows:
! The scanner updates the input data in RAYS_lib modules at the beginning of
! each iteration of the scan loop
! The "initialize" subroutine is then called with read_input = .false.
! Then rays are traced as in a stand-alone run

!   External procedures: initialize (initialize.f90), trace_rays (trace_rays.f90)
!                        finalize (finalize.f90), all in RAYS_lib

    use diagnostics_m, only : message_unit, date_v
    use scanner_m, only : initialize_scanner_m, update_scan_parameter, aggregate_run_data,&
                        & write_scan_summary, deallocate_scanner_m, n_runs, scan_id,&
                        & scan_date_v
    
    implicit none 
    logical :: read_input = .true.
    integer :: i_run
    
!   Read input file and initialize variables before starting the scan loop
    call initialize(read_input)
    call initialize_scanner_m(read_input)
    close(message_unit)
    scan_date_v = date_v ! Get scan date from first initialization
    read_input = .false.
 
! Copy messages file to log.scan_initialization so it won't get clobbered
    call system('mv messages log.scan_init.'//trim(scan_id)) 
    
    run_loop: do i_run = 1, n_runs
    
		call update_scan_parameter(i_run)

		! Initialize with updated input parameters    
		call initialize(read_input)

		! Trace the rays.
		call trace_rays

	    ! Finish up
		call finalize_run
		
		! Add run data to scan data
		call aggregate_run_data(i_run)
    
    end do run_loop
    
    call write_scan_summary
    
    call deallocate
    call deallocate_scanner_m

!********************************************************************

 end program ray_scan
