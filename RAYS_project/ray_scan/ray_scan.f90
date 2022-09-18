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

    use scanner_m, only : initialize_scanner_m, update_scan_parameter, n_iterations
    
    implicit none 
    logical :: read_input = .true.
    integer :: i_iter
    
!   Read input file and initialize variables before starting the scan loop
    call initialize(read_input)
    call initialize_scanner_m
    read_input = .false.
    
    iteration_loop: do i_iter = 1, n_iterations
    
		call update_scan_parameter(i_iter)

		! Initialize with updated input parameters    
		call initialize(read_input)

		!   Trace the rays.
		call trace_rays

	!   Finish up, do post processing if any
		call finalize
    
    end do iteration_loop

 end program ray_scan
