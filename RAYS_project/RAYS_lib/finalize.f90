 subroutine finalize
!   Finish up, do post processing if any

!   External procedures: cpu_time (intrinsic)

    use diagnostics_m, only : message_unit, message, text_message, run_label,&
           & t_start_rays, t_finish_rays, t_start_tracing, t_finish_tracing, &
           & finalize_diagnostics_m
    use constants_m, only : finalize_constants_m, rkind
    use species_m, only : finalize_species_m
    use rf_m, only : finalize_rf_m
    use damping_m, only : finalize_damping_m
    use equilibrium_m, only : finalize_equilibrium_m
    use ray_init_m, only : finalize_ray_init_m
    use ode_m, only : finalize_ode_solver_m
    use ray_results_m, only : finalize_ray_results_m, write_results_list_directed,&
                            & write_results_LD, ray_trace_time

    implicit none
    
    real(KIND=rkind) :: trace_time, code_time

!     close(94)
!     close(95)

    call cpu_time(t_finish_rays)
    ray_trace_time = t_finish_tracing - t_start_tracing 
    code_time = t_finish_rays - t_start_rays
    call message('CPU time ray tracing', ray_trace_time)
    call message('CPU time RAYS code', code_time)


    write(*,*) ' '
    write(*,*) 'CPU time ray tracing = ', ray_trace_time
    write(*,*) 'CPU time RAYS code', code_time
    write(*,*) ' '
    
    call message()
    call text_message('RAYS finished')
    call message()
    close(message_unit)
    
    if (write_results_list_directed .eqv. .true.) then
        call write_results_LD(59)
    end if

! Copy messages file to log.RAYS so it won't get clobbered by post processing
    call system('mv messages log.RAYS_'//trim(run_label)) 

! finalize all the modules that initialize, allocate, or open files
    call finalize_constants_m
    call finalize_species_m
    call finalize_rf_m  
    call finalize_damping_m
    call finalize_equilibrium_m
    call finalize_ray_init_m  
    call finalize_ode_solver_m
    call finalize_ray_results_m
    
    call finalize_diagnostics_m ! This really is the end

    write(*,*) ' '
    write(*,*) 'RAYS finished'
    write(*,*) ' '

 end subroutine finalize