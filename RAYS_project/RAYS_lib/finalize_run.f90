 subroutine finalize_run
!   Finish up, do post processing if any

!   External procedures: cpu_time (intrinsic)

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message, run_label, &
           & ray_list_unit, output_unit, &
           & t_start_rays, t_finish_rays, t_start_tracing, t_finish_tracing
    use ray_results_m, only : write_results_list_directed,&
                            & write_results_LD, run_trace_time

    implicit none
    
    real(KIND=rkind) :: trace_time, code_time

    call cpu_time(t_finish_rays)
    run_trace_time = t_finish_tracing - t_start_tracing 
    code_time = t_finish_rays - t_start_rays
    call message(1)
    call message('CPU time ray tracing', run_trace_time, 0)
    call message('CPU time RAYS code', code_time, 0)
    
    
    if (write_results_list_directed .eqv. .true.) then
        call write_results_LD
    end if
    
    close(ray_list_unit)
    close(output_unit)

    call message(1)
    call text_message('RAYS run finished', 0)
    call message(1)
    close(message_unit)

! Copy messages file to log.RAYS so it won't get clobbered by post processing
    call system('mv messages log.RAYS.'//trim(run_label)) 

 end subroutine finalize_run