 subroutine finalize_run
!   Finish up, do post processing if any

!   External procedures: cpu_time (intrinsic)

    use constants_m, only : rkind, ray_list_unit, output_unit
    use diagnostics_m, only : message_unit, message, text_message, run_label,&
           & t_start_rays, t_finish_rays, t_start_tracing, t_finish_tracing
    use ray_results_m, only : write_results_list_directed,&
                            & write_results_LD, run_trace_time

    implicit none
    
    real(KIND=rkind) :: trace_time, code_time

    call cpu_time(t_finish_rays)
    run_trace_time = t_finish_tracing - t_start_tracing 
    code_time = t_finish_rays - t_start_rays
    call message('CPU time ray tracing', run_trace_time)
    call message('CPU time RAYS code', code_time)


    write(*,*) ' '
    write(*,*) 'CPU time ray tracing = ', run_trace_time
    write(*,*) 'CPU time RAYS code', code_time
    write(*,*) ' '
    
    
    if (write_results_list_directed .eqv. .true.) then
        call write_results_LD
    end if
    
    close(ray_list_unit)
    close(output_unit)

    call message()
    call text_message('RAYS run finished')
    call message()
    close(message_unit)

! Copy messages file to log.RAYS so it won't get clobbered by post processing
    call system('mv messages log.RAYS.'//trim(run_label)) 

    write(*,*) ' '
    write(*,*) 'RAYS run finished'
    write(*,*) ' '

 end subroutine finalize_run