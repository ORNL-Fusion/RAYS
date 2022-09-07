 subroutine finalize
!   Finish up, do post processing if any

!   External procedures: cpu_time (intrinsic)

    use diagnostics_m, only : message_unit, message, text_message, run_label,&
           & t_start_rays, t_finish_rays, t_start_tracing, t_finish_tracing
    use constants_m, only : rkind
    use ray_results_m, only : write_results_list_directed, write_results_LD, ray_trace_time

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

    write(*,*) ' '
    write(*,*) 'RAYS finished'
    write(*,*) ' '

 end subroutine finalize