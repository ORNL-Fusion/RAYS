 subroutine finalize
!   Finish up, do post processing if any

    use diagnostics_m, only : message_unit, message, text_message, t_start_rays, &
           & t_finish_rays, t_start_tracing, t_finish_tracing

    implicit none

!     close(94)
!     close(95)

    call cpu_time(t_finish_rays)
    call message('CPU time ray tracing', t_finish_tracing - t_start_tracing)
    call message('CPU time RAYS code', t_finish_rays - t_start_rays)


    write(*,*) ' '
    write(*,*) 'CPU time ray tracing = ', t_finish_tracing - t_start_tracing
    write(*,*) 'CPU time RAYS code', t_finish_rays - t_start_rays
    write(*,*) ' '
    
    call message()
    call text_message('RAYS finished')
    call message()
    close(message_unit)

! Copy messages file to log.RAYS so it won't get clobbered by post processing
    call system('mv messages log.RAYS') 

    write(*,*) ' '
    write(*,*) 'RAYS finished'
    write(*,*) ' '

 end subroutine finalize