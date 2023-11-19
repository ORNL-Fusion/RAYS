 subroutine finalize_run
!   Finish up, do post processing if any

!   External procedures: cpu_time (intrinsic)

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message, run_label, &
           & ray_list_unit, output_unit, t_start_RAYS, t_finish_RAYS, day_to_seconds, &
           & date_to_julian
    use ray_results_m, only : write_results_list_directed,&
                            & write_results_LD

    implicit none

! Time and date vector - local, not the one loaded in subroutine initialize()
    integer :: date_v(8), ierr
    real(KIND=rkind) :: code_time

    if (write_results_list_directed .eqv. .true.) then
        call write_results_LD
    end if

    close(ray_list_unit)
    close(output_unit)

!   Find date and time after writing results files
    call date_and_time (values=date_v)
    !   Convert end date_v to Julian
        call date_to_julian(date_v,t_finish_RAYS,ierr)
        if (ierr .ne. 0) then
            write(*,*) 'julian finish, ierr = ', ierr
            stop
        end if

    code_time = (t_finish_RAYS - t_start_RAYS)*day_to_seconds
    call message('Wall time including file writes', code_time, 0)

    call message(1)
    call text_message('RAYS run finished', 0)
    call message(1)
    close(message_unit)

! Copy messages file to log.RAYS so it won't get clobbered by post processing
    call system('mv messages log.RAYS.'//trim(run_label))

 end subroutine finalize_run