 subroutine finalize_post_process_RAYS
!   Finish up, do post processing if any

    use diagnostics_m, only : message_unit, message, text_message, ray_list_unit, output_unit

    implicit none

    close(ray_list_unit)
    close(output_unit)
        
    call message(1)
    call text_message('Post processing RAYS finished', 1)
    call message(1)
    close(message_unit)

! Copy messages file to log.RAYS so it won't get clobbered by post processing
    call system('mv messages log.post_process_RAYS') 
    
    return
 end subroutine finalize_post_process_RAYS