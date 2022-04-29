 program post_process_RAYS

    use post_processing_m, only : post_process

!   Read input file and initialize variables.
    call initialize_post_process_RAYS

!   Do the post processing. Contained in module post_processing_m.f90
    call post_process

!   Finish up, do post processing if any
    call finalize_post_process_RAYS

 end program post_process_RAYS
