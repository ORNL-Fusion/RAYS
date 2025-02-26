 program post_process_RAYS

    use post_processing_m, only : initialize_post_processing_m, post_process_m,&
                                & finalize_post_processing_m

    implicit none
    logical :: read_input = .true.

!   Read input file and initialize variables.
    call initialize_post_processing_m(read_input)

!   Do the post processing. Contained in module post_processing_m.f90
    call post_process_m

!   Finish up, do post processing if any
    call finalize_post_processing_m

 end program post_process_RAYS
