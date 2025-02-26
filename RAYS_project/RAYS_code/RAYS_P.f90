 program rays_p

!	This version links to the post_processing_lib library and (as of 3/23) uses
!	the deposition_profiles_m module to calculate power deposition profiles after
!	ray tracing completes.
!
!   External procedures: initialize (initialize.f90), trace_rays (trace_rays.f90)
!                        finalize (finalize.f90), all in RAYS_lib

    use post_processing_m, only : initialize_post_processing_m, post_process_m, &
                                & finalize_post_processing_m

    implicit none
    logical :: read_input = .true.

!****************************************************************************************
!   Ray tracing
!****************************************************************************************

!   Read input file and initialize variables.
    call initialize(read_input)

!   Trace the rays.
    call trace_rays

!   Finish up, do post processing if any
    call finalize_run

!****************************************************************************************
!   Post Processing
!****************************************************************************************

!   Read input file and initialize variables for post_processing_m only.
    read_input = .false.
    call initialize_post_processing_m(read_input)

!   Do the post processing. Contained in module post_processing_m.f90
    call post_process_m

!   Finish up, do post processing if any
    call finalize_post_processing_m


 end program rays_p
