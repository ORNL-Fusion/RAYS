 program rays
!   simple_RAYS greatly simplified version of the MORAYS code to serve
!   as test bed for ray tracing code development.  (DBB 5/2021)
!
!   External procedures: initialize (initialize.f90), trace_rays (trace_rays.f90)
!                        finalize (finalize.f90), all in RAYS_lib
    implicit none
    logical :: read_input = .true.

!   Read input file and initialize variables.
    call initialize(read_input)

!   Trace the rays.
    call trace_rays

!   Finish up, do post processing if any
    call finalize_run

 end program rays
