 program rays
!   simple_RAYS greatly simplified version of the MORAYS code to serve
!   as test bed for ray tracing code development.  (DBB 5/2021)
!
!   MORAYS:  the Modified RAYS code.
!   The MORAYS code was written by Caiyi Wang in collaboration with 
!   D. Batchelor at Fusion Energy Division, Oak Ridge National Lab
!   between 6/95 - 9/95.  MORAYS was based on the RAYS code (1982).

    implicit none 

!   Read input file and initialize variables.
    call initialize

!   Trace the rays.
    call trace_rays

!   Finish up, do post processing if any
    call finalize

 end program rays
