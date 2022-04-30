# RAYS_project
The RAYS code and post\_process\_RAYS code in the neighboring directories have been
converted to a CMake project and CMakeLists.txt files have been written to build both codes.
More significantly the RAYS modules and subroutines are built as a static library, RAYS\_lib,
which is used by both the stand-alone RAYS code and post\_process\_RAYS.  In this way any 
code can incorporate the RAYS library without invoking the RAYS code *per se*.

The source files in the previous directories have been copied to the RAYS_project directory.
I expect to work from the project directory and that the older directories will atrophy 
and eventually disappear.  The graphics\_RAYS and examples\_RAYS directories are still active.
