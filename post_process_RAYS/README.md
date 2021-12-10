# Post Processing RAYS Output

Source codes for postprocessing output of the RAYS  code

This should build with gfortran.  At least it does on my Mac. In addition to the Makefile
there is a script "cp\_from\_RAYS\_code\_and\_make" which copies all the fortran source files from
/RAYS_code/ next door.  The RAYS code source "RAYS.f90" confuses the makefile with an extra
main program so the script removes that before doing "make"

Real soon now I will generate some example input files to run this code 
code to do rudimentary plotting of ECH case in slab geometry.



