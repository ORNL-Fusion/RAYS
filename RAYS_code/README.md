# RAYS
Plasma ray tracing code

This should build with gfortran.  At least it does on my Mac.

This is a very stripped down and simplified version of the RAYS code.  The toroidal plasma
configurations have been removed leaving just simple, perpendicularly stratified slab 
geometry.  Warm plasma routines have been removed leaving cold plasma dispersion, although
multiple species and arbitrary frequency are included.  All postprocessing and graphics
have been pulled out.  Rudimentary postprocessing and graphics routines specialized for
slab geometry are in the adjacent directories.  I expect to start adding back the more
realistic components but for now the intent is to produce a maximally simple version of
the code.

