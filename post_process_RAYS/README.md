# Post Processing RAYS Output

Source codes for postprocessing output of the RAYS code.  At present the code doesn't
actually process the rays but just generates some files that scan equilibrium and kx(x)
profiles across the x range for numerical comparison with the rays.  It also generates a
file, *graphics_description_slab.dat* that is used by the graphics code plot_RAYS_slab.py

This should build with gfortran.  At least it does on my Mac. In addition to the Makefile
there is a script "cp\_from\_RAYS\_code\_and\_make" which copies all the fortran source files from
/RAYS_code/ next door.  The RAYS code source "RAYS.f90" confuses the makefile with an extra
main program so the script removes that before doing "make"

See the examples in */examples_RAYS/ECH_90GHz_slab*


