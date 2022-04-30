# Post Processing RAYS Output

Source codes for postprocessing output of the RAYS code.  At present the code doesn't
actually process the rays since presently there are no heating or current drive profiles 
generated. It just generates some files that scan equilibrium and kx(x)
profiles across the x range for numerical comparison with the rays.  It also generates
files, *graphics_description_slab.dat* and *graphics_description_slab.dat* that are used 
by the graphics code plot\_RAYS\_slab.py to show ray trajectories.

See the examples in */examples_RAYS/ECH_90GHz_slab*

This code links to the RAYS\_lib static library.