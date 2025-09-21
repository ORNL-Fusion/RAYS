# Post Processing RAYS Output

Codes for postprocessing output of the RAYS code.  At present (9/20/2025)
it supports four different plasma geometries: axisymmetric toroid, solovev (obsolete),
multiple axisymmetric mirrors, and slab. Exactly what processing is done is dependent
on the specific geometry, see the individual processor modules for details.  The module
*post_processing_m* is a wrapper for the geometry specific processors.  It reads ray data
from a netCDF file, *run_results.<run ID>.nc* and also reads the RAYS input namelist file
from the original code run. The various processors write a graphics description file
which, along with the run results file, are used by the separate graphics codes to plot
ray trajectories etc.  The geometry specific processors optionally write netCDF files with
data for plotting ray detailed diagnostics, various deposition profiles (such as power),
and plasma equilibrium profiles. See the examples in */examples_RAYS/*

This library also links to the RAYS\_lib and math\_functions static libraries.