# mirror_magnetics\_lib

A static library containing routines to generate magnetic fields and gradients for
circular current loops and multiple current loops aligned on the Z axis as for
a multiple mirror device.  It is particularly tailored for the MPEX device at ORNL,
but should be applicable to any axisymmetric mirror concept.

The coding consists of one executable program *mirror_magnetics.f90* and two modules
*B_loop_m.f90* and *mirror_magnetics_m.f90*.  The program is essentially just a wrapper
for the routines in module *mirror_magnetics_m*.  The output of the program is
a netCDF file containing *Br, Bz* and *Aphi* evaluated on a uniform *r,z* grid.
For use in the RAYS codes the netCDF file is read by the
*mirror_magnetics_spline_interp_m* module in RAYS_lib, which produces a 2D spline
interpolation of the RZ data and outputs (*Bx,By,Bz, Aphi*) and their gradients as
a function of position (*x,y,z*).

Module B_loop_m exports several routines to calculate fields for a circular current
loop of radius a, carrying current of 1 Amp, centered at *(x,y,z) = 0*, lying in the
*z = 0* plane.  The fields are calculated using the complete elliptic integrals *E(k)*
and *K(k)*, with a series expansion near *r = 0*, where there is cancellation in the
elliptic integral form.

Module *mirror_magnetics_m* uses routines in *B_loop_m* to produce fields from multiple
current loops.

Input data comes from three namelists, in three separate files:<br>
"mirror_magnetics_list" contains general data for the particular case to be generated<br>
"coil_data_list" describes the fixed coil configuration - coil locations and radii<br>
"current_data_list" lists the current in each coil



