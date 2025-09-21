# splines\_lib

A static library containing modules and subroutines for simple 1D and 2D
cubic spline fits.  These routines were extracted from the PPPL EZ Spline library and
adapted to follow conventions of the RAYS system. The PPPL spline library has much more
generality than is kept here, but is more problematic to build, and doing this eliminates
need for external dependencies.

The spline\_lib itself does not link to other libraries, but the test executable
 *test_pspline* does link to RAYS\_lib and math\functions.