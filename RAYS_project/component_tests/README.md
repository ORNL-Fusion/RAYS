# Tests for components of RAYS_lib_

Wrapper code for tests of individual components of the RAYS\_lib static library.  It
links to RAYS\_lib and so can test components in the environment in which they will be
used.  This collection is distinct from the directory "test" which contains stand alone 
tests that don't link to RAYS\lib of simpler code pieces.

As of 11-22-2022 the only test here is is for damp\_fund\_ECH.f90

