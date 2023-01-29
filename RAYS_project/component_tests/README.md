# Tests for components of RAYS_lib

Wrapper code for tests of individual components of the RAYS\_lib static library.  It
links to RAYS\_lib and so can test components in the environment in which they will be
used.  This collection is distinct from the directory "test" which contains stand alone 
tests of simpler code pieces that don't link to RAYS\lib.

The wrapper code executable "component_tests" lives in RAYS\_project/bin.

As of 11-22-2022 the only test here is for damp\_fund\_ECH.f90

