#!/bin/bash
rm messages
rm log.*
rm ray_out.*
rm ray_list.*
rm run_results.*
rm scan_summary.*
rm fort.*
rm *.nc
rm *.pdf

# N.B. this script deletes all .nc files inclding needed file
# Brz_fields.MPEX_9_filaments_D3-6_ECH_2nd_harm.nc.  So copy saved file
# Brz_fields.MPEX_9_filaments_D3-6_ECH_2nd_harm.nc_save back
cp Brz_fields.MPEX_9_filaments_D3-6_ECH_2nd_harm.nc_save \
Brz_fields.MPEX_9_filaments_D3-6_ECH_2nd_harm.nc
