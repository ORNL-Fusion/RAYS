#! /usr/bin/env python

"""
eRAYS -> Does everything to run a RAYS case for axisymmetric totoidal geometry
Runs RAYS
Runs post_process_RAYS
Plots ray trajectories using bash function plot_rays
Plots depositon profiles using bash function plot_profs
Plots ray detailed diagnostics using bash function plot_diags
Plots equilibrium radial profiles using bash function plot_xycdf

Takes one optional command-line argument -> RAYS input namelist file. If present it
launches RAYS with that file as command line arg.  RAYS immediately copies that to its
default input file name "rays.in"

change log:
2/27/2025
 version 1.0

"""

import sys
import os
import subprocess
from simple_file_editing_functions import get_lines, read_string_var_from_nml_lines
debug = False

#----------------------------------------------------------------------------------------------
# Utility functions
#----------------------------------------------------------------------------------------------

# None so far

#----------------------------------------------------------------------------------------------
# Set up
#----------------------------------------------------------------------------------------------

# Get path to RAYS
rays_root = os.getenv('RAYS_ROOT')
rays_proj = os.getenv('RAYS_PROJECT')

# get optional RAYS input file file name -> default = rays.in
infile_name = 'rays.in'
if len(sys.argv) == 2:
    infile_name = sys.argv[1]
if len(sys.argv) > 2:
    print(' sys.argv = ', sys.argv)
    message = 'Usage: this script takes one optional argument -> RAYS input file name'
    print(message)
    raise Exception(message)

# Launch RAYS
path = os.path.join(rays_proj, 'bin/RAYS')
cmd = [path, infile_name]
print('Executing = ', cmd)
retcode = subprocess.call(cmd)
#retcode = subprocess.run([path, infile_name])
if (retcode != 0):
    logMsg = 'Error executing '.join(map(str, cmd))
    print('cmd = ', cmd)
    print('retcode = ', retcode)
    print('logMsg =  ', logMsg)
    raise Exception(logMsg)

# Launch post_process_RAYS
path = os.path.join(rays_proj, 'bin/post_process_RAYS')
cmd = [path, infile_name]
print('Executing = ', cmd)
retcode = subprocess.call(cmd)
if (retcode != 0):
    logMsg = 'Error executing '.join(map(str, cmd))
    print('cmd = ', cmd)
    print('retcode = ', retcode)
    print('logMsg =  ', logMsg)
    raise Exception(logMsg)

# Look in 'rays.in' and get run_label
input_file = 'rays.in'
lines = get_lines(input_file)

var = 'run_label'
run_label = read_string_var_from_nml_lines(lines, var, separator = ',')
print(var, ' = ', run_label)

var = 'equilib_model'
equilib_model = read_string_var_from_nml_lines(lines, var, separator = ',')
print(var, ' = ', equilib_model)

# Plot ray paths
if (equilib_model.strip("'") == 'axisym_toroid'):
    path = os.path.join(rays_root, 'graphics_RAYS/plot_RAYS_axisym_toroid.py')
    infile_name = 'run_results.' + run_label.strip("'") + '.nc'
    cmd = ['python', path, infile_name]
    print('Executing = ', cmd)
    retcode = subprocess.call(cmd)
    if (retcode != 0):
        logMsg = 'Error executing '.join(map(str, cmd))
        print('cmd = ', cmd)
        print('retcode = ', retcode)
        print('logMsg =  ', logMsg)
        raise Exception(logMsg)

# Plot ray diagnostics
path = os.path.join(rays_root, 'graphics_RAYS/plot_ray_diags.py')
infile_name = 'ray_detailed_diagnostics.' + run_label.strip("'") + '.nc'
if os.access(infile_name, os.R_OK):
    cmd = ['python', path, infile_name]
    print('Executing = ', cmd)
    retcode = subprocess.call(cmd)
    if (retcode != 0):
        logMsg = 'Error executing '.join(map(str, cmd))
        print('cmd = ', cmd)
        print('retcode = ', retcode)
        print('logMsg =  ', logMsg)
        raise Exception(logMsg)
else:
    print('File ', infile_name, ' is not readable')

# Plot deposition profiles
path = os.path.join(rays_root, 'graphics_RAYS/plot_profiles.py')
infile_name = 'deposition_profiles.' + run_label.strip("'") + '.nc'
if os.access(infile_name, os.R_OK):
    cmd = ['python', path, infile_name]
    print('Executing = ', cmd)
    retcode = subprocess.call(cmd)
    if (retcode != 0):
        logMsg = 'Error executing '.join(map(str, cmd))
        print('cmd = ', cmd)
        print('retcode = ', retcode)
        print('logMsg =  ', logMsg)
        raise Exception(logMsg)
else:
    print('File ', infile_name, ' is not readable')

# Plot equilibrium radial profiles for axisym_toroid
path = os.path.join(rays_root, 'graphics_RAYS/plot_XY_curves_netCDF.py')
infile_name = 'eq_radial_profiles.' + run_label.strip("'") + '.nc'
if os.access(infile_name, os.R_OK) and (equilib_model.strip("'") == 'axisym_toroid'):
    cmd = ['python', path, infile_name]
    print('Executing = ', cmd)
    retcode = subprocess.call(cmd)
    if (retcode != 0):
        logMsg = 'Error executing '.join(map(str, cmd))
        print('cmd = ', cmd)
        print('retcode = ', retcode)
        print('logMsg =  ', logMsg)
        raise Exception(logMsg)
else:
    print('File ', infile_name, ' is not readable')

print('Finsihed everything')