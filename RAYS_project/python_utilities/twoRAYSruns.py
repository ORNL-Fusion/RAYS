#! /usr/bin/env python

"""
twoRAYSruns -> Does everything to run multiple RAYS cases, e.g. for O-X conversion
where you do an O-mode run, then run the post processor for that case, then do an
X-mode run as a continuation, then run the post processor on that case, then plot the
results.

Works for for axisym_toroid, multiple_mirror or slab geometry

For each case:
Runs RAYS
Runs post_process_RAYS
Plots ray trajectories of first run using bash function plot_rays
Plots depositon profiles using bash function plot_profs
Plots ray detailed diagnostics using bash function plot_diags
Plots equilibrium radial profiles using bash function plot_xycdf

Takes command line arguments -> RAYS input namelist files for each case.

N.B.  The namelist variable 'run_label' in the two RAYS input files must be different, or
      the output from the second run will overwrite the files from the first run.

N.B.  To run the post processor it needs post process input files for each RAYS run.
      The post processor routine "initialize_post_processing_m" will look for file
      names of the form "post_process_<run_label>.in.  The user has to provide those
      files.

change log:
2/27/2026
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

# get RAYS input file file names
infile_names = []
pp_files = []

if len(sys.argv) > 1:
    infile_names = sys.argv[1]
    print('infiles =', infile_names)
if sys.argv = 1:
    print(' sys.argv = ', sys.argv)
    message = 'Usage: this script takes arguments -> RAYS input file names'
    print(message)
    raise Exception(message)

for infile_name in infile_names:
    print('\nrun ', infile_name)
	# Launch RAYS
	path = os.path.join(rays_proj, 'bin/RAYS')
	cmd = [path, infile_name]
	print('Executing = ', cmd)
	retcode = subprocess.call(cmd)
	if (retcode != 0):
		logMsg = 'Error executing '.join(map(str, cmd))
		print('cmd = ', cmd)
		print('retcode = ', retcode)
		print('logMsg =  ', logMsg)
		raise Exception(logMsg)

	# Look in 'rays.in' and get run_label etc.
	input_file = 'rays.in'
	lines = get_lines(input_file)

	var = 'run_label'
	run_label = read_string_var_from_nml_lines(lines, var, separator = ',')
	print(var, ' = ', run_label)

	var = 'equilib_model'
	equilib_model = read_string_var_from_nml_lines(lines, var, separator = ',')
	print(var, ' = ', equilib_model)

	# Launch post_process_RAYS
	path = os.path.join(rays_proj, 'bin/post_process_RAYS')
	pp_file = 'post_process_' + run_label.strip("'") + '.in'
	pp_files.append(pp_file)
	cmd = [path, pp_file]
	print('Executing = ', cmd)
	retcode = subprocess.call(cmd)
	if (retcode != 0):
		logMsg = 'Error executing '.join(map(str, cmd))
		print('cmd = ', cmd)
		print('retcode = ', retcode)
		print('logMsg =  ', logMsg)
		raise Exception(logMsg)

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

	if (equilib_model.strip("'") == 'multiple_mirror'):
		path = os.path.join(rays_root, 'graphics_RAYS/plot_RAYS_mirror.py')
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

	if (equilib_model.strip("'") == 'slab'):
		path = os.path.join(rays_root, 'graphics_RAYS/plot_RAYS_slab.py')
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


	# Plot equilibrium radial profiles for axisym_toroid or multiple mirror
	if (equilib_model.strip("'") == 'axisymmetric_toroid' or\
		equilib_model.strip("'") == 'multiple_mirror'):
		path = os.path.join(rays_root, 'graphics_RAYS/plot_XY_curves_netCDF.py')
		infile_name = 'eq_radial_profiles.' + run_label.strip("'") + '.nc'
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

	if (equilib_model.strip("'") == 'slab'):

		# Plot equilibrium X profiles for slab
		path = os.path.join(rays_root, 'graphics_RAYS/plot_XY_curves_netCDF.py')
		infile_name = 'eq_X_profiles.' + run_label.strip("'") + '.nc'
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
# If multiple input files plot all rays togehter

If len(infile_names) > 1:
	print('pp_files = ', pp_files)

print('Finsihed everything')