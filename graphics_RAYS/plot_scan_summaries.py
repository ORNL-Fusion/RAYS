#! /usr/bin/edim_v_vector python

"""
plot_scan_summaries.py -> Plots curves of various quantities accumulated from multiple
RAYS runs by ray_scan.f90.  The data is contained in files scan_summary.<scan_id>.
There can be multiple scans.  Paths to the summary files are specified in an input file, 
scan_paths.dat

"""
# Working notes:
#
# DBB (9/23/2022)
#

import sys
import os
import math
from netCDF4 import *
import matplotlib.pyplot as plt
import numpy as np

from simple_file_editing_functions import get_lines, list_variables_in_fortran_star_file,\
      read_var_from_fortran_star_file, lines_to_list
from plt_XY_Curves import *

debug = 2

#----------------------------------------------------------------------------------------------
# Utility functions
#----------------------------------------------------------------------------------------------

# Convert list to numpy array.  Split list into sublists of length len_sub_list
def split_flat_list(in_list, len_sub_list):

# check that len(in_list) is a multiple of len_sub_list
    if(float(int(len(in_list)/len_sub_list)) != len(in_list)/len_sub_list):
        message = 'split_flat_list: len(in_list not divoisible by len_sub_list)'
        print(message)
        raise

    d1 = int(len(in_list)/len_sub_list)
    x = np.array(in_list)
    return x.reshape(d1, len_sub_list)

#_________________________________________________________________________________________________
# Takes a flat list of numbers consisting of a concatenation of vectors each of length 
# len_sub_list.  Converts in_list to numpy array and splits into individual vector sub arrays.
# Takes a slice of each sub array with given offset and subtracts off a ref_array of the
# same length as the slice.  Calculates the 2-norm of the difference array.  Returns a
# numpy array containing the difference 2-norm for each vector in in_list
# 
def deviation_norm_sub_subarray(in_list, len_sub_list, ref_list, offset, norm = 1.):
    p_array = split_flat_list(in_list, len_sub_list)
    ref_array = np.array(ref_list)
    sub_list = []
    for sub_array in p_array:
        sub_sub_array = sub_array[offset : offset + len(ref_array)]
        diff = sub_sub_array - ref_array
        x = math.sqrt(np.sum(np.square(diff)))/norm
        sub_list.append(x)
    return np.array(sub_list)   

#----------------------------------------------------------------------------------------------
# Main program
#----------------------------------------------------------------------------------------------

# get the command line
scan_file_list = []
n_arg = len(sys.argv)
if n_arg == 1: # No arg, get paths from from scan_paths.dat
    # Get file paths from scan_paths.dat
    lines = get_lines('scan_paths.dat')
    scan_file_list = lines_to_list(lines)
    n_scan_files = len(scan_file_list)
    
    if debug > 1: print('scan_file_list = ', scan_file_list)

if n_arg > 1: # Get scan file names from command line
    n_scan_files = n_arg-1
    scan_file_list = sys.argv[1:]

# First  is path to reference run
ref_summary_file = scan_file_list[0]
print ('ref_summary_file = ', ref_summary_file)
scan_file_list = scan_file_list[1:]
print ('scan files = ', scan_file_list)

#----------------------------------------------------------------------------------------------
# Find out what variables are in files
#----------------------------------------------------------------------------------------------
lines = get_lines(scan_file_list[0]) 
scan_variables_list = list_variables_in_fortran_star_file(lines)
print('scan variables = ', scan_variables_list)

scan_id_list = []
scan_date_v_dict = {}
p_values_dict = {}
dim_v_vector_dict = {} # Should be same for all scans
trace_time_run_dict = {}
end_ray_param_run_dict = {}
end_resid_run_dict = {}
max_resid_run_dict = {}
ray_stop_flag_run_dict = {}
end_ray_vec_run_dict = {}
deviation_r_vec_run_dict = {}
deviation_k_vec_run_dict = {}

#----------------------------------------------------------------------------------------------
# Get reference r_vec,  k_vec and length of v(:) vector
#----------------------------------------------------------------------------------------------

lines = get_lines(ref_summary_file)

print(' ')
ref_start_ray_vec = read_var_from_fortran_star_file(lines, 'start_ray_vec', 'float')
ref_start_r_vec = ref_start_ray_vec[0:3]
ref_start_k_vec = ref_start_ray_vec[4:7]
dim_v_vector = read_var_from_fortran_star_file(lines, 'dim_v_vector', 'int')
print('ref_start_r_vec = ', ref_start_r_vec, '    ref_start_k_vec = ', ref_start_k_vec,\
      '    dim_v_vector = ', dim_v_vector)

print(' ')
ref_end_ray_vec = read_var_from_fortran_star_file(lines, 'end_ray_vec', 'float')
ref_end_r_vec = ref_end_ray_vec[0:3]
ref_end_k_vec = ref_end_ray_vec[4:7]
dim_v_vector = read_var_from_fortran_star_file(lines, 'dim_v_vector', 'int')
print('ref_end_r_vec = ', ref_end_r_vec, '    ref_end_k_vec = ', ref_end_k_vec,\
      '    dim_v_vector = ', dim_v_vector)

# calculate reference start k_vec norm
ref_start_k_vec_norm = math.sqrt(np.sum(np.square(ref_start_k_vec)))
print('ref_start_k_vec_norm = ', ref_start_k_vec_norm)

#----------------------------------------------------------------------------------------------
# Cycle through scan file list
#----------------------------------------------------------------------------------------------

for file in scan_file_list:

    print('\n','file = ', file)
    # scan data input file
    lines = get_lines(file) 

    scan_id = read_var_from_fortran_star_file(lines, 'scan_id', 'string')
    scan_id_list.append(scan_id)

    scan_parameter = read_var_from_fortran_star_file(lines, 'scan_parameter', 'string')

    scan_date_v = read_var_from_fortran_star_file(lines, 'scan_date_v', 'float')    
    scan_date_v_dict[scan_id] = scan_date_v

    p_values = read_var_from_fortran_star_file(lines, 'p_values', 'float')    
    p_values_dict[scan_id] = p_values 

    dim_v_vector = read_var_from_fortran_star_file(lines, 'dim_v_vector', 'int')
    dim_v_vector_dict[scan_id] = dim_v_vector     # Should be same for all scans
    
    trace_time_run = read_var_from_fortran_star_file(lines, 'trace_time_run', 'float')
    trace_time_run_dict[scan_id] = trace_time_run

    end_ray_param_run = read_var_from_fortran_star_file(lines, 'end_ray_param_run', 'float')
    end_ray_param_run_dict[scan_id] = end_ray_param_run

    end_resid_run = read_var_from_fortran_star_file(lines, 'end_resid_run', 'float')
    end_resid_run_dict[scan_id] = end_resid_run

    max_resid_run = read_var_from_fortran_star_file(lines, 'max_resid_run', 'float')
    max_resid_run_dict[scan_id] = max_resid_run
    
    ray_stop_flag_run = read_var_from_fortran_star_file(lines, 'ray_stop_flag_run', 'string')    
    ray_stop_flag_run_dict[scan_id] = ray_stop_flag_run

    end_ray_vec_run = read_var_from_fortran_star_file(lines, 'end_ray_vec_run', 'float')   
    end_ray_vec_run_dict[scan_id] = end_ray_vec_run
    
    # normalization = 1, deviation is in meters
    deviation_r_vec_run = deviation_norm_sub_subarray(end_ray_vec_run, dim_v_vector,\
                          ref_end_r_vec, 0, 1.)
    deviation_r_vec_run_dict[scan_id] = deviation_r_vec_run

    # normalization to magnitude of starting k_vec of reference run
    deviation_k_vec_run = deviation_norm_sub_subarray(end_ray_vec_run, dim_v_vector,\
                          ref_end_k_vec, 4, ref_start_k_vec_norm)
    deviation_k_vec_run_dict[scan_id] = deviation_k_vec_run

 
    if debug > 1: 
        print('\n','scan_id = ', scan_id,'\n')
        print('\n','scan_parameter = ', scan_parameter,'\n')
        print('scan_date_v = ', scan_date_v_dict[scan_id],'\n')
        print('dim_v_vector = ', dim_v_vector_dict[scan_id],'\n')
        print('trace_time_run = ', trace_time_run_dict[scan_id],'\n')
        print('end_ray_param_run = ', end_ray_param_run_dict[scan_id],'\n')
        print('end_resid_run = ', end_resid_run_dict[scan_id],'\n')
        print('max_resid_run = ', max_resid_run_dict[scan_id],'\n')
        print('ray_stop_flag_run = ', ray_stop_flag_run_dict[scan_id],'\n')
        print('end_ray_vec_run = ', end_ray_vec_run_dict[scan_id],'\n')
        print('deviation_r_vec_run [-1] = ', deviation_r_vec_run_dict[scan_id][-1],'\n')
        print('deviation_k_vec_run [-1] = ', deviation_k_vec_run_dict[scan_id][-1],'\n')

#----------------------------------------------------------------------------------------------
# Do plots 
#----------------------------------------------------------------------------------------------
 
open_file_XY_Curves_Fig('scan_summary_plots.pdf')

# Plot max_resid_run
curve_list = []
for id in scan_id_list:
    x = p_values_dict[id]
    lbl = id
    y = max_resid_run_dict[id]
    new_curve = XY_curve(x, y, label = lbl)
    curve_list.append(new_curve)

title = 'Maximum Run Residual vs ' + scan_parameter
xlabel = scan_parameter
ylabel = 'max residual'
plot1 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
plot_XY_Curves_Fig(plot1)

# Plot end_resid_run
curve_list = []
for id in scan_id_list:
    x = p_values_dict[id]
    lbl = id
    y = end_resid_run_dict[id]
    new_curve = XY_curve(x, y, label = lbl)
    curve_list.append(new_curve)

title = 'End Run Residual vs ' + scan_parameter
xlabel = scan_parameter
ylabel = 'end residual'
plot1 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
plot_XY_Curves_Fig(plot1)

# Plot trace_time_run
curve_list = []
for id in scan_id_list:
    x = p_values_dict[id]
    lbl = id
    y = trace_time_run_dict[id]
    new_curve = XY_curve(x, y, label = lbl)
    curve_list.append(new_curve)

title = 'Run Time vs ' + scan_parameter
xlabel = scan_parameter
ylabel = 'time (sec)'
plot1 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
plot_XY_Curves_Fig(plot1)

# Plot end_ray_param_run
curve_list = []
for id in scan_id_list:
    x = p_values_dict[id]
    lbl = id
    y = end_ray_param_run_dict[id]
    new_curve = XY_curve(x, y, label = lbl)
    curve_list.append(new_curve)

title = 'Ray Path Length (m), s, vs ' + scan_parameter
xlabel = scan_parameter
ylabel = 's (m)'
plot1 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
plot_XY_Curves_Fig(plot1)

# Plot end_r_vec deviation
curve_list = []
for id in scan_id_list:
    x = p_values_dict[id]
    lbl = id
    y = deviation_r_vec_run_dict[id]
    new_curve = XY_curve(x, y, label = lbl)
    curve_list.append(new_curve)

title = 'End |r_vec| Deviation (m) vs ' + scan_parameter
xlabel = scan_parameter
ylabel = 's (m)'
plot1 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
plot_XY_Curves_Fig(plot1)

# Plot end_k_vec deviation
curve_list = []
for id in scan_id_list:
    x = p_values_dict[id]
    lbl = id
    y = deviation_k_vec_run_dict[id]
    new_curve = XY_curve(x, y, label = lbl)
    curve_list.append(new_curve)

title = 'End |k_vec| Deviation from ref |k_vec| vs ' + scan_parameter
xlabel = scan_parameter
ylabel = 'Relative k end deviation'
plot1 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
plot_XY_Curves_Fig(plot1)

close_file_XY_Curves_Fig()
