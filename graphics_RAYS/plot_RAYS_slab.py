#! /usr/bin/env python

"""
plot_RAYS.py -> Plots ray trajectories and k vectors from data in file ray.out
DBB 11/19/2021

"""
# Working notes:
#
# DBB (11/23/2021)
# Changed input file from generic ray.out to run-specific ray_out.<run_label>.
# <run_label> is specified in graphics_description_slab file.
# I should add an optional command line arg to specify a choice of <run_lebel>.  That way
# can have multiple run label files in the working directory.
#

import sys
import os
import math
from netCDF4 import *
import matplotlib.pyplot as plt

from simple_file_editing_functions import get_lines, input_file_to_variable_dict
from plt_XY_Curves import *

debug = 0
k_vec_base_length = 0.001

#----------------------------------------------------------------------------------------------
# Utility functions
#----------------------------------------------------------------------------------------------
# 
def  n_evenly_spaced_integers(n, Length):

# The purpose is to generate a list of integers that can be used to index a subset of a long 
# list of length = Length, at an approximately even stride, but including the first and last
# items in the list.
# It tries to fit n evenly spaced integers in the range 0:Length-1, provided n >= Length
# If n == 0 it returns [0] i.e. it indexes the first item in the list
# If n == -1 it returns [-1] i.e. it indexes the last item in the list
# If n >= Length you are asking for more indices than would be in the list. It returns range(Length)
# If n < Length-1 returns [0, n-2 integers spaced as evenly as possible, Length-1].  Of
# course the returned list won't be exactly evenly spaced unless n divides Length.

    if (type(n) != int) or (type(Length) != int) or (Length < 1):
        print('error n_evenly_spaced_integers: arguments must be positive integers',\
              ' n = ', n, 'Length = ', Length)
              
    if n >= Length:
        return list(range(Length))
    elif n == 0:
        return [0]
    elif n == -1:
        return [-1]
    else:
        return [int((i)*float(Length-1)/(n-1)) for i in range(n)]        
#
#----------------------------------------------------------------------------------------------
# Main program
#----------------------------------------------------------------------------------------------

# Get data from graphics description input file
graphics_variable_dict = input_file_to_variable_dict('graphics_description_slab.dat')
if debug > 1: print('graphics_variable_dict = ', graphics_variable_dict)

run_description = graphics_variable_dict['run_description'] 
run_label = graphics_variable_dict['run_label'] 

# get the command line
ray_file_list = []
n_arg = len(sys.argv)
if n_arg == 1: # No arg, get run_label from graphics description file
    n_ray_files = 1
    ray_file_list.append('ray_out.' + run_label)

if n_arg > 1: # Get ray file names from command line
    n_ray_files = n_arg-1
    ray_file_list = sys.argv[1:]

print ('ray files = ', ray_file_list)

nray = 0
rays_s_list = []
rays_x_list = []
rays_y_list = []
rays_z_list = []
rays_kx_list = []
rays_ky_list = []
rays_kz_list = []
rays_knorm_list = []
rays_npoints_list =[]

#----------------------------------------------------------------------------------------------
# Cycle through ray file list
#----------------------------------------------------------------------------------------------

for file in ray_file_list:

    print('file = ', file)
    # Ray data input file
    lines = get_lines(file) 

    lines.append('0.0\n')  # Tack on an extra zero line to signal the end of file


    k_max = 0.

    for i in range(len(lines)-1):
        line = lines[i]
        split_line = line.split()
        if debug > 3: print('split_line = ', split_line)    
        num_line = [float(x) for x in split_line]
        if debug > 2: print('num_line = ', num_line)
    
        next_line = lines[i+1]
        split_next_line = next_line.split()
        num_next_line = [float(x) for x in split_next_line]

        if num_line[0] == 0:  # This is a new ray

            nray = nray+1   # Increment ray counter
            s_list = []     # Reinitialize for new ray
            x_list = []
            y_list = []
            z_list = []
            kx_list = []
            ky_list = []
            kz_list = []        
            knorm_list = []        
            npoints = 0
    
        s_list.append(num_line[0])
        x_list.append(num_line[1])
        y_list.append(num_line[2])
        z_list.append(num_line[3])
        kx_list.append(num_line[4])
        ky_list.append(num_line[5])
        kz_list.append(num_line[6])
    
        kmag = math.sqrt(pow(num_line[4],2) + pow(num_line[5],2) +pow(num_line[6],2))
        knorm_list.append(kmag)
        if kmag > k_max: k_max = kmag

        npoints = npoints +1        # Increment points counter
        
        # Test to see if this is last line of a ray
        if num_next_line[0] == 0:  # This is last line of this ray. Add lists to rays_lists
           rays_s_list.append(s_list)
           rays_x_list.append(x_list)
           rays_y_list.append(y_list)
           rays_z_list.append(z_list)
           rays_kx_list.append(kx_list)
           rays_ky_list.append(ky_list)
           rays_kz_list.append(kz_list)
           rays_knorm_list.append(knorm_list)
           rays_npoints_list.append(npoints)
 
    if debug > 1: print('\n', 'rays_x_list[0] = ', rays_x_list[0],'\n')
    if debug > 1: print('rays_z_list[1] = ', rays_z_list[0])
    if debug > 0: print('\n', 'rays_kx_list[0] = ', rays_kx_list[0],'\n')
    if debug > 0: print('rays_kz_list[0] = ', rays_kz_list[0])

    print('nray = ', nray)
    print('k_max = ', k_max)
    print ('len(rays_s_list) = ', len(rays_s_list))

#----------------------------------------------------------------------------------------------
# Generate Z-X ray plot
#----------------------------------------------------------------------------------------------


xmin = float(graphics_variable_dict['xmin'])
xmax = float(graphics_variable_dict['xmax'])
ymin = float(graphics_variable_dict['ymin'])
ymax = float(graphics_variable_dict['ymax'])
zmin = float(graphics_variable_dict['zmin'])
zmax = float(graphics_variable_dict['zmax'])

num_plot_k_vectors = int(graphics_variable_dict['num_plot_k_vectors'])
scale_k_vec = graphics_variable_dict['scale_k_vec']
set_XY_lim = graphics_variable_dict['set_XY_lim']

print('run_description = ', run_description)
print('run_label = ', run_label)
print('xmin = ', xmin, ' xmax = ', xmax, ' ymin = ',  ymin, ' ymax = ', ymax,  'zmin =',\
        zmin, ' zmax = ', zmax ) 

print('num_plot_k_vectors = ', num_plot_k_vectors)
print('scale_k_vec = ', scale_k_vec)
print('set_XY_lim = ', set_XY_lim)


#----------------------------------------------------------------------------------------------
# Generate plot using calls to plot_XY_Curves.py
#----------------------------------------------------------------------------------------------
 

# Open graphics output file.  N.B. run_label comes from graphics description file.
# Edit that if you want to customize the run label for multiple ray files.
 
open_file_XY_Curves_Fig('ray_plots.' + run_label + '.pdf')
max_size = 8.
title = run_description + '  ' + run_label

xz_ratio = (xmax-xmin)/(zmax-zmin)
z_size = max_size
x_size = z_size*xz_ratio
figsize = (z_size, x_size)
diagonal = math.sqrt(pow(x_size,2) + pow(z_size,2))

xlabel = 'z(m)'
ylabel = 'x(m)'
curve_list = []

for i_ray in range(len(rays_s_list)):
    lbl = 'ray ' + str(i_ray + 1)
    new_curve = XY_curve(rays_z_list[i_ray], rays_x_list[i_ray], label = lbl)
    curve_list.append(new_curve)

if set_XY_lim in ['True', 'true', 'T']:
    plotZX = XY_Curves_Fig(curve_list, title, xlabel, ylabel, figsize=figsize, \
             ylim = [xmin,xmax], xlim = [zmin,zmax])
else:
    plotZX = XY_Curves_Fig(curve_list, title, xlabel, ylabel, figsize=figsize)

# Add k vectors at selected points
if num_plot_k_vectors > 0:
    for i_ray in range(len(rays_s_list)):
        indices = n_evenly_spaced_integers(num_plot_k_vectors, rays_npoints_list[i_ray])
        for i in indices:
            if scale_k_vec in ['True', 'true', 'T']:
                kx = k_vec_base_length*diagonal*rays_kx_list[i_ray][i]/k_max
                kz = k_vec_base_length*diagonal*rays_kz_list[i_ray][i]/k_max
            else:
                kx = k_vec_base_length*diagonal*rays_kx_list[i_ray][i]/rays_knorm_list[i_ray][i]
                kz = k_vec_base_length*diagonal*rays_kz_list[i_ray][i]/rays_knorm_list[i_ray][i]
        
            if debug > 1: print('kx = ', kx, ' kz = ', kz)
            plt.arrow(rays_z_list[i_ray][i], rays_x_list[i_ray][i],\
              kz, kx, shape='full', head_width = 0.01)


plot_XY_Curves_Fig(plotZX)
  
#----------------------------------------------------------------------------------------------
# Finalize
#----------------------------------------------------------------------------------------------
# 
# plot_index(index, 1)
#     

close_file_XY_Curves_Fig()

# 
# if debug:
#     print('index = ', index)

