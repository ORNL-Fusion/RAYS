#! /usr/bin/env python

"""
plot_RAYS_profiles.py -> Plots profiles from deposition_profiles.<run_label>
DBB 11/20/2023

"""
# Working notes:
#
# DBB (11/23/2021)
# Changed input file from generic ray.out to run-specific ray_out.<run_label>.
# <run_label> is specified in graphics_description_solovev file.
# I should add an optional command line arg to specify a choice of <run_lebel>.  That way
# can have multiple run label files in the working directory.
#

import sys
import os
import math
from netCDF4 import *
import matplotlib.pyplot as plt
import numpy as np

from simple_file_editing_functions import get_lines, input_file_to_variable_dict,\
      dict_variable_to_list_of_floats, read_var_from_fortran_star_file
from plt_XY_Curves import *

debug = 0


#
#----------------------------------------------------------------------------------------------
# Main program
#----------------------------------------------------------------------------------------------

# Get data from graphics description input file
graphics_variable_dict = input_file_to_variable_dict('graphics_description_axisym_toroid.dat')
if debug > 1: print('graphics_variable_dict = ', graphics_variable_dict)

run_description = graphics_variable_dict['run_description']
run_label = graphics_variable_dict['run_label']

# get the command line
n_arg = len(sys.argv)
if n_arg == 1: # No arg, get run_label from graphics description file
    profile_input_filename = 'deposition_profiles.' + run_label

if n_arg =2: # Get input file name from command line
    profile_input_filename = sys.argv[2]

print ('iplot_RAYS_profiles.py: input file = ', profile_input_filename)

nray = 0
rays_s_list = []
rays_x_list = []
rays_y_list = []
rays_z_list = []
rays_r_list = []
rays_kx_list = []
rays_ky_list = []
rays_kz_list = []
rays_kr_list = []
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
            r_list = []
            kx_list = []
            ky_list = []
            kz_list = []
            kr_list = []
            knorm_list = []
            npoints = 0

        s_list.append(num_line[0])
        x_list.append(num_line[1])
        y_list.append(num_line[2])
        z_list.append(num_line[3])

        rmag = math.sqrt(pow(num_line[1],2) + pow(num_line[2],2))
        r_list.append(rmag)

        kx_list.append(num_line[4])
        ky_list.append(num_line[5])
        kz_list.append(num_line[6])

        kr_list.append((num_line[1]*num_line[4] +num_line[3]*num_line[5])/rmag)

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
           rays_r_list.append(r_list)
           rays_kx_list.append(kx_list)
           rays_ky_list.append(ky_list)
           rays_kz_list.append(kz_list)
           rays_kr_list.append(kz_list)
           rays_knorm_list.append(knorm_list)
           rays_npoints_list.append(npoints)

    if debug > 1: print('\n', 'rays_x_list[0] = ', rays_x_list[0],'\n')
    if debug > 1: print('rays_z_list[1] = ', rays_z_list[0])
    if debug > 0: print('\n', 'rays_kx_list[0] = ', rays_kx_list[0],'\n')
    if debug > 0: print('rays_kz_list[0] = ', rays_kz_list[0])

    print('nray = ', nray)
    print('k_max = ', k_max)
    print ('len(rays_s_list) = ', len(rays_s_list))


# Open graphics output file.  N.B. run_label comes from graphics description file.
# Edit that if you want to customize the run label for multiple ray files.

open_file_XY_Curves_Fig('ray_plots.' + run_label + '.pdf')
max_size = 8.
title = run_description + '  ' + run_label

#----------------------------------------------------------------------------------------------
# Generate R-Z ray plot
#----------------------------------------------------------------------------------------------


xmin = float(graphics_variable_dict['box_rmin'])
xmax = float(graphics_variable_dict['box_rmax'])
zmin = float(graphics_variable_dict['box_zmin'])
zmax = float(graphics_variable_dict['box_zmax'])

num_plot_k_vectors = int(graphics_variable_dict['num_plot_k_vectors'])
scale_k_vec = graphics_variable_dict['scale_k_vec']
set_XY_lim = graphics_variable_dict['set_XY_lim']

print('run_description = ', run_description)
print('run_label = ', run_label)
print('xmin = ', xmin, ' xmax = ', xmax, 'zmin =', zmin, ' zmax = ', zmax )

print('num_plot_k_vectors = ', num_plot_k_vectors)
print('scale_k_vec = ', scale_k_vec)
print('set_XY_lim = ', set_XY_lim)

# Generate plot using calls to plot_XY_Curves.py

xz_ratio = (xmax-xmin)/(zmax-zmin)
z_size = max_size
x_size = z_size*xz_ratio
figsize = (z_size, x_size)
diagonal = math.sqrt(pow(x_size,2) + pow(z_size,2))

xlabel = 'r(m)'
ylabel = 'z(m)'
curve_list = []

for i_ray in range(len(rays_s_list)):
    lbl = 'ray ' + str(i_ray + 1)
    new_curve = XY_curve(rays_r_list[i_ray], rays_z_list[i_ray], label = lbl)
    curve_list.append(new_curve)

if set_XY_lim in ['True', 'true', 'T']:
    plotZX = XY_Curves_Fig(curve_list, title, xlabel, ylabel, figsize=figsize, \
             ylim = [zmin,zmax], xlim = [xmin,xmax], aspect_ratio = 'equal')
else:
    plotZX = XY_Curves_Fig(curve_list, title, xlabel, ylabel, figsize=figsize, aspect_ratio = 'equal')

# Add k vectors at selected points
if num_plot_k_vectors > 0:
    for i_ray in range(len(rays_s_list)):
        indices = n_evenly_spaced_integers(num_plot_k_vectors, rays_npoints_list[i_ray])
        for i in indices:
            if scale_k_vec in ['True', 'true', 'T']:
                kr = k_vec_base_length*diagonal*rays_kr_list[i_ray][i]/k_max
                kz = k_vec_base_length*diagonal*rays_kz_list[i_ray][i]/k_max
            else:
                kr = k_vec_base_length*diagonal*rays_kr_list[i_ray][i]/rays_knorm_list[i_ray][i]
                kz = k_vec_base_length*diagonal*rays_kz_list[i_ray][i]/rays_knorm_list[i_ray][i]

            if debug > 1: print('kr = ', kr, ' kz = ', kz)
            plt.arrow(rays_r_list[i_ray][i], rays_z_list[i_ray][i],\
              kr, kz, shape='full', head_width = 0.01)

# Plot plasma boundary from R_boundary,Z_boundary

R_boundary = dict_variable_to_list_of_floats(graphics_variable_dict, 'R_boundary')
Z_boundary = dict_variable_to_list_of_floats(graphics_variable_dict, 'Z_boundary')
for i in range(len(R_boundary)):
    print('R_boundary[i] = ', R_boundary[i], '   Z_boundary[i] = ', Z_boundary[i] )

lbl = ''
new_curve = XY_curve(R_boundary, Z_boundary, label = lbl)
curve_list.append(new_curve)

n_Bpoints = len(R_boundary)
print(' ')
#print('R_boundary = ', R_boundary)
print('min(R_boundary) = ', min(R_boundary))
print('max(R_boundary) = ', max(R_boundary))
print(' ')
#print('Z_boundary = ', Z_boundary)
print('min(Z_boundary) = ', min(Z_boundary))
print('max(Z_boundary) = ', max(Z_boundary))


plot_XY_Curves_Fig(plotZX)

#----------------------------------------------------------------------------------------------
# Generate XY ray plot
#----------------------------------------------------------------------------------------------

xmax = float(graphics_variable_dict['box_rmax'])
ymax = xmax
xmin = -xmax
ymin = -xmax

num_plot_k_vectors = int(graphics_variable_dict['num_plot_k_vectors'])
scale_k_vec = graphics_variable_dict['scale_k_vec']
set_XY_lim = graphics_variable_dict['set_XY_lim']
inner_bound = float(graphics_variable_dict['inner_bound'])
outer_bound = float(graphics_variable_dict['outer_bound'])

print('run_description = ', run_description)
print('run_label = ', run_label)
print('xmin = ', xmin, ' xmax = ', xmax, 'ymin =', ymin, ' ymax = ', ymax )

print('num_plot_k_vectors = ', num_plot_k_vectors)
print('scale_k_vec = ', scale_k_vec)
print('set_XY_lim = ', set_XY_lim)

# Generate plot using calls to plot_XY_Curves.py

y_size = max_size
x_size = max_size
figsize = (y_size, x_size)
diagonal = math.sqrt(pow(x_size,2) + pow(y_size,2))

xlabel = 'x(m)'
ylabel = 'y(m)'
curve_list = []

for i_ray in range(len(rays_s_list)):
    lbl = 'ray ' + str(i_ray + 1)
    new_curve = XY_curve(rays_x_list[i_ray], rays_y_list[i_ray], label = lbl)
    curve_list.append(new_curve)

if set_XY_lim in ['True', 'true', 'T']:
    plotXY = XY_Curves_Fig(curve_list, title, xlabel, ylabel, figsize=figsize, \
             ylim = [ymin,ymax], xlim = [xmin,xmax])
else:
    plotXY = XY_Curves_Fig(curve_list, title, xlabel, ylabel, figsize=figsize)

# Add k vectors at selected points
if num_plot_k_vectors > 0:
    for i_ray in range(len(rays_s_list)):
        indices = n_evenly_spaced_integers(num_plot_k_vectors, rays_npoints_list[i_ray])
        for i in indices:
            if scale_k_vec in ['True', 'true', 'T']:
                kx = k_vec_base_length*diagonal*rays_kx_list[i_ray][i]/k_max
                ky = k_vec_base_length*diagonal*rays_ky_list[i_ray][i]/k_max
            else:
                kx = k_vec_base_length*diagonal*rays_kx_list[i_ray][i]/rays_knorm_list[i_ray][i]
                ky = k_vec_base_length*diagonal*rays_ky_list[i_ray][i]/rays_knorm_list[i_ray][i]

            if debug > 1: print('kx = ', kx, ' ky = ', ky)
            plt.arrow(rays_x_list[i_ray][i], rays_y_list[i_ray][i],\
              kx, ky, shape='full', head_width = 0.01)

fig = plt.gcf()
ax = fig.gca()
ax.set_aspect('equal')
circle1 = plt.Circle((0, 0), inner_bound, color='black', fill=False)
ax.add_patch(circle1)
circle2 = plt.Circle((0, 0), outer_bound, color='black', fill=False)
ax.add_patch(circle2)

plot_XY_Curves_Fig(plotXY)

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

