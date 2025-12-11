#! /usr/bin/env python

"""
plot_RAYS.py -> Plots ray trajectories and k vectors from data ray_results netCDF file
DBB 11/19/2021

"""
# Working notes:
#
# DBB (12/11/2025)
# Reworked the vector plots using the new version of plt_XY_curves.py, which itself has
# been reworked to use the object oriented interface of matplotlib.  The vectors are added
# through keyword args to XY_Curves_Fig() instead of being tacked on using the functional
# interface after the plot was generated.  Also added function, index_intervals(), which
# finds indices of equal arclength steps along the ray to place vectors, instead of equal
# number of ray steps.  This gives better plots when the group velocity changes a lot.

# DBB (8/26/2025)
# Changed input from reading the ASCII files 'ray_out.' + run_label to read the netCDF file
# 'ray_results.' + run_label + '.nc'. Lifted straight out of plot_RAYS_axisym_toroid.py
#
# DBB (7/29/2024)
# Adding plot of psi contours and cyclotron resonance contours to the R-Z ray plot figure.
# So far only electron cyclotron resonances are plotted.  Will do ions later.
#
# DBB (7/28/2024)
# Added optional command line args to specify multiple input ray data input netCDF files
# That makes it possible to combine ray plots from multiple runs. (Actually that was done
# a while ago). With no command line args the default file name is:
# 'run_results.' + run_label + '.nc' where run_label comes from file
# graphics_description_axisym_toroid.dat.  To set the plot title for multiple file plots
# you need to edit the graphics_description_slab.dat file.
#
# DBB (4/21/2024)
# Changed input from reading the ASCII files 'ray_out.' + run_label to read the netCDF file
# 'ray_results.' + run_label + '.nc'.  The old version is stashed in spare parts.
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
import matplotlib.cm as cm
import numpy as np
import numpy.ma as ma

from simple_file_editing_functions import get_lines, input_file_to_variable_dict,\
      dict_variable_to_list_of_floats
from plt_XY_Curves import *

debug = 0
k_vec_base_length = 0.005

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
        return

    if n >= Length:
        return list(range(Length))
    elif n == 0:
        return [0]
    elif n == -1:
        return [-1]
    else:
        return [int((i)*float(Length-1)/(n-1)) for i in range(n)]

#----------------------------------------------------------------------------------------------

def index_intervals(x, n):

# Given a list of increasing numbers, x, return n+1 indices which give x intervals
# of approximately the same length. index[0] = 0 => xmin, index[n+1] = len(n) => xmax

    if type(n) != int:
        print('\nindex_intervals error')
        print('n must be integer')
        return

    for i in range(len(x)-1):
        if (x[i+1]-x[i] <= 0.):
            print('\nindex_intervals error')
            print('x must be non-decreasing')
            return

    delta = (x[-1]-x[0])/n
    index = [i for i in range(n+1)]

    if (n >= len(x)): # As many intervals as data points
        print('\nindex_intervals error')
        print('As many intervals as data points, can not split them')
        return

    i = 0
    j = 1
    s = x[0]
    index[0] = 0
    index[-1] = len(x)

    while i < len(x)-1:  # stop before x max
        while x[i] < j*delta:
            i = i +1
        # Crossed an interval, check to see if x[i-1] or x[i] is closer to j*delta
        index[j] = i
        if abs(x[i-1] - j*delta) < abs(x[i] - j*delta):
            index[j] = i-1
        j = j +1

    return index

#
#----------------------------------------------------------------------------------------------
# Main program
#----------------------------------------------------------------------------------------------

# Get data from graphics description input file
graphics_variable_dict = input_file_to_variable_dict('graphics_description_slab.dat')
if debug > 1: print('graphics_variable_dict = ', graphics_variable_dict)

run_description = graphics_variable_dict['run_description']
run_label = graphics_variable_dict['run_label']
num_plot_k_vectors = int(graphics_variable_dict['num_plot_k_vectors'])
scale_k_vec = graphics_variable_dict['scale_k_vec']
k_vec_base_length = float(graphics_variable_dict['k_vec_base_length'])
set_XY_lim = graphics_variable_dict['set_XY_lim']

xmin = float(graphics_variable_dict['xmin'])
xmax = float(graphics_variable_dict['xmax'])
ymin = float(graphics_variable_dict['ymin'])
ymax = float(graphics_variable_dict['ymax'])
zmin = float(graphics_variable_dict['zmin'])
zmax = float(graphics_variable_dict['zmax'])

print('run_description = ', run_description)
print('run_label = ', run_label)
print('xmin = ', xmin, ' xmax = ', xmax, 'zmin =', zmin, ' zmax = ', zmax )

print('num_plot_k_vectors = ', num_plot_k_vectors)
print('scale_k_vec = ', scale_k_vec)
print('k_vec_base_length = ', k_vec_base_length)
print('set_XY_lim = ', set_XY_lim)

max_size = 8.
print('max_size = ', max_size)

# get the command line
results_file_list = []
n_arg = len(sys.argv)
if n_arg == 1: # No arg, get run_label from graphics description file
    n_results_files = 1
    results_file_list.append('run_results.' + run_label + '.nc')

if n_arg > 1: # Get ray file names from command line
    n_results_files = n_arg-1
    results_file_list = sys.argv[1:]

print ('results files = ', results_file_list)

nray = 0
ray_ymin = 0. # Keep track of range of y.  Many plots will be in Z-X plane
ray_ymax = 0.
ray_kmax = 0.

rays_s_list = []
rays_x_list = []
rays_y_list = []
rays_z_list = []
rays_kx_list = []
rays_ky_list = []
rays_kz_list = []
rays_knorm_list = []
rays_npoints_list =[]

x_draw_list = []
y_draw_list = []
z_draw_list = []

kx_draw_list = []
ky_draw_list = []
kz_draw_list = []


#----------------------------------------------------------------------------------------------
# Cycle through ray file list
#----------------------------------------------------------------------------------------------

zx_curve_list = []
zy_curve_list = []
n_all_rays = 0


x_draw_list = []; y_draw_list = []; z_draw_list = []
kx_draw_list = []; ky_draw_list = []; kz_draw_list = []

for file in results_file_list:

    print('\nProcessing CDF file ', file)
    CDF = Dataset(file, 'r', format = 'NETCDF3_CLASSIC')
    CDF_dims = CDF.dimensions

    n_rays = len(CDF_dims['number_of_rays'])
    ray_vec = ma.getdata(CDF.variables['ray_vec'])
    npoints = ma.getdata(CDF.variables['npoints']).tolist()
    print('npoints = ', npoints)
#     print('ray_vec.shape = ', ray_vec.shape)

    for i in range(n_rays):
        n_all_rays = n_all_rays +1
        x = ray_vec[i,0:npoints[i],0]
        y = ray_vec[i,0:npoints[i],1]
        z = ray_vec[i,0:npoints[i],2]
        kx = ray_vec[i,0:npoints[i],3]
        ky = ray_vec[i,0:npoints[i],4]
        kz = ray_vec[i,0:npoints[i],5]
        s = ray_vec[i,0:npoints[i],6]

# Keep track of range of y.  Many plots will be in Z-X plane
        ray_ymin = min(ray_ymin, min(y))
        ray_ymax = max(ray_ymax, max(y))

        k_norm = [math.sqrt(kx[j]**2+ky[j]**2+kz[j]**2) for j in range(len(kx))]
        ray_kmax = max(ray_kmax, max(k_norm))

        lbl = 'ray ' + str(n_all_rays)
        new_curve = XY_curve(z, x, label = lbl)
        zx_curve_list.append(new_curve)
        new_curve = XY_curve(z, y, label = lbl)
        zy_curve_list.append(new_curve)

# Get data to draw k vectors if doing that
        if num_plot_k_vectors > 0:
#             indices = n_evenly_spaced_integers(num_plot_k_vectors, npoints[i])
            indices = index_intervals(s, num_plot_k_vectors)
            x_draw = [x[j] for j in indices]
            y_draw = [y[j] for j in indices]
            z_draw = [z[j] for j in indices]
            if scale_k_vec in ['True', 'true', 'T']:
                print('Scaling k')
                kx_draw = [max_size*k_vec_base_length*kx[j]/ray_kmax for j in indices]
                ky_draw = [max_size*k_vec_base_length*ky[j]/ray_kmax for j in indices]
                kz_draw = [max_size*k_vec_base_length*kz[j]/ray_kmax for j in indices]
            else:
                print('Not Scaling k')
                kx_draw  = [max_size*k_vec_base_length*kx[j]/ \
                math.sqrt(pow(kx[j],2) + pow(ky[j],2) + pow(kz[j],2))\
                for j in indices]

                ky_draw  = [max_size*k_vec_base_length*ky[j]/ \
                math.sqrt(pow(kx[j],2) + pow(ky[j],2) + pow(kz[j],2))\
                for j in indices]

                kz_draw  = [max_size*k_vec_base_length*kz[j]/ \
                math.sqrt(pow(kx[j],2) + pow(ky[j],2) + pow(kz[j],2))\
                for j in indices]

            x_draw_list.append(x_draw)
            y_draw_list.append(y_draw)
            z_draw_list.append(z_draw)

            kx_draw_list.append(kx_draw)
            ky_draw_list.append(ky_draw)
            kz_draw_list.append(kz_draw)

print('Total number of rays = ', n_all_rays)

#----------------------------------------------------------------------------------------------
# Open graphics output file.  N.B. run_label comes from graphics description file.
# Edit graphics description that if you want to customize the run label for multiple ray files.
#----------------------------------------------------------------------------------------------

open_file_XY_Curves_Fig('ray_plots.' + run_label + '.pdf')
title = run_description + '  ' + run_label

#----------------------------------------------------------------------------------------------
# Generate Z-X ray plot
#----------------------------------------------------------------------------------------------

xz_ratio = (xmax-xmin)/(zmax-zmin)
z_size = max_size
x_size = z_size*xz_ratio
figsize = (z_size, x_size)
diagonal = math.sqrt(pow(x_size,2) + pow(z_size,2))

xlabel = 'z(m)'
ylabel = 'x(m)'

kwargs = {'figsize' : (11.,8.5)}
# kwargs = {'figsize' : figsize}
# kwargs = {}

if set_XY_lim in ['True', 'true', 'T']:
    kwargs.update({'xlim':[zmin,zmax], 'ylim':[xmin,xmax]} )

# Add k vectors at selected points
if num_plot_k_vectors > 0:
    kwargs.update({'vectors':[z_draw_list, x_draw_list, kz_draw_list, kx_draw_list],\
                   'aspect_ratio':'equal'})

plotZX = XY_Curves_Fig(zx_curve_list, title, xlabel, ylabel, **kwargs)
plot_XY_Curves_Fig(plotZX)

#----------------------------------------------------------------------------------------------
# Generate ZY ray plot
#----------------------------------------------------------------------------------------------

#Many plots will be restricted to ZX plane.  Check to see y extent at least 1% of max_size
if ray_ymax-ray_ymin > 0.01*max_size:

    zx_ratio = (zmax-zmin)/(xmax-xmin)
    z_size = max_size
    x_size = z_size*xz_ratio
    figsize = (z_size, x_size)

    xlabel = 'z(m)'
    ylabel = 'y(m)'

    kwargs = {'figsize' : figsize, 'aspect_ratio': 'equal'}

if set_XY_lim in ['True', 'true', 'T']:
    kwargs.update({'xlim':[zmin,zmax], 'ylim':[xmin,xmax] })

    # Add k vectors at selected points
    if num_plot_k_vectors > 0:
        kwargs.update({'vectors':[z_draw_list, y_draw_list, kz_draw_list, ky_draw_list],\
                   'aspect_ratio':'equal'})

    plotZY = XY_Curves_Fig(zy_curve_list, title, xlabel, ylabel, **kwargs)
    plot_XY_Curves_Fig(plotZY)

#-----------------------------------------------------------------------------------------
# Finalize
#-----------------------------------------------------------------------------------------

close_file_XY_Curves_Fig()

