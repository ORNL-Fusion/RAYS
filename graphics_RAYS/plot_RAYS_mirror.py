#! /usr/bin/env python

"""
plot_RAYS.py -> Plots ray trajectories and k vectors from data ray_results netCDF file
DBB 3/14/2026

N.B. This is taken, with only a few changes from plot_RAYS_axisym_toroid.py.  There are
     some possibly informative working notes for that code up to 7/29/2024.  See it.

Because axis of symmetry (z axis) is now horizontal rather than vertical as in toroids,
it makes sense to make z horizontal (i.e. plot x axis) and y vertical (plot y axis). Also
to add a x,z plot.


"""
# Working notes:
#
# DBB (3/16/2026)
# Modified to use the more object oriented version of plt_XY_Curves.  No changes in
# function but better control of things like vectors and added contours.
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
        while i < len(x)-1 and x[i] < j*delta:
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
graphics_variable_dict = input_file_to_variable_dict('graphics_description_mirror.dat')
if debug > 1: print('graphics_variable_dict = ', graphics_variable_dict)

run_description = graphics_variable_dict['run_description']
run_label = graphics_variable_dict['run_label']
num_plot_k_vectors = int(graphics_variable_dict['num_plot_k_vectors'])
scale_k_vec = graphics_variable_dict['scale_k_vec']
k_vec_base_length = float(graphics_variable_dict['k_vec_base_length'])
set_XY_lim = graphics_variable_dict['set_XY_lim']


xmin = -float(graphics_variable_dict['box_rmax'])
xmax = float(graphics_variable_dict['box_rmax'])
ymin = xmin
ymax = xmax
zmin = float(graphics_variable_dict['box_zmin'])
zmax = float(graphics_variable_dict['box_zmax'])

print('run_description = ', run_description)
print('run_label = ', run_label)
print('xmin = ', xmin, ' xmax = ', xmax, 'zmin =', zmin, ' zmax = ', zmax )

print('num_plot_k_vectors = ', num_plot_k_vectors)
print('scale_k_vec = ', scale_k_vec)
print('k_vec_base_length = ', k_vec_base_length)
print('set_XY_lim = ', set_XY_lim)

z_reference = float(graphics_variable_dict['z_reference'])
r_Omode_cut_at_z_ref = float(graphics_variable_dict['r_Omode_cut_at_z_ref'])
print('z_reference = ', z_reference)
print('r_Omode_cut_at_z_ref = ', r_Omode_cut_at_z_ref)

max_size = 13.
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

print ('plot_RAYS_mirror: results files = ', results_file_list)

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

x_draw_list = []
y_draw_list = []
z_draw_list = []

kx_draw_list = []
ky_draw_list = []
kz_draw_list = []


#----------------------------------------------------------------------------------------------
# Cycle through ray file list
#----------------------------------------------------------------------------------------------

xz_curve_list = []
yz_curve_list = []
xy_curve_list = []
n_all_rays = 0

for file in results_file_list:

    print('\nProcessing CDF file ', file)
    CDF = Dataset(file, 'r', format = 'NETCDF3_CLASSIC')
    CDF_dims = CDF.dimensions

    n_rays = len(CDF_dims['number_of_rays'])
    ray_vec = ma.getdata(CDF.variables['ray_vec'])
    npoints = ma.getdata(CDF.variables['npoints']).tolist()

    ray_kmax = 0.
    for i in range(n_rays):
        n_all_rays = n_all_rays +1
        x = ray_vec[i,0:npoints[i],0]
        y = ray_vec[i,0:npoints[i],1]
        z = ray_vec[i,0:npoints[i],2]
        kx = ray_vec[i,0:npoints[i],3]
        ky = ray_vec[i,0:npoints[i],4]
        kz = ray_vec[i,0:npoints[i],5]
        s = ray_vec[i,0:npoints[i],6]
        R = [math.sqrt(pow(x[j],2) + pow(y[j],2)) for j in range(npoints[i])]
        kr = [(x[j]*kx[j]+y[j]*ky[j])/R[j] for j in range(npoints[i])]

        k_norm = [math.sqrt(kx[j]**2+ky[j]**2+kz[j]**2) for j in range(len(kx))]
        ray_kmax = max(ray_kmax, max(k_norm))
        print('ray_kmax = ', ray_kmax)

        lbl = 'ray ' + str(n_all_rays)
        new_curve = XY_curve(z,x, label = lbl)
        xz_curve_list.append(new_curve)
        new_curve = XY_curve(z,y, label = lbl)
        yz_curve_list.append(new_curve)
        new_curve = XY_curve(x,y, label = lbl)
        xy_curve_list.append(new_curve)

# Get data to draw k vectors if doing that
        if num_plot_k_vectors > 0:
#             indices = n_evenly_spaced_integers(num_plot_k_vectors, npoints[i])
            indices = index_intervals(s, num_plot_k_vectors)
            x_draw = [x[j] for j in indices]
            y_draw = [y[j] for j in indices]
            z_draw = [z[j] for j in indices]
            R_draw = [R[j] for j in indices]
            if scale_k_vec in ['True', 'true', 'T']:
                print('Scaling k')
                kx_draw = [max_size*k_vec_base_length*kx[j]/ray_kmax for j in indices]
                ky_draw = [max_size*k_vec_base_length*ky[j]/ray_kmax for j in indices]
                kr_draw = [max_size*k_vec_base_length*kr[j]/ray_kmax for j in indices]
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

                kr_draw  = [max_size*k_vec_base_length*kr[j]/ \
                math.sqrt(pow(kx[j],2) + pow(ky[j],2) + pow(kz[j],2))\
                for j in indices]

            x_draw_list.append(x_draw)
            y_draw_list.append(y_draw)
            z_draw_list.append(z_draw)
#            R_draw_list.append(R_draw)

            kx_draw_list.append(kx_draw)
            ky_draw_list.append(ky_draw)
            kz_draw_list.append(kz_draw)
#            kr_draw_list.append(kr_draw)

#           if debug > 1: print('kr = ', kr, ' kz = ', kz)
#           plt.arrow(rays_r_list[i_ray][i], rays_z_list[i_ray][i],\
#             kr, kz, shape='full', head_width = 0.01)

print('Total number of rays = ', n_all_rays)

#----------------------------------------------------------------------------------------------
# Open graphics output file.  N.B. run_label comes from graphics description file.
# Edit graphics description that if you want to customize the run label for multiple ray files.
#----------------------------------------------------------------------------------------------

open_file_XY_Curves_Fig('ray_plots.' + run_label + '.pdf')
title = run_description + '  ' + run_label

#----------------------------------------------------------------------------------------------
# Generate Z,X ray plot, Z horizontal, X vertical
#----------------------------------------------------------------------------------------------

xz_ratio = (xmax-xmin)/(zmax-zmin)
z_size = max_size
x_size = 3.*z_size*xz_ratio
figsize = (z_size, x_size)
print('figsize = ', figsize)
# aspect = 3.0
aspect = 'equal'
xlabel = 'z(m)'
ylabel = 'x(m)'

kwargs = {'figsize' : (11.,8.5), 'aspect_ratio' : aspect}

if set_XY_lim in ['True', 'true', 'T']:
    kwargs.update({'xlim':[zmin,zmax], 'ylim':[xmin,xmax]} )

# Add k vectors at selected points
if num_plot_k_vectors > 0:
    kwargs.update({'vectors':[z_draw_list, x_draw_list, kz_draw_list, kx_draw_list],\
                   'aspect_ratio':'equal'})

#----------------------------------------------------------------------------------------------
# Add eq contours to  ZX ray plot
#----------------------------------------------------------------------------------------------

CDF_file_name = 'eq_contours.' + run_label + '.nc'
if os.path.exists(CDF_file_name):
    print('Processing CDF file ', CDF_file_name)
    CDF = Dataset(CDF_file_name, 'r', format = 'NETCDF3_CLASSIC')

    CDF_dim_names = list(CDF.dimensions.keys())
    CDF_dims = CDF.dimensions
    CDF_var_names = list(CDF.variables.keys())
#   print('\n***** CDF_dim_names = ' , CDF_dim_names)
#   print('\n***** CDF_dims = ' , CDF_dims)
#   print('\n***** CDF_var_names = ' , CDF_var_names)
#   print('\n***** CDF.variables = ' , CDF.variables)
#   print('***************')

    X = ma.getdata(CDF.variables['X'])
    Y = X
    Z = ma.getdata(CDF.variables['Z'])

# Stuff for AphiN plot
    AphiN = ma.getdata(CDF.variables['AphiN'])
    levels = [0.05*(i) for i in range(10)]
    levels.insert(1,0.025)

    C_plot_list = []
    new_C_plot = XY_contour(Z, X, AphiN, levels, colors='k', linewidths=0.5,  \
                linestyles='dashed')
    C_plot_list.append(new_C_plot)

    levels = [1.0]
    new_C_plot = XY_contour(Z, X, AphiN, levels, colors='green', linewidths=0.5, \
                linestyles='solid')
    C_plot_list.append(new_C_plot)

# Stuff for resonance contours
    gamma_array = ma.getdata(CDF.variables['gamma_array'])
    # Do electron resonances
    abs_gamma = np.absolute(gamma_array[0,:,:])
    gamma_min = ma.min(abs_gamma)
    gamma_max = ma.max(abs_gamma)
    levels = [0.5, 1.0]
    new_C_plot = XY_contour(Z, X, abs_gamma, levels, colors='red', linewidths=0.5, \
                linestyles=['dashed','solid'])
    C_plot_list.append(new_C_plot)

    # Do ion resonances later

# Stuff for cutoff contours
    omega_pN_array = ma.getdata(CDF.variables['omega_pN_array'])
    # Do electron plasma cutoff
    omega_pN = omega_pN_array[0,:,:]
    omega_pN_min = ma.min(omega_pN)
    omega_pN_max = ma.max(omega_pN)
    print('omega_pN_min = ', omega_pN_min, '   omega_pN_max = ', omega_pN_max)
    levels = [1.0]
    new_C_plot = XY_contour(Z, X, omega_pN, levels, colors='blue', linewidths=1.0, \
                linestyles=['solid'])
    C_plot_list.append(new_C_plot)

    kwargs.update({'contours': C_plot_list})
#----------------------------------------------------------------------------------------------
# Plot fig ZX plot
#----------------------------------------------------------------------------------------------
    plotZX = XY_Curves_Fig(xz_curve_list, title, xlabel, ylabel, **kwargs)

plot_XY_Curves_Fig(plotZX)

#----------------------------------------------------------------------------------------------
# Generate Z,Y ray plot, Z horizontal, Y vertical
#----------------------------------------------------------------------------------------------

yz_ratio = (ymax-ymin)/(zmax-zmin)
z_size = max_size
y_size = 3.*z_size*yz_ratio
figsize = (z_size, y_size)
print('figsize = ', figsize)
# aspect = 3.0
aspect = 'equal'
xlabel = 'z(m)'
ylabel = 'y(m)'

kwargs = {'figsize' : (11.,8.5), 'aspect_ratio': 'equal'}

if set_XY_lim in ['True', 'true', 'T']:
    kwargs.update({'xlim':[zmin,zmax], 'ylim':[xmin,xmax] })

# Add k vectors at selected points
if num_plot_k_vectors > 0:
    kwargs.update({'vectors':[z_draw_list, y_draw_list, kz_draw_list, ky_draw_list],\
               'aspect_ratio':'equal'})

#----------------------------------------------------------------------------------------------
# Add eq contours to  ZY ray plot.  N.B. Contours for ZY are exactly the same as for ZX
#----------------------------------------------------------------------------------------------
kwargs.update({'contours': C_plot_list})

#----------------------------------------------------------------------------------------------
# Plot fig ZY plot
#----------------------------------------------------------------------------------------------

plotZY = XY_Curves_Fig(yz_curve_list, title, xlabel, ylabel, **kwargs)

plot_XY_Curves_Fig(plotZY)

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
# inner_bound = float(graphics_variable_dict['inner_bound'])
# outer_bound = float(graphics_variable_dict['outer_bound'])

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

xlabel = 'x(m)'
ylabel = 'y(m)'


if set_XY_lim in ['True', 'true', 'T']:
    kwargs.update({'ylim':[ymin,ymax], 'xlim':[xmin,xmax]} )

# Add k vectors at selected points
if num_plot_k_vectors > 0:
    kwargs.update({'vectors':[x_draw_list, y_draw_list, kx_draw_list, ky_draw_list],\
               'aspect_ratio':'equal'})

if (r_Omode_cut_at_z_ref > 0.):
    circle_list = [XY_circle(0, 0, r_Omode_cut_at_z_ref, color='blue')]
    kwargs.update({'circles': circle_list})

plotXY = XY_Curves_Fig(xy_curve_list, title, xlabel, ylabel, **kwargs)

plot_XY_Curves_Fig(plotXY)


#----------------------------------------------------------------------------------------------
# Finalize
#-----------------------------------------------------------------------------------------

close_file_XY_Curves_Fig()

