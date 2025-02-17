#! /usr/bin/env python

"""
plot_RAYS.py -> Plots ray trajectories and k vectors from data ray_results netCDF file
DBB 11/19/2021

"""
# Working notes:
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
# graphics_description_axisym_toroid.dat.
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

#----------------------------------------------------------------------------------------------

# Solovev equilibrium ploloidal flux function defined on np.meshgrid(r, z)
# Returns psi on this grid
def psi(R, Z, rmaj, kappa, outer_bound):
 shape = np.shape(R)
 psiB = (outer_bound**2 - rmaj**2)**2/(2.0*rmaj)**2
 R0 = rmaj*np.ones(shape)
 psi = ((R*Z/(rmaj*kappa))**2 + (R**2 -R0**2 )**2/(2.*rmaj)**2)/psiB
 return psi

#
#----------------------------------------------------------------------------------------------
# Main program
#----------------------------------------------------------------------------------------------

# Get data from graphics description input file
graphics_variable_dict = input_file_to_variable_dict('graphics_description_axisym_toroid.dat')
if debug > 1: print('graphics_variable_dict = ', graphics_variable_dict)

run_description = graphics_variable_dict['run_description']
run_label = graphics_variable_dict['run_label']
num_plot_k_vectors = int(graphics_variable_dict['num_plot_k_vectors'])
scale_k_vec = graphics_variable_dict['scale_k_vec']
k_vec_base_length = float(graphics_variable_dict['k_vec_base_length'])
set_XY_lim = graphics_variable_dict['set_XY_lim']


xmin = float(graphics_variable_dict['box_rmin'])
xmax = float(graphics_variable_dict['box_rmax'])
zmin = float(graphics_variable_dict['box_zmin'])
zmax = float(graphics_variable_dict['box_zmax'])

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

Rz_curve_list = []
xy_curve_list = []
n_all_rays = 0

for file in results_file_list:

    print('Processing CDF file ', file)
    CDF = Dataset(file, 'r', format = 'NETCDF3_CLASSIC')
    CDF_dims = CDF.dimensions

    n_rays = len(CDF_dims['number_of_rays'])
    ray_vec = ma.getdata(CDF.variables['ray_vec'])
    npoints = ma.getdata(CDF.variables['npoints']).tolist()
    print('type(npoints) = ', type(npoints))
    print('npoints = ', npoints)
    print('ray_vec.shape = ', ray_vec.shape)

    x_draw_list = []; y_draw_list = []; z_draw_list = []; R_draw_list = []
    kx_draw_list = []; ky_draw_list = []; kz_draw_list = []; kr_draw_list = []

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
        knorm = [math.sqrt(pow(kx[i],2) + pow(ky[i],2) + pow(kz[i],2))\
                for j in range(npoints[i])]
        k_max = max(knorm)

#         print('npoints[i] = ', npoints[i])
#         print('len(x) = ',len(x), ' x = ', x[0:10])
#         print('len(kx) = ',len(kx), ' kx[0:10] = ', kx[0:10])
#         print('len(z) = ',len(z), ' z = ', z[0:10])
#         print('type(x) = ', type(x))
#         print('type(y) = ', type(y))
#         print('type(z) = ', type(z))
#         print('type(R) = ', type(R))
#         print('R = ', R[0:10])

        lbl = 'ray ' + str(n_all_rays)
        new_curve = XY_curve(R, z, label = lbl)
        Rz_curve_list.append(new_curve)
        new_curve = XY_curve(x, y, label = lbl)
        xy_curve_list.append(new_curve)

# Get data to draw k vectors if doing that
        if num_plot_k_vectors > 0:
            print('type(num_plot_k_vectors) = ', type(num_plot_k_vectors))
            print('type(npoints[i]) = ', type(npoints[i]))
            indices = n_evenly_spaced_integers(num_plot_k_vectors, npoints[i])
            x_draw = [x[j] for j in indices]
            y_draw = [y[j] for j in indices]
            z_draw = [z[j] for j in indices]
            R_draw = [R[j] for j in indices]
            if scale_k_vec in ['True', 'true', 'T']:
                print('Scaling k')
                kx_draw = [max_size*k_vec_base_length*kx[j]/k_max for j in indices]
                ky_draw = [max_size*k_vec_base_length*ky[j]/k_max for j in indices]
                kr_draw = [max_size*k_vec_base_length*kr[j]/k_max for j in indices]
                kz_draw = [max_size*k_vec_base_length*kz[j]/k_max for j in indices]
            else:
                print('Not Scaling k')
                kx_draw  = [max_size*k_vec_base_length*kx[j]/knorm[j] for j in indices]
                ky_draw  = [max_size*k_vec_base_length*ky[j]/knorm[j] for j in indices]
                kr_draw  = [max_size*k_vec_base_length*kr[j]/knorm[j] for j in indices]
                kz_draw  = [max_size*k_vec_base_length*kz[j]/knorm[j] for j in indices]

            x_draw_list.append(x_draw)
            y_draw_list.append(y_draw)
            z_draw_list.append(z_draw)
            R_draw_list.append(R_draw)

            kx_draw_list.append(kx_draw)
            ky_draw_list.append(ky_draw)
            kz_draw_list.append(kz_draw)
            kr_draw_list.append(kr_draw)

#           if debug > 1: print('kr = ', kr, ' kz = ', kz)
#           plt.arrow(rays_r_list[i_ray][i], rays_z_list[i_ray][i],\
#             kr, kz, shape='full', head_width = 0.01)

#----------------------------------------------------------------------------------------------
# Open graphics output file.  N.B. run_label comes from graphics description file.
# Edit graphics description that if you want to customize the run label for multiple ray files.
#----------------------------------------------------------------------------------------------

open_file_XY_Curves_Fig('ray_plots.' + run_label + '.pdf')
title = run_description + '  ' + run_label

#----------------------------------------------------------------------------------------------
# Generate R-Z ray plot
#----------------------------------------------------------------------------------------------

xz_ratio = (xmax-xmin)/(zmax-zmin)
z_size = max_size
x_size = z_size*xz_ratio
figsize = (z_size, x_size)

xlabel = 'r(m)'
ylabel = 'z(m)'

if set_XY_lim in ['True', 'true', 'T']:
    plotZX = XY_Curves_Fig(Rz_curve_list, title, xlabel, ylabel, figsize=figsize, \
             ylim = [zmin,zmax], xlim = [xmin,xmax], aspect_ratio = 'equal')
else:
    plotZX = XY_Curves_Fig(Rz_curve_list, title, xlabel, ylabel, figsize=figsize, aspect_ratio = 'equal')

# Add k vectors at selected points
if num_plot_k_vectors > 0:
    for i_ray in range(n_all_rays):
        z_draw = z_draw_list[i]
        R_draw = R_draw_list[i]
        kz_draw  = kz_draw_list[i]
        kr_draw  = kr_draw_list[i]

        for j in range(num_plot_k_vectors):
            plt.arrow(R_draw[j], z_draw[j], kr_draw[j], kz_draw[j], shape='full', head_width = 0.01)

# Plot plasma boundary from R_boundary,Z_boundary

R_boundary = dict_variable_to_list_of_floats(graphics_variable_dict, 'R_boundary')
Z_boundary = dict_variable_to_list_of_floats(graphics_variable_dict, 'Z_boundary')

# for i in range(len(R_boundary)):
#     print('R_boundary[i] = ', R_boundary[i], '   Z_boundary[i] = ', Z_boundary[i] )

lbl = ''
new_curve = XY_curve(R_boundary, Z_boundary, label = lbl, color='black')
Rz_curve_list.append(new_curve)

n_Bpoints = len(R_boundary)
print(' ')
#print('R_boundary = ', R_boundary)
print('min(R_boundary) = ', min(R_boundary))
print('max(R_boundary) = ', max(R_boundary))
print(' ')
#print('Z_boundary = ', Z_boundary)
print('min(Z_boundary) = ', min(Z_boundary))
print('max(Z_boundary) = ', max(Z_boundary))

#----------------------------------------------------------------------------------------------
# Add eq contours to  RZ ray plot
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

    R = ma.getdata(CDF.variables['R'])
    # print('R = ', R)
    Z = ma.getdata(CDF.variables['Z'])

# Stuff for psiN plot
    # print('Z= ', Z)
    psiN = ma.getdata(CDF.variables['psiN'])
#     print('ma.shape(psiN) = ', ma.shape(psiN))

#    levels = [0.25, 0.5, .75, 1.25, 1.5]
    levels = [0.05*i for i in range(21)]
    plt.contour(R, Z, psiN, levels, colors='k', linewidths=0.5, linestyles='dashed')

# Stuff for resonance contours
    gamma_array = ma.getdata(CDF.variables['gamma_array'])
#     print('ma.shape(gamma_array) = ', ma.shape(gamma_array))

    # Do electron resonances
    abs_gamma = np.absolute(gamma_array[0,:,:])
    print('ma.shape(abs_gamma) = ', ma.shape(abs_gamma))
    gamma_min = ma.min(abs_gamma)
    gamma_max = ma.max(abs_gamma)
    levels = [0.5, 1.0]
    plt.contour(R, Z, abs_gamma, levels, colors='red', linewidths=0.5, linestyles='solid')

    # Do ion resonances later

#     print('gamma_min = ', gamma_min, '  gamma_max = ', gamma_max)

#----------------------------------------------------------------------------------------------
# Plot fig RZ plot
#----------------------------------------------------------------------------------------------

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

xlabel = 'x(m)'
ylabel = 'y(m)'

if set_XY_lim in ['True', 'true', 'T']:
    plotXY = XY_Curves_Fig(xy_curve_list, title, xlabel, ylabel, figsize=figsize, \
             ylim = [ymin,ymax], xlim = [xmin,xmax])
else:
    plotXY = XY_Curves_Fig(xy_curve_list, title, xlabel, ylabel, figsize=figsize)

if num_plot_k_vectors > 0:
    for i_ray in range(n_all_rays):
        x_draw = x_draw_list[i]
        y_draw = y_draw_list[i]
        kx_draw  = kx_draw_list[i]
        ky_draw  = ky_draw_list[i]

        for j in range(num_plot_k_vectors):
            plt.arrow(x_draw[j], y_draw[j], kx_draw[j], ky_draw[j], shape='full', head_width = 0.01)

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

