#! /usr/bin/env python

"""
Plot_ray_diags.py -> Plot ray diagnostic data from RAYS. It reads a netCDF file generated
by the post_process_m module and produces plots which are written to a .pdf file.
The code takes one command-line argument -> a  CDF file path, expected extension is ".nc"

The input netCDF file must contain a set of variables to be plotted of the form:

var(max_number_of_points, nray)

The variables are plotted 1 variabale per pdf page, one curve for each ray, with point numbers
as x-axis

This script requires external modules:
    netCDF4
    plt_XY_curves - which also requires matplotlib
    numpy

change log:
 2/27/2024
 version 1.0

"""

import sys
import os
import numpy as np
import numpy.ma as ma

from plt_XY_Curves import *
from netCDF4 import *

debug = True

#----------------------------------------------------------------------------------------------
# Utility functions
#----------------------------------------------------------------------------------------------

# None so far

#----------------------------------------------------------------------------------------------
# Set up
#----------------------------------------------------------------------------------------------

# get the command line -> .config or .nc file name
if len(sys.argv) != 2:
    print(' sys.argv = ', sys.argv)
    message = 'Usage: this script takes one command line argument -> a CDF file path'
    print(message)
    raise Exception(message)

# Get CDF file name and generate plot file name
CDF_file_name = sys.argv[1]
if sys.argv[1][-3:] == '.nc':
    plot_file_name = CDF_file_name[:-3] + '.pdf'
    print('plot_file_name = ', plot_file_name)
else:
    print('Unrecognized cammand line argument.  Should be a .nc or .config file')
    raise exception('Unrecognized cammand line argument.  Should be a .nc or .config file')

open_file_XY_Curves_Fig(plot_file_name)

all_CDF_files_dict = {}
all_CDF_prof_names = []  # List of profles in all files

var_names = []
var_dict = {}
grid_names_dict ={}
grid_dict = {}
index = []
plot_count = 0


print('Processing CDF file ', CDF_file_name)
CDF = Dataset(CDF_file_name, 'r', format = 'NETCDF3_CLASSIC')

CDF_dim_names = list(CDF.dimensions.keys())
CDF_dims = CDF.dimensions
CDF_var_names = list(CDF.variables.keys())
# print('\n***** CDF_dim_names = ' , CDF_dim_names)
# print('\n***** CDF_dims = ' , CDF_dims)
# print('\n***** CDF_var_names = ' , CDF_var_names)
# print('\n***** CDF.variables = ' , CDF.variables)
# print('***************')

n_vars = len(CDF_var_names)
n_rays = len(CDF_dims['number_of_rays'])
print('n_vars = ' , n_vars, '  n_rays = ' , n_rays,)
print('Plotting variables \n')

for name in CDF_var_names:
    if name != 'date_vector' and name != 'npoints': # Skip these they are not plotable
        print('name = ', name)
        curve_list = []
        for i in range(n_rays):
            n_points = ma.getdata(CDF.variables['npoints'][i])
            y =  ma.getdata(CDF.variables[name][i, 0:n_points])
            x = [1.*i for i in range(n_points)]
#             print('n_points = ', n_points)
#             print('len(x) = ',len(x), ' x = ', x[0:10])
#             print('len(y) = ',len(y), ' y = ', y[-1])

            lbl = 'ray ' + str(i + 1)
            new_curve = XY_curve(x, y, label = lbl)
            curve_list.append(new_curve)

        xlabel = 'Step Number'
        ylabel = name
        title = name
        plot1 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
        plot_XY_Curves_Fig(plot1)


# Global attributes defined in this CDF file
# RAYS_run_label = CDF.RAYS_run_label
# date_vector = list(CDF.date_vector)
#
# global_attributes = [ ['CDF file = ', CDF_file_name],\
#                       ['RAYS_run_label label = ', RAYS_run_label]]

#     global_attributes = [ ['CDF file = ', file],\
#                           ['RAYS_run_label label = ', RAYS_run_label], \
#                           ['date_vector = ', date_vector] ]

#
# #----------------------------------------------------------------------------------------------
# # Finalize
# #----------------------------------------------------------------------------------------------
#
# plot_index(index, 1)
#
close_file_XY_Curves_Fig()
#
# if debug:
# #     print('index = ', index
