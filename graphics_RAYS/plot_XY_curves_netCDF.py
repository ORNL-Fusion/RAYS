#! /usr/bin/env python

"""
plot_XY_curves_netCDF.py -> Plot curves from a netCDF file

Takes one command-line argument -> a netCDF file path, expected extension is ".nc"

Plots are written to a .pdf file with the same name as the input .nc file with .nc
extension replaced by .pdf.  Plots one curve per page.

Expects a pdf.file of the form:

Dimensions:
n_curves = number of curves to be written
n_grid(n_CURVES) = number of grid points for each curve
str_max_len = maximum length of a name string -> fortran .nc read needs this for allocation

Variables:
grid(n_grid,n_curves) = x grid for each curve
curve(n_grid,n_curves) = y data for each curve
grid_name(n_curves)
curve_name(n_curves)

Requires external modules:
    netCDF4
    plt_XY_curves - which also requires matplotlib
    numpy
    numpy.ma

change log:
 1/29/2025
 version 1.0

"""

import sys
import os
import numpy as np
import numpy.ma as ma

from plt_XY_Curves import *
from netCDF4 import *

debug = False

#----------------------------------------------------------------------------------------------
# Utility functions
#----------------------------------------------------------------------------------------------

# None so far

#----------------------------------------------------------------------------------------------
# Set up
#----------------------------------------------------------------------------------------------

# get the command line -> .nc file name
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
    msg = 'Unrecognized cammand line argument.' + CDF_file_name + '  Should be a .nc file'
    print(msg)
    raise exception(msg)

open_file_XY_Curves_Fig(plot_file_name)

all_CDF_files_dict = {}
all_CDF_curve_names = []  # List of profles in all files

curve_names = []
curve_dict = {}
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
# print('\n***** CDF.variable n_grid = ' , CDF.variables['n_grid'][1])
# print('***************')

n_curves = CDF.dimensions['n_curves'].size
print('n_curves = ' , n_curves)

for i in list(range(n_curves)):
    n_grid = CDF.variables['n_grid'][i]
    curve_name = str(CDF.variables['curve_name'][i],'utf-8').strip()
    curve_names.append(curve_name)
    curve = ma.getdata(CDF.variables['curve'][i])[0:n_grid]
    curve_dict[curve_name] = curve
    grid_name= str(CDF.variables['grid_name'][i],'utf-8').strip()
    grid_names_dict[curve_name] = grid_name
    grid = ma.getdata(CDF.variables['grid'][i])[0:n_grid]
    grid_dict[curve_name] = grid

if debug:
    print('curve_names = ' , curve_names)
    print('curve_dict = ' , curve_dict)
    print('grid_names = ' , grid_names_dict)
    print('grid_dict = ' , grid_dict)


# # Global attributes defined in this CDF file
# RAYS_run_label = CDF.RAYS_run_label
# date_vector = list(CDF.date_vector)

# global_attributes = [ ['CDF file = ', CDF_file_name],\
#                       ['RAYS_run_label label = ', RAYS_run_label]]
#
# #     global_attributes = [ ['CDF file = ', file],\
# #                           ['RAYS_run_label label = ', RAYS_run_label], \
# #                           ['date_vector = ', date_vector] ]
#
# summary ={'file_name':CDF_file_name , 'global_attributes': global_attributes}
# print('\n************ summary = ',summary)
# plot_summary(summary)


#----------------------------------------------------------------------------------------------
# Plot all curves vs grid
#----------------------------------------------------------------------------------------------
print('Plotting curves \n')

for curve_name in curve_names:
    if debug:
        print('curve_name = ', curve_name)
#         q = [x for x in curve_dict[curve_name]]
#         print('q = ', q)

    xlabel = grid_names_dict[curve_name]
    ylabel = curve_name
    title = curve_name

    plot1 = XY_Curves_Fig(XY_curve(grid_dict[curve_name],curve_dict[curve_name]), title, xlabel, ylabel)
    plot_XY_Curves_Fig(plot1)

#----------------------------------------------------------------------------------------------
# Finalize
#----------------------------------------------------------------------------------------------

plot_index(index, 1)

close_file_XY_Curves_Fig()

if debug:
    print('index = ', index)