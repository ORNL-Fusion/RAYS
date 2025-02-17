#! /usr/bin/env python

"""
plot_psi_contourss.py -> Plot normalized psi from RAYS. It reads a netCDF file generated
by the post_process_m module ('normalized_psi.<run_label>.nc') containing normalized
poloidal flux evaluated on a uniformly spaced  R,Z grid [box_rmin, box_rmax, box_zmin,
box_zmax] and produces plots which are written to a .pdf file with name = input filename.nc

The code takes one command-line argument -> a  CDF file path, expected extension is ".nc"

The input netCDF file must contain:


This script requires external modules:
    netCDF4
    plt_XY_curves - which also requires matplotlib
    numpy

change log:
 7/27/2024
 version 1.0

"""

import sys
import os
import numpy as np
import numpy.ma as ma

from matplotlib import use
# PTB use('TkAgg')
#use('MacOSX')
use('pdf')
import matplotlib.font_manager as font_mgr
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from netCDF4 import *

debug = True

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
    print('Unrecognized cammand line argument.  Should be a .nc or .config file')
    raise exception('Unrecognized cammand line argument.  Should be a .nc or .config file')

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

xmin = ma.getdata(CDF.variables['box_rmin'])
# print('xmin = ', xmin)
xmax = ma.getdata(CDF.variables['box_rmax'])
# print('xmax = ', xmax)
zmin = ma.getdata(CDF.variables['box_zmin'])
# print('zmin = ', zmin)
zmax = ma.getdata(CDF.variables['box_zmax'])
# print('zmax = ', zmax)

R = ma.getdata(CDF.variables['R'])
# print('R = ', R)
Z = ma.getdata(CDF.variables['Z'])
# print('Z= ', Z)
psiN = ma.getdata(CDF.variables['psiN'])
# print('psiN = ', psiN)


#----------------------------------------------------------------------------------------------
# Generate R-Z psi contour plot
#----------------------------------------------------------------------------------------------

max_size = 8.
xz_ratio = (xmax-xmin)/(zmax-zmin)
z_size = max_size
x_size = z_size*xz_ratio
figsize = (z_size, x_size)

xlabel = 'r(m)'
ylabel = 'z(m)'

fig, ax = plt.subplots(figsize = figsize)

levels = [0.25, 0.5, .75, 1.0, 1.25, 1.5]
CS = ax.contour(R, Z, psiN, levels, colors='k')
#ax.clabel(CS, inline=True, fontsize=10)
ax.set_title('Normalized psi')
ax.set_xlabel('R')
ax.set_ylabel('Z')
ax.set_aspect('equal')

# if set_XY_lim in ['True', 'true', 'T']:
#     plotZX = XY_Curves_Fig(Rz_curve_list, title, xlabel, ylabel, figsize=figsize, \
#              ylim = [zmin,zmax], xlim = [xmin,xmax], aspect_ratio = 'equal')
# else:
#     plotZX = XY_Curves_Fig(Rz_curve_list, title, xlabel, ylabel, figsize=figsize, aspect_ratio = 'equal')
#
#
#
# plot_XY_Curves_Fig(plotZX)

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

plt.draw()
#plt.show()
plt.savefig(plot_file_name, format = 'pdf' )
