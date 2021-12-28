#! /usr/bin/env python

"""
plot_kx_profiles_slab.py -> kx(x) from data in file kx_profiles_slab.<run_label>
DBB 12/27/2021

"""
# Working notes:
#
# DBB (11/23/2021)
#

import sys
import os
import math
from netCDF4 import *
import matplotlib.pyplot as plt

from simple_file_editing_functions import get_lines, input_file_to_variable_dict
from plt_XY_Curves import *

debug = 4
k_vec_base_length = 0.001

#----------------------------------------------------------------------------------------------
# Get description data needed for the plot
#----------------------------------------------------------------------------------------------

# Graphics description input file
graphics_variable_dict = input_file_to_variable_dict('graphics_description_slab.dat')
if debug > 1: print('graphics_variable_dict = ', graphics_variable_dict)

run_description = graphics_variable_dict['run_description'] 
run_label = graphics_variable_dict['run_label'] 

xmin = float(graphics_variable_dict['xmin'])
xmax = float(graphics_variable_dict['xmax'])

num_plot_k_vectors = int(graphics_variable_dict['num_plot_k_vectors'])
scale_k_vec = graphics_variable_dict['scale_k_vec']
set_XY_lim = graphics_variable_dict['set_XY_lim']

print('run_description = ', run_description)
print('run_label = ', run_label)
print('xmin = ', xmin, ' xmax = ', xmax ) 


# Graphics output file
open_file_XY_Curves_Fig('kx_plots.' + run_label + '.pdf')
figsize = (6.0,4.0)

#----------------------------------------------------------------------------------------------
# Get all the kx data
#----------------------------------------------------------------------------------------------

# data input file
lines = get_lines('kx_profiles_slab.' + run_label) 

lines.append('ray\n')  # Tack on an extra "ray" line to signal the end
if debug > 4: print('len(lines) = ', len(lines))
if debug > 4: 
    for line in lines:
        print(line)

x_list = []
kx_real_plus_list = []
kx_im_plus_list = []
kx_real_minus_list = []
kx_im_minus_list = []
kx_real_fast_list =  []
kx_im_fast_list = []
kx_real_slow_list = []
kx_im_slow_list = []  

for i in range(len(lines)):
    line = lines[i]
    split_line = line.split()
    if debug > 4: print('split_line = ', split_line)
    
    if split_line[0] == 'ray':  # this is a new ray
        nray = split_line[1]
        ny = split_line[3]
        nz = split_line[5]
        if debug > 4: print('nray = ', nray, ' ny = ', ny, 'nz = ', nz)
        x_list = []
        kx_real_plus_list = []
        kx_im_plus_list = []
        kx_real_minus_list = []
        kx_im_minus_list = []
        kx_real_fast_list =  []
        kx_im_fast_list = []
        kx_real_slow_list = []
        kx_im_slow_list = []  

    elif split_line[0] == 'x':  # this is a column heading line
        pass
        
    else: # this is an (x,kx) line
        x_list.append(float(split_line[0]))
        kx_real_plus_list.append(float(split_line[1]))
        kx_im_plus_list.append(float(split_line[2]))
        kx_real_minus_list.append(float(split_line[3]))
        kx_im_minus_list.append(float(split_line[4]))
        kx_real_fast_list.append(float(split_line[5]))
        kx_im_fast_list .append(float(split_line[6]))
        kx_real_slow_list.append(float(split_line[7]))
        kx_im_slow_list.append(float(split_line[8])) 
    
    if lines[i+1].split()[0] == 'ray': # This is the last point of that ray. Do a plot.
    
    #----------------------------------------------------------------------------------------------
    # Do plot for this ray, plus/minus roots
    #----------------------------------------------------------------------------------------------

        title = 'ray = ' + str(nray) + '   ny = ' + str(ny) + '   nz = ' + str(nz)
        xlabel = 'x(m)'
        ylabel = 'kx(m-1)'
        if debug > 3: print('title = ', title, 'xlabel = ', xlabel, 'ylabel = ', ylabel)
        if debug > 4: print('x_list = ', x_list)
        if debug > 4: print('kx_real_plus_list = ', kx_real_plus_list)

        curve_list = []
        new_curve = XY_curve(x_list, kx_real_plus_list, label = 're_kx plus')
        new_curve.setColor('red')
        new_curve.setLinestyle('solid')
        curve_list.append(new_curve)
    
        new_curve = XY_curve(x_list, kx_im_plus_list, label = 'im_kx plus')
        new_curve.setColor('red')
        new_curve.setLinestyle('dashed')
        curve_list.append(new_curve)
    
        new_curve = XY_curve(x_list, kx_real_minus_list, label = 're_kx minus')
        new_curve.setColor('blue')
        new_curve.setLinestyle('solid')
        curve_list.append(new_curve)
    
        new_curve = XY_curve(x_list, kx_im_minus_list, label = 'im_kx minus')
        new_curve.setColor('blue')
        new_curve.setLinestyle('dashed')
        curve_list.append(new_curve)

        plotKX = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
        plot_XY_Curves_Fig(plotKX)
        print('len(curve_list) = ', len(curve_list))

        curve_list = []
        new_curve = XY_curve(x_list, kx_real_fast_list, label = 're_kx fast')
        new_curve.setColor('red')
        new_curve.setLinestyle('solid')
        curve_list.append(new_curve)
    
        new_curve = XY_curve(x_list, kx_im_fast_list, label = 'im_kx fast')
        new_curve.setColor('red')
        new_curve.setLinestyle('dashed')
        curve_list.append(new_curve)
    
        new_curve = XY_curve(x_list, kx_real_slow_list, label = 're_kx slow')
        new_curve.setColor('blue')
        new_curve.setLinestyle('solid')
        curve_list.append(new_curve)
    
        new_curve = XY_curve(x_list, kx_im_slow_list, label = 'im_kx slow')
        new_curve.setColor('blue')
        new_curve.setLinestyle('dashed')
        curve_list.append(new_curve)

        plotKX = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
        plot_XY_Curves_Fig(plotKX)

        print('len(curve_list) = ', len(curve_list))
    
    if i == len(lines)-2: break  # That was the last point on the last ray, exit loop   


  
#----------------------------------------------------------------------------------------------
# Finalize
#----------------------------------------------------------------------------------------------
# 

close_file_XY_Curves_Fig()


