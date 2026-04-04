#! /usr/bin/env python

"""
ls_RAYS_files.py makes soft link between specific RAYS input file to generic file 'rays.in'
Also optionally makes similar soft link to 'post_process_rays.in' from specific post-
process namelist file.

Takes 1 or 2 command line args -> specific RAYS input namelist file name and optionally
specific post process input namelist file.

Deletes generic files 'rays.in' and 'post_process_rays.in' if they exist
Makes a soft link from the specific files to the generic ones.

change log:
3/1/2026
 version 1.0

"""

import sys
import os
import subprocess
from simple_file_editing_functions import get_lines, read_string_var_from_nml_lines
debug = False

#----------------------------------------------------------------------------------------------
# Utility functions
#----------------------------------------------------------------------------------------------

# None so far

#----------------------------------------------------------------------------------------------
# Set up
#----------------------------------------------------------------------------------------------

# Get path to RAYS
rays_root = os.getenv('RAYS_ROOT')
rays_proj = os.getenv('RAYS_PROJECT')
print('sys.argv = ', sys.argv)
if len(sys.argv) > 1:
    rays_in_file = sys.argv[1]
    print('rays_in_file =', rays_in_file)
    cmd = ['rm rays.in']
    retcode = subprocess.call(cmd, shell=True)
    if (retcode != 0):
        logMsg = 'Error deleting rays.in'
        print('cmd = ', cmd)
        print('retcode = ', retcode)
        print('logMsg =  ', logMsg)

    cmd = ['ln -s '+ rays_in_file + ' rays.in']
    retcode = subprocess.call(cmd, shell=True)
    if (retcode != 0):
        logMsg = 'Error linking ' + rays_in_file + ' to rays.in'
        print('cmd = ', cmd)
        print('retcode = ', retcode)
        print('logMsg =  ', logMsg)

if len(sys.argv) > 2:
    pp_in_file = sys.argv[2]
    print('post process in file =', pp_in_file)
    cmd = ['rm post_process_rays.in']
    retcode = subprocess.call(cmd, shell=True)
    if (retcode != 0):
        logMsg = 'Error deleting post_process_rays.in'
        print('cmd = ', cmd)
        print('retcode = ', retcode)
        print('logMsg =  ', logMsg)
    cmd = ['ln -s ' + pp_in_file + ' post_process_rays.in']
    retcode = subprocess.call(cmd, shell=True)
    if (retcode != 0):
        logMsg = 'Error linking ' + pp_in_file + ' to post_process_rays.in'
        print('cmd = ', cmd)
        print('retcode = ', retcode)
        print('logMsg =  ', logMsg)

#     print('2 pp_in_file = ', pp_in_file)
#     if (pp_in_file.strip("'") == 'post_process_X-mode_ny_01.in'):
#         print('3 pp_in_file = ', pp_in_file)
#         raise

if sys.argv == 3:
    print(' sys.argv = ', sys.argv)
    message = 'Usage: this script takes 1 or 2 arguments -> RAYS input file names'
    print(message)
    raise Exception(message)
