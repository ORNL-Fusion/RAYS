#! /usr/bin/env python

"""
Code to document header information for the fortran modules in a static library.
This code assumes that the modules follow the template format -> module_template_m.f90
It does not check that or try to accomodate if it doesn't follow the template.

The header consists of:
1) the top-level comments <--> module description
2) Working notes comments
3) declaration of module data (not including variables that are input from namelists),
   definitions of derived types, and interfaces
4) declaration of variables that come in from namelists.

i.e. everything above the 'contains' line if there is one.

For more details see the module template -> module_template_m.f90
"""

# Working notes:
#

import sys
import os
import glob

debug = False
debug0 = True

key_words = ['real', 'complex', 'integer', 'character', 'logical',\
             'real(', 'complex(', 'integer(', 'character(', 'logical(']

#---------------------------------------------------------------------------------------
# classes and functions
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# Open an input file and return the lines
def get_lines(filename):
	file = open(filename, 'r')
	lines = file.readlines()
	file.close()
	return lines

#---------------------------------------------------------------------------------------
# Open an output file and write lines into it
def put_lines(filename, lines):
	file = open(filename, 'w')
	file.writelines(lines)
	file.close()

#---------------------------------------------------------------------------------------
# Convert a list of lines (ending with \n) to a string
def lines_to_string(lines_list):
	string = ''
	for line in lines_list:
		string = string + line + '\n'
	return string

#---------------------------------------------------------------------------------------
# Get next marker line (i.e. next line of form !____ ..., at least 20 underscores).
# Returns line number of found marker line, -1 if none found by end of lines list.
# Remember line numbers start with 0, so might be 1 less that you expect
def get_next_marker_line(lines, n_start_search):
	line_number = n_start_search
	for line in lines[n_start_search:]:
		if debug: print('get_next_marker_line: line_number =', line_number)
		index = line.find('!____________________')
		if index >= 0: # Found marker line
			if debug: print('next marker: return line_number = ', line_number)
			return line_number

		else: # Not a marker line
			if line_number >= len(lines)-1: # No more marker lines in lines list
#				print('next marker 2: type(line_number) = ', type(line_number), '  line_number = ', line_number)
				return -1
		line_number += 1
#		print('next marker 3: type(line_number) = ', type(line_number), '  line_number = ', line_number)

#---------------------------------------------------------------------------------------
# Get next section marker (i.e. 3 lines of the form [marker line, !_section name,
# marker line].
#
# Starts search from n_start_search
# If it finds a valid 3 line section marker it returns the number of the last line
# of marker + 1 (i.e. the first line of the data). It also returns the name of the
# section found or 'none' if it didn't find a section.

def get_next_section(lines, n_start_search):
	line_number = n_start_search
	section_name = 'none'

	while True:
#		print('before get_next_marker_line, line_number = ',line_number)
		line_1_number = get_next_marker_line(lines, line_number)
#		print('line_1_number = ',line_1_number)
		if line_1_number == -1: # No more markers to end of lines
			return (-1, section_name)
		if line_1_number >= len(lines)-1: # Marker on last line.Last line. No section.
			return (-1, section_name)

	# Found a marker line. Could be start of a section mark. Look next line.
		line_number = line_1_number + 1
		if line_number >= len(lines)-1: # Last line, can't be part of section marker.
			return (-1, section_name)
		if lines[line_number][0] == '!': # A comment, it could be a section name
			next_line = lines[line_number][1:].strip()
#			print('get_next_section: next_line = ', next_line, '  line number = ', line_number)
			if len(next_line) > 0: # Found something in the next line. Call it name
				section_name = next_line
#				print('get_next_section: section_name = ', section_name)
				# Check that line after that is a marker line
				line_number += 1
#				print('get_next_section: line_number = ', line_number)
				n_next_marker = get_next_marker_line(lines, line_number)
#				print('get_next_section: n_next_marker = ', n_next_marker)
				if n_next_marker > line_number: # This is not a section marker
					line_number = n_next_marker # Restart search from this line
#					print('really line_number = ', line_number)
					continue
				line_number += 1 # This thing was a section marker
				return (line_number, section_name)
		else: # line doesn't start with !.	Not a comment, so not a section marker
				line_number += 1 # Start search from next line

#---------------------------------------------------------------------------------------
# In a list of lines find the comment lines and strip the comment character and blank
# space from the front
def strip_comment_characters(lines):
	stripped_lines = []
	for line in lines:
		if len(line.strip()) == 0: continue # It's blank, leave it alone
		new_line = line.lstrip()
		if new_line[0] == '!': # It's a comment
			new_line = new_line[1:-1].strip()
		stripped_lines.append(new_line)
	return stripped_lines

#---------------------------------------------------------------------------------------
# In a list of lines find the variable declaration lines and make sure there is a line
# break before it.  This is used in generating the namelist data lines so that declarations
# appear in markdown as a separate line from the comment above it.  Without this
# declaration lines run together with any comments before.  Have to check the line before
# to see if it is already blank or not even there.
# def add_line_breaks(lines):
# 	fixed_lines = []
# 	line_before = ''
# 	is_first_line = True
# 	for line in lines:
# 		if is_first_line: # Just add to list and go on
# 			is_first_line = False
# 			fixed_lines.append(line)
# 			line_before = line
# 			continue
#
# 		stripped_line = line.strip()
# 		if len(stripped_line) == 0 or stripped_line[0]=='!': # Blank or comment, leave it alone
# 			print('add_line_breaks: line = ', line)
# 			fixed_lines.append(line)
# 			line_before = stripped_line
# 			continue
# 		else:
# 			temp_line = stripped_line.split()[0].lower()
# # 			print('temp_line = ', temp_line)
# 			temp_line = temp_line.split('(')[0]
# # 			print('temp_line 1 = ', temp_line)
# 			temp_line = temp_line.split(',')[0]
# # 			print('temp_line2 = ', temp_line)
# 			if temp_line in key_words: # Is a declaration
# 				if len(line_before.strip())  == 0: # Preceded by blank line, go on
# 					fixed_lines.append(stripped_line)
# 					line_before = line
# 					continue
# 				else: # Not preceded by blank line, add line break
# 					stripped_line = '<br>' + stripped_line
# 					fixed_lines.append(stripped_line)
# 					line_before = line
# 	return fixed_lines

def add_line_breaks(lines):
	fixed_lines = []
	for line in lines:

		stripped_line = line.strip()
		if len(stripped_line) == 0: # Blank , leave it alone
			fixed_lines.append(line)
			continue

		temp_line = stripped_line.split()[0].lower()
# 			print('temp_line = ', temp_line)
		temp_line = temp_line.split('(')[0]
# 			print('temp_line 1 = ', temp_line)
		temp_line = temp_line.split(',')[0]
# 			print('temp_line2 = ', temp_line)
		if temp_line in key_words: # Is a declaration, prepend a line break
			line = '<br>' + line
		fixed_lines.append(line)
	return fixed_lines

# def add_line_breaks(lines):
# 	fixed_lines = []
# 	line_before = ''
# 	is_first_line = True
# 	for line in lines:
# 		if is_first_line: # Just add to list and go on
# 			is_first_line = False
# 			fixed_lines.append(line)
# 			line_before = line
# 			continue
#
# 		stripped_line = line.strip()
# 		if len(stripped_line) == 0: # Blank, leave it alone append stripped line
# 			print('add_line_breaks: line = ', line)
# 			fixed_lines.append(stripped_line)
# 			line_before = stripped_line
# 			continue
# 		if len(stripped_line) > 0 and stripped_line[0]=='!': # Comment, leave it alone
# 			print('add_line_breaks: line = ', line)
# 			fixed_lines.append('<br>' + line)
# 			line_before = line
# 			continue
#
# 		temp_line = stripped_line.split()[0].lower()
# # 			print('temp_line = ', temp_line)
# 		temp_line = temp_line.split('(')[0]
# # 			print('temp_line 1 = ', temp_line)
# 		temp_line = temp_line.split(',')[0]
# # 			print('temp_line2 = ', temp_line)
# 		if temp_line in key_words: # Is a declaration
# 			if len(line_before.strip())  == 0: # Preceded by blank line, go on
# 				fixed_lines.append(stripped_line)
# 				line_before = line
# 				continue
# 			else: # Not preceded by blank line, add line break
# 				stripped_line = '<br>' + stripped_line
# 				fixed_lines.append(stripped_line)
# 				line_before = line
# 	return fixed_lines
#

#---------------------------------------------------------------------------------------
# Get module name.	Should be the first line
def get_mod_name(lines):

# If the file is compliant with the template the first line will be the module name
# N.B. It could be a submodule
	index = lines[0].find('module ')
	if index >= 0:
		name = lines[0][index+6: -1].strip()
		if debug: print('\nmodule name = ', name)
		return name
	else:
		if debug0: print('debug: module_name_not_found')
		return 'module_name_not_found'

#---------------------------------------------------------------------------------------
# Get top-level comments for module
# Returns tuple (line number of last data line + 1 = marker line number, top level lines
# list)
def get_top_level(lines, module_name):

	n_lines = len(lines)
	top_level_lines = []
	line_number = 0

# If the file is compliant with the template the first line will be the module name.
# Assuming that, get the next marker line.	It will mark end of top level comments
	marker_line_number = get_next_marker_line(lines, 1)

	if marker_line_number >= n_lines:
		message = 'Nothing but comments in module ' + module_name
		print(message)

	top_level_lines.extend(lines[1: marker_line_number])
	top_level_lines = strip_comment_characters(top_level_lines)

	if debug: print('\ntop_level_lines = ', top_level_lines)
	return (marker_line_number, top_level_lines)

#---------------------------------------------------------------------------------------
# Get module data
# Returns tuple (line number of last data line + 1 = marker line number, module data
# lines list)
def get_module_data(lines, module_name, n_start_search):

	n_lines = len(lines)
	module_data_lines = []

# Find module data section
	while True:
		if debug: print('get_module_data: n_start_search = ', n_start_search)
		n_data_start, section_name = get_next_section(lines, n_start_search)
		if section_name == 'Module data':
			break
		else:
			if n_data_start == -1:
				print('Did not find Module data section for ', module_name)
				return (n_start_search, module_data_lines)
		n_start_search = n_data_start

# Load up Module data lines
	# Find next marker line.  It indicates bottom of Module data section
	marker_line_number = get_next_marker_line(lines, n_data_start)
	if marker_line_number >= n_lines or marker_line_number == -1:
		message = 'Malformed header in module ' + module_name + '.	Hit bottom of file'
		print(message)
		raise Exception(message)

	module_data_lines.extend(lines[n_data_start: marker_line_number])
	module_data_lines = strip_comment_characters(module_data_lines)

	if debug: print('\n module_data_lines = ', module_data_lines)
	return (marker_line_number, module_data_lines)

#---------------------------------------------------------------------------------------
def get_namelist_list(lines, module_name):
# There could be more than one namelist.  So return a dictionary (usually with only
# one entry) containing a lines list for each namelist_lines

	if debug: print('In get_namelist_list', module_name)
	n_lines = len(lines)
	namelist_list = []
	line_number = -1

# Find namelist all namelists
	for line in lines:
		line_number = line_number+1

		if line_number == n_lines-1 and len(namelist_list) > 0:
			if debug: print('get_namelist_list: namelist_list = ', namelist_list)
			return namelist_list

		if line_number == n_lines and len(namelist_list) == 0:
			message = '\nDid not find namelists in module ' + module_name
			print(message)
			return namelist_list

		line = lines[line_number]
		index = line.find('! Namelist data')
		if index >= 0: # Found namelist. now parse out namelist name
			index_1 = line.find('/')
			index_2 = line.find('/', index_1+1)
			name = line[index_1+1: index_2]
			namelist_list.append(name)
#---------------------------------------------------------------------------------------
# Get namelist data
# There could be more than one namelist.  So return a dictionary (usually with only
# one entry) containing a lines list for each namelist_lines
def get_namelist_data(lines, module_name):
	if debug0: print('\nIn get_namelist_data', module_name)

# Don't want to include the namelists themselves, so delete everything below the first
# namelist line.  But need a marker line at the bottom.
	lines_temp = []
	for line in lines:
		if len(line) > 1 and line.lstrip().split()[0] == 'namelist': # First namelist
			break
		else:
			lines_temp.append(line)
	lines = lines_temp
	lines.append('!____________________') # Append marker line


	n_lines = len(lines)
	namelist_data_lines = []
	namelist_dict = {}
	n_start_search = 0

# Find namelist data section
	while True:
		n_data_start, section_name = get_next_section(lines, n_start_search)
		if debug: print('get_namelist_data: section name = ', section_name)

		if n_data_start == -1: # Bottom of file, return
			if debug0: print('\nget_namelist_data: namelist_dict = ',\
			      namelist_dict)
			return namelist_dict

# Check if this section is an namelist data section
		index = section_name.find('Namelist data')
		if index >= 0: # Found a namelist section. now parse out namelist name
			index_1 = section_name.find('/')
			index_2 = section_name.find('/', index_1+1)
			namelist_name = section_name[index_1+1: index_2]
			if debug: print('\nget_namelist_data: namelist_name name = ', namelist_name)

# So this is a namelist section. Get data lines for this namelist
			marker_line_number = get_next_marker_line(lines, n_data_start+1)

			if marker_line_number >= n_lines or marker_line_number == -1:
				message = 'get_namelist_data: Malformed header in module ' +\
						  module_name + '.	Hit bottom of file before marker line'
				print(message)
				raise Exception(message)

# Update namelist_dict

			namelist_data_lines = lines[n_data_start: marker_line_number]
			print('\nnamelist_data_lines = ', namelist_data_lines)
			namelist_data_lines = add_line_breaks(namelist_data_lines)
			print('\nnamelist_data_lines 2 = ', namelist_data_lines)
			namelist_data_lines = strip_comment_characters(namelist_data_lines)

			namelist_dict[namelist_name] = namelist_data_lines
			n_start_search = marker_line_number
			if debug: print('get_namelist_data: namelist_dict.keys = ',\
			      namelist_dict.keys())

		else: # Not a namelist section. Restart search from n_data_start
			n_start_search = n_data_start # Look for another namelist section

#---------------------------------------------------------------------------------------

class lib_data:

	def __init__(self, lib_name, lib_path, **kwargs):

		self.lib_name = lib_name
		if debug: print('\nlib path = ', lib_path)
		os.chdir(lib_path)
		mods = glob.glob("*_m.f90")
		mods.sort(key=str.lower)
		self.module_path_list = mods
		if debug0: print('lib_data: module_path_list = ', self.module_path_list)

		self.mod_names = []
		self.modules = {}
		self.mod_names_to_path = {}

# find modules in lib
		for path in self.module_path_list:
			module_path = os.path.join(lib_path, path)
			if debug:
				print('_________________________________________________________')
				print('\nmodule_path = ', module_path)
				print('_________________________________________________________')

			name = path[: -4].strip()
			self.modules[name] = module_data(module_path)
			mod_name = self.modules[name].module_name
			self.mod_names.append(name)
			self.mod_names_to_path[mod_name] = module_path

		if debug: print('mod_names_to_path = ', self.mod_names_to_path)

#---------------------------------------------------------------------------------------

class module_data:

	def __init__(self, module_path, **kwargs):

		if debug:
			print('_________________________________________________________')
			print('\nmodule_data: module_path = ', module_path)
			print('_________________________________________________________')

		self.module_path = module_path
		temp_lines = get_lines(module_path)
		# Delete lines below header (i.e. "contains" line if there is one)


		self.module_lines = get_lines(module_path)

# Get name of module
		self.module_name = get_mod_name(self.module_lines)
		if debug: print('module_data: self.module_name = ', self.module_name)

# Get top-level comments. i.e. comments just below module name
		end_line_number, section_lines	= \
				get_top_level(self.module_lines, self.module_name)
		self.top_level_comments = section_lines
		if debug: print('\nTop level comments =', self.top_level_comments)

# Get Module data
		n_start_search = end_line_number
		end_line_number, section_lines	= \
				get_module_data(self.module_lines, self.module_name, n_start_search)
		self.module_data = section_lines
		if debug: print('\nModule data =', self.module_data)

# Get list of all namelists in module
		self.namelist_list = get_namelist_list(self.module_lines, self.module_name)
		if debug0: print('\nnamelist list = ', self.namelist_list)

# Get namelist data

		self.namelist_dict = get_namelist_data(self.module_lines, self.module_name)
		if debug: print('\nNamelist data =', self.namelist_dict)


#_________________________________________________________________________________________________

if __name__ == '__main__':

	import os

	print('main')

	RAYS_ROOT = os.environ['RAYS_ROOT']
	print('RAYS_ROOT = ', RAYS_ROOT)

	cwd = os.getcwd()
	out_file_path = os.path.join(cwd,'doc_files')

	lib_list = ['RAYS_lib', 'post_process_lib', 'math_functions_lib', 'splines_lib', \
	            'mirror_magnetics_lib']
# 	lib_list = ['mirror_magnetics_lib']

	print('Libraries: ', lib_list)

	libs = {}

# Instantiate libraries
	for lib in lib_list:
		print('\nlib = ', lib)
		lib_path = os.path.join(RAYS_ROOT, 'RAYS_project', lib)
		if debug: print('main: lib_path = ', lib_path)
		libs[lib] = lib_data(lib, lib_path)

# Write data for each library

	for lib in lib_list:
		print('\nmain: libs[', lib, '].mod_names = ', libs[lib].mod_names)
		print('')

		lib_path = os.path.join(RAYS_ROOT, 'RAYS_project', lib)

		mod_description_lines = ['\n# Module Descriptions for ' + lib]
		mod_data_lines = ['\n# Module data for modules in ' + lib]
# 		namelist_data_lines = ['\n# Namelist Descriptions for modules in ' + lib]
		namelist_data_lines = []

		for mod_name in libs[lib].mod_names:
			mod_description_lines.append('\nDescription of module: ' + mod_name + '\n')
			lines = libs[lib].modules[mod_name].top_level_comments
			mod_description_lines.extend(lines)
			if debug:
				print(libs[lib].modules[mod_name].top_level_comments)
				for line in lines:
					print(line)

			mod_data_lines.append('\nModule data for: ' + mod_name + '\n')
			lines = libs[lib].modules[mod_name].module_data
			mod_data_lines.extend(lines)
			if debug:
				print(libs[lib].modules[mod_name].module_data)
				for line in lines:
					print(line)

			if debug:
				print('\ntype(libs[lib].modules[mod_name].namelist_dict) = ',\
					type(libs[lib].modules[mod_name].namelist_dict))
				print('libs[lib].modules[mod_name].namelist_dict = ',\
					libs[lib].modules[mod_name].namelist_dict)
				print('list(libs[lib].modules[mod_name].namelist_dict.keys()) = ',\
					list(libs[lib].modules[mod_name].namelist_dict.keys()))

			namelist_names = list(libs[lib].modules[mod_name].namelist_dict.keys())
			if len(namelist_names) == 0:
				continue
			namelist_data_lines.append('\n## Namelist data for module: ' + mod_name)
			if debug: print('\nnamelists in ', mod_name, ' = ', namelist_names)

			if len(namelist_names) > 0:

				for namelist in namelist_names:
					namelist_data_lines.append('\n### namelist /' + namelist + '/\n')
					lines = libs[lib].modules[mod_name].namelist_dict[namelist]
					namelist_data_lines.extend(lines)

				for line in namelist_data_lines:
					if debug: print(line)

			if debug:
				for line in lines:
					print(namelist_data_lines)

# Write module description file to output directory
		filename = lib + '_module_description.txt'
		file_path = os.path.join(out_file_path, filename)
		string = lines_to_string(mod_description_lines)
		put_lines(file_path, string)
		file_path = os.path.join(lib_path, 'module_description.txt')
		put_lines(file_path, string)

# Write module data file to output directory
		filename = lib + '_module_data.txt'
		file_path = os.path.join(out_file_path, filename)
		string = lines_to_string(mod_data_lines)
		put_lines(file_path, string)
		file_path = os.path.join(lib_path, 'module_data.txt')
		put_lines(file_path, string)

# Write module data file to output directory
		filename = lib + '_namelist_description.txt'
		file_path = os.path.join(out_file_path, filename)
		string = lines_to_string(namelist_data_lines)
		put_lines(file_path, string)
		file_path = os.path.join(lib_path, 'namelist_description.md')
		put_lines(file_path, string)
