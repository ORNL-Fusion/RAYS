cmake_minimum_required (VERSION 3.21)

project (RAYS_project Fortran)

# sets path for executables
set(EXECUTABLE_OUTPUT_PATH ../bin)

add_subdirectory (RAYS_lib)
add_subdirectory (RAYS_code)
add_subdirectory (post_process_RAYS)

# Add the ability to auto generate API documentation with Doxygen
# find_package (Doxygen)
# if (DOXYGEN_FOUND)
# 	configure_file (${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile @ONLY)
# 	add_custom_target (doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile
# 	                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
# 	                   COMMENT "Generating API documentation with Doxygen" VERBATIM)
# endif ()
