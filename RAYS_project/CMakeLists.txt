cmake_minimum_required (VERSION 3.21)

project (RAYS_project Fortran)

# sets path for executables
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "../../lib")

# These only make executables that go into /bin
add_subdirectory (RAYS_lib)
add_subdirectory (RAYS_code)
add_subdirectory (ray_scan)
add_subdirectory (post_process_RAYS)

# These make executables that go into both /bin and /tests
# Check CMAKE_RUNTIME_OUTPUT_DIRECTORY inside the subdirectories

add_subdirectory (solovev_2_eqdsk)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../tests")

# These only make executables that go into /tests

add_subdirectory (splines_lib)
add_subdirectory (math_functions_lib)

# Add the ability to auto generate API documentation with Doxygen
# find_package (Doxygen)
# if (DOXYGEN_FOUND)
# 	configure_file (${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile @ONLY)
# 	add_custom_target (doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile
# 	                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
# 	                   COMMENT "Generating API documentation with Doxygen" VERBATIM)
# endif ()
