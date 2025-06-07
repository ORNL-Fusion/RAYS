set (lib_sources "")
list (APPEND lib_sources
	B_loop_m.f90
	mirror_magnetics_m.f90
	../../RAYS_lib/diagnostics_m.f90
	../../RAYS_lib/constants_m.f90
)

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "../../lib")
add_library (mirror_magnetics_lib ${lib_sources})
target_include_directories (mirror_magnetics_lib
                            PUBLIC
                            $<TARGET_PROPERTY:math_functions,BINARY_DIR>
                            $<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
)
target_link_libraries (mirror_magnetics_lib
                      PUBLIC
                      RAYS_lib math_functions
)


#********************************************
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../tests")
add_executable (test_B_loop)

target_sources (test_B_loop
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_B_loop.f90>
 )

target_link_libraries (test_B_loop
                       PUBLIC
                       math_functions
                       mirror_magnetics_lib
)

#********************************************
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin")
add_executable (mirror_magnetics)

target_sources (mirror_magnetics
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mirror_magnetics.f90>
 )

target_include_directories (mirror_magnetics
                            PUBLIC
                            $<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
                            $<TARGET_PROPERTY:mirror_magnetics_lib,BINARY_DIR>
)
target_link_libraries (mirror_magnetics
                       PUBLIC
                       RAYS_lib
                       mirror_magnetics_lib
)
