set (lib_sources "")
list (APPEND lib_sources
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/zfunctions_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/quad_trapezoid_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/bin_to_uniform_grid_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/bisect_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/monotonic_function_inversion.f90>
)
add_library (math_functions ${lib_sources})

#********************************************

#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../tests")

add_executable (test_zfun)

target_sources (test_zfun
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_zfun.f90>
)

target_link_libraries (test_zfun
                       PUBLIC
                       math_functions
                       spline
)

add_executable (test_trapezoid)

target_sources (test_trapezoid
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_trapezoid.f90>
)

target_link_libraries (test_trapezoid
                       PUBLIC
                       math_functions
)


add_executable (test_uniform_grid_binner)

target_sources (test_uniform_grid_binner
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_uniform_grid_binner.f90>
)


target_link_libraries (test_uniform_grid_binner
                       PUBLIC
                       math_functions
)

add_executable (test_bisect)

target_sources (test_bisect
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_bisect.f90>
)

target_link_libraries (test_bisect
                       PUBLIC
                       math_functions
)

add_executable (test_fun_invert)

target_sources (test_fun_invert
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_fun_invert.f90>
)

target_link_libraries (test_fun_invert
                       PUBLIC
                       math_functions
)

