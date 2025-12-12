set (lib_sources "")
list (APPEND lib_sources
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/zfunctions_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/quad_trapezoid_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/bin_to_uniform_grid_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/bisect_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/monotonic_function_inversion.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/complete_elliptic_int_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vectors3_m.f90>
#                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solve_Real_Poly_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/rpoly_m.f90>
)
add_library (math_functions ${lib_sources})

target_include_directories (math_functions
                            PUBLIC
                            $<TARGET_PROPERTY:spline,BINARY_DIR>
)

target_link_libraries(math_functions
                      PUBLIC
                      spline
)

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

add_executable (test_ellip_int)

target_sources (test_ellip_int
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_ellip_int.f90>
)

target_link_libraries (test_ellip_int
                       PUBLIC
                       math_functions
)

add_executable (test_vectors3)

target_sources (test_vectors3
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_vectors3.f90>
)

target_link_libraries (test_vectors3
                       PUBLIC
                       math_functions
)

add_executable (test_rpoly)

target_sources (test_rpoly
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_rpoly.f90>
)

target_link_libraries (test_rpoly
                       PUBLIC
                       math_functions
)
