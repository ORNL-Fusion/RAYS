set (lib_sources "")
list (APPEND lib_sources
	bcspline.f90
	ibc_ck.f90
	v_spline.f90
	bcspeval.f90
	cspline.f90
	cspeval.f90
	constants_m.f90
	splinck.f90
	zonfind.f90)

add_library (spline_lib ${lib_sources})

#********************************************
add_executable (test_pspline)

target_sources (test_pspline
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_pspline.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/spline_test_fn.f90>
 )

target_include_directories (test_pspline
                            PUBLIC
                            $<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
                            $<TARGET_PROPERTY:spline_lib,BINARY_DIR>
)
target_link_libraries (test_pspline
                       PUBLIC
                       RAYS_lib
                       spline_lib
)
