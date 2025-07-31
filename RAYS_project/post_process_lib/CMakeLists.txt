set (lib_sources "")
list (APPEND lib_sources
	slab_processor_m.f90
	solovev_processor_m.f90
	axisym_toroid_processor_m.f90
	mirror_processor_m.f90
	deposition_profiles_m.f90
	OX_conv_analysis_m.f90
	post_processing_m.f90
)

add_library (post_process_lib ${lib_sources})

target_include_directories (post_process_lib
                            PUBLIC
                            $<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
                            $<TARGET_PROPERTY:math_functions,BINARY_DIR>
)
target_link_libraries (post_process_lib
                      PUBLIC
                      RAYS_lib math_functions spline
                      $<$<BOOL:${NetCDF_FOUND}>:NetCDF::NetCDF>
)
