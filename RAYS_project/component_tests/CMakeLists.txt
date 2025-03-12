add_executable (component_tests)

target_sources (component_tests
		PRIVATE
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/component_tests_RAYS.f90>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/initialize_components_in_RAYS_lib.f90>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/damping_slab_test_m.f90>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/read_results_LD_test_m.f90>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dispersion_solver_test_m.f90>
)

target_include_directories (component_tests
			PUBLIC
			$<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
)
target_link_libraries (component_tests
		    PUBLIC
		    RAYS_lib
)


add_executable (write_eq_vec)

target_sources (write_eq_vec
		PRIVATE
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/write_eq_rvec.f90>
)

target_include_directories (write_eq_vec
			PUBLIC
			$<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
)
target_link_libraries (write_eq_vec
		    PUBLIC
		    RAYS_lib
)
