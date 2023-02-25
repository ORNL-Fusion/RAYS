add_executable (component_tests)

target_sources (component_tests
		PRIVATE
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/component_tests_RAYS.f90>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/initialize_components_in_RAYS_lib.f90>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/damping_slab_test_m.f90>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/read_results_LD_test_m.f90>
)

target_include_directories (component_tests
			PUBLIC
			$<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
)
target_link_libraries (component_tests
		    PUBLIC
		    RAYS_lib
)
