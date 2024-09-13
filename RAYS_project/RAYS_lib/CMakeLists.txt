set (lib_sources "")
list (APPEND lib_sources
	get_unit_number.f90
	RK4_ode_m.f90
	SG_ode_m.f90
	check_save.f90
	constants_m.f90
	damping_m.f90
	damp_fund_ECH.f90
	deallocate.f90
	deriv_cold.f90
	deriv_num.f90
	diagnostics_m.f90
	disp_solve_cold_n1sq_vs_n3.f90
	dispersion_solvers_m.f90
	eqn_ray.f90
	equilibrium_m.f90
	finalize_run.f90
	initialize_ode_vector.f90
	intialize.f90
	ode_RAYS.f90
	ode_m.f90
	ray_init_m.f90
	ray_tracing.f90
	rf_m.f90
	simple_slab_ray_init_m.f90
	slab_eq_m.f90
	axisym_toroid_eq_m.f90
    solovev_magnetics_m.f90
    axisym_toroid_ray_init_nphi_ntheta_m.f90
    axisym_toroid_ray_init_R_Z_nphi_ntheta_m.f90
	solovev_eq_m.f90
	solovev_ray_init_nphi_ntheta_m.f90
	eqdsk_utilities_m.f90
	eqdsk_magnetics_lin_interp_m.f90
	eqdsk_magnetics_spline_interp_m.f90
	species_m.f90
	suscep_m.f90
	ray_results_m.f90
	openmp_m.f90
)

add_library (RAYS_lib ${lib_sources})

target_include_directories (RAYS_lib
                            PUBLIC
                            $<TARGET_PROPERTY:math_functions,BINARY_DIR>
                            $<TARGET_PROPERTY:spline,BINARY_DIR>
)

target_link_libraries(RAYS_lib
                      PUBLIC
                      spline math_functions
                      $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
                      $<$<BOOL:${NetCDF_FOUND}>:NetCDF::NetCDF>
)
#
# target_sources (RAYS_lib
#                 PRIVATE
#                 $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/*.f90>
# )
# target_include_directories (RAYS_lib
#                             PUBLIC
#                             $<TARGET_PROPERTY:RAYS,BINARY_DIR>
# )
