 subroutine deallocate

    use diagnostics_m, only :deallocate_diagnostics_m
    use constants_m, only : deallocate_constants_m, rkind
    use species_m, only : deallocate_species_m
    use rf_m, only : deallocate_rf_m
    use damping_m, only : deallocate_damping_m
    use equilibrium_m, only : deallocate_equilibrium_m
    use ray_init_m, only : deallocate_ray_init_m
    use ode_m, only : deallocate_ode_solver_m
    use ray_results_m, only : deallocate_ray_results_m

    implicit none
    

! deallocate all the modules that initialize, allocate, or open files
    call deallocate_constants_m
    call deallocate_species_m
    call deallocate_rf_m  
    call deallocate_damping_m
    call deallocate_equilibrium_m
    call deallocate_ray_init_m  
    call deallocate_ode_solver_m
    call deallocate_ray_results_m
    
    call deallocate_diagnostics_m ! This really is the end

    write(*,*) ' '
    write(*,*) 'RAYS deallocation finished'
    write(*,*) ' '

 end subroutine deallocate