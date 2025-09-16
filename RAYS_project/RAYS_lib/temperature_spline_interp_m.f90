module  temperature_spline_interp_m
! Implements temperature profiles based on pointwise numeerical data contained in namelist
! temperature_spline_interp_list.

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

! At present (2-22-2025) this only produces only two profiles (electron, ion).  All ions
! species have the same temperature. Profiles are normalized to one.
! This is similar to the "parabolic" profiles implemented in axisym_toroid_eq_m.
! The actual temperatures are then profile*t0s(i).
! fraction of the electron temperature given by n0s(0:nspec).  In the future, when ion dynamics
! counts (e.g. ICRF) and separate ion temperature data is provided, something more elaborate
! will be necessary.  In particular how to ensure charge neutrality and specify Zeff.
! The temperature values in the namelist are assumed to be on a uniform grid from 0.0 to 1.0
! For psi > 1.0 returns dens = d_scrape_off, dd_psi = 0.
!
! For now only one spline_temperature_model and one function is exported, spline_temperature_ne.
! So, for now, namelist variable 'spline_temperature_model' is not used.
!
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use constants_m, only : rkind
    use quick_cube_splines_m, only : cube_spline_function_1D

    implicit none

! Local data
	integer, parameter :: n_grid_max = 200 ! Max dimension of input arrays, truncated later

! Stuff for 1D spline profiles

    type(cube_spline_function_1D) :: Te_profileN  ! Te profile normalized to 1. on axis
    type(cube_spline_function_1D) :: Ti_profileN  ! Ti profile normalized to 1. on axis
 	character (len = 80) :: Te_name = 'Te_profile'
 	character (len = 80) :: Ti_name = 'Ti_profile'

!_________________________________________________________________________________________
! Namelist data for /temperature_spline_interp_list/
!_________________________________________________________________________________________

    character (len = 60) :: spline_Te_model, spline_Ti_model ! Not used, YET
	integer :: ngrid ! Actual number of points to be splined. <= n_grid_max
	real(KIND=rkind) ::  Te_in(n_grid_max) ! Values on grid (grid assumed uniform 0 to 1)
	real(KIND=rkind) ::  Ti_in(n_grid_max) ! Values on grid (grid assumed uniform 0 to 1)

  namelist /temperature_spline_interp_list/ ngrid, Te_in, Ti_in

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

  subroutine initialize_temperature_spline_interp(read_input)

    use constants_m, only : one
    use species_m, only : nspec
    use diagnostics_m, only : message, message_unit,messages_to_stdout, verbosity

    implicit none

    logical, intent(in) :: read_input

 	integer :: input_unit, get_unit_number ! External, free unit finder

	real(KIND=rkind), allocatable ::  grid(:) ! Grid points (constructed to be 0. to 1.)
	real(KIND=rkind), allocatable ::  Te_values(:) ! Values on grid
	real(KIND=rkind), allocatable ::  Ti_values(:) ! Values on grid

	integer :: i

    if (verbosity >= 0) then
		write(*,*) 'initialize_temperature_spline_interp'
    end if

    if (read_input .eqv. .true.) then
  		input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, temperature_spline_interp_list)
        close(unit=input_unit)
    end if

! Write input namelist
    if (verbosity >= 0) then
		write(message_unit, temperature_spline_interp_list)
		if (messages_to_stdout) write(*, temperature_spline_interp_list)
		call message(1)
    end if

! Allocate arrays
	allocate(grid(ngrid), Te_values(ngrid), Ti_values(ngrid))

! Construct arrays for splining
    do i = 1, ngrid
!    	grid(i) = plasma_psi_limit*(i-1)/(ngrid - 1)
    	grid(i) = one*(i-1)/(ngrid - 1)
    	Te_values(i) = Te_in(i)/Te_in(1) ! Normalized to 1.0 on axis
    	Ti_values(i) = Ti_in(i)/Ti_in(1) ! Normalized to 1.0 on axis
    end do

! Initialize spline coefficients
    call Te_profileN%cube_spline_1D_init(ngrid, grid, Te_values, Te_name)
    call Ti_profileN%cube_spline_1D_init(ngrid, grid, Ti_values, Ti_name)
    return
  end subroutine initialize_temperature_spline_interp

!****************************************************************************************

  subroutine  temperature_spline_interp(psiN, T_scrape_off, Te, dTe_psi, Ti, dTi_psi)

    use constants_m, only : one, zero
    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message

    implicit none

	real(KIND=rkind), intent(in) :: psiN, T_scrape_off
	real(KIND=rkind), intent(out) :: Te, dTe_psi, Ti, dTi_psi

   	Te = zero
  	Ti = zero
    if (psiN <= one) then
		call Te_profileN%eval_1D_fp(psiN, Te, dTe_psi)
		call Ti_profileN%eval_1D_fp(psiN, Ti, dTi_psi)
  	end if
  	if (Te < T_scrape_off) then
    	Te = T_scrape_off
    	dTe_psi = zero
    end if
    if (Ti < T_scrape_off) then
    	Ti = T_scrape_off
    	dTi_psi = zero
    end if

    return
    end subroutine  temperature_spline_interp

!********************************************************************

    subroutine deallocate_temperature_spline_interp_m
		! Nothing to deallocate
		return
    end subroutine deallocate_temperature_spline_interp_m

 !********************************************************************

end module  temperature_spline_interp_m
