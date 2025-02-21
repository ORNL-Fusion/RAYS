module  density_spline_interp_m
! Implements density profiles based on pointwise data contained in makelist
! density_spline_interp_list.
! At present (8/12/2024) this only produces only one profile (presumably electrons),
! normalized to one.  This is similar to the "parabolic" profiles implemented in the
! body of axisym_toroid_eq_m. The ion species then track the electron density with
! fraction of the electron density given by n0s(0:nspec).  In the future, when ion dynamics
! counts (e.g. ICRF) and separate ion density data is provided, something more elaborate
! will be necessary.  In particular how to ensure charge neutrality and specify Zeff.
! The density values in the namelist are assumed to be on a uniform grid from 0.0 to 1.0
! For psi > 1.0 returns dens = d_scrape_off, dd_psi = 0.
!
! For now only one spline_density_model and one function is exported, spline_density_ne.
! So, for now, namelist variable 'spline_density_model' is not used.
!

    use constants_m, only : rkind
    use quick_cube_splines_m, only : cube_spline_function_1D

    implicit none

	integer, parameter :: n_grid_max = 200 ! Max dimension of input arrays, truncated later

! namelist /density_spline_interp_list/ variables
    character (len = 60) :: spline_density_model ! Not used, YET
	integer :: ngrid ! Actual number of points to be splined
	real(KIND=rkind) ::  ne_in(n_grid_max) ! Values on grid (grid assumed uniform 0 to 1)


! Stuff for 1D spline profiles

    type(cube_spline_function_1D) :: ne_profile_N  ! ne profile normalized to 1. on axis
 	character (len = 80) :: profile_name = 'ne_profile'

  namelist /density_spline_interp_list/ ngrid, ne_in

!********************************************************************

contains

!********************************************************************

!  subroutine initialize_density_spline_interp(read_input, ne_profile_N)
  subroutine initialize_density_spline_interp(read_input)

    use constants_m, only : one
    use species_m, only : nspec
    use diagnostics_m, only : message, message_unit,messages_to_stdout, verbosity
!	use axisym_toroid_eq_m,  only : plasma_psi_limit ! N.B. This is circular dependence but
			! gfortran compiler seems to allow it.  But don't change plasma_psi_limit!
			! We need this to allow plasma outside psiN = 1.

    implicit none

    logical, intent(in) :: read_input
!    type(cube_spline_function_1D), intent(out) :: ne_profile_N

 	integer :: input_unit, get_unit_number ! External, free unit finder

	real(KIND=rkind), allocatable ::  grid(:) ! Grid points (constructed to be 0. to 1.)
	real(KIND=rkind), allocatable ::  ne_values(:) ! Values on grid
	character (len = 80) :: profile_name = 'ne_profile'

	integer :: i

    if (verbosity >= 0) then
		write(*,*) 'initialize_density_spline_interp'
    end if

    if (read_input .eqv. .true.) then
  		input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, density_spline_interp_list)
        close(unit=input_unit)
    end if

! Write input namelist
    if (verbosity >= 0) then
		write(message_unit, density_spline_interp_list)
		if (messages_to_stdout) write(*, density_spline_interp_list)
		call message(1)
    end if

! Allocate arrays
	allocate(grid(ngrid), ne_values(ngrid))

! Construct arrays for splining
    do i = 1, ngrid
!    	grid(i) = plasma_psi_limit*(i-1)/(ngrid - 1)
    	grid(i) = one*(i-1)/(ngrid - 1)
    	ne_values(i) = ne_in(i)/ne_in(1) ! Normalized to 1.0 on axis
    end do

! Initialize spline coefficients for ne_profile_N

    call ne_profile_N%cube_spline_1D_init(ngrid, grid, ne_values, profile_name)
    return
  end subroutine initialize_density_spline_interp

!****************************************************************************************

  subroutine  density_spline_interp(psi, d_scrape_off, dens, dd_psi)

    use constants_m, only : one, zero
    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message

    implicit none

	real(KIND=rkind), intent(in) :: psi, d_scrape_off
	real(KIND=rkind), intent(out) :: dens, dd_psi

    if (psi <= one) then
  	call ne_profile_N%eval_1D_fp(psi, dens, dd_psi)
    else
    	dens = d_scrape_off
    	dd_psi = zero
    end if

    return
    end subroutine  density_spline_interp

!********************************************************************

    subroutine deallocate_density_spline_interp_m
		! Nothing to deallocate
		return
    end subroutine deallocate_density_spline_interp_m

 !********************************************************************

end module  density_spline_interp_m
