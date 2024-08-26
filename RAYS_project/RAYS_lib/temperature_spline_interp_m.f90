module  temperature_spline_interp_m
! Implements temperature profiles based on pointwise data contained in makelist
! temperature_spline_interp_list.

    use constants_m, only : rkind
    use quick_cube_splines_m, only : cube_spline_function_1D

    implicit none

	integer, parameter :: n_grid_max = 200 ! Max dimension of input arrays, truncated later

! namelist /temperature_spline_interp_list/ variables
    character (len = 60) :: spline_temperature_model
	integer :: ngrid ! Actual number of points to be splined
	real(KIND=rkind) ::  ts_in(n_grid_max,0:nspec) ! Values on input for each species


! Stuff for 1D spline profiles

    type(cube_spline_function_1D) :: ts_profile_N_sp ! scalar, not array

 	character (len = 80) :: profile_name = 'Ts_profiles'

  namelist /temperature_spline_interp_list/ ngrid, ts_in

!********************************************************************

contains

!********************************************************************

  subroutine initialize_temperature_spline_interp(read_input, ts_profile_N_sp, is)

    use constants_m, only : one
    use species_m, only : nspec, spec_name
    use diagnostics_m, only : message, message_unit,messages_to_stdout, verbosity
!	use axisym_toroid_eq_m,  only : plasma_psi_limit ! N.B. This is circular dependence but
			! gfortran compiler seems to allow it.  But don't change plasma_psi_limit!
			! We need this to allow plasma outside psiN = 1.

    implicit none

    logical, intent(in) :: read_input
    type(cube_spline_function_1D), intent(out) :: ts_profiles_N_sp

 	integer :: input_unit, get_unit_number ! External, free unit finder

	real(KIND=rkind), allocatable ::  grid(:) ! Grid points (constructed to be 0. to 1.)
	real(KIND=rkind), allocatable ::  ts_values(:) ! Values on grid
	character (len = 80) :: profile_name = 'ne_profile'

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
	allocate(grid(ngrid), ts_values(ngrid))

! Construct psi grid
    do i = 1, ngrid
    	grid(i) = one*(i-1)/(ngrid - 1)
    end do

! Construct arrays for splining
    do i_spec = 0, nspec
		do i = 1, ngrid
			ts_values(i) = ts_in(i, i_spec)/ts_in(i, i_spec) ! Normalized to 1.0 on axis
		end do

		! Initialize spine function for this species
		profile_name = spec_name(i_spec)
		call ts_sp_profiles_N(i_spec)%cube_spline_1D_init(ngrid, grid, ts_values, profile_name)

    end do

! Initialize spline coefficients for ts_profiles_N_sp

    call ts_profiles_N_sp%cube_spline_1D_init(ngrid, grid, ts_values, profile_name)

    return
  end subroutine initialize_temperature_spline_interp

!****************************************************************************************

  subroutine  temperature_spline_interp(psi, dens, dd_psi)

    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message

    implicit none

	real(KIND=rkind), intent(in) :: psi
	real(KIND=rkind), intent(out) :: dens, dd_psi
 write(*,*) 'temperature_spline_interp: got to 1'

    call ts_profiles_N_sp%eval_1D_fp(psi, dens, dd_psi)
 write(*,*) 'temperature_spline_interp: got to 2'

    return
    end subroutine  temperature_spline_interp

!********************************************************************

    subroutine deallocate_temperature_spline_interp_m
		! Nothing to deallocate
		return
    end subroutine deallocate_temperature_spline_interp_m

 !********************************************************************

end module  temperature_spline_interp_m
