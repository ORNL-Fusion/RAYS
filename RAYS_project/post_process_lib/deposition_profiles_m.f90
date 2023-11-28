 module deposition_profiles_m

! Contains data for various deposition profiles (such as power versus psi summed
! over species -> Ptotal(:)) and associated routines to calculate the profiles from ray
! data in module ray_results_m.  Profile data is put into an instance of derived type
! "dep_profile".  Generating each profile requires a derived type "work_space".
!
! The profile data and the grid on which the profile is defined must be calculable from
! the ray data in ray_results_m.  For example to calculate the total power deposition
! profile on a psi grid one has at each step of each ray the spatial position, X, in the
! ray vector, v(1:3), from, which psi(X) can be calculated.  Total power deposited is
! obtained from ray vector component v(8). The calculation is done in a generic subroutine
! pointer, Q_evaluator(iray, ix, xQ, Q), where iray = ray number, ix = step number along
! the ray, and returned are xQ = grid value a point X and Q = profile quantity.
!
! Which profiles make sense depends on the equilibrium geometry.  Presently supported are:
! "slab" and "axisym_toroid" geometries.  For now assume that all profiles will be calculated
! that we know how to calculate (although there are ideas on how to let the user can select
! which profiles to calculate.  This is sitting in a previous version, in spare_parts)
!!
! Every profile needs:
! 1) An entry in the list of known profile names and the name of the grid that goes with it
! 2) An evaluator subroutine to calculate the grid and profile quantity at each ray step
! 3) A case selector block in initialize_deposition_profiles() pointing to the evaluator routine

    use constants_m, only : rkind

    implicit none

! Lists of known profiles
    character(len=20), parameter :: slab_profile_names(*) = [character(len=20) :: 'Ptotal_x']
    character(len=20), parameter :: slab_grid_names(*) = [character(len=20) :: 'x']
    character(len=20), parameter :: axisym_toroid_profile_names(*) = &
!                         & [character(len=20) :: 'Ptotal_psi','Ptotal_rho']
                         & [character(len=20) :: 'Ptotal_psi']
    character(len=20), parameter :: axisym_toroid_grid_names(*) = [character(len=20) ::  &
                         & 'psi', 'rho']

    logical :: write_results_list_directed = .true.
    integer :: n_bins = 0, default_n_bins = 100
    integer :: n_profiles

!***********************************

    abstract interface
		subroutine Q_evaluator(iray, ix, xQ, Q)
			use constants_m, only : rkind
			use ray_results_m, only : ray_vec
			implicit none
			integer, intent(in) :: iray, ix !ray number, step number along the ray
			real(KIND=rkind), intent(out) :: xQ, Q ! grid value at point X, profile quantity
		end subroutine Q_evaluator
    end interface

!****************************************************************************


	type dep_profile
		character(len=20) :: profile_name, grid_name
		integer :: n_bins
        real(KIND=rkind) :: grid_min, grid_max
        real(KIND=rkind), allocatable :: grid(:), profile(:)
        real(KIND=rkind), allocatable :: work(:,:)
        real(KIND=rkind) :: Q_sum
        procedure(Q_evaluator), pointer, nopass :: evaluator => Null()
    end type dep_profile

! Array of 1D profiles to calculate
	type(dep_profile), allocatable :: profiles_1D(:)

    namelist /deposition_profiles_list/ n_bins, write_results_list_directed

!****************************************************************************

contains

!****************************************************************************

subroutine initialize_deposition_profiles(read_input)

	use constants_m, only : rkind
	use diagnostics_m, only : message_unit, text_message, verbosity, run_label, date_v
	use ray_init_m, only : nray
	use equilibrium_m, only : equilib_model
	use slab_eq_m, only : xmin_slab => xmin, xmax_slab => xmax

	implicit none
	logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder

! Data for all profiles
	character(len=20) :: profile_name, grid_name
	real(KIND=rkind) :: grid_min, grid_max

	integer :: i, i_profile
    real(KIND=rkind) :: delta


	call text_message('initialize_deposition_profiles', 1)

	if (read_input .eqv. .true.) then
	! Read and write input namelist
  		input_unit = get_unit_number()
		open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
		read(input_unit, deposition_profiles_list)
		close(unit=input_unit)
		if (verbosity > 0) write(message_unit, deposition_profiles_list)
	end if

	if (n_bins == 0) n_bins = default_n_bins

    ! find what geometry we are in
    equilibria: select case (trim(equilib_model))
       case ('slab')
!         A 1-D slab equilibrium with stratification in x.

          n_profiles = size(slab_profile_names)
		  allocate(profiles_1D(n_profiles))

  		  grid_min = xmin_slab
		  grid_max = xmax_slab

		  ! Initialize everything inside dep_profile objects
		  do i_profile = 1, n_profiles

				profile_name = slab_profile_names(i_profile)
				grid_name = slab_grid_names(i_profile)
				profiles_1D(i_profile)%profile_name  = profile_name
				profiles_1D(i_profile)%grid_name = grid_name
				profiles_1D(i_profile)%n_bins = n_bins
				profiles_1D(i_profile)%grid_min = grid_min
				profiles_1D(i_profile)%grid_max = grid_max

				! Allocate arrays inside dep_profile objects
				allocate(profiles_1D(i_profile)%grid(n_bins+1), source = 0.0_rkind)
				allocate(profiles_1D(i_profile)%profile(n_bins), source = 0.0_rkind)
				allocate(profiles_1D(i_profile)%work(n_bins, nray), source = 0.0_rkind)

				! Generate uniform grid
				delta = (grid_max - grid_min)/real(n_bins)
				do i = 1, n_bins+1
					profiles_1D(i_profile)%grid(i) = grid_min + delta*(i-1)
				end do

				! Set pointer to proper evaluator
				select case (trim(profile_name))
					case('Ptotal_x')
						profiles_1D(i_profile)%evaluator => Ptotal_x_slab_evaluator
					case default
					  call text_message('initialize_deposition_profiles: unimplemented slab profile',&
					  & trim(profile_name), 0)
					  stop 1
				end select

		  end do

       case ('axisym_toroid')
!         A generic axisymmetric toroidal plasma model.

		  grid_min = 0.0_rkind
		  grid_max = 1.0_rkind

          n_profiles = size(axisym_toroid_profile_names)
 		  allocate(profiles_1D(n_profiles))

		  ! Initialize everything inside dep_profile objects
		  do i_profile = 1, n_profiles

				profile_name = axisym_toroid_profile_names(i_profile)
				profiles_1D(i_profile)%profile_name  = profile_name
				profiles_1D(i_profile)%grid_name  = axisym_toroid_grid_names(i_profile)
				profiles_1D(i_profile)%n_bins = n_bins
				profiles_1D(i_profile)%grid_min = grid_min
				profiles_1D(i_profile)%grid_max = grid_max

				! Allocate arrays inside dep_profile objects
				allocate(profiles_1D(i_profile)%grid(n_bins+1), source = 0.0_rkind)
				allocate(profiles_1D(i_profile)%profile(n_bins), source = 0.0_rkind)
				allocate(profiles_1D(i_profile)%work(n_bins, nray), source = 0.0_rkind)

				! Generate uniform grid, Grid range [0.0, 1.0], either psi or rho
				delta = (grid_max - grid_min)/real(n_bins)
				do i = 1, n_bins+1
					profiles_1D(i_profile)%grid(i) = grid_min + delta*(i-1)
				end do

				! Set pointer to proper evaluator
				select case (trim(profile_name))
					case('Ptotal_psi')
						profiles_1D(i_profile)%evaluator => Ptotal_axisym_psi_evaluator
! 					case('Ptotal_rho')
! 						profiles_1D(i_profile)%evaluator => Ptotal_x_slab_evaluator
					case default
					  call text_message(&
					   &'initialize_deposition_profiles: unimplemented axisym_toroid profile',&
					  & trim(profile_name), 0)
					  stop 1
				end select

		  end do

       case default
          call text_message('initialize_deposition_profiles: unimplemented equilib_model',&
          & trim(equilib_model), 0)
          stop 1
    end select equilibria

	return
	end subroutine initialize_deposition_profiles

!****************************************************************************

    subroutine calculate_deposition_profiles

		use constants_m, only : rkind
	    use diagnostics_m, only : text_message
		use ray_results_m
		use ray_init_m, only : nray
		use equilibrium_m, only : equilib_model

		implicit none

		integer :: i_profile, iray

		call text_message('Calculating deposition profiles ', 1)



		profile_loop: do i_profile = 1, n_profiles

			ray_loop: do iray = 1, nray
				call bin_a_ray(profiles_1D(i_profile), iray)
			end do ray_loop

			! sum over rays to get total profile, sum over profile to get total quantity deposited
			profiles_1D(i_profile)%profile(:) = sum(profiles_1D(i_profile)%work, 2)
			profiles_1D(i_profile)%Q_sum = sum(profiles_1D(i_profile)%profile)

			! Deallocate the work array after summing the profile over rays
			deallocate(profiles_1D(i_profile)%work)

		end do profile_loop

	return
    end subroutine calculate_deposition_profiles


!****************************************************************************

    subroutine bin_a_ray(prof, iray)

		use constants_m, only : rkind
		use ode_m, only : nstep_max
		use ray_results_m, only : npoints
		use ray_init_m, only : nray
   	    use bin_to_uniform_grid_m, only : bin_to_uniform_grid

		implicit none

        type(dep_profile) :: prof
		integer, intent(in) :: iray

		real(KIND=rkind):: Q(nstep_max+1)  ! quantity at ray point i
		real(KIND=rkind) :: xQ(nstep_max+1)  ! grid value at ray point i

		integer :: ip, ierr

		! Evaluate grid value and profile quantity at each ray point
		do ip = 1, npoints(iray)
			call prof%evaluator(iray, ip, xQ(ip), Q(ip))
		end do

    	call bin_to_uniform_grid(Q(1:npoints(iray)), xQ(1:npoints(iray)), prof%grid_min,&
    	     & prof%grid_max, prof%work(:, iray), ierr)

	return
    end subroutine bin_a_ray

!****************************************************************************

    subroutine write_deposition_profiles_list_directed

		use diagnostics_m, only : run_label
		use ode_m, only : nstep_max
		use ray_results_m, only : npoints
		use ray_init_m, only : nray

		implicit none

		integer :: dep_profile_unit, get_unit_number
		integer :: i_profile

		!  File name for message output
		character(len=80) :: out_file

		! Open file to put deposition data in
		dep_profile_unit = get_unit_number()
		out_file = 'deposition_profiles.'//trim(run_label)
		open(unit=dep_profile_unit, file=out_file, action='write', status='replace',&
		     & form='formatted')

		profile_loop: do i_profile = 1, n_profiles

			write(dep_profile_unit, *) 'profile_name = ', profiles_1D(i_profile)%profile_name
			write(dep_profile_unit, *) profiles_1D(i_profile)%profile
			write(dep_profile_unit, *) 'grid_name = ', profiles_1D(i_profile)%grid_name
			write(dep_profile_unit, *) profiles_1D(i_profile)%grid
			write(dep_profile_unit, *) 'Ptotal_total_deposition'
			write(dep_profile_unit, *) profiles_1D(i_profile)%Q_sum

		end do profile_loop

    close(unit = dep_profile_unit)

	return
    end subroutine write_deposition_profiles_list_directed

!****************************************************************************
! Evaluators
!****************************************************************************

    subroutine Ptotal_x_slab_evaluator(iray, ix, xQ, Q)

		use constants_m, only : rkind
		use ray_results_m, only : ray_vec, initial_ray_power

		implicit none

		integer, intent(in) :: iray, ix
		real(KIND=rkind), intent(out) :: xQ, Q

		xQ = ray_vec(1, ix, iray)
		Q = ray_vec(8, ix, iray)*initial_ray_power(iray)

 	return
    end subroutine Ptotal_x_slab_evaluator

!********************************************************************

    subroutine Ptotal_axisym_psi_evaluator(iray, ix, xQ, Q)

		use constants_m, only : rkind
		use ray_results_m, only : ray_vec, initial_ray_power
    	use axisym_toroid_eq_m, only: axisym_toroid_psi

		implicit none

		integer, intent(in) :: iray, ix
		real(KIND=rkind), intent(out) :: xQ, Q

    	real(KIND=rkind) :: rvec(3)
		real(KIND=rkind) :: psi, gradpsi(3), psiN, gradpsiN(3)
		rvec = ray_vec(1:3, ix, iray)
		call axisym_toroid_psi(rvec, psi, gradpsi, psiN, gradpsiN)

		xQ = psiN
		Q = ray_vec(8, ix, iray)*initial_ray_power(iray)

 	return
    end subroutine Ptotal_axisym_psi_evaluator

!********************************************************************
! Deallocate
!****************************************************************************

    subroutine deallocate_deposition_profiles_m

		implicit none
		integer :: i_profile

		profile_loop: do i_profile = 1, n_profiles
			deallocate(profiles_1D(i_profile)%grid)
			deallocate(profiles_1D(i_profile)%profile)
            ! deallocate(profiles_1D(i_profile)%work) work array already deallocated above
            ! in subroutine calculate_deposition_profiles()
		end do profile_loop

        return
    end subroutine deallocate_deposition_profiles_m

 end module deposition_profiles_m
