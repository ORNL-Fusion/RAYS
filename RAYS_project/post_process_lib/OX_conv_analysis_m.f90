 module OX_conv_analysis_m

! Analyzes ray data obtained from an instance of type run_results as defined in module
! module ray_results_m.  The presumption is that the RAYS raun consists of O-mode rays
! that approach the cutoff from low density.  The steps are
! 1) Find the point, v_max, on the ray having maximum density, which should be closest to
!    the cutoff.
! 2) Determine the point, x_cut, on the cutoff surface closest to the v_max point.
! 3) Evaluate the conversion coefficient to X mode
!
! This is written for axisymmetric mirror equilibria, but I think it will work without
! change for tokamaks

    use constants_m, only : rkind, zero, one
	use diagnostics_m, only : message_unit, text_message, verbosity, run_label, date_v, &
	                  & messages_to_stdout

    implicit none

! derived type containing data for OX_conv
	type OX_conv
		real(KIND=rkind) v_max(3) ! point on ray with maximum density, alpha < 1
		real(KIND=rkind) n_max(3) ! ray refractive index vector at v_max
		real(KIND=rkind) v_cut(3) ! point on cutoff surface closest to v_max
		real(KIND=rkind) conv_coeff ! value of conversion coefficient
	end type OX_conv

!****************************************************************************

!   namelist /OX_conv_list/

!****************************************************************************

contains

!****************************************************************************

! subroutine initialize_OX_conv_analysis(read_input)
!
! 	implicit none
! 	logical, intent(in) :: read_input
! 	integer :: input_unit, get_unit_number ! External, free unit finder
!
! ! Data for all profiles
! 	character(len=20) :: profile_name, grid_name
! 	real(KIND=rkind) :: grid_min, grid_max
!
! 	integer :: i, i_profile
!     real(KIND=rkind) :: delta
!
!
! 	call text_message('initialize_deposition_profiles', 1)
!
! 	if (read_input .eqv. .true.) then
! 	! Read and write input namelist
!   		input_unit = get_unit_number()
! 		open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
! 		read(input_unit, OX_conv_list)
! 		close(unit=input_unit)
! 	end if
!
! ! Write input namelist
!     if (verbosity >= 0) then
! 		write(message_unit, OX_conv_list)
! 		if (messages_to_stdout) write(*, OX_conv_list)
!     end if
!
!
! 	return
! 	end subroutine initialize_OX_conv_analysis

!****************************************************************************

subroutine analyze_OX_conv

	use ray_results_m, only : number_of_rays

	implicit none

	integer :: i_ray
	integer :: step_number
	real(KIND=rkind) :: alpha_max
	real(KIND=rkind) :: v_max_ray(3)
	real(KIND=rkind) :: k_max_ray(3)
	logical :: found_max

	type(OX_conv) :: conv_data_temp(number_of_rays)
	type(OX_conv), allocatable :: OX_conv_data(:)

	ray_loop: do i_ray = 1, number_of_rays

		call find_v_max_ray(i_ray, found_max, step_number, alpha_max,  v_max_ray, k_max_ray)
		conv_data_temp(i_ray)%v_max = v_max_ray
		conv_data_temp(i_ray)%n_max = k_max_ray

		write(*,*)
		write(*,*) 'ray ', i_ray,  '   found_max = ', found_max, '  alpha_max = ',&
		        &   alpha_max, '   step_number = ', step_number
		write(*,*) 'v_max = ', v_max_ray, '   k_max = ', k_max_ray

	end do ray_loop

	return
	end subroutine analyze_OX_conv

!****************************************************************************

 subroutine find_v_max_ray(i_ray, found_max, step_number, alpha_max, &
                         & x_max_ray, k_max_ray)
! Find point on ray with maximum density, i.e. maximum alpha = omega_pe**2 / omega**2
! Some refinements to consider if turns out to be needed.  Instead of taking ray pointer
! with highest ne, fit quadratic to last 3 points and interpolate to max.  Also could
! not start searching for maximum until alpha is close to 1 so as not to get fooled by
! a local maximum far from the cutoff.

	use ray_results_m, only : npoints, ray_vec
    use equilibrium_m, only : equilibrium, eq_point

	implicit none

	integer, intent(in) :: i_ray
	integer, intent(out) :: step_number
	real(KIND=rkind), intent(out) :: alpha_max
	real(KIND=rkind), intent(out) :: x_max_ray(3)
	real(KIND=rkind), intent(out) :: k_max_ray(3)
	logical, intent(out) :: found_max

	integer :: i
	real(KIND=rkind) :: alpha_low, alpha_high
	real(KIND=rkind) :: rvec(3)

	type(eq_point) :: eq

	found_max = .false.
	alpha_max = zero
	step_number = one

	rvec = ray_vec(1:3,1, i_RAY)
	call equilibrium(rvec, eq)
	alpha_low = eq%alpha(0)

	do i = 2, npoints(i_ray)
		rvec = ray_vec(1:3,i, i_RAY)
		call equilibrium(rvec, eq)
		alpha_high = eq%alpha(0)
		if (alpha_high < alpha_low) then
			found_max = .true.
			alpha_max = alpha_low
			step_number = i-1
			x_max_ray = ray_vec(1:3,i-1, i_RAY)
			k_max_ray = ray_vec(4:6,i-1, i_RAY)
			exit
		end if
		alpha_low = alpha_high
	end do

	return
 end subroutine find_v_max_ray

!****************************************************************************

!     subroutine write_deposition_profiles_LD
!
! 		use diagnostics_m, only : run_label
! 		use ode_m, only : nstep_max
! 		use ray_results_m, only : npoints
! 		use ray_init_m, only : nray
!
! 		implicit none
!
! 		integer :: dep_profile_unit, get_unit_number
! 		integer :: i_profile
!
! 		!  File name for message output
! 		character(len=80) :: out_file
!
! 		! Open file to put deposition data in
! 		dep_profile_unit = get_unit_number()
! 		out_file = 'deposition_profiles.'//trim(run_label)
! 		open(unit=dep_profile_unit, file=out_file, action='write', status='replace',&
! 		     & form='formatted')
!
! 		profile_loop: do i_profile = 1, n_profiles
!
! 			write(dep_profile_unit, *) 'profile_name = ', profiles_1D(i_profile)%profile_name
! 			write(dep_profile_unit, *) profiles_1D(i_profile)%profile
! 			write(dep_profile_unit, *) 'grid_name = ', profiles_1D(i_profile)%grid_name
! 			write(dep_profile_unit, *) profiles_1D(i_profile)%grid
! 			write(dep_profile_unit, *) 'Ptotal_total_deposition'
! 			write(dep_profile_unit, *) profiles_1D(i_profile)%Q_sum
!
! 		end do profile_loop
!  write(*,*) 'Q_sum = ', profiles_1D(1)%Q_sum
!     close(unit = dep_profile_unit)
!
! 	return
!     end subroutine write_deposition_profiles_LD
!
!
! !****************************************************************************
!
!  subroutine write_deposition_profiles_NC
!
! 	use diagnostics_m, only : date_v, run_label
! 	use ode_m, only : nstep_max
! 	use ray_results_m, only : npoints
! 	use ray_init_m, only : nray
! 	use netcdf
!
! 	implicit none
!
! 	integer :: dep_profile_unit, get_unit_number
! 	integer :: i_profile
!
! 	!  File name for output
! 	character(len=80) :: out_filename
!
! ! netCDF Declarations
! 	integer :: ncid
! 	integer :: nf90_put_var, dim_len, ierr
! 	integer :: start1(1)
! 	character(len = 20) :: dim_name
! ! Declarations: dimensions
! 	integer, parameter :: n_dims = 5, d8 = 8, d20 = 20
! !	integer :: n_profiles, n_bins, n_bins_p1
! 	integer :: n_bins_p1
! 	integer :: n_profiles_id, n_bins_id, n_bins_p1_id, d8_id, d20_id
! ! Declarations: variables
! 	integer, parameter :: n_vars =  8
! 	integer :: Q_sum_id, date_vector_id, profile_name_id, grid_name_id, grid_min_id,&
! 		         & grid_max_id, grid_id, profile_id
!
! !   Open NC file
! 	dep_profile_unit = get_unit_number()
! 	out_filename = 'deposition_profiles.'//trim(run_label)//'.nc'
! 	call check( nf90_create(trim(out_filename), nf90_clobber, ncid) )
!
! !   Define NC dimensions
!     call check( nf90_def_dim(ncid, 'n_profiles', NF90_UNLIMITED, n_profiles_id))
!     call check( nf90_def_dim(ncid, 'n_bins', n_bins, n_bins_id))
!     call check( nf90_def_dim(ncid, 'n_bins_p1', n_bins+1, n_bins_p1_id))
!     call check( nf90_def_dim(ncid, 'd20', d20, d20_id))
!
! ! Define NC variables
!     call check( nf90_def_var(ncid, 'Q_sum', NF90_DOUBLE, n_profiles_id, Q_sum_id))
!     call check( nf90_def_var(ncid, 'n_bins', NF90_INT, [n_profiles_id], n_bins_id))
!     call check( nf90_def_var(ncid, 'grid_min', NF90_DOUBLE, [n_profiles_id], grid_min_id))
!     call check( nf90_def_var(ncid, 'grid_max', NF90_DOUBLE, [n_profiles_id], grid_max_id))
!     call check( nf90_def_var(ncid, 'profile_name', NF90_CHAR, [d20_id,n_profiles_id], profile_name_id))
!     call check( nf90_def_var(ncid, 'grid_name', NF90_CHAR, [d20_id,n_profiles_id], grid_name_id))
!     call check( nf90_def_var(ncid, 'grid', NF90_DOUBLE, [n_bins_p1_id,n_profiles_id], grid_id))
!     call check( nf90_def_var(ncid, 'profile', NF90_DOUBLE, [n_bins_id,n_profiles_id], profile_id))
!
! ! Put global attributes
!     call check( nf90_put_att(ncid, NF90_GLOBAL, 'RAYS_run_label', run_label))
!     call check( nf90_put_att(ncid, NF90_GLOBAL, 'date_vector', date_v))
!
!     call check( nf90_enddef(ncid))
!
! ! Put NC variables
! 	profile_loop: do i_profile = 1, n_profiles
!
! 		call check( nf90_put_var(ncid, Q_sum_id, profiles_1D(i_profile)%Q_sum,&
! 					start=[i_profile]))
! 		call check( nf90_put_var(ncid, n_bins_id, profiles_1D(i_profile)%n_bins,&
! 					start=[i_profile]))
! 		call check( nf90_put_var(ncid, grid_min_id, profiles_1D(i_profile)%grid_min,&
! 					start=[i_profile]))
!  		call check( nf90_put_var(ncid, grid_max_id, profiles_1D(i_profile)%grid_max,&
! 					start=[i_profile]))
! 		call check( nf90_put_var(ncid, profile_name_id, profiles_1D(i_profile)%profile_name,&
! 					start=[1,i_profile], count = [d20,1]))
!  		call check( nf90_put_var(ncid, grid_name_id, profiles_1D(i_profile)%grid_name,&
! 					start=[1,i_profile], count = [d20,1]))
! 		call check( nf90_put_var(ncid, grid_id, profiles_1D(i_profile)%grid,&
! 					start=[1,i_profile], count = [n_bins+1,1]))
! 		call check( nf90_put_var(ncid, profile_id, profiles_1D(i_profile)%profile,&
! 					start=[1,i_profile], count = [n_bins,1]))
!
! 	end do profile_loop
!
! !   Close the NC file
!     call check( nf90_close(ncid) )
!
! 	return
!     end subroutine write_deposition_profiles_NC
!
! !****************************************************************************
!
!   subroutine check(status)
!     use netcdf
!     integer, intent ( in) :: status
!
!     if(status /= nf90_noerr) then
!       print *, trim(nf90_strerror(status))
!       stop 2
!     end if
!   end subroutine check
!
!
!********************************************************************
! Deallocate
!****************************************************************************

    subroutine deallocate_OX_conv_analysis_m

! Nothing to deallocate

        return
    end subroutine deallocate_OX_conv_analysis_m

 end module OX_conv_analysis_m
