 module OX_conv_analysis_m

! Analyzes ray data obtained from an instance of type run_results as defined in module
! module ray_results_m.  The presumption is that the RAYS raun consists of O-mode rays
! that approach the cutoff from low density.  The steps are
! 1) Find the point, x_max, on the ray having maximum density, which should be closest to
!    the cutoff.
! 2) Determine the point, x_cut, on the cutoff surface closest to the x_max point.
! 3) Evaluate the conversion coefficient to X mode
! 4) If don't find density maximum or don't find cutoff surface or conversion coefficient
!    is less than conversion_threshold, ignore this ray and consider that it didn't convert
! 5) Allocate and load array of type OX_conv for the rays that did convert.
!
! This is written for axisymmetric mirror equilibria, but I think it will work little
! change for tokamaks

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use constants_m, only : rkind, zero, one, two, pi
	use diagnostics_m, only : message_unit, text_message, verbosity, run_label, date_v, &
	                  & messages_to_stdout

    implicit none

	integer :: number_of_rays_converted
	real(KIND=rkind), parameter :: conversion_threshold = 0.0001

! derived type containing data for OX_conv
	type OX_conv
		real(KIND=rkind) :: x_max(3) ! point on ray with maximum density, alpha < 1
		real(KIND=rkind) :: k_max(3) ! k vector at x_max
		real(KIND=rkind) :: alpha_max ! omega_pe**2/Omga**2 at x_max
		real(KIND=rkind) :: x_cut(3) ! point on cutoff surface closest to x_max
		real(KIND=rkind) :: conv_coeff ! value of conversion coefficient
		real(KIND=rkind) :: nvecx_c(3) ! n vector at x_max in direction grad(ne)
		real(KIND=rkind) :: nvecy_c(3) ! n vector at x_max transverse to grad(ne), B
		real(KIND=rkind) :: nvecz_c(3) ! n vector at x_max perp to grad(ne) in
		                               ! grad(ne),B plane
		integer :: ray_number ! RAYS ray number for this ray
		integer :: step_number ! ray step number at x_max
	end type OX_conv

	type(OX_conv), allocatable :: OX_conv_data(:)

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

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

	integer :: i_ray ! Number of this ray from RAYS run
	integer :: step_number ! running step number along this ray
	real(KIND=rkind) :: alpha_max ! maximum alpha on this ray
	real(KIND=rkind) :: x_max_ray(3) ! position of maximum alpha on this ray
	real(KIND=rkind) :: k_max_ray(3) ! k value at found_cutoff
	real(KIND=rkind) :: x_cutoff_ray(3) ! position on O-mode cutoff surface nearest x_max_ray
	real(KIND=rkind) :: conv_coeff ! value of conversion coefficient
	real(KIND=rkind) :: nvecx_c(3)
	real(KIND=rkind) :: nvecy_c(3)
	real(KIND=rkind) :: nvecz_c(3)

	logical :: found_max, found_cutoff
	integer :: iteration, i
	integer :: converted_ray_number(number_of_rays)

	type(OX_conv) :: conv_data_temp(number_of_rays)

	number_of_rays_converted = 0
	ray_loop: do i_ray = 1, number_of_rays
		conv_data_temp(i_ray)%ray_number = i_ray

! Find x_max_ray
		call find_x_max_ray(i_ray, found_max, step_number, alpha_max,  x_max_ray, k_max_ray)
		if (found_max) then
			conv_data_temp(i_ray)%x_max = x_max_ray
			conv_data_temp(i_ray)%k_max = k_max_ray
			conv_data_temp(i_ray)%alpha_max = alpha_max
			conv_data_temp(i_ray)%step_number = step_number
		else
			conv_data_temp(i_ray)%x_max = zero
			conv_data_temp(i_ray)%k_max = zero
			conv_data_temp(i_ray)%alpha_max = zero
			conv_data_temp(i_ray)%step_number = 0
		end if

		write(*,*)
		write(*,*) 'ray ', i_ray,  '   found_max = ', found_max, '  alpha_max = ',&
		        &   alpha_max, '   step_number = ', step_number
		write(*,*) 'x_max = ', x_max_ray, '   k_max = ', k_max_ray

! Find nearest point on O-mode cutoff surface
		if (found_max) then
			call find_x_cutoff_ray(x_max_ray, found_cutoff, x_cutoff_ray, iteration)
			if(found_cutoff) then
				conv_data_temp(i_ray)%x_cut = x_cutoff_ray
			else
				conv_data_temp(i_ray)%x_cut = zero
			end if

			write(*,*) 'found_cutoff = ', found_cutoff,'  iteration = ', iteration,&
			         & '  x_cutoff = ', x_cutoff_ray
		end if

! Calculate conversion coefficient to X-mode
		if (found_max .and. found_cutoff ) then
			call OX_conv_coeff(x_max_ray, k_max_ray, x_cutoff_ray, conv_coeff, &
				& nvecx_c, nvecy_c, nvecz_c)
			if (conv_coeff > conversion_threshold) then
				number_of_rays_converted = number_of_rays_converted + 1
				converted_ray_number(number_of_rays_converted) = i_ray
				conv_data_temp(i_ray)%conv_coeff = conv_coeff
				conv_data_temp(i_ray)%nvecx_c = nvecx_c
				conv_data_temp(i_ray)%nvecy_c = nvecy_c
				conv_data_temp(i_ray)%nvecz_c = nvecz_c
				write(*,*) 'conv_coeff = ', conv_coeff
			end if
		end if

	end do ray_loop

! Now that we know how many rays actually converted, allocate OX_conv_data and load from
! conv_data_temp

 	allocate(OX_conv_data(number_of_rays_converted))
	do i = 1, number_of_rays_converted

		i_ray = converted_ray_number(i)

		OX_conv_data(i)%x_max = conv_data_temp(i_ray)%x_max
		OX_conv_data(i)%k_max = conv_data_temp(i_ray)%k_max
		OX_conv_data(i)%alpha_max = conv_data_temp(i_ray)%alpha_max
		OX_conv_data(i)%x_cut = conv_data_temp(i_ray)%x_cut
		OX_conv_data(i)%conv_coeff = conv_data_temp(i_ray)%conv_coeff
		OX_conv_data(i)%nvecx_c = conv_data_temp(i_ray)%nvecx_c
		OX_conv_data(i)%nvecy_c = conv_data_temp(i_ray)%nvecy_c
		OX_conv_data(i)%nvecz_c = conv_data_temp(i_ray)%nvecz_c
		OX_conv_data(i)%ray_number = conv_data_temp(i_ray)%ray_number
		OX_conv_data(i)%step_number = conv_data_temp(i_ray)%step_number

	end do

	if (number_of_rays_converted > 0 ) then
		call write_OX_conversion_data_LD
	end if

	return
	end subroutine analyze_OX_conv

!****************************************************************************

 subroutine find_x_max_ray(i_ray, found_max, step_number, alpha_max, &
                         & x_max_ray, k_max_ray)
! Find point on ray with maximum density, i.e. maximum alpha = omega_pe**2 / omega**2
! Some refinements to consider if turns out to be needed.  Instead of taking ray point
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
	step_number = 1

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
 end subroutine find_x_max_ray

!****************************************************************************

 subroutine find_x_cutoff_ray(x_max_ray, found_cutoff, x_cutoff_ray, iteration)
! Find point on O-mode cutoff surface closest to x_max by steepest ascent

    use equilibrium_m, only : equilibrium, eq_point

	implicit none

	real(KIND=rkind), intent(in) :: x_max_ray(3)
	real(KIND=rkind), intent(out) :: x_cutoff_ray(3)
	logical, intent(out) :: found_cutoff
	integer, intent(out) :: iteration

	real(KIND=rkind) :: alpha_temp
	real(KIND=rkind) :: grad_alpha_unit(3)
	real(KIND=rkind) :: mod_grad_alpha
	real(KIND=rkind) :: delta_x(3)
	real(KIND=rkind) :: x_temp(3)
	real(KIND=rkind) :: delta, r
	real(KIND=rkind), parameter :: alpha_tolerence = one/10.0_rkind**4
	integer, parameter :: max_iterations = 10

	type(eq_point) :: eq

	found_cutoff = .false.
	iteration = 0
	x_temp = x_max_ray
	call equilibrium(x_temp, eq)
	alpha_temp = eq%alpha(0)
! 	write(*,*) ''
! 	write(*,*) 'iteration = ', iteration
! 	write(*,*) 'x_temp = ', x_temp,  '  alpha_temp = ', alpha_temp
	do iteration = 1, max_iterations

		if (abs(alpha_temp - one) <= alpha_tolerence) then
			found_cutoff = .true.
			exit
		end if

		call equilibrium(x_temp, eq)
		alpha_temp = eq%alpha(0)
		mod_grad_alpha = norm2(eq%gradns(:, 0))*(alpha_temp/eq%ns(0))
		grad_alpha_unit = eq%gradns(:, 0)/norm2(eq%gradns(:, 0))
		r = norm2((/x_temp(1), x_temp(2)/))
		delta = (one-eq%alpha(0))/mod_grad_alpha
		! Don't take step bigger than 25% of radius to avoid stepping over axis
		delta_x = grad_alpha_unit*minval((/delta, 0.25_rkind*r/))
		x_temp = x_temp+delta_x
! 	write(*,*) ''
! 	write(*,*) 'iteration = ', iteration
! 	write(*,*) 'alpha_temp = ', alpha_temp, '   mod_grad_alpha = ', mod_grad_alpha, '   grad_alpha_unit = ', grad_alpha_unit
! 	write(*,*) 'delta_x = ', delta_x,  '   x_temp = ', x_temp
	end do
	x_cutoff_ray = x_temp

	return
 end subroutine find_x_cutoff_ray

!****************************************************************************

 subroutine OX_conv_coeff(x_max_ray, k_max_ray, x_cutoff_ray,&
          & conv_coeff, nvecx_c, nvecy_c, nvecz_c)

! Calculate O-X mode conversion coefficient according to the model of Mjolhus, 1984,
! Eq 19.  I introduce a different coordinate system and translate his notation from
! ionospheric interest to something more convenient for fusion.  In his model the
! direction of density stratification is called z (ionospheric vertical), the magnetic
! field lies in the y-z plane and makes angle alpha with the vertical, The coordinate
! orthogonal to y-z is x.  I use a coordinate system (xc, yc, zc) with xc in the direction
! of grad(ne), the magnetic field is in the zc,xc plane, and the orthogonal direction is yc.
! If the density is constant on magnetic field lines then zc is in the same direction as B.
! Generally we expect for fusion applications that B will be nearly orthogonal to grad(ne)
! then B will make a small angle to zc.  grad(ne) is evaluated at x_cutoff_ray.  Also, I
! call the angle between grad(ne) and B theta instead of alpha to avoid confusion with
! my alpha = omega_pe**2/omega_rf**2.  Also Mjolhus uses Y for omega_ce/omega at cutoff,
! whereas in RAYS equilibrium_m this is abs(gamma(0))

!

    use equilibrium_m, only : equilibrium, eq_point
    use rf_m, only : k0
    use vectors3_m, only : cross_product

	implicit none

	real(KIND=rkind), intent(in) :: x_max_ray(3)
	real(KIND=rkind), intent(in) :: k_max_ray(3)
	real(KIND=rkind), intent(in) :: x_cutoff_ray(3)
	real(KIND=rkind), intent(out) :: conv_coeff
	real(KIND=rkind), intent(out) :: nvecx_c(3)
	real(KIND=rkind), intent(out) :: nvecy_c(3)
	real(KIND=rkind), intent(out) :: nvecz_c(3)

	real(KIND=rkind) :: xc_unit(3), yc_unit(3), zc_unit(3)
	real(KIND=rkind) :: v_temp(3)
	real(KIND=rkind) :: theta, gamma, F, G
	real(KIND=rkind) :: n_parallel, ny_c, nz_c, n_vertical, n_crit
	real(KIND=rkind) :: L

	type(eq_point) :: eq

! N.B. These things are evaluated at the cutoff surface
	call equilibrium(x_cutoff_ray, eq)
	xc_unit = eq%gradns(:,0)/norm2(eq%gradns(:, 0)) ! Unit vector along grad(ne)
	v_temp = cross_product(eq%bunit, xc_unit) ! perpendicular to x and B, i.e. y direction
	yc_unit = v_temp/norm2(v_temp)
	zc_unit = cross_product(xc_unit, yc_unit)
	theta = acos(dot_product(xc_unit, eq%bunit))
	gamma = abs(eq%gamma(0))

 write(*,*) 'xc_unit = ', xc_unit
 write(*,*) 'yc_unit = ', yc_unit
 write(*,*) 'zc_unit = ', zc_unit
 write(*,*) 'theta = ', theta
 write(*,*) 'gamma = ', gamma

! N.B. k vector is evaluated on the ray at x_max_ray.  If the plasma were plane stratified,
! as in Mjohus model, n_parallel and n_transverse would be constant along the ray.  I'm
! evaluating them using bunit at the reflection point.  There might be some ambiguity
! relative to evaluating bunit at the cutoff surface, but at least they do satisfy the
! dispersion relation where they are evaluated.
	call equilibrium(x_max_ray, eq)
	n_parallel = dot_product(k_max_ray, eq%bunit)/k0
	nz_c = dot_product(k_max_ray, zc_unit)/k0
	ny_c = dot_product(k_max_ray, yc_unit)/k0
	n_vertical = dot_product(k_max_ray, xc_unit)/k0

	nvecx_c = n_vertical*xc_unit
	nvecy_c = ny_c*yc_unit
	nvecz_c = nz_c*zc_unit

	n_crit = sin(theta)*sqrt(gamma/(one + gamma))

	F = 0.5_rkind*(one+gamma)*sqrt(gamma)/((one+gamma)*cos(theta)**2 + &
	 &  sin(theta)**2/two)**1.5_rkind

	G = 0.5_rkind*sqrt(gamma)/sqrt((one+gamma)*cos(theta)**2 + sin(theta)**2/two)

	L = eq%ns(0)/norm2(eq%gradns(:,0))

    conv_coeff = exp(-pi*k0*L*(F*(abs(nz_c)-n_crit)**2 + G*abs(ny_c)**2))

 write(*,*) 'n_parallel = ', n_parallel
 write(*,*) 'ny_c = ', ny_c
 write(*,*) 'nz_c = ', nz_c
 write(*,*) 'n_vertical = ', n_vertical
 write(*,*) 'n_crit = ', n_crit
 write(*,*) 'F = ', F
 write(*,*) 'G = ', G
 write(*,*) 'L = ', L
 write(*,*) 'conv_coeff = ', conv_coeff

	return
 end subroutine OX_conv_coeff

!****************************************************************************

    subroutine write_OX_conversion_data_LD

		use diagnostics_m, only : run_label

		implicit none

		integer :: OX_conv_unit, get_unit_number
		integer :: i_ray

		!  File name for message output
		character(len=80) :: out_file

		! Open file to put deposition data in
! 		OX_conv_unit = get_unit_number()
! 		out_file = 'OX_conv_data.'//trim(run_label)
! 		open(unit=OX_conv_unit, file=out_file, action='write', status='replace',&
! 		     & form='formatted')

		do i_ray = 1, number_of_rays_converted
			write(*,*)'ray ', i_ray, 'conv coeff = ', OX_conv_data(i_ray)%conv_coeff, 'nz_c = ', &
			 & norm2(OX_conv_data(i_ray)%nvecz_c), 'ny_c = ', &
			 & norm2(OX_conv_data(i_ray)%nvecy_c)
		end do

!         close(unit = OX_conv_unit)

	return
    end subroutine write_OX_conversion_data_LD


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
