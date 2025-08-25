 module mirror_processor_m
! Post processing for mirror equilibrium
!
! Subroutine "mirror_processor" calls a number of other routines that write write data for
! further analysis or plotting:

! "write_graphics_description_file" Writes data to be read by python plotting programs

! "ray_detailed_diagnostics" Writes a netCDF file containing diagnostic data collected all
! along the rays

! "write_eq_contour_data_NC" Writes a netcdf file with data needed to plot
! contours of quantities like gamma = cyclotron freq/wave freq.  The contours can then be
! put on plots of ray trajectories

! "write_eq_XZ_grid_data_NC" Writes a netcdf file with data needed to plot
! contours of equilibrium quantites like AphiN, ne, Te, etc
! evaluated on a uniformly spaced  R,Z grid [box_rmin==0, box_rmax, box_zmin, box_zmax]

! "write_eq_radial_profile_data_NC"  Writes a netcdf file containing data needed to plot
! radial profiles for equilibrium quantites like ne, Te, etc as functions of
! AphiN and also R.  These profiles are evaluated at Z = z_reference, which is an input
! parameter in namelist group /mirror_processor_list/.  Profiles extend radially from zero
! to plasma_AphiN_limit, which is set as a namelist input to multiple_mirror_eq_m.

! N.B. This code is adapted, with very few changes from axisym_toroid_processor.f90


    use constants_m, only : rkind

    implicit none

! Number of k vectors to plot for each ray in graphics
    integer :: num_plot_k_vectors

! Scale plot k vectors to kmax on the ray, (True, False)
    character(len = 5) :: scale_k_vec = 'True'

! Base length for drawing k vectors as fraction of fig size
    real(KIND=rkind) :: k_vec_base_length

! Set plot xlim -> [xmin, xmax] and ylim -> [ymin,ymax], (True, False)
    character(len = 5) :: set_XY_lim = 'True'

! Number of plasma boundary points to calculate
	integer :: n_boundary_points
! R,Z boundary points
    real(KIND=rkind), allocatable :: R_boundary(:), Z_boundary(:)

! Number of points in X,Z grid for normalized Aphi and eq_XZ_grid netCDF files
	integer ::  N_pointsX_eq, N_pointsZ_eq

! Number of points in AphiN grid for 1D profiles
	integer ::  n_AphiN

! Number of points in rho grid for 1D profiles
	integer ::  n_rho

! Tolerance for finding R(Aphi) grid by bisection
    real(KIND=rkind) :: bisection_eps

! Reference z location at which to evaluate radial equilibrium profiles
    real(KIND=rkind) :: z_reference
	character(len=25) :: char_z_ref

! Radius of O-mode cutoff at z = z_reference
    real(KIND=rkind) :: r_Omode_cut_at_z_ref

! Radius of LUFS at z = z_reference
    real(KIND=rkind) :: r_LUFS_at_z_ref

! Flags determining what to do besides plot rays.  Can be overridden (i.e. turned off)
! in namelist file
	logical :: calculate_dep_profiles = .true. ! Calculate all depositon profiles
	logical :: write_dep_profiles = .true. ! Write depositon profile netCDF file
	logical :: calculate_ray_diag = .true. ! Write detailed ray diagnostic netCDF file
	logical :: write_contour_data = .true.  ! Write data needed to plot contours netCDF file
	logical :: write_eq_XZ_grid_data = .true.  ! Write data for equilibrium on RZ grid netCDF file
	logical :: write_eq_radial_profile_data = .true.  ! Write data for radial profiles netCDF file
	logical :: do_OX_conv_analysis = .false.  ! Special case OX conversion so, default = false

    namelist /mirror_processor_list/ num_plot_k_vectors, scale_k_vec,&
             & k_vec_base_length, set_XY_lim, &
             & calculate_dep_profiles, write_dep_profiles, calculate_ray_diag, &
             & write_contour_data, N_pointsX_eq, N_pointsZ_eq, &
             & write_eq_XZ_grid_data, write_eq_radial_profile_data, n_AphiN, &
             & bisection_eps, n_rho, z_reference, &
             & do_OX_conv_analysis

 contains

 subroutine initialize_mirror_processor(read_input)

    use constants_m, only : one
    use diagnostics_m, only : message_unit, message, text_message, verbosity
	use deposition_profiles_m, only : initialize_deposition_profiles

    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder
    if (read_input .eqv. .true.) then

 		n_AphiN = 51 ! Default value, may be reset on input
 		n_rho = n_AphiN ! Default value, may be reset on input
 		bisection_eps = one*10d-6 ! Default value, may be reset on input

    ! Read and write input namelist
  		input_unit = get_unit_number()
        open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
        read(input_unit, mirror_processor_list)
        close(unit=input_unit)
        if (verbosity > 0 )write(message_unit, mirror_processor_list)
    end if
	if (calculate_dep_profiles .eqv. .true.) call initialize_deposition_profiles(read_input)

    if (read_input .eqv. .true.) then
        call text_message('Finished initialize_mirror_processor ', 1)
    end if

    return
 end subroutine initialize_mirror_processor

!*************************************************************************

  subroutine mirror_processor

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message, verbosity
    use deposition_profiles_m, only : calculate_deposition_profiles, &
                         & write_deposition_profiles_LD, write_deposition_profiles_NC
	use OX_conv_analysis_m, only : analyze_OX_conv

    implicit none

    if (verbosity > 0) call text_message('Begin mirror_processor')

    call write_graphics_description_file

    if (calculate_dep_profiles .eqv. .true.) call calculate_deposition_profiles
    if (write_dep_profiles .eqv. .true.) call write_deposition_profiles_NC
	if (calculate_ray_diag .eqv. .true.) call ray_detailed_diagnostics
	if (write_contour_data .eqv. .true.) call write_eq_contour_data_NC
    if (write_eq_radial_profile_data .eqv. .true.) call write_eq_radial_profile_data_NC
    if (do_OX_conv_analysis .eqv. .true.) call analyze_OX_conv

    if (verbosity > 0) call text_message('Finished mirror_processor work')

    return
 end subroutine mirror_processor


!*************************************************************************

  subroutine find_plasma_boundary
    use constants_m, only : rkind
	use multiple_mirror_eq_m, only : magnetics_model
    use diagnostics_m, only : message_unit, message, text_message

	implicit none

	integer :: i
	real(KIND=rkind) :: R

	return
  end subroutine find_plasma_boundary

!*************************************************************************

  subroutine write_graphics_description_file

    use constants_m, only : rkind, zero, one
    use diagnostics_m, only : text_message, run_description, run_label, verbosity
    use multiple_mirror_eq_m, only : box_rmax, box_zmin, box_zmax, r_LUFS
    use bisect_m, only : solve_bisection

    implicit none

    integer :: graphics_descrip_unit, get_unit_number
    integer :: ierr

 !  File name for  output
    character(len=80) :: out_filename

   if (verbosity > 0) call text_message('Writing graphics_description_file')

    ! Open fortran ascii file for results output
    graphics_descrip_unit = get_unit_number()
    out_filename = 'graphics_description_mirror.dat'
    open(unit=graphics_descrip_unit, file=out_filename, &
       & action='write', status='replace', form='formatted')

   write(graphics_descrip_unit, *) 'run_description = ', run_description
   write(graphics_descrip_unit, *) 'run_label = ', run_label

   write(graphics_descrip_unit, *) 'box_rmax = ', box_rmax
   write(graphics_descrip_unit, *) 'box_zmin = ', box_zmin
   write(graphics_descrip_unit, *) 'box_zmax = ', box_zmax

   write(graphics_descrip_unit, *) 'num_plot_k_vectors = ', num_plot_k_vectors
   write(graphics_descrip_unit, *) 'scale_k_vec = ', trim(scale_k_vec)
   write(graphics_descrip_unit, *) 'k_vec_base_length = ', k_vec_base_length
   write(graphics_descrip_unit, *) 'set_XY_lim = ', trim(set_XY_lim)

! Find radius of O-mode cutoff
  call solve_bisection(f_R_alpha_e, r_Omode_cut_at_z_ref, zero, r_LUFS, one,&
       & bisection_eps, ierr)
  if (ierr == 0) r_Omode_cut_at_z_ref = zero

   write(graphics_descrip_unit, *) ' '
   write(graphics_descrip_unit, *) 'z_reference = ', z_reference
   write(graphics_descrip_unit, *) ' '
   write(graphics_descrip_unit, *) 'r_Omode_cut_at_z_ref = ', r_Omode_cut_at_z_ref

   close(unit = graphics_descrip_unit)

  end subroutine write_graphics_description_file

!*************************************************************************

  subroutine ray_detailed_diagnostics()
! Takes the data from the ray_results_m module and for each ray extracts or calculates
! a set of data values for each step along the ray (e.g. ne, omega_ce/omega_rf, Aphi...).
! That data is then written to a netcdf file for analysis or plotting by a graphics code.
! There is a lot of data, and this subroutine is probably only useful for small numbers
! of rays.  The data can be plotted using plot_ray_diags.py which is located in RAYS/graphics_RAYS

    use constants_m, only : rkind, e
    use diagnostics_m, only : integrate_eq_gradients, message, text_message, verbosity
    use species_m, only : nspec,ms
    use multiple_mirror_eq_m, only : multiple_mirror_Aphi
    use equilibrium_m, only : equilibrium, eq_point
    use rf_m, only : omgrf, k0, ray_dispersion_model
    use damping_m, only : damping_model, damping
    use ray_results_m, only : number_of_rays, max_number_of_points, dim_v_vector, npoints,&
        & ray_vec, residual_results => residual, date_vector, RAYS_run_label
    use netcdf

    implicit none

    integer :: iray, istep
    type(eq_point) :: eq

    real(KIND=rkind) :: rvec(3)
    real(KIND=rkind) :: kvec(3), k1, k3
    real(KIND=rkind) :: nvec(3), n1, n3
    real(KIND=rkind) :: v6(6)
    real(KIND=rkind) :: vth

! Needed to call multiple_mirror_Aphi. Only AphiN is used, temp_<> variables are not.
    real(KIND=rkind) :: temp_Aphi, temp_gradAphi(3), AphiN, temp_gradAphiN(3)

! Needed to call damping()
    real(KIND=rkind) :: vg(3)
    real(KIND=rkind) :: ksi(0:nspec), ki

! Needed to call deriv_cold()
    real(KIND=rkind) :: dddx(3), dddk(3), dddw

! netCDF Declarations

    integer :: ncid
! Declarations: dimensions
    integer, parameter :: n_dims = 4
    integer :: d8 = 8
    integer :: number_of_rays_id, max_number_of_points_id, dim_v_vector_id, d8_id

! Declarations: variable IDs
    integer, parameter :: n_vars =  19
    integer :: date_vector_id, npoints_id, s_id, ne_id, Te_kev_id, modB_id, alpha_e_id, &
             & gamma_e_id, Aphi_id, R_id, Z_id, n_par_id, n_perp_id, P_absorbed_id, &
             & n_imag_id, xi_0_id, xi_1_id, xi_2_id, residual_id
!   Declare local arrays
!    integer :: date_vector - from ray_results
!    integer :: npoints - from ray_results
    real(kind=rkind), allocatable :: s(:,:)
    real(kind=rkind), allocatable :: ne(:,:)
    real(kind=rkind), allocatable :: Te_kev(:,:)
    real(kind=rkind), allocatable :: modB(:,:)
    real(kind=rkind), allocatable :: alpha_e(:,:)
    real(kind=rkind), allocatable :: gamma_e(:,:)
    real(kind=rkind), allocatable :: Aphi(:,:)
    real(kind=rkind), allocatable :: R(:,:)
    real(kind=rkind), allocatable :: Z(:,:)
    real(kind=rkind), allocatable :: n_par(:,:)
    real(kind=rkind), allocatable :: n_perp(:,:)
    real(kind=rkind), allocatable :: P_absorbed(:,:)
    real(kind=rkind), allocatable :: n_imag(:,:)
    real(kind=rkind), allocatable :: xi_0(:,:)
    real(kind=rkind), allocatable :: xi_1(:,:)
    real(kind=rkind), allocatable :: xi_2(:,:)
    real(kind=rkind), allocatable :: residual(:,:)

 !  File name for  output
    character(len=128) :: out_filename

    if (verbosity > 0) call text_message('Writing ray_detailed_diagnostics')
    out_filename = 'ray_detailed_diagnostics.'//trim(RAYS_run_label)//'.nc'

!   Allocate local arrays
!    allocate(npoints(number_of_rays))
    allocate(s(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(ne(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(Te_kev(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(modB(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(alpha_e(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(gamma_e(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(Aphi(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(R(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(Z(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(n_par(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(n_perp(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(P_absorbed(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(n_imag(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(xi_0(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(xi_1(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(xi_2(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(residual(max_number_of_points,number_of_rays),source=0.0_rkind)

    ray_loop: do iray = 1, number_of_rays

        step_loop: do istep = 1, npoints(iray)

			v6 = ray_vec(1:6, istep, iray)
        	rvec(:) = v6(1:3)
        	kvec(:) = v6(4:6)
        	s(istep, iray) = ray_vec(7, istep, iray)

        	R(istep, iray) = sqrt(rvec(1)**2 + rvec(2)**2)
       		Z(istep, iray) = rvec(3)

        	call equilibrium(rvec, eq)
        	Te_kev(istep, iray) = eq%Ts(0)/e/1000.0_rkind
        	modB(istep, iray) = eq%bmag
        	alpha_e(istep, iray) = eq%alpha(0)
        	gamma_e(istep, iray) = abs(eq%gamma(0))
        	ne(istep, iray) = eq%ns(0)

        	call multiple_mirror_Aphi(rvec, temp_Aphi, temp_gradAphi, AphiN, temp_gradAphiN)
        	Aphi(istep, iray) = AphiN

			k3 = sum(kvec*eq%bunit)
			k1 = sqrt( sum((kvec-k3*eq%bunit)**2) )
			nvec = kvec/k0; n1 = k1/k0; n3 = k3/k0
			n_par(istep, iray) = n3
			n_perp(istep, iray) = n1
			residual(istep, iray) = residual_results(istep, iray)

	! Damping ki -> n_imag
			damp : if (damping_model /= 'no_damp') then

			!   Calculate the group velocity. So far this is only implemented for cold
			!   dispersion. This part is extracted from EQN_RAY.
			!   First, calculate dD/dk, dD/dx, and dD/d(omega)

					if ( ray_dispersion_model == 'cold' ) then
					   call deriv_cold(eq, nvec, dddx, dddk, dddw)
					else
					   write(*,*) 'ray_detailed_diagnostics: dispersion_model = ', &
					       & ray_dispersion_model
					   stop 'ray_detailed_diagnostics: unimplemented ray_dispersion_model'
					end if

					if ( abs(dddw) > tiny(dddw) ) then
					   vg = -dddk / dddw
					else
					   write(*,*) 'ray_detailed_diagnostics: infinite group velocity'
					   stop 'ray_detailed_diagnostics: infinite group velocity'
					end if

					call damping(eq, v6, vg, ksi, ki)
					n_imag(istep, iray) = ki/k0
					P_absorbed(istep, iray) = ray_vec(8, istep, iray)
		   end if damp

		! Zfunction arguments for electrons
			if (eq%ts(0) > 0. .and. abs(k3) > 0.) then
				! Thermal speed:
				vth = sqrt( 2.*eq%ts(0)/ms(0) )
				! Z function args for 0 to 2nd harmonic
				xi_0(istep, iray) = omgrf/(k3*vth)
				xi_1(istep, iray) = (omgrf+eq%omgc(0))/(k3*vth)
				xi_2(istep, iray) = (omgrf+2.*eq%omgc(0))/(k3*vth)
			end if

		end do step_loop

	end do ray_loop

! netCDF coding

!   Open NC file
    call check( nf90_create(trim(out_filename), nf90_clobber, ncid) )

!   Define NC dimensions
    call check( nf90_def_dim(ncid, 'number_of_rays', number_of_rays, number_of_rays_id))
    call check( nf90_def_dim(ncid, 'max_number_of_points', max_number_of_points, max_number_of_points_id))
    call check( nf90_def_dim(ncid, 'dim_v_vector', dim_v_vector, dim_v_vector_id))
    call check( nf90_def_dim(ncid, 'd8', 8, d8_id))

! Define NC variables
    call check( nf90_def_var(ncid, 'date_vector', NF90_INT, [d8_id], date_vector_id))
    call check( nf90_def_var(ncid, 'npoints', NF90_INT, [number_of_rays_id], npoints_id))
    call check( nf90_def_var(ncid, 's', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], s_id))
    call check( nf90_def_var(ncid, 'ne', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], ne_id))
    call check( nf90_def_var(ncid, 'Te_kev', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], Te_kev_id))
    call check( nf90_def_var(ncid, 'modB', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], modB_id))
    call check( nf90_def_var(ncid, 'alpha_e', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], alpha_e_id))
    call check( nf90_def_var(ncid, 'gamma_e', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], gamma_e_id))
    call check( nf90_def_var(ncid, 'Aphi', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], Aphi_id))
    call check( nf90_def_var(ncid, 'R', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], R_id))
    call check( nf90_def_var(ncid, 'Z', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], Z_id))
    call check( nf90_def_var(ncid, 'n_par', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], n_par_id))
    call check( nf90_def_var(ncid, 'n_perp', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], n_perp_id))
    call check( nf90_def_var(ncid, 'P_absorbed', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], P_absorbed_id))
    call check( nf90_def_var(ncid, 'n_imag', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], n_imag_id))
    call check( nf90_def_var(ncid, 'xi_0', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], xi_0_id))
    call check( nf90_def_var(ncid, 'xi_1', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], xi_1_id))
    call check( nf90_def_var(ncid, 'xi_2', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], xi_2_id))
    call check( nf90_def_var(ncid, 'residual', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], residual_id))

! Put global attributes
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'RAYS_run_label', RAYS_run_label))

call check( nf90_enddef(ncid))

! Put NC variables
    call check( nf90_put_var(ncid, date_vector_id, date_vector))
    call check( nf90_put_var(ncid, npoints_id, npoints))
    call check( nf90_put_var(ncid, s_id, s))
    call check( nf90_put_var(ncid, ne_id, ne))
    call check( nf90_put_var(ncid, Te_kev_id, Te_kev))
    call check( nf90_put_var(ncid, modB_id, modB))
    call check( nf90_put_var(ncid, alpha_e_id, alpha_e))
    call check( nf90_put_var(ncid, gamma_e_id, gamma_e))
    call check( nf90_put_var(ncid, Aphi_id, Aphi))
    call check( nf90_put_var(ncid, R_id, R))
    call check( nf90_put_var(ncid, Z_id, Z))
    call check( nf90_put_var(ncid, n_par_id, n_par))
    call check( nf90_put_var(ncid, n_perp_id, n_perp))
    call check( nf90_put_var(ncid, P_absorbed_id, P_absorbed))
    call check( nf90_put_var(ncid, n_imag_id, n_imag))
    call check( nf90_put_var(ncid, xi_0_id, xi_0))
    call check( nf90_put_var(ncid, xi_1_id, xi_1))
    call check( nf90_put_var(ncid, xi_2_id, xi_2))
    call check( nf90_put_var(ncid, residual_id, residual))

!   Close the NC file
    call check( nf90_close(ncid) )

  end subroutine ray_detailed_diagnostics

!****************************************************************************

subroutine write_eq_contour_data_NC
! Write a netcdf file ('eq_contours.<run_label>.nc') containing data needed to plot
! contours of quantities like gamma = cyclotron freq/wave freq.  The contours can then be
! plots of ray trajectories

	use constants_m, only : zero, one
    use diagnostics_m, only : message_unit, message, text_message, verbosity, run_label
    use species_m, only : nspec, spec_name
	use equilibrium_m, only : eq_point, equilibrium, write_eq_point
    use multiple_mirror_eq_m, only : box_rmax, box_zmin, box_zmax, &
         & plasma_AphiN_limit, multiple_mirror_Aphi, multiple_mirror_eq
    use netcdf

    implicit none

    real(KIND=rkind) :: box_xmin, box_xmax
    real(KIND=rkind) :: X(N_pointsX_eq), Z(N_pointsZ_eq), dX, dZ
    real(KIND=rkind) :: AphiN_out(N_pointsZ_eq, N_pointsX_eq)
	integer i,j
    real(KIND=rkind) :: Aphi, gradAphi(3), AphiN, gradAphiN(3)
    real(KIND=rkind) :: rvec(3)

! Stash plasma_AphiN_limit from mirror_eq_m.  Danger of circularity, equilibrium_m
! uses mirror_eq_m.
    real(KIND=rkind) ::plasma_AphiN_limit_temp

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

!	omega_cyclotron/omega
    real(KIND=rkind), allocatable :: gamma(:)
    real(KIND=rkind), allocatable :: gamma_array(:,:,:)

!	plasma_cyclotron/omega
    real(KIND=rkind), allocatable :: omega_pN(:)
    real(KIND=rkind), allocatable :: omega_pN_array(:,:,:)

! netCDF Declarations
    integer :: ncid
! Declarations: dimensions
    integer, parameter :: n_dims = 2
    integer :: n_X, n_Z, nspec_p1, d12_id
    integer :: n_X_id, n_Z_id, nspec_p1_id, spec_name_id
! Declarations: variable IDs
    integer :: X_id, Z_id, AphiN_id, omega_pN_array_id, gamma_array_id
    integer :: box_xmin_id, box_xmax_id, box_zmin_id, box_zmax_id

 !  File name for  output
    character(len=80) :: out_filename

    real(KIND=rkind) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind) :: ns(0:nspec), gradns(3,0:nspec)
    real(KIND=rkind) :: ts(0:nspec), gradts(3,0:nspec)
    character(len=60) :: equib_err

   if (verbosity > 0) call text_message('Writing equilibrium contours of Aphi')

!   Open NC file
    out_filename = 'eq_contours.'//trim(run_label)//'.nc'
    call check( nf90_create(trim(out_filename), nf90_clobber, ncid) )
    write(*,*) 'N_pointsX_eq = ', N_pointsX_eq, '  N_pointsZ_eq = ', N_pointsZ_eq

!   Define NC dimensions
    call check( nf90_def_dim(ncid, 'n_X', N_pointsX_eq, n_X_id))
    call check( nf90_def_dim(ncid, 'n_Z', N_pointsZ_eq, n_Z_id))
    nspec_p1 = nspec + 1
    call check( nf90_def_dim(ncid, 'nspec_p1', nspec_p1, nspec_p1_id))
    call check( nf90_def_dim(ncid, 'd12', 12, d12_id))

! Define NC variables
    call check( nf90_def_var(ncid, 'box_xmin', NF90_DOUBLE, box_xmin_id))
    call check( nf90_def_var(ncid, 'box_xmax', NF90_DOUBLE, box_xmax_id))
    call check( nf90_def_var(ncid, 'box_zmin', NF90_DOUBLE, box_zmin_id))
    call check( nf90_def_var(ncid, 'box_zmax', NF90_DOUBLE, box_zmax_id))
    call check( nf90_def_var(ncid, 'X', NF90_DOUBLE, [n_X_id], X_id))
    call check( nf90_def_var(ncid, 'Z', NF90_DOUBLE, [n_Z_id], Z_id))
    call check( nf90_def_var(ncid, 'AphiN', NF90_DOUBLE, [n_Z_id, n_X_id], AphiN_id))
    call check( nf90_def_var(ncid, 'omega_pN_array', NF90_DOUBLE, [n_Z_id, n_X_id,nspec_p1_id],&
                                  & omega_pN_array_id))
    call check( nf90_def_var(ncid, 'gamma_array', NF90_DOUBLE, [n_Z_id ,n_X_id,nspec_p1_id],&
                                  & gamma_array_id))
    call check( nf90_def_var(ncid, 'spec_name', NF90_CHAR, [d12_id,nspec_p1_id], spec_name_id))

! Put global attributes
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'RAYS_run_label', run_label))

! Finished NC definition
	call check( nf90_enddef(ncid))

	allocate(omega_pN(nspec_p1))
	allocate(gamma(nspec_p1))
	allocate(omega_pN_array(N_pointsZ_eq, N_pointsX_eq, nspec_p1))
	allocate(gamma_array(N_pointsZ_eq, N_pointsX_eq, nspec_p1))

	! N.B. X ranges from -rmax to rmax
	box_xmin = -box_rmax
	box_xmax = box_rmax

	dx = one/(N_pointsX_eq - 1)*(box_xmax-box_xmin)
	dz = one/(N_pointsZ_eq - 1)*(box_zmax-box_zmin)
    do i = 1, N_pointsZ_eq
 		Z(i) = box_zmin + (i-1)*dz
	end do
	Z(1) = box_zmin; Z(N_pointsZ_eq) = box_zmax
	do j = 1, N_pointsX_eq
		X(j) = box_xmin + (j-1)*dx
	end do
	X(1) = box_xmin; X(N_pointsX_eq) = box_xmax
! write(*,*) 'X = ', X
! write(*,*) 'Z = ', Z
	do i = 1, N_pointsZ_eq
	do j = 1, N_pointsX_eq
		rvec = (/X(j), zero, Z(i)/)

		call multiple_mirror_Aphi(rvec, Aphi, gradAphi, AphiN, gradAphiN)
		AphiN_out(i,j) = AphiN
!		write(*,*) 'rvec = ', rvec, '   AphiN_out(i,j) = ', AphiN_out(i,j)
!		write(*,*) 'i = ', i,'  j = ', j, '   AphiN_out(i,j) = ', AphiN_out(i,j)

		plasma_AphiN_limit_temp = plasma_AphiN_limit ! Stash plasma_AphiN_limit
		plasma_AphiN_limit = 10.0 ! Set plasma_AphiN_limit high so resonance contours
		                          ! go out to box
		call equilibrium(rvec, eq)

		plasma_AphiN_limit = plasma_AphiN_limit_temp ! Restore plasma_AphiN_limit

		omega_pN(:) = sqrt(eq%alpha(0:nspec))
!		write(*,*) 'rvec = ', rvec, '  eq%alpha(0:nspec) = ', eq%alpha(0:nspec)
		omega_pN_array(i, j, :) = omega_pN
		gamma(:) = eq%gamma(0:nspec)
		gamma_array(i, j, :) = gamma
	end do
	end do

! Put NC variables
    call check( nf90_put_var(ncid, box_xmin_id, box_xmin))
    call check( nf90_put_var(ncid, box_xmax_id, box_xmax))
    call check( nf90_put_var(ncid, box_zmin_id, box_zmin))
    call check( nf90_put_var(ncid, box_zmax_id, box_zmax))
    call check( nf90_put_var(ncid, X_id, X))
    call check( nf90_put_var(ncid, Z_id, Z))
    call check( nf90_put_var(ncid, AphiN_id, AphiN_out(:,:)))
    call check( nf90_put_var(ncid, omega_pN_array_id, omega_pN_array))
    call check( nf90_put_var(ncid, gamma_array_id, gamma_array))
    call check( nf90_put_var(ncid, spec_name_id, spec_name(0:nspec)))

!   Close the NC file
    call check( nf90_close(ncid) )
    return
 end subroutine write_eq_contour_data_NC


!*************************************************************************

 subroutine write_eq_radial_profile_data_NC
! Write a netcdf file ('eq_radial.<run_label>.nc') containing data needed to plot
! radial profiles for equilibrium quantites like ne, Te, etc as functions of
! AphiN and also R.  These profiles are evaluated at Z = z_reference, which is an input
! parameter in namelist group /mirror_processor_list/.  Profiles extend radially from zero
! to plasma_AphiN_limit, which is set as a namelist input to multiple_mirror_eq_m.

    use constants_m, only : rkind, one, zero, e
    use diagnostics_m, only : message_unit, message, text_message, verbosity, run_label
    use species_m, only : nspec, spec_name
	use equilibrium_m, only : eq_point, equilibrium, write_eq_point
    use multiple_mirror_eq_m, only : box_rmax, box_zmin, box_zmax, r_LUFS, &
         & Aphi_LUFS, plasma_AphiN_limit, multiple_mirror_Aphi, multiple_mirror_eq
    use bisect_m, only : solve_bisection
    use ray_results_m, only : date_vector, RAYS_run_label
    use netcdf
	use XY_curves_netCDF_m, only : XY_curve_netCDF, write_XY_curves_netCDF

    implicit none

! Number of profiles to generate and list of profiles
    integer, parameter :: n_profiles = 8
	type(XY_curve_netCDF) :: profile_list(n_profiles)
	integer :: n_grid(n_profiles)

!   Declare local variables
    real(KIND=rkind) :: AphiN(n_AphiN), R(n_AphiN)
    real(KIND=rkind) :: ne_AphiN(n_AphiN), Te_AphiN_ev(n_AphiN), Ti_AphiN_ev(n_AphiN)
    real(KIND=rkind) :: RBphi_AphiN(n_AphiN), rho_AphiN(n_AphiN)

    real(KIND=rkind) :: rho(n_rho), AphiN_rho(n_rho)
    real(KIND=rkind) :: ne_rho(n_rho), Te_rho_ev(n_rho), Ti_rho_ev(n_rho)

	integer :: i, ierr
    real(KIND=rkind) :: Aphi, gradAphi(3), gradAphiN(3), AphiNij, gradAphiNij(3),&
                      & gradAphiN_rho(3)
    real(KIND=rkind) :: rvec(3)

    ! Args to spline profiles, not used
    real(KIND=rkind) :: dRBphi_dAphi, drho_dAphi
    real(KIND=rkind) :: dAphi_drho, dRBphi_drho, drho_drho

! Stash plasma_AphiN_limit from mirror_eq_m.  Danger of circularity, equilibrium_m
! uses mirror_eq_m.
    real(KIND=rkind) ::plasma_AphiN_limit_temp

!    real(KIND=rkind) :: r_LUFS ! Radius of LUFS at z = z_reference, to set up r(Aphi) grid.

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

 !  File name for  output
    character(len=80) :: out_filename

    real(KIND=rkind) :: ns(0:nspec)
    real(KIND=rkind) :: ts(0:nspec)
    character(len=60) :: equib_err

   if (verbosity > 0) call text_message('Writing plasma equilibrium profiles')
   if (verbosity > 0) call message('Writing plasma z_reference', z_reference, 1)

! Get character form of z_reference
	write(char_z_ref, fmt = '(f10.3)') z_reference
	char_z_ref = adjustl(char_z_ref)

	plasma_AphiN_limit_temp = plasma_AphiN_limit ! Stash plasma_AphiN_limit
!	plasma_AphiN_limit = 10.0 ! Set plasma_AphiN_limit high so can get fields outside plasma

! Note: Want to plot ne, Te etc on a uniform Aphi grid, but they  are given by subroutine
! equilibrium(x,y,z).  So ...
! Generate Aphi grid and R grid.
! Invert AphiN(R, y=0, z_reference) -> R(AphiN, z_reference) using bisection.
! Then evaluate ne and Te at (R, 0, z_reference).
!
 write(*,*) 'r_LUFS = ', r_LUFS
	do i = 1, n_AphiN
	    AphiN(i) = plasma_AphiN_limit*(i-1)/(n_AphiN - 1)
	    call solve_bisection(f_R_AphiN, R(i), zero, 1.1_rkind*r_LUFS, AphiN(i),&
	                       & bisection_eps, ierr)

	    rvec = (/R(i), zero, z_reference/)
!    write(*,*) 'write_eq_radial_profile_data_NC, i = ', i, '  R(i) = ',  R(i), '   AphiN(i) = ',AphiN(i)
 	    if (ierr == 0) stop

		call equilibrium(rvec, eq)
		ne_AphiN(i) = eq%ns(0)
		Te_AphiN_ev(i) = eq%Ts(0)/e ! Convert from Joules to ev
		Ti_AphiN_ev(i) = eq%Ts(1)/e ! N.B. For now all ions are assumed to have the same Ti profile
	end do

!  write(*,*) " "
!  write(*,*) "AphiN", AphiN
!  write(*,*) " "
!  write(*,*) "R = ", R
!  write(*,*) " "
!  write(*,*) "ne_AphiN =  ", ne_AphiN
!  write(*,*) " "
!  write(*,*) "rho_AphiN =  ", rho_AphiN
!  write(*,*) 'bisection ierr = ', ierr
!
! !***************** Stuff for profiles versus rho *****************************

! Note: For historical reasons, and to avoid name collisions, in this section R is called
!       rho.  Want to plot on a uniform radial grid, whereas the R(i) grid constructed
!       above is uniform in AphiN. So generate a new rho grid extending from zero to
!       R(plasma_AphiN_limit) which R(n_AphiN) as calculated above at z = z_reference.
   if (verbosity > 0) call text_message('Writing plasma equilibrium profiles on R <-> rho grid')

	do i = 1, n_rho
	    rho(i) =  R(n_AphiN)*(i-1)/(n_rho - 1)
	    rvec = (/rho(i), zero, z_reference/)

		call equilibrium(rvec, eq)
		ne_rho(i) = eq%ns(0)
		Te_rho_ev(i) = eq%Ts(0)/e ! Convert from Joules to ev
		Ti_rho_ev(i) = eq%Ts(1)/e ! N.B. For now all ions are assumed to have the same Ti profile
		call multiple_mirror_Aphi(rvec, Aphi, gradAphi, AphiN_rho(i), gradAphiN_rho)
		! N.B. Aphi, gradAphi, and gradAphiN_rho are effectively dummy args. Only AphiN_rho
		!      is used below.
	end do

!  write(*,*) " "
!  write(*,*) "rho", rho
!  write(*,*) " "
!  write(*,*) "R = ", R
!  write(*,*) " "
!  write(*,*) "ne_rho =  ", ne_rho

! Load ne data into profile_list
	profile_list(1)%grid_name = 'AphiN'
	profile_list(1)%curve_name = 'ne(AphiN)'
	n_grid(1) = n_AphiN
	allocate(profile_list(1)%grid(n_grid(1)), source = 0.0_rkind)
	profile_list(1)%grid(:) = AphiN(:)
	allocate(profile_list(1)%curve(n_grid(1)), source = 0.0_rkind)
	profile_list(1)%curve(:) = ne_AphiN(:)

! Load Te data into profile_list
	profile_list(2)%grid_name = 'AphiN'
	profile_list(2)%curve_name = 'Te(AphiN)'
	n_grid(2) = n_AphiN
	allocate(profile_list(2)%grid(n_grid(2)), source = 0.0_rkind)
	profile_list(2)%grid(:) = AphiN(:)
	allocate(profile_list(2)%curve(n_grid(2)), source = 0.0_rkind)
	profile_list(2)%curve(:) = Te_AphiN_ev(:)

! Load Ti data into profile_list
	profile_list((3))%grid_name = 'AphiN'
	profile_list((3))%curve_name = 'Ti(AphiN)'
	n_grid((3)) = n_AphiN
	allocate(profile_list((3))%grid(n_grid((3))), source = 0.0_rkind)
	profile_list((3))%grid(:) = AphiN(:)
	allocate(profile_list((3))%curve(n_grid((3))), source = 0.0_rkind)
	profile_list((3))%curve(:) = Ti_AphiN_ev(:)

! Load R data into profile_list
	profile_list((4))%grid_name = 'AphiN'
	profile_list((4))%curve_name = 'R(AphiN)'
	n_grid((4)) = n_AphiN
	allocate(profile_list((4))%grid(n_grid((4))), source = 0.0_rkind)
	profile_list((4))%grid(:) = AphiN(:)
	allocate(profile_list((4))%curve(n_grid((4))), source = 0.0_rkind)
	profile_list((4))%curve(:) = R(:)

! !***************** Stuff for profiles versus rho *****************************
! Load AphiN of rho data into profile_list
	profile_list((5))%grid_name = 'R'
	profile_list((5))%curve_name = 'AphiN(R, z = '//trim(char_z_ref)//')'
	n_grid((5)) = n_rho
	allocate(profile_list((5))%grid(n_grid((5))), source = 0.0_rkind)
	profile_list((5))%grid(:) = rho(:)
	allocate(profile_list((5))%curve(n_grid((5))), source = 0.0_rkind)
	profile_list((5))%curve(:) = AphiN_rho(:)

! Load ne of rho data into profile_list
	profile_list((6))%grid_name = 'R'
	profile_list((6))%curve_name = 'ne(R, z = '//trim(char_z_ref)//')'
	n_grid((6)) = n_rho
	allocate(profile_list((6))%grid(n_grid((6))), source = 0.0_rkind)
	profile_list((6))%grid(:) = rho(:)
	allocate(profile_list((6))%curve(n_grid((6))), source = 0.0_rkind)
	profile_list((6))%curve(:) = ne_rho(:)

! Load Te of rho data into profile_list
	profile_list((7))%grid_name = 'R'
	profile_list((7))%curve_name = 'Te(R, z = '//trim(char_z_ref)//')'
	n_grid((7)) = n_rho
	allocate(profile_list((7))%grid(n_grid((7))), source = 0.0_rkind)
	profile_list((7))%grid(:) = rho(:)
	allocate(profile_list((7))%curve(n_grid((7))), source = 0.0_rkind)
	profile_list((7))%curve(:) = Te_rho_ev(:)

! Load Ti of rho data into profile_list
	profile_list(((8)))%grid_name = 'R'
	profile_list(((8)))%curve_name = 'Ti(R, z = '//trim(char_z_ref)//')'
	n_grid(((8))) = n_rho
	allocate(profile_list(((8)))%grid(n_grid(((8)))), source = 0.0_rkind)
	profile_list(((8)))%grid(:) = rho(:)
	allocate(profile_list(((8)))%curve(n_grid(((8)))), source = 0.0_rkind)
	profile_list(((8)))%curve(:) = Ti_rho_ev(:)

! Restore plasma_AphiN_limit
	plasma_AphiN_limit = plasma_AphiN_limit_temp

! Generate netCDF file
!  write(*,*) 'profile_list(1)%curve_name = ', profile_list(1)%curve_name
!  write(*,*) 'profile_list(9)%curve_name = ', profile_list(9)%curve_name
	out_filename = 'eq_radial_profiles.'//trim(run_label)
	call write_XY_curves_netCDF(profile_list, out_filename)

    return
 end subroutine write_eq_radial_profile_data_NC

!*************************************************************************

  subroutine check(status)
    use netcdf
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check

!*************************************************************************

 function f_R_alpha_e(R)
! Returns alpha_e(R,z_reference)

    use constants_m, only : rkind, zero
	use equilibrium_m, only : eq_point, equilibrium

	IMPLICIT NONE
    real(KIND=rkind) f_R_alpha_e
    real(KIND=rkind), intent(in) :: R
    real(KIND=rkind) :: rvec(3)
    type(eq_point) :: eq

	rvec = (/R, zero, z_reference/)
	    rvec = (/R, zero, z_reference/)

	call equilibrium(rvec, eq)
	f_R_alpha_e = eq%alpha(0)
!	write(*,*) 'rvec', rvec,  '  f_R_alpha_e = ', f_R_alpha_e

	return
 end function f_R_alpha_e

!*************************************************************************

 function f_R_AphiN(R)
! Returns AphiN(R,z_reference)

    use constants_m, only : rkind, zero
	use multiple_mirror_eq_m, only : multiple_mirror_Aphi

	IMPLICIT NONE
    real(KIND=rkind) f_R_AphiN
    real(KIND=rkind), intent(in) :: R
    real(KIND=rkind) :: Aphi, gradAphi(3), AphiN, gradAphiN(3)
    real(KIND=rkind) :: rvec(3)

	rvec = (/R, zero, z_reference/)
	call multiple_mirror_Aphi(rvec, Aphi, gradAphi, AphiN, gradAphiN)
	f_R_AphiN = AphiN
write(*,*) 'rvec = ', rvec,  '  f_R_Aphi = ', f_R_AphiN

	return
 end function f_R_AphiN

 end module mirror_processor_m

