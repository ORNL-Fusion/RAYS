 module axisym_toroid_processor_m
! Post processing for axisym_toroid equilibrium

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

! Number of points in R,Z grid for normalized psi netCDF file
	integer ::  N_pointsR_eq, N_pointsZ_eq

! Flags determining what to do besides plot rays.  Can be overridden (i.e. turned off)
! in namelist file
	logical :: calculate_dep_profiles = .true. ! Calculate all depositon profiles
	logical :: write_dep_profiles = .true. ! Write depositon profile netCDF file
	logical :: calculate_ray_diag = .true. ! Write detailed ray diagnostic netCDF file
	logical :: write_contour_data = .true.  ! Write data needed to plot contours netCDF file

    namelist /axisym_toroid_processor_list/ num_plot_k_vectors, scale_k_vec,&
             & k_vec_base_length, set_XY_lim, &
             & calculate_dep_profiles, write_dep_profiles, calculate_ray_diag, &
             & write_contour_data, N_pointsR_eq, N_pointsZ_eq

 contains

 subroutine initialize_axisym_toroid_processor(read_input)

    use diagnostics_m, only : message_unit, message, text_message, verbosity
	use deposition_profiles_m, only : initialize_deposition_profiles

    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder
    if (read_input .eqv. .true.) then
    ! Read and write input namelist
  		input_unit = get_unit_number()
        open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
        read(input_unit, axisym_toroid_processor_list)
        close(unit=input_unit)
        if (verbosity > 0 )write(message_unit, axisym_toroid_processor_list)
        call text_message('Finished initialize_axisym_toroid_processor ', 1)
    end if

	if (calculate_dep_profiles .eqv. .true.) call initialize_deposition_profiles(read_input)

    return
 end subroutine initialize_axisym_toroid_processor

!*************************************************************************

  subroutine axisym_toroid_processor

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message, verbosity
    use deposition_profiles_m, only : calculate_deposition_profiles, &
                         & write_deposition_profiles_LD, write_deposition_profiles_NC

    implicit none

    if (verbosity > 0) call text_message('Begin axisym_toroid_processor')

    call write_graphics_description_file

    if (calculate_dep_profiles .eqv. .true.) call calculate_deposition_profiles
    if (write_dep_profiles .eqv. .true.) call write_deposition_profiles_NC

	if (calculate_ray_diag .eqv. .true.) call ray_detailed_diagnostics

	if (write_contour_data .eqv. .true.) call write_eq_contour_data_NC

    if (verbosity > 0) call text_message('Finished axisym_toroid_processor work')

    return
 end subroutine axisym_toroid_processor


!*************************************************************************

  subroutine find_plasma_boundary
    use constants_m, only : rkind
	use axisym_toroid_eq_m, only : magnetics_model, inner_bound, outer_bound, upper_bound, lower_bound
	use solovev_magnetics_m, only : rmaj, kappa
    use diagnostics_m, only : message_unit, message, text_message
    use eqdsk_utilities_m, only : NBOUND, RBOUND, ZBOUND

	implicit none

	integer :: i
	real(KIND=rkind) :: R, dR, Zsq

    magnetics: select case (trim(magnetics_model))

    case ('solovev_magnetics')  ! N.B.  This is up-down symmetric
		n_boundary_points = 271 ! should be an odd number for this calculation
		allocate(R_boundary(n_boundary_points))
		allocate(Z_boundary(n_boundary_points))
		dR = 2.*(outer_bound - inner_bound)/(n_boundary_points)

		R_boundary(1) = inner_bound
		Z_boundary(1) = 0.
		R_boundary(n_boundary_points) = inner_bound
		Z_boundary(n_boundary_points) = 0.
		R_boundary((n_boundary_points-1)/2+1) = outer_bound
		Z_boundary(n_boundary_points) = 0.

		do i = 2, (n_boundary_points -1)/2
			R = inner_bound + i*dR
			Zsq = kappa*2/(4.*R**2)*(outer_bound**4 + 2.*(R**2 - outer_bound**2)*rmaj**2 -&
			      & R**4)
			R_boundary(i) = R
			Z_boundary(i) = sqrt(Zsq)
			R_boundary(n_boundary_points-(i-1)) = R_boundary(i)
			Z_boundary(n_boundary_points-(i-1)) = -Z_boundary(i)
		end do

    case ('eqdsk_magnetics_lin_interp')
        n_boundary_points = NBOUND ! not necessarily an odd number
		allocate(R_boundary(n_boundary_points))
		allocate(Z_boundary(n_boundary_points))
        R_boundary(:) = RBOUND(:)
        Z_boundary(:) = ZBOUND(:)

    case ('eqdsk_magnetics_spline_interp')
        n_boundary_points = NBOUND ! not necessarily an odd number
		allocate(R_boundary(n_boundary_points))
		allocate(Z_boundary(n_boundary_points))
        R_boundary(:) = RBOUND(:)
        Z_boundary(:) = ZBOUND(:)

! 	do i = 1, n_boundary_points
! 		write(*,*) 'i = ', i, '   R_boundary = ', R_boundary(i), '   Z_boundary', Z_boundary(i)
! 	end do

    case default
	  write(0,*) 'initialize_axisym_toroid_eq: unknown magnetics model =', magnetics_model
	  call text_message('initialize_axisym_toroid_eq: unknown magnetics model',&
	  & trim(magnetics_model),0)
	  stop 1
    end select magnetics

  end subroutine find_plasma_boundary

!*************************************************************************

  subroutine write_graphics_description_file

   use diagnostics_m, only : text_message, run_description, run_label, verbosity
   use axisym_toroid_eq_m, only : r_axis, z_axis, &
                          & box_rmin, box_rmax, box_zmin, box_zmax, &
                          & inner_bound, outer_bound, upper_bound, lower_bound

    implicit none

    integer :: graphics_descrip_unit, get_unit_number

 !  File name for  output
    character(len=80) :: out_filename

   if (verbosity > 0) call text_message('write_graphics_description_file')

    ! Open fortran ascii file for results output
    graphics_descrip_unit = get_unit_number()
    out_filename = 'graphics_description_axisym_toroid.dat'
    open(unit=graphics_descrip_unit, file=out_filename, &
       & action='write', status='replace', form='formatted')

   write(graphics_descrip_unit, *) 'run_description = ', run_description
   write(graphics_descrip_unit, *) 'run_label = ', run_label

   write(graphics_descrip_unit, *) 'r_axis = ', r_axis
   write(graphics_descrip_unit, *) 'z_axis = ', z_axis
   write(graphics_descrip_unit, *) 'inner_bound = ', inner_bound
   write(graphics_descrip_unit, *) 'outer_bound = ', outer_bound
   write(graphics_descrip_unit, *) 'upper_bound = ', upper_bound
   write(graphics_descrip_unit, *) 'lower_bound = ', lower_bound

   write(graphics_descrip_unit, *) 'box_rmin = ', box_rmin
   write(graphics_descrip_unit, *) 'box_rmax = ', box_rmax
   write(graphics_descrip_unit, *) 'box_zmin = ', box_zmin
   write(graphics_descrip_unit, *) 'box_zmax = ', box_zmax

   write(graphics_descrip_unit, *) 'num_plot_k_vectors = ', num_plot_k_vectors
   write(graphics_descrip_unit, *) 'scale_k_vec = ', trim(scale_k_vec)
   write(graphics_descrip_unit, *) 'k_vec_base_length = ', k_vec_base_length
   write(graphics_descrip_unit, *) 'set_XY_lim = ', trim(set_XY_lim)

   call find_plasma_boundary

   write(graphics_descrip_unit, *) ' '
   write(graphics_descrip_unit, *) 'R_boundary = ', R_boundary
   write(graphics_descrip_unit, *) ' '
   write(graphics_descrip_unit, *) 'Z_boundary = ', Z_boundary

   close(unit = graphics_descrip_unit)

  end subroutine write_graphics_description_file

!*************************************************************************

  subroutine ray_detailed_diagnostics()
! Takes the data from the ray_results_m module and for each ray extracts or calculates
! a set of data values for each step along the ray (e.g. ne, omega_ce/omega_rf, psi...).
! That data is then written to a netcdf file for analysis or plotting by a graphics code.
! There is a lot of data, and this subroutine is probably only useful for small numbers
! of rays.

    use constants_m, only : rkind
    use diagnostics_m, only : integrate_eq_gradients, message, text_message
    use species_m, only : nspec,ms
    use axisym_toroid_eq_m, only : axisym_toroid_psi
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

! Needed to call axisym_toroid_psi. Only psiN is used, temp_<> variables are not.
    real(KIND=rkind) :: temp_psi, temp_gradpsi(3), psiN, temp_gradpsiN(3)

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
    integer :: date_vector_id, npoints_id, s_id, ne_id, Te_id, modB_id, alpha_e_id, &
             & gamma_e_id, Psi_id, R_id, Z_id, n_par_id, n_perp_id, P_absorbed_id, &
             & n_imag_id, xi_0_id, xi_1_id, xi_2_id, residual_id
!   Declare local arrays
!    integer :: date_vector - from ray_results
!    integer :: npoints - from ray_results
    real(kind=rkind), allocatable :: s(:,:)
    real(kind=rkind), allocatable :: ne(:,:)
    real(kind=rkind), allocatable :: Te(:,:)
    real(kind=rkind), allocatable :: modB(:,:)
    real(kind=rkind), allocatable :: alpha_e(:,:)
    real(kind=rkind), allocatable :: gamma_e(:,:)
    real(kind=rkind), allocatable :: Psi(:,:)
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

    out_filename = 'ray_detailed_diagnostics.'//trim(RAYS_run_label)//'.nc'

!   Allocate local arrays
!    allocate(npoints(number_of_rays))
    allocate(s(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(ne(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(Te(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(modB(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(alpha_e(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(gamma_e(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(Psi(max_number_of_points,number_of_rays),source=0.0_rkind)
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
        	Te(istep, iray) = eq%Ts(0)
        	modB(istep, iray) = eq%bmag
        	alpha_e(istep, iray) = eq%alpha(0)
        	gamma_e(istep, iray) = abs(eq%gamma(0))
        	ne(istep, iray) = eq%ns(0)

        	call axisym_toroid_psi(rvec, temp_psi, temp_gradpsi, psiN, temp_gradpsiN)
        	Psi(istep, iray) = psiN

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

		! Zfunction arguments
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
    call check( nf90_def_var(ncid, 'Te', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], Te_id))
    call check( nf90_def_var(ncid, 'modB', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], modB_id))
    call check( nf90_def_var(ncid, 'alpha_e', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], alpha_e_id))
    call check( nf90_def_var(ncid, 'gamma_e', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], gamma_e_id))
    call check( nf90_def_var(ncid, 'Psi', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], Psi_id))
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
    call check( nf90_put_var(ncid, Te_id, Te))
    call check( nf90_put_var(ncid, modB_id, modB))
    call check( nf90_put_var(ncid, alpha_e_id, alpha_e))
    call check( nf90_put_var(ncid, gamma_e_id, gamma_e))
    call check( nf90_put_var(ncid, Psi_id, Psi))
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
! contours of equilibrium quantites like psi and omegaCj/omega => beta(j)
! flux evaluated on a uniformly spaced  R,Z grid [box_rmin, box_rmax, box_zmin, box_zmax]

	use constants_m, only : zero, one
    use diagnostics_m, only : message_unit, message, text_message, verbosity, run_label
    use species_m, only : nspec, spec_name
	use equilibrium_m, only : eq_point, equilibrium, write_eq_point
    use axisym_toroid_eq_m, only : box_rmin, box_rmax, box_zmin, box_zmax, &
         & axisym_toroid_psi, axisym_toroid_eq
    use netcdf

    implicit none

    real(KIND=rkind) :: R(N_pointsR_eq), Z(N_pointsZ_eq), dR, dZ
    real(KIND=rkind) :: psiN_out(N_pointsR_eq, N_pointsZ_eq)
	integer i,j
    real(KIND=rkind) :: psi, gradpsi(3), psiN, gradpsiN(3)
    real(KIND=rkind) :: rvec(3)

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

!	omega_cyclotron/omega
    real(KIND=rkind), allocatable :: gamma(:)
    real(KIND=rkind), allocatable :: gamma_array(:,:,:)

! netCDF Declarations
    integer :: ncid
! Declarations: dimensions
    integer, parameter :: n_dims = 2
    integer :: n_R, n_Z, nspec_p1, d12_id
    integer :: n_R_id, n_Z_id, nspec_p1_id, spec_name_id
! Declarations: variable IDs
    integer :: R_id, Z_id, psiN_id, gamma_array_id
    integer :: box_rmin_id, box_rmax_id, box_zmin_id, box_zmax_id

 !  File name for  output
    character(len=80) :: out_filename

    real(KIND=rkind) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind) :: ns(0:nspec), gradns(3,0:nspec)
    real(KIND=rkind) :: ts(0:nspec), gradts(3,0:nspec)
    character(len=60) :: equib_err


!   Open NC file
    out_filename = 'eq_contours.'//trim(run_label)//'.nc'
    call check( nf90_create(trim(out_filename), nf90_clobber, ncid) )

!   Define NC dimensions
    call check( nf90_def_dim(ncid, 'n_R', N_pointsR_eq, n_R_id))
    call check( nf90_def_dim(ncid, 'n_Z', N_pointsZ_eq, n_Z_id))
    nspec_p1 = nspec + 1
    call check( nf90_def_dim(ncid, 'nspec_p1', nspec_p1, nspec_p1_id))
    call check( nf90_def_dim(ncid, 'd12', 12, d12_id))

! Define NC variables
    call check( nf90_def_var(ncid, 'box_rmin', NF90_DOUBLE, box_rmin_id))
    call check( nf90_def_var(ncid, 'box_rmax', NF90_DOUBLE, box_rmax_id))
    call check( nf90_def_var(ncid, 'box_zmin', NF90_DOUBLE, box_zmin_id))
    call check( nf90_def_var(ncid, 'box_zmax', NF90_DOUBLE, box_zmax_id))
    call check( nf90_def_var(ncid, 'R', NF90_DOUBLE, [n_R_id], R_id))
    call check( nf90_def_var(ncid, 'Z', NF90_DOUBLE, [n_Z_id], Z_id))
    call check( nf90_def_var(ncid, 'psiN', NF90_DOUBLE, [n_R_id, N_Z_id], psiN_id))
    call check( nf90_def_var(ncid, 'gamma_array', NF90_DOUBLE, [n_R_id,N_Z_id,nspec_p1_id],&
                                  & gamma_array_id))
    call check( nf90_def_var(ncid, 'spec_name', NF90_CHAR, [d12_id,nspec_p1_id], spec_name_id))

! Put global attributes
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'RAYS_run_label', run_label))

! Finished NC definition
	call check( nf90_enddef(ncid))

	allocate(gamma(nspec_p1))
	allocate(gamma_array(N_pointsR_eq, N_pointsZ_eq, nspec_p1))

	dr = one/(N_pointsR_eq - 1)*(box_rmax-box_rmin)
	dz = one/(N_pointsZ_eq - 1)*(box_zmax-box_zmin)
	do i = 1, N_pointsR_eq ; do j = 1, N_pointsZ_eq
		R(i) = box_rmin + (i-1)*dr
		Z(j) = box_zmin + (j-1)*dz
		rvec = (/R(i), zero, Z(j)/)

		call axisym_toroid_psi(rvec, psi, gradpsi, psiN, gradpsiN)
		psiN_out(i,j) = psiN

		call equilibrium(rvec, eq)
		gamma(:) = eq%gamma(0:nspec)
		gamma_array(i, j, :) = gamma
	end do ; end do

! Put NC variables
    call check( nf90_put_var(ncid, box_rmin_id, box_rmin))
    call check( nf90_put_var(ncid, box_rmax_id, box_rmax))
    call check( nf90_put_var(ncid, box_zmin_id, box_zmin))
    call check( nf90_put_var(ncid, box_zmax_id, box_zmax))
    call check( nf90_put_var(ncid, R_id, R))
    call check( nf90_put_var(ncid, Z_id, Z))
    call check( nf90_put_var(ncid, psiN_id, psiN_out(:,:)))
    call check( nf90_put_var(ncid, gamma_array_id, gamma_array))
    call check( nf90_put_var(ncid, spec_name_id, spec_name(0:nspec)))

!   Close the NC file
    call check( nf90_close(ncid) )
    return
 end subroutine write_eq_contour_data_NC

!*************************************************************************

  subroutine check(status)
    use netcdf
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check


 end module axisym_toroid_processor_m

