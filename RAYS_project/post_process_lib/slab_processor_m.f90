 module slab_processor_m
! Post processing for slab equilibrium, all variation in x direction

    use constants_m, only : rkind

    implicit none

!   Number and x locations of cyclotron resonances, 2nd harmonic resonances and hybrid
!   resonances found between x_min and x_max.  Allow for multiple resonances <= n_locs
    integer, parameter :: n_locs = 5
    integer :: n_ce_res, n_2ce_res, n_hybrid_res
    real(KIND=rkind), dimension(n_locs) :: x_ce_res, x_2ce_res, x_hybrid_res

!   Number and x locations of cutoffs
    integer :: n_P_cut, n_H_cut, n_det
    real(KIND=rkind), dimension(n_locs) :: x_P_cut, x_H_cut, x_det

! Alpha at hybrid cutoffs
    real(KIND=rkind), dimension(n_locs) :: alpha_e_H_cut
    real(KIND=rkind), dimension(n_locs) :: alpha_e_det

! Number of k vectors to plot for each ray in graphics
    integer :: num_plot_k_vectors

! Scale plot k vectors to kmax on the ray, (True, False) N.B. character not logical -> python
    character(len = 5) :: scale_k_vec = 'True'

! Base length for drawing k vectors as fraction of fig size
    real(KIND=rkind) :: k_vec_base_length

! Set plot xlim -> [xmin, xmax] and ylim -> [ymin,ymax], (True, False)
    character(len = 5) :: set_XY_lim = 'True'

! number of x points in equilibrium profile plots
	integer :: n_X = 51 ! Default can be reset on input

! Flags determining what to do besides plot rays.  Can be overridden (i.e. turned off)
! in namelist file
	logical :: calculate_dep_profiles = .true. ! Calculate all depositon profiles
	logical :: write_dep_profiles = .true. ! Write depositon profile netCDF file
	logical :: calculate_ray_diag = .true. ! Write detailed ray diagnostic netCDF file
	logical :: write_eq_X_profile_data = .true.  ! Write data for radial profiles netCDF file

    namelist /slab_processor_list/ num_plot_k_vectors, scale_k_vec,&
             & k_vec_base_length, set_XY_lim, &
             & calculate_dep_profiles, write_dep_profiles, calculate_ray_diag, &
             & write_eq_X_profile_data, n_X

 contains

 subroutine initialize_slab_processor(read_input)

    use diagnostics_m, only : message_unit, message, text_message
	use deposition_profiles_m, only : initialize_deposition_profiles

    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder

    if (read_input .eqv. .true.) then
    ! Read and write input namelist
   		input_unit = get_unit_number()
        open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
        read(input_unit, slab_processor_list)
        close(unit=input_unit)
        write(message_unit, slab_processor_list)
        call text_message('Finished initialize_slab_processor')
    end if

	if (calculate_dep_profiles .eqv. .true.) call initialize_deposition_profiles(read_input)

    return
 end subroutine initialize_slab_processor

!*************************************************************************

  subroutine slab_processor

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message
    use deposition_profiles_m, only : calculate_deposition_profiles, &
                        & write_deposition_profiles_LD

    implicit none

    call write_eq_profiles

    call write_kx_profiles

    call find_res_and_cuts

    call write_graphics_description_file

	if (calculate_ray_diag .eqv. .true.) call ray_detailed_diagnostics_slab

    if (calculate_dep_profiles .eqv. .true.) call calculate_deposition_profiles

    if (write_dep_profiles .eqv. .true.) call write_deposition_profiles_LD

    if (write_eq_X_profile_data .eqv. .true.) call write_eq_X_profile_data_NC

    call text_message('Finished slab_processor work')

    return
 end subroutine slab_processor

!*************************************************************************

  subroutine ray_detailed_diagnostics_slab()
! Takes the data from the ray_results_m module and for each ray extracts or calculates
! a set of data values for each step along the ray (e.g. ne, omega_ce/omega_rf, psi...).
! That data is then written to a netcdf file for analysis or plotting by a graphics code.
! There is a lot of data, and this subroutine is probably only useful for small numbers
! of rays.  The data can be plotted using plot_ray_diags.py which is located in RAYS/graphics_RAYS

! Working notes:
! (DBB 8/28/2025)  This thing is modified from the axisymmetric version ray_detailed_diagnostics

    use constants_m, only : rkind, e
    use diagnostics_m, only : integrate_eq_gradients, message, text_message, verbosity
    use species_m, only : nspec,ms
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
             & gamma_e_id, X_id, Y_id, Z_id, n_par_id, n_perp_id, P_absorbed_id, &
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
    real(kind=rkind), allocatable :: X(:,:)
    real(kind=rkind), allocatable :: Y(:,:)
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

    if (verbosity > 0) call text_message('Writing ray_detailed_diagnostics_slab')
    out_filename = 'ray_detailed_diagnostics_slab.'//trim(RAYS_run_label)//'.nc'

!   Allocate local arrays
!    allocate(npoints(number_of_rays))
    allocate(s(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(ne(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(Te_kev(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(modB(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(alpha_e(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(gamma_e(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(X(max_number_of_points,number_of_rays),source=0.0_rkind)
    allocate(Y(max_number_of_points,number_of_rays),source=0.0_rkind)
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

        	X(istep, iray) = rvec(1)
        	Y(istep, iray) = rvec(2)
       		Z(istep, iray) = rvec(3)

        	call equilibrium(rvec, eq)
        	Te_kev(istep, iray) = eq%Ts(0)/e/1000.0_rkind
        	modB(istep, iray) = eq%bmag
        	alpha_e(istep, iray) = eq%alpha(0)
        	gamma_e(istep, iray) = abs(eq%gamma(0))
        	ne(istep, iray) = eq%ns(0)

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
					   write(*,*) 'ray_detailed_diagnostics_slab: dispersion_model = ', &
					       & ray_dispersion_model
					   stop 'ray_detailed_diagnostics_slab: unimplemented ray_dispersion_model'
					end if

					if ( abs(dddw) > tiny(dddw) ) then
					   vg = -dddk / dddw
					else
					   write(*,*) 'ray_detailed_diagnostics_slab: infinite group velocity'
					   stop 'ray_detailed_diagnostics_slab: infinite group velocity'
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
    call check( nf90_def_var(ncid, 'X', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], X_id))
    call check( nf90_def_var(ncid, 'Y', NF90_DOUBLE, [max_number_of_points_id,number_of_rays_id], Y_id))
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
    call check( nf90_put_var(ncid, X_id, X))
    call check( nf90_put_var(ncid, Y_id, Y))
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

  end subroutine ray_detailed_diagnostics_slab

!*************************************************************************

  subroutine find_res_and_cuts
! Finds x locations of cold plasma resonances and cutoffs.  Now only finds electron cyclotron
! resonances.  Later include arbitrary ion species.  Simple scan algorithm, no root finding.
! Maybe later.
!
! Nota Bene: B is assumed in y-z plane and ny, nz are constant because of x stratification.
!
! Restriction:  This assumes that the component of k in the y-x plane is parallel to B for
! all x.  Equivalent to B along z (no shear) and ky = 0.  If not then there is
! component of k_perp in the y-z plane, nz is not n_parallel and the cutoff (kx = 0) is in
! a different place.  Could be generalized.  Maybe later.

    use constants_m, only : rkind
    use diagnostics_m, only : run_label
    use equilibrium_m, only : equilibrium, eq_point
    use slab_eq_m, only : xmin, xmax
    use suscep_m, only : dielectric_cold
    use ray_init_m, only : nray, rindex_vec0

    implicit none

	integer :: res_and_cut_unit, get_unit_number ! External, free unit finder

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

    integer, parameter :: n_xpoints = 1000 ! Number of x points in scan

    real(KIND=rkind) :: rvec(3)
    real(KIND=rkind) :: v_ce_0, v_ce_1, v_2ce_0, v_2ce_1, v_hybrid_0, v_hybrid_1
    real(KIND=rkind) :: v_P_cut_0, v_P_cut_1,  v_H_cut_0, v_H_cut_1
    real(KIND=rkind) :: a, b, c, v_det_0, v_det_1 ! parameters for calculating determinant
    integer :: iray, ix
    real(KIND=rkind) :: nz, x, dx
    complex :: vH

    complex(KIND=rkind) :: eps_cold(3,3)

   res_and_cut_unit = get_unit_number()
   open(unit = res_and_cut_unit, file = 'res_and_cut.'//trim(run_label), &
            & action='write', status='replace', form='formatted')

    dx = (xmax - xmin)/(n_xpoints-1)

    ray_loop: do iray = 1, nray
        write(*,*) ' '
        write(*,*) 'ray ', iray

        write(res_and_cut_unit,*) ' '
        write(res_and_cut_unit,*) 'ray ', iray

        n_ce_res = 0
        n_2ce_res = 0
        n_hybrid_res = 0
        n_P_cut = 0
        n_H_cut = 0
        n_det = 0

        x_ce_res = 0.
        x_2ce_res = 0.
        x_hybrid_res = 0.
        x_P_cut = 0.
        x_H_cut = 0.
        x_det = 0.

        rvec = (/real(0., KIND=rkind), real(0., KIND=rkind), real(0., KIND=rkind)/)
       x_loop: do ix = 0, n_xpoints-1
            x = xmin + ix*dx
            rvec(1) = real(x, KIND=rkind)

            call equilibrium(rvec, eq )
            nz = sum(rindex_vec0(:, iray)*eq%bunit(:))
            v_ce_1 = eq%gamma(0)+1. ! remember gamma(0) is negative
            v_2ce_1 = eq%gamma(0)+0.5

            call dielectric_cold(eq, eps_cold)
            v_hybrid_1 = eps_cold(1,1)%re
            v_P_cut_1 = eps_cold(3,3)%re
            v_H_cut_1 = real(( eps_cold(1,1)**2+eps_cold(1,2)**2 - 2.*eps_cold(1,1)*nz**2 &
                    & + nz**4 ), KIND=rkind)

           a = eps_cold(1,1)%re
           b = real( -(eps_cold(1,1)**2+eps_cold(1,2)**2+eps_cold(1,1)*eps_cold(3,3)) &
              & + (eps_cold(1,1)+eps_cold(3,3))*nz**2, KIND=rkind )
           c = real( ( eps_cold(1,1)**2+eps_cold(1,2)**2 - 2.*eps_cold(1,1)*nz**2 &
              & + nz**4 ) *  eps_cold(3,3), KIND=rkind )

           v_det_1 = b**2-4.*a*c

!             vH = eps_cold(1,1)**2+eps_cold(1,2)**2 - 2.*eps_cold(1,1)*nz**2 + nz**4
!             write(*,*) 'x = ', x, ' eps_cold(1,1 ) = ', eps_cold(1,1), ' v_P_cut_1 = ',&
!              & v_P_cut_1, 'v_H = ', vH

            if (ix == 0) then ! first time initialize v_ce0 etc and skip to next point
                v_ce_0 = v_ce_1
                v_2ce_0 = v_2ce_1
                v_hybrid_0 = v_hybrid_1
                v_P_cut_0 = v_P_cut_1
                v_H_cut_0 = v_H_cut_1
                v_det_0 = v_det_1
                cycle x_loop
            end if

            ! Check for sign change between v_0 and v_1
            if (v_ce_0 * v_ce_1 <= 0.) then
              n_ce_res = n_ce_res + 1
              x_ce_res(n_ce_res) = x
            end if
            if (v_2ce_0 * v_2ce_1 <= 0.) then
              n_2ce_res = n_2ce_res + 1
              x_2ce_res(n_2ce_res) = x
            end if
            if (v_hybrid_0 * v_hybrid_1 <= 0.) then
              n_hybrid_res = n_hybrid_res + 1
              x_hybrid_res(n_hybrid_res) = x
            end if
            if (v_P_cut_0 * v_P_cut_1 <= 0.) then
              n_P_cut = n_P_cut + 1
              x_P_cut(n_P_cut) = x
            end if
            if (v_H_cut_0 * v_H_cut_1 <= 0.) then
              n_H_cut = n_H_cut + 1
              x_H_cut(n_H_cut) = x
              alpha_e_H_cut(n_H_cut) = eq%alpha(0)
            end if
            if (v_det_0 * v_det_1 <= 0.) then
              n_det = n_det + 1
              x_det(n_det) = x
              alpha_e_det(n_det) = eq%alpha(0)
            end if

            ! advance v_0
            v_ce_0 = v_ce_1
            v_2ce_0 = v_2ce_1
            v_hybrid_0 = v_hybrid_1
            v_P_cut_0 = v_P_cut_1
            v_H_cut_0 = v_H_cut_1
            v_det_0 = v_det_1

        end do x_loop

        write(*,*) 'x_ce_res = ', x_ce_res
        write(*,*) 'x_2ce_res = ', x_2ce_res
        write(*,*) 'x_hybrid_res = ', x_hybrid_res
        write(*,*) 'x_P_cut = ', x_P_cut
        write(*,*) 'x_H_cut = ', x_H_cut
        write(*,*) 'alpha_e_H_cut = ', alpha_e_H_cut
        write(*,*) 'x_det = ', x_det
        write(*,*) 'alpha_e_det = ', alpha_e_det

        write(res_and_cut_unit,*) 'x_ce_res = ', x_ce_res
        write(res_and_cut_unit,*) 'x_2ce_res = ', x_2ce_res
        write(res_and_cut_unit,*) 'x_hybrid_res = ', x_hybrid_res
        write(res_and_cut_unit,*) 'x_P_cut = ', x_P_cut
        write(res_and_cut_unit,*) 'x_H_cut = ', x_H_cut
        write(res_and_cut_unit,*) 'alpha_e_H_cut = ', alpha_e_H_cut
        write(res_and_cut_unit,*) 'x_det = ', x_det
        write(res_and_cut_unit,*) 'alpha_e_det = ', alpha_e_det

    end do ray_loop
    close(unit = res_and_cut_unit)

  end subroutine find_res_and_cuts

!*************************************************************************

  subroutine write_eq_profiles
! Writes equilibrium profiles and other stuff for plotting
! For ease of reading these are re-cast to single precision before writing

    use constants_m, only : rkind, skind, pi
    use diagnostics_m, only : run_label
    use equilibrium_m, only : equilibrium, eq_point
    use slab_eq_m, only : xmin, xmax
    implicit none

	integer :: eq_profile_unit, get_unit_number ! External, free unit finder

!   fast, slow etc
    character(len=60) :: wave_mode

    integer, parameter :: n_xpoints = 101 ! Number of x points in scan
    integer :: ix
    real(KIND=rkind) :: rvec(3)
    real(KIND=rkind) :: x, dx
    real(KIND=skind), dimension(9) :: profile_vec
    character(len = 9), dimension(9) :: prof_name
    character(len = 10), parameter  :: b10 = '          '
    character(len = 9), parameter  ::   b9 = '         '
    character(len = 8), parameter  ::   b8 = '        '
    character(len = 7), parameter  ::   b7 = '       '

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

    write(*,*)  'Start writing eq profile vectors'

    prof_name = (/'x        ', 'ne       ', 'bmag     ', 'fc_e     ', 'fp_e     ', &
                       & 'fp_e/f_rf', 'fc_e/f_rf', 'alpha    ', 'beta     '/)


    eq_profile_unit = get_unit_number()
    open(unit = eq_profile_unit, file = 'eq_profiles_slab.'//trim(run_label))

    write(eq_profile_unit,*) ' ',prof_name(1), b9, prof_name(2), b8, prof_name(3), b8,&
                     & prof_name(4), b8, prof_name(5), b7, prof_name(6), b8, prof_name(7),&
                     & b8, prof_name(8), b8, prof_name(9)


    rvec = (/real(0., KIND=rkind), real(0., KIND=rkind), real(0., KIND=rkind)/)
    dx = (xmax - xmin)/(n_xpoints-1)
    x_loop: do ix = 0, n_xpoints-1

        x = xmin + ix*dx
        rvec(1) = real(x, KIND=rkind)

        call equilibrium(rvec, eq)

        profile_vec(1) = real(x, KIND=skind)
        profile_vec(2) = real(eq%ns(0), KIND=skind) ! electron density
        profile_vec(3) = real(eq%bmag, KIND=skind) ! total B field
        profile_vec(4) = real(abs(eq%omgc(0))/(2.*pi), KIND=skind) ! electron cyclotron frequency
        profile_vec(5) = real(sqrt(eq%omgp2(0))/(2.*pi), KIND=skind) ! electron plasma frequency
        profile_vec(6) = real(sqrt(eq%alpha(0)), KIND=skind) ! f_pe/f_rf
        profile_vec(7) = real(abs(eq%gamma(0)), KIND=skind) ! f_ce/f_rf
        profile_vec(8) = real(eq%alpha(0), KIND=skind) ! (f_pe/f_rf)^2
        profile_vec(9) = real(eq%gamma(0)**2, KIND=skind) ! (f_ce/f_rf)^2

        write(eq_profile_unit,*) profile_vec

    end do x_loop

    write(*,*)  'Finished writing eq profile vectors'
    close(eq_profile_unit)

  end subroutine write_eq_profiles

!*************************************************************************

 subroutine write_eq_X_profile_data_NC
! Write a netcdf file ('eq_radial.<run_label>.nc') containing data needed to plot
! radial profiles for equilibrium quantites like psiN, ne, Te, etc as functions of
! psi (and in the near future rho = sqrt(toroidal flux))
!
! For now limit profiles to inside boundary 0 <= psiN <= 1.  Maybe later extend to outside.

    use constants_m, only : rkind, one, zero, e
    use diagnostics_m, only : message_unit, message, text_message, verbosity, run_label
    use species_m, only : nspec, spec_name
	use equilibrium_m, only : eq_point, equilibrium, write_eq_point
    use slab_eq_m, only: xmin, xmax
    use ray_results_m, only : date_vector, RAYS_run_label
    use netcdf
	use XY_curves_netCDF_m, only : XY_curve_netCDF, write_XY_curves_netCDF

    implicit none

! Number of profiles to generate and list of profiles
    integer, parameter :: n_profiles = 6
	type(XY_curve_netCDF) :: profile_list(n_profiles)
	integer :: n_grid(n_profiles)

!   Declare local variables
    real(KIND=rkind) :: X(n_X)
    real(KIND=rkind) :: ne(n_X), Te_ev(n_X), Ti_ev(n_X), alpha_e(n_X), gamma_e(n_X),&
                    &   modB_tesla(n_X)

	integer :: i, ix, ierr
    real(KIND=rkind) :: rvec(3), dx

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

 !  File name for  output
    character(len=80) :: out_filename

    real(KIND=rkind) :: ns(0:nspec)
    real(KIND=rkind) :: ts(0:nspec)
    real(KIND=rkind) :: alpha(0:nspec)
    real(KIND=rkind) :: gamma(0:nspec)
    real(KIND=rkind) :: bmag
    character(len=60) :: equib_err

   if (verbosity > 0) call text_message('Writing plasma equilibrium profiles on radial psi grid')

	rvec = (/zero, zero, zero/)
	dx = (xmax - xmin)/(n_X-1)
	do i = 1, n_X
        X(i) = xmin + (i-1)*dx
        rvec(1) = X(i)
		call equilibrium(rvec, eq)
		ne(i) = eq%ns(0)
		Te_ev(i) = eq%Ts(0)/e ! Convert from Joules to ev
		Ti_ev(i) = eq%Ts(1)/e ! N.B. For now all ions are assumed to have the same Ti profile
		modB_tesla(i) = eq%bmag
		alpha_e(i) = eq%alpha(0)
		gamma_e(i) = eq%gamma(0)
	end do

!  write(*,*) " "
!  write(*,*) "X = ", XZ
!  write(*,*) " "
!  write(*,*) "ne =  ", ne
!  write(*,*) " "
!  write(*,*) "Te =   ", Te
!  write(*,*) " "
!  write(*,*) "Ti =   ", Ti
!  write(*,*) " "
!  write(*,*) "modB_tesla =  ", modB_tesla


! Load ne data into profile_list
	profile_list(1)%grid_name = 'X'
	profile_list(1)%curve_name = 'ne(X)'
	n_grid(1) = n_X
	allocate(profile_list(1)%grid(n_grid(1)), source = zero)
	profile_list(1)%grid(:) = X(:)
	allocate(profile_list(1)%curve(n_grid(1)), source = zero)
	profile_list(1)%curve(:) = ne(:)

! Load Te data into profile_list
	profile_list(2)%grid_name = 'X'
	profile_list(2)%curve_name = 'Te(X)'
	n_grid(2) = n_X
	allocate(profile_list(2)%grid(n_grid(2)), source = zero)
	profile_list(2)%grid(:) = X(:)
	allocate(profile_list(2)%curve(n_grid(2)), source = zero)
	profile_list(2)%curve(:) = Te_ev(:)

! Load Ti data into profile_list
	profile_list((3))%grid_name = 'X'
	profile_list((3))%curve_name = 'Ti(X)'
	n_grid((3)) = n_X
	allocate(profile_list((3))%grid(n_grid((3))), source = zero)
	profile_list((3))%grid(:) = X(:)
	allocate(profile_list((3))%curve(n_grid((3))), source = zero)
	profile_list((3))%curve(:) = Ti_ev(:)

! Load modB data into profile_list
	profile_list((4))%grid_name = 'X'
	profile_list((4))%curve_name = 'midB (Tesla)'
	n_grid((4)) = n_X
	allocate(profile_list((4))%grid(n_grid((4))), source = zero)
	profile_list((4))%grid(:) = X(:)
	allocate(profile_list((4))%curve(n_grid((4))), source = zero)
	profile_list((4))%curve(:) = modB_tesla(:)

! Load alpha_e data into profile_list
	profile_list((5))%grid_name = 'X'
	profile_list((5))%curve_name = 'alpha_e(X)'
	n_grid((5)) = n_X
	allocate(profile_list((5))%grid(n_grid((5))), source = zero)
	profile_list((5))%grid(:) = X(:)
	allocate(profile_list((5))%curve(n_grid((5))), source = zero)
	profile_list((5))%curve(:) = alpha_e(:)

! Load gamma_e data into profile_list
	profile_list((6))%grid_name = 'X'
	profile_list((6))%curve_name = 'gamma_e(X)'
	n_grid((6)) = n_X
	allocate(profile_list((6))%grid(n_grid((6))), source = zero)
	profile_list((6))%grid(:) = X(:)
	allocate(profile_list((6))%curve(n_grid((6))), source = zero)
	profile_list((6))%curve(:) = gamma_e(:)


! Generate netCDF file
	out_filename = 'eq_X_profiles.'//trim(run_label)
	call write_XY_curves_netCDF(profile_list, out_filename)

    return
 end subroutine write_eq_X_profile_data_NC

!*************************************************************************

  subroutine write_kx_profiles
! Writes kx roots versus x.
! For ease of reading these are re-cast to single precision before writing

    use constants_m, only : rkind, skind, pi
    use diagnostics_m, only : run_label
    use equilibrium_m, only : equilibrium, eq_point
    use slab_eq_m, only : xmin, xmax
    use ray_init_m, only : nray, rindex_vec0
    use rf_m, only : ray_dispersion_model, k0
    use dispersion_solvers_m, only : solve_nx_vs_ny_nz_by_bz

    implicit none

	integer :: kx_profile_unit, get_unit_number ! External, free unit finder

!   fast, slow etc
    character(len=60) :: wave_mode

    integer, parameter :: n_xpoints = 101 ! Number of x points in scan
    integer :: iray, ix
    real(KIND=rkind) :: x, dx
    real(KIND=rkind) :: rvec(3)
    real(KIND=rkind) :: ny, nz
    real(KIND=skind), dimension(9) :: profile_vec
    character(len = 17), dimension(9) :: profile_name_vec
    complex(KIND=rkind) :: nx
    complex(KIND=skind) :: nx_sngl
    real(KIND=skind) :: ny_sngl, nz_sngl

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

    write(*,*)  'Start writing kx profile vectors'

    profile_name_vec = (/'x                ', 'kx_real_plus     ', 'kx_im_plus       ',&
                         'kx_real_minus    ', 'kx_im_minus      ', 'kx_real_fast     ',&
                         'kx_im_fast       ', 'kx_real_slow     ', 'kx_im_slow       '/)

    kx_profile_unit = get_unit_number()
    open(unit = kx_profile_unit, file = 'kx_profiles_slab.'//trim(run_label))

    ray_loop: do iray = 1, nray
        ny = rindex_vec0(2, iray)
        nz = rindex_vec0(3, iray)
        ny_sngl = real(ny, KIND=skind)
        nz_sngl = real(nz, KIND=skind)

        write(*,*) 'ray ', iray, ' ny ', ny_sngl, ' nz ', nz_sngl
        write(kx_profile_unit,*) 'ray ', iray, ' ny ', ny_sngl, ' nz ', nz_sngl
        write(kx_profile_unit,*) profile_name_vec

     rvec = (/real(0., KIND=rkind), real(0., KIND=rkind), real(0., KIND=rkind)/)
     dx = (xmax - xmin)/(n_xpoints-1)
     x_loop: do ix = 0, n_xpoints-1

        x = xmin + ix*dx
        rvec(1) = real(x, KIND=rkind)

        call equilibrium(rvec, eq)

            profile_vec(1) = real(x, KIND=skind)

            ! kx vs x for fast and slow cold plasma roots.  Always use +1 for k0_sign
            wave_mode = 'plus'
            call solve_nx_vs_ny_nz_by_bz(eq, ray_dispersion_model, wave_mode, +1, ny, nz, nx)
            nx_sngl = cmplx(nx, KIND=skind)
            profile_vec(2) = real(k0, KIND=skind)*real(nx_sngl)
            profile_vec(3) = real(k0, KIND=skind)*aimag(nx_sngl)

            wave_mode = 'minus'
            call solve_nx_vs_ny_nz_by_bz(eq, ray_dispersion_model, wave_mode, +1, ny, nz, nx)
            nx_sngl = cmplx(nx, KIND=skind)
            profile_vec(4) = real(k0, KIND=skind)*real(nx_sngl)
            profile_vec(5) = real(k0, KIND=skind)*aimag(nx_sngl)
            ! kx vs x for fast and slow cold plasma roots

            wave_mode = 'fast'
            call solve_nx_vs_ny_nz_by_bz(eq, ray_dispersion_model, wave_mode, +1, ny, nz, nx)
            nx_sngl = cmplx(nx, KIND=skind)
            profile_vec(6) = real(k0, KIND=skind)*real(nx_sngl)
            profile_vec(7) = real(k0, KIND=skind)*aimag(nx_sngl)

            wave_mode = 'slow'
            call solve_nx_vs_ny_nz_by_bz(eq, ray_dispersion_model, wave_mode, +1, ny, nz, nx)
            nx_sngl = cmplx(nx, KIND=skind)
            profile_vec(8) = real(k0, KIND=skind)*real(nx_sngl)
            profile_vec(9) = real(k0, KIND=skind)*aimag(nx_sngl)

            write(kx_profile_unit,*) profile_vec

        end do x_loop

    end do ray_loop

    write(*,*)  'Finished writing kx profile vectors'
    close(kx_profile_unit)

  end subroutine write_kx_profiles

!*************************************************************************

!   subroutine calculate_depositon_profiles
! ! Calculate various deposition profiles (e.g. power deposition for each species or
! ! summed over species)
!
!
!  end subroutine calculate_depositon_profiles

!*************************************************************************

  subroutine write_graphics_description_file
! Info the graphics routine needs to plot ray trajectories.

   use diagnostics_m, only : run_description, run_label
   use slab_eq_m, only : xmin, xmax, ymin, ymax, zmin, zmax

   integer :: graphics_descrip_unit, get_unit_number ! External, free unit finder

   open(unit = graphics_descrip_unit, file = 'graphics_description_slab.dat')

   write(graphics_descrip_unit, *) 'run_description = ', run_description
   write(graphics_descrip_unit, *) 'run_label = ', run_label

   write(graphics_descrip_unit, *) 'xmin = ', xmin
   write(graphics_descrip_unit, *) 'xmax = ', xmax
   write(graphics_descrip_unit, *) 'ymin = ', ymin
   write(graphics_descrip_unit, *) 'ymax = ', ymax
   write(graphics_descrip_unit, *) 'zmin = ', zmin
   write(graphics_descrip_unit, *) 'zmax = ', zmax

   write(graphics_descrip_unit, *) 'num_plot_k_vectors = ', num_plot_k_vectors
   write(graphics_descrip_unit, *) 'scale_k_vec = ', trim(scale_k_vec)
   write(graphics_descrip_unit, *) 'k_vec_base_length = ', k_vec_base_length
   write(graphics_descrip_unit, *) 'set_XY_lim = ', trim(set_XY_lim)

   close(unit = graphics_descrip_unit)

  end subroutine write_graphics_description_file


!*************************************************************************

  subroutine check(status)
    use netcdf
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check

 end module slab_processor_m

