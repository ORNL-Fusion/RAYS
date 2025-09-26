 module file_input_ray_init_m

! Generate the initial conditions for the rays
! i.e. intitial position for each ray: rvec0(1:3, iray) and
! information necessary to initialize k for each ray: kvec0(1:3, iray).
! Also generates a power weight for each ray.  It derives this data from an input
! namelist file "ray_init_<run_label>.in" where <run_label> must match the run_label
! variable in the RAYS input namelist file.
!

! This routine initializes based on initial launch position (X,Y,Z) and initial
! refractive indices (nx, ny, nz) from data in an input file.  Such data could represent
! any kind of launching structure - beam, gril, current, strap ...  The (nx, ny, nz)
! specify the ray initial direction, but don't have to constitute a unit vector, they
! are normalized later

! Not all input rays have to have valid dispersion solutions. For example if some initial
! locations are outside the allowed domain (e.g. outside the plasma) or if some rays are
! specified in a direction which is evanescent rather than propagating.  In that case
! the ray is skipped (rather than crashing) and the total number of valid initializations
! is taken for nray.

! There is a complication with respect to the power weighting.  If the user inputs a power
! weight distribution, but some don't initialize, the total weight will not sum to one. For
! now the quick-and-dirty solution is to renormalize by n_rays_in/nray.  We also allow a
! simple solution of uniform weight = 1/nray if all the weights are input as zero.


!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use constants_m, only : rkind

    implicit none

! Local data **************************************************

 !  File name for input
    character(len=80) :: in_filename

!_________________________________________________________________________________________
! Namelist data for /file_input_ray_init_list/
!_________________________________________________________________________________________

! Number of initial condition sets to be read in from namelist
    integer:: n_rays_in

! Initial positions and directions to be read in from namelist file.
! N.B. The ordering of indices in rvec_in(nray_max,3), rindex_vec_in(nray_max,3) is
!      opposite of that for rvec0(3, nray), rindex_vec0(3, nray)!!!  This is to make it
!      simpler to put data into the namelist file i.e. (X, Y, Z) on one line of the
!      namelist.

    real(KIND=rkind), allocatable :: rvec_in(:,:), rindex_vec_in(:,:)
    real(KIND=rkind), allocatable :: ray_pwr_wt_in(:)

 namelist /file_input_ray_init_list/ &
     & n_rays_in, rvec_in, rindex_vec_in, ray_pwr_wt_in

!_________________________________________________________________________________________


contains


!****************************************************************************


    subroutine file_input_ray_init(read_input, nray_max, nray, rvec0,&
             & rindex_vec0, ray_pwr_wt)

! Specify a range of initial launch positions (r,theta) and at each launch position
! an initial poloidal wave number, rindex_theta0, and toroidal wave number, rindex_phi0.
!
! ray_pwr_wt(i) = fraction of total power carried by ray i.  Should provide a ray weight
!                 subroutine as part of antenna model.  But for now al weights are just 1/nray.
!
! N.B. Some of the ray initializations may fail (e.g. initial point is outside plasma or
!      wave mode is evanescent).  This does not cause the program to stop.  It counts
!      the successful initializations and sets number of rays <=> nray, to that.
!
! N.B. Since we don't know nray until the end, allocate the output arrays rvec0,
!      rindex_vec0 at the end with proper length, nray!

! External procedures: Only from module use.

    use constants_m, only : rkind, zero, one, two
    use diagnostics_m, only: message_unit, messages_to_stdout,  message, text_message, &
                           & verbosity, run_label
    use rf_m, only : ray_dispersion_model, wave_mode, k0_sign, k0
    use one_ray_init_XYZ_k_direction_m, only : ray_init_XYZ_k_direction

    implicit none

    logical, intent(in) :: read_input

    integer, intent(in) :: nray_max
    integer, intent(out) :: nray
    real(KIND=rkind), allocatable, intent(out) :: rvec0(:, :), rindex_vec0(:, :)
    real(KIND=rkind), allocatable, intent(out) :: ray_pwr_wt(:)

    real(KIND=rkind), allocatable :: rvec_temp(:, :), rindex_vec_temp(:, :)
    real(KIND=rkind), allocatable :: ray_pwr_wt_temp(:)

    real(KIND=rkind) :: rvec(3), rindex_vec(3)
    integer :: count, i
!   Error return
	character(len=80) :: init_err = ''

 	integer :: input_unit, get_unit_number ! External, free unit finder

    call message(1)
    call text_message( 'Initializing file_input_ray_init ', 1)

! Allocate input arrays
! N.B. The ordering of indices in rvec_in(nray_max,3), rindex_vec_in(nray_max,3) is
!      opposite of that for rvec0(3, nray), rindex_vec0(3, nray)!!!
!		allocate ( rvec_in(nray_max,3), rindex_vec_in(nray_max,3), source = 0.0_rkind )
		allocate ( rvec_in(3, nray_max), rindex_vec_in(3, nray_max), source = 0.0_rkind )
		allocate ( ray_pwr_wt_in(nray_max), source = 1.0_rkind )
! Allocate temporary arrays
! N.B. The ordering of indices in rvec_temp(3, nray_max), rindex_vec_temp(3, nray_max) is
!      the same as for rvec0(3, nray), rindex_vec0(3, nray)!!!
		allocate ( rvec_temp(3, nray_max), rindex_vec_temp(3, nray_max), source = zero )
		allocate ( ray_pwr_wt_temp(nray_max), source = zero )

    if (read_input .eqv. .true.) then
! Read and write input namelist
		in_filename = 'ray_init_'//trim(run_label)//'.in'
		input_unit = get_unit_number()
		open(unit=input_unit, file=trim(in_filename),action='read', status='old', &
		   & form='formatted')
		read(input_unit, file_input_ray_init_list)
		close(unit=input_unit)
		if (verbosity >= 0) then
			write(message_unit, file_input_ray_init_list)
			if (messages_to_stdout) write(*, file_input_ray_init_list)
		end if
	end if

	if ((n_rays_in < 1) .or. (n_rays_in > nray_max)) then
		call message ('file_input_ray_init: improper number of rays  n_rays_in=', n_rays_in)
		write (*,*) 'file_input_ray_init: improper number of rays  n_rays_in=', n_rays_in
		stop 1
	end if

! N.B. Not all of these may successfully initialize because of errors.  So count successful
!      initializations.  That will be the final value of nray.
    count = 0

    ray_loop: do i = 1, n_rays_in

         rvec(:) = rvec_in(:,i)
         rindex_vec(:) = rindex_vec_in(:,i)
         write(*,*) 'rvec_in = ', rvec_in
!          write(*,*) 'rindex_vec_in = ', rindex_vec_in
        write(*,*) 'i = ', i, '   rvec_in = ', rvec_in(:,i)
        write(*,*) 'i = ', i, '   rvec = ', rvec

! Solve dispersion for magnitude of refractive index in direction rindex_vec
	   call ray_init_XYZ_k_direction(rvec, rindex_vec, init_err)

        if (trim(init_err) .ne. '') then
            if (verbosity > 0) then
				write(message_unit, *) 'ray_init_XYZ_k_direction failed for ray, i = ', i, &
				& '   init_err = ', trim(init_err), '  rvec = ', rvec, '   rindex_vec_in = ',&
				& rindex_vec_in(:,i)
            end if
            cycle ray_loop
        end if

        count = count +1
        rvec_temp(:, count) = rvec  ! N.B. The reversal of index order occurs here
        rindex_vec_temp(:, count) = rindex_vec

    end do ray_loop

    nray = count
    call message('axisym_toroid_ray_init: nray', nray, 1)

    if (nray == 0) then
        stop 'No successful ray initializations'
    end if

! Now that we know correct nray, allocate the output arrays

	allocate ( rvec0(3, nray), rindex_vec0(3, nray) )
	allocate ( ray_pwr_wt(nray) )

    rvec0 = rvec_temp(1:3, 1:nray)
    rindex_vec0 = rindex_vec_temp(1:3,1:nray)

    if (maxval(ray_pwr_wt_in) == zero) then
		ray_pwr_wt = one/nray ! Note: Q&D power model 1/nray
    else
		ray_pwr_wt = ray_pwr_wt_temp(1:nray)*n_rays_in/nray
	end if

	deallocate ( rvec_in, rindex_vec_in, ray_pwr_wt_in )
	deallocate ( rvec_temp, rindex_vec_temp, ray_pwr_wt_temp )

    end subroutine file_input_ray_init

!****************************************************************************

    subroutine deallocate_file_input_ray_init_m
! 		if (allocated(rvec0)) then
! 			deallocate ( rvec0, rindex_vec0)
! 			deallocate ( ray_pwr_wt)
! 		end if
		return ! Maybe nothing to deallocate.  rvec0 etc deallocated when
		       ! ray_init_axisym_toroid_R_Z_nphi_ntheta returns?
    end subroutine deallocate_file_input_ray_init_m

end module file_input_ray_init_m
