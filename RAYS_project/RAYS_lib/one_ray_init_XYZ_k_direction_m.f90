 module one_ray_init_XYZ_k_direction_m

! This module initializes based on initial launch position (X,Y,Z) and initial
! direction of wave vector (nX, nY, nZ))

! The main purpose of the module is as a wrapper for the contained subroutine
! ray_init_XYZ_k_direction(), which returns a refractive index vector, n_vec, that is
! a solution of the dispersion relation in the direction of k evaluated at position (X,Y,Z).
! The subroutine is used to initialize multiple rays inside some higher launcher model, such
! as for a cone or beam.  The module also contains subroutine one_ray_init_XYZ_k_direction(),
! which initializes a single ray from data in the namelist. This simple model could be
! used for testing.

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use constants_m, only : rkind, zero, one

    implicit none

!_________________________________________________________________________________________
! Namelist data for /one_ray_init_XYZ_k_direction_list/
!_________________________________________________________________________________________

! Data for initial launch position
    real(KIND=rkind) ::  X = zero, Y = zero, Z = zero
! Data for initial launch direction
    real(KIND=rkind) ::  nX = zero, nY = zero, nZ = zero
! Switch to turn off solution of disp rel.  Just use the input nX,nY,nZ as is
	logical :: use_this_n_vec = .false.

 namelist /one_ray_init_XYZ_k_direction_list/ X, Y, Z, nX, nY, nZ, use_this_n_vec

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

	subroutine one_ray_init_XYZ_n_direction(nray_max, nray, rvec0, rindex_vec0, ray_pwr_wt)

! N.B. The allocations and Ray_loop below could be avoided since we already know that
!      there is at most one successful ray initialization. But to keep this routine
!      parallel to the other ray initialization routines I'm keeping them.

    use diagnostics_m, only: message_unit, messages_to_stdout,  message, text_message, verbosity
    use equilibrium_m, only : equilibrium, eq_point, write_eq_point
    use dispersion_solvers_m, only: solve_n_vs_theta
    use rf_m, only : ray_dispersion_model, wave_mode, k0_sign, k0

    implicit none

    integer, intent(in) :: nray_max
    integer, intent(out) :: nray
    real(KIND=rkind), allocatable, intent(out) :: rvec0(:, :), rindex_vec0(:, :)
    real(KIND=rkind), allocatable, intent(out) :: ray_pwr_wt(:)
!   Error return
	character(len=80) :: init_err = ''

    type(eq_point) :: eq
 	integer :: input_unit, get_unit_number ! External, free unit finder
    integer :: i_R, count
    real(KIND=rkind) :: rvec(3), rindex_vec(3), cos_theta, theta
    complex(KIND=rkind) :: n_cmplx

    call message(1)
    call text_message( 'Initializing ray_init_axisym_toroid_nphi_ntheta ', 1)

! Read and write input namelist
  	input_unit = get_unit_number()
    open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
    read(input_unit, one_ray_init_XYZ_k_direction_list)
    close(unit=input_unit)
    if (verbosity >= 0) then
		write(message_unit, one_ray_init_XYZ_k_direction_list)
		if (messages_to_stdout) write(*, one_ray_init_XYZ_k_direction_list)
    end if
! Allocate maximum space for the initial condition vectors rvec, rindex_vec
! N.B. Not all of these may successfully initialize because of errors.  So count successful
!      initializations.  That's the final value of nray.

	nray = 1

	if ((nray < 1) .or. (nray > nray_max)) then
		call message ('one_ray_init_XYZ_n_direction: improper number of rays  nray=', nray)
		write (*,*) 'one_ray_init_XYZ_n_direction: improper number of rays  nray=', nray
		stop 1
	end if


! Now that we know correct nray, allocate the output arrays
	allocate ( rvec0(3, nray), rindex_vec0(3, nray) )
	allocate ( ray_pwr_wt(nray) )

    count = 0

	rvec = (/ X, Y, Z /)
	rindex_vec = (/ nX, nY, nZ /)

! No disp rel solution. Just use the input rvec, rindex_vec as is
	if (use_this_n_vec) then
		rvec0(:, 1) = rvec
		rindex_vec0(:, 1) = rindex_vec
		ray_pwr_wt(1) = one
		return
	end if

    Ray_loop: do i_R = 1, 1
		call ray_init_XYZ_k_direction(rvec, rindex_vec, init_err)
		if (trim(init_err) /= '') cycle Ray_loop
        count = count + 1
        rvec0( : , count) = rvec
        rindex_vec0( : , count) = rindex_vec
		ray_pwr_wt(count) = one
	end do Ray_loop

    nray = count
    call message('one_ray_init_XYZ_n_direction: nray', nray, 1)

    if (nray == 0) then
        stop 'No successful ray initializations'
    end if


    end subroutine one_ray_init_XYZ_n_direction
!****************************************************************************

    subroutine ray_init_XYZ_k_direction(rvec, nvec, init_err)

    use diagnostics_m, only: message_unit, message, text_message, verbosity
    use equilibrium_m, only : equilibrium, eq_point
    use dispersion_solvers_m, only: solve_n_vs_theta
    use rf_m, only : ray_dispersion_model, wave_mode, k0_sign

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(inout) :: nvec(3)
!   Error return
	character(len=80) :: init_err
    real(KIND=rkind) :: nvec0(3) ! Copy of input nvec.  Used for diagnostic
    type(eq_point) :: eq

    real(KIND=rkind) :: cos_theta, theta
    complex(KIND=rkind) :: n_cmplx

	nvec0 = nvec
	call equilibrium(rvec, eq)
	if (trim(eq%equib_err) /= '') then
	   call text_message('ray_init_axisym_toroid_nphi_ntheta: equib_err = ', &
						& trim(eq%equib_err), 1)
	   init_err = trim(eq%equib_err)
	   return
	end if

	! unit refractive index vector
	nvec = nvec/sqrt(dot_product(nvec,nvec))

	cos_theta = dot_product(eq%bunit, nvec) ! Parallel component
	theta = acos(cos_theta)

! Solve dispersion for complex refractive index in nvec direction
	call solve_n_vs_theta(eq, ray_dispersion_model, wave_mode, k0_sign, theta, n_cmplx)

! Test for evanescence
	if (abs(n_cmplx%im) > 10.*tiny(abs(n_cmplx))) then
		if (verbosity > 0) then
			write(message_unit, *) 'ray_init_XYZ_k_direction: non-zero Im(n_cmplx),&
			   & rvec = ',  rvec, ' nvec0 = ', nvec0, ' n_cmplx = ', n_cmplx
		end if
		init_err = 'evanescent'
		return
	end if

	nvec = n_cmplx%re*nvec  ! i.e. magnitude times unit vector

    end subroutine ray_init_XYZ_k_direction
!****************************************************************************

    subroutine deallocate_ray_init_XYZ_k_direction_m
! 		if (allocated(rvec0)) then
! 			deallocate ( rvec0, rindex_vec0)
! 			deallocate ( ray_pwr_wt)
! 		end if
		return ! Maybe nothing to deallocate.  rvec0 etc deallocated when
		       ! ray_init_XYZ_k_direction returns?
    end subroutine deallocate_ray_init_XYZ_k_direction_m

end module one_ray_init_XYZ_k_direction_m
