module  eqdsk_magnetics_spline_interp_m
! A simple  eqdsk equilibrium model. This code uses eqdsk routines in module eqdsk_utilities_m.f90
! These were adapted from similar codes by Richard Fitzpatrick.
!
! N.B. Values of Psi are shifted on initialization so that Psi is zero on axis.
! N.B  In eqdsk radial profiles (e.g. Q(psi) and RBphi(psi)) are given on a uniform grid
!      in poloidal flux.  It uses the same grid dimension as R -> NW.
!      For the splined profiles below the grid is PsiN -> [0,1]

! Working notes
! DBB (1/27/2025) Added stuff to get splined profiles of Q of psi and RBhi of psi.
! DBB (8/2/2024) Changed error return 'if psi > 1' to 'if psi > plasma_psi_limit' to
! allow rays outside last closed flux surface. plasma_psi_limit defaults to 1.0 but can be
! reset in namelist
!_________________________________________________________________________________________

    use constants_m, only : rkind
    use quick_cube_splines_m, only : cube_spline_function_1D, cube_spline_function_2D

    implicit none

! Local data **************************************************

! Stuff for 2D spline profiles i.e. psi

    type(cube_spline_function_2D) :: Psi_profile
 	character (len = 80) :: Psi_name = 'Psi_profile'

! Stuff for 1D spline profiles versus PsiN

    type(cube_spline_function_1D) :: T_profile
 	character (len = 80) :: T_name = 'T_profile'

    type(cube_spline_function_1D) :: Q_profile
 	character (len = 80) ::Q_name = 'Q_profile'

    type(cube_spline_function_1D) :: rho_profile
 	character (len = 80) ::rho_name = 'rho_profile'

    type(cube_spline_function_1D) :: Tflux_profile
 	character (len = 80) ::Tflux_name = 'Tflux_profile'

! Stuff for 1D spline profiles versus rho

    type(cube_spline_function_1D) :: psiN_profile_rho
 	character (len = 80) :: psiN_name = 'psiN_profile_rho'

    ! Flux function psi at plasma boundary
    real(KIND=rkind) :: psiB

! Namelist data for /eqdsk_magnetics_spline_interp_list/  *****************************

! Name of input eqdsk file
    character (len = 100) :: eqdsk_file_name

  namelist / eqdsk_magnetics_spline_interp_list/ eqdsk_file_name

!********************************************************************

contains

!********************************************************************

  subroutine initialize_eqdsk_magnetics_spline_interp(read_input, r_axis, z_axis, &
               & box_rmin, box_rmax, box_zmin, box_zmax, &
               & inner_bound, outer_bound, upper_bound, lower_bound)

    use constants_m, only : one
    use species_m, only : nspec
    use diagnostics_m, only : message, message_unit,messages_to_stdout, verbosity

    use eqdsk_utilities_m, only : ReadgFile, WritegFile, R_grid, Z_grid, dR, dZ, &
        & string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, &
        & PSIAXIS, PSIBOUND, B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, &
        & ZBOUND, RLIM, ZLIM

    implicit none

    logical, intent(in) :: read_input
 	integer :: input_unit, get_unit_number ! External, free unit finder

! Geometry data
    ! Magnetic axis
    real(KIND=rkind), intent(out) :: r_axis, z_axis
    ! data for bounding box of computational domain
    real(KIND=rkind), intent(out) :: box_rmin, box_rmax, box_zmin, box_zmax
    ! data for plasma boundary
    real(KIND=rkind), intent(out) :: inner_bound, outer_bound, upper_bound, lower_bound

! Local arrays
    real(KIND=rkind), dimension (:), allocatable :: PsiN_grid
    real(KIND=rkind), dimension (:), allocatable :: Tflux_on_psiNgrid
    real(KIND=rkind), dimension (:), allocatable :: rho_on_psiNgrid
    real(KIND=rkind), dimension (:), allocatable :: psiN_on_rhogrid

    integer :: i, j, nwk

    if (verbosity >= 0) then
		write(*,*) 'initialize_eqdsk_magnetics_spline_interp'
    end if

    if (read_input .eqv. .true.) then
  		input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, eqdsk_magnetics_spline_interp_list)
        close(unit=input_unit)
    end if

! Write input namelist
    if (verbosity >= 0) then
		write(message_unit, eqdsk_magnetics_spline_interp_list)
		if (messages_to_stdout) write(*, eqdsk_magnetics_spline_interp_list)
		call message(1)
    end if

    call ReadgFile(eqdsk_file_name)

	r_axis = RAXIS
	z_axis = ZAXIS

	box_rmin = RBOXLFT
	box_rmax = box_rmin + RBOXLEN
	box_zmin = ZOFF - ZBOXLEN/2.
	box_zmax = ZOFF + ZBOXLEN/2.

    inner_bound = minval(RBOUND)
    outer_bound = maxval(RBOUND)
    lower_bound = minval(ZBOUND)
    upper_bound = maxval(ZBOUND)

    if (verbosity >= 0) then
		call message('r_axis = ', r_axis)
		call message('z_axis = ', z_axis)
		call message('Inner boundary = ', inner_bound)
		call message('Outer boundary = ', outer_bound)
		call message('Upper boundary = ', upper_bound)
		call message('Lower boundary = ', lower_bound)

		write(*,*) 'r_axis = ', r_axis
		write(*,*) 'z_axis = ', z_axis
		write(*,*) 'Inner boundary = ', inner_bound
		write(*,*) 'Outer boundary = ', outer_bound
		write(*,*) 'Lower boundary = ', lower_bound
		write(*,*) 'Upper boundary = ', upper_bound
    end if

! Allocate arrays
    ! radial and Z grids that Psi is defined on
    if (.not. allocated(R_grid)) then
		allocate (R_grid(NRBOX))
		allocate (Z_grid(NZBOX))
		allocate (PsiN_grid(NRBOX))
		allocate (Tflux_on_psiNgrid(NRBOX))
		allocate (rho_on_psiNgrid(NRBOX))
		allocate (PsiN_on_rhogrid(NRBOX))
    end if

    do i = 1, NRBOX
    	R_grid(i) = box_rmin + (box_rmax - box_rmin)*(i-1)/(NRBOX - 1)
    end do

    do i = 1, NZBOX
    	Z_grid(i) = box_zmin + (box_zmax - box_zmin)*(i-1)/(NZBOX - 1)
    end do

    do i = 1, NRBOX
    	PsiN_grid(i) = one*(i-1)/(NRBOX - 1)
    end do

! shift Psi to be zero on axis
    Psi = Psi - PSIAXIS
    PSIBOUND = PSIBOUND - PSIAXIS
    psiB = PSIBOUND

! Initialize spline coefficients for psi on 2D grid
	call psi_profile%cube_spline_2D_init(NRBOX, R_grid, NZBOX, Z_grid, Psi, Psi_name)

! Initialize spline coefficients for RBphi
    call T_profile%cube_spline_1D_init(NRBOX, R_grid, T, T_name)

! Initialize spline coefficients for Q on psiN grid
    call Q_profile%cube_spline_1D_init(NRBOX, PsiN_grid, Q, Q_name)


! Initialize spline coefficients for rho on psiN grid
  ! First have to calculate rho on psiN grid
	call calculate_rho_on_psiNgrid(PsiN_grid, Tflux_on_psiNgrid, rho_on_psiNgrid)
    call rho_profile%cube_spline_1D_init(NRBOX, PsiN_grid, rho_on_psiNgrid, rho_name)
    call Tflux_profile%cube_spline_1D_init(NRBOX, PsiN_grid, Tflux_on_psiNgrid, Tflux_name)

! Initialize spline coefficients for psiN on rho grid.
! Use the same rho grid as PsiN, [0, 1]
	call calculate_PsiN_on_rhogrid(PsiN_grid, psiN_on_rhogrid)
    call psiN_profile_rho%cube_spline_1D_init(NRBOX, PsiN_grid, psiN_on_rhogrid, psiN_name)

    return
  end subroutine initialize_eqdsk_magnetics_spline_interp

!****************************************************************************************

  subroutine  eqdsk_magnetics_spline_interp(rvec, bvec, gradbtensor, psi, gradpsi, psiN, &
            & gradpsiN, equib_err)

!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message
    use eqdsk_utilities_m, only : NRBOX, NZBOX, R_grid, Z_grid, PSIBOUND

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind), intent(out) :: Psi, gradpsi(3), psiN, gradpsiN(3)
    character(len=60), intent(out) :: equib_err

    real(KIND=rkind) :: PsiR, PsiZ, PsiRR, PsiRZ, PsiZZ, RBphi, RBphiR
    real(KIND=rkind) :: x, y, z, r
    real(KIND=rkind) :: br, bz, bphi, bp0
    real(KIND=rkind) :: dd_psi, dbrdr, dbrdz, dbzdr, dbzdz, dbphidr

    equib_err = ''
    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    r = sqrt(x**2+y**2)

	call Psi_profile%eval_2D_fpp(r, z, Psi, PsiR, PsiZ, PsiRR, PsiRZ, PsiZZ)

    call T_profile%eval_1D_fp(r, RBphi, RBphiR)

!   Magnetic field
    br = PsiZ/r
    bz = -PsiR/r
    bphi = RBphi/r
    gradpsi = (/-x*bz, -y*bz, r*br/)

!   Normalized Flux function x, y, z normalized to 1.0 at last closed flux surface
    psiN = psi/PSIBOUND
    gradpsiN = gradpsi/PSIBOUND

!   Magnetic field derivatives.

!   dbrdr = d(Br)/dr, dbrdz = d(Br)/dz.
    dbrdr = -br/r + PsiRZ/r
    dbrdz = PsiZZ/r

!   dbzdr = d(Bz)/dr, dbzdz = d(Bz)/dz.
    dbzdr = -bz/r - PsiRR/r
    dbzdz = - PsiRZ/r

!   dbphidr = d(Bphi)/dr.
    dbphidr = (RBphiR-bphi)/r

!   B field in the fixed (x,y,z) coordinates.
    bvec(1) = br*x/r - bphi*y/r
    bvec(2) = br*y/r + bphi*x/r
    bvec(3) = bz

!   d(Bx)/dx, d(Bx)/dy, and d(Bx)/dz:
    gradbtensor(1,1) = ( dbrdr*x**2 + br*y**2/r + (-dbphidr+bphi/r)*x*y ) / r**2
    gradbtensor(2,1) = ( (dbrdr-br/r)*x*y - dbphidr*y**2 - bphi*x**2/r ) / r**2
    gradbtensor(3,1) = dbrdz*x/r

!   d(By)/dx, d(By)/dy, and d(By)/dz:
    gradbtensor(1,2) = ( (dbrdr-br/r)*x*y + dbphidr*x**2 + bphi*y**2/r ) / r**2
    gradbtensor(2,2) = ( dbrdr*y**2 + br*x**2/r + (dbphidr-bphi/r)*x*y ) / r**2
    gradbtensor(3,2) = dbrdz*y/r

!   d(Bz)/dx, d(Bz)/dy, and d(Bz)/dz:
    gradbtensor(1,3) = dbzdr * x/r
    gradbtensor(2,3) = dbzdr * y/r
    gradbtensor(3,3) = dbzdz

    return
    end subroutine  eqdsk_magnetics_spline_interp

!********************************************************************

  subroutine  eqdsk_magnetics_spline_interp_psi(rvec, Psi, gradpsi, psiN, gradpsiN)
!   Simple  eqdsk equilibrium model originally based on notes from 7/28/1995 by Cai-ye Wang.
!   Reworked extensively by DBB.  See notes of 2-12-2022.

    use constants_m, only : rkind
    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message
    use eqdsk_utilities_m, only :   PSIBOUND
    use eqdsk_utilities_m, only : R_grid, Z_grid, NRBOX, NZBOX, PSIBOUND

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: Psi, gradpsi(3), psiN, gradpsiN(3)

    real(KIND=rkind) :: PsiR, PsiZ
    real(KIND=rkind) :: x, y, z, R
    real(KIND=rkind) :: br, bz, bphi, bp0

    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    R = sqrt(x**2+y**2)

	call Psi_profile%eval_2D_fp(r, z, Psi, PsiR, PsiZ)

!   Magnetic field
    br = PsiZ/r
    bz = -PsiR/r
    gradpsi = (/-x*bz, -y*bz, r*br/)

!   Normalized Flux function x, y, z normalized to 1.0 at last closed flux surface
    psiN = Psi/PSIBOUND
    gradpsiN = gradpsi/PSIBOUND

    return
  end subroutine  eqdsk_magnetics_spline_interp_psi

!********************************************************************
  subroutine  eqdsk_magnetics_spline_interp_rho(rvec, rho, gradrho)

    use constants_m, only : rkind
    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message
    use eqdsk_utilities_m, only :   PSIBOUND
    use eqdsk_utilities_m, only : R_grid, Z_grid, NRBOX, NZBOX, PSIBOUND

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: rho, gradrho(3)
    real(KIND=rkind) :: psi, gradpsi(3), psiN, gradpsiN(3)

    real(KIND=rkind) :: x, y, z, R
    real(KIND=rkind) :: drho_dpsi

    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    R = sqrt(x**2+y**2)

	call eqdsk_magnetics_spline_interp_psi(rvec, psi, gradpsi, psiN, gradpsiN)
	call eqdsk_magnetics_spline_interp_rho_psiN(psiN, rho, drho_dpsi)
	gradrho = gradpsiN*drho_dpsi

    return
  end subroutine  eqdsk_magnetics_spline_interp_rho

!********************************************************************
  subroutine  eqdsk_magnetics_spline_interp_Q_psiN(PsiN, Q, dQ_dPsi)

    implicit none

    real(KIND=rkind), intent(in) :: PsiN
    real(KIND=rkind), intent(out) :: Q, dQ_dPsi

	call Q_profile%eval_1D_fp(PsiN, Q, dQ_dPsi)

    return
  end subroutine  eqdsk_magnetics_spline_interp_Q_psiN
!********************************************************************

  subroutine  eqdsk_magnetics_spline_interp_rho_PsiN(PsiN, rho, drho_dPsi)

    implicit none

    real(KIND=rkind), intent(in) :: PsiN
    real(KIND=rkind), intent(out) :: rho, drho_dPsi

	call rho_profile%eval_1D_fp(PsiN, rho, drho_dPsi)

    return
  end subroutine  eqdsk_magnetics_spline_interp_rho_PsiN
!********************************************************************

  subroutine  eqdsk_magnetics_spline_interp_PsiN_rho(rho, PsiN, dPsiN_drho)

    implicit none

    real(KIND=rkind), intent(in) :: rho
    real(KIND=rkind), intent(out) :: PsiN, dPsiN_drho

	call psiN_profile_rho%eval_1D_fp(rho,  PsiN, dPsiN_drho)

    return
  end subroutine  eqdsk_magnetics_spline_interp_PsiN_rho
!********************************************************************

  subroutine  eqdsk_magnetics_spline_interp_Q_rho(rho, Q, dQ_drho)

    implicit none

    real(KIND=rkind), intent(in) :: rho
    real(KIND=rkind), intent(out) :: Q, dQ_drho
    real(KIND=rkind) ::  PsiN, dPsiN_dPsi, dPsiN_drho, dQ_dPsi

	call eqdsk_magnetics_spline_interp_PsiN_rho(rho, PsiN, dPsiN_drho)
	call Q_profile%eval_1D_fp(PsiN,  Q, dQ_dPsi)
	dQ_drho = dQ_dPsi*dPsiN_drho

    return
  end subroutine  eqdsk_magnetics_spline_interp_Q_rho
!********************************************************************

  subroutine calculate_rho_on_psiNgrid(PsiN_grid, Tflux_on_psiNgrid, rho_on_psiNgrid)
! Calculates toroidal flux, Tflux=Phi, and rho=sqrt(Phi) on uniform psiN grid by
! integrating d(Phi) = q(psi)d(psi)

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message
    USE trapezoid_quad_m, only : trapezoid_quad_cumulative

    implicit none

    real(KIND=rkind), intent(in) :: PsiN_grid(:)
    real(KIND=rkind), intent(out) :: rho_on_psiNgrid(size(PsiN_grid))
    real(KIND=rkind), intent(out) :: Tflux_on_psiNgrid(size(PsiN_grid))

    real(KIND=rkind) :: Q(size(PsiN_grid))
    real(KIND=rkind) :: dQ_dPsi, Tflux_total
    integer :: n_grid, i

    n_grid = size(PsiN_grid)

    do i = 1, n_grid
		call eqdsk_magnetics_spline_interp_Q_psiN(PsiN_grid(i), Q(i), dQ_dPsi)
	end do

	CALL trapezoid_quad_cumulative(PsiN_grid, Q, Tflux_total, Tflux_on_psiNgrid)
	rho_on_psiNgrid = Tflux_on_psiNgrid/Tflux_total
	rho_on_psiNgrid = sqrt(rho_on_psiNgrid)

    return
  end subroutine calculate_rho_on_psiNgrid

!********************************************************************

  subroutine calculate_psiN_on_rhogrid(rho_grid, psiN_on_rhogrid)
! Invert rho(psiN) using bisection to get psiN on a uniform rho grid.  Need this for
! such things as plotting profiles versus rho.  Will use it to generate a rho(psiN) spline
! function

    use constants_m, only : rkind, one, zero
    use diagnostics_m, only : message_unit, message
    use bisect_m, only : solve_bisection

    implicit none

    real(KIND=rkind), intent(in) :: rho_grid(:)
    real(KIND=rkind), intent(out) :: psiN_on_rhogrid(size(rho_grid))

    integer :: n_grid, i, ierr

! Tolerance for finding psiN of rho by bisection
    real(KIND=rkind), parameter :: bisection_eps =  one*10d-6

    n_grid = size(rho_grid)

	do i = 1, n_grid
	    call solve_bisection(f_rho_psiN, psiN_on_rhogrid(i), rho_grid(1), rho_grid(n_grid),&
	                         &rho_grid(i), bisection_eps, ierr)
	end do

    return
  end subroutine calculate_psiN_on_rhogrid

!*************************************************************************

 function f_rho_psiN(psiN)

    use constants_m, only : rkind, zero

	IMPLICIT NONE
    real(KIND=rkind) :: f_rho_psiN
    real(KIND=rkind) :: psiN, drho_dPsi

	call eqdsk_magnetics_spline_interp_rho_PsiN(PsiN, f_rho_psiN, drho_dPsi)

	return
 end function f_rho_psiN

!********************************************************************

    subroutine deallocate_eqdsk_magnetics_spline_interp_m
		! Nothing to deallocate
		return
    end subroutine deallocate_eqdsk_magnetics_spline_interp_m

 !********************************************************************

end module  eqdsk_magnetics_spline_interp_m
