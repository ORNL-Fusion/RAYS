module  eqdsk_magnetics_spline_interp_m
! A simple  eqdsk equilibrium model. This code uses eqdsk routines in module eqdsk_utilities_m.f90
! These were adapted from similar codes by Richard Fitzpatrick.
!
! N.B. Values of Psi are shifted on initialization so that Psi is zero on axis.

! Working notes
! DBB (8/2/2024) Changed error return 'if psi > 1' to 'if psi > plasma_psi_limit' to
! allow rays outside last closed flux surface. plasma_psi_limit defaults to 1.0 but can be
! reset in namelist

    use constants_m, only : rkind
    use quick_cube_splines_m, only : cube_spline_function_1D, cube_spline_function_2D

    implicit none

    character (len = 100) :: eqdsk_file_name

! Stuff for 2D spline profiles i.e. psi

    type(cube_spline_function_2D) :: Psi_profile
 	character (len = 80) :: Psi_name = 'Psi_profile'

! Stuff for 1D spline profiles

    type(cube_spline_function_1D) :: T_profile
 	character (len = 80) :: T_name = 'T_profile'

    ! Flux function psi at plasma boundary
    real(KIND=rkind) :: psiB

  namelist / eqdsk_magnetics_spline_interp_list/ eqdsk_file_name

!********************************************************************

contains

!********************************************************************

  subroutine initialize_eqdsk_magnetics_spline_interp(read_input, r_axis, z_axis, &
               & box_rmin, box_rmax, box_zmin, box_zmax, &
               & inner_bound, outer_bound, upper_bound, lower_bound)

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
		call message('Inner boundary = ', inner_bound)
		call message('Outer boundary = ', outer_bound)
		call message('Upper boundary = ', upper_bound)
		call message('Lower boundary = ', lower_bound)

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
    end if

    do i = 1, NRBOX
    	R_grid(i) = box_rmin + (box_rmax - box_rmin)*(i-1)/(NRBOX - 1)
    end do

    do i = 1, NZBOX
    	Z_grid(i) = box_zmin + (box_zmax - box_zmin)*(i-1)/(NZBOX - 1)
    end do

! shift Psi to be zero on axis
    Psi = Psi - PSIAXIS
    PSIBOUND = PSIBOUND - PSIAXIS
    psiB = PSIBOUND

! Initialize spline coefficients for psi

	call psi_profile%cube_spline_2D_init(NRBOX, R_grid, NZBOX, Z_grid, Psi, Psi_name)

! Initialize spline coefficients for RBphi

    call T_profile%cube_spline_1D_init(NRBOX, R_grid, T, T_name)

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
    br = -PsiZ/r
    bz = PsiR/r
    bphi = RBphi/r
    gradpsi = (/x*bz, y*bz, -R*br/)

!   Normalized Flux function x, y, z normalized to 1.0 at last closed flux surface
    psiN = psi/PSIBOUND
    gradpsiN = gradpsi/PSIBOUND

!   Magnetic field derivatives.

!   dbrdr = d(Br)/dr, dbrdz = d(Br)/dz.
    dbrdr = -br/r - PsiRZ/r
    dbrdz = -PsiZZ/r

!   dbzdr = d(Bz)/dr, dbzdz = d(Bz)/dz.
    dbzdr = -bz/r + PsiRR/r
    dbzdz = PsiRZ/r

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
!
!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

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
    br = -PsiZ/r
    bz = PsiR/r
    gradpsi = (/x*bz, y*bz, -R*br/)

!   Normalized Flux function x, y, z normalized to 1.0 at last closed flux surface
    psiN = Psi/PSIBOUND
    gradpsiN = gradpsi/PSIBOUND


    return
  end subroutine  eqdsk_magnetics_spline_interp_psi
!********************************************************************

    subroutine deallocate_eqdsk_magnetics_spline_interp_m
		! Nothing to deallocate
		return
    end subroutine deallocate_eqdsk_magnetics_spline_interp_m

 !********************************************************************

end module  eqdsk_magnetics_spline_interp_m
