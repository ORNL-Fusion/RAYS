module  eqdsk_magnetics_lin_interp_m
! A simple  eqdsk equilibrium model. This code uses eqdsk routines in module eqdsk_utilities_m.f90
! These were adapted from similar codes by Richard Fitzpatrick.  It uses Fitzpatrick's functions
! for interpolating psi and its derivatives using simple 2 point or 3 point approximations.
! In that regard they are not expected to be very accurate.  This version will soon be supplanted
! by one that uses cubic splines for the interpolation.
!
! N.B. Values of Psi are shifted on initialization so that Psi is zero on axis.

    use constants_m, only : rkind
    
    implicit none

    character (len = 100) :: eqdsk_file_name

! data for magnetics
    real(KIND=rkind) :: rmaj, kappa, bphi0, iota0
    real(KIND=rkind) :: inner_bound, outer_bound, upper_bound, lower_bound

    ! Flux function psi at plasma boundary
    real(KIND=rkind) :: psiB

  namelist / eqdsk_magnetics_lin_interp_list/ eqdsk_file_name
     
!********************************************************************

contains

!********************************************************************

  subroutine initialize_eqdsk_magnetics_lin_interp(read_input, r_axis, z_axis, &
               & box_rmin, box_rmax, box_zmin, box_zmax, &
               & inner_bound, outer_bound, upper_bound, lower_bound)

    use constants_m, only : input_unit
    use species_m, only : nspec
    use diagnostics_m, only : message, message_unit, verbosity
 !   use eqdsk_utilities_m, only : ReadgFile, R_grid, Z_grid, dR, dZ

    use eqdsk_utilities_m, only : ReadgFile, WritegFile, R_grid, Z_grid, dR, dZ, &
        & string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, &
        & PSIAXIS, PSIBOUND, B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, &
        & ZBOUND, RLIM, ZLIM, &
        & GetPsi, GetRBphi, GetPsiR, GetPsiZ, GetPsiRR, GetPsiZZ, GetPsiRZ, GetRBphiR

    implicit none
    
    logical, intent(in) :: read_input

! Geometry data
    ! Magnetic axis
    real(KIND=rkind), intent(out) :: r_axis, z_axis
    ! data for bounding box of computational domain
    real(KIND=rkind), intent(out) :: box_rmin, box_rmax, box_zmin, box_zmax
    ! data for plasma boundary
    real(KIND=rkind), intent(out) :: inner_bound, outer_bound, upper_bound, lower_bound

    integer :: i
    
    write(*,*) 'initialize_eqdsk_magnetics_lin_interp'   

    if (read_input .eqv. .true.) then 
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, eqdsk_magnetics_lin_interp_list)
        close(unit=input_unit)
        write(message_unit, eqdsk_magnetics_lin_interp_list)
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

    call message('Inner boundary = ', inner_bound)
    call message('Outer boundary = ', outer_bound)
    call message('Upper boundary = ', upper_bound)
    call message('Lower boundary = ', lower_bound)

    write(*,*) 'Inner boundary = ', inner_bound
    write(*,*) 'Outer boundary = ', outer_bound
    write(*,*) 'Lower boundary = ', lower_bound
    write(*,*) 'Upper boundary = ', upper_bound
    
    ! radial and Z grids that Psi is defined on
    allocate (R_grid(NRBOX))
    allocate (Z_grid(NZBOX))

    do i = 1, NRBOX
    	R_grid(i) = box_rmin + (box_rmax - box_rmin)*(i-1)/(NRBOX - 1)
    end do

    do i = 1, NZBOX
    	Z_grid(i) = box_zmin + (box_zmax - box_zmin)*(i-1)/(NZBOX - 1)
    end do

    dR = (R_grid(2) - R_grid(1))/2.
    dZ = (Z_grid(2) - Z_grid(1))/2.

! shift Psi to be zero on axis
    Psi = Psi - PSIAXIS
    PSIBOUND = PSIBOUND - PSIAXIS
    psiB = PSIBOUND
    
    return    
  end subroutine initialize_eqdsk_magnetics_lin_interp


  subroutine  eqdsk_magnetics_lin_interp(rvec, bvec, gradbtensor, psi, gradpsi, psiN, gradpsiN, equib_err)

!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message
    use eqdsk_utilities_m, only : GetPsi, GetRBphi, GetPsiR, GetPsiZ, GetPsiRR, GetPsiZZ,&
        & GetPsiRZ, GetRBphiR, PSIBOUND

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3) 
    real(KIND=rkind), intent(out) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind), intent(out) :: psi, gradpsi(3), psiN, gradpsiN(3)
    character(len=20), intent(out) :: equib_err


    real(KIND=rkind) :: x, y, z, r
    real(KIND=rkind) :: br, bz, bphi, bp0
    real(KIND=rkind) :: dd_psi, dbrdr, dbrdz, dbzdr, dbzdz, dbphidr
    integer :: is

    equib_err = ''
    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    r = sqrt(x**2+y**2)

! Check that we are in the plasma. Set equib_err but don't stop
    if (psiN > 1.) equib_err = 'psi >1 out_of_plasma'

!   Poloidfal flux function 
    psi = GetPsi(R, z)

!   Magnetic field
    br = -GetPsiZ(R,z)/r
    bz = GetPsiR(R,z)/r
    bphi = GetRBphi(r)/r
    gradpsi = (/x*bz, y*bz, -R*br/)

!   Normalized Flux function x, y, z normalized to 1.0 at last surface
    psiN = psi/PSIBOUND
    gradpsiN = gradpsi/PSIBOUND

!   Magnetic field derivatives.

!   dbrdr = d(Br)/dr, dbrdz = d(Br)/dz.
    dbrdr = -br/r - GetPsiRZ(r,z)/r
    dbrdz = -GetPsiZZ(r, z)/r
    
!   dbzdr = d(Bz)/dr, dbzdz = d(Bz)/dz.
    dbzdr = -bz/r + GetPsiRR(r, z)/r
    dbzdz = GetPsiRZ(r,z)/r
    
!   dbphidr = d(Bphi)/dr.
    dbphidr = (GetRBphiR(r)-bphi)/r

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
    end subroutine  eqdsk_magnetics_lin_interp

!********************************************************************

  subroutine  eqdsk_magnetics_lin_interp_psi(rvec, psi, gradpsi, psiN, gradpsiN)
!   Simple  eqdsk equilibrium model originally based on notes from 7/28/1995 by Cai-ye Wang.
!   Reworked extensively by DBB.  See notes of 2-12-2022. 
!
!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

    use constants_m, only : rkind
    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message
    use eqdsk_utilities_m, only :   GetPsi, GetPsiR, GetPsiZ, PSIBOUND
   
    implicit none

    real(KIND=rkind), intent(in) :: rvec(3) 
    real(KIND=rkind), intent(out) :: psi, gradpsi(3), psiN, gradpsiN(3)


    real(KIND=rkind) :: x, y, z, R
    real(KIND=rkind) :: br, bz, bphi, bp0

    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    R = sqrt(x**2+y**2)
 
!   Poloidfal flux function 
    psi = GetPsi(R, z)


!   Magnetic field
    br = -GetPsiZ(R,z)/r
    bz = GetPsiR(R,z)/r
    gradpsi = (/x*bz, y*bz, -R*br/)

!   Normalized Flux function x, y, z normalized to 1.0 at last surface
    psiN = psi/PSIBOUND
    gradpsiN = gradpsi/PSIBOUND
    
    return 
  end subroutine  eqdsk_magnetics_lin_interp_psi
 
end module  eqdsk_magnetics_lin_interp_m
