module  eqdsk_magnetics_spline_interp_m
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

! Data for splines
    ! bcspline arguments
    real(KIND=rkind), allocatable :: fspl(:,:,:,:), wk(:)
    real(KIND=rkind), allocatable :: bcxmin(:),bcxmax(:)      ! (inth) if used
    real(KIND=rkind), allocatable :: bcthmin(:),bcthmax(:)    ! (inx) if used

    ! bcspeval arguments, if not already declared above
    real(KIND=rkind) :: fval(6)
    integer :: iselect(6)
   
! Stuff for 1D splines
! cspeval arguments

    real(KIND=rkind), allocatable ::  fspl_1D(:,:)
    real(KIND=rkind), allocatable ::  wk_1D(:)
    real(KIND=rkind) ::  fval_1D(3)
    integer :: iselect_1D(3)
    
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
    use diagnostics_m, only : message, message_unit, verbosity
 
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

! bcspline arguments
     integer ibcxmin,ibcxmax,ibcthmin,ibcthmax,ilinx, ilinth,ier

    integer :: i, j, nwk
    
    write(*,*) 'initialize_eqdsk_magnetics_spline_interp'   

    if (read_input .eqv. .true.) then 
  		input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, eqdsk_magnetics_spline_interp_list)
        close(unit=input_unit)
        write(message_unit, eqdsk_magnetics_spline_interp_list)
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
    
! Allocate arrays
    ! radial and Z grids that Psi is defined on
    if (.not. allocated(R_grid)) then
		allocate (R_grid(NRBOX))
		allocate (Z_grid(NZBOX))
		allocate (fspl(4,4,NRBOX,NZBOX))
		allocate (fspl_1D(4,NRBOX))
	    nwk = 4*NRBOX*NZBOX +5*max(NRBOX,NZBOX)
		allocate (wk(nwk))
		allocate (wk_1D(NRBOX))
		allocate (bcxmin(1),bcxmax(1),bcthmin(1),bcthmax(1)) ! not used but must be allocated	
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

 write (*,*) 'PSIAXIS = ', PSIAXIS, '   maxval(Psi) = ', maxval(Psi)
! Set f(1,1,i,j) = eqdsk Psi
    do i = 1, NRBOX
    do j = 1, NZBOX
    	fspl(1,1,i,j) = Psi(i,j)
    end do
    end do

! Set f(1,1,i) = eqdsk T
    do i = 1, NRBOX
    	fspl_1D(1,i) = T(i)
    end do

! Set up spline coefficients for Psi
    ibcxmin = 0
    bcxmin = 0.
    ibcxmax = 0
    bcxmax = 0.
    ibcthmin = 0
    bcthmin = 0.
    ibcthmax = 0
    bcthmax = 0.
    ilinx = 1
    ilinth = 1
    
    call bcspline(R_grid,NRBOX,Z_grid,NZBOX,fspl,NRBOX, &
         ibcxmin,bcxmin,ibcxmax,bcxmax, &
         ibcthmin,bcthmin,ibcthmax,bcthmax, &
         wk,nwk,ilinx,ilinth,ier)
         if (ier .ne. 0) write (*,*) 'bcspline: ier = ', ier

! Set up spline coefficients for RBphi
    
    call cspline(R_grid,NRBOX,fspl_1D,ibcxmin,bcxmin,ibcxmax,bcxmax,wk_1D,NRBOX,ilinx,ier)         
         if (ier .ne. 0) write (*,*) 'bcspline: ier = ', ier
    
    return    
  end subroutine initialize_eqdsk_magnetics_spline_interp


  subroutine  eqdsk_magnetics_spline_interp(rvec, bvec, gradbtensor, psi, gradpsi, psiN, gradpsiN, equib_err)

!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message
    use eqdsk_utilities_m, only : NRBOX, NZBOX, R_grid, Z_grid, PSIBOUND

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3) 
    real(KIND=rkind), intent(out) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind), intent(out) :: Psi, gradpsi(3), psiN, gradpsiN(3)
    character(len=20), intent(out) :: equib_err

! Spline variables
    integer :: iselect(6) = 1 ! Selector for pspline outputs. Output everthing.
    integer :: iselect_1D(3) = 1 ! Selector for pspline outputs. Output everthing.
    integer :: ilinx = 0
    integer :: ilinth = 0
    integer :: ier

    real(KIND=rkind) :: PsiR, PsiZ, PsiRR, PsiRZ, PsiZZ, RBphi, RBphiR
    real(KIND=rkind) :: x, y, z, r
    real(KIND=rkind) :: br, bz, bphi, bp0
    real(KIND=rkind) :: dd_psi, dbrdr, dbrdz, dbzdr, dbzdz, dbphidr
 
    equib_err = ''
    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    r = sqrt(x**2+y**2)

! evaluate spline fits
	call bcspeval(r,z,iselect,fval,R_grid,NRBOX,Z_grid,NZBOX,ilinx,ilinth,fspl,NRBOX,ier)            
    if (ier .ne. 0) write (*,*) 'bcspeval: ier = ', ier

    Psi = fval(1)
    PsiR = fval(2)
    PsiZ = fval(3)
    PsiRR = fval(4)
    PsiRZ = fval(5)
    PsiZZ = fval(6)

    call cspeval(r,iselect_1D,fval_1D,R_grid,NRBOX,ilinx,fspl_1D,ier)
    if (ier .ne. 0) write (*,*) 'cspeval: ier = ', ier
    
    RBphi = fval_1D(1)
    RBphiR = fval_1D(2)

!   Magnetic field
    br = -PsiZ/r
    bz = PsiR/r
    bphi = RBphi/r
    gradpsi = (/x*bz, y*bz, -R*br/)
 
!   Normalized Flux function x, y, z normalized to 1.0 at last surface
    psiN = psi/PSIBOUND
    gradpsiN = gradpsi/PSIBOUND
    ! Check that we are in the plasma. Set equib_err but don't stop
    if (psiN > 1.) equib_err = 'psi >1 out_of_plasma'

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

    integer, parameter, dimension(6) :: psi_select = (/ 1, 1, 1, 0, 0, 0 /)
    real (KIND=rkind) :: fval(6)
    real(KIND=rkind) :: PsiR, PsiZ
    real(KIND=rkind) :: x, y, z, R
    real(KIND=rkind) :: br, bz, bphi, bp0
    integer :: ier
    integer :: ilinx = 0
    integer :: ilinth = 0

    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    R = sqrt(x**2+y**2)
    fval = 0.

! evaluate spline fit
	call bcspeval(R,z,psi_select,fval,R_grid,NRBOX,Z_grid,NZBOX,ilinx,ilinth,fspl,NRBOX,ier)            
    if (ier .ne. 0) write (*,*) 'bcspeval: ier = ', ier

    Psi = fval(1) ! Poloidfal flux function 
    PsiR = fval(2)
    PsiZ = fval(3)


!   Magnetic field
    br = -PsiZ/r
    bz = PsiR/r
    gradpsi = (/x*bz, y*bz, -R*br/)

!   Normalized Flux function x, y, z normalized to 1.0 at last surface
    psiN = Psi/PSIBOUND
    gradpsiN = gradpsi/PSIBOUND

    
    return 
  end subroutine  eqdsk_magnetics_spline_interp_psi
!********************************************************************

    subroutine deallocate_eqdsk_magnetics_spline_interp_m
		if (allocated(fspl)) then
			deallocate( fspl ,wk, bcxmin, bcxmax, bcthmin, bcthmax )
			deallocate( fspl_1D, wk_1D)
		end if
		return
    end subroutine deallocate_eqdsk_magnetics_spline_interp_m

 !********************************************************************

end module  eqdsk_magnetics_spline_interp_m
