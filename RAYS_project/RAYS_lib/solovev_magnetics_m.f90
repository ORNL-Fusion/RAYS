module solovev_magnetics_m
! A simple Solovev equilibrium model based on notes from Cai-Ye Wang of 1995 and DBB 2/12/2022.
!

    use constants_m, only : rkind

    implicit none

! data for magnetics
    real(KIND=rkind) :: rmaj, kappa, bphi0, iota0
    real(KIND=rkind) :: inner_bound, outer_bound, vert_bound, r_Zmax, zmax_sq
    real(KIND=rkind) :: outer_boundary

    ! Flux function psi at plasma boundary
    real(KIND=rkind) :: psiB

! data for bounding box
    real(KIND=rkind) :: box_rmin, box_rmax, box_zmin, box_zmax

 namelist /solovev_magnetics_list/ &
     & rmaj, outer_boundary, kappa, bphi0, iota0, &
     & box_rmin, box_rmax, box_zmin, box_zmax

!********************************************************************

contains

!********************************************************************

  subroutine initialize_solovev_magnetics(read_input, r_axis, z_axis, &
	   & arg_box_rmin, arg_box_rmax, arg_box_zmin, arg_box_zmax, &
	   & inner_bound, outer_bound, upper_bound, lower_bound)

    use species_m, only : nspec
    use diagnostics_m, only : message, text_message, message_unit, messages_to_stdout, verbosity

    implicit none
    logical, intent(in) :: read_input
 	integer :: input_unit, get_unit_number ! External, free unit finder

! Geometry data
    ! Magnetic axis
    real(KIND=rkind), intent(out) :: r_axis, z_axis
    ! data for bounding box of computational domain
    real(KIND=rkind), intent(out) :: arg_box_rmin, arg_box_rmax, arg_box_zmin,arg_box_zmax
! data for plasma boundary
    real(KIND=rkind), intent(out) :: inner_bound, outer_bound, upper_bound, lower_bound

     real(KIND=rkind) :: bp0

    call message(1)
    call text_message('Initializing solovev_magnetics_m ', 1)

    if (read_input .eqv. .true.) then
  		input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, solovev_magnetics_list)
        close(unit=input_unit)
    end if

! Write input namelist
    if (verbosity >= 0) then
		write(message_unit, solovev_magnetics_list)
		if (messages_to_stdout) write(*, solovev_magnetics_list)
    end if

	arg_box_rmin = box_rmin
	arg_box_rmax = box_rmax
	arg_box_zmin = box_zmin
	arg_box_zmax = box_zmax


    outer_bound = outer_boundary

! Calculate inner and boundary
     ! Check that inner boundary is real number
     if ( outer_bound < rmaj .or. outer_bound >= Sqrt(2.)*rmaj ) then
        call message('Inner boundary complex, outer_bound >=  sqrt2*rmaj = ', outer_bound)
        write(*,*) 'Inner boundary complex, outer_bound >=  sqrt2*rmaj = ', outer_bound
        stop
    end if

!   Define
    bp0 = bphi0*iota0

!   Flux at plasma boundary
    psiB = .5*bp0 * (outer_bound**2-rmaj**2)**2/rmaj**2/4.
    inner_bound = sqrt(2.*rmaj**2 -outer_bound**2 )

    ! radius of maximum in z
    r_Zmax = (2.*outer_bound**2 * rmaj**2 - outer_bound**4)**0.25

    ! z at r_Zmax
    vert_bound = kappa/(2.*r_Zmax)*sqrt(outer_bound**4 + &
               & 2.*(r_Zmax**2 - outer_bound**2)*rmaj**2 - r_Zmax**4)

    call message('PsiB = ', psiB, 1)
    call message('Inner boundary = ', inner_bound, 1)
    call message('Outer boundary = ', outer_bound, 1)
    call message('Vertical boundary = ', vert_bound, 1)
    call message('radius of z_max = ', r_Zmax, 1)

    r_axis = rmaj
    z_axis = 0.
    upper_bound = vert_bound
    lower_bound = -vert_bound

    return
  end subroutine initialize_solovev_magnetics


  subroutine solovev_magnetics(rvec, bvec, gradbtensor, psi, gradpsi, psiN, gradpsiN, equib_err)
!   Simple solovev equilibrium model originally based on notes from 7/28/1995 .
!   by Cai-ye Wang.  Reworked extensively by DBB.  See notes of 2-12-2022.
!
!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message, text_message

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind), intent(out) :: psi, gradpsi(3), psiN, gradpsiN(3)
    character(len=60), intent(out) :: equib_err


    real(KIND=rkind) :: x, y, z, r
    real(KIND=rkind) :: br, bz, bphi, bp0
    real(KIND=rkind) :: dd_psi, dbrdr, dbrdz, dbzdr, dbzdz, dbphidr
    integer :: is

    equib_err = ''
    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    r = sqrt(x**2+y**2)


! Check that we are in the box
    if (r < box_rmin .or. r > box_rmax) equib_err = 'R out_of_bounds'
    if (z < box_zmin .or. z > box_zmax) equib_err = 'z out_of_bounds'

    if (equib_err /= '') then
        call text_message('solovev_magnetics:  equib_err ', equib_err, 1)
        return
    end if

!   Define
    bp0 = bphi0*iota0

! Get poloidal flux
    call solovev_magnetics_psi(rvec, psi, gradpsi, psiN, gradpsiN)

!   Magnetic field and its derivatives.
    br = -bp0*r*z/(rmaj*kappa)**2
    bz = bp0 * ( (z/(rmaj*kappa))**2 + .5*((r/rmaj)**2-1.) )
    bphi = bphi0*rmaj/r


!   dbrdr = d(Br)/dr, dbrdz = d(Br)/dz.
    dbrdr = br/r
    dbrdz = -bp0*r/(rmaj*kappa)**2

!   dbzdr = d(Bz)/dr, dbzdz = d(Bz)/dz.
    dbzdr = bp0*r/rmaj**2
    dbzdz = bp0*2.*z/(rmaj*kappa)**2

!   dbphidr = d(Bphi)/dr.
    dbphidr = -bphi/r

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
    end subroutine solovev_magnetics

!********************************************************************

  subroutine solovev_magnetics_psi(rvec, psi, gradpsi, psiN, gradpsiN)
!   Simple solovev equilibrium model originally based on notes from 7/28/1995 by Cai-ye Wang.
!   Reworked extensively by DBB.  See notes of 2-12-2022.
!
!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

    use constants_m, only : rkind
    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: psi, gradpsi(3), psiN, gradpsiN(3)


    real(KIND=rkind) :: x, y, z, R
    real(KIND=rkind) :: br, bz, bphi, bp0

    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    R = sqrt(x**2+y**2)

!   Define
    bp0 = bphi0*iota0

!   Flux function x, y, z normalized to one at last surface (z=0, r=outer_bound)
    psi = .5*bp0 * ( (R*z/(rmaj*kappa))**2 + ((R**2-rmaj**2)**2)/rmaj**2/4. )


!   Magnetic field
    br = -bp0*R*z/(rmaj*kappa)**2
    bz = bp0 * ( (z/(rmaj*kappa))**2 + .5*((R/rmaj)**2-1.) )
    gradpsi = (/x*bz, y*bz, -R*br/)

!   Normalized Flux function x, y, z normalized to 1.0 at last surface (z=0, R=outer_bound)
    psiN = psi/psiB
    gradpsiN = gradpsi/psiB

    return
  end subroutine solovev_magnetics_psi

!********************************************************************

    subroutine deallocate_solovev_magnetics_m
        return ! nothing to deallocate
    end subroutine deallocate_solovev_magnetics_m

end module solovev_magnetics_m
