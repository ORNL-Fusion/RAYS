module solovev_magnetics_m
! A simple Solovev equilibrium model based on notes from Cai-Ye Wang of 1995 and DBB 2/12/2022.
!

    use constants_m, only : rkind
    
    implicit none

! data for magnetics
    real(KIND=rkind) :: rmaj, kappa, bphi0, iota0
    real(KIND=rkind) :: inner_bound, outer_bound, vert_bound, r_Zmax

    ! Flux function psi at plasma boundary
    real(KIND=rkind) :: psiB
    
! data for bounding box
    real(KIND=rkind) :: box_rmin, box_rmax, box_zmin, box_zmax

 namelist /solovev_magnetics_list/ &
     & rmaj, outer_bound, kappa, bphi0, iota0, &
     & box_rmin, box_rmax, box_zmin, box_zmax
     
!********************************************************************

contains

!********************************************************************

  subroutine initialize_solovev_magnetics(read_input)

    use constants_m, only : input_unit
    use species_m, only : nspec
    use diagnostics_m, only : message, message_unit, verbosity
    
    implicit none
    logical, intent(in) :: read_input

    real(KIND=rkind) :: bp0

    if (read_input .eqv. .true.) then    
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, solovev_magnetics_list)
        close(unit=input_unit)
        write(message_unit, solovev_magnetics_list)
    end if
    
! Calculate inner and boundary
     ! Check that inner boundary is real number
     if ( outer_bound < rmaj .or. outer_bound >= Sqrt(2.)*rmaj ) then
        call message('Inner boundary complex, outer_bound >=  sqrt2) rmaj = ', outer_bound)
        write(*,*) 'Inner boundary complex, outer_bound >=  sqrt2) rmaj = ', outer_bound
    end if

!   Define
    bp0 = bphi0*iota0
  
!   Flux at plasma boundary
    psiB = .5*bp0 * (outer_bound**2-rmaj**2)**2/rmaj**2/4.
   
    inner_bound = sqrt(2.*rmaj**2 -outer_bound**2 )
    
    ! radius of maximum in z
    r_Zmax = (2.*outer_bound**2 * rmaj**2 - outer_bound**4)**0.25
    ! z at r_Zmax
    vert_bound = sqrt(kappa)/(2.*r_Zmax)*sqrt(outer_bound**4 - r_Zmax**4 +2.*r_Zmax**2*rmaj**2- &
              &  2.*outer_bound**2*rmaj**2)

    call message('PsiB = ', psiB)
    call message('Inner boundary = ', inner_bound)
    call message('Outer boundary = ', outer_bound)
    call message('Vertical boundary = ', vert_bound)
    call message('radius of z_max = ', r_Zmax)

    write(*,*) 'PsiB = ', psiB
    write(*,*) 'Inner boundary = ', inner_bound
    write(*,*) 'Outer boundary = ', outer_bound
    write(*,*) 'Vertical boundary = ', vert_bound

    if (verbosity > 2)  then
        call write_solovev_profiles
    end if
    
    return
  end subroutine initialize_solovev_magnetics


  subroutine solovev_magnetics(rvec, bvec, gradbtensor, psi, gradpsi, psiN, gradpsiN, equib_err)
!   Simple solovev equilibrium model originally based on notes from 7/28/1995 .
!   by Cai-ye Wang.  Reworked extensively by DBB.  See notes of 2-12-2022. 
!
!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message
    
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
 

! Check that we are in the box
    if (r < box_rmin .or. r > box_rmax) equib_err = 'R out_of_bounds'
    if (z < box_zmin .or. z > box_zmax) equib_err = 'z out_of_bounds'

!   Define
    bp0 = bphi0*iota0

! Get poloidal flux
    call solovev_psi(rvec, psi, gradpsi, psiN, gradpsiN)     

! Check that we are in the plasma
    if (psiN > 1.) equib_err = 'psi >1 out_of_plasma'
    
    if (equib_err /= '') then
        write (message_unit, *) 'solovev_eq:  equib_err ', equib_err
        return
    end if

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
  
  end module solovev_magnetics_m

!********************************************************************

  subroutine solovev_psi(rvec, psi, gradpsi, psiN, gradpsiN)
!   Simple solovev equilibrium model originally based on notes from 7/28/1995 by Cai-ye Wang.
!   Reworked extensively by DBB.  See notes of 2-12-2022. 
!
!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

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
  
  end subroutine solovev_psi
 
end module solovev_magnetics_m
