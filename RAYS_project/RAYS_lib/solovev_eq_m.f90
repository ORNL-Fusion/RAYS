module solovev_eq_m
! A simple Solovev equilibrium model based on notes from Cai-Ye Wang of 1995 and DBB 2/12/2022.
!

    use constants_m, only : rkind

    implicit none

! data for magnetics
    real(KIND=rkind) :: rmaj, kappa, bphi0, iota0
    real(KIND=rkind) :: inner_bound, outer_bound, vert_bound, r_Zmax

    ! Flux function psi at plasma boundary
    real(KIND=rkind) :: psiB

! data for density and temperature
    character(len=60) :: dens_prof_model
    real(KIND=rkind) :: alphan1
    real(KIND=rkind) :: alphan2
    character(len=20), allocatable :: t_prof_model(:)
    real(KIND=rkind), allocatable :: alphat1(:)
    real(KIND=rkind), allocatable :: alphat2(:)

! data for bounding box
    real(KIND=rkind) :: box_rmin, box_rmax, box_zmin, box_zmax

 namelist /solovev_eq_list/ &
     & rmaj, outer_bound, kappa, bphi0, iota0, &
     & dens_prof_model, alphan1, alphan2, t_prof_model, alphat1, alphat2, &
     & box_rmin, box_rmax, box_zmin, box_zmax

!********************************************************************

contains

!********************************************************************

  subroutine initialize_solovev_eq_m(read_input)

    use species_m, only : nspec
    use diagnostics_m, only : message, text_message, message_unit, messages_to_stdout, verbosity

    implicit none
    logical, intent(in) :: read_input
 	integer :: input_unit, get_unit_number ! External, free unit finder

    real(KIND=rkind) :: bp0

    allocate( t_prof_model(0:nspec) )
    allocate( alphat1(0:nspec), alphat2(0:nspec) )

    call message(1)
    call text_message('Initializing solovev_eq_m ', 1)

    if (read_input .eqv. .true.) then
  		input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, solovev_eq_list)
        close(unit=input_unit)
    end if

! Write input namelist
    if (verbosity >= 0) then
		write(message_unit, solovev_eq_list)
		if (messages_to_stdout) write(*, solovev_eq_list)
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
    vert_bound = kappa/(2.*r_Zmax)*sqrt(outer_bound**4 + &
               & 2.*(r_Zmax**2 - outer_bound**2)*rmaj**2 - r_Zmax**4)

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
  end subroutine initialize_solovev_eq_m

!********************************************************************

  subroutine solovev_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)
!   Simple solovev equilibrium model originally based on notes from 7/28/1995 by Cai-ye Wang.
!   Reworked extensively by DBB.  See notes of 2-12-2022.
!
!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind), intent(out) :: ns(0:nspec), gradns(3,0:nspec)
    real(KIND=rkind), intent(out) :: ts(0:nspec), gradts(3,0:nspec)
    character(len=60), intent(out) :: equib_err


    real(KIND=rkind) :: x, y, z, r
    real(KIND=rkind) :: br, bz, bphi, bp0
    real(KIND=rkind) :: psi, gradpsi(3), psiN, gradpsiN(3)
    real(KIND=rkind) :: dd_psi, dbrdr, dbrdz, dbzdr, dbzdz, dbphidr
    integer :: is

    equib_err = ''
    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    r = sqrt(x**2+y**2)


! Check that we are in the box
    if (r < box_rmin .or. r > box_rmax) equib_err = 'R out_of_box'
    if (z < box_zmin .or. z > box_zmax) equib_err = 'z out_of_box'

!   Define
    bp0 = bphi0*iota0

! Get poloidal flux
    call solovev_psi(rvec, psi, gradpsi, psiN, gradpsiN)

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

!   Density profile.

    density: select case (trim(dens_prof_model))

        case ('constant')
          ns(:nspec) = n0s(:nspec)
          gradns = 0.

        case ('parabolic')
            ns = 0.
            gradns = 0.

            if (psiN < 1.0) then
                ns(0:nspec) = n0s(0:nspec) * (1.-psiN**(alphan2))**alphan1
                dd_psi = -alphan1*alphan2*psiN**(alphan2 - 1.)*(1.-psiN**(alphan2))&
                        & **(alphan1 - 1.)
                gradns(1, 0:nspec) = n0s(0:nspec)*dd_psi*gradpsiN(1)
                gradns(2, 0:nspec) = n0s(0:nspec)*dd_psi*gradpsiN(2)
                gradns(3, 0:nspec) = n0s(0:nspec)*dd_psi*gradpsiN(3)
            end if

        case default
            write(0,*) 'solovev_eq invalid dens_prof_model =', dens_prof_model
            stop 1

    end select density


!   Temperature profile.
    do is = 0, nspec
       temperature: select case (t_prof_model(is))


        case ('zero')
          ts(is) = 0.
          gradts(:,is) = 0.

        case ('constant')
          ns(:nspec) = n0s(:nspec)
          gradns = 0.

       case ('parabolic')
!      Parabolic: N.B. psi goes something like r**2 so if alphan2 = 1 and alphan1 = 2
!      the profile is pretty much parabolic

          ts = 0.
          gradts = 0.

          if (psiN < 1.) then
              ts(is) = t0s(is) * (1.-psiN**(alphat2(is)))**alphat1(is)
              dd_psi = -alphat1(is)*alphat2(is)*psiN**(alphat2(is) - 1.)* &
                      & (1.-psiN**(alphat2(is)))**alphat1(is)
              gradts(1,is) = t0s(is)*dd_psi*gradpsiN(1)
              gradts(2,is) = t0s(is)*dd_psi*gradpsiN(2)
              gradts(3,is) = t0s(is)*dd_psi*gradpsiN(3)
          end if

       case default
          write(0,*) 'SOLOVEV: Unknown t_prof_model: ', t_prof_model(0:nspec)
          stop 1

       end select temperature
    end do

! Do some checking for invalid values

    if (minval(ns) < 0.) equib_err = 'negative_dens'
    if (minval(ts) < 0.) equib_err = 'negative_temp'

    return
 end subroutine solovev_eq

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

!   Flux function x, y, z
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

!***********************************************************************
 subroutine write_solovev_profiles
    use species_m, only : nspec
    use diagnostics_m, only : message_unit

    implicit none

    real(KIND=rkind) :: dx, x
    real(KIND=rkind) :: rvec(3)
    real(KIND=rkind) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind) :: ns(0:nspec), gradns(3,0:nspec)
    real(KIND=rkind) :: ts(0:nspec), gradts(3,0:nspec)
    real(KIND=rkind) :: psi, gradpsi(3), psiN, gradpsiN(3)
    character(len=60) :: equib_err

    integer, parameter :: nx_points = 51
    integer :: ip, i

    character (len = *), parameter :: b3 = '   '
    character (len = *), parameter :: b4 = '    '
    character (len = *), parameter :: b8 = '        '
    character (len = *), parameter :: b9 = '         '
    character (len = *), parameter :: b10 = '          '
    character (len = *), parameter :: b12 = '            '

    dx = (outer_bound-inner_bound)/(nx_points-1)

    write (message_unit,*) '    x', b9,'ne', b12, 'bx', b9, 'by', b9, 'bz', b9, 'psi', &
            & b8, 'psi/psiB', b3,  'Te',b9, 'Ti(s)'

    do ip = 1, nx_points
        x = inner_bound + (ip-1)*dx
        rvec( : ) = (/ x, real(0.,KIND=rkind), real(0.,KIND=rkind) /)
        call solovev_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)
        call solovev_psi(rvec, psi, gradpsi, psiN, gradpsiN)
        write (message_unit,'(f11.5, a, e12.5, 3f11.5, f11.5, f11.5,  7f11.5)') &
               & x,'  ', ns(0), bvec, psi, psi/psiB, (ts(i), i=0, nspec)
!        write(*,*) 'x = ', x, '  gradpsi = ', gradpsi
!        write(*,*) 'gradbtensor = ', gradbtensor
    end do

 end subroutine write_solovev_profiles

!********************************************************************

    subroutine deallocate_solovev_eq_m
		if (allocated(t_prof_model)) then
			deallocate( t_prof_model )
			deallocate( alphat1 )
		end if
		return
    end subroutine deallocate_solovev_eq_m

end module solovev_eq_m
