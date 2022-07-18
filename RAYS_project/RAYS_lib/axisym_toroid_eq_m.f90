module axisym_toroid_eq_m
! A general toroidal equilibrium module that supports a choice of magnetics, density,
! and temperature models.  It is assumed that there exists a poloidal flux function 
! on which the densities and temperatures are constant. An external magnetics module is
! required for each magnetics model which provides B, grad B tensor, and the flux function,
! psi, along with it's gradient (poloidal_flux_model).
!
! Input to axisym_toroid_eq() is rvec = (x,y,z) position
! Outputs from axisym_toroid_eq() are all quantities needed to fill an eq_point type

    use constants_m, only : rkind
    
    implicit none

! data for magnetics
    character(len=15) :: magnetics_model

<<<<<<< HEAD
! data for density and temperature
    character(len=15) :: density_prof_model
    real(KIND=rkind) :: alphan1
    real(KIND=rkind) :: alphan2
    character(len=15), allocatable :: temperature_prof_model(:)
    real(KIND=rkind), allocatable :: alphat1(:)
    real(KIND=rkind), allocatable :: alphat2(:)

 namelist /axisym_toroid_eq_list/&
     & magnetics_model, &
     & density_prof_model, &
     & alphan1, alphan2, & ! parameters for parabolic model
     & temperature_prof_model, &
     & alphat1, alphat2 ! parameters for parabolic model
=======
! data for slab density and temperature
    character(len=15) :: dens_prof_model
    character(len=15), allocatable :: t_prof_model(:)

 namelist /axisym_toroid_eq_list/&
     & poloidal_flux_model, &
     & density_prof_model, &
     & alphan1, alphan2, ! parameters for parabolic model
     & temperature_prof_model, &
     & alphat1, alphat2, & ! parameters for parabolic model
>>>>>>> 47dd208b739c6b82af6ec9c77d62cf52c1c3da99
     
!********************************************************************

contains

!********************************************************************

  subroutine initialize_axisym_toroid_eq(read_input)

    use constants_m, only : input_unit
    use species_m, only : nspec
    use diagnostics_m, only : message, message_unit, verbosity

    use solovev_magnetics_m, only : initialize_solovev_magnetics, solovev_magnetics

    implicit none
    logical, intent(in) :: read_input

    real(KIND=rkind) :: bp0
    
<<<<<<< HEAD
    allocate( temperature_prof_model(0:nspec) )
=======
    allocate( t_prof_model(0:nspec) )
>>>>>>> 47dd208b739c6b82af6ec9c77d62cf52c1c3da99
    allocate( alphat1(0:nspec), alphat2(0:nspec) )

    if (read_input .eqv. .true.) then    
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, toroid_eq_list)
        close(unit=input_unit)
        write(message_unit, toroid_eq_list)
    end if
    
<<<<<<< HEAD
    magnetics: select case (trim(magnetics_model))
=======
    magnetics: select case (trim(poloidal_flux_model))
>>>>>>> 47dd208b739c6b82af6ec9c77d62cf52c1c3da99
       case ('solovev_magnetics')
          call initialize_solovev_magnetics(read_input)

       case default
          write(0,*) 'initialize_axisym_toroid_eq: unknown magnetics model =', magnetics_model
          call text_message('initialize_axisym_toroid_eq: unknown magnetics model',&
          & trim(magnetics_model),0)
          stop 1
    end select magnetics

! Densities and temperatures.  No initialization needed beyond reading the namelist at
! this time.  If I add something to read profiles from a data file I might need to add
! an initialization routine for densities and temperatures.
  
  end subroutine initialize_axisym_toroid_eq

!********************************************************************

  subroutine axisym_toroid_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)
! Input rvec = (x,y,z) position
! Outputs all quantities needed to fill an eq_point type
!
! Calls a selected specific magnetics routine which provides all magnetics quantities 
! needed for an eq_point type + flux, grad(flux), normalized flux and grad(normalized flux)
!
! Calls a general density and temperature routine which returns densities, temperatures
! and their gradients based on density_prof_model and temperature_prof_model

    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message
<<<<<<< HEAD

    use solovev_magnetics_m, only : solovev_magnetics
=======
>>>>>>> 47dd208b739c6b82af6ec9c77d62cf52c1c3da99
    
    implicit none

    real(KIND=rkind), intent(in) :: rvec(3) 
    real(KIND=rkind), intent(out) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind), intent(out) :: ns(0:nspec), gradns(3,0:nspec)
    real(KIND=rkind), intent(out) :: ts(0:nspec), gradts(3,0:nspec)
    character(len=20), intent(out) :: equib_err


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
    if (r < box_rmin .or. r > box_rmax) equib_err = 'R out_of_bounds'
    if (z < box_zmin .or. z > box_zmax) equib_err = 'z out_of_bounds'

<<<<<<< HEAD
    magnetics: select case (trim(magnetics_model))
       case ('solovev_magnetics')
          call solovev_magnetics(rvec, bvec, gradbtensor, psi, gradpsi, psiN, gradpsiN, equib_err)
=======
    magnetics: select case (trim(poloidal_flux_model))
       case ('solovev_magnetics')
          call solovev_magnetics(rvec, bvec, gradbtensor psi, gradpsi, psiN, gradpsiN, equib_err)
>>>>>>> 47dd208b739c6b82af6ec9c77d62cf52c1c3da99
    end select magnetics


!   Density profile.

<<<<<<< HEAD
    density: select case (trim(density_prof_model))
=======
    density: select case (trim(dens_prof_model))
>>>>>>> 47dd208b739c6b82af6ec9c77d62cf52c1c3da99

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
<<<<<<< HEAD
            write(0,*) 'axisym_toriod_eq: Unknown density_prof_model =', density_prof_model
=======
            write(0,*) 'axisym_toriod_eq: Unknown dens_prof_model =', dens_prof_model
>>>>>>> 47dd208b739c6b82af6ec9c77d62cf52c1c3da99
            stop 1

    end select density


!   Temperature profile.
    do is = 0, nspec
<<<<<<< HEAD
       temperature: select case (temperature_prof_model(is))
=======
       temperature: select case (t_prof_model(is))
>>>>>>> 47dd208b739c6b82af6ec9c77d62cf52c1c3da99


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
<<<<<<< HEAD
          write(0,*) 'axisym_toroid_eq: Unknown temperature_prof_model: ', &
                     & temperature_prof_model(0:nspec)
=======
          write(0,*) 'axisym_toroid_eq: Unknown t_prof_model: ', t_prof_model(0:nspec)
>>>>>>> 47dd208b739c6b82af6ec9c77d62cf52c1c3da99
          stop 1

       end select temperature
    end do

! Do some checking for invalid values

    if (minval(ns) < 0.) equib_err = 'negative_dens'
    if (minval(ts) < 0.) equib_err = 'negative_temp'

    return
 end subroutine axisym_toroid_eq

!********************************************************************

 subroutine axisym_toroid_psi(rvec, psi, gradpsi, psiN, gradpsiN)

    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message

<<<<<<< HEAD
    use solovev_magnetics_m, only : solovev_magnetics
=======
    use solovev_magnetics_m, only : solovev_magnetics_m
>>>>>>> 47dd208b739c6b82af6ec9c77d62cf52c1c3da99
    
    implicit none

    real(KIND=rkind), intent(in) :: rvec(3) 
    real(KIND=rkind), intent(out) :: psi, gradpsi(3), psiN, gradpsiN(3)


    magnetics: select case (trim(poloidal_flux_model))
       case ('solovev_magnetics')
          call solovev_psi(rvec, psi, gradpsi, psiN, gradpsiN)
    end select magnetics

    return
  
  end subroutine axisym_toroid_psi

!***********************************************************************
 subroutine write_axisym_toroid_profiles
    use species_m, only : nspec
    use diagnostics_m, only : message_unit

    implicit none
    
    real(KIND=rkind) :: dx, x
    real(KIND=rkind) :: rvec(3) 
    real(KIND=rkind) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind) :: ns(0:nspec), gradns(3,0:nspec)
    real(KIND=rkind) :: ts(0:nspec), gradts(3,0:nspec)
    real(KIND=rkind) :: psi, gradpsi(3), psiN, gradpsiN(3)
    character(len=20) :: equib_err

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
        call axisym_toroid_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)
        call axisym_toroid_psi(rvec, psi, gradpsi, psiN, gradpsiN)  
        write (message_unit,'(f11.5, a, e12.5, 3f11.5, f11.5, f11.5,  7f11.5)') &
               & x,'  ', ns(0), bvec, psi, psi/psiB, (ts(i), i=0, nspec)
!        write(*,*) 'x = ', x, '  gradpsi = ', gradpsi
!        write(*,*) 'gradbtensor = ', gradbtensor
    end do
        
 end subroutine write_axisym_toroid_profiles
 
end module axisym_toroid_eq_m
