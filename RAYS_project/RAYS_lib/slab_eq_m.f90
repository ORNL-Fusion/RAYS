module slab_eq_m

    use constants_m, only : rkind
    
    implicit none

! data for slab magnetics
    character(len=12) :: bx_prof_model, by_prof_model, bz_prof_model
    real(KIND=rkind) :: bx0, by0, bz0
    real(KIND=rkind) :: Ln_scale, LT_scale, LBy_shear_scale, LBz_scale

! data for slab density and temperature
    character(len=60) :: dens_prof_model
    real(KIND=rkind) :: rmaj, rmin
    real(KIND=rkind) :: alphan1
    real(KIND=rkind) :: alphan2  
    character(len=20), allocatable :: t_prof_model(:)
    real(KIND=rkind), allocatable :: alphat1(:)
    real(KIND=rkind), allocatable :: alphat2(:)
    
! data for boundary    
    real(KIND=rkind) :: xmin, xmax, ymin, ymax, zmin, zmax

 namelist /slab_eq_list/ &
     & bx_prof_model, by_prof_model, bz_prof_model, bx0, by0, bz0,                    &
     & rmaj, rmin, dens_prof_model, alphan1, alphan2, t_prof_model, alphat1, alphat2, &
     & Ln_scale, LT_scale, LBy_shear_scale, LBz_scale, &
     & xmin, xmax, ymin, ymax, zmin, zmax
     
!********************************************************************

contains

!********************************************************************

  subroutine initialize_slab_eq_m(read_input)

    use constants_m, only : input_unit
    use species_m, only : nspec
    use diagnostics_m, only : message_unit, verbosity
    
    implicit none
    logical, intent(in) :: read_input

    allocate( t_prof_model(0:nspec) )
    allocate( alphat1(0:nspec), alphat2(0:nspec) )

    if (read_input .eqv. .true.) then    
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, slab_eq_list)
        close(unit=input_unit)
        write(message_unit, slab_eq_list)
    end if

    if (verbosity > 2)  then
        call write_slab_profiles
    end if
    
    return
  end subroutine initialize_slab_eq_m

!********************************************************************

 subroutine slab_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)
!   A simple slab plasma equilibrium with B(x) = By(x)*ey + Bz(x)*ez, 
!   where ey and ez are unit vectors along y and z.
!   Note that bvec = B and gradbtensor(i,j) = d[B(j)]/d[x(i)].
!
!   For now Bx is constrained to be 0.  Non-zero B component in the direction of
!   stratification (i.e, x) complicates initialization solution of the dispersion
!   relation.  Several models are given for By(x), Bz(x).
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
    character(len=20), intent(out) :: equib_err

    real(KIND=rkind) :: x, y, z
    integer :: is

    equib_err = ''
    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    gradbtensor = 0.

! Check that we are in the box
    if (x < xmin .or. x > xmax) equib_err = 'x out_of_bounds'
    if (y < ymin .or. y > ymax) equib_err = 'y out_of_bounds'
    if (z < zmin .or. z > zmax) equib_err = 'z out_of_bounds'
    
    if (equib_err /= '') then
        write (message_unit, *) 'slab_eq:  equib_err ', equib_err
        return
    end if

!   Calculate Bx, and d(Bx)/dx.
    bx: select case (trim(bx_prof_model))

       case ('zero')
          bvec(1) = 0.

       case default
          write(0,*) 'SLAB: invalid bx_prof_model = ', bx_prof_model
          stop 1

    end select bx

!   Calculate By, and d(By)/dx.
    by: select case (trim(by_prof_model))

       case ('constant')
          bvec(2) = by0

       case ('toroid')
!         Same as Tokamak like magnetic field Bz. For diagnostics.
          bvec(2) = by0 / (1.+x/rmaj)
          gradbtensor(1,2) = -bvec(2) / (rmaj+x)

       case ('linear_shear')
!         Sheared By(x), By(0) = 0
          bvec(2) = by0 * x/LBy_shear_scale
          gradbtensor(1,2) = by0 / LBy_shear_scale

       case default
          write(0,*) 'SLAB: invalid by_prof_model = ', by_prof_model
          stop 1

    end select by

!   Calculate Bz, and d(Bz)/dx, d(Bz)/dy, d(Bz)/dz.
    bz: select case (trim(bz_prof_model))

       case ('constant')
          bvec(3) = bz0

       case ('toroid')
!         Tokamak like toroidal  field.
          bvec(3) = bz0 / (1.+x/rmaj)
          gradbtensor(1,3) = -bvec(3) / (rmaj+x)

       case ('linear')
!         Linear with scale length rmin
          bvec(3) = bz0 * (1.+x/LBz_scale)
          gradbtensor(1,3) = bz0/LBz_scale

       case default
          write(0,*) 'SLAB: invalid bz_prof_model = ', bz_prof_model
          stop 1

    end select bz

!   Density profile.
    
    density: select case (trim(dens_prof_model))

        case ('constant')
          ns(:nspec) = n0s(:nspec)
          gradns = 0.

        case ('linear')
!       Linear with scale length rmin (i.e. ns = n0s at x = 0, ns = 0 at x = -rmin)
            ns = n0s(:nspec)*(1.0 + x/Ln_scale)
            gradns = 0.
            gradns(1,0:nspec) = n0s(0:nspec) * (1.0/Ln_scale)

        case ('parabolic')
!       Parabolic around x = rmin
            ns(:nspec) = n0s(:nspec) * (1.-(x/rmin-1)**alphan2)**alphan1
            gradns = 0.
            gradns(1,:) = -n0s(:nspec)/rmin*alphan1*alphan2*(x/rmin-1)**(alphan2-1)* &
                         & (1-(1-(x/rmin-1)**alphan2)**(alphan1-1))
            gradns(2:3,:) = 0. 

        case ('Gaussian')
!         Gaussian (default: alphan1=1.).
            ns(:nspec) = n0s(:nspec) * exp(-3.*alphan1*(x/rmin)**2)
            gradns = 0.
            gradns(1,:nspec) = ns * (-6.*alphan1*x/rmin**2)

        case default
            write(0,*) 'SLAB: invalid dens_prof_model =', dens_prof_model
            stop 1

    end select density

!   Temperature profile.
    
    do is = 0, nspec
       temperature: select case (t_prof_model(is))

       case ('zero')
          ts(is) = 0.
          gradts(:,is) = 0.

       case ('constant')
          ts(is) = t0s(is)
          gradts(:,is) = 0.

        case ('linear')
!       Linear with scale length rmin
            ts(is)=t0s(is)*(1.+x/LT_scale)
            gradts(1,is) = t0s(is) * (1./LT_scale)

       case ('parabolic')
!         Parabolic around x = rmin
          ts(is) = t0s(is) * (1.-(x/rmin-1)**alphat2(is))**alphat1(is)
          gradts(1,is) = -t0s(is)/rmin*(alphat1(is)*alphat2(is)*(x/rmin-1)**(alphat2(is)-1)* &
                         & (1-(x/rmin-1)**alphat2(is))**(alphat1(is)-1))
          gradts(2:3,is) = 0.

       case default
          write(0,*) 'SLAB: invalid t_prof_model = ', t_prof_model(:nspec)
          stop 1

       end select temperature
    end do

! Do some checking for invalid values

    if (minval(ns) < 0.) equib_err = 'negative_dens'
    if (minval(ts) < 0.) equib_err = 'negative_temp'

    return
 end subroutine slab_eq


 subroutine write_slab_profiles
    use species_m, only : nspec
    use diagnostics_m, only : message_unit

    implicit none
    
    real(KIND=rkind) :: dx, x
    real(KIND=rkind) :: rvec(3) 
    real(KIND=rkind) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind) :: ns(0:nspec), gradns(3,0:nspec)
    real(KIND=rkind) :: ts(0:nspec), gradts(3,0:nspec)
    character(len=60) :: equib_err

    integer, parameter :: nx_points = 51
    integer :: ip, i
    
    character (len = *), parameter :: b9 = '         '
    character (len = *), parameter :: b10 = '          '
    character (len = *), parameter :: b12 = '            '
    
    dx = 2.*rmin/(nx_points-1)
    
    write (message_unit,*) '    x', b9,'ne', b12, 'bx', b9, 'by', b9, 'bz', b9, 'Te',b9, 'Ti(s)'

    do ip = 1, nx_points
        x = -rmin + (ip-1)*dx
        rvec( : ) = (/ x, real(0.,KIND=rkind), real(0.,KIND=rkind) /)
        call slab_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)
        write (message_unit,'(f11.5, a, e12.5, 3f11.5, 7f11.5)') &
               & x,'  ', ns(0), bvec, (ts(i), i=0, nspec)
    end do
        
 end subroutine write_slab_profiles
    
!********************************************************************

    subroutine finalize_slab_eq_m
		if (allocated(t_prof_model)) then
			deallocate( t_prof_model )
			deallocate( alphat1 )
		end if
		return
    end subroutine finalize_slab_eq_m

end module slab_eq_m
