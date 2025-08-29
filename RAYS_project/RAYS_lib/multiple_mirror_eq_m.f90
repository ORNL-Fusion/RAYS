module multiple_mirror_eq_m
!
! A general axisymmetric mirror equilibrium module that supports a choice of magnetics,
! density, and temperature models.  It is assumed that there exists a flux function, Aphi,
! Aphi = Axial magnetic flux inside radius r.  Presently densities and temperatures are
! assumed constant on Aphi surfaces.  In the future axial variation of could be introduced.
! An external magnetics module is required for each specific magnetics model,
! which provides B, grad B tensor, and the flux function, Aphi. The magnetics model also
! provides certain geometry data on initialization such as the r,z box dimensions on which
! the equilibrium is defined, and the r,z position and Aphi at the limiter strike point of
! the last un-interrupted flux surface, LUFS.
!
! Input to multiple_mirror_eq() is rvec = (x,y,z) position
! Outputs from multiple_mirror_eq() are all the quantities needed to fill an eq_point type.
!
! Presently the only supported magnetics model is: mirror_magnetics_spline_interp
! although someday may add an analytic form.
!
! The radial coordinate for density and temperature profiles is AphiN = Aphi/Aphi_LUFS
! which goes like r**2 near the axis.
!
! Presently supported density models: 'constant', 'parabolic', 'density_spline_interp'
! All species presently use the same density profile model.  Otherwise charge neutrality is
! complicated.
!
! Presently supported temperature models: 'zero' 'constant', 'parabolic',
! and 'temperature_spline_interp'
! Different species can have different profile models. Although presently only two different
! spline profiles are supported: electron and ion.  All splined ion profiles are the same.
!
! Note on boundaries:
!   Outermost boundary is defined by box_rmax, box_zmin, box_zmax.  These are in the
!   initialization for the specific magnetic model implementation (which for now is only
!   splines).  The init for mirror_magnetics_spline_interp_m sets these to zero, r_max,
!   z_min, z_max respectively, which are contained in the magnetic field netCDF file.
!
!   Plasma boundary is defined in terms of normalized Aphi on the last un-interrupted
!   flux surface, == Aphi_LUFS.
!
! Working notes:

    use constants_m, only : rkind, one, zero, two
    use species_m, only : nspec
    use mirror_magnetics_spline_interp_m, only : mirror_field_NC_file, &
        & initialize_mirror_magnetics_spline_interp
    use quick_cube_splines_m, only : cube_spline_function_1D
    use density_spline_interp_m, only : initialize_density_spline_interp

    implicit none

! Local data **************************************************

! Namelist data for /multiple_mirror_eq_list/  **************************************************
    character(len=60) :: magnetics_model

! data for magnetics
! Geometry data
    real(KIND=rkind) :: box_rmax, box_zmin, box_zmax

	! Get r_LUFS, z_LUFS, Aphi_LUFS from specific mirror magnetics routine (e.g.
	! initialize_mirror_magnetics_spline_interp)
    real(KIND=rkind) :: r_LUFS, z_LUFS ! Location of scrape-off point of last flux surface
    real(KIND=rkind) :: Aphi_LUFS ! Bounding value of Aphi evaluated at r_LUFS, z_LUFS

! Maximum value of normalized Aphi considered to be inside of plasma.  Defaults to 1.0
! but can be reset in namelist.
    real(KIND=rkind) :: plasma_AphiN_limit = one

! Data for density
    character(len=60) :: density_prof_model
    real(KIND=rkind) :: alphan1, alphan2 ! parameters for parabolic model
    real(KIND=rkind) :: AphiN0_d, delta_d ! parameters for hyperbolic model
! Density outside plasma_AphiN_limit = 1 as a fraction of ne0. Defaults to 0. but can
! be reset in namelist
    real(KIND=rkind) :: d_scrape_off = zero

! Data for temperature
    character(len=60), allocatable :: temperature_prof_model(:)
    ! Parabolic model parameters
    real(KIND=rkind), allocatable :: alphat1(:), alphat2(:) ! Can be dfferent for different species
    real(KIND=rkind), allocatable :: Aphin0_t(:),delta_t(:) ! Can be dfferent for different species
! Temperature outside rho = 1 as a fraction of Te0, defaults to 0. but can be set in namelist
    real(KIND=rkind) :: T_scrape_off = zero

	integer :: i, is

 namelist /multiple_mirror_eq_list/&
	 & magnetics_model, &
     & plasma_AphiN_limit, &
     & density_prof_model, &
     & d_scrape_off, &
     & alphan1, alphan2, & ! parameters for parabolic density model
     & Aphin0_d, delta_d, & ! parameters for hyperbolic density model
     & temperature_prof_model, &
     & alphat1, alphat2, & ! parameters for parabolic temperature model.
     & Aphin0_t, delta_t, & ! parameters for hyperbolic temperature model.
     & T_scrape_off

!********************************************************************

contains

!********************************************************************

  subroutine initialize_multiple_mirror_eq_m(read_input)

    use constants_m, only : rkind
    use diagnostics_m, only : message, message_unit, text_message, messages_to_stdout, verbosity

	use species_m, only : nspec
	use mirror_magnetics_spline_interp_m, only : initialize_mirror_magnetics_spline_interp, &
	               & r_LUFS_spline, z_LUFS_spline, Aphi_LUFS_spline
    use density_spline_interp_m, only : initialize_density_spline_interp
    use temperature_spline_interp_m, only : initialize_temperature_spline_interp

    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder
	integer :: n_T_spline ! Number of species with splined temperature profiles
	integer :: i_spec_spline(0:nspec) ! Species number of any splined Ti profiles

    real(KIND=rkind) :: bp0

    call message(1)
    call text_message('initializing_multiple_mirror_eq ', 1)

	if (.not. allocated(temperature_prof_model)) then
		allocate( temperature_prof_model(0:nspec) )
		allocate( alphat1(0:nspec), alphat2(0:nspec) )
		allocate( Aphin0_t(0:nspec), delta_t(0:nspec) )
		temperature_prof_model = ' '
		alphat1 = 0.
		alphat2 = 0.
    end if

! Read input namelist
    if (read_input .eqv. .true.) then
    	input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, multiple_mirror_eq_list)
        close(unit=input_unit)
    end if

! Write input namelist
    if (verbosity > 0) then
		write(message_unit, multiple_mirror_eq_list)
		if (messages_to_stdout) write(*, multiple_mirror_eq_list)
    end if

! magnetics (For, now only have spline interpolation.  Add analytic fields someday soon.)
    magnetics: select case (trim(magnetics_model))
       case ('mirror_magnetics_spline_interp')
          call initialize_mirror_magnetics_spline_interp(read_input, &
               & box_rmax, box_zmin, box_zmax)
          r_LUFS = r_LUFS_spline
          z_LUFS = z_LUFS_spline
          Aphi_LUFS = Aphi_LUFS_spline

          case default
          write(0,*) 'initialize_multiple_mirror_eq: unknown magnetics model =', magnetics_model
          call text_message('initialize_multiple_mirror_eq: unknown magnetics model',&
          & trim(magnetics_model),0)
          stop 1
    end select magnetics

	if (verbosity > 0) then
		write(message_unit, *) 'multiple_mirror_eq: r_LUFS = ', r_LUFS,  '  z_LUFS = ', &
		                      & z_LUFS, '   Aphi_LUFS =  ', Aphi_LUFS
		write(*, *)  'multiple_mirror_eq: r_LUFS = ', r_LUFS,  '  z_LUFS = ', &
		                      & z_LUFS, '   Aphi_LUFS =  ', Aphi_LUFS
	end if

! Density and Temperature.  Check for valid model name
! Density and temperature only need initialization for spline interpolation models

    density: select case (trim(density_prof_model))
        case ('density_spline_interp')
        	call initialize_density_spline_interp(read_input)
        case ('constant')
        case ('parabolic')
        case ('hyperbolic')
        case default
				write(message_unit, *) 'multiple_mirror_eq: Unknown density_prof_model: ', &
                     & trim(density_prof_model)
				write(*, *) 'multiple_mirror_eq: Unknown density_prof_model: ', &
                     & trim(density_prof_model)
            stop 1
    end select density

! Temperature.  Check for valid model name and count number of splined profiles
	n_T_spline = 0
    do is = 0, nspec
       temperature: select case( trim(temperature_prof_model(is)) )
        case ('temperature_spline_interp')
			n_T_spline = n_T_spline + 1
			i_spec_spline(n_T_spline) = is
        case('zero')
        case('constant')
        case('parabolic')
        case default
			write(message_unit, *) 'multiple_mirror_eq: Unknown temperature_prof_model: ', &
				 & temperature_prof_model(is)
			write(*,*) 'multiple_mirror_eq: Unknown temperature_prof_model: ', &
				 & temperature_prof_model(is)
            stop 1
       end select temperature
    end do

! Initialize temperature splines if there are any
	if (n_T_spline > 0) call initialize_temperature_spline_interp(read_input)

  end subroutine initialize_multiple_mirror_eq_m

!********************************************************************

  subroutine multiple_mirror_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)
! Input rvec = (x,y,z) position
! Outputs all quantities needed to fill an eq_point type
!
! Calls a selected specific magnetics routine which provides all magnetics quantities
! needed for an eq_point type + flux, grad(flux), normalized flux and grad(normalized flux)
!
! Calls specific density and temperature routines which return densities, temperatures
! and their gradients based on density_prof_model and temperature_prof_model
!
! N.B. so far this is only implemented for mirror_magnetics_spline_interp
! although someday may add an analytic form

    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message, verbosity

	use mirror_magnetics_spline_interp_m, only : mirror_magnetics_spline_interp
    use density_spline_interp_m, only : density_spline_interp
    use temperature_spline_interp_m, only : temperature_spline_interp

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind), intent(out) :: ns(0:nspec), gradns(3,0:nspec)
    real(KIND=rkind), intent(out) :: ts(0:nspec), gradts(3,0:nspec)
    character(len=60), intent(out) :: equib_err


    real(KIND=rkind) :: x, y, z, r
    real(KIND=rkind) :: br, bz, bphi, bp0
    real(KIND=rkind) ::  Aphi, gradAphi(3), AphiN, gradAphiN(3)
    real(KIND=rkind) :: dens, dd_rho
    real(KIND=rkind) :: t_prof, dt_drho
    real(KIND=rkind) :: Te, dTe_rho, Ti, dTi_rho
    real(KIND=rkind), parameter :: Tiny = 10.0e-14_rkind
    integer :: is

!*****************************************************************************
    equib_err = ''
    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    r = sqrt(x**2+y**2)

! Check that we are in the box.  But so we don't get crash when evaluating on the
! box boundary, allow a leeway of 2*tiny(x)
    if (r > box_rmax) then
    	equib_err = 'R_out_of_box'
    	write(*,*) 'R_out_of_box: R = ', R, '   box_rmax = ', box_rmax
    end if
    if (z < box_zmin .or. z > box_zmax) then
    	equib_err = 'Z_out_of_box'
    	write(*,*) 'Z_out_of_box: Z = ', Z, '   box_zmin = ', box_zmin, '   box_zmax = ', box_zmax
    end if

    if (equib_err /= '') return

    magnetics: select case (trim(magnetics_model))
       case ('mirror_magnetics_spline_interp')
			call mirror_magnetics_spline_interp(rvec, bvec, gradbtensor, Aphi, gradAphi,&
                           & AphiN, gradAphiN, equib_err)
    end select magnetics

! Check that we are in the plasma. Set equib_err but don't stop
   if (AphiN > plasma_AphiN_limit) equib_err = 'out_of_plasma'

!   Density profile.

    density: select case (trim(density_prof_model))

        case ('constant')
          ns(:nspec) = n0s(:nspec)
          gradns = 0.

        case ('parabolic')
!      Parabolic: N.B. AphiN goes something like r**2 so if alphan2 = 1 and alphan1 = 2
!      the profile is pretty much parabolic
            call parabolic_prof(AphiN, d_scrape_off, alphan1, alphan2, dens, dd_rho)
			ns(0:nspec) = n0s(0:nspec) * dens
			gradns(1, 0:nspec) = n0s(0:nspec)*dd_rho*gradAphiN(1)
			gradns(2, 0:nspec) = n0s(0:nspec)*dd_rho*gradAphiN(2)
			gradns(3, 0:nspec) = n0s(0:nspec)*dd_rho*gradAphiN(3)

        case ('hyperbolic')
            call hyperbolic_prof(AphiN, d_scrape_off, AphiN0_d, delta_d, dens, dd_rho)
			ns(0:nspec) = n0s(0:nspec) * dens
			gradns(1, 0:nspec) = n0s(0:nspec)*dd_rho*gradAphiN(1)
			gradns(2, 0:nspec) = n0s(0:nspec)*dd_rho*gradAphiN(2)
			gradns(3, 0:nspec) = n0s(0:nspec)*dd_rho*gradAphiN(3)

        case ('density_spline_interp')
            call density_spline_interp(AphiN, d_scrape_off, dens, dd_rho)
			ns(0:nspec) = n0s(0:nspec) * dens
			gradns(1, 0:nspec) = n0s(0:nspec)*dd_rho*gradAphiN(1)
			gradns(2, 0:nspec) = n0s(0:nspec)*dd_rho*gradAphiN(2)
			gradns(3, 0:nspec) = n0s(0:nspec)*dd_rho*gradAphiN(3)

    end select density

!   Temperature profile.
   ts = 0.
   gradts = 0.
   do is = 0, nspec
       temperature: select case( trim(temperature_prof_model(is)) )

        case ('zero')
          ts(is) = 0.
          gradts(:,is) = 0.

        case ('constant')
          ts(is) = t0s(is)
          gradts = 0.

       case ('parabolic')
!      Parabolic: N.B. AphiN goes something like r**2 so if alphan2 = 1 and alphan1 = 2
!      the profile is pretty much parabolic
		  call parabolic_prof(AphiN, T_scrape_off, alphat1(is), alphat2(is), t_prof, dt_drho)
		  ts(is) = t0s(is) * t_prof
		  gradts(1,is) = t0s(is)*dt_drho*gradAphiN(1)
		  gradts(2,is) = t0s(is)*dt_drho*gradAphiN(2)
		  gradts(3,is) = t0s(is)*dt_drho*gradAphiN(3)

       case ('hyperbolic')
		  call hyperbolic_prof(AphiN, T_scrape_off, Aphin0_t(is),delta_t(is), t_prof, dt_drho)
		  ts(is) = t0s(is) * t_prof
		  gradts(1,is) = t0s(is)*dt_drho*gradAphiN(1)
		  gradts(2,is) = t0s(is)*dt_drho*gradAphiN(2)
		  gradts(3,is) = t0s(is)*dt_drho*gradAphiN(3)

       case ('temperature_spline_interp')
            call temperature_spline_interp(AphiN, T_scrape_off, Te, dTe_rho, Ti, dTi_rho)
            if (is == 0) then
            	ts(is) = t0s(is) * Te
                gradts(1,is) = t0s(is)*dTe_rho*gradAphiN(1)
                gradts(2,is) = t0s(is)*dTe_rho*gradAphiN(2)
                gradts(3,is) = t0s(is)*dTe_rho*gradAphiN(3)
        	else
            	ts(is) = t0s(is) * Ti
                gradts(1,is) = t0s(is)*dTi_rho*gradAphiN(1)
                gradts(2,is) = t0s(is)*dTi_rho*gradAphiN(2)
                gradts(3,is) = t0s(is)*dTi_rho*gradAphiN(3)
            end if

       end select temperature
    end do

! Do some checking for invalid values

    if (minval(ns) < 0.) equib_err = 'negative_dens'
    if (minval(ts) < 0.) equib_err = 'negative_temp'

    return
 end subroutine multiple_mirror_eq

!********************************************************************

 subroutine multiple_mirror_Aphi(rvec, Aphi, gradAphi, AphiN, gradAphiN)
 ! Returns Aphi (a.k.a axial flux) and normalized Aphi/Aphi_LUFS + gradients
 !
 ! N.B. so far this is only implemented for mirror_magnetics_spline_interp
 ! although someday may add an analytic form

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message

    use mirror_magnetics_spline_interp_m, only : mirror_magnetics_spline_interp_aphi

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: Aphi, gradAphi(3), AphiN, gradAphiN(3)

	call mirror_magnetics_spline_interp_Aphi(rvec, Aphi, gradAphi, AphiN, gradAphiN)

    return
  end subroutine multiple_mirror_Aphi

!***********************************************************************
!  subroutine write_multiple_mirror_profiles
!     use constants_m, only : rkind
!     use species_m, only : nspec
!     use diagnostics_m, only : message_unit, messages_to_stdout
!
!     implicit none
!
!     real(KIND=rkind) :: dx, x
!     real(KIND=rkind) :: rvec(3)
!     real(KIND=rkind) :: bvec(3), gradbtensor(3,3)
!     real(KIND=rkind) :: ns(0:nspec), gradns(3,0:nspec)
!     real(KIND=rkind) :: ts(0:nspec), gradts(3,0:nspec)
!     real(KIND=rkind) :: rho, gradrho(3), rhoN, grad_rhoN(3)
!     character(len=60) :: equib_err
!
!     integer, parameter :: nx_points = 51
!     integer :: ip, i
!
!     character (len = *), parameter :: b3 = '   '
!     character (len = *), parameter :: b4 = '    '
!     character (len = *), parameter :: b8 = '        '
!     character (len = *), parameter :: b9 = '         '
!     character (len = *), parameter :: b10 = '          '
!     character (len = *), parameter :: b12 = '            '
!
!     dx = (outer_bound-inner_bound)/(nx_points-1)
!
!     write (message_unit,*) '    x', b9,'ne', b12, 'bx', b9, 'by', b9, 'bz', b9, 'rho', &
!             & b8, 'rhoN', b8,  'Te',b9, 'Ti(s)'
!     if (messages_to_stdout) then
! 		write (*,*) '    x', b9,'ne', b12, 'bx', b9, 'by', b9, 'bz', b9, 'rho', &
! 				& b8, 'rhoN', b8,  'Te',b9, 'Ti(s)'
!     end if
!
!     do ip = 1, nx_points
!         x = inner_bound + (ip-1)*dx
!         rvec( : ) = (/ x, real(0.,KIND=rkind), real(0.,KIND=rkind) /)
!         call multiple_mirror_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)
!         call multiple_mirror_rho(rvec, rho, gradrho, rhoN, grad_rhoN)
!         write (message_unit,'(f11.5, a, e12.5, 3f11.5, f11.5, f11.5,  7f11.5)') &
!                & x,'  ', ns(0), bvec, rho, rhoN, (ts(i), i=0, nspec)
! 		if (messages_to_stdout) then
! 			write (*,'(f11.5, a, e12.5, 3f11.5, f11.5, f11.5,  7f11.5)') &
! 				   & x,'  ', ns(0), bvec, rho, rhoN, (ts(i), i=0, nspec)
! 		end if
!
! !        write(*,*) 'x = ', x, '  gradrho = ', gradrho
! !        write(*,*) 'gradbtensor = ', gradbtensor
!     end do
!
!  end subroutine write_multiple_mirror_profiles


!********************************************************************

  pure subroutine parabolic_prof(rho, f_min, alpha1, alpha2, f, fp)
! Provides a parabolic-like function f(rho) and its derivative fp(rho)
!
! Parabolic-like means:
! f = 1.0 at rho = 0
! f = (1 - rho**alpha2)**alpha1 provided f > f_min
! f = f_min if the parabola above would give f < f_min or if rho > 1.0

    implicit none

    real(KIND=rkind), intent(in) :: rho, f_min, alpha1, alpha2
    real(KIND=rkind), intent(out) :: f, fp

	f = zero
	if (rho < one) then
		f = (1.-rho**(alpha2))**alpha1
		fp = -alpha1*alpha2*rho**(alpha2 - 1.)*(1.-rho**(alpha2))&
				& **(alpha1 - 1.)
	end if

	if (f < f_min) then
		f = f_min
		fp = zero
	end if

  end subroutine parabolic_prof

!********************************************************************

  pure subroutine hyperbolic_prof(rho, f_min, rho0, delta, f, fp)
! Provides a hyperbolic-like function f(rho) and its derivative fp(rho)
!
! f_min is the floor, rho0 is the profile inflection point and delta is a gradient
! scale length

    implicit none

    real(KIND=rkind), intent(in) :: rho, f_min, rho0, delta
    real(KIND=rkind), intent(out) :: f, fp

	f = zero
	f = (tanh((rho+rho0)/delta)-tanh((rho-rho0)/delta))/two/tanh((rho0)/delta)
	fp = (one/cosh((rho+rho0)/delta)**2 - one/cosh((rho-rho0)/delta)**2)/(two*delta)/&
		& tanh((rho0)/delta)

	f = (one - f_min)*f + f_min
	fp = (one - f_min)*fp

  end subroutine hyperbolic_prof

!********************************************************************

  pure subroutine hyperbolic_prof_inside_LUFS(rho, f_min, rho0, delta, f, fp)
! Provides a hyperbolic-like function f(rho) and its derivative fp(rho)
!
! hyperbolic-like means:
! f = 1.0 at rho = 0
! f = hyperbolic provided f > f_min
! f = f_min if the hyperbolic above would give f < f_min or if rho > 1.0 (i.e. outside LUFS)
! i.e. f_min is the floor, rho0 is the profile inflection point, and delta is a gradient
!      scale length

    implicit none

    real(KIND=rkind), intent(in) :: rho, f_min, rho0, delta
    real(KIND=rkind), intent(out) :: f, fp

	f = zero
	if (rho < one) then
		f = (tanh((rho+rho0)/delta)-tanh((rho-rho0)/delta))/two/tanh((rho0)/delta)
		fp = (one/cosh((rho+rho0)/delta)**2 - one/cosh((rho-rho0)/delta)**2)/(two*delta)/&
		    & tanh((rho0)/delta)
	end if

	if (f < f_min) then
		f = f_min
		fp = zero
	end if

  end subroutine hyperbolic_prof_inside_LUFS


!********************************************************************

  subroutine deallocate_multiple_mirror_eq_m
	use mirror_magnetics_spline_interp_m, only : deallocate_mirror_magnetics_spline_interp_m

	if (allocated(temperature_prof_model)) then
		deallocate( temperature_prof_model )
		deallocate( alphat1 )
		deallocate( alphat2 )
	end if

	call deallocate_mirror_magnetics_spline_interp_m
	return

  end subroutine deallocate_multiple_mirror_eq_m

!********************************************************************

end module multiple_mirror_eq_m
