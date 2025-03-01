module axisym_toroid_eq_m
! A general toroidal equilibrium module that supports a choice of magnetics, density,
! and temperature models.  It is assumed that there exists a poloidal flux function
! on which the densities and temperatures are constant. An external magnetics module is
! required for each magnetics model which provides B, grad B tensor, and the flux function,
! psi, along with it's gradient (poloidal_flux_model).  The magnetics model also provides
! certain geometry data on initialization.
!
! Input to axisym_toroid_eq() is rvec = (x,y,z) position
! Outputs from axisym_toroid_eq() are all quantities needed to fill an eq_point type
!
! Presently supported magnetics models are: 'solovev_magnetics', 'eqdsk_magnetics_lin_interp'
! and 'eqdsk_magnetics_spline_interp'
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
! Working notes:

    use constants_m, only : rkind, one, zero
    use species_m, only : nspec
    use quick_cube_splines_m, only : cube_spline_function_1D
    use density_spline_interp_m, only : initialize_density_spline_interp

    implicit none

! Local data **************************************************

! Namelist data for /axisym_toroid_eq_list/  **************************************************
! data for magnetics
    character(len=60) :: magnetics_model

! Geometry data
! Magnetic axis
    real(KIND=rkind) :: r_axis, z_axis
! data for bounding box
    real(KIND=rkind) :: box_rmin, box_rmax, box_zmin, box_zmax
! data for plasma boundary
    real(KIND=rkind) :: inner_bound, outer_bound, upper_bound, lower_bound
! Maximum value of psi considered to be inside of plasma.  Defaults to 1.0 but can be
! reset in namelist.
    real(KIND=rkind) :: plasma_psi_limit = one

! Data for density
    character(len=60) :: density_prof_model
    real(KIND=rkind) :: alphan1 ! parameters for parabolic model
    real(KIND=rkind) :: alphan2
! Density outside psi = 1 as a fraction of ne0, defaults to 0. but can be set in namelist
    real(KIND=rkind) :: d_scrape_off = zero

! Data for temperature
    character(len=60), allocatable :: temperature_prof_model(:)
    ! Parabolic model parameters
    real(KIND=rkind), allocatable :: alphat1(:) ! Can be dfferent for different species
    real(KIND=rkind), allocatable :: alphat2(:) ! Can be dfferent for different species
! Temperature outside psi = 1 as a fraction of Te0, defaults to 0. but can be set in namelist
    real(KIND=rkind) :: T_scrape_off = zero

	integer :: i, is

 namelist /axisym_toroid_eq_list/&
     & magnetics_model, &
     & plasma_psi_limit, &
     & density_prof_model, &
     & d_scrape_off, &
     & alphan1, alphan2, & ! parameters for parabolic density model
     & temperature_prof_model, &
     & alphat1, alphat2, & ! parameters for parabolic temperature model.  Same for all species now
     & T_scrape_off

!********************************************************************

contains

!********************************************************************

  subroutine initialize_axisym_toroid_eq_m(read_input)

    use constants_m, only : rkind
    use diagnostics_m, only : message, message_unit, text_message, messages_to_stdout, verbosity

	use species_m, only : nspec
    use solovev_magnetics_m, only : initialize_solovev_magnetics
    use eqdsk_magnetics_lin_interp_m, only : initialize_eqdsk_magnetics_lin_interp
    use eqdsk_magnetics_spline_interp_m, only : initialize_eqdsk_magnetics_spline_interp
    use density_spline_interp_m, only : initialize_density_spline_interp
    use temperature_spline_interp_m, only : initialize_temperature_spline_interp

    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder
	integer :: n_T_spline ! Number of species with splined temperature profiles
	integer :: i_spec_spline(0:nspec) ! Species number of any splined Ti profiles

    real(KIND=rkind) :: bp0

    call message(1)
    call text_message('initializing_axisym_toroid_eq ', 1)

	if (.not. allocated(temperature_prof_model)) then
		allocate( temperature_prof_model(0:nspec) )
		allocate( alphat1(0:nspec), alphat2(0:nspec) )
		temperature_prof_model = ' '
		alphat1 = 0.
		alphat2 = 0.
    end if

! Read input namelist
    if (read_input .eqv. .true.) then
    	input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, axisym_toroid_eq_list)
        close(unit=input_unit)
    end if

! Write input namelist
    if (verbosity > 0) then
		write(message_unit, axisym_toroid_eq_list)
		if (messages_to_stdout) write(*, axisym_toroid_eq_list)
    end if

    magnetics: select case (trim(magnetics_model))
       case ('solovev_magnetics')
          call initialize_solovev_magnetics(read_input, r_axis, z_axis, &
               & box_rmin, box_rmax, box_zmin, box_zmax, &
               & inner_bound, outer_bound, upper_bound, lower_bound)

        case ('eqdsk_magnetics_lin_interp')
           call initialize_eqdsk_magnetics_lin_interp(read_input, r_axis, z_axis, &
               & box_rmin, box_rmax, box_zmin, box_zmax, &
               & inner_bound, outer_bound, upper_bound, lower_bound)

        case ('eqdsk_magnetics_spline_interp')
           call initialize_eqdsk_magnetics_spline_interp(read_input, r_axis, z_axis, &
               & box_rmin, box_rmax, box_zmin, box_zmax, &
               & inner_bound, outer_bound, upper_bound, lower_bound)

          case default
          write(0,*) 'initialize_axisym_toroid_eq: unknown magnetics model =', magnetics_model
          call text_message('initialize_axisym_toroid_eq: unknown magnetics model',&
          & trim(magnetics_model),0)
          stop 1
    end select magnetics

! Density and Temperature.  Check for valid model name
! Density and temperature only need initialization for spline interpolation models

    density: select case (trim(density_prof_model))
        case ('density_spline_interp')
        	call initialize_density_spline_interp(read_input)
        case ('constant')
        case ('parabolic')
        case default
			if (verbosity > 0) then
				write(message_unit, *) 'axisym_toroid_eq: Unknown density_prof_model: ', &
                     & trim(density_prof_model)
				write(*, *) 'axisym_toroid_eq: Unknown density_prof_model: ', &
                     & trim(density_prof_model)
			end if
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
			if (verbosity > 0) then
				write(message_unit, *) 'axisym_toroid_eq: Unknown temperature_prof_model: ', &
                     & temperature_prof_model(is)
				write(*,*) 'axisym_toroid_eq: Unknown temperature_prof_model: ', &
                     & temperature_prof_model(is)
			end if
            stop 1
       end select temperature
    end do

! Initialize temperature splines if there are any
	if (n_T_spline > 0) call initialize_temperature_spline_interp(read_input)

  end subroutine initialize_axisym_toroid_eq_m

!********************************************************************

  subroutine axisym_toroid_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)
! Input rvec = (x,y,z) position
! Outputs all quantities needed to fill an eq_point type
!
! Calls a selected specific magnetics routine which provides all magnetics quantities
! needed for an eq_point type + flux, grad(flux), normalized flux and grad(normalized flux)
!
! Calls specific density and temperature routines which return densities, temperatures
! and their gradients based on density_prof_model and temperature_prof_model

    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message, verbosity

    use solovev_magnetics_m, only : solovev_magnetics
    use eqdsk_magnetics_lin_interp_m, only : eqdsk_magnetics_lin_interp
    use eqdsk_magnetics_spline_interp_m, only : eqdsk_magnetics_spline_interp
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
    real(KIND=rkind) :: psi, gradpsi(3), psiN, grad_psiN(3)
    real(KIND=rkind) :: dens, dd_psi
    real(KIND=rkind) :: t_prof, dt_dpsi
    real(KIND=rkind) :: Te, dTe_psi, Ti, dTi_psi
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
    if (r < box_rmin-Tiny .or. r > box_rmax+Tiny) then
    	equib_err = 'R_out_of_box'
    	write(*,*) 'R_out_of_box: R = ', R, '   box_rmin = ', box_rmin, '   box_rmax = ', box_rmax
    end if
    if (z < box_zmin-Tiny .or. z > box_zmax+Tiny) then
    	equib_err = 'Z_out_of_box'
    	write(*,*) 'Z_out_of_box: Z = ', Z, '   box_zmin = ', box_zmin, '   box_zmax = ', box_zmax
    end if

    if (equib_err /= '') return

    magnetics: select case (trim(magnetics_model))
       case ('solovev_magnetics')
          call solovev_magnetics(rvec, bvec, gradbtensor, psi, gradpsi, psiN, grad_psiN,&
            &  equib_err)

       case ('eqdsk_magnetics_lin_interp')
          call eqdsk_magnetics_lin_interp(rvec, bvec, gradbtensor, psi, gradpsi, psiN,&
              & grad_psiN, equib_err)

       case ('eqdsk_magnetics_spline_interp')
          call eqdsk_magnetics_spline_interp(rvec, bvec, gradbtensor, psi, gradpsi, psiN,&
             & grad_psiN, equib_err)

    end select magnetics

    ! Check that we are in the plasma. Set equib_err but don't stop
    if (psiN > plasma_psi_limit) equib_err = 'out_of_plasma'

!   Density profile.

    density: select case (trim(density_prof_model))

        case ('constant')
          ns(:nspec) = n0s(:nspec)
          gradns = 0.

        case ('parabolic')
            call parabolic_prof(psiN, d_scrape_off, alphan1, alphan2, dens, dd_psi)
			ns(0:nspec) = n0s(0:nspec) * dens
			gradns(1, 0:nspec) = n0s(0:nspec)*dd_psi*grad_psiN(1)
			gradns(2, 0:nspec) = n0s(0:nspec)*dd_psi*grad_psiN(2)
			gradns(3, 0:nspec) = n0s(0:nspec)*dd_psi*grad_psiN(3)

        case ('density_spline_interp')
            call density_spline_interp(psiN, d_scrape_off, dens, dd_psi)
			ns(0:nspec) = n0s(0:nspec) * dens
			gradns(1, 0:nspec) = n0s(0:nspec)*dd_psi*grad_psiN(1)
			gradns(2, 0:nspec) = n0s(0:nspec)*dd_psi*grad_psiN(2)
			gradns(3, 0:nspec) = n0s(0:nspec)*dd_psi*grad_psiN(3)

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
!      Parabolic: N.B. psi goes something like r**2 so if alphan2 = 1 and alphan1 = 2
!      the profile is pretty much parabolic
		  call parabolic_prof(psiN, T_scrape_off, alphat1(is), alphat2(is), t_prof, dt_dpsi)
		  ts(is) = t0s(is) * t_prof
		  gradts(1,is) = t0s(is)*dt_dpsi*grad_psiN(1)
		  gradts(2,is) = t0s(is)*dt_dpsi*grad_psiN(2)
		  gradts(3,is) = t0s(is)*dt_dpsi*grad_psiN(3)

       case ('temperature_spline_interp')
            call temperature_spline_interp(psiN, T_scrape_off, Te, dTe_psi, Ti, dTi_psi)
            if (is == 0) then
            	ts(is) = t0s(is) * Te
                gradts(1,is) = t0s(is)*dTe_psi*grad_psiN(1)
                gradts(2,is) = t0s(is)*dTe_psi*grad_psiN(2)
                gradts(3,is) = t0s(is)*dTe_psi*grad_psiN(3)
        	else
            	ts(is) = t0s(is) * Ti
                gradts(1,is) = t0s(is)*dTi_psi*grad_psiN(1)
                gradts(2,is) = t0s(is)*dTi_psi*grad_psiN(2)
                gradts(3,is) = t0s(is)*dTi_psi*grad_psiN(3)
            end if

       end select temperature
    end do

! Do some checking for invalid values

    if (minval(ns) < 0.) equib_err = 'negative_dens'
    if (minval(ts) < 0.) equib_err = 'negative_temp'

    return
 end subroutine axisym_toroid_eq

!********************************************************************

 subroutine axisym_toroid_psi(rvec, psi, gradpsi, psiN, grad_psiN)
 ! Returns poloidal flux -> psi(x,y,z) etc

    use constants_m, only : rkind
    use species_m, only : nspec, n0s, t0s
    use diagnostics_m, only : message_unit, message

    use solovev_magnetics_m, only : solovev_magnetics_psi
    use eqdsk_magnetics_lin_interp_m, only : eqdsk_magnetics_lin_interp_psi
    use eqdsk_magnetics_spline_interp_m, only : eqdsk_magnetics_spline_interp_psi

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: psi, gradpsi(3), psiN, grad_psiN(3)


    magnetics: select case (trim(magnetics_model))
       case ('solovev_magnetics')
          call solovev_magnetics_psi(rvec, psi, gradpsi, psiN, grad_psiN)

       case ('eqdsk_magnetics_lin_interp')
         call  eqdsk_magnetics_lin_interp_psi(rvec, psi, gradpsi, psiN, grad_psiN)

       case ('eqdsk_magnetics_spline_interp')
         call  eqdsk_magnetics_spline_interp_psi(rvec, psi, gradpsi, psiN, grad_psiN)
    end select magnetics

    return

  end subroutine axisym_toroid_psi
!********************************************************************

 subroutine axisym_toroid_rho(rvec, rho, gradrho)
 ! Returns normalized square root toroidal flux -> rho(x,y,z).
 ! N.B. Unlike axisym_toroid_psi rho is already normalized so there rare no rhoN or
 ! grad_rhoN arguments
 ! N.B. so far this is only implemented in eqdsk_magnetics_spline_interp

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message

!     use solovev_magnetics_m, only : solovev_magnetics_rho
!     use eqdsk_magnetics_lin_interp_m, only : eqdsk_magnetics_lin_interp_rho
    use eqdsk_magnetics_spline_interp_m, only : eqdsk_magnetics_spline_interp_rho

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: rho, gradrho(3)


    magnetics: select case (trim(magnetics_model))
!        case ('solovev_magnetics')
!           call solovev_magnetics_rho(rvec, rho, gradrho, rhoN, grad_rhoN)
!
!        case ('eqdsk_magnetics_lin_interp')
!          call  eqdsk_magnetics_lin_interp_rho(rvec, rho, gradrho, rhoN, grad_rhoN)

       case ('eqdsk_magnetics_spline_interp')
         call  eqdsk_magnetics_spline_interp_rho(rvec, rho, gradrho)

        case default
          write(0,*) 'axisym_toroid_rho: unknown magnetics model =', magnetics_model
          call text_message('axisym_toroid_rho: unknown magnetics model',&
          & trim(magnetics_model),0)
          stop 1

    end select magnetics

    return
  end subroutine axisym_toroid_rho

!***********************************************************************
 subroutine write_axisym_toroid_profiles
    use constants_m, only : rkind
    use species_m, only : nspec
    use diagnostics_m, only : message_unit, messages_to_stdout

    implicit none

    real(KIND=rkind) :: dx, x
    real(KIND=rkind) :: rvec(3)
    real(KIND=rkind) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind) :: ns(0:nspec), gradns(3,0:nspec)
    real(KIND=rkind) :: ts(0:nspec), gradts(3,0:nspec)
    real(KIND=rkind) :: psi, gradpsi(3), psiN, grad_psiN(3)
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
            & b8, 'psiN', b8,  'Te',b9, 'Ti(s)'
    if (messages_to_stdout) then
		write (*,*) '    x', b9,'ne', b12, 'bx', b9, 'by', b9, 'bz', b9, 'psi', &
				& b8, 'psiN', b8,  'Te',b9, 'Ti(s)'
    end if

    do ip = 1, nx_points
        x = inner_bound + (ip-1)*dx
        rvec( : ) = (/ x, real(0.,KIND=rkind), real(0.,KIND=rkind) /)
        call axisym_toroid_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)
        call axisym_toroid_psi(rvec, psi, gradpsi, psiN, grad_psiN)
        write (message_unit,'(f11.5, a, e12.5, 3f11.5, f11.5, f11.5,  7f11.5)') &
               & x,'  ', ns(0), bvec, psi, psiN, (ts(i), i=0, nspec)
		if (messages_to_stdout) then
			write (*,'(f11.5, a, e12.5, 3f11.5, f11.5, f11.5,  7f11.5)') &
				   & x,'  ', ns(0), bvec, psi, psiN, (ts(i), i=0, nspec)
		end if

!        write(*,*) 'x = ', x, '  gradpsi = ', gradpsi
!        write(*,*) 'gradbtensor = ', gradbtensor
    end do

 end subroutine write_axisym_toroid_profiles


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

    subroutine deallocate_axisym_toroid_eq_m
		use solovev_magnetics_m, only : deallocate_solovev_magnetics_m
		use eqdsk_magnetics_lin_interp_m, only : deallocate_eqdsk_magnetics_lin_interp_m
		use eqdsk_magnetics_spline_interp_m, only : deallocate_eqdsk_magnetics_spline_interp_m

		if (allocated(temperature_prof_model)) then
			deallocate( temperature_prof_model )
			deallocate( alphat1 )
			deallocate( alphat2 )
		end if

		call deallocate_solovev_magnetics_m
		call deallocate_eqdsk_magnetics_lin_interp_m
		call deallocate_eqdsk_magnetics_spline_interp_m
		return

    end subroutine deallocate_axisym_toroid_eq_m

!********************************************************************

end module axisym_toroid_eq_m
