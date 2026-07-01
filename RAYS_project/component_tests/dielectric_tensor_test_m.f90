module dielectric_tensor_test_m
! Tests various dielectric tensors.  Evaluates susceptibility tensor for each species
! and dielectric tensor (i.e. susceptibility tensor summed over species and with
! k x k term added).  These are obtained from module suscep_m.
!
! There are two kinds of test:
! 1) one_point -> Evaluates dielectric tensor at one equilibrium point for one pair of
!    (nx, nz) for one or more dispersion models
! 2) Scans through a series of position x points and gets the equilibrium at those points.
!    Then cycles through a series of nx, nz values and evaluates eps(x, nx, nz)
! N.B. Scan stuff is commented out.  I'll fix it if it is ever needed.
!
! Dispersion tensor models supported are cold and warm_bessel
!
! N.B. For convenience refractive indices are input rather than k values

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use diagnostics_m, only : message_unit, message, text_message, verbosity
    use constants_m, only : rkind, zero, one, pi
	use write_matrix_m, only : write_matrix
     implicit none

! Switches to select test kind
    logical :: one_point = .true.
    logical :: scan = .false.  ! default, can be reset on input

! Switches to select dispersion models
    logical :: test_cold = .true.
    logical :: test_warm_bessel = .false.  ! default, can be reset on input

! Data for one_point test
	integer :: number_of__models ! As of now 1 or 2
    real(KIND=rkind) :: x_in, y_in, z_in ! Eq point to evaluate
    real(KIND=rkind) :: n1_in, n3_in ! n perp and n|| to evaluate
    complex(KIND=rkind) :: n1_c, n3_c ! n perp and n|| to evaluate

! Data for x scan and nz scan
    integer :: n_xpoints ! Number of x points in scan
    real(KIND=rkind) :: xmin, xmax ! min and max of x scan
    integer :: n_nz_points ! Number of points in nz scan
    real(KIND=rkind) :: nz_min, nz_max ! min and max of nz scan

! Switches to select specific parameter scans
    logical :: scan_x = .false.
    logical :: scan_nz = .false.

    namelist /dielectric_tensor_test_list/ one_point, scan, &
            & test_cold, test_warm_bessel, &
            & x_in, y_in, z_in, n1_in, n3_in, &
            & scan_x, scan_nz, &
            & n_xpoints, xmin, xmax, n_nz_points, nz_min, nz_max

 !*************************************************************************

 contains

 subroutine init_dielectric_tensor_test(read_input)

    use diagnostics_m, only : message_unit, message, text_message, run_label

    implicit none

    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder

    if (read_input .eqv. .true.) then
    ! Read and write input namelist
   		input_unit = get_unit_number()
        open(unit=input_unit, file='component_test_rays.in',action='read', status='old',&
           & form='formatted')
        read(input_unit, dielectric_tensor_test_list)
        close(unit=input_unit)
        write(message_unit, dielectric_tensor_test_list)
    end if

    call text_message('Finished initialize_dielectric_tensor_test ')

    return
 end subroutine init_dielectric_tensor_test

 !*************************************************************************

  subroutine dielectric_tensor_test

	use constants_m, only : e
    use diagnostics_m, only : message_unit, message, text_message, messages_to_stdout
    use equilibrium_m, only : equilibrium, eq_point
    use rf_m, only : ray_dispersion_model, wave_mode, k0
    use species_m, only : nspec, spec_name, spec_model
    use suscep_m, only : suscep_cold, dielectric_cold, &
              & suscep_bessel_n1_n3_real, dielectric_general_n1_n3_real, &
              & suscep_bessel_n1_cmplx_n3_real, dielectric_general_n1_cmplx_n3_real, &
              & disp_fun_cold_real, disp_fun_general_real

    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

    complex(KIND=rkind) :: chi(3,3), eps(3,3)

! Results
    real(KIND=rkind) :: x, x_vec(n_xpoints), nz_vec(n_nz_points)
    complex(KIND=rkind) :: nxsq_vec(n_nz_points), nsq_vec(n_nz_points)

    real(KIND=rkind) :: rvec(3)
    real(KIND=rkind) :: dx, dnz, dtheta, nxr, nzr
    integer :: get_unit_number, out_unit
    external get_unit_number
    integer :: i, is, ix, inz, itheta
    integer :: k_sign = 1 ! Choose positive sign for k_x
    character(len=80) :: label

 !*************************************************************************

    ! Open file to put write(*,*) data in
! 	out_unit = message_unit
	out_unit = get_unit_number()
	open (unit = out_unit, file = 'dielectric_tensor_test_m.out', status = 'unknown')


    if (one_point) then

      	rvec = (/x_in, y_in, z_in/)
      	call equilibrium(rvec, eq)

		write(out_unit, *) ' '
		write(out_unit, *) 'rvec = ', rvec
		write(out_unit, *) 'n1_in = ', n1_in, '  n3_in = ', n3_in
		write(out_unit, *) 'species ', spec_name(0:nspec)
		write(out_unit, *) '|B| = ', eq%bmag
		write(out_unit, *) 'n(s) = ', eq%ns
		write(out_unit, *) 'T(s) ev = ', eq%ts/e
		write(*, *) 'rvec = ', rvec
		write(*, *) 'species ', spec_name(0:nspec)
		write(*, *) '|B| = ', eq%bmag
		write(*, *) 'n(s) = ', eq%ns
		write(*, *) 'T(s) ev = ', eq%ts/e

      	if (test_cold .eqv. .true.) then
			do is = 0, nspec
				call suscep_cold(eq, is, chi)
				label = 'cold suscep '//trim(spec_name(is))
				write(out_unit, '(a," ")')
				call write_matrix(label, chi, 3, 3, out_unit)
			end do

			call dielectric_cold(eq, eps)
			label = 'cold dielectric'
			call write_matrix(label, eps, 3, 3, out_unit)
		end if

      	if (test_warm_bessel .eqv. .true.) then
			do is = 0, nspec
				call suscep_bessel_n1_n3_real(eq, n1_in, n3_in, is, chi)
				label = 'suscep_bessel_n1_n3_real '//trim(spec_name(is))
				write(out_unit, '(a," ")')
				call write_matrix(label, chi, 3, 3, out_unit)
			end do

			spec_model(:) = 'bessel_n1_n3_real'
			write(out_unit,*) 'spec_model = ', trim(spec_model(0))
			call dielectric_general_n1_n3_real(eq, n1_in, n3_in, eps)
			label = 'dielectric_general_n1_n3_real'
			write(out_unit, '(a," ")')
			call write_matrix(label, eps, 3, 3, out_unit)

			n1_c = cmplx(n1_in, zero)
			do is = 0, nspec
				call suscep_bessel_n1_cmplx_n3_real(eq, n1_c, n3_in, is, chi)
				label = 'suscep_bessel_n1_cmplx_n3_real '//trim(spec_name(is))
				write(out_unit, '(a," ")')
				call write_matrix(label, chi, 3, 3, out_unit)
			end do

			spec_model(:) = 'bessel_n1_n3_real'
			write(out_unit,*) 'spec_model(0) = ', trim(spec_model(0))
			call dielectric_general_n1_n3_real(eq, n1_in, n3_in, eps)
			label = 'dielectric_general_n1_n3_real'
			write(out_unit, '(a," ")')
			call write_matrix(label, eps, 3, 3, out_unit)
		end if

! Write dispersion functions

      	if (test_cold .eqv. .true.) then
			label = 'disp_fun_cold_real'
			write(out_unit, '(a," ")')
			write(out_unit,*) label, ' = ',&
			  &   disp_fun_cold_real(eq, n1_in, zero, n3_in)
		end if

      	if (test_warm_bessel .eqv. .true.) then
			spec_model(:) = 'cold'
			write(out_unit, '(a," ")')
			write(out_unit,*) 'spec_model(0) = ', trim(spec_model(0))
			label = 'disp_fun_general_real'
			write(out_unit,*) label, ' = ',&
			  &   disp_fun_general_real(eq, n1_in, n3_in)
		end if

      	if (test_warm_bessel .eqv. .true.) then
			spec_model(:) = 'bessel_n1_n3_real'
			write(out_unit, '(a," ")')
			write(out_unit,*) 'spec_model(0) = ', trim(spec_model(0))
			label = 'disp_fun_general_real'
			write(out_unit,*) label, ' = ',&
			  &   disp_fun_general_real(eq, n1_in, n3_in)
		end if

		write(out_unit, '(a," ")')

	end if

!     if (scan) then
!
! 		if (n_nz_points > 1) then
! 			dnz = (nz_max - nz_min)/(n_nz_points-1)
! 		else
! 			dnz = zero
! 		end if
!
! 		! get x vector
! 		if (n_xpoints > 1) then
! 			dx = (xmax - xmin)/(n_xpoints-1)
! 		else
! 			dx = zero
! 		end if
!
! 		do ix = 0, n_xpoints-1
! 			x = xmin + ix*dx
! 			x_vec(ix+1) = x
! 		end do
!
! 		write(out_unit, *) 'x_vec'
! 		write(out_unit, *) x_vec
!
! 	!*************************************************************************
! 	! nz scan_nz
! 	!*************************************************************************
! 		nz_scan: if (scan_nz) then
! 		! get nz vector
! 			if (n_nz_points > 1) then
! 				dnz = (nz_max - nz_min)/(n_nz_points-1)
! 			else
! 				dnz = 0.
! 			end if
!
! 			do inz = 0, n_xpoints-1
! 				nz = nz_min + ix*dnz
! 				nz_vec(inz+1) = nz
! 			end do
! 			write(out_unit, *) 'nz_vec'
! 			write(out_unit, *) nz_vec
!
! 			x_loop: do ix = 1, n_xpoints
! 				rvec = (/ x_vec(ix), 0.0_rkind,  0.0_rkind /)
! 				call equilibrium(rvec, eq )
!
! 				nz_loop: do inz = 0, n_nz_points-1
! 					nz = nz_min + inz*dx
!
! 					call solve_cold_n1sq_vs_n3(eq, nz, nxsq)
! 					write (*,*) 'nxsq = ', nxsq
! 					write (*,*) 'nsq = ', nxsq + cmplx(nz**2,zero)
! 		! 			nxsq_vec(inz+1) = nxsq
!
! 					do i = 1,4
! 						nxr = real(sqrt(nxsq(i)), KIND=rkind)
! 						if (abs(nxr) > zero) then
! 							theta = atan(nxr/nz)
! 						else
! 							theta = pi/2.
! 						end if
! 						call solve_cold_nsq_vs_theta(eq, theta, nsq)
! 						write (*,*) 'theta = ', theta, 'AH nsq = ', nsq
! 					end do
! 		!			nsq_vec(inz+1) = nsq
! 		! 			write(out_unit,*) 'x = ', x_vec(ix), '  nz = ', nz
! 		! 			write(out_unit,*) 'nxsq = ', nxsq, '  nxsq + nz**2 = ', nxsq+nz**2
! 		! 			write(out_unit,*) 'theta = ', theta, '  nsq = ', nsq
!
! 				end do nz_loop
! 		end do x_loop
!
! 		end if
!
!     end if

!*************************************************************************

	close (unit = out_unit)

    return
 end subroutine dielectric_tensor_test

! !*************************************************************************
!
!     subroutine write_cmplx_matrix(out_unit, mess, value, m_dim, n_dim)
! ! Adapted from cmatrixdbl_message() in diagnostics_m.f90
!
!     implicit none
!     character (len=*), intent (in) :: mess
!     integer, intent (in) :: out_unit, m_dim, n_dim
!     complex(kind=rkind), intent (in) :: value(m_dim, n_dim)
!
!     real(kind=rkind) :: v_min, v_max
!     integer :: i, j
!
! !******************************
!         v_min = huge(v_min)
!         v_max = tiny(v_max)
!
!         do i = 1, m_dim
!             do j = 1, n_dim
!
!                 if ((abs(value(i,j)) > 0.).and.(abs(value(i,j)) < v_min)) &
!                     & v_min = abs(value(i,j))
!                 if (abs(value(i,j)) > v_max ) v_max = abs(value(i,j))
!             end do
!         end do
!
!
!         if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then
!
!             write(out_unit, '(a," ")')
!             write(out_unit, '(a, " = ")' ) trim(mess)
! 			write(*, '(a, " = ")' ) trim(mess)
!
!             do i = 1, m_dim
!                 write(out_unit, '(i3,2x,8(3x,2f14.6))' )  i, (value(i, j), j=1, n_dim)
! 				write(*, '(i3,2x,8(3x,2f14.6))' )  i, (value(i, j), j=1, n_dim)
!             end do
!
!         else
!
!             write(out_unit, '(a," ")')
!             write(out_unit, '(a, " = ")' ) trim(mess)
! 			write(*, '(a, " = ")' ) trim(mess)
!             do i = 1, m_dim
!                 write(out_unit, '(i3,2x,8(3x,2(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
! 				write(*, '(i3,2x,8(3x,2(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
!             end do
!
!         end if
!
! !******************************
!     return
!     end subroutine write_cmplx_matrix

!*************************************************************************


end module dielectric_tensor_test_m