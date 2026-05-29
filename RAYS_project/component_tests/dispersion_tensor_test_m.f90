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
! N.B. For convenience refractive indices are input.  These are converted to k values
!      and put into module

    use constants_m, only : rkind, zero, one, pi

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

! Data for x scan and nz scan
    integer :: n_xpoints ! Number of x points in scan
    real(KIND=rkind) :: xmin, xmax ! min and max of x scan
    integer :: n_nz_points ! Number of points in nz scan
    real(KIND=rkind) :: nz_min, nz_max ! min and max of nz scan

! Switches to select specific parameter scans
    logical :: scan_x = .true.
    logical :: scan_nz = .true.

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
        open(unit=input_unit, file='component_test_rays.in',action='read', status='old', form='formatted')
        read(input_unit, dielectric_tensor_test_list)
        close(unit=input_unit)
        write(message_unit, dielectric_tensor_test_list)
    end if

    call text_message('Finished initialize_dielectric_tensor_test ')

    return
 end subroutine init_dielectric_tensor_test

 !*************************************************************************

  subroutine dielectric_tensor_test

    use diagnostics_m, only : message_unit, message, text_message, run_label
    use equilibrium_m, only : equilibrium, eq_point
    use rf_m, only : ray_dispersion_model, wave_mode, k0
    use species_m, only : nspec, spec_name
    use suscep_m, only : dielectric_cold

    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

    complex(KIND=rkind) :: eps(3,3)

! Args for solve_cold_n1sq_vs_n3
   complex(KIND=rkind) :: nxsq(4), nx1(4)
   real(KIND=rkind) :: nx1r(4)
   real(KIND=rkind) :: nz
   real(KIND=rkind) :: disp_fun_cold_n1_n3

! Args for solve_cold_nsq_vs_theta
   real(KIND=rkind) :: nsq(4)
   complex(KIND=rkind) :: n1(4)
   real(KIND=rkind) :: n1z(4)
   complex(KIND=rkind) :: n1x(4)
   real(KIND=rkind) :: theta
   real(KIND=rkind) :: disp_fun_cold_n_theta

! Results
    real(KIND=rkind) :: x, x_vec(n_xpoints), nz_vec(n_nz_points)
    complex(KIND=rkind) :: nxsq_vec(n_nz_points), nsq_vec(n_nz_points)

    real(KIND=rkind) :: rvec(3)
    real(KIND=rkind) :: dx, dnz, dtheta, nxr, nzr
    integer :: get_unit_number, out_unit
    external get_unit_number
    integer :: i, ix, inz, itheta
    integer :: k_sign = 1 ! Choose positive sign for k_x

 !*************************************************************************

    ! Open file to put write(*,*) data in
	out_unit = get_unit_number()
	open (unit = out_unit, file = 'dielectric_tensor_test_m.out', status = 'unknown')


    if (one_point) then

      	rvec = (/x_in, y_in, z_in/)
      	call equilibrium(rvec, eq)

		write(out_unit, *) ' '
		write(out_unit, *) 'cold plasma: rvec = ', rvec, '  species ', spec_name(0:nspec)
		write(out_unit, *) '|B| = ', eq%bmag


      	if (test_cold .eqv. .true.) then
			call dielectric_cold(eq, eps)
		end if

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
    call text_message('Finished dielectric_tensor_test work')

    return
 end subroutine dielectric_tensor_test


!*************************************************************************

end module dielectric_tensor_test_m