module dispersion_solver_test_m
! Tests various dispersion solvers.  At present this cycles through a series of
! position x points and gets the equilibrium at those points.  Then cycles through
! a series of nz values and evaluates nx(x,nz)
! As of 3/3/2025 the only things evaluated are solve_cold_n1sq_vs_n3 and
! solve_cold_nsq_vs_theta.
! Also the equilibrium position is taken as rvec = (x,0,0)

    use constants_m, only : rkind

    implicit none


! Data for x scan and nz scan
    integer :: n_xpoints ! Number of x points in scan
    real(KIND=rkind) :: xmin, xmax ! min and max of x scan
    integer :: n_nz_points ! Number of points in nz scan
    real(KIND=rkind) :: nz_min, nz_max ! min and max of nz scan
    integer :: n_theta_points ! Number of points in theta scan
    real(KIND=rkind) :: theta_min, theta_max ! min and max of nz scan

! Switches to select specific parameter scans
    logical :: scan_nz = .true.
    logical :: scan_theta = .true.

    namelist /dispersion_solver_test_list/ scan_nz, scan_nz, &
            & n_xpoints, xmin, xmax, n_nz_points, nz_min, nz_max, &
            & n_theta_points, theta_min, theta_max

 !*************************************************************************

 contains

 subroutine init_dispersion_solver_test(read_input)

    use diagnostics_m, only : message_unit, message, text_message, run_label

    implicit none

    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder

    if (read_input .eqv. .true.) then
    ! Read and write input namelist
   		input_unit = get_unit_number()
        open(unit=input_unit, file='component_test_rays.in',action='read', status='old', form='formatted')
        read(input_unit, dispersion_solver_test_list)
        close(unit=input_unit)
        write(message_unit, dispersion_solver_test_list)
    end if

    call text_message('Finished initialize_dispersion_solver_test ')

    return
 end subroutine init_dispersion_solver_test

 !*************************************************************************

  subroutine dispersion_solver_test

    use constants_m, only : rkind, zero, one, pi
    use diagnostics_m, only : message_unit, message, text_message, run_label
    use equilibrium_m, only : equilibrium, eq_point
    use rf_m, only : ray_dispersion_model, wave_mode, k0
    use species_m, only : nspec
    use dispersion_solvers_m, only : solve_n1_vs_n2_n3

    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

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
    real(KIND=rkind) :: x, x_vec(n_xpoints), nz_vec(n_nz_points), theta_vec(n_theta_points)
    complex(KIND=rkind) :: nxsq_vec(n_nz_points), nsq_vec(n_nz_points)
    complex(KIND=rkind) :: nsq_vec_theta(n_theta_points)

    real(KIND=rkind) :: rvec(3)
    real(KIND=rkind) :: dx, dnz, dtheta, nxr, nzr
    integer :: get_unit_number, out_unit
    external get_unit_number
    integer :: i, ix, inz, itheta
    integer :: k_sign = 1 ! Choose positive sign for k_x

 !*************************************************************************

    ! Open file to put write(*,*) data in
	out_unit = get_unit_number()
	open (unit = out_unit, file = 'dispersion_solver_test_m.out', status = 'unknown')

    if (n_nz_points > 1) then
    	dnz = (nz_max - nz_min)/(n_nz_points-1)
    else
    	dnz = zero
    end if

! get x vector
    if (n_xpoints > 1) then
    	dx = (xmax - xmin)/(n_xpoints-1)
    else
    	dx = zero
    end if

	do ix = 0, n_xpoints-1
		x = xmin + ix*dx
		x_vec(ix+1) = x
	end do
	write(out_unit, *) 'x_vec'
	write(out_unit, *) x_vec

!*************************************************************************
! nz scan_nz
!*************************************************************************
	nz_scan: if (scan_nz) then
	! get nz vector
		if (n_nz_points > 1) then
			dnz = (nz_max - nz_min)/(n_nz_points-1)
		else
			dnz = 0.
		end if

		do inz = 0, n_xpoints-1
			nz = nz_min + ix*dnz
			nz_vec(inz+1) = nz
		end do
		write(out_unit, *) 'nz_vec'
		write(out_unit, *) nz_vec

		x_loop: do ix = 1, n_xpoints
			rvec = (/ x_vec(ix), 0.0_rkind,  0.0_rkind /)
			call equilibrium(rvec, eq )

			nz_loop: do inz = 0, n_nz_points-1
				nz = nz_min + inz*dx

				call solve_cold_n1sq_vs_n3(eq, nz, nxsq)
				write (*,*) 'nxsq = ', nxsq
				write (*,*) 'nsq = ', nxsq + cmplx(nz**2,zero)
	! 			nxsq_vec(inz+1) = nxsq

				do i = 1,4
					nxr = real(sqrt(nxsq(i)), KIND=rkind)
					if (abs(nxr) > zero) then
						theta = atan(nxr/nz)
					else
						theta = pi/2.
					end if
					call solve_cold_nsq_vs_theta(eq, theta, nsq)
					write (*,*) 'theta = ', theta, 'AH nsq = ', nsq
				end do
	!			nsq_vec(inz+1) = nsq
	! 			write(out_unit,*) 'x = ', x_vec(ix), '  nz = ', nz
	! 			write(out_unit,*) 'nxsq = ', nxsq, '  nxsq + nz**2 = ', nxsq+nz**2
	! 			write(out_unit,*) 'theta = ', theta, '  nsq = ', nsq

			end do nz_loop
	end do x_loop

	end if nz_scan

!*************************************************************************
! scan theta
!*************************************************************************
	theta_scan: if (scan_theta) then
	! get theta vector
		if (n_theta_points > 1) then
			dtheta = (theta_max - theta_min)/(n_theta_points-1)
		else
			dtheta = zero
		end if

		do itheta = 0, n_xpoints-1
			theta = theta_min + itheta*dtheta
			theta_vec(itheta+1) = theta
		end do
!		write(*, *) 'theta_vec'
!		write(*, *) theta_vec

		x_loop_theta: do ix = 1, n_xpoints
			rvec = (/ x_vec(ix), zero,  zero /)
			call equilibrium(rvec, eq )

			theta_loop: do itheta = 0, n_theta_points-1
				write(*,*) ''
				write(*,*) ''
				write(*,*) ''
				theta = theta_min + itheta*dtheta
				call solve_cold_nsq_vs_theta(eq, theta, nsq)
				write (*,*) 'theta = ', theta,'  nsq = ', nsq
				n1 = sqrt(cmplx(nsq, zero))
				write (*,*) 'n1 = ', n1, '   residual = ', disp_fun_cold_n_theta(eq, theta, n1)

				do i = 1,4
					! Comparison to nx(nz) meaningful only for n(theta) real
					if (abs(aimag(n1(i))) < tiny(one)) then
						n1z(i) = n1(i)*cos(theta)
						n1x(i) = n1(i)*sin(theta)

						write(*,*) ''
						write (*,*) 'n1z(i) = ', n1z(i)
						write (*,*) 'n1x(i) = ', n1x(i)
						call solve_cold_n1sq_vs_n3(eq, n1z(i), nxsq)
						write (*,*) 'nxsq(i) = ', nxsq(i)
						nx1(i) = Sqrt(nxsq(i))
						write (*,*) 'nx1(i) = ', nx1(i)
						nx1r = real(nx1(i), KIND=rkind)
						write (*,*) 'nx1r = ', nx1r(i), '   residual = ', &
						       & disp_fun_cold_n1_n3(eq, n1z(i), nx1(i))
						write (*,*) 'cross residual = ', &
						       & disp_fun_cold_n1_n3(eq, n1z(i), n1x(i))
					end if
				end do

			end do theta_loop
	end do x_loop_theta

	end if theta_scan
!*************************************************************************

	close (unit = out_unit)
    call text_message('Finished dispersion_solver_test work')

    return
 end subroutine dispersion_solver_test


!*************************************************************************

end module dispersion_solver_test_m