module ray_deriv_test_m
! Tests various deriv routines.
! There is one kind of test:
! 1) one_point -> Evaluates ray derivatives at one equilibrium point for
!    one pair of (nx, nz) for one or more deriv routines
!
! deriv routines supported are cold, general, and num
!
! N.B. For convenience refractive indices are input rather than k values although
!      deriv_num uses k from ray vector v.

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

! Switches to select deriv models. Defaults can be reset on input
    logical :: test_cold = .true.
    logical :: test_general = .true.
    logical :: test_num = .true.

! Switches to select specific parameter scans
    logical :: scan_nx = .false.

! Data for one_point test
    integer :: number_of__models ! As of now 1 or 2
    real(KIND=rkind) :: x_in, y_in, z_in ! Eq point to evaluate
    real(KIND=rkind) :: n1_in, n3_in ! real part of n perp and n|| to evaluate
    real(KIND=rkind) :: n1_im = zero, n3_im = zero ! imaginary part of n perp and n||
    complex(KIND=rkind) :: n1_c, n3_c ! complex n perp and n|| to evaluate

! Data for x scan and nx scan
    integer :: n_xpoints ! Number of x points in scan
    real(KIND=rkind) :: xmin, xmax ! min and max of x scan
    integer :: n_nx_points ! Number of points in nx scan
    real(KIND=rkind) :: nx_min, nx_max ! min and max of nx scan

    namelist /ray_deriv_test_list/ &
            & test_cold, test_general, test_num, scan_nx, &
            & x_in, y_in, z_in, n1_in, n3_in, &
            & n_xpoints, xmin, xmax, n_nx_points, nx_min, nx_max

 !*************************************************************************

 contains

 subroutine init_ray_deriv_test(read_input)

    use diagnostics_m, only : message_unit, message, text_message, run_label

    implicit none

    logical, intent(in) :: read_input
    integer :: input_unit, get_unit_number ! External, free unit finder

    if (read_input .eqv. .true.) then
    ! Read and write input namelist
        input_unit = get_unit_number()
        open(unit=input_unit, file='component_test_rays.in',action='read', status='old',&
           & form='formatted')
        read(input_unit, ray_deriv_test_list)
        close(unit=input_unit)
        write(message_unit, ray_deriv_test_list)
    end if

    call text_message('Finished initialize_ray_deriv_test ')

    return
 end subroutine init_ray_deriv_test

 !*************************************************************************

  subroutine ray_deriv_test

    use constants_m, only : e
    use diagnostics_m, only : message, text_message, messages_to_stdout
    use equilibrium_m, only : equilibrium, eq_point
    use ode_m, only : nv, ray_deriv_name
    use rf_m, only : k0, frf
    use species_m, only : nspec, spec_name, spec_model, nmins, nmaxs
    use suscep_m, only : suscep_cold, dielectric_cold, &
              & suscep_bessel, dielectric_general, &
              & disp_fun_cold, disp_fun_general, disp_fun_general_Hermitian

    implicit none

    real(KIND=rkind) :: rvec(3), nvec(3), kvec(3)
    integer :: get_unit_number, out_unit
    external get_unit_number
    character(len=80) :: label
    real(KIND=rkind) :: v(nv)

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

! Results
    real(KIND=rkind) :: x, x_vec(n_xpoints), nx_vec(n_nx_points)
    real(KIND=rkind) :: dx, dnx, nx, nz
    integer :: i, ix, inx
    real(KIND=rkind):: dddx(3), dddk(3), dddw

 !*************************************************************************

    ! Open file to put write(*,*) data in
!   out_unit = message_unit
    out_unit = get_unit_number()
    open (unit = out_unit, file = 'ray_deriv_test_m.out', status = 'unknown')

    rvec = (/x_in, y_in, z_in/)
    call equilibrium(rvec, eq)

    nvec = (/n1_in, zero, n3_in/)
    kvec = k0*nvec

    v = zero
    v(1:3) = rvec
    v(4:6) = kvec

    write(out_unit, *) ' '
    write(out_unit, *) 'rvec = ', rvec
    write(out_unit, *) 'frf = ', frf
    write(out_unit, *) 'n1_in = ', n1_in, '  n3_in = ', n3_in
    write(out_unit, *) 'species ', spec_name(0:nspec)
    write(out_unit, *) '|B| = ', eq%bmag
    write(out_unit, *) 'n(s) = ', eq%ns
    write(out_unit, *) 'T(s) ev = ', eq%ts/e
    write(out_unit, *) 'nmin(s) = ', nmins
    write(out_unit, *) 'nmax(s) = ', nmaxs
    write(*, *) 'rvec = ', rvec
    write(*, *) 'species ', spec_name(0:nspec)
    write(*, *) '|B| = ', eq%bmag
    write(*, *) 'n(s) = ', eq%ns
    write(*, *) 'T(s) ev = ', eq%ts/e
    write(*, *) 'nmin(s) = ', nmins
    write(*, *) 'nmax(s) = ', nmaxs

    write(*, *) ' '
    write(*, *) 'alpha(0) = ', eq%alpha(0)
    write(*, *) 'gamma(0) = ', eq%gamma(0)
    write(out_unit, *) ' '
    write(out_unit, *) 'alpha(0) = ', eq%alpha(0)
    write(out_unit, *) 'gamma(0) = ', eq%gamma(0)

    if (test_cold .eqv. .true.) then
        ray_deriv_name = 'cold'
        call deriv_cold(eq, v, dddx, dddk, dddw)
        label = 'deriv_cold '
        write(*,*) ''
        write(*, '(a)') label
        write(*,'(a,1p3e12.4)') 'dddx =', dddx
        write(*,'(a,1p3e12.4)') 'dddk =', dddk
        write(*,'(a,1p1e12.4)') 'dddw =', dddw
        write(*,*) ''
        write(out_unit,*) ''
        write(out_unit, '(a)') label
        write(out_unit,'(a,1p3e12.4)') 'dddx =', dddx
        write(out_unit,'(a,1p3e12.4)') 'dddk =', dddk
        write(out_unit,'(a,1p1e12.4)') 'dddw =', dddw
        write(out_unit,*) ''
    end if

    if (test_general .eqv. .true.) then
        ray_deriv_name = 'general'
        call deriv_general(eq, v, dddx, dddk, dddw)
        label = 'deriv_general'
        write(*,*) ''
        write(*, '(a)') label
        write(*,'(a,1p3e12.4)') 'dddx =', dddx
        write(*,'(a,1p3e12.4)') 'dddk =', dddk
        write(*,'(a,1p1e12.4)') 'dddw =', dddw
        write(*,*) ''
        write(out_unit,*) ''
        write(out_unit, '(a)') label
        write(out_unit,'(a,1p3e12.4)') 'dddx =', dddx
        write(out_unit,'(a,1p3e12.4)') 'dddk =', dddk
        write(out_unit,'(a,1p1e12.4)') 'dddw =', dddw
        write(out_unit,*) ''
    end if

    if (test_num .eqv. .true.) then
    !     Numerical differentiation.
    !     N.B. Must be called with v(), not just nvec.  Evaluates eq at other positions
    !     so v(1:3) is needed

        ray_deriv_name = 'general'
        call deriv_num(eq, v, dddx, dddk, dddw)
        label = 'deriv_num '
        write(*,*) ''
        write(*, '(a)') label
        write(*,'(a,1p3e12.4)') 'dddx =', dddx
        write(*,'(a,1p3e12.4)') 'dddk =', dddk
        write(*,'(a,1p1e12.4)') 'dddw =', dddw
        write(*,*) ''
        write(out_unit,*) ''
        write(out_unit, '(a)') label
        write(out_unit,'(a,1p3e12.4)') 'dddx =', dddx
        write(out_unit,'(a,1p3e12.4)') 'dddk =', dddk
        write(out_unit,'(a,1p1e12.4)') 'dddw =', dddw
        write(out_unit,*) ''
    end if

!*************************************************************************
! nx, x scan_nx
!*************************************************************************
	nx_scan: if (scan_nx) then
	! generate nx vector
		if (n_nx_points > 1) then
			dnx = (nx_max - nx_min)/(n_nx_points-1)
		else
			dnx = 0.
		end if

		do inx = 0, n_xpoints-1
			nx = nx_min + ix*dnx
			nx_vec(inx+1) = nx
		end do
		write(out_unit, *) 'nx_vec'
		write(out_unit, *) nx_vec

		x_loop: do ix = 1, n_xpoints
			rvec = (/ x_vec(ix), 0.0_rkind,  0.0_rkind /)
			call equilibrium(rvec, eq )

			nx_loop: do inx = 0, n_nx_points-1
				nx = nx_min + inx*dx

				if (test_cold .eqv. .true.) then
					ray_deriv_name = 'cold'
					call deriv_cold(eq, v, dddx, dddk, dddw)
					label = 'deriv_cold '
					write(out_unit,*) ''
					write(out_unit, '(a)') label
					write(out_unit,'(a,1p3e12.4)') 'dddx =', dddx
					write(out_unit,'(a,1p3e12.4)') 'dddk =', dddk
					write(out_unit,'(a,1p1e12.4)') 'dddw =', dddw
					write(out_unit,*) ''
				end if

				if (test_general .eqv. .true.) then
					ray_deriv_name = 'general'
					call deriv_general(eq, v, dddx, dddk, dddw)
					label = 'deriv_general'
					write(out_unit,*) ''
					write(out_unit, '(a)') label
					write(out_unit,'(a,1p3e12.4)') 'dddx =', dddx
					write(out_unit,'(a,1p3e12.4)') 'dddk =', dddk
					write(out_unit,'(a,1p1e12.4)') 'dddw =', dddw
					write(out_unit,*) ''
				end if

				if (test_num .eqv. .true.) then
				!     Numerical differentiation.
				!     N.B. Must be called with v(), not just nvec.  Evaluates eq at
				!          other positions so v(1:3) is needed

					ray_deriv_name = 'general'
					call deriv_num(eq, v, dddx, dddk, dddw)
					label = 'deriv_num '
					write(out_unit,*) ''
					write(out_unit, '(a)') label
					write(out_unit,'(a,1p3e12.4)') 'dddx =', dddx
					write(out_unit,'(a,1p3e12.4)') 'dddk =', dddk
					write(out_unit,'(a,1p1e12.4)') 'dddw =', dddw
					write(out_unit,*) ''
				end if

			end do nx_loop
		end do x_loop

	end if nx_scan

!*************************************************************************

    close (unit = out_unit)

    return
 end subroutine ray_deriv_test

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
!           write(*, '(a, " = ")' ) trim(mess)
!
!             do i = 1, m_dim
!                 write(out_unit, '(i3,2x,8(3x,2f14.6))' )  i, (value(i, j), j=1, n_dim)
!               write(*, '(i3,2x,8(3x,2f14.6))' )  i, (value(i, j), j=1, n_dim)
!             end do
!
!         else
!
!             write(out_unit, '(a," ")')
!             write(out_unit, '(a, " = ")' ) trim(mess)
!           write(*, '(a, " = ")' ) trim(mess)
!             do i = 1, m_dim
!                 write(out_unit, '(i3,2x,8(3x,2(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
!               write(*, '(i3,2x,8(3x,2(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
!             end do
!
!         end if
!
! !******************************
!     return
!     end subroutine write_cmplx_matrix

!*************************************************************************


end module ray_deriv_test_m