		call restart_X_mode(OX_conv_data(i)%x_cut, OX_conv_data(i)%ny_c, &
		   & OX_conv_data(i)%nz_c, OX_conv_data(i)%x_restart, OX_conv_data(i)%k_restart)

!****************************************************************************

 subroutine restart_X_mode(rvecO_cutoff, ny_c, nz_c, rvecX_restart, nvec_Xrestart)

! Calculate the initial conditions for a the x-mode on high density side of the
! conversion layer.
!
! Subject to the slab geometry approximation inherent in the conversion model we assume
! that the components of n perpendicular to grad(n) (i.e ny_c and nz_c) are
! constant along the path from the ray turning point to the X-mode propagation location
! (these were also assumed the same at the O-mode cutoff surface to calculate the
! conversion coefficient).  For more discussion see notes "O-X conversion",11/21/2025.
!
! First we find the density for the X-mode cutoff by solving for the root of the nx**0
! coefficient of the Booker quartic at or above the alpha = 1 surface.

! First find the point on the high density side where the X-mode
! starts to propagate. We take the restart location to be the in the direction of grad(n)
! from the O-mode cutoff position, x_cut, found above.  We search a distane, s, along the
! grad(n) vector until we find a root of the dispersion function. That will be the X-mode
! cutoff.  Then extend it a little further so the X-mode is propagating.

    use constants_m, only : rkind, zero, one, two
    use equilibrium_m, only : equilibrium, eq_point
    use rf_m, only : k0
    use vectors3_m, only : cross_product

	implicit none

	real(KIND=rkind), intent(in) :: rvecO_cutoff(3)
	real(KIND=rkind), intent(in) :: ny_c
	real(KIND=rkind), intent(in) :: nvz_c
	real(KIND=rkind), intent(out) :: rvecX_restart(3)
	real(KIND=rkind), intent(out) :: nvecX_restart(3)

	real(KIND=rkind), intent(out) :: h_vec(3)
	real(KIND=rkind) :: xc_unit(3), yc_unit(3), zc_unit(3)
	real(KIND=rkind) :: v_temp(3)
	real(KIND=rkind) :: n_parallel, ny_c, nz_c, n_vertical, n_crit

! Distance along grad(n) for X-mode cutoff
	real(KIND=rkind) :: s_cut

! Maximum distance along grad(n) to search for X-mode cutoff with bisection
    real(KIND=rkind), parameter :: s_max = 0.02_rkind

! Tolerance for finding s, distance along grad(n) for X-mode cutoff
    real(KIND=rkind), parameter :: bisection_eps = one*10d-6

! Distance along grad(n) to extend the restart point past the X-mode cutoff
	real(KIND=rkind), parameter :: s_extend = 0.001_rkind ! one millimeter

! Additional data to feed to solve_bisection
	real(KIND=rkind) :: data(4)

	type(eq_point) :: eq

! N.B. These things are evaluated at the O-mode cutoff surface
	call equilibrium(rvecO_cutoff, eq)
	xc_unit = eq%gradns(:,0)/norm2(eq%gradns(:, 0)) ! Unit vector along grad(ne)
	v_temp = cross_product(eq%bunit, xc_unit) ! perpendicular to x and B, i.e. y direction
	theta = acos(dot_product(xc_unit, eq%bunit))
	gamma = abs(eq%gamma(0))

! Find density alphaX = (omega plasma/omega rf)**2 of X-mode cutoff (i.e. root of Booker
! quartic)
	data(1) = gamma
	data(2) = ny_c
	data(3) = nz_c
	data(4) = theta
    call solve_bisection(f_Booker0, alphaX, one, two, zero,&
       & bisection_eps, ierr, data)


	call equilibrium(x_cutoff_ray, eq)
	xc_unit = eq%gradns(:,0)/norm2(eq%gradns(:, 0)) ! Unit vector along grad(ne)

	h_vec = x_cutoff_ray - x_max_ray
	norm_h_vec = norm2(h_vec)
	k_restart = k_max_ray

! N.B. These things are evaluated at the cutoff surface
	call equilibrium(x_cutoff_ray, eq)
	xc_unit = eq%gradns(:,0)/norm2(eq%gradns(:, 0)) ! Unit vector along grad(ne)
	v_temp = cross_product(eq%bunit, xc_unit) ! perpendicular to x and B, i.e. y direction
	yc_unit = v_temp/norm2(v_temp)
	zc_unit = cross_product(xc_unit, yc_unit)
	theta = acos(dot_product(xc_unit, eq%bunit))
	gamma = abs(eq%gamma(0))

 write(*,*) 'x_cutoff_ray = ', x_cutoff_ray
 write(*,*) 'x_restart = ', x_restart
 write(*,*) 'k_max_ray = ', k_max_ray
 write(*,*) 'k_restart = ', k_restart

	return
 end subroutine restart_X_mode


 function f_Booker0(X, data) result(d)
! Evaluates the nx**0 coefficient of the Booker quartic versus X=(plasma freq/rf freq)**2
! This is the magnetoionic notation.  Gives the density of cutoff. Components
! of n_vec perpendicular to grad(in) are taken to be constant, and the same as at the
! O-mode turning point, and also at the O-mode cutoff.

	implicit none
	real(KIND=rkind) :: X
	real(KIND=rkind) :: data(4)
	real(KIND=rkind) :: d


! Local coordinate unit vectors at location s
	real(KIND=rkind) :: xc_unit(3), yc_unit(3), zc_unit(3)
	real(KIND=rkind) :: v_temp(3)
	real(KIND=rkind) :: n_parallel, ny_c, nz_c, n_vertical, n_crit

	Y = data(1)
	ny = data(2)
	nz = data(3)
	theta = data(4)

	d = &
     & ((1 - X)*((1 - X)**2 - Y**2))/(1 - Y**2) + ny**4*(1 - X/(1 - Y**2)) - &
     & ny**2*(-((X*Y**2)/(1 - Y**2)) + 2*(1 - X)*(1 - X/(1 - Y**2))) + &
     & ny**2*nz**2*(2*(1 - X/(1 - Y**2)) + (X*Y**2*Cos(theta)**2)/(1 - Y**2)) + &
     & nz**4*(1 - (X*(1 - Y**2*Cos(theta)**2))/(1 - Y**2)) - &
     & nz**2*(2*(1 - X)*(1 - X/(1 - Y**2)) - (X*Y**2*Sin(theta)**2)/(1 - Y**2))

    return
 end function f_Booker0