module B_loop_m

! Routines to calculate magnetic field from a circular loop of radius a, carrying
! current of 1 Amp, centered at (x,y,z) = 0, lying in the z = 0 plane

!   use constants_m, only : rkind, pi, mu0, zero, one

    use complete_elliptic_int_m, only : elliptic_Em, elliptic_Km

    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind),parameter :: zero = 0.0_rkind
    real(KIND=rkind),parameter :: one = 1.0_rkind
    real(KIND=rkind),parameter :: n2 = 2.0_rkind, n3 = 3.0_rkind, n4 = 4.0_rkind
    real(KIND=rkind),parameter :: n5 = 5.0_rkind, n6 = 6.0_rkind, n8 = 8.0_rkind
    real(KIND=rkind),parameter :: n10 = 10.0_rkind,  n12 = 12.0_rkind
    real(KIND=rkind),parameter :: n14 = 14.0_rkind, n15 = 15.0_rkind, n16 = 16.0_rkind
    real(KIND=rkind),parameter :: n32 = 32.0_rkind, n45 = 45.0_rkind, n128 = 128.0_rkind
    real(KIND=rkind),parameter :: n256 = 256.0_rkind, n484 = 484.0_rkind

    real(KIND=rkind),parameter :: r0 = 0.001_rkind ! Threshold for small arg expansion

    real(KIND=rkind),parameter :: pi = 3.1415926535897932385

    real(KIND=rkind),parameter :: mu0 = pi * 4.e-7
!    real(KIND=rkind),parameter :: mu0 = 1.0
    real(KIND=rkind),parameter :: c0 = mu0/(n2*pi)

!******************************

contains

!******************************

subroutine Brz_loop(r, z, a, Br, Bz, Aphi)

    implicit none

    real(KIND=rkind), intent(in) :: r, z, a
    real(KIND=rkind), intent(out) :: Br, Bz, Aphi
    real(KIND=rkind) :: m0, m, alpha, beta, gamma, f
    real(KIND=rkind) :: a2, a3, a4, a5, a6, r2, r3, r4, z2, z4
    real(KIND=rkind) :: Em, Km

    r2 = r**2
    a2 = a**2
    z2 = z**2
    m0=(a+r)**2 + z2

! On axis
    if (r < n2*tiny(one)) then
        Br = zero
        Bz = mu0*a2/n2/(a2+z2)**1.5_rkind
        Aphi = zero
        return
    end if

! Near axis
    if (r < r0) then
        a3 = a*a2
        a4 = a2*a2
        a5 = a*a4
        a6 = a3*a3
        z4 = z2*z2
        r3 = r*r2
        r4 = r2*r2
        f = a2+z2
        Br = n3*a2*z*r/(n4*f**2.5_rkind)
        Br = Br - n15*z*r3*(-n3*a4+n4*a2*z2)/(n32*f**4.5_rkind)
        Br = mu0*Br
        Bz = a2/n2/f**1.5_rkind + n3/n8*(a4-n4*a2*z2)*r2/f**3.5_rkind + &
           & n45/n128*(a6-n12*a4*z2+n8*a2*z4)*r4/f**5.5_rkind
        Bz = mu0*Bz
        Aphi = (a2*r2)/(n4*pi*f**1.5_rkind)+ &
           & n3*(a4-n4*a2*z2)*r4/(n32*pi*f**3.5_rkind)
        return
    end if

! Otherwise
    m = n4*a*r/m0
    alpha = a2 + r2 + z2
    beta = a2 - r2 -z2
    gamma = (a-r)**2 + z2
    Em = elliptic_Em(m)
    Km = elliptic_Km(m)
    Br = c0*z/(r*sqrt(m0))*(alpha/gamma*Em - Km)
    Bz = c0/sqrt(m0)*(beta/gamma*Em + Km)
    Aphi = n4*c0*a*r/(m*sqrt(m0))*((one-m/n2)*Km-Em)
    return
end subroutine Brz_loop

!******************************

subroutine Aphi_rz_loop(r, z, a, Aphi)

    implicit none

    real(KIND=rkind), intent(in) :: r, z, a
    real(KIND=rkind), intent(out) :: Aphi
    real(KIND=rkind) :: m0, m, alpha, beta, gamma, f
    real(KIND=rkind) :: a2, a3, a4, a5, a6, r2, r3, r4, z2, z4
    real(KIND=rkind) :: Em, Km

    r2 = r**2
    a2 = a**2
    z2 = z*z
    m0=(a+r)**2 + z2

! On axis
    if (r < n2*tiny(one)) then
        Aphi = zero
        return
    end if

! Near axis
    if (r < r0) then
        a4 = a2*a2
        r3 = r*r2
        f = a2+z2
        Aphi = (a2*r2)/(n4*pi*f**1.5_rkind)+ &
           & n3*(a4-n4*a2*z2)*r4/(n32*pi*f**3.5_rkind)
        return
    end if

! Otherwise
    m = n4*a*r/m0
    alpha = a2 + r2 + z2
    beta = a2 - r2 -z2
    gamma = (a-r)**2 + z2
    Em = elliptic_Em(m)
    Km = elliptic_Km(m)
    Aphi = n4*c0*a*r/(m*sqrt(m0))*((one-m/n2)*Km-Em)

    return
end subroutine Aphi_rz_loop

!******************************

subroutine Bxyz_loop(x, y, z, a, Bx, By, Bz)

    implicit none

    real(KIND=rkind), intent(in) :: x, y, z, a
    real(KIND=rkind), intent(out) :: Bx, By, Bz
    real(KIND=rkind) :: r, r2, Br, a2, z2
    real(KIND=rkind) :: Aphi

    r2 = x**2+y**2
    r = sqrt(r2)

! On axis
    if (r < n2*tiny(one)) then
		Bx = zero
		By = zero
		a2 = a*a
		z2 = z*z
        Bz = mu0*a2/n2/(a2+z2)**1.5_rkind
        return
    end if

	call Brz_loop(r, z, a, Br, Bz, Aphi)
	Bx = x*Br/r
	By = y*Br/r

    return
end subroutine Bxyz_loop

!******************************

subroutine Aphi_xyz_loop(x, y, z, a, Aphi)

    implicit none

    real(KIND=rkind), intent(in) :: x, y, z, a
    real(KIND=rkind), intent(out) :: Aphi
    real(KIND=rkind) :: r, r2

    r2 = x**2+y**2
    r = sqrt(r2)
	call Aphi_rz_loop(r, z, a, Aphi)
    return
end subroutine Aphi_xyz_loop

!******************************

subroutine Brz_loop_scaled(r, z, Br, Bz, Aphi)

    implicit none

    real(KIND=rkind), intent(in) :: r, z
    real(KIND=rkind), intent(out) :: Br, Bz, Aphi
    real(KIND=rkind) :: m0, m, alpha, beta, gamma, f
    real(KIND=rkind) :: r2, r3, r4, z2, z4
    real(KIND=rkind) :: Em, Km

    r2 = r**2
    z2 = z**2

! On axis
    if (r < n2*tiny(one)) then
        Br = zero
        Bz = mu0/(n2*(one + z2)**1.5_rkind)
        Aphi = zero
       return
    end if

! Near axis
    if (r < r0) then
        r3 = r*r2
        r4 = r2*r2
        z4 = z2*z2
        f = one+z2
        Br = n3*z*r/(n4*f**2.5_rkind)
        Br = Br - n15*z*r3*(-n3+n4*z2)/(n32*f**4.5_rkind)
        Br = mu0*Br

        Bz = one/n2/f**1.5_rkind + n3/n8*(one-n4*z2)*r2/f**3.5_rkind + &
           & n45/n128*(one-n12*z2+n8*z4)*r4/f**5.5_rkind
        Bz = mu0*Bz
        Aphi = r2/(n4*pi*f**1.5_rkind)+ &
           & n3*(one-n4*z2)*r4/(n32*pi*f**3.5_rkind)
!            write(*,*) 'r = ', r, 'Aphi = ', Aphi
!            stop
       return
    end if

! Otherwise
    m0=(one + r)**2 + z2
    m = n4*r/m0
    alpha = one + r2 + z2
    beta = one - r2 -z2
    gamma = (one-r)**2 + z2
    Em = elliptic_Em(m)
    Km = elliptic_Km(m)
    Br = c0*z/(r*sqrt(m0))*(alpha/gamma*Em - Km)
    Bz = c0/sqrt(m0)*(beta/gamma*Em + Km)
    Aphi = -sqrt(m0)*Em + alpha/sqrt(m0)*Km
    Aphi = c0*Aphi

    return
end subroutine Brz_loop_scaled

!******************************

subroutine Aphi_rz_loop_scaled(r, z, Aphi)

    implicit none

    real(KIND=rkind), intent(in) :: r, z
    real(KIND=rkind), intent(out) :: Aphi
    real(KIND=rkind) :: m0, m, alpha, beta, gamma, f
    real(KIND=rkind) :: r2, r3, r4, z2, z4
    real(KIND=rkind) :: Em, Km

    r2 = r**2
    z2 = z**2

! On axis
    if (r < n2*tiny(one)) then
        Aphi = zero
        return
    end if

! Near axis
    if (r < r0) then
        r4 = r2*r2
        f = one+z2
        Aphi = r2/(n4*pi*f**1.5_rkind)+ &
           & n3*(one-n4*z2)*r4/(n32*pi*f**3.5_rkind)
        return
    end if

! Otherwise
    m0=(one + r)**2 + z2
    m = n4*r/m0
    alpha = one + r2 + z2
    beta = one - r2 -z2
    gamma = (one-r)**2 + z2
    Em = elliptic_Em(m)
    Km = elliptic_Km(m)
    Aphi = -sqrt(m0)*Em + alpha/sqrt(m0)*Km
    Aphi = c0*Aphi

    return
end subroutine Aphi_rz_loop_scaled

!******************************

subroutine Bxyz_loop_scaled(x, y, z, Bx, By, Bz)

    implicit none

    real(KIND=rkind), intent(in) :: x, y, z
    real(KIND=rkind), intent(out) :: Bx, By, Bz
    real(KIND=rkind) :: r, r2, z2, Br
    real(KIND=rkind) :: Aphi

    r2 = x**2+y**2
    r = sqrt(r2)

! On axis
    if (r < n2*tiny(one)) then
		Bx = zero
		By = zero
		z2 = z*z
        Bz = mu0/(n2*(one + z2)**1.5_rkind)
        return
    end if

	call Brz_loop_scaled(r, z, Br, Bz, Aphi)
	Bx = x*Br/r
	By = y*Br/r

    return
end subroutine Bxyz_loop_scaled

!******************************

subroutine Aphi_xyz_loop_scaled(x, y, z, Aphi)

    implicit none

    real(KIND=rkind), intent(in) :: x, y, z
    real(KIND=rkind), intent(out) :: Aphi
    real(KIND=rkind) :: r, r2, Br

    r2 = x**2+y**2
    r = sqrt(r2)
	call Aphi_rz_loop_scaled(r, z, Aphi)

    return
end subroutine Aphi_xyz_loop_scaled


end module B_loop_m
