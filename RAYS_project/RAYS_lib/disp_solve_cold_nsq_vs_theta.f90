 subroutine solve_cold_nsq_vs_theta(eq, theta, nsq)
! calculates cold plasma refractive index squared versus theta = angle between n and B.
! In cold plasma nsq is real, but may be negative
!
! The dispersion relation is quadratic in n^2.  The roots are classified according to
! whether the sign in front of the discriminant is + or -.  Roots are also classified as
! to which |n| is larger -> slow mode, or smaller -> fast mode.  This returns both choices:
! nsq(1) -> plus
! nsq(2) -> minus
! nsq(3) -> fast
! nsq(4) -> slow

    use constants_m, only : rkind, one
    use diagnostics_m, only : message
    use equilibrium_m, only : eq_point
!    use suscep_m, only :  dielectric_cold
    use suscep_m, only :  RLSDP_cold

       implicit none

!      Derived type containing equilibrium data for a spatial point in the plasma
       type(eq_point), intent(in) :: eq

       real(KIND=rkind), intent(in) :: theta
!       complex(KIND=rkind), intent(out) :: nsq(4)
       real(KIND=rkind), intent(out) :: nsq(4)

       real(KIND=rkind) :: S ,D , P,  R, L
       real(KIND=rkind) :: sin2, cos2, a, b, c, sgn_b
       real(KIND=rkind) :: discr, sqrt_discr
       real, parameter :: two = 2.0_rkind, four = 4.0_rkind

       call RLSDP_cold(eq, S ,D , P, R, L)

	   cos2 = cos(theta)**2
	   sin2 = one - cos2

!      Coefficients for A(n3)*(nsq)^2 + B(n3)*nsq + C(n3) = 0.
       a = S*sin2 + P*cos2
       b = -R*L*sin2 - P*S*(one + cos2)
       c = P*R*L
       sgn_b = sign(one,b)   ! SIGN(A,B)=|A| if B>=0 and SIGN(A,B)=-|A| if B<0.
       discr = b**2-four*a*c

       if ( discr < 0. ) then
          call message('disp_solve_cold_nsq_vs_theta:  evanescent root', 0)
          call message('disp_solve_cold_nsq_vs_theta: b**2-4.*a*c', discr, 0)
          write(*,*) 'disp_solve_cold_nsq_vs_theta: b**2-4.*a*c', discr
          return ! This should never happen since in cold plasma discr is real and non-negative.
       end if
       sqrt_discr = sqrt(discr)
       if (sgn_b < 0.) then
           nsq(1) = ( -b + sqrt_discr) / (two*a) ! plus root
           nsq(2) = two*c/(-b + sqrt_discr)      ! minus root
       else
           nsq(2) = ( -b - sqrt_discr) / (two*a) ! minus root
           nsq(1) = two*c/(-b - sqrt_discr)      ! plus root
       end if

!   write(*,*) 'a = ', a, ' b = ', b, ' c = ', c
!   write(*,*) 'nsq = ', nsq, '    sqrt_d = ', sqrt_d

!      Store smaller root (fast wave) into nsq(3) and larger root (slow wave) into nsq(4).
       if ( abs(nsq(1)) <= abs(nsq(2)) ) then
           nsq(3) = nsq(1) ! fast root
           nsq(4) = nsq(2) ! slow root
       else
          nsq(3) = nsq(2) ! fast root
          nsq(4) = nsq(1) ! slow root
       end if

       return
 end subroutine solve_cold_nsq_vs_theta
!****************************************************************************

complex function disp_fun_cold_n_theta(eq, theta, n)
! calculates the residual of cold plasma dispersion relation versus refractive index, n,
! and angle theta
! N.B. n is complex, theta is real

    use constants_m, only : rkind, one
    use diagnostics_m, only : message
    use equilibrium_m, only : eq_point
    use suscep_m, only :  RLSDP_cold

       implicit none

!      Derived type containing equilibrium data for a spatial point in the plasma
       type(eq_point), intent(in) :: eq

       complex(KIND=rkind), intent(in) :: n
       real(KIND=rkind), intent(in) :: theta

       real(KIND=rkind) :: S ,D , P,  R, L
       real(KIND=rkind) :: sin2, cos2, a, b, c
       complex(KIND=rkind) :: nsq
       real, parameter :: two = 2.0_rkind, four = 4.0_rkind

       call RLSDP_cold(eq, S ,D , P, R, L)

	   cos2 = cos(theta)**2
	   sin2 = one - cos2

!      Coefficients for A(n3)*(nsq)^2 + B(n3)*nsq + C(n3) = 0.
       a = S*sin2 + P*cos2
       b = -R*L*sin2 - P*S*(one + cos2)
       c = P*R*L

       nsq = n**2
       disp_fun_cold_n_theta = a*nsq**2 + b*nsq + c


       return
 end function disp_fun_cold_n_theta

