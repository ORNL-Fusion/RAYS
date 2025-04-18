 subroutine solve_cold_n1sq_vs_n3(eq, n3, n1sq)
! calculates cold plasma perpendicular refractive index squared versus parallel
! refractive index <==> n3. n1sq may be < 0 in which case nx is complex(KIND=rkind).  It's
!  the users responsibility to handle that properly.
!
! The dispersion relation is quadratic in nx^2.  The roots are classified according to
! whether the sign in front of the discriminant is + or -.  Roots are also classified as
! to which |nx| is larger -> slow mode, or smaller -> fast mode.  This returns both choices:
! n1sq(1) -> plus
! n1sq(2) -> minus
! n1sq(3) -> fast
! n1sq(4) -> slow

! N.B. refractive index convention
!   n1 = component of nvec perpendicular to B to be solved for (e.g. x or radial component)
!   n2 = transverse component of nvec i.e. perpendicular to both B and n1
!   n3 = parallel component of kvec (nvec)
!   nperp_sq = square of perpendicular component of nvec = n1**2 +n2**2

    use constants_m, only : rkind
    use diagnostics_m, only : message
    use equilibrium_m, only : eq_point
!    use suscep_m, only :  dielectric_cold
    use suscep_m, only :  RLSDP_cold

       implicit none

!      Derived type containing equilibrium data for a spatial point in the plasma
       type(eq_point), intent(in) :: eq

       real(KIND=rkind), intent(in) :: n3
       complex(KIND=rkind), intent(out) :: n1sq(4)

!       complex(KIND=rkind) :: eps(3,3)
       real(KIND=rkind) :: S ,D , P,  R, L

       real(KIND=rkind) :: a, b, c, sgn_b, discr
       complex(KIND=rkind) :: sqrt_d
       real, parameter :: two = 2.0_rkind, four = 4.0_rkind


! !      Calculate cold plasma dielectric tensor.
!        call dielectric_cold(eq, eps)
!
!  !    write(*,*) 'eps = ', eps
!
! !      Coefficients for A(n3)*(n1sq)^2 + B(n3)*n1sq + C(n3) = 0.
!        a = real(eps(1,1), KIND=rkind)
!        b = real(-(eps(1,1)**2+eps(1,2)**2+eps(1,1)*eps(3,3)) &
!           & + (eps(1,1)+eps(3,3))*n3**2, KIND=rkind)
!        c = real(( eps(1,1)**2+eps(1,2)**2 - two*eps(1,1)*n3**2 &
!           & + n3**4 ) *  eps(3,3), KIND=rkind)
       call RLSDP_cold(eq, S ,D , P, R, L)

!      Coefficients for A(n3)*(n1sq)^2 + B(n3)*n1sq + C(n3) = 0.
       a = S
       b = -R*L - P*S +n3**2*(P+S)
       c = P*(n3**2 - R)*(n3**2 - L)
       discr = b**2-four*a*c

       sgn_b = sign(real(1., KIND=rkind),b)   ! SIGN(A,B)=|A| if B>=0 and SIGN(A,B)=-|A| if B<0.
       sqrt_d = sqrt(cmplx(discr, KIND=rkind))

       if ( (b**2-4.*a*c) < 0. ) then
          call message('disp_solve_cold_n1sq_vs_n3:  evanescent root', 0)
          call message('disp_solve_cold_n1sq_vs_n3: b**2-4.*a*c', discr, 0)
       end if

       if (sgn_b < 0.) then
           n1sq(1) = ( -b + sqrt_d) / (two*a) ! plus root
           n1sq(2) = two*c/(-b + sqrt_d)      ! minus root
       else
           n1sq(2) = ( -b - sqrt_d) / (two*a) ! minus root
           n1sq(1) = two*c/(-b - sqrt_d)      ! plus root
       end if

!   write(*,*) 'a = ', a, ' b = ', b, ' c = ', c
!   write(*,*) 'n1sq = ', n1sq, '    sqrt_d = ', sqrt_d

!      Store smaller root in (fast wave) into n1sq(3) and larger root (slow wave) into n1sq(4).
       if ( abs(n1sq(1)) <= abs(n1sq(2)) ) then
           n1sq(3) = n1sq(1) ! fast root
           n1sq(4) = n1sq(2) ! slow root
       else
          n1sq(3) = n1sq(2) ! fast root
          n1sq(4) = n1sq(1) ! slow root
       end if

       return
 end subroutine solve_cold_n1sq_vs_n3
!****************************************************************************

complex function residual_cold_n1_n3(eq, n1, n3)
! calculates the residual of cold plasma dispersion relation versus perpendicular, n1 and
! parallel, n3, refractive index indices
! N.B. n1 is complex, n3 is real

    use constants_m, only : rkind, one
    use diagnostics_m, only : message
    use equilibrium_m, only : eq_point
    use suscep_m, only :  RLSDP_cold

       implicit none

!      Derived type containing equilibrium data for a spatial point in the plasma
       type(eq_point), intent(in) :: eq

       complex(KIND=rkind), intent(in) :: n1
       real(KIND=rkind), intent(in) :: n3

       real(KIND=rkind) :: S ,D , P,  R, L
       real(KIND=rkind) :: a, b, c
       complex(KIND=rkind) :: n1sq

       call RLSDP_cold(eq, S ,D , P, R, L)

!      Coefficients for A(n3)*(n1sq)^2 + B(n3)*n1sq + C(n3) = 0.
       a = S
       b = -R*L - P*S +n3**2*(P+S)
       c = P*(n3**2 - R)*(n3**2 - L)

       n1sq = n1**2
       residual_cold_n1_n3 = a*n1sq**2 + b*n1sq + c

       return
 end function residual_cold_n1_n3

