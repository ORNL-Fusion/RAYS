 subroutine disp_solve_cold_nxsq_vs_nz(nz, nxsq)
! calculates cold plasma perpendicular refractive index squared versus parallel 
! refractive index. It may be <0 in which case nx is complex.  It's the users responsibility
! to handle that properly.
!
! The dispersion relation is quadratic in nx^2.	 The roots are classified according to
! whether the sign in front of the discriminant is + or -.  Roots are also classified as
! to which |nx| is larger -> slow mode, or smaller -> fast mode.  This returns both choices:
! nxsq(1) -> plus
! nxsq(2) -> minus
! nxsq(3) -> fast
! nxsq(4) -> slow

	use diagnostics_m, only : message
	use suscep_m, only : dielectric_cold, eps_cold

	   implicit none 

	   real, intent(in) :: nz
	   complex, intent(out) :: nxsq(4)

	   complex :: eps(3,3)
	   real :: a, b, c, sgn_b
	   complex :: sqrt_d


!	   Calculate cold plasma dielectric tensor.
	   call dielectric_cold
	   eps = eps_cold

!	   Coefficients for A(nz)*(nxsq)^2 + B(nz)*nxsq + C(nz) = 0.
	   a = real( eps(1,1) )
	   b = real( -(eps(1,1)**2+eps(1,2)**2+eps(1,1)*eps(3,3)) &
		  & + (eps(1,1)+eps(3,3))*nz**2 )
	   c = real( ( eps(1,1)**2+eps(1,2)**2 - 2.*eps(1,1)*nz**2 &
		  & + nz**4 ) *	 eps(3,3) )
		  
	   sgn_b = sign(1.,b)	! SIGN(A,B)=|A| if B>=0 and SIGN(A,B)=-|A| if B<0.
	   sqrt_d = csqrt(cmplx(b**2-4.*a*c))

	   if ( (b**2-4.*a*c) < 0. ) then
		  call message('disp_solve_cold_nxsq_vs_nz:	 evanescent root', 0)
		  call message('disp_solve_cold_nxsq_vs_nz: b**2-4.*a*c', b**2-4.*a*c, 0)
	   end if

	   if (sgn_b < 0.) then
		   nxsq(1) = ( -b + sqrt_d) / (2.*a) ! plus root
		   nxsq(2) = 2.*c/(-b + sqrt_d)		 ! minus root
	   else
		   nxsq(2) = ( -b - sqrt_d) / (2.*a) ! minus root
		   nxsq(1) = 2.*c/(-b - sqrt_d)		 ! plus root
	   end if
	   
!   write(*,*) 'a = ', a, ' b = ', b, ' c = ', c
!   write(*,*) 'nxsq = ', nxsq, '	 sqrt_d = ', sqrt_d

!	   Store smaller root in (fast wave) into nxsq(3) and larger root (slow wave) into nxsq(4).
	   if ( abs(nxsq(1)) <= abs(nxsq(2)) ) then
		   nxsq(3) = nxsq(1) ! fast root
		   nxsq(4) = nxsq(2) ! slow root
	   else
		  nxsq(3) = nxsq(2) ! fast root
		  nxsq(4) = nxsq(1) ! slow root
	   end if

	   return
 end subroutine disp_solve_cold_nxsq_vs_nz
	
