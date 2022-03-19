 subroutine deriv_num(eq0, v, dddx, dddk, dddw)
!   numerically calculates the derivatives of D with respect to k, r, omega.
!   v(1:3) = (x,y,z); v(4:6) = (kx, ky, kz),
!   dddx = dD/dx, dddk = dD/dk, dddw = dD/d(omega).
!
!   N.B. omgrf and k0 are module variables and therefore global. They get changed
!        and reset in this subroutine.  Therefore this subroutine, and omgrf, k0 
!        are not presently thread safe!!  In particular equilibrium uses omgrf to
!        calculate alpha and beta so equilibrium() must be called after omgrf is
!        changed.

    use constants_m, only : rkind, clight
    use equilibrium_m, only : equilibrium, eq_point, write_eq_point
    use rf_m, only : omgrf, k0, ray_dispersion_model
    use ode_m, only : nv
    use diagnostics_m, only : message_unit, verbosity
 
    implicit none

    type(eq_point), intent(in) :: eq0
    real(KIND=rkind), intent(in) :: v(nv)
    real(KIND=rkind), intent(out) :: dddx(3), dddk(3), dddw

    real(KIND=rkind) :: rvec(3), kvec(3)
    real(KIND=rkind) :: rvec0(3), kvec0(3), omgrf0
    type(eq_point) :: eq_plus, eq_minus
    real(KIND=rkind) :: det_plus, det_minus, delta, change, k_plus, k_minus

    integer :: i

! write(*,*) 'v = ', v
!   Save parameters.
    rvec0 = v(1:3); kvec0 = v(4:6); omgrf0 = omgrf
    kvec = kvec0

!   Step for computing the derivaties.
    delta = 1.e-6

!   Derivatives of D with respect to r.
    do i = 1, 3
       rvec = rvec0
!       change = max(delta, abs(delta*rvec(i)))/2.
       change = delta
       rvec(i) = rvec0(i) + change
       call equilibrium(rvec, eq_plus)
!        write(*,*) ' '
!        write(*,*) 'deriv_num: rvec = ', rvec
!        call write_eq_point(eq_plus)
       rvec(i) = rvec0(i) - change
       call equilibrium(rvec, eq_minus)
!        write(*,*) ' '
!        write(*,*) 'deriv_num: rvec = ', rvec
!        call write_eq_point(eq_minus)
       det_plus = determ(eq_plus)
       det_minus = determ(eq_minus)
       dddx(i) = (det_plus-det_minus)/(2.*change)
    end do

!   Derivatives of D with respect to k.
    do i = 1, 3
       kvec = kvec0
       change = max(delta, abs(delta*kvec(i)))/2.
       kvec(i) = kvec0(i) + change
       det_plus = determ(eq0)
       kvec(i) = kvec0(i) - change
       det_minus = determ(eq0)
       dddk(i) = (det_plus-det_minus)/(2.*change)
    end do

!   Derivative of D with respect to omega.
    kvec = kvec0
    omgrf = omgrf0 * (1.+delta/2.)
    k0 = omgrf/clight
    call equilibrium(rvec0, eq_plus)
    det_plus = determ(eq_plus)
    omgrf = omgrf0 * (1.-delta/2.)
    k0 = omgrf/clight
    call equilibrium(rvec0, eq_minus)
    det_minus = determ(eq_minus)
    dddw = (det_plus-det_minus)/(omgrf0*delta)

! Reset omgrf and k0 to module values
    omgrf = omgrf0 
    k0 = omgrf/clight
    
    verb: if (verbosity > 3) then
       write(*,*) ''
       write(*,'(a,1p3e12.4)') 'deriv_num: dddx =', dddx
       write(*,'(a,1p3e12.4)') 'deriv_num: dddk =', dddk
       write(*,'(a,1p1e12.4)') 'deriv_num: dddw =', dddw
       write(*,*) ''   

    end if verb

    return

 contains
    
    real(KIND=rkind) function determ(eq)
!      calculates the determinant for epsn = eps + nn -n^2I.

       use equilibrium_m, only : eq_point
       use rf_m, only : ray_dispersion_model
       use suscep_m, only : dielectric_cold

       implicit none 

       type(eq_point), intent(in) :: eq       
       complex(KIND=rkind) :: eps(3,3), eps_h(3,3), epsn(3,3) 
       complex(KIND=rkind) :: ctmp
       real(KIND=rkind) :: k1, k3, n(3)
       integer :: i, j

       k3 = dot_product(kvec, eq%bunit)
       k1 = sqrt( sum((kvec-k3*eq%bunit)**2) )
!      Refractive index.
       n(1) = k1/k0; n(2) = 0.; n(3) = k3/k0
       call dielectric_cold(eq,eps)

!      Hermitian part.
       eps_h = .5 * (eps + conjg(transpose(eps)))

!      epsn = eps + nn -n^2I:
!      epsn.E = eps.E + n x n x E = (eps + nn -n^2I).E,
!      where E = (Ex,Ey,Ez)^T and I is the unit 3X3 tensor.
       do i = 1, 3; do j = 1, 3
          epsn(i,j) = eps_h(i,j) + n(i)*n(j) - int(i/j)*int(j/i)*sum(n**2)
       end do; end do

!      Determinant for 3X3 epsn:
       ctmp = &
          &   epsn(3,3)*(epsn(1,1)*epsn(2,2)-epsn(2,1)*epsn(1,2)) &
          & - epsn(3,2)*(epsn(1,1)*epsn(2,3)-epsn(2,1)*epsn(1,3)) &
          & + epsn(3,1)*(epsn(1,2)*epsn(2,3)-epsn(2,2)*epsn(1,3))

!      For a Hermitian matrix, the imaginary part of its determinant vanishes.
       if ( abs(ctmp%im) > 1.e-7 ) then
          write(0,'(a,1p1e12.4)') 'RESIDUAL: Im(det) =', ctmp%im
          stop 1
       end if

       determ = 0.
       if ( ray_dispersion_model == 'cold' ) then
!         For a cold plasma.  Note that the factor here corresponds to that
!         used in DERIV_COLD.
          determ = ctmp%re * product(1.-eq%gamma**2)
       else if ( ray_dispersion_model == 'warm' ) then
!         For a warm plasma.
          determ = ctmp%re
       end if

       return
    end function determ

 end subroutine deriv_num
 
