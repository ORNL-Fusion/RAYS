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

    use constants_m, only : rkind, zero, one, two, ten, clight
    use equilibrium_m, only : equilibrium, eq_point, write_eq_point
    use rf_m, only : omgrf, k0, ray_dispersion_model
    use ode_m, only : nv
    use diagnostics_m, only : verbosity

    implicit none

    type(eq_point), intent(in) :: eq0
    real(KIND=rkind), intent(in) :: v(nv)
    real(KIND=rkind), intent(out) :: dddx(3), dddk(3), dddw

    real(KIND=rkind) :: rvec(3), kvec(3)
    real(KIND=rkind) :: rvec0(3), kvec0(3), omgrf0
    type(eq_point) :: eq_plus, eq_minus
    real(KIND=rkind) :: det_plus, det_minus, delta, change

    integer :: i

!  write(*,*) 'deriv_num: v = ', v
!   Save parameters.
    rvec0 = v(1:3); kvec0 = v(4:6); omgrf0 = omgrf
    kvec = kvec0

!   Step for computing the derivaties.
    delta = 1.0d-8

!   Derivatives of D with respect to r.
    do i = 1, 3
       rvec = rvec0
!        write(*,*) ' '
!        write(*,*) 'deriv_num: rvec0 = ', rvec0
!       change = max(delta, abs(delta*rvec(i)))/2.
       change = delta
       rvec(i) = rvec0(i) + change
!        write(*,*) ' '
!        write(*,*) 'deriv_num: rvec(+) = ', rvec
       call equilibrium(rvec, eq_plus)
!        call write_eq_point(eq_plus)
       rvec(i) = rvec0(i) - change
!        write(*,*) ' '
!        write(*,*) 'deriv_num: rvec(-) = ', rvec
       call equilibrium(rvec, eq_minus)
!        call write_eq_point(eq_minus)
       det_plus = determ(eq_plus)
       det_minus = determ(eq_minus)
       dddx(i) = (det_plus-det_minus)/(2.*change)
    end do

!   Derivatives of D with respect to k.
    do i = 1, 3
       kvec = kvec0
       change = max(delta, abs(delta*kvec(i)))/two
       kvec(i) = kvec0(i) + change
       det_plus = determ(eq0)
       kvec(i) = kvec0(i) - change
       det_minus = determ(eq0)
       dddk(i) = (det_plus-det_minus)/(two*change)
    end do

!   Derivative of D with respect to omega.
    kvec = kvec0
    omgrf = omgrf0 * (one+delta/two)
    k0 = omgrf/clight
    call equilibrium(rvec0, eq_plus)
    det_plus = determ(eq_plus)
    omgrf = omgrf0 * (one-delta/two)
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
! calculates the determinant for epsn = eps + nn -n^2I.
! N.B. Gets arg kvec from above by host association

       use constants_m, only : rkind, zero, ten, epsmach
       use equilibrium_m, only : eq_point
       use rf_m, only : ray_dispersion_model
       use suscep_m, only : dielectric_cold, dielectric_general
       use matrix3x3_m, only : hermitian3x3, anti_hermitian3x3, determinant3x3

       implicit none

       type(eq_point), intent(in) :: eq
       complex(KIND=rkind) :: eps(3,3), eps_h(3,3), epsn(3,3)
       complex(KIND=rkind) :: ctmp
       real(KIND=rkind) :: k1, k3, n(3)
       complex(KIND=rkind) :: nc(3)
       integer :: i, j

       eps = zero
       eps_h = zero
       epsn = zero

       k3 = dot_product(kvec, eq%bunit)
       k1 = sqrt( sum((kvec-k3*eq%bunit)**2) )
!      Refractive index.
       n(1) = k1/k0; n(2) = zero; n(3) = k3/k0
       nc(1) = cmplx(n(1), zero, rkind)
       nc(2) = zero
       nc(3) = cmplx(n(3), zero, rkind)
    deriv_name: select case (trim(ray_dispersion_model))
       case ('cold')
           call dielectric_cold(eq, eps)
       case ('general')
           call dielectric_general(eq, nc(1), nc(3), eps)
   end select deriv_name

!      Hermitian part.
       eps_h = hermitian3x3(eps)
!      epsn = eps + nn -n^2I:
!      epsn.E = eps.E + n x n x E = (eps + nn -n^2I).E,
!      where E = (Ex,Ey,Ez)^T and I is the unit 3X3 tensor.

       do i = 1, 3; do j = 1, 3
          epsn(i,j) = eps_h(i,j) + n(i)*n(j) - int(i/j)*int(j/i)*sum(n**2)
       end do; end do

!      Determinant for 3X3 epsn:
       ctmp = determinant3x3(epsn)

!      For a Hermitian matrix, the imaginary part of its determinant vanishes.
       if ( abs(aimag(ctmp)) > ten*epsmach ) then
          write(0,'(a,1p1e12.4)') 'RESIDUAL: Im(det) = ', aimag(ctmp)
          stop 1
       end if

       determ = zero
       if ( ray_dispersion_model == 'cold' ) then
!         To compare with ray_deriv_name = 'cold' multiply by the factor
!         product(one-eq%gamma**2), which is  used in DERIV_COLD to eliminate
!         resonant denominators, uncomment next line.
!           determ = ctmp%re * product(one-eq%gamma**2)
          determ = ctmp%re
       else if ( trim(ray_dispersion_model) == 'general' ) then
!         For a warm plasma.
          determ = ctmp%re
       end if

       return
    end function determ

 end subroutine deriv_num

