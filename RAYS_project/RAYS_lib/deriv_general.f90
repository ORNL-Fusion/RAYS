 subroutine deriv_general(eq, v, dddx, dddk, dddw)
!   calculates the derivatives of D with respect to k, r, omega.
!   v(1:3) = (x,y,z); v(4:6) = (kx, ky, kz),
!   dddx = dD/dx
!   dddk = dD/dk
!   dddw = dD/d(omega)

!   This uses the convention for storing the six independent values of
!   eps(i,j) as a six vector
!   eps(i) -> [ (eps(1,1), eps(2,2), eps(3,3), eps(1,2), eps(1,3), eps(2,3) ]

!   N./B. Send depsdx_h, depsdk_h, depsdw_h to depsdq_h by host association.

    use constants_m, only : rkind, zero
    use equilibrium_m, only : eq_point, write_eq_point
    use rf_m, only : omgrf, k0
    use ode_m, only : nv
    use write_matrix_m, only : write_matrix, write_vector

    implicit none

    type(eq_point), intent(in) :: eq
    real(KIND=rkind), intent(in) :: v(nv)
    real(KIND=rkind), intent(out) :: dddx(3), dddk(3), dddw

    real(KIND=rkind) :: rvec(3), kvec(3), nvec(3)
    real(KIND=rkind) :: n1, n3
    real(KIND=rkind) :: dn3dk(3), dn1dk(3), dn3dx(3), dn1dx(3), dn3dw, dn1dw
    complex(KIND=rkind) :: depsdk_h(6,3), depsdx_h(6,3), depsdw_h(6)

    complex(KIND=rkind) :: g(6), h1, h3
    complex(KIND=rkind) :: ddxc(3), ddkc(3), ddwc

    integer :: i

    rvec = v(1:3)
    kvec = v(4:6)
    nvec = kvec/k0

    n3 = dot_product(nvec, eq%bunit)
    n1 = sqrt( sum((nvec-n3*eq%bunit)**2) )

    dn1dk = (nvec-n3*eq%bunit)/(k0*n1)
    dn3dk = eq%bunit/k0

    dn3dx = matmul(eq%gradbunit,nvec)
    dn1dx = -n3/n1*dn3dx

    dn1dw = -n1/omgrf
    dn3dw = -n3/omgrf


! *******************************
    depsdw_h = zero
    depsdk_h = zero
    depsdx_h = zero

! Send args depsdx_h, depsdk_h, depsdw_h by host association
    call depsdq_h(eq, kvec)
!  write(*,*) ' '
!  call write_vector('deriv_general: depsdw_h', depsdw_h, 6)
!  call write_matrix('deriv_general: depsdx_h', depsdx_h, 6, 3)
!  call write_matrix('deriv_general: depsdk_h', depsdk_h, 6, 3)

! ********************************

    call g_and_h(eq, nvec, g, h1, h3)

! ********************************

    do i = 1,3
        ddkc(i) = sum( g(:)*depsdk_h(:,i) ) + h1*dn1dk(i) + h3*dn3dk(i)
        ddxc(i) = sum(g(:)*depsdx_h(:,i)) + h1*dn1dx(i) + h3*dn3dx(i)

        dddk(i) = ddkc(i)%re
        dddx(i) = ddxc(i)%re
    end do

    ddwc = sum(depsdw_h*g) + h1*dn1dw + h3*dn3dw
    dddw = ddwc%re

!  write(*,*) ' '
!  call write_vector('deriv_general: ddxc', ddxc, 3)
!  call write_vector('deriv_general: ddkc', ddkc, 3)
!  write(*,*) 'deriv_general: ddwc = ', ddwc

  return

! **************************************************************************

  contains

! **************************************************************************


! get args from host association above

 subroutine depsdq_h(eq, kvec)

!   Calculates the partial of derivatives of eps_h with respect to k, x, omega.
!   sums over  species is.
!   Get depsdx_h, depsdk_h, depsdw_h from above by host association.

    use constants_m, only : rkind, one, zero
    use species_m, only : nspec, spec_model

    implicit none

    type(eq_point), intent(in) :: eq
    real(KIND=rkind), intent(in) :: kvec(3)

    complex(KIND=rkind) :: depsdw_hs(6), depsdk_hs(6,3), depsdx_hs(6,3)

    integer :: is

    do is = 0, nspec
      plasma_model: select case (spec_model(is) )

          case ('cold')
            call depsdq_cold(eq,is, depsdw_hs, depsdk_hs, depsdx_hs)

          case('bessel')
           call depsdq_bessel(eq, is, kvec, depsdw_hs, depsdk_hs, depsdx_hs)

          case default
            write (0,*) 'depsdq_h: unimplemented species model =', spec_model(is)

      end select plasma_model

      depsdw_h = depsdw_h + depsdw_hs
      depsdk_h = depsdk_h + depsdk_hs
      depsdx_h = depsdx_h + depsdx_hs

    end do

    return
    end subroutine depsdq_h

! *************************************************************

subroutine depsdq_cold(eq, is, depsdw_hs, depsdk_hs, depsdx_hs)

!   calculates the partial of derivatives of the dielectric tensor with respect
!   to k, r, omega for species is using the cold plasma susceptibility
!   Includes collision frequency/omgrf  => nu_collision to avoid singularity at
!   fundamental resonance.

    use constants_m, only : rkind, zero, one, two, zi=>i
    use equilibrium_m, only : eq_point
    use species_m, only : nus
    use write_matrix_m, only : write_matrix, write_vector

    implicit none

    type(eq_point), intent(in) :: eq
    integer, intent(in) :: is

    complex(KIND=rkind), intent(out) :: depsdw_hs(6), depsdk_hs(6,3), depsdx_hs(6,3)

    complex(KIND=rkind) :: alpha_c, gamma_c, chis_6v(6)
    complex(KIND=rkind) :: depsdw(6), depsdk(6,3), depsdx(6,3)
    complex(KIND=rkind) :: d_chi_d_alpha_c(6), d_chi_d_gamma_c(6), d_alpha_c_dw, d_gamma_c_dw

    integer :: ivec

     alpha_c = eq%alpha(is)/cmplx(one,nus(is),rkind)**2
     gamma_c = eq%gamma(is)/cmplx(one,nus(is),rkind)

     chis_6v(1) = -alpha_c/(one-gamma_c**2)
     chis_6v(2) = chis_6v(1)
     chis_6v(3) = -alpha_c
     chis_6v(4) = -zi*gamma_c*alpha_c/(one-gamma_c**2)
     chis_6v(5) = zero
     chis_6v(6) = zero

!   Generate derivatives of dielectric tensor for this species ********

    d_chi_d_alpha_c(1) = -one/(one-gamma_c**2)
    d_chi_d_alpha_c(2) = d_chi_d_alpha_c(1)
    d_chi_d_alpha_c(3) = -one
    d_chi_d_alpha_c(4) = -zi*gamma_c/(one-gamma_c**2)
    d_chi_d_alpha_c(5) = zero
    d_chi_d_alpha_c(6) = zero

    d_chi_d_gamma_c(1) = -two*gamma_c*alpha_c/(one-gamma_c**2)**2
    d_chi_d_gamma_c(2) = d_chi_d_gamma_c(1)
    d_chi_d_gamma_c(3) = zero
    d_chi_d_gamma_c(4) = -zi*alpha_c*(one+gamma_c**2)/(one-gamma_c**2)**2
    d_chi_d_gamma_c(5) = zero
    d_chi_d_gamma_c(6) = zero

! Derivatives with respect to omgrf

    d_alpha_c_dw = -two*alpha_c/omgrf*cmplx(one,-nus(is),rkind)
    d_gamma_c_dw = -gamma_c/omgrf/cmplx(one,nus(is),rkind)

    depsdw(1) = d_chi_d_alpha_c(1)*d_alpha_c_dw + d_chi_d_gamma_c(1)*d_gamma_c_dw
    depsdw(2) = d_chi_d_alpha_c(2)*d_alpha_c_dw + d_chi_d_gamma_c(2)*d_gamma_c_dw
    depsdw(3) = d_chi_d_alpha_c(3)*d_alpha_c_dw + d_chi_d_gamma_c(3)*d_gamma_c_dw
    depsdw(4) = d_chi_d_alpha_c(4)*d_alpha_c_dw + d_chi_d_gamma_c(4)*d_gamma_c_dw
    depsdw(5) = zero
    depsdw(6) = zero

! Get Hermitian part

    call v6_Hermitian(depsdw, depsdw_hs)

!   Derivatives with respect to k.
!   For cold plasma all derivatives of chi w.r.t. k are zero

    do ivec = 1, 3
       depsdk(:,ivec) = zero
       depsdk_hs(:,ivec) = zero
    end do

!   dddx = dD/dx:
!   Derivatives with respect to space coordinates.
!  write(*,*) 'depsdq_cold: dne/dx /ne = ', eq%gradns(1,is)/eq%ns(is)
!  write(*,*) 'depsdq_cold: d|B|/dx /B = ', eq%gradbmag(1)/eq%bmag

    depsdx = zero
    do ivec = 1, 3
       depsdx(:,ivec) = d_chi_d_alpha_c(:)*alpha_c*eq%gradns(ivec,is)/eq%ns(is)+ &
                      & d_chi_d_gamma_c(:)*gamma_c*eq%gradbmag(ivec)/eq%bmag
       call v6_Hermitian( depsdx(:,ivec), depsdx_hs(:,ivec) )
    end do

    end subroutine depsdq_cold


! *************************************************************



! *************************************************************

  subroutine depsdq_bessel(eq, is, kvec, depsdw_hs, depsdk_hs, depsdx_hs)

!   calculates the partial of derivatives of the dielectric tensor with respect
!   to k, r, omega for species is using the Bessel function susceptibility

    use constants_m, only : rkind, zero, one, two, ten, zi=>i
    use equilibrium_m, only : eq_point
    use rf_m, only : omgrf
    use species_m, only : ms, nmins, nmaxs, n_limit
    use zfunctions_m, only : zfun0_real_arg
    use suscep_m, only : suscep_bessel
    use write_matrix_m, only : write_matrix, write_vector

    implicit none

    type(eq_point), intent(in) :: eq
    integer, intent(in) :: is
    real(KIND=rkind), intent(in) :: kvec(3)
    complex(KIND=rkind), intent(out) :: depsdw_hs(6), depsdk_hs(6,3), depsdx_hs(6,3)

    real(KIND=rkind) :: k1, k3, n1, n3
    real(KIND=rkind) :: bmag, ns, ts, omgp2, omgc
    real(KIND=rkind) :: bunit(3), gradbunit(3,3), gradbmag(3), gradns(3), gradts(3)

    real(KIND=rkind) :: iomgc
    real(KIND=rkind) :: vth, beta, lambda
    real(KIND=rkind), dimension(-n_limit:n_limit) :: xi
    complex(KIND=rkind), dimension(-n_limit:n_limit) :: zf, zfp, zfpp, ei, eip, eipp
    complex(KIND=rkind) :: a, chin(6,-n_limit:n_limit), chis_6v(6)

    complex(KIND=rkind) :: depsdw(6), depsdk(6,3), depsdx(6,3)
    complex(KIND=rkind), dimension(6) :: depsdb, depsdl
    complex(KIND=rkind), dimension(6,-n_limit:n_limit) :: depsdxi, ctmp
    complex(KIND=rkind) :: de1, de2, de3
    real(KIND=rkind) :: dbdk(3), dldk(3), dxidk(3,-n_limit:n_limit)
    real(KIND=rkind) :: dbdx(3), dldx(3), dxidx(3,-n_limit:n_limit)
    real(KIND=rkind) :: dbdw, dxidw

    complex(KIND=rkind) :: depsdk_xi(6,3), depsdx_xi(6,3)

    integer :: nmin, nmax, n, ivec
    real(KIND=rkind) :: half = one/two

    interface
      subroutine ebessel_dbb(z, nmin, nmax, ein, einp)
!        calculates exp(-z)*I_n(z) and exp(-z)*I'_n(z).
         use constants_m, only : rkind
         implicit none
         complex(KIND=rkind), intent(in) :: z
         integer, intent(in) :: nmin, nmax
         complex(KIND=rkind), intent(out) :: ein(nmin:nmax), einp(nmin:nmax)
      end subroutine ebessel_dbb
    end interface

    chin = zero
    chis_6v = zero

    bmag = eq%bmag
    ns = eq%ns(is)
    ts = eq%ts(is)
    omgp2 = eq%omgp2(is)
    omgc = eq%omgc(is)
    bunit = eq%bunit
    gradbunit = eq%gradbunit
    gradbmag = eq%gradbmag
    gradns(:) = eq%gradns(:,is)
    gradts(:) = eq%gradts(:,is)

!   nmin and nmax for this species
    nmin = nmins(is)
    nmax = nmaxs(is)

!   Sign of omgc.
    iomgc = sign(one, omgc)

!   Thermal speed.
    vth = sqrt( two*ts/ms(is) )

!   k_perp and k_parallel
    k3 = dot_product(kvec, eq%bunit)
    k1 = sqrt( sum((kvec-k3*eq%bunit)**2) )
    n1 = k1/k0
    n3 = k3/k0

!   Argument for Bessel functions.
    lambda = half * (k1*vth/omgc)**2

!   Arguments for plasma Z functions.
    xi = zero
    do n = nmin, nmax
       xi(n) = (omgrf-n*omgc) / (k3*vth)
    end do

!   Generate Bessel functions ******************************************
!   ei = exp(-lambda)*I_n(lambda), eip = exp(-lambda)*I'_n(lambda).

    call ebessel_dbb(cmplx(lambda, zero, rkind), nmin, nmax, ei(nmin:nmax), eip(nmin:nmax))

!   eipp = exp(-lambda)*I''_n(lambda).
!   Use the Bessel'e equation: y'' + (1/x)y' - (1+n^2/x^2)y = 0.
    do n = nmin, nmax
       eipp(n) = (one+(n/lambda)**2)*ei(n) - eip(n)/lambda
    end do

! write(*,*) 'depsdq_bessel: n1 = ', n1, '  n3 = ', n3, '  lambda = ', lambda
! call write_vector('depsdq_bessel: xi', xi(:), 2*nmax+1)
! call write_vector('depsdq_bessel: ei', ei, 2*nmax+1)
! call write_vector('depsdq_bessel: eip', eip, 2*nmax+1)
! call write_vector('depsdq_bessel: eipp', eipp, 2*nmax+1)

! Generate susceptibility

  if ( k3 /= zero ) then

       beta = omgp2/(omgrf*k3*vth)

!    Generate Plasma Z functions ***************************************

   do n = nmin, nmax
      zf(n) = zfun0_real_arg(xi(n), k3)
      zfp(n) = -two * ( one + xi(n)*zf(n) )
      zfpp(n) = -two * ( zf(n) + xi(n)*zfp(n) )
    end do

   do n = nmin, nmax
      chin(1,n) = n**2 * (ei(n)/lambda) * zf(n)
      chin(2,n) = chin(1,n) + two*lambda*(ei(n)-eip(n)) * zf(n)
      chin(3,n) = -ei(n) * xi(n) * zfp(n)

      chin(4,n) = zi * n * (eip(n)-ei(n)) * zf(n)
      chin(5,n) = -iomgc * sqrt(half/lambda) * n * ei(n) * zfp(n)
      chin(6,n) = iomgc * zi * sqrt(half*lambda) * (eip(n)-ei(n)) * zfp(n)
   end do

  else
!      For k3=0.
!      See Eq.(11-32). Here, an=A_n, but bn=B_n/k3.
       do n = nmin, nmax
          a = -one / (omgrf-n*eq%omgc(is))

          chin(1,n) = n**2 * (ei(n)/lambda) * a
          chin(2,n) = chin(1,n) + two*lambda*(ei(n)-eip(n)) * a
          chin(3,n) = ei(n) * a
          chin(4,n) = zi * n * (eip(n)-ei(n)) * a
          chin(5,n) = zero
          chin(6,n) = zero
       end do

        beta = eq%omgp2(is)/omgrf

  end if

!   Sum suscep over harmonic number

     chis_6v(1) = beta*sum(chin(1,nmin:nmax))
     chis_6v(2) = beta*sum(chin(2,nmin:nmax))
     chis_6v(3) = beta*sum(chin(3,nmin:nmax))
     chis_6v(4) = beta*sum(chin(4,nmin:nmax))
     chis_6v(5) = beta*sum(chin(5,nmin:nmax))
     chis_6v(6) = beta*sum(chin(6,nmin:nmax))

!   Derivatives of dielectric tensor eps with respect to beta:
!   depsdb = d(eps)/d(beta).

    depsdb = chis_6v/beta

!   Derivatives of dielectric tensor eps with respect to xi:
!   depsdxi = d(eps)/d(xi) = d(chi)/d(xi)
    do n = nmin, nmax
       depsdxi(1,n) = beta * n**2 / lambda * ei(n) * zfp(n)
       depsdxi(2,n) = depsdxi(1,n) + two*beta * lambda &
          & * (ei(n)-eip(n)) * zfp(n)
       depsdxi(3,n) = -beta * ei(n) * (xi(n)*zfpp(n)+zfp(n))
       depsdxi(4,n) = zi * beta * n * (eip(n)-ei(n)) * zfp(n)
       depsdxi(5,n) = -iomgc * beta * sqrt(half/lambda) * n &
          & * ei(n) * zfpp(n)
       depsdxi(6,n) = iomgc * zi * beta * sqrt(half*lambda) &
          & * (eip(n)-ei(n)) * zfpp(n)
    end do

! Derivatives with respect to w
!   dbdw = d(beta)/d(omega); dxidw = d(xi)/d(omega); d(lambda)/d(omega) = 0.

    dbdw = -beta/omgrf
    dxidw = xi(0)/omgrf

    do ivec = 1, 6
       depsdw(ivec) = depsdb(ivec)*dbdw &
       & + sum( depsdxi(ivec,nmin:nmax) )*dxidw
    end do

! Get Hermitian part

    call v6_Hermitian(depsdw, depsdw_hs)

!   Derivatives with respect to k.
!   dbdk = d(beta)/dk; dldk = d(lambda)/dk; dxidk = d(xi)/dk

!   Derivatives of dielectric tensor eps with respect to lambda:
!   depsdl = d(eps)/d(lambda) = d(chi)/d(lambda)

    do n = nmin, nmax
       de1 = -ei(n) + eip(n)
       de2 = -eip(n) + eipp(n)
       de3 = -ei(n) + two*eip(n) - eipp(n)
       ctmp(1,n) = n**2 * zf(n) *(-ei(n)/lambda**2 + de1/lambda)
       ctmp(2,n) = ctmp(1,n) + two * zf(n) * (-de1 + lambda * de3)
       ctmp(3,n) = -xi(n)*zfp(n) * de1
       ctmp(4,n) = -zi * n * zf(n) * de3
       ctmp(5,n) = -iomgc/sqrt(two) * n * zfp(n)*(-half/lambda**(1.5) * ei(n) &
                 & + de1/sqrt(lambda))
        ctmp(6,n) = iomgc *sqrt(half)* zi * zfp(n)* (half/sqrt(lambda) * de1  &
                 &  -sqrt(lambda) * de3 )
    end do

    ctmp = ctmp * beta

!   Sum over harmonics.
    depsdl = sum(ctmp(:,nmin:nmax),dim=2)

    dbdk(:) = -beta/k3 * bunit
    dldk(:) = two*lambda/k1**2 * (kvec-k3*bunit)

    depsdk_xi = zero
    do n = nmin, nmax
       dxidk(:,n) = -xi(n)/k3 * bunit
       do i = 1,6
           depsdk_xi(i,:) = depsdk_xi(i,:) + depsdxi(i,n)*dxidk(:,n)
       end do
    end do

    do ivec = 1, 3
       depsdk(:,ivec) = depsdb*dbdk(ivec) + depsdl*dldk(ivec) &
          & + depsdk_xi(:,ivec)

        call v6_Hermitian( depsdk(:,ivec), depsdk_hs(:,ivec) )
    end do

!   dddx = dD/dx:
!   Derivatives with respect to space coordinates.
!   dbdx = d(beta)/dx; dldx = d(lambda)/dx; dxidx = d(xi)/dx

       dbdx(:) = beta * ( gradns(:)/ns  &
        & - matmul(gradbunit,kvec)/k3 - gradts(:)/(two*ts) )
!  write(*,*) 'depsdq_bessel: dbdx = ', dbdx

       dldx(:) = lambda &
       & * ( -two*k3*matmul(gradbunit,kvec)/k1**2 &
       &      + gradts(:)/ts - two*gradbmag/bmag )
!  write(*,*) 'depsdq_bessel: dldx = ', dldx

    depsdx_xi = zero
    do n = nmin, nmax

       dxidx(:,n) = -xi(n) &
       & * ( matmul(gradbunit,kvec)/k3 + gradts(:)/(two*ts) ) &
       & - n*omgc/(k3*vth) * gradbmag/bmag

       do i = 1,6
            depsdx_xi(i,:) = depsdx_xi(i,:) + depsdxi(i,n)*dxidx(:,n)
        end do
    end do
!  write(*,*) 'depsdq_bessel: dxidx = ', dxidx

    do ivec = 1, 3

       depsdx(:,ivec) = depsdb*dbdx(ivec) + depsdl*dldx(ivec) &
          & + depsdx_xi(:,ivec)

        call v6_Hermitian( depsdx(:,ivec), depsdx_hs(:,ivec) )

    end do
!  write(*,*) 'depsdq_bessel: depsdx = ', depsdx
!  write(*,*) 'depsdq_bessel: depsdx_hs = ', depsdx_hs

  end subroutine depsdq_bessel

! *************************************************************

  subroutine g_and_h(eq, nvec, g, h1, h3)
! Generates the quantities g(1:6), h1 and h3 as defined in the notes from
! 11-4-99.  These are elements in the chain rule for dD/dq.

    use constants_m, only : rkind, zero, one, two, three, four
    use equilibrium_m, only : eq_point
    use rf_m, only : ray_dispersion_model
    use suscep_m, only : dielectric_cold, dielectric_general
    USE matrix3x3_m, only : hermitian3x3
    use write_matrix_m, only : write_matrix, write_vector

    implicit none
    type(eq_point), intent(in) :: eq
    real(KIND=rkind), intent(in) :: nvec(3)
    complex(KIND=rkind), intent(out) :: g(6), h1, h3

    real(KIND=rkind) :: n1, n3
    complex(KIND=rkind) :: n1c, n3c
    complex(KIND=rkind) :: eps(3,3)
    complex(KIND=rkind) :: eps_h(3,3)   ! Hermitian part of eps_6v
    complex(KIND=rkind) :: e1,e2,e3,e4,e5,e6

    n3 = dot_product(nvec, eq%bunit)
    n1 = sqrt( sum((nvec-n3*eq%bunit)**2) )
    n3c = cmplx(n3, zero, rkind)
    n1c = cmplx(n1, zero, rkind)

    dispersion_model: select case (ray_dispersion_model )
      case ('cold')
        call dielectric_cold(eq, eps)
      case('general')
!      N.B.  Dielectric routines take complex args n, but in ray tracing nvec is real
       call dielectric_general(eq, n1c, n3c, eps)
    end select dispersion_model

    eps_h = hermitian3x3(eps)

    e1 = eps_h(1,1); e2 = eps_h(2,2); e3 = eps_h(3,3); e4 = eps_h(1,2)
    e5 = eps_h(1,3); e6 = eps_h(2,3);

    g(1) = e6**2 + e2*(e3 - n1**2) - (e3 - n1**2)*(n1**2 + n3**2)
    g(2) = -e5**2 + e1*(e3 - n1**2) - two*e5*n1*n3 - e3*n3**2
    g(3) = e1*e2 + e4**2 - e1*n1**2 - e1*n3**2 - e2*n3**2 + n1**2*n3**2 + n3**4
    g(4) = two*(e3*e4 + e5*e6 + n1*(-e4*n1 + e6*n3))
    g(5) = two*(e4*e6 - e2*(e5 + n1*n3) + (e5 + n1*n3)*(n1**2 + n3**2))
    g(6) = two*(e4*(e5 + n1*n3) + e6*(e1 - n3**2))

    h1 = two*(-(e4**2*n1) + e5**2*n1 - e2*e5*n3 + e4*e6*n3 + three*e5*n1**2*n3 + &
          & e3*n1*n3**2 + e5*n3**3 + e1*n1*(-e2 - e3 + two*n1**2 + n3**2))

    h3 = two*(e4*e6*n1 + e5*n1**3 - e1*e3*n3 + e5**2*n3 - e6**2*n3 + e1*n1**2*n3 + &
         & e3*n1**2*n3 + three*e5*n1*n3**2 + two*e3*n3**3 - e2*(e5*n1 + e3*n3))

    return
  end subroutine g_and_h

  end subroutine deriv_general


! ********************************************************************************

    subroutine v6_Hermitian(v,v_h)

    use constants_m, only : rkind, one, two

    implicit none

    complex(KIND=rkind), intent(in) :: v(6)
    complex(KIND=rkind), intent(out) :: v_h(6)

    v_h(1) = v(1) + conjg( v(1) )
    v_h(2) = v(2) + conjg( v(2) )
    v_h(3) = v(3) + conjg( v(3) )
    v_h(4) = v(4) - conjg( v(4) )
    v_h(5) = v(5) + conjg( v(5) )
    v_h(6) = v(6) - conjg( v(6) )

    v_h = v_h/two

    return
    end subroutine v6_Hermitian

! ********************************************************************************

    subroutine v6_3x3(x_6v, x_3x3)

    use constants_m, only : rkind, one, two

    implicit none

    complex(KIND=rkind), intent(in) :: x_6v(6)
    complex(KIND=rkind), intent(out) :: x_3x3(3,3)


     x_3x3(1,1) = x_6v(1)
     x_3x3(2,2) = x_6v(2)
     x_3x3(3,3) = x_6v(3)
     x_3x3(1,2) = x_6v(4)
     x_3x3(1,3) = x_6v(5)
     x_3x3(2,3) = x_6v(6)

     x_3x3(2,1) = -x_3x3(1,2)
     x_3x3(3,1) = x_3x3(1,3)
     x_3x3(3,2) = -x_3x3(2,3)

    return
    end subroutine v6_3x3

