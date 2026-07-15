 subroutine deriv_general(eq, v, dddx, dddk, dddw, depsdw_h)
!   calculates the derivatives of D with respect to k, r, omega.
!   v(1:3) = (x,y,z); v(4:6) = (kx, ky, kz),
!   dddx = dD/dx
!   dddk = dD/dk
!   dddw = dD/d(omega)

!   depsdw_h = Hermitian part of d(eps)/d(omega), and warm plasma eps are
!   used in routine POYNTING().  A side effect of routine depsdq_ below is
!   to calculate these and store them in the appropriate module.

!   This uses the convention for storing the six independent values of
!   eps(i,j) as a six vector
!   eps(i) -> [ (eps(1,1), eps(2,2), eps(3,3), eps(1,2), eps(1,3), eps(2,3) ]

    use constants_m, only : rkind
    use equilibrium_m, only : eq_point
    use rf_m, only : omgrf, k0
    use ode_m, only : nv

    implicit none

    type(eq_point), intent(in) :: eq
    real(KIND=rkind), intent(in) :: v(nv)
    real(KIND=rkind), intent(out), optional :: dddx(3), dddk(3), dddw
    complex(KIND=rkind), intent(out), optional :: depsdw_h(6)

    real(KIND=rkind) :: rvec(3), kvec(3), nvec(3)
    complex(KIND=rkind) :: n1, n3
    real(KIND=rkind) :: dn3dk(3), dn1dk(3), dn3dx(3), dn1dx(3), dn3dw, dn1dw
    complex(KIND=rkind) :: depsdk_h(6,3), depsdx_h(6,3)

    complex(KIND=rkind) :: g(6), h1, h3

    integer :: i

    logical :: only_depsdw

!   If optional output variables dddx, dddk, dddw are absent then just need to
!   calculate depsdw_h.  So set flag

    only_depsdw = .false.
    if ( .not. ( present(dddx) .or. present(dddk) .or. present(dddw) ) ) &
        &   only_depsdw = .true.

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

    call depsdq_h(eq, kvec, depsdx_h, depsdk_h, depsdw_h)

    if (only_depsdw) return

! ********************************

    call g_and_h(eq, nvec, g, h1, h3)

! ********************************

    do i = 1,3

    dddk(i) = real( sum( g(:)*depsdk_h(:,i) ) + h1*dn1dk(i) + h3*dn3dk(i) )
    dddx(i) = real( sum( g(:)*depsdx_h(:,i) ) + h1*dn1dx(i) + h3*dn3dx(i) )

    end do

    dddw = real( sum(depsdw_h*g) + h1*dn1dw + h3*dn3dw )


    return

! **************************************************************************

    contains

! **************************************************************************



 subroutine depsdq_h(eq, kvec, depsdx_h, depsdk_h, depsdw_h)

!   Calculates the partial of derivatives of eps_h with respect to k, x, omega.
!   sums over  species is.

    use constants_m, only : rkind, one, zero
    use species_m, only : nspec, spec_model

    implicit none

    type(eq_point), intent(in) :: eq
    real(KIND=rkind), intent(in) :: kvec(3)
    complex(KIND=rkind), intent(out) :: depsdw_h(6), depsdk_h(6,3), depsdx_h(6,3)

    complex(KIND=rkind) :: depsdw_hs(6), depsdk_hs(6,3), depsdx_hs(6,3)

    integer :: is

    depsdw_h = zero
    depsdk_h = zero
    depsdx_h = zero

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
!   to k, r, omega for  species is using the cold plasma susceptibility
!   Also generates and stores the contribution of species is to chis(:,:,:)
!   Includes collision frequency/omgrf  => nu_collision to avoid singularity at
!   fundamental resonance.

    use constants_m, only : rkind, zero, one, two, zi=>i
    use equilibrium_m, only : eq_point
    use species_m, only : nus

    implicit none

    type(eq_point), intent(in) :: eq
    integer, intent(in) :: is

    real(KIND=rkind) :: iomgc
    complex(KIND=rkind) :: alpha_c, gamma_c, chis_6v(6), chis(3,3)
    complex(KIND=rkind), intent(out) :: depsdw_hs(6), depsdk_hs(6,3), depsdx_hs(6,3)

    complex(KIND=rkind) :: depsdw(6), depsdk(6,3), depsdx(6,3)
    complex(KIND=rkind) :: d_chi_d_alpha_c(6), d_chi_d_gamma_c(6), d_alpha_c_dw, d_gamma_c_dw

    integer :: ivec

!   Sign of omgc.
    iomgc = sign(one,eq%omgc(is))
!
     alpha_c = eq%alpha(is)/cmplx(one,nus(is))**2
     gamma_c = eq%gamma(is)/cmplx(one,nus(is))

     chis_6v(1) = -alpha_c/(1.-gamma_c**2)
     chis_6v(2) = chis_6v(1)
     chis_6v(3) = -alpha_c
     chis_6v(4) = -zi*gamma_c*chis_6v(1)
     chis_6v(5) = zero
     chis_6v(6) = zero

     chis(1,1) = chis_6v(1)
     chis(2,2) = chis_6v(2)
     chis(3,3) = chis_6v(3)
     chis(1,2) = chis_6v(4)
     chis(1,3) = chis_6v(5)
     chis(2,3) = chis_6v(6)

     chis(2,1) = -chis(1,2)
     chis(3,1) = chis(1,3)
     chis(3,2) = -chis(2,3)

!    call v6_3x3( chis_6v, chis(:,:,is) )

!   Generate derivatives of dielectric tensor for this species ********

    d_chi_d_alpha_c(1) = -one/(one-gamma_c**2)
    d_chi_d_alpha_c(2) = d_chi_d_alpha_c(1)
    d_chi_d_alpha_c(3) = -one
    d_chi_d_alpha_c(4) = -zi*gamma_c/(one-gamma_c**2)
    d_chi_d_alpha_c(5) = zero
    d_chi_d_alpha_c(6) = zero

    d_chi_d_gamma_c(1) = -two*chis_6v(1)*gamma_c/(one-gamma_c**2)
    d_chi_d_gamma_c(2) = d_chi_d_gamma_c(1)
    d_chi_d_gamma_c(3) = zero
    d_chi_d_gamma_c(4) = -zi*chis_6v(1)*(1.-gamma_c**2/(one-gamma_c**2))
    d_chi_d_gamma_c(5) = zero
    d_chi_d_gamma_c(6) = zero

! Derivatives with respect to omgrf
!

    d_alpha_c_dw = -two*alpha_c/omgrf*cmplx(1.,-nus(is))
    d_gamma_c_dw = -gamma_c/omgrf/cmplx(1.,nus(is))

    depsdw(1) = d_chi_d_alpha_c(1)*d_alpha_c_dw + d_chi_d_gamma_c(1)*d_gamma_c_dw
    depsdw(2) = d_chi_d_alpha_c(2)*d_alpha_c_dw + d_chi_d_gamma_c(2)*d_gamma_c_dw
    depsdw(3) = d_chi_d_alpha_c(3)*d_alpha_c_dw + d_chi_d_gamma_c(3)*d_gamma_c_dw
    depsdw(4) = d_chi_d_alpha_c(4)*d_alpha_c_dw + d_chi_d_gamma_c(4)*d_gamma_c_dw
    depsdw(5) = zero
    depsdw(6) = zero

! Get Hermitian part

    call v6_Hermitian(depsdw, depsdw_hs)

! ****************************************************************************

!   If we only need depsdw (for Poynting) exit routine here.
    if (only_depsdw) return

! ****************************************************************************


!   Derivatives with respect to k.
!   For cold plasma all derivatives of chi w.r.t. k are zero



    do ivec = 1, 3
       depsdk(:,ivec) = zero
       call v6_Hermitian( depsdk(:,ivec), depsdk_hs(:,ivec) )
    end do



!   dddx = dD/dx:
!   Derivatives with respect to space coordinates.

    do ivec = 1, 3
       depsdx(:,ivec) = d_chi_d_alpha_c*eq%gradns(ivec,is)+d_chi_d_gamma_c*eq%gradbmag(ivec)
       call v6_Hermitian( depsdx(:,ivec), depsdx_hs(:,ivec) )
    end do

    end subroutine depsdq_cold


! *************************************************************



! *************************************************************

subroutine depsdq_bessel(eq, is, kvec, depsdw_hs, depsdk_hs, depsdx_hs)

!   calculates the partial of derivatives of the dielectric tensor with respect
!   to k, r, omega for species is using the Bessel function susceptibility

    use constants_m, only : rkind, zero, one, two, zi=>i
    use equilibrium_m, only : eq_point
    use rf_m, only : omgrf
    use species_m, only : ms, nmins, nmaxs, n_limit
    use zfunctions_m, only : zfun0_real_arg

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
    complex(KIND=rkind) :: a, b, chin(6,-n_limit:n_limit), chis_6v(6), chis(3,3)

    complex(KIND=rkind) :: depsdw(6), depsdk(6,3), depsdx(6,3)
    complex(KIND=rkind), dimension(6) :: depsdb, depsdl
    complex(KIND=rkind), dimension(6,-n_limit:n_limit) :: depsdxi, ctmp

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

! 	k_perp and k_parallel
    k3 = dot_product(kvec, eq%bunit)
    k1 = sqrt( sum((kvec-k3*eq%bunit)**2) )
	n1 = k1/k0
	n3 = k3/k0

!   Argument for Bessel's function.
    lambda = half * (k1*vth/omgc)**2

!   Argument for plasma Z function.
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

! Generate susceptibility


    if ( k3 /= zero ) then

       beta = omgp2/(omgrf*k3*vth)


!    Generate Plasma Z functions ***************************************

    do n = nmin, nmax
       zf(n) = zfun0_real_arg(xi(n), k3)

!     Z' = dZ/d(xi) & Z'' = d^2Z/d(xi)^2.
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
          a = -one / (omgrf-n*omgc)
          b = -half * (vth/(omgrf-n*omgc))**2

!         Eq.(10-57):
          chin(1,n) = n**2 * (ei(n)/lambda) * a
          chin(2,n) = chin(1,n) + two*lambda*(ei(n)-eip(n)) * a
          chin(3,n) = two*(omgrf-n*omgc)/vth**2 * ei(n) * b

          chin(4,n) = zi * n * (eip(n)-ei(n)) * a
          chin(5,n) = zero
          chin(6,n) = zero

       end do

    end if

!   Sum suscep over harmonic number

     chis_6v(1) = beta*sum(chin(1,nmin:nmax))
     chis_6v(2) = beta*sum(chin(2,nmin:nmax))
     chis_6v(3) = beta*sum(chin(3,nmin:nmax))
     chis_6v(4) = beta*sum(chin(4,nmin:nmax))
     chis_6v(5) = beta*sum(chin(5,nmin:nmax))
     chis_6v(6) = beta*sum(chin(6,nmin:nmax))

     chis(1,1) = chis_6v(1)
     chis(2,2) = chis_6v(2)
     chis(3,3) = chis_6v(3)
     chis(1,2) = chis_6v(4)
     chis(1,3) = chis_6v(5)
     chis(2,3) = chis_6v(6)

     chis(2,1) = -chis(1,2)
     chis(3,1) = chis(1,3)
     chis(3,2) = -chis(2,3)

!    call v6_3x3( chis_6v, chis(:,:,is) )


!   Generate derivatives of dielectric tensor for this species ********

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

! ****************************************************************************

!   If we only need depsdw (for Poynting) exit routine here.
    if (only_depsdw) return

! ****************************************************************************


!   Derivatives with respect to k.
!   dbdk = d(beta)/dk; dldk = d(lambda)/dk; dxidk = d(xi)/dk


!   Derivatives of dielectric tensor eps with respect to lambda:
!   depsdl = d(eps)/d(lambda) = d(chi)/d(lambda)

    do n = nmin, nmax
       ctmp(1,n) = beta * n**2 / lambda &
          & * (eip(n)-(one+one/lambda)*ei(n)) * zf(n)
       ctmp(2,n) = ctmp(1,n) + two * beta &
          & * (ei(n)*(one-lambda)+eip(n)*(two*lambda-one)-lambda*eipp(n)) &
          & * zf(n)
       ctmp(3,n) = beta * (ei(n)-eip(n)) * xi(n)*zfp(n)

       ctmp(4,n) = zi * beta * n * (eipp(n)-two*eip(n)+ei(n)) &
          & * zf(n)
       ctmp(5,n) = iomgc * beta * sqrt(half/lambda) * n &
           & * (ei(n)*(one+half/lambda)-eip(n)) * zfp(n)
       ctmp(6,n) = iomgc * zi * beta * sqrt(half*lambda) &
          & * (eipp(n)+eip(n)*(half/lambda-two)+ei(n)*(one-half/lambda)) &
          & * zfp(n)
    end do

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

       dldx(:) = lambda &
       & * ( -two*k3*matmul(gradbunit,kvec)/k1**2 &
       &      + gradts(:)/ts - two*gradbmag/bmag )

    depsdx_xi = zero

    do n = nmin, nmax

       dxidx(:,n) = -xi(n) &
       & * ( matmul(gradbunit,kvec)/k3 + gradts(:)/(two*ts) ) &
       & - n*omgc/(k3*vth) * gradbmag/bmag

       do i = 1,6
            depsdx_xi(i,:) = depsdx_xi(i,:) + depsdxi(i,n)*dxidx(:,n)
        end do
    end do

    do ivec = 1, 3

       depsdx(:,ivec) = depsdb*dbdx(ivec) + depsdl*dldx(ivec) &
          & + depsdx_xi(:,ivec)

        call v6_Hermitian( depsdx(:,ivec), depsdx_hs(:,ivec) )

    end do

    end subroutine depsdq_bessel


! *************************************************************

    subroutine g_and_h(eq, nvec, g, h1, h3)
! Generates the quantities g(1:6), h1 and h3 as defined in the notes from
! 11-4-99.  These are elements in the chain rule for dD/dq.

    use constants_m, only : rkind, zero, one, two
    use equilibrium_m, only : eq_point
    use suscep_m, only : dielectric_general
    USE matrix3x3_m, only : hermitian3x3

    implicit none
    type(eq_point), intent(in) :: eq
    real(KIND=rkind), intent(in) :: nvec(3)
    complex(KIND=rkind), intent(out) :: g(6), h1, h3

    real(KIND=rkind) :: n1, n3
    complex(KIND=rkind) :: eps(3,3)
    complex(KIND=rkind) :: eps_h(3,3)   ! Hermitian part of eps_6v
    real(KIND=rkind) :: half = one/two

    n3 = dot_product(nvec, eq%bunit)
    n1 = sqrt( sum((nvec-n3*eq%bunit)**2) )

! N.B.  Dielectric routines take complex args, but in ray tracing nvec is real
	call dielectric_general(eq, cmplx(n1, zero, rkind), cmplx(n3, zero, rkind), eps)
    eps_h = hermitian3x3(eps)

    g(1) = n1**4 + n1**2*n3**2 + n1**2*(-eps_h(2,2) - eps_h(3,3)) - n3**2*eps_h(3,3) +&
            & eps_h(2,2)*eps_h(3,3) + eps_h(2,3)**2
    g(2) = -(n1**2*eps_h(1,1)) - n3**2*eps_h(3,3) + eps_h(1,1)*eps_h(3,3) - two*n1*n3*eps_h(1,3) -&
            & eps_h(1,3)**2
    g(3) = n1**2*n3**2 + n3**4 - n1**2*eps_h(1,1) + n3**2*(-eps_h(1,1) - eps_h(2,2)) +&
            & eps_h(1,1)*eps_h(2,2) + eps_h(1,2)**2
    g(4) = -two*n1**2*eps_h(1,2) + two*eps_h(3,3)*eps_h(1,2) + two*n1*n3*eps_h(2,3) + two*eps_h(1,3)*eps_h(2,3)
    g(5) = two*n1**3*n3 + two*n1*n3**3 - two*n1*n3*eps_h(2,2) + two*n1**2*eps_h(1,3) + &
            &2*n3**2*eps_h(1,3) - two*eps_h(2,2)*eps_h(1,3) + two*eps_h(1,2)*eps_h(2,3)
    g(6) = two*n1*n3*eps_h(1,2) + two*eps_h(1,2)*eps_h(1,3) - two*n3**2*eps_h(2,3) + two*eps_h(1,1)*eps_h(2,3)

    h1 = 4*n1**3*eps_h(1,1) + 6.0_rkind*n1**2*n3*eps_h(1,3) + two*n3**3*eps_h(1,3) + &
        & n1*(2*n3**2*(eps_h(1,1) + eps_h(3,3)) + two*(-(eps_h(1,1)*eps_h(2,2)) - eps_h(1,1)*eps_h(3,3)&
         & - eps_h(1,2)**2 + eps_h(1,3)**2)) + n3*(-two*eps_h(2,2)*eps_h(1,3) + two*eps_h(1,2)*eps_h(2,3))

    h3 = 4.0_rkind*n3**3*eps_h(3,3) + two*n1**2*n3*(eps_h(1,1) + eps_h(3,3)) + two*n1**3*eps_h(1,3) +&
        & 6.0_rkind*n1*n3**2*eps_h(1,3) + n1*(-two*eps_h(2,2)*eps_h(1,3) + two*eps_h(1,2)*eps_h(2,3)) + &
        & two*n3*(-(eps_h(1,1)*eps_h(3,3)) - eps_h(2,2)*eps_h(3,3) + eps_h(1,3)**2 - eps_h(2,3)**2)

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

    v_h = v_h/2.

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

