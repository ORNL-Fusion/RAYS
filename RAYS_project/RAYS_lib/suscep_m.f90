 module suscep_m
!   Contains routines to calculate susceptibility tensor, chi, for a single species,
!   dielectric tensor, eps, and dispersion functions D = k x k + eps
!
!   N.B. is = species number (0:nspec), and eps = I + sum(chis(is)). The dispersion
!             relation is generally written in terms of eps.
!
!   There are several forms for the dielectric tensor and dispersion functions.  For
!   ray tracing, refractive indices (n vectors, also k vectors), are all real.  But for
!   root finding the refractive indices are in general complex.  So included are routines
!   with n arguments real and also with n-parallel real and n-perp complex.  I'll deal
!   with the general case of complex n-parallel later.
!   1) A general routine, dielectric_tensor, where the susceptibility model for each
!      species is selected by setting "spec_model" in module species_m to various values
!      e.g. "bessel" or "cold".  Different suscep routines must be written for each
!      spec_model.
!
!   2) A specific dielectric_cold tensor where all species are cold
!
!   3) A specific dielectric_bessel tensor where all species are warm plasma full bessel
!
! Definitions:
!    susceptibility tensor, chi, in this module means 4*pi*i*/omega*sigma where sigma
!    is the conductivity tensor, i.e. plasma current J = sigma dot E.  Also note that total
!    chi is the sum over individual species chi(is).
!
!    dielectric tensor, eps, is I + sum[chi(is)] where I is the identity.

! N.B. Plasma quantities come in from type eq_point defined in equilibrium_m.
! An equilibrium routine must have been called previously to generate eq

! Routines exported:
!   suscep_cold(is)
!   suscep_bessel(is)
!	RLSDP_cold

!   dielectric_tensor
!   dielectric_cold
!   dielectric_bessel

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use constants_m, only : rkind

    implicit none

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

! ********************************************************************************
!       SUSCEPTIBILITY ROUTINES
! ********************************************************************************

 subroutine suscep_cold(eq, is, chi)
!   calculates the cold plasma susceptibility tensor chi for a single species species, s.

    use constants_m, only : rkind, zi=>i
    use equilibrium_m, only : eq_point

    implicit none

    type(eq_point), intent(in) :: eq
    integer, intent(in) :: is
    complex(KIND=rkind), intent(out) :: chi(3,3)

    real(KIND=rkind) :: alphas, gammas

!   alpha = (omgp/omgrf)^2, gamma = (omgc/omgrf).

    alphas= eq%alpha(is)
    gammas= eq%gamma(is)

    chi(1,1) = cmplx(-alphas / (1.-gammas**2), kind=rkind)
    chi(2,2) = cmplx(chi(1,1), kind=rkind)
    chi(3,3) = cmplx(-alphas, kind=rkind)
    chi(1,2) = cmplx(-zi*alphas*gammas / (1.-gammas**2), kind=rkind)

    chi(1,3) = cmplx(0., kind=rkind)
    chi(2,3) = cmplx(0., kind=rkind)

    chi(2,1) = -chi(1,2)
    chi(3,1) = chi(1,3)
    chi(3,2) = -chi(2,3)


    return
 end subroutine suscep_cold

! ********************************************************************************

 subroutine suscep_bessel_n1_n3_real(eq, n1, n3, is, chi)
!   Calculates the warm plasma susceptibility for a single species "is".
!   N.B. Args n1, n3 are real
!   Notations in Stix's book are used.

    use constants_m, only : zi=>i, zero, one, two
    use species_m, only : ms, n_limit, nmins, nmaxs
    use equilibrium_m, only : eq_point
    use rf_m, only : omgrf, k0
    use zfunctions_m, only : zfun, zfun0, zfun0_real_arg

    implicit none

    type(eq_point), intent(in) :: eq
    real(KIND=rkind), intent(in) :: n1
    real(KIND=rkind), intent(in) :: n3
    integer, intent(in) :: is
    complex(KIND=rkind), intent(out) :: chi(3,3)

    real(KIND=rkind) :: alphas, gammas

    real(KIND=rkind) :: a, b
    complex(KIND=rkind) :: lambda, ei(-n_limit:n_limit), eip(-n_limit:n_limit)
    real(KIND=rkind) :: xi(-n_limit:n_limit)
    complex(KIND=rkind) :: zf(-n_limit:n_limit), zfp(-n_limit:n_limit)
    complex(KIND=rkind) :: chin(6,-n_limit:n_limit)
    real(KIND=rkind) :: vth, beta, iomgc

    integer :: n, nmin, nmax

!   alpha = (omgp/omgrf)^2, gamma = (omgc/omgrf).
    alphas= eq%alpha(is)
    gammas= eq%gamma(is)
    nmin = nmins(is)
    nmax = nmaxs(is)

!   Sign of omgc:
    iomgc = sign(one,eq%omgc(is))

!   Thermal speed:
    vth = sqrt( two*eq%ts(is)/ms(is) )

!   Eq.(10-55) for lambda.  N.B. eq%omgc(is) carries the sign of the charge
    b= k0*n1*vth/eq%omgc(is)/sqrt(two)
    lambda = cmplx(b**2, zero)

!   ei = exp(-lambda)*I_n(lambda), eip = exp(-lambda)*I'_n(lambda).

       call ebessel_dbb(lambda, nmin, nmax, ei(nmin:nmax), eip(nmin:nmax))
!  write(*,*) 'suscep_bessel: vth = ', vth, '  lambda = ', lambda
!  write(*,*) 'suscep_bessel: ei = ', ei
!  write(*,*) 'suscep_bessel: eip = ', eip

    if ( abs(n3) > zero ) then
       beta = eq%omgp2(is)/(omgrf*k0*n3*vth)
!  write(*,*) 'suscep_bessel: species ', is, '   beta = ', beta

!      Z function.
       do n = nmin, nmax
         xi(n) = (omgrf-n*eq%omgc(is)) / (k0*n3*vth)
		 zf(n) = zfun0_real_arg(xi(n), n3)
		 zfp(n) = -two*(1+xi(n)*zf(n)) ! Z'
       end do
!  write(*,*) 'suscep_bessel: species ', is, '  xi = ', xi
!  write(*,*) 'suscep_bessel: species ', is, '  zf = ', zf
!  write(*,*) 'suscep_bessel: species ', is, '  zfp = ', zfp

       do n = nmin, nmax
          chin(1,n) = n**2 * (ei(n)/lambda) * zf(n)
          chin(2,n) = chin(1,n) + two*lambda*(ei(n)-eip(n)) * zf(n)
          chin(3,n) = -ei(n) * xi(n) * zfp(n)

          chin(4,n) = zi * n * (eip(n)-ei(n)) * zf(n)
          chin(5,n) = -one/sqrt(two)/b * n * ei(n) * zfp(n)
          chin(6,n) = zi/sqrt(two) * b * (eip(n)-ei(n)) * zfp(n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(1,n) = ', chin(1,n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(2,n) = ', chin(2,n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(3,n) = ', chin(3,n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(4,n) = ', chin(4,n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(5,n) = ', chin(5,n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(6,n) = ', chin(6,n)
      end do

    else ! n3 = 0
!      See Eq.(11-32). Here, an=A_n, but bn=B_n/k3.

       do n = nmin, nmax
          a = -one / (omgrf-n*eq%omgc(is))
          b = -one/two * (vth/(omgrf-n*eq%omgc(is)))**2

!         Eq.(10-57):
          chin(1,n) = n**2 * (ei(n)/lambda) * a
          chin(4,n) = zi * n * (eip(n)-ei(n)) * a
          chin(5,n) = zero

          chin(2,n) = chin(1,n) + two*lambda*(ei(n)-eip(n)) * a
          chin(6,n) = zero

          chin(3,n) = two*(omgrf-n*eq%omgc(is))/vth**2 * ei(n) * b

       end do

        beta = eq%omgp2(is)/omgrf

    end if

     chi(1,1) = beta*sum(chin(1,nmin:nmax))
     chi(2,2) = beta*sum(chin(2,nmin:nmax))
     chi(3,3) = beta*sum(chin(3,nmin:nmax))
     chi(1,2) = beta*sum(chin(4,nmin:nmax))
     chi(1,3) = beta*sum(chin(5,nmin:nmax))
     chi(2,3) = beta*sum(chin(6,nmin:nmax))

     chi(2,1) = -chi(1,2)
     chi(3,1) = chi(1,3)
     chi(3,2) = -chi(2,3)

    return
 end subroutine suscep_bessel_n1_n3_real

! ********************************************************************************

 subroutine suscep_bessel_n1_cmplx_n3_real(eq, n1, n3, is, chi)
!   Calculates the warm plasma susceptibility for a single species "is".
!   N.B. Args n1 is complex, n3 is real
!   Notations in Stix's book are used.

    use constants_m, only : zi=>i, zero, one, two
    use species_m, only : ms, n_limit, nmins, nmaxs
    use equilibrium_m, only : eq_point
    use rf_m, only : omgrf, k0
    use zfunctions_m, only : zfun, zfun0, zfun0_real_arg

    implicit none

    type(eq_point), intent(in) :: eq
    complex(KIND=rkind), intent(in) :: n1
    real(KIND=rkind), intent(in) :: n3
    integer, intent(in) :: is
    complex(KIND=rkind), intent(out) :: chi(3,3)

    real(KIND=rkind) :: alphas, gammas

    real(KIND=rkind) :: a, b
    complex(KIND=rkind) :: lambda, ei(-n_limit:n_limit), eip(-n_limit:n_limit)
    real(KIND=rkind) :: xi(-n_limit:n_limit)
    complex(KIND=rkind) :: zf(-n_limit:n_limit), zfp(-n_limit:n_limit)
    complex(KIND=rkind) :: chin(6,-n_limit:n_limit)
    real(KIND=rkind) :: vth, beta, iomgc

    integer :: n, nmin, nmax

!   alpha = (omgp/omgrf)^2, gamma = (omgc/omgrf).
    alphas= eq%alpha(is)
    gammas= eq%gamma(is)
    nmin = nmins(is)
    nmax = nmaxs(is)

!   Sign of omgc:
    iomgc = sign(one,eq%omgc(is))

!   Thermal speed:
    vth = sqrt( two*eq%ts(is)/ms(is) )

!   Eq.(10-55) for lambda.  N.B. eq%omgc(is) carries the sign of the charge
    b= k0*n1*vth/eq%omgc(is)/sqrt(two)
    lambda = b**2

!   ei = exp(-lambda)*I_n(lambda), eip = exp(-lambda)*I'_n(lambda).

       call ebessel_dbb(lambda, nmin, nmax, ei(nmin:nmax), eip(nmin:nmax))
!  write(*,*) 'suscep_bessel: vth = ', vth, '  lambda = ', lambda
!  write(*,*) 'suscep_bessel: ei = ', ei
!  write(*,*) 'suscep_bessel: eip = ', eip

    if ( abs(n3) > zero ) then
       beta = eq%omgp2(is)/(omgrf*k0*n3*vth)
!  write(*,*) 'suscep_bessel: species ', is, '   beta = ', beta

!      Z function.
       do n = nmin, nmax
         xi(n) = (omgrf-n*eq%omgc(is)) / (k0*n3*vth)
		 zf(n) = zfun0_real_arg(xi(n), n3)
		 zfp(n) = -two*(1+xi(n)*zf(n)) ! Z'
       end do
!  write(*,*) 'suscep_bessel: species ', is, '  xi = ', xi
!  write(*,*) 'suscep_bessel: species ', is, '  zf = ', zf
!  write(*,*) 'suscep_bessel: species ', is, '  zfp = ', zfp

       do n = nmin, nmax
          chin(1,n) = n**2 * (ei(n)/lambda) * zf(n)
          chin(2,n) = chin(1,n) + two*lambda*(ei(n)-eip(n)) * zf(n)
          chin(3,n) = -ei(n) * xi(n) * zfp(n)

          chin(4,n) = zi * n * (eip(n)-ei(n)) * zf(n)
          chin(5,n) = -one/sqrt(two)/b * n * ei(n) * zfp(n)
          chin(6,n) = zi/sqrt(two) * b * (eip(n)-ei(n)) * zfp(n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(1,n) = ', chin(1,n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(2,n) = ', chin(2,n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(3,n) = ', chin(3,n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(4,n) = ', chin(4,n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(5,n) = ', chin(5,n)
!   write(*,*) 'suscep_bessel: ','  n = ', n, '  chin(6,n) = ', chin(6,n)
      end do

    else ! n3 = 0
!      See Eq.(11-32). Here, an=A_n, but bn=B_n/k3.

       do n = nmin, nmax
          a = -one / (omgrf-n*eq%omgc(is))
          b = -one/two * (vth/(omgrf-n*eq%omgc(is)))**2

!         Eq.(10-57):
          chin(1,n) = n**2 * (ei(n)/lambda) * a
          chin(4,n) = zi * n * (eip(n)-ei(n)) * a
          chin(5,n) = zero

          chin(2,n) = chin(1,n) + two*lambda*(ei(n)-eip(n)) * a
          chin(6,n) = zero

          chin(3,n) = two*(omgrf-n*eq%omgc(is))/vth**2 * ei(n) * b

       end do

        beta = eq%omgp2(is)/omgrf

    end if

     chi(1,1) = beta*sum(chin(1,nmin:nmax))
     chi(2,2) = beta*sum(chin(2,nmin:nmax))
     chi(3,3) = beta*sum(chin(3,nmin:nmax))
     chi(1,2) = beta*sum(chin(4,nmin:nmax))
     chi(1,3) = beta*sum(chin(5,nmin:nmax))
     chi(2,3) = beta*sum(chin(6,nmin:nmax))

     chi(2,1) = -chi(1,2)
     chi(3,1) = chi(1,3)
     chi(3,2) = -chi(2,3)

    return
 end subroutine suscep_bessel_n1_cmplx_n3_real

! ********************************************************************************

 subroutine RLSDP_cold(eq, S ,D , P, R, L)
! Calculates S,D,P, R,L from cold plaama theory. e.g. Stix Eq 1.19 - 1.22
! N.B. Collisions not included so these are all real.

    use constants_m, only : rkind, zero, one, two
    use species_m, only : nspec, spec_model
    use equilibrium_m, only : eq_point
    use rf_m, only : omgrf

    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point), intent(in) :: eq

    real(KIND=rkind), intent(out) :: S ,D , P,  R, L
    real(KIND=rkind) :: alphas, gammas 	! alpha = (omgp/omgrf)^2, gamma = (omgc/omgrf).


    integer :: is, i

    R= zero ; L = zero; S = zero; D = zero; P = zero;
    do is = 0, nspec
		alphas= eq%alpha(is)
		gammas= eq%gamma(is)
		R = R - alphas/(one + gammas)
		L = L - alphas/(one - gammas)
		S = S - alphas/(1.-gammas**2)
		D = D - alphas*gammas/(1.-gammas**2)
		P = P - alphas
    end do

	R = one + R
	L = one + L
!	S = one + S
	S = (R + L)/two
	D = (R-L)/two
	P = one + P

    return
 end subroutine RLSDP_cold

! ********************************************************************************
!       Dielectric tensor ROUTINES
! ********************************************************************************

 subroutine dielectric_general_n1_cmplx_n3_real(eq, n1, n3, eps)

!   General form: Plasma dielectric tensor eps summed each species susceptibility, chi.
!   N.B. Args n1 is complex, n3 is real
!   Different species can have different susceptibility models.

    use constants_m, only : rkind, zero, one, two
    use species_m, only : nspec, spec_model
    use equilibrium_m, only : eq_point

    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point(nspec=nspec)), intent(in) :: eq

!   Refractive indices, n1 -> n perp, n3 -> n parallel
    complex(KIND=rkind), intent(in) :: n1
    real(KIND=rkind), intent(in) :: n3

    complex(KIND=rkind), intent(out) :: eps(3,3)
    complex(KIND=rkind) :: chi(3,3)

    integer :: is, i

    eps = zero

!   Get susceptibility tensor for each species.
    do is = 0, nspec

	  plasma_model: select case (spec_model(is) )

		  case ('cold')
			call suscep_cold(eq, is, chi)

		  case ('bessel_n1_cmplx_n3_real')
			call suscep_bessel_n1_cmplx_n3_real(eq, n1, n3, is, chi)

		  case default
			write (0,*) 'bessel_n1_cmplx_n3_real:&
			  &  unimplemented suscep model = ', spec_model(is)

	  end select plasma_model

	  eps = eps + chi

    end do

!   Dielectric tensor.

    do i =1,3
        eps(i,i) = eps(i,i) + one
    end do

    return
 end subroutine dielectric_general_n1_cmplx_n3_real

! ********************************************************************************

 subroutine dielectric_general_n1_n3_real(eq, n1, n3, eps)

!   General form: Plasma dielectric tensor eps summed each species susceptibility, chi.
!   N.B. Args n1 , n3 are real
!   Different species can have different susceptibility models.

    use constants_m, only : rkind, zero, one, two
    use species_m, only : nspec, spec_model
    use equilibrium_m, only : eq_point

    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point(nspec=nspec)), intent(in) :: eq

!   Refractive indices, n1 -> n perp, n3 -> n parallel
    real(KIND=rkind), intent(in) :: n1
    real(KIND=rkind), intent(in) :: n3

    complex(KIND=rkind), intent(out) :: eps(3,3)
    complex(KIND=rkind) :: chi(3,3)

    integer :: is, i

    eps = zero

!   Get susceptibility tensor for each species.
    do is = 0, nspec

	  plasma_model: select case (spec_model(is) )

		  case ('cold')
			call suscep_cold(eq, is, chi)

		  case ('bessel_n1_n3_real')
			call suscep_bessel_n1_n3_real(eq, n1, n3, is, chi)
		  case default
			write (0,*) 'dielectric_general_n1_n3_real: unimplemented suscep model = ',&
			           & trim(spec_model(is))

	  end select plasma_model

	  eps = eps + chi

    end do

!   Dielectric tensor.

    do i =1,3
        eps(i,i) = eps(i,i) + one
    end do

    return
 end subroutine dielectric_general_n1_n3_real

! ********************************************************************************

 subroutine dielectric_cold(eq, eps)
!   calculates the cold plasma dielectric tensor eps for each species using suscep_cold().
!   Output eps is derived type dielectric_tensor defined above.

    use constants_m, only : rkind, zero, one, two
    use species_m, only : nspec, spec_model
    use equilibrium_m, only : eq_point

    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point), intent(in) :: eq

    complex(KIND=rkind), intent(out) :: eps(3,3)

    complex(KIND=rkind) :: chi(3,3)

    integer :: is, i

    eps = zero

!   Get susceptibility tensor for each species.
    do is = 0, nspec
        call suscep_cold(eq, is, chi)
        eps = eps + chi

    end do

!   Dielectric tensor.

    do i =1,3
        eps(i,i) = eps(i,i) + one
    end do

    return
 end subroutine dielectric_cold

! ********************************************************************************

 subroutine dielectric_bessel_n1_cmplx_n3_real(eq, n1, n3, eps)
!   Calculates the dielectric tensor eps using suscep_bessel() for all species

    use constants_m, only : rkind, zero, one, two
    use species_m, only : nspec, spec_model
    use equilibrium_m, only : eq_point

    implicit none

    type(eq_point), intent(in) :: eq
    complex(KIND=rkind), intent(in) :: n1
    real(KIND=rkind), intent(in) :: n3
    complex(KIND=rkind), intent(out) :: eps(3,3)

    complex(KIND=rkind) :: chi(3,3)

    integer :: is, i

    chi = cmplx(zero, zero)
    eps = cmplx(zero, zero)

!   Get susceptibility tensor for each species.
    do is = 0, nspec
        call suscep_bessel_n1_cmplx_n3_real(eq, n1, n3, is, chi)
        eps = eps + chi

    end do

!   Dielectric tensor.

    do i =1,3
        eps(i,i) = eps(i,i) + cmplx(one, zero)
    end do

    return
 end subroutine dielectric_bessel_n1_cmplx_n3_real

!****************************************************************************

complex function disp_fun_cold_cmplx(eq, n1, n2, n3)
! calculates the cold plasma dispersion function versus the components of n perpendicular
! to B (i.e. n1, n2), and the component parallel to B (i.e. n3)

! N.B. The components of n, and the return value, disp_fun_cold_cmplx, are complex.

    use constants_m, only : rkind, one
    use equilibrium_m, only : eq_point

       implicit none

!      Derived type containing equilibrium data for a spatial point in the plasma
       type(eq_point), intent(in) :: eq

       complex(KIND=rkind), intent(in) :: n1, n2, n3

       real(KIND=rkind) :: S ,D , P,  R, L
       real(KIND=rkind) :: a, b, c
       complex(KIND=rkind) :: n_perp_sq

       call RLSDP_cold(eq, S ,D , P, R, L)

!      Coefficients for A(n3)*(n1sq)^2 + B(n3)*n1sq + C(n3) = 0.
       a = S
       b = -R*L - P*S +n3**2*(P+S)
       c = P*(n3**2 - R)*(n3**2 - L)

       n_perp_sq = n1**2 + n2**2
       disp_fun_cold_cmplx = a*n_perp_sq**2 + b*n_perp_sq + c

       return
 end function disp_fun_cold_cmplx
!****************************************************************************

real function disp_fun_cold_real(eq, n1, n2, n3)
! calculates the cold plasma dispersion function versus the components of n perpendicular
! to B (i.e. n1, n2), and the component parallel to B (i.e. n3)

! N.B. The components of n, and the return value, disp_fun_cold_real, are real.

    use constants_m, only : rkind, one
    use equilibrium_m, only : eq_point

       implicit none

!      Derived type containing equilibrium data for a spatial point in the plasma
       type(eq_point), intent(in) :: eq

       real(KIND=rkind), intent(in) :: n1, n2, n3

       real(KIND=rkind) :: S ,D , P,  R, L
       real(KIND=rkind) :: a, b, c
       real(KIND=rkind) :: n_perp_sq

       call RLSDP_cold(eq, S ,D , P, R, L)

!      Coefficients for A(n3)*(n1sq)^2 + B(n3)*n1sq + C(n3) = 0.
       a = S
       b = -R*L - P*S +n3**2*(P+S)
       c = P*(n3**2 - R)*(n3**2 - L)

       n_perp_sq = n1**2 + n2**2
       disp_fun_cold_real = a*n_perp_sq**2 + b*n_perp_sq + c

       return
 end function disp_fun_cold_real



! ********************************************************************************

 end module suscep_m
