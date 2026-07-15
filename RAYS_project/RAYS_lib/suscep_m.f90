 module suscep_m
!   Contains routines to calculate susceptibility tensor, chi, for a single species,
!   dielectric tensor, eps, and dispersion functions, disp_fun
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
!
!    dispersion function, disp_fun, is det(nn-n.nI + eps)

! N.B. Plasma quantities come in from type eq_point defined in equilibrium_m.
! An equilibrium routine must have been called previously to generate eq

! Routines exported:
!   suscep_cold(is)
!   suscep_bessel(is)
!   RLSDP_cold

!   dielectric_tensor
!   dielectric_cold
!   dielectric_bessel

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use constants_m, only : rkind, zi=>i, zero, one, two, ten, tinyr
    use equilibrium_m, only : eq_point

    implicit none

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

! ********************************************************************************
!       SUSCEPTIBILITY ROUTINES
! ********************************************************************************

 subroutine suscep_cold(eq, is, chi)
!   calculates the cold plasma susceptibility tensor chi for a single species species, s.

    implicit none

    type(eq_point), intent(in) :: eq
    integer, intent(in) :: is
    complex(KIND=rkind), intent(out) :: chi(3,3)

    real(KIND=rkind) :: alphas, gammas

!   alpha = (omgp/omgrf)^2, gamma = (omgc/omgrf).

    alphas= eq%alpha(is)
    gammas= eq%gamma(is)

    chi(1,1) = cmplx(-alphas / (one-gammas**2), kind=rkind)
    chi(2,2) = cmplx(chi(1,1), kind=rkind)
    chi(3,3) = cmplx(-alphas, kind=rkind)
    chi(1,2) = cmplx(-zi*alphas*gammas / (one-gammas**2), kind=rkind)

    chi(1,3) = zero
    chi(2,3) = zero

    chi(2,1) = -chi(1,2)
    chi(3,1) = chi(1,3)
    chi(3,2) = -chi(2,3)


    return
 end subroutine suscep_cold

! ********************************************************************************

 subroutine suscep_bessel(eq, n1, n3, is, chi)
!   Calculates the warm plasma susceptibility for a single species "is".
!   N.B. Args n1, n3 are real
!   Notations in Stix's book are used.

    use species_m, only : ms, n_limit, nmins, nmaxs
    use rf_m, only : omgrf, k0
    use zfunctions_m, only : zfun, zfun0, zfun0_real_arg

    implicit none

    type(eq_point), intent(in) :: eq
    complex(KIND=rkind), intent(in) :: n1, n3
    integer, intent(in) :: is
    complex(KIND=rkind), intent(out) :: chi(3,3)

    real(KIND=rkind) :: alphas, gammas

    real(KIND=rkind) :: an
    complex(KIND=rkind) :: b, lambda, ei(-n_limit:n_limit), eip(-n_limit:n_limit)
    complex(KIND=rkind) :: xi(-n_limit:n_limit)
    complex(KIND=rkind) :: zf(-n_limit:n_limit), zfp(-n_limit:n_limit)
    complex(KIND=rkind) :: chin(6,-n_limit:n_limit)
    real(KIND=rkind) :: vth, iomgc
    complex(KIND=rkind) :: beta

    integer :: n, nmin, nmax

   chi = zero

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
 write(*,*) ' '
 write(*,*) 'suscep_bessel: vth = ', vth, '  lambda = ', lambda
 do n = nmin, nmax
	 write(*,*) 'suscep_bessel: n = ', n,' ei = ', ei(n)
	 write(*,*) 'suscep_bessel: n = ', n,' eip = ', eip(n)
 end do

    if ( abs(n3) > zero ) then
       beta = eq%omgp2(is)/(omgrf*k0*n3*vth)
 write(*,*) 'suscep_bessel: species ', is, '   beta = ', beta

!      Z function.
       do n = nmin, nmax
         xi(n) = (omgrf-n*eq%omgc(is)) / (k0*n3*vth)
         if (abs(xi(n)%im) < ten*tinyr) then ! check for real arg
             zf(n) = zfun0_real_arg(xi(n)%re, n3%re)
         else
             zf(n) = zfun0(xi(n), n3%re)
         end if

         zfp(n) = -two*(1+xi(n)*zf(n)) ! Z'
write(*,*) 'suscep_bessel: n = ', n, ' xi(n) = ', xi(n)
write(*,*) 'suscep_bessel: zf(n) = ', zf(n)
write(*,*) 'suscep_bessel: zfp(n) = ', zfp(n)
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
  write(*,*) 'suscep_bessel: ','  n = ', n, '  beta*chin(1,n) = ', beta*chin(1,n)
  write(*,*) 'suscep_bessel: ','  n = ', n, '  beta*chin(2,n) = ', beta*chin(2,n)
  write(*,*) 'suscep_bessel: ','  n = ', n, '  beta*chin(3,n) = ', beta*chin(3,n)
  write(*,*) 'suscep_bessel: ','  n = ', n, '  beta*chin(4,n) = ', beta*chin(4,n)
  write(*,*) 'suscep_bessel: ','  n = ', n, '  beta*chin(5,n) = ', beta*chin(5,n)
  write(*,*) 'suscep_bessel: ','  n = ', n, '  beta*chin(6,n) = ', beta*chin(6,n)
      end do

    else ! n3 = 0
!      See Eq.(11-32, 11-33). Here, an=A_n, but bn=B_n/k3.

       do n = nmin, nmax
          an = -one / (omgrf-n*eq%omgc(is))

          chin(1,n) = n**2 * (ei(n)/lambda) * an
          chin(2,n) = chin(1,n) + two*lambda*(ei(n)-eip(n)) * an
          chin(3,n) = ei(n) * an
          chin(4,n) = zi * n * (eip(n)-ei(n)) * an
          chin(5,n) = zero
          chin(6,n) = zero
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
 end subroutine suscep_bessel


! ********************************************************************************

 subroutine RLSDP_cold(eq, S ,D , P, R, L)
! Calculates S,D,P, R,L from cold plaama theory. e.g. Stix Eq 1.19 - 1.22
! N.B. Collisions not included so these are all real.

    use species_m, only : nspec, spec_model
    use rf_m, only : omgrf

    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point), intent(in) :: eq

    real(KIND=rkind), intent(out) :: S ,D , P,  R, L
    real(KIND=rkind) :: alphas, gammas  ! alpha = (omgp/omgrf)^2, gamma = (omgc/omgrf).


    integer :: is, i

    R= zero ; L = zero; S = zero; D = zero; P = zero;
    do is = 0, nspec
        alphas= eq%alpha(is)
        gammas= eq%gamma(is)
        R = R - alphas/(one + gammas)
        L = L - alphas/(one - gammas)
        S = S - alphas/(one-gammas**2)
        D = D - alphas*gammas/(one-gammas**2)
        P = P - alphas
    end do

    R = one + R
    L = one + L
!   S = one + S
    S = (R + L)/two
    D = (R-L)/two
    P = one + P

    return
 end subroutine RLSDP_cold

! ********************************************************************************
!       Dielectric tensor ROUTINES
! ********************************************************************************

 subroutine dielectric_general(eq, n1, n3, eps)

!   General form: Plasma dielectric tensor eps summed each species susceptibility, chi.
!   N.B. Args n1 is complex, n3 is real
!   Different species can have different susceptibility models.

    use species_m, only : nspec, spec_model

    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point(nspec=nspec)), intent(in) :: eq

!   Refractive indices, n1 -> n perp, n3 -> n parallel
    complex(KIND=rkind), intent(in) :: n1, n3

    complex(KIND=rkind), intent(out) :: eps(3,3)
    complex(KIND=rkind) :: chi(3,3)

    integer :: is, i

    eps = zero
    chi = zero

!   Get susceptibility tensor for each species.
    do is = 0, nspec

      plasma_model: select case (spec_model(is) )

          case ('cold')
            call suscep_cold(eq, is, chi)

          case ('bessel')
            call suscep_bessel(eq, n1, n3, is, chi)

          case default
            write (0,*) 'dielectric_general:&
              &  unimplemented suscep model = ', spec_model(is)

      end select plasma_model

      eps = eps + chi

    end do

!   Dielectric tensor.

    do i =1,3
        eps(i,i) = eps(i,i) + cmplx(one, zero)
    end do

    return
    end subroutine dielectric_general

! ********************************************************************************

 subroutine dielectric_cold(eq, eps)
!   calculates the cold plasma dielectric tensor eps for each species using suscep_cold().
!   Output eps is derived type dielectric_tensor defined above.

    use constants_m, only : unitmat3
    use species_m, only : nspec, spec_model

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
        eps(i,i) = eps(i,i) + cmplx(one,zero)
    end do

    return
 end subroutine dielectric_cold

! ********************************************************************************

 subroutine dielectric_bessel(eq, n1, n3, eps)
!   Calculates the dielectric tensor eps using suscep_bessel() for all species

    use species_m, only : nspec, spec_model

    implicit none

    type(eq_point), intent(in) :: eq
    complex(KIND=rkind), intent(in) :: n1, n3
    complex(KIND=rkind), intent(out) :: eps(3,3)

    complex(KIND=rkind) :: chi(3,3)

    integer :: is, i

    chi = cmplx(zero, zero)
    eps = cmplx(zero, zero)

!   Get susceptibility tensor for each species.
    do is = 0, nspec
        call suscep_bessel(eq, n1, n3, is, chi)
        eps = eps + chi

    end do

!   Dielectric tensor.

    do i =1,3
        eps(i,i) = eps(i,i) + cmplx(one, zero)
    end do

    return
 end subroutine dielectric_bessel

!****************************************************************************

complex function disp_fun_cold(eq, n1, n2, n3)
! calculates the cold plasma dispersion function versus the components of n perpendicular
! to B (i.e. n1, n2), and the component parallel to B (i.e. n3)

! N.B. The components of n, and the return value, disp_fun_cold_cmplx, are complex.

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
       disp_fun_cold = a*n_perp_sq**2 + b*n_perp_sq + c

       return
 end function disp_fun_cold

! ********************************************************************************

 complex function disp_fun_general(eq, n1, n3)
! Calculates the general plasma dispersion function versus the components of n perpendicular
! to B (i.e. n1), and the component parallel to B (i.e. n3).
! N.B. n1 and n2 are complex.
! The dispersion function is the determinant of the dispersion tensor

    USE matrix3x3_m, only : hermitian3x3, determinant3x3

    implicit none

    ! Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point), intent(in) :: eq
    complex(KIND=rkind), intent(in) :: n1, n3

    complex(KIND=rkind) :: eps(3,3), disp_tensor(3,3)

    call dielectric_general(eq, n1, n3, eps)
    disp_tensor = hermitian3x3(eps)
    disp_tensor(1,1) = disp_tensor(1,1)-n3**2
    disp_tensor(1,3) = disp_tensor(1,3)+n1*n3
    disp_tensor(2,2) = disp_tensor(2,2)-n1**2-n3**2
    disp_tensor(3,1) = disp_tensor(3,1)+n1*n3
    disp_tensor(3,3) = disp_tensor(3,3)-n1**2

    disp_fun_general = determinant3x3(disp_tensor)

       return
 end function disp_fun_general

! ********************************************************************************

 function disp_fun_general_Hermitian(eq, n1, n3)
! Calculates the general plasma dispersion function versus the components of n perpendicular
! to B (i.e. n1), and the component parallel to B (i.e. n3).
! N.B. n1 and n2 are complex.
! The dispersion function is the determinant of the Hermitian part of dispersion tensor
! and is therefore real.  Used for root finding to initialize k for ray tracing.

    USE matrix3x3_m, only : hermitian3x3, determinant3x3

    implicit none

    real(KIND=rkind) :: disp_fun_general_Hermitian

    ! Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point), intent(in) :: eq
    complex(KIND=rkind), intent(in) :: n1, n3

    complex(KIND=rkind) :: eps(3,3), disp_tensor(3,3)

    call dielectric_general(eq, n1, n3, eps)
    disp_tensor = hermitian3x3(eps)
    disp_tensor(1,1) = disp_tensor(1,1)-n3**2
    disp_tensor(1,3) = disp_tensor(1,3)+n1*n3
    disp_tensor(2,2) = disp_tensor(2,2)-n1**2-n3**2
    disp_tensor(3,1) = disp_tensor(3,1)+n1*n3
    disp_tensor(3,3) = disp_tensor(3,3)-n1**2

    disp_fun_general_Hermitian = determinant3x3(disp_tensor)

       return
 end function disp_fun_general_Hermitian

! ********************************************************************************

 end module suscep_m
