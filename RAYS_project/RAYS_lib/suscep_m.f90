 module suscep_m
!   Contains routines to calculate susceptibility tensor, chi, for a single species and
!   dielectric tensor, eps. eps is of derived type dielectric_tensor.
!
!   N.B. is = species number (0:nspec), and eps = I + sum(chis(is))
!
!   There are 2 forms for dielectric tensor:
!   1) A general routine, dielectric_tensor, where the susceptibility model for each
!      species is selected by setting "spec_model" in module species_m to various values
!      e.g. "bessel" or "cold".  different suscep routines must be written for each spec_model
!      As of right now only "suscep_cold" is implemented.
!
!   2) A specific cold dielectric tensor where all species are cold
!
! N.B. Plasma quantities come in from equilibrium_m.  An equilibrium routine must have been
! called previously.

!   Routines:
!   suscep_cold(is)
!   suscep_bessel(is) -> To be implemented now in the spare parts box

!   dielectric_tensor
!   dielectric_cold

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
!       Dielectric tensor ROUTINES
! ********************************************************************************

 subroutine dielectric(eq, eps)

!   calculates the plasma dielectric tensor eps from each species susceptibility, chi.
!   Output eps is derived type dielectric_tensor defined above.  Different species can have
!   different susceptibility models.

    use species_m, only : nspec, spec_model
    use equilibrium_m, only : eq_point

    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point(nspec=nspec)), intent(in) :: eq

    complex(KIND=rkind), intent(out) :: eps(3,3)
    complex(KIND=rkind) :: chi(3,3)

    integer :: is, i

    eps = 0.0_rkind

!   Get susceptibility tensor for each species.
    do is = 0, nspec

          plasma_model: select case (spec_model(is) )

              case ('cold')
                call suscep_cold(eq, is, chi)

              case default
                write (0,*) 'dielectric_tensor: unimplemented species model =', spec_model(is)

          end select plasma_model

          eps = eps + chi

    end do

!   Dielectric tensor.

    do i =1,3
        eps(i,i) = eps(i,i) + 1.0_rkind
    end do

    return
 end subroutine dielectric

! ********************************************************************************

 subroutine dielectric_cold(eq, eps)
!   calculates the cold plasma dielectric tensor eps for each species using suscep_cold().
!   Output eps is derived type dielectric_tensor defined above.

    use species_m, only : nspec, spec_model
    use equilibrium_m, only : eq_point

    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point), intent(in) :: eq

    complex(KIND=rkind), intent(out) :: eps(3,3)

    complex(KIND=rkind) :: chi(3,3)

    integer :: is, i

    eps = 0.0_rkind

!   Get susceptibility tensor for each species.
    do is = 0, nspec
        call suscep_cold(eq, is, chi)
        eps = eps + chi

    end do

!   Dielectric tensor.

    do i =1,3
        eps(i,i) = eps(i,i) + 1.0_rkind
    end do

    return
 end subroutine dielectric_cold

! ********************************************************************************

 subroutine RLSDP_cold(eq, S ,D , P, R, L)
! Calculates S,D,P, R,L from cold plaama theory. e.g. Stix Eq 1.19 - 1.22
! N.B. Collisions not included so these are all real.

    use constants_m, only : rkind, zero, one
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
	S = (R + L)/2._rkind
	D = (R-L)/2._rkind
	P = one + P

    return
 end subroutine RLSDP_cold

! ********************************************************************************

 end module suscep_m
