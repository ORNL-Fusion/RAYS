 module suscep_m
!   Contains susceptibility tensor chis for each species and dielectric tensor eps
!   and routines to calculate these tensors
!
! N.B. is = species number (0:nspec), and eps = I + sum(chis(is))
!
!   There are 2 forms:
!   1) A general form chis(i,j,is), eps(i,j) where the susceptibility model for each 
!      species is selected by setting spec_model to various values e.g. "bessel"
!      or "cold"
!   2) A specific cold model chis_cold(i,j,is), eps_cold(i,j) where all species are cold
!   It is therefore possible to simultaneously have a general result and cold result
!
! N.B. Plasma quantities come in from equilibrium_m.  An equilibrium routine must have been
! called previously.

!   Routines:
!   suscep_bessel(is)
!   suscep_cold_chis(is)
!   suscep_cold(is)

!   dielectric_tensor
!   dielectric_cold

    use constants_m, only : rkind
    
    implicit none   

!   Susceptibility tensor for each species [complex(KIND=rkind) :: chis_cold(3,3,0:nspec), chis(3,3,0:nspec)].

    complex(KIND=rkind), allocatable :: chis(:,:,:), chis_cold(:,:,:)
       
!   Dielectric tensor eps.

    complex(KIND=rkind) :: eps(3,3), eps_cold(3,3)
    
! ********************************************************************************

contains

! ********************************************************************************

  subroutine initialize_suscep
 
    use species_m, only : nspec
     
    implicit none

    allocate(chis(3,3,0:nspec), chis_cold(3,3,0:nspec))

    return
  end subroutine initialize_suscep

! ********************************************************************************
!       SUSCEPTIBILITY ROUTINES
! ********************************************************************************


! ********************************************************************************

 subroutine suscep_cold_chis(is)
 
!   calculates the cold plasma susceptibility tensor 
!   chi for a single cold species species.
!   Output, chis, is stored in module suscep_m above


    use constants_m, only : zi=>i
    use equilibrium_m, only : alpha, gamma

    implicit none

    integer, intent(in) :: is
    real(KIND=rkind) :: alphas, gammas

!   alpha = (omgp/omgrf)^2, gamma = (omgc/omgrf).

    alphas= alpha(is)
    gammas= gamma(is)

    chis(1,1,is) = -alphas / (1.-gammas**2)
    chis(2,2,is) = chis(1,1,is)
    chis(3,3,is) = -alphas
    chis(1,2,is) = -zi*alphas*gammas / (1.-gammas**2)

    chis(1,3,is) = 0.
    chis(2,3,is) = 0.
     
    chis(2,1,is) = -chis(1,2,is)
    chis(3,1,is) = chis(1,3,is)
    chis(3,2,is) = -chis(2,3,is)


    return
 end subroutine suscep_cold_chis

! ********************************************************************************

 subroutine suscep_cold(is)
 
! Same as suscep_cold_chis but stores the result in chis_cold instead of in chis.  This
! is in case susceptibility with pure cold and mixed cold/warm species are both needed


    use constants_m, only : zi=>i
    use equilibrium_m, only : alpha, gamma

    implicit none

    integer, intent(in) :: is
    real(KIND=rkind) :: alphas, gammas

!   alpha = (omgp/omgrf)^2, gamma = (omgc/omgrf).

    alphas= alpha(is)
    gammas= gamma(is)

    chis_cold(1,1,is) = -alphas / (1.-gammas**2)
    chis_cold(2,2,is) = chis_cold(1,1,is)
    chis_cold(3,3,is) = -alphas
    chis_cold(1,2,is) = -zi*alphas*gammas / (1.-gammas**2)

    chis_cold(1,3,is) = 0.
    chis_cold(2,3,is) = 0.
     
    chis_cold(2,1,is) = -chis_cold(1,2,is)
    chis_cold(3,1,is) = chis_cold(1,3,is)
    chis_cold(3,2,is) = -chis_cold(2,3,is)

    return
 end subroutine suscep_cold
 
 
! ********************************************************************************
!       Dielectric tensor ROUTINES
! ********************************************************************************

 subroutine dielectric_tensor
 
!   calculates the plasma dielectric tensor eps from each
!   species' susceptibility, chi.
!   Output, "chis, and eps" are stored in module suscep_m

!   The routines suscep_... store the susceptibility tensor for
!   each species, "is", using specified plasma model, "spec_model"

    use species_m, only : nspec, spec_model

    implicit none

    integer :: is, i
    
!   Susceptibility tensor for each species.

    do is = 0, nspec
    
  plasma_model: select case (spec_model(is) )
    
      case ('cold')
        call suscep_cold_chis(is)
        
      case default
        write (0,*) 'dielectric_tensor: unimplemented species model =', spec_model(is)
        
  end select plasma_model
        
    end do

!   Dielectric tensor.

    eps = sum(chis,dim=3)
    
    do i =1,3
        eps(i,i) = eps(i,i) + 1.
    end do

    return
 end subroutine dielectric_tensor

! ********************************************************************************

 subroutine dielectric_cold
 
!   calculates the COLD plasma dielectric tensor eps_cold for each species.  This
!   is in case susceptibility with pure cold and mixed cold/warm species are both needed

!   The routine suscep_cold... stores the susceptibility tensor for
!   each species, "is", in chis_cold

    use species_m, only : nspec, spec_model

    implicit none

    integer :: is, i
    
!   Susceptibility tensor for each species.

    do is = 0, nspec
    
    call suscep_cold(is)
        
    end do


!   Dielectric tensor.

    eps_cold = sum(chis_cold,dim=3)
    
    do i =1,3
        eps_cold(i,i) = eps_cold(i,i) + 1.
    end do

    return
 end subroutine dielectric_cold
  

! ********************************************************************************

 end module suscep_m
