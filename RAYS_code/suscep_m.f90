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

    implicit none   

!   Susceptibility tensor for each species [complex :: chis_cold(3,3,0:nspec), chis(3,3,0:nspec)].

    complex, allocatable :: chis(:,:,:), chis_cold(:,:,:)
       
!   Dielectric tensor eps.

    complex :: eps(3,3), eps_cold(3,3)
    
    

!   expand_z0: if the Z function argument > "expand_z0", use large argument
!              expansion of Z in SUSCEP.
!   expand_z1: if the Z function argument > "expand_z1", use large argument
!              expansion of Z in DERIV_WARM.

    real :: expand_z0, expand_z1
    


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

 subroutine suscep_bessel(is)
 
!   calculates the warm plasma susceptibility for  a single species "is".
!   Notations in Stix's book are used.
!   Output, chis, is stored in module suscep_m

    use constants_m, only : zi=>i
    use equilibrium_m, only : ts, b0, omgc, omgp2
    use rf_m, only : omgrf, k1, k3
    use species_m, only : qs, ms, nmins, nmaxs, n_limit

    implicit none

    integer, intent(in) :: is

    complex :: lambda, a, b, ei(-n_limit:n_limit), eip(-n_limit:n_limit)
    complex :: zf(-n_limit:n_limit), zfp(-n_limit:n_limit)
    complex :: chin(6,-n_limit:n_limit)
    real :: vth, beta, xi(-n_limit:n_limit), iomgc

    complex, external :: zfun0


    integer :: n, nmin, nmax
    
    nmin = nmins(is)
    nmax = nmaxs(is)
    
!   Sign of omgc:
    iomgc = sign(1.,omgc(is))

!   Thermal speed:
    vth = sqrt( 2.*ts(is)/ms(is) )

!   Eq.(10-55) for lambda:
    lambda = k1**2*ts(is) / (ms(is)*omgc(is)**2)

!   ei = exp(-lambda)*I_n(lambda), eip = exp(-lambda)*I'_n(lambda).

       call ebessel_dbb(lambda, nmin, nmax, ei(nmin:nmax), eip(nmin:nmax))


    if ( k3 /= 0. ) then
!      beta = (Omga_p/Omega)^2 * [Omega/(k3*vth)].
       beta = omgp2(is)/(omgrf*k3*vth)

!      Z function.
       do n = nmin, nmax
          xi(n) = (omgrf-n*omgc(is)) / (k3*vth)

          if ( abs(xi(n)) < expand_z0 ) then
             zf(n) = zfun0(cmplx(xi(n)), k3)
             zfp(n) = -2.*(1+xi(n)*zf(n)) ! Z'
          else
             zf(n) = -1./xi(n) * (1.+.5/xi(n)**2+.75/xi(n)**4)
             zfp(n) = 1./xi(n)**2 * (1.+1.5/xi(n)**2+3.75/xi(n)**4)
          end if
       end do

       do n = nmin, nmax
          chin(1,n) = n**2 * (ei(n)/lambda) * zf(n)
          chin(2,n) = chin(1,n) + 2.*lambda*(ei(n)-eip(n)) * zf(n)
          chin(3,n) = -ei(n) * xi(n) * zfp(n)

          chin(4,n) = zi * n * (eip(n)-ei(n)) * zf(n)
          chin(5,n) = -iomgc * sqrt(.5/lambda) * n * ei(n) * zfp(n)
          chin(6,n) = iomgc * zi * sqrt(.5*lambda) * (eip(n)-ei(n)) * zfp(n)

       end do
       
    else
!      For k3=0.
!      See Eq.(11-32). Here, an=A_n, but bn=B_n/k3.
       do n = nmin, nmax
          a = -1. / (omgrf-n*omgc(is))
          b = -.5 * (vth/(omgrf-n*omgc(is)))**2
         
!         Eq.(10-57):
          chin(1,n) = n**2 * (ei(n)/lambda) * a 
          chin(4,n) = zi * n * (eip(n)-ei(n)) * a
          chin(5,n) = 0.

          chin(2,n) = chin(1,n) + 2.*lambda*(ei(n)-eip(n)) * a
          chin(6,n) = 0.

          chin(3,n) = 2.*(omgrf-n*omgc(is))/vth**2 * ei(n) * b
          
       end do
       
        beta = omgp2(is)/omgrf

    end if
     
     chis(1,1,is) = beta*sum(chin(1,nmin:nmax))
     chis(2,2,is) = beta*sum(chin(2,nmin:nmax))
     chis(3,3,is) = beta*sum(chin(3,nmin:nmax))
     chis(1,2,is) = beta*sum(chin(4,nmin:nmax))
     chis(1,3,is) = beta*sum(chin(5,nmin:nmax))
     chis(2,3,is) = beta*sum(chin(6,nmin:nmax))
     
     chis(2,1,is) = -chis(1,2,is)
     chis(3,1,is) = chis(1,3,is)
     chis(3,2,is) = -chis(2,3,is)

    return
 end subroutine suscep_bessel

! ********************************************************************************

 subroutine suscep_cold_chis(is)
 
!   calculates the cold plasma susceptibility tensor 
!   chi for a single cold species species.
!   Output, chis, is stored in module suscep_m above


    use constants_m, only : zi=>i
    use equilibrium_m, only : alpha, gamma

    implicit none

    integer, intent(in) :: is
    real :: alphas, gammas

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
    real :: alphas, gammas

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
        
      case('bessel')
        call suscep_bessel(is)
        
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
