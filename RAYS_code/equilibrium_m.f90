 module equilibrium_m
 
!   contains equilibrium quantities.

    implicit none

!   Switch to select equilibrium.
    character(len = 15) :: equilib_model
    
!   B field at reference point (e.g. magnetic axis)    
    real :: b0

!   B field. Note that bvec = B, bmag = |B|, gradbmag(i) = d(bmag)/d[x(i)],
!   bunit = B/|B|, and gradbunit(i,j) = d[B(j)/bmag]/d[x(i)],
!   gradbtensor(i,j) = d[B(j)]/d[x(i)].
    real :: bvec(3), bmag, gradbmag(3), bunit(3), gradbunit(3,3), gradbtensor(3,3)
    
! Flux function psi
    real :: psi, gradpsi(3)

!   Density.
    real, allocatable :: ns(:), gradns(:,:)

!   Temperature.
    real, allocatable :: ts(:), gradts(:,:)

!   Some often used plasma parameters.
    real, allocatable :: omgc(:), omgp2(:)
    real, allocatable :: alpha(:), gamma(:)
      
    namelist /equilibrium_list/ equilib_model, b0
    
!********************************************************************

contains

!********************************************************************

  subroutine initialize_equilibrium

    use constants_m, only : input_unit    
    use diagnostics_m, only : message_unit, message, text_message
    use species_m, only : nspec
    use slab_eq_m, only : initialize_slab_eq

    implicit none

! Allocate plasma quantities

    allocate(ns(0:nspec), gradns(3,0:nspec))
    allocate(ts(0:nspec), gradts(3,0:nspec))
    allocate(omgc(0:nspec), omgp2(0:nspec))
    allocate(alpha(0:nspec), gamma(0:nspec))

! Read and write input namelist
    open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
    read(input_unit, equilibrium_list)
    close(unit=input_unit)
    write(message_unit, equilibrium_list)
    
    equilibria: select case (trim(equilib_model))

       case ('slab')
!         A 1-D slab equilibrium.
          call initialize_slab_eq

       case default
          write(0,*) 'initialize_equilibrium: improper equilib_model =', equilib_model
          call text_message('initialize_equilibrium: improper equilib_model',&
          & trim(equilib_model),0)
          stop 1

    end select equilibria


    return
  end subroutine initialize_equilibrium

!********************************************************************

 subroutine equilibrium(rvec)
!   Calculates the plasma equilibrium.
!   Note that bvec = B, bmag = |B|, gradbmag(i) = d(bmag)/d[x(i)],
!   bunit = B/|B|, and gradbunit(i,j) = d[B(j)/bmag]/d[x(i)],
!   gradbtensor(i,j) = d[B(j)]/d[x(i)].
!
!   Selects from multiple equilibrium model subroutines.  As of now there is only one.

    use constants_m, only : eps0
    use rf_m, only : omgrf
    use species_m, only : nspec, ms, qs 
    use diagnostics_m, only : message, equib_err
    
    use slab_eq_m, only : slab_eq
 
    implicit none

    real, intent(in) :: rvec(3)

    integer :: ivec, ivec1, ivec2

    equilibria: select case (trim(equilib_model))

       case ('slab')
!         A 1-D slab equilibrium.
          call slab_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts)

       case default
          continue

    end select equilibria
                           
!   If equilibrium subroutine has set equib_err return for outside handling. Does not crash.
    if (equib_err /= '') return

!   bmag = |B|, and bunit = B/|B|.
    bmag = sqrt( sum(bvec**2) )
    bunit = bvec / bmag

!   Use d[B(j)]/d[x(i)] to calculate d(bmag)/d[x(i)].
    do ivec = 1, 3
       gradbmag(ivec) = sum( gradbtensor(ivec,:)*bunit )
    end do

!   Use d[B(j)]/d[x(i)] and d(bmag)/d[x(i)] to calculate d[B(j)/bmag]/d[x(i)].
    do ivec1 = 1, 3; do ivec2 = 1, 3
       gradbunit(ivec1, ivec2) &
       & = ( gradbtensor(ivec1,ivec2) - gradbmag(ivec1)*bunit(ivec2) ) / bmag
    end do; end do 

!   Calculate some often used plasma parameters.

    omgc = qs(:nspec)*bmag/ms(:nspec)
    omgp2 = ns(:nspec)*qs(:nspec)**2/(eps0*ms(:nspec))
    alpha = omgp2/omgrf**2
    gamma = omgc/omgrf

    return
 end subroutine equilibrium


 end module equilibrium_m
