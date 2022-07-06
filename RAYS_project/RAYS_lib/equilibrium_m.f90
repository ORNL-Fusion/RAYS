! Generic wrapper for various equilibrium routines.  Selects specific equilibrium model from
! namelist /equilibrium_list/ equilib_model.
!
! Working notes:
! 2/21/2022 (DBB) Organizational decision: The purpose of this module is to provide the 
! equilibrium data needed for ray tracing for any plasma geometry.  So I have decided to 
! eliminate any data that is geometry specific, most notably flux functions.  Specific
! equilibrium modules are required to provide all data for type eq_point, but are free
! to provide other data and services in their own module structure such as flux functions
! or boundary data that make sense for that particular geometry.  The hope is that this
! module will therefore never require modification other than addition of new equilibrium
! models to the case constructs below. 
! N.B. A consequence is that the specific geometries must be provided by modules, not 
! submodules, so that their entities are accessible outside the equilibrium_m module.
!
! 1/10/2022 (DBB) converted magnetic and species data at a spatial pint to derived 
! type -> eq_point so we can have multiple instances of equilibrium data in memory at
! the same time.

 module equilibrium_m
 
!   contains equilibrium quantities.

    use constants_m, only : rkind
    
    implicit none

!   Switch to select specific equilibrium model.
    character(len = 15) :: equilib_model

! Derived type containing equilibrium data for a spatial point in the plasma
    type eq_point
        
    !   B field. Note that bvec = B, bmag = |B|, gradbmag(i) = d(bmag)/d[x(i)],
    !   bunit = B/|B|, and gradbunit(i,j) = d[B(j)/bmag]/d[x(i)],
    !   gradbtensor(i,j) = d[B(j)]/d[x(i)].
        real(KIND=rkind) :: bvec(3), bmag, gradbmag(3), bunit(3), gradbunit(3,3), gradbtensor(3,3)

    !   Density.
        real(KIND=rkind), allocatable :: ns(:), gradns(:,:)

    !   Temperature.
        real(KIND=rkind), allocatable :: ts(:), gradts(:,:)

    !   Some often used plasma parameters.
        real(KIND=rkind), allocatable :: omgc(:), omgp2(:)
        real(KIND=rkind), allocatable :: alpha(:), gamma(:)

    !   Error returns
        character(len=20) :: equib_err = ''

    end type eq_point
          
    namelist /equilibrium_list/ equilib_model
    
!********************************************************************

contains

!********************************************************************

  subroutine initialize_equilibrium(read_input)

    use constants_m, only : input_unit    
    use diagnostics_m, only : message_unit, message, text_message
    use slab_eq_m, only : initialize_slab_eq
    use solovev_eq_m, only : initialize_solovev_eq

    implicit none
    logical, intent(in) :: read_input

    if (read_input .eqv. .true.) then    
    ! Read and write input namelist
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, equilibrium_list)
        close(unit=input_unit)
        write(message_unit, equilibrium_list)
    end if
    
    equilibria: select case (trim(equilib_model))

       case ('slab')
!         A 1-D slab equilibrium with stratification in x
          call initialize_slab_eq(read_input)

       case ('solovev')
!         A 1-D slab equilibrium with stratification in x
          call initialize_solovev_eq(read_input)


       case ('axisym_toroid')
!         A 1-D slab equilibrium with stratification in x
          call initialize_solovev_eq(read_input)

       case default
          write(0,*) 'initialize_equilibrium: improper equilib_model =', equilib_model
          call text_message('initialize_equilibrium: improper equilib_model',&
          & trim(equilib_model),0)
          stop 1

    end select equilibria


    return
  end subroutine initialize_equilibrium

!********************************************************************

 subroutine equilibrium(rvec, eq)
!   Calculates the plasma equilibrium.
!   Note that bvec = B, bmag = |B|, gradbmag(i) = d(bmag)/d[x(i)],
!   bunit = B/|B|, and gradbunit(i,j) = d[B(j)/bmag]/d[x(i)],
!   gradbtensor(i,j) = d[B(j)]/d[x(i)].
!
!   Selects from multiple equilibrium model subroutines.  As of now there is only one.

    use constants_m, only : eps0
    use rf_m, only : omgrf
    use species_m, only : nspec, ms, qs 
    use diagnostics_m, only : message
    
    use slab_eq_m, only : slab_eq
    use solovev_eq_m, only : solovev_eq
 
    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    type(eq_point), intent(out) :: eq

! Local variables
!   B field. Note that bvec = B, bmag = |B|, gradbmag(i) = d(bmag)/d[x(i)],
!   bunit = B/|B|, and gradbunit(i,j) = d[B(j)/bmag]/d[x(i)],
!   gradbtensor(i,j) = d[B(j)]/d[x(i)].
    real(KIND=rkind) :: bvec(3), bmag, gradbmag(3), bunit(3), gradbunit(3,3), gradbtensor(3,3)

!   Density.
    real(KIND=rkind) :: ns(0:nspec), gradns(3, 0:nspec)

!   Temperature.
    real(KIND=rkind)  :: ts(0:nspec), gradts(3, 0:nspec)

!   Some often used plasma parameters.
    real(KIND=rkind)  :: omgc(0:nspec), omgp2(0:nspec)
    real(KIND=rkind)  :: alpha(0:nspec), gamma(0:nspec)

    character(len=20) :: equib_err

    integer :: ivec, ivec1, ivec2

    equilibria: select case (trim(equilib_model))

       case ('slab')
!         A 1-D slab equilibrium with stratification in x
          call slab_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)

       case ('solovev')
          call solovev_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)

       case ('axisym_toroid')
          call solovev_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)

       case default
          write(0,*) 'equilibrium_m: invalid equilibrium model = ', trim(equilib_model)
          stop 1

    end select equilibria
                           
!   If equilibrium subroutine has set equib_err then return for outside handling. Do not crash.
    if (equib_err /= '') then
        eq%equib_err = equib_err        
        return
    end if

! Load up eq values
    eq%bvec = bvec
    eq%gradbtensor = gradbtensor
    eq%ns = ns
    eq%gradns = gradns
    eq%ts = ts
    eq%gradts = gradts

!   bmag = |B|, and bunit = B/|B|.
    bmag = sqrt( sum(bvec**2) )
    bunit = bvec/bmag
    eq%bmag = bmag
    eq%bunit = bunit

!   Use d[B(j)]/d[x(i)] to calculate d(bmag)/d[x(i)].
    do ivec = 1, 3
       gradbmag(ivec) = sum( gradbtensor(ivec,:)*bunit )
    end do
    eq%gradbmag = gradbmag
    
!     write(*,*) 'gradbtensor = ', gradbtensor
!     write(*,*) 'gradbmag = ', gradbmag
!     write(*,*) 'eq%gradbmag = ', eq%gradbmag

!   Use d[B(j)]/d[x(i)] and d(bmag)/d[x(i)] to calculate d[B(j)/bmag]/d[x(i)].
    do ivec1 = 1, 3; do ivec2 = 1, 3
       gradbunit(ivec1, ivec2) &
       & = ( gradbtensor(ivec1,ivec2) - gradbmag(ivec1)*bunit(ivec2) ) / bmag
    end do; end do 
    eq%gradbunit = gradbunit

!   Calculate some often used plasma parameters.

    omgc = qs(:nspec)*bmag/ms(:nspec)
    omgp2 = ns(:nspec)*qs(:nspec)**2/(eps0*ms(:nspec))
    alpha = omgp2/omgrf**2
    gamma = omgc/omgrf
    eq%omgc = omgc
    eq%omgp2 = omgp2
    eq%alpha = alpha
    eq%gamma = gamma

    return
 end subroutine equilibrium

!********************************************************************

 subroutine write_eq_point(eq, unit)
! Writes all element of eq_point type to stdout, unless optional arg 'unit' is
! present, in which case it writes to 'unit'

    use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
    
    implicit none
    
    type(eq_point), intent(in) :: eq
    integer, intent(in), optional :: unit
    integer :: out_unit 
    
    out_unit = stdout
    if (present(unit)) out_unit = unit
    
    write(out_unit, *) ' '
    write(out_unit, *) 'eq_point = '
    write(out_unit, *) 'bvec = ', eq%bvec
    write(out_unit, *) 'gradbtensor = ', eq%gradbtensor
    write(out_unit, *) 'ns = ', eq%ns
    write(out_unit, *) 'gradns = ', eq%gradns
    write(out_unit, *) 'ts = ', eq%ts
    write(out_unit, *) 'gradts = ', eq%gradts
    write(out_unit, *) 'bmag = ', eq%bmag
    write(out_unit, *) 'bunit = ', eq%bunit
    write(out_unit, *) 'gradbmag = ', eq%gradbmag
    write(out_unit, *) 'radbunit = ', eq%gradbunit
    write(out_unit, *) 'omgc = ', eq%omgc
    write(out_unit, *) 'omgp2 = ', eq%omgp2
    write(out_unit, *) 'alpha = ', eq%omgp2
    write(out_unit, *) 'gamma = ', eq%gamma
    write(out_unit, *) 'omgc = ', eq%omgc
    write(out_unit, *) 'omgp2 = ', eq%omgp2
    write(out_unit, *) 'alpha = ', eq%alpha
    write(out_unit, *) 'gamma = ', eq%gamma

 end  subroutine write_eq_point
    
!********************************************************************

 end module equilibrium_m
