! Generic wrapper for various equilibrium routines.  Selects specific equilibrium model from
! namelist /equilibrium_list/ equilib_model.
!
! Working notes:
!
! 10/3/2023 (DBB) Changed plasma species data to be dimensioned 0:nspec0 rather than
! allocated.  The implicit allocations in loading data into the eq_point type may be
! slowing down parallel execution.
!
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
    use species_m, only : nspec0

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
		real(KIND=rkind) :: ns(0:nspec0), gradns(3, 0:nspec0)

	!   Temperature.
		real(KIND=rkind)  :: ts(0:nspec0), gradts(3, 0:nspec0)

	!   Some often used plasma parameters.
		real(KIND=rkind)  :: omgc(0:nspec0), omgp2(0:nspec0)
		real(KIND=rkind)  :: alpha(0:nspec0), gamma(0:nspec0)

    !   Error returns
        character(len=60) :: equib_err = ''

    end type eq_point

    namelist /equilibrium_list/ equilib_model

!********************************************************************

contains

!********************************************************************

  subroutine initialize_equilibrium_m(read_input)

    use diagnostics_m, only : message_unit, message, text_message, messages_to_stdout, verbosity
    use species_m, only : nspec0
    use slab_eq_m, only : initialize_slab_eq_m
    use solovev_eq_m, only : initialize_solovev_eq_m
    use axisym_toroid_eq_m, only : initialize_axisym_toroid_eq_m

    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder

    if (read_input .eqv. .true.) then
    ! Read and write input namelist
    	input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, equilibrium_list)
        close(unit=input_unit)
    end if

! Write input namelist
    if (verbosity >= 0) then
		write(message_unit, equilibrium_list)
		if (messages_to_stdout) write(*, equilibrium_list)
    end if

    equilibria: select case (trim(equilib_model))

       case ('slab')
!         A 1-D slab equilibrium with stratification in x
          call initialize_slab_eq_m(read_input)

       case ('solovev')
!         A simple analytic tokamak model
          call initialize_solovev_eq_m(read_input)


       case ('axisym_toroid')

!         A generic axisymmetric toroidal plasma model
          call initialize_axisym_toroid_eq_m(read_input)

       case default
          write(0,*) 'initialize_equilibrium: improper equilib_model =', equilib_model
          call text_message('initialize_equilibrium: improper equilib_model',&
          & trim(equilib_model),0)
          stop 1

    end select equilibria

    return
  end subroutine initialize_equilibrium_m

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
    use diagnostics_m, only : message, text_message

    use slab_eq_m, only : slab_eq
    use solovev_eq_m, only : solovev_eq
    use axisym_toroid_eq_m, only : axisym_toroid_eq

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

    character(len=60) :: equib_err

    integer :: ivec, ivec1, ivec2
    equilibria: select case (trim(equilib_model))
       case ('slab')
!         A 1-D slab equilibrium with stratification in x
          call slab_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)

       case ('solovev')
          call solovev_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)

       case ('axisym_toroid')
          call axisym_toroid_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)

       case default
          write(0,*) 'equilibrium_m: invalid equilibrium model = ', trim(equilib_model)
          stop 1

    end select equilibria

!   If equilibrium subroutine has set equib_err then return for outside handling. Do not crash.
    if (equib_err /= '') then
        eq%equib_err = equib_err
        call text_message('equilibrium:  equib_err = ', equib_err,1)
        return
    end if

! Initialize eq values
    eq%bvec = 0.
    eq%gradbtensor = 0.
    eq%ns = 0.
    eq%gradns = 0.
    eq%ts = 0.
    eq%gradts = 0.

!   bmag = |B|, and bunit = B/|B|.
    bmag = 0.
    bunit = 0.
    eq%bmag = 0.
    eq%bunit = 0.
    eq%gradbmag = 0.
    eq%gradbunit = 0.

    omgc = 0.
    omgp2 = 0.
    alpha = 0.
    gamma = 0.
    eq%omgc = 0.
    eq%omgp2 = 0.
    eq%alpha = 0.
    eq%gamma = 0.

! Load up eq values
    eq%bvec = bvec
    eq%gradbtensor = gradbtensor
    eq%ns(0:nspec) = ns(0:nspec)
    eq%gradns(:,0:nspec) = gradns(:,0:nspec)
    eq%ts(0:nspec) = ts(0:nspec)
    eq%gradts(:,0:nspec) = gradts(:,0:nspec)

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

    omgc(:nspec) = qs(:nspec)*bmag/ms(:nspec)
    omgp2(:nspec) = ns(:nspec)*qs(:nspec)**2/(eps0*ms(:nspec))
    alpha(:nspec) = omgp2(:nspec)/omgrf**2
    gamma(:nspec) = omgc(:nspec)/omgrf
    eq%omgc(:nspec) = omgc(:nspec)
    eq%omgp2(:nspec) = omgp2(:nspec)
    eq%alpha(:nspec) = alpha(:nspec)
    eq%gamma(:nspec) = gamma(:nspec)

    return
 end subroutine equilibrium

!********************************************************************

 subroutine write_eq_point(eq, unit)
! Writes all element of eq_point type to stdout, unless optional arg 'unit' is
! present, in which case it writes to 'unit'

    use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
    use diagnostics_m, only : message, message_unit, text_message, verbosity
    use species_m, only : nspec

    implicit none

    type(eq_point), intent(in) :: eq
    integer, intent(in), optional :: unit
    integer :: write_unit

    write_unit = 6
    if (present(unit)) write_unit = unit

    write(write_unit, *) ' '
    write(write_unit, *) 'eq_point = '
    write(write_unit, *) 'bvec = ', eq%bvec
    write(write_unit, *) 'gradbtensor = ', eq%gradbtensor
    write(write_unit, *) 'ns = ', eq%ns
    write(write_unit, *) 'gradns = ', eq%gradns
    write(write_unit, *) 'ts = ', eq%ts
    write(write_unit, *) 'gradts = ', eq%gradts
    write(write_unit, *) 'bmag = ', eq%bmag
    write(write_unit, *) 'bunit = ', eq%bunit
    write(write_unit, *) 'gradbmag = ', eq%gradbmag
    write(write_unit, *) 'gradbunit = ', eq%gradbunit
    write(write_unit, *) 'omgc = ', eq%omgc
    write(write_unit, *) 'omgp2 = ', eq%omgp2
    write(write_unit, *) 'alpha = ', eq%omgp2
    write(write_unit, *) 'gamma = ', eq%gamma

!     call message(1)
!     call text_message('eq_point')
!     call message('bvec = ', eq%bvec)
!     call message('gradbtensor = ', eq%gradbtensor, 3, 3)
!     call message('ns = ', eq%ns)
!     call message('gradns = ', eq%gradns, 3, nspec+1)
!     call message('ts = ', eq%ts)
!     call message('gradts = ', eq%gradts, 3, nspec+1)
!     call message('bmag = ', eq%bmag)
!     call message('bunit = ', eq%bunit)
!     call message('gradbmag = ', eq%gradbmag)
!     call message('gradbunit = ', eq%gradbunit, 3, 3)
!     call message('omgc = ', eq%omgc)
!     call message('omgp2 = ', eq%omgp2)
!     call message('alpha = ', eq%alpha)
!     call message('gamma = ', eq%gamma)

 end  subroutine write_eq_point

!********************************************************************

    subroutine deallocate_equilibrium_m
		use slab_eq_m, only : deallocate_slab_eq_m
		use solovev_eq_m, only : deallocate_solovev_eq_m
		use axisym_toroid_eq_m, only : deallocate_axisym_toroid_eq_m

		call deallocate_slab_eq_m
		call deallocate_solovev_eq_m
		call deallocate_axisym_toroid_eq_m
    end subroutine deallocate_equilibrium_m

 end module equilibrium_m
