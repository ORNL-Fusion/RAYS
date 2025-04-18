 module species_m
!   contains species data.

    use constants_m, only : rkind

    implicit none

!   Electron density at reference point (e.g. at magnetic axis or peak electron density)
!   The profiles generated in the various equilibrium modules are normalized to one at
!   the reference location (i.e. where ne = n0s(0) = n0).  The ion densities are specified
!   as a fraction of electron density, eta(i).
    real(KIND=rkind) :: n0

!   Maximum No. of ion species: nspec0;
!   Actual No. of ion species: nspec

    integer, parameter :: nspec0 = 5
    integer :: nspec  ! Calculated from species contained in namelist species_list

!   is = species number
!   is=0 is reserved for electrons and the rest for ions.
!   qs: charge          ms: mass
!   eta: concentration as fraction of electron density
!   n0s: number density = n0*eta
!   t0s_eV: temperature in eV.  Entered in eV in namelist file for convenience
!   t0s: temperature in MKS (Joules) converted below from t0s_eV
!   tseps_eV: temperature in eV.  Entered in eV in namelist file for convenience
!   tseps: temperature in MKS (Joules) converted below from tseps_eV
!	alfas: T_perp/T_paral
!   v0s: parallel drift velocity
!   nus: collision frequency i.e. nu/omega

! Working notes:
! (2-22-2025 DBB) As of now, none of the equilibrium models use the edge density and
! temperature variables (nseps, tseps).  Instead there are variables in the equilibrium
! modules to specify ne and Te at the edge (d_scrapeoff and T_scrapeoff) as fractions of
! the peak electron values (n0s(0), t0s(0)).

!   Static declarations and initializations

    character(len=9), dimension(0:5) :: spec_name0 = &
      & (/ 'electron ', 'hydrogen ', 'deuterium', 'tritium  ', '3He      ', 'alpha    '/)
    real(KIND=rkind), dimension(0:5) :: qs0 = (/-1., 1., 1., 1., 2., 2. /) ! units of electron charge
    real(KIND=rkind), dimension(0:5) :: ms0 = (/1., 1836., 3670., 5497., 5496., 7294. /) ! in units of me

!   Criterion for checking charge neutrality.
    real(KIND=rkind) :: neutrality = 1.e-10

!   Actual species parameters to be obtained from datain()

    character(len=12), dimension(0:nspec0) :: spec_name = ''
    real(KIND=rkind), dimension(0:nspec0) :: qs = 0.
    real(KIND=rkind), dimension(0:nspec0) :: ms = 0.
    real(KIND=rkind), dimension(0:nspec0) :: eta = 0.
    real(KIND=rkind), dimension(0:nspec0) :: n0s = 0.
    real(KIND=rkind), dimension(0:nspec0) :: nseps = 0.
    real(KIND=rkind), dimension(0:nspec0) :: t0s_eV = 0.
    real(KIND=rkind), dimension(0:nspec0) :: t0s = 0.
    real(KIND=rkind), dimension(0:nspec0) :: tseps_eV = 0.
    real(KIND=rkind), dimension(0:nspec0) :: tseps = 0.
    real(KIND=rkind), dimension(0:nspec0) :: alfas = 0.
    real(KIND=rkind), dimension(0:nspec0) :: v0s = 0.
    real(KIND=rkind), dimension(0:nspec0) :: nus = 0.

!   An array indicating which plasma model is to be used for each species
!   spec_model(is) = 'cold' susceptibility model is cold plasma
!   spec_model(is) = 'bessel' susceptibility model is full besssel function

    character(len=12) :: spec_model(0:nspec0) = ''


!   nmins, nmaxs: minimum and maximum harmonics kept in susceptibility tensor for
!   for species is.  +/- n_limit is size of nimns,nmax arrays

    integer, parameter :: n_limit = 5
    integer, dimension(-n_limit:n_limit) :: nmins, nmaxs


    namelist /species_list/ &
      & n0, nseps, spec_name, spec_model, qs, ms, t0s_eV, tseps_eV, eta, neutrality


!********************************************************************

contains

!********************************************************************

    subroutine initialize_species_m(read_input)
!   Loads charge and mass values for common plasma species from species names in namelist file
!   Called from initialize()

        use constants_m, only : e, me
        use diagnostics_m, only : message_unit, message, verbosity

        implicit none
        logical, intent(in) :: read_input

	 	integer :: input_unit, get_unit_number ! External, free unit finder
        integer :: is, j
        real(KIND=rkind) :: charge

! Read and write input namelist
        if (read_input .eqv. .true.) then
  		  	input_unit = get_unit_number()
            open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
            read(input_unit, species_list)
            close(unit=input_unit)
        end if
        if (verbosity >= 0) write(message_unit, species_list)

! Electrons:
        spec_name(0) = 'electron'
        ms(0) = 1.
        qs(0) = -1.
        eta(0) = 1.

! load up arrays with values from non-zero eta

        nspec=0
        do is = 1, nspec0
          if (eta(is) > 0.) then
            nspec = nspec+1
            do j = 1, nspec0
                if ( trim(spec_name(nspec)) == trim(spec_name0(j)) ) then
                    ms(nspec) = ms0(j)
                    qs(nspec) = qs0(j)
                end if
            end do
          end if
        end do
!   Check neutrality:
        charge = dot_product(qs(:nspec), eta(:nspec))
        if ( abs(charge) > neutrality ) then
           write(message_unit, *) 'initialize_species: charge neutrality violated, charge =',&
            &  charge
           call message('initialize_species: charge neutrality violated, charge', charge, 0)
           stop 1
        end if

!   Convert mass unit to kg:
        ms = me*ms
!   Convert charge unit from atomic number to C:
        qs = e * qs

!   Species number density relative to reference electron density:
        n0s = eta * n0

!   Convert temperature unit from eV to Joules:
        t0s = e * t0s_eV

        if (verbosity > 0) then
			write(message_unit,*) ' is      qs         ms         eta        t0s(eV)      n0s'
			do is = 0, nspec
			   write(message_unit,'(1x,i2,1p5e12.4)')                   &
			   & is, qs(is), ms(is), eta(is), t0s(is)/e, n0s(is)
        end do
        end if

    end subroutine initialize_species_m

    subroutine deallocate_species_m
        return ! nothing to deallocate
    end subroutine deallocate_species_m

 end module species_m
