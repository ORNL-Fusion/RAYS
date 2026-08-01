module uniform_eq_m
! Simple model for uniform plasma, useful for testing

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use constants_m, only : rkind, zero
    implicit none

!_________________________________________________________________________________________
! Namelist data for /uniform_eq_list/
!_________________________________________________________________________________________

! Geometry data
    ! data for bounding box (meters)
    real(KIND=rkind) :: xmin, xmax, ymin, ymax, zmin, zmax

! data for uniform magnetics
    ! Magnetic field in Tesla
    real(KIND=rkind) :: B0
    ! Direction of magnetic field.  Doesn't have to be unit vector, gets normalized
    real(KIND=rkind) :: B_direction(3)

! data for density
    ! No namelist data, n0 is set in /species_list/

! data for slab temperature
    ! No namelist data, t0s_eV is set in /species_list/

 namelist /uniform_eq_list/  B0, B_direction, xmin, xmax, ymin, ymax, zmin, zmax

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

  subroutine initialize_uniform_eq_m(read_input)

    use diagnostics_m, only : message_unit, messages_to_stdout, verbosity

    implicit none
    logical, intent(in) :: read_input
    integer :: input_unit, get_unit_number ! External, free unit finder

    if (read_input .eqv. .true.) then
        input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, uniform_eq_list)
        close(unit=input_unit)
    end if

! Write input namelist
    if (verbosity >= 0) then
        write(message_unit, uniform_eq_list)
        if (messages_to_stdout) write(*, uniform_eq_list)
    end if

    return
  end subroutine initialize_uniform_eq_m

!********************************************************************

 subroutine uniform_eq(rvec, bvec, gradbtensor, ns, gradns, ts, gradts, equib_err)
! Simple model for uniform plasma, useful for testing

    use species_m, only : nspec, n0s, t0s, eta
    use diagnostics_m, only : message_unit, message

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind), intent(out) :: ns(0:nspec), gradns(3,0:nspec)
    real(KIND=rkind), intent(out) :: ts(0:nspec), gradts(3,0:nspec)
    character(len=60), intent(out) :: equib_err

    real(KIND=rkind) :: x, y, z

    equib_err = ''
    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    gradbtensor = zero
    ns = zero
    gradns = zero
    ts = zero
    gradts = zero

! Check that we are in the box
    if (x < xmin .or. x > xmax) equib_err = 'x out_of_bounds'
    if (y < ymin .or. y > ymax) equib_err = 'y out_of_bounds'
    if (z < zmin .or. z > zmax) equib_err = 'z out_of_bounds'
    if (equib_err /= '') then
        write (message_unit, *) 'slab_eq:  equib_err ', trim(equib_err)
        return
    end if

!   Magnetic field
    bvec(:) = B0*B_direction/norm2(B_direction)

!   Density
    ns(:nspec) = n0s(:nspec)

!   Temperature
    ts(:nspec) = t0s(:nspec)


    return
 end subroutine uniform_eq

end module uniform_eq_m
