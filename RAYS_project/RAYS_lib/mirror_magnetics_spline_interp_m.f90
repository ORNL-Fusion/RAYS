module  mirror_magnetics_spline_interp_m

! Calculates standard magnetic quantities for multiple magnetic mirrors aligned on the
! z axis.  It reads a netCDF file containing Br,Bz, Aphi on an r,z grid and exports 2D
! spline function derived types (Br_spline, Br_spline, Aphi_spline ).  The netCDF file is
! written by stand-alone program mirror_magnetics.f90 which resides in mirror_magnetics_lib.
!
! A flux-function-like quantity, Aphi, is provided to serve as a radial coordinate that is
! constant along field lines.  This is normalized to be unity on the last un-interrupted
! flux surface (or field line), LUFS.
!
! N.B. This module uses the mirror_magnetics_m module in mirror_magnetics_lib to read the
!      netCDF file and gets the data for initialization from that module.


! Working notes:
!_________________________________________________________________________________________

    use constants_m, only : rkind, zero, one, two
    use quick_cube_splines_m, only : cube_spline_function_1D, cube_spline_function_2D

    implicit none

! Local data **************************************************

! Introduce r_LUFS_spline, z_LUFS_spline, Aphi_LUFS_spline as local module variables so
! don't have to use multiple_mirror_eq_m in subroutines below  -> avoid circularity

    real(KIND=rkind) ::  r_LUFS_spline, z_LUFS_spline, Aphi_LUFS_spline

! Stuff for 2D spline profiles

    type(cube_spline_function_2D) :: Br_spline
 	character (len = 80) :: Br_spline_name = 'Br_spline'
    type(cube_spline_function_2D) :: Bz_spline
 	character (len = 80) :: Bz_spline_name = 'Bz_spline'
    type(cube_spline_function_2D) :: Aphi_spline
 	character (len = 80) :: Aphi_spline_name = 'Aphi_spline'

! Namelist data for /mirror_magnetics_spline_interp_list/  *****************************

    character (len = 100) :: mirror_field_NC_file ! Input in multiple_mirror_eq_m namelist

	namelist / mirror_magnetics_spline_interp_list/ mirror_field_NC_file

!_________________________________________________________________________________________

contains
!_________________________________________________________________________________________

  subroutine initialize_mirror_magnetics_spline_interp(read_input, &
               & box_rmax, box_zmin, box_zmax)


    use constants_m, only : one
    use species_m, only : nspec
    use diagnostics_m, only : message, text_message, message_unit,messages_to_stdout,&
                            & verbosity
    use mirror_magnetics_m, only : n_r, n_z, r_min, r_max, z_min, z_max, &
                                  & r_grid, z_grid, Br, Bz, Aphi, &
                                  & r_LUFS, z_LUFS, &
                                   & read_mirror_fields_Brz_NC
! N.B. For generality module mirror_magnetics_m allows for a non-zero r_min, although for
!      normal application to magnetic mirror devices with coils centered on the z axis
!      r_min should always be 0.  Here we check for that and crash if it is non-zero.

   implicit none

    logical, intent(in) :: read_input
 	integer :: input_unit, get_unit_number ! External, free unit finder

! Geometry data
    ! data for bounding box of computational domain
    real(KIND=rkind), intent(out) :: box_rmax, box_zmin, box_zmax

    if (verbosity >= 0) then
		write(*,*) 'initialize_mirror_magnetics_spline_interp'
    end if
! Read input namelist
    if (read_input .eqv. .true.) then
    	input_unit = get_unit_number()
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, mirror_magnetics_spline_interp_list)
        close(unit=input_unit)
    end if

! Write input namelist
    if (verbosity > 0) then
		write(message_unit, mirror_magnetics_spline_interp_list)
		if (messages_to_stdout) write(*, mirror_magnetics_spline_interp_list)
    end if

! Allocate the grids and field array and load the arrays
    call read_mirror_fields_Brz_NC(trim(mirror_field_NC_file))

    if (r_min /= zero) then
    	call text_message('initialize_mirror_magnetics_spline_interp: non-zero r_min')
    	write (*,*) 'initialize_mirror_magnetics_spline_interp: non-zero r_min'
    	stop
    end if

    box_rmax = r_max
    box_zmin = z_min
    box_zmax = z_max

    r_LUFS_spline = r_LUFS
    z_LUFS_spline = z_LUFS

! Initialize spline coefficients for Br_spline on 2D grid
	call Br_spline%cube_spline_2D_init(n_r, r_grid, n_z, z_grid, Br, Br_spline_name)

! Initialize spline coefficients for Bz_spline on 2D grid
	call Bz_spline%cube_spline_2D_init(n_r, r_grid, n_z, z_grid, Bz, Bz_spline_name)

! Initialize spline coefficients for Aphi_spline on 2D grid
	call Aphi_spline%cube_spline_2D_init(n_r, r_grid, n_z, z_grid, Aphi, Aphi_spline_name)

! Calculate Aphi at strike point of LUFS for normalization purposes.
	call Aphi_spline%eval_2D_f(r_LUFS_spline, z_LUFS_spline, Aphi_LUFS_spline)
	write(*,*) 'r_LUFS_spline = ', r_LUFS_spline,'  z_LUFS_spline = ', z_LUFS_spline,'  Aphi_LUFS_spline = ', Aphi_LUFS_spline

    return
  end subroutine initialize_mirror_magnetics_spline_interp

!****************************************************************************************

  subroutine  mirror_magnetics_spline_interp(rvec, bvec, gradbtensor, Aphi, gradAphi,&
                           AphiN, gradAphiN, equib_err)

!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

    use diagnostics_m, only : message_unit, message

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: bvec(3), gradbtensor(3,3)
    real(KIND=rkind), intent(out) :: Aphi, gradAphi(3), AphiN, gradAphiN(3)
    character(len=60), intent(out) :: equib_err

    real(KIND=rkind) :: x, y, z, r
    real(KIND=rkind) :: br, bz, bphi, bp0
    real(KIND=rkind) :: dd_rho, dbrdr, dbrdz, dbzdr, dbzdz, dbphidr
    real(KIND=rkind) :: dAphidr, dAphidz

    equib_err = ''
    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    r = sqrt(x**2+y**2)

	call Br_spline%eval_2D_fp(r, z, br, dbrdr, dbrdz)
	call Bz_spline%eval_2D_fp(r, z, bz, dbzdr, dbzdz)
	call Aphi_spline%eval_2D_fp(r, z, Aphi, dAphidr, dAphidz)

    if (r < two*tiny(one)) then ! On axis
		bvec = (/ zero, zero, bz /)

		gradbtensor = zero
		gradbtensor(1,1) = -dbzdz/two
		gradbtensor(2,2) = -dbzdz/two
		gradbtensor(3,3) = dbzdz
		Aphi = zero
		gradAphi = zero
		AphiN = zero
		gradAphiN = zero

	else ! Off axis
		bvec = (/ x*br/r, y*br/r, bz /)

!   d(Bx)/dx, d(Bx)/dy, and d(Bx)/dz:
		gradbtensor(1,1) = (one - (x/r)**2)*br/r + (x/r)**2*dbrdr
		gradbtensor(2,1) =  x*y/r**2*(dbrdr - br/r)
		gradbtensor(3,1) = dbrdz*x/r

!   d(By)/dx, d(By)/dy, and d(By)/dz:
		gradbtensor(1,2) = x*y/r**2*(dbrdr - br/r)
		gradbtensor(2,2) = (one - (y/r)**2)*br/r + (y/r)**2*dbrdr
		gradbtensor(3,2) = dbrdz*y/r

!   d(Bz)/dx, d(Bz)/dy, and d(Bz)/dz:
		gradbtensor(1,3) = dbzdr * x/r
		gradbtensor(2,3) = dbzdr * y/r
		gradbtensor(3,3) = dbzdz

! Aphi
		gradAphi(1) = dAphidr*x/r
		gradAphi(2) = dAphidr*y/r
		gradAphi(3) = dAphidz

    end if

!   Normalized Flux function x, y, z normalized to 1.0 at last un-interupted flux surface
    AphiN = Aphi/Aphi_LUFS_spline
    gradAphiN = gradAphi/Aphi_LUFS_spline

    return
    end subroutine  mirror_magnetics_spline_interp
!****************************************************************************************

  subroutine  mirror_magnetics_spline_interp_Aphi(rvec, Aphi, gradAphi, AphiN, gradAphiN)

!   Checks for some error conditions and sets equib_err for outside handling.  Does not
!   stop.

    use diagnostics_m, only : message_unit, message

    implicit none

    real(KIND=rkind), intent(in) :: rvec(3)
    real(KIND=rkind), intent(out) :: Aphi, gradAphi(3), AphiN, gradAphiN(3)

    real(KIND=rkind) :: x, y, z, r
    real(KIND=rkind) :: br, bz, bphi, bp0
    real(KIND=rkind) :: dd_rho, dbrdr, dbrdz, dbzdr, dbzdz, dbphidr
    real(KIND=rkind) :: dAphidr, dAphidz

    x = rvec(1)
    y = rvec(2)
    z = rvec(3)
    r = sqrt(x**2+y**2)

	call Aphi_spline%eval_2D_fp(r, z, Aphi, dAphidr, dAphidz)

    if (r < two*tiny(one)) then ! On axis
		Aphi = zero
		gradAphi = zero
		AphiN = zero
		gradAphiN = zero

	else ! Off axis
		gradAphi(1) = dAphidr*x/r
		gradAphi(2) = dAphidr*y/r
		gradAphi(3) = dAphidz

    end if

!   Normalized Flux function x, y, z normalized to 1.0 at last un-interupted flux surface
    AphiN = Aphi/Aphi_LUFS_spline
    gradAphiN = gradAphi/Aphi_LUFS_spline

    return
    end subroutine  mirror_magnetics_spline_interp_Aphi

!********************************************************************

    subroutine deallocate_mirror_magnetics_spline_interp_m
		! Nothing to deallocate
		return
    end subroutine deallocate_mirror_magnetics_spline_interp_m

 !********************************************************************

end module  mirror_magnetics_spline_interp_m
