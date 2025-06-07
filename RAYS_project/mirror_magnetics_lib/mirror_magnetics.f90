program mirror_magnetics
! Code to generate magnetic fields for multiple mirrors aligned along the z-axis.
! The primary output is a netCDF file containing Br, Bz and Aphi evaluated on
! a uniform r,z grid.  It (optionally) prints field data for diagnostic purposes.
! Most of the actual work is done in the module mirror_magnetics_m.

! The intent is to use the field data from the netCDF output in RAYS_lib modules
! to generate spline fits for magnetic equilibria in ray tracing.  In the future
! I may add analytic derivatives so that these routines can be used directly for
! ray tracing, eliminating the spline interpolation.

! Input data comes from three namelists, in three separate files:
! "mirror_magnetics_list" contains general data for the particular case to be generated
! "coil_data_list" contains data on the fixed coil configuration e.g. coil locations
! "current_data_list" lists the current in each coil
! See "mirror_magnetics_m" for details.

! These codes link with RAYS_lib and other RAYS_project libraries

    use constants_m, only : initialize_constants_m, rkind, zero, one
    use diagnostics_m, only : initialize_diagnostics,  message, text_message, &
                       & message_unit, messages_to_stdout, verbosity

	use mirror_magnetics_m, only : initialize_mirror_magnetics, n_coils, coil, &
	                & calculate_B_on_rz_grid, n_r, n_z, r_grid, z_grid,  &
	                & Br_2D => Br, Bz_2D => Bz, Aphi_2D=> Aphi, &
	                & write_mirror_fields_Brz_NC,read_mirror_fields_Brz_NC, &
	                & NC_file_name

	implicit none
    real(KIND=rkind) :: rvec(3)
    real(KIND=rkind)  :: Bx, By, Bz, Br, Aphi
    real(KIND=rkind)  :: a, r, z, zc


    integer :: i,j
    logical :: read_input = .true.
    logical :: run_tests = .false.

!   Error returns
    character(len=60) :: equib_err = ''

	write(*,*) 'Starting mirror_magnetics'

! ****** Initialize the RAYS_lib modules, they read input data from module namelists   ****
    call initialize_diagnostics(.false.)
    call message(1)
    call initialize_constants_m
    call message(1)
	call initialize_mirror_magnetics(read_input)

! ****** Test subroutine coil_rz_field_unit, no sum over coils  ****

 tests: if(run_tests) then
	write(*,*)
	write(*,*) '****** Test subroutine coil_rz_field_unit, no sum over coils  ****'
	write(message_unit,*) '****** Test subroutine coil_rz_field_unit, no sum over coils  ****'

	rvec = (/zero,zero,zero/)
	r = sqrt(rvec(1)**2 + rvec(2)**2)
	z = rvec(3)
	write(*,*) ' '
	write(*,*) 'Evaluate field from coil "i" at location rvec = ', rvec
	do i = 1, n_coils
		zc =  coil(i)%z_center
		write(*,*) 'coil number = ', i, 'z_center =  ', zc
		write(message_unit,*) 'coil number = ', i, 'z_center =  ', zc

		call coil(i)%coil_Brz_field_1Amp(r, z, Br, Bz, Aphi, equib_err)
		write(*,*) 'r = ', r, '  z-zc = ', z-zc,'  Br = ', Br, '  Bz = ', Bz, &
				   & ' Aphi = ', Aphi
		write(message_unit,*) 'r = ', r, '  z-zc = ', z-zc,'  Br = ', Br, '  Bz = ', Bz, &
				   & ' Aphi = ', Aphi
	end do

	rvec = (/zero, zero, one/)
	r = sqrt(rvec(1)**2 + rvec(2)**2)
	z = rvec(3)
	write(*,*) ' '
	write(*,*) 'Evaluate field from coil "i" at location rvec = ', rvec
	do i = 1, n_coils
		zc =  coil(i)%z_center
		write(*,*) 'coil number = ', i, 'z_center =  ', zc
		write(message_unit,*) 'coil number = ', i, 'z_center =  ', zc

		call coil(i)%coil_Brz_field_1Amp(r, z, Br, Bz, Aphi, equib_err)
		write(*,*) 'r = ', r, '  z-zc = ', z-zc,'  Br = ', Br, '  Bz = ', Bz, &
				   & ' Aphi = ', Aphi
		write(message_unit,*) 'r = ', r, '  z-zc = ', z-zc,'  Br = ', Br, '  Bz = ', Bz, &
				   & ' Aphi = ', Aphi
	end do

	rvec = (/zero, 0.5*one, one/)
	r = sqrt(rvec(1)**2 + rvec(2)**2)
	z = rvec(3)
	write(*,*) ' '
	write(*,*) 'Evaluate field from coil "i" at location rvec = ', rvec
	do i = 1, n_coils
		zc =  coil(i)%z_center
		write(*,*) 'coil number = ', i, 'z_center =  ', zc
		write(message_unit,*) 'coil number = ', i, 'z_center =  ', zc

		call coil(i)%coil_Brz_field_1Amp(r, z, Br, Bz, Aphi, equib_err)
		write(*,*) 'r = ', r, '  z-zc = ', z-zc,'  Br = ', Br, '  Bz = ', Bz, &
				   & ' Aphi = ', Aphi
		write(message_unit,*) 'r = ', r, '  z-zc = ', z-zc,'  Br = ', Br, '  Bz = ', Bz, &
				   & ' Aphi = ', Aphi
	end do
 end if tests

! ****** Calculate mirror_Brz_field on r/z grid, sum over coils  ****

	write(*,*)
	write(*,*)
	write(*,*) 'Calculate mirror_Brz_field on r/z grid, sum over coils  ****'
	write(message_unit,*) 'Calculate mirror_Brz_field on r/z grid, sum over coils  ****'

	call calculate_B_on_rz_grid
! 	write(*,*) 'n_r = ', n_r
! 	do i = 1, n_r
! 		write(*,*)   ' '
! 		write(*,*)  " r = ", r_grid(i)
! 		do j = 1, n_z
! 			write(*,*) '  z = ', z_grid(j), '  Br = ', Br_2D(2,j), '  Bz = ', Bz_2D(2,j),&
! 			         & '  Aphi = ', Aphi_2D(2,j)
! 		end do
! 	end do

	call write_mirror_fields_Brz_NC

	call read_mirror_fields_Brz_NC(NC_file_name)



end program mirror_magnetics


