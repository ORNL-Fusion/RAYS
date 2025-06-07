module mirror_magnetics_m
! Routines and derived types to generate magnetic fields for multiple miorrors aligned
! along the z-axis.  The present intent is to be able to generate the fields, then to use
! these to do spline fits in various r,z domains.  The spline fits and the spline
! derivatives would then be used in RAYS rather than using these fields directly.  We'll
! see how it turns out.

! We define a derived type "coil" which has the data defining the location of conductors
! in the coil.  Inner and outer coil radius and coil z center are self explanatory. There
! can be thousands of turns in the coil, so we discretize the arrangement of turns into
! a small number of filaments arranged in a (n_r_layers X n_z_slices) grid such that
! n_turns = n_r_layers * n_z_slices.

! A type-bound procedure 'coil_rz_field_1Amp' is provided which gives (Br, Bz, Aphi)
! per amp of current as functions of (r,z).

! Input data comes from three namelists, in three separate files:
! "mirror_magnetics_list" contains general data for the particular case to be generated
! "coil_data_list" contains data on the fixed coil configuration e.g. coil locations
! "current_data_list" lists the current in each coil

!    use constants_m, only : rkind, zero

    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind),parameter :: zero = 0.0_rkind

    character(len=80) :: coil_data_file = ''
    character(len=80) :: current_data_file = ''
    character(len=80) :: case_name = ''
    character(len=80) :: base_file_name = ''
    character(len=100) :: NC_file_name = ''
	integer :: n_coils

! Data for coils
    character(len=80) :: coil_set_name = ''
    real(KIND=rkind), allocatable :: inner_radius(:) ! Inner radius of coil conductors
    real(KIND=rkind), allocatable :: outer_radius(:) ! Outer radius of coil conductors
    real(KIND=rkind), allocatable :: z_width(:) ! z width of coils
    real(KIND=rkind), allocatable :: z_center(:) ! z of coil center
    integer(KIND=rkind), allocatable :: n_turns(:) ! number of turns in coil
	integer, allocatable :: n_r_layers(:) ! number of turn layers in radius
	integer, allocatable :: n_z_slices(:) ! number of turn slices in z
	integer, allocatable :: n_filaments(:) ! number of sub-groups of turns = n_r_layers*n_z_slices

! Data for currents
    character(len=80) :: current_set_name = ''
    real(KIND=rkind), allocatable :: I_coil(:) ! current in coil

! Derived type with data for a mirror coil consisting of one or more current loops
	type coil_type
	    real(KIND=rkind) :: inner_radius ! Inner radius of coil conductors
	    real(KIND=rkind) :: outer_radius ! Outer radius of coil conductors
	    real(KIND=rkind) :: z_center ! z of coil center
	    real(KIND=rkind) :: z_width ! coil width
	    integer :: n_turns ! number of turns
	    integer :: n_filaments ! number of sub-groups of the turns
	    integer :: n_r_layers ! number of turn layers in radius
	    integer :: n_z_slices ! number of turn slices in z
	    real(KIND=rkind), allocatable :: r_filament(:) ! r positon of filament, size n_r_layers
	    real(KIND=rkind), allocatable :: z_filament(:) ! z positon of filament, size n_z_slices
	contains
		procedure :: coil_Brz_field_1Amp ! Br(r,z), Bz(r,z for one amp in coil)
	end type coil_type

	type(coil_type), allocatable :: coil(:)

! 2D field data
	integer :: n_r, n_z
    real(KIND=rkind)  :: r_min, r_max, z_min, z_max
    real(KIND=rkind), allocatable  :: r_grid(:), z_grid(:)
    real(KIND=rkind), allocatable  :: Br(:,:), Bz(:,:), Aphi(:,:)

! Location of scrape-off point of last flux surface
    real(KIND=rkind) :: r_LUFS, z_LUFS

 namelist /mirror_magnetics_list/ case_name, coil_data_file, current_data_file, n_coils,&
         & n_r, n_z, r_min, r_max, z_min, z_max, r_LUFS, z_LUFS

 namelist /coil_data_list/ coil_set_name, inner_radius, outer_radius, &
          & z_width, z_center, n_turns, n_r_layers, n_z_slices

 namelist /current_data_list/ current_set_name, I_coil
!********************************************************************

contains

!********************************************************************

  subroutine initialize_mirror_magnetics(read_input)

     use diagnostics_m, only : message, text_message, message_unit, messages_to_stdout, verbosity

    implicit none
    logical, intent(in) :: read_input
    integer :: n_args ! number of command line arguments
    character(len=80) :: namelist_file = ''

 	integer :: input_unit, get_unit_number ! External, free unit finder

 	integer :: i, j
	real(KIND=rkind) :: delta_r, delta_z

 	verbosity = 1
    messages_to_stdout = .true.
    call text_message('Initializing mirror_magnetics_m ', 1)


! Default filename is 'mirror_magnetics.nml'.  Optionally get input file name from command
! line. then copy that file to 'mirror_magnetics.nml'
    n_args = command_argument_count()
    if(n_args > 1) then
        write(*,*) 'Takes zero or one command line argument -> namelist filename'
        stop 'incorrect command line arguments'
    else if (n_args == 1) then
        call get_command_argument(1,namelist_file)
        if (trim(namelist_file) /= 'mirror_magnetics.nml') then ! Don't copy if input
                                                                !filename already mirror_magnetics.nml
 	       call system('cp '//trim(namelist_file)//' mirror_magnetics.nml')
        end if
    end if

! Read mirror magnetics data
    if (read_input .eqv. .true.) then
  	   input_unit = get_unit_number()
       open(unit=input_unit, file='mirror_magnetics.nml',action='read', status='old', form='formatted')
       read(input_unit, mirror_magnetics_list)
        close(unit=input_unit)
    end if
! Write input namelist
    if (verbosity >= 0) then
		write(message_unit, mirror_magnetics_list)
		if (messages_to_stdout) write(*, mirror_magnetics_list)
    end if
	allocate(inner_radius(n_coils))
	allocate(outer_radius(n_coils))
	allocate(z_width(n_coils))
	allocate(z_center(n_coils))
	allocate(n_turns(n_coils))
	allocate(n_r_layers(n_coils))
	allocate(n_z_slices(n_coils))
	allocate(n_filaments(n_coils))
	allocate(I_coil(n_coils))

! Read coil data
    if (read_input .eqv. .true.) then
  	   input_unit = get_unit_number()
       open(unit=input_unit, file=trim(coil_data_file),action='read', status='old', form='formatted')
       read(input_unit, coil_data_list)
        close(unit=input_unit)
    end if
! Write input namelist
    if (verbosity >= 0) then
		write(message_unit, coil_data_list)
		if (messages_to_stdout) write(*, coil_data_list)
    end if

! Read current data
    if (read_input .eqv. .true.) then
  	   input_unit = get_unit_number()
       open(unit=input_unit, file=trim(current_data_file),action='read', status='old', form='formatted')
       read(input_unit, current_data_list)
       close(unit=input_unit)
    end if
! Write input namelist
    if (verbosity >= 0) then
		write(message_unit, current_data_list)
		if (messages_to_stdout) write(*, current_data_list)
    end if

	base_file_name = trim(coil_set_name)//'_'//trim(current_set_name)//'_'//trim(case_name)
	write(*,*) 'base_file_name = ', trim(base_file_name)

	allocate(coil(n_coils))

	do i = 1, n_coils
		coil(i)%inner_radius = inner_radius(i)
		coil(i)%outer_radius = outer_radius(i)
		coil(i)%z_width = z_width(i)
		coil(i)%z_center = z_center(i)
		coil(i)%n_turns = n_turns(i)
		coil(i)%n_r_layers = n_r_layers(i)
		coil(i)%n_z_slices = n_z_slices(i)

		allocate(coil(i)%r_filament(n_r_layers(i)))
		allocate(coil(i)%z_filament(n_z_slices(i)))

! Calculate positions of filaments
		delta_r = (outer_radius(i) - inner_radius(i))/(n_r_layers(i)+1)
		do j = 1,n_r_layers(i)
			coil(i)%r_filament(j) = inner_radius(i)+j*delta_r
		end do
		delta_z = z_width(i)/(n_z_slices(i)+1)
		do j = 1,n_z_slices(i)
			coil(i)%z_filament(j) = z_center(i) - z_width(i)/2._rkind + j*delta_z
		end do
	end do

  end subroutine initialize_mirror_magnetics
!********************************************************************

  subroutine coil_Brz_field_1Amp(this, r, z, Br, Bz, Aphi, equib_err)
! Calculates the field at location (r,z) for this coil instance.  Sums over multiple
! filaments of the coil.

    use diagnostics_m, only : message_unit, message, text_message
    use B_loop_m, only : Brz_loop_scaled

    implicit none

    class(coil_type) :: this

    real(KIND=rkind), intent(in) :: r, z
    real(KIND=rkind), intent(out) :: Br, Bz, Aphi

    !   Error returns
    character(len=60) :: equib_err

    real(KIND=rkind) ::  a ! Filament loop radius
    real(KIND=rkind) ::  z_relative ! z coordiate of field point relative to filament z
    real(KIND=rkind) ::  Br_f, Bz_f, Aphi_f
    integer :: i,j

    equib_err = ''

! Check that we are in the box
!     if (r < box_rmin .or. r > box_rmax) equib_err = 'R out_of_bounds'
!     if (z < box_zmin .or. z > box_zmax) equib_err = 'z out_of_bounds'
!
!     if (equib_err /= '') then
!         call text_message('mirror_magnetics:  equib_err ', equib_err, 1)
!         return
!     end if

	Br = zero
	Bz = zero
	Aphi = zero

	do i = 1,this%n_r_layers
		a = this%r_filament(i)
		do j = 1, this%n_z_slices
			z_relative = z - this%z_filament(j)
!			write(*,*) 'a = ', a, 'this%z_filament(j) = ', this%z_filament(j), '   z_relative = ',z_relative
			call Brz_loop_scaled(r/a, z_relative/a, Br_f, Bz_f, Aphi_f)
			Br = Br + Br_f/a
			Bz = Bz + Bz_f/a
			Aphi = Aphi + Aphi_f*a
		end do
	end do

! normalize to number of filaments for total one amp coil current
	Br = Br/(this%n_r_layers*this%n_z_slices)
	Bz = Bz/(this%n_r_layers*this%n_z_slices)
	Aphi = Aphi/(this%n_r_layers*this%n_z_slices)

	return
 end subroutine coil_Brz_field_1Amp

!********************************************************************

  subroutine mirror_Brz_field(r, z, Br, Bz, Aphi, equib_err)
! Sums over multiple coils to get total field at field point (r,z)

    use diagnostics_m, only : message_unit, message, text_message
    use B_loop_m, only : Brz_loop_scaled

    implicit none

    real(KIND=rkind), intent(in) :: r, z
    real(KIND=rkind), intent(out) :: Br, Bz, Aphi

    !   Error returns
    character(len=60) :: equib_err

    real(KIND=rkind) ::  Br_i, Bz_i, Aphi_i ! Fields for i'th coil
    integer :: i,j

    equib_err = ''

! Check that we are in the box
!     if (r < box_rmin .or. r > box_rmax) equib_err = 'R out_of_bounds'
!     if (z < box_zmin .or. z > box_zmax) equib_err = 'z out_of_bounds'
!
!     if (equib_err /= '') then
!         call text_message('mirror_magnetics:  equib_err ', equib_err, 1)
!         return
!     end if

	Br = zero
	Bz = zero
	Aphi = zero

	do i = 1, n_coils
		call coil(i)%coil_Brz_field_1Amp(r, z, Br_i, Bz_i, Aphi_i, equib_err)
		Br = Br + Br_i*coil(i)%n_turns*I_coil(i)
		Bz = Bz + Bz_i*coil(i)%n_turns*I_coil(i)
		Aphi = Aphi + Aphi_i*coil(i)%n_turns*I_coil(i)
	end do

	return
 end subroutine mirror_Brz_field

!********************************************************************

  subroutine calculate_B_on_rz_grid
! Calculates total (Br,Bz) field from all coils on an evenly spaced grid in (r,z)

	implicit none

    real(KIND=rkind) :: rvec(3)

    !   Error returns
    character(len=60) :: equib_err

    integer :: i,j

	allocate(r_grid(n_r))
	allocate(z_grid(n_z))
	allocate(Br(n_r, n_z))
	allocate(Bz(n_r, n_z))
	allocate(Aphi(n_r, n_z))

	if (n_r == 1) then
		r_grid(1) = r_min
	else
		do i = 1, n_r
			r_grid(i) = r_min +(r_max - r_min)*(i-1)/(n_r - 1)
		end do
	end if
	write(*,*) ' r_min = ', r_min, ' r_max = ', r_max
! 	write(*,*) ' r_grid = ', r_grid

	if (n_z == 1) then
		z_grid(1) = z_min
	else
		do j = 1, n_z
		z_grid(j) = z_min +(z_max - z_min)*(j-1)/(n_z - 1)
		end do
	end if
	write(*,*) '   z_min = ', z_min, '   z_max = ', z_max
!	write(*,*) ' z_grid = ', z_grid

	do i = 1, n_r
		do j = 1, n_z
			call mirror_Brz_field(r_grid(i), z_grid(j), Br(i,j), Bz(i,j), Aphi(i,j), equib_err)
		end do
	end do
!	write(*,*) 'Aphi(2,:) = ', Aphi(2,:)

  end subroutine calculate_B_on_rz_grid
!****************************************************************************

  subroutine write_mirror_fields_Brz_NC
! writes a netCDF file with grid and field data as generated by calculate_B_on_rz_grid
! subroutine above.

    use diagnostics_m, only : run_label
    use netcdf

    implicit none

    integer :: j

! netCDF Declarations
    integer :: ncid

! Declarations: dimensions
    integer, parameter :: n_dims = 2
    integer :: n_r_id, n_z_id

! Declarations: variable IDs
    integer, parameter :: n_vars =  9
   integer :: r_min_id, r_max_id, z_min_id, z_max_id, r_LUFS_id, z_LUFS_id, &
            & r_grid_id, z_grid_id, Br_id, Bz_id, Aphi_id

 !  File name for  output
    character(len=100) :: out_filename

    NC_file_name = 'Brz_fields.'//trim(base_file_name)//'.nc'

!   Open NC file
    call check( nf90_create(trim(NC_file_name), nf90_clobber, ncid) )
    write(*,*) 'writing file ', trim(NC_file_name)


!   Define NC dimensions
    call check( nf90_def_dim(ncid, 'n_r', n_r, n_r_id))
    call check( nf90_def_dim(ncid, 'n_z', n_z, n_z_id))

! Define NC variables
    call check( nf90_def_var(ncid, 'r_min', NF90_DOUBLE, r_min_id))
    call check( nf90_def_var(ncid, 'r_max', NF90_DOUBLE, r_max_id))
    call check( nf90_def_var(ncid, 'z_min', NF90_DOUBLE, z_min_id))
    call check( nf90_def_var(ncid, 'z_max', NF90_DOUBLE, z_max_id))
    call check( nf90_def_var(ncid, 'r_LUFS', NF90_DOUBLE, r_LUFS_id))
    call check( nf90_def_var(ncid, 'z_LUFS', NF90_DOUBLE, z_LUFS_id))
    call check( nf90_def_var(ncid, 'r_grid', NF90_DOUBLE, [n_r_id], r_grid_id))
    call check( nf90_def_var(ncid, 'z_grid', NF90_DOUBLE, [n_z_id], z_grid_id))
    call check( nf90_def_var(ncid, 'Br', NF90_DOUBLE, [n_r_id,n_z_id], Br_id))
    call check( nf90_def_var(ncid, 'Bz', NF90_DOUBLE, [n_r_id,n_z_id], Bz_id))
    call check( nf90_def_var(ncid, 'Aphi', NF90_DOUBLE, [n_r_id,n_z_id], Aphi_id))

! Put global attributes
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'NC_file_name', NC_file_name))

	call check( nf90_enddef(ncid))

! Put NC variables
    call check( nf90_put_var(ncid, r_min_id, r_min))
    call check( nf90_put_var(ncid, r_max_id, r_max))
    call check( nf90_put_var(ncid, z_min_id, z_min))
    call check( nf90_put_var(ncid, z_max_id, z_max))
    call check( nf90_put_var(ncid, r_LUFS_id, r_LUFS))
    call check( nf90_put_var(ncid, z_LUFS_id, z_LUFS))
    call check( nf90_put_var(ncid, r_grid_id, r_grid))
    call check( nf90_put_var(ncid, z_grid_id, z_grid))
    call check( nf90_put_var(ncid, Br_id, Br))
    call check( nf90_put_var(ncid, Bz_id, Bz))
    call check( nf90_put_var(ncid, Aphi_id, Aphi))

!   Close the NC file
    call check( nf90_close(ncid) )

! 		do j = 1, n_z
! 			write(*,*) '  z = ', z_grid(j), '  Br = ', Br(2,j), '  Bz = ', Bz(2,j),&
! 			         & '  Aphi = ', Aphi(2,j)
! 		end do

  end subroutine write_mirror_fields_Brz_NC

!****************************************************************************

  subroutine read_mirror_fields_Brz_NC(in_filename)
! Reads a netCDF file as written by write routine above

    use diagnostics_m, only : run_label
    use netcdf

    implicit none

    integer j
 !  File name for netCDF input
    character(len=*) :: in_filename

! netCDF Declarations
    integer :: ncid

! Declarations: dimensions
    integer, parameter :: n_dims = 2
    integer :: n_r_id, n_z_id

! Declarations: variable IDs
    integer, parameter :: n_vars =  9
   integer :: r_min_id, r_max_id, z_min_id, z_max_id, r_LUFS_id, z_LUFS_id, &
            & r_grid_id, z_grid_id, Br_id, Bz_id, Aphi_id

!   Open NC file
    call check( nf90_open(in_filename, nf90_nowrite, ncid) )
    write(*,*) 'reading file ', trim(in_filename)


!   Inquire NC dimensions
    call check( nf90_inq_dimid(ncid, 'n_r', n_r_id))
    call check(nf90_inquire_dimension(ncid, n_r_id, len=n_r))
    call check( nf90_inq_dimid(ncid, 'n_z', n_z_id))
    call check(nf90_inquire_dimension(ncid, n_z_id, len=n_z))

! If arrays are already allocated deallocate arrays
    if (allocated(r_grid)) deallocate(r_grid)
    if (allocated(z_grid)) deallocate(z_grid)
    if (allocated(Br)) deallocate(Br)
    if (allocated(Bz)) deallocate(Bz)
    if (allocated(Aphi)) deallocate(Aphi)

!   Allocate arrays
    allocate(r_grid(n_r))
    allocate(z_grid(n_z))
    allocate(Br(n_r,n_z))
    allocate(Bz(n_r,n_z))
    allocate(Aphi(n_r,n_z))

!   Inquire and get NC variables
    call check(nf90_inq_varid(ncid, 'r_min', r_min_id))
    call check( nf90_get_var(ncid, r_min_id, r_min))
    call check(nf90_inq_varid(ncid, 'r_max', r_max_id))
    call check( nf90_get_var(ncid, r_max_id, r_max))
    call check(nf90_inq_varid(ncid, 'z_min', z_min_id))
    call check( nf90_get_var(ncid, z_min_id, z_min))
    call check(nf90_inq_varid(ncid, 'z_max', z_max_id))
    call check( nf90_get_var(ncid, z_max_id, z_max))
    call check(nf90_inq_varid(ncid, 'r_LUFS', r_LUFS_id))
    call check( nf90_get_var(ncid, r_LUFS_id, r_LUFS))
    call check(nf90_inq_varid(ncid, 'z_LUFS', z_LUFS_id))
    call check( nf90_get_var(ncid, z_LUFS_id, z_LUFS))
    call check(nf90_inq_varid(ncid, 'r_grid', r_grid_id))
    call check( nf90_get_var(ncid, r_grid_id, r_grid))
    call check(nf90_inq_varid(ncid, 'z_grid', z_grid_id))
    call check( nf90_get_var(ncid, z_grid_id, z_grid))
    call check(nf90_inq_varid(ncid, 'Br', Br_id))
    call check( nf90_get_var(ncid, Br_id, Br))
    call check(nf90_inq_varid(ncid, 'Bz', Bz_id))
    call check( nf90_get_var(ncid, Bz_id, Bz))
    call check(nf90_inq_varid(ncid, 'Aphi', Aphi_id))
    call check( nf90_get_var(ncid, Aphi_id, Aphi))

! Get global attributes
    call check( nf90_get_att(ncid, NF90_GLOBAL, 'NC_file_name', NC_file_name))

!   Close the NC file
    call check( nf90_close(ncid) )

! 		do j = 1, n_z
! 			write(*,*) '  z = ', z_grid(j), '  Br = ', Br(2,j), '  Bz = ', Bz(2,j),&
! 			         & '  Aphi = ', Aphi(2,j)
! 		end do

  end subroutine read_mirror_fields_Brz_NC

!****************************************************************************

  subroutine check(status)
    use netcdf
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check

!****************************************************************************

!********************************************************************

    subroutine deallocate_mirror_magnetics_m
        return ! nothing to deallocate
    end subroutine deallocate_mirror_magnetics_m

end module mirror_magnetics_m
