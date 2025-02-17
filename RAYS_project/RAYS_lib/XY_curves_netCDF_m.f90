 module XY_curves_netCDF_m

! Utilities to conveniently generate a netCDF file containing lists of X/Y vectors
! suitable for plotting, therefore called curves.  It provides a derived type,
! XY_curve_netCDF, containing the data needed for a single curve. It provides a subroutine,
! write_XY_curves_netCDF, which accepts a list of type XY_curve_netCDF instances and writes
! a netCDF file contining the curve data plus some metadata.

	implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    integer, parameter :: str_max_len = 60

! Derived type containing data needed to define an XY curve
	type XY_curve_netCDF
		character(len=:), allocatable :: grid_name
		character(len=:), allocatable :: curve_name
		real(kind=rkind), allocatable :: grid(:)
		real(kind=rkind), allocatable :: curve(:)
	end type XY_curve_netCDF

!********************************************************************

contains

!********************************************************************

  subroutine write_XY_curves_netCDF(curve_list, out_filename)

	use netcdf

    implicit none

! List of curves
	type(XY_curve_netCDF), intent(in) :: curve_list(:)
 !  File name for  output
    character(len=*), intent(in) :: out_filename

	integer :: i, itemp
    integer :: ncid

! Declarations: dimensions
!     integer :: n_dims
	integer, allocatable :: n_grid(:)
    integer :: n_curves_id, grid_max_len_id, curve_max_len_id
    integer :: n_curves, grid_max_len, curve_max_len
    integer :: grid_name_max_len_id, curve_name_max_len_id, str_max_len_id
    integer :: grid_name_max_len, curve_name_max_len

! Declarations: variables
! Declarations: variable IDs
    integer :: n_vars ! =  n_curves + 2 (grid_name and curve_name)
    integer :: curve_name_id, n_grid_id, grid_name_id, grid_id, curve_id

!   Declare netCDF variables
	integer, allocatable :: n_grid_nc(:)
    real(kind=rkind), allocatable :: grid_nc(:,:)
    real(kind=rkind), allocatable :: curve_nc(:,:)
    character(len=:), allocatable :: curve_name_nc(:)
    character(len=:), allocatable :: grid_name_nc(:)

    character(len=NF90_MAX_NAME)   :: NC_file

	write(*,*) " "
	write(*,*) "write_XY_curves_netCDF, file = ", out_filename//'.nc'
	n_curves = size(curve_list)
! 	do i = 1, n_curves
! 	 write(*,*) " "
! 	 write(*,*) "i = ", i
! 	 write(*,*) "curve_list(i)%grid_name = ", curve_list(i)%grid_name
! 	 write(*,*) "curve_list(i)%grid = ", curve_list(i)%grid
! 	 write(*,*) "curve_list(i)%curve_name = ", curve_list(i)%curve_name
! 	 write(*,*) "curve_list(i)%curve = ", curve_list(i)%curve
! 	end do

    allocate(n_grid(n_curves),source=0)

! Find Maximum sizes of stuff
	grid_max_len = 0
	grid_name_max_len = 0
	curve_name_max_len = 0
	do i = 1, n_curves
		n_grid(i) = size(curve_list(i)%grid)
		if (n_grid(i) > grid_max_len) grid_max_len = n_grid(i)
		itemp = len(trim(curve_list(i)%grid_name))
		if (itemp > grid_name_max_len) grid_name_max_len = itemp
		itemp = len(trim(curve_list(i)%curve_name))
		if (itemp > curve_name_max_len) curve_name_max_len = itemp
	end do

! If nc arrays already allocated deallocate them
    if (allocated(grid_nc)) deallocate(n_grid_nc)
    if (allocated(grid_nc)) deallocate(grid_nc)
    if (allocated(curve_nc)) deallocate(curve_nc)
    if (allocated(grid_name_nc)) deallocate(grid_name_nc)
    if (allocated(curve_name_nc)) deallocate(curve_name_nc)

!   Allocate nc arrays
    allocate(grid_nc(grid_max_len,n_curves),source=0.0_rkind)
    allocate(curve_nc(grid_max_len,n_curves),source=0.0_rkind)
    allocate(character(len=grid_name_max_len) :: grid_name_nc(n_curves))
    allocate(character(len=curve_name_max_len) :: curve_name_nc(n_curves))

! Open NC file
	NC_file = trim(out_filename)//'.nc'
	call check( nf90_create(trim(NC_file), nf90_clobber, ncid) )

!   Define NC dimensions
    call check( nf90_def_dim(ncid, 'n_curves', n_curves, n_curves_id))
    call check( nf90_def_dim(ncid, 'grid_max_len', grid_max_len, grid_max_len_id))
!     call check( nf90_def_dim(ncid, 'str_max_len', str_max_len, str_max_len_id))
    call check( nf90_def_dim(ncid, 'grid_name_max_len_id', grid_name_max_len, grid_name_max_len_id))
    call check( nf90_def_dim(ncid, 'curve_name_max_len_id', curve_name_max_len, curve_name_max_len_id))

! Define NC variables
    call check( nf90_def_var(ncid, 'n_grid', NF90_INT, [n_curves_id], n_grid_id))
    call check( nf90_def_var(ncid, 'grid', NF90_DOUBLE, [grid_max_len_id,n_curves_id], grid_id))
    call check( nf90_def_var(ncid, 'curve', NF90_DOUBLE, [grid_max_len_id,n_curves_id], curve_id))
    call check( nf90_def_var(ncid, 'grid_name', NF90_CHAR, [grid_name_max_len_id,n_curves_id],&
                                   & grid_name_id))
    call check( nf90_def_var(ncid, 'curve_name', NF90_CHAR, [curve_name_max_len_id,n_curves_id],&
                                   & curve_name_id))


! Put global attributes
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'out_filename', out_filename))
    call check( nf90_put_att(ncid, NF90_GLOBAL, 'date_vector', '12345'))
!  write(*,*) "Got to 5"

	call check( nf90_enddef(ncid))

!Load up NC variables
	do i = 1, n_curves
		grid_nc(:,i) = curve_list(i)%grid(:)
		curve_nc(:,i) = curve_list(i)%curve(:)
		grid_name_nc(i) = trim(curve_list(i)%grid_name)
		curve_name_nc(i) = trim(curve_list(i)%curve_name)
	end do

! 	write(*,*) " "
! 	write(*,*) "NC variables loaded"
! 	do i = 1, n_curves
! 	 write(*,*) " "
! 	 write(*,*) "i = ", i
! 	 write(*,*) "grid_name_nc(i) = ", grid_name_nc(i)
! 	 write(*,*) "grid_nc(:,i) = ", grid_nc(:,i)
! 	 write(*,*) "curve_name_nc(i) = ", curve_name_nc(i)
! 	 write(*,*) "curve_nc(:,i) = ", curve_nc(:,i)
! 	end do

! Put NC variables
    call check( nf90_put_var(ncid, n_grid_id, n_grid))
    call check( nf90_put_var(ncid, grid_id, grid_nc))
    call check( nf90_put_var(ncid, curve_id, curve_nc))
    call check( nf90_put_var(ncid, grid_name_id, grid_name_nc))
    call check( nf90_put_var(ncid, curve_name_id, curve_name_nc))
!  write(*,*) "Got to 4"
!  write(*,*) "grid name = ", grid_name_nc

! Close the NC file
    call check( nf90_close(ncid) )

 end subroutine write_XY_curves_netCDF

!*************************************************************************

  subroutine check(status)
    use netcdf
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check

!*************************************************************************

 end module XY_curves_netCDF_m
