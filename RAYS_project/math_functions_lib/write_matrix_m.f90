 MODULE write_matrix_m

! Generic interface to functions that write out m x n matrix to file or std out
! write_matrix_S : single precision, real
! write_matrix_rkind : precision rkind, real
! write_Cmatrix_S : single precision, complex
! write_Cmatrix_rkind : precision rkind, complex

! N.B. Default unit is ISO_FORTRAN_ENV OUTPUT_UNIT. Can be changed in call by using
! optional argument 'out_unit'

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

	use iso_fortran_env, only: output_unit

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind),parameter :: two = 2.0_rkind

    interface write_matrix
        module procedure write_matrix_S, write_matrix_rkind, write_Cmatrix_S,&
              & write_Cmatrix_rkind
    end interface

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

    subroutine write_matrix_S(mess, value, m_dim, n_dim, unit)

    implicit none

    character (len=*), intent (in) :: mess
    integer, intent (in) :: m_dim, n_dim
    real, intent (in) :: value(m_dim, n_dim)
    integer, optional, intent (in) :: unit
    integer :: out_unit

    real(kind=rkind) :: v_min, v_max
    integer :: i, j

    if(present(unit)) then
    	out_unit = unit
    else
    	out_unit = output_unit
    end if

	v_min = huge(v_min)
	v_max = tiny(v_max)

	do i = 1, m_dim
		do j = 1, n_dim

			if ((abs(value(i,j)) > 0.).and.(abs(value(i,j)) < v_min)) &
				& v_min = abs(value(i,j))
			if (abs(value(i,j)) > v_max ) v_max = abs(value(i,j))
		end do
	end do

	if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then ! Write in floating point

		write(out_unit, '(a," ")')
		write(out_unit, '(a, " = ")' ) trim(mess)
!		write(*, '(a, " = ")' ) trim(mess)

		do i = 1, m_dim
			write(out_unit, '(i3,2x,8(3x,f14.8))' )  i, (value(i, j), j=1, n_dim)
!			write(*, '(i3,2x,8(3x,f14.8))' )  i, (value(i, j), j=1, n_dim)
		end do

	else ! Write in scientific notation

		write(out_unit, '(a," ")')
		write(out_unit, '(a, " = ")' ) trim(mess)
!		write(*, '(a, " = ")' ) trim(mess)
		do i = 1, m_dim
			write(out_unit, '(i3,2x,8(3x,(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
!			write(*, '(i3,2x,8(3x,(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
		end do

	end if

    return
    end subroutine write_matrix_S

    subroutine write_matrix_rkind(mess, value, m_dim, n_dim, unit)

    implicit none

    character (len=*), intent (in) :: mess
    integer, intent (in) :: m_dim, n_dim
    real(kind=rkind), intent (in) :: value(m_dim, n_dim)
    integer, optional, intent (in) :: unit
    integer :: out_unit

    real(kind=rkind) :: v_min, v_max
    integer :: i, j

    if(present(unit)) then
    	out_unit = unit
    else
    	out_unit = output_unit
    end if

	v_min = huge(v_min)
	v_max = tiny(v_max)

	do i = 1, m_dim
		do j = 1, n_dim

			if ((abs(value(i,j)) > 0.).and.(abs(value(i,j)) < v_min)) &
				& v_min = abs(value(i,j))
			if (abs(value(i,j)) > v_max ) v_max = abs(value(i,j))
		end do
	end do

	if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then ! Write in floating point

		write(out_unit, '(a," ")')
		write(out_unit, '(a, " = ")' ) trim(mess)
!		write(*, '(a, " = ")' ) trim(mess)

		do i = 1, m_dim
			write(out_unit, '(i3,2x,8(3x,f14.8))' )  i, (value(i, j), j=1, n_dim)
!			write(*, '(i3,2x,8(3x,f14.8))' )  i, (value(i, j), j=1, n_dim)
		end do

	else ! Write in scientific notation

		write(out_unit, '(a," ")')
		write(out_unit, '(a, " = ")' ) trim(mess)
!		write(*, '(a, " = ")' ) trim(mess)
		do i = 1, m_dim
			write(out_unit, '(i3,2x,8(3x,(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
!			write(*, '(i3,2x,8(3x,(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
		end do

	end if

    return
    end subroutine write_matrix_rkind

    subroutine write_Cmatrix_S(mess, value, m_dim, n_dim, unit)

    implicit none

    character (len=*), intent (in) :: mess
    integer, intent (in) :: m_dim, n_dim
    complex, intent (in) :: value(m_dim, n_dim)
    integer, optional, intent (in) :: unit
    integer :: out_unit

    real(kind=rkind) :: v_min, v_max
    integer :: i, j

    if(present(unit)) then
    	out_unit = unit
    else
    	out_unit = output_unit
    end if

	v_min = huge(v_min)
	v_max = tiny(v_max)

	do i = 1, m_dim
		do j = 1, n_dim

			if ((abs(value(i,j)) > 0.).and.(abs(value(i,j)) < v_min)) &
				& v_min = abs(value(i,j))
			if (abs(value(i,j)) > v_max ) v_max = abs(value(i,j))
		end do
	end do

	if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then ! Write in floating point

		write(out_unit, '(a," ")')
		write(out_unit, '(a, " = ")' ) trim(mess)
!		write(*, '(a, " = ")' ) trim(mess)

		do i = 1, m_dim
			write(out_unit, '(i3,2x,8(3x,2f14.8))' )  i, (value(i, j), j=1, n_dim)
!			write(*, '(i3,2x,8(3x,2f14.8))' )  i, (value(i, j), j=1, n_dim)
		end do

	else ! Write in scientific notation

		write(out_unit, '(a," ")')
		write(out_unit, '(a, " = ")' ) trim(mess)
!		write(*, '(a, " = ")' ) trim(mess)
		do i = 1, m_dim
			write(out_unit, '(i3,2x,8(3x,2(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
!			write(*, '(i3,2x,8(3x,2(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
		end do

	end if

    return
    end subroutine write_Cmatrix_S

    subroutine write_Cmatrix_rkind(mess, value, m_dim, n_dim, unit)

    implicit none

    character (len=*), intent (in) :: mess
    integer, intent (in) :: m_dim, n_dim
    complex(kind=rkind), intent (in) :: value(m_dim, n_dim)
    integer, optional, intent (in) :: unit
    integer :: out_unit

    real(kind=rkind) :: v_min, v_max
    integer :: i, j

    if(present(unit)) then
    	out_unit = unit
    else
    	out_unit = output_unit
    end if

	v_min = huge(v_min)
	v_max = tiny(v_max)

	do i = 1, m_dim
		do j = 1, n_dim

			if ((abs(value(i,j)) > 0.).and.(abs(value(i,j)) < v_min)) &
				& v_min = abs(value(i,j))
			if (abs(value(i,j)) > v_max ) v_max = abs(value(i,j))
		end do
	end do

	if ((v_max < 1.e4) .and. (v_min > 1.e-3)) then ! Write in floating point

		write(out_unit, '(a," ")')
		write(out_unit, '(a, " = ")' ) trim(mess)
!		write(*, '(a, " = ")' ) trim(mess)

		do i = 1, m_dim
			write(out_unit, '(i3,2x,8(3x,2f14.8))' )  i, (value(i, j), j=1, n_dim)
		end do

	else ! Write in scientific notation

		write(out_unit, '(a," ")')
		write(out_unit, '(a, " = ")' ) trim(mess)
!		write(*, '(a, " = ")' ) trim(mess)
		do i = 1, m_dim
			write(out_unit, '(i3,2x,8(3x,2(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
!			write(*, '(i3,2x,8(3x,2(1pe14.6)))' )  i, (value(i, j), j=1, n_dim)
		end do

	end if

    return
    end subroutine write_Cmatrix_rkind

 END MODULE write_matrix_m