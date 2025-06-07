program test_B_loop

  use B_loop_m, only : Brz_loop, Brz_loop_scaled, Bxyz_loop, Bxyz_loop_scaled, &
                    &  Aphi_xyz_loop, Aphi_xyz_loop_scaled

  implicit none
    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
	real(KIND=rkind),parameter :: zero = 0.0_rkind
	real(KIND=rkind),parameter :: one = 1.0_rkind

  integer :: i, j

  real(KIND=rkind)  :: a, x, y, z, r
  real(KIND=rkind)  :: Bx, By, Bz, Br, Aphi
  real(KIND=rkind)  :: Bx_scaled, By_scaled, Bz_scaled, Br_scaled, Aphi_scaled
  real(KIND=rkind)  :: rN(3) = (/zero, .009_rkind, .5_rkind/)
  real(KIND=rkind)  :: zN(2) = (/zero, 0.5_rkind/)

	y = 0._rkind

 write(*,*) ' '
 write(*,*) 'Br_loop and scaled using Br_loop_scaled'

	a = 1.0_rkind
	write(*,*) 'a =  ', a
	do i = 1, 3
    	r = rN(i)
        do j = 1,2
            write(*,*) ' '
            z = zN(j)
            call Brz_loop(r, z, a, Br, Bz, Aphi)
            write(*,*) 'r = ', r, '  z = ', z,'  Br = ', Br, '  Bz = ', Bz, &
                       & ' Aphi = ', Aphi
            call Brz_loop_scaled(r/a, z/a, Br_scaled, Bz_scaled, Aphi_scaled)
            Br_scaled = Br_scaled/a
            Bz_scaled = Bz_scaled/a
            Aphi_scaled = Aphi_scaled*a
            write(*,*) 'Scaled'
            write(*,*) 'r = ', r, '  z = ', z,'  Br = ', Br_scaled,&
                    & '  Bz = ', Bz_scaled, ' Aphi = ', Aphi_scaled
        end do
    end do

 write(*,*) ' '
 write(*,*) 'Br_loop and scaled using Br_loop'

	a = 1.0_rkind
	write(*,*) 'a =  ', a
	do i = 1, 3
    	r = rN(i)
        do j = 1,2
            write(*,*) ' '
            z = zN(j)
            call Brz_loop(r, z, a, Br, Bz, Aphi)
            write(*,*) 'r = ', r, '  z = ', z,'  Br = ', Br, '  Bz = ', Bz, &
                       & ' Aphi = ', Aphi
            call Brz_loop(r/a, z/a, 1.0_rkind, Br_scaled, Bz_scaled, Aphi_scaled)
            Br_scaled = Br_scaled/a
            Bz_scaled = Bz_scaled/a
            Aphi_scaled = Aphi_scaled*a
            write(*,*) 'Scaled'
            write(*,*) 'r = ', r, '  z = ', z,'  Br = ', Br_scaled,&
                    & '  Bz = ', Bz_scaled, ' Aphi = ', Aphi_scaled
        end do
    end do


 write(*,*) ' '
 write(*,*) ' '
 write(*,*) '******************************************************* '
 write(*,*) 'Br_loop and scaled using Br_loop_scaled'

	a = 1.5_rkind
	write(*,*) 'a =  ', a
	do i = 1, 3
    	r = rN(i)
        do j = 1,2
            write(*,*) ' '
            z = zN(j)
            call Brz_loop(r, z, a, Br, Bz, Aphi)
            write(*,*) 'r = ', r, '  z = ', z,'  Br = ', Br, '  Bz = ', Bz, &
                       & ' Aphi = ', Aphi
            call Brz_loop_scaled(r/a, z/a, Br_scaled, Bz_scaled, Aphi)
            Br_scaled = Br_scaled/a
            Bz_scaled = Bz_scaled/a
            Aphi_scaled = Aphi_scaled*a
            write(*,*) 'Scaled'
            write(*,*) 'r = ', r, '  z = ', z,'  Br = ', Br_scaled,&
                    & '  Bz = ', Bz_scaled, ' Aphi = ', Aphi_scaled
        end do
    end do

 write(*,*) ' '
 write(*,*) 'Br_loop and scaled using Br_loop'

	a = 1.5_rkind
	write(*,*) 'a =  ', a
	do i = 1, 3
    	r = rN(i)
        do j = 1,2
            write(*,*) ' '
            z = zN(j)
            call Brz_loop(r, z, a, Br, Bz, Aphi)
            write(*,*) 'r = ', r, '  z = ', z,'  Br = ', Br, '  Bz = ', Bz, &
                       & ' Aphi = ', Aphi
            call Brz_loop(r/a, z/a, 1.0_rkind, Br_scaled, Bz_scaled, Aphi_scaled)
            Br_scaled = Br_scaled/a
            Bz_scaled = Bz_scaled/a
            Aphi_scaled = Aphi_scaled*a
            write(*,*) 'Scaled'
            write(*,*) 'r = ', r, '  z = ', z,'  Br = ', Br_scaled,&
                    & '  Bz = ', Bz_scaled, ' Aphi = ', Aphi_scaled
        end do
    end do


 write(*,*) ' '
 write(*,*) ' '
 write(*,*) '******************************************************* '
 write(*,*) 'Bxyz'

	a = 1.5_rkind
	write(*,*) 'a =  ', a
	do i = 1, 3
    	x = rN(i)
        do j = 1,2
            write(*,*) ' '
            z = zN(j)
            call Bxyz_loop(x, y, z, a, Bx, By, Bz)
            write(*,*) 'x = ', x, '  y = ', y, '  z = ', z,'  Bx = ', Bx, '  By = ', By, &
               & '  Bz = ', Bz
            call Bxyz_loop_scaled(x/a, y/a, z/a, Bx_scaled, By_scaled, Bz_scaled)
            Bx_scaled = Bx_scaled/a
            By_scaled = By_scaled/a
            Bz_scaled = Bz_scaled/a
            write(*,*) 'Scaled'
            write(*,*) 'x = ', x, '  y = ', y, '  z = ', z,'  Bx = ', Bx_scaled,&
              & '  By = ', By_scaled, '  Bz = ', Bz_scaled
        end do
    end do

 write(*,*) ' '
 write(*,*) ' '
 write(*,*) '******************************************************* '
 write(*,*) 'Aphi'

	a = 1.5_rkind
	write(*,*) 'a =  ', a
	do i = 1, 3
    	x = rN(i)
        do j = 1,2
            write(*,*) ' '
            z = zN(j)
            call Aphi_xyz_loop(x, y, z, a, Aphi)
            call Aphi_xyz_loop_scaled(x/a, y/a, z/a, Aphi_scaled)
            Aphi_scaled = a*Aphi_scaled
            write(*,*) 'x = ', x, '  y = ', y, '  z = ', z,'  Aphi = ', Aphi,&
              & '  Aphi_scaled = ', Aphi_scaled
        end do
    end do


end program test_B_loop


