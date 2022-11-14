program test_zfun

	use zfunctions_m, only : zfun, zfun0

	implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals

	complex z
	complex zf

	complex(kind = rkind) :: zD
	complex(kind = rkind) ::  zfD

	integer i,j
	character*15 fname
	

! Single precision

5	format (' z = ',2(f12.5),'  ', a,'  = ',2(1pe18.10))
6	format (' ')

	do j = -1, 1
	do i = -10, 10
		
	z = cmplx(i,j)
	zD = cmplx(i,j)
	
	fname = 'zfun0 kx > 0   '
	zf = ZFUN0(z, 1.)
	write (6,5) z, fname, zf

	fname = 'zfun0_D, kx > 0   '
	zfD = ZFUN0(zD, 1._rkind)
	write (6,5) zD, fname, zfD
	
	fname = 'zfun0 kx < 0   '
	zf = ZFUN0(z, -1.)
	write (6,5) z, fname, zf

 	fname = 'zfun_D      '
 	zfD = zfun(zD)
 	write (6,5) zD, fname, zfD
	 	
 	fname = 'zfun        '
 	zf = zfun(z)
 	write (6,5) z, fname, zf
	 	
 	fname = 'zfun_D'
 	zfD = zfun(zD)
 	write (6,5) zD, fname, zfD
		
   	write (6,6)

    end do
	write (6,30)
30	format(/)
    end do

end program test_zfun
