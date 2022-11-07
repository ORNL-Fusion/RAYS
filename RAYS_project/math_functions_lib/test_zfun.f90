program test_zfun

! Note: fzeta agrees with ZFUN0 for z real or Im(z) > 0 but differs for 
! Im(z) < 0

	implicit none
	complex z
	complex zf, ZFUN0, fzeta
	integer i,j
	character*15 fname
	
	go to 100
	
	do i= -4, 4, 2
	do j = -4, 4, 2
	
	
	z = cmplx(i,j)
	
	fname = 'zfun0'
	zf = ZFUN0(z, 1.)
	write (6,5) z, fname, zf
	
	fname = 'fzeta'
	zf = fzeta(z)
	write (6,5) z, fname, zf
	
5	format ('z = ',2(f12.5),'  ', a,' = ',2(1pe18.8))
	write (6,30)
30	format(/)

    end do
    end do

100	continue

	do j = -1, 1
	do i = -10, 10
	
	
	z = cmplx(i,j)
	
	fname = 'zfun0 kx > 0'
	write (6,*) ' '
	zf = ZFUN0(z, 1.)
	write (6,5) z, fname, zf
	
	fname = 'zfun0 kx < 0'
	zf = ZFUN0(z, -1.)
	write (6,5) z, fname, zf
	 	
 	fname = 'fzeta'
 	zf = fzeta(z)
 	write (6,5) z, fname, zf
		
	!	write (6,30)

    end do
    end do

end program test_zfun
