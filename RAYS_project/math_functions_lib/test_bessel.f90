  implicit none
  INTEGER, PARAMETER  :: rkind = selected_real_kind(15,307)
  real(KIND=rkind),parameter :: zero = 0.0_rkind, one = 1.0_rkind, two = 2.0_rkind
  integer, parameter :: nmin=-2, nmax=2
  complex(KIND=rkind)::z, b(nmin:nmax), bp(nmin:nmax)
  integer::i,j

!   do i=-2, 2; do j=-2, 2
  do i=0,1; do j=0, 0
  z = cmplx(i*one, j*one)
  call ebessel_dbb(z, nmin, nmax, b, bp)
  write (6, "(/,' z=',2(f12.4))" ) z
  write (6, "( 'b=', 10(2(1pe15.7),5x) )" ) b
  write(6, "( 'bp=', 10(2(1pe15.7),5x) )" ) bp
  end do; end do

  end

