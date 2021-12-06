!   implicit none
!   integer, parameter :: nmin=-3, nmax=5
!   complex::z, b(nmin:nmax), bp(nmin:nmax)
!   integer::i,j
!   
!   do i=-2, 2; do j=-2, 2
!   z = cmplx(i,j)
!   call ebessel_dbb(z, nmin, nmax, b, bp)
!   write (6, "(/,' z=',2(f12.4))" ) z
!   write (6, "( 'b=', 10(2(1pe15.7),5x) )" ) b
!   write(6, "( 'bp=', 10(2(1pe15.7),5x) )" ) bp
!   end do; end do
!   
!   pause
!   end
 
 
 subroutine ebessel_dbb(z, nmin, nmax, ein, einp)
!   calculates exp(-z)*I_n(z) and exp(-z)*I'_n(z) for nmin ² n ² nmax
!   Agrees with Mathematica calculation to 7 places over the range -3 ² n ² 3, 
!   -2 ² re{z}, Im{z} ² 2

    USE Complex_Bessel, only : dp, cbesi

    implicit none
    
    integer, parameter :: n_limit_ebessel = 50

    complex, intent(in) :: z
    integer, intent(in) :: nmin, nmax
    complex, intent(out) :: ein(nmin:nmax), einp(nmin:nmax)

!    complex :: b(n_limit_ebessel+1)
    integer :: ier, ncalc, i, nz

    COMPLEX(dp) :: zdp
    complex(dp) :: b(n_limit_ebessel+1)
!    complex :: b(n_limit_ebessel+1)

    zdp = z 
    ncalc = max(abs(nmin), abs(nmax))+1
    
    if (ncalc > n_limit_ebessel+1) then
     write (0, "( 'ebessel_dbb: order exceeds n_limit_ebessel,  ncalc=', i8 )" ) ncalc
         return
    end if
    
! Test for zero argument

    if (cabs(z) == 0. ) then
    
        do i=0, ncalc
        ein(i)=cmplx(0., 0.)
        einp(i)=cmplx(0., 0.)
    end do
    
    ein(0)=cmplx(1., 0.)
    einp(1)=cmplx(0.5, 0.)
    
    else
    
!  for  nonzero arg call EXPBESIC to get exp(-z)*I_n(z) for 0 < n < calc

        call cbesi(zdp, 0._dp ,2 , ncalc, b, nz, ier)
!        call EXPBESIC (z, ncalc, b, ier)
    
    do i= max(0,nmin), nmax
        ein(i)=b(i+1)
        einp(i)=b(i+2) + real(i)*b(i+1)/z
    end do
    end if
    
! for negative orders fill in by symmetry

    if (nmin < 0) then
    
        do i=nmin, -1
        ein(i) = ein(-i)
        einp(i) = einp(-i)
    end do
    end if
        
    return 
 end subroutine ebessel_dbb
