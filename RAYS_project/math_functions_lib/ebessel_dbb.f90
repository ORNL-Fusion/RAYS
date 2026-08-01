 subroutine ebessel_dbb(z, nmin, nmax, ein, einp)
!   calculates exp(-z)*I_n(z) and exp(-z)*I'_n(z) for nmin ² n ² nmax
!   Agrees with Mathematica calculation to 7 places over the range -3 ² n ² 3,
!   -2 ² re{z}, Im{z} ² 2

    USE Complex_Bessel, only : rkind, cbesi
!     use write_matrix_m, only : write_matrix, write_vector

    implicit none

    integer, parameter :: n_limit_ebessel = 50
    real(KIND=rkind),parameter :: zero = 0.0_rkind, one = 1.0_rkind, two = 2.0_rkind

    complex(rkind), intent(in) :: z
    integer, intent(in) :: nmin, nmax
    complex(rkind), intent(out) :: ein(nmin:nmax), einp(nmin:nmax)

    integer :: ier, ncalc, i, nz

    complex(rkind) :: b(n_limit_ebessel+1)

    b = zero; ein = zero; einp = zero
    ncalc = max(abs(nmin), abs(nmax))+2
    if (ncalc > n_limit_ebessel+1) then
     write (0, "( 'ebessel_dbb: order exceeds n_limit_ebessel,  ncalc=', i8 )" ) ncalc
         return
    end if

! Test for zero argument

    if (abs(z) == zero ) then

!         do i=0, ncalc
!             ein(i)=cmplx(zero, zero)
!             einp(i)=cmplx(zero, zero)
!         end do

        ein(0)=cmplx(one, zero, rkind)
        einp(1)=cmplx(one/two, zero, rkind)

    else ! non-zero

!  for  nonzero arg call EXPBESIC to get exp(-z)*I_n(z) for 0 < n < calc

        call cbesi(z, zero ,2 , ncalc, b, nz, ier)
!        call EXPBESIC (z, ncalc, b, ier)
!  call write_vector('ebessel_dbb: b', b, n_limit_ebessel+1)

        do i= max(0,nmin), nmax
            ein(i)=b(i+1)
            einp(i)=b(i+2) + i*b(i+1)/z
        end do
!  call write_vector('ebessel_dbb: ein', ein, size(ein))
!  call write_vector('ebessel_dbb: einp', einp, size(ein))

    end if

! for negative orders fill in by symmetry

    if (nmin < 0) then
        do i=nmin, -1
        ein(i) = ein(-i)
        einp(i) = einp(-i)
        end do
    end if
!  write(*,*)
!  write(*,*) 'ebessel_dbb: z = ', z
!  write(*,*) 'ebessel_dbb: ein = ', ein

    return
 end subroutine ebessel_dbb
