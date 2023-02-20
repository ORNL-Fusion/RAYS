program test_uniform_grid_binner_m

    USE bin_to_uniform_grid_m, only : bin_to_uniform_grid

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

    REAL(kind = skind) :: pis = 3.1415926
    REAL(kind = rkind) :: pir = 4._rkind*atan(1._rkind)
    REAL(kind = rkind) :: pi = 3.1415926535897931


    INTEGER, PARAMETER :: nx = 103, n_bins = 5
    REAL(kind = skind) :: Q(nx), x(nx), binned_Q(n_bins)
    REAL(kind = skind) ::  xmin, xmax, delta_x

! Coarse Q grid test
    INTEGER, PARAMETER :: nxc = 10, n_binsc = 30
    REAL(kind = skind) :: Qc(nxc), xc(nxc), binned_Qc(n_binsc)
    REAL(kind = skind) ::  xminc, xmaxc, delta_xc

! Test rkind
    INTEGER, PARAMETER :: nxr = 103, n_binsr = 5
    REAL(kind = rkind) :: Qr(nx), xr(nx), binned_Qr(n_bins)
    REAL(kind = rkind) ::  xminr, xmaxr, delta_xr

    INTEGER :: ix, ib, ierr
    
    xmin = 0.
    xmax = 1.
    delta_x = (xmax - xmin)/real(nx-1)
	write(*,*)
	write(*,*) 'delta_x = ', delta_x
	write(*,*)
    
    do ix = 1, nx
    	x(ix) = delta_x*(ix-1)
    	Q(ix) = x(ix)
    end do
 
!  	write(*,*) 'Q'
! 	write(*,*) Q
   
	call bin_to_uniform_grid(Q, x, xmin, xmax, binned_Q, ierr)

	write(*,*)
	write(*,*) 'binned_Q linear'
	write(*,*) binned_Q
	write(*,*)
	write(*,*) 'sum(binned_Q)'
	write(*,*) sum(binned_Q)
	write(*,*)

!*****************************************************************************
    do ix = 1, nx
    	x(ix) = delta_x*(ix-1)
    	Q(ix) = sin(pis*x(ix))
    end do
 
!  	write(*,*) 'Q'
! 	write(*,*) Q
   
	call bin_to_uniform_grid(Q, x, xmin, xmax, binned_Q, ierr)

	write(*,*)
	write(*,*) 'binned_Q'
	write(*,*) binned_Q
	write(*,*)
	write(*,*) 'sum(binned_Q sin)'
	write(*,*) sum(binned_Q)
	write(*,*)

!*****************************************************************************
! Test coarse Q grid so that a Q step is binned into more than 2 bins
!*****************************************************************************

	write(*,*) '************************************************************************'
	write(*,*) 'Test coarse Q grid'
	write(*,*) '************************************************************************'
	write(*,*)
    
    xminc = 0.
    xmaxc = 1.
    delta_xc = (xmaxc - xminc)/real(nxc-1)
	write(*,*)
	write(*,*) 'delta_xc = ', delta_xc
	write(*,*)
    
    do ix = 1, nxc
    	xc(ix) = delta_xc*(ix-1)
    	Qc(ix) = xc(ix)
    end do
 
  	write(*,*) 'Qc'
 	write(*,*) Qc
   
	call bin_to_uniform_grid(Qc, xc, xminc, xmaxc, binned_Qc, ierr)

	write(*,*)
	write(*,*) 'binned_Qc linear'
	write(*,*) binned_Qc
	write(*,*)
	write(*,*) 'sum(binned_Qc)'
	write(*,*) sum(binned_Qc)
	write(*,*)

!*****************************************************************************

    do ix = 1, nxc
    	xc(ix) = delta_xc*(ix-1)
    	Qc(ix) = sin(pis*xc(ix))
    end do
 
!  	write(*,*) 'Q'
! 	write(*,*) Q
   
	call bin_to_uniform_grid(Qc, xc, xminc, xmaxc, binned_Qc, ierr)

	write(*,*)
	write(*,*) 'binned_Qc sin'
	write(*,*) binned_Qc
	write(*,*)
	write(*,*) 'sum(binned_Qc sin)'
	write(*,*) sum(binned_Qc)
	write(*,*)

!*****************************************************************************
! Test kind = selected_real_kind(15,307)
!*****************************************************************************

	write(*,*) '************************************************************************'
	write(*,*) 'Test rkind = selected_real_kind(15,307)'
	write(*,*) '************************************************************************'
	write(*,*)
    
    xminr = 0._rkind
    xmaxr = 1._rkind
    delta_xr = (xmaxr - xminr)/real(nxr-1)
	write(*,*)
	write(*,*) 'delta_xr = ', delta_xr
	write(*,*)
    
    do ix = 1, nxr
    	xr(ix) = delta_xr*(ix-1)
    	Qr(ix) = sin(pir*xr(ix))
    end do
 
  	write(*,*) 'Qr'
 	write(*,*) Qr
   
	call bin_to_uniform_grid(Qr, xr, xminr, xmaxr, binned_Qr, ierr)

	write(*,*)
	write(*,*) 'binned_Qr sin'
	write(*,*) binned_Qr
	write(*,*)
	write(*,*) 'sum(binned_Qr)'
	write(*,*) sum(binned_Qr)
	write(*,*)

 write (*,*) 'sin(pis) = ', sin(pis), '    4.*atan(1.) = ', 4.*atan(1.)
 write (*,*) 'sin(pi) = ', sin(pi), '    pi = ', pi
 write (*,*) 'sin(pir) = ', sin(pir), '    4._rkind*atan(1._rkind) = ', 4._rkind*atan(1._rkind)
 write (*,*) 'sin(4._rkind*atan(1._rkind)) = ', sin(4._rkind*atan(1._rkind))
end program test_uniform_grid_binner_m
