program solovev_2_eqdsk

! A simple code to generate an eqdsk file from a solovev equilibrium.  It links to RAYS\_lib
! uses the solovev\_magnetics routine there, and uses the same solovev\_magnetics\_list
! namelist to define the equilibrium.  It also directly uses eqdsk\_utilities\_m to write
! the eqdsk file.  This routine only provides the eqdsk variables needed to calculate the
! magnetic field quantities used in ray tracing, i.e. those provided by solovev_magnetics().
! Unused eqdsk variables, such as current, pressure, FF', P', and q are set to zero.
!
! N.B. To be compatible with initialize_solovev_magnetics(), which expects a generic input
!      file name to be 'rays.in', that is also the namelist file name for this code, 
!      although a more descriptive name would be 'solovev_2_eqdsk.in'

  use constants_m, only : rkind, input_unit, output_unit
  
  use solovev_magnetics_m, only : initialize_solovev_magnetics, solovev_magnetics, &
        & solovev_magnetics_psi, rmaj, kappa, bphi0, iota0, psiB

  use eqdsk_utilities_m, only : WritegFile, R_grid, Z_grid, dR, dZ, &
        & string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, &
        & PSIAXIS, PSIBOUND, B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, &
        & ZBOUND, RLIM, ZLIM

  
  implicit none

  character (len = 100) :: eqdsk_file_name

! N.B. The variables "box_..."" also appear in module solovev_magnetics_m.  They are
!      declared here because they are needed as arguments in initialize_solovev_magnetics().
!      So we don't access them from use association in the module
  
    ! data for bounding box of computational domain
    real(KIND=rkind) :: box_rmin, box_rmax, box_zmin, box_zmax
    ! data for plasma boundary
    real(KIND=rkind) :: inner_bound, outer_bound, upper_bound, lower_bound

  real(KIND=rkind) :: zero = 0.
  integer :: i,j
  real(KIND=rkind) :: rvec(3) 
  real(KIND=rkind)  :: gradpsi(3), psiN, gradpsiN(3)

  logical :: read_input = .true.

  real(KIND=rkind) :: R, Zsq, dR_bound

  namelist /solovev_2_eqdsk_list/ &
     & eqdsk_file_name, string, &
     & NRBOX, NZBOX, NBOUND


 open(unit=input_unit, file='rays.in', action='read', status='old', form='formatted')
 read(input_unit, solovev_2_eqdsk_list)
 close(unit=input_unit)
 write(*, solovev_2_eqdsk_list)


 call initialize_solovev_magnetics(read_input, RAXIS, ZAXIS, &
	   & box_rmin, box_rmax, box_zmin,box_zmax, &
	   & inner_bound, outer_bound, upper_bound, lower_bound)

 
 write(*,*) 'Inner boundary = ', inner_bound
 write(*,*) 'Outer boundary = ', outer_bound
 write(*,*) 'Lower boundary = ', lower_bound
 write(*,*) 'Upper boundary = ', upper_bound

 RBOXLEN = box_rmax - box_rmin
 ZBOXLEN = box_zmax - box_zmin
 RBOXLFT = box_rmin
 ZOFF = 0.
 PSIAXIS = 0.
 PSIBOUND = psiB
 B0 = bphi0
 
! Unused eqdsk variables
 CURRENT = 0.
 
 allocate (P(NRBOX))
 allocate (TTp(NRBOX))
 allocate (Pp(NRBOX))
 allocate (Q(NRBOX))
 
 P = 0.
 TTp = 0.
 Pp = 0.
 Q = 0. 

 NLIM = 1  ! There is no limiter model so set to 1 to avoid allocation problems
 allocate (RLIM(NLIM))
 allocate (ZLIM(NLIM))
 RLIM = 0.
 ZLIM = 0. 
    
    ! radial and Z grids
    allocate (R_grid(NRBOX))
    allocate (Z_grid(NZBOX))
    ! R*Bphi and psi
    allocate (T(NRBOX))
    allocate (Psi(NRBOX, NZBOX)) 

    do i = 1, NRBOX
        R_grid(i) = box_rmin + (box_rmax - box_rmin)*(i-1)/(NRBOX - 1)
        T(i) = bphi0*rmaj
    end do 

    do j = 1, NZBOX
        Z_grid(j) = box_zmin + (box_zmax - box_zmin)*(j-1)/(NZBOX - 1)
    end do 
 
    do i = 1, NRBOX
       do j = 1, NZBOX
           rvec = (/ R_grid(i), zero, Z_grid(j) /)
           call solovev_magnetics_psi(rvec, Psi(i,j) , gradpsi, psiN, gradpsiN)
       end do
    end do 
    
! Calculate boundary points.  This is cloned from axisym_toroid_processor_m
! N.B.  This is up-down symmetric, and NBOUND needs to be an odd number

 allocate(RBOUND(NBOUND), ZBOUND(NBOUND))

		RBOUND(1) = inner_bound
		ZBOUND(1) = 0.
		RBOUND(NBOUND) = inner_bound
		ZBOUND(NBOUND) = 0.
		RBOUND((NBOUND+1)/2) = outer_bound
		ZBOUND((NBOUND+1)/2) = 0. 

        dR_bound = 2.*(outer_bound - inner_bound)/(NBOUND-1)
		do i = 2, (NBOUND-1)/2
			R = inner_bound + (i-1)*dR_bound
			Zsq = kappa/(4.*R**2)*(outer_bound**4 + 2.*(R**2 - outer_bound**2)*rmaj**2 -&
			      & R**4)
			RBOUND(i) = R
			ZBOUND(i) = sqrt(Zsq)
			RBOUND(NBOUND-(i-1)) = RBOUND(i)
			ZBOUND(NBOUND-(i-1)) = -ZBOUND(i)
		end do 
		
! 		do i = 1, NBOUND
! 			write(*,*) 'i = ', i, '   RBOUND = ', RBOUND(i), '   ZBOUND', ZBOUND(i)
! 		end do
    
  
  write (*, *) ' '
  write (*, *) 'string,  i3,      NRBOX,    NZBOX'
  write (*, '(a48, 3i4)') string,  i3,      NRBOX,    NZBOX
  write (*, *) ' '
  write (*, *) 'RBOXLEN, ZBOXLEN, R0,       RBOXLFT,  ZOFF'
  write (*, '(5e16.9  )') RBOXLEN, ZBOXLEN, R0,       RBOXLFT,  ZOFF
  write (*, *) ' '
  write (*, *) 'RAXIS,   ZAXIS,   PSIAXIS,  PSIBOUND, B0'
  write (*, '(5e16.9  )') RAXIS,   ZAXIS,   PSIAXIS,  PSIBOUND, B0
  write (*, *) ' '
  write (*, *) 'CURRENT, PSIAXIS, zero,     RAXIS,    zero'
  write (*, '(5e16.9  )') CURRENT, PSIAXIS, zero,     RAXIS,    zero
  write (*, *) ' '
  write (*, *) 'ZAXIS,   zero,    PSIBOUND, zero,     zero'
  write (*, '(5e16.9  )') ZAXIS,   zero,    PSIBOUND, zero,     zero
  write (*, *) ' '
  write (*, *)  '(T(i), i = 1, NRBOX)'
  write (*, '(5e16.9)')  (T   (i), i = 1, NRBOX)
  write (*, *) ' '
  write (*, *)  '(P(i), i = 1, NRBOX)'
  write (*, '(5e16.9)')  (P   (i), i = 1, NRBOX)
  write (*, *) ' '
  write (*, *)  '(TTp(i), i = 1, NRBOX)'
  write (*, '(5e16.9)')  (TTp (i), i = 1, NRBOX)
  write (*, *) ' '
  write (*, *)  '(Pp(i), i = 1, NRBOX)'
  write (*, '(5e16.9)')  (Pp  (i), i = 1, NRBOX)
  write (*, *) ' '
  write (*, *)  '((Psi(i, j), i = 1, NRBOX), j = 1, NZBOX)'
  write (*, '(5e16.9)')  ((Psi  (i, j), i = 1, NRBOX), j = 1, NZBOX)
  write (*, *) ' '
  write (*, *)  '(Q (i), i = 1, NRBOX)'
  write (*, '(5e16.9)')  (Q (i), i = 1, NRBOX)
  write (*, *) ' '
  write (*, *)      'NBOUND, NLIM'
  write (*, '(2i5)')      NBOUND, NLIM
  write (*, *) ' '
  write (*, *)  '(RBOUND (i), ZBOUND (i), i = 1, NBOUND)'
  write (*, '(5e16.9)')  (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  write (*, *) ' '
  write (*, *)  'minval(RBOUND)= ', minval(RBOUND), '  maxval(RBOUND)= ', maxval(RBOUND)
  write (*, *)  'minval(ZBOUND)= ', minval(ZBOUND), '  maxval(ZBOUND)= ', maxval(ZBOUND)
  write (*, *) ' '
  write (*, *)  '(RLIM(i), ZLIM(i), i = 1, NLIM)'
  write (*, '(5e16.9)')  (RLIM   (i), ZLIM   (i), i = 1, NLIM)
        
  ! ..................
  ! Write output gFile
  ! ..................

  call WritegFile(eqdsk_file_name)
   
  ! ........
  ! Clean up
  ! ........
  
  deallocate (T)
  deallocate (P)
  deallocate (TTp)
  deallocate (Pp)
  deallocate (Q)
  deallocate (Psi)
  deallocate (RBOUND)
  deallocate (ZBOUND)
  deallocate (RLIM)
  deallocate (ZLIM)
  
end program solovev_2_eqdsk
