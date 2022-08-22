program test_gFileReadWrite_Fitzpatrick

  use constants_m, only : rkind
  use eqdsk_utilities_m, only : ReadgFile, WritegFile, R_grid, Z_grid, dR, dZ, &
        & string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, &
        & PSIAXIS, PSIBOUND, B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, &
        & ZBOUND, RLIM, ZLIM, &
        & GetPsi, GetRBphi, GetPsiR, GetPsiZ, GetPsiRR, GetPsiZZ, GetPsiRZ, GetRBphiR, &
        & psi_derivs_lin_interp

  
  implicit none

  character (len = 100) :: eqdsk_in = 'eqdsk_in.dat'
  character (len = 100) :: eqdsk_out = 'eqdsk_out.dat'
  
    ! data for bounding box of computational domain
    real(KIND=rkind) :: box_rmin, box_rmax, box_zmin, box_zmax
    ! data for plasma boundary
    real(KIND=rkind) :: inner_bound, outer_bound, upper_bound, lower_bound

  real(KIND=rkind) :: zero = 0.
  integer :: i,j
  real(KIND=rkind) :: R, Z
  real(KIND=rkind) :: PsiR, PsiZ, PsiRR, PsiZZ, PsiRZ, RBphiR 
  
  ! ................
  ! Read input gFile
  ! ................

  call ReadgFile(eqdsk_in)
  
  write (*, '(a48, 3i4)') string,  i3,      NRBOX,    NZBOX
  write (*, '(5e16.9  )') RBOXLEN, ZBOXLEN, R0,       RBOXLFT,  ZOFF
  write (*, '(5e16.9  )') RAXIS,   ZAXIS,   PSIAXIS,  PSIBOUND, B0
  write (*, '(5e16.9  )') CURRENT, PSIAXIS, zero,     RAXIS,    zero
  write (*, '(5e16.9  )') ZAXIS,   zero,    PSIBOUND, zero,     zero
  write (*, '(5e16.9)')  (T   (i), i = 1, NRBOX)
  write (*, '(5e16.9)')  (P   (i), i = 1, NRBOX)
  write (*, '(5e16.9)')  (TTp (i), i = 1, NRBOX)
  write (*, '(5e16.9)')  (Pp  (i), i = 1, NRBOX)
  write (*, '(5e16.9)')  ((Psi  (i, j), i = 1, NRBOX), j = 1, NZBOX)
  write (*, '(5e16.9)')  (Q (i), i = 1, NRBOX)
  write (*, '(2i5)')      NBOUND, NLIM
  write (*, '(5e16.9)')  (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  write (*, '(5e16.9)')  (RLIM   (i), ZLIM   (i), i = 1, NLIM)
 
    box_rmin = RBOXLFT
    box_rmax = box_rmin + RBOXLEN
    box_zmin = ZOFF - ZBOXLEN/2.
    box_zmax = ZOFF + ZBOXLEN/2.

    inner_bound = minval(RBOUND)
    outer_bound = maxval(RBOUND)
    lower_bound = minval(ZBOUND)
    upper_bound = maxval(ZBOUND)

    write(*,*) 'Inner boundary = ', inner_bound
    write(*,*) 'Outer boundary = ', outer_bound
    write(*,*) 'Lower boundary = ', lower_bound
    write(*,*) 'Upper boundary = ', upper_bound
    
    ! radial and Z grids
    allocate (R_grid(NRBOX))
    allocate (Z_grid(NZBOX))

    do i = 1, NRBOX
        R_grid(i) = inner_bound + (outer_bound - inner_bound)*(i-1)/(NRBOX - 1)
    end do

    do i = 1, NZBOX
        Z_grid(i) = lower_bound + (upper_bound - lower_bound)*(i-1)/(NZBOX - 1)
    end do
    
    dr = (R_grid(2) - R_grid(1))/2.
    dz = (Z_grid(2) - Z_grid(1))/2.
    R = (outer_bound - inner_bound)/2.
    Z = (upper_bound - abs(lower_bound))/2.

    write(*,*) 'dr = ', dr, '   dz = ', dz, '   R = ', R, '   Z = ', Z
    write(*,*) 'GetPsi(R, Z) = ', GetPsi(R, Z), '   GetRBphi(R) = ', GetRBphi(R)
    write(*,*) 'GetPsiR(R, Z) = ', GetPsiR(R, Z)
    write(*,*) 'GetPsiZ(R, Z) = ', GetPsiZ(R, Z)
    write(*,*) 'GetPsiRR(R, Z) = ', GetPsiRR(R, Z)
    write(*,*) 'GetPsiZZ(R, Z) = ', GetPsiZZ(R, Z)
    write(*,*) 'GetPsiRZ(R, Z) = ', GetPsiRZ(R, Z)
    write(*,*) 'GetRBphiR(R) = ', GetRBphiR(R)
    
    write(*,*) ' '
    write(*,*) 'Derivatives from psi_derivs_lin_interp() '
    
    call psi_derivs_lin_interp(R, Z, PsiR, PsiZ, PsiRR, PsiZZ, PsiRZ, RBphiR)
    
    write(*,*) 'PsiR = ', PsiR
    write(*,*) 'PsiZ = ', PsiZ
    write(*,*) 'PsiRR = ', PsiRR
    write(*,*) 'PsiZZ = ', PsiZZ
    write(*,*) 'PsiRZ = ', PsiRZ
    write(*,*) 'RBphiR(R) = ', RBphiR
    
    
        
  ! ..................
  ! Write output gFile
  ! ..................

  call WritegFile(eqdsk_out)
   
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
  
end program test_gFileReadWrite_Fitzpatrick
