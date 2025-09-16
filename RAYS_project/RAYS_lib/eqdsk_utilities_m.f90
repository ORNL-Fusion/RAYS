module eqdsk_utilities_m

! Routines for reading, writing, and interpolating geqdsk files.  These were adapted from
! routines written by Richard Fitzpatrick,found in his EPEC github repo, and made into a
! module.  Fitzpatrick's variable names are different from those found in Lang Lau's
! description "G EQSDK FORMAT", 2/7/1997.  Most are obvious, but some are not completely:
! Fitzpatrick's -> Lau's
! T -> fpol == R*Btoroidal
! P -> pres
! TTp -> ffprim
! Pp -> pprime
!
! Fitzpatrick also provides some interpolation routines GetPsi, and derivatives that rely
! on linear or 3 point interpolations on the R and Z grids.  I have adapted these routines
! although I don't expect them to be as accurate as a cubic spline interpolation.  Note,
! I expect the optimum values for dR, dZ to be 1/2 the grid spacing for first derivatives
! and equal to the grid spacing for 2nd derivatives. so I have coded in spacing of 2*dR
! for 2nd derivatives, and anticipate that the module initialization will set dR, dZ to
! the grid spacing.

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

  use constants_m, only : rkind

  implicit none

! eqdsk data
  character (len = 100) :: string
  integer :: i, j, i3, NRBOX, NZBOX, NBOUND, NLIM
  real(KIND=rkind) :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
  real(KIND=rkind), dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM
  real(KIND=rkind), dimension (:, :), allocatable :: Psi

! data needed for interpolation of psi etc
  real(KIND=rkind), dimension (:), allocatable :: R_grid, Z_grid
  real(KIND=rkind) :: dR, dZ

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

! %%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to read gFile
! %%%%%%%%%%%%%%%%%%%%%%%%

subroutine ReadgFile (eqdsk_file_name)


  implicit none

  character (len = *), intent(in) :: eqdsk_file_name

  open (unit = 100, file = trim(eqdsk_file_name), status = 'old')

  read (100, '(a48, 3i4)') string,  i3,      NRBOX,   NZBOX
  read (100, '(5e16.9  )') RBOXLEN, ZBOXLEN, R0,      RBOXLFT,  ZOFF
  read (100, '(5e16.9  )') RAXIS,   ZAXIS,   PSIAXIS, PSIBOUND, B0
  read (100, '(5e16.9  )') CURRENT, zero,    zero,    zero,     zero
  read (100, '(5e16.9  )') zero,    zero,    zero,    zero,     zero

  if (.not. allocated(T)) then
	  allocate (T   (NRBOX))
	  allocate (P   (NRBOX))
	  allocate (TTp (NRBOX))
	  allocate (Pp  (NRBOX))
	  allocate (Q   (NRBOX))
	  allocate (Psi (NRBOX, NZBOX))
  end if

  read (100, '(5e16.9)') (T   (i), i = 1, NRBOX)
  read (100, '(5e16.9)') (P   (i), i = 1, NRBOX)
  read (100, '(5e16.9)') (TTp (i), i = 1, NRBOX)
  read (100, '(5e16.9)') (Pp  (i), i = 1, NRBOX)

  read (100, '(5e16.9)') ((Psi  (i, j), i = 1, NRBOX), j = 1, NZBOX)

  read (100, '(5e16.9)') (Q (i), i = 1, NRBOX)

  read (100, '(2i5)') NBOUND, NLIM

  if (.not. allocated(RBOUND)) then
	  allocate (RBOUND (NBOUND))
	  allocate (ZBOUND (NBOUND))
	  allocate (RLIM   (NLIM))
	  allocate (ZLIM   (NLIM))
  end if

!   write(*,*) 'NBOUND = ', NBOUND, '   NLIM = ', NLIM
!   write(*,*) 'shape(RBOUND) = ', shape(RBOUND)

  read (100, '(5e16.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
!   write(*,*) 'RBOUND = ', RBOUND
!   write(*,*) 'ZBOUND = ', RBOUND

  read (100, '(5e16.9)') (RLIM   (i), ZLIM   (i), i = 1, NLIM)

  close (unit = 100)

end subroutine ReadgFile

! %%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to write gFile
! %%%%%%%%%%%%%%%%%%%%%%%%%

subroutine WritegFile (eqdsk_out)

  implicit none

  character (len = *), intent(in) :: eqdsk_out

  open (unit = 100, file = trim(eqdsk_out), status = 'unknown')

  zero = 0.

  write (100, '(a48, 3i4)') string,  i3,      NRBOX,    NZBOX
  write (100, '(5e16.9  )') RBOXLEN, ZBOXLEN, R0,       RBOXLFT,  ZOFF
  write (100, '(5e16.9  )') RAXIS,   ZAXIS,   PSIAXIS,  PSIBOUND, B0
  write (100, '(5e16.9  )') CURRENT, PSIAXIS, zero,     RAXIS,    zero
  write (100, '(5e16.9  )') ZAXIS,   zero,    PSIBOUND, zero,     zero
  write (100, '(5e16.9)')  (T   (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  (P   (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  (TTp (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  (Pp  (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  ((Psi  (i, j), i = 1, NRBOX), j = 1, NZBOX)
  write (100, '(5e16.9)')  (Q (i), i = 1, NRBOX)
  write (100, '(2i5)')      NBOUND, NLIM
  write (100, '(5e16.9)')  (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  write (100, '(5e16.9)')  (RLIM   (i), ZLIM   (i), i = 1, NLIM)

  close (unit = 100)

end subroutine WritegFile

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated Psi
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetPsi(R, Z)

  implicit none

  real(KIND=rkind)  :: GetPsi

  integer :: i, j

  real(KIND=rkind) :: R, Z, x, y

  i = 1 + int ((R - R_grid (1)) /(R_grid (2) - R_grid (1)))
  j = 1 + int ((Z - Z_grid (1)) /(Z_grid (2) - Z_grid (1)))

  x = (R - R_grid (i)) /(R_grid (2) - R_grid (1))
  y = (Z - Z_grid (j)) /(Z_grid (2) - Z_grid (1))

  GetPsi = Psi (i, j) * (1.-x)*(1.-y) + Psi (i+1, j) * x*(1.-y) + Psi (i, j+1) * (1.-x)*y + Psi(i+1, j+1) * x*y

end function GetPsi

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated T == R*Bphi
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetRBphi(R)

  implicit none

  real(KIND=rkind)  :: GetRBphi

  integer :: i

  real(KIND=rkind) :: R, x

  i = 1 + int ((R - R_grid (1)) /(R_grid (2) - R_grid (1)))

  x = (R - R_grid (i)) /(R_grid (2) - R_grid (1))

  GetRBphi = T(i) * (1.-x) + T(i+1)*x

end function GetRBphi

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated dPsi/dR
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

real(KIND=rkind) function GetPsiR(R, Z)

  use constants_m, only : rkind

!  use Function_Defs_1
  implicit none

  real(KIND=rkind) :: R, Z, R1, R2

  R1 = R - dR
  R2 = R + dR

  GetPsiR = (GetPsi(R2, Z) - GetPsi(R1, Z))/2./dR

end function GetPsiR

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated dPsi/dZ
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

real(KIND=rkind) function GetPsiZ(R, Z)

  use constants_m, only : rkind

  implicit none

  real(KIND=rkind) :: R, Z, Z1, Z2

  Z1 = Z - dZ
  Z2 = Z + dZ

  GetPsiZ = (GetPsi(R, Z2) - GetPsi(R, Z1)) /2./dZ

end function GetPsiZ

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated d^2Psi/dR^2
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

real(KIND=rkind) function GetPsiRR(R, Z)

  use constants_m, only : rkind

!  use Function_Defs_1

  implicit none

  real(KIND=rkind) :: R, Z, R1, R2

  R1 = R - 2.*dR
  R2 = R + 2.*dR

  GetPsiRR = (GetPsi(R2, Z) - 2. * GetPsi(R, Z)&
       + GetPsi(R1, Z)) /dR/dR

end function GetPsiRR

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated d^2Psi/dZ^2
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

real(KIND=rkind) function GetPsiZZ(R, Z)

  use constants_m, only : rkind

  implicit none

  real(KIND=rkind) :: R, Z, Z1, Z2

  Z1 = Z - 2.*dZ
  Z2 = Z + 2.*dZ

  GetPsiZZ = (GetPsi(R, Z2) - 2. * GetPsi(R, Z)&
       +  GetPsi(R, Z1)) /dZ/dZ

end function GetPsiZZ

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated d^2Psi/dRdZ
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

real(KIND=rkind) function GetPsiRZ(R, Z)

  use constants_m, only : rkind

  implicit none

  real(KIND=rkind) :: R, Z, R1, R2, Z1, Z2

  R1 = R - dR
  R2 = R + dR
  Z1 = Z - dZ
  Z2 = Z + dZ

  GetPsiRZ = (GetPsi(R2, Z2) - GetPsi(R1, Z2)&
       -  GetPsi(R2, Z1) + GetPsi(R1, Z1)) /4./dR/dZ

end function GetPsiRZ

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated d(R Bphi)/dR
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

real(KIND=rkind) function GetRBphiR(R)

  use constants_m, only : rkind

  implicit none

  real(KIND=rkind) :: R, R1, R2

  R1 = R - dR
  R2 = R + dR

  GetRBphiR = (GetRBphi(R2) - GetRBphi(R1))/2./dR

end function GetRBphiR

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to find O- and X-points
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine FindOXPoint(R, Z, p)

  use constants_m, only : rkind

  implicit none

  integer :: i

  real(KIND=rkind) :: R, Z, p, pr, pz, prr, pzz, prz, det

  do i = 1, 50

     pr  = GetPsiR  (R, Z)
     pz  = GetPsiZ  (R, Z)
     prr = GetPsiRR (R, Z)
     pzz = GetPsiZZ (R, Z)
     prz = GetPsiRZ (R, Z)

     det = prr*pzz - prz*prz

     R = R + (prz*pz - pzz*pr) /det
     Z = Z + (prz*pr - prr*pz) /det

     p = GetPsi (R, Z)

  enddo

end subroutine FindOXPoint

!********************************************************************

    subroutine deallocate_eqdsk_utilities_m
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
	  deallocate (R_grid)
	  deallocate (Z_grid)
	  return
    end subroutine deallocate_eqdsk_utilities_m

end module eqdsk_utilities_m
