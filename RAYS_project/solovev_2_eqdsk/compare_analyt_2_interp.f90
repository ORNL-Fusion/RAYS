program compare_analyt_2_interp

! A code to compare flux function and its derivatives from an analytic model to ones 
! calculated by interpolation of an eqdsk file
!
! N.B. To be compatible with initialize_solovev_magnetics(), which expects a generic input
!      file name to be 'rays.in', that is also the namelist file name for this code, 
!      although a more descriptive name would be 'solovev_2_eqdsk.in'

  use constants_m, only : rkind, input_unit, output_unit
  
  use solovev_magnetics_m, only : initialize_solovev_magnetics, solovev_magnetics, &
        & solovev_magnetics_psi, rmaj, kappa, bphi0, iota0, psiB

  use eqdsk_magnetics_lin_interp_m, only : initialize_eqdsk_magnetics_lin_interp, &
       & eqdsk_magnetics_lin_interp_psi
  
  implicit none

  character (len = 100) :: eqdsk_file_name

    ! Magnetic axis
   real(KIND=rkind) :: r_axis, z_axis
   real(KIND=rkind) :: eqd_r_axis, eqd_z_axis ! from eqdsk

! N.B. The variables "box_..."" also appear in module solovev_magnetics_m.  They are
!      declared here because they are needed as arguments in initialize_solovev_magnetics().
!      So we don't access them from use association in the module
  
    ! data for bounding box of computational domain
    real(KIND=rkind) :: box_rmin, box_rmax, box_zmin, box_zmax
    real(KIND=rkind) :: eqd_box_rmin, eqd_box_rmax, eqd_box_zmin, eqd_box_zmax ! from eqdsk
    ! data for plasma boundary
    real(KIND=rkind) :: inner_bound, outer_bound, upper_bound, lower_bound
    real(KIND=rkind) :: eqd_inner_bound, eqd_outer_bound, eqd_upper_bound, eqd_lower_bound !from eqdsk

  integer :: i,j
  real(KIND=rkind) :: rvec(3) 
  real(KIND=rkind)  :: psi, gradpsi(3), psiN, gradpsiN(3)
  real(KIND=rkind)  :: eqd_psi, eqd_gradpsi(3), eqd_psiN, eqd_gradpsiN(3) ! from eqdsk

  logical :: read_input = .true.
  logical :: inconsitent = .false.

! Test grid for comparison
  integer :: n_Rgrid, n_Zgrid
  real(KIND=rkind), dimension (:), allocatable :: R_grid, Z_grid 
  real(KIND=rkind) :: dR, dZ    

  real(KIND=rkind) :: abserr, relerr
  real(KIND=rkind) :: geom_tol = 1.e-4

  namelist /compare_analyt_2_interp_list/ &
     & n_Rgrid, n_Zgrid


 open(unit=input_unit, file='rays.in', action='read', status='old', form='formatted')
 read(input_unit, compare_analyt_2_interp_list)
 close(unit=input_unit)
 write(*, compare_analyt_2_interp_list)


 call initialize_solovev_magnetics(read_input, r_axis, z_axis, &
	   & box_rmin, box_rmax, box_zmin,box_zmax, &
	   & inner_bound, outer_bound, upper_bound, lower_bound)

 call initialize_eqdsk_magnetics_lin_interp(read_input, eqd_r_axis, eqd_z_axis, &
               & eqd_box_rmin, eqd_box_rmax, eqd_box_zmin, eqd_box_zmax, &
               & eqd_inner_bound, eqd_outer_bound, eqd_upper_bound, eqd_lower_bound)
 
! Check consistency of analytic versus eqdsk inputs
 
 if (abs(r_axis - eqd_r_axis) > geom_tol) then
     write (*,*) 'r_axis inconsistent, r_axis = ', r_axis, '  eqd_r_axis = ', eqd_r_axis
     write(*,*) 'r_axis = ', r_axis
     inconsitent = .true.
 end if
 if (abs(z_axis - eqd_z_axis) > geom_tol) then
     write (*,*) 'z_axis inconsistent, axis = ', z_axis, '  eqd_z_axis = ', eqd_z_axis
     inconsitent = .true.
 end if
 if (abs(box_rmin - eqd_box_rmin) > geom_tol) then
     write (*,*) 'box_rmin inconsistent, box_rmin = ', box_rmin, '  eqd_box_rmin = ', eqd_box_rmin
     inconsitent = .true.
 end if
 if (abs(box_rmax - eqd_box_rmax) > geom_tol) then
     write (*,*) 'box_rmax inconsistent, box_rmax = ', box_rmax, '  eqd_box_rmax = ', eqd_box_rmax
     inconsitent = .true.
 end if
 if (abs(box_zmin - eqd_box_zmin) > geom_tol) then
     write (*,*) 'box_zmin inconsistent, box_zmin = ', box_zmin, '  eqd_box_zmin = ', eqd_box_zmin
     inconsitent = .true.
 end if
 if (abs(box_zmax - eqd_box_zmax) > geom_tol) then
     write (*,*) 'box_rmax inconsistent, box_zmax = ', box_zmax, '  eqd_box_zmax = ', eqd_box_zmax
     inconsitent = .true.     
 end if
 if (abs(inner_bound - eqd_inner_bound) > geom_tol) then
     write (*,*) 'inner_bound inconsistent, nner_bound = ', inner_bound, '  eqd_inner_bound = ', eqd_inner_bound
     inconsitent = .true.
 end if
 if (abs(outer_bound - eqd_outer_bound) > geom_tol) then
     write (*,*) 'outer_bound inconsistent, outer_bound = ', outer_bound, '  eqd_outer_bound = ', eqd_outer_bound
     inconsitent = .true.
 end if
 if (abs(upper_bound - eqd_upper_bound) > geom_tol) then
     write (*,*) 'upper_bound inconsistent, upper_bound = ', upper_bound, '  eqd_upper_bound = ', eqd_upper_bound
     inconsitent = .true.
 end if
 if (abs(lower_bound - eqd_lower_bound) > geom_tol) then
     write (*,*) 'lower_bound inconsistent, lower_bound = ', lower_bound, '  eqd_lower_bound = ', eqd_lower_bound
     inconsitent = .true.
 end if
 if (inconsitent .eqv. .true.) stop
     
 write(*,*) 'Inner boundary = ', inner_bound
 write(*,*) 'Outer boundary = ', outer_bound
 write(*,*) 'Lower boundary = ', lower_bound
 write(*,*) 'Upper boundary = ', upper_bound

! set up Rgrid, Zgrid
 allocate(R_grid(n_Rgrid))
 
 dr = (box_rmax - box_rmin)/(n_Rgrid-1)
 do i = 1, n_Rgrid
     R_grid(i) = box_rmin + (i-1)*dR
 end do
 
 ! For now just have to Z values - zero and halfway up plasma
 n_Zgrid = 2
 allocate(Z_grid(n_Rgrid))
 Z_grid(1) = 0.
 Z_grid(2) = upper_bound/2.
 
! Come back and allocate
 write (*,*) ' '

 do j = 1, n_Zgrid
     write(*,*) 'Z = ', Z_grid(j)
     write(*,*) 'R(i)            psi        eqd_psi        abserr      relerr'
     do i = 1, n_Rgrid
         rvec = (/ R_grid(i), 0._rkind, Z_grid(j) /)
         call solovev_magnetics_psi(rvec, psi, gradpsi, psiN, gradpsiN)
         call eqdsk_magnetics_lin_interp_psi(rvec, eqd_psi, eqd_gradpsi, eqd_psiN, eqd_gradpsiN)
         call abserr_relerr(psi, eqd_psi, abserr, relerr)
         write(*,*) R_grid(i), psi, eqd_psi, abserr, relerr
     end do
     write (*,*) ' '
    write (*,*) '***************************************************************** '
    write (*,*) ' '
 end do    

 do j = 1, n_Zgrid
     write(*,*) 'Z = ', Z_grid(j)
     write(*,*) 'R(i)            psiN        eqd_psiN        abserr      relerr'
     do i = 1, n_Rgrid
         rvec = (/ R_grid(i), 0._rkind, Z_grid(j) /)
         call solovev_magnetics_psi(rvec, psi, gradpsi, psiN, gradpsiN)
         call eqdsk_magnetics_lin_interp_psi(rvec, eqd_psi, eqd_gradpsi, eqd_psiN, eqd_gradpsiN)
         call abserr_relerr(psiN, eqd_psiN, abserr, relerr)
         write(*,*) R_grid(i), psiN, eqd_psiN, abserr, relerr
     end do
     write (*,*) ' '
    write (*,*) '***************************************************************** '
    write (*,*) ' '
 end do    

 do j = 1, n_Zgrid
     write(*,*) 'Z = ', Z_grid(j)
     write(*,*) 'R(i)            gradpsi(1)    eqd_gradpsi(1)        abserr      relerr'
     do i = 1, n_Rgrid
         rvec = (/ R_grid(i), 0._rkind, Z_grid(j) /)
         call solovev_magnetics_psi(rvec, psi, gradpsi, psiN, gradpsiN)
         call eqdsk_magnetics_lin_interp_psi(rvec, eqd_psi, eqd_gradpsi, eqd_psiN, eqd_gradpsiN)
         call abserr_relerr(psi, eqd_psi, abserr, relerr)
         write(*,*) R_grid(i), gradpsi(1), eqd_gradpsi(1), abserr, relerr
     end do
     write (*,*) ' '
    write (*,*) '***************************************************************** '
    write (*,*) ' '
 end do    

 do j = 1, n_Zgrid
     write(*,*) 'Z = ', Z_grid(j)
     write(*,*) 'R(i)            gradpsi(2)    eqd_gradpsi(2)        abserr      relerr'
     do i = 1, n_Rgrid
         rvec = (/ R_grid(i), 0._rkind, Z_grid(j) /)
         call solovev_magnetics_psi(rvec, psi, gradpsi, psiN, gradpsiN)
         call eqdsk_magnetics_lin_interp_psi(rvec, eqd_psi, eqd_gradpsi, eqd_psiN, eqd_gradpsiN)
         call abserr_relerr(psi, eqd_psi, abserr, relerr)
         write(*,*) R_grid(i), gradpsi(2), eqd_gradpsi(2), abserr, relerr
     end do
     write (*,*) ' '
    write (*,*) '***************************************************************** '
    write (*,*) ' '
 end do    
    
 do j = 1, n_Zgrid
     write(*,*) 'Z = ', Z_grid(j)
     write(*,*) 'R(i)            gradpsi(3)    eqd_gradpsi(3)        abserr      relerr'
     do i = 1, n_Rgrid
         rvec = (/ R_grid(i), 0._rkind, Z_grid(j) /)
         call solovev_magnetics_psi(rvec, psi, gradpsi, psiN, gradpsiN)
         call eqdsk_magnetics_lin_interp_psi(rvec, eqd_psi, eqd_gradpsi, eqd_psiN, eqd_gradpsiN)
         call abserr_relerr(psi, eqd_psi, abserr, relerr)
         write(*,*) R_grid(i), gradpsi(3), eqd_gradpsi(3), abserr, relerr
     end do
     write (*,*) ' '
    write (*,*) '***************************************************************** '
    write (*,*) ' '
 end do
 
 contains
 
  subroutine abserr_relerr(x, x1, abserr, relerr)
	
	  use constants_m, only : rkind

	  implicit none
  
	  real(KIND=rkind), intent(in) :: x, x1
	  real(KIND=rkind), intent(out) :: abserr, relerr
  
	  abserr = x1 - x
  
	  if (abs(x) > 10.*tiny(x)) then
		relerr = abserr/abs(x)
	  else
		relerr = huge(x)
	  end if
  end subroutine abserr_relerr
     
end program compare_analyt_2_interp
