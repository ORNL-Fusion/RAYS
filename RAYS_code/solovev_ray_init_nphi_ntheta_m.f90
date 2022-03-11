 module solovev_ray_init_nphi_ntheta_m
 
! Parameters to generate the initial conditions for the rays
! i.e. intitial position for each ray: rvec0(1:3, iray) and 
! information necessary to initialize k for each ray: kvec0(1:3, iray)
! Note: initialization of k requires solution of the dispersion relation
!   and depends on magnetic field and density information.  The actual
!   initialization of k occurs in subroutine 'initialize' which is called from 
!   subroutine 'ray_tracing'.   Here the information necessary to do this
!   calculation (e.g. toroidal, poloidal mode number, toroidal, poloidal wave
!   number) is generated from data on the number of launch angles etc. and put
!   into rvec0, kvec0.  The values for rvec0, kvec0 are generated in the module
!   routines 'ray_init_...'

!   There are several reasonable ways to initialize k values,
!   particularly in toroidal geometry.  At low frequencies toroidal and poloidal mode
!   numbers are meaningful.  At higher frequencies the launch angle is more useful.
!   For ECH one launches a beam with a distribution of rays around a central beam axis. 
!   The variable ray_init_model determines which of several k intialization schemes are used

    use constants_m, only : rkind
    
    implicit none   

! Data for initial launch position
! N.B. r_launch <--> minor radius, theta is relative to equatorial plane   
    integer:: n_r_launch = 1
    real(KIND=rkind) ::  r_launch0 = 0., dr_launch = 0.
    integer:: n_theta_launch = 1
    real(KIND=rkind) ::  theta_launch0 = 0., dtheta_launch = 0.

! Data initial launch refractive index: theta = poloidal direction, phi = toroidal direction
    integer:: n_rindex_theta = 1
    real(KIND=rkind) ::  rindex_theta0 = 0., delta_rindex_theta = 0.
    integer:: n_rindex_phi = 1
    real(KIND=rkind) ::  rindex_phi0 = 0., delta_rindex_phi = 0.
     
 namelist /solovev_ray_init_nphi_ktheta_list/ &
     & n_r_launch,  r_launch0, dr_launch, dr_launch, &
     & n_theta_launch, theta_launch0, dtheta_launch, &
     & n_rindex_theta, rindex_theta0, delta_rindex_theta, &
     & n_rindex_phi, rindex_phi0, delta_rindex_phi

! Set up the initial conditions (before root finding of the dispersion relation)
! for each ray 

!****************************************************************************


contains


!****************************************************************************

                                             
    subroutine ray_init_solovev_nphi_ntheta(nray_max, nray, rvec0, rindex_vec0)
    
! Specify a range of initial launch positions (r,theta) and at each launch position
! an initial poloidal wave number, rindex_theta0, and toroidal wave number, rindex_phi0.
!
! N.B. Some of the ray initializations may fail (e.g. initial point is outside plasma or 
!      wave mode is evanescent).  This does not cause the program to stop.  It counts
!      the successful initializations and sets number of rays, nray, to that.

    use constants_m, only : input_unit
    use diagnostics_m, only: message_unit, message, text_message
    use species_m, only : nspec
    use equilibrium_m, only : equilibrium, eq_point
    use dispersion_solvers_m, only: solve_n1_vs_n2_n3
    use rf_m, only : ray_dispersion_model, wave_mode, k0_sign, k0
    use solovev_eq_m, only: rmaj, solovev_psi

    implicit none
    
    integer, intent(in) :: nray_max
    integer, intent(out) :: nray
    type(eq_point(nspec=nspec)) :: eq
    real(KIND=rkind), allocatable, intent(out) :: rvec0(:, :), rindex_vec0(:, :)
    
    integer :: iray, itheta, i_ntheta, i_nphi, count
    real(KIND=rkind) :: x, z, rmin_launch, theta, rindex_theta, rindex_phi, n2, n3
    real(KIND=rkind) :: psi, gradpsi(3), psiN, gradpsiN(3)
    real(KIND=rkind) :: rvec(3), phi_unit(3), psi_unit(3), theta_unit(3), trans_unit(3)
    real(KIND=rkind) :: rindex_vec(3)
    complex(KIND=rkind) :: npsi_cmplx
    real(KIND=rkind) :: npsi
    real(KIND=rkind) :: nperp

! Read and write input namelist
    open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
    read(input_unit, solovev_ray_init_nphi_ktheta_list)
    close(unit=input_unit)
    write(message_unit, solovev_ray_init_nphi_ktheta_list)

! Allocate maximum space for the initial condition vectors rvec0, rindex_vec0
! N.B. Not all of these may successfully initialize because of errors.  So count successful
!      initializations.  That's the final value of nray.
    
        nray = n_r_launch * n_theta_launch * n_rindex_theta * n_rindex_theta

        if ((nray > 0) .and. (nray <= nray_max)) then
        allocate ( rvec0(3, nray), rindex_vec0(3, nray) )
        else
        write (*,*) 'solovev ray init: improper number of rays  nray=', nray
        stop 1
        end if  

    count = 0

    ray_loop: do iray = 1, n_r_launch

    theta_loop: do itheta = 1, n_theta_launch
        theta = theta_launch0 + (itheta-1) * dtheta_launch

        rmin_launch = r_launch0+(iray-1)*dr_launch 
        x = rmaj + rmin_launch*cos(theta)
        z = rmin_launch*sin(theta)
        rvec(:) = (/ x, 0._rkind, z /)  ! Without loss of generality take y = 0

    rindex_theta_loop: do i_ntheta = 1, n_rindex_theta 
        rindex_theta = rindex_theta0 + (i_ntheta-1) * delta_rindex_theta 

    rindex_phi_loop: do i_nphi = 1, n_rindex_phi
        rindex_phi = rindex_phi0 + (i_nphi-1) * delta_rindex_phi

        call equilibrium(rvec, eq)
           if (trim(eq%equib_err) /= '') cycle rindex_phi_loop

        call solovev_psi(rvec, psi, gradpsi, psiN, gradpsiN)  
                     
! Calculate a bunch of unit vectors, parallel and transverse components of k
        psi_unit = gradpsi/sqrt(dot_product(gradpsi,gradpsi))

        phi_unit = (/ 0._rkind, 1._rkind, 0._rkind/)
        
        theta_unit = (/-gradpsi(3), 0._rkind, gradpsi(1) /) ! psi_unit X phi_unit
        theta_unit = theta_unit/sqrt(dot_product(theta_unit,theta_unit))        

        ! bunit X psi_unit
        trans_unit = (/eq%bunit(2)*psi_unit(3)-eq%bunit(3)*psi_unit(2), &
                    & eq%bunit(3)*psi_unit(1)-eq%bunit(1)*psi_unit(3), &
                    & eq%bunit(1)*psi_unit(2)-eq%bunit(2)*psi_unit(1)  /)
        
        ! Refreactive index vector projected onto flux surface (i.e. no psi component)
        rindex_vec = rindex_phi*phi_unit + rindex_theta*theta_unit

        n3 = dot_product(eq%bunit, rindex_vec) ! Parallel component
        n2 = dot_product(trans_unit, rindex_vec) ! Transverse component

! Solve dispersion for complex refractive index in psi direction then cast as real        
        call solve_n1_vs_n2_n3(eq, ray_dispersion_model, wave_mode, k0_sign, &
             &  n2, n3, npsi_cmplx)

        if (npsi_cmplx%im /= 0.) then
            write(message_unit, *) 'toroid_ray_init: evanescent ray, rvec = ', rvec, &
            & ' n2 = ', n2, ' n3 = ', n3
            cycle rindex_phi_loop
        end if
        npsi = npsi_cmplx%re

        count = count +1
        rvec0( : , count) = rvec
        
        ! Take psi component to point inward i.e. direction -grad_psi
        rindex_vec0( : , count) = rindex_vec - npsi*psi_unit

 write(*,*) 'wave_mode = ', wave_mode
 write(*,*) 'eq%bunit = ', eq%bunit
 write(*,*) 'psi_unit = ', psi_unit, '  phi_unit = ', phi_unit
 write(*,*) 'theta_unit = ', theta_unit, '  trans_unit = ', trans_unit
 write(*,*) 'rindex_vec = ', rindex_vec, '  n3 = ', n3, '  n2 = ', n2
 write(*,*) 'npsi = ', npsi, '  n3 = ', n3, '  n2 = ', n2
 write(*,*) 'rindex_vec0( : , count) = ', rindex_vec0( : , count)

 nperp = sqrt(npsi**2+n2**2)
 write(*,*) 'nperp = ', nperp
! write(*,*) 'residual = ', residual(eq, k0*nperp, k0*n3)
            
    end do rindex_phi_loop
    end do rindex_theta_loop
    end do theta_loop
    end do ray_loop

    nray = count
    call message('toroid_ray_init: nray', nray)
    if (nray == 0) stop 'No successful ray initializations' 

    end subroutine ray_init_solovev_nphi_ntheta

!****************************************************************************
    real(KIND=rkind) function residual_2(eq, k1, k3)
!      calculates the residual for given k1 and k3.
!      get dielectric tensor from module suscep_m

       use species_m, only : nspec
       use suscep_m, only :  dielectric_cold

       use equilibrium_m, only : eq_point
       use rf_m, only : ray_dispersion_model, k0

       implicit none 
       
       type(eq_point(nspec=nspec)), intent(in) :: eq
 
       real(KIND=rkind) :: k1, k3

       complex(KIND=rkind) :: eps(3,3), eps_h(3,3), epsn(3,3), ctmp
       complex(KIND=rkind) :: eps_norm(3,3)
       real(KIND=rkind) :: n(3)

       integer :: i, j

 write(*,*) 'inside function residual'

!   Need dielectric tensor.

    if (ray_dispersion_model == "cold") then
 write(*,*) 'call dielectric_cold'
        call dielectric_cold(eq, eps)
 write(*,*) ' return from call dielectric_cold'
    end if
    
!    write(*,*) 'eps = ', eps

!      Hermitian part.
       eps_h = .5 * ( eps + conjg(transpose(eps)) )

!      Refractive index.
       n(1) = k1/k0; n(2) = 0.; n(3) = k3/k0

!      epsn = eps + nn -n^2I:
!      epsn.E = eps.E + n x n x E = (eps + nn -n^2I).E,
!      where E = (Ex,Ey,Ez)^T and I is the unit 3X3 tensor.

       do i = 1, 3; do j = 1, 3
          epsn(i,j) = eps_h(i,j) + n(i)*n(j) - int(i/j)*int(j/i)*sum(n**2)
          eps_norm(i,j) = abs( eps_h(i,j) ) + abs(n(i)*n(j))
       end do; end do

!      Determinant for 3X3 epsn:
       ctmp = &
          &   epsn(3,3)*(epsn(1,1)*epsn(2,2)-epsn(2,1)*epsn(1,2)) &
          & - epsn(3,2)*(epsn(1,1)*epsn(2,3)-epsn(2,1)*epsn(1,3)) &
          & + epsn(3,1)*(epsn(1,2)*epsn(2,3)-epsn(2,2)*epsn(1,3))

!      For a Hermitian matrix, the imaginary part of its determinant vanishes.
       if ( abs(aimag(ctmp)) > 1.e-6 ) then
          write(0,'(a,1p1e12.4)') 'RESIDUAL: Im(det) = ', aimag(ctmp)
          stop 1
       end if

       residual_2 = abs(ctmp) / &
          & ( eps_norm(3,3)*(eps_norm(1,1)*eps_norm(2,2)) &
          & + eps_norm(3,3)*(eps_norm(2,1)*eps_norm(1,2)) &
          & + eps_norm(3,2)*(eps_norm(1,1)*eps_norm(2,3)) &
          & + eps_norm(3,2)*(eps_norm(2,1)*eps_norm(1,3)) &
          & + eps_norm(3,1)*(eps_norm(1,2)*eps_norm(2,3)) &
          & + eps_norm(3,1)*(eps_norm(2,2)*eps_norm(1,3)) ) 
      

       return
    end function residual_2


end module solovev_ray_init_nphi_ntheta_m
