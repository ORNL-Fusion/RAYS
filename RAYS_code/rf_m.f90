 module rf_m
 
    implicit none

!   Name of dispersion model used in ray tracing
    character(len=15) :: ray_dispersion_model

!   RF wave frequency: omgrf = 2*pi*frf, frf in Hz, and wave number
!   in vacuum: k0 = clight / omgrf
    real :: omgrf, frf, k0
        
!   kvec (nvec) = k (k/k0) in xyz coordinates (vector):
!   k1 (n1) = perpedicular component of kvec (nvec) (magnitude).
!   k3 (n3) = parallel component of kvec (nvec) (magnitude).
    real :: kvec(3), k1, k3
    real :: nvec(3), n1, n3  
 
!   A switch to select which root from the dispersion to be used for
!   k initialization: i.e., fast wave, slow wave, IBW, or KAW, etc. 
!   Cold plasma gives quadratic for nx^2(nz).  So can choose wave mode = plus, minus
!   fast or slow, where fast is the smaller of |nx|
    character(len=15) :: wave_mode

!   A switch to determine whether the rays are initialized to propagate into
!   or out of the plasma.  Typically k_init_sign = 1 for waves to propagate inward
!   (unless the wave is a backward wave).
    integer:: k0_sign

!   ray_param = 'arcl': default. Integrate with respect to the arclength along the ray.
!   ray_param = 'time': Integrate with respect to time along the ray.
    character(len = 4) :: ray_param = 'arcl'

!   Maximum allowable residual of dispersion function. Used in checksave()
    real :: dispersion_resid_limit
    
    namelist /rf_list/ ray_dispersion_model, frf, wave_mode, k0_sign, ray_param, &
                     & dispersion_resid_limit
    
!********************************************************************

contains

!********************************************************************

  subroutine initialize_rf
 
    use constants_m, only : input_unit, pi, clight
    use diagnostics_m, only : message_unit, message
    
    implicit none

! Read and write input namelist
    open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
    read(input_unit, rf_list)
    close(unit=input_unit)
    write(message_unit, rf_list)

    if ( frf <= 0. ) then
       write (0,*) 'initialize_rf: frf =', frf
       stop 1
    else
       omgrf = 2.*pi * frf
    end if

!   Wave number in vacuum.
    k0 = omgrf/clight 
    call message ('initialize_rf: k0', k0, 1)

    return
  end subroutine initialize_rf

 end module rf_m
