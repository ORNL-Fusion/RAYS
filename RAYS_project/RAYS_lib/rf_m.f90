 module rf_m

    use constants_m, only : rkind
     
    implicit none

!   Name of dispersion model used in ray tracing
    character(len=60) :: ray_dispersion_model

!   RF wave frequency: omgrf = 2*pi*frf, frf in Hz, and wave number
!   in vacuum: k0 = clight / omgrf
    real(KIND=rkind) :: omgrf, frf, k0
 
!   A switch to select which root from the dispersion to be used for
!   k initialization: i.e., fast wave, slow wave, IBW, or KAW, etc. 
!   Cold plasma gives quadratic for nx^2(nz).  So can choose wave mode = plus, minus
!   fast or slow, where fast is the smaller of |nx|
    character(len=60) :: wave_mode

!   A switch to determine whether the rays are initialized to propagate into
!   or out of the plasma.  Typically k_init_sign = 1 for waves to propagate inward
!   (unless the wave is a backward wave).
    integer:: k0_sign

!   ray_param = 'arcl': default. Integrate with respect to the arclength along the ray.
!   ray_param = 'time': Integrate with respect to time along the ray.
    character(len = 4) :: ray_param = 'arcl'

!   Maximum allowable residual of dispersion function. Used in checksave()
    real(KIND=rkind) :: dispersion_resid_limit
    
    namelist /rf_list/ ray_dispersion_model, frf, wave_mode, k0_sign, ray_param, &
                     & dispersion_resid_limit
    
!********************************************************************

contains

!********************************************************************

  subroutine initialize_rf_m(read_input)
 
    use constants_m, only : pi, clight
    use diagnostics_m, only : message_unit, message
    
    implicit none
    logical, intent(in) :: read_input
 	integer :: input_unit, get_unit_number ! External, free unit finder   
 	
    if (read_input .eqv. .true.) then    
		! Read and write input namelist
  		  	input_unit = get_unit_number()
			open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
			read(input_unit, rf_list)
			close(unit=input_unit)
			write(message_unit, rf_list)
    end if

    if ( frf <= 0. ) then
       write (0,*) 'initialize_rf: frf =', frf
       stop 1
    else
       omgrf = 2.*pi * frf
    end if

!   Wave number in vacuum.
    k0 = omgrf/clight 
    call message ('initialize_rf: k0', k0, 1)
    
    write(*,*) 'wave_mode = ', wave_mode

    return
  end subroutine initialize_rf_m

    subroutine deallocate_rf_m
        continue ! nothing to deallocate
    end subroutine deallocate_rf_m

 end module rf_m
