 module damping_m
!   contains parameters and routines to calculate damping

    use constants_m, only : rkind
    
    implicit none
     
!   Switch to select which model to use for damping calculation
!   damping_model = 'no_damp" do not calculate damping
!   damping_model = 'poynting' use multicomponent Poyntings theorem
!   damping_model = 'fund_ECH' use simple weak damping approx. for fundamental ECH
    character(len=60) :: damping_model
    
!   Multi species damping.  Only meaningful if damping_model /= 'no_damp"
!   If .true. integrate damping by individual species as well as total damping
!   If .false. integrate damping only total damping
    logical :: multi_spec_damping = .false.

!   Ray is considered totally damped if damping > total_damping_limit
    real(KIND=rkind) :: total_damping_limit = 0.99
    
! The warm plasma quantities are needed in routine poynting - generated in routine deriv()
    
    complex(KIND=rkind) :: depsdw_h3x3(3,3)
    
    namelist /damping_list/ damping_model, multi_spec_damping, total_damping_limit
    

!********************************************************************

 contains

!********************************************************************

  subroutine initialize_damping_m(read_input)
 
    use constants_m, only : input_unit
    use diagnostics_m, only : message_unit
    
    implicit none
    logical, intent(in) :: read_input

    if (read_input .eqv. .true.) then
    ! Read and write input namelist
        open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
        read(input_unit, damping_list)
        close(unit=input_unit)
        write(message_unit, damping_list)
    end if

    return
  end subroutine initialize_damping_m


!********************************************************************


 subroutine damping(v, k1, k3, ksi, ki, vg)
 
!   Wrapper subroutine to call one of several damping models based on the
!   value of damping_model


!   k1: component of k perpendicular to B (magnitude).
!   k3: component of k parallel to B (magnitude).

    use diagnostics_m, only : message, text_message
    use species_m, only : nspec

    implicit none

    real(KIND=rkind), intent(in) :: v(:), k1, k3
    real(KIND=rkind), intent(out) :: ksi(0:nspec), ki

    real(KIND=rkind), intent(in) :: vg(3)    
    
    model: select case (trim(damping_model))

        case ('no_damp')    ! do not calculate damping
    
            ksi(0:nspec) = 0.
            ki=sum(ksi(0:nspec))
    
!         case ('poynting')   ! Poyntings theorem
!     
!             call poynting(v, k1, k3, ksi, ki, vg)
!         
!         case ('fund_ECH')   ! simple weak damping approximation for fundamental ECH
!     
!             call damp_fund_ECH (ksi, ki, vg)
    
        case default
    
            write (0, *) 'damping: Unimplemented damping model', trim(damping_model)
            call text_message('damping: Unimplemented damping model', trim(damping_model), 0)
            stop 1
    
    end select model
    
    
    return
    
 end subroutine damping
    
!********************************************************************

    subroutine finalize_damping_m
        return ! nothing to clean up
    end subroutine finalize_damping_m

 end module damping_m

