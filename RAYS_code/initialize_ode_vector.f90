 subroutine initialize_ode_vector(iray, nv, v)

! Initializes ode vector v before integrating each ray. 
! Called from subroutine ray_tracing

    use constants_m, only : rkind
    use diagnostics_m, only : integrate_eq_gradients, message
    use species_m, only : nspec, n0s
    use equilibrium_m, only : equilib_model, equilibrium, eq_point
    use rf_m, only : k0
    use ray_init_m, only : rvec0, rindex_vec0
    use damping_m, only : damping_model, multi_spec_damping

    implicit none

    integer, intent(in) :: iray, nv
    real(KIND=rkind), intent(inout) :: v(nv)

    type(eq_point(nspec=nspec)) :: eq

    integer :: nv0

!   Initialize v = (1:6) 

    v(1:3) = rvec0(1:3, iray)
    v(4:6) = k0*rindex_vec0(1:3, iray)
    v(7) = 0.

! Start counting length of v(:) i.e. nv0
    nv0 = 7
    
!   Initialize power deposition for each species.
    damping : if (damping_model /= 'no_damp') then    
        v(8) = 0.  ! Total absorbed power
        nv0 = 8
        if (multi_spec_damping .eqv. .true.) then
          v(9:9+nspec) = 0. ! Absorption by species
          nv0 = nv0 + 1 + nspec
        end if
    end if damping
    
! Initialize diagnostic checks on consistency of equilibrium quantities and
!  their gradients

    gradients : if ( integrate_eq_gradients .eqv. .true. ) then
       call equilibrium(v(1:3),eq)
!      For checking if grad(B) is consistent with B.
       v(nv0+1:nv0+3) = eq%bvec(1:3)
 
!      For checking if grad(ne) and grad(Te) are consistent with ne and Te.
       v(nv0+4) = eq%ns(0)/n0s(0)
       v(nv0+5) = eq%ts(0)
     
    end if gradients
    
    return 
 end subroutine initialize_ode_vector
