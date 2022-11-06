 subroutine eqn_ray(s, v, dvds, ray_stop)

!   Calculates the RHS of the ray's equation.
!   By default there are seven equations; three for r(s), three for k(s) and one for s
!   v(1:3) = r, v(4:6) = k,  v(7) = s.
!
!   The v(7) equation integrates v_group along the ray and serves as a diagnostic.  
!   If ray_param = 'arcl' it serves as a check on the consistency of v_group with s.
!   If ray_param = 'time' it gives arc length s(t) along the ray trajectory.
!
!   Optionally v(8) integrates total absorption, ki(s)
!   Optionally v(9:9+nspec)integrates absorption by species, kis(s)
!   Optionally v(5 more elements)integrates gradients of the equilibrium quantities to check 
!   consistency
!
!   Selects between different ray dispersion models to calculate derivatives of D(r,x,omega)
!
!   N.B. In integrating v(7), damping, or diagnostic gradients the integral is with respect
!        to arc length, so introduce parameter, dsd_ray_param = 1 for ray_param = 'arcl',
!        dsd_ray_param = |v group| == vg0 for ray_param = 'time'
 
!   External procedures: deriv_cold (deriv_cold)

    use constants_m, only : rkind
    use diagnostics_m, only : integrate_eq_gradients, verbosity, message_unit, message
    use species_m, only : nspec
    use equilibrium_m, only : equilibrium, eq_point
    use ode_m, only : nv, ray_deriv_name, ode_stop
    use rf_m, only : k0, ray_dispersion_model, ray_param
    use damping_m, only : damping_model, damping, multi_spec_damping
    
     implicit none
     
    real(KIND=rkind), intent(in) :: s
    real(KIND=rkind), intent(in) :: v(nv)
    real(KIND=rkind), intent(out) :: dvds(nv) 
    type(ode_stop), intent(out)  :: ray_stop

    real(KIND=rkind) :: rvec(3)
    real(KIND=rkind) :: kvec(3), k1, k3
    real(KIND=rkind) :: nvec(3) 

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

    real(KIND=rkind) :: dddx(3), dddk(3), dddw, vg(3), vg0, vg_unit(3)
    real(KIND=rkind) :: ksi(0:nspec), ki
    real(KIND=rkind) :: dsd_ray_param
    integer :: nv0

    integer :: is

!***********************************
 
    interface deriv_cold
       subroutine deriv_cold(eq, nvec, dddx, dddk, dddw)
          use constants_m, only : rkind
          use equilibrium_m, only : equilibrium, eq_point
          type(eq_point), intent(in) :: eq
          real(KIND=rkind), intent(in) :: nvec(3)
          real(KIND=rkind), intent(out) :: dddx(3), dddk(3), dddw
       end subroutine deriv_cold
    end interface deriv_cold
 
    interface deriv_num
       subroutine deriv_num(eq, v, dddx, dddk, dddw)
            use constants_m, only : rkind, clight
            use equilibrium_m, only : equilibrium, eq_point, write_eq_point
            use rf_m, only : omgrf, k0, ray_dispersion_model
            use ode_m, only : nv
            use diagnostics_m, only : message_unit, verbosity
            type(eq_point), intent(in) :: eq
            real(KIND=rkind), intent(in) :: v(nv)
            real(KIND=rkind), intent(out) :: dddx(3), dddk(3), dddw
       end subroutine deriv_num
    end interface deriv_num

!   kvec (nvec) = k (k/k0) in xyz coordinates (vector)
    rvec = v(1:3); kvec = v(4:6)
    nvec = kvec/k0

!   Calculate the plasma equilibrium.
    call equilibrium(rvec, eq)

! Check if equib_err has been set in equilibrium
    if (trim(eq%equib_err) /= '') then
        ray_stop%stop_ode = .true.
        write (*,*) 'eqn_ray: s = ', s, '  equib_err = ', eq%equib_err
        write (message_unit,'(a,g12.4,3a,3g12.4)') 'eqn_ray: s = ', s, '  equib_err = ',&
              & trim(eq%equib_err), '  r_end = ', rvec
        ray_stop%ode_stop_flag = eq%equib_err
        return
    end if

!   Calculate dD/dk, dD/dx, and dD/d(omega) 

    dispersion_model: select case (trim(ray_deriv_name))

        case ('cold' )
    !      Derivatives of D for a cold plasma.
         call deriv_cold(eq, nvec, dddx, dddk, dddw)

        case ('numerical' )
    !      Numerical differentiation.
    !      N.B. Must be called with v(), not just nvec.  Evaluates eq at other positions so v(1:3)
    !           is needed

         call deriv_num(eq, v, dddx, dddk, dddw)

        case default
           write(*,*) 'EQN_RAY: invalid value, ray_deriv_name = ', ray_deriv_name
           stop 1
    end select dispersion_model
       
   if (verbosity > 3) then  ! Debugging diagnostics
       write (*,*) 'eqn_ray: rvec = ', rvec, '  kvec = ', kvec
       write (message_unit,*) 'eqn_ray: dddx = ', dddx, 'dddk = ', dddk, 'dddw = ', dddw
       write (*,*) 'eqn_ray: dddx = ', dddx, 'dddk = ', dddk, 'dddw = ', dddw
   end if


!   Group velocity. 

    if ( dddw /= 0. ) then
       vg = -dddk / dddw
       vg0 = sqrt(sum(vg**2))
!      Unit vector along the group velocity.
       vg_unit = vg / vg0
!      write(6,'(a,1p3e12.4)') 'EQN_RAY: vg/|vg| =', vg_unit
    else
       write(*,*) 'EQN_RAY: infinite group velocity, dddw = ', dddw
       ray_stop%stop_ode = .true.
       ray_stop%ode_stop_flag = 'infinite Vg'
       return
    end if

!   Geometrical optics equations for a ray.

    parameter: select case (ray_param)

       case ('arcl')
!         Integrate with respect to the arclength along the ray.
          if ( any(dddk/=0.) ) then
!            dr/ds:
             dvds(1:3) = -sign(real(1.,KIND=rkind),dddw) * dddk / sqrt(sum(dddk**2)) 

!            dk/ds:
             dvds(4:6) =  sign(real(1.,KIND=rkind),dddw) * dddx / sqrt(sum(dddk**2))

!            ds/d(ray parameter)
             dsd_ray_param = 1.
             
!            Unit vector along the group velocity.
             vg_unit = vg / vg0
!            write(6,'(a,1p3e12.4)') 'EQN_RAY: vg/|vg| =', vg_unit
          else
             write(0,*) 'EQN_RAY: ray stalled, dddk = ', dddk
             write(message_unit,*) 'EQN_RAY: ray stalled, dddk = ', dddk
             ray_stop%stop_ode = .true.
             ray_stop%ode_stop_flag = 'ray stalled'
             return
          end if

       case ('time')
!         Integrate with respect to time along the ray.
!         dr/dt:
          dvds(1:3) = -dddk / dddw

!         dk/dt:
          dvds(4:6) =  dddx / dddw

!         ds/d(ray parameter)
          dsd_ray_param = vg0
          
!         Calculate ksi and ki.
          call damping(eq, v, vg, ksi, ki)


       case default
          write(0,*) 'EQN_RAY: invalid ray parameter = ', ray_param
          stop 1

    end select parameter

! Integrated arc length
    dvds(7) = dsd_ray_param

! Start counting length of v(:) i.e. nv0
    nv0 = 7

! Damping Calculate ksi and ki.
    damp : if (damping_model /= 'no_damp') then    
    
        call damping(eq, v, vg, ksi, ki)

!       Differential equation for total power absorption.
        nv0 = nv0 + 1
        dvds(nv0) = dsd_ray_param * 2.*ki*(1.-v(nv0))        


!       Differential equation for power absorption for each species.
        if (multi_spec_damping .eqv. .true.) then
          do is = 0, nspec
             dvds(nv0+1+is) = dsd_ray_param * 2.*ksi(is)*(1.-v(nv0))
          end do
          nv0 = nv0 +1 +nspec
        end if
    end if damp
    
! Optionally integrate additional diagnostic equations along ray

    if (integrate_eq_gradients .eqv..true.) then

!       Check if grad(B) is consistent with B.
        dvds(nv0+1)  = sum(dsd_ray_param*vg_unit*eq%gradbtensor(:,1))
        dvds(nv0+2) = sum(dsd_ray_param*vg_unit*eq%gradbtensor(:,2))
        dvds(nv0+3) = sum(dsd_ray_param*vg_unit*eq%gradbtensor(:,3))

!       Check if grad(ne) and grad(Te) are consistent with ne and Te.
!       ne normalized to peak electron density 
        dvds(nv0+4) = sum(dsd_ray_param*vg_unit*eq%gradns(:,0))
        dvds(nv0+5) = sum(dsd_ray_param*vg_unit*eq%gradts(:,0))
     
    end if

    if (verbosity > 3) then  ! Debugging diagnostics
        write (*, *) 'eqn_ray: dvds = ', dvds
    end if
                
    return
 end subroutine eqn_ray
