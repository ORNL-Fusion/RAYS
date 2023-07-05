 subroutine check_save(s, nv, v, resid, ray_stop)
!   Does checking for ray stop criteria and saves output after each step.
!   External routines: deriv_cold (deriv_cold.f90)

    use constants_m, only : rkind
    use diagnostics_m, only : integrate_eq_gradients, message, text_message, &
                            & write_formatted_ray_files, output_unit
    use species_m, only : nspec, n0s
    use equilibrium_m, only : equilibrium, eq_point
    use ode_m, only : ds, ode_stop
    use rf_m, only : ray_dispersion_model, ray_param, k0, dispersion_resid_limit
    use damping_m, only : damping_model, damping, multi_spec_damping, total_damping_limit

    implicit none

    real(KIND=rkind), intent(in) :: s
    integer, intent(in) :: nv
    real(KIND=rkind), intent(in) :: v(nv)
    real(KIND=rkind), intent(out) :: resid
    type(ode_stop), intent(out)  :: ray_stop

    type(eq_point(nspec=nspec)) :: eq
   
    real(KIND=rkind) :: kvec(3), k1, k3
    real(KIND=rkind) :: nvec(3), n1, n3  

    integer :: nv0
    real(KIND=rkind) :: dddx(3), dddk(3), dddw, vg(3), vg0
    real(KIND=rkind) :: ksi(0:nspec), ki
    real(KIND=rkind) :: total_absorption, total_species_absorption
    real(KIND=rkind) :: bmod, diff_vec(3)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 
!   Calculate the plasma equilibrium.
    call equilibrium(v(1:3), eq)
    if (eq%equib_err .ne. '') then
        ray_stop%ode_stop_flag = eq%equib_err
    end if

    call message ('check_save: ne', eq%ns(0), 1)
    call message ('check_save: B(tesla)', eq%bmag, 1)
    call message ('check_save: alpha s', eq%alpha(0:nspec), nspec+1, 1)
    call message ('check_save: gamma s', eq%gamma(0:nspec), nspec+1, 1)


!   Parallel and perpendicular components of nvec.
    kvec = v(4:6)
    k3 = sum(kvec*eq%bunit)
    k1 = sqrt( sum((kvec-k3*eq%bunit)**2) )

    nvec = kvec/k0; n1 = k1/k0; n3 = k3/k0

    call message ('check_save: ray k-parallel', k3, 1)
    call message ('check_save: ray k-perp', k1, 1)
    call message ('check_save: ray n-parallel', n3, 1)
    call message ('check_save: ray n-perp', n1, 1)

!   Calculate the residual.
    
    resid = residual(eq, k1,k3)
    call message ('check_save: residual', resid, 1)
    if (resid > dispersion_resid_limit) then
        ray_stop%stop_ode = .true.
        ray_stop%ode_stop_flag = 'dispersion_residual'
    end if
        



!   Calculate the group velocity. It's a piece of useful information,
!   especially, when iray_param = 'time', i.e., integrate the ray equations with 
!   respect to time.  This part is extracted from EQN_RAY.
!   First, calculate dD/dk, dD/dx, and dD/d(omega)


       
    if ( ray_dispersion_model == 'cold' ) then

!      Derivatives of D for a cold plasma.
       call deriv_cold(eq, nvec, dddx, dddk, dddw)      
    else
       call text_message('CHECK_SAVE: ray_dispersion_model = ', ray_dispersion_model)
       write(*,*) 'CHECK_SAVE: ray_dispersion_model = ', ray_dispersion_model
       stop 'check_save: unimplemented ray_dispersion_model'
    end if

!   Group velocity.
    if ( abs(dddw) > tiny(dddw) ) then
       vg = -dddk / dddw
       vg0 = sqrt(sum(vg**2))

       if ( trim(ray_param) == 'time' ) then
!         Note that when ray_param = 'time', "ds" in the code is "dt". Time steps are a
!         fraction of a nonosecond.
          if ( vg0*ds > 1. ) then
             call message('CHECK_SAVE: time step too big, Vg*dt (m) =', vg0 * ds, 1)
          end if
       end if
       
       call message('check_save: Vg',vg, 3, 1)
       
    else
       call message('CHECK_SAVE: dddw = ', dddw)
       call text_message('check_save: infinite group velocity, dddw = 0', 1)
       ray_stop%stop_ode = .true.
       ray_stop%ode_stop_flag = 'infinite_Vg'       
    end if

! Start counting length of v(:) i.e. nv0
    nv0 = 7

    damp : if (damping_model /= 'no_damp') then    
        call damping(eq, v, vg, ksi, ki)
        nv0 = nv0 + 1

        total_absorption = v(nv0)
        call message ('check_save: ni', ki/k0, 1)
        call message ('check_save: Total absorption', total_absorption, 1)
        if (total_absorption > total_damping_limit) then
 	       ray_stop%stop_ode = .true.
           ray_stop%ode_stop_flag = 'total_absorption'
           call text_message('ray_stop%ode_stop_flag = ', ray_stop%ode_stop_flag, 1)
        end if

!       Check if the sum of all by species absorption is equal to the total.        
        if (multi_spec_damping .eqv. .true.) then
           total_species_absorption =  sum(v(nv0+1:nv0+1+nspec))
           call message ('check_save: Total species absorption', total_species_absorption, 1)
           nv0 = nv0+1+nspec
        end if
    end if damp

!   Check the accuracy of integrated gradients of equilibrium compared to quantities
!   provided by the equilirium routines.

    integrate_gradients : if (integrate_eq_gradients .eqv. .true.) then
!      Check if grad(B) is consistent with B.
        bmod = sqrt(sum(eq%bvec(1:3))**2)
        diff_vec= (v(nv0+1:nv0+3) - eq%bvec(1:3))/bmod
        call message ('check_save: B relative error', diff_vec, 3, 1)

!       Check if grad(Te) and grad(ne) are consistent with Te and ne.
        call message('check_save: ne relative error', (v(nv0+4)-eq%ns(0))/n0s(0), 1)
        call message('check_save: Te difference',  v(nv0+5)-eq%ts(0),1)
    end if integrate_gradients

    
!****** write ray vector for this step to output file *********************************    
    
    if (write_formatted_ray_files) then
		write(output_unit, *) s, v
    end if

    return


!*************************************************************************     

contains
 
    real(KIND=rkind) function residual(eq, k1, k3)
!      calculates the residual for given k1 and k3.
!      get dielectric tensor from module suscep_m

       use species_m, only : nspec
       use suscep_m, only :  dielectric_cold

       implicit none 
       
       type(eq_point(nspec=nspec)), intent(in) :: eq
 
       real(KIND=rkind) :: k1, k3

       complex(KIND=rkind) :: eps(3,3), eps_h(3,3), epsn(3,3), ctmp
       complex(KIND=rkind) :: eps_norm(3,3)
       real(KIND=rkind) :: n(3)

       integer :: i, j

!   Need dielectric tensor.  If ray_dispersion_model = "warm" or "numeric" eps will
!   have already been calculated.  If ray_model = "cold" must calculate eps.

    if (ray_dispersion_model == "cold") then
        call dielectric_cold(eq, eps)
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

       residual = abs(ctmp) / &
          & ( eps_norm(3,3)*(eps_norm(1,1)*eps_norm(2,2)) &
          & + eps_norm(3,3)*(eps_norm(2,1)*eps_norm(1,2)) &
          & + eps_norm(3,2)*(eps_norm(1,1)*eps_norm(2,3)) &
          & + eps_norm(3,2)*(eps_norm(2,1)*eps_norm(1,3)) &
          & + eps_norm(3,1)*(eps_norm(1,2)*eps_norm(2,3)) &
          & + eps_norm(3,1)*(eps_norm(2,2)*eps_norm(1,3)) ) 
      

       return
    end function residual
 
 end subroutine check_save
