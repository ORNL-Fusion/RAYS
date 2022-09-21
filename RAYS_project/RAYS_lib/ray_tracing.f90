 subroutine trace_rays
! Does all ray tracing.  Outer loop iterates over different rays.  Inner loop iterates
! over ray parameter step, s, to calculate ray trajectory.  Output is saved after each step/

! External procedures: cpu_time (intrinsic), check_save (check_save.f90)

    use constants_m, only : rkind, ray_list_unit
    use diagnostics_m, only : message_unit, message, verbosity, &
                               t_start_tracing, t_finish_tracing
    use ode_m, only : ode_solver, ray_init_ode_solver, nv, ds, s_max, nstep_max, ode_stop
    use ray_init_m, only : nray
    use damping_m, only : damping_model, multi_spec_damping
    use species_m, only : nspec
    use ray_results_m, only : ray_stop_flag, ray_vec, residual, npoints, end_residuals,&
                            & max_residuals, end_ray_parameter, end_ray_vec, ray_trace_time

    implicit none

    integer :: iray, nstep
    real(KIND=rkind) :: s, sout, resid
 
!   v: Vector to be integrated by ode solver
    real(KIND=rkind) :: v(nv)
    type(ode_stop)  :: ray_stop

    interface
       subroutine eqn_ray(s, v, dvds)
          use constants_m, only : rkind
          use ode_m, only : nv
          real ( kind = rkind ), intent(in) :: s
          real ( kind = rkind ), intent(in) :: v(:)
          real ( kind = rkind ), intent(out) :: dvds(:)
       end subroutine eqn_ray
    end interface

    call cpu_time(t_start_tracing)
       
    ray_loop: do iray = 1, nray  

         write (*,'(/,a,i4)') 'ray #', iray
         call message()
         call message ('trace_rays: ray #', iray, 0)

         nstep=0
         s = 0. 
         sout = 0.
         resid = 0. 
         ray_stop%stop_ode = .false.
         ray_stop%ode_stop_flag = ''
         ray_stop_flag(iray) = ray_stop%ode_stop_flag

    !    Reset ode solver for beginning of ray
         call ray_init_ode_solver
       
    !    Initialization of ray vector v.
         call initialize_ode_vector(iray, nv, v)
         
    !    Save in ray_results_m 
		 ray_vec(:,1,iray) = v(:)
 
         call message()
         call message ('trace_rays: initial (x,y,z)', v(1:3), 3, 1)
         call message ('trace_rays: initial (kx,ky,kz)', v(4:6), 3, 1)

    !    Do some checking and save initial values.
         call check_save(sout, nv, v, resid, ray_stop)
         if (ray_stop%stop_ode .eqv. .true.) then  ! Ray didn't start, initial conditions bad
            ray_stop_flag(iray) = ray_stop%ode_stop_flag
            npoints(iray) = 1 ! i.e. initial point
            write(message_unit, *) 'ray ', iray, ' did not start. ', ray_stop%ode_stop_flag
            cycle ray_loop
         end if
         
      
    !***********************************    

         trajectory: do
       
            s = sout
            sout = sout + ds
            nstep=nstep+1

            call message()
            call message ('trace_rays: start step', nstep, 1)
            call message ('trace_rays: s', s, 1)
            call message ('trace_rays: sout', sout, 1)
          
! check limits on s
            if(sout > s_max) then
                call message ('trace_rays: stop ray, sout > s_max, s',s,0)
                write (*, *) 'trace_rays: stop ray, sout > s_max, s = ',s
               
                write (*, '( "ray ",i3, "  s=", g12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
                & v(1:3), v(4:6)
                
                write (message_unit, '( "ray ",i3, "  s=", g12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
                & v(1:3), v(4:6)
                
                ray_stop%stop_ode = .true.
                ray_stop%ode_stop_flag = 'sout > s_max'
                ray_stop_flag(iray) = ray_stop%ode_stop_flag
                exit trajectory 
            end if
            
! check limits on nstep
            if(nstep > nstep_max) then
                call message ('trace_rays: terminate ray, nstep > nstep_max, nstep', nstep, 0)
                write (*, *) 'trace_rays: terminate ray, nstep > nstep_max, nstep = ',nstep
                ray_stop%ode_stop_flag =  ' nstep > nstep_max'
                
                nstep = nstep_max ! Last valid step was nstep_max
               
                write (*, '( "ray ",i3, " stopped, last valid step  s=", g12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
                & v(1:3), v(4:6)
                
                write (message_unit, '( "ray ",i3, " stopped, last valid step   s=", g12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
                & v(1:3), v(4:6)
                
                ray_stop%stop_ode = .true.
                ray_stop_flag(iray) = ray_stop%ode_stop_flag
               exit trajectory
            end if
          
! Integrate from s to sout
            call ode_solver(eqn_ray, nv, v, s, sout, ray_stop)

! check for stop condition inside ode solver. Step failed.
            if (ray_stop%stop_ode .eqv. .true.) then 
                ray_stop_flag(iray) = ray_stop%ode_stop_flag
                write(message_unit, *) 'ray ', iray, ' stopped in ODE solver. ', &
                    & ray_stop%ode_stop_flag
                write(*, *) 'ray ', iray, ' stopped in ODE solver. ', &
                    & ray_stop%ode_stop_flag
                
                write (*, '( "ray ",i3, " stopped, last valid step   s=", g12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
                & v(1:3), v(4:6)
                
                write (message_unit, '( "ray ",i3, " stopped, last valid step   s=", g12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
                & v(1:3), v(4:6)
                exit trajectory
            end if    

            call message ('trace_rays: (x,y,z)', v(1:3), 3, 1)
            call message ('trace_rays: (kx,ky,kz)', v(4:6), 3, 1)
            call message ('trace_rays: integrated path length s', v(7), 1)

            damping : if (damping_model /= 'no_damp') then    
              call message ('trace_rays: Total P-abs', v(8), 1)

              if ((multi_spec_damping .eqv. .true.) .and. (verbosity >= 1)) then
                write (message_unit,*) 'trace_rays: P_abs(species)', v(9:9+nspec)
              end if  
            end if damping
 
! Do some more checking checking and save output from step
            call check_save(s, nv, v, resid, ray_stop)

            if (ray_stop%stop_ode .eqv. .true.) then
                write(message_unit, *) 'ray ', iray, ' stopped, last valid step  ', ray_stop%ode_stop_flag
                write(*, *) 'ray ', iray, ' stopped, last valid step  ', ray_stop%ode_stop_flag,&
                  & '  residual  ', resid
                
                write (*, '( " s=", g12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') s, nstep, &
                & v(1:3), v(4:6)
                
                write (message_unit, '( " s=", g12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') s, nstep, &
                & v(1:3), v(4:6)
                exit trajectory
            end if

! Save in ray_results_m
            ray_vec(:,nstep+1,iray) = v(:)
            residual(nstep,iray) = resid
            
        end do trajectory
        
! Save in ray_results_m

! Tricky point: The number of points in the ray = point 0 + number of valid steps.   
! nstep is incremented at the top of the loop, but "exit trajectory" above means last 
! step failed, and so is 1 too big, so actual npoints = nstep, not nstep + 1
        npoints(iray) = nstep
        end_residuals(iray) = residual(nstep,iray)
        max_residuals(iray) = maxval(abs(residual(1:nstep,iray)))
        end_ray_parameter(iray) = s
        ray_stop_flag(iray) = ray_stop%ode_stop_flag
        end_ray_vec(:, iray) = v(:)

    end do ray_loop

!   Signature.
!    write(94) -12345

!   Write ray file description
     write(ray_list_unit, *) nray
     write(ray_list_unit, *) nv
     write(ray_list_unit, *) npoints
     write(ray_list_unit, *) end_residuals
     write(ray_list_unit, *) ray_stop_flag

! Below are writes for binary file.  For now use formatted
!      write (95) nray
!      write (95) npoints
!      write (95) nv
!      write (95) ray_stop

	call cpu_time(t_finish_tracing)
	ray_trace_time = t_finish_tracing - t_start_tracing

    return
 end subroutine trace_rays
