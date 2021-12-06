 subroutine trace_rays

    use constants_m, only : ray_list_unit
    use diagnostics_m, only : message_unit, message, stop_ode, verbosity, ray_stop_flag, &
                               t_start_tracing, t_finish_tracing
    use ode_m, only : ode_solver, ray_init_ode_solver, nv, v, ds, s_max, nstep_max
    use ray_init_m, only : nray
    use damping_m, only : damping_model, multi_spec_damping
    use species_m, only : nspec

    implicit none

    integer :: iray, nstep, npoints(nray)
    character(len = 20) :: ray_stop(nray)
    real :: s, sout

    interface
       subroutine eqn_ray(s, v, dvds)
          use ode_m, only : nv
          real ( kind = 4 ), intent(in) :: s
          real ( kind = 4 ), intent(in) :: v(:)
          real ( kind = 4 ), intent(out) :: dvds(:)
       end subroutine eqn_ray
    end interface

    call cpu_time(t_start_tracing)
       
    ray_loop: do iray = 1, nray  

         write (*,'(/,a,i4)') 'ray #', iray
         call message()
         call message ('ray_tracing: ray #', iray, 0)

         nstep=0
         s = 0. 
         sout = 0. 
         stop_ode = .false.
         ray_stop_flag = ''

    !    Reset ode solver for beginning of ray
         call ray_init_ode_solver
       
    !    Initialization of ray vector v.
         call initialize_ode_vector(iray, nv, v)
 
         call message()
         call message ('ray_tracing: initial (x,y,z)', v(1:3), 3, 1)
         call message ('ray_tracing: initial (kx,ky,kz)', v(4:6), 3, 1)

    !    Do some checking and save initial values.
         call check_save(sout, nv, v)
         if (ray_stop_flag .ne. '') then  ! Ray didn't get started, initial conditions bad
            ray_stop(iray) = ray_stop_flag
            npoints(iray) = 1 ! i.e. initial point
            write(message_unit, *) 'ray ', iray, ' did not start. ', ray_stop_flag
            cycle ray_loop
         end if
       
    !***********************************    

         trajectory: do
       
            s = sout
            sout = sout + ds
            nstep=nstep+1
            
! check limits on s
            if(sout > s_max) then
                call message ('trace_rays: terminate ray, sout > s_max, s',s,0)
                write (*, *) 'trace_rays: terminate ray, sout > s_max, s = ',s
                ray_stop_flag = 'sout > s_max'
                ray_stop(iray) = ray_stop_flag
                exit trajectory 
            end if
            
! check limits on nstep
            if(nstep > nstep_max) then
                call message ('trace_rays: terminate ray, nstep > nstep_max, nstep', nstep, 0)
                write (*, *) 'trace_rays: terminate ray, nstep > nstep_max, nstep = ',nstep
                ray_stop_flag =  ' nstep > nstep_max'
                ray_stop(iray) = ray_stop_flag
               exit trajectory
            end if
          
            call message ('trace_rays: nstep', nstep, 2)
            call message ('trace_rays: sout', sout, 2)

! Integrate from s to sout
            call ode_solver(eqn_ray, nv, v, s, sout)

! check for stop condition inside ode solver. Step failed.
            if (stop_ode .eqv. .true.) then 
                ray_stop(iray) = ray_stop_flag
                write(message_unit, *) 'ray ', iray, ' stoped in ODE solver. ', ray_stop_flag
                
                write (*, '( "ray ",i3, " stopped  s=", f12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
                & v(1:3), v(4:6)
                
                write (message_unit, '( "ray ",i3, " stopped  s=", f12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
                & v(1:3), v(4:6)
                exit trajectory
            end if    

            call message()
            call message ('trace_rays: nstep', nstep, 1)
            call message ('trace_rays: s', s, 1)
            call message ('trace_rays: (x,y,z)', v(1:3), 3, 1)
            call message ('trace_rays: (kx,ky,kz)', v(4:6), 3, 1)
            call message ('trace_rays: integrated path length s', v(7), 1)

            damping : if (damping_model /= 'no_damp') then    
              call message ('trace_rays: Total P-abs', v(8), 1)

              if ((multi_spec_damping .eqv. .true.) .and. (verbosity >= 1)) then
                write (message_unit,*) 'trace_rays: P_abs(species)', v(9:9+nspec)
              end if  
            end if damping
 
! Do some checking and save output from step
            call check_save(s, nv, v)

            if (ray_stop_flag .ne. '') then
                ray_stop(iray) = ray_stop_flag
                write(message_unit, *) 'ray ', iray, ' stoped. ', ray_stop_flag
                
                write (*, '( "ray ",i3, " stopped  s=", f12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
                & v(1:3), v(4:6)
                
                write (message_unit, '( "ray ",i3, " stopped  s=", f12.4, "   nstep=", i4, /, &
                &  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
                & v(1:3), v(4:6)
                exit trajectory
            end if

        end do trajectory 

! Tricky point: The number of points in the ray = point 0 + number of valid steps.   
! nstep is incremented at the top of the loop, but "exit trajectory" above means last 
! step failed, and so is 1 too big, so actual npoints = nstep not nstep + 1
        npoints(iray)=nstep

    end do ray_loop

!   Signature.
!    write(94) -12345

!   Write ray file description
     write(ray_list_unit, *) nray
     write(ray_list_unit, *) npoints
     write(ray_list_unit, *) nv
     write(ray_list_unit, *) ray_stop

! Below are writes for binary file.  For now use formatted
!      write (95) nray
!      write (95) npoints
!      write (95) nv
!      write (95) ray_stop

    call cpu_time(t_finish_tracing)

    return
 end subroutine trace_rays
