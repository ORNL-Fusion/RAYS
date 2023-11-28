 subroutine trace_rays
! Does all ray tracing.  Outer loop iterates over different rays.  Inner loop iterates
! over ray parameter step, s, to calculate ray trajectory.  Output is saved after each step/

! External procedures: cpu_time (intrinsic), check_save (check_save.f90)

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, verbosity, write_formatted_ray_files, &
                             & messages_to_stdout,ray_list_unit, day_to_seconds, &
                             & date_to_julian
    use ode_m, only : ode_solver, ray_init_ode_solver, nv, ds, s_max, nstep_max, ode_stop
    use ray_init_m, only : nray, ray_pwr_wt
    use damping_m, only : damping_model, multi_spec_damping
    use species_m, only : nspec
    use ray_results_m, only : ray_stop_flag, ray_vec, residual, npoints, end_residuals,&
                            & max_residuals, end_ray_parameter, start_ray_vec, end_ray_vec,&
                            & initial_ray_power, ray_trace_time, total_trace_time
    use openmp_m, only : num_threads
    use omp_lib
    implicit none

! Time and date vector - local, not the one loaded in subroutine initialize()
    integer :: date_v(8), ierr
    real(KIND=rkind) :: trace_time, code_time

    integer :: iray, nstep
    real(KIND=rkind) :: s, sout, resid, t_start_ray, t_finish_ray
    real(KIND=rkind) :: t_start_tracing, t_finish_tracing

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

!   Get date and time i.e. before ray loop, convert to Julian -> t_start_tracing
    call date_and_time (values=date_v)
	call date_to_julian(date_v,t_start_tracing,ierr)
	if (ierr .ne. 0) then
		write(*,*) 'julian t_start_tracing, ierr = ', ierr
		stop
	end if

!!$    t_start_tracing = omp_get_wtime()

!$OMP parallel do schedule(static) DEFAULT(FIRSTPRIVATE) &
!$OMP& SHARED(ray_stop_flag, ray_vec, npoints, residual, ray_trace_time, initial_ray_power, &
!$OMP& end_residuals, max_residuals, end_ray_parameter, start_ray_vec, end_ray_vec)

!!$OMP DO
    ray_loop: do iray = 1, nray
!$  write(12,*) 'ray_tracing begin: ray# = ',iray,'  omp_get_thread_num = ', omp_get_thread_num()

         call message(1)
         call message ('ray_tracing: ray #', iray, 1)
         call message ('ray_tracing: omp_thread_num = ', omp_get_thread_num(), 1)

!         call cpu_time(t_start_ray)
!$      t_start_ray = omp_get_wtime()

         nstep=0
         s = 0.
         sout = 0.
         resid = 0.
         ray_stop%stop_ode = .false.
         ray_stop%ode_stop_flag = ''
         ray_stop_flag(iray) = ray_stop%ode_stop_flag

    !    Reset ode solver for beginning of ray.
         call ray_init_ode_solver(ray_stop)

    !    Initialization of ray vector v.
         call initialize_ode_vector(iray, nv, v)

    !    Save in ray_results_m
		 ray_vec(:,1,iray) = v(:)

         call message(1)
         call message ('trace_rays: initial (x,y,z)', v(1:3), 3, 1)
         call message ('trace_rays: initial (kx,ky,kz)', v(4:6), 3, 1)

    !    Do some checking and save initial values.
         call check_save(sout, nv, v, resid, ray_stop)
         if (ray_stop%stop_ode .eqv. .true.) then  ! Ray didn't start, initial conditions bad
            ray_stop_flag(iray) = ray_stop%ode_stop_flag
            npoints(iray) = 1 ! i.e. initial point

 			if (verbosity > 0) then
				write(message_unit, *) 'ray ', iray, ' did not start. ', ray_stop%ode_stop_flag
				if (messages_to_stdout) write(*, *) 'ray ', iray, ' did not start. ', &
				              & ray_stop%ode_stop_flag
			end if

            cycle ray_loop
 		  end if

    !***********************************

         trajectory: do

            s = sout
            sout = sout + ds

            call message(1)
            call message ('trace_rays: ray #', iray, 1)
            call message ('trace_rays: start step', nstep + 1, 1)
            call message ('trace_rays: s', s, 1)
            call message ('trace_rays: sout', sout, 1)

! check limits on s
            if(sout > s_max) then
                call message ('trace_rays: stop ray, sout > s_max, s',s, 1)

                if (verbosity > 0) then
					write (message_unit, '( "ray ",i3, "  s=", g12.4, "   nstep=", i4, /, &
					&  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
					& v(1:3), v(4:6)

					if (messages_to_stdout) then
						 write (*, '( "ray ",i3, "  s=", g12.4, "   nstep=", i4, /, &
						&  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') iray, s, nstep, &
						& v(1:3), v(4:6)
					end if
				end if

                ray_stop%stop_ode = .true.
                ray_stop%ode_stop_flag = 'sout > s_max'
                ray_stop_flag(iray) = ray_stop%ode_stop_flag
                exit trajectory
            end if

! check limits on nstep
            if(nstep + 1 > nstep_max) then
                ray_stop%stop_ode = .true.
                ray_stop%ode_stop_flag =  ' nstep > nstep_max'
                ray_stop_flag(iray) = ray_stop%ode_stop_flag
                nstep = nstep_max ! Last valid step was nstep_max

                call message ('trace_rays: terminate ray, nstep+1 > nstep_max, nstep', nstep, 1)

                if (verbosity > 0) then

					write (message_unit, '( "ray ",i3, " stopped, last valid step   s=", g12.4, &
					& "   nstep=", i4, /, "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') &
					&  iray, s, nstep, v(1:3), v(4:6)

					if (messages_to_stdout) then
						write (message_unit, '( "ray ",i3, " stopped, last valid step   s=", g12.4, &
						& "   nstep=", i4, /, "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') &
						&  iray, s, nstep, v(1:3), v(4:6)
					end if
				end if

                exit trajectory
            end if

! Integrate from s to sout
            call ode_solver(eqn_ray, nv, v, s, sout, ray_stop)

! check for stop condition inside ode solver. Step failed.
            if (ray_stop%stop_ode .eqv. .true.) then
                ray_stop_flag(iray) = ray_stop%ode_stop_flag

                if (verbosity > 0) then
					write(message_unit, *) 'ray ', iray, ' stopped in ODE solver on step #',&
					  & nstep+1, ',  s = ', s, ', ode_stop_flag =  ', trim(ray_stop%ode_stop_flag)
					write (message_unit, '( "last valid step   s=", g12.4, "   nstep=", i4, /, &
					&  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') s, nstep, &
					& v(1:3), v(4:6)

					if (messages_to_stdout) then
						write(*, *) 'ray ', iray, ' stopped in ODE solver on step #',&
						  & nstep+1, ',  s = ', s, ', ode_stop_flag =  ', trim(ray_stop%ode_stop_flag)
						write (message_unit, '( "last valid step   s=", g12.4, "   nstep=", i4, /, &
						&  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') s, nstep, &
						& v(1:3), v(4:6)
					end if
				end if

                exit trajectory
            end if

            call message ('trace_rays: (x,y,z)', v(1:3), 3, 1)
            call message ('trace_rays: (kx,ky,kz)', v(4:6), 3, 1)
            call message ('trace_rays: integrated path length s', v(7), 1)

            damping : if (damping_model /= 'no_damp') then
              call message ('trace_rays: Total P-abs', v(8), 1)

              if (multi_spec_damping .eqv. .true.) then
                call message('trace_rays: P_abs(species)', v(9:9+nspec), 1)
              end if
            end if damping

! Do some more checking checking and save output from step
            call check_save(s, nv, v, resid, ray_stop)

            if (ray_stop%stop_ode .eqv. .true.) then
                if (verbosity > 0) then
					write(message_unit, *) 'ray ', iray, ' stopped in check_save  ', &
										& ray_stop%ode_stop_flag, '  residual =', resid

					 write (message_unit, '( " s=", g12.4, "   nstep=", i4, /, &
					&  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') s, nstep, &
					& v(1:3), v(4:6)

					if (messages_to_stdout) then
						write(message_unit, *) 'ray ', iray, ' stopped in check_save  ', &
											& ray_stop%ode_stop_flag, '  residual =', resid

						 write (message_unit, '( " s=", g12.4, "   nstep=", i4, /, &
						&  "  (x,y,z)=", 3(f10.4),/,  "  (kx,ky,kz)=", 3(f10.4) )') s, nstep, &
						& v(1:3), v(4:6)
					end if
				end if

                exit trajectory
            end if

! increment nstep here if that was a valid step
            nstep=nstep+1

! Save in ray_results_m
! Tricky point: The number of points in the ray = point 0 + number of valid steps.
! So  actual npoints = nstep + 1
            ray_vec(:,nstep+1,iray) = v(:)
            residual(nstep+1,iray) = resid

        end do trajectory

	 !   call cpu_time(t_finish_ray)
!$      t_finish_ray = omp_get_wtime()

! Save in ray_results_m

        npoints(iray) = nstep + 1
	    initial_ray_power(iray) = ray_pwr_wt(iray)
	    ray_trace_time(iray) = t_finish_ray - t_start_ray
        end_residuals(iray) = residual(nstep,iray)
        max_residuals(iray) = maxval(abs(residual(1:nstep,iray)))
        end_ray_parameter(iray) = v(7)
        ray_stop_flag(iray) = ray_stop%ode_stop_flag
        start_ray_vec(:,iray) = ray_vec(:,1,iray)
        end_ray_vec(:, iray) = v(:)

!!$  write(12,*) 'ray_tracing end: ray# = ',iray,'  omp_get_thread_num = ', omp_get_thread_num()

    end do ray_loop

!$omp end parallel do
!$      t_finish_tracing = omp_get_wtime()

!   Get date and time i.e. after ray loop, convert to Julian -> t_finish_tracing
    call date_and_time (values=date_v)
	call date_to_julian(date_v,t_finish_tracing,ierr)
	if (ierr .ne. 0) then
		write(*,*) 'julian t_finish_tracing, ierr = ', ierr
		stop
	end if
	total_trace_time = (t_finish_tracing - t_start_tracing)*day_to_seconds
    call message('Wall time ray tracing', total_trace_time, 0)

!   Write ray file description
    if (write_formatted_ray_files) then
       write(ray_list_unit, *) nray
       write(ray_list_unit, *) nv
       write(ray_list_unit, *) npoints
       write(ray_list_unit, *) end_residuals
       write(ray_list_unit, *) ray_stop_flag
    end if

! Below are writes for binary file.  For now use formatted
!      write (95) nray
!      write (95) npoints
!      write (95) nv
!      write (95) ray_stop

!	call cpu_time(t_finish_tracing)

    return
 end subroutine trace_rays
