 submodule (ode_m) SG_ode_m
!   contains parameters specific to SG_ode ODE solver.

    use constants_m, only : rkind

!   Access type definition ode_stop should be available by host association from ode_m.
!   But it confuses make not to use ode_m since it thinks SG_ode_m doesn't depend on it.
!   Therefore:
    use ode_m, only : ode_stop
    
    implicit none

!   Target relative and absolute errors for the ODE solver.
    real(KIND=rkind) :: rel_err0, abs_err0

!   Evolving relative and absolute errors from the ODE solver.
    real(KIND=rkind) :: rel_err, abs_err

!   Total ODE error limit abs(rel_err)+abs(abs_err) above which to bail.
    real(KIND=rkind) :: SG_error_limit = 0.1  ! Default

!   Return status flag
    integer :: iflag


 namelist /SG_ode_list/ rel_err0, abs_err0, SG_error_limit

!********************************************************************

contains

!********************************************************************

  module subroutine initialize_SG_ode(read_input)

     use diagnostics_m, only : message_unit, message, text_message, messages_to_stdout, verbosity
    
    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder   

    call message(1)
    call text_message('Initializing SG_ode_m ', 1)

	if (read_input .eqv. .true.) then        
	! Read and write input namelist
  		input_unit = get_unit_number()
		open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
		read(input_unit, SG_ode_list)
		close(unit=input_unit)
	end if

! Write input namelist
    if (verbosity > 0) then
		write(message_unit, SG_ode_list)
		if (messages_to_stdout) write(*, SG_ode_list)
    end if

!   Check error criterion for the ODE solver.  Put a limit on how small target error can be.
    if ( rel_err0 < 1.e-10 .or. abs_err0 < 1.e-10 ) then
       write(0,*) 'initialize_SG_ode: rel_err0, abs_err0 too small =', rel_err0, abs_err0
       stop 1
    end if

    return
  end subroutine initialize_SG_ode

!********************************************************************

  module subroutine ray_init_SG_ode
  ! Resets SG_ode solver to appropriate state for start of ray.
  ! This is nothing but resetting rel_err and abs_err

    rel_err = rel_err0
    abs_err = abs_err0

    return
  end subroutine ray_init_SG_ode

!********************************************************************

  module subroutine SG_ode(eqn_ray, nv, v, s, sout, ray_stop)
  ! Implements generic interface for SG_ode solver (i.e. Shampine & Gordon ode.f90)
  ! Takes one step from s to sout.  SG_ode may take intermediate steps to get to sout
  ! as indicated by iflag=3.  So call to ode is in a do loop.  See comments in ode()
  ! for details.

    use diagnostics_m, only : message_unit, message, text_message, verbosity
    use constants_m, only : rkind

  ! Arguments of ODE
    external eqn_ray
    integer, intent(in) :: nv
    real(KIND=rkind), intent(inout) :: v(nv)
    real(KIND=rkind), intent(inout) :: s, sout
    type(ode_stop), intent(out)  :: ray_stop

    real(KIND=rkind) :: work(100+21*nv)
    integer :: iwork(5)
    
    real(KIND=rkind) :: total_error

    odeloop: do

       iflag = 1
 
       call ode(eqn_ray, nv, v, s, sout, rel_err, abs_err, &
            & iflag, work, iwork, ray_stop)

       if (verbosity > 2) write(message_unit,'(/,1(a,i4),2(a,g12.6),2(a,1pe10.4))') &
            & ' iflag =', iflag,'  s=', s, '  sout=',sout,  &
            & '  rel_err= ', rel_err,' abs_err = ', abs_err

         if (ray_stop%stop_ode) then  ! stop_ode has been set in eqn_ray or ode itself
            sout=s
            exit odeloop
         end if

         if (iflag == 2) then ! Normal return with s = sout
             exit odeloop
         else if (iflag == 3) then ! rel_err/abs_err have been adjusted try again
         
             total_error = abs(rel_err) + abs(abs_err)
             if (total_error > SG_error_limit) then
               call message('SG_ode: Total error too big', total_error, 0)
               ray_stop%ode_stop_flag = 'ODE total error'
               ray_stop%stop_ode = .true.
               exit odeloop
             end if
             
             cycle odeloop
         else ! error return
              ray_stop%stop_ode = .true.
              call message('SG_ode: Error return: iflag', iflag, 0)
              if (iflag == 5) call text_message('Equation is stiff', 0)
              ray_stop%ode_stop_flag = 'ODE iflag error'
              exit odeloop
         end if

    end do odeloop
    return
  end subroutine SG_ode

 end submodule SG_ode_m

