  subroutine classify_mode_OX(eq, dispersion_model, n_vec, mode)
! Checks whether input n_vec satisfies dispersion relation for O or X mode. Can be
! either, or both, or neither.  Output mode is a logical 2-vector
! It does this by calculating 2-norm of n_vec and comparing that to dispersion solution.
! Tolerance for comparison is set in parameter below.

    use constants_m, only : rkind, one
    use diagnostics_m, only : message, text_message
    use species_m, only : nspec
    use equilibrium_m, only : eq_point

    implicit none

! Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point(nspec=nspec)), intent(in) :: eq

!   Name of dispersion model (can be different from that used in ray tracing)
    character(len=20), intent(in) :: dispersion_model

! Input n vector
    real(KIND=rkind), intent(in) :: n_vec(3)

! Output mode vector
	logical :: mode(2)

    real(KIND=rkind) :: nunit(3), nmag, nmag2, theta
    real(KIND=rkind) :: nsq(4) !n perp square

	real(KIND=rkind), parameter :: tolerance = one/10.0_rkind**6

	mode =(/.false., .false./)

! normalize n_vec
    nmag = norm2(n_vec)
    nmag2 = nmag**2
    nunit = n_vec/nmag
	theta = acos(dot_product(nunit,eq%bunit))

    dispersion_relation: select case (trim(dispersion_model))

       case ('cold')
          ! Solve for n squared
          call solve_cold_nsq_vs_theta(eq, theta, nsq)
          if (abs(nmag - sqrt(abs(nsq(1))))  <= tolerance) mode(1) = .true.
          if (abs(nmag - sqrt(abs(nsq(2))))  <= tolerance) mode(2) = .true.
!          if (sqrt(abs(nmag2 - nsq(2)))  <= tolerance) mode(2) = .true.
!  write(*,*) 'classify_mode_OX: n_vec(1) = ', n_vec(1),&
!   & '  sin(theta)*sqrt(nsq(1)) = ', sin(theta)*sqrt(nsq(1)),&
!   & '  sin(theta)*sqrt(nsq(2)) = ', sin(theta)*sqrt(nsq(2))
       case default
          write(0,*) 'classify_mode_OX: unimplemented dispersion_model =', trim(dispersion_model)
          call text_message('classify_mode_OX: unimplemented dispersion_model ',trim(dispersion_model),0)
          stop 1

    end select dispersion_relation

    return
  end  subroutine classify_mode_OX
