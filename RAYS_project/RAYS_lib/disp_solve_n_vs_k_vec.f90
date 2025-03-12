  subroutine solve_n_vs_k_vec(eq, dispersion_model, wave_mode, k_sign, k_vec, n_out)
! Solves for n in the direction of k_vec
! N.B. n_out is complex(KIND=rkind).  Calling program must account for that.

    use constants_m, only : rkind, one
    use diagnostics_m, only : message, text_message
    use species_m, only : nspec
    use equilibrium_m, only : eq_point

    implicit none

! Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point(nspec=nspec)), intent(in) :: eq

!   Name of dispersion model (can be different from that used in ray tracing)
    character(len=20), intent(in) :: dispersion_model

!   fast, slow etc
    character(len=20), intent(in) :: wave_mode

    integer, intent(in) :: k_sign
    real(KIND=rkind), intent(in) :: k_vec(3) ! Input k vector
    complex(KIND=rkind), intent(out) :: n_out

    real(KIND=rkind) :: kunit(3), kmag, theta
    complex(KIND=rkind) :: nsq(4) !n perp square
    integer :: i_mode  ! plus mode -> 1, minus mode -> 2, fast mode -> 3, slow mode -> 3

! normalize k_vec
    kmag = sqrt( sum(kunit**2) )
    kunit = kunit/kmag

	theta = acos(dot_product(kunit,eq%bunit))
    mode: select case (trim(wave_mode))

       case ('plus')
          i_mode = 1

       case ('minus')
          i_mode = 2

       case ('fast')
          i_mode = 3

       case ('slow')
          i_mode = 4

       case default
          write(0,*) 'solve_disp: improper wave_mode = ', trim(wave_mode)
          call text_message('solve_disp: improper wave_mode = ', trim(wave_mode),0)
          stop 1

     end select mode

     dispersion_relation: select case (trim(dispersion_model))

       case ('cold')
          ! Solve for n squared
          call solve_cold_nsq_vs_theta(eq,theta, nsq)
          n_out = real(k_sign, kind=rkind) * sqrt(nsq(i_mode))

       case default
          write(0,*) 'solve_disp: unimplemented dispersion_model =', trim(dispersion_model)
          call text_message('solve_disp: unimplemented dispersion_model ',trim(dispersion_model),0)
          stop 1

    end select dispersion_relation

    return
  end  subroutine solve_n_vs_k_vec
