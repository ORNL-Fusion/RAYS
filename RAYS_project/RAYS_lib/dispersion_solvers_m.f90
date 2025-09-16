 module dispersion_solvers_m
! Contains routines to solve various forms of the plasma dispersion relation
!
! The particular model to be solved is specified by argument "dispersion_model"
! Presently one model is supported: cold, (warm real soon now)
!
! Presently two forms for this model are supported:
! 1) solve_n1_vs_n2_n3
!   n1 = component of nvec perpendicular to B to be solved for (e.g. x or radial component)
!   n2 = transverse component of nvec i.e. perpendicular to both B and n1
!   n3 = parallel component of kvec (nvec)
!
!2) solve_n_vs_theta
!   n = scalar refractive index
!	theta = angle between vector refractive index and B
!   This is essentially Appleton_Hartree
!
! The user also specifies which mode to select (wave_mode).
! Cold plasma presents 4 possibilites -> plus/minus in solution of quadratic, or
! fast/slow meaning smaller/larger magnitude of refractive index.
! For warm plasma there will be other modes (e.g. Bernstein).  Haven't figured that out yet.
!
! N.B. The plasma magnetic and species quantities come in through argument 'eq' which
! is a derived type 'eq_point'

! N.B. refractive index convention
!   nperp_sq = square of perpendicular component of nvec = n1**2 +n2**2

! External procedures:
! disp_solve_cold_nxsq_vs_nz in disp_solve_cold_nxsq_vs_nz.f90
! solve_cold_nsq_vs_theta in solve_cold_nsq_vs_theta.f90

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use constants_m, only : rkind

    implicit none

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

  subroutine solve_n1_vs_n2_n3(eq, dispersion_model, wave_mode, k_sign, n2, n3, n1)

! N.B. n1 output is complex(KIND=rkind).  Calling program must account for that.

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
    real(KIND=rkind), intent(in) :: n2 ! n transverse
    real(KIND=rkind), intent(in) :: n3 ! n parallel
    complex(KIND=rkind), intent(out) :: n1

    complex(KIND=rkind) :: nperp_sq(4) !n perp square
    integer :: i_mode  ! plus mode -> 1, minus mode -> 2, fast mode -> 3, slow mode -> 3

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
          ! Solve for n-perp squared
          call solve_cold_n1sq_vs_n3(eq, n3, nperp_sq)
          n1 = real(k_sign, kind=rkind) * sqrt(nperp_sq(i_mode)-n2**2)

       case default
          write(0,*) 'solve_disp: unimplemented dispersion_model =', trim(dispersion_model)
          call text_message('solve_disp: unimplemented dispersion_model ',trim(dispersion_model),0)
          stop 1

    end select dispersion_relation

    return
  end  subroutine solve_n1_vs_n2_n3

!****************************************************************************

  subroutine solve_nx_vs_ny_nz_by_bz(eq, dispersion_model, wave_mode, k_sign, ny, nz, nx)

! Convenience routine that resolves ny, nz into n2 and n3 (n transverse and n parallel) and
! calls subroutine above.

! N.B. This requires that the magnetic field be in the y-z plane, otherwise n_parallel
!      has a component along nx.  That would make dispersion relation 4th order instead of
!      bi-quadratic

! N.B. nx output is complex(KIND=rkind).  Calling program must account for that.

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
    real(KIND=rkind), intent(in) :: ny, nz
    complex(KIND=rkind), intent(out) :: nx

    real(KIND=rkind) :: n2  ! n transverse
    real(KIND=rkind) :: n3  ! n parallel

    n2 = ny*eq%bunit(3) - nz*eq%bunit(2)
    n3 = ny*eq%bunit(2) + nz*eq%bunit(3)

    dispersion_relation: select case (trim(dispersion_model))

       case ('cold')
          ! Solve for n-perp squared
          call solve_n1_vs_n2_n3(eq, dispersion_model, wave_mode, k_sign, n2, n3, nx)

       case default
          write(0,*) 'solve_disp: unimplemented dispersion_model =', trim(dispersion_model)
          call text_message('solve_disp: unimplemented dispersion_model ',trim(dispersion_model),0)
          stop 1

    end select dispersion_relation

    return
  end  subroutine solve_nx_vs_ny_nz_by_bz
!****************************************************************************

  subroutine solve_n_vs_theta(eq, dispersion_model, wave_mode, k_sign, theta, n_out)
! Appleton_Hartree like solver
! N.B. n_out output is complex(KIND=rkind).  Calling program must account for that.

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
    real(KIND=rkind), intent(in) :: theta ! Angle between vector n and B
    complex(KIND=rkind), intent(out) :: n_out

    real(KIND=rkind) :: nsq(4) ! N.B. n square, might be negative
    integer :: i_mode  ! plus mode -> 1, minus mode -> 2, fast mode -> 3, slow mode -> 3

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
  end  subroutine solve_n_vs_theta

!********************************************************************

    subroutine deallocate_dispersion_solvers_m
        return ! nothing to deallocate
    end subroutine deallocate_dispersion_solvers_m

end module dispersion_solvers_m
