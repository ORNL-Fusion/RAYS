 module dispersion_solvers_m
! Contains routines to solve various forms of the plasma dispersion relation
! N.B. These routines require that the plasma magnetic and species quantities
! have been set by a previous call to subroutine equilibrium(rvec)
! N.B. remember convention
!   n1 = perpedicular component of nvec.
!   n12 = square of perpendicular component of nvec
!   n3 = parallel component of kvec (nvec).

! External procedures: disp_solve_cold_nxsq_vs_nz (disp_solve_cold_nxsq_vs_nz.f90)

    use constants_m, only : rkind
    
    implicit none

!****************************************************************************

contains

!****************************************************************************

    
  subroutine solve_disp_nx_vs_ny_nz(eq, dispersion_model, wave_mode, k_sign, ny, nz, nx)    
! N.B. This requires that the magnetic field be in the y-z plane, otherwise n_parallel 
!      has a component along nx.  Complicates dispersion relation.
! N.B. nx is complex(KIND=rkind)

    use diagnostics_m, only : message, text_message
    use species_m, only : nspec
    use equilibrium_m, only : eq_point

    implicit none

! Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point(nspec=nspec)), intent(in) :: eq

!   Name of dispersion model (can be different from that used in ray tracing)
    character(len=15), intent(in) :: dispersion_model

!   fast, slow etc
    character(len=15), intent(in) :: wave_mode
    
    integer, intent(in) :: k_sign
    real(KIND=rkind), intent(in) :: ny, nz
    complex(KIND=rkind), intent(out) :: nx

    real(KIND=rkind) :: n3  ! n parallel
    complex(KIND=rkind) :: n12(4) !n perp square
    integer :: i_mode  ! plus mode -> 1, minus mode -> 2, fast mode -> 3, slow mode -> 3

    n3 = ny*eq%bunit(2) + nz*eq%bunit(3)
 
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
          call disp_solve_cold_nxsq_vs_nz(eq, n3, n12) 
          nx = k_sign*sqrt(n12(i_mode))

       case default
          write(0,*) 'solve_disp: unimplemented dispersion_model =', trim(dispersion_model)
          call text_message('solve_disp: unimplemented dispersion_model ',trim(dispersion_model),0)
          stop 1

    end select dispersion_relation

    return
  end  subroutine solve_disp_nx_vs_ny_nz 
    

end module dispersion_solvers_m
