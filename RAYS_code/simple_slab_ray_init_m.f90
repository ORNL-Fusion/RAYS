 module simple_slab_ray_init_m
! Contains data to initialize rays for a simple slab equilibrium (i.e. all variations are
! in the x coordinate), 
! Contains a subroutine to do the initialization: simple_slab_ray_init,
! which generates initial positions, rvec0 = (x0, 0., 0 : nray), and initial refractive
! index vector, rindex_vec0 = (nx0, ny0, nz0 : nray).

! External procedures: solve_disp_nx_vs_ny_nz (solve_disp_nx_vs_ny_nz.f90)
    use constants_m, only : rkind
    
    implicit none

    integer:: n_x_launch = 1
    real(KIND=rkind) ::  x_launch0 = 0., dx_launch = 0.
    integer:: n_y_launch = 1
    real(KIND=rkind) ::  y_launch0 = 0., dy_launch = 0.
    integer:: n_z_launch = 1
    real(KIND=rkind) ::  z_launch0 = 0., dz_launch = 0.

    integer:: n_ky_launch, n_kz_launch
    real(KIND=rkind) ::  rindex_y0, delta_rindex_y0, rindex_z0, delta_rindex_z0
     
 namelist /simple_slab_ray_init_list/ &
     & n_x_launch, x_launch0, dx_launch, n_y_launch, y_launch0, dy_launch, &
     & n_z_launch, z_launch0, dz_launch, n_ky_launch, rindex_y0,           &
     & delta_rindex_y0, n_kz_launch, rindex_z0, delta_rindex_z0

!****************************************************************************

contains

!****************************************************************************

    
    subroutine simple_slab_ray_init(nray_max, nray, rvec0, rindex_vec0)     
! initializer for slab geometry

! Choose initial y and z refractive indices (rindex_y0, delta_rindex_y0) and
! solve dispersion relation for rindex_x0)
!
! N.B. Some of the ray initializations may fail (e.g. initial point is outside plasma or 
!      wave mode is evanescent).  This does not cause the program to stop.  It counts
!      the successful initializations and sets number of rays, nray, to that.

    use constants_m, only : input_unit
    use diagnostics_m, only: message_unit, message, text_message
    use species_m, only : nspec
    use equilibrium_m, only : equilibrium, eq_point
    use dispersion_solvers_m, only: solve_disp_nx_vs_ny_nz
    use rf_m, only : ray_dispersion_model,wave_mode, k0_sign

    implicit none
    
    integer, intent(in) :: nray_max
    integer, intent(out) :: nray
    type(eq_point(nspec=nspec)) :: eq
    real(KIND=rkind), allocatable, intent(out) :: rvec0(:, :), rindex_vec0(:, :)
    
    integer :: ix, iy, iz, iky, ikz, count
    real(KIND=rkind) :: x, y, z, rindex_y, rindex_z
    real(KIND=rkind) :: rvec(3)
    complex(KIND=rkind) :: rindex_x

! Read and write input namelist
    open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
    read(input_unit, simple_slab_ray_init_list)
    close(unit=input_unit)
    write(message_unit, simple_slab_ray_init_list)

! allocate space for the initial condition vectors rvec0, rindex_vec0
    nray = n_x_launch * n_ky_launch * n_kz_launch

        if (nray > 0) then
            allocate (rvec0(3, nray), rindex_vec0(3, nray))
        else
            write (6,*) 'ray_init_slab: invalid number of rays  nray=', nray
            stop 1
        end if  

! Load up initial position and k vectors for each ray.  Count successful initializations.
    count = 0
    zloop: do iz = 1, n_z_launch
        z = z_launch0  + (iz-1) * dy_launch   
    yloop: do iy = 1, n_y_launch
        y = y_launch0  + (iy-1) * dy_launch   
    xloop: do ix = 1, n_x_launch
        x = x_launch0  + (ix-1) * dx_launch   

        rvec( : ) = (/ x, y, z /)  

        kyloop: do iky = 1, n_ky_launch
            rindex_y = rindex_y0 + (iky-1) * delta_rindex_y0

            kzloop: do ikz = 1, n_kz_launch   
            
                if (count > nray_max) then
                  write(0,*) 'simple_slab_ray_init: ray count exceeds nray_max'
                  call text_message('simple_slab_ray_init: ray count exceeds nray_max',0)
                  stop 1
                end if
                                  
                rindex_z = rindex_z0 + (ikz-1) * delta_rindex_z0  

                call equilibrium(rvec, eq)
                   if (trim(eq%equib_err) /= '') cycle kzloop
                call solve_disp_nx_vs_ny_nz(eq, ray_dispersion_model, wave_mode, k0_sign,&
                     &  rindex_y, rindex_z, rindex_x)
                if (aimag(rindex_x) /= 0.) then
                    write(message_unit, *) 'slab_init: evanescent ray x = ', x, &
                    & ' ny = ', rindex_y, ' nz = ', rindex_z
               
                    cycle kzloop
                end if

                count = count +1
                rvec0( : , count) = rvec
                rindex_vec0( : , count) = (/ real(rindex_x, KIND=rkind), rindex_y, rindex_z /)

            end do kzloop 
        end do kyloop
    end do xloop
    end do yloop
    end do zloop

    nray = count
    call message('simple_slab_ray_init: nray', nray)
    if (nray == 0) stop 'No successful ray initializations' 
    end  subroutine simple_slab_ray_init 

end module simple_slab_ray_init_m
