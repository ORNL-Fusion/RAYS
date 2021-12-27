 module slab_processor_m
! Post processing for slab equilibrium, all variation in x direction

    use constants_m, only : rkind

    implicit none

    character(len=80) :: processor = ''
    integer, parameter :: graphics_descrip_unit = 96
    integer, parameter :: eq_profile_unit = 97
    integer, parameter :: kx_profile_unit = 98
    integer, parameter :: res_and_cut_unit = 99

!   Number and x locations of cyclotron resonances, 2nd harmonic resonances and hybrid
!   resonances found between x_min and x_max.  Allow for multiple resonances <= n_locs
    integer, parameter :: n_locs = 5
    integer :: n_ce_res, n_2ce_res, n_hybrid_res
    real(KIND=rkind), dimension(n_locs) :: x_ce_res, x_2ce_res, x_hybrid_res

!   Number and x locations of cutoffs
    integer :: n_P_cut, n_H_cut, n_det
    real(KIND=rkind), dimension(n_locs) :: x_P_cut, x_H_cut, x_det

! Alpha at hybrid cutoffs
    real(KIND=rkind), dimension(n_locs) :: alpha_e_H_cut
    real(KIND=rkind), dimension(n_locs) :: alpha_e_det

! Number of k vectors to plot for each ray in graphics
    integer :: num_plot_k_vectors

! Scale plot k vectors to kmax on the ray, (True, False)
    character(len = 5) :: scale_k_vec = 'True'

! Set plot xlim -> [xmin, xmax] and ylim -> [ymin,ymax], (True, False)
    character(len = 5) :: set_XY_lim = 'True'

    namelist /slab_processor_list/ processor, num_plot_k_vectors, scale_k_vec, set_XY_lim

 contains

 subroutine initialize_slab_processor

    use diagnostics_m, only : message_unit, message, text_message
    use constants_m, only : input_unit

    implicit none

! Read and write input namelist
    open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
    read(input_unit, slab_processor_list)
    close(unit=input_unit)
    write(message_unit, slab_processor_list)
    call text_message('Finished initialize_slab_processor ', processor)
    
    return
 end subroutine initialize_slab_processor

!*************************************************************************     

  subroutine slab_processor

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message

    implicit none
    
    call write_eq_profiles
    
    call write_kx_profiles
    
    call find_res_and_cuts
    
    call write_graphics_description_file
    
    call text_message('Finished slab_processor work')

    return
 end subroutine slab_processor


!*************************************************************************     

  subroutine find_res_and_cuts
! Finds x locations of cold plasma resonances and cutoffs.  Now only finds electron cyclotron
! resonances.  Later include arbitrary ion species.  Simple scan algorithm, no root finding.  
! Maybe later.
!
! Nota Bene: B is assumed in y-z plane and ny, nz are constant because of x stratification.
!
! Restriction:  This assumes that the component of k in the y-x plane is parallel to B for
! all x.  Equivalent to B along z (no shear) and ky = 0.  If not then there is 
! component of k_perp in the y-z plane, nz is not n_parallel and the cutoff (kx = 0) is in
! a different place.  Could be generalized.  Maybe later.

    use constants_m, only : rkind
    use diagnostics_m, only : run_label
    use equilibrium_m, only : equilibrium, bunit, alpha, gamma
    use slab_eq_m, only : xmin, xmax
    use suscep_m, only : dielectric_cold, eps_cold
    use ray_init_m, only : nray, rindex_vec0
    
    implicit none

    integer, parameter :: n_xpoints = 1000 ! Number of x points in scan
    real(KIND=rkind) :: v_ce_0, v_ce_1, v_2ce_0, v_2ce_1, v_hybrid_0, v_hybrid_1
    real(KIND=rkind) :: v_P_cut_0, v_P_cut_1,  v_H_cut_0, v_H_cut_1
    real(KIND=rkind) :: a, b, c, v_det_0, v_det_1 ! parameters for calculating determinant
    integer :: iray, ix
    real(KIND=rkind) :: nz, x, dx
    complex :: vH


	open(unit = res_and_cut_unit, file = 'res_and_cut.'//trim(run_label), &
			& action='write', status='replace', form='formatted')

    dx = (xmax - xmin)/(n_xpoints-1)

    ray_loop: do iray = 1, nray
        write(*,*) ' '
        write(*,*) 'ray ', iray

        write(res_and_cut_unit,*) ' '
        write(res_and_cut_unit,*) 'ray ', iray

        n_ce_res = 0
        n_2ce_res = 0
        n_hybrid_res = 0
        n_P_cut = 0
        n_H_cut = 0
        n_det = 0
        
        x_ce_res = 0.
        x_2ce_res = 0.
        x_hybrid_res = 0.
        x_P_cut = 0.
        x_H_cut = 0.
        x_det = 0.
    
        nz = sum(rindex_vec0(:, iray)*bunit(:))
        x_loop: do ix = 0, n_xpoints-1
            x = xmin + ix*dx
        
            call equilibrium( (/real(x, KIND=rkind), real(0., KIND=rkind), real(0., KIND=rkind)/) )
            call dielectric_cold
            v_ce_1 = gamma(0)+1. ! remember gamma(0) is negative
            v_2ce_1 = gamma(0)+0.5
            v_hybrid_1 = real(eps_cold(1,1), KIND=rkind)
            v_P_cut_1 = real(eps_cold(3,3), KIND=rkind)
            v_H_cut_1 = real(( eps_cold(1,1)**2+eps_cold(1,2)**2 - 2.*eps_cold(1,1)*nz**2 &
                    & + nz**4 ), KIND=rkind)
 
		   a = real( eps_cold(1,1), KIND=rkind )
		   b = real( -(eps_cold(1,1)**2+eps_cold(1,2)**2+eps_cold(1,1)*eps_cold(3,3)) &
			  & + (eps_cold(1,1)+eps_cold(3,3))*nz**2, KIND=rkind )
		   c = real( ( eps_cold(1,1)**2+eps_cold(1,2)**2 - 2.*eps_cold(1,1)*nz**2 &
			  & + nz**4 ) *	 eps_cold(3,3), KIND=rkind )
		  
		   v_det_1 = b**2-4.*a*c
           
!             vH = eps_cold(1,1)**2+eps_cold(1,2)**2 - 2.*eps_cold(1,1)*nz**2 + nz**4
!             write(*,*) 'x = ', x, ' eps_cold(1,1 ) = ', eps_cold(1,1), ' v_P_cut_1 = ',&
!              & v_P_cut_1, 'v_H = ', vH
                   
            if (ix == 0) then ! first time initialize v_ce0 etc and skip to next point
                v_ce_0 = v_ce_1
                v_2ce_0 = v_2ce_1
                v_hybrid_0 = v_hybrid_1
                v_P_cut_0 = v_P_cut_1
                v_H_cut_0 = v_H_cut_1
                v_det_0 = v_det_1
                cycle x_loop
            end if
            
            ! Check for sign change between v_0 and v_1
            if (v_ce_0 * v_ce_1 <= 0.) then
              n_ce_res = n_ce_res + 1
              x_ce_res(n_ce_res) = x
            end if
            if (v_2ce_0 * v_2ce_1 <= 0.) then
              n_2ce_res = n_2ce_res + 1
              x_2ce_res(n_2ce_res) = x
            end if
            if (v_hybrid_0 * v_hybrid_1 <= 0.) then
              n_hybrid_res = n_hybrid_res + 1
              x_hybrid_res(n_hybrid_res) = x
            end if
            if (v_P_cut_0 * v_P_cut_1 <= 0.) then
              n_P_cut = n_P_cut + 1
              x_P_cut(n_P_cut) = x
            end if
            if (v_H_cut_0 * v_H_cut_1 <= 0.) then
              n_H_cut = n_H_cut + 1
              x_H_cut(n_H_cut) = x
              alpha_e_H_cut(n_H_cut) = alpha(0)
            end if
            if (v_det_0 * v_det_1 <= 0.) then
              n_det = n_det + 1
              x_det(n_det) = x
              alpha_e_det(n_det) = alpha(0)
            end if
            
            ! advance v_0
            v_ce_0 = v_ce_1
            v_2ce_0 = v_2ce_1
            v_hybrid_0 = v_hybrid_1
            v_P_cut_0 = v_P_cut_1
            v_H_cut_0 = v_H_cut_1
            v_det_0 = v_det_1
            
        end do x_loop
   
        write(*,*) 'x_ce_res = ', x_ce_res
        write(*,*) 'x_2ce_res = ', x_2ce_res
        write(*,*) 'x_hybrid_res = ', x_hybrid_res
        write(*,*) 'x_P_cut = ', x_P_cut
        write(*,*) 'x_H_cut = ', x_H_cut
        write(*,*) 'alpha_e_H_cut = ', alpha_e_H_cut
        write(*,*) 'x_det = ', x_det
        write(*,*) 'alpha_e_det = ', alpha_e_det

        write(res_and_cut_unit,*) 'x_ce_res = ', x_ce_res
        write(res_and_cut_unit,*) 'x_2ce_res = ', x_2ce_res
        write(res_and_cut_unit,*) 'x_hybrid_res = ', x_hybrid_res
        write(res_and_cut_unit,*) 'x_P_cut = ', x_P_cut
        write(res_and_cut_unit,*) 'x_H_cut = ', x_H_cut
        write(res_and_cut_unit,*) 'alpha_e_H_cut = ', alpha_e_H_cut
        write(res_and_cut_unit,*) 'x_det = ', x_det
        write(res_and_cut_unit,*) 'alpha_e_det = ', alpha_e_det
     
    end do ray_loop
    close(unit = res_and_cut_unit)
    
  end subroutine find_res_and_cuts

!*************************************************************************     

  subroutine write_eq_profiles
! Writes equilibrium profiles and other stuff for plotting
! For ease of reading these are re-cast to single precision before writing

    use constants_m, only : rkind, skind, pi
    use diagnostics_m, only : run_label
    use equilibrium_m, only : equilibrium, bmag, bunit, ns, omgc, omgp2, alpha, gamma
    use slab_eq_m, only : xmin, xmax
    use suscep_m, only : dielectric_cold, eps_cold
    use dispersion_solvers_m, only : solve_disp_nx_vs_ny_nz    
    implicit none

!   fast, slow etc
    character(len=15) :: wave_mode

    integer, parameter :: n_xpoints = 101 ! Number of x points in scan
    integer :: iray, ix
    real(KIND=rkind) :: x, dx
    real(KIND=skind), dimension(9) :: profile_vec
    character(len = 9), dimension(9) :: prof_name
    character(len = 10), parameter  :: b10 = '          '
    character(len = 9), parameter  ::   b9 = '         '
    character(len = 8), parameter  ::   b8 = '        '
    character(len = 7), parameter  ::   b7 = '       '
    
    write(*,*)  'Start writing eq profile vectors'  

    prof_name = (/'x        ', 'ne       ', 'bmag     ', 'fc_e     ', 'fp_e     ', &
                       & 'fp_e/f_rf', 'fc_e/f_rf', 'alpha    ', 'beta     '/)
        
    open(unit = eq_profile_unit, file = 'eq_profiles_slab.'//trim(run_label))
   
    write(eq_profile_unit,*) ' ',prof_name(1), b9, prof_name(2), b8, prof_name(3), b8,&
                     & prof_name(4), b8, prof_name(5), b7, prof_name(6), b8, prof_name(7),&
                     & b8, prof_name(8), b8, prof_name(9)

    dx = (xmax - xmin)/(n_xpoints-1)
    x_loop: do ix = 0, n_xpoints-1

        x = xmin + ix*dx
        
        call equilibrium( (/real(x, KIND=rkind), real(0., KIND=rkind), real(0., KIND=rkind)/) )
        
        profile_vec(1) = real(x, KIND=skind)
        profile_vec(2) = real(ns(0), KIND=skind) ! electron density
        profile_vec(3) = real(bmag, KIND=skind) ! total B field
        profile_vec(4) = real(abs(omgc(0))/(2.*pi), KIND=skind) ! electron cyclotron frequency
        profile_vec(5) = real(sqrt(omgp2(0))/(2.*pi), KIND=skind) ! electron plasma frequency
        profile_vec(6) = real(sqrt(alpha(0)), KIND=skind) ! f_pe/f_rf
        profile_vec(7) = real(abs(gamma(0)), KIND=skind) ! f_ce/f_rf
        profile_vec(8) = real(alpha(0), KIND=skind) ! (f_pe/f_rf)^2
        profile_vec(9) = real(gamma(0)**2, KIND=skind) ! (f_ce/f_rf)^2

        write(eq_profile_unit,*) profile_vec

    end do x_loop
    
    write(*,*)  'Finished writing eq profile vectors'  
    close(eq_profile_unit)   
    
  end subroutine write_eq_profiles

!*************************************************************************     

  subroutine write_kx_profiles
! Writes kx roots versus x.
! For ease of reading these are re-cast to single precision before writing

    use constants_m, only : rkind, skind, pi
    use diagnostics_m, only : run_label
    use equilibrium_m, only : equilibrium, bmag, bunit, ns, omgc, omgp2, alpha, gamma
    use slab_eq_m, only : xmin, xmax
    use suscep_m, only : dielectric_cold, eps_cold
    use ray_init_m, only : nray, rindex_vec0
    use rf_m, only : ray_dispersion_model, k0, k0_sign
    use dispersion_solvers_m, only : solve_disp_nx_vs_ny_nz    
    implicit none

!   fast, slow etc
    character(len=15) :: wave_mode

    integer, parameter :: n_xpoints = 101 ! Number of x points in scan
    integer :: iray, ix
    real(KIND=rkind) :: x, dx
    real(KIND=rkind) :: ny, nz
    real(KIND=skind), dimension(9) :: profile_vec
    character(len = 17), dimension(9) :: profile_name_vec
    complex(KIND=rkind) :: nx
    complex(KIND=skind) :: nx_sngl

    write(*,*)  'Start writing kx profile vectors'  

    profile_name_vec = (/'x                ', 'kx_real_plus     ', 'kx_im_plus       ',&
                         'kx_real_minus    ', 'kx_im_minus      ', 'kx_real_fast     ',&
                         'kx_im_fast       ', 'kx_real_slow     ', 'kx_im_slow       '/)
                              
    open(unit = kx_profile_unit, file = 'kx_profiles_slab.'//trim(run_label))
   

    dx = (xmax - xmin)/(n_xpoints-1)

    ray_loop: do iray = 1, nray
        write(*,*) 'ray ', iray
        write(kx_profile_unit,*) 'ray ', iray
        write(kx_profile_unit,*) profile_name_vec

    
        ny = rindex_vec0(2, iray)
        nz = rindex_vec0(3, iray)
 
        x_loop: do ix = 0, n_xpoints-1
            x = xmin + ix*dx
             
            call equilibrium( (/real(x, KIND=rkind), real(0., KIND=rkind), real(0., KIND=rkind)/) )

            profile_vec(1) = real(x, KIND=skind)
           
            ! kx vs x for fast and slow cold plasma roots
            wave_mode = 'plus'
            call solve_disp_nx_vs_ny_nz(ray_dispersion_model, wave_mode, k0_sign, ny, nz, nx)            
            nx_sngl = cmplx(nx, KIND=skind)
            profile_vec(2) = k0*real(nx_sngl)
            profile_vec(3) = k0*aimag(nx_sngl)
            
            wave_mode = 'minus'
            call solve_disp_nx_vs_ny_nz(ray_dispersion_model, wave_mode, k0_sign, ny, nz, nx)            
            nx_sngl = cmplx(nx, KIND=skind)
            profile_vec(4) = k0*real(nx_sngl)
            profile_vec(5) = k0*aimag(nx_sngl)
            ! kx vs x for fast and slow cold plasma roots

            wave_mode = 'fast'
            call solve_disp_nx_vs_ny_nz(ray_dispersion_model, wave_mode, k0_sign, ny, nz, nx)            
            nx_sngl = cmplx(nx, KIND=skind)
            profile_vec(6) = k0*real(nx_sngl)
            profile_vec(7) = k0*aimag(nx_sngl)
            
            wave_mode = 'slow'
            call solve_disp_nx_vs_ny_nz(ray_dispersion_model, wave_mode, k0_sign, ny, nz, nx)            
            nx_sngl = cmplx(nx, KIND=skind)
            profile_vec(8) = k0*real(nx_sngl)
            profile_vec(9) = k0*aimag(nx_sngl)

            write(kx_profile_unit,*) profile_vec

        end do x_loop
     
    end do ray_loop
    
    write(*,*)  'Finished writing kx profile vectors'  
    close(kx_profile_unit)   
    
  end subroutine write_kx_profiles

!*************************************************************************     

  subroutine write_graphics_description_file
  
   use diagnostics_m, only : run_description, run_label
   use slab_eq_m, only : xmin, xmax, ymin, ymax, zmin, zmax
   
   open(unit = graphics_descrip_unit, file = 'graphics_description_slab.dat')
  
   write(graphics_descrip_unit, *) 'run_description = ', run_description
   write(graphics_descrip_unit, *) 'run_label = ', run_label
  
   write(graphics_descrip_unit, *) 'xmin = ', xmin
   write(graphics_descrip_unit, *) 'xmax = ', xmax
   write(graphics_descrip_unit, *) 'ymin = ', ymin
   write(graphics_descrip_unit, *) 'ymax = ', ymax
   write(graphics_descrip_unit, *) 'zmin = ', zmin
   write(graphics_descrip_unit, *) 'zmax = ', zmax
  
   write(graphics_descrip_unit, *) 'num_plot_k_vectors = ', num_plot_k_vectors
   write(graphics_descrip_unit, *) 'scale_k_vec = ', trim(scale_k_vec)
   write(graphics_descrip_unit, *) 'set_XY_lim = ', trim(set_XY_lim)

   close(unit = graphics_descrip_unit)

  end subroutine write_graphics_description_file


 end module slab_processor_m

