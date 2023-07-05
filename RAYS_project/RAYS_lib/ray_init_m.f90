 module ray_init_m
 
! Generates intitial position for each ray: rvec0(1:3, iray = 1:nray) and 
! initial refractive index vector for each ray: rindex_vec0(1:3, iray = 1:nray)
!
! How the starting points for the rays are set depends on the plasma geometry e.g
! slab, tokamak, mirror machine, ionosphere, etc and on the geometry of the antenna 
! There are several reasonable ways to initialize k values. In toroidal geometry at low  
! frequencies toroidal and poloidal mode numbers are meaningful.  At higher frequencies the 
! launch angle is more useful e.g. for ECH one launches a beam with a distribution of rays  
! around a central beam axis. Also initialization of k requires solution of the dispersion  
! relation, of which several different ones might be supported.  So initialization depends 
! geometry, antenna model, and dispersion model, at least.  These details are delegated to
! more specialized ray_init modules, this subroutine merely serves as a selector for
! those more specialized functions.  The dispersion solvers typically work with refractive
! index rather than k, so these routines initialize refractive index.  These are converted
! to k later.
!
! So far the ray initialization modules implemented are:
!
! simple_slab_init_m
!
!The action of the ray_init modules is to fill nray, rvec0(:,:) and rindex_vec0(:,:)
!
! ray_pwr_wt(i) = fraction of total power carried by ray i.  Should provide a ray weight
!                 subroutine as part of antenna model.  But for now al wights are just 1/nray.
!
! This module is called from main program.

    use constants_m, only : rkind

    implicit none   
    
    character(len=60) :: ray_init_model

!   Initial position and wavenumber of the ray.  The number of rays to be traced, nray,
!   is calculated in ray_init_m from input values of initial loacations and initial k.
!   space for rvec0 and rindex_vec0 is allocated in ray_init_m

!   Maximum number of rays allowed.
    integer :: nray_max

!   Number of rays.
    integer :: nray

    real(KIND=rkind), allocatable :: rvec0(:,:)
    real(KIND=rkind), allocatable :: rindex_vec0(:,:)
    real(KIND=rkind), allocatable :: ray_pwr_wt(:)

    namelist /ray_init_list/ ray_init_model, nray_max

!****************************************************************************

contains

!****************************************************************************

    subroutine initialize_ray_init_m(read_input)

        use diagnostics_m, only : message_unit, message, text_message
        use simple_slab_ray_init_m, only : simple_slab_ray_init
        use solovev_ray_init_nphi_ntheta_m, only : ray_init_solovev_nphi_ntheta
        use axisym_toroid_ray_init_nphi_ntheta_m, only : ray_init_axisym_toroid_nphi_ntheta
        use axisym_toroid_ray_init_R_Z_nphi_ntheta_m, only : ray_init_axisym_toroid_R_Z_nphi_ntheta
 
        implicit none
        logical, intent(in) :: read_input
 		integer :: input_unit, get_unit_number ! External, free unit finder   
     
        if (read_input .eqv. .true.) then    
        ! Read and write input namelist
  		  	input_unit = get_unit_number()
            open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
            read(input_unit, ray_init_list)
            close(unit=input_unit)
        end if
        write(message_unit, ray_init_list)

        init_model: select case (trim(ray_init_model))
    
            case ('simple_slab')
                call simple_slab_ray_init(nray_max, nray, rvec0, rindex_vec0, ray_pwr_wt)

             case ('solovev')
                call ray_init_solovev_nphi_ntheta(nray_max, nray, rvec0, rindex_vec0, ray_pwr_wt)

             case ('axisym_toroid_nphi_ntheta')
                call ray_init_axisym_toroid_nphi_ntheta(nray_max, nray, rvec0, rindex_vec0,&
                   & ray_pwr_wt)

             case ('axisym_toroid_ray_init_R_Z_nphi_ntheta')
                call ray_init_axisym_toroid_R_Z_nphi_ntheta(nray_max, nray, rvec0,&
                   & rindex_vec0, ray_pwr_wt)

            case default
                write(0,*) 'initialize_ray_init: invalid ray_init_model = ', trim(ray_init_model)
                call text_message('initialize_ray_init: invalid ray_init_model = ', trim(ray_init_model),0)
                stop 1

        end select init_model
 
    return
    end subroutine initialize_ray_init_m

!********************************************************************

    subroutine deallocate_ray_init_m
        use simple_slab_ray_init_m, only : deallocate_simple_slab_ray_init_m
        use solovev_ray_init_nphi_ntheta_m, only : deallocate_solovev_ray_init_nphi_ntheta_m
        use axisym_toroid_ray_init_nphi_ntheta_m, only :&
          & deallocate_axisym_toroid_ray_init_nphi_ntheta_m
        use axisym_toroid_ray_init_R_Z_nphi_ntheta_m, only :&
          & deallocate_axisym_toroid_ray_init_R_Z_nphi_ntheta_m

        call deallocate_simple_slab_ray_init_m
        call deallocate_solovev_ray_init_nphi_ntheta_m
        call deallocate_axisym_toroid_ray_init_nphi_ntheta_m
        call deallocate_axisym_toroid_ray_init_R_Z_nphi_ntheta_m
        
        return
    end subroutine deallocate_ray_init_m

 end module ray_init_m
