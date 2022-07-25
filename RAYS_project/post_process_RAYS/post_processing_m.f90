 module post_processing_m

    use constants_m, only : rkind

    implicit none
    
! Calculated below from data in input files
    integer :: npoints_max
    integer, allocatable :: npoints(:)
    real(KIND=rkind), allocatable :: s_vec(:,:), v_vec(:,:,:)

! Switch to select specific post processor
    character(len=80) :: processor = ''

    namelist /post_process_list/ processor

 contains

 subroutine initialize_post_processing_m(read_input)

    use diagnostics_m, only : message_unit, message, text_message, verbosity
    use constants_m, only : input_unit, output_unit, ray_list_unit
    use slab_processor_m, only : initialize_slab_processor
    use solovev_processor_m, only : initialize_solovev_processor
    use axisym_toroid_processor_m, only : initialize_axisym_toroid_processor
    use ray_init_m, only : nray_ray_init => nray

    implicit none
    logical, intent(in) :: read_input

    integer :: nray, nv, iray, ipoint
    real(KIND=rkind) :: s
    real(KIND=rkind), allocatable :: v(:)
    real(KIND=rkind), allocatable :: residuals(:) 
    character(len = 20), allocatable :: ray_stop(:)

    if (read_input .eqv. .true.) then    
	! Read and write input namelist
        open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
        read(input_unit, post_process_list)
        close(unit=input_unit)
        write(message_unit, post_process_list)
	end if

    select case (trim(processor))

       case ('slab')
          call initialize_slab_processor(read_input)

       case ('solovev')
          call initialize_solovev_processor(read_input)

       case ('axisym_toroid')
          call initialize_axisym_toroid_processor(read_input)

       case default
          write(*,*) 'post_process_rays: unimplemented post_processor =', trim(processor)
          call text_message('post_process_rays: unimplemented post_processor', trim(processor))
          stop 1

       end select

!****** Read in all ray data *********************************    

!  Read number of rays and number of points points on each ray from binary file ray_list.bin
!  Then read the ray data from binary file rays.bin
!  Note: files ray_list.bin = unit 95 and rays.bin = unit 94 are opened in subroutine 
!  initialize()

    call text_message('initialize_post_process_rays', 3)

!    read (95) nray this is a read for binary file.  For now use formatted file instead

    read(ray_list_unit, *) nray
    allocate(npoints(nray))
    read(ray_list_unit, *) nv
    allocate(v(nv))    

! Check for consistency between nray and the value obtained in ray_init.  This is a weak 
! check that the rays.in and the ray_list.bin file are from the same run.
    If (nray .ne. nray_ray_init) then
        write (*,*) 'initialize_post_processing: nray = ', nray,&
        & ' inconsistent with nray_ray_init = ', nray_ray_init
        write (message_unit,*) 'initialize_post_processing: nray = ', nray,&
        & ' inconsistent with nray_ray_init = ', nray_ray_init
        !stop 1
    end if

    read(ray_list_unit, *) npoints
    npoints_max = maxval(npoints)
    allocate(residuals(nray))
    read(ray_list_unit, *) residuals
    allocate(ray_stop(nray))
    read(ray_list_unit, *) ray_stop

!     read (95) npoints
!     npoints_max = maxval(npoints)
!     read (95) nv
!     allocate(v(nv))    
!     allocate(ray_stop(nray))
!     read (95) ray_stop

    allocate(s_vec(nray, npoints_max))
    allocate(v_vec(nray, npoints_max, nv ))

    ray_loop : do iray = 1, nray
        do ipoint = 1, npoints(iray)

            read (output_unit, *) s, v
!            read (94) s, v
            s_vec(iray, ipoint) = s
            v_vec(iray, ipoint, :) = v
        end do
    end do ray_loop

!     close(94)
!     close(95)

! Write data to an ASCII file for diagnostic    
    if (verbosity >= 3) then
       call message()
       call text_message('Ray data')    
       do iray = 1, nray
            call message('ray number = ', iray)
            do ipoint = 1, npoints(iray)
                write(message_unit, *) s_vec(iray, ipoint), v_vec(iray, ipoint, 1:nv)
            end do
            call text_message('stopped ',ray_stop(iray))
        end do
    end if

    call text_message('Finished initialize_post_process_rays')
    
    return
 end subroutine initialize_post_processing_m

!*************************************************************************     

  subroutine post_process

    use diagnostics_m, only : message_unit, message, text_message, verbosity
    use slab_processor_m, only : slab_processor
    use solovev_processor_m, only : solovev_processor
    use axisym_toroid_processor_m, only : axisym_toroid_processor
!	use solovev_magnetics_m, only : inner_bound, outer_bound, rmaj, kappa
    use axisym_toroid_eq_m, only : r_axis, z_axis, &
                          & box_rmin, box_rmax, box_zmin, box_zmax, &
                          & inner_bound, outer_bound, upper_bound, lower_bound

    implicit none

   write(*, *) 'r_axis = ', r_axis
   write(*, *) 'z_axis = ', z_axis
   write(*, *) 'inner_bound = ', inner_bound
   write(*, *) 'outer_bound = ', outer_bound
   write(*, *) 'upper_bound = ', upper_bound
   write(*, *) 'lower_bound = ', lower_bound
   
    select case (trim(processor))

       case ('slab')
          write(*,*) 'calling slab_processor'
          call slab_processor

       case ('solovev')
          write(*,*) 'calling solovev_processor'
          call solovev_processor

 
       case ('axisym_toroid')
          write(*,*) 'calling axisym_toroid_processor'
          call axisym_toroid_processor

      case default
          write(0,*) 'post_process_rays: unimplemented post_processor =', trim(processor)
          call text_message('post_process_rays: unimplemented post_processor', trim(processor))
          stop 1

       end select

    return
 end subroutine post_process

!*************************************************************************     

 
 end module post_processing_m

