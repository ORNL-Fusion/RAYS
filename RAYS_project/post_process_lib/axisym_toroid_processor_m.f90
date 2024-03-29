 module axisym_toroid_processor_m
! Post processing for axisym_toroid equilibrium

    use constants_m, only : rkind

    implicit none

    character(len=80) :: processor = ''
    integer, parameter :: graphics_descrip_unit = 96

! Number of k vectors to plot for each ray in graphics
    integer :: num_plot_k_vectors

! Scale plot k vectors to kmax on the ray, (True, False)
    character(len = 5) :: scale_k_vec = 'True'

! Set plot xlim -> [xmin, xmax] and ylim -> [ymin,ymax], (True, False)
    character(len = 5) :: set_XY_lim = 'True'

! Number of plasma boundary points to calculate
	integer :: n_boundary_points
! R,Z boundary points
    real(KIND=rkind), allocatable :: R_boundary(:), Z_boundary(:)

	logical :: calculate_dep_profiles, write_dep_profiles
    
    namelist /axisym_toroid_processor_list/ processor, num_plot_k_vectors, scale_k_vec,&
             & set_XY_lim, calculate_dep_profiles, write_dep_profiles

 contains

 subroutine initialize_axisym_toroid_processor(read_input)

    use diagnostics_m, only : message_unit, message, text_message, verbosity
	use deposition_profiles_m, only : initialize_deposition_profiles

    implicit none
    logical, intent(in) :: read_input
	 integer :: input_unit, get_unit_number ! External, free unit finder   
 
    if (read_input .eqv. .true.) then    
    ! Read and write input namelist
  		input_unit = get_unit_number()
        open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
        read(input_unit, axisym_toroid_processor_list)
        close(unit=input_unit)
        if (verbosity > 0 )write(message_unit, axisym_toroid_processor_list)
        call text_message('Finished initialize_axisym_toroid_processor ', processor, 1)
    end if

	if (calculate_dep_profiles .eqv. .true.) call initialize_deposition_profiles(read_input)
    
    return
 end subroutine initialize_axisym_toroid_processor

!*************************************************************************     

  subroutine axisym_toroid_processor

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message, verbosity
    use deposition_profiles_m, only : calculate_deposition_profiles, write_deposition_profiles

    implicit none
    
    call write_graphics_description_file

    if (calculate_dep_profiles .eqv. .true.) call calculate_deposition_profiles

    if (write_dep_profiles .eqv. .true.) call write_deposition_profiles
    
    if (verbosity > 0) call text_message('Finished axisym_toroid_processor work')

    return
 end subroutine axisym_toroid_processor


!*************************************************************************     

  subroutine find_plasma_boundary
    use constants_m, only : rkind
	use axisym_toroid_eq_m, only : magnetics_model, inner_bound, outer_bound, upper_bound, lower_bound
	use solovev_magnetics_m, only : rmaj, kappa
    use diagnostics_m, only : message_unit, message, text_message
    use eqdsk_utilities_m, only : NBOUND, RBOUND, ZBOUND
	
	implicit none
	
	integer :: i
	real(KIND=rkind) :: R, dR, Zsq
	 
    magnetics: select case (trim(magnetics_model))

    case ('solovev_magnetics')  ! N.B.  This is up-down symmetric
		n_boundary_points = 271 ! should be an odd number for this calculation
		allocate(R_boundary(n_boundary_points))
		allocate(Z_boundary(n_boundary_points))
		dR = 2.*(outer_bound - inner_bound)/(n_boundary_points)

		R_boundary(1) = inner_bound
		Z_boundary(1) = 0.
		R_boundary(n_boundary_points) = inner_bound
		Z_boundary(n_boundary_points) = 0.
		R_boundary((n_boundary_points-1)/2+1) = outer_bound
		Z_boundary(n_boundary_points) = 0.

		do i = 2, (n_boundary_points -1)/2
			R = inner_bound + i*dR
			Zsq = kappa*2/(4.*R**2)*(outer_bound**4 + 2.*(R**2 - outer_bound**2)*rmaj**2 -&
			      & R**4)
			R_boundary(i) = R
			Z_boundary(i) = sqrt(Zsq)
			R_boundary(n_boundary_points-(i-1)) = R_boundary(i)
			Z_boundary(n_boundary_points-(i-1)) = -Z_boundary(i)
		end do

    case ('eqdsk_magnetics_lin_interp')
        n_boundary_points = NBOUND ! not necessarily an odd number
		allocate(R_boundary(n_boundary_points))
		allocate(Z_boundary(n_boundary_points))
        R_boundary(:) = RBOUND(:)
        Z_boundary(:) = ZBOUND(:)

    case ('eqdsk_magnetics_spline_interp')
        n_boundary_points = NBOUND ! not necessarily an odd number
		allocate(R_boundary(n_boundary_points))
		allocate(Z_boundary(n_boundary_points))
        R_boundary(:) = RBOUND(:)
        Z_boundary(:) = ZBOUND(:)
        
        
	do i = 1, n_boundary_points
		write(*,*) 'i = ', i, '   R_boundary = ', R_boundary(i), '   Z_boundary', Z_boundary(i)
	end do
  
    case default
	  write(0,*) 'initialize_axisym_toroid_eq: unknown magnetics model =', magnetics_model
	  call text_message('initialize_axisym_toroid_eq: unknown magnetics model',&
	  & trim(magnetics_model),0)
	  stop 1
    end select magnetics

  end subroutine find_plasma_boundary
 
!*************************************************************************     

  subroutine write_graphics_description_file
  
   use diagnostics_m, only : run_description, run_label
   use axisym_toroid_eq_m, only : r_axis, z_axis, &
                          & box_rmin, box_rmax, box_zmin, box_zmax, &
                          & inner_bound, outer_bound, upper_bound, lower_bound
  
   open(unit = graphics_descrip_unit, file = 'graphics_description_axisym_toroid.dat')
  
   write(graphics_descrip_unit, *) 'run_description = ', run_description
   write(graphics_descrip_unit, *) 'run_label = ', run_label
  
   write(graphics_descrip_unit, *) 'r_axis = ', r_axis
   write(graphics_descrip_unit, *) 'z_axis = ', z_axis
   write(graphics_descrip_unit, *) 'inner_bound = ', inner_bound
   write(graphics_descrip_unit, *) 'outer_bound = ', outer_bound
   write(graphics_descrip_unit, *) 'upper_bound = ', upper_bound
   write(graphics_descrip_unit, *) 'lower_bound = ', lower_bound
     
   write(graphics_descrip_unit, *) 'box_rmin = ', box_rmin
   write(graphics_descrip_unit, *) 'box_rmax = ', box_rmax
   write(graphics_descrip_unit, *) 'box_zmin = ', box_zmin
   write(graphics_descrip_unit, *) 'box_zmax = ', box_zmax
  
   write(graphics_descrip_unit, *) 'num_plot_k_vectors = ', num_plot_k_vectors
   write(graphics_descrip_unit, *) 'scale_k_vec = ', trim(scale_k_vec)
   write(graphics_descrip_unit, *) 'set_XY_lim = ', trim(set_XY_lim)
   
   call find_plasma_boundary

   write(graphics_descrip_unit, *) ' '
   write(graphics_descrip_unit, *) 'R_boundary = ', R_boundary
   write(graphics_descrip_unit, *) ' '
   write(graphics_descrip_unit, *) 'Z_boundary = ', Z_boundary

   close(unit = graphics_descrip_unit)

  end subroutine write_graphics_description_file


 end module axisym_toroid_processor_m

