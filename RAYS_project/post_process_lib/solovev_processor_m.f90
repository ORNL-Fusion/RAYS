 module solovev_processor_m
! Post processing for solovev equilibrium

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message, messages_to_stdout, verbosity

    implicit none

!_________________________________________________________________________________________
! Namelist data for /solovev_processor_list/
!_________________________________________________________________________________________

    character(len=80) :: processor = ''

! Number of k vectors to plot for each ray in graphics
    integer :: num_plot_k_vectors

! Scale plot k vectors to kmax on the ray, (True, False)
    character(len = 5) :: scale_k_vec = 'True'

! Set plot xlim -> [xmin, xmax] and ylim -> [ymin,ymax], (True, False)
    character(len = 5) :: set_XY_lim = 'True'

    namelist /solovev_processor_list/ processor, num_plot_k_vectors, scale_k_vec, set_XY_lim

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

 subroutine initialize_solovev_processor(read_input)

    use diagnostics_m, only : message_unit, message, text_message

    implicit none
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder

    if (read_input .eqv. .true.) then
    ! Read and write input namelist
   		input_unit = get_unit_number()
        open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
        read(input_unit, solovev_processor_list)
        close(unit=input_unit)
        if (verbosity > 0) write(message_unit, solovev_processor_list)
        call text_message('Finished initialize_solovev_processor ', processor, 1)
    end if

    return
 end subroutine initialize_solovev_processor

!*************************************************************************

  subroutine solovev_processor

    use constants_m, only : rkind
    use diagnostics_m, only : text_message

    implicit none

    call write_graphics_description_file

    call text_message('Finished solovev_processor work', 1)

    return
 end subroutine solovev_processor


!*************************************************************************


!*************************************************************************

  subroutine write_graphics_description_file

   use diagnostics_m, only : run_description, run_label
   use solovev_eq_m, only : rmaj, kappa, bphi0, iota0, &
                          & box_rmin, box_rmax, box_zmin, box_zmax, &
                          & inner_bound, outer_bound, vert_bound, r_Zmax, &
                          &  psiB

   implicit none
   integer :: graphics_descrip_unit, get_unit_number

   graphics_descrip_unit = get_unit_number()
   open(unit = graphics_descrip_unit, file = 'graphics_description_solovev.dat')

   write(graphics_descrip_unit, *) 'run_description = ', run_description
   write(graphics_descrip_unit, *) 'run_label = ', run_label

   write(graphics_descrip_unit, *) 'rmaj = ', rmaj
   write(graphics_descrip_unit, *) 'kappa = ', kappa
   write(graphics_descrip_unit, *) 'bphi0 = ', bphi0
   write(graphics_descrip_unit, *) 'iota0 = ', iota0
   write(graphics_descrip_unit, *) 'psiB = ', psiB
   write(graphics_descrip_unit, *) 'inner_bound = ', inner_bound
   write(graphics_descrip_unit, *) 'outer_bound = ', outer_bound
   write(graphics_descrip_unit, *) 'vert_bound = ', vert_bound


   write(graphics_descrip_unit, *) 'box_rmin = ', box_rmin
   write(graphics_descrip_unit, *) 'box_rmax = ', box_rmax
   write(graphics_descrip_unit, *) 'box_zmin = ', box_zmin
   write(graphics_descrip_unit, *) 'box_zmax = ', box_zmax

   write(graphics_descrip_unit, *) 'num_plot_k_vectors = ', num_plot_k_vectors
   write(graphics_descrip_unit, *) 'scale_k_vec = ', trim(scale_k_vec)
   write(graphics_descrip_unit, *) 'set_XY_lim = ', trim(set_XY_lim)

   close(unit = graphics_descrip_unit)

  end subroutine write_graphics_description_file


 end module solovev_processor_m

