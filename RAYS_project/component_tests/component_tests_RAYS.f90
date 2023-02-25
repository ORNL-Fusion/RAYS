 program component_test_RAYS

    use diagnostics_m, only : message_unit, message, text_message, verbosity
    use constants_m, only : input_unit, output_unit, ray_list_unit
    
    use damping_slab_test_m, only : init_damping_slab_test, damping_slab_test
    use read_results_LD_test_m, only : init_read_results_LD_test, read_results_LD_test

    implicit none 
    logical :: read_input = .true.
    character(len=60) :: test_name

    namelist /component_test_list/ test_name

!   Read input files and initialize variables needed from RAYS_lib
    call initialize_components_in_RAYS_lib(read_input)

    if (read_input .eqv. .true.) then    
	! Read and write input namelist
        open(unit=input_unit, file='component_test_rays.in',action='read', status='old',&
                                  & form='formatted')
        read(input_unit, component_test_list)
        close(unit=input_unit)
        write(message_unit, component_test_list)
	end if

    select case (trim(test_name))

       case ('damping_slab')
          call init_damping_slab_test(read_input)
          call damping_slab_test

       case ('read_results_LD')
          call init_read_results_LD_test(read_input)
          call read_results_LD_test

       case default
          write(*,*) 'Unimplemented component tester =', trim(test_name)
          call text_message('Unimplemented component testor', trim(test_name))
          stop 1

       end select

 end program component_test_RAYS
