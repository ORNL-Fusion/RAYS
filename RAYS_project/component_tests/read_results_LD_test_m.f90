module read_results_LD_test_m

    use constants_m, only : rkind    
    
    implicit none
    
 !  File name for input
    character(len=80) :: in_filename

    namelist /read_results_LD_test_list/ in_filename

 !*************************************************************************     

 contains

 subroutine init_read_results_LD_test(read_input)

    use diagnostics_m, only : message_unit, message, text_message, run_label

    implicit none
    
    logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder   

    if (read_input .eqv. .true.) then    
    ! Read and write input namelist
   		input_unit = get_unit_number()
        open(unit=input_unit, file='component_test_rays.in',action='read', status='old', form='formatted')
        read(input_unit, read_results_LD_test_list)
        close(unit=input_unit)
        write(message_unit, read_results_LD_test_list)
    end if


    call text_message('Finished read_results_LD_test_list ',1)
            
    return
 end subroutine init_read_results_LD_test
 
 !*************************************************************************     

  subroutine read_results_LD_test
   
    use diagnostics_m, only : text_message
    use ray_results_m, only : read_results_LD
    
    implicit none
    
    call read_results_LD(in_filename)

    call text_message('Finished read_results_LD_test work')

    return
 end subroutine read_results_LD_test


!*************************************************************************     

end module read_results_LD_test_m