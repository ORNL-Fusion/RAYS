 module openmp_m
! Initializes the open mp setup, i.e. the number of threads.  There is no namelist data.
! In this case the 'read_input' variable tells the init subroutine whether to set num_threads
! to the minimum of (num_procs,nray) or to use the module variable num_threads which has been
! set directly in a host program.  This make it possible to do scans of thread number scaling

    use ray_init_m, only : nray
    use omp_lib
    use diagnostics_m, only : message

    implicit none

    integer :: num_threads, num_procs

 namelist /openmp_list/ num_threads

 contains
    subroutine initialize_openmp_m(read_input)

       implicit none
       logical, intent(in) :: read_input
 	   integer :: input_unit, get_unit_number ! External, free unit finder

 write(*,*) 'num_threads = ',num_threads
	   num_procs = OMP_get_num_procs()
       if (read_input .eqv. .true.) then  !See if there is num_threads != 0 in namelist
            num_threads = 0
 			! Read num_threads from namelist.
			input_unit = get_unit_number()
			open(unit=input_unit, file='rays.in',action='read', status='old', form='formatted')
			read(input_unit, openmp_list, end=1, err=1)
			close(unit=input_unit)
1       	continue ! error return, no namelist
            if (num_threads == 0) then ! Calculate num_threads
!			   num_threads = min(num_procs,nray)
			   num_threads = min(8,nray)
!			   call omp_set_num_threads(num_threads)
			end if
       end if
       ! read_input = false, use num_threads value set from outside

    write(12,*) 'initialize_openmp_m: num_procs = ',num_procs,'  num threads = ', num_threads
        call message ('initialize_openmp_m: num_procs', num_procs, 0)
        call message ('initialize_openmp_m: num threads', num_threads, 0)

    end subroutine initialize_openmp_m

!********************************************************************

    subroutine deallocate_openmp_m
        return ! nothing to deallocate
    end subroutine deallocate_openmp_m

 end module openmp_m
