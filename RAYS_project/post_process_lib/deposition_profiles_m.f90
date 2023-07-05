 module deposition_profiles_m
 
! Contains data for various deposition profiles (such as power versus psi summed
! over species -> Ptotal(:)) and associated routines to calculate the profiles from ray 
! data in module ray_results_m.
!
! For now take all profile grids to have the same grid_min and grid_max, i.e. covering the
! whole plasma. Different profiles can have different number of points.
! grid_min, grid_max, and grid_name are picked up from specific equilibrium.  For slab
! they are [xmin, xmax], tor toroidal equilibria they are [0.0, 1.0]

    use constants_m, only : rkind

    implicit none
    
    logical :: write_results_list_directed = .true.

	type dep_profile
		character(len=20) :: profile_name,  grid_name
		integer :: n_bins
        real(KIND=rkind), allocatable :: grid(:), profile(:)
        real(KIND=rkind) :: Q_sum
    end type dep_profile

	type work_space
        real(KIND=rkind), allocatable :: work(:,:)
    end type work_space

! Data for all profiles
	character(len=20) :: grid_name
	real(KIND=rkind) :: grid_min, grid_max

! Data for total power deposition, Ptotal
	logical :: calc_Ptotal
	character(len=20) :: profile_name,  grid_name_Ptotal
	integer :: n_bins_Ptotal
	real(KIND=rkind) :: grid_min_Ptotal, grid_max_Ptotal
	type(work_space) :: Ptotal_work
	type(dep_profile) :: Ptotal_prof

! Data for other deposition profiles, someday

!***********************************

    abstract interface
		subroutine Q_evaluator(iray, ix, xQ, Q)   
			use constants_m, only : rkind
			use ray_results_m, only : ray_vec	  
			implicit none	
			integer, intent(in) :: iray, ix
			real(KIND=rkind), intent(out) :: xQ, Q
		end subroutine Q_evaluator
    end interface

    procedure(Q_evaluator), pointer :: evaluator => Null()

    namelist /deposition_profiles_list/ n_bins_Ptotal
    
!****************************************************************************

contains

!****************************************************************************

subroutine initialize_deposition_profiles(read_input)

	use constants_m, only : rkind
	use diagnostics_m, only : message_unit, text_message, run_label, date_v
	use ray_init_m, only : nray
	use equilibrium_m, only : equilib_model
	use slab_eq_m, only : xmin_slab => xmin, xmax_slab => xmax

	implicit none
	logical, intent(in) :: read_input
	integer :: input_unit, get_unit_number ! External, free unit finder   

	integer :: i
    real(KIND=rkind) :: delta
	

	write(*,*) 'initialize_deposition_profiles'     
	call text_message('initialize_deposition_profiles', 3)
	
	calc_Ptotal = .false.
	n_bins_Ptotal = 0

	if (read_input .eqv. .true.) then    
	! Read and write input namelist
  		input_unit = get_unit_number()
		open(unit=input_unit, file='post_process_rays.in',action='read', status='old', form='formatted')
		read(input_unit, deposition_profiles_list)
		close(unit=input_unit)
		write(message_unit, deposition_profiles_list)
	end if
		
    ! These hold for all profiles
    equilibria: select case (trim(equilib_model))
       case ('slab')
!         A 1-D slab equilibrium with stratification in x. Get grid range from slab xmin, xmax
          grid_min = xmin_slab
          grid_max = xmax_slab
          grid_name  = 'x'

       case ('solovev')
!         A simple analytic tokamak model. Grid range [0.0, 1.0]
          grid_min = 0.0_rkind
          grid_max = 1.0_rkind
          grid_name  = 'psi'

       case ('axisym_toroid')    
!         A generic axisymmetric toroidal plasma model.  Grid range [0.0, 1.0]
          grid_min = 0.0_rkind
          grid_max = 1.0_rkind
          grid_name  = 'psi'

       case default
          write(0,*) 'initialize_deposition_profiles: improper equilib_model =',&
                     & trim(equilib_model)
          call text_message('initialize_deposition_profiles: improper equilib_model',&
          & trim(equilib_model), 0)
          stop 1

    end select equilibria

		
	! initialize stuff for Ptotal
		if (n_bins_Ptotal > 1) then ! Otherwise there is no profile
			calc_Ptotal = .true.
			Ptotal_prof%profile_name = 'Ptotal'
			Ptotal_prof%grid_name = grid_name
			Ptotal_prof%n_bins = n_bins_Ptotal
			
			! Allocate type arrays for Ptotal
			allocate(Ptotal_prof%grid(n_bins_Ptotal+1), source = 0.0_rkind)
			allocate(Ptotal_prof%profile(n_bins_Ptotal), source = 0.0_rkind)
			allocate(Ptotal_work%work(n_bins_Ptotal, nray), source = 0.0_rkind)
			
			! Generate uniform grid for Ptotal	
			delta = (grid_max - grid_min)/real(n_bins_Ptotal)    
			do i = 1, n_bins_Ptotal+1
				Ptotal_prof%grid(i) = grid_min + delta*(i-1)
			end do
		end if

	! initialize stuff for other profiles just like above, someday

	   
	return
	end subroutine initialize_deposition_profiles

!****************************************************************************

    subroutine calculate_deposition_profiles
    
		use constants_m, only : rkind
		use ray_results_m
		use ray_init_m, only : nray
		use equilibrium_m, only : equilib_model
	  
		implicit none
	
		integer :: iray
		write(*,*) 'Calculating deposition profiles '
	
		ray_loop: do iray = 1, nray
		
			! Ptotal profile
			if (calc_Ptotal .eqv. .true.) then
				equilibria: select case (trim(equilib_model))
				   case ('slab')
						evaluator => Ptotal_slab_evaluator
				   case ('axisym_toroid')    
						evaluator => Ptotal_axisym_evaluator
				end select equilibria
				call bin_a_ray(iray, n_bins_Ptotal, grid_min, grid_max,&
				& Ptotal_work%work(:,iray))
			end if
	
		end do ray_loop
	
		! sum over rays to get total profile, sum over profile to get total quantity deposited
		if (calc_Ptotal .eqv. .true.) then
			Ptotal_prof%profile(:) = sum(Ptotal_work%work, 2)
			Ptotal_prof%Q_sum = sum(Ptotal_prof%profile)
    	end if

		! sum over rays to get other total profiles just like above, someday
    	
	return
    end subroutine calculate_deposition_profiles


!****************************************************************************

    subroutine bin_a_ray(iray, n_bins, xmin, xmax, binned_Q)
    
		use constants_m, only : rkind
		use ode_m, only : nstep_max
		use ray_results_m, only : npoints
		use ray_init_m, only : nray
   	    use bin_to_uniform_grid_m, only : bin_to_uniform_grid
	  
		implicit none
	
		integer, intent(in) :: iray, n_bins
		real(KIND=rkind), intent(in) :: xmin, xmax
		real(KIND=rkind), intent(out) :: binned_Q(n_bins)
		
		real(KIND=rkind):: Q(nstep_max+1)  ! quantity at ray point i
		real(KIND=rkind) :: xQ(nstep_max+1)  ! grid value at ray point i

		integer :: ip, ierr
		
		! Evaluate grid value and profile quantity at each ray point
		do ip = 1, npoints(iray)
			call evaluator(iray, ip, xQ(ip), Q(ip))
		end do
		
    	call bin_to_uniform_grid(Q(1:npoints(iray)), xQ(1:npoints(iray)), xmin, xmax,&
    	        &  binned_Q, ierr)
    	
	return
    end subroutine bin_a_ray

!****************************************************************************

    subroutine write_deposition_profiles
    
		use diagnostics_m, only : run_label
		use ode_m, only : nstep_max
		use ray_results_m, only : npoints
		use ray_init_m, only : nray
	  
		implicit none
		
		integer :: dep_profile_unit, get_unit_number
		
		!  File name for message output
		character(len=80) :: out_file
	
		! Open file to put deposition data in
		dep_profile_unit = get_unit_number()
		out_file = 'deposition_profiles.'//trim(run_label)
		open(unit=dep_profile_unit, file=out_file, action='write', status='replace',&
		     & form='formatted') 

	! Stuff for Ptotal
		
		if (calc_Ptotal .eqv. .true.) then
			write(dep_profile_unit, *) Ptotal_prof%grid_name
			write(dep_profile_unit, *) Ptotal_prof%grid
			write(dep_profile_unit, *) Ptotal_prof%profile_name
			write(dep_profile_unit, *) Ptotal_prof%profile
			write(dep_profile_unit, *) 'Ptotal_total_deposition'
			write(dep_profile_unit, *) Ptotal_prof%Q_sum
		end if

	! Stuff for other profiles, someday
		
    close(unit = dep_profile_unit)

	return
    end subroutine write_deposition_profiles

!****************************************************************************
! Evaluators
!****************************************************************************

    subroutine Ptotal_slab_evaluator(iray, ix, xQ, Q)
    
		use constants_m, only : rkind
		use ray_results_m, only : ray_vec
	  
		implicit none
	
		integer, intent(in) :: iray, ix
		real(KIND=rkind), intent(out) :: xQ, Q
		
		xQ = ray_vec(1, ix, iray)
		Q = ray_vec(8, ix, iray)

 	return
    end subroutine Ptotal_slab_evaluator

!********************************************************************

    subroutine Ptotal_axisym_evaluator(iray, ix, xQ, Q)
    
		use constants_m, only : rkind
		use ray_results_m, only : ray_vec
    	use axisym_toroid_eq_m, only: axisym_toroid_psi
	  
		implicit none
	
		integer, intent(in) :: iray, ix
		real(KIND=rkind), intent(out) :: xQ, Q

    	real(KIND=rkind) :: rvec(3) 
		real(KIND=rkind) :: psi, gradpsi(3), psiN, gradpsiN(3)
		rvec = ray_vec(1:3, ix, iray)
		call axisym_toroid_psi(rvec, psi, gradpsi, psiN, gradpsiN)
		
		xQ = psiN
		Q = ray_vec(8, ix, iray)

 	return
    end subroutine Ptotal_axisym_evaluator

!********************************************************************
! Deallocate
!****************************************************************************

    subroutine deallocate_deposition_profiles_m
        
        if (calc_Ptotal .eqv. .true.) then      	
			deallocate(Ptotal_prof%grid)
			deallocate(Ptotal_prof%profile)
			deallocate(Ptotal_work%work)
		end if
        
        return
    end subroutine deallocate_deposition_profiles_m
     
 end module deposition_profiles_m
