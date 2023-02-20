module damping_slab_test_m

    use constants_m, only : rkind    
    
    implicit none

	abstract interface
		 subroutine damping_subroutine(eq, v, vg, ksi, ki)
			use constants_m, only : clight, rkind
			use rf_m, only : omgrf, k0
			use species_m, only : nspec, qs, ms
			use ode_m, only : nv
			use equilibrium_m, only : eq_point
			use zfunctions_m, only : zfun0

			implicit none
	
			type(eq_point), intent(in) :: eq
			real(KIND=rkind), intent(in) :: v(6)
			real(KIND=rkind), intent(in) :: vg(3)
            real(KIND=rkind), intent(out) :: ksi(0:nspec), ki
		end subroutine damping_subroutine
	end interface
     
    procedure (damping_subroutine), pointer :: damping => NULL()
    external damp_fund_ECH

    integer :: test_out_unit, formatted_out_unit
    
! Data for x scan and nz scan
    integer :: n_xpoints ! Number of x points in scan
    real(KIND=rkind) :: xmin, xmax ! min and max of x scan
    integer :: n_nz_points ! Number of points in nz scan
    real(KIND=rkind) :: nz_min, nz_max ! min and max of nz scan
    
! Switch to select specific damping model
    character(len=60) :: damping_model = ''

    namelist /damping_slab_test_list/ damping_model,&
            & n_xpoints, xmin, xmax, n_nz_points, nz_min, nz_max

 !*************************************************************************     

 contains

 subroutine init_damping_slab_test(read_input)

    use diagnostics_m, only : message_unit, message, text_message, run_label
    use constants_m, only : input_unit

    implicit none
    
    logical, intent(in) :: read_input

    if (read_input .eqv. .true.) then    
    ! Read and write input namelist
        open(unit=input_unit, file='component_test_rays.in',action='read', status='old', form='formatted')
        read(input_unit, damping_slab_test_list)
        close(unit=input_unit)
        write(message_unit, damping_slab_test_list)
    end if

    select case (trim(damping_model))

       case ('damp_fund_ECH')
          damping => damp_fund_ECH

       case default
          write(*,*) 'Unimplemented damping model =', trim(damping_model)
          call text_message('Unimplemented damping model', trim(damping_model))
          stop 1

       end select

    call text_message('Finished initialize_damping_slab_test ', damping_model)
            
    return
 end subroutine init_damping_slab_test
 
 !*************************************************************************     

  subroutine damping_slab_test

    use constants_m, only : rkind
    use diagnostics_m, only : message_unit, message, text_message, run_label
    use equilibrium_m, only : equilibrium, eq_point
    use rf_m, only : ray_dispersion_model, wave_mode, k0
    use species_m, only : nspec
    use dispersion_solvers_m, only : solve_n1_vs_n2_n3
    use trapezoid_quad_m, only : trapezoid_quad
   
    implicit none

!   Derived type containing equilibrium data for a spatial point in the plasma
    type(eq_point) :: eq

! Args for dispersion subroutine 
   complex(KIND=rkind) :: nx
   real(KIND=rkind) :: ny = 0.
   real(KIND=rkind) :: nz
   integer :: k_sign = 1 ! Choose positive sign for k_x
 

! Inputs for damping subroutine 
    real(KIND=rkind) :: v(6)
    real(KIND=rkind) :: vg(3)
    real(KIND=rkind) :: ksi(0:nspec), ki

    integer :: ix, inz
    real(KIND=rkind) :: x, dx, dnz
    real(KIND=rkind) :: rvec(3), kvec(3)

! Results
    real(KIND=rkind) :: x_vec(n_xpoints), ki_vec(n_xpoints), ni_vec(n_xpoints)
    real(KIND=rkind) :: alpha ! Optical depth

    integer :: get_unit_number
    external get_unit_number

 !*************************************************************************     

    ! Open file to put write(*,*) data in
	test_out_unit = get_unit_number()
	open (unit = test_out_unit, file = 'slab_damping.out', status = 'unknown')

    ! Open file to put x, k_Im in
	formatted_out_unit = get_unit_number()
	open(unit=formatted_out_unit, file='x_k_Im_'//trim(run_label),action='write',&
			 &  status='unknown', form='formatted')
    
    dx = (xmax - xmin)/(n_xpoints-1)
    if (n_nz_points > 1) then
    	dnz = (nz_max - nz_min)/(n_nz_points-1)
    else
    	dnz = 0.
    end if
    ny = 0.
    
    ! Set vg = (1,0,0) to calculate ki along x direction
    vg = (/ 1.0_rkind, 0.0_rkind,  0.0_rkind /)

! get x vector
	do ix = 0, n_xpoints-1
		x = xmin + ix*dx
		x_vec(ix+1) = x
	end do
	write(test_out_unit, *) 'x_vec' 
	write(test_out_unit, *) x_vec

    nz_loop: do inz = 0, n_nz_points-1
		nz = nz_min + inz*dx
		
		x_loop: do ix = 1, n_xpoints
			rvec = (/ x_vec(ix), 0.0_rkind,  0.0_rkind /)	
			call equilibrium(rvec, eq )
            call solve_n1_vs_n2_n3(eq, ray_dispersion_model, wave_mode, k_sign,&
                                     & ny, nz, nx)
			kvec = k0*(/ real(nx), ny, nz /)
			v(1:3) = rvec
			v(4:6) = kvec		
 			call damping(eq, v, vg, ksi, ki)
			ki_vec(ix) = ki
			ni_vec(ix) = ki/k0
			write(formatted_out_unit,*) x_vec(ix), '  ', ni_vec(ix)
		end do x_loop

! get optical depth
        call trapezoid_quad(x_vec, ki_vec, alpha)    
		
		write(test_out_unit, *) 'nz'
		write(test_out_unit, *) nz
		write(test_out_unit, *) 'alpha'
		write(test_out_unit, *) alpha
		write(test_out_unit, *) 'ki_vec'
		write(test_out_unit, *) ki_vec
		write(test_out_unit, *) 'ni_vec'
		write(test_out_unit, *) ni_vec

	end do nz_loop
       
	close (unit = test_out_unit)
	close (unit = formatted_out_unit)
    call text_message('Finished damping_slab_test work')

    return
 end subroutine damping_slab_test


!*************************************************************************     

end module damping_slab_test_m