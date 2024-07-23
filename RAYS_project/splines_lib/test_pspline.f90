    program test_psplines

    use quick_cube_splines_m, only : cube_spline_function_1D, cube_spline_function_2D

    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind), parameter :: zero = 0.0_rkind, one = 1.0_rkind
    real(KIND=rkind), parameter :: pi = atan(zero, -one), sqrt_pi = sqrt(pi)

! Setup

    logical :: do_2D = .false., do_1D = .false., do_C = .false., do_Z = .false.
    logical :: do_1D_module = .false., do_2D_module = .true.
    logical :: write_details = .true.

! Stuff for 2D splines

    integer, parameter :: nx = 21, ny = 21
    real(KIND=rkind) ::  x_grid_min = -10., x_grid_max = 10., y_grid_min = -10., y_grid_max = 10.

    integer, parameter :: nx_get = 21, ny_get = 21
    real(KIND=rkind) ::  x_get_min = -10., x_get_max = 10., y_get_min = -10., y_get_max = 10.
    integer :: nr = 5, nTheta = 10
    real(KIND=rkind) :: r, theta, rget

    integer :: i, j
    real(KIND=rkind) ::  f, fx, fy, fxx, fxy, fyy

    real :: start_time_coef, end_time_coef, start_time_fn, end_time_fn, &
                            & start_time_spline, end_time_spline

! bcspline arguments

  integer, parameter :: inx = nx, inth = ny, inf3 = nx, nwk = 4*inx*inth +5*max(inx,inth)

  integer ibcxmin,ibcxmax,ibcthmin,ibcthmax,ilinx, ilinth,ier

  !============
  real(KIND=rkind) :: x(inx),th(inth),fspl_2D(4,4,inf3,inth),wk(nwk)
  real(KIND=rkind) :: bcxmin(inth),bcxmax(inth)      ! (inth) if used
  real(KIND=rkind) :: bcthmin(inx),bcthmax(inx)    ! (inx) if used

! bcspeval arguments, if not already declared above

    real(KIND=rkind) :: xget, yget, fget, y(ny), fval(6)
    integer :: iselect(6)

! Stuff for 1D splines

    integer, parameter :: nx1D =1001
    real(KIND=rkind) ::  x1D_grid_min = -5., x1D_grid_max = 5., x1D_grid(nx1D)

    integer, parameter :: nx1D_get = 101
    real(KIND=rkind) ::  x1D_get_min = -5., x1D_get_max = 5., x1D(nx1D_get)
    real(KIND=rkind) ::  f1D(nx1D_get), fx1D(nx1D_get), fxx1D(nx1D_get)

! bcspeval arguments, if not already declared above

    real(KIND=rkind) ::  fspl1D(4, nx1D), fval1D(3)
    integer :: iselect1D(3)
    integer, parameter :: iwk = nx1D
    real(KIND=rkind) ::  wk1D(iwk)

! Stuff for C splines

    integer, parameter :: nxC = 2001
    real(KIND=rkind) ::  xC_grid_min = -5., xC_grid_max = 5., xC_grid(nxC)

    integer, parameter :: nxC_get = 1000
    real(KIND=rkind) ::  xC_get_min = -5., xC_get_max = 5., xC(nxC_get)
    complex(KIND=rkind) ::  fC, fxC, fxxC
    real(KIND=rkind) ::  fRe(nxC_get), fxRe(nxC_get), fxxRe(nxC_get)
    real(KIND=rkind) ::  fIm(nxC_get), fxIm(nxC_get), fxxIm(nxC_get)

! bcspeval arguments, if not already declared above

    real(KIND=rkind) ::  fsplRe(4, nxC)=0., fvalRe(3)=0., fsplIm(4, nxC)=0., fvalIm(3)=0.
    integer :: iselectC(3)
    integer, parameter :: iwkC = nxC
    real(KIND=rkind) ::  wkC(iwk)

! Stuff for cube_spline_function_1D

    type(cube_spline_function_1D) :: spl_func
	character (len = 80) :: func_name
    real(KIND=rkind) :: fsp, fxsp, fxxsp

! Stuff for cube_spline_function_2D

    type(cube_spline_function_2D) :: spl_func_2D
    real(KIND=rkind) :: f2sp, f2xsp, f2ysp, f2xxsp, f2xysp, f2yysp

! Statistics

    real(KIND=rkind) ::  random, rel_err, abs_err, s_abs, s_rel

!***************************************************************************************
! 2D splines

  D2: if (do_2D .eqv. .true.) then

    write (*,*)
    write (*,*) '2D results'

    call cpu_time(start_time_coef)

! Load grids and function array
    do j = 1, ny
        y(j) = y_grid_min + (j-1)*(y_grid_max - y_grid_min)/(ny-1)
    end do
    do i = 1, nx
        x(i) = x_grid_min + (i-1)*(x_grid_max - x_grid_min)/(nx-1)
    end do

    do j = 1, ny
        do i = 1, nx
            call spline_test_fn(x(i), y(j), f, fx, fy, fxx, fxy, fyy)
            fspl_2D(1,1,i,j) = f
         end do
    end do

! Get spline coefficients
    th = y
    ibcxmin = 0
    bcxmin = 0.
    ibcxmax = 0
    bcxmax = 0.
    ibcthmin = 0
    bcthmin = 0.
    ibcthmax = 0
    bcthmax = 0.
    ilinx = 0
    ilinth = 0

    call bcspline(x,inx,th,inth,fspl_2D,inf3, &
         ibcxmin,bcxmin,ibcxmax,bcxmax, &
         ibcthmin,bcthmin,ibcthmax,bcthmax, &
         wk,nwk,ilinx,ilinth,ier)
         if (ier .ne. 0) write (*,*) 'bcspline: ier = ', ier

    call cpu_time(end_time_coef)

!    write (*,*) f(1,1,:,:)

! Evaluate splined function on data

    call cpu_time(start_time_spline)

    s_abs = 0._rkind; s_rel = 0._rkind
    iselect = 1 ! Selector for outputs

    do i = 1, nr
        do j = 1, nTheta
            r = real(i)
            theta = (j-1)*pi/nTheta

            xget = r*cos(theta)
            yget = r*sin(theta)
            rget = sqrt(xget**2 + yget**2)

            call bcspeval(xget,yget,iselect,fval, &
                  x,inx,th,inth,ilinx,ilinth,fspl_2D,inf3,ier)
            if (ier .ne. 0) write (*,*) 'bcspeval: ier = ', ier

            call spline_test_fn(xget, yget, f, fx, fy, fxx, fxy, fyy)

			abs_err = f - fval(1)
			rel_err = abs_err/fval(1)
			s_abs = s_abs + abs(abs_err)
			s_rel = s_rel + abs(rel_err)

            if (write_details) then
				write (*,*)
				write (*,*) 'xget = ', xget, '     yget = ', yget
				write (*,*) 'fget = ', f, '  fval(1) = ', fval(1)
				write (*,*) 'fx = ', fx, '  fval(2) = ', fval(2)
				write (*,*) 'fy = ', fy, '  fval(3) = ', fval(3)
				write (*,*) 'fxx = ', fxx, '  fval(4) = ', fval(4),'fxy = ', fxy, &
				&  '  fval(6) = ', fval(6), '  fyy = ', fyy, '  fval(5) = ', fval(5)
            end if

        end do
    end do

    call cpu_time(end_time_spline)

 		write(*,*) ' '
		write(*,*) 'sqrt(x**2 + y**2)'
		write(*,*) 'x grid points = ', nx, ' y grid points =  = ', ny
		write(*,*) 'CPU time coeff eval = ', end_time_coef - start_time_coef
		write(*,*) 'CPU time spline eval = ', end_time_spline - start_time_spline
		write(*,*) 'Average absolute error = ', s_abs/(nr*nTheta)
		write(*,*) 'Average relative error = ', s_rel/(nr*nTheta)

  end if D2

!***************************************************************************************
! 1D splines

  D1: if (do_1D .eqv. .true.) then

    write (*,*)
    write (*,*) '1D results'

    call cpu_time(start_time_coef)

! 1D Grid
    do i = 1, nx1D
        x1D_grid(i) = x1D_grid_min + (i-1)*(x1D_grid_max - x1D_grid_min)/(nx1D-1)
    end do

! Grid is generated above, evaluate function array
    do i = 1, nx1D
        call spline_test_fn_1D(x1D_grid(i), f, fx, fxx)
        fspl1D(1,i) = f
    end do

! Get spline coefficients
    ibcxmin = 0
    bcxmin = 0.
    ibcxmax = 0
    bcxmax = 0.
    ibcthmin = 0
    bcthmin = 0.
    ibcthmax = 0
    bcthmax = 0.
    ilinx = 0

    call cspline(x1D_grid,nx1D,fspl1D,ibcxmin,bcxmin,ibcxmax,bcxmax,wk1D,iwk,ilinx,ier)
    write (*,*) 'ilinx = ', ilinx

    call cpu_time(end_time_coef)

! Evaluate splined function on data

    iselect1D = (/ 1, 1, 1 /)

    do i = 1, nx1D_get
        call random_number(random)
        x1D(i) = x1D_get_min + (x1D_get_max - x1D_get_min)*random
        x1D(i) = x1D_get_min + (i-1)*(x1D_get_max - x1D_get_min)/(nx1D_get-1)/sqrt(2.)
    end do

! Evaluate the function analytically for timing
    call cpu_time(start_time_fn)
    do i = 1, nx1D_get
        call spline_test_fn_1D(x1D(i), f1D(i), fx1D(i), fxx1D(i))
    end do
    call cpu_time(end_time_fn)

! Evaluate using splines
    call cpu_time(start_time_spline)

    s_abs = 0._rkind; s_rel = 0._rkind
    do i = 1, nx1D_get
        call cspeval(x1D(i),iselect1D,fval1D,x1D_grid,nx1D,ilinx,fspl1D,ier)
        if (ier .ne. 0) write (*,*) 'cspeval: ier = ', ier

        abs_err = fval1D(1) - f1D(i)
        rel_err = abs_err/f1D(i)
        s_abs = s_abs + abs(abs_err)
        s_rel = s_rel + abs(rel_err)

   if (write_details) then
        write (*,*) ' '
        write (*,*) 'x = ', x1D(i), '  f1D(i) = ', f1D(i), '  f spline = ', fval1D(1)
        write (*,*) 'abs_err = ', abs_err, '  rel_err = ', rel_err
        write (*,*) 'fx = ', fx, '  fx spline = ', fval1D(2),'  fxx = ', fxx,&
               &  '  fxx spline = ', fval1D(3)
    end if

    end do
    call cpu_time(end_time_spline)

 		write(*,*) ' '
		write(*,*) '1 + SIN(k x)  grid points = ', nx1D, '  sample points = ', nx1D_get
		write(*,*) 'CPU time coeff eval = ', end_time_coef - start_time_coef
		write(*,*) 'CPU time function eval = ', end_time_fn - start_time_fn
		write(*,*) 'CPU time spline eval = ', end_time_spline - start_time_spline
		write(*,*) 'Average absolute error = ', s_abs/nx1D_get
		write(*,*) 'Average relative error = ', s_rel/nx1D_get

  end if D1


!***************************************************************************************
! Complex splines

  C: if (do_C .eqv. .true.) then

    write (*,*)
    write (*,*) 'Complex results'

    call cpu_time(start_time_coef)

! C Grid
    do i = 1, nxC
        xC_grid(i) = xC_grid_min + (i-1)*(xC_grid_max - xC_grid_min)/(nxC-1)
    end do

! Grid is generated above, evaluate function array
    do i = 1, nxC
        call spline_test_fn_C(xC_grid(i), fC, fxC, fxxC)
        fsplRe(1,i) = real(fC, rkind); fsplIm(1,i) = aimag(fC)
    end do

! Get spline coefficients
    ibcxmin = 0
    bcxmin = 0.
    ibcxmax = 0
    bcxmax = 0.
    ibcthmin = 0
    bcthmin = 0.
    ibcthmax = 0
    bcthmax = 0.
    ilinx = 0

    call cspline(xC_grid,nxC,fsplRe,ibcxmin,bcxmin,ibcxmax,bcxmax,wkC,iwkC,ilinx,ier)
    call cspline(xC_grid,nxC,fsplIm,ibcxmin,bcxmin,ibcxmax,bcxmax,wkC,iwkC,ilinx,ier)
    write (*,*) 'ilinx = ', ilinx

    call cpu_time(end_time_coef)

! Evaluate splined function on data

    iselectC = (/ 1, 1, 1 /)

    do i = 1, nxC_get
        call random_number(random)
        xC(i) = xC_get_min + (xC_get_max - xC_get_min)*random
        xC(i) = xC_get_min + (i-1)*(xC_get_max - xC_get_min)/(nxC_get-1)/sqrt(2.)
    end do

! Evaluate the function analytically for timing
    call cpu_time(start_time_fn)
    do i = 1, nxC_get
        call spline_test_fn_C(xC(i), fC, fxC, fxxC)
        fRe(i) = real(fC,rkind); fIm(i) = aimag(fC)
        fxRe(i) = real(fxC,rkind); fxIm(i) = aimag(fxC)
        fxxRe(i) = real(fxxc,rkind); fxxIm(i) = aimag(fxxC)
    end do
    call cpu_time(end_time_fn)

! Evaluate using splines
    call cpu_time(start_time_spline)

    s_abs = 0.; s_rel = 0.
    do i = 1, nxC_get
        call cspeval(xC(i),iselectC,fvalRe,xC_grid,nxC,ilinx,fsplRe,ier)
        call cspeval(xC(i),iselectC,fvalIm,xC_grid,nxC,ilinx,fsplIm,ier)
        if (ier .ne. 0) write (*,*) 'cspeval: ier = ', ier

        abs_err = abs(cmplx(fvalRe(1) - fRe(i), fvalIm(1) - fIm(i), kind=rkind))
        rel_err = abs_err/abs(cmplx(fRe(i),fIm(i), kind=rkind))
        s_abs = s_abs + abs_err
        s_rel = s_rel + rel_err

       if (write_details) then
			write (*,*) ' '
			write (*,*) 'x = ', xC(i)
			write (*,*) '  fRe(i) = ', fRe(i), '  Re(f spline) = ', fvalRe(1)
			write (*,*) '  fIm(i) = ', fIm(i), '  Im(f spline) = ', fvalIm(1)
			write (*,*) '  fxRe(i) = ', fxRe(i), '  Re(fx spline) = ', fvalRe(2)
			write (*,*) '  fxIm(i) = ', fxIm(i), '  Im(fx spline) = ', fvalIm(2)
			write (*,*) '  fxxRe(i) = ', fxxRe(i), '  Re(fxx spline) = ', fvalRe(3)
			write (*,*) '  fxxIm(i) = ', fxxIm(i), '  Im(fxx spline) = ', fvalIm(3)
			write (*,*) 'abs_err = ', abs_err, '  rel_err = ', rel_err
			write (*,*) 'diff Re = ', fvalRe(1) - fRe(i), '  Diff Im = ', fvalIm(1) - fIm(i)
        end if

    end do
    call cpu_time(end_time_spline)

    write(*,*) ' '
    write(*,*) 'Complex function  grid points = ', nxC, '  sample points = ', nxC_get
    write(*,*) 'CPU time coeff eval = ', end_time_coef - start_time_coef
    write(*,*) 'CPU time function eval = ', end_time_fn - start_time_fn
    write(*,*) 'CPU time spline eval = ', end_time_spline - start_time_spline
    write(*,*) 'Average absolute error = ', s_abs/nxC_get
    write(*,*) 'Average relative error = ', s_rel/nxC_get

  end if C
!***************************************************************************************
! Zfunction splines

  Z: if (do_Z .eqv. .true.) then

    write (*,*)
    write (*,*) 'Splined Zfunction results'

    call cpu_time(start_time_coef)

! C Grid
    do i = 1, nxC
        xC_grid(i) = xC_grid_min + (i-1)*(xC_grid_max - xC_grid_min)/(nxC-1)
     end do

! Grid is generated above, evaluate function array
    do i = 1, nxC
        call spline_test_fn_Zf(xC_grid(i), fC, fxC, fxxC)
        fsplRe(1,i) = real(fC, kind=rkind); fsplIm(1,i) = aimag(fC)
    end do

! Get spline coefficients
    ibcxmin = 0
    bcxmin = 0.
    ibcxmax = 0
    bcxmax = 0.
    ibcthmin = 0
    bcthmin = 0.
    ibcthmax = 0
    bcthmax = 0.
    ilinx = 0

    call cspline(xC_grid,nxC,fsplRe,ibcxmin,bcxmin,ibcxmax,bcxmax,wkC,iwkC,ilinx,ier)
    call cspline(xC_grid,nxC,fsplIm,ibcxmin,bcxmin,ibcxmax,bcxmax,wkC,iwkC,ilinx,ier)
    write (*,*) 'ilinx = ', ilinx

    call cpu_time(end_time_coef)

! Evaluate splined function on data

    iselectC = (/ 1, 0, 0 /)

    do i = 1, nxC_get
        call random_number(random)
        xC(i) = xC_get_min + (xC_get_max - xC_get_min)*random
!        xC(i) = xC_get_min + (i-1)*(xC_get_max - xC_get_min)/(nxC_get-1)/sqrt(2.)
!        xC(i) = xC_get_min + (i-1)*(xC_get_max - xC_get_min)/(nxC_get-1)+ 10.0e-7_rkind
!        xC(i) = xC_grid(i) + 10.0e-4
    end do

! Evaluate the function analytically for timing
    call cpu_time(start_time_fn)
    do i = 1, nxC_get
        call spline_test_fn_Zf(xC(i), fC, fxC, fxxC)
        fRe(i) = real(fC, kind=rkind); fIm(i) = aimag(fC)
        fxRe(i) = real(fxC, kind=rkind); fxIm(i) = aimag(fxC)
        fxxRe(i) = real(fxxc, kind=rkind); fxxIm(i) = aimag(fxxC)
    end do
    call cpu_time(end_time_fn)

! Evaluate using splines
    call cpu_time(start_time_spline)

    s_abs = 0.; s_rel = 0.
    do i = 1, nxC_get
        call cspeval(xC(i),iselectC,fvalRe,xC_grid,nxC,ilinx,fsplRe,ier)
        call cspeval(xC(i),iselectC,fvalIm,xC_grid,nxC,ilinx,fsplIm,ier)
        if (ier .ne. 0) write (*,*) 'cspeval: ier = ', ier

        abs_err = abs(cmplx(fvalRe(1) - fRe(i), fvalIm(1) - fIm(i), kind=rkind))
        rel_err = abs_err/abs(cmplx(fRe(i), fIm(i), kind=rkind))
        s_abs = s_abs + abs_err
        s_rel = s_rel + rel_err

       if (write_details) then
			write (*,*) ' '
			write (*,*) 'x = ', xC(i)
			write (*,*) '  fRe(i) = ', fRe(i), '  Re(f spline) = ', fvalRe(1)
			write (*,*) '  fIm(i) = ', fIm(i), '  Im(f spline) = ', fvalIm(1)
			write (*,*) '  fxRe(i) = ', fxRe(i), '  Re(fx spline) = ', fvalRe(2)
			write (*,*) '  fxIm(i) = ', fxIm(i), '  Im(fx spline) = ', fvalIm(2)
			write (*,*) '  fxxRe(i) = ', fxxRe(i), '  Re(fxx spline) = ', fvalRe(3)
			write (*,*) '  fxxIm(i) = ', fxxIm(i), '  Im(fxx spline) = ', fvalIm(3)
			write (*,*) 'abs_err = ', abs_err, '  rel_err = ', rel_err
        end if

    end do
    call cpu_time(end_time_spline)

    write(*,*) ' '
    write(*,*) 'Splined Z function  grid points = ', nxC, '  sample points = ', nxC_get
    write(*,*) 'CPU time coeff eval = ', end_time_coef - start_time_coef
    write(*,*) 'CPU time function eval = ', end_time_fn - start_time_fn
    write(*,*) 'CPU time spline eval = ', end_time_spline - start_time_spline
    write(*,*) 'Average absolute error = ', s_abs/nxC_get
    write(*,*) 'Average relative error = ', s_rel/nxC_get


    write(*,*) ' '
    write (*,*) ' '
    write (*,*) 'Zfunction real part spline, imaginary part analytic'

! Evaluate using splines
    call cpu_time(start_time_spline)

    s_abs = 0.; s_rel = 0.
    fvalRe(2:3) = 0.; fvalIm(2:3) = 0.
    do i = 1, nxC_get
        call cspeval(xC(i),iselectC,fvalRe,xC_grid,nxC,ilinx,fsplRe,ier)
        if (ier .ne. 0) write (*,*) 'cspeval: ier = ', ier
        fvalIm(1) = sqrt_pi*exp(-xC(i)**2)

        abs_err = abs(cmplx(fvalRe(1) - fRe(i), fvalIm(1) - fIm(i), kind=rkind))
        rel_err = abs_err/abs(cmplx(fRe(i), fIm(i), kind=rkind))
        s_abs = s_abs + abs_err
        s_rel = s_rel + rel_err

       if (write_details) then
			write (*,*) ' '
			write (*,*) 'x = ', xC(i)
			write (*,*) '  fRe(i) = ', fRe(i), '  Re(f spline) = ', fvalRe(1)
			write (*,*) '  fIm(i) = ', fIm(i), '  Im(f spline) = ', fvalIm(1)
			write (*,*) '  fxRe(i) = ', fxRe(i), '  Re(fx spline) = ', fvalRe(2)
			write (*,*) '  fxIm(i) = ', fxIm(i), '  Im(fx spline) = ', fvalIm(2)
			write (*,*) '  fxxRe(i) = ', fxxRe(i), '  Re(fxx spline) = ', fvalRe(3)
			write (*,*) '  fxxIm(i) = ', fxxIm(i), '  Im(fxx spline) = ', fvalIm(3)
			write (*,*) 'abs_err = ', abs_err, '  rel_err = ', rel_err
        end if

    end do
    call cpu_time(end_time_spline)

    write(*,*) ' '
    write(*,*) 'Z (spline/analytic)  grid points = ', nxC, '  sample points = ', nxC_get
    write(*,*) 'CPU time coeff eval = ', end_time_coef - start_time_coef
    write(*,*) 'CPU time function eval = ', end_time_fn - start_time_fn
    write(*,*) 'CPU time spline eval = ', end_time_spline - start_time_spline
    write(*,*) 'Average absolute error = ', s_abs/nxC_get
    write(*,*) 'Average relative error = ', s_rel/nxC_get

  end if Z

!***************************************************************************************
! 1D_module tests

  D1_module: if (do_1D_module .eqv. .true.) then

    write (*,*)
    write (*,*) '1D_module results'

    call cpu_time(start_time_coef)

! 1D Grid
    do i = 1, nx1D
        x1D_grid(i) = x1D_grid_min + (i-1)*(x1D_grid_max - x1D_grid_min)/(nx1D-1)
    end do

! Grid is generated above, evaluate function array
    do i = 1, nx1D
        call spline_test_fn_1D(x1D_grid(i), f, fx, fxx)
        fspl1D(1,i) = f
    end do

    func_name = 'cos(k*x)'
    call spl_func%cube_spline_1D_init(nx1D, x1D_grid, fspl1D(1,:), func_name)

    call cpu_time(end_time_coef)

! Evaluate splined function on data

    do i = 1, nx1D_get
        call random_number(random)
        x1D(i) = x1D_get_min + (x1D_get_max - x1D_get_min)*random
        x1D(i) = x1D_get_min + (i-1)*(x1D_get_max - x1D_get_min)/(nx1D_get-1)/sqrt(2.)
    end do

! Evaluate the function analytically for timing
    call cpu_time(start_time_fn)
    do i = 1, nx1D_get
        call spline_test_fn_1D(x1D(i), f1D(i), fx1D(i), fxx1D(i))
    end do
    call cpu_time(end_time_fn)

! Evaluate using splines
    call cpu_time(start_time_spline)

    s_abs = 0._rkind; s_rel = 0._rkind
    do i = 1, nx1D_get
!         call spl_func%eval_1D_f(x1D(i), fsp)
!         call spl_func%eval_1D_fp(x1D(i), fsp, fxsp)
        call spl_func%eval_1D_fpp(x1D(i), fsp, fxsp, fxxsp)

        abs_err = fsp - f1D(i)
        rel_err = abs_err/f1D(i)
        s_abs = s_abs + abs(abs_err)
        s_rel = s_rel + abs(rel_err)

   if (write_details) then
        write (*,*) ' '
        write (*,*) 'x = ', x1D(i), '  f1D(i) = ', f1D(i), '  f spline = ', fsp
        write (*,*) 'abs_err = ', abs_err, '  rel_err = ', rel_err
        write (*,*) 'fx = ', fx1D(i), '  fx spline = ', fxsp,'  fxx = ', fxx1D(i),&
               &  '  fxx spline = ', fxxsp
    end if

    end do
    call cpu_time(end_time_spline)

 		write(*,*) ' '
		write(*,*) '1 + SIN(k x)  grid points = ', nx1D, '  sample points = ', nx1D_get
		write(*,*) 'CPU time coeff eval = ', end_time_coef - start_time_coef
		write(*,*) 'CPU time function eval = ', end_time_fn - start_time_fn
		write(*,*) 'CPU time spline eval = ', end_time_spline - start_time_spline
		write(*,*) 'Average absolute error = ', s_abs/nx1D_get
		write(*,*) 'Average relative error = ', s_rel/nx1D_get

  end if D1_module

!***************************************************************************************
! 2D_module tests

  D2_module: if (do_2D_module .eqv. .true.) then

    write (*,*)
    write (*,*) '2D_module results'

    call cpu_time(start_time_coef)

! Load grids and function array
    do j = 1, ny
        y(j) = y_grid_min + (j-1)*(y_grid_max - y_grid_min)/(ny-1)
    end do
    do i = 1, nx
        x(i) = x_grid_min + (i-1)*(x_grid_max - x_grid_min)/(nx-1)
    end do

    do j = 1, ny
        do i = 1, nx
            call spline_test_fn(x(i), y(j), f, fx, fy, fxx, fxy, fyy)
            fspl_2D(1,1,i,j) = f
         end do
    end do

    func_name = 'sqrt(x**2 + y**2)'

	call spl_func_2D%cube_spline_2D_init(nx, x, ny, y,fspl_2D(1,1,:,:), func_name)

    call cpu_time(end_time_coef)

! Evaluate splined function on data

! Evaluate the function analytically for timing
    call cpu_time(start_time_fn)
    do i = 1, nr
        do j = 1, nTheta
            r = real(i)
            theta = (j-1)*pi/nTheta
            xget = r*cos(theta)
            yget = r*sin(theta)
            rget = sqrt(xget**2 + yget**2)
            call spline_test_fn(xget, yget, f, fx, fy, fxx, fxy, fyy)
        end do
    end do
    call cpu_time(end_time_fn)

! Evaluate using splines
    call cpu_time(start_time_spline)

    s_abs = 0._rkind; s_rel = 0._rkind

    do i = 1, nr
        do j = 1, nTheta
            r = real(i)
            theta = (j-1)*pi/nTheta

            xget = r*cos(theta)
            yget = r*sin(theta)
            rget = sqrt(xget**2 + yget**2)

            call spl_func_2D%eval_2D_fpp(xget, yget, f2sp, f2xsp, f2ysp, f2xxsp, f2xysp, f2yysp)

            call spline_test_fn(xget, yget, f, fx, fy, fxx, fxy, fyy)

			abs_err = f2sp - f
			rel_err = abs_err/f
			s_abs = s_abs + abs(abs_err)
			s_rel = s_rel + abs(rel_err)

            if (write_details) then
				write (*,*)
				write (*,*) 'xget = ', xget, '     yget = ', yget
				write (*,*) 'fget = ', f, '  f2sp = ', f2sp
				write (*,*) 'fx = ', fx, '  f2xsp = ', f2xsp
				write (*,*) 'fy = ', fy, '  f2ysp = ', f2ysp
				write (*,*) 'fxx = ', fxx, '  f2xxsp = ', f2xxsp
				write (*,*) 'fxy = ', fxy, '  f2xysp = ', f2xysp
				write (*,*) 'fyy = ', fyy, '  f2yysp = ', f2yysp
            end if
        end do
    end do

    call cpu_time(end_time_spline)

 		write(*,*) ' '
		write(*,*) 'sqrt(x**2 + y**2)'
		write(*,*) 'x grid points = ', nx, ' y grid points =  = ', ny
		write(*,*) 'CPU time coeff eval = ', end_time_coef - start_time_coef
		write(*,*) 'CPU time function eval = ', end_time_fn - start_time_fn
		write(*,*) 'CPU time spline eval = ', end_time_spline - start_time_spline
		write(*,*) 'Average absolute error = ', s_abs/(nr*nTheta)
		write(*,*) 'Average relative error = ', s_rel/(nr*nTheta)

  end if D2_module

    end program test_psplines