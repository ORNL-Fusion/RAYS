    program test_psplines

    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals

! setup

    integer, parameter :: nxMax = 101, nyMax = 101
    integer, parameter :: nx = nxMax, ny = nyMax, nr = 5, nTheta = 10, nx1D = 11

    integer :: i, j
    real(KIND=rkind) ::  ff, fx, fy, fxx, fxy, fyy
    real(KIND=rkind) :: r, theta, rget
    real(KIND=rkind) :: pi = 3.1415926

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
! bcspeval arguments, if not already declared above

    real(KIND=rkind) ::  fspl_1D(4, nx), fval1D(3)
    integer :: iselect1D(3), iwk

!***************************************************************************************
! 2D splines

! Load grids and function array
    do j = 1, ny
        y(j) = real(j-11)/10.
    end do
    do i = 1, nx
        x(i) = real(i-11)/10.
    end do

    do j = 1, ny
        do i = 1, nx
            call spline_test_fn(x(i), y(j), ff, fx, fy, fxx, fxy, fyy)
            fspl_2D(1,1,i,j) = ff
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


    write (*,*) 'ilinx = ', ilinx, '  ilinth = ', ilinth
!    write (*,*) f(1,1,:,:)

! Evaluate splined function on data

    iselect = 1 ! Selector for outputs
 !   iselect(1) = 1 ! Output function

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

            call spline_test_fn(xget, yget, fget, fx, fy, fxx, fxy, fyy)
            write (*,*)
            write (*,*) 'xget = ',xget, '  yget = ', yget, '  fget = ', fget, '  fval(1) = ', fval(1)
            write (*,*) 'fx = ', fx, '  fval(2) = ', fval(2),'fy = ', fy, '  fval(3) = ', fval(3)
            write (*,*) 'fxx = ', fxx, '  fval(4) = ', fval(4),'fxy = ', fxy, &
            &  '  fval(6) = ', fval(6), '  fyy = ', fyy, '  fval(5) = ', fval(5)


        end do
    end do


!***************************************************************************************
! 1D splines
    write (*,*)
    write (*,*) '1D results'

! Grid is generated above, evaluate function array

    do i = 1, nx
        call spline_test_fn_1D(x(i), ff, fx, fxx)
        fspl_1D(1,i) = ff
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

    call cspline(x,nx,fspl_1D,ibcxmin,bcxmin,ibcxmax,bcxmax,wk,iwk,ilinx,ier)
    write (*,*) 'ilinx = ', ilinx

! Evaluate splined function on data

    iselect1D = (/ 1, 1, 1 /)

    do i = 1, nx1D

        xget = 0.4*(i - int(real(nx1D-1)/2.))

        call cspeval(xget,iselect1D,fval1D,x,nx,ilinx,fspl_1D,ier)
        if (ier .ne. 0) write (*,*) 'cspeval: ier = ', ier

        call spline_test_fn_1D(xget,fget, fx, fxx)
        write (*,*) 'x =   ',xget
        write (*,*) 'f =   ', fget, '  fval(1) = ', fval1D(1), '  err = ', fval1D(1) - fget
        write (*,*) 'fx =  ', fx, '  fval(2) = ', fval1D(2), '  err = ', fval1D(2) - fx
        write (*,*) 'fxx = ', fxx, '  fval(3) = ', fval1D(3), '  err = ', fval1D(3) - fxx

    end do

    end program test_psplines