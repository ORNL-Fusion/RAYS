!  cspline -- dmc 15 Feb 1999
!
!  a standard interface to the 1d spline setup routine
!    modified dmc 3 Mar 2000 -- to use Wayne Houlberg's v_spline code.
!    new BC options added.
!
subroutine cspline(x,nx,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax, &
     wk,iwk,ilinx,ier)
  use iso_c_binding, only: fp => c_double
  !
  !============
  implicit none
  integer iwk,nx,ierx,inum,i
  !============
  real(fp) :: half,sixth
  !============
  real(fp) :: x(nx)                        ! x axis (in)
  real(fp) :: fspl(4,nx)                   ! spline data (in/out)
  integer ibcxmin                   ! x(1) BC flag (in, see comments)
  real(fp) :: bcxmin                       ! x(1) BC data (in, see comments)
  integer ibcxmax                   ! x(nx) BC flag (in, see comments)
  real(fp) :: bcxmax                       ! x(nx) BC data (in, see comments)
  real(fp) :: wk(iwk)                      ! workspace of size at least nx
  integer ilinx                     ! even spacing flag (out)
  integer ier                       ! output, =0 means OK
  !
  !  ** note wk(...) array is not used unless ibcxmin=-1 (periodic spline
  !  evaluation)
  !
  !  this routine computes spline coefficients for a 1d spline --
  !  evaluation of the spline can be done by cspeval.f90 subroutines
  !  or directly by inline code.
  !
  !  the input x axis x(1...nx) must be strictly ascending, i.e.
  !  x(i+1).gt.x(i) is required for i=1 to nx-1.  This is checked and
  !  ier=1 is set and the routine exits if the test is not satisfied.
  !
  !  on output, ilinx=1 is set if, to a reasonably close tolerance,
  !  all grid spacings x(i+1)-x(i) are equal.  This allows a speedier
  !  grid lookup algorithm on evaluation of the spline.  If on output
  !  ilinx=2, this means the spline x axis is not evenly spaced.
  !
  !  the input data for the spline are given in f[j] = fspl(1,j).  The
  !  output data are the spline coefficients fspl(2,j),fspl(3,j), and
  !  fspl(4,j), j=1 to nx.  The result is a spline s(x) satisfying the
  !  boundary conditions and with the properties
  !
  !     s(x(j)) = fspl(1,j)
  !     s'(x) is continuous even at the grid points x(j)
  !     s''(x) is continuous even at the grid points x(j)
  !
  !  the formula for evaluation of s(x) is:
  !
  !     let dx = x-x(i), where x(i).le.x.le.x(i+1).  Then,
  !     s(x)=fspl(1,i) + dx*(fspl(2,i) +dx*(fspl(3,i) + dx*fspl(4,i)))
  !
  !  ==>boundary conditions.  Complete specification of a 1d spline
  !  requires specification of boundary conditions at x(1) and x(nx).
  !
  !  this routine provides 4 options:
  !
  ! -1 ***** PERIODIC BC
  !  ibcxmin=-1  --  periodic boundary condition.  This means the
  !    boundary conditions s'(x(1))=s'(x(nx)) and s''(x(1))=s''(x(nx))
  !    are imposed.  Note that s(x(1))=s(x(nx)) (i.e. fspl(1,1)=fspl(1,nx))
  !    is not required -- that is determined by the fspl array input data.
  !    The periodic boundary condition is to be preferred for periodic
  !    data.  When splining periodic data f(x) with period P, the relation
  !    x(nx)=x(1)+n*P, n = the number of periods (usually 1), should hold.
  !    (ibcxmax, bcxmin, bcxmax are ignored).
  !
  !  if a periodic boundary condition is set, this covers both boundaries.
  !  for the other types of boundary conditions, the type of condition
  !  chosen for the x(1) boundary need not be the same as the type chosen
  !  for the x(nx) boundary.
  !
  !  0 ***** NOT A KNOT BC
  !  ibcxmin=0 | ibcxmax=0 -- this specifies a "not a knot" boundary
  !    condition -- see cubsplb.f90.  This is a common way for inferring
  !    a "good" spline boundary condition automatically from data in the
  !    vicinity of the boundary.  (bcxmin | bcxmax are ignored).
  !
  !  1 ***** BC:  SPECIFIED SLOPE
  !  ibcxmin=1 | ibcxmax=1 -- boundary condition is to have s'(x(1)) |
  !    s'(x(nx)) match the passed value (bcxmin | bcxmax).
  !
  !  2 ***** BC:  SPECIFIED 2nd DERIVATIVE
  !  ibcxmin=2 | ibcxmax=2 -- boundary condition is to have s''(x(1)) |
  !    s''(x(nx)) match the passed value (bcxmin | bcxmax).
  !
  !  3 ***** BC:  SPECIFIED SLOPE = 0.0
  !  ibcxmin=3 | ibcxmax=3 -- boundary condition is to have s'(x(1)) |
  !    s'(x(nx)) equal to ZERO.
  !
  !  4 ***** BC:  SPECIFIED 2nd DERIVATIVE = 0.0
  !  ibcxmin=4 | ibcxmax=4 -- boundary condition is to have s''(x(1)) |
  !    s''(x(nx)) equal to ZERO.
  !
  !  5 ***** BC:  1st DIVIDED DIFFERENCE
  !  ibcxmin=5 | ibcxmax=5 -- boundary condition is to have s'(x(1)) |
  !    s'(x(nx)) equal to the slope from the 1st|last 2 points
  !
  !  6 ***** BC:  2nd DIVIDED DIFFERENCE
  !  ibcxmin=6 | ibcxmax=6 -- boundary condition is to have s''(x(1)) |
  !    s''(x(nx)) equal to the 2nd derivative from the 1st|last 3 points
  !
  !  7 ***** BC:  3rd DIVIDED DIFFERENCE
  !  ibcxmin=7 | ibcxmax=7 -- boundary condition is to have s'''(x(1)) |
  !    s'''(x(nx)) equal to the 3rd derivative from the 1st|last 4 points
  !
  !---------------------------------------------------------------------
  data half/0.5_fp/
  data sixth/0.166666666666666667_fp/
  !
  !  error checks
  !
  ier = 0
  if(nx.lt.2) then
     write(6,'('' ?cspline:  at least 2 x points required.'')')
     ier=1
  end if
  call ibc_ck(ibcxmin,'cspline','xmin',-1,7,ier)
  if(ibcxmin.ge.0) call ibc_ck(ibcxmax,'cspline','xmax',0,7,ier)
  !
  !  x axis check
  !
  call splinck(x,nx,ilinx,1.0E-3_fp,ierx)
  if(ierx.ne.0) ier=2
  !
  if(ier.eq.2) then
     write(6,'('' ?cspline:  x axis not strict ascending'')')
  end if
  !
  if(ibcxmin.eq.-1) then
     inum=nx
     if(iwk.lt.inum) then
        write(6,1009) inum,iwk,nx
1009    format( &
             ' ?cspline:  workspace too small.  need:  ',i6,' got:  ',i6/ &
             '  (need = nx, nx=',i6)
        ier=3
     end if
  end if
  !
  if(ier.ne.0) return
  !
  !  OK -- evaluate spline
  !
  if(ibcxmin.eq.1) then
     fspl(2,1)=bcxmin
  else if(ibcxmin.eq.2) then
     fspl(3,1)=bcxmin
  end if
  !
  if(ibcxmax.eq.1) then
     fspl(2,nx)=bcxmax
  else if(ibcxmax.eq.2) then
     fspl(3,nx)=bcxmax
  end if
  !
  call v_spline(ibcxmin,ibcxmax,nx,x,fspl,wk)
  !
  do i=1,nx
     fspl(3,i)=half*fspl(3,i)
     fspl(4,i)=sixth*fspl(4,i)
  end do
  !
  return
end subroutine cspline
