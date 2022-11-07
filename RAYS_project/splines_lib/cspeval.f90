!  cspeval -- eval cubic spline function and/or derivatives
!
subroutine cspeval(xget,iselect,fval,x,nx,ilinx,f,ier)
  use iso_c_binding, only: fp => c_double
  !
  !============
  implicit none
  integer ier,nx
  !============
  real(fp) :: xget                         ! interpolation target
  real(fp) :: fval(3)                      ! output values
  real(fp) :: x(nx),f(4,nx)                ! spline data
  !
  integer iselect(3)                ! output selector
  integer ilinx                     ! =1 if x(...) is evenly spaced
  !
  !  modification -- dmc 11 Jan 1999 -- remove SAVE stmts; break routine
  !    into these parts:
  !
  !    cspevx -- find grid cell of target pt.
  !    cspevfn -- evaluate function using output of bcpsevxy
  !
  !    in cases where multiple functions are defined on the same grid,
  !    time can be saved by using cspevx once and then cspevfn
  !    multiple times.
  !
  !  input:
  !     (xget)        location where interpolated value is desired
  !                   x(1).le.xget.le.x(nx) expected
  !
  !     iselect       select desired output
  !
  !                     iselect(1)=1 -- want function value (f) itself
  !                     iselect(2)=1 -- want  df/dx
  !                     iselect(3)=1 -- want  d2f/dx2
  !
  !              example:  iselect(1)=iselect(2)=iselect(3)=1
  !                            f, df/dx, and d2f/dx2 are all evaluated
  !                            and returned in order in fval(1), fval(2),
  !                            and fval(3)
  !                        iselect(1)=0, iselect(2)=1, iselect(3)=0
  !                            only the 1st derivative is evaluated
  !                            and returned in fval(1).
  !
  !                     set iselect(1)=3 to get d3f/dx3, 1 value only.
  !
  !                   see fval (output) description.
  !
  !     x(1...nx)     independent coordinate x, strict ascending
  !
  !     ilinx  --  =1: flag that x is linearly spaced (avoid search for speed)
  !
  !  **CAUTION** actual even spacing of x, is NOT CHECKED HERE!
  !
  !     f             the function values (at grid points) and spline coefs
  !
  !  evaluation formula:  for point x btw x(i) and x(i+1), dx=x-x(i)
  !
  !      spline value =
  !        f(1,i) + dx*f(2,i) + dx**2*f(3,i) + dx**3*f(4,i)
  !
  !  output:
  !      up to 3 elements of fval, ordered as follows:
  !        fval(1)=function value or lowest order derivative requested
  !        fval(2)=next order derivative
  !             etc
  !        the ordering is a subset of the sequence given under the "iselect"
  !        description.
  !
  !      ier = 0 -- successful completion; = 1 -- an error occurred.
  !
  !-------------------------------------------------------------------
  !  local
  !
  integer :: ia(1) = (/ 0 /)
  real(fp) :: dxa(1)
  !
  !--------------------------
  !
  call cspevx(xget,x,nx,ilinx,ia(1),dxa(1),ier)
  if(ier.ne.0) return
  !
  call cspevfn(iselect,1,1,fval,ia,dxa,f,nx)
  !
  return
end subroutine cspeval
!
!-------------------------------------------------------------------------
!  cspevx -- look up x zone
!
!  this is the "first part" of cspeval, see comments, above.
!
subroutine cspevx(xget,x,nx,ilinx,i,dx,ier)
  use iso_c_binding, only: fp => c_double
  !
  !============
  implicit none
  integer nxm,ii
  !============
  real(fp) :: zxget,zxtol
  !============
  integer nx                        ! x array dimension
  !
  real(fp) :: xget                         ! target point
  real(fp) :: x(nx)                        ! independent coord. array
  !
  integer ilinx                     ! =1:  assume x evenly spaced
  !
  !  output of cspevx
  !
  integer i                         ! index to cell containing target pt
  real(fp) :: dx                           ! displacement of target pt w/in cell
  ! dx = x-x(i)
  !
  !  the input argument range is checked...
  !
  integer ier                       ! return ier.ne.0 on error
  !
  !------------------------------------
  !
  ier=0
  !
  !  range check
  !
  zxget=xget
  if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
     zxtol=4.0E-7_fp*max(abs(x(1)),abs(x(nx)))
     if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
        ier=1
        write(6,1001) xget,x(1),x(nx)
1001    format(' ?cspeval:  xget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((xget.lt.x(1)-0.5_fp*zxtol).or. &
             (xget.gt.x(nx)+0.5_fp*zxtol)) &
             write(6,1011) xget,x(1),x(nx)
1011    format(' %cspeval:  xget=',1pe15.8,' beyond range ', &
             1pe15.8,' to ',1pe15.8,' (fixup applied)')
        if(xget.lt.x(1)) then
           zxget=x(1)
        else
           zxget=x(nx)
        end if
     end if
  end if
  if(ier.ne.0) return
  !
  !  now find interval in which target point lies..
  !
  nxm=nx-1
  !
  if(ilinx.eq.1) then
     ii=1+nxm*(zxget-x(1))/(x(nx)-x(1))
     i=min(nxm, ii)
     if(zxget.lt.x(i)) then
        i=i-1
     else if(zxget.gt.x(i+1)) then
        i=i+1
     end if
  else
     if((1.le.i).and.(i.lt.nxm)) then
        if((x(i).le.zxget).and.(zxget.le.x(i+1))) then
           continue  ! already have the zone
        else
           call zonfind(x,nx,zxget,i)
        end if
     else
        i=nx/2
        call zonfind(x,nx,zxget,i)
     end if
  end if
  !
  dx=zxget-x(i)
  !
  return
end subroutine cspevx
!------------------------------------------------------------------------
!  cspevfn -- OK now evaluate the cubic spline
!
subroutine cspevfn(ict,ivec,ivd,fval,iv,dxv,f,nx)
  use iso_c_binding, only: fp => c_double
  !
  !  input:
  !============
  implicit none
  integer nx,iaval,i
  !============
  real(fp) :: dx
  !============
  integer ict(3)                    ! selector:
  !        ict(1)=1 for f      (don't evaluate f if ict(1)=0)
  !        ict(2)=1 for df/dx   ""
  !        ict(3)=1 for d2f/dx2
  !
  !        set ict(1)=3 to get d3f/dx3 (only)
  !
  integer ivec,ivd                  ! vector dimensioning
  !
  !    ivec-- number of vector pts (spline values to look up)
  !    ivd -- 1st dimension of fval, .ge.ivec
  !
  ! output:
  real(fp) :: fval(ivd,*)                 ! output array
  !
  !    v = index to element in vector;
  !  fval(v,1) = first item requested by ict(...),
  !  fval(v,2) = 2nd item requested,  ...etc...
  !
  !  input:
  integer iv(ivec)                  ! grid cell indices -- vectors
  real(fp) :: dxv(ivec)                    ! displacements w/in cell -- vectors
  !
  real(fp) :: f(4,nx)                      ! cubic fcn spline coeffs array
  !
  !  usage example:
  !
  !  1.  for each element xx(v) in a vector of x values:
  !    find the x zone index and displacement with respect to the
  !    lower end point of the zone; store thes in vectors iv and dxv.
  !
  !  2.  set ict(1)=0, ict(2)=1, ict(3)=0 -- get only the 1st derivative
  !
  !  3.  ivec is the length of the vector; ivd is the 1st dimension of the
  !      array fval to receive the output
  !
  !      real fval(ivd,3)
  !      real xv(ivd)
  !      real dxv(ivd)
  !      integer iv(ivd)
  !      integer ict(3)
  !
  !      real fspline(4,nx)  ! spline coeffs
  !      data ict/0,1,0/     ! this call:  want 1st derivative only
  !                          !  this will be output to fval(*,1).
  !      ...
  !      do iv=1,ivec
  !        ...               ! comput indices & displacements
  !      end do
  !      call cspevfn(ict,ivec,ivd,fval,iv,dxv,fspline,nx)
  !
  !--------------------------------------------------------------------
  !  local:
  !
  integer v                         ! vector element index
  !
  !  OK can now do evaluations
  !
  iaval=0  ! fval addressing
  !
  if(ict(1).eq.3) then
     !
     !  fxxx = d3f/dx3
     !
     iaval=iaval+1
     do v=1,ivec
        i=iv(v)
        fval(v,iaval)=6.0_fp*f(4,i)
     end do
  else
     !
     !  normal call...
     !
     if(ict(1).gt.0) then
        !  evaluate f
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           dx=dxv(v)
           fval(v,iaval)=f(1,i)+dx*(f(2,i)+dx*(f(3,i)+dx*f(4,i)))
        end do
     end if
     !
     if(ict(2).gt.0) then
        !  evaluate df/dx
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           dx=dxv(v)
           fval(v,iaval)=f(2,i)+dx*(2.0_fp*f(3,i)+dx*3.0_fp*f(4,i))
        end do
     end if
     !
     if(ict(3).gt.0) then
        !  evaluate d2f/dx2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           dx=dxv(v)
           fval(v,iaval)=2.0_fp*f(3,i)+dx*6.0_fp*f(4,i)
        end do
     end if
  end if
  !
  return
end subroutine cspevfn
!----------------------
