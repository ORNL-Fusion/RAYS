!
!  bcspeval -- eval bicubic spline function and/or derivatives
!
subroutine bcspeval(xget,yget,iselect,fval, &
     x,nx,y,ny,ilinx,iliny,f,inf3,ier)
  use iso_c_binding, only: fp => c_double
  !
  implicit none
  integer iselect(6)
  integer ilinx,iliny,nx,ny,inf3,ier
  !
  real(fp) :: xget,yget
  real(fp) :: fval(*)
  real(fp) :: x(nx),y(ny),f(4,4,inf3,ny)
  !
  !  modification -- dmc 11 Jan 1999 -- remove SAVE stmts; break routine
  !    into these parts:
  !
  !    bcspevxy -- find grid cell of target pt.
  !    bcspevfn -- evaluate function using output of bcpsevxy
  !
  !    in cases where multiple functions are defined on the same grid,
  !    time can be saved by using bcspevxy once and then bcspevfn
  !    multiple times.
  !
  !  input:
  !     (xget,yget)   location where interpolated value is desired
  !                   x(1).le.xget.le.x(nx) expected
  !                   y(1).le.yget.le.y(ny) expected
  !
  !     iselect       select desired output
  !
  !                     iselect(1)=1 -- want function value (f) itself
  !                     iselect(2)=1 -- want  df/dx
  !                     iselect(3)=1 -- want  df/dy
  !                     iselect(4)=1 -- want  d2f/dx2
  !                     iselect(5)=1 -- want  d2f/dy2
  !                     iselect(6)=1 -- want  d2f/dxdy
  !
  !              example:  iselect(1)=iselect(2)=iselect(3)=1
  !                            f, df/dx, and df/dy all evaluated
  !                        iselect(4)=iselect(5)=iselect(6)=0
  !                            2nd derivatives not evaluated.
  !
  !                   the number of non zero values iselect(1:6)
  !                   determines the number of outputs...
  !                   see fval (output) description.
  !
  !  new dmc December 2005 -- access to higher derivatives (even if not
  !  continuous-- but can only go up to 3rd derivatives on any one coordinate.
  !     if iselect(1)=3 -- want 3rd derivatives
  !          iselect(2)=1 for d3f/dx3
  !          iselect(3)=1 for d3f/dx2dy
  !          iselect(4)=1 for d3f/dxdy2
  !          iselect(5)=1 for d3f/dy3
  !               number of non-zero values iselect(2:5) gives no. of outputs
  !     if iselect(1)=4 -- want 4th derivatives
  !          iselect(2)=1 for d4f/dx3dy
  !          iselect(3)=1 for d4f/dx2dy2
  !          iselect(4)=1 for d4f/dxdy3
  !               number of non-zero values iselect(2:4) gives no. of outputs
  !     if iselect(1)=5 -- want 5th derivatives
  !          iselect(2)=1 for d5f/dx3dy2
  !          iselect(3)=1 for d5f/dx2dy3
  !               number of non-zero values iselect(2:3) gives no. of outputs
  !     if iselect(1)=6 -- want 6th derivatives
  !          d6f/dx3dy3 -- one value is returned.
  !
  !     x(1...nx)     independent coordinate x, strict ascending
  !     y(1...ny)     independent coordinate y, strict ascending
  !
  !     ilinx  --  =1: flag that x is linearly spaced (avoid search for speed)
  !     iliny  --  =1: flag that y is linearly spaced (avoid search for speed)
  !
  !  **CAUTION** actual even spacing of x, y is NOT CHECKED HERE!
  !
  !
  !     f             the function values (at grid points) and spline coefs
  !
  !  evaluation formula:  for point x btw x(i) and x(i+1), dx=x-x(i)
  !                             and y btw y(j) and y(j+1), dy=y-y(j),
  !
  !      spline value =
  !        f(1,1,i,j) + dx*f(2,1,i,j) + dx**2*f(3,1,i,j) + dx**3*f(4,1,i,j)
  !   +dy*(f(1,2,i,j) + dx*f(2,2,i,j) + dx**2*f(3,2,i,j) + dx**3*f(4,2,i,j))
  !   +d2*(f(1,3,i,j) + dx*f(2,3,i,j) + dx**2*f(3,3,i,j) + dx**3*f(4,3,i,j))
  !   +d3*(f(1,4,i,j) + dx*f(2,4,i,j) + dx**2*f(3,4,i,j) + dx**3*f(4,4,i,j))
  !
  !      where d2=dy**2 and d3=dy**3.
  !
  !  output:
  !      up to 6 elements of fval, ordered as follows:
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
  integer :: i(1)
  integer :: j(1)
  !
  real(fp) :: dx(1),dy(1)
  !
  !--------------------------
  !
  i(1)=0
  j(1)=0
  !
  call bcspevxy(xget,yget,x,nx,y,ny,ilinx,iliny, &
       i(1),j(1),dx(1),dy(1),ier)
  if(ier.ne.0) return
  !
  call bcspevfn(iselect,1,1,fval,i,j,dx,dy,f,inf3,ny)
  !
  return
end subroutine bcspeval
!
!-------------------------------------------------------------------------
!  bcspevxy -- look up x-y zone
!
!  this is the "first part" of bcspeval, see comments, above.
!
subroutine bcspevxy(xget,yget,x,nx,y,ny,ilinx,iliny, &
     i,j,dx,dy,ier)
  use iso_c_binding, only: fp => c_double
  !
  !============
  implicit none
  integer nxm,nym,ii,jj
  !============
  real(fp) :: zxget,zyget,zxtol,zytol
  !============
  integer nx,ny                     ! array dimensions
  !
  real(fp) :: xget,yget                    ! target point
  real(fp) :: x(nx),y(ny)                  ! indep. coords.
  !
  integer ilinx                     ! =1:  assume x evenly spaced
  integer iliny                     ! =1:  assume y evenly spaced
  !
  !  output of bcspevxy
  !
  integer i,j                       ! index to cell containing target pt
  real(fp) :: dx,dy                        ! displacement of target pt w/in cell
  ! dx=x-x(i)  dy=y-y(j)
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
  zyget=yget

  if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
     zxtol=4.0E-7_fp*max(abs(x(1)),abs(x(nx)))
     if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
        ier=1
        write(6,1001) xget,x(1),x(nx)
1001    format(' ?bcspeval:  xget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((xget.lt.x(1)-0.5_fp*zxtol).or. &
             (xget.gt.x(nx)+0.5_fp*zxtol)) &
             write(6,1011) xget,x(1),x(nx)
1011    format(' %bcspeval:  xget=',1pe15.8,' beyond range ', &
             1pe15.8,' to ',1pe15.8,' (fixup applied)')
        if(xget.lt.x(1)) then
           zxget=x(1)
        else
           zxget=x(nx)
        end if
     end if
  end if
  if((yget.lt.y(1)).or.(yget.gt.y(ny))) then
     zytol=4.0E-7_fp*max(abs(y(1)),abs(y(ny)))
     if((yget.lt.y(1)-zytol).or.(yget.gt.y(ny)+zytol)) then
        ier=1
        write(6,1002) yget,y(1),y(ny)
1002    format(' ?bcspeval:  yget=',1pe11.4,' out of range ', &
             1pe11.4,' to ',1pe11.4)
     else
        if((yget.lt.y(1)-0.5_fp*zytol).or.(yget.gt.y(ny)+0.5_fp*zytol)) &
             write(6,1012) yget,y(1),y(ny)
1012    format(' %bcspeval:  yget=',1pe15.8,' beyond range ', &
             1pe15.8,' to ',1pe15.8,' (fixup applied)')
        if(yget.lt.y(1)) then
           zyget=y(1)
        else
           zyget=y(ny)
        end if
     end if
  end if
  if(ier.ne.0) return
  !
  !  now find interval in which target point lies..
  !
  nxm=nx-1
  nym=ny-1
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
  if(iliny.eq.1) then
     jj=1+nym*(zyget-y(1))/(y(ny)-y(1))
     j=min(nym, jj)
     if(zyget.lt.y(j)) then
        j=j-1
     else if(zyget.gt.y(j+1)) then
        j=j+1
     end if
  else
     if((1.le.j).and.(j.lt.nym)) then
        if((y(j).le.zyget).and.(zyget.le.y(j+1))) then
           continue  ! already have the zone
        else
           call zonfind(y,ny,zyget,j)
        end if
     else
        j=ny/2
        call zonfind(y,ny,zyget,j)
     end if
  end if
  !
  dx=zxget-x(i)
  dy=zyget-y(j)
  !
  return
end subroutine bcspevxy
!------------------------------------------------------------------------
!  bcspevfn -- OK now evaluate the bicubic spline
!
subroutine bcspevfn(ict,ivec,ivd,fval,iv,jv,dxv,dyv,f,inf3,ny)
  use iso_c_binding, only: fp => c_double
  !
  !  input:
  !============
  implicit none
  integer ny,iaval,i,j
  !============
  real(fp) :: dx,dy
  !============
  integer ict(6)                    ! selector:
  !        ict(1)=1 for f      (don't evaluate f if ict(1)=0)
  !        ict(2)=1 for df/dx   ""
  !        ict(3)=1 for df/dy   ""
  !        ict(4)=1 for d2f/dx2
  !        ict(5)=1 for d2f/dy2
  !        ict(6)=1 for d2f/dxdy
  !
  !    note:  if ict(1)=-1, evaluate f,d2f/dx2,d2f/dy2,d4f/dx2dy2
  !
  !                   the number of non zero values ict(1:6)
  !                   determines the number of outputs...
  !
  !  new dmc December 2005 -- access to higher derivatives (even if not
  !  continuous-- but can only go up to 3rd derivatives on any one coordinate.
  !     if ict(1)=3 -- want 3rd derivatives
  !          ict(2)=1 for d3f/dx3
  !          ict(3)=1 for d3f/dx2dy
  !          ict(4)=1 for d3f/dxdy2
  !          ict(5)=1 for d3f/dy3
  !               number of non-zero values ict(2:5) gives no. of outputs
  !     if ict(1)=4 -- want 4th derivatives
  !          ict(2)=1 for d4f/dx3dy
  !          ict(3)=1 for d4f/dx2dy2
  !          ict(4)=1 for d4f/dxdy3
  !               number of non-zero values ict(2:4) gives no. of outputs
  !     if ict(1)=5 -- want 5th derivatives
  !          ict(2)=1 for d5f/dx3dy2
  !          ict(3)=1 for d5f/dx2dy3
  !               number of non-zero values ict(2:3) gives no. of outputs
  !     if ict(1)=6 -- want 6th derivatives
  !          d6f/dx3dy3 -- one value is returned.
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
  integer iv(ivec),jv(ivec)         ! grid cell indices -- vectors
  real(fp) :: dxv(ivec),dyv(ivec)          ! displacements w/in cell -- vectors
  !
  integer inf3                      ! 3rd dimension of f -- .ge. nx
  real(fp) :: f(4,4,inf3,ny)               ! bicubic fcn spline coeffs array
  !
  !  usage example:
  !
  !  1.  for each element (xx(v),yy(v)) in a vector of (x,y) pairs,
  !    find the x and y zone indices and displacements with respect
  !    to the "lower left corner" of the zone; store these in vectors
  !    iv,jv and dxv,dyv.
  !
  !  2.  set ict(1)=0, ict(2)=1, ict(3)=1, the rest zero -- get only
  !      the 1st derivatives.
  !
  !  3.  ivec is the length of the vector; ivd is the 1st dimension
  !      of the array fval to receive the output
  !
  !      real fval(ivd,6)
  !      real xv(ivd),yv(ivd)
  !      integer iv(ivd),jv(ivd)
  !      real dxv(ivd),dyv(ivd)
  !      integer ict(6)
  !
  !      real fspline(4,4,nx,ny)  ! spline coeffs
  !      data ict/0,1,1,0,0,0/    ! this call:  want 1st derivatives
  !                               ! only ... these will be output to
  !                               ! fval(*,1) fval(*,2)
  !      ...
  !      do iv=1,ivec
  !        ...                    ! find indices and displacements
  !      end do
  !      call bcspevfn(ict,ivec,ivd,fval,iv,jv,dxv,dyv,fspline,nx,ny)
  !
  !-------------------------------------------------------------------
  !  local:
  !
  integer v                         ! vector element index
  !
  !  OK can now do evaluations
  !
  iaval=0  ! fval addressing
  !
  if(ict(1).le.2) then
     if((ict(1).gt.0).or.(ict(1).eq.-1)) then
        !  evaluate f
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                f(1,1,i,j)+dy*(f(1,2,i,j)+dy*(f(1,3,i,j)+dy*f(1,4,i,j))) &
                +dx*(f(2,1,i,j)+dy*(f(2,2,i,j)+dy*(f(2,3,i,j)+dy*f(2,4,i,j))) &
                +dx*(f(3,1,i,j)+dy*(f(3,2,i,j)+dy*(f(3,3,i,j)+dy*f(3,4,i,j))) &
                +dx*(f(4,1,i,j)+dy*(f(4,2,i,j)+dy*(f(4,3,i,j)+dy*f(4,4,i,j))))))
        end do
     end if
     !
     if((ict(2).gt.0).and.(ict(1).ne.-1)) then
        !  evaluate df/dx
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                f(2,1,i,j)+dy*(f(2,2,i,j)+dy*(f(2,3,i,j)+dy*f(2,4,i,j))) &
                +2.0_fp*dx*( &
                f(3,1,i,j)+dy*(f(3,2,i,j)+dy*(f(3,3,i,j)+dy*f(3,4,i,j))) &
                +1.5_fp*dx*( &
                f(4,1,i,j)+dy*(f(4,2,i,j)+dy*(f(4,3,i,j)+dy*f(4,4,i,j))) &
                ))
        end do
     end if
     !
     if((ict(3).gt.0).and.(ict(1).ne.-1)) then
        !  evaluate df/dy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                f(1,2,i,j)+dy*(2.0_fp*f(1,3,i,j)+dy*3.0_fp*f(1,4,i,j)) &
                +dx*(f(2,2,i,j)+dy*(2.0_fp*f(2,3,i,j)+dy*3.0_fp*f(2,4,i,j)) &
                +dx*(f(3,2,i,j)+dy*(2.0_fp*f(3,3,i,j)+dy*3.0_fp*f(3,4,i,j)) &
                +dx*(f(4,2,i,j)+dy*(2.0_fp*f(4,3,i,j)+dy*3.0_fp*f(4,4,i,j)) &
                )))
        end do
     end if
     !
     if((ict(4).gt.0).or.(ict(1).eq.-1)) then
        !  evaluate d2f/dx2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                2.0_fp*( &
                f(3,1,i,j)+dy*(f(3,2,i,j)+dy*(f(3,3,i,j)+dy*f(3,4,i,j)))) &
                +6.0_fp*dx*( &
                f(4,1,i,j)+dy*(f(4,2,i,j)+dy*(f(4,3,i,j)+dy*f(4,4,i,j))))
        end do
     end if
     !
     if((ict(5).gt.0).or.(ict(1).eq.-1)) then
        !  evaluate d2f/dy2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                2.0_fp*f(1,3,i,j)+6.0_fp*dy*f(1,4,i,j) &
                +dx*(2.0_fp*f(2,3,i,j)+6.0_fp*dy*f(2,4,i,j) &
                +dx*(2.0_fp*f(3,3,i,j)+6.0_fp*dy*f(3,4,i,j) &
                +dx*(2.0_fp*f(4,3,i,j)+6.0_fp*dy*f(4,4,i,j))))
        end do
     end if
     !
     if((ict(6).gt.0).and.(ict(1).ne.-1)) then
        !  evaluate d2f/dxdy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                f(2,2,i,j)+dy*(2.0_fp*f(2,3,i,j)+dy*3.0_fp*f(2,4,i,j)) &
                +2._fp*dx*(f(3,2,i,j)+dy*(2.0_fp*f(3,3,i,j)+dy*3.0_fp*f(3,4,i,j)) &
                +1.5_fp*dx*(f(4,2,i,j)+dy*(2.0_fp*f(4,3,i,j)+dy*3.0_fp*f(4,4,i,j)) &
                ))
        end do
     end if
     !
     if(ict(1).eq.-1) then
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                4.0_fp*f(3,3,i,j)+12.0_fp*dy*f(3,4,i,j) &
                +dx*(12.0_fp*f(4,3,i,j)+36.0_fp*dy*f(4,4,i,j))
        end do
     end if
     !
     !-----------------------------------
     !  access to 3rd derivatives
     !
  else if(ict(1).eq.3) then
     if(ict(2).eq.1) then
        !  evaluate d3f/dx3 (not continuous)
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                +6.0_fp*( &
                f(4,1,i,j)+dy*(f(4,2,i,j)+dy*(f(4,3,i,j)+dy*f(4,4,i,j))))
        end do
     end if
     !
     if(ict(3).eq.1) then
        !  evaluate d3f/dx2dy
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                2.0_fp*( &
                f(3,2,i,j)+dy*(2.0_fp*f(3,3,i,j)+dy*3.0_fp*f(3,4,i,j))) &
                +6.0_fp*dx*( &
                f(4,2,i,j)+dy*(2.0_fp*f(4,3,i,j)+dy*3.0_fp*f(4,4,i,j)))
        end do
     end if
     !
     if(ict(4).eq.1) then
        !  evaluate d3f/dxdy2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                (2.0_fp*f(2,3,i,j)+6.0_fp*dy*f(2,4,i,j) &
                +2.0_fp*dx*(2.0_fp*f(3,3,i,j)+6.0_fp*dy*f(3,4,i,j) &
                +1.5_fp*dx*(2.0_fp*f(4,3,i,j)+6.0_fp*dy*f(4,4,i,j)) &
                ))
        end do
     end if

     if(ict(5).eq.1) then
        !  evaluate d3f/dy3 (not continuous)
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           fval(v,iaval)=6.0_fp*(f(1,4,i,j)+ &
                dx*(f(2,4,i,j)+dx*(f(3,4,i,j)+dx*f(4,4,i,j))))
        end do
     end if
     !
     !-----------------------------------
     !  access to 4th derivatives
     !
  else if(ict(1).eq.4) then
     if(ict(2).eq.1) then
        !  evaluate d4f/dx3dy (not continuous)
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                +6.0_fp*( &
                f(4,2,i,j)+dy*2.0_fp*(f(4,3,i,j)+dy*1.5_fp*f(4,4,i,j)))
        end do
     end if
     !
     if(ict(3).eq.1) then
        !  evaluate d4f/dx2dy2
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                4.0_fp*f(3,3,i,j)+12.0_fp*dy*f(3,4,i,j) &
                +dx*(12.0_fp*f(4,3,i,j)+36.0_fp*dy*f(4,4,i,j))
        end do
     end if
     !
     if(ict(4).eq.1) then
        !  evaluate d4f/dxdy3 (not continuous)
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           fval(v,iaval)= &
                6.0_fp*(f(2,4,i,j) &
                +2.0_fp*dx*(f(3,4,i,j)+1.5_fp*dx*f(4,4,i,j)))
        end do
     end if
     !
     !-----------------------------------
     !  access to 5th derivatives
     !
  else if(ict(1).eq.5) then
     if(ict(2).eq.1) then
        !  evaluate d5f/dx3dy2 (not continuous)
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dy=dyv(v)
           fval(v,iaval)= &
                +12.0_fp*(f(4,3,i,j)+dy*3.0_fp*f(4,4,i,j))
        end do
     end if
     !
     if(ict(3).eq.1) then
        !  evaluate d5f/dx3dy2 (not continuous)
        iaval=iaval+1
        do v=1,ivec
           i=iv(v)
           j=jv(v)
           dx=dxv(v)
           fval(v,iaval)= &
                12.0_fp*(f(3,4,i,j)+dx*3.0_fp*f(4,4,i,j))
        end do
     end if
     !
     !-----------------------------------
     !  access to 6th derivatives
     !
  else if(ict(1).eq.6) then
     !  evaluate d6f/dx3dy3 (not continuous)
     iaval=iaval+1
     do v=1,ivec
        i=iv(v)
        j=jv(v)
        fval(v,iaval)= &
             36.0_fp*f(4,4,i,j)
     end do
  end if
  !
  return
end subroutine bcspevfn
!----------------------
