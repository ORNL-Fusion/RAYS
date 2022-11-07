subroutine zonfind(x,nx,zxget,i)
  use iso_c_binding, only: fp => c_double
  implicit none
  integer nx,nxm,i1,i2,ij,ii
  real(fp) :: dx
  real(fp) :: x(nx),zxget
  integer i
  !
  !  find index i such that x(i).le.zxget.le.x(i+1)
  !
  !  x(1...nx) is strict increasing and x(1).le.zxget.le.x(nx)
  !  (this is assumed to already have been checked -- no check here!)
  !
  nxm=nx-1
  if((i.lt.1).or.(i.gt.nxm)) then
     i1=1
     i2=nx-1
     go to 10
  end if

  if(x(i).gt.zxget) then
     ! look down
     dx=x(i+1)-x(i)
     if((x(i)-zxget).gt.4*dx) then
        i1=1
        i2=i-1
        go to 10
     else
        i2=i-1
        do ij=i2,1,-1
           if((x(ij).le.zxget).and.(zxget.le.x(ij+1))) then
              i=ij
              return
           end if
        end do
        i=1
        return
     end if
  else if(x(i+1).lt.zxget) then
     ! look up
     dx=x(i+1)-x(i)
     if((zxget-x(i+1)).gt.4*dx) then
        i1=i+1
        i2=nxm
        go to 10
     else
        i2=i+1
        do ij=i2,nxm
           if((x(ij).le.zxget).and.(zxget.le.x(ij+1))) then
              i=ij
              return
           end if
        end do
        ij=nxm
        return
     end if
  else
     ! already there...
     return
  end if

  !---------------------------
  !  binary search
  !
10 continue

  if(i1.eq.i2) then
     ! found by proc. of elimination
     i=i1
     return
  end if

  ii=(i1+i2)/2

  if(zxget.lt.x(ii)) then
     i2=ii-1
  else if(zxget.gt.x(ii+1)) then
     i1=ii+1
  else
     ! found
     i=ii
     return
  end if

  go to 10

  return
end subroutine zonfind

