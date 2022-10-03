subroutine splinck(x,inx,ilinx,ztol,ier)
  use iso_c_binding, only: fp => c_double
  !
  !  check if a grid is strictly ascending and if it is evenly spaced
  !  to w/in ztol
  !
  !============
  implicit none
  integer inx,ix
  !============
  real(fp) :: dxavg,zeps,zdiffx,zdiff
  !============
  real(fp) :: x(inx)                       ! input -- grid to check
  !
  integer ilinx                     ! output -- =1 if evenly spaced =2 O.W.
  !
  real(fp) :: ztol                         ! input -- spacing check tolerance
  !
  integer ier                       ! output -- =0 if OK
  !
  !  ier=1 is returned if x(1...inx) is NOT STRICTLY ASCENDING...
  !
  !-------------------------------
  !
  ier=0
  ilinx=1
  if(inx.le.1) return
  !
  dxavg=(x(inx)-x(1))/(inx-1)
  zeps=abs(ztol*dxavg)
  !
  do ix=2,inx
     zdiffx=(x(ix)-x(ix-1))
     if(zdiffx.le.0.0_fp) ier=2
     zdiff=zdiffx-dxavg
     if(abs(zdiff).gt.zeps) then
        ilinx=2
     end if
  end do
10 continue
  !
  return
end subroutine splinck
