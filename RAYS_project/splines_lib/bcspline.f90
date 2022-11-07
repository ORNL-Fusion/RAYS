!  bcspline -- dmc 30 May 1996
!
!  set up coefficients for bicubic spline with following BC's:
!  FULL BC CONTROL at all bdys
!
!  inhomogeneous explicit BCs -- this means setting of 1st or 2nd
!  derivative at boundary to a non-zero value.
!
!  periodic, not-a-knot, zero derivative, and divided-difference based
!  BCs are "homogeneous"-- i.e. if splines s & t satisfy the BC then
!  the spline (c*s + t) formed as a linear combination of these two
!  splines, also satisfies the BC.
!
!  algorithm note -- handling of inhomogeneous explicit BC's while
!  maintaining full C2 differentiability is delicate.  Basic method:  use
!  a fully C2 method based on the "not-a-knot" BC, and then, correct to
!  meet each user BC by calculating a C2 spline that is zero at all grid
!  points but satisfies a BC which is the difference btw the user spec
!  and the not-a-knot result; add the coeffs of this into the original.
!
!  for this more workspace is needed: nwk .ge. 4*inx*inth +5*max(inx,inth)
!
subroutine bcspline(x,inx,th,inth,fspl,inf3, &
     ibcxmin,bcxmin,ibcxmax,bcxmax, &
     ibcthmin,bcthmin,ibcthmax,bcthmax, &
     wk,nwk,ilinx,ilinth,ier)
  use iso_c_binding, only: fp => c_double
  !
  !============
  implicit none
  integer inth,inf3,ibcxmin,ibcxmax,ibcthmin,ibcthmax,nwk,ilinx
  integer ilinth,ier,inx,iflg2,ix,itest,ierx,ierth,inxo,ith
  integer intho,ic,ibcthmina,ibcthmaxa,iasc,iinc,iawk,jx,jth,ii
  integer iadr,ia5w,iaspl
  !============
  real(fp) :: xo2,xo6,zhxn,zhth,zdiff1,zdiff2
  !============
  real(fp) :: x(inx),th(inth),fspl(4,4,inf3,inth),wk(nwk)
  real(fp) :: bcxmin(*),bcxmax(*)      ! (inth) if used
  real(fp) :: bcthmin(*),bcthmax(*)    ! (inx) if used
  !
  !  input:
  !    x(1...inx) -- abscissae, first dimension of data
  !   th(1...inth) -- abscissae, second dimension of data  f(x,th)
  !   fspl(1,1,1..inx,1..inth) -- function values
  !   inf3 -- fspl dimensioning, inf3.ge.inx required.
  !
  !  boundary conditions input:
  !   ibcxmin -- indicator for boundary condition at x(1):
  !    bcxmin(...) -- boundary condition data
  !     =-1 -- periodic boundary condition
  !     =0 -- use "not a knot", bcxmin(...) ignored
  !     =1 -- match slope, specified at x(1),th(ith) by bcxmin(ith)
  !     =2 -- match 2nd derivative, specified at x(1),th(ith) by bcxmin(ith)
  !     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all th(j)
  !     =4 -- boundary condition is d2f/dx2=0 at x(1), all th(j)
  !     =5 -- match 1st derivative df/dx to 1st divided difference
  !     =6 -- match 2nd derivative d2f/dx2 to 2nd divided difference
  !     =7 -- match 3rd derivative d3f/dx3 3rd divided difference
  !           (for more detailed definition of BCs 5-7, see the
  !           comments of subroutine mkspline)
  !   NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
  !
  !   ibcxmax -- indicator for boundary condition at x(nx):
  !    bcxmax(...) -- boundary condition data
  !     (interpretation as with ibcxmin, bcxmin)
  !   NOTE:  if ibcxmin=-1, ibcxmax is ignored! ...and the BC is periodic.
  !
  !   ibcthmin -- indicator for boundary condition at th(1):
  !    bcthmin(...) -- boundary condition data
  !     (interpretation as with ibcxmin, bcxmin)
  !   ibcthmax -- indicator for boundary condition at th(inth):
  !    bcthmax(...) -- boundary condition data
  !     (interpretation as with ibcxmin, bcxmin)
  !   NOTE:  if ibcthmin=-1, ibcthmax is ignored! ...and the BC is periodic.
  !
  !   NOTE the bcxmin,bcxmax,bcthmin,bcthmax arrays are only used if the
  !     corresponding boundary condition flags are set to 1 or 2.
  !     Carefully note the dimensioning of these arrays!
  !
  !  output:
  !   fspl(*,*,1..inx,1..inth) -- bicubic spline coeffs (4x4)
  !   ...fspl(1,1,*,*) is not replaced.
  !
  !   ilinx -- =1 on output if x(inx) pts are nearly evenly spaced (tol=1e-3)
  !   ilinth-- =1 on output if th(inth) evenly spaced (tol=1e-3)
  !
  !   ier -- completion code, 0 for normal
  !
  !  workspace:
  !   wk -- must be at least 5*max(inx,inth) large
  !                          5*max(inx,inth) + 4*inx*inth large
  !                          if explicit non-zero d/dth or d2/dth2 BC's
  !                          are supplied.
  !  nwk -- size of workspace of workspace provided
  !
  !---------------------------------
  !  in what follows, "f" is an abbreviation for "fspl"...
  !
  !  compute bicubic spline of 2d function, given values at the grid
  !  grid crossing points, f(1,1,i,j)=f(x(i),th(j)).
  !
  !  on evaluation:  for point x btw x(i) and x(i+1), dx=x-x(i)
  !                       and th btw th(j) and th(j+1), dt=th-th(j),
  !
  !      spline =
  !        f(1,1,i,j) + dx*f(2,1,i,j) + dx**2*f(3,1,i,j) + dx**3*f(4,1,i,j)
  !   +dt*(f(1,2,i,j) + dx*f(2,2,i,j) + dx**2*f(3,2,i,j) + dx**3*f(4,2,i,j))
  !   +d2*(f(1,3,i,j) + dx*f(2,3,i,j) + dx**2*f(3,3,i,j) + dx**3*f(4,3,i,j))
  !   +d3*(f(1,4,i,j) + dx*f(2,4,i,j) + dx**2*f(3,4,i,j) + dx**3*f(4,4,i,j))
  !
  !      where d2=dt**2 and d3=dt**3.
  !
  integer iselect1(10)
  integer iselect2(10)
  !
  real(fp) :: zcur(1)
  !---------------------------------
  !
  !  see if 2nd pass is needed due to "non-linear" d/dth bdy cond.
  !
  iflg2=0
  if(ibcthmin.ne.-1) then
     if((ibcthmin.eq.1).or.(ibcthmin.eq.2)) then
        do ix=1,inx
           if (bcthmin(ix).ne.0.0_fp) iflg2=1
        end do
     end if
     if((ibcthmax.eq.1).or.(ibcthmax.eq.2)) then
        do ix=1,inx
           if (bcthmax(ix).ne.0.0_fp) iflg2=1
        end do
     end if
  end if
  !
  ier=0
  itest=5*max(inx,inth)
  if(iflg2.eq.1) then
     itest=itest +4*inx*inth
  end if
  !
  if(nwk.lt.itest) then
     write(6,9901) nwk,itest
9901 format(' ?bcspline:  workspace too small:'/ &
          '  user supplied:  nwk=',i6,'; need at least:  ',i6/ &
          '  nwk=4*nx*ny +5*max(nx,ny) will work for any user'/ &
          '  choice of bdy conditions.')
     ier=1
  end if
  if(inx.lt.2) then
     write(6,'('' ?bcspline:  at least 2 x points required.'')')
     ier=1
  end if
  if(inth.lt.2) then
     write(6,'('' ?bcspline:  need at least 2 theta points.'')')
     ier=1
  end if
  !
  call ibc_ck(ibcxmin,'bcspline','xmin',-1,7,ier)
  if(ibcxmin.ge.0) call ibc_ck(ibcxmax,'bcspline','xmax',0,7,ier)
  call ibc_ck(ibcthmin,'bcspline','thmin',-1,7,ier)
  if(ibcthmin.ge.0) call ibc_ck(ibcthmax,'bcspline','thmax',0,7,ier)
  !
  !  check ilinx & x vector
  !
  call splinck(x,inx,ilinx,1.0E-3_fp,ierx)
  if(ierx.ne.0) ier=2
  !
  if(ier.eq.2) then
     write(6,'('' ?bcspline:  x axis not strict ascending'')')
  end if
  !
  !  check ilinth & th vector
  !
  call splinck(th,inth,ilinth,1.0E-3_fp,ierth)
  if(ierth.ne.0) ier=3
  !
  if(ier.eq.3) then
     write(6,'('' ?bcspline:  th axis not strict ascending'')')
  end if
  !
  if(ier.ne.0) return
  !
  !------------------------------------
  !
  xo2=0.5_fp
  xo6=1.0_fp/6.0_fp
  !
  !  spline in x direction
  !
  inxo=4*(inx-1)
  do ith=1,inth
     !
     !  copy the function in
     !
     do ix=1,inx
        wk(4*(ix-1)+1)=fspl(1,1,ix,ith)
     end do
     !
     if(ibcxmin.eq.1) then
        wk(2)=bcxmin(ith)
     else if(ibcxmin.eq.2) then
        wk(3)=bcxmin(ith)
     end if
     !
     if(ibcxmax.eq.1) then
        wk(inxo+2)=bcxmax(ith)
     else if(ibcxmax.eq.2) then
        wk(inxo+3)=bcxmax(ith)
     end if
     !
     !  use Wayne's routine
     !
     call v_spline(ibcxmin,ibcxmax,inx,x,wk,wk(4*inx+1))
     !
     !  copy the coefficients out
     !
     do ix=1,inx
        fspl(2,1,ix,ith)=wk(4*(ix-1)+2)
        fspl(3,1,ix,ith)=wk(4*(ix-1)+3)*xo2
        fspl(4,1,ix,ith)=wk(4*(ix-1)+4)*xo6
     end do
     !
  end do
  !
  !-----------------------------------
  !
  !  spline in theta direction
  !
  intho=4*(inth-1)
  do ix=1,inx
     !
     !  spline each x coeff
     !
     do ic=1,4
        !
        !  copy ordinates in
        !
        do ith=1,inth
           wk(4*(ith-1)+1)=fspl(ic,1,ix,ith)
        end do
        !
        !  first pass:  use a linear BC -- if flag indicates BC correction
        !  will be needed, it will be done later
        !
        wk(2)=0.0_fp
        wk(3)=0.0_fp
        wk(intho+2)=0.0_fp
        wk(intho+3)=0.0_fp
        !
        ibcthmina=ibcthmin
        ibcthmaxa=ibcthmax
        if(iflg2.eq.1) then
           if((ibcthmin.eq.1).or.(ibcthmin.eq.2)) ibcthmina=0
           if((ibcthmax.eq.1).or.(ibcthmax.eq.2)) ibcthmaxa=0
        end if
        !
        call v_spline(ibcthmina,ibcthmaxa,inth,th,wk,wk(4*inth+1))
        !
        !  copy coeffs out
        !
        do ith=1,inth
           fspl(ic,2,ix,ith)=wk(4*(ith-1)+2)
           fspl(ic,3,ix,ith)=wk(4*(ith-1)+3)*xo2
           fspl(ic,4,ix,ith)=wk(4*(ith-1)+4)*xo6
        end do
        !
     end do
     !
  end do
  !
  !  now make correction for user BC's if needed
  !
  if(iflg2.eq.1) then
     !
     iasc=1                         ! wk addr for correction splines
     iinc=4*inth                    ! spacing btw correction splines
     iawk=iasc+4*inth*inx
     !
     !  last grid zone widths
     !
     zhxn=x(inx)-x(inx-1)
     jx=inx-1
     zhth=th(inth)-th(inth-1)
     jth=inth-1
     !
     do ii=1,10
        iselect1(ii)=0
        iselect2(ii)=0
     end do
     if(ibcthmin.eq.1) iselect1(3)=1
     if(ibcthmin.eq.2) iselect1(5)=1
     if(ibcthmax.eq.1) iselect2(3)=1
     if(ibcthmax.eq.2) iselect2(5)=1
     !
     !  loop over BC's
     !
     do ix=1,inx
        !
        !  (a) d/dth @ th(1) difference btw current BC and user request
        !
        if(ibcthmin.eq.1) then
           if(ix.lt.inx) then
              zcur(1)=fspl(1,2,ix,1)   ! 1st deriv.
           else
              zcur(1)=fspl(1,2,jx,1)+zhxn*(fspl(2,2,jx,1)+zhxn* &
                   (fspl(3,2,jx,1)+zhxn*fspl(4,2,jx,1)))
           end if
           zdiff1=bcthmin(ix)-zcur(1)
        else if(ibcthmin.eq.2) then
           if(ix.lt.inx) then
              zcur(1)=2.0_fp*fspl(1,3,ix,1) ! 2nd deriv.
           else
              zcur(1)=2.0_fp* &
                   (fspl(1,3,jx,1)+zhxn*(fspl(2,3,jx,1)+zhxn* &
                   (fspl(3,3,jx,1)+zhxn*fspl(4,3,jx,1))))
           end if
           zdiff1=bcthmin(ix)-zcur(1)
        else
           zdiff1=0.0_fp
        end if
        !
        !  (b) d/dth @ th(inth) difference btw current BC and user request
        !
        if(ibcthmax.eq.1) then
           if(ix.lt.inx) then
              !  1st deriv.
              zcur(1)= &
                   fspl(1,2,ix,jth)+zhth*(2.0_fp*fspl(1,3,ix,jth)+ &
                   zhth*3.0_fp*fspl(1,4,ix,jth))
           else
              call bcspeval(x(inx),th(inth),iselect2,  zcur, &
                   x,inx,th,inth,ilinx,ilinth,fspl,inf3,ier)
              if(ier.ne.0) return
           end if
           zdiff2=bcthmax(ix)-zcur(1)
        else if(ibcthmax.eq.2) then
           if(ix.lt.inx) then
              !  2nd deriv.
              zcur(1)=2.0_fp*fspl(1,3,ix,jth)+ &
                   6.0_fp*zhth*fspl(1,4,ix,jth)
           else
              call bcspeval(x(inx),th(inth),iselect2,  zcur(1), &
                   x,inx,th,inth,ilinx,ilinth,fspl,inf3,ier)
              if(ier.ne.0) return
           end if
           zdiff2=bcthmax(ix)-zcur(1)
        else
           zdiff2=0.0_fp
        end if
        !
        !  ok compute the theta spline with BC's to span the difference(s)
        !  these theta "correction splines" are zero at all the grid points
        !  but have at least one non-zero 1st or 2nd derivative BC
        !
        iadr=iasc+(ix-1)*iinc
        do ith=1,inth
           wk(iadr+4*(ith-1))=0.0_fp
        end do
        !
        wk(iadr+1)=0.0_fp
        wk(iadr+2)=0.0_fp
        wk(iadr+intho+1)=0.0_fp
        wk(iadr+intho+2)=0.0_fp
        !
        if(ibcthmin.eq.1) then
           wk(iadr+1)=zdiff1
        else if(ibcthmin.eq.2) then
           wk(iadr+2)=zdiff1
        end if
        !
        if(ibcthmax.eq.1) then
           wk(iadr+intho+1)=zdiff2
        else if(ibcthmax.eq.2) then
           wk(iadr+intho+2)=zdiff2
        end if
        !
        call v_spline(ibcthmin,ibcthmax,inth,th,wk(iadr),wk(iawk))
     end do
     !
     !  add in results to main array -- th spline coef corrections
     !
     do ix=1,inx
        iadr=iasc+(ix-1)*iinc
        do ith=1,inth-1
           wk(iadr+4*(ith-1)+2)=wk(iadr+4*(ith-1)+2)*xo2
           wk(iadr+4*(ith-1)+3)=wk(iadr+4*(ith-1)+3)*xo6
           if(ix.lt.inx) then
              fspl(1,2,ix,ith)=fspl(1,2,ix,ith)+wk(iadr+4*(ith-1)+1)
              fspl(1,3,ix,ith)=fspl(1,3,ix,ith)+wk(iadr+4*(ith-1)+2)
              fspl(1,4,ix,ith)=fspl(1,4,ix,ith)+wk(iadr+4*(ith-1)+3)
           end if
        end do
     end do
     !
     !  compute the x splines of the th spline correction coeffs
     !
     ia5w=iawk+4*inx
     !
     do ith=1,inth-1
        do ic=2,4
           do ix=1,inx
              iaspl=iasc+iinc*(ix-1)
              wk(iawk+4*(ix-1))=wk(iaspl+4*(ith-1)+(ic-1))
           end do
           !
           !  use zero BCs for this correction spline
           !
           wk(iawk+1)=0.0_fp
           wk(iawk+2)=0.0_fp
           wk(iawk+inxo+1)=0.0_fp
           wk(iawk+inxo+2)=0.0_fp
           !
           !  periodic spline of correction spline higher coeffs (1st coeffs are
           !  all zero by defn of the correction spline
           !
           call v_spline(ibcxmin,ibcxmax,inx,x,wk(iawk),wk(ia5w))
           !
           do ix=1,inx-1
              fspl(2,ic,ix,ith)=fspl(2,ic,ix,ith)+ &
                   wk(iawk+4*(ix-1)+1)
              fspl(3,ic,ix,ith)=fspl(3,ic,ix,ith)+ &
                   wk(iawk+4*(ix-1)+2)*xo2
              fspl(4,ic,ix,ith)=fspl(4,ic,ix,ith)+ &
                   wk(iawk+4*(ix-1)+3)*xo6
           end do
           !
        end do
     end do                          ! ith
     !
  end if                             ! BC correction needs test
  !
  return
end subroutine bcspline
