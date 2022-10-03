subroutine ibc_ck(ibc,slbl,xlbl,imin,imax,ier)
  use iso_c_binding, only: fp => c_double
  implicit none
  ! Check that spline routine ibc flag is in range
  integer, intent(in)  :: ibc          ! flag value
  character(len=*), intent(in) :: slbl ! subroutine name
  character(len=*), intent(in) :: xlbl ! axis label
  integer, intent(in)  :: imin         ! min allowed value
  integer, intent(in)  :: imax         ! max allowed value
  integer, intent(out) :: ier          ! output -- set =1 if error detected

  if((ibc.lt.imin).or.(ibc.gt.imax)) then
     ier=1
     write(6,1001) slbl,xlbl,ibc,imin,imax
1001 format(' ?',a,' -- ibc',a,' = ',i9,' out of range ',i2,' to ',i2)
  end if

  return
end subroutine ibc_ck
