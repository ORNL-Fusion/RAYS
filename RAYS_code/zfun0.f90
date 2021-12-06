 complex function zfun0(z, kz)
 
!   calculates plasma function defined by Eq.(8-82), Stix.

    implicit none

    complex :: z
    real :: kz

    integer :: izfunc
    complex :: fu
!    complex, external :: DISP

    interface disp
      complex function disp(z)
        complex z
      end function disp
    end interface disp
    
    izfunc = 1

!   See Eq.(8-82).
    if ( kz > 0. ) then
       if ( izfunc == 1 ) then
          zfun0 = disp(z)
       else if ( izfunc == 2 ) then
          call zfun(z, fu)
          zfun0 = fu
       else
          write(0,*) 'ZFUN0: izfunc=', izfunc
          stop 11
       end if

    else if ( kz < 0. ) then
       if ( izfunc == 1 ) then
          zfun0 = -disp(-z)
       else if ( izfunc == 2 ) then
          call zfun(z, fu)
          zfun0 = -fu
       else
          write(0,*) 'ZFUN0: izfunc=', izfunc
          stop 12
       end if

    else
       write(0, *) 'ZFUN0: kz=', kz
       stop 1
    end if

!   write (*,*) 'ZFUN0: zfun0=', zfun0
    return
    
 end function zfun0
 