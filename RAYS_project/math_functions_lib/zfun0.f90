 complex function zfun0(z, kz)
!   calculates plasma function defined by Eq.(8-82), Stix.
! N.B. This assumes that kz is real and non-zero

    implicit none

    complex :: z
    real :: kz

    complex, external :: DISP

!   See Eq.(8-82).
    if ( kz > 0. ) then
          zfun0 = DISP(z)
    else if ( kz < 0. ) then
          zfun0 = -disp(-z)
    else
       write(0, *) 'ZFUN0: Error, kz must be real and non-zero, kz= ', kz
       stop 1
    end if
    return
 end function zfun0
