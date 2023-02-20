 complex function zfun0(z, kz)
! Calculates plasma function defined by Eq.(8-82), Stix.
! N.B. This assumes that kz is real and non-zero

    implicit none

    complex :: z
    real :: kz

    complex, external :: zfun

!   See Eq.(8-82).
    if ( kz > 0. ) then
          zfun0 = zfun(z)
    else if ( kz < 0. ) then
          zfun0 = -zfun(-z)
    else
       write(0, *) 'ZFUN0: Error, kz must be real and non-zero, kz= ', kz
       stop 1
    end if
    return
 end function zfun0
