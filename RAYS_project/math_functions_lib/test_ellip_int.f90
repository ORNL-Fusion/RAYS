PROGRAM test_ellip_int
   USE complete_elliptic_int_m, only : elliptic_Em, elliptic_Km

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
!    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

    REAL(kind = rkind) :: x, E, K, E_math, K_math

	x= 0.0_rkind
	E = elliptic_Em(x)
	K = elliptic_Km(x)
	E_math = 1.570796326794897_rkind
	K_math = 1.570796326794897_rkind

    write(*,*) ' '
    write(*,*) 'x = ', x
    WRITE (*, *) ' E  = ', E, '  E Mathematica = ', E_math
    WRITE (*, *) ' K  = ', K, '  K Mathematica = ', K_math

	x= 0.1_rkind
	E = elliptic_Em(x)
	K = elliptic_Km(x)
	E_math = 1.530757636897763_rkind
	K_math = 1.612441348720219_rkind

    write(*,*) ' '
    write(*,*) 'x = ', x
    WRITE (*, *) ' E  = ', E, '  E Mathematica = ', E_math
    WRITE (*, *) ' K  = ', K, '  K Mathematica = ', K_math

	x= 0.5_rkind
	E = elliptic_Em(x)
	K = elliptic_Km(x)
	E_math = 1.350643881047676_rkind
	K_math = 1.854074677301372_rkind

    write(*,*) ' '
    write(*,*) 'x = ', x
    WRITE (*, *) ' E  = ', E, '  E Mathematica = ', E_math
    WRITE (*, *) ' K  = ', K, '  K Mathematica = ', K_math

	x= 0.9_rkind
	E = elliptic_Em(x)
	K = elliptic_Km(x)
	E_math = 1.104774732704073_rkind
	K_math = 2.578092113348173_rkind

    write(*,*) ' '
    write(*,*) 'x = ', x
    WRITE (*, *) ' E  = ', E, '  E Mathematica = ', E_math
    WRITE (*, *) ' K  = ', K, '  K Mathematica = ', K_math

END PROGRAM test_ellip_int
