 module constants_m
!   contains useful contants. Note that units are MKS.

    implicit none

    integer, parameter :: input_unit = 70 ! Unit number for namelist input file
    integer, parameter :: output_unit = 75 ! Unit number for formatted ray data output
    integer, parameter :: ray_list_unit = 74 ! Unit number for formatted ray list output
    complex :: i                ! positive square root of -1
    real :: pi                  ! usual pi= 3.14159.....
    real :: clight, mu0, eps0   ! consistent electro-magnetic constants
    real :: me, mp, e           ! electron mass, proton mass and unit charge
    real :: epsmach             ! measure of machine precision
    real :: tinyr               ! measure of smallest number represented
    real :: unitmat2(2,2)   ! 2X2 unit matrix
    real :: unitmat3(3,3)   ! 3X3 unit matrix
    real :: unlike      ! a real number unlikely set by input file
    real :: nlike       ! an integer number unlikely set by input file

 contains
    subroutine initialize_constants_m

       i = (0.,1.)
       pi = 3.1415926535897932385

       clight = 2.997930e8
       mu0 = pi * 4.e-7
       eps0=1./(mu0*clight**2)   ! chosen to make c=1/sqrt(eps0*mu0)

       me = 9.1094e-31
       mp = 1.6726e-27
       e = 1.6022e-19

       epsmach = epsilon(clight)
       tinyr = tiny(clight)

       unitmat2 = reshape( (/ 1., 0., 0., 1. /), (/ 2, 2 /) )
       unitmat3 = reshape( (/ 1., 0., 0., 0., 1., 0., 0., 0., 1. /),    &
          & (/ 3, 3 /) )

       unlike = - sqrt(1.5) * 1.e25
       nlike = -123454321 

    end subroutine initialize_constants_m

 end module constants_m
