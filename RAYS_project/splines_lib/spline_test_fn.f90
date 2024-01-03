! simple function f = rho(x,y) and derivatives. N.B. Derivatives singular at origin.
    subroutine spline_test_fn(x, y, f, fx, fy, fxx, fxy, fyy)

    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind), intent(in) :: x, y
    real(KIND=rkind), intent(out) ::  f, fx, fy, fxx, fxy, fyy

    f = sqrt(x**2 + y**2)

    fx = 0.; fy = 0.; fxx = 0.; fxy = 0.; fyy = 0.;

    if (f > 0. ) then
		fx = x/f
		fy = y/f
		fxx = 1./f - x**2/f**3
		fxy = -x*y/f**3
		fyy = 1./f - y**2/f**3
    end if

    return
    end subroutine spline_test_fn

!****************************************************************************************
! simple real function f(x) and derivatives.
    subroutine spline_test_fn_1D(x, f, fx, fxx)

    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind), intent(in) :: x
    real(KIND=rkind), intent(out) ::  f, fx, fxx

    real(KIND=rkind) ::  k = 10.

    f = 1.0_rkind + cos(k*x)
    fx = -k*sin(k*x)
    fxx = -k**2*cos(k*x)

    return
    end subroutine spline_test_fn_1D

!****************************************************************************************
! simple complex function f(x) and derivatives.
    subroutine spline_test_fn_c(x, f, fx, fxx)

    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind), intent(in) :: x
    complex(KIND=rkind), intent(out) ::  f, fx, fxx

    real(KIND=rkind) ::  k = 10.

    f = 1. + exp(cmplx(0.,k*x))
    f = 1.0_rkind + cmplx(cos(k*x), sin(k*x), kind=rkind)
    fx = cmplx(0._rkind,1._rkind)*k*exp(cmplx(0._rkind,k*x, kind=rkind))
    fxx = -k**2*exp(cmplx(0._rkind,k*x, kind=rkind))

    return
    end subroutine spline_test_fn_C

!****************************************************************************************
! simple Plasma dispersion function function f(x) and derivatives.
    subroutine spline_test_fn_Zf(x, f, fx, fxx)

    use zfunctions_m, only : zfun, zfun0

    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind), intent(in) :: x
    complex(KIND=rkind), intent(out) ::  f, fx, fxx

    f = zfun0(cmplx(x, 0.0_rkind, kind=rkind), 1.0_rkind)
    fx = 0.
    fxx = 0.

    return
    end subroutine spline_test_fn_Zf

