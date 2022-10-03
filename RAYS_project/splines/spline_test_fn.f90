! simple function f = rho(x,y) and derivatives. N.B. Derivatives singular at origin.
    subroutine spline_test_fn(x, y, f, fx, fy, fxx, fxy, fyy)
    
    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind), intent(in) :: x, y
    real(KIND=rkind), intent(out) ::  f, fx, fy, fxx, fxy, fyy
    
    f = sqrt(x**2 + y**2)
    fx = x/f
    fy = y/f
    fxx = 1./f - x**2/f**3
    fxy = -x*y/f**3
    fyy = 1./f - y**2/f**3
    
    return
    end subroutine spline_test_fn
    
! simple function f(x) and derivatives.
    subroutine spline_test_fn_1D(x, f, fx, fxx)
    
    implicit none

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind), intent(in) :: x
    real(KIND=rkind), intent(out) ::  f, fx, fxx
    
    f = 1. + sin(0.3*x)
    fx = 0.3*cos(0.3*x)
    fxx = -0.3**2*sin(0.3*x)
    
    return
    end subroutine spline_test_fn_1D
    
