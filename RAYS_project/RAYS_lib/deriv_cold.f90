 subroutine deriv_cold(eq, nvec, dddx, dddk, dddw)
!   calculates the derivatives of D with respect to k, r, omega.
!   v(1:3) = (x,y,z); v(4:6) = (kx, ky, kz),
!   dddx = dD/dx, dddk = dD/dk, dddw = dD/d(omega).

    use constants_m, only : rkind
    use equilibrium_m, only : eq_point, write_eq_point
    use rf_m, only : omgrf, k0
    use species_m, only : nspec0, nspec
    use diagnostics_m, only : message_unit, verbosity

    implicit none

    type(eq_point), intent(in) :: eq
    real(KIND=rkind), intent(in) :: nvec(3)
    real(KIND=rkind), intent(out) :: dddx(3), dddk(3), dddw

    real(KIND=rkind) :: alpha(0:nspec), gamma(0:nspec)
    real(KIND=rkind) :: n1, n3
    real(KIND=rkind) :: p, t, dtdg(0:nspec)
    real(KIND=rkind) :: q,  dqda(0:nspec),  dqdg(0:nspec)
    real(KIND=rkind) :: q1, dq1da(0:nspec), dq1dg(0:nspec)
    real(KIND=rkind) :: q2, dq2da(0:nspec), dq2dg(0:nspec)
    real(KIND=rkind) :: u,  duda(0:nspec),  dudg(0:nspec)

    real(KIND=rkind) :: dddn12, dddn3, ddda(0:nspec), dddg(0:nspec)
    real(KIND=rkind) :: dn3dk(3), dn12dk(3)
    real(KIND=rkind) :: dn3dx(3), dn12dx(3), dadx(3,0:nspec), dgdx(3,0:nspec)
    real(KIND=rkind) :: dn3dw, dn12dw, dadw(0:nspec), dgdw(0:nspec)

    real(KIND=rkind), dimension(0:nspec,0:nspec) :: gp, gm, gpm

    integer :: is, is1, is2, ivec

	do is = 0, nspec
		alpha(is) = eq%alpha(is)
		gamma(is) = eq%gamma(is)
	end do

    n3 = dot_product(nvec, eq%bunit)
    n1 = sqrt( sum((nvec-n3*eq%bunit)**2) )

!   Derivatives with respect to k.
!   dn3dk = d(n3)/dk; dn12dk = d(n1^2)/dk
    dn3dk = eq%bunit / k0
    dn12dk = (2./k0) * (nvec-n3*eq%bunit)

!   Derivatives with respect to space coordinates.
    do ivec = 1, 3
       dn3dx(ivec) = sum( eq%gradbunit(ivec,:) * nvec )
    end do
    dn12dx = -2.*n3*dn3dx

    dadx = 0.
    dgdx = 0.
    do ivec = 1, 3
        do is = 0, nspec
       if (verbosity > 4) call write_eq_point(eq)
           dadx(ivec,is) = eq%alpha(is) * eq%gradns(ivec,is)/eq%ns(is)
           dgdx(ivec,is) = gamma(is) * eq%gradbmag(ivec)/eq%bmag
        end do
    end do

!   Derivatives with respect to omega.
!   dn3dw = d(n3)/d(omega); dn12dw = d(n1^2)/d(omega)
!   dadw = d(alpha)/d(omega); dgdw = d(gamma)/d(omega)
    dn3dw = -n3/omgrf
    dn12dw = (-2./omgrf) * n1**2
    dadw = -2./omgrf * alpha
    dgdw = -1./omgrf * gamma

!   Over species.
    p = 1. - sum(alpha)
    t = product(1.-gamma**2)

!   Derivatives with respect to alpha.
!   dQ1/d(alpha) & dQ2/d(alpha):
    dq1da = 1; dq2da = 1
    do is1 = 0, nspec
       do is = 0, nspec
          if ( is /= is1 ) then
             dq1da(is1) = dq1da(is1) * (1.+gamma(is))
             dq2da(is1) = dq2da(is1) * (1.-gamma(is))
          end if
       end do
    end do

!   Q1 & Q2:
    q1 = sum(alpha * dq1da)
    q2 = sum(alpha * dq2da)

!   U = TS = T - ...
    u = t - sum(alpha*dq1da*dq2da)

!   Q = 2U - T + Q1*Q2.
    q = 2.*u - t + q1*q2

!   dU/d(alpha):
    duda = -dq1da * dq2da

!   dQ/d(alpha):
    dqda = 2.*duda + dq1da*q2 + q1*dq2da

!   dD/d(alpha):
    ddda = -t * n3**4                                                &
      & + ( 2.*(u-p*duda) + (-t+duda)*n1**2 ) * n3**2                &
      & - q + p*dqda - (dqda-u+p*duda)*n1**2 + duda*n1**4

!   Derivatives with respect to gamma.
!   Kernel gp, gm, and gpm
    gp = 1.; gm = 1.
    do is1 = 0, nspec; do is2 = 0, nspec
       do is = 0, nspec
          if ( is /= is1 .and. is /= is2 ) then
             gp(is1,is2) = gp(is1,is2) * (1.+gamma(is))
             gm(is1,is2) = gm(is1,is2) * (1.-gamma(is))
          end if
       end do
    end do; end do
    gpm = gp * gm

!   dT/d(gamma):
    dtdg = 2.*gamma*duda

!   dU/d(gamma):
    do is = 0, nspec
       dudg(is) = sum(alpha*gpm(:,is))
    end do
    dudg = dtdg + 2.*gamma*(dudg+alpha*duda)

!   dQ1/d(gamma):
    do is = 0, nspec
       dq1dg(is) = sum(alpha*gp(:,is))
    end do
    dq1dg = dq1dg - alpha*dq1da

!   dQ2/d(gamma):
    do is = 0, nspec
       dq2dg(is) = sum(alpha*gm(:,is))
    end do
    dq2dg = -dq2dg + alpha*dq2da

!   dQd(gamma):
    dqdg = 2.*dudg - dtdg + dq1dg*q2 + q1*dq2dg

!   dD/d(gamma):
    dddg = dtdg*p * n3**4                                              &
      & + ( -2.*p*dudg + (dtdg*p+dudg)*n1**2 ) * n3**2                 &
      & + p*dqdg - (dqdg+p*dudg)*n1**2 + dudg*n1**4

!   dddn3 = dD/d(n3) and dddn12 = dD/d(n1^2).
    dddn3 = ( 4.*t*p*n3**2 + 2.*(-2.*p*u+(t*p+u)*n1**2) ) * n3
    dddn12= (t*p+u)*n3**2 - (q+p*u) + 2.*u*n1**2

!   Finally, dD/dk, dD/dx, dD/d(omega)
!   dddk = dD/dk:
    dddk = dddn3 * dn3dk + dddn12 * dn12dk

!   dddx = dD/dx:
    do ivec = 1, 3
       dddx(ivec) = sum( ddda*dadx(ivec,:) + dddg*dgdx(ivec,:) )
    end do
    dddx = dddx + dddn3 * dn3dx + dddn12 * dn12dx

!   dddw = dD/d(omega):
    dddw = sum(ddda*dadw + dddg*dgdw) + dddn3 * dn3dw + dddn12 * dn12dw

    verb: if (verbosity > 4) then
      write(*,*) ''
      write (*,*) 'dn3dk = ',  dn3dk
      write (*,*) 'dn12dk = ',  dn12dk
      write (*,*) 'dn3dx = ',  dn3dx
      write (*,*) 'dadx = ',  dadx
      write (*,*) 'dgdx = ',  dgdx

       write (*,*) 'dn3dw = ',  dn3dw
       write (*,*) 'dn12dw = ',  dn12dw
       write (*,*) 'dn12dw = ',  dn12dw
       write (*,*) 'dn12dw = ',  dn12dw

       write (*,*) 'p = ',  p
       write (*,*) 't = ',  t

       write (*,*) 'dq1da = ',  dq1da
       write (*,*) 'dq2da = ',  dq2da

       write (*,*) 'q1 = ',  q1
       write (*,*) 'q2 = ',  q2
       write (*,*) 'u = ',  u
       write (*,*) 'q = ',  q

       write (*,*) 'duda = ',  duda
       write (*,*) 'dqda = ',  dqda
       write (*,*) 'ddda = ',  ddda
       write (*,*) 'gp = ',  gp
       write (*,*) 'gm = ',  gm
       write (*,*) 'gpm = ',  gpm
       write (*,*) 'dtdg = ',  dtdg
       write (*,*) 'dudg = ',  dudg
       write (*,*) 'dq1dg = ',  dq1dg
       write (*,*) 'dq2dg = ',  dq2dg
       write (*,*) 'dqdg = ',  dqdg

       write(*,'(a,1p3e12.4)') 'dddn3 =', dddn3
       write(*,'(a,1p3e12.4)') 'dn3dw =', dn3dw
       write(*,'(a,1p3e12.4)') 'dddn12 =', dddn12
       write(*,'(a,1p3e12.4)') 'dn12dw =', dn12dw
       write(*,'(a,1p3e12.4)') 'ddda =', ddda
       write(*,'(a,1p3e12.4)') 'dddg =', dddg
       write(*,'(a,1p3e12.4)') 'dadw =', dadw
       write(*,'(a,1p3e12.4)') 'dgdw =', dgdw

       write(*,'(a,1p3e12.4)') 'dn3dk =', dn3dk
       write(*,'(a,1p3e12.4)') 'dn12dk =', dn12dk

      write(*,'(a,1p3e12.4)') 'dddx =', dddx
      write(*,'(a,1p3e12.4)') 'dddk =', dddk
      write(*,'(a,1p1e12.4)') 'dddw =', dddw
      write(*,*) ''

    end if verb
    return
 end subroutine deriv_cold
