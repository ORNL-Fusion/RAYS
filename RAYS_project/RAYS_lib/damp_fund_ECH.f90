! A simple weak damping routine for fundamental electron
! cyclotron absorption based on notes from 4-24-81
! This routine neglects ion dynamics altogether
!
! Working notes
! (DBB 11-1-2022) Adapted to new RAYS_lib environment
!

    SUBROUTINE DAMP_FUND_ECH(eq, v_kx, vg, ksi, ki)

    use constants_m, only : clight, rkind
    use rf_m, only : omgrf, k0
    use species_m, only : nspec, qs, ms
    use ode_m, only : nv
    use equilibrium_m, only : eq_point, write_eq_point
    use zfunctions_m, only : zfun, zfun0, zfun0_real_arg

    implicit none

    type(eq_point), intent(in) :: eq
    real(KIND=rkind), intent(in) :: v_kx(6) ! N.B. need x and k from v(:), not the rest of it.
    real(KIND=rkind), intent(in) :: vg(3)
    real(KIND=rkind), intent(out) :: ksi(0:nspec), ki

    real(KIND=rkind) :: kvec(3), nvec(3), vth, vg_unit(3), expand_z0, xi

    complex(KIND=rkind) :: zf, zfun0, zfun0_real_arg

	REAL(KIND=rkind)  :: ALFAE, BETAE, VT
	REAL(KIND=rkind)  :: B1, P, Q, R1, R3, k1, k3, RDX, R1S, R3S, RS
    REAL(KIND=rkind)  :: LAMBDA1,LAMBDA2, LAMBDA5, A, B, C
	REAL(KIND=rkind)  :: DDNX2, DDNZ, DNparallel(3), DNperp2(3), DDN(3)

	COMPLEX D_WARM, DELTA

	ksi(0:nspec) = 0.
	ki = ksi(0)

!   kvec (nvec) = k (k/k0) in xyz coordinates (vector)
    kvec = v_kx(4:6)
    nvec = kvec/k0

! 1 and 3 below refer to perp and parallel to B
	k3 = dot_product(kvec, eq%bunit)
	k1 = sqrt( sum((kvec-k3*eq%bunit)**2) )
	R3 = k3/k0
	R1 = k1/k0
      R1S=R1**2
      R3S=R3**2
      RS=R1S+R3S

      B1=eq%gamma(0)
      BETAE=B1**2

! check if k||=0, if so there is no damping
	if (R3 == 0.) return

! Get warm plasma terms (note: B1 as defined here carries the sign of the
! electron charge (i.e. is negative)  Omega-sub-e in the notes does not )

!   Thermal speed:
    vth = sqrt( 2.*eq%ts(0)/ms(0) )
    VT = vth/clight

!      Z function.

    xi = (omgrf+eq%omgc(0))/(k3*vth)

! Check if arg too large to produce damping
    if (abs(xi)> 5.) return
!    zf = zfun0(cmplx(xi,kind=rkind), real(k3,kind=rkind))
    zf = zfun0_real_arg(xi, real(k3,kind=rkind))

	P=eq%alpha(0)
	Q=P/2./(1-B1)

	LAMBDA1=(1.-Q)*RS*R1S + (1.-P)*RS*R3S - (1.-Q)*(1.-P)*(RS+R3S)		&
  &          -(1-2.*Q)*R1S + (1-2*Q)*(1-P)

        LAMBDA2=-P/B1 * ( RS*R1S - (1.-2.*Q)*R1S ) +  				&
  &		P**2/4./BETAE * R1S/R3S * ( RS+R3S-2.*(1.-2.*Q) )

     	LAMBDA5= P*( RS*R3S-(1.-Q)*(RS+R3S) + (1.-2.*Q) )

	D_WARM = -(1.-B1)*R3*VT* 								&
   &		( LAMBDA1 + LAMBDA2 + R1S/2./R3/BETAE*VT*xi*LAMBDA5 ) * 		&
   &		(xi+1./zf)

! Get cold plasma terms

	A = 1.- P - BETAE

	B = -( (1.-P)*A + (1.-P)**2 - BETAE ) +						&
   &		( A + (1.-P)*(1.-BETAE) )*R3S

     	C = (1.-P) * ( (1.-P)**2-BETAE+(1.-BETAE)*R3S**2-2.*A*R3S )

	DDNX2 = 2.*A*R1S + B

	DDNZ = 2.*R3 * (  ( A +(1.-P)*(1.-BETAE) )*R1S + 				&
   &		(1-P) * ( 2.*(1.-BETAE)*R3S - 2.*A )  )

   	DNparallel = eq%bunit
	DNperp2 = 2*(nvec - R3*eq%bunit)

	DDN=DDNX2 * DNperp2 + DDNZ * DNparallel


! Choose direction for unit vector for k_imaginary (e.g. along x or
! along Vg)

	 vg_unit=vg/sqrt(sum(vg**2))

! calculate correction to refractive index

	DELTA = - D_WARM / dot_product(DDN, vg_unit)

    ksi(0) = k0 * aimag(DELTA)
	ksi(1:nspec) = 0.
	ki = ksi(0)

!		  write(*, *) 'x, eq%bmag, DELTA, ki = ', x, eq%bmag, DELTA, ki

      RETURN
      END
