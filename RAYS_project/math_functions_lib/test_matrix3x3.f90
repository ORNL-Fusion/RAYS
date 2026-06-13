PROGRAM test_matrix3x3
   USE matrix3x3_m, only : transpose3x3, adjoint3x3, hermitian3x3, anti_hermitian3x3,&
                         & inverse3x3, determinant3x3
    use write_matrix_m, only : write_matrix
    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind),parameter :: zero = 0.0_rkind
    real(KIND=rkind),parameter :: one = 1.0_rkind
    real(KIND=rkind),parameter :: two = 2.0_rkind
    real(KIND=rkind),parameter :: ten = 10.0_rkind
    complex, parameter :: c0 = cmplx(0., 0.), c1 = cmplx(1., 0.), ci = cmplx(0., 1.0)
	complex(kind=rkind), parameter :: c0_rk = cmplx(zero, zero), c1_rk = cmplx(one, zero)
	complex(kind=rkind), parameter :: ci_rk = cmplx(zero,one)

	real :: A_S_re(3,3), U_S(3,3), L_S(3,3), D_S(3,3)
	complex :: A_S_cmplx(3,3)
	real(kind=rkind) :: A_rk_re(3,3), U_rk(3,3), L_rk(3,3), D_rk(3,3)
	complex(kind=rkind) :: A_rk_cmplx(3,3)

	complex :: cId_S(3,3)
	complex(kind=rkind) :: cId_rk(3,3)

	cId_S = c0
	cId_S(1,1) = c1; cId_S(2,2) = c1; cId_S(3,3) = c1
	cId_rk = c0_rk
	cId_rk(1,1) = c1_rk; cId_rk(2,2) = c1_rk; cId_rk(3,3) = c1_rk;

	A_S_re = 0.; A_S_cmplx = c0; U_S = 0.; L_S = 0.; D_S = 0.;
	A_rk_re = zero; A_rk_cmplx = c0; U_rk = zero ; L_rk = zero; D_rk = zero


	U_S(1,2) = 1.0; U_S(1,3) = 1.0; U_S(2,3) = 1.0
	U_rk(1,2) = c1_rk; U_rk(1,3) = c1_rk; U_rk(2,3) = c1_rk
	L_S(2,1) = 1.0; L_S(3,1) = 1.0; L_S(3,2) = 1.0
	L_rk(2,1) = c1_rk; L_rk(3,1) = c1_rk; L_rk(3,2) = c1_rk
	D_S(1,1) = 1.0; 	D_S(2,2) = 2.0; 	D_S(3,3) = 3.0;
	D_rk(1,1) = one; 	D_rk(2,2) = two; 	D_rk(3,3) = one+two;

	call write_matrix('U_S', U_S, 3,3)
	call write_matrix('L_S', L_S, 3,3)
	call write_matrix('D_S', D_S, 3,3)
	write(*,*) ' '

	A_S_re = U_S + D_S
	call write_matrix('A', A_S_re, 3,3)
	call write_matrix('transpose(A)', transpose3x3(A_S_re), 3,3)
	write(*,*) 'determinant(A)', determinant3x3(A_S_re)
	write(*,*) ' '

	A_rk_re = U_rk + D_rk
	call write_matrix('A', A_rk_re, 3,3)
	call write_matrix('transpose(A)', transpose3x3(A_rk_re), 3,3)
	write(*,*) 'determinant(A)', determinant3x3(A_rk_re)
	write(*,*) ' '

	A_S_cmplx = ci*U_S + D_S
	call write_matrix('A', A_S_cmplx, 3,3)
	call write_matrix('transpose(A)', transpose3x3(A_S_cmplx), 3,3)
	call write_matrix('hermitian(A)', hermitian3x3(A_S_cmplx), 3,3)
	call write_matrix('anti_hermitian(A)', anti_hermitian3x3(A_S_cmplx), 3,3)
	call write_matrix('adjoint(A)', adjoint3x3(A_S_cmplx), 3,3)
	write(*,*) 'determinant(A)', determinant3x3(A_S_cmplx)
	call write_matrix('inverse(A)', inverse3x3(A_S_cmplx), 3,3)
	call write_matrix('A x Ainv', matmul(A_S_cmplx, inverse3x3(A_S_cmplx)), 3,3)
	call write_matrix('A x Ainv -I', matmul(A_S_cmplx, inverse3x3(A_S_cmplx))-cId_s, 3,3)
	write(*,*) ' '

	A_rk_cmplx = ci*U_rk + D_rk
	call write_matrix('A', A_rk_cmplx, 3,3)
	call write_matrix('transpose(A)', transpose3x3(A_rk_cmplx), 3,3)
	call write_matrix('hermitian(A)', hermitian3x3(A_rk_cmplx), 3,3)
	call write_matrix('anti_hermitian(A)', anti_hermitian3x3(A_rk_cmplx), 3,3)
	call write_matrix('adjoint(A)', adjoint3x3(A_rk_cmplx), 3,3)
	write(*,*) 'determinant(A)', determinant3x3(A_rk_cmplx)
	call write_matrix('inverse(A)', inverse3x3(A_rk_cmplx), 3,3)
	call write_matrix('A x Ainv', matmul(A_rk_cmplx, inverse3x3(A_rk_cmplx)), 3,3)
	call write_matrix('A x Ainv - I', matmul(A_rk_cmplx, inverse3x3(A_rk_cmplx)) - cID_rk&
	     &, 3,3)
	write(*,*) ' '

	A_rk_cmplx = ci*U_rk - ci*L_rk + D_rk
	call write_matrix('A', A_rk_cmplx, 3,3)
	call write_matrix('transpose(A)', transpose3x3(A_rk_cmplx), 3,3)
	call write_matrix('hermitian(A)', hermitian3x3(A_rk_cmplx), 3,3)
	call write_matrix('anti_hermitian(A)', anti_hermitian3x3(A_rk_cmplx), 3,3)
	call write_matrix('adjoint(A)', adjoint3x3(A_rk_cmplx), 3,3)
	write(*,*) 'determinant(A)', determinant3x3(A_rk_cmplx)
	call write_matrix('inverse(A)', inverse3x3(A_rk_cmplx), 3,3)
	call write_matrix('A x Ainv', matmul(A_rk_cmplx, inverse3x3(A_rk_cmplx)), 3,3)
	call write_matrix('A x Ainv - I', matmul(A_rk_cmplx, inverse3x3(A_rk_cmplx)) - cID_rk&
	     &, 3,3)
	write(*,*) ' '

	A_rk_cmplx = ci*U_rk + ci*L_rk
	call write_matrix('A', A_rk_cmplx, 3,3)
	call write_matrix('transpose(A)', transpose3x3(A_rk_cmplx), 3,3)
	call write_matrix('hermitian(A)', hermitian3x3(A_rk_cmplx), 3,3)
	call write_matrix('anti_hermitian(A)', anti_hermitian3x3(A_rk_cmplx), 3,3)
	call write_matrix('adjoint(A)', adjoint3x3(A_rk_cmplx), 3,3)
	write(*,*) 'determinant(A)', determinant3x3(A_rk_cmplx)
	call write_matrix('inverse(A)', inverse3x3(A_rk_cmplx), 3,3)
	call write_matrix('A x Ainv', matmul(A_rk_cmplx, inverse3x3(A_rk_cmplx)), 3,3)
	call write_matrix('A x Ainv - I', matmul(A_rk_cmplx, inverse3x3(A_rk_cmplx)) - cID_rk&
	     &, 3,3)
	write(*,*) ' '



END PROGRAM test_matrix3x3
