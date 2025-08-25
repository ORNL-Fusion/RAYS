PROGRAM test_vectors3
   USE vectors3_m, only : cross_product, triple_product, triple_vector_product

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

	real(kind=skind) :: A_S(3), B_S(3), C_S(3), cross_S(3), trip_S, vec_S(3)
	real(kind=rkind) :: A_D(3), B_D(3), C_D(3), cross_D(3), trip_D, vec_D(3)
	real(kind=rkind) :: one = 1.0_rkind, zero = 0.0_rkind


! test single
	A_S = (/1.0, 0.0, 0.0/)
	B_S = (/0.0, 1., 0.0/)
	C_S = (/0.0, 0.0, 1.0/)

	vec_S = triple_vector_product(A_S,B_S,C_S)

    write(*,*) ' '
    write(*,*) 'single precision'
    WRITE (*, *) ' A x B  = ',cross_product(A_S,B_S)
    WRITE (*, *) ' B x C  = ',cross_product(B_S,C_S)
    WRITE (*, *) ' C x A  = ',cross_product(C_S,A_S)
    WRITE (*, *) ' A x A  = ',cross_product(A_S,A_S)
    WRITE (*, *) ' B x B  = ',cross_product(B_S,B_S)
    WRITE (*, *) ' C x C  = ',cross_product(C_S,C_S)
    WRITE (*, *) ' A dot B x C  = ',triple_product(A_S,B_S,C_S)
! 	vec_S = triple_vector_product(A_S,B_S,C_S)
!     WRITE (*, *) ' A x B x C  = ',vec_S
!     WRITE (*, *) ' A x B x A  = ',triple_vector_product(A_S,B_S,A_S)
    WRITE (*, *) ' A x B x C  = ',triple_vector_product(A_S,B_S,C_S)
    WRITE (*, *) ' A x B x A  = ',triple_vector_product(A_S,B_S,A_S)
    WRITE (*, *) ' A x C x A  = ',triple_vector_product(A_S,C_S,A_S)

! test double
	A_D = (/one, zero, zero/)
	B_D = (/zero, one, zero/)
	C_D = (/zero, zero, one/)

    write(*,*) ' '
    write(*,*) 'double precision'
    WRITE (*, *) ' A x B  = ',cross_product(A_D,B_D)
    WRITE (*, *) ' B x C  = ',cross_product(B_D,C_D)
    WRITE (*, *) ' C x A  = ',cross_product(C_D,A_D)
    WRITE (*, *) ' A x A  = ',cross_product(A_D,A_D)
    WRITE (*, *) ' B x B  = ',cross_product(B_D,B_D)
    WRITE (*, *) ' C x C  = ',cross_product(C_D,C_D)
    WRITE (*, *) ' A dot B x C  = ',triple_product(A_D,B_D,C_D)
    WRITE (*, *) ' A x B x C  = ',triple_vector_product(A_D,B_D,C_D)
    WRITE (*, *) ' A x B x A  = ',triple_vector_product(A_D,B_D,A_D)
    WRITE (*, *) ' A x C x A  = ',triple_vector_product(A_D,C_D,A_D)

END PROGRAM test_vectors3
