 MODULE matrix3x3_m

! Standard functions of 3x3 input matrix 'A'
! generic procedure: transpose3x3, for completeness, redundant with intrinsic
! generic procedure: adjoint3x3, conjugate transpose
! generic procedure: hermitian3x3
! generic procedure: anti_hermitian3x3
! generic procedure: determinant3x3
! N.B. transpose and determinant make sense for both real and complex matrices.
!      the others, adjoint, hermitian, anti_hermitian only apply to complex matrices

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    real(KIND=rkind),parameter :: zero = 0.0_rkind
    real(KIND=rkind),parameter :: one = 1.0_rkind
    real(KIND=rkind),parameter :: two = 2.0_rkind
    real(KIND=rkind),parameter :: ten = 10.0_rkind
    complex, parameter :: c0 = cmplx(0., 0.), c1 = cmplx(1., 0.), ci = cmplx(0., 1.0)
	complex(kind=rkind), parameter :: c0_rk = cmplx(zero, zero), c1_rk = cmplx(one, zero)
	complex(kind=rkind), parameter :: ci_rk = cmplx(zero,one)

    interface transpose3x3
        module procedure transpose_rkind, transpose_S, cmplx_transpose_rkind,&
              & cmplx_transpose_S
    end interface

    interface adjoint3x3
        module procedure adjoint_rkind, adjoint_S
    end interface

    interface hermitian3x3
        module procedure hermitian_rkind, hermitian_S
    end interface

    interface anti_hermitian3x3
        module procedure anti_hermitian_rkind, anti_hermitian_S
    end interface

    interface determinant3x3
        module procedure determinant_rkind, determinant_S, cmplx_determinant_rkind,&
             &  cmplx_determinant_S
    end interface

    interface inverse3x3
        module procedure inverse_rkind, inverse_S
    end interface

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

    pure function transpose_S(A)
        IMPLICIT NONE
        REAL, INTENT(IN) :: A(3,3)
        REAL :: transpose_S(3,3)
        transpose_S = transpose(A)
    return
    end function transpose_S

    pure function cmplx_transpose_S(A)
        IMPLICIT NONE
        COMPLEX, INTENT(IN) :: A(3,3)
        COMPLEX :: cmplx_transpose_S(3,3)
        cmplx_transpose_S = transpose(A)
    return
    end function cmplx_transpose_S

    pure function transpose_rkind(A)
        IMPLICIT NONE
        REAL(KIND=rkind), INTENT(IN) :: A(3,3)
        REAL(KIND=rkind) :: transpose_rkind(3,3)
        transpose_rkind = transpose(A)
    return
    end function transpose_rkind

    pure function cmplx_transpose_rkind(A)
        IMPLICIT NONE
        COMPLEX(KIND=rkind), INTENT(IN) :: A(3,3)
        COMPLEX(KIND=rkind) :: cmplx_transpose_rkind(3,3)
        cmplx_transpose_rkind = transpose(A)
    return
    end function cmplx_transpose_rkind
 !*********************************************************************************

    pure function adjoint_S(A)
        IMPLICIT NONE
        COMPLEX, INTENT(IN) :: A(3,3)
        COMPLEX :: adjoint_S(3,3)
        adjoint_S(1,1) = -(A(2,3)*A(3,2)) + A(2,2)*A(3,3)
        adjoint_S(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
        adjoint_S(1,3) = -(A(1,3)*A(2,2)) + A(1,2)*A(2,3)
        adjoint_S(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
        adjoint_S(2,2) = -(A(1,3)*A(3,1)) + A(1,1)*A(3,3)
        adjoint_S(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
        adjoint_S(3,1) = -(A(2,2)*A(3,1)) + A(2,1)*A(3,2)
        adjoint_S(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
        adjoint_S(3,3) = -(A(1,2)*A(2,1)) + A(1,1)*A(2,2)
    return
    end function adjoint_S

    pure function adjoint_rkind(A)
        IMPLICIT NONE
        COMPLEX(KIND=rkind), INTENT(IN) :: A(3,3)
        COMPLEX(KIND=rkind) :: adjoint_rkind(3,3)
        adjoint_rkind(1,1) = -(A(2,3)*A(3,2)) + A(2,2)*A(3,3)
        adjoint_rkind(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
        adjoint_rkind(1,3) = -(A(1,3)*A(2,2)) + A(1,2)*A(2,3)
        adjoint_rkind(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
        adjoint_rkind(2,2) = -(A(1,3)*A(3,1)) + A(1,1)*A(3,3)
        adjoint_rkind(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
        adjoint_rkind(3,1) = -(A(2,2)*A(3,1)) + A(2,1)*A(3,2)
        adjoint_rkind(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
        adjoint_rkind(3,3) = -(A(1,2)*A(2,1)) + A(1,1)*A(2,2)
    return
    end function adjoint_rkind
 !*********************************************************************************

    pure function hermitian_S(A)
        IMPLICIT NONE
        COMPLEX, INTENT(IN) :: A(3,3)
        COMPLEX :: hermitian_S(3,3)
        hermitian_S = ( A + conjg(transpose(A)))/2.
    return
    end function hermitian_S

    pure function hermitian_rkind(A)
        IMPLICIT NONE
        COMPLEX(KIND=rkind), INTENT(IN) :: A(3,3)
        COMPLEX(KIND=rkind) :: hermitian_rkind(3,3)
        hermitian_rkind = ( A + conjg(transpose(A)))/two
    return
    end function hermitian_rkind
 !*********************************************************************************

    pure function anti_hermitian_S(A)
        IMPLICIT NONE
        COMPLEX, INTENT(IN) :: A(3,3)
        COMPLEX :: anti_hermitian_S(3,3)
        anti_hermitian_S = ( A - conjg(transpose(A)))/2.
    return
    end function anti_hermitian_S

    pure function anti_hermitian_rkind(A)
        IMPLICIT NONE
        COMPLEX(KIND=rkind), INTENT(IN) :: A(3,3)
        COMPLEX(KIND=rkind) :: anti_hermitian_rkind(3,3)
        anti_hermitian_rkind = ( A - conjg(transpose(A)))/two
    return
    end function anti_hermitian_rkind
 !*********************************************************************************

    pure function determinant_S(A)
        IMPLICIT NONE
        REAL, INTENT(IN) :: A(3,3)
        REAL :: determinant_S
        determinant_S = &
          &   A(3,3)*(A(1,1)*A(2,2)-A(2,1)*A(1,2)) &
          & - A(3,2)*(A(1,1)*A(2,3)-A(2,1)*A(1,3)) &
          & + A(3,1)*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
    return
    end function determinant_S

    pure function cmplx_determinant_S(A)
        IMPLICIT NONE
        COMPLEX, INTENT(IN) :: A(3,3)
        COMPLEX :: cmplx_determinant_S
        cmplx_determinant_S = &
          &   A(3,3)*(A(1,1)*A(2,2)-A(2,1)*A(1,2)) &
          & - A(3,2)*(A(1,1)*A(2,3)-A(2,1)*A(1,3)) &
          & + A(3,1)*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
    return
    end function cmplx_determinant_S

    pure function determinant_rkind(A)
        IMPLICIT NONE
        REAL(kind = rkind), INTENT(IN) :: A(3,3)
        REAL(kind = rkind) :: determinant_rkind
        determinant_rkind = &
          &   A(3,3)*(A(1,1)*A(2,2)-A(2,1)*A(1,2)) &
          & - A(3,2)*(A(1,1)*A(2,3)-A(2,1)*A(1,3)) &
          & + A(3,1)*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
    return
    end function determinant_rkind

    pure function cmplx_determinant_rkind(A)
        IMPLICIT NONE
        COMPLEX(kind = rkind), INTENT(IN) :: A(3,3)
        COMPLEX(kind = rkind) :: cmplx_determinant_rkind
        cmplx_determinant_rkind = &
          &   A(3,3)*(A(1,1)*A(2,2)-A(2,1)*A(1,2)) &
          & - A(3,2)*(A(1,1)*A(2,3)-A(2,1)*A(1,3)) &
          & + A(3,1)*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
    return
    end function cmplx_determinant_rkind

!*********************************************************************************
! Inverse functions.  Arg matrix A is complex.  The inverse of a real matrix is in
! general complex, so just make arg complex too.

    pure function inverse_S(A)
        IMPLICIT NONE
        COMPLEX, INTENT(IN) :: A(3,3)
        COMPLEX :: inverse_S(3,3)
 	    complex :: det
        det = determinant3x3(A)
        inverse_S = adjoint_S(A)/det
    return
    end function inverse_S

    pure function inverse_rkind(A)
        IMPLICIT NONE
        COMPLEX(KIND=rkind), INTENT(IN) :: A(3,3)
        COMPLEX(KIND=rkind) :: inverse_rkind(3,3)
 	    complex(kind=rkind) :: det
        det = determinant3x3(A)
        inverse_rkind = adjoint_rkind(A)/det
    return
    end function inverse_rkind

 END MODULE matrix3x3_m