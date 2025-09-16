 MODULE vectors3_m

! generic procedure: cross_product, calculates usual A x B
! generic procedure: triple_product, calculates usual A dot B x C
! generic procedure: triple_vector_product, calculates usual A x B x C

!_________________________________________________________________________________________
! Working notes:
!_________________________________________________________________________________________

!_________________________________________________________________________________________
! Module data
!_________________________________________________________________________________________

    IMPLICIT NONE

    integer, parameter :: rkind = selected_real_kind(15,307) ! kind parameter for reals
    integer, parameter :: skind = selected_real_kind(6,37) ! kind parameter single precision

    interface cross_product
        module procedure cross_product_D, cross_product_S
    end interface

    interface triple_product
        module procedure triple_product_D, triple_product_S
    end interface

    interface triple_vector_product
        module procedure triple_vector_product_D, triple_vector_product_S
    end interface

!_________________________________________________________________________________________
contains
!_________________________________________________________________________________________

    function cross_product_S(A, B)
        IMPLICIT NONE
        REAL(kind = skind), INTENT(IN) :: A(3), B(3)
        REAL(kind = skind) :: cross_product_S(3)
        cross_product_S = (/A(2)*B(3)- A(3)*B(2), A(3)*B(1) - A(1)*B(3), &
                            A(1)*B(2) - A(2)*B(1)/)
    return
    end function cross_product_S
 !*********************************************************************************

    function triple_product_S(A, B, C)
        IMPLICIT NONE
        REAL(kind = skind), INTENT(IN) :: A(3), B(3), C(3)
        REAL(kind = skind) :: triple_product_S
        triple_product_S = dot_product(A, cross_product_S(B,C))
    return
    end function triple_product_S
 !*********************************************************************************

    function triple_vector_product_S(A, B, C)
        IMPLICIT NONE
        REAL(kind = skind), INTENT(IN) :: A(3), B(3), C(3)
        REAL(kind = skind) :: triple_vector_product_S(3)
        triple_vector_product_S = cross_product(A, cross_product_S(B,C))
    return
    end function triple_vector_product_S
 !*********************************************************************************

    function cross_product_D(A, B)
        IMPLICIT NONE
        REAL(kind = rkind), INTENT(IN) :: A(3), B(3)
        REAL(kind = rkind) :: cross_product_D(3)
        cross_product_D = (/A(2)*B(3)- A(3)*B(2), A(3)*B(1) - A(1)*B(3), &
                            A(1)*B(2) - A(2)*B(1)/)
    return
    end function cross_product_D
 !*********************************************************************************

    function triple_product_D(A, B, C)
        IMPLICIT NONE
        REAL(kind = rkind), INTENT(IN) :: A(3), B(3), C(3)
        REAL(kind = rkind) :: triple_product_D
        triple_product_D = dot_product(A, cross_product_D(B,C))
    return
    end function triple_product_D
 !*********************************************************************************

    function triple_vector_product_D(A, B, C)
        IMPLICIT NONE
        REAL(kind = rkind), INTENT(IN) :: A(3), B(3), C(3)
        REAL(kind = rkind) :: triple_vector_product_D(3)
        triple_vector_product_D = cross_product(A, cross_product_D(B,C))
    return
    end function triple_vector_product_D

 END MODULE vectors3_m