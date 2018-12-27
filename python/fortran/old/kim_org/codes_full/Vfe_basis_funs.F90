MODULE VFE_BASIS_FUNS   

    USE KiapsBase,  ONLY : iulog => KIM_IU_LOG, r8 => KIM_REAL8_KIND, i4 => KIM_INT_KIND
    USE Dimensions, ONLY : plev => nlev, plevp => nlevp
    USE VFE_MOD,    ONLY : debug_vfe, basis_degree, deri_order, nquad, basis_order, nlevbc, nknot, nbasis
    
    IMPLICIT NONE

    PUBLIC :: COMPUTE_BASIS_FUNS_ETA,  &
              COMPUTE_BASIS_FUNS_QUAD, &
              BASIS_FUNS, &
              BASIS_DERIS, &
              FIND_KNOT_SPAN, &
              DEFINE_QUADRATURE_SET

CONTAINS
!
!---------------------------------------------------------------------------------!
! SUBROUTINE COMPUTE_BASIS_FUNS_ETA: To compute basis functions at eta full levels
! 
!---------------------------------------------------------------------------------!
!
SUBROUTINE COMPUTE_BASIS_FUNS_ETA(etaf,knots,basis_fun)    
!   IMPLICIT NONE    
    REAL(r8), INTENT(IN),  DIMENSION(0:plev+1)      :: etaf
    REAL(r8), INTENT(IN),  DIMENSION(nknot)         :: knots
    REAL(r8), INTENT(OUT), DIMENSION(nlevbc,nbasis) :: basis_fun
    
    INTEGER(i4)                                     :: ip
    INTEGER(i4)                                     :: ileft
    REAL(r8)                                        :: eta_point
    REAL(r8), DIMENSION(nbasis)                     :: basisf
    !
    !----------------------------------------------------------------------------------------------!
    ! B-spline basis functions at eta full level points
    !----------------------------------------------------------------------------------------------!    
    !
    IF (debug_vfe) THEN
        WRITE (*, '(a)') ''
        WRITE (*, '(a)') '  B-spline basis functions at eta full level points'
        WRITE (*, '(a)') '  i  left   eta           B1(X)         B2(X)         B3(X), .....'
        WRITE (*, '(a)') '  ----------------------------------------------------------------'
    END IF
    
    basisf(:) = 0.0_r8
    basis_fun(:,:) = 0.0_r8
    DO ip = 1, nlevbc
        eta_point = etaf(ip-1) 
        
        ! Find left-most points with x-point from the knot
        CALL FIND_KNOT_SPAN(nbasis,basis_order,knots,eta_point,ileft)
        
        ! Get B(I,K)(X) in VALUES(1:N): B(LEFT-K+1,K)(X), ..., B(LEFT,K)(X). All the others are zero.
        CALL BASIS_FUNS(basis_order,ileft,knots,eta_point,basisf(ileft+1-basis_order))
        basis_fun(ip,:) = basisf(:)
        
        IF (debug_vfe) WRITE (*, '(2I4, F15.10, 100F15.10)') ip, ileft, eta_point, basisf(1:nbasis)  
        
        ! Zero out the values just computed in preparation for the next evaluation point
        basisf(:) = 0.0_r8
    END DO   
END SUBROUTINE COMPUTE_BASIS_FUNS_ETA
!
!---------------------------------------------------------------------------------!
! SUBROUTINE COMPUTE_BASIS_FUNS_QUAD: To compute basis functions at eta quadrature points 
! 
!---------------------------------------------------------------------------------!
!
SUBROUTINE COMPUTE_BASIS_FUNS_QUAD(etaf,knots,basis_fun_quad,basis_der_quad)    
!   IMPLICIT NONE    
    
    REAL(r8), INTENT(IN),  DIMENSION(0:plev+1)              :: etaf
    REAL(r8), INTENT(IN),  DIMENSION(nknot)                 :: knots
    REAL(r8), INTENT(OUT), DIMENSION(nquad*(plev+1),nbasis) :: basis_fun_quad
    REAL(r8), INTENT(OUT), DIMENSION(nquad*(plev+1),nbasis) :: basis_der_quad
    
    INTEGER(i4)                                             :: ip, icount, iquad, ileft
    REAL(r8)                                                :: eta_point
    REAL(r8), DIMENSION(nquad)                              :: quadx
    REAL(r8), DIMENSION(nquad)                              :: quadw
    REAL(r8), DIMENSION(nbasis)                             :: basisf
    REAL(r8), DIMENSION(nbasis)                             :: basisd
    !
    !----------------------------------------------------------------------------------------------!
    ! B-spline basis functions at eta full level points
    !----------------------------------------------------------------------------------------------!    
    !
    CALL DEFINE_QUADRATURE_SET(nquad,quadw,quadx)

    IF (debug_vfe) THEN 
        WRITE (*, '(a)') ''
        WRITE (*, '(a)') '  B-spline basis functions at gauss quadrature points'
        WRITE (*, '(a)') '  i  left   eta           B1(X)         B2(X)         B3(X), .....'
        WRITE (*, '(a)') '  ----------------------------------------------------------------'
    END IF
    
    icount = 1  
    basisf(:) = 0.0_r8
    basisd(:) = 0.0_r8
    DO ip = 1, nlevbc-1        
        DO iquad = 1, nquad
            eta_point = ((1.0_r8-quadx(iquad))*etaf(ip-1)+(quadx(iquad)-(-1.0_r8))*etaf(ip))/(1.0_r8-(-1.0_r8))

            ! Find left-most points with x-point from the knot
            CALL FIND_KNOT_SPAN(nbasis,basis_order,knots,eta_point,ileft)
            
            CALL BASIS_DERIS(basis_order,deri_order,ileft,knots,eta_point,basisf(ileft+1-basis_order),basisd(ileft+1-basis_order))
            IF (debug_vfe) WRITE (*, '(2I4, F10.5, 100F12.7)') icount, ileft, eta_point, basisf(1:nbasis)  
            basis_fun_quad(icount,:) = basisf(:)
            basis_der_quad(icount,:) = basisd(:)
            
            icount = icount + 1
            basisf(:) = 0.0_r8
            basisd(:) = 0.0_r8
        END DO
    END DO
END SUBROUTINE COMPUTE_BASIS_FUNS_QUAD
!
!---------------------------------------------------------------------------------!
! SUBROUTINE BASIS_FUNS
! 
!---------------------------------------------------------------------------------!
!
SUBROUTINE BASIS_FUNS(basis_order,ileft,knots,pts,basisf)
    IMPLICIT NONE
    
    INTEGER(i4), INTENT(IN)                              :: basis_order
    INTEGER(i4), INTENT(IN)                              :: ileft
    
    REAL(r8), INTENT(IN)                                 :: pts
    REAL(r8), INTENT(IN),DIMENSION(ileft+basis_order)    :: knots        
    REAL(r8), INTENT(OUT),DIMENSION(basis_order)         :: basisf
    
    INTEGER(i4)                                          :: i, j, k
    INTEGER(i4)                                          :: basis_degree
    
    REAL(r8)                                             :: saved, temp
    REAL(r8), DIMENSION(basis_order-1)                   :: delta_left
    REAL(r8), DIMENSION(basis_order-1)                   :: delta_rght
    REAL(r8), DIMENSION(0:basis_order-1,0:basis_order-1) :: NDU

    basis_degree = basis_order - 1
    !
    !---------------------------------------------------------------------------------!
    ! To compute B-spline bassis function
    !---------------------------------------------------------------------------------!
    !    
    NDU(0,0) = 1.0_r8
    
    DO j = 1, basis_degree
        delta_left(j) = pts - knots(ileft+1-j)
        delta_rght(j) = knots(ileft+j) - pts
        
        saved = 0.0_r8        
        DO k = 0, j-1
            NDU(j,k) = delta_rght(k+1) + delta_left(j-k)
            temp = NDU(k,j-1)/NDU(j,k)
            NDU(k,j) = saved + delta_rght(k+1)*temp
            saved = delta_left(j-k)*temp
        END DO
        NDU(j,j) = saved
    END DO
    
    ! B-spline basis function 
    DO i = 1, basis_order
        basisf(i) = NDU(i-1,basis_degree) 
    END DO            
END SUBROUTINE basis_funs
!
!---------------------------------------------------------------------------------!
! SUBROUTINE BASIS_DERIS
! 
!---------------------------------------------------------------------------------!
!
SUBROUTINE BASIS_DERIS(basis_order,deri_order,ileft,knots,pts,basisf,basisd)
    IMPLICIT NONE
    
    INTEGER(i4), INTENT(IN)                              :: basis_order
    INTEGER(i4), INTENT(IN)                              :: deri_order
    INTEGER(i4), INTENT(IN)                              :: ileft

    REAL(r8), INTENT(IN)                                 :: pts
    REAL(r8), INTENT(IN),  DIMENSION(ileft+basis_order)  :: knots
    REAL(r8), INTENT(OUT), DIMENSION(basis_order)        :: basisf
    REAL(r8), INTENT(OUT), DIMENSION(basis_order)        :: basisd
    
    INTEGER(i4)                                          :: i, j, k, L
    INTEGER(i4)                                          :: ROW1, ROW2, DIFF_LK, DIFF_PK, LOW, HIGH
    INTEGER(i4)                                          :: basis_degree
    
    REAL(r8)                                             :: saved, temp    
    REAL(r8), DIMENSION(basis_order-1)                   :: delta_left
    REAL(r8), DIMENSION(basis_order-1)                   :: delta_rght
    REAL(r8), DIMENSION(0:deri_order,0:basis_order)      :: deris
    REAL(r8), DIMENSION(0:1,0:basis_order)               :: A
    REAL(r8), DIMENSION(0:basis_order-1,0:basis_order-1) :: NDU
    !
    basis_degree = basis_order - 1
    !
    !---------------------------------------------------------------------------------!
    ! To compute B-spline bassis function
    !---------------------------------------------------------------------------------!
    !    
    NDU(0,0) = 1.0_r8
    
    DO j = 1, basis_order-1
        delta_left(j) = pts - knots(ileft+1-j)
        delta_rght(j) = knots(ileft+j) - pts
        
        saved = 0.0_r8        
        DO k = 0, j-1
            NDU(j,k) = delta_rght(k+1) + delta_left(j-k)
            temp = NDU(k,j-1)/NDU(j,k)
            NDU(k,j) = saved + delta_rght(k+1)*temp
            saved = delta_left(j-k)*temp
        END DO
        NDU(j,j) = saved
    END DO
    !
    !---------------------------------------------------------------------------------!
    ! To compute the k-th order derivative of B-spline bassis function
    !---------------------------------------------------------------------------------!
    !
    DO j = 0, basis_degree
        deris(0,j) = NDU(j,basis_degree)
    END DO
    
    DO L = 0, basis_degree
        ROW1 = 0
        ROW2 = 1
        A(0,0) = 1.0_r8
        
        DO K = 1, deri_order
            temp = 0.0_r8
            DIFF_LK = L - K
            DIFF_PK = basis_degree - K
            
            IF (DIFF_LK>=0) THEN
                A(ROW2,0) = A(ROW1,0) / NDU(DIFF_PK+1,DIFF_LK)
                temp = A(ROW2,0) * NDU(DIFF_LK,DIFF_PK)
            END IF
            
            IF (DIFF_LK>=-1) THEN
                LOW = 1
            ELSE
                LOW = -DIFF_LK
            END IF
            
            IF (L-1<=DIFF_PK) THEN
                HIGH = K - 1
            ELSE
                HIGH = basis_degree - L
            END IF
            
            DO J = LOW, HIGH
                A(ROW2,J) = (A(ROW1,J)-A(ROW1,J-1)) / NDU(DIFF_PK+1,DIFF_LK+J)
                temp = temp + A(ROW2,J) * NDU(DIFF_LK+J,DIFF_PK)
            END DO
            
            IF (L<=DIFF_PK) THEN
                A(ROW2,K) = -A(ROW1,K-1) / NDU(DIFF_PK+1,L)
                temp = temp + A(ROW2,K) * NDU(L,DIFF_PK)
            END IF
            
            deris(K,L) = temp
            J = ROW1
            ROW1 = ROW2
            ROW2 = J
        END DO
    END DO
    
    L = basis_degree
    DO K = 1, deri_order
        DO J = 0, basis_degree
            deris(K,J) = deris(K,J) * L
        END DO    
        L = L*(basis_degree-K)
    END DO
    
    DO i = 1, basis_order
        basisf(i) = deris(0,i-1) ! Fisrt derivative of B-spline basis 
        basisd(i) = deris(1,i-1) ! Fisrt derivative of B-spline basis 
    END DO            
END SUBROUTINE BASIS_DERIS
!
!---------------------------------------------------------------------------------!
! SUBROUTINE FIND_KNOT_SPAN
! 
!---------------------------------------------------------------------------------!
!
SUBROUTINE FIND_KNOT_SPAN(nbasis,basis_order,knots,pts,ileft)
    IMPLICIT NONE
    
    REAL(r8), PARAMETER                                    :: TOLERANCE=5.D-15 
    
    INTEGER(i4), INTENT(IN)                                :: nbasis
    INTEGER(i4), INTENT(IN)                                :: basis_order
    REAL(r8),    INTENT(IN)                                :: pts
    REAL(r8),    INTENT(IN), DIMENSION(nbasis+basis_order) :: knots
 
    INTEGER(i4)                                            :: ileft
    INTEGER(i4)                                            :: basis_degree, nknot, LOW, HIGH, MID

    nknot = nbasis + basis_order
    basis_degree = basis_order - 1
      
    IF ((pts<knots(1)-TOLERANCE) .OR. (pts>knots(nknot)+TOLERANCE)) THEN
        WRITE (*, '(3X, A, F10.5)') 'STOP! The point is out of range:', pts
        STOP
        ileft = -1
    ELSE IF (DABS(pts-knots(nknot-basis_degree))<TOLERANCE) THEN
        ileft = nknot - basis_order
    ELSE
        LOW = basis_degree
        HIGH = nknot - basis_degree
        MID = (LOW+HIGH)/2
        DO WHILE (pts<knots(MID) .OR. pts>=knots(MID+1))
            IF (pts<knots(MID)) THEN
                HIGH = MID
            ELSE
                LOW = MID
            END IF
            MID = (LOW+HIGH)/2
        END DO
        ileft = MID
    END IF
END SUBROUTINE FIND_KNOT_SPAN
!
!---------------------------------------------------------------------------------!
! SUBROUTINE DEFINE_QUADRATURE_SET
! 
!---------------------------------------------------------------------------------!
!
SUBROUTINE DEFINE_QUADRATURE_SET(nquad,quadw,quadx)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)                     :: nquad
    REAL(r8), DIMENSION(nquad), INTENT(OUT) :: quadw, quadx
    
    IF (nquad == 1) THEN
        quadx(1) =  0.0_r8
        quadw(1) =  2.0_r8
    ELSE IF (nquad == 2) THEN
        quadx(1) = -0.5773502691896257645091487805020_r8
        quadx(2) =  0.5773502691896257645091487805020_r8

        quadw(1) =  1.0_r8
        quadw(2) =  1.0_r8
    ELSE IF (nquad == 3) THEN
        quadx(1) = -0.7745966692414833770358530799560_r8
        quadx(2) =  0.0_r8
        quadx(3) =  0.7745966692414833770358530799560_r8

        quadw(1) =  5.0_r8/9.0_r8
        quadw(2) =  8.0_r8/9.0_r8
        quadw(3) =  5.0_r8/9.0_r8
    ELSE IF (nquad == 4) THEN
        quadx(1) = -0.8611363115940525752239464888930_r8
        quadx(2) = -0.3399810435848562648026657591030_r8
        quadx(3) =  0.3399810435848562648026657591030_r8
        quadx(4) =  0.8611363115940525752239464888930_r8

        quadw(1) =  0.3478548451374538573730639492220_r8
        quadw(2) =  0.6521451548625461426269360507780_r8
        quadw(3) =  0.6521451548625461426269360507780_r8
        quadw(4) =  0.3478548451374538573730639492220_r8
    ELSE IF (nquad == 5) THEN
        quadx(1) = -0.9061798459386639927976268782990_r8
        quadx(2) = -0.5384693101056830910363144207000_r8
        quadx(3) =  0.0_r8
        quadx(4) =  0.5384693101056830910363144207000_r8
        quadx(5) =  0.9061798459386639927976268782990_r8

        quadw(1) =  0.2369268850561890875142640407200_r8
        quadw(2) =  0.4786286704993664680412915148360_r8
        quadw(3) =  0.5688888888888888888888888888890_r8
        quadw(4) =  0.4786286704993664680412915148360_r8
        quadw(5) =  0.2369268850561890875142640407200_r8
    ELSE IF (nquad == 6) THEN
        quadx(1) = -0.9324695142031520278123015544940_r8
        quadx(2) = -0.6612093864662645136613995950200_r8
        quadx(3) = -0.2386191860831969086305017216810_r8
        quadx(4) =  0.2386191860831969086305017216810_r8
        quadx(5) =  0.6612093864662645136613995950200_r8
        quadx(6) =  0.9324695142031520278123015544940_r8

        quadw(1) =  0.1713244923791703450402961421730_r8
        quadw(2) =  0.3607615730481386075698335138380_r8
        quadw(3) =  0.4679139345726910473898703439900_r8
        quadw(4) =  0.4679139345726910473898703439900_r8
        quadw(5) =  0.3607615730481386075698335138380_r8
        quadw(6) =  0.1713244923791703450402961421730_r8
    END IF 
END SUBROUTINE DEFINE_QUADRATURE_SET

END MODULE VFE_BASIS_FUNS
