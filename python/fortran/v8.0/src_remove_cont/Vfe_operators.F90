MODULE VFE_OPERATORS

    USE KiapsBase,      ONLY : iulog => KIM_IU_LOG, r8 => KIM_REAL8_KIND, i4 => KIM_INT_KIND
    USE Dimensions,     ONLY : plev => nlev, plevp => nlevp
    USE VFE_MOD,        ONLY : debug_vfe, basis_degree, deri_order, nquad, basis_order, nlevbc, nknot, nbasis,MINV
    USE VFE_BASIS_FUNS, ONLY : COMPUTE_BASIS_FUNS_ETA, COMPUTE_BASIS_FUNS_QUAD, DEFINE_QUADRATURE_SET 
    USE KiapsParallel,  ONLY : KIM_Par
        
    IMPLICIT NONE

    PUBLIC :: SUB_DERI_OPERATOR, SUB_INTE_OPERATOR

CONTAINS
!
!---------------------------------------------------------------------------------!
! SUBROUTINE SUB_DERI_OPERATOR: Vertical derivative operator (VOPDERI)
!---------------------------------------------------------------------------------!
!
SUBROUTINE SUB_DERI_OPERATOR( debug_vfe,etaf,knots,VOPDERI )
    IMPLICIT NONE

    ! Variables for hybrid coordinates    
    LOGICAL,  INTENT(IN)                        :: debug_vfe
    REAL(r8), INTENT(IN),  DIMENSION(0:plev+1)  :: etaf
    REAL(r8), INTENT(IN),  DIMENSION(nknot)     :: knots
    REAL(r8), INTENT(OUT), DIMENSION(plev,plev) :: VOPDERI
        
    ! Variables for b-spline basis function 
    REAL(r8), DIMENSION(nlevbc,nbasis)          :: bspline_funs
    REAL(r8), DIMENSION(nquad*(plev+1),nbasis)  :: ibspline_funs_quad
    REAL(r8), DIMENSION(nquad*(plev+1),nbasis)  :: kbspline_funs_quad
    REAL(r8), DIMENSION(nquad*(plev+1),nbasis)  :: jbspline_deri_quad

    ! Variables for VFE
    REAL(r8)                                    :: ZDET
    REAL(r8), DIMENSION(nlevbc,nbasis)          :: ZP
    REAL(r8), DIMENSION(nlevbc,nbasis)          :: ZPINV
    REAL(r8), DIMENSION(nlevbc,nbasis)          :: ZS

    REAL(r8), DIMENSION(nlevbc,nbasis)          :: ZA
    REAL(r8), DIMENSION(nlevbc,nbasis)          :: ZAINV
    REAL(r8), DIMENSION(nlevbc,nbasis)          :: ZB    
    REAL(r8), DIMENSION(nlevbc,nbasis)          :: ZDERI
    REAL(r8), DIMENSION(nlevbc,nbasis)          :: ZTEMP

    INTEGER(i4)                                 :: i, j, k, ilev, icount, iquad
    REAL(r8)                                    :: deta, ZW1, ZW2
    REAL(r8), DIMENSION(nquad)                  :: quadx, quadw

    !
    !--------------------------------------------------------------------------------------------------------!
    ! STEP 1: Matrix P^-1
    ! ZP(i,j) = BASISFUNC_j(eta_i)) from FE space of the function to be differenciated into physical space
    ! i: row->eta points, j: column->basis function
    !--------------------------------------------------------------------------------------------------------!
    !    
    CALL COMPUTE_BASIS_FUNS_ETA(etaf,knots,bspline_funs)
    ZP(:,:) = bspline_funs(:,:)

    ! Invert SF: ZSFINV projects the function to be integrated from physical space onto FE space
    ZDET = 0.0_r8
    CALL MINV(ZP,nlevbc,nbasis,ZDET,ZPINV)
    !
    !--------------------------------------------------------------------------------------------------------!
    ! STEP 2: Matrix S
    ! ZS(i,j)=BASISDER_j(eta_i) to project the derivative from FE space to physical space 
    !--------------------------------------------------------------------------------------------------------!    
    !
    ZS(:,:) = ZP(:,:)
  
    IF (debug_vfe) THEN 
        WRITE (*, '(A)') ''
        WRITE (*, '(A)') '  Matrix S'
        WRITE (*, '(A)') '  ----------------------------------------------------------------'
        DO i = 1, nlevbc
            WRITE (*, '(50F6.1)') (100.*ZS(i,j), j = 1, nbasis)
        END DO
    END IF
    !
    ! Compute basis functions E_i(eta), E_k(eta), and dE_j(eta)/deta at eta quadrature points (quadx) 
    !
    CALL DEFINE_QUADRATURE_SET(nquad,quadw,quadx)
    CALL COMPUTE_BASIS_FUNS_QUAD(etaf,knots,ibspline_funs_quad,jbspline_deri_quad)
    kbspline_funs_quad(:,:) = ibspline_funs_quad(:,:)   
    !
    !--------------------------------------------------------------------------------------------------------!    
    ! STEP 3: Mass matrix A
    ! Compute A_ik(eta): integrate the product of basis functions at all quadrature points
    !--------------------------------------------------------------------------------------------------------!
    !
    DO i = 1, nbasis
        DO k = 1, nbasis        
            icount = 1
            ZA(i,k) = 0.0_r8
            DO ilev = 1, nlevbc-1
                deta = etaf(ilev) - etaf(ilev-1)
                DO iquad = 1, nquad
                    ZA(i,k) = ZA(i,k) + 0.50_r8*deta*quadw(iquad)*ibspline_funs_quad(icount,i)*kbspline_funs_quad(icount,k)
                    icount = icount + 1
                END DO                
            END DO  
         END DO      
    END DO  
    !
    ! Invert mass matrix: A^-1 
    ZDET = 0.0_r8
    CALL MINV(ZA,nbasis,nbasis,ZDET,ZAINV)
    
    IF (debug_vfe) THEN 
        WRITE (*, '(a)') ''
        WRITE (*, '(a)') '  Mass matrix: A'
        WRITE (*, '(a)') '  ----------------------------------------------------------------'
        DO i = 1, nbasis
            WRITE (*, '(50F6.1)') (100.*ZA(i,j), j = 1, nbasis)
        END DO
    END IF
    !
    !--------------------------------------------------------------------------------------------------------!    
    ! STEP 4: Operator matrix B
    ! Compute B_kj(eta): integrate the product of basis function and its derivative at all quadrature points
    !--------------------------------------------------------------------------------------------------------!    
    !
    DO j = 1, nbasis
        DO k = 1, nbasis       
            icount = 1
            ZB(k,j) = 0.0_r8
            DO ilev = 1, nlevbc-1
                deta = etaf(ilev) - etaf(ilev-1)
                DO iquad = 1, nquad
                    ZB(k,j) = ZB(k,j) + 0.50_r8*deta*quadw(iquad)*jbspline_deri_quad(icount,j)*kbspline_funs_quad(icount,k)
                    icount = icount + 1
                END DO                
            END DO  
         END DO      
    END DO
    
    IF (debug_vfe) THEN 
        WRITE (*, '(a)') ''
        WRITE (*, '(a)') '  Operator matrix: B'
        WRITE (*, '(a)') '  ----------------------------------------------------------------'
        DO k = 1, nbasis
            WRITE (*, '(50F6.1)') (100.*ZB(k,j), j = 1, nbasis)
        END DO
    END IF
    !
    !--------------------------------------------------------------------------------------------------------!    
    ! STEP 5: ZDERI
    ! Compute the product of S*A^-1*B*P^-1
    !--------------------------------------------------------------------------------------------------------!    
    !
    ! B*P^-1
    ZDERI(:,:) = 0.0_r8
    DO j = 1, nbasis
        DO i = 1, nbasis
            DO k = 1, nbasis
                ZDERI(i,j) = ZDERI(i,j) + ZB(i,k)*ZPINV(k,j)        
            END DO
        END DO
    END DO
    !
    ! A^-1*B*P^-1
    ZTEMP(:,:) = 0.0_r8
    DO j = 1, nbasis
        DO i = 1, nbasis
            DO k = 1, nbasis
                ZTEMP(i,j) = ZTEMP(i,j) + ZAINV(i,k)*ZDERI(k,j)
            END DO
        END DO
    END DO
    !
    ! ZDERI=S*A^-1*B*P^-1
    ZDERI(:,:) = 0.0_r8
    DO j = 1, nbasis
        DO i = 1, nbasis
            DO k = 1, nbasis
                ZDERI(i,j) = ZDERI(i,j) + ZS(i,k)*ZTEMP(k,j)
            END DO
        END DO
    END DO
    
    IF (debug_vfe) THEN 
        WRITE (*, '(a)') ''
        WRITE (*, '(a)') '  Product matrix: ZDERI'
        WRITE (*, '(a)') '  ----------------------------------------------------------------'
        DO k = 1, nbasis
            WRITE (*, '(50F6.1)') (100.*ZDERI(k,j), j = 1, nbasis)
        END DO
    END IF
    !
    !--------------------------------------------------------------------------------------------------------!    
    ! STEP 6: VOPDERI
    ! Reconstruction: ZDERI(N+2,N+2) -> VOPDERI(N,N)
    ! 
    ! Boundary conditions: linear extrapolation to nodes 0 and plev+1 derivative at nodes 1 and plev
    !--------------------------------------------------------------------------------------------------------!
    !
    ZW1 = 1.0_r8
    ZW2 = 0.0_r8

    DO i = 1, plev
        DO j = 2, plev-1
            VOPDERI(i,j) = ZDERI(i+1,j+1)
        END DO
        VOPDERI(i,1) = ZDERI(i+1,2) + ZW1*ZDERI(i+1,1)
        VOPDERI(i,2) = ZDERI(i+1,3) + ZW2*ZDERI(i+1,1)
        VOPDERI(i,plev-1) = ZDERI(i+1,plev  ) + ZW2*ZDERI(i+1,plev+2)
        VOPDERI(i,plev  ) = ZDERI(i+1,plev+1) + ZW1*ZDERI(i+1,plev+2)
    END DO
    
    IF (debug_vfe) THEN 
        WRITE (*, '(a)') ''
        WRITE (*, '(a)') '  Derivative operator matrix: VOPDERI(N,N)'
        WRITE (*, '(a)') '  ----------------------------------------------------------------'
        DO i = 1, plev
            WRITE (*, '(50F6.2)') (VOPDERI(i,j), j = 1, plev)
        END DO
    END IF

END SUBROUTINE SUB_DERI_OPERATOR
!
!---------------------------------------------------------------------------------!
! SUBROUTINE SUB_INTE_OPERATOR: Vertical integral operator (VOPINTE)
!---------------------------------------------------------------------------------!
!
SUBROUTINE SUB_INTE_OPERATOR( debug_vfe,etaf,knots,VOPINTE )
    IMPLICIT NONE

    ! Variables for hybrid coordinates    
    LOGICAL,  INTENT(IN)                          :: debug_vfe
    REAL(r8), INTENT(IN),  DIMENSION(0:plev+1)    :: etaf
    REAL(r8), INTENT(IN),  DIMENSION(nknot)       :: knots
    REAL(r8), INTENT(OUT), DIMENSION(plev+1,plev) :: VOPINTE
        
    ! Variables for b-spline basis function
    REAL(r8), DIMENSION(nlevbc,nbasis)            :: bspline_funs
    REAL(r8), DIMENSION(nquad*(plev+1),nbasis)    :: ebspline_funs_quad
    REAL(r8), DIMENSION(nquad*(plev+1),nbasis)    :: dbspline_funs_quad
    REAL(r8), DIMENSION(nquad*(plev+1),nbasis)    :: jbspline_deri_quad
    
    ! Variables for VFE
    REAL(r8)                                      :: ZDET  
    REAL(r8), DIMENSION(nlevbc-1,nbasis-1)        :: ZP
    REAL(r8), DIMENSION(nlevbc-1,nbasis-1)        :: ZS
    REAL(r8), DIMENSION(nlevbc-1,nbasis-1)        :: ZSINV
    
    REAL(r8), DIMENSION(nbasis-1,nbasis-1)        :: ZA
    REAL(r8), DIMENSION(nbasis-1,nbasis-1)        :: ZAINV
    REAL(r8), DIMENSION(nbasis-1,nbasis-1)        :: ZB
    REAL(r8), DIMENSION(nlevbc-1,nbasis-1)        :: ZINTE
    REAL(r8), DIMENSION(nlevbc-1,nbasis-1)        :: ZTEMP 
    
    INTEGER(i4)                                   :: i, j, k, ilev, icount, iquad
    REAL(r8)                                      :: deta, sum1, sum2
    REAL(r8), DIMENSION(nquad)                    :: quadx, quadw
    REAL(r8), DIMENSION(0:plev+1)                 :: ZD

    !
    !--------------------------------------------------------------------------------------------------------!
    ! STEP 1: Matrix P: ZP
    ! ZP(i,j) = BASISFUNC_j(eta_i)) from FE space of the function to be differenciated into physical space
    ! i: row->eta points, j: column->basis function
    !--------------------------------------------------------------------------------------------------------!
    ! 
    CALL COMPUTE_BASIS_FUNS_ETA(etaf,knots,bspline_funs)
    
    ZP(:,:)=0.0_r8
    DO i = 2, nlevbc
        DO j = 2, nbasis
            ZP(i-1,j-1) = bspline_funs(i,j)
        END DO        
    END DO
    !
    !--------------------------------------------------------------------------------------------------------!    
    ! STEP 2: Matrix S^-1: ZSINV
    ! ZS(i,j)=BASISDER_j(eta_i) to project the derivative from FE space to physical space 
    !--------------------------------------------------------------------------------------------------------!    
    !
    ZS(:,:) = ZP(:,:)
    
    ZDET = 0.0_r8
    CALL MINV(ZS,nlevbc-1,nbasis-1,ZDET,ZSINV)
    
    IF (debug_vfe) THEN 
        WRITE (*, '(A)') ''
        WRITE (*, '(A)') '  Matrix ZP'
        WRITE (10, '(A)') '  Matrix ZP'
        WRITE (*, '(A)') '  ----------------------------------------------------------------'
        DO i = 1, nlevbc-1
            WRITE (*, '(50F6.1)') (100.*ZP(i,j), j = 1, nbasis-1)
        END DO
    END IF
    !
    !--------------------------------------------------------------------------------------------------------!
    ! STEP 3: Mass matrix A^-1
    ! Compute A_ik(eta): integrate the product of basis functions at all quadrature points
    !--------------------------------------------------------------------------------------------------------!
    !
    ! Compute basis functions E_i(eta), E_k(eta), and dE_j(eta)/deta at eta quadrature points (quadx) 
    !
    ZA(:,:)=0.0_r8

    CALL DEFINE_QUADRATURE_SET(nquad,quadw,quadx)
    CALL COMPUTE_BASIS_FUNS_QUAD(etaf,knots,ebspline_funs_quad,jbspline_deri_quad)
    dbspline_funs_quad(:,:) = ebspline_funs_quad(:,:)
    
    dbspline_funs_quad(:,1) = 0.0_r8
    
    DO i = 2, nbasis
        DO k = 2, nbasis
            icount = 1
            ZA(i,k) = 0.0_r8
            DO ilev = 1, nlevbc-1
                deta = etaf(ilev) - etaf(ilev-1)                
                DO iquad = 1, nquad
                    ZA(i-1,k-1) = ZA(i-1,k-1) + 0.5d0*deta*quadw(iquad)*dbspline_funs_quad(icount,i)*dbspline_funs_quad(icount,k)
                    icount = icount + 1
                END DO                
            END DO
         END DO      
    END DO
    !
    ! Invert mass matrix: A^-1 
    !
    ZDET = 0.0_r8
    CALL MINV(ZA,nbasis-1,nbasis-1,ZDET,ZAINV)
    !ZAINV(:,:) = ZAINV(:,:)*NLEVEL
      
    IF (debug_vfe) THEN 
        WRITE (*, *    ) '  Determinant of A*NLEV=', ZDET
        WRITE (*, '(a)') ''
        WRITE (*, '(a)') '  Mass matrix: A'
        WRITE (*, '(a)') '  ----------------------------------------------------------------'
        DO i = 1, nbasis-1
         !  WRITE (*, '(50F20.10)') (ZA(i,j), j = 1, nbasis-1)
            WRITE (*, '(50F6.1)') (100.*ZA(i,j), j = 1, nbasis-1)
        END DO
        WRITE (*, '(a)') ''
        WRITE (*, '(a)') '  Mass matrix: AINV'
        WRITE (*, '(a)') '  ----------------------------------------------------------------'
        DO i = 1, nbasis-1
            WRITE (*, '(50F6.0)') (ZAINV(i,j), j = 1, nbasis-1)
        END DO
    END IF
    !
    !--------------------------------------------------------------------------------------------------------!
    ! STEP 4: Integral matrix B
    ! Compute B_ik(eta): integrate the product of basis functions at all quadrature points
    !--------------------------------------------------------------------------------------------------------! 
    !    
    ZB(:,:) = 0.0_r8
    DO j = 2, nbasis
        DO i = 1, nbasis-1      
            icount = 1
            sum1 = 0.0_r8 
            DO ilev = 1, j
                deta = etaf(ilev) - etaf(ilev-1)                
                DO iquad = 1, nquad
                    sum1 = sum1 + 0.5D0*deta*quadw(iquad)*ebspline_funs_quad(icount,i)
                    icount = icount + 1
                END DO
            END DO

            icount = 1
            sum2 = 0.0_r8        
            DO ilev = 1, nlevbc-1
                deta = etaf(ilev) - etaf(ilev-1)                
                DO iquad = 1, nquad
                    sum2 = sum2 + 0.5D0*deta*quadw(iquad)*dbspline_funs_quad(icount,j)*sum1
                    icount = icount + 1
                END DO
            END DO
            ZB(j-1,i) = sum2
        END DO
    END DO
    !  
    !
    ! Temporarily adopted a part of ZB from verfe1i.f90 for comparison ---------------------------------------------! 
    !
    ZD(:) = 0.0_r8
    DO i = 1, plev+1
        IF (i==1) THEN
            ZD(i) = etaf(1)
        ELSE IF(i==plev+1) THEN
            ZD(i) = 1.0_r8 - etaf(plev)
        ELSE
            ZD(i) = etaf(i) - etaf(i-1)
        END IF
    END DO
    
    DO i = 1, plev+1
        DO j = 1, plev+1
            IF (i==j+2) THEN
                ZB(j,i) = ZD(j+1)**2/24.0_r8
            ELSE IF(i==j+1) THEN
                ZB(j,i) = ZD(j)**2/8.0_r8 + ZD(j)*ZD(j+1)/4.0_r8 + ZD(j+1)**2/8.0_r8
            ELSE IF(i==j) THEN
                ZB(j,i) = ZD(j-1)*ZD(j)/4.0_r8 + 5.0_r8*ZD(j)**2/24.0_r8 + ZD(j+1)*(ZD(j-1)+ZD(j))/4.0_r8                
            END IF
        END DO
    END DO
    ZB(plev+1,1) = ZD(1)*ZD(plev+1)/4.0_r8
    ZB(plev  ,plev+1) = ZD(plev)**2/8.0_r8 + ZD(plev)*ZD(plev+1)/4.0_r8 + ZD(plev+1)**2/6.0_r8
    ZB(plev+1,plev+1) = ZD(plev)*ZD(plev+1)/4.0_r8 + ZD(plev+1)**2/3.0_r8
    !
    ! Temporarily adopted a part of ZA from verfe1i.f90 for comparison ---------------------------------------------! 
    !
    !
    IF (debug_vfe) THEN 
        WRITE (*, '(a)') ''
        WRITE (*, '(a)') '  Integral matrix: ZB*1E6'
        WRITE (*, '(a)') '  ----------------------------------------------------------------'
        DO i = 1, nbasis-1
            WRITE(*, '(50F6.1)') (1.E6*ZB(i,j), j = 1, nbasis-1)
        END DO
    END IF    
    !
    !--------------------------------------------------------------------------------------------------------!
    ! STEP 5:
    ! Compute integral operator ZINTE = P*A^-1*B*S^-1
    !--------------------------------------------------------------------------------------------------------!   
    !
    ! B*S^-1
    ZINTE(:,:) = 0.0_r8
    DO j = 1, nbasis-1
        DO i = 1, nbasis-1
            DO k = 1, nbasis-1
                ZINTE(i,j) = ZINTE(i,j) + ZB(i,k)*ZSINV(k,j)        
            END DO
        END DO
    END DO
    !
    ! A^-1*B*P^-1
    ZTEMP(:,:) = 0.0_r8
    DO j = 1, nbasis-1
        DO i = 1, nbasis-1
            DO k = 1, nbasis-1
                ZTEMP(i,j) = ZTEMP(i,j) + ZAINV(i,k)*ZINTE(k,j)
            END DO
        END DO
    END DO
    !
    ! ZDERI=P*A^-1*B*S^-1
    ZINTE(:,:) = 0.0D0
    DO j = 1, nbasis-1
        DO i = 1, nbasis-1
            DO k = 1, nbasis-1
                ZINTE(i,j) = ZINTE(i,j) + ZP(i,k)*ZTEMP(k,j)
            END DO
        END DO
    END DO
    !
    !--------------------------------------------------------------------------------------------------------!
    ! STEP 6:
    ! Reconstruction: ZINTE(N+1,N+1) -> VOPINTE(N+1,N)    
    !--------------------------------------------------------------------------------------------------------!  
    !
    VOPINTE(:,:) = 0.0_r8
    DO i = 1, plev
        DO j = 1, plev+1
            IF(i==1) THEN
                VOPINTE(j,i) = ZINTE(j,i) + ZINTE(j,i+1)
            ELSE
                VOPINTE(j,i) = ZINTE(j,i+1)
            END IF
        END DO
    END DO
    
    IF (debug_vfe) THEN 
        WRITE (*, '(a)') ''
        WRITE (*, '(a)') '  Integral operator matrix: VOPINTE(N+1,N)'
        WRITE (*, '(a)') '  ----------------------------------------------------------------'
        DO j = 1, plev+1
            WRITE (*, '(50F6.1)') (VOPINTE(j,i), i = 1, plev)
        END DO

        ZDET=0.0
        CALL MINV(ZINTE,nbasis-1,nbasis-1,ZDET,ZTEMP)
        WRITE (*, *    ) '  Determintant of ZINTE =', ZDET
        WRITE (*, '(a)') '  Inverse of matrix for integrals :       '
        WRITE (*, '(a)') '  ----------------------------------------------------------------'
        DO j = 1, plev+1
           WRITE (*, '(50F6.1)') (ZTEMP(j,i), i = 1, plev+1)
        ENDDO 

        WRITE(*,*) ' DERIVATIVES OF FUNCTION f(x)=1'
        DO j=1,plev
           sum1=0.0d0
           DO i=1,plev+1
              sum1=sum1+ZTEMP(j+1,i)
           ENDDO
           WRITE(*,*) ' LEVEL ',j,' DERIVATIVE=',sum1
        ENDDO
    END IF

    
END SUBROUTINE SUB_INTE_OPERATOR

END MODULE VFE_OPERATORS
