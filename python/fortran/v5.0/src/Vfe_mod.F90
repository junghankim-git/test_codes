#include "KIM.h"
!
!-----------------------------------------------------------------------------------------!
! VFE_MOD: main routine for VFE
!
!-----------------------------------------------------------------------------------------!
!
MODULE VFE_MOD

    USE KiapsBase,      ONLY : iulog => KIM_IU_LOG, r8 => KIM_REAL8_KIND, i4 => KIM_INT_KIND
    USE Dimensions,     ONLY : plev => nlev, plevp => nlevp
    USE KiapsParallel,  ONLY : KIM_Par

    IMPLICIT NONE

    PRIVATE

    LOGICAL               , PUBLIC                      :: debug_vfe=.FALSE.

    ! Variables for the De Boor algorithm
    INTEGER(i4), PARAMETER, PUBLIC                      :: basis_degree = 1      ! B-spline degree: 1=linear, 3=cubic
    INTEGER(i4), PARAMETER, PUBLIC                      :: deri_order = 1        ! First-order derivatives
    INTEGER(i4), PARAMETER, PUBLIC                      :: nquad = 5             ! Gaussian quadrature order

    INTEGER(i4), PARAMETER, PUBLIC                      :: basis_order = basis_degree + 1
    INTEGER(i4), PARAMETER, PUBLIC                      :: nlevbc = 50 + 2                  ! Nmber of eta levels (=nlev+2 for linear)
    INTEGER(i4), PARAMETER, PUBLIC                      :: nknot = nlevbc + 2*basis_degree  ! Number of knot components
    INTEGER(i4), PARAMETER, PUBLIC                      :: nbasis = nknot - basis_order     ! Number of basis functions = nknot - basis_order
    
    REAL(r8), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: VOPDERI
    REAL(r8), PUBLIC, DIMENSION(:,:), ALLOCATABLE, SAVE :: VOPINTE

    PUBLIC :: INI_VFE, FIN_VFE, SUB_VERINT, SUB_VERDER
    PUBLIC :: VFE_MASS_CORRECTION
    PUBLIC :: MINV, LUDCMP, LUBKSB
CONTAINS
!
!-----------------------------------------------------------------------------------------!
! SUBROUTINE INI_VFE
!
!-----------------------------------------------------------------------------------------!
!
SUBROUTINE INI_VFE(N)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N

    ALLOCATE(VOPDERI(N,N),VOPINTE(N+1,N))
END SUBROUTINE INI_VFE
!
!-----------------------------------------------------------------------------------------!
! SUBROUTINE FIN_VFE
!
!-----------------------------------------------------------------------------------------!
!
SUBROUTINE FIN_VFE
    IMPLICIT NONE

    DEALLOCATE(VOPDERI)
    DEALLOCATE(VOPINTE)
END SUBROUTINE FIN_VFE
!
!-----------------------------------------------------------------------------------------!
! SUBROUTINE VERIN
! Vertical Integration using matrix multiplication
!-----------------------------------------------------------------------------------------!
!
SUBROUTINE SUB_VERINT(N, VOPINTE, ZIN, ZOUT)
    IMPLICIT NONE
    
    INTEGER(i4), INTENT(IN)  :: N
    REAL(r8),    INTENT(IN)  :: VOPINTE(N+1,N), ZIN(N)
    REAL(r8),    INTENT(OUT) :: ZOUT(N+1)

    INTEGER(i4) :: I, J
    REAL(r8)    :: TOTAL
 
    ZOUT(:) = 0.0_r8
    DO J = 1, N
        DO I = 1, N+1
            ZOUT(I) = ZOUT(I) + VOPINTE(I,J)*ZIN(J)
        END DO
    END DO
END SUBROUTINE SUB_VERINT
!
!-----------------------------------------------------------------------------------------!
! SUBROUTINE SUB_VERDER
! Vertical Derivative using matrix multiplication
!-----------------------------------------------------------------------------------------!
!
SUBROUTINE SUB_VERDER(N, VOPDERI, ZIN, ZOUT)
    IMPLICIT NONE
    
    INTEGER(i4), INTENT(IN)  :: N
    REAL(r8),    INTENT(IN)  :: VOPDERI(N,N), ZIN(N)
    REAL(r8),    INTENT(OUT) :: ZOUT(N)

    INTEGER :: I,J
 
    ZOUT(:) = 0.0_r8
    DO J = 1, N
        DO I = 1, N
            ZOUT(I) = ZOUT(I) + VOPDERI(I,J)*ZIN(J)
        END DO
    END DO
END SUBROUTINE SUB_VERDER
!
!-----------------------------------------------------------------------------------------!
! SUBROUTINE VFE_MASS_CORRECTION
! mass(pressure) correction in VFE scheme
!-----------------------------------------------------------------------------------------!
!
SUBROUTINE VFE_MASS_CORRECTION (etai, hyai, hybi, hyam, hybm, dhyam, dhybm)

    IMPLICIT NONE

    REAL(r8),DIMENSION(plev+1), INTENT(IN ) :: etai, hyai, hybi
    REAL(r8),DIMENSION(plev),   INTENT(OUT) :: hyam, hybm, dhyam, dhybm

    REAL(r8)                   :: ZFAC
    REAL(r8),DIMENSION(plev)   :: DEL_ETAF, DEL_AF, DEL_BF, DVAF, DVBF
    REAL(r8),DIMENSION(plev+1) :: VAFC, VBFC, VAF, VBF
    INTEGER                    :: IREF,ITER, i,j,k


    IF (KIM_Par%isMasterProc) write(*,*)
    IF (KIM_Par%isMasterProc) write(*,*) ' *** Mass correction for VFE *** '

    VAFC(:) = 0.0_r8
    VBFC(:) = 0.0_r8
    DVAF(:) = 0.0_r8
    DVBF(:) = 0.0_r8
     VAF(:) = 0.0_r8
     VBF(:) = 0.0_r8

    ! 1. Define dB/dn in the usual FD way at middle levels
    DO k = 1, plev
        DEL_ETAF(k) = etai(k+1) -etai(k) ! dn
        DEL_AF(k) = hyai(k+1) - hyai(k) ! dA
        DEL_BF(k) = hybi(k+1) - hybi(k) ! dB
        DVBF(k) = DEL_BF(k)/DEL_ETAF(k) ! dB/dn
    END DO

    ! 2. Integrate dB/dn through the whole column, VBFC(plev+1) has total integral
    IF (KIM_Par%isMasterProc) write(*,*) 'Vcoord coefficients B : '

    CALL SUB_VERINT(plev, VOPINTE, DVBF, VBFC)
    IF (KIM_Par%isMasterProc) WRITE(*,'(A,ES20.10)') "  BEFORE-CORRECTION: Total_Int B =", VBFC(plev+1)

    ! 3. To obtain a corrected vector, B' and dB'/dn
    DO k = 1, plev
       DVBF(k)= DVBF(k)/VBFC(plev+1) ! dB'/dn
    END DO

    CALL SUB_VERINT(plev, VOPINTE, DVBF, VBF) ! B'
    DO k = 1, plev
       IF( VBF(k) < 0.0_r8 ) VBF(k) = 0.0_r8 ! to correct negative value
    END DO
    IF (KIM_Par%isMasterProc) WRITE(*,'(A,ES20.10)') "   AFTER-CORRECTION: Total_Int B =", VBF(plev+1)

    ! 4. To obtain a corrected vector, A' and dA'/dn: Use an iterative process to correct the original vector
    IF (KIM_Par%isMasterProc) write(*,*) 'Vcoord coefficients A : '
    ITER = 0
    DO
        ITER = ITER + 1

        DO k = 1, 1
            DVAF(k) = DEL_AF(k)/DEL_ETAF(k)
            IF (DEL_AF(k)<0.0_r8) IREF = k
        END DO

        DO k = 2, plev
            DVAF(k) = DEL_AF(k)/DEL_ETAF(k)
            IF (DEL_AF(k)<0.0_r8 .AND. DEL_AF(k-1)>=0.0_r8) IREF = k
        END DO

        CALL SUB_VERINT(plev, VOPINTE, DVAF, VAFC) ! A'

        ZFAC = -1.0_r8*VAFC(IREF)/(VAFC(plev+1)-VAFC(IREF))
        DO k = 1, IREF
            DEL_AF(k) = DEL_AF(k)/ZFAC
            DVAF(k) = DEL_AF(k)/DEL_ETAF(k) ! dA'/dn
        END DO
        IF (DABS(VAFC(plev+1))<1.D-14) THEN
            GOTO 111
        ELSE
            IF (KIM_Par%isMasterProc)  WRITE(*,'(A,2I5,ES20.10)') "  BEFORE-CORRECTION: IREF,ITER,Total_Int A =", IREF, ITER, VAFC(plev+1)
        END IF

        IF (ITER>20) THEN
           IF (KIM_Par%isMasterProc)  WRITE (*,*) ' STOP: The number of iteration exceeds ITER=20: HybridMassVFE.F90 '
            STOP
        END IF
    END DO

111 CALL SUB_VERINT(plev, VOPINTE, DVAF, VAF) ! A'    
    DO k = 1, plev
       IF( VAF(k) < 0.0_r8 ) VAF(k) = 0.0_r8 ! to correct negative value
    END DO
    IF (KIM_Par%isMasterProc) WRITE(*,'(A,2I5,ES20.10)') "   AFTER-CORRECTION: IREF,ITER,Total_Int A =", IREF, ITER, VAF(plev+1)
    IF (KIM_Par%isMasterProc) WRITE(*,*)

    hyam(:) = VAF(1:plev)
    hybm(:) = VBF(1:plev)
    dhyam(:)= DVAF(:)
    dhybm(:)= DVBF(:)

 
END SUBROUTINE VFE_MASS_CORRECTION
 
!
!---------------------------------------------------------------------!
!     Subroutines for matrix inversion using LU Decomposition
!---------------------------------------------------------------------!

  !===================================================================!
  !  !---- Further Comments ------------------------------------------!
  !
  !     A*X= B
  !
  !   CALL LUDCMP(A,N,M,INDX,D)
  !   CALL LUBKSB(A,N,M,INDX,B) ! The answer X is returned in B.
  !                              ! Original Matrix A is destroyed.
  !                              ! For further calculation with
  !                              ! different B : CALL LUBKSB once more.
  !
  !---- Further Comments ----------------------------------------------!
  !====================================================================!
  SUBROUTINE MINV( A, N, M, D, Y )
  !--------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER , INTENT(IN ) :: N,M
    REAL(r8), INTENT(IN ) :: A(M,M)
    REAL(r8), INTENT(OUT) :: Y(M,M)
    REAL(r8), INTENT(OUT) :: D
    REAL(r8) :: AA(M,M)
    INTEGER  :: INDX(M),i,j
  !--------------------------------------------------------------------!

      AA(:,:)= A(:,:)  ! to save A matrix

        DO i=1,N
           DO j=1,N
              Y(i,j)= 0.
           ENDDO
              Y(i,i)= 1.D0
        ENDDO

           CALL LUDCMP ( AA, N, M, INDX, D )
        DO j=1,N
           D= D*AA(j,j) ! determinant of the original matrix A
        ENDDO

        DO j=1,N
           CALL LUBKSB ( AA, N, M, INDX, Y(1,j) )  ! Y= Inverse Matrix of A
        ENDDO                                       ! A is destroyed. 

  !--------------------------------------------------------------------!
  END SUBROUTINE MINV
  !====================================================================!
  !
  !====================================================================!
  SUBROUTINE LUDCMP ( A, N, M, INDX, D )
  !--------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER , INTENT(IN   ):: N,M
    REAL(r8), INTENT(INOUT):: A(M,M)
    REAL(r8), INTENT(INOUT):: D
    REAL(r8) :: aamax, sum, dum, tiny, VV(999)
    INTEGER  :: INDX(M),i,j,k,imax
    DATA TINY /1.0D-20/
  !--------------------------------------------------------------------!
  !
       IF(N.GT.999) STOP ' Too small work-array size ! '
       D= 1.D0
       DO i=1,N
             aamax= 0.
          DO j=1,N
             IF(dabs(a(i,j)).gt.aamax) aamax= dabs(a(i,j))
          ENDDO
             IF(aamax.eq.0) STOP 'Singular Matrix in ludcmp'
             vv(i)=1.D0/aamax
       END DO
!
       DO j=1,N ! 19
          DO i=1,j-1
                sum= A(i,j)
             DO k=1,i-1
                sum= sum-A(i,k)*A(k,j)
             ENDDO
             A(i,j)= sum
          ENDDO
          aamax= 0.
          DO i=j,N  ! 16
                sum= A(i,j)
             DO k=1,j-1
                sum= sum-A(i,k)*A(k,j)
             ENDDO
             A(i,j)= sum
             dum = vv(i)*dabs(sum)
             IF(dum.ge.aamax)THEN
               imax = i
               aamax= dum
             ENDIF
          ENDDO  ! 16

!            
          IF(j.ne.imax)THEN
            DO k=1,N
               dum=A(imax,k)
               A(imax,k)= A(j,k)
               A(j,k)= dum
            ENDDO
            D= -D
            vv(imax)= vv(j)
            ENDIF
            indx(j)= imax
            IF(a(j,j).eq.0) A(j,j)= TINY
            IF(j.ne.N) THEN
               dum= 1.d0/A(j,j)
               DO i=j+1,N
                  A(i,j)= A(i,j)*dum
               ENDDO
            ENDIF
       ENDDO   !---19

  !--------------------------------------------------------------------!
  END SUBROUTINE LUDCMP
  !====================================================================!
  !
  !
  !====================================================================!
  SUBROUTINE LUBKSB ( A, N, M, INDX, B )
  !--------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER , INTENT(IN   ):: N,M
    REAL(r8), INTENT(INOUT):: A(M,M)
    REAL(r8), INTENT(INOUT):: B(N)
    REAL(r8) :: sum,dum,tiny
    INTEGER  :: INDX(M),i,j,k,ii,LL
  !--------------------------------------------------------------------!

       ii=0
       DO i=1,N
             LL= indx(i)
            sum= B(LL)
          b(LL)= B(i)
          IF(ii.ne.0) THEN
             DO j=ii,i-1
                sum= sum - A(i,j)*B(j)
             ENDDO
          ELSEIF(sum.ne.0)THEN
             ii= i
          ENDIF
             B(i)= sum
       ENDDO

       DO i=N,1,-1
             sum= B(i)
          DO j=i+1,N
             sum= sum-A(i,j)*B(j)
          ENDDO
             B(i)= sum/A(i,i)
       ENDDO

  !--------------------------------------------------------------------!
  END SUBROUTINE LUBKSB
  !====================================================================!


END MODULE VFE_MOD
