PROGRAM Test
  USE Kinds,         ONLY : i4, r8
  USE GMRES_Solver
  USE ReproSum,          ONLY : SCS, DDPDD
  IMPLICIT NONE

  ! local
    INTEGER(i4), PARAMETER             :: n = 6
    REAL(r8), DIMENSION(n)             :: eta_m
    REAL(r8), DIMENSION(n,n)           :: LU_m, LUinv_m
    REAL(r8), DIMENSION(n,n)           :: Ainv_m, tmp_m
    REAL(r8), DIMENSION(n)             :: ipivot_m
    REAL(r8), DIMENSION(n)             :: work_m
! START
    REAL(r8), DIMENSION(2*n*n)         :: wwork_m
    REAL(r8), DIMENSION(n)             :: polyArg_m, y_m, check_m
    REAL(r8)                           :: tsum, errval
! END
    INTEGER(i4)                        :: ip, ie, ilev, info, i, j
    !

    PRINT *, '########  Start    FEM      Setting  ########'
    PRINT *, ' '


    DO i = 1, n
      eta_m(i)       = 1.0_r8/DBLE(n-1)*DBLE(i-1)
      PRINT *, 'eta = ', eta_m(i)
    END DO
    DO j = 1, n
      DO i = 1, n
        LU_m(i, j)   = eta_m(i)**DBLE(j-1)
        !CALL RANDOM_NUMBER(LU_m(i, j))
      END DO
    END DO
    LUinv_m(:, :) = LU_m(:, :)

    ! Factorization & Inverse Matrix
    CALL DGETRF(n, n, LUinv_m, n, ipivot_m, info)
    IF(info .NE. 0) STOP
    CALL DGETRI(n, LUinv_m, n, ipivot_m, work_m, n, info)
    IF(info .NE. 0) STOP


!    DO j = 1, n
!      DO i = 1, n
!        !LUinv_m(i, j) = Ainv_m(i, j)
!        LUinv_m(i, j) = LU_m(i, j)
!      END DO
!    END DO


    tmp_m(:,:) = 0.0_r8
    errval = 0.0_r8
    DO j = 1, n
     DO i = 1, n
      DO ip = 1, n
!        tmp_m(i,j) = tmp_m(i,j) + LU_m(i,ip)*LUinv_m(ip,j)
        CALL DDPDD(tmp_m(i,j), errval, LU_m(i,ip)*LUinv_m(ip,j))
      END DO
      IF(i == j) tmp_m(i,j) = tmp_m(i,j) - 1.0_r8
     END DO
    END DO
    PRINT *, '### Check DGETRF + DGETRI : '
    PRINT *, ' '
    PRINT *, 'A = '
    DO j = 1, n
      WRITE(*,'(10(D13.6,X))') LU_m(j, :)
    END DO
    PRINT *, ' '
    PRINT *, 'A^-1 = '
    DO j = 1, n
      WRITE(*,'(10(D13.6,X))') LUinv_m(j, :)
    END DO
    PRINT *, ' '
    PRINT *, 'AA^-1 - I = '
    DO j = 1, n
      WRITE(*,'(10(D13.6,X))') tmp_m(j, :)
    END DO
    PRINT *, ' '


!!! Check DGESV
    DO j = 1, n
      DO i = 1, n
        LUinv_m(i, j) = LU_m(i, j)
      END DO
      !polyArg_m(j) = eta_m(j)**2.0_r8
      CALL RANDOM_NUMBER(polyArg_m(j))
    END DO
    y_m(:) = polyArg_m(:)
    CALL DGELS('n', n, n, 1, LUinv_m, n, polyArg_m, n, wwork_m, 2*n*n, info)
!    CALL DGESV(n, 1, LUinv_m, n, ipivot_m, polyArg_m, n, info)
    IF(info .NE. 0) STOP
    check_m(:) = 0.0_r8
    errval = 0.0_r8
    DO j = 1, n
     DO i = 1, n
!        check_m(i) = check_m(i) + LU_m(i,j)*polyArg_m(j)
        CALL DDPDD(check_m(i), errval, LU_m(i,j)*polyArg_m(j))
     END DO
    END DO
    check_m(:) = check_m(:) - y_m(:)
    PRINT *, '### Check DGESV : '
    PRINT *, ' '
    PRINT *, 'y = '
      WRITE(*,'(10(D13.6,X))') y_m(:)
    PRINT *, ' '
    PRINT *, 'x = '
      WRITE(*,'(10(D13.6,X))') polyArg_m(:)
    PRINT *, ' '
    PRINT *, 'Ax - y = '
      WRITE(*,'(10(D13.6,X))') check_m(:)
    PRINT *, ' '


!!! Check GMRES
    DO j = 1, n
      DO i = 1, n
        LUinv_m(i, j) = LU_m(i, j)
      END DO
      !polyArg_m(j) = eta_m(j)**2.0_r8
!      CALL RANDOM_NUMBER(polyArg_m(j))
    END DO
!    y_m(:) = polyArg_m(:)
!    CALL DGELS('n', n, n, 1, LUinv_m, n, polyArg_m, n, wwork_m, 2*n*n, info)
!    CALL DGESV(n, 1, LUinv_m, n, ipivot_m, polyArg_m, n, info)
    PRINT *, 'Start GMRES...'
    CALL SolveGMRES(n, LUinv_m, y_m, polyArg_m)
    PRINT *, 'End   GMRES...'
    IF(info .NE. 0) STOP
    check_m(:) = 0.0_r8
    errval = 0.0_r8
    DO j = 1, n
     DO i = 1, n
!        check_m(i) = check_m(i) + LU_m(i,j)*polyArg_m(j)
        CALL DDPDD(check_m(i), errval, LU_m(i,j)*polyArg_m(j))
     END DO
    END DO
    check_m(:) = check_m(:) - y_m(:)
    PRINT *, '### Check GMRES : '
    PRINT *, ' '
    PRINT *, 'y = '
      WRITE(*,'(10(D13.6,X))') y_m(:)
    PRINT *, ' '
    PRINT *, 'x = '
      WRITE(*,'(10(D13.6,X))') polyArg_m(:)
    PRINT *, ' '
    PRINT *, 'Ax - y = '
      WRITE(*,'(10(D13.6,X))') check_m(:)
    PRINT *, ' '



END PROGRAM Test
