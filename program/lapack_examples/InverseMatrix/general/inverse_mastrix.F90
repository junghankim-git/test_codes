PROGRAM LinearEqu
IMPLICIT NONE
INTEGER i,j,k,n,INFO
PARAMETER (n=2)
INTEGER IPIVOT(n)
REAL b(n), A(n,n), P(n,n), tmpP(n), WORK(n)

  b(1)=1.0;b(2)=1.0
  A(1,1)=1.0;A(1,2)=1.0;A(2,1)=2.0;A(2,2)=3.0
  P(1,1)=1.0;P(1,2)=0.0;P(2,1)=0.0;P(2,2)=1.0
  PRINT *, 'A ='
  DO i=1,n
    PRINT *, (A(i,j), j=1,n)
  ENDDO
  
  ! Factorization
  CALL SGETRF(n, n, A, n, IPIVOT, INFO)
  
  PRINT *, 'A(L-I+U) ='
  DO i=1,n
    PRINT *, (A(i,j), j=1,n)
  ENDDO
  PRINT *, 'PIVOT(1) = ', IPIVOT(1)
  PRINT *, 'PIVOT(2) = ', IPIVOT(2)
  DO i=1,n
    IF(IPIVOT(i) .ne. i) THEN
      DO j=1,n
        tmpP(j) = P(i,j)
        P(i,j)=P(IPIVOT(i),j)
        P(IPIVOT(i),j)=tmpP(j)
      ENDDO
    ENDIF
  ENDDO
  PRINT *, 'P = ', P
  
  ! Solve AX = B
!  CALL SGETRS('N', n, 1, A, n, IPIVOT, B, n, INFO)
  
  PRINT *, 'x = ', b
  
  ! Solve A^{-1}
  CALL SGETRI(n, A, n, IPIVOT, WORK, n, INFO)
  PRINT *, 'Inverse of A = '
  DO i=1,n
    PRINT *, (A(i,j), j=1,n)
  ENDDO
STOP
END
