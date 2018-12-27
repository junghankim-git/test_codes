PROGRAM LinearEqu
IMPLICIT NONE
INTEGER i,j,k,n,INFO
PARAMETER (n=2)
INTEGER IPIVOT(n)
REAL b(n,3), A(n,n), P(n,n), tmpP(n)



!b(1,1)=1.0;b(2,1)=1.0
b(1,1)=3.0;b(2,1)=2.0
b(1,2)=-1.0;b(2,2)=-1.0
b(1,3)=2.0;b(2,3)=2.0
A(1,1)=1.0;A(1,2)=1.0
A(2,1)=2.0;A(2,2)=3.0
!
!P(1,1)=1.0;P(1,2)=0.0;P(2,1)=0.0;P(2,2)=1.0


PRINT *, 'b = ', b(:,1)
PRINT *, 'b = ', b(:,2)
PRINT *, 'b = ', b(:,3)
PRINT *, 'A='
DO i=1,n
  PRINT *, (A(i,j), j=1,n)
ENDDO

CALL SGESV(n, 3, A, n, IPIVOT, b, n, INFO)

PRINT *, 'x = ', b(:,1)
PRINT *, 'x = ', b(:,2)
PRINT *, 'x = ', b(:,3)



!PRINT *, 'L = '
!PRINT *, '[    1.000000     ,    0.000000     ]'
!PRINT *, '[',A(2,1), ',    1.000000     ]'
!PRINT *, 'U = '
!PRINT *, '[',A(1,1), ',',A(1,2),']'
!PRINT *, '[    0.000000     , ', A(2,2),']'
!PRINT *, 'PIVOT(1) = ', IPIVOT(1)
!PRINT *, 'PIVOT(2) = ', IPIVOT(2)
!DO i=1,n
!  IF(IPIVOT(i) .ne. i) THEN
!    DO j=1,n 
!      tmpP(j) = P(i,j)
!      P(i,j)=P(IPIVOT(i),j)
!      P(IPIVOT(i),j)=tmpP(j)
!    ENDDO
!  ENDIF
!ENDDO
!PRINT *, 'P = ', P
!STOP
END
