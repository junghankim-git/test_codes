PROGRAM LeastSquare
IMPLICIT NONE
INTEGER   :: i,j,k,m,n,INFO
!PARAMETER (m=3,n=3)
PARAMETER (m=2,n=3)
REAL b(m), y(n), A(m,n), WORK(2*m*n)


b(1)=1.0;b(2)=0;b(3)=1.0
A(1,1)=1.0;A(1,2)=-1.0;A(1,3)=1.0
A(2,1)=1.0;A(2,2)= 0.0;A(2,3)=0.0
A(3,1)=1.0;A(3,2)= 1.0;A(3,3)=1.0

CALL SGELS('N', m, n, 1, A, m, b, m, WORK,2*m*n,INFO)
     
PRINT *, 'b1 = ', b(1)
PRINT *, 'b2 = ', b(2)
PRINT *, 'b3 = ', b(3)

PRINT *, 'A='
DO i=1,m
  PRINT *, (A(i,j), j=1,n)
ENDDO
STOP
END
