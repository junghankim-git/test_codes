PROGRAM test
 USE Kinds,   ONLY: i4, l4, r4, r8
 USE Matrix,  ONLY: InvMatrix, CheckIdentity
 IMPLICIT NONE

 REAL(r8), PARAMETER      :: pi = ACOS(-1.0_r8)
 INTEGER(i4), PARAMETER   :: n = 2
 REAL(r8), DIMENSION(2,2) :: A, invA
 
 INTEGER(i4) :: i


  A(1,1) =  2.3_r8
  A(2,1) =  1.3_r8
  A(1,2) =  6.3_r8
  A(2,2) = -4.3_r8

  CALL InvMatrix(n, A, invA)
  PRINT *, CheckIdentity(n, A, invA)


END PROGRAM test
