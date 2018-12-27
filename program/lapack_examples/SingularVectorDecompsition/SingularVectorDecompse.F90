PROGRAM Test
!  USE Kinds,         ONLY : i4, r8
!  USE ControlMatrix, ONLY : PrintMatrix
  IMPLICIT NONE
  integer, parameter :: i4 = 4, r8 = 8

!  ! KINDs
!    INTEGER, PARAMETER                 :: i4 = 4, r8 = 8
  ! Variables
    INTEGER(i4), PARAMETER             :: m=8, n=8
    REAL(r8), DIMENSION(m,n)           :: A, Res, Sigma, tmp
    REAL(r8), DIMENSION(m,m)           :: U, IdenM
    REAL(r8), DIMENSION(n,n)           :: VT, IdenN
    REAL(r8), DIMENSION(MIN(n,m))      :: S
    INTEGER(i4), PARAMETER             :: LWORK=5*m*n
    REAL(r8), DIMENSION(LWORK)         :: WORK


    INTEGER(i4)                        :: INFO
    INTEGER(i4)                        :: i, j, k
    !

    PRINT *, ' '
    PRINT *, '########  START:   Singular Vector Decomposition ########'
    PRINT *, ' '

    ! INPUT
!    DATA A / &
!            1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
!            0.0_r8, 0.0_r8, 0.0_r8, 4.0_r8, &
!            0.0_r8, 3.0_r8, 0.0_r8, 0.0_r8, &
!            0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
!            2.0_r8, 0.0_r8, 0.0_r8, 0.0_r8 /

    PRINT *, '# START: Solve A matrix..'
    CALL GetMatrix(A)
    PRINT *, '# END  : Solve A matrix..'

    CALL PrintMatrix('# Check Matrix (A)           = ', A, m, n)

    ! Singular Vector Decomposition
    CALL DGESVD('A', 'A', m, n, A, m, S, U, m, VT, n, WORK, LWORK, INFO)
    IF(INFO .NE. 0) THEN
      PRINT *, 'DGESVD: Check arguments.... Program will be stoped.'
      STOP
    END IF


    CALL PrintMatrix('# Singular Values (Sigma)    = ', S, 1, MIN(m,n))

    CALL PrintMatrix('# Left  Singular Vector (U)  = ', U, m, m)

    CALL PrintMatrix('# Right Singular Vector (VT) = ', VT, n, n)

    IdenM = 0.0_r8
    DO j = 1, m
      DO i = 1, m
        DO k = 1, m
         IdenM(i,j) = IdenM(i,j) + U(i,k)*U(j,k)
        ENDDO
      ENDDO
    ENDDO

    CALL PrintMatrix('# Check Identity (U     )    = ', IdenM, m, m)

    IdenN = 0.0_r8
    DO j = 1, n
      DO i = 1, n
        DO k = 1, n
         IdenN(i,j) = IdenN(i,j) + VT(i,k)*VT(j,k)
        ENDDO
      ENDDO
    ENDDO

    CALL PrintMatrix('# Check Identity (VT)        = ', IdenN, n, n)

    DO j = 1, n
      DO i = 1, m
        IF(i == j) THEN
           Sigma(i,j) = S(i)
        ELSE
           Sigma(i,j) = 0.0_r8
        END IF
      ENDDO
    ENDDO
    tmp = 0.0_r8
    DO j = 1, n
      DO i = 1, m
        DO k = 1, m
         tmp(i,j) = tmp(i,j) + U(i,k)*Sigma(k,j)
        ENDDO
      ENDDO
    ENDDO
    Res = 0.0_r8
    DO j = 1, n
      DO i = 1, m
        DO k = 1, n
         Res(i,j) = Res(i,j) + tmp(i,k)*VT(k,j)
        ENDDO
      ENDDO
    ENDDO
    CALL PrintMatrix('# Check Matrix (U*S*VT)      = ', Res, m, n)

    PRINT *, ' '
    PRINT *, '########  END  :   Singular Vector Decomposition ########'
    PRINT *, ' '

END PROGRAM Test




SUBROUTINE PrintMatrix(String_in, A_in, m_in, n_in)
   IMPLICIT NONE
   CHARACTER(*), INTENT(IN)                   :: String_in
   INTEGER(4), INTENT(IN)                     :: m_in, n_in
   REAL(8), DIMENSION(m_in, n_in), INTENT(IN) :: A_in
 
   ! local 
   INTEGER(4)      :: i, j
   CHARACTER(20)   :: frmt
 
   WRITE(frmt, '(A,I02,A)') '(', n_in, '(F10.5,X))'
 
   PRINT *, String_in
   DO i = 1, m_in
      WRITE(*,frmt) (A_in(i,j), j=1,n_in)
   END DO
   PRINT *, ' '

END SUBROUTINE PrintMatrix




SUBROUTINE GetMatrix(mat)

USE RealSPH

IMPLICIT NONE

integer, parameter        :: n=8
real(r8), dimension(n,n), INTENT(INOUT)  :: mat
!INTEGER, PARAMETER                 :: i4 = 4, r8 = 8
integer                   :: l,m,i,j
integer, dimension(3)     :: ll, mm
real(r8), dimension(n)    :: lat, long
real(r8), dimension(n)    :: cylm,sylm
real(r8), dimension(2,n)  :: CY
real(r8), dimension(3,n)  :: SY
real(r8), dimension(n)    :: Y00, Y10, Y20

!set lat, long


 lat(1)  = asin(1._r8/(sqrt(3._r8)))
 long(1) = 0._r8 
 lat(2)  = -asin(1._r8/(sqrt(3._r8)))
 long(2) = 0._r8

 lat(3)  = asin(1._r8/(sqrt(3._r8)))
 long(3) = pi/2._r8

 lat(4)  = -asin(1._r8/(sqrt(3._r8)))
 long(4) = pi/2._r8

 lat(5)  = asin(1._r8/(sqrt(3._r8)))
 long(5) = pi

 lat(6)  = -asin(1._r8/(sqrt(3._r8)))
 long(6) = pi

 lat(7)  = asin(1._r8/(sqrt(3._r8)))
 long(7) = 3._r8*pi/2._r8

 lat(8)  = -asin(1._r8/(sqrt(3._r8)))
 long(8) = 3._r8*pi/2._r8

!compute Y00, Y10, Y20

do i=1,n
 Y00(i) = 1._r8/(2._r8*(sqrt(pi)))
 Y10(i) = (1._r8/2._r8)*(sqrt(3._r8/pi))*(sin(lat(i)))
 Y20(i) = (1._r8/4._r8)*(sqrt(5._r8/pi))*((3._r8*(sin(lat(i)))**2)-1._r8)
enddo

!read grid

open(10,file='l_m.txt')
do i=1,3
  read(10,*) ll(i),mm(i)
  
  !compute SPH
  
  l=ll(i)
  m=mm(i)
  CALL IniRecurCoef
  CALL ComputSPH (l, m, n, lat, long, cylm, sylm)
  
  if (l==1.and.m==1) then
  CY(1,:)= cylm(:)
  SY(1,:)= sylm(:)
  elseif(l==2.and.m==1) then
  CY(2,:)= cylm(:)
  SY(2,:)= sylm(:)
  elseif(l==2.and.m==2) then
  SY(3,:)= sylm(:)
  endif
enddo

!write matrix

do j=1,n
 do i=1,n
 if (j==1)then
  mat(i,j)= Y00(i)
 elseif (j==2) then
  mat(i,j)= SY(1,i)
 elseif (j==3) then
  mat(i,j)= Y10(i)
 elseif (j==4) then
  mat(i,j)= CY(1,i)
 elseif (j==5) then
  mat(i,j)= SY(3,i)
 elseif (j==6) then
  mat(i,j)= SY(2,i)
 elseif (j==7) then
  mat(i,j)= Y20(i)
 elseif (j==8) then
  mat(i,j)= CY(2,i)
 endif
 enddo
enddo

OPEN(unit=21, file='out.txt', form='formatted', status='unknown')

do j = 1, n
do i = 1, n
WRITE(21,*) mat(i,j)
enddo
enddo

CLOSE(21)

END SUBROUTINE GetMatrix
