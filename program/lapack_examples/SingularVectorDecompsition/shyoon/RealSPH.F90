!-------------------------------------------------------------------------------
!
!  MODULE RealSPH
!
!> @brief
!>  - Fortran module to compute real spherical harmonics at node points
!>  - given on the unit sphere
!>
!> @date 01MAY2014
!>  - IS SONG : First crack
!
!-------------------------------------------------------------------------------

MODULE RealSPH

!  USE Base,              ONLY: r8 => KGM_REAL8_KIND  ! KGM_REAL8_KIND = SELECTED_REAL_KIND(12)
!  USE PhysicalConstants, ONLY: pi => DD_PI           ! DD_PI = 3.141592653589793238462643383279_r8

  IMPLICIT NONE

  PRIVATE

  INTEGER,  PUBLIC,  PARAMETER :: r8 = SELECTED_REAL_KIND(12)
  REAL(r8), PUBLIC,  PARAMETER :: pi = 3.141592653589793238462643383279_r8

  INTEGER,  PUBLIC,  PARAMETER :: M_MAX = 8000
  INTEGER,  PUBLIC,  PARAMETER :: L_MAX = M_MAX
  REAL(r8), PRIVATE, PARAMETER :: MISSING_VAL = 1.e+100_r8

  REAL(r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: alm
  REAL(r8), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: blm


  PUBLIC :: IniRecurCoef
  PUBLIC :: ComputSPH
! PUBLIC :: ConstructSPHMatrix

CONTAINS

  SUBROUTINE IniRecurCoef

    IMPLICIT NONE

    REAL(r8), DIMENSION(-M_MAX:M_MAX) :: amm
    REAL(r8), DIMENSION(-M_MAX:M_MAX) :: prod

    INTEGER :: l, m, k

!   Begins
!   ------

    IF (.NOT. ALLOCATED(alm)) ALLOCATE(alm(0:L_MAX,-M_MAX:M_MAX))
    IF (.NOT. ALLOCATED(blm)) ALLOCATE(blm(0:L_MAX,-M_MAX:M_MAX))

!   amm

    DO m = -M_MAX, M_MAX
      prod(m) = 1._r8
    END DO
    DO m = 0, M_MAX
      DO k = 1, ABS(m)
        prod(m) = prod(m)*(2._r8*DBLE(k)+1._r8)/(2._r8*DBLE(k))
      END DO
    END DO
    DO m = -M_MAX, -1
      prod(m) = prod(-m)
    END DO
    DO m = -M_MAX, M_MAX
      amm(m) = ((-1)**ABS(m))*SQRT(1._r8/(4._r8*pi)*prod(m))
    END DO

!   alm

    DO l = 1, L_MAX
      alm(l,0) = SQRT((4._r8*DBLE(l)**2-1)/DBLE(l)**2)
    END DO
    DO m = 0, M_MAX
      alm(m,m) = amm(m)
    END DO
    DO m = 1, M_MAX
      DO l = 0, m-1
        alm(l,m) = MISSING_VAL
      END DO
      DO l = m+1, L_MAX
        alm(l,m) = SQRT((4._r8*DBLE(l)**2-1)/(DBLE(l)**2-DBLE(m)**2))
      END DO
    END DO
    DO m = -M_MAX,-1
      DO l = 0, L_MAX
        alm(l,m) = alm(l,-m)
      END DO
    END DO

!   blm

    DO l = 0, 1
      blm(l,0) = MISSING_VAL
    END DO
    DO l = 2, L_MAX
      blm(l,0) = -SQRT((2._r8*DBLE(l)+1._r8)/(2._r8*DBLE(l)-3._r8)*  &
                       (DBLE(l-1)**2/DBLE(l)**2))
    END DO
    DO m = 1, M_MAX
      DO l = 0, m+1
        IF (l <= L_MAX) blm(l,m) = MISSING_VAL
      END DO
      DO l = m+2, L_MAX
        blm(l,m) = -SQRT((2._r8*DBLE(l)+1._r8)/(2._r8*DBLE(l)-3._r8)*  &
                         (DBLE(l-1)**2-DBLE(m)**2)/(DBLE(l)**2-DBLE(m)**2))
      END DO
    END DO 
    DO m = -M_MAX, -1
      DO l = 0, L_MAX
        blm(l,m) = blm(l,-m)
      END DO
    END DO

    RETURN
  END SUBROUTINE IniRecurCoef

  SUBROUTINE ComputSPH(l, m, n, lat, long, cylm, sylm)

    IMPLICIT NONE

    INTEGER               , INTENT(IN)    :: l, m, n
    REAL(r8), DIMENSION(n), INTENT(IN)    :: lat, long  ! in radians
    REAL(r8), DIMENSION(n), INTENT(INOUT) :: cylm, sylm

    REAL(r8), DIMENSION(n,3) :: lgdrp
    REAL(r8), DIMENSION(n) :: cosl, sinl
    INTEGER :: i, j, mp1, mp2

!   Begins
!   ------

    IF (l < 0) THEN
      WRITE(6,'(A)') '(ComputeSPH): Negative degree (l)'
      STOP
    END IF
    IF (l > L_MAX) THEN
      WRITE(6,'(A)') '(ComputeSPH): Degree (l) larger than L_MAX'
      STOP
    END IF
    IF (m < 0) THEN
      WRITE(6,'(A)') '(ComputeSPH): Negative order (m)'
      STOP
    END IF
    IF (m > M_MAX ) THEN
      WRITE(6,'(A)') '(ComputeSPH): Order (m) larger than M_MAX'
      STOP
    END IF
    IF (alm(l,m) == MISSING_VAL) THEN
      WRITE(6,'(A)') '(ComputeSPH): Order (m) and degree (l) are '//  &
                     'inconsistent with M_MAX'
      STOP
    END IF

    DO i = 1, n
      sinl (i) = SIN(lat(i))
    END DO
    DO i = 1, n
      cosl(i) = (1._r8-sinl(i)*sinl(i))**(ABS(m)*0.5_r8)
    END DO

    IF (m <= M_MAX-2) THEN
      DO i = 1, n
        lgdrp(i,1) = alm(m  ,m)*cosl(i)
      END DO
      DO i = 1, n
        lgdrp(i,2) = alm(m+1,m)*sinl(i)*lgdrp(i,1)
      END DO
      DO i = 1, n
        lgdrp(i,3) = alm(m+2,m)*sinl(i)*lgdrp(i,2) + blm(m+2,m)*lgdrp(i,1)
      END DO 
      j = m+3
      DO WHILE (j <= l)
        DO i = 1, n
          lgdrp(i,1) = lgdrp(i,2)
        END DO
        DO i = 1, n
          lgdrp(i,2) = lgdrp(i,3)
        END DO
        DO i = 1, n
          lgdrp(i,3) = alm(j,m)*sinl(i)*lgdrp(i,2) + blm(j,m)*lgdrp(i,1)
        END DO
        j = j+1
      END DO
    ELSE IF (m == M_MAX-1) THEN
      DO i = 1, n
        lgdrp(i,1) = alm(m  ,m)*cosl(i)
      END DO
      DO i = 1, n
        lgdrp(i,2) = alm(m+1,m)*sinl (i)*lgdrp(i,1)
      END DO
      DO i = 1, n
        lgdrp(i,3) = lgdrp(i,2)
      END DO
    ELSE IF (m == M_MAX) THEN
      DO i = 1, n
        lgdrp(i,1) = alm(m  ,m)*cosl(i)
      END DO
      DO i = 1, n
        lgdrp(i,2) = lgdrp(i,1)
      END DO
      DO i = 1, n
        lgdrp(i,3) = lgdrp(i,1)
      END DO
    END IF

    IF (l == m) THEN
      DO i = 1, n
        lgdrp(i,3) = lgdrp(i,1)
      END DO
    ELSE IF (l == m+1) THEN
      DO i = 1, n
        lgdrp(i,3) = lgdrp(i,2)
      END DO
    END IF

    IF (m == 0) THEN
      DO i = 1, n
        cylm(i) = lgdrp(i,3)
      END DO
      DO i = 1, n
        sylm(i) = MISSING_VAL
      END DO
    ELSE
      DO i = 1, n
        cylm(i) = SQRT(2._r8)*((-1)**m)*lgdrp(i,3)*COS(DBLE(m)*long(i))
      END DO
      DO i = 1, n
        sylm(i) = SQRT(2._r8)*((-1)**m)*lgdrp(i,3)*SIN(DBLE(m)*long(i))
      END DO
    END IF

    RETURN
  END SUBROUTINE ComputSPH

END MODULE RealSPH
