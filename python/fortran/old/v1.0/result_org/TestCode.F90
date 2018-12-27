MODULE mymod
END MODULE mymod
!>  - JH KIM  : Added to "ghost_exchangeV" subroutine. (CAM-SE 5.3)
!-------------------------------------------------------------------------------
#include "KIM.h"
MODULE testCOde
 USE MyModules1, ONLY : var11, var12, var13
 USE MyModules2, ONLY : var11, var12, var13
 USE mymodules3, ONLY : fjdkl, fjdkla, gja,fjdla,fdjalk, fjklda
 IMPLICIT NONE
 PRIVATE
 TYPE my_type
  INTEGER(i4) :: a, b
  REAL(r8) :: fjdk, fjdkla
 END TYPE
 PUBLIC :: ghost_exchangeV_new
 INTERFACE my_interface
  MODULE PROCEDURE my_sub1
  MODULE PROCEDURE my_sub2
  MODULE PROCEDURE my_sub3
 END INTERFACE myinterface
CONTAINS





 SUBROUTINE my_sub()
   a = 2
 END SUBROUTINE





 SUBROUTINE my_sub0()
 END SUBROUTINE





 SUBROUTINE my_sub1(a, b, c)
  USE mymodule3, ONLY : j, fj, jfkd
  IMPLICIT NONE
  TYPE my_sub_t
   REAL :: abc
  END TYPE my_sub_t
  TYPE(my_sub_t) :: mytt
  REAL(r8), DIMENSION(:, :), ALLOCATABLE :: a
  INTEGER(i4) :: b, c
    ii = 10
    DO i = 1, 100
      jfdk = fjdk
    END DO
    DO j = 1, 100
      DO i = 1, 100
        b = abd
      END DO
    END DO
    DO k = 1, 100
      DO j = 1, 100
        DO i = 1, 100
          b = abd
        END DO
        END DO
    END DO
    ii = 1
    IF(ii .EQ. 1) some = some1
    IF(ii .NE. 1) THEN
      a = b
      DO i = 1, n
        a = b
        DO j = 1, m
          f = 1
        ENDDO
        a = b
      ENDDO
      a = b
    ENDIF
    IF(a == 1) THEN
      b = 2
    ELSE IF(a == 2) THEN
      b = 3
    ELSE
      b = 0
    ENDIF
    RETURN
 END SUBROUTINE my_sub1





 FUNCTION my_fun(a, b) RESULT(d)
  USE mymod1, ONLY : dj, fj
  USE mymod2
  IMPLICIT NONE
  INCLUDE 'jfkdla.inc'
  INTEGER(i4), DIMENSION(:), INTENT(INOUT) :: a, b
!local
  INTEGER(i4) :: i, n
    n = SIZE(a)
    d = 0
    DO i = 1, n
      d = a(i) + b(i)
    END DO
 END FUNCTION





END MODULE
