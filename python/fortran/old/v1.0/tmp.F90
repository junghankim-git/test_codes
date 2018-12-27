!-------------------------------------------------------------------------------
!>
!> @brief
!>  - bndry module.
!>
!> @date ?????2012
!>  - JH KIM  : First written from the HOMME and was modified for KIAPSGM framework.
!> @date 30JAN2015
!>  - JH KIM  : Added to "ghost_exchangeV" subroutine. (CAM-SE 5.3)
!-------------------------------------------------------------------------------
#include "KIM.h"
MODULE testCOde
 USE MyModules1, ONLY : var11, var12, var13
 USE MyModules2, ONLY : var11, var12, var13
 IMPLICIT NONE
 PRIVATE (3)
 TYPE my_type (4)
  INTEGER(i4) :: a, b
  REAL(r8) :: fjdk, fjdkla (3)
 END TYPE (3)

 PUBLIC :: ghost_exchangeV_new
 INTERFACE my_interface
   MODULE PROCEDURE my_func1
   MODULE PROCEDURE my_func2
   MODULE PROCEDURE my_func3
 END INTERFACE myinterface
CONTAINS

 SUBROUTINE my_func1(a, b, c)
  USE mymodule3, ONLY : j, fj, jfkd
  IMPLICIT NONE
  REAL(r8), DIMENSION(:, :), ALLOCATABLE :: a
  INTEGER(i4) :: b, c
    jj = 100 (5)
    DO i = 1, 100 (4)
      jfdk = fjdk (5)
    END DO (5)
    ii = 1
    IF(ii .EQ. 1) some = some1 (5)
    IF(ii .NE. 1) THEN (4)
      DO i = 1, n (4)
        DO j = 1, m (4)
          f = 1 (5)
      ENDDO (4)
    ENDDO (5)

    ENDIF

    RETURN
  END SUBROUTINE my_func1
