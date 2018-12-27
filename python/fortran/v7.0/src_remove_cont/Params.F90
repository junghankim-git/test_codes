#include <KIM.h>
!
MODULE Params
!
  INTEGER, PUBLIC, PARAMETER :: INTERNAL_EDGE  = 0
  INTEGER, PUBLIC, PARAMETER :: EXTERNAL_EDGE  = 1
!
! Type of partitioning methods
!
  INTEGER, PUBLIC, PARAMETER :: RECURSIVE  = 0,  &
                                KWAY       = 1,  &
                                VOLUME     = 2,  &
                                WRECURSIVE = 3,  &
                                SFCURVE    = 4
!
END MODULE Params
