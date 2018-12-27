#include <KIM.h>
!
!-------------------------------------------------------------------------------
   module params
!
!
   integer, public, parameter :: internal_edge = 0
   integer, public, parameter :: external_edge = 1
!
! Type of partitioning methods
!
   integer, public, parameter :: recursive = 0, &
                                                                     kway = 1, &
                                                                   volume = 2, &
                                                               wrecursive = 3, &
                                                                    sfcurve = 4
!
!
   end module params
!-------------------------------------------------------------------------------
