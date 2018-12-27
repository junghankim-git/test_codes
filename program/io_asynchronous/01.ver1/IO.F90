MODULE IO
 USE Infra, ONLY: i4, l4, r4, r8
 USE Infra, ONLY: Communicator_t
 IMPLICIT NONE
 INCLUDE 'mpif.h'
 PRIVATE


 CONTAINS




 SUBROUTINE IniIO(com, io, grp)
  IMPLICIT NONE
  TYPE(Communicator_t), INTENT(IN) :: com, io, grp
  ! local


 END SUBROUTINE IniIO




 SUBROUTINE WriteIO(com, io, grp)
  IMPLICIT NONE
  TYPE(Communicator_t), INTENT(IN) :: com, io, grp
  ! local


 END SUBROUTINE WriteIO


END MODULE IO
