MODULE Model
 USE Infra, ONLY: i4, l4, r4, r8
 USE KIO,   ONLY: Communicator_t
 IMPLICIT NONE
 INCLUDE 'mpif.h'
 PRIVATE


  PUBLIC :: RunModel

 CONTAINS



 SUBROUTINE RunModel(com)
  USE Computation, ONLY: gnn, lnn, ista, iend, var
  USE Computation, ONLY: IniCom, FinCom, RunCom
  IMPLICIT NONE
  TYPE(Communicator_t), INTENT(IN) :: com
  ! local

    CALL IniCom(com)

    CALL RunCom(com)
!    CALL CreateFile(io_t, 'data.txt')
!    CALL WriteVariable(io_t, lnn, var)
!    CALL CloseFile(io_t)

    CALL FinCom(com)

 END SUBROUTINE RunModel


END MODULE Model
