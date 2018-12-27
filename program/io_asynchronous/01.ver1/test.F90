PROGRAM test
 USE Infra,       ONLY: i4, r4, r8
 USE Infra,       ONLY: Communicator_t
 USE Infra,       ONLY: nprocs, rank, comm
 USE Infra,       ONLY: nios, useAsynchronous
 USE Infra,       ONLY: IniInfra, FinInfra, PrintStatus
 USE Computation, ONLY: gnn, lnn, ista, iend, var
 USE Computation, ONLY: IniCom, FinCom, RunCom
 USE IO
 IMPLICIT NONE
 INCLUDE 'mpif.h'

 ! local 
 INTEGER(i4) :: i, ierr

  TYPE(Communicator_t) :: glb, com, io, grp

  CALL IniInfra(glb, com, io, grp)

  CALL PrintStatus(glb,      'Global Communicator', glb)
  CALL PrintStatus(com, 'Computation Communicator', glb)
  CALL PrintStatus(io,          'I/O Communicator', glb)
  CALL PrintStatus(grp,       'Group Communicator', glb)
 
  CALL MPI_Barrier(glb%comm, ierr)


  IF (com%isuse) THEN

    CALL IniCom(com)

    CALL RunCom(com)

    CALL FinCom(com)

  END IF

  IF (io%isuse) THEN


  END IF


 
  CALL FinInfra()

END PROGRAM test
