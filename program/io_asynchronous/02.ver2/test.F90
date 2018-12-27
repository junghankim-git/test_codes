PROGRAM test
 USE Infra,       ONLY: i4, r4, r8
 USE Infra,       ONLY: nprocs, rank, comm
 USE Infra,       ONLY: nios, useIODecomposition, useAsynchronous
 USE Infra,       ONLY: IniInfra, FinInfra
 USE KIO,         ONLY: KIO_t
 USE KIO,         ONLY: IniIO, FinIO, PrintCommunicator
 USE KIO,         ONLY: CreateFile, CloseFile, WriteVariable
 USE Model,       ONLY: RunModel
 IMPLICIT NONE
 INCLUDE 'mpif.h'

 ! local 
 INTEGER(i4) :: i, ierr

  TYPE(KIO_t) :: io_t

  CALL IniInfra()
  CALL IniIO(io_t, comm, nprocs, rank, nios, useIODecomposition, useAsynchronous)

  CALL PrintCommunicator(io_t%glb,      'Global Communicator', io_t%glb)
  CALL PrintCommunicator(io_t%com, 'Computation Communicator', io_t%glb)
  CALL PrintCommunicator(io_t%io,          'I/O Communicator', io_t%glb)
  CALL PrintCommunicator(io_t%grp,       'Group Communicator', io_t%glb)
 
  CALL MPI_Barrier(io_t%glb%comm, ierr)




  IF (io_t%com%isuse) THEN
    ! Model run

    CALL RunModel(io_t%com)

  END IF

  IF (io_t%io%isuse) THEN
    !  I/O  run

!    CALL CreateFile(io_t, 'data.txt')
!    CALL WriteVariable(io_t, lnn, var)
!    CALL CloseFile(io_t)

  END IF


 
  CALL FinIO(io_t)
  CALL FinInfra()

END PROGRAM test
