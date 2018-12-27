MODULE Infra
 IMPLICIT NONE
 INCLUDE 'mpif.h'
 PRIVATE

 INTEGER, PARAMETER :: i4 = 4
 INTEGER, PARAMETER :: l4 = 4
 INTEGER, PARAMETER :: r4 = 4
 INTEGER, PARAMETER :: r8 = 8

 INTEGER(i4) :: comm, nprocs, rank

 INTEGER(i4) :: nios = 2
 LOGICAL(l4) :: useIODecomposition = .FALSE.
 LOGICAL(l4) :: useAsynchronous    = .FALSE.

 PUBLIC :: i4, l4, r4, r8
 PUBLIC :: nprocs, rank, comm
 PUBLIC :: nios, useIODecomposition, useAsynchronous
 PUBLIC :: IniInfra, FinInfra


 CONTAINS




 SUBROUTINE IniInfra()
  IMPLICIT NONE
  ! local
  INTEGER(i4) :: ierr

   CALL MPI_Init(ierr)
   comm = MPI_COMM_WORLD
   CALL MPI_COMM_SIZE(comm, nprocs, ierr)
   CALL MPI_COMM_RANK(comm,   rank, ierr)

   CALL ReadNL()


 END SUBROUTINE IniInfra
 
 


 SUBROUTINE ReadNL()
  IMPLICIT NONE
  ! local
  INTEGER(i4) :: ierr
  NAMELIST /inputs/ nios, useAsynchronous, useIODecomposition

  IF(rank == 0) THEN
    OPEN(21, FILE='inputs.nl', STATUS='old')
    READ(21, NML=inputs)
    CLOSE(21)
  END IF

  CALL MPI_BCAST(nios,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(useAsynchronous,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(useIODecomposition, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

 END SUBROUTINE ReadNL




 SUBROUTINE FinInfra()
  IMPLICIT NONE
  ! local
  INTEGER(i4) :: ierr
 
   CALL MPI_Finalize(ierr)

 END SUBROUTINE FinInfra





END MODULE Infra
