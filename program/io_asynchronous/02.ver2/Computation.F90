MODULE Computation
 USE Infra, ONLY: i4, l4, r4, r8
 USE KIO,   ONLY: Communicator_t
 USE KIO,   ONLY: Decompose1D
 IMPLICIT NONE
 INCLUDE 'mpif.h'
 PRIVATE


  INTEGER(i4), PARAMETER    :: gnn    = 100
  INTEGER(i4), PARAMETER    :: nsteps = 5
  INTEGER(i4), DIMENSION(:), ALLOCATABLE :: var

  INTEGER(i4)               :: lnn
  INTEGER(i4)               :: ista, iend


  PUBLIC :: gnn, lnn, ista, iend
  PUBLIC :: var
  PUBLIC :: IniCom, RunCom, FinCom

 CONTAINS



 SUBROUTINE IniCom(com)
  IMPLICIT NONE
  TYPE(Communicator_t), INTENT(IN) :: com
  ! local
  INTEGER(i4) :: i

   lnn = Decompose1D(1, gnn, com%nprocs, com%rank, ista, iend)

   ALLOCATE(var(ista:iend))

   DO i = ista, iend
     var(i) = 0
   END DO

 END SUBROUTINE IniCom



 SUBROUTINE RunCom(com)
  IMPLICIT NONE
  TYPE(Communicator_t), INTENT(IN) :: com
  ! local
  INTEGER(i4) :: istep, i, ierr
  INTEGER(i4) :: gsum

   DO istep = 1, nsteps

     var  = GenVar(ista, iend, istep)

     gsum = GlobalSum(com, ista, iend, var)
     IF (com%ismaster) THEN
       PRINT '(A,I3,A,I8)', 'istep = ', istep, ',  global sum = ', gsum
     END IF

   END DO

 END SUBROUTINE RunCom



 FUNCTION GenVar(is, ie, step) RESULT(buf)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN) :: is, ie, step
  INTEGER(i4)             :: buf(is:ie)
  ! local
  INTEGER(i4) :: i, j

   DO i = ista, iend
     buf(i) = 100*step + i
   END DO

 END FUNCTION GenVar



 FUNCTION GlobalSum(com, is, ie, buf) RESULT(glbsum)
  IMPLICIT NONE
  TYPE(Communicator_t), INTENT(IN) :: com
  INTEGER(i4), INTENT(IN) :: is, ie
  INTEGER(i4), INTENT(IN) :: buf(is:ie)
  INTEGER(i4)             :: glbsum
  ! local
  INTEGER(i4) :: i, ierr, lsum

   lsum = 0
   DO i = ista, iend
     lsum = lsum + var(i)
   END DO
   CALL MPI_REDUCE(lsum, glbsum, 1, MPI_INTEGER, MPI_SUM, com%root, com%comm, ierr)

 END FUNCTION GlobalSum


 SUBROUTINE FinCom(com)
  IMPLICIT NONE
  TYPE(Communicator_t), INTENT(IN) :: com
  ! local

   DEALLOCATE(var)

 END SUBROUTINE FinCom



END MODULE Computation
