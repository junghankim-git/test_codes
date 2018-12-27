MODULE KIO
 IMPLICIT NONE
 INCLUDE 'mpif.h'
 PRIVATE

 INTEGER, PARAMETER :: i4 = 4
 INTEGER, PARAMETER :: l4 = 4
 INTEGER, PARAMETER :: r4 = 4
 INTEGER, PARAMETER :: r8 = 8

 TYPE Communicator_t
   INTEGER(i4) :: comm, nprocs, rank, root
   LOGICAL(l4) :: ismaster, isuse
 END TYPE Communicator_t


 TYPE KIO_t
   LOGICAL(l4) :: useIODecomposition
   LOGICAL(l4) :: useAsynchronous

   INTEGER(i4) :: ncoms
   INTEGER(i4) :: nios
   INTEGER(i4) :: ngrps

   TYPE(Communicator_t) :: glb
   TYPE(Communicator_t) :: com
   TYPE(Communicator_t) :: io
   TYPE(Communicator_t) :: grp
 END TYPE KIO_t


 LOGICAL(l4) :: isopen = .FALSE. 

 PUBLIC :: Communicator_t, KIO_t
 PUBLIC :: PrintCommunicator, Decompose1D
 PUBLIC :: IniIO, FinIO
 PUBLIC :: CreateFile, CloseFile, WriteVariable

 CONTAINS



 SUBROUTINE IniIO(io_t, comm, nprocs, rank, nios, useIODecomposition, useAsynchronous)
  IMPLICIT NONE
  TYPE(KIO_t), INTENT(INOUT) :: io_t
  INTEGER(i4), INTENT(IN) :: comm, nprocs, rank
  INTEGER(i4), INTENT(INOUT) :: nios
  LOGICAL(l4), INTENT(INOUT) :: useIODecomposition, useAsynchronous
  ! local

  INTEGER(i4) :: wrld_grp, com_grp, io_grp, grp_grp
  INTEGER(i4) :: ncoms, ngrps
  INTEGER(i4), DIMENSION(:), ALLOCATABLE :: comranks, ioranks, grpranks
  INTEGER(i4) :: ierr

   io_t%glb%comm     = comm
   io_t%glb%nprocs   = nprocs
   io_t%glb%rank     = rank
   io_t%glb%root     = 0
   io_t%glb%isuse    = .TRUE.
   io_t%glb%ismaster = .FALSE.
   IF (io_t%glb%rank == io_t%glb%root) io_t%glb%ismaster = .TRUE.

   io_t%com%root   = 0
   io_t%io%root    = 0
   io_t%grp%root   = 0


   !-----------------------------------
   ! Asynchronous: .TRUE.
   !-----------------------------------
   ! case 1 (Decomposition: .TRUE.)
   ! * * * * * * * * * * * *    : world
   ! @ @ @ @ @ @ ! 1/2(max)     : I/O
   !-----------------------------------

   !-----------------------------------
   ! Asynchronous: .FALSE.
   !-----------------------------------
   ! case 1 (Decomposition: .FALSE.)
   ! * * * * * * * * * * * *    : World
   ! @ @ @ @ @ @ @ @ @ @ @ @    :  I/O
   ! case 2 (Decomposition: .TRUE.)
   ! * * * * * * * * * * * *    : World
   ! @ @ @ @ @ @ ! 1/2(max)     :  I/O
   !-----------------------------------

   IF (useAsynchronous .AND. .NOT. useIODecomposition) THEN
     useIODecomposition = .TRUE.
     IF (io_t%glb%ismaster) PRINT *, 'useAsynchronous is .TRUE., then set useIODecomposition is .TRUE.'
   END IF

   IF (.NOT. useIODecomposition .AND. nios .LE. nprocs) THEN
     nios = nprocs
     IF (io_t%glb%ismaster) PRINT *, 'useIODecomposition is .FALSE., then set nios = nprocs', nios
   END IF

   IF (useIODecomposition .AND. nios .GT. nprocs/2) THEN
     nios = nprocs/2
     IF (io_t%glb%ismaster) PRINT *, 'useIODecomposition is .TRUE., then set nios = nprocs/2', nios
     IF (nios .EQ. 0) THEN
       nios = 1
       useAsynchronous  = .FALSE.
       IF (io_t%glb%ismaster) PRINT *, 'nios = 0, then set nios = 1, useAsynchronous = .FALSE.', nios
     END IF
   END IF


   IF (nprocs .LT. nios) THEN
     nios = nprocs
     useAsynchronous    = .FALSE.
     useIODecomposition = .FALSE.
     IF (io_t%glb%ismaster) PRINT *, 'nprocs < nios, then set nios = nprocs, useAsynchronous = .FALSE., useIODecomposition = .FALSE.'
   END IF


   IF (useAsynchronous) THEN
     ncoms = nprocs - nios
   ELSE
     ncoms = nprocs
   END IF



   ALLOCATE(comranks(ncoms))
   ALLOCATE( ioranks( nios))
   CALL GetComIoGrpRanks(io_t%glb, useAsynchronous, ncoms, nios, comranks, ioranks, ngrps, grpranks)
   io_t%com%isuse = .FALSE.
   IF (MINVAL(ABS(comranks-rank)) .EQ. 0) THEN
     io_t%com%isuse = .TRUE.
   END IF
   io_t%io%isuse = .FALSE.
   IF (MINVAL(ABS(ioranks-rank)) .EQ. 0) THEN
     io_t%io%isuse = .TRUE.
   END IF
   io_t%grp%isuse = .FALSE.
   IF (MINVAL(ABS(grpranks-rank)) .EQ. 0) THEN
     io_t%grp%isuse = .TRUE.
   END IF

   IF (io_t%glb%ismaster) THEN
     PRINT *, '* Final Setting '
     PRINT *, ' - useAsynchronous    = ', useAsynchronous
     PRINT *, ' - useIODecomposition = ', useIODecomposition
     PRINT *, ' - ncoms              = ', ncoms
     PRINT *, ' - nios               = ', nios
   END IF
     PRINT *, ' - ngrps              = ', ngrps
   io_t%useAsynchronous = useAsynchronous
   io_t%useIODecomposition = useIODecomposition
   io_t%ncoms = ncoms
   io_t%nios  = nios
   io_t%ngrps = ngrps



   CALL MPI_COMM_GROUP(comm, wrld_grp, ierr)
   CALL MPI_GROUP_INCL(wrld_grp, ncoms, comranks, com_grp, ierr)
   CALL MPI_GROUP_INCL(wrld_grp,  nios,  ioranks,  io_grp, ierr)
   CALL MPI_GROUP_INCL(wrld_grp, ngrps, grpranks, grp_grp, ierr)
   io_t%com%comm = -1
    io_t%io%comm = -1
   io_t%grp%comm = -1
   CALL MPI_COMM_CREATE(comm, com_grp, io_t%com%comm, ierr)
   CALL MPI_COMM_CREATE(comm,  io_grp,  io_t%io%comm, ierr)
   CALL MPI_COMM_CREATE(comm, grp_grp, io_t%grp%comm, ierr)
   io_t%com%nprocs = -1
    io_t%io%nprocs = -1
   io_t%grp%nprocs = -1
   io_t%com%rank   = -1
    io_t%io%rank   = -1
   io_t%grp%rank   = -1
   IF (io_t%com%isuse) THEN
     CALL MPI_COMM_SIZE(io_t%com%comm, io_t%com%nprocs, ierr)
     CALL MPI_COMM_RANK(io_t%com%comm,   io_t%com%rank, ierr)
   END IF
   IF (io_t%io%isuse) THEN
     CALL MPI_COMM_SIZE( io_t%io%comm,  io_t%io%nprocs, ierr)
     CALL MPI_COMM_RANK( io_t%io%comm,    io_t%io%rank, ierr)
   END IF
   IF (io_t%grp%isuse) THEN
     CALL MPI_COMM_SIZE(io_t%grp%comm, io_t%grp%nprocs, ierr)
     CALL MPI_COMM_RANK(io_t%grp%comm,   io_t%grp%rank, ierr)
   END IF

   io_t%com%ismaster = .FALSE.
   io_t%io%ismaster  = .FALSE.
   io_t%grp%ismaster = .FALSE.
   IF (io_t%com%rank == io_t%com%root) io_t%com%ismaster = .TRUE.
   IF (io_t%io%rank  ==  io_t%io%root)  io_t%io%ismaster = .TRUE.
   IF (io_t%grp%rank == io_t%grp%root) io_t%grp%ismaster = .TRUE.

   DEALLOCATE(comranks)
   DEALLOCATE( ioranks)
   DEALLOCATE(grpranks)


 END SUBROUTINE IniIO





 SUBROUTINE FinIO(io_t)
  IMPLICIT NONE
  TYPE(KIO_t), INTENT(IN) :: io_t

 END SUBROUTINE FinIO






 SUBROUTINE GetComIoGrpRanks(glb_t, useAsynchronous, ncoms, nios, comranks, ioranks, ngrps, grpranks)
  IMPLICIT NONE
  TYPE(Communicator_t), INTENT(IN) :: glb_t
  LOGICAL(l4), INTENT(IN)  :: useAsynchronous
  INTEGER(i4), INTENT(IN)  :: ncoms, nios
  INTEGER(i4), INTENT(OUT) :: comranks(ncoms), ioranks(nios)
  INTEGER(i4), INTENT(OUT) :: ngrps
  INTEGER(i4), ALLOCATABLE, INTENT(OUT) :: grpranks(:)
  ! local
  INTEGER(i4) :: irank, i, j
  INTEGER(i4) :: ista, iend, nranks
  INTEGER(i4) :: icom, iio, igrp

   iio  = 0
   DO i = 0, nios-1

     nranks = Decompose1D(0, glb_t%nprocs-1, nios, i, ista, iend)

     ! Set IO ranks
     IF (nranks == 1) THEN
       ioranks(i+1) = ista
     ELSE
       ioranks(i+1) = ista + (iend-ista+1)/nranks
     END IF

     ! Set Grp ranks
     IF (glb_t%rank .GE. ista .AND. glb_t%rank .LE. iend) THEN
       ngrps = nranks
       ALLOCATE(grpranks(ngrps))
!----------------------------------
! org: no sorting
!       DO j = 1, nranks
!         grpranks(j) = ista+j-1
!       END DO
!----------------------------------
!----------------------------------
! new: I/O rank is the master
       grpranks(1) = ioranks(i+1)
       igrp = 1
       DO j = 1, nranks
         irank = ista+j-1
         IF (irank .NE. ioranks(i+1)) THEN
           igrp  = igrp + 1
           grpranks(igrp) = irank
         END IF
       END DO
!----------------------------------
     END IF

   END DO

   ! Set Com ranks
   icom = 0
   DO irank = 0, glb_t%nprocs - 1

     IF (useAsynchronous) THEN

       IF (MINVAL(ABS(ioranks-irank)) .NE. 0) THEN
         icom = icom + 1
         comranks(icom)   = irank
       END IF

     ELSE

       icom = icom + 1
       comranks(icom) = irank

     END IF

   END DO

 END SUBROUTINE GetComIoGrpRanks




 SUBROUTINE PrintCommunicator(cstr, string, glb)
  IMPLICIT NONE
  TYPE(Communicator_t), INTENT(IN) :: cstr, glb
  CHARACTER(LEN=*), INTENT(IN)     :: string
  ! local
  INTEGER(i4) :: i, j, ierr  
  INTEGER(i4) :: icomm, inprocs, irank, iroot, stat(MPI_STATUS_SIZE)
  LOGICAL(l4) :: iuse

   IF (glb%ismaster) THEN
     PRINT '(A)', ' '
     PRINT '(A)', string
     PRINT '(A)', '==========================================================='
     PRINT '(A)', '|  rank  |    comm    | nprocs |  rank  |  root  |  isuse |'
     PRINT '(A)', '-----------------------------------------------------------'
     PRINT '(A,I8,A,I12,A,I8,A,I8,A,I8,A,L8,A)', '|', cstr%root, '|', cstr%comm, '|', cstr%nprocs, '|', cstr%rank, '|', cstr%root, '|', cstr%isuse, '|'
     PRINT '(A)', '-----------------------------------------------------------'
     DO i = 1, glb%nprocs-1
       CALL MPI_RECV(icomm,   1, MPI_INTEGER, i, 111, glb%comm, stat, ierr)
       CALL MPI_RECV(inprocs, 1, MPI_INTEGER, i, 112, glb%comm, stat, ierr)
       CALL MPI_RECV(irank,   1, MPI_INTEGER, i, 113, glb%comm, stat, ierr)
       CALL MPI_RECV(iroot,   1, MPI_INTEGER, i, 114, glb%comm, stat, ierr)
       CALL MPI_RECV(iuse,    1, MPI_LOGICAL, i, 115, glb%comm, stat, ierr)
     PRINT '(A,I8,A,I12,A,I8,A,I8,A,I8,A,L8,A)', '|', i, '|', icomm, '|', inprocs, '|', irank, '|', iroot, '|', iuse, '|'
     PRINT '(A)', '-----------------------------------------------------------'
     END DO
     !PRINT '(A)', ' '
   ELSE
     CALL MPI_SEND(cstr%comm,   1, MPI_INTEGER, glb%root, 111, glb%comm, ierr)
     CALL MPI_SEND(cstr%nprocs, 1, MPI_INTEGER, glb%root, 112, glb%comm, ierr)
     CALL MPI_SEND(cstr%rank,   1, MPI_INTEGER, glb%root, 113, glb%comm, ierr)
     CALL MPI_SEND(cstr%root  , 1, MPI_INTEGER, glb%root, 114, glb%comm, ierr)
     CALL MPI_SEND(cstr%isuse , 1, MPI_LOGICAL, glb%root, 115, glb%comm, ierr)
   END IF

 END SUBROUTINE PrintCommunicator



 FUNCTION Decompose1D(n1, n2, nprocs, rank, ista, iend) RESULT(nd)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN)            :: n1, n2
  INTEGER(i4), INTENT(IN)            :: nprocs, rank
  INTEGER(i4), OPTIONAL, INTENT(OUT) :: ista, iend
  INTEGER(i4)                        :: nd

  ! local
  INTEGER(i4)  :: l_sta, l_end, domain, extra

   domain = n2-n1+1
   nd     = domain/nprocs
   extra  = MOD(domain, nprocs)
   l_sta  = nd*rank+n1 + MIN(rank, extra)
   l_end  = l_sta+nd-1
   IF (rank <= extra-1) THEN
     l_end = l_end + 1
   ENDIF
   nd     = l_end-l_sta+1

 
   IF (PRESENT(ista)) THEN
     ista = l_sta
   END IF
   IF (PRESENT(iend)) THEN
     iend = l_end
   END IF

   RETURN

 END FUNCTION Decompose1D








 SUBROUTINE CreateFile_com(io_t, filename)
  IMPLICIT NONE
  TYPE(KIO_t), INTENT(IN)      :: io_t
  CHARACTER(LEN=*), INTENT(IN) :: filename
  ! local


  IF (io_t%io%ismaster) THEN

    OPEN(21, FILE=filename, STATUS='UNKNOWN', FORM='FORMATTED')

  END IF
  isopen = .TRUE.

 END SUBROUTINE CreateFile_com








 SUBROUTINE CreateFile(io_t, filename)
  IMPLICIT NONE
  TYPE(KIO_t), INTENT(IN)      :: io_t
  CHARACTER(LEN=*), INTENT(IN) :: filename
  ! local


  IF (io_t%io%ismaster) THEN

    OPEN(21, FILE=filename, STATUS='UNKNOWN', FORM='FORMATTED')

  END IF
  isopen = .TRUE.

 END SUBROUTINE CreateFile





 SUBROUTINE CloseFile(io_t)
  IMPLICIT NONE
  TYPE(KIO_t), INTENT(IN)      :: io_t
  ! local


  IF (io_t%io%ismaster) THEN

    CLOSE(21)

  END IF
  isopen = .FALSE.

 END SUBROUTINE CloseFile




 SUBROUTINE WriteVariable(io_t, n, var)
  IMPLICIT NONE
  TYPE(KIO_t), INTENT(IN)      :: io_t
  INTEGER(i4), INTENT(IN) :: n
  INTEGER(i4), INTENT(IN) :: var(n)
  ! local
  INTEGER(i4) :: i

   IF (isopen) THEN
     IF (io_t%io%ismaster) THEN
       DO i = 1, n
         WRITE(21,*) var(i)
       END DO
     END IF
   ELSE
     IF (io_t%io%ismaster) THEN
       PRINT *, 'could not writing, check file create...'
     END IF
   END IF


 END SUBROUTINE WriteVariable


END MODULE KIO
