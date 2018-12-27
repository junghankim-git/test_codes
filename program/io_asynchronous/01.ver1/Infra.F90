MODULE Infra
 IMPLICIT NONE
 INCLUDE 'mpif.h'
 PRIVATE

 INTEGER, PARAMETER :: i4 = 4
 INTEGER, PARAMETER :: l4 = 4
 INTEGER, PARAMETER :: r4 = 4
 INTEGER, PARAMETER :: r8 = 8

 INTEGER(i4) :: nios = 2
 LOGICAL(l4) :: useIODecomposition = .FALSE.
 LOGICAL(l4) :: useAsynchronous    = .FALSE.
 INTEGER(i4) :: nprocs, rank, comm

 TYPE Communicator_t
   INTEGER(i4) :: comm, nprocs, rank, root
   LOGICAL(l4) :: ismaster, isuse
 END TYPE Communicator_t

 INTEGER(i4) :: ierr

 PUBLIC :: Communicator_t
 PUBLIC :: i4, l4, r4, r8
 PUBLIC :: nios, useAsynchronous
 PUBLIC :: nprocs, rank, comm
 PUBLIC :: IniInfra, FinInfra, PrintStatus, Decompose1D


 CONTAINS




 SUBROUTINE IniInfra(glb, com, io, grp)
  IMPLICIT NONE
  TYPE(Communicator_t), INTENT(OUT) :: glb, com, io, grp
  ! local
  INTEGER(i4) :: wrld_grp, com_grp, io_grp, grp_grp
  INTEGER(i4) :: ncoms, ngrps
  INTEGER(i4), DIMENSION(:), ALLOCATABLE :: comranks, ioranks, grpranks

   CALL MPI_Init(ierr)
   comm = MPI_COMM_WORLD
   CALL MPI_COMM_SIZE(comm, nprocs, ierr)
   CALL MPI_COMM_RANK(comm,   rank, ierr)

   CALL ReadNL()

   glb%comm   = comm
   glb%nprocs = nprocs
   glb%rank   = rank
   glb%root   = 0
   glb%isuse  = .TRUE.
   glb%ismaster = .FALSE.
   IF (glb%rank == glb%root) glb%ismaster = .TRUE.

   com%root   = 0
   io%root    = 0
   grp%root   = 0


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
     IF (glb%ismaster) PRINT *, 'useAsynchronous is .TRUE., then set useIODecomposition is .TRUE.'
   END IF

   IF (.NOT. useIODecomposition .AND. nios .LE. nprocs) THEN
     nios = nprocs
     IF (glb%ismaster) PRINT *, 'useIODecomposition is .FALSE., then set nios = nprocs', nios
   END IF

   IF (useIODecomposition .AND. nios .GT. nprocs/2) THEN
     nios = nprocs/2
     IF (glb%ismaster) PRINT *, 'useIODecomposition is .TRUE., then set nios = nprocs/2', nios
     IF (nios .EQ. 0) THEN
       nios = 1
       useAsynchronous  = .FALSE.
       IF (glb%ismaster) PRINT *, 'nios = 0, then set nios = 1, useAsynchronous = .FALSE.', nios
     END IF
   END IF

   IF (nprocs .LT. nios) THEN
     nios = nprocs
     useAsynchronous    = .FALSE.
     useIODecomposition = .FALSE.
     IF (glb%ismaster) PRINT *, 'nprocs < nios, then set nios = nprocs, useAsynchronous = .FALSE., useIODecomposition = .FALSE.'
   END IF


   IF (useAsynchronous) THEN
     ncoms = nprocs - nios
   ELSE
     ncoms = nprocs
   END IF

!old
#if 0
   IF (useAsynchronous) THEN
     IF (IsIoRank(rank)) THEN
       com%isuse = .FALSE.
       io%isuse  = .TRUE.
     ELSE
       com%isuse = .TRUE.
       io%isuse  = .FALSE.
     END IF
   ELSE
     com%isuse = .TRUE.
     IF (IsIoRank(rank)) THEN
       io%isuse  = .TRUE.
     ELSE
       io%isuse  = .FALSE.
     END IF
   END IF
#endif

   ALLOCATE(comranks(ncoms))
   ALLOCATE( ioranks( nios))
   CALL GetComIoGrpRanks(ncoms, nios, comranks, ioranks, ngrps, grpranks)
   com%isuse = .FALSE.
   IF (MINVAL(ABS(comranks-rank)) .EQ. 0) THEN
     com%isuse = .TRUE.
   END IF
   io%isuse = .FALSE.
   IF (MINVAL(ABS(ioranks-rank)) .EQ. 0) THEN
     io%isuse = .TRUE.
   END IF
   grp%isuse = .FALSE.
   IF (MINVAL(ABS(grpranks-rank)) .EQ. 0) THEN
     grp%isuse = .TRUE.
   END IF

   IF (glb%ismaster) THEN
     PRINT *, '* Final Setting '
     PRINT *, ' - useAsynchronous    = ', useAsynchronous
     PRINT *, ' - useIODecomposition = ', useIODecomposition
     PRINT *, ' - ncoms              = ', ncoms
     PRINT *, ' - nios               = ', nios
   END IF
     PRINT *, ' - ngrps              = ', ngrps



   CALL MPI_COMM_GROUP(comm, wrld_grp, ierr)
   CALL MPI_GROUP_INCL(wrld_grp, ncoms, comranks, com_grp, ierr)
   CALL MPI_GROUP_INCL(wrld_grp,  nios,  ioranks,  io_grp, ierr)
   CALL MPI_GROUP_INCL(wrld_grp, ngrps, grpranks, grp_grp, ierr)
   com%comm = -1
    io%comm = -1
   grp%comm = -1
   CALL MPI_COMM_CREATE(comm, com_grp, com%comm, ierr)
   CALL MPI_COMM_CREATE(comm,  io_grp,  io%comm, ierr)
   CALL MPI_COMM_CREATE(comm, grp_grp, grp%comm, ierr)
   com%nprocs = -1
    io%nprocs = -1
   grp%nprocs = -1
   com%rank   = -1
    io%rank   = -1
   grp%rank   = -1
   IF (com%isuse) THEN
     CALL MPI_COMM_SIZE(com%comm, com%nprocs, ierr)
     CALL MPI_COMM_RANK(com%comm,   com%rank, ierr)
   END IF
   IF (io%isuse) THEN
     CALL MPI_COMM_SIZE( io%comm,  io%nprocs, ierr)
     CALL MPI_COMM_RANK( io%comm,    io%rank, ierr)
   END IF
   IF (grp%isuse) THEN
     CALL MPI_COMM_SIZE(grp%comm, grp%nprocs, ierr)
     CALL MPI_COMM_RANK(grp%comm,   grp%rank, ierr)
   END IF

   com%ismaster = .FALSE.
   io%ismaster  = .FALSE.
   grp%ismaster  = .FALSE.
   IF (com%rank == com%root) com%ismaster = .TRUE.
   IF (io%rank  ==  io%root)  io%ismaster = .TRUE.
   IF (grp%rank == grp%root) grp%ismaster = .TRUE.

   DEALLOCATE(comranks)
   DEALLOCATE( ioranks)
   DEALLOCATE(grpranks)

 END SUBROUTINE IniInfra
 
 


 SUBROUTINE ReadNL()
  IMPLICIT NONE
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
 
   CALL MPI_Finalize(ierr)

 END SUBROUTINE FinInfra




! Check here
 FUNCTION IsIoRank(irank) RESULT(isIO)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN)  :: irank
  LOGICAL(l4) :: isIO
  !local
  INTEGER(i4) :: i
  INTEGER(i4) :: base, stride

   isIO = .FALSE.
   stride = nprocs/nios
   IF (stride .GE. 2) THEN
     base   = stride/2
   ELSE
     base   = 0
   END IF

   IF (nprocs .GT. 1) THEN

     DO i = 0, nios-1
       IF (irank == base+i*stride) THEN
         isIO = .TRUE.
         EXIT
       END IF
     END DO

   ELSE

     isIO = .TRUE.

   END IF


 END FUNCTION IsIoRank






 SUBROUTINE GetComIoGrpRanks(ncoms, nios, comranks, ioranks, ngrps, grpranks)
  IMPLICIT NONE
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

     nranks = Decompose1D(0, nprocs-1, nios, i, ista, iend)

     ! Set IO ranks
     IF (nranks == 1) THEN
       ioranks(i+1) = ista
     ELSE
       ioranks(i+1) = ista + (iend-ista+1)/nranks
     END IF

     ! Set Grp ranks
     IF (rank .GE. ista .AND. rank .LE. iend) THEN
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
   DO irank = 0, nprocs - 1

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




 SUBROUTINE GetComIoGrpRanks_old(ncoms, comranks, nios, ioranks)
  IMPLICIT NONE
  INTEGER(i4)              :: ncoms, nios
  INTEGER(i4), INTENT(OUT) :: comranks(ncoms), ioranks(nios)
  ! local
  INTEGER(i4) :: irank
  INTEGER(i4) :: icom, iio

   icom = 0
   iio  = 0


   DO irank = 0, nprocs - 1

     IF (IsIoRank(irank)) THEN
       iio  = iio  + 1
       ioranks(iio)   = irank
     END IF

     IF (useAsynchronous) THEN

       IF (.NOT. IsIoRank(irank)) THEN
         icom = icom + 1
         comranks(icom)   = irank
       END IF

     ELSE

       icom = icom + 1
       comranks(icom) = irank

     END IF

   END DO

 END SUBROUTINE GetComIoGrpRanks_old




 SUBROUTINE PrintStatus(cstr, string, glb)
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

 END SUBROUTINE PrintStatus



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



END MODULE Infra
