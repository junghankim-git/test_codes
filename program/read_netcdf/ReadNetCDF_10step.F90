PROGRAM Test
 USE Kinds, ONLY: i4, r4, r8
 USE NetCDF_Mod, ONLY: OpenFile, CloseFile, IniFile, FinFile, GetDimSize, ReadVariable
 IMPLICIT NONE

 INTEGER(i4), PARAMETER :: nfiles = 11
 CHARACTER(LEN=64)      :: basepath
 CHARACTER(LEN=2)       :: dyn
 CHARACTER(LEN=64)      :: tmp
 CHARACTER(LEN=128)     :: filename

 INTEGER(i4) :: nhoriz, nlevs
 REAL(r8), DIMENSION(:,:), ALLOCATABLE  :: var_dyn, var_phy
 INTEGER(i4) :: ierr, ifile, i, ncount
 LOGICAL :: isFirst


 WRITE(dyn, '(A2)') 'SW'


 isFirst = .TRUE.
 WRITE(basepath, '(A)') '/scratch/jhkim/TestBed/Data/advec/10h/'
 DO ifile = 1, nfiles
     print '(I2)', i+2
     WRITE(filename,'(A,A2,A,I2.2,A,I2.2,A)') TRIM(basepath), dyn, '_G/101/UP-2011072512', ifile-1, '00-0000', ifile, '.nc'
    ! 72 character
    !CALL OpenFile("/scratch/jhkim/TestBed/Data/0.25.06/10h/SH_G/101/UP-20110725130000-000002.nc", ierr)
    CALL OpenFile(filename, ierr)

    CALL IniFile(ierr)
    IF (isFirst) THEN
       nhoriz = GetDimSize("ncol", ierr)
       nlevs   = GetDimSize("lev", ierr)
       ALLOCATE(var_dyn(nhoriz, nlevs))
       ALLOCATE(var_phy(nhoriz, nlevs))
       isFirst = .FALSE.
    END IF
   
   
    CALL ReadVariable("qr_dyn", var_dyn, ierr)
    CALL ReadVariable("qr", var_phy, ierr)
   
    !PRINT *, var(:,1)
   
!    CALL CheckSmall(nhoriz, nlevs, var_dyn, ncount, .TRUE., .FALSE.) 
PRINT *, 'dynamics'
    CALL CheckSmall(nhoriz, nlevs, var_dyn, ncount, .TRUE., .FALSE.) 
    CALL CheckSmall(nhoriz, nlevs, var_dyn, ncount, .FALSE., .TRUE.) 
PRINT *, 'physics'
    CALL CheckSmall(nhoriz, nlevs, var_phy, ncount, .TRUE., .FALSE.) 
    CALL CheckSmall(nhoriz, nlevs, var_phy, ncount, .FALSE., .TRUE.) 
PRINT *, 'd = epsilon, p=0'
    CALL CheckTwo(nhoriz, nlevs, var_dyn, var_phy, ncount)
    CALL CheckTwo(nhoriz, nlevs, var_phy, var_dyn, ncount)
!    CALL CheckSmall(nhoriz, nlevs, var, ncount, .TRUE.,  .TRUE.) 
   
    CALL FinFile(ierr)
    CALL CloseFile(ierr)

 END DO

 DEALLOCATE(var_dyn)
 DEALLOCATE(var_phy)

 CONTAINS


 SUBROUTINE CheckSmall(nx, ny, var_in, ncount_in, isMinus, isSmall)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN)                :: nx, ny
  REAL(r8), DIMENSION(nx,ny), INTENT(IN) :: var_in
  INTEGER(i4), INTENT(OUT)               :: ncount_in
  LOGICAL, INTENT(IN)                    :: isMinus, isSmall

  ! local
  INTEGER(i4) :: i, j
  REAL(r8)    :: small, zero

  ncount_in = 0
  zero      = 0.0D0
  small     = 1.0D-150

  IF (isMinus .AND. .NOT. isSmall) THEN

    DO j = 1, ny
      DO i = 1, nx
        IF (var_in(i,j) .LT. zero) THEN
          ncount_in = ncount_in + 1
        END IF
      END DO
    END DO

  ELSE IF (.NOT. isMinus .AND. isSmall) THEN

    DO j = 1, ny
      DO i = 1, nx
        IF ((var_in(i,j) .LT. small) .AND. (var_in(i,j) .GT. -small) .AND. var_in(i,j) .NE. 0.0D0) THEN
          ncount_in = ncount_in + 1
 !         print *, var_in(i,j)
        END IF
      END DO
    END DO

  ELSE IF (isMinus .AND. isSmall) THEN

    DO j = 1, ny
      DO i = 1, nx
        IF ((var_in(i,j) .LT. zero) .AND. ((var_in(i,j) .LT. small) .AND. (var_in(i,j) .GT. -small) .AND. var_in(i,j) .NE. 0.0D0)) THEN
          ncount_in = ncount_in + 1
 !         print *, var_in(i,j)
        END IF
      END DO
    END DO

  ELSE

    WRITE(*,*) 'Check check small option...'
    STOP

  END IF

  IF (isMinus .AND. .NOT. isSmall) THEN
    WRITE(*,'(A,I8,X,I8,X,F10.7,X,A)') 'minus check: ', ncount_in, nx*ny, DBLE(ncount_in)/DBLE(nx*ny)*100.0D0, '%'
  ELSE IF (.NOT. isMinus .AND. isSmall) THEN
    WRITE(*,'(A,I8,X,I8,X,F10.7,X,A)') 'small check: ', ncount_in, nx*ny, DBLE(ncount_in)/DBLE(nx*ny)*100.0D0, '%'
!    WRITE(*,*) 'small value = ', small, TINY(small), HUGE(small)
  ELSE IF (isMinus .AND. isSmall) THEN
    WRITE(*,'(A,I8,X,I8,X,F10.7,X,A)') 'two   check: ', ncount_in, nx*ny, DBLE(ncount_in)/DBLE(nx*ny)*100.0D0, '%'
!    WRITE(*,*) 'small value = ', small, TINY(small)
  END IF

 END SUBROUTINE CheckSmall


 SUBROUTINE CheckTwo(nx, ny, var1_in, var2_in, ncount_in)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN)                :: nx, ny
  REAL(r8), DIMENSION(nx,ny), INTENT(IN) :: var1_in, var2_in
  INTEGER(i4), INTENT(OUT)               :: ncount_in

  ! local
  INTEGER(i4) :: i, j
  REAL(r8)    :: small, zero

  ncount_in = 0
  zero      = 0.0D0
  small     = 1.0D-150


    DO j = 1, ny
      DO i = 1, nx
        IF (((var1_in(i,j) .LT. small) .AND. (var1_in(i,j) .GT. -small) .AND. (var1_in(i,j) .NE. 0.0D0)) .AND. (var2_in(i,j) .EQ. 0.0D0)) THEN
          ncount_in = ncount_in + 1
 !         print *, var_in(i,j)
        END IF
      END DO
    END DO


    WRITE(*,'(A,I8,X,I8,X,F10.7,X,A)') ' two check: ', ncount_in, nx*ny, DBLE(ncount_in)/DBLE(nx*ny)*100.0D0, '%'

 END SUBROUTINE CheckTwo


END PROGRAM Test
