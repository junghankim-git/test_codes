PROGRAM TestMain
 USE Kinds,      ONLY: i4, l4, r8
 USE Logger,     ONLY: PrintLog, toString
 USE IoASCII,   ONLY: OpenFile, CloseFile, WriteData_1Log1, WriteData_1Int1, WriteData_1Real, WriteData_String
 USE Quadrature, ONLY: quadrature_t, gausslobatto
 IMPLICIT NONE
 INTEGER(i4)                         :: unitnum
 INTEGER(i4)                         :: nn
 INTEGER(i4)                         :: elev, llev
 REAL(r8), DIMENSION(:), ALLOCATABLE :: x
 REAL(r8), DIMENSION(:), ALLOCATABLE :: cx
 TYPE(quadrature_t)                  :: qp

   elev = 0
   llev = 2

   nn = 4

   CALL ReadNameList('inputs.nl')
   ALLOCATE(x(nn))
   ALLOCATE(cx(nn+1))
   x(:)  = 0.0_r8
   cx(:) = 0.0_r8
   
!   qp4=gausslobatto(4)
!   qp6=gausslobatto(6)
!   print *, qp4%points(:)
!   print *, ((qp4%points(2)-qp4%points(1))-0.5_r8*(qp4%points(3)-qp4%points(2)))
!   print *, qp5%points(:)
!   print *, (0.5_r8*(qp5%points(3)-qp5%points(2)) - (qp5%points(2)-qp5%points(1)))
!   print *, qp6%points(:)
!   print *, (0.5_r8*(qp6%points(3)-qp6%points(2)) - (qp6%points(2)-qp6%points(1)))


   qp=gausslobatto(nn)
   x(:) = qp%points(:)
   CALL GenControlVolume(nn, x, cx)
   CALL CheckGridCell(nn, x, cx)

   CALL OpenFile('Result.dat', unitnum)
   CALL WriteData_1Int1(unitnum, nn)
   CALL WriteData_1Real(unitnum, nn, x)
   CALL WriteData_1Real(unitnum, nn+1, cx)
   CALL CloseFile(unitnum)


   DEALLOCATE(x)
   DEALLOCATE(cx)
 CONTAINS


 SUBROUTINE ReadNameList(nlFileName)
  IMPLICIT NONE
   CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: nlFileName
   NAMELIST /inputs/ nn

   CALL PrintLog('START: Read Namelist...', elev, llev)
   IF(PRESENT(nlFileNAme)) THEN
     OPEN(21, FILE=TRIM(nlFileName), STATUS='old')
     READ(21, NML=inputs)
     CLOSE(21)
   ELSE
     READ(*, NML=inputs)
   END IF

   CALL PrintLog('nn = '//toString(nn), elev, llev+1)
   CALL PrintLog('END  : Read Namelist...', elev, llev)

 END SUBROUTINE ReadNameList



 SUBROUTINE GenControlVolume(n, grid, cell)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN)               :: n
  REAL(r8), DIMENSION(n), INTENT(IN)    :: grid
  REAL(r8), DIMENSION(0:n), INTENT(OUT) :: cell

  !local
  INTEGER(i4)          :: i
  LOGICAL(l4)          :: leven
  REAL(r8)             :: dx0, dxc

   IF (n < 2) THEN
     PRINT *, 'must be greater then 2 ...'
   END IF

   IF ((n/2)*2 == n) THEN
     leven = .TRUE.
   ELSE
     leven = .FALSE.
   END IF

   IF (leven) THEN
 
     dxc = 0.5_r8*(grid(n/2+1)-grid(n/2))

     ! center
     cell(n/2) = grid(n/2)+dxc
 
     DO i = n/2-1, 0, -1
       cell(  i) = grid(i+1) - (cell(i+1)-grid(i+1))
       cell(n-i) = grid(n-i) - (cell(n-i-1)-grid(n-i))
     END DO

   ELSE

     IF (n == 2) THEN
       dx0 = 0.5_r8*(grid(2)-grid(1))
     ELSE
       dx0 = 0.5_r8*(grid(2)-grid(1))**2.0_r8/(grid(3)-grid(2))
     END IF
     cell(  0) = grid( 1)-dx0
     cell(  1) = grid( 1)+dx0
     cell(n-1) = grid( n)-dx0
     cell(  n) = grid( n)+dx0
 
     DO i = 2, n/2
       !2-3
       cell(  i) = grid(  i) + (grid(  i) - cell(  i-1))
       !n-2, n-3
       !n-2      = grid(n-2) + (grid(n-2) - cell(n-3))
       cell(n-i) = grid(n-i+1) - (cell(n-i+1) - grid(n-i+1))
     END DO

   END IF


 END SUBROUTINE GenControlVolume




 SUBROUTINE CheckGridCell(n, grid, cell)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN)              :: n
  REAL(r8), DIMENSION(n), INTENT(IN)   :: grid
  REAL(r8), DIMENSION(0:n), INTENT(IN) :: cell

  !local
  INTEGER(i4)              :: i
  INTEGER(i4), PARAMETER   :: width = 9              ! width of decimal number (include plus, minus, point), (must be less than 10)
  !INTEGER(i4), PARAMETER   :: bel_dec_width = 3      ! width of below decimal number (must be less than or equal to width-3)
  INTEGER(i4), PARAMETER   :: up_dec_width = 4       ! width of upper decimal number (must be greater than or equal to 3)
  CHARACTER(LEN=width)     :: str_g, str_c
  CHARACTER(LEN=width*2+1) :: two_nums
  CHARACTER(LEN=(n+2)*(width*2+1)) :: str_pos, str_val
  ! formats
  CHARACTER(LEN=4)         :: grid_idx_fmt
  CHARACTER(LEN=9)         :: cell_idx_fmt
  CHARACTER(LEN=6)         :: each_val_fmt
  CHARACTER(LEN=7)         :: two_nums_fmt


   ! formats
   WRITE(grid_idx_fmt,'(A,I1,A)')      '(I', width, ')'
   WRITE(cell_idx_fmt,'(A,I1,A)')      '(I', width-2, '.1,A2)'
   WRITE(each_val_fmt,'(A,I1,A,I1,A)') '(F', width, '.', width-up_dec_width, ')'
   WRITE(two_nums_fmt,'(A)')           '(A,X,A)'
#if 0
   PRINT *, 'grid_idx_fmt (eg.         1) = ', grid_idx_fmt
   PRINT *, 'cell_idx_fmt (eg.       3/2) = ', cell_idx_fmt
   PRINT *, 'each_val_fmt (eg.      3.12) = ', each_val_fmt
   PRINT *, 'two_nums_fmt (eg. 3.12 2.22) = ', two_nums_fmt
   PRINT *, 'grid_idx_fmt (eg.         1) = ', grid_idx_fmt
   PRINT *, 'cell_idx_fmt (eg.       3/2) = ', cell_idx_fmt
   PRINT *, 'each_val_fmt (eg.      3.12) = ', each_val_fmt
   PRINT *, 'two_nums_fmt (eg. 3.12 2.22) = ', two_nums_fmt
#endif
   !
 
 
   ! position
   WRITE(str_c,cell_idx_fmt) (0*2+1), '/2'
   WRITE(two_nums, two_nums_fmt) '      ', str_c
   str_pos = two_nums
   DO i = 1, n
     WRITE(str_g,grid_idx_fmt) i
     WRITE(str_c,cell_idx_fmt) (i*2+1), '/2'
     WRITE(two_nums, two_nums_fmt) (str_g), (str_c)
     str_pos = TRIM(str_pos)//' '//two_nums
   END DO
 
   ! value
   WRITE(str_c,each_val_fmt) cell(0)
   WRITE(two_nums, two_nums_fmt) '      ', str_c
   str_val = two_nums
   DO i = 1, n
     WRITE(str_g, each_val_fmt) grid(i)
     WRITE(str_c, each_val_fmt) cell(i)
     WRITE(two_nums, two_nums_fmt) (str_g), (str_c)
     str_val = TRIM(str_val)//' '//two_nums
   END DO
 
   ! check grid, cell
   CALL PrintLog(' ', elev)
   CALL PrintLog('START: Check position of points...', elev, llev)
   CALL PrintLog('Positions of points...', elev, llev+1)
   CALL PrintLog(str_pos, elev, llev+1)
   CALL PrintLog(str_val, elev, llev+1)
   !CALL PrintLog(' ', elev)
 
 
   DO i = 2, n
     IF (grid(i) <= grid(i-1)) THEN
       CALL PrintLog('check grids..., index:'//toString(i), -2, llev+1)
     END IF
   END DO
   DO i = 1, n
     IF (cell(i) <= cell(i-1)) THEN
       CALL PrintLog('check cells..., index:'//toString(i), -2, llev+1)
     END IF
   END DO

   DO i = 1, n/2
     IF (cell(i)-cell(i-1) /= cell(n+1-i)-cell(n+1-i-1)) THEN
       CALL PrintLog('check delta cells..., index:'//toString(i), -1, llev+1)
       PRINT *, cell(i)-cell(i-1), cell(n+1-i)-cell(n+1-i-1)
     END IF
   END DO

   CALL PrintLog('END  : Check position of points...', elev, llev)
   CALL PrintLog(' ', elev)


 END SUBROUTINE CheckGridCell


END PROGRAM TestMain
