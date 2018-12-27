PROGRAM TestMain
 USE Kinds,          ONLY: i4, l4, r8
 USE Logger,         ONLY: PrintLog, toString
 USE IoASCII,       ONLY: OpenFile, CloseFile, WriteData_1Log1, WriteData_1Int1, WriteData_1Real, WriteData_String
 USE Controls,       ONLY: np, ne, ntimelevs, dt, velocity
 USE Element,        ONLY: element_t, nUPs
 USE Element,        ONLY: IniElement, FinElement, PrintElement, DiffElement, SetPosition, SetState
 USE Element,        ONLY: ElementsToUPs_Pos, ElementsToUPs_State, UPsToElements_Pos, UPsToElements_State
 USE SemiLagrangian, ONLY: SemiLagrangian_t
 USE SemiLagrangian, ONLY: IniSemiLagrangian, FinSemiLagrangian, GenControlVolume, RunSemiLagrangian
 IMPLICIT NONE
 INTEGER(i4)  :: unitnum
 INTEGER(i4)  :: elev, llev
 REAL(r8)     :: xmin, xmax
 TYPE(element_t), DIMENSION(:), ALLOCATABLE :: elem, elem_out
 TYPE(SemiLagrangian_t) :: sl
! REAL(r8), DIMENSION(:), ALLOCATABLE :: xgrid
! REAL(r8), DIMENSION(:), ALLOCATABLE :: xcell
! REAL(r8), DIMENSION(:), ALLOCATABLE :: v, T

  ! logger setting
  elev = 0
  llev = 2

  ! domain
  xmin = -1.0_r8
  xmax =  1.0_r8

  ! np, ne
  ne   = 1
  np   = 4

  CALL ReadNameList('inputs.nl')
  CALL IniElement(elem)
  CALL IniSemiLagrangian(sl)
!  ALLOCATE(xgrid(nUPs))
!  ALLOCATE(xcell(0:nUPs))
!  ALLOCATE(v(nUPs))
!  ALLOCATE(T(nUPs))

  CALL SetPosition(elem, xmin, xmax)
  CALL SetState(elem)
  CALL PrintElement(elem)
  CALL GenControlVolume(sl, elem)


  CALL RunSemiLagrangian(sl, elem, dt)


  CALL OpenFile('Result.dat', unitnum)
  CALL WriteData_1Int1(unitnum, np)
  CALL WriteData_1Int1(unitnum, ne)
  CALL WriteData_1Int1(unitnum, nUPs)
  CALL WriteData_1Real(unitnum, nUPs, sl%grid)
  CALL WriteData_1Real(unitnum, nUPs+1, sl%cell)
  CALL CloseFile(unitnum)


! START: back test
!  CALL IniElement(elem_out)
!  CALL ElementsToUPs_State(elem, v, 1)
!  CALL ElementsToUPs_State(elem, T, 2)
!  CALL UPsToElements_Pos(elem_out, xgrid, xcell)
!  CALL UPsToElements_State(elem_out, v, 1)
!  CALL UPsToElements_State(elem_out, T, 2)
!  CALL DiffElement(elem, elem_out)
!  CALL FinElement(elem_out)
! END  : back test


!  DEALLOCATE(xgrid)
!  DEALLOCATE(xcell)
!  DEALLOCATE(v)
!  DEALLOCATE(T)
  CALL FinSemiLagrangian(sl)
  CALL FinElement(elem)

 CONTAINS




 SUBROUTINE ReadNameList(nlFileName)
  IMPLICIT NONE
   CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: nlFileName
   NAMELIST /inputs/ np, ne, dt, velocity

   CALL PrintLog('START: Read Namelist...', elev, llev)
   IF(PRESENT(nlFileNAme)) THEN
     OPEN(21, FILE=TRIM(nlFileName), STATUS='old')
     READ(21, NML=inputs)
     CLOSE(21)
   ELSE
     READ(*, NML=inputs)
   END IF

   CALL PrintLog('ne        = '//toString(ne), elev, llev+1)
   CALL PrintLog('np        = '//toString(np), elev, llev+1)
   CALL PrintLog('ntimelevs = '//toString(ntimelevs), elev, llev+1)
   CALL PrintLog('dt        = '//toString(dt), elev, llev+1)
   CALL PrintLog('velocity  = '//toString(velocity), elev, llev+1)
   CALL PrintLog('END  : Read Namelist...', elev, llev)

 END SUBROUTINE ReadNameList




END PROGRAM TestMain
