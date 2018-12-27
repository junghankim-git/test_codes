!-------------------------------------------------------------------------------
! MODULE WaveComp
!
!> @brief
!>  - Wave model component's driver
!>
!> @date 18MAR2016
!>  - Junghan Kim : First written (KIM-CPL v0.25)
!> @date 15APR2016
!>  - Junghan Kim : The version for application to the KIM (KIM-CPL v1.0)
!-------------------------------------------------------------------------------
#undef DEBUG

MODULE wav_comp
 ! Kinds
 USE KindsCPL,     ONLY: i4, r8, i4
 ! Shared
 USE Controls,     ONLY: lunits
 USE ParallelCPL,  ONLY: Decompose1D, SyncPar, GetStartCount
 USE LoggerCPL,    ONLY: elev, llev, PrintLog, toString, CreateLogFile, CloseLogFile
 ! Coupler
 USE Coupler,      ONLY: id_atm, id_wav, id_ocn, id_cpl, nComps
 USE Coupler,      ONLY: wav_nlfile
 USE Coupler,      ONLY: useStridedMAP
 USE Coupler,      ONLY: Component_t, CopyComponent
 USE Coupler,      ONLY: Coupler_t
 USE Coupler,      ONLY: SetGlobalSize,     SetMap
 USE Coupler,      ONLY: CreateCoupler,     DeleteCoupler
 USE Coupler,      ONLY: IniCouplerWorld,   FinCouplerWorld
 USE Coupler,      ONLY: IniEnvironment,    FinEnvironment
 USE Coupler,      ONLY: ExchangeVariables, WaitExchangeVariables
 USE Coupler,      ONLY: PutVariable,       GetVariable
 ! Model
 USE wav_main,     ONLY: lvarsize, gvarsize
 USE wav_main,     ONLY: DOFs
 USE wav_main,     ONLY: SetWaveMain => Set, IniWaveMain => Ini, RunWaveMain => Run, FinWaveMain => Fin
  ! recv variables
 USE wav_main,     ONLY: wav_u10m, wav_v10m
 IMPLICIT NONE

 PRIVATE

 ! Component & Coupler
 TYPE(Component_t)                          :: comp
 TYPE(Coupler_t), DIMENSION(:), ALLOCATABLE :: cpls
 INTEGER(i4)                                :: myid

 ! LoggerCPL
 INTEGER(i4) :: log_unit

 ! APIs
 PUBLIC :: Set, Ini, IniCPL, ExchangeCPL, Run, Fin


 CONTAINS



 SUBROUTINE Set(comp_in, id_in)
  IMPLICIT NONE
  TYPE(Component_t), INTENT(IN) :: comp_in
  INTEGER(i4),       INTENT(IN) :: id_in
  CHARACTER(LEN=512) :: logfilename

   !---------------------------------------------------------------------------
   ! Component
   CALL CopyComponent(comp_in, comp)
   myid = id_in
   log_unit = lunits(myid)
   WRITE(logfilename, '(A)'), 'wav.log'

   !---------------------------------------------------------------------------
   ! Coupler
   IF (comp%iscoupling) THEN

     CALL IniCouplerWorld(comp)

   END IF
   !---------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! LoggerCPL
   CALL CreateLogFile(comp%par, TRIM(logfilename), log_unit)
   CALL PrintLog(comp%par, 'Set WaveComp', elev, llev, un=log_unit)

   !---------------------------------------------------------------------------
   ! WaveMain
   CALL SetWaveMain(comp%par, wav_nlfile)

 END SUBROUTINE Set





 SUBROUTINE Ini()
  IMPLICIT NONE
  ! local

   CALL PrintLog(comp%par, 'Ini WaveComp', elev, llev, un=log_unit)

   !---------------------------------------------------------------------------
   ! WaveMain
   CALL IniWaveMain()
   !---------------------------------------------------------------------------

   IF (comp%iscoupling) THEN
     ! Send global size to Coupler     (User)
     CALL SetGlobalSize(comp, myid, gvarsize)
   END IF

 END SUBROUTINE Ini




 SUBROUTINE IniCPL()
  IMPLICIT NONE
  ! local
  INTEGER(i4), DIMENSION(:), ALLOCATABLE :: starts, counts

   !---------------------------------------------------------------------------
   ! Coupler
   IF (comp%iscoupling) THEN

     CALL PrintLog(comp%par, 'Ini Coupling', elev, llev, un=log_unit)

     ! Allocate Coupler derived type
     CALL CreateCoupler(comp, cpls)

     ! DOF's map (grid impormations)   (User)
     CALL GetStartCount(DOFs, starts, counts, useStridedMAP)
     CALL SetMap(comp, cpls, 1, starts, counts)
     DEALLOCATE(starts, counts)

     ! Initialize CPL environment
     CALL IniEnvironment(comp, cpls)

   END IF
   !---------------------------------------------------------------------------


 END SUBROUTINE IniCPL






 SUBROUTINE ExchangeCPL(itime)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN) :: itime

   !---------------------------------------------------------------------------
   ! Coupler
   IF (comp%iscoupling) THEN

     CALL PrintLog(comp%par, 'Do Coupling  ('//toString(itime)//')', elev, llev, un=log_unit)

     !---------------------------------------------------------------------------
     != User's variables (Put -> Exchange(send, recv) -> Get)
     !---------------------------------------------------------------------------
     ! STEP1: Put (user)
!     CALL PutVariable(comp, cpls, 1, id_atm, 'du_wav', wav_du, itime)
!     CALL PutVariable(comp, cpls, 1, id_atm, 'dv_wav', wav_dv, itime)
  
     ! STEP2: Exchange Variables
     CALL ExchangeVariables(comp, cpls, itime)

   END IF
   !---------------------------------------------------------------------------

 END SUBROUTINE ExchangeCPL







 SUBROUTINE Run(itime)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN) :: itime

   CALL PrintLog(comp%par, 'Run WaveComp  ('//toString(itime)//')', elev, llev, un=log_unit)

   IF (comp%iscoupling) THEN
     CALL WaitExchangeVariables(comp, cpls, itime)

     ! STEP3: Get (user)
     CALL GetVariable(comp, cpls, 1, id_atm, 'u10m', wav_u10m, itime)
     CALL GetVariable(comp, cpls, 1, id_atm, 'v10m', wav_v10m, itime)
   END IF

   !---------------------------------------------------------------------------
   ! WaveMain
   CALL RunWaveMain()
   !---------------------------------------------------------------------------

 END SUBROUTINE Run






 SUBROUTINE Fin()
  IMPLICIT NONE

   CALL PrintLog(comp%par, 'Fin WaveComp', elev, llev, un=log_unit)

   !---------------------------------------------------------------------------
   ! WaveMain
   CALL FinWaveMain()
   !---------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! Coupler
   IF (comp%iscoupling) THEN

     CALL FinEnvironment(comp, cpls)
     CALL DeleteCoupler(comp, cpls)
!     CALL FinCouplerWorld(comp) => FinEnvironment

   END IF
   !---------------------------------------------------------------------------

   CALL CloseLogFile(comp%par, log_unit)

 END SUBROUTINE Fin




END MODULE wav_comp
