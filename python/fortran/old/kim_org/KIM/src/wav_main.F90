!-------------------------------------------------------------------------------
! MODULE wav_main
!
!> @brief
!>  - Wave model's main (test)
!>
!> @date 18MAR2016
!>  - Junghan Kim : First written (KIM-CPL v0.25)
!> @date 12APR2016
!>  - Junghan Kim : The version for application to the KIM (KIM-CPL v1.0)
!-------------------------------------------------------------------------------

MODULE wav_main
 USE KindsCPL,     ONLY: i4, r8, l4
 USE ParallelCPL,  ONLY: Parallel_t, CopyCommunicator, Decompose1D, CreateDOF
 USE Coupler,      ONLY: nComps
 USE Auxiliary,    ONLY: CreateX, CreateVar, PrintGlobalVar
 USE NetCDF_CPL,   ONLY: netcdf_t, CreateNCFile, CloseNCFile, AddDimension, DefVariable, AddVariable
 IMPLICIT NONE

 PRIVATE

 ! Component
 TYPE(Parallel_t)  :: par

 ! Grid
 INTEGER(i4) :: gvarsize, lvarsize
 LOGICAL(l4) :: usecontinuous
 INTEGER(i4), DIMENSION(:), ALLOCATABLE :: DOFs
 REAL(r8), DIMENSION(:), ALLOCATABLE    :: wav_u10m, wav_v10m
! REAL(r8), DIMENSION(:), ALLOCATABLE    :: wav_x, wav_u, wav_v, wav_du, wav_dv

 ! namelist
 CHARACTER(LEN=512)                     :: nlfilename

 ! test
 TYPE(netcdf_t) :: my_nc

 PUBLIC :: Set, Ini, Run, Fin
 PUBLIC :: lvarsize, gvarsize, DOFs
 PUBLIC :: wav_u10m, wav_v10m


 CONTAINS



 SUBROUTINE Set(par_in, nlfilename_in)
  IMPLICIT NONE
  TYPE(Parallel_t), INTENT(IN) :: par_in
  CHARACTER(LEN=*), INTENT(IN) :: nlfilename_in

   CALL CopyCommunicator(par_in, par)
   nlfilename = nlfilename_in

 END SUBROUTINE Set




 SUBROUTINE ReadNL()
  IMPLICIT NONE
  NAMELIST /wav/ gvarsize, usecontinuous

   OPEN(21, FILE=TRIM(nlfilename), STATUS='old')
   READ(21, NML=wav)
   CLOSE(21)

 END SUBROUTINE ReadNL





 SUBROUTINE Ini()
  IMPLICIT NONE
  INTEGER(i4) :: ista, iend, i

   ! Domain decomposition and creat vars
   CALL ReadNL()
   lvarsize = Decompose1D(1, gvarsize, par%nprocs, par%rank, ista, iend)

   IF (usecontinuous) THEN
     DOFs    = CreateDOF(ista, iend)
   ELSE
     DOFs    = CreateDOF(gvarsize, par%rank, par%nprocs)
   END IF
   ALLOCATE(wav_u10m(lvarsize))
   ALLOCATE(wav_v10m(lvarsize))
   wav_u10m(:) = 0.0_r8
   wav_v10m(:) = 0.0_r8

   !test
   CALL CreateNCFile(my_nc, par, "file_wav.nc")
   CALL AddDimension(my_nc, "ncols", gvarsize)
   CALL DefVariable(my_nc, "wav_u10m")
!   CALL AddVariable(my_nc, "wav_u10m", wav_u10m, DOFs)

 END SUBROUTINE Ini






 SUBROUTINE Run()
  IMPLICIT NONE
  INTEGER(i4) :: i

   CALL AddVariable(my_nc, "wav_u10m", wav_u10m, DOFs)
!   DO i = 1, lvarsize
!     wav_u(i)  = wav_u(i) + wav_du(i)
!     wav_v(i)  = wav_v(i) + wav_dv(i)
!     wav_du(i) = -1.0_r8*wav_du(i)
!     wav_dv(i) = -1.0_r8*wav_dv(i)
!   END DO

 END SUBROUTINE Run






 SUBROUTINE Fin()
  IMPLICIT NONE

   DEALLOCATE(DOFs)
   DEALLOCATE(wav_u10m)
   DEALLOCATE(wav_v10m)

   ! test
   CALL CloseNCFile(my_nc)

 END SUBROUTINE Fin






END MODULE wav_main
