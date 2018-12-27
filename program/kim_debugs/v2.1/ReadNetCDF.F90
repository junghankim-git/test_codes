PROGRAM Test
 USE Kinds, ONLY: i4, r4, r8
 USE NetCDF
 USE NetCDF_Mod, ONLY: netcdf_file_t, OpenFile, CloseFile, GetDimSize, ReadVariable
 IMPLICIT NONE

 TYPE(netcdf_file_t)    :: ncfile
 CHARACTER(LEN=128)     :: filename_conf, filename1, filename2
 CHARACTER(LEN=128)     :: outfilename
 INTEGER(i4) :: ncid
 INTEGER(i4) :: ncol_dimid, one_id, nlev_dimid
 INTEGER(i4) :: lon_varid, lat_varid, t1_varid, t2_varid, diff_T_varid, diff_ps_varid, points_varid

 INTEGER(i4) :: nhoriz, nlev
 REAL(r8), DIMENSION(:), ALLOCATABLE   :: lat, lon
 REAL(r8), DIMENSION(:), ALLOCATABLE   :: ps1, ps2, diff_ps
 REAL(r8), DIMENSION(:,:), ALLOCATABLE :: T1, T2, diff_T
 REAL(r8), DIMENSION(:,:), ALLOCATABLE :: theta_c, theta_w
 REAL(r8), DIMENSION(:), ALLOCATABLE   :: points
 INTEGER(i4) :: ierr, i, ilev, ncount
 REAL(r8)    :: p0, R, cp, err, sum1, sum2, ave1, ave2, mx1, mx2, mn1, mn2
 LOGICAL     :: isIni1, isIni2

 p0   = 100000.0_r8
 R    = 287.04_r8
 cp   = 1005.0_r8
 nlev = 50

 WRITE(filename_conf, '(A)') '/data/jpark/2.1.01/KIM/exp_ne30np4/run.ctl/UP-20151004120000-000001.nc'
 !WRITE(filename_conf, '(A)') '/home/yslee/warm_debugging/UP-20151004000000-000002.nc'

 !  ATM1: /home/yslee/warm_debugging/20151003_18UTC_SW_GFS_50lev_ne120np4.nc
 !  ATM2: /data/KIM/inputdata/ne120/vgecore/atm/Init_2015100400_SW_50lev_ne120np4.nc
 !  ATM3: /data/KIM/inputdata/ne120/vgecore/atm/20151004_00UTC_SW_GFS_50lev_ne120np4.nc
 ! out00: /home/yslee/warm_debugging/UP-20151004000000-000002.nc


#if 0
 ! check 1
 ! ATM2
 WRITE(filename1, '(A)') '/data/KIM/inputdata/ne120/vgecore/atm/Init_2015100400_SW_50lev_ne120np4.nc'
 isIni1 = .TRUE.
 ! ATM3
 WRITE(filename2, '(A)') '/data/KIM/inputdata/ne120/vgecore/atm/20151004_00UTC_SW_GFS_50lev_ne120np4.nc'
 isIni2 = .TRUE.

 ! check 2
 ! ATM2
 WRITE(filename1, '(A)') '/data/KIM/inputdata/ne120/vgecore/atm/Init_2015100400_SW_50lev_ne120np4.nc'
 isIni1 = .TRUE.
 ! out00
 WRITE(filename2, '(A)') '/home/yslee/warm_debugging/UP-20151004000000-000002.nc'
 isIni2 = .FALSE.

 ! check 3
 ! out00
 WRITE(filename1, '(A)') '/home/yslee/warm_debugging/UP-20151004000000-000002.nc'
 isIni1 = .FALSE.
 ! ATM1
 WRITE(filename2, '(A)') '/home/yslee/warm_debugging/20151003_18UTC_SW_GFS_50lev_ne120np4.nc'
 isIni2 = .TRUE.

 ! check 4-1: v0.26.04
 ! out18
 WRITE(filename1, '(A)') '/home/yslee/warm_debugging/0.26.04/KIM/exp_ne120np4/run.ctl/UP-20151004000000-000003.nc'
 isIni1 = .FALSE.
 ! out00_1
 !WRITE(filename2, '(A)') '/home/yslee/warm_debugging/0.26.04/KIM/exp_ne120np4/run.ctl/UP-20151003180000-000001.nc'
 !isIni2 = .FALSE.
 WRITE(filename2, '(A)') '/home/yslee/warm_debugging/20151003_18UTC_SW_GFS_50lev_ne120np4.nc'
 isIni2 = .TRUE.

 ! check 4-2: v0.26.01
 ! out18
 WRITE(filename1, '(A)') '/home/yslee/warm_debugging/0.26.01/KIM/exp_ne120np4/history/UP-20151004000000-000004.nc'
 isIni1 = .FALSE.
 ! out00_1
 !WRITE(filename2, '(A)') '/home/yslee/warm_debugging/0.26.04/KIM/exp_ne120np4/run.ctl/UP-20151003180000-000001.nc'
 !isIni2 = .FALSE.
 WRITE(filename2, '(A)') '/home/yslee/warm_debugging/20151003_18UTC_SW_GFS_50lev_ne120np4.nc'
 isIni2 = .TRUE.

 ! check 4-3: v0.25.01
 ! out18
 WRITE(filename1, '(A)') '/home/yslee/warm_debugging/0.25.01/KIM/exp_ne120np4/history/UP-20151004000000-000002.nc'
 isIni1 = .FALSE.
 ! out00_1
 !WRITE(filename2, '(A)') '/home/yslee/warm_debugging/0.26.04/KIM/exp_ne120np4/run.ctl/UP-20151003180000-000001.nc'
 !isIni2 = .FALSE.
 WRITE(filename2, '(A)') '/home/yslee/warm_debugging/20151003_18UTC_SW_GFS_50lev_ne120np4.nc'
 isIni2 = .TRUE.
#endif

 ! check 5: 100312 - 100318
 WRITE(filename1, '(A)') '/home/yslee/warm_debugging/100312.cold/UP-20151003120000-000001.nc'
 isIni1 = .FALSE.
 ! out00_1
 WRITE(filename2, '(A)') '/home/yslee/warm_debugging/100312.cold/UP-20151003180000-000002.nc'
 !WRITE(filename2, '(A)') '/home/yslee/warm_debugging/100312.cold/UP-20151004060000-000004.nc'
 isIni2 = .FALSE.

 ! check 6: 100318 - 100400 (sfco: 1004 00UTC)
 ! out00
 WRITE(filename1, '(A)') '/scratch/jhkim/TestBed/Data/2.1_bug/UP-20151004000000-000002.nc'
 isIni1 = .FALSE.
 ! ATM1
 WRITE(filename2, '(A)') '/home/yslee/warm_debugging/20151003_18UTC_SW_GFS_50lev_ne120np4.nc'
 isIni2 = .TRUE.


#if 0
 ! check 3
 ! out00
 WRITE(filename1, '(A)') '/home/yslee/warm_debugging/20151003_18UTC_SW_GFS_50lev_ne120np4.nc'
 isIni1 = .TRUE.
 WRITE(filename2, '(A)') '/home/yslee/warm_debugging/UP-20151004000000-000002.nc'
 isIni2 = .FALSE.
 ! ATM1

 ! check 4-1: v0.26.04
 ! out18
 WRITE(filename1, '(A)') '/home/yslee/warm_debugging/20151003_18UTC_SW_GFS_50lev_ne120np4.nc'
 isIni1 = .TRUE.
 WRITE(filename2, '(A)') '/home/yslee/warm_debugging/0.26.04/KIM/exp_ne120np4/run.ctl/UP-20151004000000-000003.nc'
 isIni2 = .FALSE.

 ! new sst, seaice
 ! out18
 WRITE(filename1, '(A)') '/home/yslee/warm_debugging/20151003_18UTC_SW_GFS_50lev_ne120np4.nc'
 isIni1 = .TRUE.
 WRITE(filename2, '(A)') '/home/yslee/warm_debugging/2.1/KIM/exp_ne120np4/run.ctl/UP-20151004000000-000002.nc'
 isIni2 = .FALSE.
#endif


 ! Read Configuration: START
 CALL OpenFile(ncfile, filename_conf, ierr)
! CALL IniFile(ncfile, ierr)
 nhoriz = GetDimSize(ncfile, "ncol", ierr)
 nlev   = GetDimSize(ncfile, "lev", ierr)
 ALLOCATE(lat(nhoriz))
 ALLOCATE(lon(nhoriz))
 CALL ReadVariable(ncfile, "lat", lat, ierr)
 CALL ReadVariable(ncfile, "lon", lon, ierr)
! CALL FinFile(ncfile, ierr)
 CALL CloseFile(ncfile, ierr)
 ! Read Configuration: END

 ALLOCATE(ps1(nhoriz))
 ALLOCATE(ps2(nhoriz))
 ALLOCATE(diff_ps(nhoriz))
 ALLOCATE(T1(nhoriz,nlev))
 ALLOCATE(T2(nhoriz,nlev))
 ALLOCATE(diff_T(nhoriz,nlev))
 ALLOCATE(theta_c(nhoriz, nlev))
 ALLOCATE(theta_w(nhoriz, nlev))
 ALLOCATE(points(nhoriz))
   

 ! Read Cold: START
 CALL OpenFile(ncfile, filename1, ierr)
! CALL IniFile(ncfile, ierr)
 CALL RestoreVar1D(ncfile, "ps", nhoriz, isIni1, ps1)
 CALL RestoreVar(ncfile, "T",  nhoriz, nlev, isIni1, ps1, T1)
! CALL FinFile(ncfile, ierr)
 CALL CloseFile(ncfile, ierr)
 ! Read Cold: END

 ! Read Warm: START
 CALL OpenFile(ncfile, filename2, ierr)
! CALL IniFile(ncfile, ierr)
 CALL RestoreVar1D(ncfile, "ps", nhoriz, isIni2, ps2)
 CALL RestoreVar(ncfile, "T",  nhoriz, nlev, isIni2, ps2, T2)
! CALL FinFile(ncfile, ierr)
 CALL CloseFile(ncfile, ierr)
 ! Read Warm: END

 ncount = 0
 sum1 = 0.0_r8
 sum2 = 0.0_r8
 mx1  = -1.0D-30
 mx2  = -1.0D-30
 mn1  = 1.0D+30
 mn2  = 1.0D+30
 points(:) = 0.0_r8
 DO i = 1, nhoriz
   diff_ps(i) = ps2(i) - ps1(i)
   diff_T(i,1)  = T2(i,1)  - T1(i,1)
   IF (lat(i) < -65.0_r8) THEN
     ncount = ncount + 1
     points(i) = 1.0_r8
     sum1 = sum1 + T1(i,1)
     sum2 = sum2 + T2(i,1)
     mx1 = MAX(mx1, T1(i,1))
     mx2 = MAX(mx2, T2(i,1))
     mn1 = MIN(mn1, T1(i,1))
     mn2 = MIN(mn2, T2(i,1))
   END IF
 END DO
 ave1 = sum1/DBLE(ncount)
 ave2 = sum2/DBLE(ncount)

 PRINT *, ' ============================== '
 PRINT *, ' T1 = ', ave1, mn1, mx1
 PRINT *, ' T2 = ', ave2, mn2, mx2
 PRINT *, ' ============================== '



 outfilename = './anal.nc'

 ierr = NF90_CREATE(outfilename, NF90_CLOBBER, ncid)
 ierr = NF90_DEF_DIM(ncid, "ncol", nhoriz, ncol_dimid)
 ierr = NF90_DEF_DIM(ncid, "nlev", nlev,   nlev_dimid)
 ierr = NF90_DEF_DIM(ncid, "one", 1, one_id)
 !ierr = NF90_CLOSE(ncid)

 ierr = NF90_DEF_VAR(ncid, 'lon',     NF90_REAL8, ncol_dimid, lon_varid    )
 ierr = NF90_DEF_VAR(ncid, 'lat',     NF90_REAL8, ncol_dimid, lat_varid    )
 ierr = NF90_DEF_VAR(ncid, 'T1',      NF90_REAL8, ncol_dimid, t1_varid     )
 ierr = NF90_DEF_VAR(ncid, 'T2',      NF90_REAL8, ncol_dimid, t2_varid     )
 ierr = NF90_DEF_VAR(ncid, 'diff_T',  NF90_REAL8, ncol_dimid, diff_T_varid )
 ierr = NF90_DEF_VAR(ncid, 'diff_ps', NF90_REAL8, ncol_dimid, diff_ps_varid)
 ierr = NF90_DEF_VAR(ncid, 'points',  NF90_REAL8, ncol_dimid, points_varid )
 ierr = NF90_ENDDEF(ncid)

 !ierr = NF90_CREATE(outfilename, NF90_CLOBBER, ncid)
 ierr = NF90_PUT_VAR(ncid, lon_varid, lon)
 ierr = NF90_PUT_VAR(ncid, lat_varid, lat)
 ierr = NF90_PUT_VAR(ncid, t1_varid, T1(:,1))
 ierr = NF90_PUT_VAR(ncid, t2_varid, T2(:,1))
 ierr = NF90_PUT_VAR(ncid, diff_T_varid, diff_T(:,1))
 ierr = NF90_PUT_VAR(ncid, diff_ps_varid, diff_ps(:))
 ierr = NF90_PUT_VAR(ncid, points_varid, points(:))

 ierr = NF90_CLOSE(ncid)

 DEALLOCATE(lat)
 DEALLOCATE(lon)
 DEALLOCATE(ps1)
 DEALLOCATE(ps2)
 DEALLOCATE(diff_ps)
 DEALLOCATE(T1)
 DEALLOCATE(T2)
 DEALLOCATE(diff_T)
 DEALLOCATE(theta_c)
 DEALLOCATE(theta_w)
 DEALLOCATE(points)

 CONTAINS




 SUBROUTINE RestoreVar1D(nctype, varname, nx, isIni, var)
  IMPLICIT NONE
  TYPE(netcdf_file_t), INTENT(INOUT) :: nctype
  CHARACTER(LEN=*), INTENT(IN)       :: varname
  INTEGER(i4), INTENT(IN)            :: nx
  LOGICAL, INTENT(IN)                :: isIni
  REAL(r8), DIMENSION(nx)            :: var

  ! local
  INTEGER(i4) :: l_ierr

    IF (isIni) THEN
      CALL ReadVariable(nctype, "psfc", var, l_ierr)
    ELSE
      CALL ReadVariable(nctype, "ps",   var, l_ierr)
    END IF

 END SUBROUTINE RestoreVar1D




 SUBROUTINE RestoreVar(nctype, varname, nx, ny, isIni, ps_in, var, istheta)
  IMPLICIT NONE
  TYPE(netcdf_file_t), INTENT(INOUT) :: nctype
  CHARACTER(LEN=*), INTENT(IN)       :: varname
  INTEGER(i4), INTENT(IN)            :: nx, ny
  LOGICAL, INTENT(IN)                :: isIni
  REAL(r8), DIMENSION(nx, 1), INTENT(IN)     :: ps_in
  REAL(r8), DIMENSION(nx, ny), INTENT(INOUT) :: var
  LOGICAL, OPTIONAL, INTENT(IN)      :: istheta

  ! local
  INTEGER(i4) :: i, j, l_ierr
  REAL(r8), DIMENSION(nx, ny)   :: var_tmp

    IF (isIni) THEN
      CALL ReadVariable(nctype, "theta", var_tmp, l_ierr)
      DO j = 1, ny
        DO i = 1, nx
          var(i,j) = var_tmp(i, ny-j+1)
        END DO
      END DO
      DO i = 1, nx
        var(i,1)    = var(i,1)*(p0/ps_in(i,1))**(-R/cp)
      END DO
    ELSE
      CALL ReadVariable(nctype, "T", var, l_ierr)
    END IF

 END SUBROUTINE RestoreVar


END PROGRAM Test
