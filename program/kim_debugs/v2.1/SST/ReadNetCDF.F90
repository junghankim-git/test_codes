PROGRAM Test
 USE Kinds, ONLY: i4, r4, r8
 USE NetCDF_Mod
 !USE NetCDF_Mod, ONLY: netcdf_file_t, OpenFile, CloseFile, GetDimSize, ReadVariable
 !USE NetCDF_Mod, ONLY: CreateFile, AddDimension, AddVariable, WriteVariable
 IMPLICIT NONE

 TYPE(netcdf_file_t)    :: ncfile, out_ncfile
 CHARACTER(LEN=128)     :: infile_conf, infile_sfco, infile_sst, infile_lsmsk, infile_seaice
 CHARACTER(LEN=128)     :: outfilename
 INTEGER(i4) :: ncid
 INTEGER(i4) :: ncol_dimid, one_id, nlev_dimid
 INTEGER(i4) :: lon_varid, lat_varid, sfco_varid, sst_varid, asst_varid, ini_varid, seaice_varid, lsmsk_varid, points_varid

 INTEGER(i4) :: nhoriz, nlev
 REAL(r8), DIMENSION(:), ALLOCATABLE    :: lat, lon
 INTEGER(i4), DIMENSION(:), ALLOCATABLE :: lsmsk
 REAL(r8), DIMENSION(:), ALLOCATABLE    :: sfco, sst, asst, seaice, ini
 REAL(r8), DIMENSION(:), ALLOCATABLE    :: points
 INTEGER(i4) :: ierr, i, ilev, ncount
 REAL(r8)    :: p0, R, cp, err, sum1, sum2, ave1, ave2, mx1, mx2, mn1, mn2
 LOGICAL     :: isIni1, isIni2

 p0   = 100000.0_r8
 R    = 287.04_r8
 cp   = 1005.0_r8
 nlev = 50

 WRITE(infile_conf,   '(A)') '/data/jpark/2.1.01/KIM/exp_ne30np4/run.ctl/UP-20151004120000-000001.nc'

 WRITE(infile_sfco,   '(A)') '/data/KIM/inputdata/ne120/vgecore/sfco/sfco_2015100318_ne120np4.nc'
 WRITE(infile_lsmsk,  '(A)') '/data/KIM/inputdata/ne120/vgecore/lsmsk.nc'
 ! org
 !WRITE(infile_sst,    '(A)') '/data/KIM/inputdata/ne120/vgecore/sst/UM_OSTIA_sst_20151003_ne120np4.nc.org'
 !WRITE(infile_seaice, '(A)') '/data/KIM/inputdata/ne120/vgecore/seaice/UM_OSTIA_seaice_20151003_ne120np4.nc.org'
 ! new
 WRITE(infile_sst,    '(A)') '/data/KIM/inputdata/ne120/vgecore/sst/UM_OSTIA_sst_20151003_ne120np4.nc'
 WRITE(infile_seaice, '(A)') '/data/KIM/inputdata/ne120/vgecore/seaice/UM_OSTIA_seaice_20151003_ne120np4.nc'



 ! Read Configuration: START
 CALL OpenFile(ncfile, infile_conf, ierr)
 nhoriz = GetDimSize(ncfile, "ncol", ierr)
 nlev   = GetDimSize(ncfile, "lev", ierr)
 ALLOCATE(lat(nhoriz))
 ALLOCATE(lon(nhoriz))
 CALL ReadVariable(ncfile, "lat", lat, ierr)
 CALL ReadVariable(ncfile, "lon", lon, ierr)
 CALL CloseFile(ncfile, ierr)
 ! Read Configuration: END

 ALLOCATE(sfco(nhoriz))
 ALLOCATE(sst(nhoriz))
 ALLOCATE(asst(nhoriz))
 ALLOCATE(ini(nhoriz))
 ALLOCATE(lsmsk(nhoriz))
 ALLOCATE(seaice(nhoriz))
 ALLOCATE(points(nhoriz))
   

 ! SFCO
 CALL OpenFile(ncfile, infile_sfco, ierr)
 CALL ReadVariable(ncfile, "tsfc", sfco, ierr)
 CALL CloseFile(ncfile, ierr)


 ! SST
 CALL OpenFile(ncfile, infile_sst, ierr)
 CALL ReadVariable(ncfile, "sst", sst, ierr)
 CALL CloseFile(ncfile, ierr)

 ! lsmsk
 CALL OpenFile(ncfile, infile_lsmsk, ierr)
 CALL ReadVariable(ncfile, "lsmsk", lsmsk, ierr)
 CALL CloseFile(ncfile, ierr)


 ! seaice
 CALL OpenFile(ncfile, infile_seaice, ierr)
 CALL ReadVariable(ncfile, "icec", seaice, ierr)
 CALL CloseFile(ncfile, ierr)



  points(:) = 0.0_r8
  asst(:)   = 0.0_r8
  ini(:)    = sfco(:)
  DO i =  1, nhoriz
    !IF (slimsk(i,1) == 0.0_r8) THEN
    IF (lsmsk(i) == 0._r8 .AND. seaice(i) < 0.5_r8) THEN
      ini(i)    = sst(i)
      points(i) = 1.0_r8
      asst(i)   = sst(i)
    END IF
  END DO
! ncount = 0
! sum1 = 0.0_r8
! sum2 = 0.0_r8
! mx1  = -1.0D-30
! mx2  = -1.0D-30
! mn1  = 1.0D+30
! mn2  = 1.0D+30
! points(:,:) = 0.0_r8
! DO i = 1, nhoriz
!   lsmsk(i,1) = sst(i,1) - sfco(i,1)
!   diff_T(i,1)  = T2(i,1)  - seaice(i,1)
!   IF (lat(i,1) < -65.0_r8) THEN
!     ncount = ncount + 1
!     points(i,1) = 1.0_r8
!     sum1 = sum1 + seaice(i,1)
!     sum2 = sum2 + T2(i,1)
!     mx1 = MAX(mx1, seaice(i,1))
!     mx2 = MAX(mx2, T2(i,1))
!     mn1 = MIN(mn1, seaice(i,1))
!     mn2 = MIN(mn2, T2(i,1))
!   END IF
! END DO
! ave1 = sum1/DBLE(ncount)
! ave2 = sum2/DBLE(ncount)
!
! PRINT *, ' ============================== '
! PRINT *, ' seaice = ', ave1, mn1, mx1
! PRINT *, ' T2 = ', ave2, mn2, mx2
! PRINT *, ' ============================== '
!
!
!


#if 0
 outfilename = './anal.nc'

 ierr = NF90_CREATE(outfilename, NF90_CLOBBER, ncid)
 ierr = NF90_DEF_DIM(ncid, "ncol", nhoriz, ncol_dimid)
 ierr = NF90_DEF_DIM(ncid, "nlev", nlev,   nlev_dimid)
 ierr = NF90_DEF_DIM(ncid, "one", 1, one_id)
 !ierr = NF90_CLOSE(ncid)

 ierr = NF90_DEF_VAR(ncid, 'lon',    NF90_REAL8, ncol_dimid, lon_varid)
 ierr = NF90_DEF_VAR(ncid, 'lat',    NF90_REAL8, ncol_dimid, lat_varid)
 ierr = NF90_DEF_VAR(ncid, 'sfco',   NF90_REAL8, ncol_dimid, sfco_varid)
 ierr = NF90_DEF_VAR(ncid, 'sst',    NF90_REAL8, ncol_dimid, sst_varid)
 ierr = NF90_DEF_VAR(ncid, 'asst',   NF90_REAL8, ncol_dimid, asst_varid)
 ierr = NF90_DEF_VAR(ncid, 'ini',    NF90_REAL8, ncol_dimid, ini_varid)
 ierr = NF90_DEF_VAR(ncid, 'lsmsk',  NF90_INT, ncol_dimid, lsmsk_varid)
 ierr = NF90_DEF_VAR(ncid, 'seaice', NF90_REAL8, ncol_dimid, seaice_varid)
 ierr = NF90_DEF_VAR(ncid, 'points',  NF90_REAL8, ncol_dimid, points_varid )
 ierr = NF90_ENDDEF(ncid)

 !ierr = NF90_CREATE(outfilename, NF90_CLOBBER, ncid)
 ierr = NF90_PUT_VAR(ncid, lon_varid, lon)
 ierr = NF90_PUT_VAR(ncid, lat_varid, lat)
 ierr = NF90_PUT_VAR(ncid, sfco_varid, sfco(:))
 ierr = NF90_PUT_VAR(ncid, sst_varid, sst(:))
 ierr = NF90_PUT_VAR(ncid, asst_varid, asst(:))
 ierr = NF90_PUT_VAR(ncid, ini_varid, ini(:))
 ierr = NF90_PUT_VAR(ncid, seaice_varid, seaice(:))
 ierr = NF90_PUT_VAR(ncid, lsmsk_varid, lsmsk(:))
 ierr = NF90_PUT_VAR(ncid, points_varid, points(:))

 ierr = NF90_CLOSE(ncid)
#endif

 outfilename = './anal.nc'
 CALL CreateFile(out_ncfile, outfilename, ierr)

 CALL AddDimension(out_ncfile, 'ncol', nhoriz, ierr)
 CALL AddDimension(out_ncfile, 'nlev',   nlev, ierr)
 CALL AddDimension(out_ncfile,  'one',      1, ierr)

 CALL AddVariable(out_ncfile, 'lon',    IO_REAL8, (/'ncol'/), ierr)
 CALL AddVariable(out_ncfile, 'lat',    IO_REAL8, (/'ncol'/), ierr)
 CALL AddVariable(out_ncfile, 'sfco',   IO_REAL8, (/'ncol'/), ierr)
 CALL AddVariable(out_ncfile, 'sst',    IO_REAL8, (/'ncol'/), ierr)
 CALL AddVariable(out_ncfile, 'asst',   IO_REAL8, (/'ncol'/), ierr)
 CALL AddVariable(out_ncfile, 'ini',    IO_REAL8, (/'ncol'/), ierr)
 CALL AddVariable(out_ncfile, 'lsmsk',  IO_INT4,  (/'ncol'/), ierr)
 CALL AddVariable(out_ncfile, 'seaice', IO_REAL8, (/'ncol'/), ierr)
 CALL AddVariable(out_ncfile, 'points', IO_REAL8, (/'ncol'/), ierr)

 CALL WriteVariable(out_ncfile, 'lon',    lon,    ierr)
 CALL WriteVariable(out_ncfile, 'lat',    lat,    ierr)
 CALL WriteVariable(out_ncfile, 'sfco',   sfco,   ierr)
 CALL WriteVariable(out_ncfile, 'sst',    sst,    ierr)
 CALL WriteVariable(out_ncfile, 'asst',   asst,   ierr)
 CALL WriteVariable(out_ncfile, 'ini',    ini,    ierr)
 CALL WriteVariable(out_ncfile, 'lsmsk',  seaice, ierr)
 CALL WriteVariable(out_ncfile, 'seaice', lsmsk , ierr)
 CALL WriteVariable(out_ncfile, 'points', points, ierr)


 CALL CloseFile(out_ncfile, ierr)


 DEALLOCATE(lat)
 DEALLOCATE(lon)
 DEALLOCATE(sfco)
 DEALLOCATE(asst)
 DEALLOCATE(ini)
 DEALLOCATE(lsmsk)
 DEALLOCATE(seaice)
 DEALLOCATE(points)

 CONTAINS




 SUBROUTINE RestoreVar1D(nctype, varname, nx, ny, isIni, var)
  IMPLICIT NONE
  TYPE(netcdf_file_t), INTENT(INOUT) :: nctype
  CHARACTER(LEN=*), INTENT(IN)       :: varname
  INTEGER(i4), INTENT(IN)            :: nx, ny
  LOGICAL, INTENT(IN)                :: isIni
  REAL(r8), DIMENSION(nx, ny)        :: var

  ! local
  INTEGER(i4) :: l_ierr

   CALL ReadVariable(nctype, TRIM(varname), var, l_ierr)

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
