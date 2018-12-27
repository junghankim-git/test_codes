MODULE CheckRestart
 USE NetCDF
 IMPLICIT NONE

  INTEGER, PARAMETER :: i4 = 4, r4 = 4, r8 = 8, l4 = 4

  TYPE diminfo_t
    INTEGER(i4), DIMENSION(:), ALLOCATABLE       :: sizes
  END TYPE diminfo_t

  
  INTEGER(i4), PARAMETER                     :: nfiles = 2
  CHARACTER(LEN=512), DIMENSION(nfiles)      :: filenames

  !INTEGER(i4), PARAMETER                     :: nvars = 117
  INTEGER(i4), PARAMETER                     :: nvars = 51
  CHARACTER(LEN=64), DIMENSION(nvars,nfiles) :: varnames
  INTEGER(i4), DIMENSION(nvars)              :: ndims
!  INTEGER(i4), DIMENSION(nvars)              :: diminfo
  TYPE(diminfo_t), DIMENSION(nvars)          :: diminfo

  ! dimensions
  INTEGER(i4), PARAMETER :: ncol      = 48602
  INTEGER(i4), PARAMETER :: nlev      = 50
  INTEGER(i4), PARAMETER :: nilev     = 51
  INTEGER(i4), PARAMETER :: ndate     = 4
  INTEGER(i4), PARAMETER :: natype    = 5
  INTEGER(i4), PARAMETER :: mvar1d    = 11
  INTEGER(i4), PARAMETER :: laz       = 16
  INTEGER(i4), PARAMETER :: loz_gmao  = 58
  INTEGER(i4), PARAMETER :: numsfcsp9 = 69

  INTERFACE IsSame
    MODULE PROCEDURE IsSame_1D
  END INTERFACE IsSame

CONTAINS

 SUBROUTINE Ini(filename1, filename2)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: filename1
  CHARACTER(LEN=*), INTENT(IN) :: filename2
  ! local
  INTEGER(i4) :: ivar

   filenames(1) = TRIM(filename1)
   filenames(2) = TRIM(filename2)

   varnames(001,2) = 'lats';           ndims(001) = 1
   varnames(002,2) = 'lons';           ndims(002) = 1
   varnames(003,2) = 'levs';           ndims(003) = 1
   varnames(004,2) = 'hyam';           ndims(004) = 1
   varnames(005,2) = 'hybm';           ndims(005) = 1
   varnames(006,2) = 'ilevs';          ndims(006) = 1
   varnames(007,2) = 'hyai';           ndims(007) = 1
   varnames(008,2) = 'hybi';           ndims(008) = 1
   varnames(009,2) = 'pdst';           ndims(009) = 1
   varnames(010,2) = 'ps';             ndims(010) = 1
   varnames(011,2) = 'ht';             ndims(011) = 1
   varnames(012,2) = 'dht_1';          ndims(012) = 1
   varnames(013,2) = 'dht_2';          ndims(013) = 1
   varnames(014,2) = 'pds_2';          ndims(014) = 1
   varnames(015,2) = 'pdsb';           ndims(015) = 1
   varnames(016,2) = 'u';              ndims(016) = 2
   varnames(017,2) = 'v';              ndims(017) = 2
   varnames(018,2) = 'theta';          ndims(018) = 2
   varnames(019,2) = 'omega';          ndims(019) = 2
   varnames(020,2) = 'mub';            ndims(020) = 2
   varnames(021,2) = 'pb';             ndims(021) = 2
   varnames(022,2) = 't_init';         ndims(022) = 2
   varnames(023,2) = 'alb';            ndims(023) = 2
   varnames(024,2) = 'mu_2';           ndims(024) = 2
   varnames(025,2) = 'mu_2_w';         ndims(025) = 2
   varnames(026,2) = 'p';              ndims(026) = 2
   varnames(027,2) = 'alt';            ndims(027) = 2
   varnames(028,2) = 'al';             ndims(028) = 2
   varnames(029,2) = 'p8w';            ndims(029) = 2
   varnames(030,2) = 'ph0';            ndims(030) = 2
   varnames(031,2) = 'mut_w';          ndims(031) = 2
   varnames(032,2) = 'rv_tendf1';      ndims(032) = 2
   varnames(033,2) = 'rv_tendf2';      ndims(033) = 2
   varnames(034,2) = 'rw_tendf';       ndims(034) = 2
   varnames(035,2) = 't_tendf';        ndims(035) = 2
   varnames(036,2) = 'T';              ndims(036) = 2
   varnames(037,2) = 'w_2';            ndims(037) = 2
   varnames(038,2) = 'mub_w';          ndims(038) = 2
   varnames(039,2) = 'phb';            ndims(039) = 2
   varnames(040,2) = 'p_phy';          ndims(040) = 2
   varnames(041,2) = 'p_hyd_w';        ndims(041) = 2
   varnames(042,2) = 'ph_2';           ndims(042) = 2
   varnames(043,2) = 'mut';            ndims(043) = 2
   varnames(044,2) = 'q';              ndims(044) = 2
   varnames(045,2) = 'qc';             ndims(045) = 2
   varnames(046,2) = 'qi';             ndims(046) = 2
   varnames(047,2) = 'qr';             ndims(047) = 2
   varnames(048,2) = 'qs';             ndims(048) = 2
   varnames(049,2) = 'oz';             ndims(049) = 2
   varnames(050,2) = 'cld';            ndims(050) = 2
   varnames(051,2) = 'cld_f';          ndims(051) = 2


!   varnames(061,2) = 'tsea_oml_init';  ndims(061) = 1
!   varnames(062,2) = 'tsea_oml_init';  ndims(062) = 1
!   varnames(063,2) = 'tsea_oml_init';  ndims(063) = 1
!   varnames(064,2) = 'tsea_oml_init';  ndims(064) = 1
!   varnames(065,2) = 'tsea_oml_init';  ndims(065) = 1
!   varnames(066,2) = 'tsea_oml_init';  ndims(066) = 1
!   varnames(067,2) = 'tsea_oml_init';  ndims(067) = 1
!   varnames(068,2) = 'tsea_oml_init';  ndims(068) = 1
!   varnames(069,2) = 'tsea_oml_init';  ndims(069) = 1
!   varnames(070,2) = 'tsea_oml_init';  ndims(070) = 1

   DO ivar = 1 ,nvars
     IF      (ndims(ivar) .EQ. 1) THEN
       ALLOCATE(diminfo(ivar)%sizes(1))
     ELSE IF (ndims(ivar) .EQ. 2) THEN
       ALLOCATE(diminfo(ivar)%sizes(2))
     ELSE IF (ndims(ivar) .EQ. 3) THEN
       ALLOCATE(diminfo(ivar)%sizes(3))
     ELSE IF (ndims(ivar) .EQ. 4) THEN
       ALLOCATE(diminfo(ivar)%sizes(4))
     END IF
   END DO
   diminfo(001)%sizes(:) = (/ncol/)
   diminfo(002)%sizes(:) = (/ncol/)
   diminfo(003)%sizes(:) = (/nlev/)
   diminfo(004)%sizes(:) = (/nlev/)
   diminfo(005)%sizes(:) = (/nlev/)
   diminfo(006)%sizes(:) = (/nilev/)
   diminfo(007)%sizes(:) = (/nilev/)
   diminfo(008)%sizes(:) = (/nilev/)
   diminfo(009)%sizes(:) = (/ncol/)
   diminfo(010)%sizes(:) = (/ncol/)
   diminfo(011)%sizes(:) = (/ncol/)
   diminfo(012)%sizes(:) = (/ncol/)
   diminfo(013)%sizes(:) = (/ncol/)
   diminfo(014)%sizes(:) = (/ncol/)
   diminfo(015)%sizes(:) = (/ncol/)
   diminfo(016)%sizes(:) = (/ncol,nlev/)
   diminfo(017)%sizes(:) = (/ncol,nlev/)
   diminfo(018)%sizes(:) = (/ncol,nlev/)
   diminfo(019)%sizes(:) = (/ncol,nlev/)
   diminfo(020)%sizes(:) = (/ncol,nlev/)
   diminfo(021)%sizes(:) = (/ncol,nlev/)
   diminfo(022)%sizes(:) = (/ncol,nlev/)
   diminfo(023)%sizes(:) = (/ncol,nlev/)
   diminfo(024)%sizes(:) = (/ncol,nlev/)
   diminfo(025)%sizes(:) = (/ncol,nlev/)
   diminfo(026)%sizes(:) = (/ncol,nlev/)
   diminfo(027)%sizes(:) = (/ncol,nlev/)
   diminfo(028)%sizes(:) = (/ncol,nlev/)
   diminfo(029)%sizes(:) = (/ncol,nlev/)
   diminfo(030)%sizes(:) = (/ncol,nlev/)
   diminfo(031)%sizes(:) = (/ncol,nlev/)
   diminfo(032)%sizes(:) = (/ncol,nlev/)
   diminfo(033)%sizes(:) = (/ncol,nlev/)
   diminfo(034)%sizes(:) = (/ncol,nlev/)
   diminfo(035)%sizes(:) = (/ncol,nlev/)
   diminfo(036)%sizes(:) = (/ncol,nlev/)
   diminfo(037)%sizes(:) = (/ncol,nilev/)
   diminfo(038)%sizes(:) = (/ncol,nilev/)
   diminfo(039)%sizes(:) = (/ncol,nilev/)
   diminfo(040)%sizes(:) = (/ncol,nilev/)
   diminfo(041)%sizes(:) = (/ncol,nilev/)
   diminfo(042)%sizes(:) = (/ncol,nilev/)
   diminfo(043)%sizes(:) = (/ncol,nilev/)
   diminfo(044)%sizes(:) = (/ncol,nlev/)
   diminfo(045)%sizes(:) = (/ncol,nlev/)
   diminfo(046)%sizes(:) = (/ncol,nlev/)
   diminfo(047)%sizes(:) = (/ncol,nlev/)
   diminfo(048)%sizes(:) = (/ncol,nlev/)
   diminfo(049)%sizes(:) = (/ncol,nlev/)
   diminfo(050)%sizes(:) = (/ncol,nlev/)
   diminfo(051)%sizes(:) = (/ncol,nlev/)

   ! change varname
   varnames(:,1)  = varnames(:,2)

   varnames(001,1) = 'lat'
   varnames(002,1) = 'lon'
   varnames(003,1) = 'lev'
   varnames(006,1) = 'ilev'
   varnames(011,1) = 'z'
   varnames(044,1) = 'Q_001'
   varnames(045,1) = 'Q_002'
   varnames(046,1) = 'Q_003'
   varnames(047,1) = 'Q_004'
   varnames(048,1) = 'Q_005'
   varnames(049,1) = 'Q_006'
   varnames(050,1) = 'Q_007'
   varnames(051,1) = 'cld_forcing'

 END SUBROUTINE Ini





 SUBROUTINE Fin()
  IMPLICIT NONE
  ! local
  INTEGER(i4) :: ivar

   DO ivar = 1 ,nvars
     DEALLOCATE(diminfo(ivar)%sizes)
   END DO

 END SUBROUTINE Fin





 SUBROUTINE Run()
  IMPLICIT NONE
  ! local
  INTEGER(i4) :: fid(2), vid(2)
  INTEGER(i4) :: ifile, ivar, ierr, verr(nfiles)
  LOGICAL(l4) :: same
  REAL(r8), DIMENSION(:,:),       ALLOCATABLE :: buf_1d
  REAL(r8), DIMENSION(:,:,:),     ALLOCATABLE :: buf_2d
  REAL(r8), DIMENSION(:,:,:,:),   ALLOCATABLE :: buf_3d
  REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: buf_4d

   PRINT *, '# Open 2 files'
   DO ifile = 1, nfiles
     ierr = NF90_OPEN(TRIM(filenames(ifile)), NF90_NOWRITE, fid(ifile))
     PRINT *, '  - stat: ', ierr
   END DO
   PRINT *, ' '

   DO ivar = 1 ,nvars

     PRINT *, '# Read Variable: ', TRIM(varnames(ivar,2))
     vid = -10000000
     IF      (ndims(ivar) .EQ. 1) THEN

       ALLOCATE(buf_1d(diminfo(ivar)%sizes(1), nfiles))

       ! start
       buf_1d(:,1) = -099.0_r8
       buf_1d(:,2) = -999.0_r8
       DO ifile = 1, nfiles
       verr(ifile) = NF90_INQ_VARID(fid(ifile), TRIM(varnames(ivar,ifile)), vid(ifile))
       ierr = NF90_GET_VAR(fid(ifile), vid(ifile), buf_1d(:,ifile))
       END DO
       same = IsSame(buf_1d(:,1), buf_1d(:,2))
       ! stop

       DEALLOCATE(buf_1d)

     ELSE IF (ndims(ivar) .EQ. 2) THEN

       ALLOCATE(buf_2d(diminfo(ivar)%sizes(1), diminfo(ivar)%sizes(2), nfiles))

       ! start
       buf_2d(:,:,1) = -099.0_r8
       buf_2d(:,:,2) = -999.0_r8
       DO ifile = 1, nfiles
       verr(ifile) = NF90_INQ_VARID(fid(ifile), TRIM(varnames(ivar,ifile)), vid(ifile))
       !IF (ierr .NE. 0) PRINT *, 'check var:', TRIM(varnames(ivar,ifile)), ifile
       ierr = NF90_GET_VAR(fid(ifile), vid(ifile), buf_2d(:,:,ifile))
       END DO
       same = IsSame(buf_2d(:,1,1), buf_2d(:,1,2))
       ! stop

       DEALLOCATE(buf_2d)

     ELSE IF (ndims(ivar) .EQ. 3) THEN

       ALLOCATE(buf_3d(diminfo(ivar)%sizes(1), diminfo(ivar)%sizes(2), diminfo(ivar)%sizes(3), nfiles))

       ! start
       buf_3d(:,:,:,1) = -099.0_r8
       buf_3d(:,:,:,2) = -999.0_r8
       DO ifile = 1, nfiles
       verr(ifile) = NF90_INQ_VARID(fid(ifile), TRIM(varnames(ivar,ifile)), vid(ifile))
       !IF (ierr .NE. 0) PRINT *, 'check var:', TRIM(varnames(ivar,ifile)), ifile
       ierr = NF90_GET_VAR(fid(ifile), vid(ifile), buf_3d(:,:,:,ifile))
       END DO
       same = IsSame(buf_3d(:,1,1,1), buf_3d(:,1,1,2))
       ! stop

       DEALLOCATE(buf_3d)

     ELSE IF (ndims(ivar) .EQ. 4) THEN

       ALLOCATE(buf_4d(diminfo(ivar)%sizes(1), diminfo(ivar)%sizes(2), diminfo(ivar)%sizes(3), diminfo(ivar)%sizes(4), nfiles))

       ! start
       buf_4d(:,:,:,:,1) = -099.0_r8
       buf_4d(:,:,:,:,2) = -999.0_r8
       DO ifile = 1, nfiles
       verr(ifile) = NF90_INQ_VARID(fid(ifile), TRIM(varnames(ivar,ifile)), vid(ifile))
       !IF (ierr .NE. 0) PRINT *, 'check var:', TRIM(varnames(ivar,ifile)), ifile
       ierr = NF90_GET_VAR(fid(ifile), vid(ifile), buf_4d(:,:,:,:,ifile))
       END DO
       same = IsSame(buf_4d(:,1,1,1,1), buf_4d(:,1,1,1,2))
       ! stop

       DEALLOCATE(buf_4d)

     END IF

     IF (verr(1) .NE. 0 .OR. verr(2) .NE. 0) THEN
       PRINT *, '   => check varname', verr(1), verr(2)
     ELSE
     IF (same) THEN
       PRINT *, '   => is same.'
     ELSE
       PRINT *, '   => is diff.'
     END IF
     END IF

   END DO ! DO ivar = 


   PRINT *, '# Close 2 files'
   DO ifile = 1, nfiles
     ierr = NF90_CLOSE(fid(ifile))
     PRINT *, '  - stat: ', ierr
   END DO
   PRINT *, ' '

 END SUBROUTINE Run



 FUNCTION IsSame_1D(var1, var2) RESULT(issame)
  IMPLICIT NONE
  REAL(r8), DIMENSION(:), INTENT(IN) :: var1, var2
  LOGICAL(l4)                        :: issame
  ! local
  INTEGER(i4) :: dim1, i1

   issame = .TRUE.
   dim1   = SIZE(var1, DIM=1)

   DO i1 = 1, dim1
     IF (var1(i1) .NE. var2(i1)) THEN
       issame = .FALSE.
       EXIT
     END IF
   END DO

 END FUNCTION IsSame_1D

END MODULE CheckRestart
