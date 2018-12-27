MODULE ResHis2Post
 USE NetCDF
 PRIVATE

 INTEGER, PARAMETER :: i4=4, r4=4, r8=8, l4=4
! CHARACTER(LEN=512) :: resfilename
! CHARACTER(LEN=512) :: hisfilename
 CHARACTER(LEN=512) :: outfilename

 !res: ncol, nlev // hyam, hybm, hyai, hybi, u, v, q
 !his: ps, topo, p, pint, T

 ! Variables
 INTEGER(i4) :: ncol, nlev, nilev
 REAL(r8), DIMENSION(:),   ALLOCATABLE :: hyam, hybm, hyai, hybi, levs, ilevs
 REAL(r8), DIMENSION(:),   ALLOCATABLE :: ps, topo
 REAL(r8), DIMENSION(:,:), ALLOCATABLE :: u, v, T, p, pint, q

 ! NetCDF
 INTEGER(i4) :: resid, hisid


 PUBLIC :: Ini, Run, Fin

CONTAINS

 SUBROUTINE Ini(resfile, hisfile, outfile)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: resfile
  CHARACTER(LEN=*), INTENT(IN) :: hisfile
  CHARACTER(LEN=*), INTENT(IN) :: outfile
  ! local
  INTEGER(i4) :: ierr, ncol_tmp, nlev_tmp
  INTEGER(i4) :: ncol_id, nlev_id

!   resfilename = TRIM(resfile)
!   hisfilename = TRIM(hisfile)
   outfilename = TRIM(outfile)

   !=== restart file
   ierr = NF90_OPEN(TRIM(resfile), NF90_NOWRITE, resid)
   IF (ierr .NE. 0) THEN
     PRINT *, 'Can not open file : ', TRIM(resfile)
     STOP
   END IF
   ! ncol
   ierr = NF90_INQ_DIMID(resid, 'ncol', ncol_id)
   IF (ierr .NE. 0) THEN
     PRINT *, 'Can not find dimension (ncol) in file : ', TRIM(resfile)
     STOP
   END IF
   ierr = NF90_INQUIRE_DIMENSION(resid, ncol_id, LEN=ncol)
   ! nlev
   ierr = NF90_INQ_DIMID(resid, 'lev', nlev_id)
   IF (ierr .NE. 0) THEN
     PRINT *, 'Can not find dimension (lev) in file : ', TRIM(resfile)
     STOP
   END IF
   ierr = NF90_INQUIRE_DIMENSION(resid, nlev_id, LEN=nlev)

   !=== history file
   ierr = NF90_OPEN(TRIM(hisfile), NF90_NOWRITE, hisid)
   IF (ierr .NE. 0) THEN
     PRINT *, 'Can not open file : ', TRIM(hisfile)
     STOP
   END IF
   ! ncol
   ierr = NF90_INQ_DIMID(hisid, 'ncol', ncol_id)
   IF (ierr .NE. 0) THEN
     PRINT *, 'Can not find dimension (ncol) in file : ', TRIM(hisfile)
     STOP
   END IF
   ierr = NF90_INQUIRE_DIMENSION(hisid, ncol_id, LEN=ncol_tmp)
   ! nlev
   ierr = NF90_INQ_DIMID(hisid, 'lev', nlev_id)
   IF (ierr .NE. 0) THEN
     PRINT *, 'Can not find dimension (lev) in file : ', TRIM(hisfile)
     STOP
   END IF
   ierr = NF90_INQUIRE_DIMENSION(hisid, nlev_id, LEN=nlev_tmp)


   !=== Check dimension
   IF (ncol .NE. ncol_tmp) THEN
     PRINT *, 'Check ncol sizes in two files...'
     STOP
   END IF
   IF (nlev .NE. nlev_tmp) THEN
     PRINT *, 'Check nlev sizes in two files...'
     STOP
   END IF


   !=== Initialize Variables
   nilev = nlev+1
   ALLOCATE(hyam(nlev))
   ALLOCATE(hybm(nlev))
   ALLOCATE(hyai(nilev))
   ALLOCATE(hybi(nilev))
   ALLOCATE(levs(nlev))
   ALLOCATE(ilevs(nilev))
   ALLOCATE(ps(ncol))
   ALLOCATE(topo(ncol))
   ALLOCATE(u(ncol,nlev))
   ALLOCATE(v(ncol,nlev))
   ALLOCATE(T(ncol,nlev))
   ALLOCATE(p(ncol,nlev))
   ALLOCATE(q(ncol,nlev))
   ALLOCATE(pint(ncol,nilev))

 END SUBROUTINE Ini



 
 SUBROUTINE Run()
  IMPLICIT NONE
  INTEGER(i4) :: ierr
               !  res,      res,      res,    res,    res,      res,    his,    his,     res, res,  his,  his,   res,  his
  INTEGER(i4) :: hyam_id, hybm_id, hyai_id, hybi_id, levs_id, ilevs_id, ps_id, topo_id, u_id, v_id, T_id, p_id, q_id, pint_id
 !res: ncol, nlev // hyam, hybm, hyai, hybi, u, v, q
 !his: ps, topo, p, pint, T


   ! restart
   ierr = NF90_INQ_VARID(resid, 'hyam', hyam_id)
   CALL CheckNC(ierr, '   . inq. hyam in res. file')
   ierr = NF90_GET_VAR(resid, hyam_id, hyam)
   CALL CheckNC(ierr, '   . read hyam in res. file')
   ierr = NF90_INQ_VARID(resid, 'hybm', hybm_id)
   CALL CheckNC(ierr, '   . inq. hybm in res. file')
   ierr = NF90_GET_VAR(resid, hybm_id, hybm)
   CALL CheckNC(ierr, '   . read hybm in res. file')
   ierr = NF90_INQ_VARID(resid, 'hyai', hyai_id)
   CALL CheckNC(ierr, '   . inq. hyai in res. file')
   ierr = NF90_GET_VAR(resid, hyai_id, hyai)
   CALL CheckNC(ierr, '   . read hyai in res. file')
   ierr = NF90_INQ_VARID(resid, 'hybi', hybi_id)
   CALL CheckNC(ierr, '   . inq. hybi in res. file')
   ierr = NF90_GET_VAR(resid, hybi_id, hybi)
   CALL CheckNC(ierr, '   . read hybi in res. file')
   ierr = NF90_INQ_VARID(resid, 'lev', levs_id)
   CALL CheckNC(ierr, '   . inq. levs in res. file')
   ierr = NF90_GET_VAR(resid, levs_id, levs)
   CALL CheckNC(ierr, '   . read levs in res. file')
   ierr = NF90_INQ_VARID(resid, 'ilev', ilevs_id)
   CALL CheckNC(ierr, '   . inq. ilevs in res. file')
   ierr = NF90_GET_VAR(resid, ilevs_id, ilevs)
   CALL CheckNC(ierr, '   . read ilevs in res. file')


   ierr = NF90_INQ_VARID(resid, 'u', u_id)
   CALL CheckNC(ierr, '   . inq. u in res. file')
   ierr = NF90_GET_VAR(resid, u_id, u)
   CALL CheckNC(ierr, '   . read u in res. file')
   ierr = NF90_INQ_VARID(resid, 'v', v_id)
   CALL CheckNC(ierr, '   . inq. v in res. file')
   ierr = NF90_GET_VAR(resid, v_id, v)
   CALL CheckNC(ierr, '   . read v in res. file')
   ierr = NF90_INQ_VARID(resid, 'Q_001', q_id)
   CALL CheckNC(ierr, '   . inq. q in res. file')
   ierr = NF90_GET_VAR(resid, q_id, q)
   CALL CheckNC(ierr, '   . read q in res. file')

   ! history
   ierr = NF90_INQ_VARID(hisid, 'ps', ps_id)
   CALL CheckNC(ierr, '   . inq. ps in his. file')
   ierr = NF90_GET_VAR(hisid, ps_id, ps)
   CALL CheckNC(ierr, '   . read ps in his. file')
   ierr = NF90_INQ_VARID(hisid, 'topo', topo_id)
   CALL CheckNC(ierr, '   . inq. topo in his. file')
   ierr = NF90_GET_VAR(hisid, topo_id, topo)
   CALL CheckNC(ierr, '   . read topo in his. file')
   ierr = NF90_INQ_VARID(hisid, 'T', T_id)
   CALL CheckNC(ierr, '   . inq. T in his. file')
   ierr = NF90_GET_VAR(hisid, T_id, T)
   CALL CheckNC(ierr, '   . read T in his. file')
   ierr = NF90_INQ_VARID(hisid, 'p', p_id)
   CALL CheckNC(ierr, '   . inq. p in his. file')
   ierr = NF90_GET_VAR(hisid, p_id, p)
   CALL CheckNC(ierr, '   . read p in his. file')
   ierr = NF90_INQ_VARID(hisid, 'pint', pint_id)
   CALL CheckNC(ierr, '   . inq. pint in his. file')
   ierr = NF90_GET_VAR(hisid, pint_id, pint)
   CALL CheckNC(ierr, '   . read pint in his. file')

   ! Close file
   ierr = NF90_CLOSE(resid)
   ierr = NF90_CLOSE(hisid)


   ! Write Result
   CALL WriteResult()

 END SUBROUTINE Run



 
 SUBROUTINE Fin()
  IMPLICIT NONE

   DEALLOCATE(hyam)
   DEALLOCATE(hybm)
   DEALLOCATE(hyai)
   DEALLOCATE(hybi)
   DEALLOCATE(levs)
   DEALLOCATE(ilevs)
   DEALLOCATE(ps)
   DEALLOCATE(topo)
   DEALLOCATE(u)
   DEALLOCATE(v)
   DEALLOCATE(T)
   DEALLOCATE(p)
   DEALLOCATE(q)
   DEALLOCATE(pint)

 END SUBROUTINE Fin



 SUBROUTINE WriteResult()
  IMPLICIT NONE
  INTEGER(i4) :: ierr
  INTEGER(i4) :: wmod, file_id
  INTEGER(i4) :: ncol_id, nlev_id, nilev_id
  INTEGER(i4) :: hyam_id, hybm_id, levs_id, hyai_id, hybi_id, ilevs_id
  INTEGER(i4) :: ps_id, topo_id, u_id, v_id, T_id, p_id, q_id, pint_id
  
   ! Open file
   wmod = IOR(NF90_CLOBBER, NF90_64BIT_OFFSET)
   ierr = NF90_CREATE(TRIM(outfilename), wmod, file_id)

   ! Define dimensions
   ierr = NF90_DEF_DIM(file_id, 'ncol',  ncol,  ncol_id)
   ierr = NF90_DEF_DIM(file_id, 'lev',   nlev,  nlev_id)
   ierr = NF90_DEF_DIM(file_id, 'ilev',  nilev, nilev_id)

   ! Define variables
   ierr = NF90_DEF_VAR(file_id, 'hyam', NF90_REAL8,             nlev_id,  hyam_id)
   ierr = NF90_DEF_VAR(file_id, 'hybm', NF90_REAL8,             nlev_id,  hybm_id)
   ierr = NF90_DEF_VAR(file_id, 'hyai', NF90_REAL8,            nilev_id,  hyai_id)
   ierr = NF90_DEF_VAR(file_id, 'hybi', NF90_REAL8,            nilev_id,  hybi_id)
   ierr = NF90_DEF_VAR(file_id, 'lev',  NF90_REAL8,             nlev_id,  levs_id)
   ierr = NF90_DEF_VAR(file_id, 'ilev', NF90_REAL8,            nilev_id, ilevs_id)
   ierr = NF90_DEF_VAR(file_id, 'ps',   NF90_REAL8,             ncol_id,    ps_id)
   ierr = NF90_DEF_VAR(file_id, 'topo', NF90_REAL8,             ncol_id,  topo_id)
   ierr = NF90_DEF_VAR(file_id, 'u',    NF90_REAL8, (/ncol_id,nlev_id/),     u_id)
   ierr = NF90_DEF_VAR(file_id, 'v',    NF90_REAL8, (/ncol_id,nlev_id/),     v_id)
   ierr = NF90_DEF_VAR(file_id, 'T',    NF90_REAL8, (/ncol_id,nlev_id/),     T_id)
   ierr = NF90_DEF_VAR(file_id, 'p',    NF90_REAL8, (/ncol_id,nlev_id/),     p_id)
   ierr = NF90_DEF_VAR(file_id, 'q',    NF90_REAL8, (/ncol_id,nlev_id/),     q_id)
   ierr = NF90_DEF_VAR(file_id, 'pint', NF90_REAL8, (/ncol_id,nilev_id/), pint_id)

   ! End Define
   ierr = NF90_ENDDEF(file_id)

   ! Write Variable
   ierr = NF90_PUT_VAR(file_id, hyam_id,  hyam)
   ierr = NF90_PUT_VAR(file_id, hybm_id,  hybm)
   ierr = NF90_PUT_VAR(file_id, hyai_id,  hyai)
   ierr = NF90_PUT_VAR(file_id, hybi_id,  hybi)
   ierr = NF90_PUT_VAR(file_id, levs_id,  levs)
   ierr = NF90_PUT_VAR(file_id,ilevs_id, ilevs)

   ierr = NF90_PUT_VAR(file_id,   ps_id,  ps)
   ierr = NF90_PUT_VAR(file_id, topo_id,  topo)
   ierr = NF90_PUT_VAR(file_id,    u_id,  u)
   ierr = NF90_PUT_VAR(file_id,    v_id,  v)
   ierr = NF90_PUT_VAR(file_id,    T_id,  T)
   ierr = NF90_PUT_VAR(file_id,    p_id,  p)
   ierr = NF90_PUT_VAR(file_id,    q_id,  q)
   ierr = NF90_PUT_VAR(file_id, pint_id,  pint)

   ! Close file
   ierr = NF90_CLOSE(file_id)

 END SUBROUTINE WriteResult



 SUBROUTINE CheckNC(err, message)
  IMPLICIT NONE
  INTEGER(i4),      INTENT(IN) :: err
  CHARACTER(LEN=*), INTENT(IN) :: message

   IF (err .EQ. 0) THEN
     PRINT *, TRIM(message)
   ELSE
     PRINT *, TRIM(message)
     STOP
   END IF

 END SUBROUTINE CheckNC


END MODULE ResHis2Post
