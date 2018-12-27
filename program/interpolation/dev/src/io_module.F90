!-------------------------------------------------------------------------------
   module io_module
!
   use kinds, only: i4, l4, r8
   use netcdf
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4) :: nfun, norder, n, m
   integer(i4) :: fid
   integer(i4) :: nfid, noid, nid, mid
   integer(i4) :: sxid, syid, txid
   integer(i4) :: ayid, ay1id, ay2id, ayiid
   integer(i4) :: yid, y1id, y2id, yiid
!
   public :: create_ncfile
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine create_ncfile(filename,testcase,isgll,nf,no,nn,mm,sx,sy,tx,ty_a,typ_a,typp_a,tint_a, &
                                                           ty,typ,typp,tint)
!-------------------------------------------------------------------------------
   implicit none
   character(len=*),              intent(in   ) :: filename
   integer(i4),                   intent(in   ) :: testcase, nf, no, nn, mm
   logical(l4),                   intent(in   ) :: isgll
   real(r8)   , dimension(nn)   , intent(in   ) :: sx, sy
   real(r8)   , dimension(mm)   , intent(in   ) :: tx, ty_a, typ_a, typp_a, tint_a
   real(r8)   , dimension(mm,nf), intent(in   ) :: ty, typ, typp, tint
! local
   integer(i4) :: err, gll
!
   nfun = nf; norder = no; n = nn; m = mm
   if (isgll) then
     gll = 1
   else
     gll = 0
   endif
!
   err = nf90_create(trim(filename),ior(nf90_clobber,nf90_netcdf4),fid)
   err = nf90_put_att(fid,nf90_global,'testcase',testcase)
   err = nf90_put_att(fid,nf90_global,'isgll',gll)
   err = nf90_def_dim(fid,'nfun',nf,nfid)
   err = nf90_def_dim(fid,'norder',no,noid)
   err = nf90_def_dim(fid,'n',n,nid)
   err = nf90_def_dim(fid,'m',m,mid)
   err = nf90_def_var(fid,'sx',nf90_double,nid,sxid)
   err = nf90_def_var(fid,'sy',nf90_double,nid,syid)
   err = nf90_def_var(fid,'tx',nf90_double,mid,txid)
   err = nf90_def_var(fid,'ty_a',nf90_double,mid,ayid)
   err = nf90_def_var(fid,'typ_a',nf90_double,mid,ay1id)
   err = nf90_def_var(fid,'typp_a',nf90_double,mid,ay2id)
   err = nf90_def_var(fid,'tint_a',nf90_double,mid,ayiid)
   err = nf90_def_var(fid,'ty',nf90_double,(/mid,nfid/),yid)
   err = nf90_def_var(fid,'typ',nf90_double,(/mid,nfid/),y1id)
   err = nf90_def_var(fid,'typp',nf90_double,(/mid,nfid/),y2id)
   err = nf90_def_var(fid,'tint',nf90_double,(/mid,nfid/),yiid)
   err = nf90_enddef(fid)
   err = nf90_put_var(fid,sxid,sx)
   err = nf90_put_var(fid,syid,sy)
   err = nf90_put_var(fid,txid,tx)
   err = nf90_put_var(fid,ayid,ty_a)
   err = nf90_put_var(fid,ay1id,typ_a)
   err = nf90_put_var(fid,ay2id,typp_a)
   err = nf90_put_var(fid,ayiid,tint_a)
   err = nf90_put_var(fid,yid,ty)
   err = nf90_put_var(fid,y1id,typ)
   err = nf90_put_var(fid,y2id,typp)
   err = nf90_put_var(fid,yiid,tint)
   err = nf90_close(fid)
!
   return
   end subroutine create_ncfile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module io_module
!-------------------------------------------------------------------------------
