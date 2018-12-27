!-------------------------------------------------------------------------------
   module io_module
!
   use kinds, only: i4, l4, r8
   use netcdf
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4) :: n, m, l
   integer(i4) :: fid
   integer(i4) :: nid, np1id, mid, mp1id, lid
   integer(i4) :: sxid, scxid, syid, siid
   integer(i4) :: txid, tcxid, tyid, tiid, tyaid
   integer(i4) :: axid, ayid, tyfid, tyiid
!
   public :: create_ncfile, write_ncfile
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine create_ncfile(filename,use_spline,order,bndry,ismono,testcase,nn,mm,ll)
!-------------------------------------------------------------------------------
   implicit none
   character(len=*), intent(in   ) :: filename
   logical(l4),      intent(in   ) :: use_spline
   integer(i4),      intent(in   ) :: order, bndry
   logical(l4),      intent(in   ) :: ismono
   integer(i4),      intent(in   ) :: testcase, nn, mm, ll
! local
   integer(i4) :: err, spline, mono
!
   n = nn; m = mm; l = ll
   spline = 0; mono = 0
   if (use_spline) spline = 1
   if (ismono) mono = 1
!
   err = nf90_create(trim(filename),ior(nf90_clobber,nf90_netcdf4),fid)
   err = nf90_def_dim(fid,'n',n,nid)
   err = nf90_def_dim(fid,'np1',n+1,np1id)
   err = nf90_def_dim(fid,'m',m,mid)
   err = nf90_def_dim(fid,'mp1',m+1,mp1id)
   err = nf90_def_dim(fid,'l',l,lid)
   err = nf90_def_var(fid,'sx',nf90_double,nid,sxid)
   err = nf90_def_var(fid,'scx',nf90_double,np1id,scxid)
   err = nf90_def_var(fid,'sy',nf90_double,nid,syid)
   err = nf90_def_var(fid,'siy',nf90_double,np1id,siid)
   err = nf90_def_var(fid,'tx',nf90_double,mid,txid)
   err = nf90_def_var(fid,'tcx',nf90_double,mp1id,tcxid)
   err = nf90_def_var(fid,'ty',nf90_double,mid,tyid)
   err = nf90_def_var(fid,'tiy',nf90_double,mp1id,tiid)
   err = nf90_def_var(fid,'ty_ana',nf90_double,mid,tyaid)
   err = nf90_def_var(fid,'x_ana',nf90_double,lid,axid)
   err = nf90_def_var(fid,'y_ana',nf90_double,lid,ayid)
   err = nf90_def_var(fid,'ty_fun',nf90_double,lid,tyfid)
   err = nf90_def_var(fid,'ty_int',nf90_double,lid,tyiid)
   err = nf90_put_att(fid,nf90_global,'use_spline',spline)
   err = nf90_put_att(fid,nf90_global,'order',order)
   err = nf90_put_att(fid,nf90_global,'ismono',mono)
   err = nf90_put_att(fid,nf90_global,'bndry',bndry)
   err = nf90_put_att(fid,nf90_global,'testcase',testcase)
!
   return
   end subroutine create_ncfile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ncfile(sx,scx,sy,siy,tx,tcx,ty,tiy,ty_ana,x_ana,y_ana,ty_fun,ty_int)
!-------------------------------------------------------------------------------
   implicit none
   real(r8), dimension(n)   :: sx, sy
   real(r8), dimension(n+1) :: scx, siy
   real(r8), dimension(m)   :: tx, ty, ty_ana
   real(r8), dimension(n+1) :: tcx, tiy
   real(r8), dimension(l)   :: x_ana, y_ana, ty_fun, ty_int
! local
   integer(i4) :: err
!
   err = nf90_enddef(fid)
   err = nf90_put_var(fid,sxid,sx)
   err = nf90_put_var(fid,scxid,scx)
   err = nf90_put_var(fid,syid,sy)
   err = nf90_put_var(fid,siid,siy)
   err = nf90_put_var(fid,txid,tx)
   err = nf90_put_var(fid,tcxid,tcx)
   err = nf90_put_var(fid,tyid,ty)
   err = nf90_put_var(fid,tiid,tiy)
   err = nf90_put_var(fid,tyaid,ty_ana)
   err = nf90_put_var(fid,axid,x_ana)
   err = nf90_put_var(fid,ayid,y_ana)
   err = nf90_put_var(fid,tyfid,ty_fun)
   err = nf90_put_var(fid,tyiid,ty_int)
   err = nf90_close(fid)
!
   return
   end subroutine write_ncfile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module io_module
!-------------------------------------------------------------------------------
