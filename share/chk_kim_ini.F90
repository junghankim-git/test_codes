!-------------------------------------------------------------------------------
   module chk_kim_ini
!-------------------------------------------------------------------------------
!
!  abstract : basic parameters and types
!
!  history log :
!    2016-07-26   junghan kim     first written
!
!-------------------------------------------------------------------------------
   use kinds, only: i4, r4, r8, l4
   use netcdf
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   type kim_ini_t
     character(len=512) :: filename
     logical(l4)        :: initial
     integer(i4)        :: ncol, nlev
     real(r4), dimension(:)  , allocatable :: ps
     real(r4), dimension(:,:), allocatable :: u, v, qv, qc, qr, qs, qi
   end type kim_ini_t
!
   public :: kim_ini_t, initialize, finalize, analysis
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine initialize(kim_t, filename)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_ini_t), intent(inout) :: kim_t
   character(len=*)               :: filename
! local variables
   integer(i4) :: fid, did, vid, ndims, nvars, ncol, nlev
   logical(l4) :: initial
   character(len=16), dimension(8) :: varnames
   real(r8), dimension(:,:), allocatable :: buf
!
   print *, '# loading file: '//trim(filename)
   print *, ' '
   kim_t%filename = trim(filename)
!
   call check_nc(nf90_open(trim(filename),nf90_nowrite,fid),'open')
   call check_nc(nf90_inquire(fid,ndims,nvars),'ndims, nvars')
   if (ndims.eq.2) then
     initial = .true.
   else
     initial = .false.
   endif
   kim_t%initial = initial
!
! dimension
   call check_nc(nf90_inq_dimid(fid,'ncol',did),'ncol')
   call check_nc(nf90_inquire_dimension(fid,did,len=ncol),'ncol')
   if (initial) then
     call check_nc(nf90_inq_dimid(fid,'level',did),'nlev')
     call check_nc(nf90_inquire_dimension(fid,did,len=nlev),'nlev')
   else
     call check_nc(nf90_inq_dimid(fid,'nlev',did),'nlev')
     call check_nc(nf90_inquire_dimension(fid,did,len=nlev),'nlev')
   endif
   kim_t%ncol = ncol
   kim_t%nlev = nlev
!
! allocate
   allocate(kim_t%ps(ncol))
   allocate(kim_t%u(ncol,nlev))
   allocate(kim_t%v(ncol,nlev))
   allocate(kim_t%qv(ncol,nlev))
   allocate(kim_t%qc(ncol,nlev))
   allocate(kim_t%qr(ncol,nlev))
   allocate(kim_t%qs(ncol,nlev))
   allocate(kim_t%qi(ncol,nlev))
   kim_t%ps = 0.0
   kim_t%u  = 0.0
   kim_t%v  = 0.0
   kim_t%qv = 0.0
   kim_t%qc = 0.0
   kim_t%qr = 0.0
   kim_t%qs = 0.0
   kim_t%qi = 0.0
!
! variables
   varnames(1) = 'ps'
   varnames(2) = 'u'
   varnames(3) = 'v'
   varnames(4) = 'q'
   varnames(5) = 'qc'
   varnames(6) = 'qr'
   varnames(7) = 'qs'
   varnames(8) = 'qi'
   if (initial) varnames(1) = 'psfc'
!
   call check_nc(nf90_inq_varid(fid,trim(varnames(1)),vid),'inq_var(1)')
   call check_nc(nf90_get_var(fid,vid,kim_t%ps),'get_var(1)')
   if (initial) then
     allocate(buf(ncol,nlev))
     call check_nc(nf90_inq_varid(fid,trim(varnames(2)),vid),'inq_var(2)')
     call check_nc(nf90_get_var(fid,vid,buf),'get_var(2)')
     kim_t%u = buf(:,nlev:1:-1)
     call check_nc(nf90_inq_varid(fid,trim(varnames(3)),vid),'inq_var(3)')
     call check_nc(nf90_get_var(fid,vid,buf),'get_var(3)')
     kim_t%v = buf(:,nlev:1:-1)
     call check_nc(nf90_inq_varid(fid,trim(varnames(4)),vid),'inq_var(4)')
     call check_nc(nf90_get_var(fid,vid,buf),'get_var(4)')
     kim_t%qv = buf(:,nlev:1:-1)
     call check_nc(nf90_inq_varid(fid,trim(varnames(5)),vid),'inq_var(5)')
     call check_nc(nf90_get_var(fid,vid,buf),'get_var(5)')
     kim_t%qc = buf(:,nlev:1:-1)
     call check_nc(nf90_inq_varid(fid,trim(varnames(6)),vid),'inq_var(6)')
     call check_nc(nf90_get_var(fid,vid,buf),'get_var(6)')
     kim_t%qr = buf(:,nlev:1:-1)
     call check_nc(nf90_inq_varid(fid,trim(varnames(7)),vid),'inq_var(7)')
     call check_nc(nf90_get_var(fid,vid,buf),'get_var(7)')
     kim_t%qs = buf(:,nlev:1:-1)
     call check_nc(nf90_inq_varid(fid,trim(varnames(8)),vid),'inq_var(8)')
     call check_nc(nf90_get_var(fid,vid,buf),'get_var(8)')
     kim_t%qi = buf(:,nlev:1:-1)
     deallocate(buf)
   else
     call check_nc(nf90_inq_varid(fid,trim(varnames(2)),vid),'inq_var(2)')
     call check_nc(nf90_get_var(fid,vid,kim_t%u),'get_var(2)')
     call check_nc(nf90_inq_varid(fid,trim(varnames(3)),vid),'inq_var(3)')
     call check_nc(nf90_get_var(fid,vid,kim_t%v),'get_var(3)')
     call check_nc(nf90_inq_varid(fid,trim(varnames(4)),vid),'inq_var(4)')
     call check_nc(nf90_get_var(fid,vid,kim_t%qv),'get_var(4)')
     call check_nc(nf90_inq_varid(fid,trim(varnames(5)),vid),'inq_var(5)')
     call check_nc(nf90_get_var(fid,vid,kim_t%qc),'get_var(5)')
     call check_nc(nf90_inq_varid(fid,trim(varnames(6)),vid),'inq_var(6)')
     call check_nc(nf90_get_var(fid,vid,kim_t%qr),'get_var(6)')
     call check_nc(nf90_inq_varid(fid,trim(varnames(7)),vid),'inq_var(7)')
     call check_nc(nf90_get_var(fid,vid,kim_t%qs),'get_var(7)')
     call check_nc(nf90_inq_varid(fid,trim(varnames(8)),vid),'inq_var(8)')
     call check_nc(nf90_get_var(fid,vid,kim_t%qi),'get_var(8)')
   endif
!
   call check_nc(nf90_close(fid),'close')
!
#if 0
! print min/max
   print *, 'min/max (ps) = ',minval(kim_t%ps),maxval(kim_t%ps)
   print *, 'min/max ( u) = ',minval(kim_t%u) ,maxval(kim_t%u)
   print *, 'min/max ( v) = ',minval(kim_t%v) ,maxval(kim_t%v)
   print *, 'min/max (qv) = ',minval(kim_t%qv),maxval(kim_t%qv)
   print *, 'min/max (qc) = ',minval(kim_t%qc),maxval(kim_t%qc)
   print *, 'min/max (qr) = ',minval(kim_t%qr),maxval(kim_t%qr)
   print *, 'min/max (qs) = ',minval(kim_t%qs),maxval(kim_t%qs)
   print *, 'min/max (qi) = ',minval(kim_t%qi),maxval(kim_t%qi)
   print *, ' '
#endif
!
   return
   end subroutine initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine finalize(kim_t)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_ini_t), intent(inout) :: kim_t
! local variables
!
! deallocate
   deallocate(kim_t%ps)
   deallocate(kim_t%u)
   deallocate(kim_t%v)
   deallocate(kim_t%qv)
   deallocate(kim_t%qc)
   deallocate(kim_t%qr)
   deallocate(kim_t%qs)
   deallocate(kim_t%qi)
!
   return
   end subroutine finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine analysis(kim1_t, kim2_t)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_ini_t), intent(in   ) :: kim1_t, kim2_t
! local variables
   integer(i4) :: ncol, nlev, k
   real(r4)    :: minv, maxv
!
   ncol = kim1_t%ncol
   nlev = kim1_t%nlev
   if (kim2_t%ncol.ne.ncol) then
     print *, 'diff. ncol...', ncol, kim2_t%ncol
   endif
   if (kim2_t%nlev.ne.nlev) then
     print *, 'diff. nlev...', nlev, kim2_t%nlev
   endif
!
   call check_var_1d('ps',kim1_t%ps,kim2_t%ps)
   call check_var_2d('u',kim1_t%u,kim2_t%u)
   call check_var_2d('v',kim1_t%v,kim2_t%v)
   call check_var_2d('qv',kim1_t%qv,kim2_t%qv)
   call check_var_2d('qc',kim1_t%qc,kim2_t%qc)
   call check_var_2d('qr',kim1_t%qr,kim2_t%qr)
   call check_var_2d('qs',kim1_t%qs,kim2_t%qs)
   call check_var_2d('qi',kim1_t%qi,kim2_t%qi)
!
#if 0
   print *, 'start u'
   do k=1,nlev
     print *, kim1_t%u(1,k), kim2_t%u(1,k)
   enddo
   print *, 'end'
   print *, 'start v'
   do k=1,nlev
     print *, kim1_t%v(1,k), kim2_t%v(1,k)
   enddo
   print *, 'end'
#endif
!
   return
   end subroutine analysis
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_var_1d(varname, var1, var2)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),       intent(in   ) :: varname
   real(r4), dimension(:), intent(in   ) :: var1, var2
! local variables
   real(r4)    :: minv1, maxv1, minv2, maxv2
   logical(l4) :: diff
   character(len=128) :: string
!
   minv1 = minval(var1)
   maxv1 = maxval(var1)
   minv2 = minval(var2)
   maxv2 = maxval(var2)
   diff  = .false.
   if (minv1.ne.minv2.or.maxv1.ne.maxv2) then
     diff = .true.
   else
     diff = diff_var_1d(var1,var2)
   endif
!
   write(string,'(a8,a)') trim(varname),' :'
   if (diff) then
     print *,trim(string),' diff min/max = ',minv1,minv2,maxv1,maxv2
   else
     print *,trim(string),' same min/max = ',minv1,maxv1
   endif
!
   return
   end subroutine check_var_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_var_2d(varname, var1, var2)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),         intent(in   ) :: varname
   real(r4), dimension(:,:), intent(in   ) :: var1, var2
! local variables
   real(r4)    :: minv1, maxv1, minv2, maxv2
   logical(l4) :: diff
   character(len=128) :: string
!
   minv1 = minval(var1)
   maxv1 = maxval(var1)
   minv2 = minval(var2)
   maxv2 = maxval(var2)
   diff  = .false.
   if (minv1.ne.minv2.or.maxv1.ne.maxv2) then
     diff = .true.
   else
     diff = diff_var_2d(var1,var2)
   endif
!
   write(string,'(a8,a)') trim(varname),' :'
   if (diff) then
     print *,trim(string),' diff min/max = ',minv1,minv2,maxv1,maxv2
   else
     print *,trim(string),' same min/max = ',minv1,maxv1
   endif
!
   return
   end subroutine check_var_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function diff_var_1d(var1, var2) result(diff)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:), intent(in   ) :: var1, var2
   logical(l4)                           :: diff
! local variables
   integer(i4) :: ncol, i
!
   ncol = size(var1)
   diff = .false.
   do i = 1,ncol
     if (var1(i).ne.var2(i)) then
       diff = .true.
       exit
     endif
   enddo
!
   end function diff_var_1d 
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function diff_var_2d(var1, var2) result(diff)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), dimension(:,:), intent(in   ) :: var1, var2
   logical(l4)                             :: diff
! local variables
   integer(i4) :: ncol, nlev, i, j
!
   ncol = size(var1,dim=1)
   nlev = size(var1,dim=2)
   diff = .false.
   do j = 1,nlev
     do i = 1,ncol
       if (var1(i,j).ne.var2(i,j)) then
         diff = .true.
         exit
!         print *, var1(i,j), var2(i,j)
       endif
     enddo
   enddo
!
   end function diff_var_2d 
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_nc(err, message)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4)     , intent(in   ) :: err
   character(len=*), intent(in   ) :: message
! local variables
!
   if (err.ne.nf90_noerr) then
     print *, trim(message), ': ', trim(nf90_strerror(err))
     stop
   endif
!
   return
   end subroutine check_nc
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module chk_kim_ini
!-------------------------------------------------------------------------------
