!-------------------------------------------------------------------------------
   module check_var_in_netcdf
!-------------------------------------------------------------------------------
!
!  abstract : converter between binary and decimal
!
!  histroy log :
!    2018-09-08   junghan kim   initial setup
!
!  variable :
!
!-------------------------------------------------------------------------------
   use kinds, only: l4, i4, r4, r8=>r8d
   use netcdf
   private
!
!
   public :: do_check
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine do_check(fname1, fname2, vname1, vname2)
   implicit none
   character(len=*), intent(in   ) :: fname1, fname2, vname1, vname2
! local variables
   integer(i4) :: i, k, ndiff
   integer(i4) :: err, fid, did, vid
   integer(i4) :: ncol, nlev
   real(r8) :: verr, vsum
   real(r8), dimension(:,:), allocatable :: v1, v2, diff
!
!  file 1
   err = nf90_open(trim(fname1),nf90_nowrite,fid)
   call nc_check(err)
   ! dim
   err = nf90_inq_dimid(fid,'ncol',did)
   call nc_check(err)
   err = nf90_inquire_dimension(fid,did,len=ncol)
   call nc_check(err)
   err = nf90_inq_dimid(fid,'nlev',did)
   call nc_check(err)
   err = nf90_inquire_dimension(fid,did,len=nlev)
   call nc_check(err)
   ! allocate
   allocate(v1(ncol,nlev))
   allocate(v2(ncol,nlev))
   allocate(diff(ncol,nlev))
   ! var
   err = nf90_inq_varid(fid,trim(vname1),vid)
   call nc_check(err)
   err = nf90_get_var(fid,vid,v1)
   call nc_check(err)
   err = nf90_close(fid)
   call nc_check(err)
!
!  file 2
   err = nf90_open(trim(fname2),nf90_nowrite,fid)
   call nc_check(err)
   err = nf90_inq_varid(fid,trim(vname2),vid)
   call nc_check(err)
   err = nf90_get_var(fid,vid,v2)
   call nc_check(err)
   err = nf90_close(fid)
   call nc_check(err)
!
! algorithm
   diff = abs(v2-v1)
   do k = 1,nlev
     print '(i2,x,f7.4,x,f7.4,x,f7.4)',k,minval(diff(:,k)),maxval(diff(:,k)),sum(diff(:,k))/real(ncol,r8)
   enddo
#if 0
   ndiff = 0
   do k = 1,nlev
   !do k = 1,1
     do i = 1,ncol
       if (v1(i,k).ne.v2(i,k)) then
         ndiff = ndiff+1
         print *,i,k,v1(i,k)
         print *,i,k,v2(i,k)
         print *,' '
       endif
     enddo
   enddo
   print*,ncol,ncol*nlev,ndiff
   vsum = 0.0_r8
   verr = 0.0_r8
   do k=1,nlev
   do i=1,ncol
     call ddpdd(vsum,verr,v1(i,k))
   enddo
   enddo
   print*,minval(v1),maxval(v1),vsum
#endif
!
   deallocate(v1)
   deallocate(v2)
   deallocate(diff)
!
   return
   end subroutine do_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_check(ncstat_in)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: ncstat_in
!
   if (ncstat_in.ne.nf90_noerr) then
     print *, trim(nf90_strerror(ncstat_in))
     stop
   endif
!
   return
   end subroutine nc_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine scs(add_a, err, add_b)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(inout) :: add_a, err
   real(r8), intent(in   ) :: add_b
! local variables
   real(r8) :: tsum, new_err, tmp
!-------------------------------------------------------------------------------
!
   tmp = err+add_b
   tsum = add_a+tmp
   new_err = tmp-(tsum-add_a)
   add_a = tsum
   err = new_err
!
   return
   end subroutine scs
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ddpdd(add_a, err, add_b)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(inout) :: add_a, err
   real(r8), intent(in   ) :: add_b
! local variables
   real(r8) :: tsum, new_err
!-------------------------------------------------------------------------------
!
   tsum = add_a+add_b
   new_err = (add_b-(tsum-add_a))+(add_a-(tsum-(tsum-add_a)))
   call scs(tsum,err,new_err)
   add_a = tsum
!
   return
   end subroutine ddpdd
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module check_var_in_netcdf
!-------------------------------------------------------------------------------
