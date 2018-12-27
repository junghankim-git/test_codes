!-------------------------------------------------------------------------------
   module nc_control
   use kinds,    only: i4, r4, r8, l4
   use netcdf
!-------------------------------------------------------------------------------
!
   private
!
   type netcdf_t
     integer(i4)        :: in_nc, ou_nc
     integer(i4)        :: ndims, nvars, natts
     integer(i4)        :: compress
     character(len=512) :: infilename, oufilename
   end type netcdf_t
!
   public :: netcdf_t, nc_initialize, nc_finalize
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_initialize(nc, infilename, oufilename, compress)
   implicit none
   type(netcdf_t),   intent(inout) :: nc
   character(len=*), intent(in   ) :: infilename
   character(len=*), intent(in   ) :: oufilename
   integer(i4),      intent(in   ) :: compress
! local variables
   integer(i4), parameter       :: nmax = 100, maxndims=4
   ! dimension
   integer(i4), dimension(nmax) :: dids, dimsize
   ! attributes
   integer(i4), dimension(nmax) :: aids, atype
   integer(i4)                  :: iatt
   real(r4)                     :: ratt
   character(len=64)            :: satt
   ! variables
   integer(i4), dimension(nmax) :: vids, ndims, natts, vtype
   integer(i4), dimension(maxndims,nmax) :: dimids
   integer(i4)                  :: ierr, id, iv, ia
   character(len=256)           :: nam
!
   nc%infilename = trim(infilename)
   nc%oufilename = trim(oufilename)
   nc%compress   = compress
!
   ierr = nf90_create(trim(oufilename),ior(nf90_clobber,nf90_netcdf4),nc%ou_nc)
!
   ierr = nf90_open(trim(infilename),nf90_nowrite,nc%in_nc)
   call nc_check(ierr,'nc_initialize::cannot open...')
   ierr = nf90_inquire(nc%in_nc,nc%ndims,nc%nvars,nc%natts)
   call nc_check(ierr,'nc_initialize::nf90_inquire')
   print *, 'infile : ', trim(infilename)
   print *, ' - ndims ', nc%ndims
   print *, ' - nvars ', nc%nvars
   print *, ' - natts ', nc%natts
!
! dimensions
   do id = 1,nc%ndims
     ierr = nf90_inquire_dimension(nc%in_nc,id,name=nam,len=dimsize(id))
     ierr = nf90_def_dim(nc%ou_nc,trim(nam),dimsize(id),dids(id))
   enddo
! attributes
   do ia = 1,nc%natts
     ierr = nf90_inq_attname(nc%in_nc,nf90_global,ia,nam)
     ierr = nf90_inquire_attribute(nc%in_nc,nf90_global,trim(nam),xtype=atype(ia))
     if     (atype(ia).eq.nf90_int ) then
       ierr = nf90_get_att(nc%in_nc,nf90_global,trim(nam),iatt)
       ierr = nf90_put_att(nc%ou_nc,nf90_global,trim(nam),iatt)
     elseif (atype(ia).eq.nf90_real) then
       ierr = nf90_get_att(nc%in_nc,nf90_global,trim(nam),ratt)
       ierr = nf90_put_att(nc%ou_nc,nf90_global,trim(nam),ratt)
     elseif (atype(ia).eq.nf90_char) then
       ierr = nf90_get_att(nc%in_nc,nf90_global,trim(nam),satt)
       ierr = nf90_put_att(nc%ou_nc,nf90_global,trim(nam),trim(satt))
     else
       print *, 'undefined attributes...', atype(ia)
     endif
   enddo
   ierr = nf90_put_att(nc%ou_nc,nf90_global,'compress',nc%compress)
! variables
   do iv = 1,nc%nvars
     ierr = nf90_inquire_variable(nc%in_nc,iv,name=nam,xtype=vtype(iv),        &
                           ndims=ndims(iv),dimids=dimids(:,iv),natts=natts(iv))
     ierr = nf90_def_var(nc%ou_nc,trim(nam),vtype(iv),dimids(1:ndims(iv),iv),  &
                           vids(iv),deflate_level=nc%compress)
     ! atts
     do ia = 1,natts(iv)
       ierr = nf90_inq_attname(nc%in_nc,vids(iv),ia,nam)
       ierr = nf90_inquire_attribute(nc%in_nc,vids(iv),trim(nam),xtype=atype(ia))
       if     (atype(ia).eq.nf90_int ) then
         ierr = nf90_get_att(nc%in_nc,vids(iv),trim(nam),iatt)
         ierr = nf90_put_att(nc%ou_nc,vids(iv),trim(nam),iatt)
       elseif (atype(ia).eq.nf90_real) then
         ierr = nf90_get_att(nc%in_nc,vids(iv),trim(nam),ratt)
         ierr = nf90_put_att(nc%ou_nc,vids(iv),trim(nam),ratt)
       elseif (atype(ia).eq.nf90_char) then
         ierr = nf90_get_att(nc%in_nc,vids(iv),trim(nam),satt)
         ierr = nf90_put_att(nc%ou_nc,vids(iv),trim(nam),trim(satt))
       else
         print *, 'undefined attributes...', atype(ia)
       endif
     enddo
   enddo
   ierr = nf90_enddef(nc%ou_nc)
!
   end subroutine nc_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_finalize(nc)
   implicit none
   type(netcdf_t),   intent(inout) :: nc
! local variables
   integer(i4)  :: ierr
!
   ierr = nf90_close(nc%ou_nc)
!
   ierr = nf90_close(nc%in_nc)
!
   end subroutine nc_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_check(ierr,message)
   implicit none
   integer(i4),                intent(in   ) :: ierr
   character(len=*), optional, intent(in   ) :: message
! local variables
!
   if (ierr.ne.nf90_noerr) then
     if (present(message)) print *, trim(message)
     stop
   endif
!
   end subroutine nc_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module nc_control
!-------------------------------------------------------------------------------
