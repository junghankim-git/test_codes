!-------------------------------------------------------------------------------
   module nc_control
   use kinds,    only: i4, r4, r8, l4
   use netcdf
!-------------------------------------------------------------------------------
!
   private
!
   integer(i4), parameter  :: maxnvars = 200, maxndims=4
!
   type dimension_t
     integer(i4)                      :: id
     integer(i4)                      :: dimsize
   end type dimension_t
   type variable_t
     integer(i4)                      :: id, id_in, typ
     integer(i4)                      :: ndims, natts
     integer(i4), dimension(maxndims) :: dimids
   end type variable_t
!
   type netcdf_t
     integer(i4)        :: in_nc, ou_nc
     integer(i4)        :: ndims, nvars, natts
     integer(i4)        :: compress
     character(len=512) :: infilename, oufilename
     type(dimension_t), dimension(maxnvars) :: dims
     type(variable_t),  dimension(maxnvars) :: vars
   end type netcdf_t
!
   public :: netcdf_t, nc_initialize, nc_define, nc_compress, nc_finalize
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
   integer(i4)                     :: ierr
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
   end subroutine nc_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_define(nc)
   implicit none
   type(netcdf_t), intent(inout) :: nc
! local variables
   ! attributes
   integer(i4), dimension(maxnvars)  :: aids, atype
   integer(i4)                   :: iatt
   real(r4)                      :: ratt
   character(len=64)             :: satt
   ! variables
   integer(i4)                   :: ierr, ii, jj
   integer(i4)                   :: ou_vid, did, aid
   integer(i4)                   :: vtype, ndims, natts
   integer(i4), dimension(maxndims) :: dimids
   character(len=256)            :: nam
!
! dimensions
   do ii = 1,nc%ndims
     ierr = nf90_inquire_dimension(nc%in_nc,ii,name=nam,len=nc%dims(ii)%dimsize)
     ierr = nf90_def_dim(nc%ou_nc,trim(nam),nc%dims(ii)%dimsize,nc%dims(ii)%id)
   enddo
! attributes
   do ii = 1,nc%natts
     ierr = nf90_inq_attname(nc%in_nc,nf90_global,ii,nam)
     ierr = nf90_inquire_attribute(nc%in_nc,nf90_global,trim(nam),xtype=atype(ii))
     if     (atype(ii).eq.nf90_int ) then
       ierr = nf90_get_att(nc%in_nc,nf90_global,trim(nam),iatt)
       ierr = nf90_put_att(nc%ou_nc,nf90_global,trim(nam),iatt)
     elseif (atype(ii).eq.nf90_real) then
       ierr = nf90_get_att(nc%in_nc,nf90_global,trim(nam),ratt)
       ierr = nf90_put_att(nc%ou_nc,nf90_global,trim(nam),ratt)
     elseif (atype(ii).eq.nf90_char) then
       ierr = nf90_get_att(nc%in_nc,nf90_global,trim(nam),satt)
       ierr = nf90_put_att(nc%ou_nc,nf90_global,trim(nam),trim(satt))
     else
       print *, 'undefined attributes...', atype(ii)
     endif
   enddo
   ierr = nf90_put_att(nc%ou_nc,nf90_global,'compress',nc%compress)
! variables
   do ii = 1,nc%nvars
     ierr = nf90_inquire_variable(nc%in_nc,ii,name=nam,xtype=vtype,ndims=ndims,&
                                  dimids=dimids,natts=natts)
     ierr = nf90_def_var(nc%ou_nc,trim(nam),vtype,dimids(1:ndims),ou_vid,         &
                         deflate_level=nc%compress)
     nc%vars(ii)%id_in  = ii
     nc%vars(ii)%id     = ou_vid
     nc%vars(ii)%typ    = vtype
     nc%vars(ii)%ndims  = ndims
     nc%vars(ii)%natts  = natts
     nc%vars(ii)%dimids = dimids
     ! atts
     do jj = 1,natts
       ierr = nf90_inq_attname(nc%in_nc,ii,jj,nam)
       ierr = nf90_inquire_attribute(nc%in_nc,ii,trim(nam),xtype=atype(jj))
       if     (atype(jj).eq.nf90_int ) then
         ierr = nf90_get_att(nc%in_nc,    ii,trim(nam),iatt)
         ierr = nf90_put_att(nc%ou_nc,ou_vid,trim(nam),iatt)
       elseif (atype(jj).eq.nf90_real) then
         ierr = nf90_get_att(nc%in_nc,    ii,trim(nam),ratt)
         ierr = nf90_put_att(nc%ou_nc,ou_vid,trim(nam),ratt)
       elseif (atype(jj).eq.nf90_char) then
         ierr = nf90_get_att(nc%in_nc,    ii,trim(nam),satt)
         ierr = nf90_put_att(nc%ou_nc,ou_vid,trim(nam),trim(satt))
       else
         print *, 'undefined attributes...', atype(jj)
         stop
       endif
     enddo
   enddo
   ierr = nf90_enddef(nc%ou_nc)
!
   end subroutine nc_define
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_compress(nc)
   implicit none
   type(netcdf_t), intent(inout) :: nc
! local variables
   ! variables
!   integer(i4)                   :: ierr, ii, jj
!   integer(i4)                   :: vid, did, aid
!   integer(i4)                   :: vtype, ndims, natts
!   integer(i4), dimension(maxndims) :: dimids
!   character(len=256)            :: nam
!
   integer(i4)                   :: ierr, ii, jj, vsize, vid, vtype, ndims
   integer(i4), dimension(:), allocatable :: dsizes
   integer(i4), dimension(:), allocatable :: buf_i4
   real(r4),    dimension(:), allocatable :: buf_r4
   real(r8),    dimension(:), allocatable :: buf_r8
!
! variables
   do ii = 1,nc%nvars
     vid   = nc%vars(ii)%id
     vtype = nc%vars(ii)%typ
     ndims = nc%vars(ii)%ndims
     vsize = 1
     allocate(dsizes(ndims))
     do jj =1,ndims
       dsizes(jj) = nc_get_dimsize(nc,nc%vars(ii)%dimids(jj))
       vsize      = vsize*dsizes(jj)
     enddo
     if (vsize.eq.-1.or.vsize.eq.1) then
       print *, 'cannot find dimid...'
       stop
     endif
     if     (vtype.eq.nf90_int)   then
       allocate(buf_i4(vsize))
       ierr = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_i4);print *, ierr
       if     (ndims.eq.1) then
         ierr = nf90_put_var(nc%ou_nc,vid,buf_i4)
       elseif (ndims.eq.2) then
         ierr = nf90_put_var(nc%ou_nc,vid,reshape(buf_i4,shape=(/dsizes(1),dsizes(2)/)))
       elseif (ndims.eq.3) then
         ierr = nf90_put_var(nc%ou_nc,vid,reshape(buf_i4,shape=(/dsizes(1),dsizes(2),dsizes(3)/)))
       endif
       deallocate(buf_i4)
     elseif (vtype.eq.nf90_real)  then
       allocate(buf_r4(vsize))
       ierr = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_r4);print *, ierr
       if     (ndims.eq.1) then
         ierr = nf90_put_var(nc%ou_nc,vid,buf_r4)
       elseif (ndims.eq.2) then
         ierr = nf90_put_var(nc%ou_nc,vid,reshape(buf_r4,shape=(/dsizes(1),dsizes(2)/)))
       elseif (ndims.eq.3) then
         ierr = nf90_put_var(nc%ou_nc,vid,reshape(buf_r4,shape=(/dsizes(1),dsizes(2),dsizes(3)/)))
       endif
       deallocate(buf_r4)
     elseif (vtype.eq.nf90_real8) then
       allocate(buf_r8(vsize))
       ierr = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_r8);print *, ierr
       if     (ndims.eq.1) then
         ierr = nf90_put_var(nc%ou_nc,vid,buf_r8)
       elseif (ndims.eq.2) then
         ierr = nf90_put_var(nc%ou_nc,vid,reshape(buf_r8,shape=(/dsizes(1),dsizes(2)/)))
       elseif (ndims.eq.3) then
         ierr = nf90_put_var(nc%ou_nc,vid,reshape(buf_r8,shape=(/dsizes(1),dsizes(2),dsizes(3)/)))
       endif
       deallocate(buf_r8)
     else
       print *, 'check variable type...'
       stop
     endif
     call nc_check(ierr,'check variable')
     deallocate(dsizes)
   enddo
!
   end subroutine nc_compress
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function nc_get_dimsize(nc, id) result(dimsize)
   implicit none
   type(netcdf_t),   intent(in   ) :: nc
   integer(i4),      intent(in   ) :: id
   integer(i4)                     :: dimsize
! local variables
   integer(i4) :: i
!
   dimsize = -1
   do i = 1,nc%ndims
     if (nc%dims(i)%id.eq.id) then
       dimsize = nc%dims(i)%dimsize
       exit
     endif
   enddo
!
   end function nc_get_dimsize
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
