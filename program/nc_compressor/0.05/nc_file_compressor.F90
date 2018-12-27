!-------------------------------------------------------------------------------
   module nf_file_compressor
   use kinds,    only: i4, r4, r8, l4
   use netcdf
!-------------------------------------------------------------------------------
!
   private
!
   integer(i4), parameter  :: max_ndims = 200
   integer(i4), parameter  :: max_nvars = 200
   integer(i4), parameter  :: max_ndims_var=4
!
   type dimension_t
     integer(i4)        :: id
     integer(i4)        :: dimsize
   end type dimension_t
   type variable_t
     character(len=16)  :: name
     integer(i4)        :: id, id_in, type
     integer(i4)        :: ndims, natts
     integer(i4), dimension(max_ndims_var) :: dimids
   end type variable_t
!
   type nfc_t
     integer(i4)        :: in_nc, ou_nc
     integer(i4)        :: ndims, nvars, natts
     integer(i4)        :: compress
     character(len=512) :: infilename, oufilename
     type(dimension_t), dimension(max_ndims) :: dims
     type(variable_t),  dimension(max_nvars) :: vars
   end type nfc_t
!
   public :: nfc_t, nfc_initialize, nfc_compress, nfc_finalize
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nfc_initialize(nc, infilename, oufilename, compress)
   implicit none
   type(nfc_t),      intent(inout) :: nc
   character(len=*), intent(in   ) :: infilename
   character(len=*), intent(in   ) :: oufilename
   integer(i4),      intent(in   ) :: compress
! local variables
   integer(i4)                     :: err
!
   nc%infilename = trim(infilename)
   nc%oufilename = trim(oufilename)
   nc%compress   = compress
!
   err = nf90_create(trim(oufilename),ior(nf90_clobber,nf90_netcdf4),nc%ou_nc)
!
   err = nf90_open(trim(infilename),nf90_nowrite,nc%in_nc)
   call nfc_check(err,'nfc_initialize::cannot open...')
   err = nf90_inquire(nc%in_nc,nc%ndims,nc%nvars,nc%natts)
   call nfc_check(err,'nfc_initialize::nf90_inquire')
   print *, 'infile : ', trim(infilename)
   print *, ' - ndims ', nc%ndims
   print *, ' - nvars ', nc%nvars
   print *, ' - natts ', nc%natts
!
   end subroutine nfc_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nfc_finalize(nc)
   implicit none
   type(nfc_t),   intent(inout) :: nc
! local variables
   integer(i4)  :: err
!
   err = nf90_close(nc%ou_nc)
!
   err = nf90_close(nc%in_nc)
!
   end subroutine nfc_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nfc_compress(nc)
   implicit none
   type(nfc_t),   intent(inout) :: nc
!
   call make_outfile_with_infile(nc)
!
   call put_compressed_var(nc)
!
   end subroutine nfc_compress
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_outfile_with_infile(nc)
   implicit none
   type(nfc_t), intent(inout) :: nc
! local variables
   ! attributes
   integer(i4), dimension(max_nvars) :: aids, atype
   integer(i4)                   :: iatt
   real(r4)                      :: ratt
   character(len=64)             :: satt
   ! variables
   integer(i4)                   :: err1, err2, ii, jj
   integer(i4)                   :: ou_vid, did, aid
   integer(i4)                   :: vtype, ndims, natts
   integer(i4), dimension(max_ndims_var) :: dimids
   character(len=256)            :: nam
!
! dimensions
   do ii = 1,nc%ndims
     err1 = nf90_inquire_dimension(nc%in_nc,ii,name=nam,len=nc%dims(ii)%dimsize)
     err2 = nf90_def_dim(nc%ou_nc,trim(nam),nc%dims(ii)%dimsize,nc%dims(ii)%id)
   enddo
! attributes
   do ii = 1,nc%natts
     err1 = nf90_inq_attname(nc%in_nc,nf90_global,ii,nam)
     err1 = nf90_inquire_attribute(nc%in_nc,nf90_global,trim(nam),xtype=atype(ii))
     if     (atype(ii).eq.nf90_int ) then
       err1 = nf90_get_att(nc%in_nc,nf90_global,trim(nam),iatt)
       err2 = nf90_put_att(nc%ou_nc,nf90_global,trim(nam),iatt)
     elseif (atype(ii).eq.nf90_real) then
       err1 = nf90_get_att(nc%in_nc,nf90_global,trim(nam),ratt)
       err2 = nf90_put_att(nc%ou_nc,nf90_global,trim(nam),ratt)
     elseif (atype(ii).eq.nf90_char) then
       err1 = nf90_get_att(nc%in_nc,nf90_global,trim(nam),satt)
       err2 = nf90_put_att(nc%ou_nc,nf90_global,trim(nam),trim(satt))
     else
       print *, 'undefined attributes...', atype(ii)
     endif
   enddo
   err1 = nf90_put_att(nc%ou_nc,nf90_global,'compress',nc%compress)
! variables
   do ii = 1,nc%nvars
     err1 = nf90_inquire_variable(nc%in_nc,ii,name=nam,xtype=vtype,ndims=ndims,&
                                                     dimids=dimids,natts=natts)
     err2 = nf90_def_var(nc%ou_nc,trim(nam),vtype,dimids(1:ndims),ou_vid,      &
                                                     deflate_level=nc%compress)
     nc%vars(ii)%id_in  = ii
     nc%vars(ii)%id     = ou_vid
     nc%vars(ii)%name   = trim(nam)
     nc%vars(ii)%type   = vtype
     nc%vars(ii)%ndims  = ndims
     nc%vars(ii)%natts  = natts
     nc%vars(ii)%dimids = dimids
     ! atts
     do jj = 1,natts
       err1 = nf90_inq_attname(nc%in_nc,ii,jj,nam)
       err1 = nf90_inquire_attribute(nc%in_nc,ii,trim(nam),xtype=atype(jj))
       if     (atype(jj).eq.nf90_int ) then
         err1 = nf90_get_att(nc%in_nc,    ii,trim(nam),iatt)
         err2 = nf90_put_att(nc%ou_nc,ou_vid,trim(nam),iatt)
       elseif (atype(jj).eq.nf90_real) then
         err1 = nf90_get_att(nc%in_nc,    ii,trim(nam),ratt)
         err2 = nf90_put_att(nc%ou_nc,ou_vid,trim(nam),ratt)
       elseif (atype(jj).eq.nf90_char) then
         err1 = nf90_get_att(nc%in_nc,    ii,trim(nam),satt)
         err2 = nf90_put_att(nc%ou_nc,ou_vid,trim(nam),trim(satt))
       else
         print *, 'undefined attributes...', atype(jj)
         stop
       endif
     enddo
   enddo
   err2 = nf90_enddef(nc%ou_nc)
!
   end subroutine make_outfile_with_infile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine put_compressed_var(nc)
   implicit none
   type(nfc_t), intent(inout) :: nc
! local variables
   integer(i4) :: err1, err2, ii, jj, vsize, vid, vtype, ndims
   integer(i4), dimension(:), allocatable :: dsizes
   ! buffer (i4)
   integer(i4), dimension(:),       allocatable :: buf_i4_1d
   integer(i4), dimension(:,:),     allocatable :: buf_i4_2d
   integer(i4), dimension(:,:,:),   allocatable :: buf_i4_3d
   integer(i4), dimension(:,:,:,:), allocatable :: buf_i4_4d
   real(r4),    dimension(:),       allocatable :: buf_r4_1d
   real(r4),    dimension(:,:),     allocatable :: buf_r4_2d
   real(r4),    dimension(:,:,:),   allocatable :: buf_r4_3d
   real(r4),    dimension(:,:,:,:), allocatable :: buf_r4_4d
   real(r8),    dimension(:),       allocatable :: buf_r8_1d
   real(r8),    dimension(:,:),     allocatable :: buf_r8_2d
   real(r8),    dimension(:,:,:),   allocatable :: buf_r8_3d
   real(r8),    dimension(:,:,:,:), allocatable :: buf_r8_4d
!
! variables
   do ii = 1,nc%nvars
     vid   = nc%vars(ii)%id
     vtype = nc%vars(ii)%type
     ndims = nc%vars(ii)%ndims
     vsize = 1
     allocate(dsizes(ndims))
     do jj =1,ndims
       dsizes(jj) = get_dim_size(nc,nc%vars(ii)%dimids(jj))
       vsize      = vsize*dsizes(jj)
     enddo
     if (vsize.eq.-1.or.vsize.eq.1) then
       print *, 'cannot find dimid...'
       stop
     endif
     if     (vtype.eq.nf90_int)   then
       if     (ndims.eq.1) then
         allocate(buf_i4_1d(dsizes(1)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_i4_1d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_i4_1d)
         deallocate(buf_i4_1d)
       elseif (ndims.eq.2) then
         allocate(buf_i4_2d(dsizes(1),dsizes(2)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_i4_2d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_i4_2d)
         deallocate(buf_i4_2d)
       elseif (ndims.eq.3) then
         allocate(buf_i4_3d(dsizes(1),dsizes(2),dsizes(3)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_i4_3d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_i4_3d)
         deallocate(buf_i4_3d)
       elseif (ndims.eq.4) then
         allocate(buf_i4_4d(dsizes(1),dsizes(2),dsizes(3),dsizes(4)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_i4_4d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_i4_4d)
         deallocate(buf_i4_4d)
       endif
     elseif (vtype.eq.nf90_real)  then
       if     (ndims.eq.1) then
         allocate(buf_r4_1d(dsizes(1)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_r4_1d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_r4_1d)
         deallocate(buf_r4_1d)
       elseif (ndims.eq.2) then
         allocate(buf_r4_2d(dsizes(1),dsizes(2)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_r4_2d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_r4_2d)
         deallocate(buf_r4_2d)
       elseif (ndims.eq.3) then
         allocate(buf_r4_3d(dsizes(1),dsizes(2),dsizes(3)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_r4_3d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_r4_3d)
         deallocate(buf_r4_3d)
       elseif (ndims.eq.4) then
         allocate(buf_r4_4d(dsizes(1),dsizes(2),dsizes(3),dsizes(4)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_r4_4d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_r4_4d)
         deallocate(buf_r4_4d)
       endif
     elseif (vtype.eq.nf90_real8) then
       if     (ndims.eq.1) then
         allocate(buf_r8_1d(dsizes(1)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_r8_1d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_r8_1d)
         deallocate(buf_r8_1d)
       elseif (ndims.eq.2) then
         allocate(buf_r8_2d(dsizes(1),dsizes(2)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_r8_2d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_r8_2d)
         deallocate(buf_r8_2d)
       elseif (ndims.eq.3) then
         allocate(buf_r8_3d(dsizes(1),dsizes(2),dsizes(3)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_r8_3d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_r8_3d)
         deallocate(buf_r8_3d)
       elseif (ndims.eq.4) then
         allocate(buf_r8_4d(dsizes(1),dsizes(2),dsizes(3),dsizes(4)))
         err1 = nf90_get_var(nc%in_nc,nc%vars(ii)%id_in,buf_r8_4d)
         err2 = nf90_put_var(nc%ou_nc,vid,buf_r8_4d)
         deallocate(buf_r8_4d)
       endif
     else
       print *, 'check variable type...'
       stop
     endif
     print *, 'processing... var :', trim(nc%vars(ii)%name)
     call nfc_check(err1,'check get variable')
     call nfc_check(err2,'check put variable')
     deallocate(dsizes)
   enddo
!
   end subroutine put_compressed_var
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_dim_size(nc, id) result(dimsize)
   implicit none
   type(nfc_t), intent(in   ) :: nc
   integer(i4), intent(in   ) :: id
   integer(i4)                :: dimsize
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
   end function get_dim_size
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nfc_check(err,message)
   implicit none
   integer(i4),                intent(in   ) :: err
   character(len=*), optional, intent(in   ) :: message
! local variables
!
   if (err.ne.nf90_noerr) then
     if (present(message)) print *, trim(message)
     stop
   endif
!
   end subroutine nfc_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module nf_file_compressor
!-------------------------------------------------------------------------------
