!-------------------------------------------------------------------------------
   module netcdf_utility
!-------------------------------------------------------------------------------
!
!  abstract :  NetCDF file compressor utility
!
!  history log :
!    2017-02-28  junghan kim    initial setup
!
!  structure :
!
!-------------------------------------------------------------------------------
   use kinds,    only: i4, i8, r4, r8, l4
   use netcdf
   use pnetcdf
   use pnetcdf,  only: nf90mpi_clobber=>nf90_clobber,nf90mpi_noerr=>nf90_noerr,&
                       nf90mpi_64bit_offset=>nf90_64bit_offset,                &
                       nf90mpi_nowrite=>nf90_nowrite,                          &
                       nf90mpi_global=>nf90_global,                            &
                       nf90mpi_int=>nf90_int, nf90mpi_int64=>nf90_int64,       &
                       nf90mpi_real=>nf90_real,nf90mpi_real8=>nf90_real8,      &
                       nf90mpi_char=>nf90_char
!-------------------------------------------------------------------------------
!
   private
!
   include 'mpif.h'
!
   integer(i4), parameter  :: max_ndims =  20
   integer(i4), parameter  :: max_natts =  40
   integer(i4), parameter  :: max_nvars = 100
   integer(i4), parameter  :: max_ndims_var =  4
   integer(i4), parameter  :: max_natts_var = 10
   integer(i4), parameter  :: len_att_nam   = 32
   integer(i4), parameter  :: len_att_val   = 256
   integer(i4), parameter  :: len_dim_nam   = 16
   integer(i4), parameter  :: len_var_nam   = 32
   integer(i4), parameter  :: len_filenam   = 512
!
   type dimension_t
     character(len_dim_nam) :: name
     integer(i4)            :: id
     integer(i4)            :: size
   end type dimension_t
   type attribute_t
     character(len_att_nam) :: name
     integer(i4)            :: id
     integer(i4)            :: type
     integer(i4)            :: iatt
     integer(i8)            :: latt
     real(r4)               :: ratt
     real(r8)               :: datt
     character(len_att_val) :: satt
   end type attribute_t
   type variable_t
     character(len_var_nam) :: name
     logical(l4)            :: selected, large
     integer(i4)            :: id, id_in, type
     integer(i4)            :: ndims, natts
     integer(i4)            :: vsize
     type(dimension_t), dimension(max_ndims_var) :: dims
     type(attribute_t), dimension(max_natts_var) :: atts
   end type variable_t
!
   type ncu_t
     integer(i4)            :: comm, nprocs, rank
     character(len_filenam) :: infile, oufile
     integer(i4)            :: fin, fou
     integer(i4)            :: ndims, nvars, natts
     integer(i4)            :: compress
     integer(i4)            :: cvt        ! 0(off), 1(data->nc4), 2(nc4->data)
     integer(i4)            :: tin, tou   ! type: 1(nc4), 2(pnc)
     character(len_var_nam), dimension(:), allocatable :: varnames
     type(dimension_t), dimension(max_ndims) :: dims
     type(attribute_t), dimension(max_natts) :: atts
     type(variable_t) , dimension(max_nvars) :: vars
     integer(i4), dimension(:),       allocatable :: buf_i4_1d
     integer(i4), dimension(:,:),     allocatable :: buf_i4_2d
     integer(i4), dimension(:,:,:),   allocatable :: buf_i4_3d
     integer(i4), dimension(:,:,:,:), allocatable :: buf_i4_4d
     integer(i8), dimension(:),       allocatable :: buf_i8_1d
     integer(i8), dimension(:,:),     allocatable :: buf_i8_2d
     integer(i8), dimension(:,:,:),   allocatable :: buf_i8_3d
     integer(i8), dimension(:,:,:,:), allocatable :: buf_i8_4d
     real(r4),    dimension(:),       allocatable :: buf_r4_1d
     real(r4),    dimension(:,:),     allocatable :: buf_r4_2d
     real(r4),    dimension(:,:,:),   allocatable :: buf_r4_3d
     real(r4),    dimension(:,:,:,:), allocatable :: buf_r4_4d
     real(r8),    dimension(:),       allocatable :: buf_r8_1d
     real(r8),    dimension(:,:),     allocatable :: buf_r8_2d
     real(r8),    dimension(:,:,:),   allocatable :: buf_r8_3d
     real(r8),    dimension(:,:,:,:), allocatable :: buf_r8_4d
   end type ncu_t
!
   public :: len_filenam, len_att_nam, len_var_nam
   public :: ncu_t, ncu_initialize, ncu_run, ncu_finalize
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ncu_initialize(nc, infile, oufile, compress, convert, varnames)
   implicit none
   type(ncu_t)     , intent(inout) :: nc
   character(len=*), intent(in   ) :: infile
   character(len=*), intent(in   ) :: oufile
   integer(i4)     , intent(in   ) :: compress, convert
   character(len_var_nam), dimension(:), optional, intent(in   ) :: varnames
! local variables
   integer(i4) :: err
!
   call mpi_init(err)
   nc%comm = mpi_comm_world
   call mpi_comm_size(nc%comm,nc%nprocs,err)
   call mpi_comm_rank(nc%comm,  nc%rank,err)
!
   if (nc%nprocs.ne.1) then
     if (nc%rank.eq.0) then
       write(*,'(a)')'this program must be excuted by only 1 process'
     endif
     call mpi_abort(nc%comm)
   endif
!
   nc%infile = trim(infile)
   nc%oufile = trim(oufile)
   nc%compress   = compress
   nc%cvt        = convert
   if     (nc%cvt.eq.1) then
     nc%tin = 4
     nc%tou = 3
   elseif (nc%cvt.eq.2) then
     nc%tin = 3
     nc%tou = 4
   else
     nc%tin = 3
     nc%tou = 3
   endif
   if (present(varnames)) then
     allocate(nc%varnames(size(varnames)))
     nc%varnames = varnames
   else
     allocate(nc%varnames(1))
     nc%varnames(1) = 'all'
   endif
!
   call nc_create(nc)
   call nc_open(nc)
   call nc_inquire(nc)
   call ncu_initialize_log(nc)
!
   end subroutine ncu_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ncu_initialize_log(nc)
   implicit none
   type(ncu_t), intent(in   ) :: nc
! local variables
!
   write(*,'(a,a)')   '# initialize: file name: ',trim(nc%infile)
   write(*,'(a,x,i3)')' * ndims: ',nc%ndims
   write(*,'(a,x,i3)')' * natts: ',nc%natts
   write(*,'(a,x,i3)')' * nvars: ',nc%nvars
!
   end subroutine ncu_initialize_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ncu_finalize(nc)
   implicit none
   type(ncu_t), intent(inout) :: nc
! local variables
   integer(i4)  :: err
!
   if     (nc%tin.eq.4) then
     err = nf90mpi_close(nc%fin)
   else
     err = nf90_close(nc%fin)
   endif
   if     (nc%tou.eq.4) then
     err = nf90mpi_close(nc%fou)
   else
     err = nf90_close(nc%fou)
   endif
!
   call mpi_finalize(err)
!
   end subroutine ncu_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ncu_run(nc)
   implicit none
   type(ncu_t), intent(inout) :: nc
!
   call put_metadata(nc)
   call put_vars(nc)
!
   end subroutine ncu_run
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine put_metadata(nc)
   implicit none
   type(ncu_t), intent(inout) :: nc
! local variables
   integer(i4) :: i, j
!
   write(*,'(a)')   '# read/write metadata'
!
! dimensions (read)
   write(*,'(a)')' * read  dimensions.'
   do i = 1,nc%ndims
     call nc_inquire_dimension(nc,i)
   enddo
! dimensions (write)
   write(*,'(a)')' * write dimensions.'
   do i = 1,nc%ndims
     call nc_def_dim(nc,i)
   enddo
!
! attributes (read)
   write(*,'(a)')' * read  global attributes.'
! gnu bug
!#if 0
   do j = 1,nc%natts
     nc%atts(i)%id = j
     call nc_inq_attname(nc,nf90_global,j)
     call nc_inquire_attribute(nc,nf90_global,j)
     !print *,'DEBUG: before j = ',j
     call nc_get_att(nc,nf90_global,j)
     !print *,'DEBUG: after  j = ',j
   enddo
! attributes (write)
   write(*,'(a)')' * write global attributes.'
   do i = 1,nc%natts
     call nc_put_att(nc,nf90_global,i)
   enddo
   call nc_put_att_compress(nc)
!#endif
!
! variables (read)
   write(*,'(a)')' * read  metadata in variables.'
   do i = 1,nc%nvars
     nc%vars(i)%id_in = i
     call nc_inquire_variable(nc,i)
     nc%vars(i)%selected = selected(nc,trim(nc%vars(i)%name))
     !
     if (nc%vars(i)%selected) then
       ! atts
       do j = 1,nc%vars(i)%natts
         call nc_inq_attname(nc,i,j)
         call nc_inquire_attribute(nc,i,j)
         call nc_get_att(nc,i,j)
       enddo
     endif
     nc%vars(i)%large = .true.
     !
   enddo
! variables (write)
   write(*,'(a)')' * write metadata in variables.'
!   write(*,*) ' '
   do i = 1,nc%nvars
     if (nc%vars(i)%selected) then
       call nc_def_var(nc,i)
       ! atts
       do j = 1,nc%vars(i)%natts
         call nc_put_att(nc,i,j)
       enddo
     endif
   enddo
   call nc_enddef(nc)
   call ncu_metadata_log(nc)
!
   end subroutine put_metadata
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ncu_metadata_log(nc)
   implicit none
   type(ncu_t), intent(in   ) :: nc
! local variables
   integer(i4) :: i, ndims
   character(len=32) :: sfmt
!
   write(*,'(a)')' * file information'
   write(*,'(a,x,i3)')'  - ndims: ',nc%ndims
   do i = 1,nc%ndims
     write(*,'(a10,i3,i15)') trim(nc%dims(i)%name),nc%dims(i)%id,nc%dims(i)%size
   enddo
   write(*,'(a,x,i3)')'  - natts: ',nc%natts
   do i = 1,nc%natts
     if     (nc%atts(i)%type.eq.nf90_int) then
       write(*,'(a10,a10,x,i16)')   trim(nc%atts(i)%name),'int4 ',nc%atts(i)%iatt
     elseif (nc%atts(i)%type.eq.nf90_int64) then
       write(*,'(a10,a10,x,i16)')   trim(nc%atts(i)%name),'int8 ',nc%atts(i)%latt
     elseif (nc%atts(i)%type.eq.nf90_real) then
       write(*,'(a10,a10,x,f16.2)') trim(nc%atts(i)%name),'real4',nc%atts(i)%ratt
     elseif (nc%atts(i)%type.eq.nf90_real8) then
       write(*,'(a10,a10,x,f16.2)') trim(nc%atts(i)%name),'real8',nc%atts(i)%datt
     elseif (nc%atts(i)%type.eq.nf90_char) then
       write(*,'(a10,a10,x,a16)')   trim(nc%atts(i)%name),'char ',nc%atts(i)%satt
     endif
   enddo
   write(*,'(a,x,i3)')'  - nvars: ',nc%nvars
   do i = 1,nc%nvars
     ndims = nc%vars(i)%ndims
     if     (ndims.eq.1) then
       sfmt = '(a10,i4,2x,a,a6,a,i4,a)'
     elseif (ndims.eq.2) then
       sfmt = '(a10,i4,2x,a,2a6,a,2i4,a)'
     elseif (ndims.eq.3) then
       sfmt = '(a10,i4,2x,a,3a6,a,3i4,a)'
     elseif (ndims.eq.4) then
       sfmt = '(a10,i4,2x,a,4a6,a,4i4,a)'
     endif
     write(*,trim(sfmt)) trim(nc%vars(i)%name),nc%vars(i)%ndims,               &
       '(',nc%vars(i)%dims(1:ndims)%name,'), (',nc%vars(i)%dims(1:ndims)%size,')'
   enddo
!
!
   end subroutine ncu_metadata_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine put_vars(nc)
   implicit none
   type(ncu_t), intent(inout) :: nc
! local variables
   integer(i4) :: err, i, j
!
! variables
   write(*,'(a)')'# read/write variables'
   do i = 1,nc%nvars
     if (nc%vars(i)%selected) then
       write(*,'(a,a)')'  - processing variable: ',trim(nc%vars(i)%name)
       call nc_get_var(nc,i)
       call nc_put_var(nc,i)
     endif
   enddo
!
   end subroutine put_vars
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function selected(nc, varname) result(res)
   implicit none
   type(ncu_t),      intent(in   ) :: nc
   character(len=*), intent(in   ) :: varname
   logical(l4)                     :: res
! local variables
   integer(i4) :: i, nvars
!
   res   = .false.
   nvars = size(nc%varnames)
!
   if (nvars.eq.1.and.trim(nc%varnames(1)).eq.'all') then
     res = .true.
   else
     do i = 1,nvars
       if (trim(varname).eq.trim(nc%varnames(i))) then
         res = .true.
       endif
     enddo
   endif
!
   end function selected
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_create(nc)
   implicit none
   type(ncu_t), intent(inout) :: nc
! local variables
   integer(i4) :: err
!
   if (nc%tou.eq.4) then
     err = nf90mpi_create(nc%comm,trim(nc%oufile),ior(nf90_clobber,nf90_64bit_data), &
                          mpi_info_null,nc%fou)
   else
     err = nf90_create(trim(nc%oufile),ior(nf90_clobber,nf90_netcdf4),nc%fou)
   endif
   call nc_check(err,'nc_create: cannot create file...')
!
   end subroutine nc_create
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_open(nc)
   implicit none
   type(ncu_t), intent(inout) :: nc
! local variables
   integer(i4) :: err
!
   if (nc%tin.eq.4) then
     err = nf90mpi_open(nc%comm,trim(nc%infile),nf90_nowrite,mpi_info_null,nc%fin)
   else
     err = nf90_open(trim(nc%infile),nf90_nowrite,nc%fin)
   endif
   call nc_check(err,'nc_open: cannot open file...')
!
   end subroutine nc_open
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_inquire(nc)
   implicit none
   type(ncu_t), intent(inout) :: nc
! local variables
   integer(i4) :: err
!
   if (nc%tin.eq.4) then
     err = nf90mpi_inquire(nc%fin,nc%ndims,nc%nvars,nc%natts)
   else
     err = nf90_inquire(nc%fin,nc%ndims,nc%nvars,nc%natts)
   endif
   call nc_check(err,'nc_inquire')
!
   end subroutine nc_inquire
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_inquire_dimension(nc, i)
   implicit none
   type(ncu_t), intent(inout) :: nc
   integer(i4), intent(in   ) :: i
! local variables
   integer(i4) :: err
   integer( 8) :: tmp
!
   if (nc%tin.eq.4) then
     err = nf90mpi_inquire_dimension(nc%fin,i,name=nc%dims(i)%name,len=tmp)
     nc%dims(i)%size = tmp
   else
     err = nf90_inquire_dimension(nc%fin,i,name=nc%dims(i)%name,len=nc%dims(i)%size)
   endif
   call nc_check(err,'nc_inquire_dimension')
!
   end subroutine nc_inquire_dimension
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_def_dim(nc,i)
   implicit none
   type(ncu_t), intent(inout) :: nc
   integer(i4), intent(in   ) :: i
! local variables
   integer(i4) :: err
!
   if (nc%tou.eq.4) then
     err = nf90mpi_def_dim(nc%fou,trim(nc%dims(i)%name),int(nc%dims(i)%size,8),nc%dims(i)%id)
   else
     err = nf90_def_dim(nc%fou,trim(nc%dims(i)%name),nc%dims(i)%size,nc%dims(i)%id)
   endif
   call nc_check(err,'nc_def_dim')
!
   end subroutine nc_def_dim
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_inq_attname(nc,i,j)
   implicit none
   type(ncu_t), intent(inout) :: nc
   integer(i4), intent(in   ) :: i, j
! local variables
   integer(i4) :: err
!
   if (nc%tin.eq.4) then
     if (i.eq.nf90_global) then
       err = nf90mpi_inq_attname(nc%fin,i,j,nc%atts(j)%name)
     else
       err = nf90mpi_inq_attname(nc%fin,i,j,nc%vars(i)%atts(j)%name)
     endif
   else
     if (i.eq.nf90_global) then
       err = nf90_inq_attname(nc%fin,i,j,nc%atts(j)%name)
     else
       err = nf90_inq_attname(nc%fin,i,j,nc%vars(i)%atts(j)%name)
     endif
   endif
   call nc_check(err,'nc_inq_attname')
!
   end subroutine nc_inq_attname
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_inquire_attribute(nc,i,j)
   implicit none
   type(ncu_t), intent(inout) :: nc
   integer(i4), intent(in   ) :: i,j
! local variables
   integer(i4) :: err
!
   if (nc%tin.eq.4) then
     if (i.eq.nf90_global) then
       err = nf90mpi_inquire_attribute(nc%fin,i,trim(nc%atts(j)%name),xtype=nc%atts(j)%type)
     else
       err = nf90mpi_inquire_attribute(nc%fin,i,trim(nc%vars(i)%atts(j)%name),xtype=nc%vars(i)%atts(j)%type)
     endif
   else
     if (i.eq.nf90_global) then
       err = nf90_inquire_attribute(nc%fin,i,trim(nc%atts(j)%name),xtype=nc%atts(j)%type)
     else
       err = nf90_inquire_attribute(nc%fin,i,trim(nc%vars(i)%atts(j)%name),xtype=nc%vars(i)%atts(j)%type)
     endif
   endif
   call nc_check(err,'')
!
   end subroutine nc_inquire_attribute
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_get_att(nc,i,j)
   implicit none
   type(ncu_t), intent(inout) :: nc
   integer(i4), intent(in   ) :: i, j
! local variables
   integer(i4) :: err, atype
   character(len_att_nam) :: attname
   integer(i4)            :: iatt
   integer(i8)            :: latt
   real(r4)               :: ratt
   real(r8)               :: datt
   character(len_att_val) :: satt
!
   if (i.eq.nf90_global) then
     attname = nc%atts(j)%name
     atype   = nc%atts(j)%type
   else
     attname = nc%vars(i)%atts(j)%name
     atype   = nc%vars(i)%atts(j)%type
   endif
   if (atype.ne.nf90_int.and.atype.ne.nf90_int64.and.atype.ne.nf90_real.and.   &
                               atype.ne.nf90_real8.and.atype.ne.nf90_char) then
     write(*,'(a,i2)') 'undefined attributes in nc_get_att: ',atype
     call mpi_abort(nc%comm)
   endif
!
   if (nc%tin.eq.4) then
     if     (atype.eq.nf90_int ) then
       err = nf90mpi_get_att(nc%fin,i,trim(attname),iatt)
     elseif (atype.eq.nf90_int64) then
       err = nf90mpi_get_att(nc%fin,i,trim(attname),latt)
     elseif (atype.eq.nf90_real) then
       err = nf90mpi_get_att(nc%fin,i,trim(attname),ratt)
     elseif (atype.eq.nf90_real8) then
       err = nf90mpi_get_att(nc%fin,i,trim(attname),datt)
     elseif (atype.eq.nf90_char) then
       err = nf90mpi_get_att(nc%fin,i,trim(attname),satt)
     endif     
   else
     if     (atype.eq.nf90_int ) then
       err = nf90_get_att(nc%fin,i,trim(attname),iatt)
     elseif (atype.eq.nf90_int64) then
       err = nf90_get_att(nc%fin,i,trim(attname),latt)
     elseif (atype.eq.nf90_real) then
       err = nf90_get_att(nc%fin,i,trim(attname),ratt)
     elseif (atype.eq.nf90_real8) then
       err = nf90_get_att(nc%fin,i,trim(attname),datt)
     elseif (atype.eq.nf90_char) then
       err = nf90_get_att(nc%fin,i,trim(attname),satt)
     endif     
   endif
   call nc_check(err,'nc_get_att')
!
   if (i.eq.nf90_global) then
     if     (atype.eq.nf90_int ) then
       nc%atts(j)%iatt = iatt
     elseif (atype.eq.nf90_int64) then
       nc%atts(j)%latt = latt
     elseif (atype.eq.nf90_real) then
       nc%atts(j)%ratt = ratt
     elseif (atype.eq.nf90_real8) then
       nc%atts(j)%datt = datt
     elseif (atype.eq.nf90_char) then
       nc%atts(j)%satt = satt
     endif
   else
     if     (atype.eq.nf90_int ) then
       nc%vars(i)%atts(j)%iatt = iatt
     elseif (atype.eq.nf90_int64) then
       nc%vars(i)%atts(j)%latt = latt
     elseif (atype.eq.nf90_real) then
       nc%vars(i)%atts(j)%ratt = ratt
     elseif (atype.eq.nf90_real8) then
       nc%vars(i)%atts(j)%datt = datt
     elseif (atype.eq.nf90_char) then
       nc%vars(i)%atts(j)%satt = satt
     endif
   endif
!
   end subroutine nc_get_att
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_put_att(nc,i,j)
   implicit none
   type(ncu_t), intent(inout) :: nc
   integer(i4), intent(in   ) :: i, j
! local variables
   integer(i4) :: err, atype
   character(len_att_nam) :: attname
   integer(i4)            :: iatt
   integer(i8)            :: latt
   real(r4)               :: ratt
   real(r8)               :: datt
   character(len_att_val) :: satt
!
   if (i.eq.nf90_global) then
     attname = trim(nc%atts(j)%name)
     atype   = nc%atts(j)%type
   else
     attname = nc%vars(i)%atts(j)%name
     atype   = nc%vars(i)%atts(j)%type
   endif
   if (atype.ne.nf90_int.and.atype.ne.nf90_int64.and.atype.ne.nf90_real.and.   &
                               atype.ne.nf90_real8.and.atype.ne.nf90_char) then
     write(*,'(a,i2)') 'undefined attributes in nc_put_att: ',atype
     call mpi_abort(nc%comm)
   endif
!
   if (i.eq.nf90_global) then
     if     (atype.eq.nf90_int ) then
       iatt = nc%atts(j)%iatt
     elseif (atype.eq.nf90_int64) then
       latt = nc%atts(j)%latt
     elseif (atype.eq.nf90_real) then
       ratt = nc%atts(j)%ratt
     elseif (atype.eq.nf90_real8) then
       datt = nc%atts(j)%datt
     elseif (atype.eq.nf90_char) then
       satt = nc%atts(j)%satt
     endif
   else
     if     (atype.eq.nf90_int ) then
       iatt = nc%vars(i)%atts(j)%iatt
     elseif (atype.eq.nf90_int64) then
       latt = nc%vars(i)%atts(j)%latt
     elseif (atype.eq.nf90_real) then
       ratt = nc%vars(i)%atts(j)%ratt
     elseif (atype.eq.nf90_real8) then
       datt = nc%vars(i)%atts(j)%datt
     elseif (atype.eq.nf90_char) then
       satt = nc%vars(i)%atts(j)%satt
     endif
   endif
!
   if (nc%tou.eq.4) then
     if     (atype.eq.nf90_int ) then
       err = nf90mpi_put_att(nc%fou,i,trim(attname),iatt)
     elseif (atype.eq.nf90_int64) then
       err = nf90mpi_put_att(nc%fou,i,trim(attname),latt)
     elseif (atype.eq.nf90_real) then
       err = nf90mpi_put_att(nc%fou,i,trim(attname),ratt)
     elseif (atype.eq.nf90_real8) then
       err = nf90mpi_put_att(nc%fou,i,trim(attname),datt)
     elseif (atype.eq.nf90_char) then
       err = nf90mpi_put_att(nc%fou,i,trim(attname),satt)
     endif     
   else
     if     (atype.eq.nf90_int ) then
       err = nf90_put_att(nc%fou,i,trim(attname),iatt)
     elseif (atype.eq.nf90_int64) then
       err = nf90_put_att(nc%fou,i,trim(attname),latt)
     elseif (atype.eq.nf90_real) then
       err = nf90_put_att(nc%fou,i,trim(attname),ratt)
     elseif (atype.eq.nf90_real8) then
       err = nf90_put_att(nc%fou,i,trim(attname),datt)
     elseif (atype.eq.nf90_char) then
       err = nf90_put_att(nc%fou,i,trim(attname),satt)
     endif     
   endif
   call nc_check(err,'nc_put_att')
!
   end subroutine nc_put_att
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_put_att_compress(nc)
   implicit none
   type(ncu_t), intent(inout) :: nc
! local variables
   integer(i4) :: err
!
   if (nc%tou.eq.4) then
     err = nf90mpi_put_att(nc%fou,nf90_global,'compress',nc%compress)
   else
     err = nf90_put_att(nc%fou,nf90_global,'compress',nc%compress)
   endif
   call nc_check(err,'nc_put_att_compress')
!
   end subroutine nc_put_att_compress
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_enddef(nc)
   implicit none
   type(ncu_t), intent(inout) :: nc
! local variables
   integer(i4) :: err
!
   if (nc%tou.eq.4) then
      err = nf90mpi_enddef(nc%fou)
   else
      err = nf90_enddef(nc%fou)
   endif
   call nc_check(err,'')
!
   end subroutine nc_enddef
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_inquire_variable(nc,i)
   implicit none
   type(ncu_t), intent(inout) :: nc
   integer(i4), intent(in   ) :: i
! local variables
   integer(i4) :: err, j, dimids(max_ndims_var)
!
   if (nc%tin.eq.4) then
     err = nf90mpi_inquire_variable(nc%fin,i,name=nc%vars(i)%name,           &
                             xtype=nc%vars(i)%type,ndims=nc%vars(i)%ndims,   &
                             dimids=dimids,natts=nc%vars(i)%natts)
   else
     err = nf90_inquire_variable(nc%fin,i,name=nc%vars(i)%name,              &
                             xtype=nc%vars(i)%type,ndims=nc%vars(i)%ndims,   &
                             dimids=dimids,natts=nc%vars(i)%natts)
   endif
   call nc_check(err,'nc_inquire_variable')
!
   nc%vars(i)%vsize = 1
   do j = 1,nc%vars(i)%ndims
     nc%vars(i)%dims(j)%id = dimids(j)
     nc%vars(i)%dims(j)%size = nc%dims(nc%vars(i)%dims(j)%id)%size
     nc%vars(i)%dims(j)%name = nc%dims(nc%vars(i)%dims(j)%id)%name
     nc%vars(i)%vsize = nc%vars(i)%vsize*nc%vars(i)%dims(j)%size
   enddo
!
   
!
   end subroutine nc_inquire_variable
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_def_var(nc,i)
   implicit none
   type(ncu_t), intent(inout) :: nc
   integer(i4), intent(in   ) :: i
! local variables
   integer(i4) :: err, j, dimids(max_ndims_var)
!
   do j = 1,nc%vars(i)%ndims
     dimids(j) = nc%vars(i)%dims(j)%id
   enddo
!
   if (nc%tou.eq.4) then
     err = nf90mpi_def_var(nc%fou,nc%vars(i)%name,nc%vars(i)%type,        &
                        dimids(1:nc%vars(i)%ndims), nc%vars(i)%id)
   else
     err = nf90_def_var(nc%fou,nc%vars(i)%name,nc%vars(i)%type,           &
                        dimids(1:nc%vars(i)%ndims),                       &
                        nc%vars(i)%id,deflate_level=nc%compress)
   endif
   call nc_check(err,'nc_def_var')
!
   end subroutine nc_def_var
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_get_var(nc,i)
   implicit none
   type(ncu_t), intent(inout) :: nc
   integer(i4), intent(in   ) :: i
! local variables
   integer(i4) :: err, vtype, ndims, dim1, dim2, dim3, dim4
!
   vtype = nc%vars(i)%type
   ndims = nc%vars(i)%ndims
   dim1 = nc%vars(i)%dims(1)%size
   dim2 = nc%vars(i)%dims(2)%size
   dim3 = nc%vars(i)%dims(3)%size
   dim4 = nc%vars(i)%dims(4)%size
!
   if (nc%tin.eq.4) then
     if     (vtype.eq.nf90_int) then
       if     (ndims.eq.1) then
         allocate(nc%buf_i4_1d(dim1))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_i4_1d)
       elseif (ndims.eq.2) then
         allocate(nc%buf_i4_2d(dim1,dim2))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_i4_2d)
       elseif (ndims.eq.3) then
         allocate(nc%buf_i4_3d(dim1,dim2,dim3))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_i4_3d)
       elseif (ndims.eq.4) then
         allocate(nc%buf_i4_4d(dim1,dim2,dim3,dim4))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_i4_4d)
       endif
     elseif (vtype.eq.nf90_real) then
       if     (ndims.eq.1) then
         allocate(nc%buf_r4_1d(dim1))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_r4_1d)
       elseif (ndims.eq.2) then
         allocate(nc%buf_r4_2d(dim1,dim2))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_r4_2d)
       elseif (ndims.eq.3) then
         allocate(nc%buf_r4_3d(dim1,dim2,dim3))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_r4_3d)
       elseif (ndims.eq.4) then
         allocate(nc%buf_r4_4d(dim1,dim2,dim3,dim4))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_r4_4d)
       endif
     elseif (vtype.eq.nf90_real8) then
       if     (ndims.eq.1) then
         allocate(nc%buf_r8_1d(dim1))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_r8_1d)
       elseif (ndims.eq.2) then
         allocate(nc%buf_r8_2d(dim1,dim2))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_r8_2d)
       elseif (ndims.eq.3) then
         allocate(nc%buf_r8_3d(dim1,dim2,dim3))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_r8_3d)
       elseif (ndims.eq.4) then
         allocate(nc%buf_r8_4d(dim1,dim2,dim3,dim4))
         err = nf90mpi_get_var_all(nc%fin,nc%vars(i)%id_in,nc%buf_r8_4d)
       endif
     endif
   else ! not pnetcdf
     if     (vtype.eq.nf90_int) then
       if     (ndims.eq.1) then
         allocate(nc%buf_i4_1d(dim1))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_i4_1d)
       elseif (ndims.eq.2) then
         allocate(nc%buf_i4_2d(dim1,dim2))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_i4_2d)
       elseif (ndims.eq.3) then
         allocate(nc%buf_i4_3d(dim1,dim2,dim3))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_i4_3d)
       elseif (ndims.eq.4) then
         allocate(nc%buf_i4_4d(dim1,dim2,dim3,dim4))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_i4_4d)
       endif
     elseif (vtype.eq.nf90_real) then
       if     (ndims.eq.1) then
         allocate(nc%buf_r4_1d(dim1))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_r4_1d)
       elseif (ndims.eq.2) then
         allocate(nc%buf_r4_2d(dim1,dim2))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_r4_2d)
       elseif (ndims.eq.3) then
         allocate(nc%buf_r4_3d(dim1,dim2,dim3))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_r4_3d)
       elseif (ndims.eq.4) then
         allocate(nc%buf_r4_4d(dim1,dim2,dim3,dim4))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_r4_4d)
       endif
     elseif (vtype.eq.nf90_real8) then
       if     (ndims.eq.1) then
         allocate(nc%buf_r8_1d(dim1))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_r8_1d)
       elseif (ndims.eq.2) then
         allocate(nc%buf_r8_2d(dim1,dim2))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_r8_2d)
       elseif (ndims.eq.3) then
         allocate(nc%buf_r8_3d(dim1,dim2,dim3))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_r8_3d)
       elseif (ndims.eq.4) then
         allocate(nc%buf_r8_4d(dim1,dim2,dim3,dim4))
         err = nf90_get_var(nc%fin,nc%vars(i)%id_in,nc%buf_r8_4d)
       endif
     endif
   endif
   call nc_check(err,'nc_get_var')
!
   end subroutine nc_get_var
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_put_var(nc,i)
   implicit none
   type(ncu_t), intent(inout) :: nc
   integer(i4), intent(in   ) :: i
! local variables
   integer(i4) :: err, vtype, ndims
   integer(i4) :: dim1, dim2, dim3, dim4
   integer(i4) :: k
   integer(i8) :: bsiz
   integer(i8), allocatable :: sta(:,:), cnt(:,:)
!
   vtype = nc%vars(i)%type
   ndims = nc%vars(i)%ndims
!
   if (nc%tou.eq.4) then
     if (nc%vars(i)%large) then
       dim1 = nc%vars(i)%dims(1)%size
       dim2 = nc%vars(i)%dims(2)%size
       dim3 = nc%vars(i)%dims(3)%size
       dim4 = nc%vars(i)%dims(4)%size
       if     (ndims.eq.1) then
         allocate(sta(ndims,1))
         allocate(cnt(ndims,1))
         sta(:,1) = (/1/)
         cnt(:,1) = (/dim1/)
         bsiz = dim1
       elseif (ndims.eq.2) then
         allocate(sta(ndims,dim2))
         allocate(cnt(ndims,dim2))
         do k = 1,dim2
           sta(:,k) = (/1,k/)
           cnt(:,k) = (/dim1,1/)
         enddo
         bsiz = dim1
       elseif (ndims.eq.3) then
         allocate(sta(ndims,dim3))
         allocate(cnt(ndims,dim3))
         do k = 1,dim3
           sta(:,k) = (/1,1,k/)
           cnt(:,k) = (/dim1,dim2,1/)
         enddo
         bsiz = dim1*dim2
       elseif (ndims.eq.4) then
         allocate(sta(ndims,dim4))
         allocate(cnt(ndims,dim4))
         do k = 1,dim4
           sta(:,k) = (/1,1,1,k/)
           cnt(:,k) = (/dim1,dim2,dim3,1/)
         enddo
         bsiz = dim1*dim2*dim3
       endif
     endif
     if     (vtype.eq.nf90_int) then
       if     (ndims.eq.1) then
         !err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_i4_1d)
         err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,1),cnt(:,1),      &
                                  nc%buf_i4_1d(:),bsiz,mpi_int)
         deallocate(nc%buf_i4_1d)
       elseif (ndims.eq.2) then
         if (nc%vars(i)%large) then
           do k = 1,dim2
             err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,k),cnt(:,k),&
                                        nc%buf_i4_2d(:,k),bsiz,mpi_int)
           enddo
         else
           err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_i4_2d)
         endif
         deallocate(nc%buf_i4_2d)
       elseif (ndims.eq.3) then
         if (nc%vars(i)%large) then
           do k = 1,dim3
             err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,k),cnt(:,k),&
                                        nc%buf_i4_3d(:,:,k),bsiz,mpi_int)
           enddo
         else
           err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_i4_3d)
         endif
         deallocate(nc%buf_i4_3d)
       elseif (ndims.eq.4) then
         if (nc%vars(i)%large) then
           do k = 1,dim4
             err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,k),cnt(:,k),&
                                        nc%buf_i4_4d(:,:,:,k),bsiz,mpi_int)
           enddo
         else
           err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_i4_4d)
         endif
         deallocate(nc%buf_i4_4d)
       endif
     elseif (vtype.eq.nf90_real) then
       if     (ndims.eq.1) then
         !err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_r4_1d)
         err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,1),cnt(:,1),      &
                                  nc%buf_r4_1d(:),bsiz,mpi_real)
         deallocate(nc%buf_r4_1d)
       elseif (ndims.eq.2) then
         if (nc%vars(i)%large) then
           do k = 1,dim2
             err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,k),cnt(:,k),&
                                        nc%buf_r4_2d(:,k),bsiz,mpi_real)
           enddo
         else
           err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_r4_2d)
         endif
         deallocate(nc%buf_r4_2d)
       elseif (ndims.eq.3) then
         if (nc%vars(i)%large) then
           do k = 1,dim3
             err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,k),cnt(:,k),&
                                        nc%buf_r4_3d(:,:,k),bsiz,mpi_real)
           enddo
         else
           err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_r4_3d)
         endif
         deallocate(nc%buf_r4_3d)
       elseif (ndims.eq.4) then
         if (nc%vars(i)%large) then
           do k = 1,dim4
             err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,k),cnt(:,k),&
                                        nc%buf_r4_4d(:,:,:,k),bsiz,mpi_real)
           enddo
         else
           err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_r4_4d)
         endif
         deallocate(nc%buf_r4_4d)
       endif
     elseif (vtype.eq.nf90_real8) then
       if     (ndims.eq.1) then
         !err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_r8_1d)
         err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,1),cnt(:,1),      &
                                  nc%buf_r8_1d(:),bsiz,mpi_real8)
         deallocate(nc%buf_r8_1d)
       elseif (ndims.eq.2) then
         if (nc%vars(i)%large) then
           do k = 1,dim2
             err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,k),cnt(:,k),&
                                        nc%buf_r8_2d(:,k),bsiz,mpi_real8)
           enddo
         else
           err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_r8_2d)
         endif
         deallocate(nc%buf_r8_2d)
       elseif (ndims.eq.3) then
         if (nc%vars(i)%large) then
           do k = 1,dim3
             err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,k),cnt(:,k),&
                                        nc%buf_r8_3d(:,:,k),bsiz,mpi_real8)
           enddo
         else
           err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_r8_3d)
         endif
         deallocate(nc%buf_r8_3d)
       elseif (ndims.eq.4) then
         if (nc%vars(i)%large) then
           do k = 1,dim4
             err = nfmpi_put_vara_all(nc%fou,nc%vars(i)%id,sta(:,k),cnt(:,k),&
                                        nc%buf_r8_4d(:,:,:,k),bsiz,mpi_real8)
           enddo
         else
           err = nf90mpi_put_var_all(nc%fou,nc%vars(i)%id,nc%buf_r8_4d)
         endif
         deallocate(nc%buf_r8_4d)
       endif
     endif
     if (nc%vars(i)%large) then
       deallocate(sta)
       deallocate(cnt)
     endif
   else ! not pnetcdf
     if     (vtype.eq.nf90_int) then
       if     (ndims.eq.1) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_i4_1d)
         deallocate(nc%buf_i4_1d)
       elseif (ndims.eq.2) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_i4_2d)
         deallocate(nc%buf_i4_2d)
       elseif (ndims.eq.3) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_i4_3d)
         deallocate(nc%buf_i4_3d)
       elseif (ndims.eq.4) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_i4_4d)
         deallocate(nc%buf_i4_4d)
       endif
     elseif (vtype.eq.nf90_real) then
       if     (ndims.eq.1) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_r4_1d)
         deallocate(nc%buf_r4_1d)
       elseif (ndims.eq.2) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_r4_2d)
         deallocate(nc%buf_r4_2d)
       elseif (ndims.eq.3) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_r4_3d)
         deallocate(nc%buf_r4_3d)
       elseif (ndims.eq.4) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_r4_4d)
         deallocate(nc%buf_r4_4d)
       endif
     elseif (vtype.eq.nf90_real8) then
       if     (ndims.eq.1) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_r8_1d)
         deallocate(nc%buf_r8_1d)
       elseif (ndims.eq.2) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_r8_2d)
         deallocate(nc%buf_r8_2d)
       elseif (ndims.eq.3) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_r8_3d)
         deallocate(nc%buf_r8_3d)
       elseif (ndims.eq.4) then
         err = nf90_put_var(nc%fou,nc%vars(i)%id,nc%buf_r8_4d)
         deallocate(nc%buf_r8_4d)
       endif
     endif
   endif
   call nc_check(err,'nc_put_var',.true.)
!
   end subroutine nc_put_var
#if 0
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_temp(nc,i)
   implicit none
   type(ncu_t), intent(inout) :: nc
   integer(i4), intent(in   ) :: i
! local variables
   integer(i4) :: err
!
   if (nc%tin.eq.4) then
   else
   endif
   call nc_check(err,'')
!
   end subroutine nc_temp
#endif
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_check(err,message,pnc)
   implicit none
   integer(i4),                intent(in   ) :: err
   character(len=*), optional, intent(in   ) :: message
   logical(l4),      optional, intent(in   ) :: pnc
! local variables
   logical(l4) :: is_pnc
!
   if (present(pnc)) then
     is_pnc = pnc
   else
     is_pnc = .false.
   endif
!
   if (is_pnc) then
     if (err.ne.nf90mpi_noerr) then
       if (present(message)) print*,trim(message)//': '//trim(nf90mpi_strerror(err))
       stop
     endif
   else
     if (err.ne.nf90_noerr) then
       if (present(message)) print*,trim(message)//': '//trim(nf90_strerror(err))
       stop
     endif
   endif
!
   end subroutine nc_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module netcdf_utility
!-------------------------------------------------------------------------------
