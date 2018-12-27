!-------------------------------------------------------------------------------
   module remap_matrix
!-------------------------------------------------------------------------------
!
!  abstract : 
!
!  history log :
!    2016-03-18  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!    2017-06-02  junghan kim    apply the remapping matrix of kim
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only: i4, r4, r8, l4
   use netcdf
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
! matrix sizes
   !    ne030:   48602
   !    ne060:  194402
   !    ne120:  777602
   !    ne240: 3110402
   !    ne320: 
   !  360x180:   64800
   !  720x330:  237600
   !  720x360:  259200
   ! 1024x768:  786432
   integer(i4), dimension(4), parameter :: cs_size  = (/48602,194402,777602,3110402/)
   integer(i4), dimension(4), parameter :: ne_def   = (/   30,    60,   120,    240/)
   integer(i4), dimension(3), parameter :: ll_size  = (/64800,259200,786432/)
   integer(i4), dimension(3), parameter :: nlon_def = (/  360,   720,  1024/)
   integer(i4), dimension(3), parameter :: nlat_def = (/  180,   360,   768/)
   character(len=16), parameter :: dpath = '/data/KIM3.0/'
   character(len=32), parameter :: ll_file = 'latlon.nc'
!
   type rm_t
     character(len=512) :: filename
     logical(l4)        :: is_scrip, has_llinfo         ! for khkim
     integer(i4)        :: src_size, dst_size
     integer(i4)        :: num_links
     character(len=16)  :: unit_lonlat
     integer(i4), dimension(:)  , allocatable :: src_adrs, dst_adrs
     real(r8)   , dimension(:,:), allocatable :: matrix
     real(r8)   , dimension(:)  , allocatable :: src_lats, src_lons, dst_lats, dst_lons
!
     logical(l4)        :: dst_ll
     integer(i4)        :: dst_nlat, dst_nlon
   end type rm_t
!
   interface rm_remap
     module procedure rm_remap_1d_i4
     module procedure rm_remap_2d_i4
     module procedure rm_remap_3d_i4
     module procedure rm_remap_4d_i4
     module procedure rm_remap_5d_i4
     module procedure rm_remap_1d_r4
     module procedure rm_remap_2d_r4
     module procedure rm_remap_3d_r4
     module procedure rm_remap_4d_r4
     module procedure rm_remap_5d_r4
     module procedure rm_remap_1d_r8
     module procedure rm_remap_2d_r8
     module procedure rm_remap_3d_r8
     module procedure rm_remap_4d_r8
     module procedure rm_remap_5d_r8
   end interface rm_remap
!
   public :: rm_t
   public :: rm_initialize, rm_finalize, rm_write_ll_info, rm_check, rm_rewrite, rm_remap
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_initialize(rmat, remap_file)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),       intent(  out) :: rmat
   character(len=*), intent(in   ) :: remap_file
! local variables
   integer(i4)        :: idx
   logical(l4)        :: has_src_info, has_dst_info
   logical(l4)        :: src_use_cs, dst_use_cs
   character(len=512) :: src_ll_file, dst_ll_file
!
   rmat%is_scrip = check_scrip_file(trim(remap_file))
!
   rmat%has_llinfo = .false.
!
   rmat%filename = trim(remap_file)
   if (rmat%is_scrip) then
     rmat%has_llinfo = .true.
     call read_scrip_matrix(trim(rmat%filename),rmat%src_size,rmat%dst_size,     &
                          rmat%num_links,rmat%src_adrs,rmat%dst_adrs,rmat%matrix)
     call read_scrip_latlon(trim(rmat%filename),rmat%src_size,rmat%dst_size,     &
                        rmat%src_lats,rmat%src_lons,rmat%dst_lats,rmat%dst_lons, &
                        rmat%unit_lonlat)
     rmat%dst_ll = .false.
   else
     call read_kim_matrix(trim(rmat%filename),rmat%src_size,rmat%dst_size,       &
                          rmat%num_links,rmat%src_adrs,rmat%dst_adrs,rmat%matrix)
     has_src_info = grid_info(rmat%src_size,src_ll_file,src_use_cs,idx)
     has_dst_info = grid_info(rmat%dst_size,dst_ll_file,dst_use_cs,idx)
     if (has_src_info.and.has_dst_info) then
       rmat%has_llinfo = .true.
       call read_kim_latlon(trim(rmat%filename),rmat%src_size,src_ll_file,src_use_cs,&
                                                rmat%dst_size,dst_ll_file,dst_use_cs,&
            rmat%src_lats,rmat%src_lons,rmat%dst_lats,rmat%dst_lons,rmat%unit_lonlat)
     endif
     if (has_dst_info.and..not.dst_use_cs) then
       rmat%dst_ll    = .true.
       rmat%dst_nlon  = nlon_def(idx)
       rmat%dst_nlat  = nlat_def(idx)
     else
       rmat%dst_ll = .false.
     endif
   endif
   print*,' '
   print*,' * loading remap file = '//trim(remap_file)
   print*,'   - matrix (min/max) = ',minval(rmat%matrix),maxval(rmat%matrix)
   print *, size(rmat%matrix(:,:))
!
   return
   end subroutine rm_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function check_scrip_file(remap_file) result(isscrip)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*), intent(in   ) :: remap_file
   logical(l4)                     :: isscrip
! local variables
   integer(i4) :: err, fid, did
!
   isscrip = .true.
   err = nf90_open(trim(remap_file),nf90_nowrite,fid)
   call check_nc(err,'cannot open file : '//trim(remap_file))
   err = nf90_inq_dimid(fid,'src_grid_size',did)
   if (err.ne.nf90_noerr) isscrip = .false.
   err = nf90_inq_dimid(fid,'dst_grid_size',did)
   if (err.ne.nf90_noerr) isscrip = .false.
   err = nf90_inq_dimid(fid,'num_links',did)
   if (err.ne.nf90_noerr) isscrip = .false.
   err = nf90_inq_dimid(fid,'num_wgts',did)
   if (err.ne.nf90_noerr) isscrip = .false.
   err = nf90_close(fid)
!
   return
   end function check_scrip_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_scrip_matrix(filename, src_size, dst_size, num_links, src_adrs, dst_adrs, matrix)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),                         intent(in   ) :: filename
   integer(i4)     ,                         intent(inout) :: src_size, dst_size, num_links
   integer(i4), dimension(:),   allocatable, intent(inout) :: src_adrs, dst_adrs
   real(r8),    dimension(:,:), allocatable, intent(inout) :: matrix
! local variables
   integer(i4) :: err, fid, did, vid
!
   err = nf90_open(trim(filename),nf90_nowrite,fid)
   call check_nc(err,'cannot open file : '//trim(filename))
   err = nf90_inq_dimid(fid,'src_grid_size',did)
   call check_nc(err,'check dimname : src_grid_size')
   err = nf90_inquire_dimension(fid,did,len=src_size)
   err = nf90_inq_dimid(fid,'dst_grid_size',did)
   call check_nc(err,'check dimname : dst_grid_size')
   err = nf90_inquire_dimension(fid,did,len=dst_size)
   err = nf90_inq_dimid(fid,'num_links',did)
   call check_nc(err,'check dimname : num_links')
   err = nf90_inquire_dimension(fid,did,len=num_links)
! rmat info.
   allocate(src_adrs(num_links))
   allocate(dst_adrs(num_links))
   allocate(matrix(3,num_links))
   matrix = 0.0_r8
! rmat info.
   err = nf90_inq_varid(fid,'src_address',vid)
   call check_nc(err,'check varname : src_address')
   err = nf90_get_var(fid,vid,src_adrs)
   err = nf90_inq_varid(fid,'dst_address',vid)
   call check_nc(err,'check varname : dst_address')
   err = nf90_get_var(fid,vid,dst_adrs)
   err = nf90_inq_varid(fid,'remap_matrix',vid)
   call check_nc(err,'check varname : remap_matrix')
   err = nf90_get_var(fid,vid,matrix)
!
   err = nf90_close(fid)
!
   return
   end subroutine read_scrip_matrix
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_scrip_latlon(filename, src_size, dst_size, src_lats, src_lons, dst_lats, dst_lons, unit_lonlat)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),                    intent(in   ) :: filename
   integer(i4)     ,                    intent(in   ) :: src_size, dst_size
   real(r8), dimension(:), allocatable, intent(inout) :: src_lats, src_lons
   real(r8), dimension(:), allocatable, intent(inout) :: dst_lats, dst_lons
   character(len=*),                    intent(inout) :: unit_lonlat
! local variables
   integer(i4) :: err, fid, did, vid
!
! lats,lons
   allocate(src_lats(src_size))
   allocate(src_lons(src_size))
   allocate(dst_lats(dst_size))
   allocate(dst_lons(dst_size))
!
   err = nf90_open(trim(filename),nf90_nowrite,fid)
   call check_nc(err,'cannot open file : '//trim(filename))
! lat,lon (src)
   err = nf90_inq_varid(fid,'src_grid_center_lat',vid)
   call check_nc(err,'check varname : src_grid_center_lat')
   err = nf90_get_var(fid,vid,src_lats)
   err = nf90_inq_varid(fid,'src_grid_center_lon',vid)
   call check_nc(err,'check varname : src_grid_center_lon')
   err = nf90_get_var(fid,vid,src_lons)
   err = nf90_get_att(fid,vid,'units',unit_lonlat)
   call check_nc(err,'check attribute : units')
! lat,lon (dst)
   err = nf90_inq_varid(fid,'dst_grid_center_lat',vid)
   call check_nc(err,'check varname : dst_grid_center_lat')
   err = nf90_get_var(fid,vid,dst_lats)
   err = nf90_inq_varid(fid,'dst_grid_center_lon',vid)
   call check_nc(err,'check varname : dst_grid_center_lon')
   err = nf90_get_var(fid,vid,dst_lons)
!
   err = nf90_close(fid)
!
   return
   end subroutine read_scrip_latlon
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_kim_matrix(filename, src_size, dst_size, num_links, src_adrs, dst_adrs, matrix)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),                         intent(in   ) :: filename
   integer(i4)     ,                         intent(inout) :: src_size, dst_size, num_links
   integer(i4), dimension(:),   allocatable, intent(inout) :: src_adrs, dst_adrs
   real(r8),    dimension(:,:), allocatable, intent(inout) :: matrix
! local variables
   integer(i4)        :: i, d, err, fid, sid, did, vid
   !-khkim
   integer(i4) :: mat_size
   integer(i4), allocatable :: src_adrs_kim(:,:)
   real(r8)   , allocatable :: matrix_kim(:,:)
!
   err = nf90_open(trim(filename),nf90_nowrite,fid)
   call check_nc(err,'cannot open file : '//trim(filename))
   ! dim
   err = nf90_inq_dimid(fid,'src_size',sid)
   call check_nc(err,'check dimname : src_size')
   err = nf90_inquire_dimension(fid,sid,len=src_size)
   call check_nc(err,'check dimension : src_size')
   err = nf90_inq_dimid(fid,'dst_size',did)
   call check_nc(err,'check dimname : dst_size')
   err = nf90_inquire_dimension(fid,did,len=dst_size)
   call check_nc(err,'check dimension : dst_size')
   err = nf90_inq_dimid(fid,'mat_size',did)
   call check_nc(err,'check dimname : mat_size')
   err = nf90_inquire_dimension(fid,did,len=mat_size)
   call check_nc(err,'check dimension : mat_size')
   !
   allocate(src_adrs_kim(mat_size,dst_size))
   allocate(matrix_kim(mat_size,dst_size))
   matrix_kim = 0.0_r8
   ! var
   err = nf90_inq_varid(fid,'src_address',vid)
   call check_nc(err,'check varname : src_address')
   err = nf90_get_var(fid,vid,src_adrs_kim)
   call check_nc(err,'check variable : src_address')
   err = nf90_inq_varid(fid,'remap_matrix',vid)
   call check_nc(err,'check varname : remap_matrix')
   err = nf90_get_var(fid,vid, matrix_kim)
   call check_nc(err,'check variable : remap_matrix')
!
   err = nf90_close(fid)
!
   num_links = dst_size*mat_size
   ! rmat info.
   allocate(src_adrs(num_links))
   allocate(dst_adrs(num_links))
   allocate(matrix(3,num_links))
   matrix = 0.0_r8
   src_adrs(:) = reshape(src_adrs_kim+1,(/num_links/))
   do d = 1,dst_size
     do i = 1,mat_size
       dst_adrs((d-1)*mat_size+i) = d
     enddo
   enddo
   matrix(1,:) = reshape(matrix_kim,(/num_links/))
   deallocate(src_adrs_kim)
   deallocate(matrix_kim)
!
   return
   end subroutine read_kim_matrix
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_kim_latlon(filename, src_size, src_ll_file, src_use_cs, dst_size, dst_ll_file, dst_use_cs, src_lats, src_lons, dst_lats, dst_lons, unit_lonlat)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),                    intent(in   ) :: filename
   integer(i4)     ,                    intent(in   ) :: src_size, dst_size
   character(len=*),                    intent(in   ) :: src_ll_file, dst_ll_file
   logical(l4)     ,                    intent(in   ) :: src_use_cs, dst_use_cs
   real(r8), dimension(:), allocatable, intent(inout) :: src_lats, src_lons
   real(r8), dimension(:), allocatable, intent(inout) :: dst_lats, dst_lons
   character(len=*),                    intent(inout) :: unit_lonlat
! local variables
!   logical(l4)        :: has_src_info, has_dst_info
!
   allocate(src_lats(src_size))
   allocate(src_lons(src_size))
   allocate(dst_lats(dst_size))
   allocate(dst_lons(dst_size))
   call get_latlon_kim(trim(src_ll_file),src_use_cs,src_lats,src_lons)
   call get_latlon_kim(trim(dst_ll_file),dst_use_cs,dst_lats,dst_lons)
!
   return
   end subroutine read_kim_latlon
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function grid_info(gsize, ll_file, use_cs, idx) result(has_grid)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4)     ,               intent(in   ) :: gsize
   character(len=*),               intent(  out) :: ll_file
   logical(l4)     ,               intent(  out) :: use_cs
   integer(i4)     ,               intent(  out) :: idx
   logical(l4)                                   :: has_grid
! local variables
!
! matrix sizes
   !    ne030:   48602
   !    ne060:  194402
   !    ne120:  777602
   !    ne240: 3110402
   !    ne320: 
   !  360x180:   64800
   !  720x330:  237600
   !  720x360:  259200
   ! 1024x768:  786432
   has_grid = .false.
   if      (minval(abs(cs_size(:)-gsize)).eq.0) then
     has_grid = .true.
     use_cs   = .true.
     idx = minloc(abs(cs_size(:)-gsize),dim=1)
     write(ll_file,'(a,a,i3.3,a)') trim(dpath),'cs_grid/cs_grid_ne',ne_def(idx),'np4_rotated.nc'
   else if (minval(abs(ll_size(:)-gsize)).eq.0) then
     has_grid = .true.
     use_cs   = .false.
     idx = minloc(abs(ll_size(:)-gsize),dim=1)
     if (idx.ge.3) then
       write(ll_file,'(a,a,i4.4,a,i3.3,a)') trim(dpath),                       &
                      'll_grid/ll_coord_',nlon_def(idx),'x',nlat_def(idx),'.nc'
     else
       write(ll_file,'(a,a,i3.3,a,i3.3,a)') trim(dpath),                       &
                      'll_grid/ll_coord_',nlon_def(idx),'x',nlat_def(idx),'.nc'
     endif
   else
     has_grid = .false.
     use_cs   = .false.
   endif
!
   return
   end function grid_info
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_latlon_kim(ll_file, use_cs, lats, lons)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),               intent(in   ) :: ll_file
   logical(l4)     ,               intent(in   ) :: use_cs
   real(r8)        , dimension(:), intent(inout) :: lats, lons
! local variables
   integer(i4) :: i, j, idx, err, fid, did, vid
   integer(i4) :: lsize, nlat, nlon
   real(r8)   , allocatable :: latlons(:,:), tlats(:), tlons(:)
!
   err = nf90_open(trim(ll_file),nf90_nowrite,fid)
   call check_nc(err,'cannot open file : '//trim(ll_file))
!
   if (use_cs) then
     ! dim
     err = nf90_inq_dimid(fid,'up_size',did)
     call check_nc(err,'check dimname : up_size')
     err = nf90_inquire_dimension(fid,did,len=lsize)
     call check_nc(err,'check dimension : dst_size')
     !
     allocate(latlons(2,lsize))
     ! var
     err = nf90_inq_varid(fid,'latlons',vid)
     call check_nc(err,'check varname : latlons')
     err = nf90_get_var(fid,vid,latlons)
     call check_nc(err,'check variable : latlons')
     do i = 1,lsize
       lons(i) = latlons(2,i)
       lats(i) = latlons(1,i)
     enddo
     deallocate(latlons)
   else
     ! dim
     err = nf90_inq_dimid(fid,'nlon',did)
     call check_nc(err,'check dimname : nlon')
     err = nf90_inquire_dimension(fid,did,len=nlon)
     call check_nc(err,'check dimension : nlon')
     err = nf90_inq_dimid(fid,'nlat',did)
     call check_nc(err,'check dimname : nlat')
     err = nf90_inquire_dimension(fid,did,len=nlat)
     call check_nc(err,'check dimension : nlat')
     !
     allocate(tlons(nlon),tlats(nlat))
     !
     ! var
     err = nf90_inq_varid(fid,'lons',vid)
     call check_nc(err,'check varname : lons')
     err = nf90_get_var(fid,vid,tlons)
     call check_nc(err,'check variable : lons')
     err = nf90_inq_varid(fid,'lats',vid)
     call check_nc(err,'check varname : lats')
     err = nf90_get_var(fid,vid,tlats)
     call check_nc(err,'check variable : lats')
     !
     do j = 1,nlat
       do i = 1,nlon
         idx = (j-1)*nlon+i
         lons(idx) = tlons(i)
         lats(idx) = tlats(j)
       enddo
     enddo
     !
     deallocate(tlons,tlats)
   endif
!
   err = nf90_close(fid)
!
!
   return
   end subroutine get_latlon_kim
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_finalize(rmat)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t), intent(inout) :: rmat
! rmat info.
   deallocate(rmat%src_adrs)
   deallocate(rmat%dst_adrs)
   deallocate(rmat%matrix)
   if (rmat%has_llinfo) then
     ! lats, lons
     deallocate(rmat%src_lats)
     deallocate(rmat%src_lons)
     deallocate(rmat%dst_lats)
     deallocate(rmat%dst_lons)
   endif
!
   return
   end subroutine rm_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_check(rmat, error)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),         intent(in   ) :: rmat
   real(r8), optional, intent(in   ) :: error
! local variables
   integer(i4) :: i, sadr, dadr, idst
   integer(i4) :: nthrs, ith
   real(r8)    :: lerror
   real(r8), dimension(:), allocatable :: weights
!
   if (present(error)) then
     lerror = error
   else
     lerror = 1.0d-13
   endif
!
   allocate(weights(rmat%dst_size))
!
   print*,''
   print*,'File = ',trim(rmat%filename)
   print*,'total links = ',rmat%num_links
   print*,' inum,src,dst,weight'
!
   weights(:) = 0.0_r8
   ith = 0
   do i = 1,rmat%num_links
     sadr = rmat%src_adrs(i)
     dadr = rmat%dst_adrs(i)
     weights(dadr) = weights(dadr)+rmat%matrix(1,i)
     if (rmat%matrix(1,i).lt.0.00_r8.or.rmat%matrix(1,i).gt.1.00_r8) then
       ith = ith+1
       print *,ith,sadr,dadr,rmat%matrix(1,i) !,rmat%matrix(2,i),rmat%matrix(3,i)
     endif
   enddo
!
   nthrs = 0
   do idst = 1,rmat%dst_size
     if ((weights(idst).gt.1.0_r8+lerror).or.(weights(idst).lt.0.0_r8-lerror)) then
       print *,'check weight : ',idst,weights(idst)
       nthrs = nthrs+1
     endif
   enddo
!
   write(*,'(i9,a,i9,a,e10.2,a)') nthrs,' points in ',rmat%dst_size,          &
                                             ' : greater than error(',error,')'
!
   deallocate(weights)
!
   return
   end subroutine rm_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_1d_i4(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),                intent(in   ) :: rm
   integer(i4), dimension(:), intent(in   ) :: src
   integer(i4), dimension(:), intent(  out) :: dst
! local variables
   integer(i4) :: i, sadr, dadr
   real(r8), allocatable :: tmp(:)
!
   if (rm%src_size.ne.size(src)) then
     print*,'error: rm_remap_1d(check src size...)'
     stop
   endif
   if (rm%dst_size.ne.size(dst)) then
     print*,'error: rm_remap_1d(check dst size...)'
     stop
   endif
!
   allocate(tmp(rm%dst_size))
   dst(:) = 0_i4
   tmp(:) = 0.0_r8
   do i = 1,rm%num_links
     sadr = rm%src_adrs(i)
     dadr = rm%dst_adrs(i)
     tmp(dadr) = tmp(dadr)+rm%matrix(1,i)*src(sadr)
   enddo
   dst(:) = tmp(:)
   deallocate(tmp)
!
   return
   end subroutine rm_remap_1d_i4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_2d_i4(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),                  intent(in   ) :: rm
   integer(i4), dimension(:,:), intent(in   ) :: src
   integer(i4), dimension(:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=2)
   if (size(dst,dim=2).ne.n) then
     print*,'error: rm_remap_2d_i4(check n2: src, dst size...)',size(dst,dim=2),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_1d_i4(rm,src(:,i),dst(:,i))
   enddo
!
   return
   end subroutine rm_remap_2d_i4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_3d_i4(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),                    intent(in   ) :: rm
   integer(i4), dimension(:,:,:), intent(in   ) :: src
   integer(i4), dimension(:,:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=3)
   if (size(dst,dim=3).ne.n) then
     print*,'error: rm_remap_3d_i4(check n: src, dst size...)',size(dst,dim=3),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_2d_i4(rm,src(:,:,i),dst(:,:,i))
   enddo
!
   return
   end subroutine rm_remap_3d_i4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_4d_i4(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),                      intent(in   ) :: rm
   integer(i4), dimension(:,:,:,:), intent(in   ) :: src
   integer(i4), dimension(:,:,:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=4)
   if (size(dst,dim=4).ne.n) then
     print*,'error: rm_remap_4d_i4(check n: src, dst size...)',size(dst,dim=4),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_3d_i4(rm,src(:,:,:,i),dst(:,:,:,i))
   enddo
!
   return
   end subroutine rm_remap_4d_i4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_5d_i4(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),                        intent(in   ) :: rm
   integer(i4), dimension(:,:,:,:,:), intent(in   ) :: src
   integer(i4), dimension(:,:,:,:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=5)
   if (size(dst,dim=5).ne.n) then
     print*,'error: rm_remap_5d_i4(check n: src, dst size...)',size(dst,dim=5),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_4d_i4(rm,src(:,:,:,:,i),dst(:,:,:,:,i))
   enddo
!
   return
   end subroutine rm_remap_5d_i4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_1d_r4(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),             intent(in   ) :: rm
   real(r4), dimension(:), intent(in   ) :: src
   real(r4), dimension(:), intent(  out) :: dst
! local variables
   integer(i4) :: i, sadr, dadr
!
   if (rm%src_size.ne.size(src)) then
     print*,'error: rm_remap_1d(check src size...)'
     stop
   endif
   if (rm%dst_size.ne.size(dst)) then
     print*,'error: rm_remap_1d(check dst size...)'
     stop
   endif
!
   dst(:) = 0.0_r4
   do i = 1,rm%num_links
     sadr = rm%src_adrs(i)
     dadr = rm%dst_adrs(i)
     dst(dadr) = dst(dadr)+rm%matrix(1,i)*src(sadr)
   enddo
!
   return
   end subroutine rm_remap_1d_r4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_2d_r4(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),               intent(in   ) :: rm
   real(r4), dimension(:,:), intent(in   ) :: src
   real(r4), dimension(:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=2)
   if (size(dst,dim=2).ne.n) then
     print*,'error: rm_remap_2d_r4(check n2: src, dst size...)',size(dst,dim=2),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_1d_r4(rm,src(:,i),dst(:,i))
   enddo
!
   return
   end subroutine rm_remap_2d_r4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_3d_r4(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),                 intent(in   ) :: rm
   real(r4), dimension(:,:,:), intent(in   ) :: src
   real(r4), dimension(:,:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=3)
   if (size(dst,dim=3).ne.n) then
     print*,'error: rm_remap_3d_r4(check n: src, dst size...)',size(dst,dim=3),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_2d_r4(rm,src(:,:,i),dst(:,:,i))
   enddo
!
   return
   end subroutine rm_remap_3d_r4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_4d_r4(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),                   intent(in   ) :: rm
   real(r4), dimension(:,:,:,:), intent(in   ) :: src
   real(r4), dimension(:,:,:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=4)
   if (size(dst,dim=4).ne.n) then
     print*,'error: rm_remap_4d_r4(check n: src, dst size...)',size(dst,dim=4),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_3d_r4(rm,src(:,:,:,i),dst(:,:,:,i))
   enddo
!
   return
   end subroutine rm_remap_4d_r4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_5d_r4(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),                     intent(in   ) :: rm
   real(r4), dimension(:,:,:,:,:), intent(in   ) :: src
   real(r4), dimension(:,:,:,:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=5)
   if (size(dst,dim=5).ne.n) then
     print*,'error: rm_remap_5d_r4(check n: src, dst size...)',size(dst,dim=5),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_4d_r4(rm,src(:,:,:,:,i),dst(:,:,:,:,i))
   enddo
!
   return
   end subroutine rm_remap_5d_r4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_1d_r8(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),             intent(in   ) :: rm
   real(r8), dimension(:), intent(in   ) :: src
   real(r8), dimension(:), intent(  out) :: dst
! local variables
   integer(i4) :: i, sadr, dadr
!
   if (rm%src_size.ne.size(src)) then
     print*,'error: rm_remap_1d(check src size...)'
     stop
   endif
   if (rm%dst_size.ne.size(dst)) then
     print*,'error: rm_remap_1d(check dst size...)'
     stop
   endif
!
   dst(:) = 0.0_r8
   do i = 1,rm%num_links
     sadr = rm%src_adrs(i)
     dadr = rm%dst_adrs(i)
     !if (sadr.gt.48602) print *, 'DEBUG: ',sadr
     dst(dadr) = dst(dadr)+rm%matrix(1,i)*src(sadr)
   enddo
!
   return
   end subroutine rm_remap_1d_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_2d_r8(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),               intent(in   ) :: rm
   real(r8), dimension(:,:), intent(in   ) :: src
   real(r8), dimension(:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=2)
   if (size(dst,dim=2).ne.n) then
     print*,'error: rm_remap_2d_r8(check n2: src, dst size...)',size(dst,dim=2),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_1d_r8(rm,src(:,i),dst(:,i))
   enddo
!
   return
   end subroutine rm_remap_2d_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_3d_r8(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),                 intent(in   ) :: rm
   real(r8), dimension(:,:,:), intent(in   ) :: src
   real(r8), dimension(:,:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=3)
   if (size(dst,dim=3).ne.n) then
     print*,'error: rm_remap_3d_r8(check n: src, dst size...)',size(dst,dim=3),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_2d_r8(rm,src(:,:,i),dst(:,:,i))
   enddo
!
   return
   end subroutine rm_remap_3d_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_4d_r8(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),                   intent(in   ) :: rm
   real(r8), dimension(:,:,:,:), intent(in   ) :: src
   real(r8), dimension(:,:,:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=4)
   if (size(dst,dim=4).ne.n) then
     print*,'error: rm_remap_4d_r8(check n: src, dst size...)',size(dst,dim=4),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_3d_r8(rm,src(:,:,:,i),dst(:,:,:,i))
   enddo
!
   return
   end subroutine rm_remap_4d_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_remap_5d_r8(rm, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),                     intent(in   ) :: rm
   real(r8), dimension(:,:,:,:,:), intent(in   ) :: src
   real(r8), dimension(:,:,:,:,:), intent(  out) :: dst
! local variables
   integer(i4) :: i, n
!
   n = size(src,dim=5)
   if (size(dst,dim=5).ne.n) then
     print*,'error: rm_remap_5d_r8(check n: src, dst size...)',size(dst,dim=5),n
     stop
   endif
!
   do i = 1,n
     call rm_remap_4d_r8(rm,src(:,:,:,:,i),dst(:,:,:,:,i))
   enddo
!
   return
   end subroutine rm_remap_5d_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_rewrite(rmat)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t), intent(inout) :: rmat
! local variables
   integer(i4) :: num_links, src_size, dst_size
   integer(i4) :: err, fid, vid
!
   err = nf90_open(trim(rmat%filename),nf90_write,fid)
   call check_nc(err,'cannot open file : '//trim(rmat%filename))
!
   err = nf90_inq_varid(fid,'remap_matrix',vid)
   call check_nc(err,'check varname in rewrite : remap_matrix')
   err = nf90_put_var(fid,vid,rmat%matrix)
   call check_nc(err,'check put_var in rewrite')
!
   err = nf90_close(fid)
!
   return
   end subroutine rm_rewrite
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rm_write_ll_info(rmat, oufile)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t),       intent(inout) :: rmat
   character(len=*), intent(in   ) :: oufile
! local variables
   integer(i4) :: wmod, info, err
   integer(i4) :: fid
   integer(i4) :: src_did, src_lat_id, src_lon_id
   integer(i4) :: dst_did, dst_lat_id, dst_lon_id
!
! file
   wmod = ior(nf90_clobber,nf90_64bit_offset)
   err = nf90_create(trim(oufile),wmod,fid)
   call check_nc(err,'rm_write_ll_info: 1')
!
! dimensions
!PRINT *, 'DEBUG: ', rmat%src_size, rmat%dst_size
   err = nf90_def_dim(fid,'src_size',rmat%src_size,src_did)
   call check_nc(err,'rm_write_ll_info: 2')
   err = nf90_def_dim(fid,'dst_size',rmat%dst_size,dst_did)
   call check_nc(err,'rm_write_ll_info: 3')
!
! variables
   err = nf90_def_var(fid,'src_lats',nf90_real8,src_did,src_lat_id)
   call check_nc(err,'rm_write_ll_info: 4')
   err = nf90_put_att(fid,src_lat_id,'units',trim(rmat%unit_lonlat))
   call check_nc(err,'rm_write_ll_info: 5')
   err = nf90_def_var(fid,'src_lons',nf90_real8,src_did,src_lon_id)
   call check_nc(err,'rm_write_ll_info: 6')
   err = nf90_put_att(fid,src_lon_id,'units',trim(rmat%unit_lonlat))
   call check_nc(err,'rm_write_ll_info: 7')
   err = nf90_def_var(fid,'dst_lats',nf90_real8,dst_did,dst_lat_id)
   call check_nc(err,'rm_write_ll_info: 8')
   err = nf90_put_att(fid,dst_lat_id,'units',trim(rmat%unit_lonlat))
   call check_nc(err,'rm_write_ll_info: 9')
   err = nf90_def_var(fid,'dst_lons',nf90_real8,dst_did,dst_lon_id)
   call check_nc(err,'rm_write_ll_info: 10')
   err = nf90_put_att(fid,dst_lon_id,'units',trim(rmat%unit_lonlat))
   call check_nc(err,'rm_write_ll_info: 11')
!
   err = nf90_enddef(fid)
!
! write
   err = nf90_put_var(fid,src_lat_id,rmat%src_lats)
   call check_nc(err,'rm_write_ll_info: 12')
   err = nf90_put_var(fid,src_lon_id,rmat%src_lons)
   call check_nc(err,'rm_write_ll_info: 13')
   err = nf90_put_var(fid,dst_lat_id,rmat%dst_lats)
   call check_nc(err,'rm_write_ll_info: 14')
   err = nf90_put_var(fid,dst_lon_id,rmat%dst_lons)
   call check_nc(err,'rm_write_ll_info: 15')
!
   err = nf90_close(fid)
   call check_nc(err,'rm_write_ll_info: 16')
!
   return
   end subroutine rm_write_ll_info
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_nc(err, message)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),      intent(in   ) :: err
   character(len=*), intent(in   ) :: message
!
   if (err.ne.0) then
     print*, trim(message)
     stop
   endif
!
   return
   end subroutine check_nc
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module remap_matrix
