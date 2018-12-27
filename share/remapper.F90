!-------------------------------------------------------------------------------
   module remapper
!-------------------------------------------------------------------------------
!
!  abstract : 
!
!  history log :
!    2016-03-18  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!    2017-06-02  junghan kim    apply the remapping matrix of kim
!    2017-08-29  junghan kim    lat-lon grid
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use netcdf
   use kinds,        only: i4, l4, r8
   use buffer,       only: buffer_t, buffer_initialize, buffer_finalize,       &
                           buffer_min, buffer_max
   use remap_matrix, only: rm_t, rm_initialize, rm_finalize, rm_write_ll_info, rm_remap
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   character(len=32), parameter :: ll_file = 'latlon.nc'
!
   integer, parameter :: max_ndims   = 5
   integer, parameter :: max_strname = 32
!
   type dim_t
     character(len=max_strname) :: name         ! ini
     integer(i4)                :: id           ! netcdf
     integer(i4)                :: size         ! netcdf
   end type dim_t
!
   type var_t
     character(len=max_strname) :: name                ! ini
     integer(i4) :: sndims, dndims
     type(dim_t) :: sdims(max_ndims), ddims(max_ndims)
     integer(i4) :: sid, did                           ! netcdf
     integer(i4) :: xtype
     logical(l4) :: isremap
   end type var_t
!
   type remap_t
     type(rm_t)         :: mat
     character(len=512) :: infile, oufile
     integer(i4)        :: sid, did ! source and destination netcdf id
     integer(i4)        :: sndims, dndims
     integer(i4)        :: nvars
     integer(i4)        :: lon_id, lat_id
     logical(l4)        :: src_ll, dst_ll
     type(dim_t)   , dimension(:), allocatable :: sdims, ddims
     type(var_t)   , dimension(:), allocatable :: vars
     type(buffer_t), dimension(:), allocatable :: src
     type(buffer_t), dimension(:), allocatable :: dst
   end type remap_t
!
   public :: remap_t
   public :: remap_initialize, remap_finalize, remap_driver, write_ll_info
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine remap_initialize(remap, remapfile, infile, oufile,               &
                   ndims, dimnames, nvars, varnames, ndims_vars, dimnames_vars)
!-------------------------------------------------------------------------------
   implicit none
!
   type(remap_t)   ,                   intent(inout) :: remap
   character(len=*),                   intent(in   ) :: remapfile
   character(len=*),                   intent(in   ) :: infile
   character(len=*),                   intent(in   ) :: oufile
   integer(i4)     ,                   intent(in   ) :: ndims
   character(len=*), dimension(ndims), intent(in   ) :: dimnames
   integer(i4)     ,                   intent(in   ) :: nvars
   character(len=*), dimension(nvars), intent(in   ) :: varnames
   integer(i4)     , dimension(nvars), intent(in   ) :: ndims_vars
   character(len=*), dimension(:,:)  , intent(in   ) :: dimnames_vars
! local variables
   integer(i4) :: err, dd, vv, ix
   integer(i4), dimension(:), allocatable :: ids
!
   call rm_initialize(remap%mat,trim(remapfile))
   remap%dst_ll = remap%mat%dst_ll
!
   remap%infile = trim(infile)
   remap%oufile = trim(oufile)
!
!===================================================================================
! INPUT
!===================================================================================
! dimensions
   remap%sndims = ndims
   allocate(remap%sdims(ndims))
   allocate(remap%ddims(ndims+1))
   !
   err = nf90_open(trim(infile),nf90_nowrite,remap%sid)
   call check_nc(err,'cannot open file : '//trim(infile))
   ! src (name, size, id)
   do dd = 1,ndims
     remap%sdims(dd)%name = trim(dimnames(dd))
     err = nf90_inq_dimid(remap%sid,trim(dimnames(dd)),remap%sdims(dd)%id)
     call check_nc(err,'check dimension name : '//trim(dimnames(dd)))
     err = nf90_inquire_dimension(remap%sid,remap%sdims(dd)%id,len=remap%sdims(dd)%size)
     call check_nc(err,'cannot read dimsize ')
   enddo
   ! dst (ndims, name, size)
   if (remap%sdims(1)%size.eq.remap%mat%src_size) then
     remap%src_ll = .false.
     if (remap%dst_ll) then ! START
       remap%dndims = ndims+1
       remap%ddims(1)%size = remap%mat%dst_nlon
       remap%ddims(1)%name = 'nlon'
       remap%ddims(2)%size = remap%mat%dst_nlat
       remap%ddims(2)%name = 'nlat'
       do dd = 2,ndims
         remap%ddims(dd+1)%size = remap%sdims(dd)%size
         remap%ddims(dd+1)%name = trim(dimnames(dd))
       enddo
     else
       remap%dndims = ndims
       remap%ddims(1)%size = remap%mat%dst_size
       remap%ddims(1)%name = trim(dimnames(1))
       do dd = 2,ndims
         remap%ddims(dd)%size = remap%sdims(dd)%size
         remap%ddims(dd)%name = trim(dimnames(dd))
       enddo
     endif
   elseif (remap%sdims(1)%size*remap%sdims(2)%size.eq.remap%mat%src_size) then
     remap%src_ll = .true.
     if (remap%dst_ll) then
       remap%dndims = ndims
       remap%ddims(1)%size = remap%mat%dst_nlon
       remap%ddims(1)%name = 'nlon'
       remap%ddims(2)%size = remap%mat%dst_nlat
       remap%ddims(2)%name = 'nlat'
       do dd = 3,ndims
         remap%ddims(dd)%size = remap%sdims(dd)%size
         remap%ddims(dd)%name = trim(dimnames(dd))
       enddo
     else
       remap%dndims = ndims-1
       remap%ddims(1)%size = remap%mat%dst_size
       remap%ddims(1)%name = 'ncol'
       do dd = 3,ndims
         remap%ddims(dd-1)%size = remap%sdims(dd)%size
         remap%ddims(dd-1)%name = trim(dimnames(dd))
       enddo
     endif
   else
     print *, 'check dimension for remap...'
     stop
   endif
!
! variables
   remap%nvars = nvars
   allocate(remap%src(nvars))
   allocate(remap%dst(nvars))
   allocate(remap%vars(nvars))
   !
   ! dimensions of var
   do vv = 1,nvars
     remap%vars(vv)%name = trim(varnames(vv))
     ! src (sndims, name, size)
     remap%vars(vv)%sndims = ndims_vars(vv)
     do dd = 1,ndims_vars(vv)
       ix = get_dimix(remap%sdims,trim(dimnames_vars(dd,vv)))
       remap%vars(vv)%sdims(dd)%name = dimnames_vars(dd,vv)
       remap%vars(vv)%sdims(dd)%size = remap%sdims(ix)%size
     enddo
     remap%vars(vv)%isremap = .false.
     do dd = 1,ndims_vars(vv)
       ix = get_dimix(remap%sdims,trim(dimnames_vars(dd,vv)))
       if (dd.eq.1.and.ix.eq.1) then
         remap%vars(vv)%isremap = .true.
       endif
     enddo
     ! dst (dndims, name, size)
     if (remap%vars(vv)%isremap) then
       if (.not.remap%src_ll) then
         if (.not.remap%dst_ll) then
           remap%vars(vv)%dndims = ndims_vars(vv)
           remap%vars(vv)%ddims(1)%name = remap%ddims(1)%name
           remap%vars(vv)%ddims(1)%size = remap%ddims(1)%size
           do dd = 2,ndims_vars(vv)
             ix = get_dimix(remap%sdims,trim(dimnames_vars(dd,vv)))
             remap%vars(vv)%ddims(dd)%name = dimnames_vars(dd,vv)
             remap%vars(vv)%ddims(dd)%size = remap%sdims(ix)%size
           enddo
         else
           remap%vars(vv)%dndims = ndims_vars(vv)+1
           remap%vars(vv)%ddims(1)%name = remap%ddims(1)%name
           remap%vars(vv)%ddims(1)%size = remap%ddims(1)%size
           remap%vars(vv)%ddims(2)%name = remap%ddims(2)%name
           remap%vars(vv)%ddims(2)%size = remap%ddims(2)%size
           do dd = 2,ndims_vars(vv)
             ix = get_dimix(remap%sdims,trim(dimnames_vars(dd,vv)))
             remap%vars(vv)%ddims(dd+1)%name = dimnames_vars(dd,vv)
             remap%vars(vv)%ddims(dd+1)%size = remap%sdims(ix)%size
           enddo
         endif
       else
         remap%vars(vv)%dndims = ndims_vars(vv)-1
         remap%vars(vv)%ddims(1)%name = 'ncol'
         remap%vars(vv)%ddims(1)%size = remap%ddims(1)%size
         do dd = 3,ndims_vars(vv)
           ix = get_dimix(remap%ddims,trim(dimnames_vars(dd,vv)))
           remap%vars(vv)%ddims(dd-1)%name = dimnames_vars(dd,vv)
           remap%vars(vv)%ddims(dd-1)%size = remap%ddims(ix)%size
         enddo
       endif
     else
       remap%vars(vv)%dndims = ndims_vars(vv)
       do dd = 1,ndims_vars(vv)
         ix = get_dimix(remap%sdims,trim(dimnames_vars(dd,vv)))
         remap%vars(vv)%ddims(dd)%name = dimnames_vars(dd,vv)
         remap%vars(vv)%ddims(dd)%size = remap%ddims(ix)%size
       enddo
     endif
   enddo
!
   ! read variable info
   do vv = 1,nvars
     err = nf90_inq_varid(remap%sid,trim(remap%vars(vv)%name),remap%vars(vv)%sid)
     call check_nc(err,'check variable name : '//trim(remap%vars(vv)%name))
     err = nf90_inquire_variable(remap%sid,remap%vars(vv)%sid,xtype=remap%vars(vv)%xtype)
     call check_nc(err,'cannot read xtype ')
     remap%vars(vv)%xtype   = nf90_double
   enddo
!
!===================================================================================
! OUTPUT
!===================================================================================
!
   err = nf90_create(trim(oufile),ior(nf90_clobber,nf90_64bit_offset),remap%did)
   call check_nc(err,'cannot create file : '//trim(oufile))
!
! define dimenisons
   do dd = 1,remap%dndims
     err = nf90_def_dim(remap%did,trim(remap%ddims(dd)%name),remap%ddims(dd)%size,remap%ddims(dd)%id)
     call check_nc(err,'cannot define dimension : '//trim(remap%ddims(dd)%name))
   enddo
!
! define variables
   do vv = 1,nvars
     allocate(ids(remap%vars(vv)%dndims))
     do dd = 1,remap%vars(vv)%dndims
       ids(dd) = get_dimid(remap%ddims,trim(remap%vars(vv)%ddims(dd)%name))
       ix      = get_dimix(remap%ddims,trim(remap%vars(vv)%ddims(dd)%name))
     enddo
     err = nf90_def_var(remap%did,trim(remap%vars(vv)%name),remap%vars(vv)%xtype,&
                        ids(1:remap%vars(vv)%dndims),remap%vars(vv)%did)
     call check_nc(err,'cannot define variable : '//trim(remap%vars(vv)%name))
     deallocate(ids)
   enddo
!
   err = nf90_enddef(remap%did)
!
   return
   end subroutine remap_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine remap_finalize(remap)
!-------------------------------------------------------------------------------
   implicit none
!
   type(remap_t), intent(inout) :: remap
! local variables
   integer(i4) :: vv, err
!
   if (remap%mat%has_llinfo) then
     call rm_write_ll_info(remap%mat,trim(ll_file))
     call write_ll_info(remap,trim(ll_file))
   endif
!
   err = nf90_close(remap%sid)
   call check_nc(err,'close src file')
   err = nf90_close(remap%did)
   call check_nc(err,'close dst file')
!
   call rm_finalize(remap%mat)
   do vv = 1,remap%nvars
     call buffer_finalize(remap%src(vv))
     call buffer_finalize(remap%dst(vv))
   enddo
   deallocate(remap%sdims)
   deallocate(remap%ddims)
   deallocate(remap%vars)
   deallocate(remap%src)
   deallocate(remap%dst)
!
   return
   end subroutine remap_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine remap_driver(remap)
!-------------------------------------------------------------------------------
   implicit none
!
   type(remap_t), intent(inout) :: remap
! local variables
   integer(i4) :: nvars, vv, dd
   integer(i4), dimension(:), allocatable :: src_dimsizes, dst_dimsizes
!
   nvars = remap%nvars
   print*,' '
   print*,'## start : processing file : ',trim(remap%infile),' ##'
!
   do vv = 1,nvars
!
     print*,' '
     print*,' 0. processing variable : '//trim(remap%vars(vv)%name)
     allocate(src_dimsizes(remap%vars(vv)%sndims))
     allocate(dst_dimsizes(remap%vars(vv)%dndims))
     src_dimsizes(:) = remap%vars(vv)%sdims(1:remap%vars(vv)%sndims)%size
     dst_dimsizes(:) = remap%vars(vv)%ddims(1:remap%vars(vv)%dndims)%size
     print*,' 1. initialize buffer(src) '
     call buffer_initialize(remap%src(vv),remap%vars(vv)%sndims,src_dimsizes)
     print*,' 2. initialize buffer(dst) '
     call buffer_initialize(remap%dst(vv),remap%vars(vv)%dndims,dst_dimsizes)
     deallocate(src_dimsizes)
     deallocate(dst_dimsizes)
 ! 
     print*,' 3. get source variable : '//trim(remap%vars(vv)%name)
     call get_src_file(remap,vv,remap%src(vv))
     if (remap%vars(vv)%isremap) then
       print*,' 4. remapping!(need remap : true) '
     else
       print*,' 4. remapping!(need remap : false) '
     endif
     call do_remap(remap,remap%vars(vv)%isremap,remap%src(vv),remap%dst(vv))
     print*,'  - src (min/max) = ',buffer_min(remap%src(vv)),buffer_max(remap%src(vv))
     print*,'  - dst (min/max) = ',buffer_min(remap%dst(vv)),buffer_max(remap%dst(vv))
     print*,' 5. put destination variable : '//trim(remap%vars(vv)%name)
     call put_dst_file(remap,vv,remap%dst(vv))
     print*,' '
!
   enddo ! do vv = 1,nvars
!
   print*,'## end : processing file : ',trim(remap%infile),' ##'
   print*,' '
   print*,' '
!
   return
   end subroutine remap_driver
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_dimid(dimt, dimname) result(id)
!-------------------------------------------------------------------------------
   implicit none
!
   type(dim_t),      intent(in   ) :: dimt(:)
   character(len=*), intent(in   ) :: dimname
   integer(i4)                     :: id
! local variables
   integer(i4) :: ndims, i
!
   id = -1
   ndims = size(dimt)
!
   do i = 1,ndims
     if (trim(dimt(i)%name).eq.trim(dimname)) then
       id = dimt(i)%id
       exit
     endif
   enddo
!
   if (id.eq.-1) then
     print*,'check dimname in get_dimid...',trim(dimname)
     stop
   endif
!
   end function
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_dimix(dimt, dimname) result(idx)
!-------------------------------------------------------------------------------
   implicit none
!
   type(dim_t),      intent(in   ) :: dimt(:)
   character(len=*), intent(in   ) :: dimname
   integer(i4) :: idx
! local variables
   integer(i4) :: ndims, i
!
   idx = -1
   ndims = size(dimt)
!
   do i = 1,ndims
     if (trim(dimt(i)%name).eq.trim(dimname)) then
       idx = i
       exit
     endif
   enddo
!
   if (idx.eq.-1) then
     print*,'check dimname in get_dimix...',trim(dimname)
     stop
   endif
!
   end function
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine do_remap(remap, hasremap, src, dst)
!-------------------------------------------------------------------------------
   implicit none
!
   type(remap_t),  intent(inout) :: remap
   logical(l4),    intent(in   ) :: hasremap
   type(buffer_t), intent(in   ) :: src
   type(buffer_t), intent(inout) :: dst
! local variables
   integer(i4) :: sndims, dndims, dimsizes(1), dst_size, i1, i2, i3, i4, i5
   real(r8), allocatable :: tmp(:)
!
   sndims = src%ndims
   dndims = dst%ndims
!
! remap
   if (.not.hasremap) then
     if (sndims.eq.1) then
       dst%d1 = src%d1
     elseif (sndims.eq.2) then
       dst%d2 = src%d2
     elseif (sndims.eq.3) then
       dst%d3 = src%d3
     elseif (sndims.eq.4) then
       dst%d4 = src%d4
     elseif (sndims.eq.5) then
       dst%d5 = src%d5
     endif
   else
     if (.not.remap%src_ll.and..not.remap%dst_ll) then
       if (sndims.eq.1) then
         call rm_remap(remap%mat,src%d1(:),dst%d1(:))
       elseif (sndims.eq.2) then
         do i2 = 1,src%dimsizes(2)
           call rm_remap(remap%mat,src%d2(:,i2),dst%d2(:,i2))
         enddo
       elseif (sndims.eq.3) then
         do i3 = 1,src%dimsizes(3)
         do i2 = 1,src%dimsizes(2)
           call rm_remap(remap%mat,src%d3(:,i2,i3),dst%d3(:,i2,i3))
         enddo
         enddo
       elseif (sndims.eq.4) then
         do i4 = 1,src%dimsizes(4)
         do i3 = 1,src%dimsizes(3)
         do i2 = 1,src%dimsizes(2)
           call rm_remap(remap%mat,src%d4(:,i2,i3,i4),dst%d4(:,i2,i3,i4))
         enddo
         enddo
         enddo
       elseif (sndims.eq.5) then
         do i5 = 1,src%dimsizes(5)
         do i4 = 1,src%dimsizes(4)
         do i3 = 1,src%dimsizes(3)
         do i2 = 1,src%dimsizes(2)
           call rm_remap(remap%mat,src%d5(:,i2,i3,i4,i5),dst%d5(:,i2,i3,i4,i5))
         enddo
         enddo
         enddo
         enddo
       endif
     elseif (.not.remap%src_ll.and.remap%dst_ll) then
       if (sndims.eq.1) then
         dst_size = size(dst%d2(:,:))
         allocate(tmp(dst_size))
         call rm_remap(remap%mat,src%d1(:),tmp)
         dst%d2 = reshape(tmp,shape(dst%d2))
       elseif (sndims.eq.2) then
         dst_size = size(dst%d3(:,:,1))
         allocate(tmp(dst_size))
         do i2 = 1,src%dimsizes(2)
           call rm_remap(remap%mat,src%d2(:,i2),tmp)
           dst%d3(:,:,i2) = reshape(tmp,shape(dst%d3(:,:,i2)))
         enddo
       elseif (sndims.eq.3) then
         dst_size = size(dst%d4(:,:,1,1))
         allocate(tmp(dst_size))
         do i3 = 1,src%dimsizes(3)
         do i2 = 1,src%dimsizes(2)
           call rm_remap(remap%mat,src%d3(:,i2,i3),tmp)
           dst%d4(:,:,i2,i3) = reshape(tmp,shape(dst%d4(:,:,i2,i3)))
         enddo
         enddo
       elseif (sndims.eq.4) then
         dst_size = size(dst%d5(:,:,1,1,1))
         allocate(tmp(dst_size))
         do i4 = 1,src%dimsizes(4)
         do i3 = 1,src%dimsizes(3)
         do i2 = 1,src%dimsizes(2)
           call rm_remap(remap%mat,src%d4(:,i2,i3,i4),tmp)
           dst%d5(:,:,i2,i3,i4) = reshape(tmp,shape(dst%d5(:,:,i2,i3,i4)))
         enddo
         enddo
         enddo
       endif
       deallocate(tmp)
     elseif (remap%src_ll.and..not.remap%dst_ll) then
       dimsizes(1) = src%dimsizes(1)*src%dimsizes(2)
       if (sndims.eq.1) then
         print *, 'error'
         stop
       elseif (sndims.eq.2) then
         call rm_remap(remap%mat,reshape(src%d2(:,:),dimsizes),dst%d1(:))
       elseif (sndims.eq.3) then
         do i3 = 1,src%dimsizes(3)
           call rm_remap(remap%mat,reshape(src%d3(:,:,i3),dimsizes),dst%d2(:,i3))
         enddo
       elseif (sndims.eq.4) then
         do i4 = 1,src%dimsizes(4)
         do i3 = 1,src%dimsizes(3)
           call rm_remap(remap%mat,reshape(src%d4(:,:,i3,i4),dimsizes),dst%d3(:,i3,i4))
         enddo
         enddo
       elseif (sndims.eq.5) then
         do i5 = 1,src%dimsizes(5)
         do i4 = 1,src%dimsizes(4)
         do i3 = 1,src%dimsizes(3)
           call rm_remap(remap%mat,reshape(src%d5(:,:,i3,i4,i5),dimsizes),dst%d4(:,i3,i4,i5))
         enddo
         enddo
         enddo
       endif
     elseif (remap%src_ll.and.remap%dst_ll) then
       dimsizes(1) = src%dimsizes(1)*src%dimsizes(2)
       if (sndims.eq.1) then
         print *, 'error 1'
         stop
       elseif (sndims.eq.3) then
         dst_size = size(dst%d3(:,:,1))
         allocate(tmp(dst_size))
         do i3 = 1,src%dimsizes(2)
           call rm_remap(remap%mat,reshape(src%d3(:,:,i3),dimsizes),tmp)
           dst%d3(:,:,i3) = reshape(tmp,shape(dst%d3(:,:,i3)))
         enddo
       elseif (sndims.eq.4) then
         dst_size = size(dst%d4(:,:,1,1))
         allocate(tmp(dst_size))
         do i4 = 1,src%dimsizes(3)
         do i3 = 1,src%dimsizes(2)
           call rm_remap(remap%mat,reshape(src%d4(:,:,i3,i4),dimsizes),tmp)
           dst%d4(:,:,i3,i4) = reshape(tmp,shape(dst%d4(:,:,i3,i4)))
         enddo
         enddo
       elseif (sndims.eq.5) then
         dst_size = size(dst%d5(:,:,1,1,1))
         allocate(tmp(dst_size))
         do i5 = 1,src%dimsizes(4)
         do i4 = 1,src%dimsizes(3)
         do i3 = 1,src%dimsizes(2)
           call rm_remap(remap%mat,reshape(src%d5(:,:,i3,i4,i5),dimsizes),tmp)
           dst%d5(:,:,i3,i4,i5) = reshape(tmp,shape(dst%d5(:,:,i3,i4,i5)))
         enddo
         enddo
         enddo
       endif
       deallocate(tmp)
     endif
   endif
!
   return
   end subroutine do_remap
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_src_file(remap, vv, buf)
!-------------------------------------------------------------------------------
   implicit none
!
   type(remap_t),  intent(in   ) :: remap
   integer(i4),    intent(in   ) :: vv
   type(buffer_t), intent(inout) :: buf
! local variables
   integer(i4) :: ndims
   integer(i4) :: err
!
   ndims = buf%ndims
!
   if     (ndims.eq.1) then
     err = nf90_get_var(remap%sid,remap%vars(vv)%sid,buf%d1)
   elseif (ndims.eq.2) then
     err = nf90_get_var(remap%sid,remap%vars(vv)%sid,buf%d2)
   elseif (ndims.eq.3) then
     err = nf90_get_var(remap%sid,remap%vars(vv)%sid,buf%d3)
   elseif (ndims.eq.4) then
     err = nf90_get_var(remap%sid,remap%vars(vv)%sid,buf%d4)
   elseif (ndims.eq.5) then
     err = nf90_get_var(remap%sid,remap%vars(vv)%sid,buf%d5)
   endif
!
   return
   end subroutine get_src_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine put_dst_file(remap, vv, buf)
!-------------------------------------------------------------------------------
   implicit none
!
   type(remap_t),  intent(in   ) :: remap
   integer(i4),    intent(in   ) :: vv
   type(buffer_t), intent(inout) :: buf
! local variables
   integer(i4) :: ndims, err
   integer(i4) :: d1, d2, d3, d4, d5, nlon, nlat
!
   ndims = buf%ndims
!
   if     (ndims.eq.1) then
     err = nf90_put_var(remap%did,remap%vars(vv)%did,buf%d1)
   elseif (ndims.eq.2) then
     err = nf90_put_var(remap%did,remap%vars(vv)%did,buf%d2)
   elseif (ndims.eq.3) then
     err = nf90_put_var(remap%did,remap%vars(vv)%did,buf%d3)
   elseif (ndims.eq.4) then
     err = nf90_put_var(remap%did,remap%vars(vv)%did,buf%d4)
   elseif (ndims.eq.5) then
     err = nf90_put_var(remap%did,remap%vars(vv)%did,buf%d5)
   endif
   call check_nc(err,'put_dst_file')
!
   return
   end subroutine put_dst_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_ll_info(remap, oufile)
!-------------------------------------------------------------------------------
   implicit none
!
   type(remap_t),    intent(inout) :: remap
   character(len=*), intent(in   ) :: oufile
! local variables
   integer(i4) :: wmod, info, err
   integer(i4) :: fid, iv
   integer(i4) :: src_did, src_vid
   integer(i4) :: dst_did, dst_vid
!
   iv = 1
! file
   wmod = ior(nf90_clobber,nf90_64bit_offset)
   err = nf90_open(trim(oufile),nf90_write,fid)
   call check_nc(err,'write_ll_info: open')
!
! dimensions
   err = nf90_inq_dimid(fid,'src_size',src_did)
   call check_nc(err,'write_ll_info: src_size')
   err = nf90_inq_dimid(fid,'dst_size',dst_did)
   call check_nc(err,'write_ll_info: dst_size')
!
! variables
   err = nf90_redef(fid)
   err = nf90_def_var(fid,'src_var',nf90_real8,src_did,src_vid)
   call check_nc(err,'write_ll_info: src_var')
   err = nf90_def_var(fid,'dst_var',nf90_real8,dst_did,dst_vid)
   call check_nc(err,'write_ll_info: dst_var')
   err = nf90_enddef(fid)
   call check_nc(err,'write_ll_info: enddef')
!
! write
   if     (remap%src(iv)%ndims.eq.1) then
     err = nf90_put_var(fid,src_vid,remap%src(iv)%d1)
   elseif (remap%src(iv)%ndims.eq.2) then
     err = nf90_put_var(fid,src_vid,remap%src(iv)%d2(:,1))
   elseif (remap%src(iv)%ndims.eq.3) then
     err = nf90_put_var(fid,src_vid,remap%src(iv)%d3(:,1,1))
   elseif (remap%src(iv)%ndims.eq.4) then
     err = nf90_put_var(fid,src_vid,remap%src(iv)%d4(:,1,1,1))
   endif
   if     (remap%dst(iv)%ndims.eq.1) then
     err = nf90_put_var(fid,dst_vid,remap%dst(iv)%d1)
   elseif (remap%dst(iv)%ndims.eq.2) then
     err = nf90_put_var(fid,dst_vid,remap%dst(iv)%d2(:,1))
   elseif (remap%dst(iv)%ndims.eq.3) then
     err = nf90_put_var(fid,dst_vid,remap%dst(iv)%d3(:,1,1))
   elseif (remap%dst(iv)%ndims.eq.4) then
     err = nf90_put_var(fid,dst_vid,remap%dst(iv)%d4(:,1,1,1))
   endif
   call check_nc(err,'write_ll_info: put_var')
!
   err = nf90_close(fid)
   call check_nc(err,'write_ll_info: close')
!
   return
   end subroutine write_ll_info
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
   if (err.ne.nf90_noerr) then
     print*,trim(message)//': '//trim(nf90_strerror(err))
     stop
   endif
!
   return
   end subroutine check_nc
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine debug()
!-------------------------------------------------------------------------------
   implicit none
! local variables
   integer(i4), save :: num = 0
!
   num = num+1
   print '(a,i3)', 'DEBUG: step ',num
!
   return
   end subroutine debug
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module remapper
