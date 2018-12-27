!-------------------------------------------------------------------------------
   module remap_netcdf
!-------------------------------------------------------------------------------
!
!  abstract : 
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds,  only: i4, r8, l4
   use buffer, only: buffer_t
   use netcdf
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   integer, parameter :: max_ndims   = 5
   integer, parameter :: max_strname = 32
!
   type dim_t
     character(len=max_strname) :: name         ! ini
     integer(i4)                :: id           ! netcdf
     integer(i4)                :: ssize, dsize ! netcdf
   end type dim_t

   type var_t
     character(len=max_strname) :: name                ! ini
     integer(i4)                :: ndims               ! ini
     character(len=max_strname) :: dimnames(max_ndims) ! ini
     integer(i4) :: sid, did                           ! netcdf
     integer(i4) :: xtype
     integer(i4) :: dimsizes(max_ndims)                ! netcdf
     integer(i4) :: dimids(max_ndims)                  ! netcdf
     logical(l4) :: isremap
   end type var_t
!
   type rnc_t
     character(len=512) :: infile, oufile
     integer(i4)        :: sid, did ! source and destination netcdf id
     integer(i4)        :: ndims
     integer(i4)        :: r_dim    ! dimension's index for remapping
     integer(i4)        :: nvars
     type(dim_t), dimension(:), allocatable :: dims
     type(var_t), dimension(:), allocatable :: vars
     logical(l4)        :: dst_ll
     integer(i4)        :: nlat, nlon, x_id, y_id
   end type rnc_t
!
   public :: rnc_t
   public :: rnc_initialize, rnc_finalize, rnc_get_src_file, rnc_put_dst_file
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_dimid(io_t, dimname) result(id)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rnc_t),      intent(in   ) :: io_t
   character(len=*), intent(in   ) :: dimname
   integer(i4)                     :: id
! local variables
   integer(i4) :: i
!
   id = -1
!
   do i = 1,io_t%ndims
     if (trim(io_t%dims(i)%name).eq.trim(dimname)) then
       id = io_t%dims(i)%id
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
   function get_dimidx(io_t, dimname) result(idx)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rnc_t),      intent(in   ) :: io_t
   character(len=*), intent(in   ) :: dimname
   integer(i4) :: idx
! local variables
   integer(i4) :: i
!
   idx =-1
!
   do i = 1,io_t%ndims
     if (trim(io_t%dims(i)%name).eq.trim(dimname)) then
       idx = i
       exit
     endif
   enddo
!
   if (idx.eq.-1) then
     print*,'check dimname in get_dimidx...',trim(dimname)
     stop
   endif
!
   end function
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rnc_initialize(io_t, infile, oufile, ssize, dsize,               &
       ndims, dimnames, nvars, varnames, ndims_vars, dimnames_vars, nlat, nlon)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rnc_t)     ,                   intent(  out) :: io_t
   character(len=*),                   intent(in   ) :: infile, oufile
   integer(i4)     ,                   intent(in   ) :: ssize, dsize ! source size, dest size 
   integer(i4)     ,                   intent(in   ) :: ndims
   character(len=*), dimension(ndims), intent(in   ) :: dimnames
   integer(i4)     ,                   intent(in   ) :: nvars
   character(len=*), dimension(nvars), intent(in   ) :: varnames
   integer(i4)     , dimension(nvars), intent(in   ) :: ndims_vars
   character(len=*), dimension(:,:),   intent(in   ) :: dimnames_vars
   integer(i4)     , optional,         intent(in   ) :: nlat, nlon
! local variables
   integer(i4) :: dd, vv, dimidx, ierr, dimids(10)
!
   io_t%infile = trim(infile)
   io_t%oufile = trim(oufile)
!
!===================================================================================
! INPUT
!===================================================================================
! dimensions
   io_t%ndims = ndims
   allocate(io_t%dims(ndims))
   do dd = 1,ndims
     io_t%dims(dd)%name = trim(dimnames(dd))
   enddo
!
! variables
   io_t%nvars = nvars
   allocate(io_t%vars(nvars))
   do vv = 1,nvars
     io_t%vars(vv)%name = trim(varnames(vv))
     io_t%vars(vv)%ndims = ndims_vars(vv)
     do dd = 1,ndims_vars(vv)
       io_t%vars(vv)%dimnames(dd) = dimnames_vars(dd,vv)
     enddo
   enddo
!
! read dimensions
   ierr = nf90_open(trim(infile),nf90_nowrite,io_t%sid)
   call fail_message(ierr,'cannot open file : '//trim(infile))
   io_t%r_dim = -1
   do dd = 1,ndims
     ierr = nf90_inq_dimid(io_t%sid,trim(io_t%dims(dd)%name),io_t%dims(dd)%id)
     call fail_message(ierr,'check dimension name : '//trim(io_t%dims(dd)%name))
     ierr = nf90_inquire_dimension(io_t%sid,io_t%dims(dd)%id,len=io_t%dims(dd)%ssize)
     call fail_message(ierr,'cannot read dimsize ')
     if (io_t%dims(dd)%ssize.eq.ssize) then
       io_t%r_dim          = dd
       io_t%dims(dd)%dsize = dsize
     else
       io_t%dims(dd)%dsize = io_t%dims(dd)%ssize
     endif
   enddo
!
! read variable info
   do vv = 1,nvars
     ierr = nf90_inq_varid(io_t%sid,trim(io_t%vars(vv)%name),io_t%vars(vv)%sid)
     call fail_message(ierr,'check variable name : '//trim(io_t%vars(vv)%name))
     ierr = nf90_inquire_variable(io_t%sid,io_t%vars(vv)%sid,xtype=io_t%vars(vv)%xtype)
     call fail_message(ierr,'cannot read xtype ')
     io_t%vars(vv)%xtype=nf90_double
     io_t%vars(vv)%isremap = .false.
     do dd = 1,ndims_vars(vv)
       dimidx = get_dimidx(io_t,trim(io_t%vars(vv)%dimnames(dd)))
       io_t%vars(vv)%dimsizes(dd) = io_t%dims(dimidx)%ssize
       io_t%vars(vv)%dimids(dd)   = get_dimid(io_t,trim(io_t%vars(vv)%dimnames(dd)))
       if (dd.eq.1.and.dimidx.eq.io_t%r_dim) then
         io_t%vars(vv)%isremap = .true.
       endif
     enddo
   enddo
!
!===================================================================================
! OUTPUT
!===================================================================================
   if (present(nlat).and.present(nlon)) then
     io_t%dst_ll = .true.
     io_t%nlon   = nlon
     io_t%nlat   = nlat
   else
     io_t%dst_ll = .false.
   endif
!
   ierr = nf90_create(trim(oufile),ior(nf90_clobber,nf90_64bit_offset),io_t%did)
   call fail_message(ierr,'cannot create file : '//trim(oufile))
!
! define dimenisons
   do dd = 1,ndims
     ierr = nf90_def_dim(io_t%did,trim(io_t%dims(dd)%name),io_t%dims(dd)%dsize,io_t%dims(dd)%id)
     call fail_message(ierr,'cannot define dimension : '//trim(io_t%dims(dd)%name))
   enddo
   if (io_t%dst_ll) then
     ierr = nf90_def_dim(io_t%did,'nlon',nlon,io_t%x_id)
     ierr = nf90_def_dim(io_t%did,'nlat',nlat,io_t%y_id)
   endif
!
! define variables
   do vv = 1,nvars
     if (io_t%vars(vv)%isremap.and.io_t%dst_ll) then
       io_t%vars(vv)%dimids(1) = io_t%x_id
       io_t%vars(vv)%dimids(2) = io_t%y_id
       do dd = 2,ndims_vars(vv)
         io_t%vars(vv)%dimids(dd+1) = get_dimid(io_t,trim(io_t%vars(vv)%dimnames(dd)))
       enddo
       ierr = nf90_def_var(io_t%did,trim(io_t%vars(vv)%name),io_t%vars(vv)%xtype,&
                          io_t%vars(vv)%dimids(1:ndims_vars(vv)+1),io_t%vars(vv)%did)
       call fail_message(ierr,'cannot define variable : '//trim(io_t%vars(vv)%name))
     else
       do dd = 1,ndims_vars(vv)
         io_t%vars(vv)%dimids(dd) = get_dimid(io_t,trim(io_t%vars(vv)%dimnames(dd)))
       enddo
       ierr = nf90_def_var(io_t%did,trim(io_t%vars(vv)%name),io_t%vars(vv)%xtype,&
                          io_t%vars(vv)%dimids(1:ndims_vars(vv)),io_t%vars(vv)%did)
       call fail_message(ierr,'cannot define variable : '//trim(io_t%vars(vv)%name))
     endif
   enddo
!
   ierr = nf90_enddef(io_t%did)
!
   return
   end subroutine rnc_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rnc_finalize(io_t)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rnc_t), intent(inout) :: io_t
! local variables
   integer(i4) :: ierr
!
! close
   ierr = nf90_close(io_t%sid)
   ierr = nf90_close(io_t%did)
!
   return
   end subroutine rnc_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rnc_get_src_file(io_t, vv, buf)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rnc_t),    intent(in   ) :: io_t
   integer(i4),    intent(in   ) :: vv
   type(buffer_t), intent(inout) :: buf
! local variables
   integer(i4) :: ndims
   integer(i4) :: ierr
!
   ndims = buf%ndims
!
   if     (ndims.eq.1) then
     ierr = nf90_get_var(io_t%sid,io_t%vars(vv)%sid,buf%d1)
   elseif (ndims.eq.2) then
     ierr = nf90_get_var(io_t%sid,io_t%vars(vv)%sid,buf%d2)
   elseif (ndims.eq.3) then
     ierr = nf90_get_var(io_t%sid,io_t%vars(vv)%sid,buf%d3)
   elseif (ndims.eq.4) then
     ierr = nf90_get_var(io_t%sid,io_t%vars(vv)%sid,buf%d4)
   elseif (ndims.eq.5) then
     ierr = nf90_get_var(io_t%sid,io_t%vars(vv)%sid,buf%d5)
   endif
!
   return
   end subroutine rnc_get_src_file
#if 0
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rnc_put_dst_file(io_t, vv, buf)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rnc_t),    intent(in   ) :: io_t
   integer(i4),    intent(in   ) :: vv
   type(buffer_t), intent(inout) :: buf
! local variables
   integer(i4) :: ndims
   integer(i4) :: ierr
!
   ndims = buf%ndims
!
   if     (ndims.eq.1) then
     ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,buf%d1)
   elseif (ndims.eq.2) then
     ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,buf%d2)
   elseif (ndims.eq.3) then
     ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,buf%d3)
   elseif (ndims.eq.4) then
     ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,buf%d4)
   elseif (ndims.eq.5) then
     ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,buf%d5)
   endif
!
   return
   end subroutine rnc_put_dst_file
#endif
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rnc_put_dst_file(io_t, vv, buf)
!-------------------------------------------------------------------------------
   implicit none
!
   type(rnc_t),    intent(in   ) :: io_t
   integer(i4),    intent(in   ) :: vv
   type(buffer_t), intent(inout) :: buf
! local variables
   integer(i4) :: ndims, ierr
   integer(i4) :: d1, d2, d3, d4, d5, d6, nlon, nlat
!
   ndims = buf%ndims
!
   if (io_t%vars(vv)%isremap.and.io_t%dst_ll) then
     nlon = io_t%nlon
     nlat = io_t%nlat
     if     (ndims.eq.1) then
       ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,reshape(buf%d1,shape=(/nlon,nlat/)))
     elseif (ndims.eq.2) then
       d2 = size(buf%d2,dim=2)
       ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,reshape(buf%d2,shape=(/nlon,nlat,d2/)))
     elseif (ndims.eq.3) then
       d2 = size(buf%d3,dim=2)
       d3 = size(buf%d3,dim=3)
       ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,reshape(buf%d3,shape=(/nlon,nlat,d2,d3/)))
     elseif (ndims.eq.4) then
       d2 = size(buf%d4,dim=2)
       d3 = size(buf%d4,dim=3)
       d4 = size(buf%d4,dim=4)
       ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,reshape(buf%d4,shape=(/nlon,nlat,d2,d3,d4/)))
     elseif (ndims.eq.5) then
       d2 = size(buf%d5,dim=2)
       d3 = size(buf%d5,dim=3)
       d4 = size(buf%d5,dim=4)
       d5 = size(buf%d5,dim=5)
       ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,reshape(buf%d5,shape=(/nlon,nlat,d2,d3,d4,d5/)))
     endif
   else
     if     (ndims.eq.1) then
       ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,buf%d1)
     elseif (ndims.eq.2) then
       ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,buf%d2)
     elseif (ndims.eq.3) then
       ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,buf%d3)
     elseif (ndims.eq.4) then
       ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,buf%d4)
     elseif (ndims.eq.5) then
       ierr = nf90_put_var(io_t%did,io_t%vars(vv)%did,buf%d5)
     endif
   endif

!
   return
   end subroutine rnc_put_dst_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine fail_message(ierr, message)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),      intent(in   ) :: ierr
   character(len=*), intent(in   ) :: message
!
   if (ierr.ne.0) then
     print*,trim(message)
     stop
   endif
!
   return
   end subroutine fail_message
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module remap_netcdf
