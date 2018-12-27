!-------------------------------------------------------------------------------
   module scrip_input
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
   use kinds, only : i4, l4, r4, r8
   use coordinates
   use netcdf
!
   type scrip_t
     ! filename
     character(len=512) :: filename
     ! dimensions
     integer(i4) :: grid_rank
     integer(i4) :: grid_size
     integer(i4) :: grid_max_corners
     ! global attribute
     character(len=512) :: description
     ! variables
     integer(i4), dimension(:),   allocatable :: grid_dims       ! grid_rank
     integer(i4), dimension(:),   allocatable :: grid_imask      ! grid_size
     real(r8)   , dimension(:),   allocatable :: grid_center_lat ! grid_size
     real(r8)   , dimension(:),   allocatable :: grid_center_lon ! grid_size
     real(r8)   , dimension(:,:), allocatable :: grid_corner_lat ! grid_max_corners, grid_size
     real(r8)   , dimension(:,:), allocatable :: grid_corner_lon ! grid_max_corners, grid_size
   end type scrip_t
!
   real(r8), parameter :: pi = 3.141592653589793238462643383279_r8
   real(r8), parameter :: none_ = -999._r8
   integer(i4) :: nchks = 0
   integer(i4) :: ulog  = 71
!
   public :: scrip_t, none_
   public :: scrip_initialize, scrip_initialize_read, scrip_write, scrip_finalize, scrip_check
   public :: check_corners, radians2degrees, degrees2radians
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine scrip_initialize(input, filename, rank, ncenter, ncorner, description)
!-------------------------------------------------------------------------------
   implicit none
!
   type(scrip_t), intent(inout) :: input
   character(len=*),  intent(in   ) :: filename
   integer(i4),       intent(in   ) :: rank, ncenter, ncorner
   character(len=*),  intent(in   ) :: description
!
   input%filename         = trim(filename)
   input%grid_rank        = rank
   input%grid_size        = ncenter
   input%grid_max_corners = ncorner
   input%description      = trim(description)
!
   allocate(input%grid_dims(rank))
   allocate(input%grid_imask(ncenter))
   allocate(input%grid_center_lat(ncenter))
   allocate(input%grid_center_lon(ncenter))
   allocate(input%grid_corner_lat(ncorner,ncenter))
   allocate(input%grid_corner_lon(ncorner,ncenter))
!
   return
   end subroutine scrip_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine scrip_finalize(input)
!-------------------------------------------------------------------------------
   implicit none
!
   type(scrip_t), intent(inout) :: input
!
   deallocate(input%grid_dims)
   deallocate(input%grid_imask)
   deallocate(input%grid_center_lat)
   deallocate(input%grid_center_lon)
   deallocate(input%grid_corner_lat)
   deallocate(input%grid_corner_lon)
!
   return
   end subroutine scrip_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine scrip_initialize_read(input, filename)
!-------------------------------------------------------------------------------
   implicit none
!
   type(scrip_t), intent(inout) :: input
   character(len=*),  intent(in   ) :: filename
! local variables
   integer(i4) :: fid, rank_did, size_did, corners_did
   integer(i4) :: vid
   integer(i4) :: rank, ncenter, ncorner
   character(len=512) :: description
   integer(i4) :: ierr
! open
   ierr = nf90_open(trim(filename),nf90_nowrite,fid)
   call nc_check(ierr)
! get dimensions
   ierr = nf90_inq_dimid(fid,'grid_rank',rank_did)
   call nc_check(ierr)
   ierr = nf90_inquire_dimension(fid,rank_did,len=rank)
   call nc_check(ierr)
   ierr = nf90_inq_dimid(fid,'grid_size',size_did)
   call nc_check(ierr)
   ierr = nf90_inquire_dimension(fid,size_did,len=ncenter)
   call nc_check(ierr)
   ierr = nf90_inq_dimid(fid,'grid_corners',corners_did)
   call nc_check(ierr)
   ierr = nf90_inquire_dimension(fid,corners_did,len=ncorner)
   call nc_check(ierr)
! get global attribute
   ierr = nf90_get_att(fid,nf90_global,'title',description)
   call nc_check(ierr)
!
! initialize
   call scrip_initialize(input,trim(filename),rank,ncenter,ncorner,trim(description))
!
! get variables
   ierr = nf90_inq_varid(fid,'grid_dims',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%grid_dims)
   call nc_check(ierr)
   ierr = nf90_inq_varid(fid,'grid_imask',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%grid_imask)
   call nc_check(ierr)
   ierr = nf90_inq_varid(fid,'grid_center_lat',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%grid_center_lat)
   call nc_check(ierr)
   ierr = nf90_inq_varid(fid,'grid_center_lon',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%grid_center_lon)
   call nc_check(ierr)
   ierr = nf90_inq_varid(fid,'grid_corner_lat',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%grid_corner_lat)
   call nc_check(ierr)
   ierr = nf90_inq_varid(fid,'grid_corner_lon',vid)
   call nc_check(ierr)
   ierr = nf90_get_var(fid,vid,input%grid_corner_lon)
   call nc_check(ierr)
!
! close
   ierr = nf90_close(fid)
!
   return
!-------------------------------------------------------------------------------
   end subroutine scrip_initialize_read
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine scrip_write(input)
!-------------------------------------------------------------------------------
   implicit none
!
   type(scrip_t), intent(inout) :: input
! local variables
   integer(i4) :: fid, rank_did, size_did, corners_did
   integer(i4) :: dims_vid, imask_vid, clat_vid, clon_vid, cclat_vid, cclon_vid
   integer(i4) :: ierr
!
! open
   print*,'start: write file (',trim(input%filename),')'
   ierr = nf90_create(trim(input%filename),ior(nf90_clobber,nf90_64bit_offset),fid)
   call nc_check(ierr)
! define dimensions
   ierr = nf90_def_dim(fid,'grid_rank',input%grid_rank,rank_did)
   ierr = nf90_def_dim(fid,'grid_size',input%grid_size,size_did)
   ierr = nf90_def_dim(fid,'grid_corners',input%grid_max_corners,corners_did)
! define variables (and attributes)
   ierr = nf90_def_var(fid,'grid_dims',nf90_int,rank_did,dims_vid)
!   ierr = nf90_put_att(fid,dims_vid, 'long_name','grid dimensions')
!   ierr = nf90_put_att(fid,dims_vid, 'units',    '___')
   ierr = nf90_def_var(fid,'grid_imask',nf90_int,size_did,imask_vid)
!   ierr = nf90_put_att(fid,imask_vid,'long_name','grid masking')
   ierr = nf90_put_att(fid,imask_vid,'units','unitless')
   ierr = nf90_def_var(fid,'grid_center_lat',nf90_real8,size_did,clat_vid)
   ierr = nf90_put_att(fid,clat_vid,'units','radians')
   ierr = nf90_def_var(fid,'grid_center_lon',nf90_real8,size_did,clon_vid)
   ierr = nf90_put_att(fid,clon_vid,'units','radians')
   ierr = nf90_def_var(fid,'grid_corner_lat',nf90_real8,(/corners_did,size_did/),cclat_vid)
   ierr = nf90_put_att(fid,cclat_vid,'units','radians')
   ierr = nf90_def_var(fid,'grid_corner_lon',nf90_real8,(/corners_did,size_did/),cclon_vid)
   ierr = nf90_put_att(fid,cclon_vid,'units','radians')
! put global attribute
   ierr = nf90_put_att(fid,nf90_global,'title',input%description)
   call nc_check(ierr)
! end define
   ierr = nf90_enddef(fid)
! put variables
   ierr = nf90_put_var(fid,dims_vid,input%grid_dims)
   ierr = nf90_put_var(fid,imask_vid,input%grid_imask)
   ierr = nf90_put_var(fid,clat_vid,input%grid_center_lat)
   ierr = nf90_put_var(fid,clon_vid,input%grid_center_lon)
   ierr = nf90_put_var(fid,cclat_vid,input%grid_corner_lat)
   ierr = nf90_put_var(fid,cclon_vid,input%grid_corner_lon)
! close
   ierr = nf90_close(fid)
   call nc_check(ierr)
   print*,'end  : write file (',trim(input%filename),')'
!
   return
!-------------------------------------------------------------------------------
   end subroutine scrip_write
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine scrip_check(input)
!-------------------------------------------------------------------------------
   implicit none
!
   type(scrip_t), intent(inout) :: input
! local variables
   integer(i4) :: iup, ncorners
!
   open(ulog,file='points.txt',form='FORMATTED',status='UNKNOWN')
! corners
   !call check_centers(input%grid_size,input%grid_center_lon,input%grid_center_lat)
   ncorners = input%grid_max_corners
   do iup = 1,input%grid_size
     call check_corners(iup,input%grid_center_lon(iup),input%grid_center_lat(iup), &
            ncorners,input%grid_corner_lon(:,iup),input%grid_corner_lat(:,iup))
   enddo
!
   close(ulog)
!
   return
   end subroutine scrip_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_centers(npts, lons, lats)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: npts
   real(r8)   , dimension(npts), intent(inout) :: lons, lats ! radian
! local variables
   integer(i4) :: i, j
!
! 1) remove redundant points
   do i = 1,npts
     do j = i+1,npts
       if (lons(i).eq.lons(j).and.lats(i).eq.lats(j)) then
         print*,'same points :',i,j,lons(i)*180_r8/pi,lats(i)*180_r8/pi
       endif
     enddo
   enddo
!
   return
!-------------------------------------------------------------------------------
   end subroutine check_centers
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_corners(ix, lon0, lat0, npts, lons, lats)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: ix
   real(r8)   ,                  intent(in   ) :: lon0, lat0 ! radian
   integer(i4),                  intent(in   ) :: npts
   real(r8)   , dimension(npts), intent(inout) :: lons, lats ! radian
! local variables
   integer(i4) :: i, j
   real(r8)    :: tlon, tlat
   logical(l4) :: iscclock, isinter
!
! check min/max
   if ((minval(lons).lt.0.0_r8).or.(maxval(lons).gt.2.0_r8*pi)) then
     call warn_message(ix,'min/max lons',lon0,lat0,npts,lons,lats)
   endif
   if ((minval(lats).lt.-0.5_r8*pi).or.(maxval(lats).gt.0.5_r8*pi)) then
     call warn_message(ix,'min/max lats',lon0,lat0,npts,lons,lats)
   endif
!
! 1) clockwise at pole
   if (abs(lat0-pi/2.0_r8).le.1.d0-17) then
     tlon    = lons(3)
     tlat    = lats(3)
     lons(3) = lons(2)
     lats(3) = lats(2)
     lons(2) = tlon
     lats(2) = tlat
   endif
!
! 2) remove redundant points (move to last index)
   call redundant_corners(npts,lons,lats)
!
! 3) remove pole
   call remove_pole(ix,npts,lons,lats)
!
! 4) check counter-clockwise
   iscclock = iscounterclock(npts,lons,lats)
!
! 5) check interior
   isinter  = isinterior(lon0,lat0,npts,lons,lats)
!
   ! write points.txt
   write(ulog,*) lons(:)
   write(ulog,*) lats(:)
!
   if (.not.iscclock) call warn_message(ix,'c-clockwise',lon0,lat0,npts,lons,lats)
   if (.not.isinter)  call warn_message(ix,'interior',lon0,lat0,npts,lons,lats)
!
   return
!-------------------------------------------------------------------------------
   end subroutine check_corners
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine remove_pole(ix,npts, lons, lats)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: ix, npts
   real(r8)   , dimension(npts), intent(inout) :: lons, lats ! radian
! local variables
   integer(i4) :: i, j
!
   do i = 1,npts
     if (abs(pi/2.0_r8-lats(i)).le.1.0d-15) lats(i) = lats(i)-1.0d-8
     if (abs(pi/2.0_r8+lats(i)).le.1.0d-15) lats(i) = lats(i)+1.0d-8
   enddo
!
   return
!-------------------------------------------------------------------------------
   end subroutine remove_pole
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine redundant_corners(npts, lons, lats)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: npts
   real(r8)   , dimension(npts), intent(inout) :: lons, lats ! radian
! local variables
   integer(i4) :: i, j
!
! processing double points
   do i = 1,npts-1
     if (issame(lons(i),lats(i),lons(i+1),lats(i+1))) then
       do j = i,npts-1
         lons(j) = lons(j+1)
         lats(j) = lats(j+1)
       enddo
       lons(npts) = lons(npts-1)
       lats(npts) = lats(npts-1)
     endif
   enddo
   if (issame(lons(npts),lats(npts),lons(1),lats(1))) then
     lons(npts) = lons(npts-1)
     lats(npts) = lats(npts-1)
   endif
!
   return
!-------------------------------------------------------------------------------
   end subroutine redundant_corners
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function issame(a1, a2, b1, b2) result(same)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: a1, a2, b1, b2
   logical(l4)             :: same
! local variables
   real(r8) :: small
!
   small = 0.0
!
   same = .false.
   if (abs(b1-a1).le.small.and.abs(b2-a2).le.small) then
     same = .true.
   endif
!
   return
   end function issame
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine radians2degrees(input)
!-------------------------------------------------------------------------------
   implicit none
!
   type(scrip_t), intent(inout) :: input
! local variables
   integer(i4) :: i, j, ncenter, ncorner
!
   ncenter = input%grid_size
   ncorner = input%grid_max_corners
!
   do i = 1,ncenter
     input%grid_center_lat(i) = input%grid_center_lat(i)*180_r8/pi
     input%grid_center_lon(i) = input%grid_center_lon(i)*180_r8/pi
   enddo
!
   do j = 1,ncenter
     do i = 1,ncorner
       input%grid_corner_lat(i,j) = input%grid_corner_lat(i,j)*180_r8/pi
       input%grid_corner_lon(i,j) = input%grid_corner_lon(i,j)*180_r8/pi
     enddo
   enddo
!
   return
   end subroutine radians2degrees
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine degrees2radians(input)
!-------------------------------------------------------------------------------
   implicit none
!
   type(scrip_t), intent(inout) :: input
! local variables
   integer(i4) :: i, j, ncenter, ncorner
!
   ncenter = input%grid_size
   ncorner = input%grid_max_corners
!
   do i = 1,ncenter
     input%grid_center_lat(i) = input%grid_center_lat(i)*pi/180_r8
     input%grid_center_lon(i) = input%grid_center_lon(i)*pi/180_r8
   enddo
   do j = 1,ncenter
     do i = 1,ncorner
       input%grid_corner_lat(i,j) = input%grid_corner_lat(i,j)*pi/180_r8
       input%grid_corner_lon(i,j) = input%grid_corner_lon(i,j)*pi/180_r8
     enddo
   enddo
!
   return
   end subroutine degrees2radians
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine warn_message(ix,message,lon0,lat0,n,lons,lats)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),               intent(in   ) :: ix
   character(len=*),          intent(in   ) :: message
   real(r8)   ,               intent(in   ) :: lon0, lat0
   integer(i4),               intent(in   ) :: n
   real(r8)   , dimension(n), intent(in   ) :: lons, lats
!
   integer(i4)        :: i, xx, yy
   character(len=256) :: string
!
   !xx = mod(ix,240); if (xx.eq.0) xx = 240
   !yy = ix/240     ; if (xx.eq.0) yy = yy+1
   xx = mod(ix-1,240)+1
   yy = (ix-1)/240+1
   print *, 'xx, yy = ', xx, yy
   write(string,'(a,i7,a,a12,a,f17.13,a,f17.13,a)') '* ',ix,' : check (',trim(message),'), (lon0,lat0) = (',lon0*180.0_r8/pi,', ',lat0*180.0_r8/pi,')'
   write(*,*) trim(string)
!
   do i = 1,n
     if (i.eq.1) then
       write(string,'(a22,a,f17.13,a,f17.13,a)') 'corners = ',' (',lons(i)*180.0_r8/pi,',',lats(i)*180.0_r8/pi,')'
     else
       string = trim(string)//', '
       write(string,'(a,a,f17.13,a,f17.13,a)') trim(string),' (',lons(i)*180.0_r8/pi,',',lats(i)*180.0_r8/pi,')'
     endif
   enddo
   write(*,*) trim(string)
   print*,' '
!
   return
!-------------------------------------------------------------------------------
   end subroutine warn_message
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
   end module scrip_input
!-------------------------------------------------------------------------------
