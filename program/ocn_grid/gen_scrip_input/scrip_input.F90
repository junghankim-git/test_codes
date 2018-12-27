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
   print*,'start : writting file : ',trim(input%filename)
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
   print*,'end : writting file : ',trim(input%filename)
!
   return
!-------------------------------------------------------------------------------
   end subroutine scrip_write
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
   subroutine check_corners(idx, lon0, lat0, npts, lons, lats)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: idx
   real(r8)   ,                  intent(in   ) :: lon0, lat0 ! radian
   integer(i4),                  intent(in   ) :: npts
   real(r8)   , dimension(npts), intent(inout) :: lons, lats ! radian
! local variables
   integer(i4) :: i, j, nzeros
   logical(l4) :: iszero, issameval, iscclock, isin
   real(r8), dimension(npts) :: rlons, rlats ! radian
!
! check min/max
   if ((minval(lons).lt.0.0_r8).or.(maxval(lons).gt.2.0_r8*pi)) then
     call warn_message(idx,'min/max lons',lon0,lat0,npts,lons,lats)
   endif
   if ((minval(lats).lt.-0.5_r8*pi).or.(maxval(lats).gt.0.5_r8*pi)) then
     call warn_message(idx,'min/max lats',lon0,lat0,npts,lons,lats)
   endif
!
   iszero    = .false.
   nzeros    = 0
   issameval = .false.
!
! fill blank (0,0)->(lon,lat)
! copy lons,lats -> rlons,rlat
   do i = 1,npts
     if ((lons(i).eq.0.0_r8.and.lats(i).eq.0.0_r8).or.(lons(i).eq.none_.and.lats(i).eq.none_)) then
       lons(i) = lons(i-1)
       lats(i) = lats(i-1)
       iszero = .true.
       nzeros = nzeros+1
     endif
     rlons(i) = lons(i)
     rlats(i) = lats(i)
   enddo
!
! 1) check same values
   do i = 1,npts
     do j = 1,npts
       if ((i.ne.j).and.(lons(i).eq.lons(j)).and.(lats(i).eq.lats(j))) then
         issameval = .true.
       endif
     enddo
   enddo
!
! 2) check counter-clockwise
   iscclock = iscounterclock(npts-nzeros,rlons,rlats)
!
! 3) check interior (fail)
!   isin     = isinterior(lon0,lat0,npts,rlons,rlats)
!
   ! write points.txt
   write(ulog,*) rlons(:)
   write(ulog,*) rlats(:)
!
   !if (issameval)     call warn_message(idx,'has double',npts,rlons,rlats)
   if (.not.iscclock) call warn_message(idx,'c-clockwise',lon0,lat0,npts,rlons,rlats)
!   if (.not.isin)     call warn_message(idx,'interior',lon0,lat0,npts,rlons,rlats)
!
   return
!-------------------------------------------------------------------------------
   end subroutine check_corners
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine warn_message(idx,message,lon0,lat0,n,lons,lats)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),               intent(in   ) :: idx
   character(len=*),          intent(in   ) :: message
   real(r8)   ,               intent(in   ) :: lon0, lat0
   integer(i4),               intent(in   ) :: n
   real(r8)   , dimension(n), intent(inout) :: lons, lats
!
   integer(i4)        :: i
   character(len=256) :: string
!
   write(string,'(a,i7,a,a12,a,f7.4,a,f7.4,a)') '* ',idx,' : check (',trim(message),'), (lon0,lat0) = (',lon0,', ',lat0,')'
   write(*,*) trim(string)
!
   do i = 1,n
     if (i.eq.1) then
       write(string,'(a22,a,f7.4,a,f7.4,a)') 'corners = ',' (',lons(i),',',lats(i),')'
     else
       string = trim(string)//', '
       write(string,'(a,a,f7.4,a,f7.4,a)') trim(string),' (',lons(i),',',lats(i),')'
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
   function iscounterclock(npts, rlons, rlats) result(iscclock)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: npts
   real(r8), dimension(npts), intent(inout) :: rlons, rlats ! radian
   logical(l4)                              :: iscclock
! local variables
   integer(i4) :: i, i0, ip
   integer(i4),dimension(npts) :: quadrants
   real(r8)   ,dimension(npts) :: thetas
! get theta
   thetas(:)    = get_thetas(npts,1.0_r8,rlons(:),rlats(:))
   quadrants(:) = get_quadrants(npts,thetas)
   !print *, thetas(:)
! check counter-clockwise
   iscclock = .true.
   do i = 1,npts
     ip = i; i0 = i-1
     if (ip.eq.1) i0 = npts
     if (quadrants(ip).ge.quadrants(i0)) then
       if (thetas(ip).lt.thetas(i0)) then
         iscclock = .false.
         exit
       endif
     else
       if (thetas(ip).lt.(thetas(i0)-2.0_r8*pi)) then
         iscclock = .false.
         exit
       endif
     endif
   enddo
!
!-------------------------------------------------------------------------------
   end function iscounterclock
!-------------------------------------------------------------------------------
!
! fail
!-------------------------------------------------------------------------------
   function isinterior(lon0, lat0, npts, lons, lats) result(isin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)   ,                  intent(in   ) :: lon0, lat0 ! radian
   integer(i4),                  intent(in   ) :: npts
   real(r8)   , dimension(npts), intent(inout) :: lons, lats ! radian
   logical(l4)                                 :: isin
! local variables
   integer(i4) :: i
   real(r8)    :: minv, maxv
!
   isin = .true.
!
   minv = minval(lons)
   maxv = maxval(lons)
   if (lon0.lt.minv.or.lon0.gt.maxv) isin = .false.
   minv = minval(lats)
   maxv = maxval(lats)
   if (lat0.lt.minv.or.lat0.gt.maxv) isin = .false.
!
   return
!-------------------------------------------------------------------------------
   end function isinterior
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_thetas(npts, r, rlons, rlats) result(thetas)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: npts
   real(r8),                  intent(in   ) :: r
   real(r8), dimension(npts), intent(inout) :: rlons, rlats ! radian
   real(r8), dimension(npts) :: thetas
! local variables
   integer(i4) :: i
   real(r8),dimension(npts) :: x,y,z,rx,ry,rz
   real(r8) :: x0,y0,z0,r0,rx0,ry0,rz0
   real(r8) :: rlon0,rlat0,min4q
!
! lon,lat -> x,y,z
!print '(4(f8.4,x))', rlons(:)
!print '(4(f8.4,x))', rlats(:)
   do i = 1,npts
     x(i) = r*cos(rlats(i))*cos(rlons(i))
     y(i) = r*cos(rlats(i))*sin(rlons(i))
     z(i) = r*sin(rlats(i))
     !print '(5(f8.4,x))', rlons(i),rlats(i),x(i),y(i),z(i)
   enddo
! normalized center (x0,y0,z0)
   x0 = sum(x)/real(npts,8)
   y0 = sum(y)/real(npts,8)
   z0 = sum(z)/real(npts,8)
   r0 = sqrt(x0**2.0_r8+y0**2.0_r8+z0**2.0_r8)
   x0 = x0*r/r0
   y0 = y0*r/r0
   z0 = z0*r/r0
   r0 = sqrt(x0**2.0_r8+y0**2.0_r8+z0**2.0_r8)
!
   rlon0 = xy2angle(x0,y0) ! atan(y0/x0)
   rlat0 = asin(z0/r0)
!
! rotation  z(-lon0),y(lat0-pi/2)
   call rotations('z',-rlon0,npts,x,y,z,rx,ry,rz)
   x = rx; y = ry; z = rz
   call rotations('y',rlat0-pi/2.0_r8,npts,x,y,z,rx,ry,rz)
!
!print *,sqrt(rx(1)**2.0_r8+ry(1)**2.0_r8+rz(1)**2.0_r8)
   do i = 1,npts
     thetas(i) = xy2angle(rx(i),ry(i))
   enddo
!
   end function get_thetas
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rotations(axis, angle, npts, x, y, z, rx, ry, rz)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),          intent(in   ) :: axis
   real(r8),                  intent(in   ) :: angle
   integer(i4),               intent(in   ) :: npts
   real(r8), dimension(npts), intent(in   ) :: x, y, z
   real(r8), dimension(npts), intent(inout) :: rx, ry, rz
! local variables
   integer(i4) :: i
!
   do i = 1,npts
     call rotation(trim(axis),angle,x(i),y(i),z(i),rx(i),ry(i),rz(i))
   enddo
!
   return
   end subroutine rotations
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rotation(axis, angle, x, y, z, rx, ry, rz)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*), intent(in   ) :: axis
   real(r8),         intent(in   ) :: angle
   real(r8),         intent(in   ) :: x, y, z
   real(r8),         intent(inout) :: rx, ry, rz
! local variables
   integer(i4) :: i, j
   real(r8) :: sn, cs
   real(r8), dimension(3, 3) :: mat
!
   sn = sin(angle)
   cs = cos(angle)
!
   if (trim(axis)=='x') then
     mat(1,1) = 1.0_r8
     mat(2,1) = 0.0_r8
     mat(3,1) = 0.0_r8
     mat(1,2) = 0.0_r8
     mat(2,2) = cs
     mat(3,2) =-sn
     mat(1,3) = 0.0_r8
     mat(2,3) = sn
     mat(3,3) = cs
   elseif (trim(axis)=='y') then
     mat(1,1) = cs
     mat(2,1) = 0.0_r8
     mat(3,1) = sn
     mat(1,2) = 0.0_r8
     mat(2,2) = 1.0_r8
     mat(3,2) = 0.0_r8
     mat(1,3) =-sn
     mat(2,3) = 0.0_r8
     mat(3,3) = cs
   elseif (trim(axis)=='z') then
     mat(1,1) = cs
     mat(2,1) =-sn
     mat(3,1) = 0.0_r8
     mat(1,2) = sn
     mat(2,2) = cs
     mat(3,2) = 0.0_r8
     mat(1,3) = 0.0_r8
     mat(2,3) = 0.0_r8
     mat(3,3) = 1.0_r8
   else
     print*,'check axis...'
     stop
   endif
!
   rx = mat(1,1)*x+mat(2,1)*y+mat(3,1)*z
   ry = mat(1,2)*x+mat(2,2)*y+mat(3,2)*z
   rz = mat(1,3)*x+mat(2,3)*y+mat(3,3)*z
!
   return
   end subroutine rotation
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_thetas_old(npts, r, rlons, rlats, quadrants) result(thetas)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: npts
   real(r8),                     intent(in   ) :: r
   real(r8), dimension(npts),    intent(inout) :: rlons, rlats ! radian
   integer(i4), dimension(npts), intent(inout) :: quadrants
   real(r8), dimension(npts) :: thetas
! local variables
   integer(i4) :: i
   real(r8), dimension(npts) :: x, y, z
   real(r8) :: x0, y0, z0, r0
   real(r8) :: rlon0, rlat0, min4q
   real(r8) :: vx, vy
   integer(i4) :: n1quad, n4quad, quad0
!
! check quadrants and shift lon, lat
   n1quad = 0
   n4quad = 0
   quadrants(:) = get_quadrants(npts,rlons)
   min4q = huge(1.0_r8)
   do i = 1,npts
     if (quadrants(i)==1) then
       n1quad = n1quad+1
     elseif (quadrants(i)==4) then
       n4quad = n4quad+1
       min4q = min(rlons(i),min4q)
     endif
   enddo
!
! move 4 quadrant to 1 quadrant
!   if ((n1quad+n4quad==npts).and.(n1quad>0).and.(n4quad>0)) then
   if ((n1quad>0).and.(n4quad>0)) then
     do i = 1,npts
       if (quadrants(i)==4) then
         rlons(i) = rlons(i)-min4q
       else
         rlons(i) = rlons(i)+(2.0_r8*pi-min4q)
       endif
     enddo
   endif
!
! lon,lat -> x,y,z
   do i = 1,npts
     x(i) = r*cos(rlats(i))*cos(rlons(i))
     y(i) = r*cos(rlats(i))*sin(rlons(i))
     z(i) = r*sin(rlats(i))
   enddo
!
! center (x0,y0,z0)
   x0 = sum(x)/real(npts,8)
   y0 = sum(y)/real(npts,8)
   z0 = sum(z)/real(npts,8)
   r0 = sqrt(x0**2.0_r8+y0**2.0_r8+z0**2.0_r8)
   x0 = x0*r/r0
   y0 = y0*r/r0
   z0 = z0*r/r0
   r0 = sqrt(x0**2.0_r8+y0**2.0_r8+z0**2.0_r8)
!
   rlon0 = xy2angle(x0,y0) ! atan(y0/x0)
   rlat0 = asin(z0/r0)
!
   do i = 1,npts
     vx = cos(rlat0)*(rlons(i)-rlon0)
     vy = rlats(i)-rlat0
     thetas(i) = xy2angle(vx,vy)
   enddo
!
   end function get_thetas_old
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_quadrants(npts, radians) result(quads)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: npts
   real(r8), dimension(npts), intent(in   ) :: radians
   integer(i4), dimension(npts) :: quads
! local variables
   integer(i4) :: i
!
   do i = 1,npts
     quads(i) = get_quadrant(radians(i))
   enddo
!
   end function get_quadrants
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_quadrant(radian) result(quad)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: radian
   integer(i4) :: quad
!
   if (radian>=0.0_r8.and.radian<0.5_r8*pi) then
     quad = 1
   elseif (radian>=0.5_r8*pi.and.radian<pi) then
     quad = 2
   elseif (radian>=pi.and.radian<1.5_r8*pi) then
     quad = 3
   elseif (radian>=1.5_r8*pi.and.radian<2.0_r8*pi) then
     quad = 4
   endif
!
   end function get_quadrant
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function xy2angle(x, y) result(theta)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: x, y
   real(r8) :: theta
!
   if (x==0.0_r8.and.y==0.0_r8) then
     !print*, 'x, y is 0...'
     !stop
     theta = 0.0_r8
   ! 1
   elseif (x>0.0_r8.and.y>=0.0_r8.or.x==0.0_r8.and.y>0.0_r8) then
     theta = atan(y/x)
   ! 2
   elseif (x<0.0_r8.and.y>0.0_r8) then ! atan(y/x) =-theta
     theta = pi+atan(y/x)
   ! 3
   elseif (x<0.0_r8.and.y<=0.0_r8) then ! atan(y/x)
     theta = pi+atan(y/x)
   ! 4
   elseif (x>=0.0_r8.and.y<0.0_r8) then
     theta = 2.0_r8*pi+atan(y/x)
   ! other ?
   else
     print*,'check xy2angle...'
     stop
   endif
!
   end function xy2angle
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
