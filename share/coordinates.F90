!-------------------------------------------------------------------------------
   module coordinates
!-------------------------------------------------------------------------------
!
!  abstract : coordinates control module
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only : i4, l4, r4, r8
!
   private
!
   type coordinate_t
     real(r8) :: lon, lat
     real(r8) :: x, y, z
   end type coordinate_t
!
   real(r8), parameter :: pi       = acos(-1.0_16)
   real(r8), parameter :: rad2deg  = 180.0_r8/pi
   real(r8), parameter :: deg2rad  = pi/180.0_r8
   real(r8), parameter :: zero_thr = 1.0d-9
!
   public :: rad2deg, deg2rad, zero_thr
   public :: coordinate_t, adjust_corners
   public :: cartesian2lonlat, lonlat2cartesian, rotation, rotations
   public :: get_lons_twovector_to_xyplane, get_lonlat_twovector_to_xyplane
   public :: rotation_lonlat_ll0_to_axis !, inv_rotation_lonlat_ll0_to_axis
   public :: rotation_lonslats_ll0_to_axis, inv_rotation_lonslats_ll0_to_axis
   public :: lonlat_to_great_circle_angle, great_circle_angle_to_lonlat
   public :: rotation_lonlat_xaxis_to_ll0
   public :: get_quadrants, get_quadrant, xy2angle, get_thetas, get_theta
   public :: iscounterclock, isinterior
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cartesian2lonlat(x, y, z, r, lon, lat)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: x, y, z
   real(r8), intent(  out) :: r, lon, lat
!
   r = sqrt(x*x+y*y+z*z)
   lon = get_lon(x, y)
   lat = asin(z/r)
!
   return
   end subroutine cartesian2lonlat
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine lonlat2cartesian(r, lon, lat, x, y, z)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: r, lon, lat
   real(r8), intent(inout) :: x, y, z
!
   if (lon.lt.0.0_r8.or.lon.gt.2.0_r8*pi) then
     print*,'check the range of lon... in lonlat2cartesian',lon
   endif
   x = r*cos(lat)*cos(lon)
   y = r*cos(lat)*sin(lon)
   z = r*sin(lat)
   !x = r*cos(real(lat,16))*cos(real(lon,16))
   !y = r*cos(real(lat,16))*sin(real(lon,16))
   !z = r*sin(real(lat,16))
!
   return
   end subroutine lonlat2cartesian
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rotation_lonlat(axis, angle, r1, lon1, lat1, r2, lon2, lat2)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*), intent(in   ) :: axis
   real(r8),         intent(in   ) :: angle
   real(r8),         intent(in   ) :: r1, lon1, lat1
   real(r8),         intent(  out) :: r2, lon2, lat2
! local variables
   integer(i4) :: i
   real(r8)    :: x1, y1, z1, x2, y2, z2
!
   call lonlat2cartesian(r1,lon1,lat1,x1,y1,z1)
   call rotation(trim(axis),angle,x1,y1,z1,x2,y2,z2)
   call cartesian2lonlat(x2,y2,z2,r2,lon2,lat2)
!
   return
   end subroutine rotation_lonlat
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
   real(r8),         intent(  out) :: rx, ry, rz
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
   subroutine rotation_lonlat_ll0_to_axis(lon0, lat0, axis, rr, lon, lat, rrs, rlon, rlat)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),         intent(in   ) :: lon0, lat0
   character(len=*), intent(in   ) :: axis
   real(r8),         intent(in   ) :: rr, lon, lat
   real(r8),         intent(inout) :: rrs, rlon, rlat
! local variables
   real(r8)    ::  xx,   yy,   zz
   real(r8)    :: rxx1, ryy1, rzz1, rxx2, ryy2, rzz2
   character(len=1), dimension(2) :: axises
   real(r8),         dimension(2) :: angles
!
   if     (trim(axis)=='x') then
     axises(1) = 'z'; angles(1) =-lon0
     axises(2) = 'y'; angles(2) = lat0
   elseif (trim(axis)=='y') then
     axises(1) = 'z'; angles(1) = pi/2.0_r8-lon0
     axises(2) = 'x'; angles(2) =-lat0
   elseif (trim(axis)=='z') then
     axises(1) = 'z'; angles(1) =-lon0
     axises(2) = 'y'; angles(2) = lat0-(pi/2.0_r8)
   else
     print*,'check axis...(in rotation_lonlat...) '
     stop
   endif
!
   call lonlat2cartesian(rr,lon,lat,xx,yy,zz)
   call rotation(trim(axises(1)),angles(1), xx , yy , zz ,rxx1,ryy1,rzz1)
   call rotation(trim(axises(2)),angles(2),rxx1,ryy1,rzz1,rxx2,ryy2,rzz2)
   call cartesian2lonlat(rxx2,ryy2,rzz2,rrs,rlon,rlat)
!
   if ( (abs(sqrt( xx **2+ yy **2+ zz **2)-rr).gt.zero_thr).or.                &
        (abs(sqrt(rxx1**2+ryy1**2+rzz1**2)-rr).gt.zero_thr).or.                &
        (abs(sqrt(rxx2**2+ryy2**2+rzz2**2)-rr).gt.zero_thr) ) then
     print '(a,a,x,2(f7.4,x))','* check rotation (axis: ',trim(axis)//' ) = ',lon0*rad2deg,lat0*rad2deg
     print '(a,4(f7.4,x))','  - before: r, x, y, z : ',sqrt(xx**2+yy**2+zz**2),xx,yy,zz
     print '(a,4(f7.4,x))','  - interm: r, x, y, z : ',sqrt(rxx1**2+ryy1**2+rzz1**2),rxx1,ryy1,rzz1
     print '(a,4(f7.4,x))','  - after : r, x, y, z : ',sqrt(rxx2**2+ryy2**2+rzz2**2),rxx2,ryy2,rzz2
     print *,' '
   endif
!
   return
   end subroutine rotation_lonlat_ll0_to_axis
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine inv_rotation_lonlat_ll0_to_axis(lon0, lat0, axis, rs, lons, lats, rrs, rlons, rlats)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),         intent(in   ) :: lon0, lat0
   character(len=*), intent(in   ) :: axis
   real(r8),         intent(in   ) :: rs, lons, lats
   real(r8),         intent(  out) :: rrs, rlons, rlats
! local variables
   real(r8) :: xx, yy, zz
   real(r8) :: rxx1, ryy1, rzz1, rxx2, ryy2, rzz2
   character(len=1), dimension(2) :: axises
   real(r8),         dimension(2) :: angles
!
   if (trim(axis)=='x') then
     axises(1) = 'y'; angles(1) =-lat0
     axises(2) = 'z'; angles(2) = lon0
   elseif (trim(axis)=='y') then
     axises(1) = 'x'; angles(1) = lat0
     axises(2) = 'z'; angles(2) =-pi/2.0_r8+lon0
   elseif (trim(axis)=='z') then
     axises(1) = 'y'; angles(1) =-lat0+(pi/2.0_r8)
     axises(2) = 'z'; angles(2) = lon0
   else
     print*,'check axis...(in inv_rotation_lonslats...) '
     stop
   endif
!
   call lonlat2cartesian(rs,lons,lats,xx,yy,zz)
   call rotation(trim(axises(1)),angles(1),xx,yy,zz,rxx1,ryy1,rzz1)
   call rotation(trim(axises(2)),angles(2),rxx1,ryy1,rzz1,rxx2,ryy2,rzz2)
   call cartesian2lonlat(rxx2,ryy2,rzz2,rrs,rlons,rlats)
!
   return
   end subroutine inv_rotation_lonlat_ll0_to_axis
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine rotation_lonslats_ll0_to_axis(lon0, lat0, axis, npts, rs, lons, lats, rrs, rlons, rlats)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),                  intent(in   ) :: lon0, lat0
   character(len=*),          intent(in   ) :: axis
   integer(i4),               intent(in   ) :: npts
   real(r8), dimension(npts), intent(in   ) :: rs, lons, lats
   real(r8), dimension(npts), intent(inout) :: rrs, rlons, rlats
! local variables
   integer(i4) :: i
!
   do i = 1,npts
     call rotation_lonlat_ll0_to_axis(lon0,lat0,trim(axis),rs(i),lons(i),lats(i),rrs(i),rlons(i),rlats(i))
   enddo
!
   return
   end subroutine rotation_lonslats_ll0_to_axis
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine inv_rotation_lonslats_ll0_to_axis(lon0, lat0, axis, npts, rs, lons, lats, rrs, rlons, rlats)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),                  intent(in   ) :: lon0, lat0
   character(len=*),          intent(in   ) :: axis
   integer(i4),               intent(in   ) :: npts
   real(r8), dimension(npts), intent(in   ) :: rs, lons, lats
   real(r8), dimension(npts), intent(inout) :: rrs, rlons, rlats
! local variables
   integer(i4) :: i
!
   do i = 1,npts
     call inv_rotation_lonlat_ll0_to_axis(lon0,lat0,axis,rs(i),lons(i),lats(i),rrs(i),rlons(i),rlats(i))
   enddo
!
   return
   end subroutine inv_rotation_lonslats_ll0_to_axis
!-------------------------------------------------------------------------------
!
! fail!! what problems!?
!-------------------------------------------------------------------------------
   subroutine rotation_lonlat_xaxis_to_ll0(lon0, lat0, rr, lon, lat, rrr, rlon, rlat)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),         intent(in   ) :: lon0, lat0
   real(r8),         intent(in   ) :: rr, lon, lat
   real(r8),         intent(inout) :: rrr, rlon, rlat
! local variables
   real(r8)    ::  xx,  yy,  zz ! initial vector
   real(r8)    :: rxx, ryy, rzz
!
   call lonlat2cartesian(rr,lon,lat,xx,yy,zz)
!   rxx = cos(lat0)*cos(lon0)*xx+cos(lat0)*sin(lon0)*yy+sin(lat0)*zz
!   ryy =-sin(lon0)*xx+cos(lon0)*yy
!   rzz =-sin(lat0)*cos(lon0)*xx-sin(lat0)*sin(lon0)*yy+cos(lat0)*zz
   rxx = cos(lat0)*cos(lon0)*xx-sin(lon0)*yy-sin(lat0)*cos(lon0)*zz
   ryy = cos(lat0)*sin(lon0)*xx+cos(lon0)*yy-sin(lat0)*sin(lon0)*zz
   rzz = sin(lat0)*xx+cos(lat0)*zz
   call cartesian2lonlat(rxx,ryy,rzz,rrr,rlon,rlat)
!
   return
   end subroutine rotation_lonlat_xaxis_to_ll0
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine lonlat_to_great_circle_angle(lon, lat, alpha, beta)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: lon, lat
   real(r8), intent(inout) :: alpha, beta
! local variables
   real(r8) :: x, y, z, tmp
!
! lon -> alpha
   if (lon.ge.0.0_r8.and.lon.le.pi/2.0_r8) then
     alpha = lon
   elseif (lon.ge.1.5_r8*pi.and.lon.le.2.0_r8*pi) then
     alpha = 2.0_r8*pi-lon
   else
     print*,'check range of lon...(in lonlattogreatecircleangles) ...'
   endif
!
! lon,lat -> beta
   call lonlat2cartesian(1.0_r8,lon,lat,x,y,z)
   tmp = get_lon(x,z)
!
   if (tmp.ge.0.0_r8.and.tmp.le.pi/2.0_r8) then
     beta = tmp
   elseif (tmp.ge.1.5_r8*pi.and.tmp.le.2.0_r8*pi) then
     beta = 2.0_r8*pi-tmp
   else
     print*,'check range of tmp...(in lonlattogreatecircleangles) ...'
   endif
!
   return
   end subroutine lonlat_to_great_circle_angle
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine great_circle_angle_to_lonlat(alpha, beta, lon, lat)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: alpha, beta
   real(r8), intent(inout) :: lon, lat
! local variables
   real(r8) :: x, y, z, tmp
!
! alpha -> lon
   if (alpha.ge.0.0_r8.and.alpha.le.pi/2.0_r8) then
     lon = alpha
   elseif (alpha.ge.-pi/2.0_r8.and.alpha.lt.0.0_r8) then
     lon = 2.0_r8*pi+alpha
   else
     print*,'check range of alpha...(in greatecircleanglestolonlat) ...'
   endif
! alpha,beta -> lat
   if (beta.ge.0.0_r8.and.beta.le.pi/2.0_r8) then
   elseif (beta.ge.-pi/2.0_r8.and.beta.lt.0.0_r8) then
   else
     print*,'check range of beta...(in greatecircleanglestolonlat) ...'
   endif
!
   x = sqrt(1.0_r8/(1.0_r8+tan(alpha)**2.0_r8+tan(beta)**2.0_r8))
   y = x*tan(alpha)
   z = x*tan(beta)
   call cartesian2lonlat(x,y,z,tmp,lon,lat)
!
   return
   end subroutine great_circle_angle_to_lonlat
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function getquadrants(npts, radians) result(quads)
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
     quads(i) = getquadrant(radians(i))
   enddo
!
   end function getquadrants
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function getquadrant(radian) result(quad)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: radian
   integer(i4) :: quad

   if (radian>= 0.0_r8.and.radian<0.5_r8*pi) then
     quad = 1
   elseif (radian>= 0.5_r8*pi.and.radian<pi) then
     quad = 2
   elseif (radian>= pi.and.radian<1.5_r8*pi) then
     quad = 3
   elseif (radian>= 1.5_r8*pi.and.radian<2.0_r8*pi) then
     quad = 4
   endif
!
   end function getquadrant
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_lon(x, y) result(theta)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: x, y
   real(r8) :: theta

   if (x==0.0_r8.and.y==0.0_r8) then
     theta = 0.0_r8
   ! 1
   elseif (x>0.0_r8.and.y>= 0.0_r8.or.x==0.0_r8.and.y>0.0_r8) then
     theta = atan(y/x)
   ! 2
   elseif (x<0.0_r8.and.y>0.0_r8) then ! atan(y/x) =-theta
     theta = pi+atan(y/x)
   ! 3
   elseif (x<0.0_r8.and.y<=0.0_r8) then ! atan(y/x)
     theta = pi+atan(y/x)
   ! 4
   elseif (x>= 0.0_r8.and.y<0.0_r8) then
     theta = 2.0_r8*pi+atan(y/x)
   ! other ?
   else
     print*, 'check get_lon...'
     stop
   endif
!
   if (theta.ge.2.0_r8*pi) theta = theta-2.0_r8*pi
!
   end function get_lon
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_lons_twovector_to_xyplane_tmp(r1, lon1, lat1, r2, lon2, lat2, theta1, theta2)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: r1, lon1, lat1, r2, lon2, lat2
   real(r8), intent(inout) :: theta1, theta2
! local variables
   real(r8) :: x1, y1, z1, x2, y2, z2 ! two vectors
   real(r8) :: rx1, ry1, rz1, rx2, ry2, rz2 ! two rotated vectors
   real(r8) :: nx, ny, nz, nr, nlon, nlat ! normal vector
!
   call lonlat2cartesian(r1,lon1,lat1,x1,y1,z1)
   call lonlat2cartesian(r2,lon2,lat2,x2,y2,z2)
!
   nx = y1*z2-z1*y2
   ny = z1*x2-x1*z2
   nz = x1*y2-y1*x2
   call cartesian2lonlat(nx,ny,nz,nr,nlon,nlat)
!
   call rotation('z',-nlon,x1,y1,z1,rx1,ry1,rz1)
   call rotation('y',nlat-pi/2.0_r8,rx1,ry1,rz1,rx1,ry1,rz1)
   call rotation('z',-nlon,x2,y2,z2,rx2,ry2,rz2)
   call rotation('y',nlat-pi/2.0_r8,rx2,ry2,rz2,rx2,ry2,rz2)
!
   if (abs(rz1).gt.1.0d-15.or.abs(rz2).gt.1.0d-15) then
     print*,'some error... ',rz1,rz2
   endif
!
   theta1 = get_lon(rx1,ry1)
   theta2 = get_lon(rx2,ry2)
!
   return
   end subroutine get_lons_twovector_to_xyplane_tmp
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_lons_twovector_to_xyplane(r1, lon1, lat1, r2, lon2, lat2, theta1, theta2)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: r1, lon1, lat1, r2, lon2, lat2
   real(r8), intent(inout) :: theta1, theta2
! local variables
   real(r8) :: x1, y1, z1, x2, y2, z2 ! two vectors
   real(r8) :: rx1, ry1, rz1, rx2, ry2, rz2 ! two rotated vectors
   real(r8) :: nx, ny, nz, nr, nlon, nlat ! normal vector
!
   call lonlat2cartesian(r1,lon1,lat1,x1,y1,z1)
   call lonlat2cartesian(r2,lon2,lat2,x2,y2,z2)
!
   nx = y1*z2-z1*y2
   ny = z1*x2-x1*z2
   nz = x1*y2-y1*x2
   call cartesian2lonlat(nx,ny,nz,nr,nlon,nlat)
!
   call rotation('z',-nlon,x1,y1,z1,rx1,ry1,rz1)
   call rotation('y',nlat-pi/2.0_r8,rx1,ry1,rz1,rx1,ry1,rz1)
   call rotation('z',-nlon,x2,y2,z2,rx2,ry2,rz2)
   call rotation('y',nlat-pi/2.0_r8,rx2,ry2,rz2,rx2,ry2,rz2)
!
   if (abs(rz1).gt.1.0d-15.or.abs(rz2).gt.1.0d-15) then
     print*,'some error... ',rz1,rz2
   endif
!
   theta1 = get_lon(rx1,ry1)
   theta2 = get_lon(rx2,ry2)
!
   return
   end subroutine get_lons_twovector_to_xyplane
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_lonlat_twovector_to_xyplane(r1, lon1, lat1, r2, lon2, lat2, r0, lon0, lat0, r, lon, lat)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: r1, lon1, lat1, r2, lon2, lat2, r0, lon0, lat0
   real(r8), intent(inout) :: r, lon, lat
! local variables
   real(r8) :: x0, y0, z0, x, y, z, rr
   real(r8) :: x1, y1, z1, x2, y2, z2 ! two vectors
   real(r8) :: rx1, ry1, rz1, rx2, ry2, rz2 ! two rotated vectors
   real(r8) :: nx, ny, nz, nr, nlon, nlat ! normal vector
!
   call lonlat2cartesian(r1,lon1,lat1,x1,y1,z1)
   call lonlat2cartesian(r2,lon2,lat2,x2,y2,z2)
!
   nx = y1*z2-z1*y2
   ny = z1*x2-x1*z2
   nz = x1*y2-y1*x2
   call cartesian2lonlat(nx,ny,nz,nr,nlon,nlat)
!
   call lonlat2cartesian(r0,lon0,lat0,x0,y0,z0)
!
   call rotation('y',-nlat+pi/2.0_r8,x0,y0,z0,x,y,z)
   call rotation('z',+nlon,x,y,z,x,y,z)
!
   call cartesian2lonlat(x,y,z,r,lon,lat)
!
   return
   end subroutine get_lonlat_twovector_to_xyplane
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
!
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
   real(r8)    :: x0, y0, z0, x1, y1, z1, x2, y2, z2, cx, cy, cz, area(npts)
!
   isin = .true.
!
   do i = 1,npts
     call lonlat2cartesian(1.0_r8,lon0,lat0,x0,y0,z0)
     call lonlat2cartesian(1.0_r8,lons(i),lats(i),x1,y1,z1)
     if (i.ne.npts) then
       call lonlat2cartesian(1.0_r8,lons(i+1),lats(i+1),x2,y2,z2)
     else
       call lonlat2cartesian(1.0_r8,lons(1),lats(1),x2,y2,z2)
     endif
     ! area = (x0,y0,z0)dot(cx,cy,cz) = (x0,y0,z0)dot{(x1,y1,z1)x(x2,y2,z2)}
     cx = y1*z2-z1*y2
     cy = z1*x2-x1*z2
     cz = x1*y2-y1*x2
     area(i) = x0*cx+y0*cy+z0*cz
     if (area(i).lt.0.0_r8) then
       isin = .false.
     endif
   enddo
!
   if (.not.isin) then
     print *,'area = ',area
     print *,'x1,y1,z1 = ',x1,y1,z1
     print *,'x2,y2,z2 = ',x2,y2,z2
     print *,'cx,cy,cz = ',cx,cy,cz
   endif
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
   real(r8), dimension(npts) :: x, y, z, rx, ry, rz
   real(r8) :: x0, y0, z0, r0, rx0, ry0, rz0
   real(r8) :: rlon0, rlat0, min4q
!
! lon,lat -> x,y,z
   do i = 1,npts
     x(i) = r*cos(rlats(i))*cos(rlons(i))
     y(i) = r*cos(rlats(i))*sin(rlons(i))
     z(i) = r*sin(rlats(i))
   enddo
!
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
   rlon0 = get_lon(x0,y0) ! atan(y0/x0)
   rlat0 = asin(z0/r0)
!
! rotation  z(-lon0),y(lat0-pi/2)
   call rotations('z',-rlon0,npts,x,y,z,rx,ry,rz)
   x = rx; y = ry; z = rz
   call rotations('y',rlat0-pi/2.0_r8,npts,x,y,z,rx,ry,rz)
!
   do i = 1,npts
     thetas(i) = get_lon(rx(i),ry(i))
   enddo
!
   end function get_thetas
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
   end module coordinates
