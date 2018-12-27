!-------------------------------------------------------------------------------
   program test
   use kinds, only : i4, l4, r4, r8
   use coordinates, only : lonlat2cartesian, cartesian2lonlat
   implicit none
!
   real(r8), parameter :: pi = acos(-1.0_r8)
!
   real(r8) :: r0, lon0, lat0, r1, lon1, lat1
   real(r8) :: x0, y0, z0, x1, y1, z1
   real(r8) :: slon, elon, slat, elat, dlon, dlat
!
   integer(i4) :: n
   integer(i4) :: i, j
!
#if 0
   n = 10000
   slon = 0.0_r8 ; elon = 2.0_r8*pi
   slat =-0.5_r8*pi ; elat = 0.5_r8*pi
!
   dlon =(elon-slon)/real(n,8)
   dlat =(elat-slat)/real(n,8)
!
   r0 = 1.0_r8
   lat0 = slat
   do j = 1,n
     lon0 = slon
     do i = 1,n
       call lonlat2cartesian(r0,lon0,lat0,x0,y0,z0)
       call cartesian2lonlat(x0,y0,z0,r1,lon1,lat1)
       call lonlat2cartesian(r1,lon1,lat1,x1,y1,z1)
       call check_status(r0,lon0,lat0,r1,lon1,lat1,x0,y0,z0,x1,y1,z1)
       lon0 = lon0+dlon
     enddo
     lat0 = lat0+dlat
   enddo
#endif
   x0 = 0.10374337955698959d0
   y0 =-0.52157501394867634d0
   z0 =-3.9617518417372022d0
   print*,'x0,y0,z0 = ',x0,y0,z0
   call cartesian2lonlat(x0,y0,z0,r0,lon0,lat0)
   print*,'r0,lon0,lat0 = ',r0,lon0,lat0
   call lonlat2cartesian(r0,lon0,lat0,x1,y1,z1)
   print*,'x1,y1,z1 = ',x1,y1,z1
   call cartesian2lonlat(x1,y1,z1,r1,lon1,lat1)
   call check_status(r0,lon0,lat0,r1,lon1,lat1,x0,y0,z0,x1,y1,z1)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_status(sr, slon, slat, tr, tlon, tlat, sx, sy, sz, tx, ty, tz)
   implicit none
   real(r8), intent(in   ) :: sr, slon, slat, tr, tlon, tlat
   real(r8), intent(in   ) :: sx, sy, sz, tx, ty, tz
!
   if ((abs(sr-tr).gt.1.0d-10).or.(abs(slon-tlon).gt.1.0d-10).or.(abs(slat-tlat).gt.1.0d-10)) then
     print*,'check lon,lat...'
     print*,sr,slon,slat
     print*,tr,tlon,tlat
     print*,' '
   endif
   if ((abs(sx-tx).gt.1.0d-10).or.(abs(sy-ty).gt.1.0d-10).or.(abs(sz-tz).gt.1.0d-10)) then
     print*,'check x,y,z...'
     print*,sx,sy,sz
     print*,tx,ty,tz
     print*,' '
   endif
!
   end subroutine check_status
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program test
!-------------------------------------------------------------------------------
