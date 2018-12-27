!-------------------------------------------------------------------------------
   program test
!
   use kinds, only : i4, l4, r4, r8
   use coordinates, only : cartesian2lonlat, lonlat2cartesian, rotation
   use coordinates, only : get_lons_twovectortoxyplane, get_lonlat_twovectortoxyplane
!
   implicit none
!
   integer(i4) :: i
   real(r8), parameter :: pi = acos(-1.0_r8)
   real(r8) :: r1, lon1, lat1, r2, lon2, lat2
   real(r8) :: r, delta, frac
   real(r8), dimension(4) :: u, theta, lons, lats
!
   r1 = 1.0_r8
   r2 = 1.0_r8
!
   lon1 = pi/3.0_r8
   lat1 = 0.0_r8
   lon2 = 2.0_r8*pi/3.0_r8
   lat2 = 0.0_r8
!
#if 0
   lon1 = 0.0_r8
   lat1 = 0.0_r8
   lon2 = pi/2.0_r8
   lat2 = pi/4.0_r8
#endif
!
   u(:) =(/-1.0_r8,-0.6_r8,0.6_r8,1.0_r8/)
!
   call get_lons_twovectortoxyplane(r1,lon1,lat1,r2,lon2,lat2,theta(1),theta(4))
!
   delta = 0.25_r8*(theta(4)-theta(1))
   theta(2) = theta(1)+delta
   theta(3) = theta(2)+delta
   frac = 0.5_r8*(theta(4)-theta(1))
!
   print*,'start = ',lon1,lat1
   print*,'end = ',lon2,lat2
   do i = 1,4
     call get_lonlat_twovectortoxyplane(r1,lon1,lat1,r2,lon2,lat2,1.0_r8,theta(i),0.0_r8,r,lons(i),lats(i))
     print*,theta(i),lons(i),lats(i)
   enddo
!
   end program test
!-------------------------------------------------------------------------------
