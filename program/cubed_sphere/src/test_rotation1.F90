!-------------------------------------------------------------------------------
   program test
!
   use kinds, only : i4, l4, r4, r8
   use coordinates, only : cartesian2lonlat, lonlat2cartesian, rotation
!
   implicit none
!
   integer(i4) :: i
   real(r8), parameter :: pi = acos(-1.0_r8)
   real(r8) :: r, lon, lat, rr, rlon, rlat
   real(r8) :: x, y, z, rx, ry, rz
!
   r = 1.0_r8
   lon = pi/2.0_r8
   lat = 0.0_r8
!
   call lonlat2cartesian(r,lon,lat,x,y,z)
   print*,'sphere : ',r,lon,lat
   print*,'cartesian : ',x,y,z
!
   call rotation('z',-lon,x,y,z,rx,ry,rz)
   call rotation('y',lat-pi/2.0_r8,rx,ry,rz,rx,ry,rz)
   call cartesian2lonlat(rx,ry,rz,rr,rlon,rlat)
   print*,'sphere : ',rr,rlon,rlat
   print*,'cartesian : ',rx,ry,rz
!
   end program test
!-------------------------------------------------------------------------------
