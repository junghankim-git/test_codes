!-------------------------------------------------------------------------------
   program test
!
   use kinds,       only: i4, l4, r4, r8
   use coordinates, only: rotation_lonlat_xaxis_to_ll0
!
   implicit none
!
   real(r8), parameter :: pi = acos(-1.0_r8)
   real(r8) :: r0, lon0, lat0
   real(r8) :: rr, lon, lat, rrr, rlon, rlat
   integer(i4) :: i
!
   r0    = 1.0_r8
   lon0  = pi/4.0_r8
   lat0  = pi/4.0_r8
   rr    = 1.0_r8
   lon   = pi/2.0_r8
   lat   = 0.0_r8
!
   call print_status('original',rr,lon,lat)
   call rotation_lonlat_xaxis_to_ll0(lon0,lat0,rr,lon,lat,rrr,rlon,rlat)
   call print_status('new     ',rrr,rlon,rlat)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_status(message, rr, lon, lat)
   implicit none
   character(len=*), intent(in   ) :: message
   real(r8),         intent(in   ) :: rr, lon, lat
!
   print '(x,a) ',trim(message)
   print '(x,x,f6.2,x,f6.2,a,x,f6.2,a) ',rr,lon/pi,' pi',lat/pi,' pi'
!
   end subroutine print_status
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program test
!-------------------------------------------------------------------------------
