!-------------------------------------------------------------------------------
   program test
!
   use kinds, only : i4, l4, r4, r8
   use coordinates, only : rotation_lonslats_ll0_to_axis
!
   implicit none
!
   real(r8), parameter :: pi = acos(-1.0_r8)
   real(r8) :: r0, lon0, lat0
   integer(i4), parameter :: npts = 2
   real(r8), dimension(npts) :: rs, lons, lats
   real(r8), dimension(npts) :: rrs, rlons, rlats
   integer(i4) :: i
!
   r0 = 1.0_r8
   lon0 = pi/4.0_r8
   lat0 = pi/4.0_r8
   rs(:) = 1.0_r8
! test1
  !lons(1) = pi/4.0_r8; lats(1) = pi/4.0_r8
  !lons(2) =    0.0_r8; lats(2) = pi/2.0_r8
! test2
   lons(1) = 7.0_r8*pi/4.0_r8; lats(1) =-pi/4.0_r8
   lons(2) = 0.0_r8; lats(2) = 0.0_r8
!
   call printstatus('original',npts,rs,lons,lats)
!
   call rotation_lonslats_ll0_to_axis(lon0,lat0,'x',npts,rs,lons,lats,rrs,rlons,rlats)
   call printstatus('x axis',npts,rrs,rlons,rlats)
   call rotation_lonslats_ll0_to_axis(lon0,lat0,'y',npts,rs,lons,lats,rrs,rlons,rlats)
   call printstatus('y axis',npts,rrs,rlons,rlats)
   call rotation_lonslats_ll0_to_axis(lon0,lat0,'z',npts,rs,lons,lats,rrs,rlons,rlats)
   call printstatus('z axis',npts,rrs,rlons,rlats)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine printstatus(message, n, rs, lons, lats)
   implicit none
   character(len=*), intent(in   ) :: message
   integer(i4), intent(in   ) :: n
   real(r8), dimension(n), intent(in   ) :: rs, lons, lats
!
   print '(x,a) ',trim(message)
   do i = 1,npts
     print '(x,x,i2,x,f6.2,x,f6.2,a,x,f6.2,a) ',i,rs(i),lons(i)/pi,' pi',lats(i)/pi,' pi'
   enddo
!
   end subroutine printstatus
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program test
!-------------------------------------------------------------------------------
