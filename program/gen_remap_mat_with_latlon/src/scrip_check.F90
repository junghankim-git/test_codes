!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  your   name    code comment
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds,        only: i4, l4, r4, r8
   use cubed_sphere, only: cubed_sphere_t, cs_do_check
   use remap_matrix, only: rm_t, rm_initialize, rm_finalize, rm_check, rm_rewrite
!
   type(cubed_sphere_t) :: cs
   type(rm_t)           :: rmtrx
   integer(i4) :: np1, ne1, nprocs1
   integer(i4) :: np2, ne2, nprocs2
   logical(l4) :: rotated1, rotated2
   logical(l4) :: dogenremap = .false.
   character(len=512) :: filename1, filename2
   character(len=16) :: lonname1, latname1, lonname2, latname2
   character(len=512) :: scripfile1, scripfile2
   character(len=512) :: remapfile12, remapfile21
!
   call readnl()
   call cs_do_check(cs,np1,ne1,nprocs1,'result/scrip_1.nc','scrip_1.nc')
!
   call rm_initialize(rmtrx,trim(remapfile12))
   call rm_check(rmtrx,1.0d-11)
   call rm_finalize(rmtrx)
!
   call rm_initialize(rmtrx,trim(remapfile21))
   call rm_check(rmtrx,1.0d-11)
   call rm_finalize(rmtrx)
!
   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine readnl()
!-------------------------------------------------------------------------------
   implicit none
!
   namelist /inputs/ np1,ne1,nprocs1,rotated1,filename1,lonname1,latname1,     &
                     np2,ne2,nprocs2,rotated2,filename2,lonname2,latname2,     &
                     dogenremap,scripfile1,scripfile2,remapfile12,remapfile21
   open(31,file='inputs.nl',status='old')
   read(31,inputs)
   close(31)
!
   end subroutine readnl
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end program test
