!-------------------------------------------------------------------------------
   program main
   use kinds,      only: i4, r4, r8, l4
   use nc_control, only: netcdf_t, nc_initialize, nc_define, nc_compress, nc_finalize
!
   type(netcdf_t) :: nc
!
   call nc_initialize(nc,'../data/UP-20110725220000-000010.nc','./out.nc',1) 
!
   call nc_define(nc)
   call nc_compress(nc)
!
   call nc_finalize(nc)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine dummy
   implicit none
!
!
   end subroutine dummy
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program main
!-------------------------------------------------------------------------------
