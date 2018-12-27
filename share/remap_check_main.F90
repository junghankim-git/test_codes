!
!-------------------------------------------------------------------------------
   program remap_check_main
   use kinds,        only: i4, l4, r8
   use remap_matrix, only: rm_t, rm_initialize, rm_finalize, rm_check, rm_rewrite
!-------------------------------------------------------------------------------
   implicit none
!
   type(rm_t) :: rmtrx
   character(len=512) :: remapfile
!
   remapfile = '/home/jhkim/work/program/gen_remap_mat_with_latlon/remap_12.nc'
   call rm_initialize(rmtrx,trim(remapfile))
   call rm_check(rmtrx,1.0d-12)
   !call rm_rewrite(rmtrx)
   call rm_finalize(rmtrx)
!
   remapfile = '/home/jhkim/work/program/gen_remap_mat_with_latlon/remap_21.nc'
   call rm_initialize(rmtrx,trim(remapfile))
   call rm_check(rmtrx,1.0d-12)
   call rm_finalize(rmtrx)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine readnl()
!-------------------------------------------------------------------------------
   implicit none
! local variables
   integer(i4) :: id, iv
!
!  read(*,nml=inputs)
!
   return
   end subroutine readnl
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program remap_check_main
