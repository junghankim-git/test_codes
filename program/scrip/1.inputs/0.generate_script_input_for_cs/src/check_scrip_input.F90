!
   program checker
   use kinds, only : i4, r4, r8
   use scrip_input, only : scrip_t
   use scrip_input, only : scrip_initialize, scrip_initialize_read, scrip_write, scrip_finalize, scrip_check
!
   type(scrip_t) :: scrip_in
!
print *, 'step 1'
   !call scrip_initialize_read(scrip_in,'/home/jhkim/study/library/main/scrip/0.inputs/CS/syjung_ne030_unrotated.nc')
   call scrip_initialize_read(scrip_in,'/home/jhkim/work/program/gen_remap_mat_with_latlon/scrip_2.nc')
print *, 'step 2'
!
   call scrip_check(scrip_in)
!
!  scrip_in%filename = '../CS/syjung_ne030_unrotated.nc'
!  call scrip_write(scrip_in)
!
   call scrip_finalize(scrip_in)
!
   end program checker
