!
   program degrees2radians_main
   use kinds, only : i4, r4, r8
   use scrip_input, only : scrip_t
   use scrip_input, only : scrip_initialize, scrip_initialize_read, scrip_write, scrip_finalize, adjust_corners
   use scrip_input, only : degrees2radians, radians2degrees
!
   type(scrip_t) :: scrip_in
!
   call scrip_initialize_read(scrip_in, '../cs/syjung_ne030_unrotated_degrees.nc')
!
   call degrees2radians(scrip_in)
!
   scrip_in%filename = '../cs/syjung_ne030_unrotated.nc'
   call scrip_write(scrip_in)
!
   call scrip_finalize(scrip_in)
!
   end program degrees2radians_main
