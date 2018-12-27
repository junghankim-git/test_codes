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
   use binary_decimal, only: r4, r8, binary64_t,                               &
                          dec_to_bin64, cut_64_to_32, bin_to_dec64, print_bin64
   implicit none
!
   type(binary64_t) :: bin64
   real(r8)         :: dec
!
   dec = 0.999_r8
   call dec_to_bin64(dec,bin64)
   call print_bin64(bin64)
!
   call cut_64_to_32(bin64)
   call bin_to_dec64(bin64,dec)
   print *,dec
!
   dec = 0.999
   print *,dec
   call dec_to_bin64(dec,bin64)
   call print_bin64(bin64)
!
   end program test
!-------------------------------------------------------------------------------
