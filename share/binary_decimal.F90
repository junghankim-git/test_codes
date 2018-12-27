#define ROUNDOFF_BIT
!#undef  ROUNDOFF_BIT
!#define ZERO_EXPO -32
#define ZERO_EXPO -40
!-------------------------------------------------------------------------------
   module binary_decimal
!-------------------------------------------------------------------------------
!
!  abstract : converter between binary and decimal
!
!  histroy log :
!    2015-08-20   junghan kim   initial setup
!
!  variable :
!
!-------------------------------------------------------------------------------
   use kinds, only: l4, i4, r4, r8=>r8d
   private
!
!   integer, parameter, public :: r4 = 4
!   integer, parameter, public :: r8 = 8
   integer, parameter, public :: nbit32 = 24
   integer, parameter, public :: nbit64 = 53
   integer, parameter, public :: nexp32 = 10
   integer, parameter, public :: nexp64 = 18
!
   type :: binary32_t
     logical(l4)                   :: sig
     integer(i4)                   :: exp
     real(r4)                      :: fra
     integer(1), dimension(nbit32) :: bit
   end type
!
   type :: binary64_t
     logical(l4)                   :: sig
     integer(i4)                   :: exp, n
     real(r8)                      :: fra
     integer(1), dimension(nbit64) :: bit
   end type
!
   interface bin_copy
     module procedure bin32_copy
     module procedure bin64_copy
   end interface bin_copy
   interface bin_truncate
     module procedure bin32_truncate
     module procedure bin64_truncate
   end interface bin_truncate
   interface bin_add_bit
     module procedure bin32_add_bit
     module procedure bin64_add_bit
   end interface bin_add_bit
   interface bin_remove_bit
!     module procedure bin32_remove_bit
     module procedure bin64_remove_bit
   end interface bin_remove_bit
   interface bin_roundoff
     module procedure bin32_roundoff
     module procedure bin64_roundoff
   end interface bin_roundoff
   interface dec_truncate
     module procedure dec32_truncate
     module procedure dec64_truncate
   end interface dec_truncate
   interface dec_roundoff
     module procedure dec32_roundoff
     module procedure dec64_roundoff
   end interface dec_roundoff
   interface dec_to_bin
     module procedure dec_to_bin32
     module procedure dec_to_bin64
   end interface dec_to_bin
   interface bin_to_dec
     module procedure bin32_to_dec
     module procedure bin64_to_dec
   end interface bin_to_dec
   interface bin_print
     module procedure bin32_print
     module procedure bin64_print
   end interface bin_print
!
   public :: binary32_t, binary64_t,                                           &
             bin_copy, bin_truncate, bin_add_bit, bin_remove_bit,              &
             dec_truncate, bin_roundoff, dec_roundoff,                         &
             dec_to_bin, bin_to_dec, bin_print
!
   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin32_copy(bin1, bin2)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary32_t), intent(in   ) :: bin1
   type(binary32_t), intent(  out) :: bin2
! local
!
   bin2%sig    = bin1%sig
   bin2%exp    = bin1%exp
   bin2%bit(:) = bin1%bit(:)
!
   return
   end subroutine bin32_copy
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin64_copy(bin1, bin2)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary64_t), intent(in   ) :: bin1
   type(binary64_t), intent(  out) :: bin2
! local
!
   bin2%sig    = bin1%sig
   bin2%exp    = bin1%exp
   bin2%bit(:) = bin1%bit(:)
!
   return
   end subroutine bin64_copy
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin32_truncate(bin,n)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary32_t), intent(inout) :: bin
   integer(i4)     , intent(in   ) :: n
! local
   integer :: i
!
   do i = nbit32-n+1,nbit32
     bin%bit(i) = 0
   enddo
!
   return
   end subroutine bin32_truncate
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin64_truncate(bin,n)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary64_t), intent(inout) :: bin
   integer(i4)     , intent(in   ) :: n
! local
   integer :: i
!
   do i = nbit64-n+1,nbit64
     bin%bit(i) = 0
   enddo
!
   return
   end subroutine bin64_truncate
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin32_add_bit(bin, n)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary32_t), intent(inout) :: bin
   integer(i4)     , intent(in   ) :: n
! local
   integer  :: i, j
!
   bin%bit(nbit32-n+1) = bin%bit(nbit32-n+1)+1
!
   do i = nbit32-n+1,1,-1
     if (bin%bit(i).eq.2) then
       bin%bit(i) = 0
       if (i.ne.1) then
         bin%bit(i-1) = bin%bit(i-1)+1
       else
         bin%exp = bin%exp+1
         bin%bit(1:nbit32-1) = bin%bit(2:nbit32)
         bin%bit(1) = 1
       endif
     else
       exit
     endif
   enddo
!
   return
   end subroutine bin32_add_bit
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin64_add_bit(bin, n)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary64_t), intent(inout) :: bin
   integer(i4)     , intent(in   ) :: n
! local
   integer  :: i, j
!
   bin%bit(nbit64-n+1) = bin%bit(nbit64-n+1)+1
!
   do i = nbit64-n+1,1,-1
     if (bin%bit(i).eq.2) then
       bin%bit(i) = 0
       if (i.ne.1) then
         bin%bit(i-1) = bin%bit(i-1)+1
       else
         bin%exp = bin%exp+1
         bin%bit(1:nbit64-1) = bin%bit(2:nbit64)
         bin%bit(1) = 1
       endif
     else
       exit
     endif
   enddo
!
   return
   end subroutine bin64_add_bit
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin64_remove_bit(bin, n)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary64_t), intent(inout) :: bin
   integer(i4)     , intent(in   ) :: n
! local
   integer  :: i, j
!
   bin%bit(nbit64-n+1) = bin%bit(nbit64-n+1)-1
!
   do i = nbit64-n+1,1,-1
     if (bin%bit(i).eq.-1) then
       bin%bit(i) = 1
       if (i.ne.1) then
         bin%bit(i-1) = bin%bit(i-1)-1
       else
         bin%exp = bin%exp-1
         bin%bit(1:nbit64-1) = bin%bit(2:nbit64)
         bin%bit(1) = 1
       endif
     else
       exit
     endif
   enddo
!
   return
   end subroutine bin64_remove_bit
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin32_roundoff(bin,n)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary32_t), intent(inout) :: bin
   integer(i4)     , intent(in   ) :: n
! local
   real(r4) :: vsca, vint
!
   if (n.eq.0) return
   if (bin%exp.le.ZERO_EXPO) then
     bin%sig    = .true.
     bin%exp    = 0
     bin%fra    = 0.0_r4
     bin%bit(:) = 0
     return
   endif
!
#ifdef ROUNDOFF_BIT
   call bin32_truncate(bin,n-1)
   if (bin%bit(nbit32-n+1).eq.0) then
     return
   else
     call bin32_add_bit(bin,n)
   endif
#else
   vsca = bin%fra*10_r4**n
   vint = real(int(vsca,r4),r4)
   if ((vsca-vint).ge.0.5_r4) vint = vint+1.0_r4
   bin%fra = vint/10_r4**n
   call fraction_to_32bit(bin%fra,bin%bit)
#endif
!
   return
   end subroutine bin32_roundoff
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin64_roundoff(bin,n)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary64_t), intent(inout) :: bin
   integer(i4)     , intent(in   ) :: n
! local
   integer, parameter :: nc = 4
   real(r8) :: vsca, vint
!
   if (n.eq.0) return
   if (bin%exp.le.ZERO_EXPO) then
     bin%sig    = .true.
     bin%exp    = 0
     bin%fra    = 0.0_r8
     bin%bit(:) = 0
     return
   endif
!
#ifdef ROUNDOFF_BIT
   call bin64_truncate(bin,n-1)
   if (bin%bit(nbit64-n+1).eq.0) then
     return
   else
     call bin64_add_bit(bin,n)
   endif
#else
   vsca = bin%fra*10_r8**(nexp64-n+1)
   vint = real(int(vsca,r8),r8)
   if ((vsca-vint).ge.0.5_r8) vint = vint+1.0_r8
   bin%fra = vint/10_r8**(nexp64-n+1)
   call fraction_to_64bit(bin%fra,bin%bit)
#endif
!
   return
   end subroutine bin64_roundoff
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine dec32_truncate(dec,n)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4)   , intent(inout) :: dec
   integer(i4), intent(in   ) :: n
! local
   type(binary32_t) :: bin32
!
   call dec_to_bin32(dec,bin32)
   call bin32_truncate(bin32,n)
   call bin32_to_dec(bin32,dec)
!
   return
   end subroutine dec32_truncate
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine dec64_truncate(dec,n)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)   , intent(inout) :: dec
   integer(i4), intent(in   ) :: n
! local
   type(binary64_t) :: bin64
!
   call dec_to_bin64(dec,bin64)
   call bin64_truncate(bin64,n)
   call bin64_to_dec(bin64,dec)
!
   return
   end subroutine dec64_truncate
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine dec32_roundoff(dec,n)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4)   , intent(inout) :: dec
   integer(i4), intent(in   ) :: n
! local
   type(binary32_t) :: bin32
!
   call dec_to_bin32(dec,bin32)
   call bin32_roundoff(bin32,n)
   call bin32_to_dec(bin32,dec)
!
   return
   end subroutine dec32_roundoff
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine dec64_roundoff(dec,n)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)   , intent(inout) :: dec
   integer(i4), intent(in   ) :: n
! local
   type(binary64_t) :: bin64
!
   call dec_to_bin64(dec,bin64)
   call bin64_roundoff(bin64,n)
   call bin64_to_dec(bin64,dec)
!
   return
   end subroutine dec64_roundoff
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine integer_to_32bit(dec, bin, n)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4)        , intent(in   ) :: dec
   type(binary32_t), intent(  out) :: bin
   integer(i4)     , intent(  out) :: n
! local
   integer  :: i, num
   integer(1), dimension(nbit32) :: tmp
!
   tmp(:) = 0
   num = int(dec)
   if (num.eq.0) then
     n = 0
   else
     do i = 1,nbit32
       if (mod(num,2).eq.0) then
         tmp(i) = 0
       else
         tmp(i) = 1
       endif
       if (num/2.eq.0) then
         exit
       endif
       num = num/2
     enddo
     n = i
   endif
!
   bin%bit(1:n) = tmp(n:1:-1)
!
   return
   end subroutine integer_to_32bit
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine integer_to_64bit(dec, bin, n)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)        , intent(in   ) :: dec
   type(binary64_t), intent(  out) :: bin
   integer(i4)     , intent(  out) :: n
! local
   integer  :: i, num
   integer(1), dimension(nbit64) :: tmp
!
   tmp(:) = 0
   num = int(dec)
   if (num.eq.0) then
     n = 0
   else
     do i = 1,nbit64
       if (mod(num,2).eq.0) then
         tmp(i) = 0
       else
         tmp(i) = 1
       endif
       if (num/2.eq.0) then
         exit
       endif
       num = num/2
     enddo
     n = i
   endif
!
   bin%bit(1:n) = tmp(n:1:-1)
!
   return
   end subroutine integer_to_64bit
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine fraction_to_32bit(dec, bit, n)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4)                      , intent(in   ) :: dec
   !type(binary32_t),           intent(  out) :: bin
   integer(1) , dimension(nbit32), intent(  out) :: bit
   integer(i4), optional         , intent(in   ) :: n
! local
   integer(i4) :: i, m
   real(r4)    :: a, dec_tmp
!
   if (present(n)) then
     m = n
   else
     m = 0
   endif
!
   dec_tmp = dec
   bit(m+1:) = 0
   dec_tmp = dec_tmp-int(dec_tmp)
   do i = m+1,nbit32
     a = 2.0_r4*dec_tmp
     if (a.ge.1.0_r4) then
       bit(i) = 1
       dec_tmp = a-1.0_r4
     else
       bit(i) = 0
       dec_tmp = a
     endif
   enddo
!
   return
   end subroutine fraction_to_32bit
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine fraction_to_64bit(dec, bit, n)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)                      , intent(in   ) :: dec
   !type(binary64_t),           intent(  out) :: bin
   integer(1) , dimension(nbit64), intent(  out) :: bit
   integer(i4), optional         , intent(in   ) :: n
! local
   integer(i4) :: i, m
   real(r8)    :: a, dec_tmp
!
   if (present(n)) then
     m = n
   else
     m = 0
   endif
!
   dec_tmp = dec
   bit(m+1:) = 0
   dec_tmp = dec_tmp-int(dec_tmp)
   do i = m+1,nbit64
     a = 2.0_r8*dec_tmp
     if (a.ge.1.0_r8) then
       bit(i) = 1
       dec_tmp = a-1.0_r8
     else
       bit(i) = 0
       dec_tmp = a
     endif
   enddo
!
   return
   end subroutine fraction_to_64bit
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine dec_to_bin32(dec, bin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4)        , intent(in   ) :: dec
   type(binary32_t), intent(  out) :: bin
! local
   real(r4) :: a, dec_tmp
!
   if (dec.ge.0.0)then
     bin%sig = .true.
     dec_tmp  = dec
   else
     bin%sig = .false.
     dec_tmp  = -dec
   endif
   bin%bit(:) = 0
!
   bin%fra = fraction(dec_tmp)
   bin%exp = exponent(dec_tmp)
!
   call fraction_to_32bit(bin%fra,bin%bit)
!
   return
   end subroutine dec_to_bin32
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine dec_to_bin64(dec, bin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)        , intent(in   ) :: dec
   type(binary64_t), intent(  out) :: bin
! local
   real(r8) :: a, dec_tmp
!
   if (dec.ge.0.0)then
     bin%sig = .true.
     dec_tmp  = dec
   else
     bin%sig = .false.
     dec_tmp  = -dec
   endif
   bin%bit(:) = 0
!
   bin%fra = fraction(dec_tmp)
   bin%exp = exponent(dec_tmp)
!
   call fraction_to_64bit(bin%fra,bin%bit)
!
   return
   end subroutine dec_to_bin64
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin32_to_dec(bin, dec)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary32_t), intent(in   ) :: bin
   real(r4)        , intent(  out) :: dec
! local
   integer  :: i, j
!
   dec = 0.0_r4
!
   do i = 1,nbit32
     if (bin%bit(i).ne.0) then
       dec = dec+1.0_r4/(2.0_r4**real(i,r4))
     endif
   enddo
!
   dec = dec*2.0**real(bin%exp,r4)
   if (.not.bin%sig) dec = -dec
!
   return
   end subroutine bin32_to_dec
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin64_to_dec(bin, dec)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary64_t), intent(in   ) :: bin
   real(r8)        , intent(  out) :: dec
! local
   integer  :: i, j
!
   dec = 0.0_r8
!
   do i = 1,nbit64
     if (bin%bit(i).ne.0) then
       dec = dec+1.0_r8/(2.0_r8**real(i,r8))
     endif
   enddo
!
   dec = dec*2.0_r8**real(bin%exp,r8)
   if (.not.bin%sig) dec = -dec
!
   return
   end subroutine bin64_to_dec
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine dec_to_bin64_special(dec, bin)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)        , intent(in   ) :: dec
   type(binary64_t), intent(  out) :: bin
! local
   integer(i4) :: n
   real(r8)    :: a, dec_tmp
!
   if (dec.ge.0.0)then
     bin%sig = .true.
     dec_tmp  = dec
   else
     bin%sig = .false.
     dec_tmp  = -dec
   endif
   bin%bit(:) = 0
!
   call integer_to_64bit(dec_tmp,bin,n)
   call fraction_to_64bit(dec_tmp,bin%bit,n)
!
   return
   end subroutine dec_to_bin64_special
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin64_to_dec_special(bin, dec)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary64_t), intent(in   ) :: bin
   real(r8)        , intent(  out) :: dec
! local
   integer  :: i, j
!
   dec = 0.0_r8
!
! integer to binary
   j = 0
   do i = bin%n,1,-1
     if (bin%bit(i).ne.0) then
       dec = dec+2.0_r8**real(j,r8)
     endif
     j = j+1
   enddo
!
! decimal to binary
   j = 1
   do i = bin%n+1,nbit64
     if (bin%bit(i).ne.0) then
       dec = dec+1.0_r8/(2.0_r8**real(j,r8))
     endif
     j = j+1
   enddo
!
   if (.not.bin%sig) dec = -dec
!
   return
   end subroutine bin64_to_dec_special
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine cut_64_to_32(bin)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary64_t), intent(inout) :: bin
! local
   integer :: i
!
   do i = nbit32+1, nbit64
     bin%bit(i) = 0
   enddo
!
   bin%bit(nbit32) = 1
   return
   end subroutine cut_64_to_32
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine cut_64bit(bin,digit)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary64_t), intent(inout) :: bin
   integer(i4)     , intent(in   ) :: digit
! local
   integer :: i
!
   do i = nbit64-digit,nbit64
     bin%bit(i) = 0
   enddo
!
   bin%bit(nbit32) = 1
!
   return
   end subroutine cut_64bit
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin32_print(bin)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary32_t), intent(in   ) :: bin
! local
   integer :: i, j
   real(r8) :: a, dec_tmp
!
   if (bin%sig) then
     print 111,'  +0.',bin%bit(1:nbit32),'(2^',bin%exp,')'
   else
     print 111,'  -0.',bin%bit(1:nbit32),'(2^',bin%exp,')'
   endif
   111 format(a,24(i1),x,a,i4,a)
!
   return
   end subroutine bin32_print
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bin64_print(bin)
!-------------------------------------------------------------------------------
   implicit none
!
   type(binary64_t), intent(in   ) :: bin
! local
   integer :: i, j
   real(r8) :: a, dec_tmp
!
   if (bin%sig) then
     print 111,'  +0.',bin%bit(1:nbit64),'(2^',bin%exp,')'
   else
     print 111,'  -0.',bin%bit(1:nbit64),'(2^',bin%exp,')'
   endif
   111 format(a,24(i1),x,29(i1),x,a,i4,a)
!
   return
   end subroutine bin64_print
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module binary_decimal
!-------------------------------------------------------------------------------
