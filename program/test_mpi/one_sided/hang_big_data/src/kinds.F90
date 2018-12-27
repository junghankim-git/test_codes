!-------------------------------------------------------------------------------
!
!  abstract :
!    Module for appllication-wide base parameters
!
!  history log :
!    2012-12-27   junghan kim    HOMME (CAM-SE) kinds
!
!-------------------------------------------------------------------------------
   module kinds
!
   private
!
   integer(kind=4), parameter, public :: i4  = selected_int_kind( 6)
   integer(kind=4), parameter, public :: i8  = selected_int_kind(13)
   integer(kind=4), parameter, public :: l4  = kind(.true.)
   integer(kind=4), parameter, public :: r4  = selected_real_kind( 6)
   integer(kind=4), parameter, public :: r8  = selected_real_kind(12)
   integer(kind=4), parameter, public :: r8d = selected_real_kind(12)
   integer(kind=4), parameter, public :: r16 = r8d
   integer(kind=4), parameter, public :: max_str_len = 80
   integer(kind=4), parameter, public :: log_unit = 6
   integer(kind=4), parameter, public :: std_in   = 5
   integer(kind=4), parameter, public :: std_out  = 6
   integer(kind=4), parameter, public :: std_err  = 0
!-------------------------------------------------------------------------------
   end module kinds
!-------------------------------------------------------------------------------
