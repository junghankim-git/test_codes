   module buffer
   implicit none
!
   integer, parameter :: buf_i4 = 1
   integer, parameter :: buf_i8 = 2
   integer, parameter :: buf_r4 = 3
   integer, parameter :: buf_r8 = 4
!
   type buffer_t
     integer    :: typ
     integer    :: n
     integer(4), dimension(:), allocatable :: i4
     integer(8), dimension(:), allocatable :: i8
     real(4),    dimension(:), allocatable :: r4
     real(8),    dimension(:), allocatable :: r8
   end type buffer_t
!
   interface buf_add
     module procedure buf_add_i8
   end interface buf_add
!
   public :: buffer_t, buf_i4, buf_i8, buf_r4, buf_r8
   public :: buf_initialize, buf_finalize, buf_add, buf_add_i8
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   contains
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine buf_initialize(buf, typ)
   implicit none
   type(buffer_t), intent(inout) :: buf
   integer,        intent(in   ) :: typ
!
   buf%typ = typ
   buf%n   = 0
!
   end subroutine buf_initialize
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine buf_finalize(buf)
   implicit none
   type(buffer_t), intent(inout) :: buf
!
   buf%n = 0
   if     (buf%typ.eq.buf_i4) then
     if (allocated(buf%i4)) deallocate(buf%i4)
   elseif (buf%typ.eq.buf_i8) then
     if (allocated(buf%i8)) deallocate(buf%i8)
   elseif (buf%typ.eq.buf_r4) then
     if (allocated(buf%r4)) deallocate(buf%r4)
   elseif (buf%typ.eq.buf_r8) then
     if (allocated(buf%r8)) deallocate(buf%r8)
   endif
!
   end subroutine buf_finalize
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine buf_add_i8(buf, num)
   implicit none
   type(buffer_t), intent(inout) :: buf
   integer(8),     intent(in   ) :: num
!
   integer(8), dimension(:), allocatable :: tmp
!
   if (buf%typ.ne.buf_i8) then
     print *, 'check buffer type... '
     stop
   endif
!
   if (buf%n.eq.0) then
     buf%n = 1
     allocate(buf%i8(1))
     buf%i8(1) = num
   else
     buf%n = buf%n+1
     allocate(tmp(buf%n))
     tmp(1:buf%n-1) = buf%i8(:)
     tmp(buf%n)     = num
     deallocate(buf%i8)
     allocate(buf%i8(buf%n))
     buf%i8(:) = tmp(:)
     deallocate(tmp)
   endif
!
   end subroutine buf_add_i8
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   end module buffer
