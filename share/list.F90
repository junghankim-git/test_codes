!------------------------------------------------------------------------------
   module list
!-------------------------------------------------------------------------------
!
!  abstract : 
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer, parameter :: maxlen_str = 128
!
   integer, parameter :: list_l4 = 1
   integer, parameter :: list_i4 = 2
   integer, parameter :: list_i8 = 3
   integer, parameter :: list_r4 = 4
   integer, parameter :: list_r8 = 5
   integer, parameter :: list_st = 6
!
   type list_t
     integer    :: typ
     integer    :: n
     logical(4), dimension(:), allocatable :: l4
     integer(4), dimension(:), allocatable :: i4
     integer(8), dimension(:), allocatable :: i8
     real(4),    dimension(:), allocatable :: r4
     real(8),    dimension(:), allocatable :: r8
     character(len=maxlen_str), dimension(:), allocatable :: st
   end type list_t
!
   interface list_add
     module procedure list_add_l4
     module procedure list_add_i4
     module procedure list_add_i8
     module procedure list_add_r4
     module procedure list_add_r8
     module procedure list_add_st
   end interface list_add
   interface list_change
     module procedure list_change_l4
     module procedure list_change_i4
     module procedure list_change_i8
     module procedure list_change_r4
     module procedure list_change_r8
     module procedure list_change_st
   end interface list_change
!
   public :: list_t, list_l4, list_i4, list_i8, list_r4, list_r8, list_st
   public :: list_initialize, list_finalize, list_add, list_change
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   contains
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_initialize(list, typ)
   implicit none
   type(list_t), intent(inout) :: list
   integer,      intent(in   ) :: typ
! local variables
   character(len=128) :: message
!
   if (.not.check_type(typ)) then
     write(message,'(a,i2)') 'list_initialize: check type : ', typ
   endif
   list%typ = typ
   list%n   = 0
!
   end subroutine list_initialize
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_finalize(list)
   implicit none
   type(list_t), intent(inout) :: list
!
   list%n = 0
   if     (list%typ.eq.list_i4) then
     if (allocated(list%i4)) deallocate(list%i4)
   elseif (list%typ.eq.list_i8) then
     if (allocated(list%i8)) deallocate(list%i8)
   elseif (list%typ.eq.list_r4) then
     if (allocated(list%r4)) deallocate(list%r4)
   elseif (list%typ.eq.list_r8) then
     if (allocated(list%r8)) deallocate(list%r8)
   elseif (list%typ.eq.list_st) then
     if (allocated(list%st)) deallocate(list%st)
   endif
!
   end subroutine list_finalize
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_add_l4(list, val)
   implicit none
   type(list_t), intent(inout) :: list
   logical(4),   intent(in   ) :: val
!
   logical(4), dimension(:), allocatable :: tmp
!
   if (list%typ.ne.list_l4) then
     print *, 'check list type... ', list%typ, list_l4
     stop
   endif
!
   if (list%n.eq.0) then
     list%n = 1
     allocate(list%l4(1))
     list%l4(1) = val
   else
     list%n = list%n+1
     allocate(tmp(list%n))
     tmp(1:list%n-1) = list%l4(:)
     tmp(list%n)     = val
     deallocate(list%l4)
     allocate(list%l4(list%n))
     list%l4(:) = tmp(:)
     deallocate(tmp)
   endif
!
   end subroutine list_add_l4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_add_i4(list, val)
   implicit none
   type(list_t), intent(inout) :: list
   integer(4),   intent(in   ) :: val
!
   integer(4), dimension(:), allocatable :: tmp
!
   if (list%typ.ne.list_i4) then
     print *, 'check list type... ', list%typ, list_i4
     stop
   endif
!
   if (list%n.eq.0) then
     list%n = 1
     allocate(list%i4(1))
     list%i4(1) = val
   else
     list%n = list%n+1
     allocate(tmp(list%n))
     tmp(1:list%n-1) = list%i4(:)
     tmp(list%n)     = val
     deallocate(list%i4)
     allocate(list%i4(list%n))
     list%i4(:) = tmp(:)
     deallocate(tmp)
   endif
!
   end subroutine list_add_i4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_add_i8(list, val)
   implicit none
   type(list_t), intent(inout) :: list
   integer(8),   intent(in   ) :: val
!
   integer(8), dimension(:), allocatable :: tmp
!
   if (list%typ.ne.list_i8) then
     print *, 'check list type... ', list%typ, list_i8
     stop
   endif
!
   if (list%n.eq.0) then
     list%n = 1
     allocate(list%i8(1))
     list%i8(1) = val
   else
     list%n = list%n+1
     allocate(tmp(list%n))
     tmp(1:list%n-1) = list%i8(:)
     tmp(list%n)     = val
     deallocate(list%i8)
     allocate(list%i8(list%n))
     list%i8(:) = tmp(:)
     deallocate(tmp)
   endif
!
   end subroutine list_add_i8
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_add_r4(list, val)
   implicit none
   type(list_t), intent(inout) :: list
   real(4),      intent(in   ) :: val
!
   real(4), dimension(:), allocatable :: tmp
!
   if (list%typ.ne.list_r4) then
     print *, 'check list type... ', list%typ, list_r4
     stop
   endif
!
   if (list%n.eq.0) then
     list%n = 1
     allocate(list%r4(1))
     list%r4(1) = val
   else
     list%n = list%n+1
     allocate(tmp(list%n))
     tmp(1:list%n-1) = list%r4(:)
     tmp(list%n)     = val
     deallocate(list%r4)
     allocate(list%r4(list%n))
     list%r4(:) = tmp(:)
     deallocate(tmp)
   endif
!
   end subroutine list_add_r4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_add_r8(list, val)
   implicit none
   type(list_t), intent(inout) :: list
   real(8),      intent(in   ) :: val
!
   real(8), dimension(:), allocatable :: tmp
!
   if (list%typ.ne.list_r8) then
     print *, 'check list type... ', list%typ, list_r8
     stop
   endif
!
   if (list%n.eq.0) then
     list%n = 1
     allocate(list%r8(1))
     list%r8(1) = val
   else
     list%n = list%n+1
     allocate(tmp(list%n))
     tmp(1:list%n-1) = list%r8(:)
     tmp(list%n)     = val
     deallocate(list%r8)
     allocate(list%r8(list%n))
     list%r8(:) = tmp(:)
     deallocate(tmp)
   endif
!
   end subroutine list_add_r8
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_add_st(list, val)
   implicit none
   type(list_t),     intent(inout) :: list
   character(len=*), intent(in   ) :: val
!
   integer :: i
   character(len=maxlen_str), dimension(:), allocatable :: tmp
!
   if (list%typ.ne.list_st) then
     print *, 'check list type... ', list%typ, list_st
     stop
   endif
!
   if (list%n.eq.0) then
     list%n = 1
     allocate(list%st(1))
     list%st(1) = trim(val)
   else
     list%n = list%n+1
     allocate(tmp(list%n))
     do i=1,list%n-1
       tmp(i) = trim(list%st(i))
     enddo
     tmp(list%n) = trim(val)
!
     deallocate(list%st)
     allocate(list%st(list%n))
     do i=1,list%n
     list%st(i) = tmp(i)
     enddo
     deallocate(tmp)
print *, '1 = ', val
print *, '2 = ', tmp(:)
print *, '3 = ', list%st(:)
   endif
!
   end subroutine list_add_st
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_change_l4(list, il, val)
   implicit none
   type(list_t), intent(inout) :: list
   integer(4),   intent(in   ) :: il
   logical(4),   intent(in   ) :: val
! local variable
!
   if (list%typ.ne.list_l4) then
     print *, 'list_change_l4: check list type... '
     stop
   endif
   if (il.gt.list%n) then
     print *, 'list_change_l4: check il... '
     stop
   endif
!
   list%l4(il) = val
!
   end subroutine list_change_l4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_change_i4(list, il, val)
   implicit none
   type(list_t), intent(inout) :: list
   integer(4),   intent(in   ) :: il
   integer(4),   intent(in   ) :: val
! local variable
!
   if (list%typ.ne.list_i4) then
     print *, 'list_change_i4: check list type... '
     stop
   endif
   if (il.gt.list%n) then
     print *, 'list_change_i4: check il... '
     stop
   endif
!
   list%i4(il) = val
!
   end subroutine list_change_i4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_change_i8(list, il, val)
   implicit none
   type(list_t), intent(inout) :: list
   integer(4),   intent(in   ) :: il
   integer(8),   intent(in   ) :: val
! local variable
!
   if (list%typ.ne.list_i8) then
     print *, 'list_change_i8: check list type... '
     stop
   endif
   if (il.gt.list%n) then
     print *, 'list_change_i8: check il... '
     stop
   endif
!
   list%i8(il) = val
!
   end subroutine list_change_i8
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_change_r4(list, il, val)
   implicit none
   type(list_t), intent(inout) :: list
   integer(4),   intent(in   ) :: il
   real(4),      intent(in   ) :: val
! local variable
!
   if (list%typ.ne.list_r4) then
     print *, 'list_change_r4: check list type... '
     stop
   endif
   if (il.gt.list%n) then
     print *, 'list_change_r4: check il... '
     stop
   endif
!
   list%r4(il) = val
!
   end subroutine list_change_r4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_change_r8(list, il, val)
   implicit none
   type(list_t), intent(inout) :: list
   integer(4),   intent(in   ) :: il
   real(8),      intent(in   ) :: val
! local variable
!
   if (list%typ.ne.list_r8) then
     print *, 'list_change_r8: check list type... '
     stop
   endif
   if (il.gt.list%n) then
     print *, 'list_change_r8: check il... '
     stop
   endif
!
   list%r8(il) = val
!
   end subroutine list_change_r8
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_change_st(list, il, val)
   implicit none
   type(list_t), intent(inout) :: list
   integer(4),   intent(in   ) :: il
   character(len=maxlen_str), intent(in   ) :: val
! local variable
!
   if (list%typ.ne.list_st) then
     print *, 'list_change_st: check list type... '
     stop
   endif
   if (il.gt.list%n) then
     print *, 'list_change_st: check il... '
     stop
   endif
!
   list%st(il) = val
!
   end subroutine list_change_st
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   function check_type(typ) result(res)
   implicit none
   integer, intent(in   ) :: typ
   logical                :: res
! local variables
!
   if (typ.ne.list_l4.and.typ.ne.list_i4.and.typ.ne.list_r4.and.typ.ne.list_r8.and.typ.ne.list_st) then
     res = .true.
   else
     res = .false.
   endif
!
   return
   end function check_type
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine list_die(message)
   implicit none
   character(len=*), intent(in   ) :: message
! local variables
!
   write(*,*) trim(message)
   stop
!
   return
   end subroutine list_die
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   end module list
