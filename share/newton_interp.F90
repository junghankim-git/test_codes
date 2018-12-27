!-------------------------------------------------------------------------------
   module newton_interp
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    2017-12-11  junghan kim    initial setup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer, parameter :: i4 = 4, r4 = 4, r8 = 8
!
   type :: newton_t
     ! n: # of points
     integer(i4) :: n
     ! x, y : source
     real(r8), dimension(:), allocatable :: x, y
     ! a0, a1, ..., an
     real(r8), dimension(:), allocatable :: a
   end type
!
   public :: newton_t
   public :: newton_initialize, newton_finalize, newton_set_x, newton_set
   public :: newton_get
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine newton_initialize(np, n)
!-------------------------------------------------------------------------------
   implicit none
!
   type(newton_t), intent(inout) :: np
   integer(i4),    intent(in   ) :: n
! local variables
   integer(i4) :: i,j
!
   np%n = n
   allocate(np%x(n))
   allocate(np%y(n))
   allocate(np%a(n))
!
   return
   end subroutine newton_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine newton_finalize(np)
!-------------------------------------------------------------------------------
   implicit none
!
   type(newton_t), intent(inout) :: np
!
   deallocate(np%x)
   deallocate(np%y)
   deallocate(np%a)
!
   return
   end subroutine newton_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine newton_set_x(np, x)
!-------------------------------------------------------------------------------
   implicit none
!
   type(newton_t),         intent(inout) :: np
   real(r8), dimension(:), intent(in   ) :: x
! local variables
   integer(i4) :: i, j
!
   do i = 1,np%n
     np%x(i) = x(i)
   enddo
!
   return
   end subroutine newton_set_x
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine newton_set(np, y)
!-------------------------------------------------------------------------------
   implicit none
!
   type(newton_t),         intent(inout) :: np
   real(r8), dimension(:), intent(in   ) :: y
! local variables
   integer(i4) :: i, j
!
   do i = 1,np%n
     np%y(i) = y(i)
   enddo
!
   do i = 1,np%n
     np%a(i) = divided_difference(np%n,np%x,np%y,1,i)
   enddo
!
   return
   end subroutine newton_set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   recursive function divided_difference(n, x, y, is, ie) result(res)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),            intent(in   ) :: n
   real(r8), dimension(n), intent(in   ) :: x, y
   integer(i4),            intent(in   ) :: is, ie
   real(r8)                              :: res
! local variables
   integer(i4) :: i, j
!
   if ((is.lt.1.or.is.gt.n).or.(ie.lt.1.or.ie.gt.n)) then
     print *,'error: check is or ie in divided_difference...'
     stop
   endif
!
   if (is.eq.ie) then
     res = y(is)
   else
     res = (divided_difference(n,x,y,is+1,ie)-divided_difference(n,x,y,is,ie-1))/(x(ie)-x(is))
   endif
!
   end function divided_difference
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine newton_get(np, m, x, y, order)
!-------------------------------------------------------------------------------
   implicit none
!
   type(newton_t),         intent(inout) :: np
   integer(i4),            intent(in   ) :: m
   real(r8), dimension(m), intent(in   ) :: x
   real(r8), dimension(m), intent(  out) :: y
   integer(i4), optional,  intent(in   ) :: order
! local variables
   integer(i4) :: i, j, k, l_order
   real(r8)    :: tmp
!
   if (present(order)) then
     l_order = order
   else
     l_order = 0
   endif
   if (l_order>2) l_order = 0
!
   ! check boundary
   if (x(1)<np%x(1).or.x(m)>np%x(np%n)) then
     print*,'check boundary.... in newton_get'
     stop
   endif
!
   if (l_order.eq.0) then
!
     do i = 1,m
       y(i) = 0.0_r8
       do j = 1,np%n
         tmp = np%a(j)
         do k = 1,j-1
           tmp = tmp*(x(i)-np%x(k))
         enddo
         y(i) = y(i)+tmp
       enddo
     enddo
!
   elseif (l_order.eq.1) then
!
!
   elseif (l_order.eq.2) then
!
!
   endif
!
   return
   end subroutine newton_get
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module newton_interp
