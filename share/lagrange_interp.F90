!-------------------------------------------------------------------------------
   module lagrange_interp
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
   implicit none
!
   integer, parameter :: i4 = 4, r4 = 4, r8 = 8
!
   type :: lagrange_t
     ! nn: # of points
     integer(i4) :: nn, ndx
     ! x, y, dx, dy: source
     real(r8), dimension(:),   allocatable :: sx, sy
     real(r8), dimension(:),   allocatable :: sdx, sdy
     ! 1/(xj-xi)
     real(r8), dimension(:,:), allocatable :: coef
   end type
!
   public :: lagrange_t
   public :: lagrange_initialize, lagrange_finalize, lagrange_set_x, lagrange_set
   public :: lagrangian, lagrange_difference, lagrange_get, lagrange_integral
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine lagrange_initialize(lp, n)
!-------------------------------------------------------------------------------
   implicit none
!
   type(lagrange_t), intent(inout) :: lp
   integer(i4),      intent(in   ) :: n
! local variables
   integer(i4) :: i,j
!
   lp%nn = n
   lp%ndx = 1000
   allocate(lp%sx(n))
   allocate(lp%sy(n))
   allocate(lp%sdx(n-1))
   allocate(lp%sdy(n-1))
   allocate(lp%coef(n,n))
!
   return
   end subroutine lagrange_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine lagrange_finalize(lp)
!-------------------------------------------------------------------------------
   implicit none
!
   type(lagrange_t), intent(inout) :: lp
!
   deallocate(lp%sx)
   deallocate(lp%sy)
   deallocate(lp%sdx)
   deallocate(lp%sdy)
   deallocate(lp%coef)
!
   return
   end subroutine lagrange_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine lagrange_set_x(lp, xin)
!-------------------------------------------------------------------------------
   implicit none
!
   type(lagrange_t),       intent(inout) :: lp
   real(r8), dimension(:), intent(in   ) :: xin
! local variables
   integer(i4) :: i, j
!
   do i = 1,lp%nn
     lp%sx(i) = xin(i)
   enddo
!
   do i = 1,lp%nn-1
     lp%sdx(i) =(lp%sx(i+1)-lp%sx(i))/dble(lp%nn)
   enddo
!
   call setlagcoef(lp)
!
   return
   end subroutine lagrange_set_x
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine lagrange_set(lp, yin)
!-------------------------------------------------------------------------------
   implicit none
!
   type(lagrange_t),       intent(inout) :: lp
   real(r8), dimension(:), intent(in   ) :: yin
! local variables
   integer(i4) :: i, j
!
   do i = 1,lp%nn
     lp%sy(i) = yin(i)
   enddo
!
   do i = 1,lp%nn-1
     lp%sdy(i) =(lp%sy(i+1)-lp%sy(i))/dble(lp%nn)
   enddo
!
   return
   end subroutine lagrange_set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
! set c_ij
   subroutine setlagcoef(lp)
!-------------------------------------------------------------------------------
   implicit none
!
   type(lagrange_t), intent(inout) :: lp
! local variables
   integer(i4) :: i, j
!
   do j = 1,lp%nn
   do i = 1,lp%nn
     if (i.ne.j) then
       lp%coef(i,j) = 1.0_r8/(lp%sx(j)-lp%sx(i))
     else
       lp%coef(i,j) =-9999.0_r8
     endif
   enddo
   enddo
!
   return
   end subroutine setlagcoef
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
! return L_m(x)
   function lagrangian(lp, m, x)
!-------------------------------------------------------------------------------
   implicit none
!
   type(lagrange_t), intent(inout) :: lp
   integer(i4),      intent(in   ) :: m
   real(r8),         intent(in   ) :: x
   real(r8) :: lagrangian
! local variables
   integer(i4) :: i
!
   lagrangian = 1.0_r8
   do i = 1,lp%nn
     if (i.ne.m) lagrangian = lagrangian*lp%coef(i,m)*(x-lp%sx(i))
   enddo
!
   end function lagrangian
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
! return dL_m(x)/dx
   function lagrange_difference(lp, m, x)
!-------------------------------------------------------------------------------
   implicit none
!
   type(lagrange_t), intent(inout) :: lp
   integer(i4),      intent(in   ) :: m
   real(r8),         intent(in   ) :: x
   real(r8) :: lagrange_difference
! local variables
   integer(i4) :: i, j
   real(r8) :: tmp
!
   lagrange_difference = 0.0_r8
!
   do i = 1,lp%nn
     tmp = 1.0_r8
     do j = 1,lp%nn
       if (j.ne.i.and.j.ne.m) then
         tmp = tmp*lp%coef(j,m)*(x-lp%sx(j))
       endif
     enddo
     if (i.ne.m) lagrange_difference = lagrange_difference+lp%coef(i,m)*tmp
   enddo
!
   end function lagrange_difference
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine lagrange_get(lp, m, x, y, order)
!-------------------------------------------------------------------------------
   implicit none
!
   type(lagrange_t),       intent(inout) :: lp
   integer(i4),            intent(in   ) :: m
   real(r8), dimension(m), intent(in   ) :: x
   real(r8), dimension(m), intent(  out) :: y
   integer(i4), optional,  intent(in   ) :: order
! local variables
   integer(i4) :: k, i, l_order
!
   if (present(order)) then
     l_order = order
   else
     l_order = 0
   endif
   if (l_order>2) l_order = 0
!
   ! check boundary
   if (x(1)<lp%sx(1).or.x(m)>lp%sx(lp%nn)) then
     print*,'Check boundary.... in lagrange_get'
     stop
   endif
!
   if (l_order.eq.0) then
!
     y = 0.0_r8
     do k = 1,m
       y(k) = 0.0_r8
       do i = 1,lp%nn
         y(k) = y(k)+lp%sy(i)*lagrangian(lp,i,x(k))
       enddo
     enddo
!
   elseif (l_order.eq.1) then
!
     y = 0.0_r8
     do k = 1,m
       y(k) = 0.0_r8
       do i = 1,lp%nn
         y(k) = y(k)+lp%sy(i)*lagrange_difference(lp,i,x(k))
       enddo
     enddo
!
   elseif (l_order.eq.2) then
!
     y =-999.0_r8
!
   endif
!
   return
   end subroutine lagrange_get
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
! return \int_0^x_k {L_m(x)} dx  (x_k : full or half level index)
   function lagrange_integral(lp, n, m, k)
!-------------------------------------------------------------------------------
   implicit none
!
   type(lagrange_t), intent(inout) :: lp
   integer(i4),      intent(in   ) :: n, m, k
   real(r8) :: lagrange_integral
! local variables
   integer(i4) :: i, j, np
   real(r8) :: tmp, x
!
   lagrange_integral = 0.0_r8
!
   do j = 1,k
     if (j.ne.k) then
       do i = 1,lp%ndx
         x = lp%sx(j)+lp%sdx(j)*(dble(i)-0.5)
         lagrange_integral = lagrange_integral+lagrangian(lp,m,x)*lp%sdx(j)
       enddo
     else
       do i = 1,lp%ndx/2
         x = lp%sx(j)+lp%sdx(j)*(dble(i)-0.5)
         lagrange_integral = lagrange_integral+lagrangian(lp,m,x)*lp%sdx(j)
       enddo
     endif
   enddo
!
   end function lagrange_integral
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module lagrange_interp
