!-------------------------------------------------------------------------------
   module spline_interp
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    2015-07-27  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer, parameter :: i4 = 4, l4 = 4, r4 = 4, r8 = 8
!
   type :: spline_t
     ! order
     integer(i4) :: order
     ! size
     integer(i4) :: n
     ! array order
     integer(i4) :: array_order
     ! linear    spline
     real(r8) :: y1_l = 0.0_r8, yn_l = 0.0_r8
     ! quadratic spline
     real(r8) :: yp1_q = 0.0_r8
     ! cubic     spline
     real(r8) :: ypp1_c = 0.0_r8, yppn_c = 0.0_r8
     real(r8), dimension(:), allocatable :: b_c
     real(r8), dimension(:,:), allocatable :: a_c, inva_c
     ! values, derivative
     real(r8), dimension(:), allocatable :: x, y, dx, dy
     real(r8), dimension(:), allocatable :: yp, ypp
     ! integral values
     real(r8), dimension(:), allocatable :: inty
   end type
!
   public :: spline_initialize, spline_finalize, spline_set_x, spline_set, spline_get
   public ::  spline_integral, spline_interface
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_initialize(sp, n, order)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),        intent(inout) :: sp
   integer(i4),           intent(in   ) :: n
   integer(i4), optional, intent(in   ) :: order
! local variables
   integer(i4) :: i, k
!
   sp%n = n
!
   if (present(order)) then
     sp%order = order
   else
     sp%order = 3 ! cubic
   endif
!
   if (order.gt.3) then
     print*,'not support order ',order
     sp%order = 3
   endif
!
   allocate(sp%x(n))
   allocate(sp%y(n))
   allocate(sp%dx(n-1))
   allocate(sp%dy(n-1))
!
   if (order.eq.2) then
     allocate(sp%yp(n))
   elseif (order.eq.3) then
     allocate(sp%yp(n))
     allocate(sp%ypp(n))
     allocate(sp%a_c(2:n-1,2:n-1))
     allocate(sp%inva_c(2:n-1,2:n-1))
     allocate(sp%b_c(2:n-1))
   endif
   allocate(sp%inty(n-1))
!
   return
   end subroutine spline_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_finalize(sp)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t), intent(inout) :: sp
!
   deallocate(sp%x)
   deallocate(sp%y)
   deallocate(sp%dx)
   deallocate(sp%dy)
!
   if (sp%order.eq.2) then
     deallocate(sp%yp)
   elseif (sp%order.eq.3) then
     deallocate(sp%yp)
     deallocate(sp%ypp)
     deallocate(sp%a_c)
     deallocate(sp%inva_c)
     deallocate(sp%b_c)
   endif
   deallocate(sp%inty)
!
   return
   end subroutine spline_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_set_x(sp, x_s)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),         intent(inout) :: sp
   real(r8), dimension(:), intent(in   ) :: x_s
! local variables
   integer(i4) :: i, j, k, l_ierr, n
   real(r8) :: tmp
   real(r8), dimension(:), allocatable :: ipivot, work
!
   n = sp%n
!
! check array order
   sp%array_order = check_array_order(x_s)
   if (sp%array_order.eq.0) then
     print *,'error: check array order (source)'
     stop
   endif
!
   allocate(ipivot(n-2))
   allocate(work(n-2))
!
   ! set x
   sp%x(:) = x_s(:)
   do k = 1,n-1
     sp%dx(k) = sp%x(k+1)-sp%x(k)
   enddo
!
   if (sp%order.eq.1) then


   elseif (sp%order.eq.2) then


   elseif (sp%order.eq.3) then
!  
     ! set A matrix
     do j = 2,n-1
     do i = 2,n-1
       if (i.eq.j) then
         sp%a_c(i,j) =(sp%x(i+1)-sp%x(i-1))/3.0_r8
       elseif (i.eq.j-1) then
         sp%a_c(i,j) = sp%dx(i)/6.0_r8
       elseif (i.eq.j+1) then
         sp%a_c(i,j) = sp%dx(j)/6.0_r8
       else
       endif
     enddo
     enddo
     !
     sp%inva_c(:,:) = sp%a_c(:,:)
     !
     call dgetrf(n-2,n-2,sp%inva_c,n-2,ipivot,l_ierr)
     if (l_ierr.ne.0) stop
     call dgetri(n-2,sp%inva_c,n-2,ipivot,work,n-2,l_ierr)
     if (l_ierr.ne.0) stop
!
   endif
!
   deallocate(ipivot)
   deallocate(work)
!
   return
   end subroutine spline_set_x
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_set(sp, y_s, yp1, ypp1, yppn)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),         intent(inout) :: sp
   real(r8), dimension(:), intent(in   ) :: y_s
   real(r8), optional,     intent(in   ) :: yp1, ypp1, yppn
! local variables
   integer(i4) :: i, j, k, n, l_ierr
   real(r8) :: tmp
!
   n = sp%n
!
   ! quadratic
   if (sp%order.eq.2) then
     if (present(yp1)) then
       sp%yp1_q = yp1
     else
       sp%yp1_q = 0.0_r8
     endif
   endif
   ! cubic
   if (sp%order.eq.3) then
     if (present(ypp1)) then
       sp%ypp1_c = ypp1
     else
       sp%ypp1_c = 0.0_r8
     endif
     if (present(yppn)) then
       sp%yppn_c = yppn
     else
       sp%yppn_c = 0.0_r8
     endif
   endif
   ! set y
   sp%y(:) = y_s(:)
   do k = 1,n-1
     sp%dy(k) = sp%y(k+1)-sp%y(k)
   enddo
!
   if (sp%order.eq.1) then


   elseif (sp%order.eq.2) then
!
     ! boundary
     !sp%yp(1) = sp%yp1_q
     sp%yp(1) = sp%dy(1)/sp%dx(1)
     do k = 1,n-1
       sp%yp(k+1) =-sp%yp(k)+2.0_r8*sp%dy(k)/sp%dx(k)
     enddo
!
   elseif (sp%order.eq.3) then
!
     ! 2rd derivative
     do k = 1,n
       sp%ypp(k) = 0.0
     enddo
     ! boundary
     sp%ypp(1) = sp%ypp1_c
     sp%ypp(n) = sp%yppn_c
     ! set b
     do k = 2,n-1
       sp%b_c(k) = sp%dy(k)/sp%dx(k)-sp%dy(k-1)/sp%dx(k-1)
     enddo
     sp%b_c(2) = sp%b_c(2)-sp%dx(2-1)*sp%ypp(2-1)/6.0_r8
     sp%b_c(n-1) = sp%b_c(n-1)-sp%dx(n-1)*sp%ypp(n-1)/6.0_r8
     ! set y''
     do j = 2,n-1
     do i = 2,n-1
       sp%ypp(i) = sp%ypp(i)+sp%inva_c(i,j)*sp%b_c(j)
     enddo
     enddo
!
#if DEBUG
     print*,'source x = '
     print*,sx
     print*,'source y = '
     print*,sy
     print*,'source dx = '
     print*,sdx
     print*,'source dy = '
     print*,sdy
     print*,'source b = '
     print*,b_c
     call printmatrix('source a = ',a_c,n-2,n-2)
     call printmatrix('source inva = ',inva_c,n-2,n-2)
     print*,'identity = '
     do j = 2,n-1
     do i = 2,n-1
     tmp = 0.0_r8
     do k = 2,n-1
     tmp = tmp+a_c(i,k)*inva_c(k,j)
     enddo
     print*,tmp
     enddo
     enddo
     print*,'source sypp = '
     print*,sypp
#endif
!
   endif
!
   ! integral
   if (sp%order.eq.1) then
!
     do k = 1,n-1
       sp%inty(k) = 0.5_r8*sp%dx(k)*(sp%y(k)+sp%y(k+1))
     enddo
!
   elseif (sp%order.eq.2) then
!
     do k = 1,n-1
       sp%inty(k) =(sp%yp(k+1)+2.0_r8*sp%yp(k))/6.0_r8*sp%dx(k)**2.0_r8+sp%y(k)*sp%dx(k)
     enddo
!
   elseif (sp%order.eq.3) then
!
     do k = 1,n-1
       sp%inty(k) = 0.5_r8*sp%dx(k)*(sp%y(k+1)+sp%y(k))-sp%dx(k)**3.0_r8/48.0_r8*(sp%ypp(k)+sp%ypp(k+1))
     enddo
!
   endif
!
   return
   end subroutine spline_set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function spline_value_in_piece(sp, kp, x, order) result(value)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),        intent(inout) :: sp
   integer(i4),           intent(in   ) :: kp
   real(r8),              intent(in   ) :: x
   integer(i4), optional, intent(in   ) :: order ! 0:function, 1:first derivative, 2:second derivative
   real(r8)                             :: value
! local variables
   integer(i4) :: n, l_order
   real(r8) :: a, b
!
   n = sp%n
!
   if (present(order)) then
     l_order = order
   else
     l_order = 0
   endif
   if (l_order>2) l_order = 0
!
     ! integration
     if (l_order.eq.-1) then
  
     if (sp%order.eq.1) then
     !
     a =(sp%x(kp+1)-x)/sp%dx(kp)
     b =(x-sp%x(kp))/sp%dx(kp)
     value =-0.5_r8*sp%dx(kp)*((a*a-1.0_r8)*sp%y(kp)-b*b*sp%y(kp+1))
     !
     elseif (sp%order.eq.2) then
     !
     a =(sp%x(kp+1)-x)/sp%dx(kp)
     b =(x-sp%x(kp))/sp%dx(kp)
     value = sp%dx(kp)*b*((sp%yp(kp+1)-sp%yp(kp))*sp%dx(kp)*b*b+0.5_r8*sp%yp(kp)*sp%dx(kp)*b+sp%y(kp))
     !
     elseif (sp%order.eq.3) then
     !
     a =(sp%x(kp+1)-x)/sp%dx(kp)
     b =(x-sp%x(kp))/sp%dx(kp)
     value =-0.5_r8*sp%dx(kp)*(sp%y(kp)*a*a-sp%y(kp+1)*b*b) -sp%dx(kp)**3.0_r8/12.0_r8*((a**4.0_r8/4.0_r8-a**2.0_r8/2.0_r8)*sp%ypp(kp)-(b**4.0_r8/4.0_r8-b**2.0_r8/2.0_r8)*sp%ypp(kp+1)) +0.5_r8*sp%dx(kp)*sp%y(kp)-sp%dx(kp)**3.0_r8/48.0_r8*(sp%ypp(kp)+sp%ypp(kp+1))
     !
     endif
!
   ! 0th value: interpolation value
   elseif (l_order.eq.0) then
!
     if (sp%order.eq.1) then
       value =(sp%x(kp+1)-x)/sp%dx(kp)*sp%y(kp)+(x-sp%x(kp))/sp%dx(kp)*sp%y(kp+1)
     elseif (sp%order.eq.2) then
       a =(sp%x(kp+1)-x)/sp%dx(kp)
       b =(x-sp%x(kp))/sp%dx(kp)
       value = 0.5_r8*(sp%yp(kp+1)-sp%yp(kp))*sp%dx(kp)*b**2.0_r8+sp%yp(kp)*sp%dx(kp)*b+sp%y(kp)
     elseif (sp%order.eq.3) then
       a =(sp%x(kp+1)-x)/sp%dx(kp)
       b =(x-sp%x(kp))/sp%dx(kp)
       value =(a*sp%y(kp)+b*sp%y(kp+1))+((a**3.0_r8-a)*sp%ypp(kp)+(b**3.0_r8-b)*sp%ypp(kp+1))*sp%dx(kp)*sp%dx(kp)/6.0_r8
     endif
!
   ! 1st derivative value
   elseif (l_order.eq.1) then
!
     if (sp%order.eq.1) then
       value = sp%dy(kp)/sp%dx(kp)
     elseif (sp%order.eq.2) then
       a =(sp%x(kp+1)-x)/sp%dx(kp)
       b =(x-sp%x(kp))/sp%dx(kp)
       value =(sp%yp(kp+1)-sp%yp(kp))*b+sp%yp(kp)
     elseif (sp%order.eq.3) then
       a =(sp%x(kp+1)-x)/sp%dx(kp)
       b =(x-sp%x(kp))/sp%dx(kp)
       value = sp%dy(kp)/sp%dx(kp)-sp%dx(kp)/6.0_r8*((3.0_r8*(a**2.0_r8)-1.0_r8)*sp%ypp(kp)-(3.0_r8*(b**2.0_r8)-1.0_r8)*sp%ypp(kp+1))
     endif
!
   ! 2nd derivative value
   elseif (l_order.eq.2) then
!
     if (sp%order.eq.1) then
       value = 0.0_r8
     elseif (sp%order.eq.2) then
       value =(sp%yp(kp+1)-sp%yp(kp))/sp%dx(kp)
     elseif (sp%order.eq.3) then
       a =(sp%x(kp+1)-x)/sp%dx(kp)
       b =(x-sp%x(kp))/sp%dx(kp)
       value = a*sp%ypp(kp)+b*sp%ypp(kp+1)
     endif
!
   endif
!
   end function spline_value_in_piece
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_get(sp, m, x, y, order)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),         intent(inout) :: sp
   integer(i4),            intent(in   ) :: m
   real(r8), dimension(m), intent(in   ) :: x
   real(r8), dimension(m), intent(  out) :: y
   ! order = 0:function, 1:first derivative, 2:second derivative
   integer(i4), optional,  intent(in   ) :: order
! local variables
   integer(i4) :: kp, k, i, n
   integer(i4) :: l_order, array_order
   logical(l4) :: in_src
!
   n = sp%n
!
   if (present(order)) then
     l_order = order
   else
     l_order = 0
   endif
   if (l_order>2) l_order = 0
!
! check array order
   array_order = check_array_order(x)
   if (sp%array_order.ne.array_order) then
     print *,'error: check array order (target)',sp%array_order,array_order
     stop
   endif
! check boundary
   in_src = check_array_boundary(sp%x,x,array_order)
   if (.not.in_src) then
     print*,'error: check boundary.... in spline_get',sp%x(1),x(1),sp%x(n),x(m)
     stop
   endif
!
   do i = 1,m
     kp   = get_idx_in_array(n,sp%x,x(i),sp%array_order)
     y(i) = spline_value_in_piece(sp,kp,x(i),l_order)
   enddo
!
   return
   end subroutine spline_get
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_integral(sp, m, x, y, order)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),         intent(inout) :: sp
   integer(i4),            intent(in   ) :: m
   real(r8), dimension(m), intent(in   ) :: x
   real(r8), dimension(m), intent(  out) :: y
   ! order = 0:function, 1:first derivative, 2:second derivative
   integer(i4), optional,  intent(in   ) :: order
! local variables
   integer(i4) :: kp(m)
   integer(i4) :: nk, k, i, n
   integer(i4) :: l_order, array_order
   logical(l4) :: in_src
   real(r8) :: areak, areakp1, areai
!
   n = sp%n
   l_order =-1
!
! check array order
   array_order = check_array_order(x)
   if (sp%array_order.ne.array_order) then
     print *,'error: check array order (target)',sp%array_order,array_order
     stop
   endif
! check boundary
   in_src = check_array_boundary(sp%x,x,array_order)
   if (.not.in_src) then
     print*,'error: check boundary.... in piecewise_get_fun_integral'
     stop
   endif
!
   do i = 1,m
     kp(i) = get_idx_in_array(n,sp%x,x(i),sp%array_order)
   enddo
!
   y(1) = 0.0_r8
   do i = 1,m-1
!
     if (kp(i).eq.kp(i+1)) then
       areak   = spline_value_in_piece(sp,kp(i),x(i),l_order)
       areakp1 = spline_value_in_piece(sp,kp(i+1),x(i+1),l_order)
       y(i+1)  = areakp1-areak
     else
       areai = 0.0_r8
       do k = kp(i),kp(i+1)-1
         areai = areai+sp%inty(k)
       enddo
       areak   = spline_value_in_piece(sp,kp(i),x(i),l_order)
       areakp1 = spline_value_in_piece(sp,kp(i+1),x(i+1),l_order)
       y(i+1)  = areai-areak+areakp1
     endif
!
   enddo
!
   return
   end subroutine spline_integral
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_interface(sp, n, sx, sy, m, tx, ty, method, yp, ypp1, yppn)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),         intent(inout) :: sp
   integer(i4),            intent(in   ) :: n
   real(r8), dimension(n), intent(in   ) :: sx
   real(r8), dimension(n), intent(in   ) :: sy
   integer(i4),            intent(in   ) :: m
   real(r8), dimension(m), intent(in   ) :: tx
   real(r8), dimension(m), intent(  out) :: ty
   ! method = 1:linear, 2:quadratic, 3: cubic
   integer(i4),            intent(in   ) :: method
   real(r8),               intent(in   ) :: yp, ypp1, yppn
! local variables
   integer(i4) :: i
   real(r8)    :: vmin
!
   call spline_initialize(sp,n,method)
   call spline_set_x(sp,sx)
   call spline_set(sp,sy,yp,ypp1,yppn)
   call spline_get(sp,n,tx,ty,0)
   call spline_finalize(sp)
!
   vmin = 1.0d-100
   do i=1,m
     if (abs(ty(i)).le.vmin) ty(i) = 0.0_r8
   enddo
!
   return
   end subroutine spline_interface
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function check_array_order(array) result(order) ! order=0(mix), order=1(ascend), order=-1(descend)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:), intent(in   ) :: array
   integer(i4)                           :: order
! local variables
   integer(i4) :: n, i
!
   n = size(array)
   order = 0
!
   do i = 2, n
     if (i.eq.2) then
       if (array(i).ge.array(i-1)) then
         order = 1
       else
         order = -1
       endif
     else
       if ((array(i).gt.array(i-1).and.order.eq.-1).or.(array(i).lt.array(i-1).and.order.eq.1)) then
         order = 0
         exit
       endif
     endif
   enddo
!
   end function check_array_order
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function check_array_boundary(src,tgt,order) result(in_src)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:), intent(in   ) :: src, tgt
   integer(i4),            intent(in   ) :: order
   logical(l4)                           :: in_src
! local variables
   integer(i4) :: n, m
!
   n = size(src)
   m = size(tgt)
   in_src = .true.
   if (order.eq.1) then
     if (tgt(1).lt.src(1).or.tgt(m).gt.src(n)) then
       in_src = .false.
     endif
   elseif (order.eq.-1) then
     if (tgt(1).gt.src(1).or.tgt(m).lt.src(n)) then
       in_src = .false.
     endif
   else
     print *,'error: check_array_boundary'
     stop
   endif
!
   end function check_array_boundary
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_idx_in_array(n, array, val, array_order) result(idx)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),            intent(in   ) :: n
   real(r8), dimension(n), intent(in   ) :: array
   real(r8),               intent(in   ) :: val
   integer(i4)                           :: array_order
   integer(i4)                           :: idx
! local variables
   integer(i4) :: i
   logical(l4) :: ascending_order
!
   idx = -1
   if (val.eq.array(n)) then
       idx = n-1
   else
     if (array_order.eq.1) then
       do i = 2,n
         if (val.ge.array(i-1).and.val.lt.array(i)) then
           idx = i-1
           exit
         endif
       enddo
     else    ! if (array_order.eq.-1)
       do i = 2,n
           if (val.le.array(i-1).and.val.gt.array(i)) then
           idx = i-1
           exit
         endif
       enddo
     endif
   endif
   if (idx.eq.-1) then
     print *,'error: check val range...'
   endif
!
   end function get_idx_in_array
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module spline_interp
!-------------------------------------------------------------------------------
