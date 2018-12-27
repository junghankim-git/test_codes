!#define USE_LAGRANGIAN
#undef USE_LAGRANGIAN
!-------------------------------------------------------------------------------
   module piecewise_remap
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    2015-09-12  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure:
!
!-------------------------------------------------------------------------------
#ifdef USE_LAGRANGIAN
   use lagrange_interp, only: lagrange_t, lagrange_initialize, lagrange_finalize
   use lagrange_interp, only: lagrange_set_x, lagrange_set, lagrange_get, lagrange_integral
#endif
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   integer, parameter :: i4 = 4, l4 = 4, r4 = 4, r8 = 8
!
   type, public :: piecewise_t
     ! order (input)
     integer(i4) :: order
     ! monotonic (input)
     logical(l4) :: mono
     ! size
     integer(i4) :: n
     ! bndry condition (input)
     integer(i4) :: bndry, ha, lha, rha ! 0:mirror, 1:mirror+derivative, 2:periodic(need base)
     ! array order
     integer(i4) :: array_order
     ! values, derivative
     real(r8), dimension(:),   allocatable :: x, y, cx, dx, dy
     real(r8), dimension(:),   allocatable :: cy, cinty
     real(r8), dimension(:,:), allocatable :: mcy
     real(r8), dimension(:),   allocatable :: dely
     real(r8), dimension(:),   allocatable :: inty
   end type
!
   public :: piecewise_initialize, piecewise_finalize
   public :: piecewise_set_x, piecewise_set
   public :: piecewise_get_fun, piecewise_get_fun_integral
   public :: piecewise_get, piecewise_interface, piecewise_cyclic_interface
   public :: piecewise_print_status
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine piecewise_initialize(pw, n, order, mono, bndry, lha, rha)
!-------------------------------------------------------------------------------
   implicit none
!
   type(piecewise_t),     intent(inout) :: pw
   integer(i4),           intent(in   ) :: n
   integer(i4), optional, intent(in   ) :: order
   logical(i4), optional, intent(in   ) :: mono
   integer(i4), optional, intent(in   ) :: bndry
   integer(i4), optional, intent(in   ) :: lha, rha
! local variables
   integer(i4) :: i, ha
   logical :: mono_l
!
   pw%n = n
!
   if (present(mono)) then
     pw%mono = mono
   else
     pw%mono = .false. ! cubic
   endif
!
   if (present(order)) then
     pw%order = order
   else
     pw%order = 2 ! cubic
   endif
!
   if (present(bndry)) then
     pw%bndry = bndry
   else
     pw%bndry = 0
   endif
!
   if (pw%bndry==2) then
     if (present(lha)) then
       pw%lha = lha
     else
       pw%lha = 0
     endif
     if (present(rha)) then
       pw%rha = rha
     else
       pw%rha = 0
     endif
     if (pw%lha>0.and.pw%rha>0) then
       print *,'error: check lha or rha ',pw%lha,pw%rha
       stop
     endif
   endif
!
   if (order.gt.3) then
     print *,'error: not support order ',order
     stop
   endif
!
   pw%ha = 2
   ha = pw%ha
!
   allocate(pw%x(1-ha:n+ha))       ! x at grid point
   allocate(pw%y(1-ha:n+ha))       ! y at grid point
   allocate(pw%cx(1-ha-1:n+ha))    ! x at cell point
   allocate(pw%dx(1-ha:n+ha))      ! dx between cells
   allocate(pw%cinty(1-ha-1:n+ha)) ! integral value at cell point:dx_i*y_i
!
   allocate(pw%cy(1-1:n))
   allocate(pw%mcy(2,0:n))
   allocate(pw%dy(1-ha:n+ha))      ! dx at grid point
   allocate(pw%dely(1-ha:n+ha))    ! 
   allocate(pw%inty(0:n))
!
   return
   end subroutine piecewise_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine piecewise_finalize(pw)
!-------------------------------------------------------------------------------
   implicit none
!
   type(piecewise_t), intent(inout) :: pw
!
   deallocate(pw%x)
   deallocate(pw%y)
   deallocate(pw%cx)
   deallocate(pw%dx)
   deallocate(pw%cinty)
!
   deallocate(pw%cy)
   deallocate(pw%mcy)
   deallocate(pw%dy)
   deallocate(pw%dely)
   deallocate(pw%inty)
!
   return
   end subroutine piecewise_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine piecewise_set_x(pw, cx_s)
!-------------------------------------------------------------------------------
   implicit none
!
   type(piecewise_t),      intent(inout) :: pw
   real(r8), dimension(:), intent(in   ) :: cx_s
! local variables
   integer(i4) :: i, n
!
! check array order
   pw%array_order = check_array_order(cx_s)
   if (pw%array_order.eq.0) then
     print *,'error: check array order (source)'
     stop
   endif
!
   n = pw%n
! set cell x
   pw%cx(0:n) = cx_s(:)
   ! halo
   if (pw%bndry==-1) then
   elseif (pw%bndry==0.or.pw%bndry==1) then ! mirror
     do i = 1,pw%ha
       pw%cx(0-i) = pw%cx(0)-(pw%cx(i)-pw%cx(0))
       pw%cx(n+i) = pw%cx(n)+(pw%cx(n)-pw%cx(n-i))
     enddo
   elseif (pw%bndry==2) then ! periodic
     do i = 1,pw%ha
       !pw%cx(0-i) = pw%cx(0)-(pw%cx(n-pw%base)-pw%cx(n-pw%base-i))
       !pw%cx(n+i) = pw%cx(n)+(pw%cx(pw%base+i)-pw%cx(pw%base+0))
       pw%cx(0-i) = pw%cx(0)-(pw%cx(n-1-pw%lha)-pw%cx(n-1-pw%lha-i))
       pw%cx(n+i) = pw%cx(n)+(pw%cx(pw%rha+1+i)-pw%cx(pw%rha+1+0))
     enddo
   else
     print *,'error: not support bndry condition...'
     stop
   endif
!
! set grid x
! grid x
   do i = 1-pw%ha,n+pw%ha
     pw%x(i) = 0.5_r8*(pw%cx(i-1)+pw%cx(i))
   enddo
!
! set dx
! dx
   do i = 1-pw%ha,n+pw%ha
     pw%dx(i) = pw%cx(i)-pw%cx(i-1)
   enddo
!
   return
   end subroutine piecewise_set_x
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine piecewise_set(pw, y_s)
!-------------------------------------------------------------------------------
   implicit none
!
   type(piecewise_t),      intent(inout) :: pw
   real(r8), dimension(:), intent(in   ) :: y_s
! local variables
   integer(i4) :: i, n, ha
#ifdef USE_LAGRANGIAN
   type(lagrange_t) :: lag
#else
   real(r8) :: tmp
   real(r8), allocatable, dimension(:,:) :: arg
#endif
!
   n  = pw%n
   ha = pw%ha
!
! set grid y and halo
   pw%y(1:n) = y_s(:)
   ! halo
   if (pw%bndry==-1) then
   elseif (pw%bndry==0) then ! mirror
     do i = 1,ha
       pw%y(1-i) = pw%y(i)
       pw%y(n+i) = pw%y(n-i+1)
     enddo
   elseif (pw%bndry==1) then ! mirror+derivative
     do i = 1,ha
       pw%y(1-i) = pw%y(1)-(pw%y(1+i)-pw%y(1))
       pw%y(n+i) = pw%y(n)+(pw%y(n)-pw%y(n-i))
     enddo
   elseif (pw%bndry==2) then ! periodic
     if (pw%lha>0) then
       do i = 1,ha
         pw%y(1-i) = pw%y(n+1-pw%lha-i)
         pw%y(n+i) = pw%y(pw%lha+i)
       enddo
     else
       do i = 1,ha
         pw%y(1-i) = pw%y(n+1-pw%rha-i)
         pw%y(n+i) = pw%y(pw%rha+i)
       enddo
     endif
   else
     print *,'error: not support bndry condition...'
     stop
   endif
!
! calculate the y(i-1/2),y(i+1/2)
   if (pw%order.eq.0) then
!
   ! don't need
!
   else
!
     ! arguments for cell values
#ifdef USE_LAGRANGIAN
     pw%cinty(0-ha) = 0.0_r8
     do i = 1-ha,n+ha
     pw%cinty(i) = pw%cinty(i-1)+pw%y(i)*pw%dx(i)
     enddo
     call lagrange_initialize(lag,5)
     do i = 0,n
       call lagrange_set_x(lag,pw%cx(i-ha:i+ha))
       call lagrange_set(lag,pw%cinty(i-ha:i+ha))
       call lagrange_get(lag,1,pw%cx(i),pw%cy(i),1) ! y(i+1/2) = d lag/d x
     enddo
     call lagrange_finalize(lag)
#else
     allocate(arg(10,0-(ha):n+(ha)))
     do i = 0,n
       arg(1,i) = pw%dx(i)/(pw%dx(i)+pw%dx(i+1))
       arg(2,i) = 1.0_r8/(pw%dx(i-1)+pw%dx(i)+pw%dx(i+1)+pw%dx(i+2))
       arg(3,i) =(2.0_r8*pw%dx(i+1)*pw%dx(i))/(pw%dx(i)+pw%dx(i+1))
       arg(4,i) =(pw%dx(i-1)+pw%dx(i))/(2.0_r8*pw%dx(i)+pw%dx(i+1))
       arg(5,i) =(pw%dx(i+2)+pw%dx(i+1))/(2.0_r8*pw%dx(i+1)+pw%dx(i))
       arg(6,i) = pw%dx(i)*(pw%dx(i-1)+pw%dx(i))/(2.0_r8*pw%dx(i)+pw%dx(i+1))
       arg(7,i) = pw%dx(i+1)*(pw%dx(i+1)+pw%dx(i+2))/(pw%dx(i)+2.0_r8*pw%dx(i+1))
     enddo
     do i = 0,n+1
       arg(8,i) = pw%dx(i)/(pw%dx(i-1)+pw%dx(i)+pw%dx(i+1))
       arg(9,i) =(2.0_r8*pw%dx(i-1)+pw%dx(i))/(pw%dx(i+1)+pw%dx(i))
       arg(10,i) =(pw%dx(i)+2.0_r8*pw%dx(i+1))/(pw%dx(i-1)+pw%dx(i))
     enddo
     ! average slop in the i'th cell
     do i = 0,n+1
       pw%dely(i) = arg(8,i)*(arg(9,i)*(pw%y(i+1)-pw%y(i))+arg(10,i)*(pw%y(i)-pw%y(i-1)))
     enddo
     ! apply monotinic for delta y
     if (pw%mono) then
       do i = 0,n+1
         if (((pw%y(i+1)-pw%y(i))*(pw%y(i)-pw%y(i-1)))>0.0_r8) then
           tmp = min(abs(pw%dely(i)),abs(2.0_r8*(pw%y(i+1)-pw%y(i))),abs(2.0_r8*(pw%y(i)-pw%y(i-1))))
           pw%dely(i) = sign(tmp,pw%dely(i))
         else
           pw%dely(i) = 0.0_r8
         endif
       enddo
     endif
     ! set cell y:  y(i-1/2) y(i+1/2)
     do i = 0,n
       pw%cy(i) = pw%y(i)+arg(1,i)*(pw%y(i+1)-pw%y(i))+arg(2,i)*(arg(3,i)*(arg(4,i)-arg(5,i))*(pw%y(i+1)-pw%y(i))-arg(6,i)*pw%dely(i+1)+arg(7,i)*pw%dely(i))
     enddo
     deallocate(arg)
#endif
     !  
     ! modified cell y
     do i = 1,n
       pw%mcy(1,i) = pw%cy(i-1)
       pw%mcy(2,i) = pw%cy(i)
     enddo
     !
     if (pw%mono) then
       do i = 1,n
         if (((pw%cy(i)-pw%y(i))*(pw%y(i)-pw%cy(i-1))).le.0.0_r8) then
           pw%mcy(1,i) = pw%y(i)
           pw%mcy(2,i) = pw%y(i)
         elseif (((pw%cy(i)-pw%cy(i-1))*(pw%y(i)-0.5_r8*(pw%cy(i-1)+pw%cy(i)))).gt.(pw%cy(i)-pw%cy(i-1))**2.0_r8/6.0_r8) then
           pw%mcy(1,i) = 3.0_r8*pw%y(i)-2.0_r8*pw%mcy(2,i)
         elseif (-(pw%cy(i)-pw%cy(i-1))**2.0_r8/6.0_r8.gt.((pw%cy(i)-pw%cy(i-1))*(pw%y(i)-0.5_r8*(pw%cy(i-1)+pw%cy(i))))) then
           pw%mcy(2,i) = 3.0_r8*pw%y(i)-2.0_r8*pw%mcy(1,i)
         endif
       enddo
     endif
!
   endif ! pw%order.ne.0
!
! set grid dy
   do i = 1,n
     pw%dy(i) = pw%mcy(2,i)-pw%mcy(1,i)
! DEBUG
!print *,pw%mcy(1,i),pw%mcy(2,i)
   enddo
!
! integral
   pw%inty(0) = 0.0
   do i = 1,n
     pw%inty(i) = pw%y(i)*pw%dx(i)
   enddo
!
   return
   end subroutine piecewise_set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine piecewise_print_status(pw)
!-------------------------------------------------------------------------------
   implicit none
!
   type(piecewise_t), intent(inout) :: pw
! local variables
   integer :: k
   character(len=20) :: form_idx, form_n, form_n_p1, form_n_halo, form_n_halo_p1, form_n_halo_m1
   real(r8) :: mass
!
   write(form_idx,'(a,i0.2,a)') '(a,x,',pw%n+2*pw%ha,'(i7.3)) '
   write(form_n,'(a,i0.2,a)') '(a,x,',pw%n,'(f7.3)) '
   write(form_n_p1,'(a,i0.2,a)') '(a,x,',pw%n+1,'(f7.3)) '
   write(form_n_halo_p1,'(a,i0.2,a)') '(a,x,',pw%n+2*pw%ha+1,'(f7.3)) '
   write(form_n_halo_m1,'(a,i0.2,a)') '(a,x,',pw%n+2*pw%ha-1,'(f7.3)) '
   write(form_n_halo,'(a,i0.2,a)') '(a,x,',pw%n+2*pw%ha,'(f7.3)) '
!
   print form_idx,'idx = ',(/-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21/)
   print form_n_halo,' x = ',pw%x(:)
   print form_n_halo_p1,'cx = ',pw%cx(:)
   print form_n_halo,'dx = ',pw%dx(:)
   print form_n_halo,' y = ',pw%y(:)
   print form_n_p1,'cy = ',pw%cy(:)
   print form_n_halo,'dy = ',pw%dy(:)
   print form_n_halo,'dely = ',pw%dely(:)
!
   return
   end subroutine piecewise_print_status
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_spline_value_in_piece(pw, kp, x, vorder) result(value)
!-------------------------------------------------------------------------------
   implicit none
!
   type(piecewise_t),     intent(inout) :: pw
   integer(i4),           intent(in   ) :: kp
   real(r8),              intent(in   ) :: x
   integer(i4), optional, intent(in   ) :: vorder ! 0:function, 1:first derivative, 2:second derivative
   real(r8)                             :: value
! local variables
   integer(i4) :: n
   real(r8) :: a, b, f6
!
   ! integration
   if (vorder.eq.-1) then
   !
     if (pw%order.eq.0) then
       value =(x-pw%cx(kp-1))*pw%y(kp)
     elseif (pw%order.eq.1) then
       b =(x-pw%cx(kp-1))/pw%dx(kp)
       value = pw%dx(kp)*((pw%y(kp)-pw%mcy(1,kp))*b**2.0_r8+pw%mcy(1,kp)*b)
     elseif (pw%order.eq.2) then
       !a  = (pw%x(kp+1)-x)/pw%dx(kp)
       b =(x-pw%cx(kp-1))/pw%dx(kp)
       f6 = 6.0_r8*(pw%y(kp)-0.5_r8*(pw%mcy(1,kp)+pw%mcy(2,kp)))
       value = pw%dx(kp)*(pw%mcy(1,kp)*b+0.5_r8*(pw%dy(kp)+f6)*b**2.0_r8-f6/3.0_r8*b**3.0_r8)
     elseif (pw%order.eq.3) then
  
     endif
   !
   ! 0th value: interpolation value
   elseif (vorder.eq.0) then
   !
     if (pw%order.eq.0) then
       value = pw%y(kp)
     elseif (pw%order.eq.1) then
       b =(x-pw%cx(kp-1))/pw%dx(kp)
       value = 2.0_r8*(pw%y(kp)-pw%mcy(1,kp))*b+pw%mcy(1,kp)
     elseif (pw%order.eq.2) then
       b =(x-pw%cx(kp-1))/pw%dx(kp)
       f6 = 6.0_r8*(pw%y(kp)-0.5_r8*(pw%mcy(1,kp)+pw%mcy(2,kp)))
       value = pw%mcy(1,kp)+b*(pw%dy(kp)+f6*(1.0_r8-b))
     elseif (pw%order.eq.3) then
     !
     endif
   !
   ! 1th value: y'
   elseif (vorder.eq.1) then
   !
     if (pw%order.eq.0) then
       value = 0.0
     elseif (pw%order.eq.1) then
       value = 2.0_r8*(pw%y(kp)-pw%mcy(1,kp))/pw%dx(kp)
     elseif (pw%order.eq.2) then
       b =(x-pw%cx(kp-1))/pw%dx(kp)
       f6 = 6.0_r8*(pw%y(kp)-0.5_r8*(pw%mcy(1,kp)+pw%mcy(2,kp)))
       value = (pw%dy(kp)+f6-2.0*f6*b)/pw%dx(kp)
     elseif (pw%order.eq.3) then
     !
     endif
   !
   ! 2th value: y''
   elseif (vorder.eq.1) then
   !
     if (pw%order.eq.0) then
       value = 0.0
     elseif (pw%order.eq.1) then
       value = 0.0
     elseif (pw%order.eq.2) then
       f6 = 6.0_r8*(pw%y(kp)-0.5_r8*(pw%mcy(1,kp)+pw%mcy(2,kp)))
       value = -2.0*f6/pw%dx(kp)**2.0
     elseif (pw%order.eq.3) then
     !
     endif
!
   endif ! vorder.eq.0
!
   end function get_spline_value_in_piece
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine piecewise_get_fun(pw, m, x, y, vorder)
!-------------------------------------------------------------------------------
   implicit none
!
   type(piecewise_t),      intent(inout) :: pw
   integer(i4),            intent(in   ) :: m
   real(r8), dimension(m), intent(in   ) :: x
   real(r8), dimension(m), intent(  out) :: y
   integer(i4), optional,  intent(in   ) :: vorder ! 0:function, 1:first derivative, 2:second derivative
! local variables
   integer(i4) :: kp, i, n, l_vorder
   logical(l4) :: in_src
!
   if (present(vorder)) then
     l_vorder = vorder
   else
     l_vorder = 0
   endif
   if (l_vorder.ne.-1.and.l_vorder.ne.0.and.l_vorder.ne.1.and.l_vorder.ne.2.and.l_vorder.ne.3) then
     print *,'vorder must be -1,0,1,2,3 in piecewise_get_fun...'
     stop
   endif
!
   n = pw%n
!
! check boundary
   in_src = check_array_boundary(pw%cx(0:n),x,pw%array_order)
   if (.not.in_src) then
     print*,'error: check boundary.... in piecewise_get_fun'
     stop
   endif
!
   do i = 1,m
     kp   = get_idx_in_array(n+1,pw%cx(0:n),x(i),pw%array_order)
     y(i) = get_spline_value_in_piece(pw,kp,x(i),l_vorder)
   enddo
!
   return
   end subroutine piecewise_get_fun
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine piecewise_get_fun_integral(pw, m, x, y)
!-------------------------------------------------------------------------------
   implicit none
!
   type(piecewise_t),      intent(inout) :: pw
   integer(i4),            intent(in   ) :: m
   real(r8), dimension(m), intent(in   ) :: x
   real(r8), dimension(m), intent(  out) :: y
! local variables
   integer(i4) :: kp(m)
   integer(i4) :: nk, k, i, n
   real(r8) :: areak, areakp1, areai
!
   n = pw%n
!
   do i = 1,m
     kp(i) = get_idx_in_array(n+1,pw%cx(0:n),x(i),pw%array_order)
   enddo
!
   y(1) = 0.0_r8
   do i = 1,m-1
!
     if (kp(i).eq.kp(i+1)) then
       areak   = get_spline_value_in_piece(pw,kp(i),x(i),-1)
       areakp1 = get_spline_value_in_piece(pw,kp(i+1),x(i+1),-1)
       y(i+1)  = areakp1-areak
     else ! if (kp(i).ne.kp(i+1)) then
       areai = 0.0_r8
       do k = kp(i),kp(i+1)-1
         areai = areai+pw%inty(k)
       enddo
       areak   = get_spline_value_in_piece(pw,kp(i),x(i),-1)
       areakp1 = get_spline_value_in_piece(pw,kp(i+1),x(i+1),-1)
       y(i+1)  = areakp1+areai-areak
     endif ! if (kp(i).ne.kp(i+1)) then
!
   enddo ! do i=1,m-1
!
   return
   end subroutine piecewise_get_fun_integral
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine piecewise_get(pw, m, cx, y)
!-------------------------------------------------------------------------------
   implicit none
!
   type(piecewise_t),        intent(inout) :: pw
   integer(i4),              intent(in   ) :: m
   real(r8), dimension(0:m), intent(in   ) :: cx
   real(r8), dimension(m),   intent(  out) :: y
! local variables
   real(r8), dimension(0:m) :: icy
   integer(i4) :: k, n, array_order
   logical(l4) :: in_src
!
   n = pw%n
!
! check array order
   array_order = check_array_order(cx)
   if (pw%array_order.ne.array_order) then
     print *,'error: check array order (target)',pw%array_order,array_order
     stop
   endif
! check boundary
   in_src = check_array_boundary(pw%cx(0:n),cx,array_order)
   if (.not.in_src) then
     print*,'error: check boundary.... in piecewise_get'
     stop
   endif
!
! integration
   call piecewise_get_fun_integral(pw,m+1,cx,icy)
   do k = 1,m
     y(k) = icy(k)/(cx(k)-cx(k-1))
   enddo
!
   return
   end subroutine piecewise_get
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine piecewise_interface(n, scx, sy, m, tcx, ty, method, ismono, bndry)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(0:n), intent(in   ) :: scx
   real(r8), dimension(n),   intent(in   ) :: sy
   integer(i4),              intent(in   ) :: m
   real(r8), dimension(0:m), intent(in   ) :: tcx
   real(r8), dimension(m),   intent(  out) :: ty
   ! method = 0:constant, 1:linear, 2:parabolic
   integer(i4), optional,    intent(in   ) :: method
   logical(i4), optional,    intent(in   ) :: ismono
   integer(i4), optional,    intent(in   ) :: bndry
! local variables
   type(piecewise_t) :: pw
   real(r8)          :: vmin
!
   call piecewise_initialize(pw,n,method,ismono,bndry)
   call piecewise_set_x(pw,scx)
   call piecewise_set(pw,sy)
   call piecewise_get(pw,m,tcx,ty)
   call piecewise_finalize(pw)
!
   vmin = 1.0d-100
   ty   = merge(ty,0.0_r8,abs(ty).ge.vmin)
!
   return
   end subroutine piecewise_interface
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine piecewise_cyclic_interface(n, scx, sy, m, tcx, ty, method, ismono, bndry)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(0:n), intent(in   ) :: scx
   real(r8), dimension(n),   intent(in   ) :: sy
   integer(i4),              intent(in   ) :: m
   real(r8), dimension(0:m), intent(in   ) :: tcx
   real(r8), dimension(m),   intent(  out) :: ty
   ! method = 0:constant, 1:linear, 2:parabolic
   integer(i4), optional,    intent(in   ) :: method
   logical(i4), optional,    intent(in   ) :: ismono
   integer(i4), optional,    intent(in   ) :: bndry
! local variables
   type(piecewise_t) :: pw
   real(r8)       :: vmin, lbase, rbase
   integer(i4)    :: i, nn, array_order
   real(r8), allocatable, dimension(:) :: escx, esy
!
   nn = 3*n
   allocate(escx(0:nn),esy(nn))
!
! expand x (source)
   array_order = check_array_order(scx)
   if (array_order.eq.0) then
     print*,'check array order of source in piecewise_cyclic_interface.'
     stop
   endif
   rbase           = scx(n)-scx(0)
   lbase           = scx(0)-scx(n)
   escx(n:2*n)     = scx(:)
   escx(2*n+1:3*n) = scx(1:n)+rbase   ! right
   escx(0:n-1)     = scx(0:n-1)+lbase ! left
!
! expand y (source)
   esy(n+1:2*n)   = sy(:)
   esy(1:n)       = sy(1:n)
   esy(2*n+1:3*n) = sy(1:n)
!
   call piecewise_initialize(pw,nn,method,ismono,2)
   call piecewise_set_x(pw,escx)
   call piecewise_set(pw,esy)
   call piecewise_get(pw,m,tcx,ty)
   call piecewise_finalize(pw)
!
   vmin = 1.0d-100
   ty   = merge(ty,0.0_r8,abs(ty).ge.vmin)
!
   deallocate(escx,esy)
!
   return
   end subroutine piecewise_cyclic_interface
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
   end module piecewise_remap
!-------------------------------------------------------------------------------
