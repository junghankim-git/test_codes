!-------------------------------------------------------------------------------
   module gen_test_func
!-------------------------------------------------------------------------------
!
!  abstract :  bi-linear interpolation module
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-15  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds,    only: i4, l4, r4, r8
!-------------------------------------------------------------------------------
   implicit none
!
   private
   real(r8), parameter :: pi = 3.14159265358979323846264_r8
!
   public :: get_func, get_func_area
   public :: get_poly1, get_poly2, get_poly3
   public :: get_harmonic, get_gaussian, get_step
   public :: get_poly1_area, get_poly2_area, get_poly3_area
   public :: get_harmonic_area, get_gaussian_area, get_step_area
   public :: get_random
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function func_poly1(xx, xmin, xmax, mu, order) result(yy)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)   , intent(in   ) :: xx
   real(r8)   , intent(in   ) :: xmin, xmax, mu
   integer(i4), intent(in   ) :: order
   real(r8)                   :: yy
! local variables
   real(r8) :: a
!
!  y = 2(x-mu)
!
   a = 2.0
!
   if (order.eq.0) then
     yy = a*(xx-mu)
   elseif (order.eq.1) then
     yy = a
   elseif (order.eq.2) then
     yy = 0.0
   elseif (order.eq.-1) then
     yy = 0.5*a*(xx-mu)**2.0-0.5*a*(xmin-mu)**2.0
   else
     print*, 'check order... '
   endif
!
   end function func_poly1
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function func_poly2(xx, xmin, xmax, mu, order) result(yy)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)   , intent(in   ) :: xx
   real(r8)   , intent(in   ) :: xmin, xmax, mu
   integer(i4), intent(in   ) :: order
   real(r8)                   :: yy
! local variables
!
!  y = (x-mu)^2 - 1
!
   if (order.eq.0) then
     yy = 1.0*(xx-mu)**2.0-1.0
   elseif (order.eq.1) then
     yy = 2.0*(xx-mu)
   elseif (order.eq.2) then
     yy = 2.0
   elseif (order.eq.-1) then
     yy = 1.0/3.0*(xx-mu)**3.0-xx-(1.0/3.0*(xmin-mu)**3.0-xmin)
   else
     print*, 'check order... '
   endif
!
   end function func_poly2
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function func_poly3(xx, xmin, xmax, mu, order) result(yy)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)   , intent(in   ) :: xx
   real(r8)   , intent(in   ) :: xmin, xmax, mu
   integer(i4), intent(in   ) :: order
   real(r8)                   :: yy
! local variables
!
!  y = (x-mu)^3 - 1
!
   if (order.eq.0) then
     yy =(xx-mu)**3.0-1.0
   elseif (order.eq.1) then
     yy = 3.0*(xx-mu)**2.0
   elseif (order.eq.2) then
     yy = 6.0*(xx-mu)
   elseif (order.eq.-1) then
     yy = 1.0/4.0*(xx-mu)**4.0-xx-(1.0/4.0*(xmin-mu)**4.0-xmin)
   else
     print*, 'check order... '
   endif
!
   end function func_poly3
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function func_harmonic(xx, xmin, xmax, mu, order) result(yy)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)   , intent(in   ) :: xx
   real(r8)   , intent(in   ) :: xmin, xmax, mu
   integer(i4), intent(in   ) :: order
   real(r8)                   :: yy
! local variables
   real(r8)    :: a, x
   integer(i4) :: i
!
!  y = sin( 4pi(x-xmin)/(xmax-xmin) )
!
   a = 4.0*pi/(xmax-xmin)
   x = a*(xx-xmin)
   if (order.eq.0) then
     yy = sin(x)
   elseif (order.eq.1) then
     yy = a*cos(x)
   elseif (order.eq.2) then
     yy =-a*a*sin(x)
   elseif (order.eq.-1) then
     yy =-1.0/a*(cos(x)-1.0)
   else
     print*, 'check order... '
   endif
!
   end function func_harmonic
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function func_gaussian(xx, xmin, xmax, mu, order) result(yy)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)   , intent(in   ) :: xx
   real(r8)   , intent(in   ) :: xmin, xmax, mu
   integer(i4), intent(in   ) :: order
   real(r8)                   :: yy
! local variables
   real(r8)    :: a, sig, x, alpha, tint, iint
   integer(i4) :: i
!
!  y = 
!
   a     = 1.0
   sig   = 1.0/sqrt(2.0)*(xmax-xmin)/(4.0*pi)
   alpha = 0.5/sig**2.0
!
   tint = a*sqrt(pi/alpha)
!
   x = xx-mu
   if (order.eq.0) then
     yy = a*exp(-alpha*x**2)
   elseif (order.eq.1) then
     yy =-2.0*a*alpha*x*exp(-alpha*x**2)
   elseif (order.eq.2) then
     yy = a*(2.0*alpha*x**2.0-1.0)/sig**2.0*exp(-alpha*x**2)
   elseif (order.eq.-1) then
     iint = tint*sqrt(1.0-exp(-sqrt(2.0)*alpha*x**2))
     if (x.le.0) then
       yy = 0.5*(tint-iint)
     else
       yy = 0.5*(tint+iint)
     endif
   else
     print*, 'check order... '
   endif
!
   end function func_gaussian
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function func_step(xx, xmin, xmax, mu, order) result(yy)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)   , intent(in   ) :: xx
   real(r8)   , intent(in   ) :: xmin, xmax, mu
   integer(i4), intent(in   ) :: order
   real(r8)                   :: yy
! local variables
!
   if (order.eq.0) then
     if (xx.lt.mu) then
       yy = 0.0
     else
       yy = 1.0
     endif
   elseif (order.eq.1) then
     yy = 0.0
   elseif (order.eq.2) then
     yy = 0.0
   elseif (order.eq.-1) then
     if (xx.lt.mu) then
       yy = 0.0
     else
       yy =(xx-mu)
     endif
   else
     print*,'check order... '
   endif
!
   end function func_step
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function func_some(xx, xmin, xmax, mu, order) result(yy)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8)   , intent(in   ) :: xx
   real(r8)   , intent(in   ) :: xmin, xmax, mu
   integer(i4), intent(in   ) :: order
   real(r8)                   :: yy
! local variables
   integer(i4), parameter :: n = 50
   integer(i4) :: i, ik
   real(r8)    :: dx, cx(0:n), x(n), f(n), area(n)
!
   dx = (xmax-xmin)/real(n,r8)
   cx(0) = xmin
   do i=1,n
     cx(i) = cx(i-1)+dx
     x(i)  = 0.5*(cx(i)+cx(i-1))
   enddo
   f(01)=0.00; f(02)=0.00; f(03)=0.00; f(04)=0.00; f(05)=0.15 ! 0-1
   f(06)=0.55; f(07)=0.75; f(08)=0.76; f(09)=0.74; f(10)=0.78 ! 1-2
   f(11)=0.82; f(12)=0.88; f(13)=0.92; f(14)=0.93; f(15)=0.81 ! 2-3
   f(16)=0.60; f(17)=0.45; f(18)=0.40; f(19)=0.42; f(20)=0.50 ! 3-4
   f(21)=0.60; f(22)=0.68; f(23)=0.68; f(24)=0.58; f(25)=0.51 ! 4-5
   f(26)=0.56; f(27)=0.71; f(28)=0.70; f(29)=0.50; f(30)=0.25 ! 5-6
   f(31)=0.15; f(32)=0.10; f(33)=0.02; f(34)=0.00; f(35)=0.00 ! 6-7
   f(36)=0.00; f(37)=0.00; f(38)=0.00; f(39)=0.00; f(40)=0.00 ! 7-8
   f(41)=0.00; f(42)=0.00; f(43)=0.00; f(44)=0.00; f(45)=0.00 ! 8-9
   f(46)=0.00; f(47)=0.00; f(48)=0.00; f(49)=0.00; f(50)=0.00 ! 9-10
   ik = get_idx_in_array(n+1,cx,xx,1)
   do i=1,n
     area(i) = f(i)*dx
   enddo
!
   if (order.eq.0) then
     yy = f(ik)
   elseif (order.eq.1) then
     yy = 0.0
   elseif (order.eq.2) then
     yy = 0.0
   elseif (order.eq.-1) then
     yy = 0.0
     do i=1,ik-1
       yy = yy+area(i)
     enddo
     yy = yy+f(ik)*(xx-cx(i-1))
   else
     print*,'check order... '
   endif
!
   end function func_some
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_func(ifun, nn, xx, yy, xmin, xmax, mu, order)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),             intent(in   ) :: ifun
   integer(i4),             intent(in   ) :: nn
   real(r8), dimension(nn), intent(in   ) :: xx
   real(r8), dimension(nn), intent(  out) :: yy
   real(r8),                intent(in   ) :: xmin, xmax, mu
   integer(i4),             intent(in   ) :: order
! local variables
   integer(i4) :: i
!
   do i = 1, nn
     if (ifun.eq.1) then
       yy(i) = func_poly1(xx(i),xmin,xmax,mu,order)
     elseif (ifun.eq.2) then
       yy(i) = func_poly2(xx(i),xmin,xmax,mu,order)
     elseif (ifun.eq.3) then
       yy(i) = func_poly3(xx(i),xmin,xmax,mu,order)
     elseif (ifun.eq.4) then
       yy(i) = func_harmonic(xx(i),xmin,xmax,mu,order)
     elseif (ifun.eq.5) then
       yy(i) = func_gaussian(xx(i),xmin,xmax,mu,order)
     elseif (ifun.eq.6) then
       yy(i) = func_step(xx(i),xmin,xmax,mu,order)
     elseif (ifun.eq.7) then
       yy(i) = func_some(xx(i),xmin,xmax,mu,order)
     else
       print *,'check ifun...'
       stop
     endif
   enddo
!
   return
   end subroutine get_func
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_func_area(ifun, nn, cxx, yy, xmin, xmax, mu)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: ifun
   integer(i4),               intent(in   ) :: nn
   real(r8), dimension(0:nn), intent(in   ) :: cxx
   real(r8), dimension(nn),   intent(  out) :: yy
   real(r8),                  intent(in   ) :: xmin, xmax, mu
! local variables
   integer(i4) :: i
   real(r8)    :: a2, a1
!
   do i = 1,nn
     if (ifun.eq.1) then
       a2 = func_poly1(cxx(i),xmin,xmax,mu,-1)
       a1 = func_poly1(cxx(i-1),xmin,xmax,mu,-1)
     elseif (ifun.eq.2) then
       a2 = func_poly2(cxx(i),xmin,xmax,mu,-1)
       a1 = func_poly2(cxx(i-1),xmin,xmax,mu,-1)
     elseif (ifun.eq.3) then
       a2 = func_poly3(cxx(i),xmin,xmax,mu,-1)
       a1 = func_poly3(cxx(i-1),xmin,xmax,mu,-1)
     elseif (ifun.eq.4) then
       a2 = func_harmonic(cxx(i),xmin,xmax,mu,-1)
       a1 = func_harmonic(cxx(i-1),xmin,xmax,mu,-1)
     elseif (ifun.eq.5) then
       a2 = func_gaussian(cxx(i),xmin,xmax,mu,-1)
       a1 = func_gaussian(cxx(i-1),xmin,xmax,mu,-1)
     elseif (ifun.eq.6) then
       a2 = func_step(cxx(i),xmin,xmax,mu,-1)
       a1 = func_step(cxx(i-1),xmin,xmax,mu,-1)
     elseif (ifun.eq.7) then
       a2 = func_some(cxx(i),xmin,xmax,mu,-1)
       a1 = func_some(cxx(i-1),xmin,xmax,mu,-1)
     else
       print *,'check ifun...'
       stop
     endif
     yy(i) = (a2-a1)/(cxx(i)-cxx(i-1))
   enddo
!
   return
   end subroutine get_func_area
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_poly1(nn, xx, yy, xmin, xmax, mu, order)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),             intent(in   ) :: nn
   real(r8), dimension(nn), intent(in   ) :: xx
   real(r8), dimension(nn), intent(  out) :: yy
   real(r8),                intent(in   ) :: xmin, xmax, mu
   integer(i4),             intent(in   ) :: order
! local variables
   integer(i4) :: i
!
   do i = 1, nn
     yy(i) = func_poly1(xx(i),xmin,xmax,mu,order)
   enddo
!
   return
   end subroutine get_poly1
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_poly2(nn, xx, yy, xmin, xmax, mu, order)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),             intent(in   ) :: nn
   real(r8), dimension(nn), intent(in   ) :: xx
   real(r8), dimension(nn), intent(  out) :: yy
   real(r8),                intent(in   ) :: xmin, xmax, mu
   integer(i4),             intent(in   ) :: order
! local variables
   integer(i4) :: i
!
   do i = 1, nn
     yy(i) = func_poly2(xx(i),xmin,xmax,mu,order)
   enddo
!
   return
   end subroutine get_poly2
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_poly3(nn, xx, yy, xmin, xmax, mu, order)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),             intent(in   ) :: nn
   real(r8), dimension(nn), intent(in   ) :: xx
   real(r8), dimension(nn), intent(  out) :: yy
   real(r8),                intent(in   ) :: xmin, xmax, mu
   integer(i4),             intent(in   ) :: order
! local variables
   integer(i4) :: i
!
   do i = 1, nn
     yy(i) = func_poly3(xx(i),xmin,xmax,mu,order)
   enddo
!
   return
   end subroutine get_poly3
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_harmonic(nn, xx, yy, xmin, xmax, mu, order)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),             intent(in   ) :: nn
   real(r8), dimension(nn), intent(in   ) :: xx
   real(r8), dimension(nn), intent(  out) :: yy
   real(r8),                intent(in   ) :: xmin, xmax, mu
   integer(i4),             intent(in   ) :: order
! local variables
   integer(i4) :: i
!
   do i = 1, nn
     yy(i) = func_harmonic(xx(i),xmin,xmax,mu,order)
   enddo
!
   return
   end subroutine get_harmonic
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_gaussian(nn, xx, yy, xmin, xmax, mu, order)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),             intent(in   ) :: nn
   real(r8), dimension(nn), intent(in   ) :: xx
   real(r8), dimension(nn), intent(  out) :: yy
   real(r8),                intent(in   ) :: xmin, xmax, mu
   integer(i4),             intent(in   ) :: order
! local variables
   integer(i4) :: i
!
   do i = 1, nn
     yy(i) = func_gaussian(xx(i),xmin,xmax,mu,order)
   enddo
!
   return
   end subroutine get_gaussian
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_step(nn, xx, yy, xmin, xmax, mu, order)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),             intent(in   ) :: nn
   real(r8), dimension(nn), intent(in   ) :: xx
   real(r8), dimension(nn), intent(  out) :: yy
   real(r8),                intent(in   ) :: xmin, xmax, mu
   integer(i4),             intent(in   ) :: order
! local variables
   integer(i4) :: i
!
   do i = 1,nn
     yy(i) = func_step(xx(i),xmin,xmax,mu,order)
   enddo
!
   return
   end subroutine get_step
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_poly1_area(nn, cxx, yy, xmin, xmax, mu)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: nn
   real(r8), dimension(0:nn), intent(in   ) :: cxx
   real(r8), dimension(nn),   intent(  out) :: yy
   real(r8),                  intent(in   ) :: xmin, xmax, mu
! local variables
   integer(i4) :: i
   real(r8)    :: a2, a1
!
   do i = 1,nn
     a2    = func_poly1(cxx(i),xmin,xmax,mu,-1)
     a1    = func_poly1(cxx(i-1),xmin,xmax,mu,-1)
     yy(i) = (a2-a1)/(cxx(i)-cxx(i-1))
   enddo
!
   return
   end subroutine get_poly1_area
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_poly2_area(nn, cxx, yy, xmin, xmax, mu)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: nn
   real(r8), dimension(0:nn), intent(in   ) :: cxx
   real(r8), dimension(nn),   intent(  out) :: yy
   real(r8),                  intent(in   ) :: xmin, xmax, mu
! local variables
   integer(i4) :: i
   real(r8)    :: a2, a1
!
   do i = 1,nn
     a2    = func_poly2(cxx(i),xmin,xmax,mu,-1)
     a1    = func_poly2(cxx(i-1),xmin,xmax,mu,-1)
     yy(i) = (a2-a1)/(cxx(i)-cxx(i-1))
   enddo
!
   return
   end subroutine get_poly2_area
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_poly3_area(nn, cxx, yy, xmin, xmax, mu)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: nn
   real(r8), dimension(0:nn), intent(in   ) :: cxx
   real(r8), dimension(nn),   intent(  out) :: yy
   real(r8),                  intent(in   ) :: xmin, xmax, mu
! local variables
   integer(i4) :: i
   real(r8)    :: a2, a1
!
   do i = 1,nn
     a2    = func_poly3(cxx(i),xmin,xmax,mu,-1)
     a1    = func_poly3(cxx(i-1),xmin,xmax,mu,-1)
     yy(i) = (a2-a1)/(cxx(i)-cxx(i-1))
   enddo
!
   return
   end subroutine get_poly3_area
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_harmonic_area(nn, cxx, yy, xmin, xmax, mu)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: nn
   real(r8), dimension(0:nn), intent(in   ) :: cxx
   real(r8), dimension(nn),   intent(  out) :: yy
   real(r8),                  intent(in   ) :: xmin, xmax, mu
! local variables
   integer(i4) :: i
   real(r8)    :: a2, a1
!
   do i = 1,nn
     a2    = func_harmonic(cxx(i),xmin,xmax,mu,-1)
     a1    = func_harmonic(cxx(i-1),xmin,xmax,mu,-1)
     yy(i) = (a2-a1)/(cxx(i)-cxx(i-1))
   enddo
!
   return
   end subroutine get_harmonic_area
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_gaussian_area(nn, cxx, yy, xmin, xmax, mu)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: nn
   real(r8), dimension(0:nn), intent(in   ) :: cxx
   real(r8), dimension(nn),   intent(  out) :: yy
   real(r8),                  intent(in   ) :: xmin, xmax, mu
! local variables
   integer(i4) :: i
   real(r8)    :: a2, a1
!
   do i = 1,nn
     a2    = func_gaussian(cxx(i),xmin,xmax,mu,-1)
     a1    = func_gaussian(cxx(i-1),xmin,xmax,mu,-1)
     yy(i) = (a2-a1)/(cxx(i)-cxx(i-1))
   enddo
!
   return
   end subroutine get_gaussian_area
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_step_area(nn, cxx, yy, xmin, xmax, mu)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: nn
   real(r8), dimension(0:nn), intent(in   ) :: cxx
   real(r8), dimension(nn),   intent(  out) :: yy
   real(r8),                  intent(in   ) :: xmin, xmax, mu
! local variables
   integer(i4) :: i
   real(r8)    :: a2, a1
!
   do i = 1,nn
     a2    = func_step(cxx(i),xmin,xmax,mu,-1)
     a1    = func_step(cxx(i-1),xmin,xmax,mu,-1)
     yy(i) = (a2-a1)/(cxx(i)-cxx(i-1))
   enddo
!
   return
   end subroutine get_step_area
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_random(nn, xx, yy, xmin, xmax, mu, order)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),             intent(in   ) :: nn
   real(r8), dimension(nn), intent(in   ) :: xx
   real(r8), dimension(nn), intent(  out) :: yy
   real(r8),                intent(in   ) :: xmin, xmax, mu
   integer(i4),             intent(in   ) :: order
! local variables
   integer(i4) :: i, n
!
   call random_seed()
   do i = 1, nn
     if (order.eq.0) then
       call random_number(yy(i))
     elseif (order.eq.1) then
       yy(i) = 0.0_r8
     elseif (order.eq.2) then
       yy(i) = 0.0_r8
     else
       print*, 'check order... '
     endif
   enddo
!
   return
   end subroutine get_random
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
   end module gen_test_func
!-------------------------------------------------------------------------------
