!-------------------------------------------------------------------------------
   module cubic_spline_interp
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    2015-08-02  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only : i4, r4, r8
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4) :: n
!
   real(r8), dimension(:), allocatable :: xx, yy, y2
!
   public :: cs_interp_initialize, cs_interp_finalize, cs_interp_get
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_interp_initialize(n_in, x_in, y_in, yp1, ypn)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: n_in
   real(r8), dimension(n_in), intent(in   ) :: x_in, y_in
   real(r8),                  intent(in   ) :: yp1, ypn
! local variables
   integer(i4) :: i, k
!
   n = n_in
!
   allocate(xx(n))
   allocate(yy(n))
   allocate(y2(n))
!
   call setcubicspline(x_in,y_in,yp1,ypn)
!
   return
   end subroutine cs_interp_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine setcubicspline(x_in, y_in, yp1, ypn)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(n), intent(in   ) :: x_in, y_in
   real(r8),               intent(in   ) :: yp1, ypn
! local variables
   integer(i4) :: i, k
   real(r8) :: p, qn, sig, un, u(200)
!
   do i = 1,n
     xx(i) = x_in(i)
     yy(i) = y_in(i)
   enddo
!
   if (yp1.gt..99e30) then
     y2(1) = 0.0_r8
     u(1)  = 0.0_r8
   else
     y2(1) =-0.5_r8
     u(1)  = (3.0_r8/(xx(2)-xx(1)))*((yy(2)-yy(1))/(xx(2)-xx(1))-yp1)
   endif
!
   do i = 2,n-1
     sig =(xx(i)-xx(i-1))/(xx(i+1)-xx(i-1))
     p = sig*y2(i-1)+2.0_r8
     y2(i) =(sig-1.0_r8)/p
     u(i) =(6.0_r8*((yy(i+1)-yy(i))/(xx(i+1)-xx(i))-(yy(i)-yy(i-1))/(xx(i)-xx(i-1)))/(xx(i+1)-xx(i-1))-sig*u(i-1))/p
   enddo
!
   if (ypn.gt..99e30) then
     qn = 0.0_r8
     un=0.0_r8
   else
     qn = 0.5_r8
     un =(3.0_r8/(xx(n)-xx(n-1)))*(ypn-(yy(n)-yy(n-1))/(xx(n)-xx(n-1)))
   endif
!
   y2(n) =(un-qn*u(n-1))/(qn*y2(n-1)+1.0_r8)
!
   do k = n-1,1,-1
     y2(k) = y2(k)*y2(k+1)+u(k)
   enddo
!
   return
   end subroutine setcubicspline
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cs_interp_finalize()
!
   deallocate(xx)
   deallocate(yy)
   deallocate(y2)
!
   return
   end subroutine cs_interp_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function cs_interp_get(x_in)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: x_in
   real(r8) :: cs_interp_get
! local variables
   integer(i4) :: i, k, khi, klo
   real(r8) :: a, b, h
!
   klo = 1
   khi = n
1  if (khi-klo.gt.1) then
     k =(khi+klo)/2
     if (xx(k).gt.x_in) then
       khi = k
     else
       klo = k
     endif
     goto 1
   endif
!
   h = xx(khi)-xx(klo)
!
   if (h.eq.0.0_r8) then
     print*,'Bad input....'
   endif
!
   a = (xx(khi)-x_in)/h
   b = (x_in-xx(klo))/h
   cs_interp_get = a*yy(klo)+b*yy(khi)+((a**3-a)*y2(klo)+(b**3-b)*y2(khi))*(h**2)/6.0_r8
!
   end function cs_interp_get
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module cubic_spline_interp
