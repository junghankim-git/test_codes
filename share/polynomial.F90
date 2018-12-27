!-------------------------------------------------------------------------------
   module polynomial
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
   use kinds, only : i4, r4, r8
!-------------------------------------------------------------------------------
   implicit none
!
   ! nl : # of points
   ! ndx: # of points between interfaces (must be even number)
   integer(i4) :: nl, ndx
!
   ! \eta_i, \eta_m
   real(r8), dimension(:), allocatable :: x_i, x_m
   real(r8), dimension(:), allocatable :: arg_i, arg_m
   ! for LU Factorization
   real(r8), dimension(:,:), allocatable, public :: lu_i ! n+1, n+1
   real(r8), dimension(:,:), allocatable, public :: lu_m ! n, n
   real(r8), dimension(:,:), allocatable, public :: luinv_i ! n+1, n+1
   real(r8), dimension(:,:), allocatable, public :: luinv_m ! n, n
!
   public :: inipoly, finpoly, poly, diffpoly, intpoly
!
   contains
!
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine inipoly(n, xin_i, xin_m)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n+1), intent(in   ) :: xin_i
   real(r8), dimension(n),   intent(in   ) :: xin_m
! local variables
   real(r8), dimension(n+1, n+1) :: ainv_i, tmp_i
   real(r8), dimension(n, n) :: ainv_m, tmp_m
   real(r8), dimension(n+1) :: ipivot_i
   real(r8), dimension(n) :: ipivot_m
   real(r8), dimension(n+1) :: work_i
   real(r8), dimension(n) :: work_m
   integer(i4) :: i, j, info
!
   nl = n
!
   allocate(x_i(n+1))
   allocate(x_m(n))
   allocate(arg_i(n+1))
   allocate(arg_m(n))
   ! for LU factorization
   allocate(lu_i(n+1,n+1))
   allocate(lu_m(n,n))
   allocate(luinv_i(n+1,n+1))
   allocate(luinv_m(n,n))
!
   do i = 1,n
     x_i(i) = xin_i(i)
     x_m(i) = xin_m(i)
   enddo
   x_i(n+1) = xin_i(n+1)
!
   do j = 1,n+1
   do i = 1,n+1
     lu_i(i,j) = xin_i(i)**dble(j-1)
     ainv_i(i,j) = lu_i(i,j)
   enddo
   enddo
!
   do j = 1,n
   do i = 1,n
     lu_m(i,j) = xin_m(i)**dble(j-1)
     ainv_m(i,j) = lu_m(i,j)
   enddo
   enddo
! factorization & inverse matrix
     call dgetrf(n+1,n+1,ainv_i,n+1,ipivot_i,info)
   if (info.ne.0) stop
     call dgetri(n+1,ainv_i,n+1,ipivot_i,work_i,n+1,info)
   if (info.ne.0) stop
     call dgetrf(n,n,ainv_m,n,ipivot_m,info)
   if (info.ne.0) stop
     call dgetri(n,ainv_m,n,ipivot_m,work_m,n,info)
   if (info.ne.0) stop
!
   do j = 1,n+1
   do i = 1,n+1
     luinv_i(i,j) = ainv_i(i,j)
   enddo
   enddo
!
   do j = 1,n
   do i = 1,n
     luinv_m(i,j) = ainv_m(i,j)
   enddo
   enddo
!
   return
   end subroutine inipoly
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine finpoly()
!
   deallocate(x_i)
   deallocate(x_m)
   deallocate(arg_i)
   deallocate(arg_m)
   deallocate(lu_i)
   deallocate(lu_m)
   deallocate(luinv_i)
   deallocate(luinv_m)
!
   return
   end subroutine finpoly
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine argpoly(n, val) !, arg)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),            intent(in   ) :: n
   real(r8), dimension(n), intent(in   ) :: val
! local variables
   integer(i4) :: i, k
!
   if (n.eq.nl+1) then
     do i = 1,n
     do k = 1,n
       arg_i(i) = arg_i(i)+luinv_i(i,k)*val(k)
     enddo
     enddo
   elseif (n.eq.nl) then
     do i = 1,n
     do k = 1,n
       arg_m(i) = arg_m(i)+luinv_m(i,k)*val(k)
     enddo
     enddo
   endif
!
   return
   end subroutine argpoly
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function poly(n, x)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: n
   real(r8),    intent(in   ) :: x
   real(r8) :: poly
! local variables
   integer(i4) :: i, k
!
   poly = 0.0_r8
!
   if (n.eq.nl+1) then
     do i = 1,n
       poly = poly+arg_i(i)*(x**dble(i-1))
     enddo
   elseif (n.eq.nl) then
     do i = 1,n
       poly = poly+arg_m(i)*(x**dble(i-1))
     enddo
   endif
!
   end function poly
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!  Return dL_m(x)/dx
   function diffpoly(n, x)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: n
   real(r8),    intent(in   ) :: x
   real(r8) :: diffpoly
! local variables
   integer(i4) :: i, k
!
   diffpoly = 0.0_r8
!
   if (n.eq.nl+1) then
     do i = 1,n
       diffpoly = diffpoly+dble(i-1)*arg_i(i)*(x**dble(i-2))
     enddo
   elseif (n.eq.nl) then
     do i = 1,n
       diffpoly = diffpoly+dble(i-1)*arg_m(i)*(x**dble(i-2))
     enddo
   endif
!
   end function diffpoly
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!  return \int_0^x_k {l_m(x)} dx  (x_k : full or half level index)
   function intpoly(n, x)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: n
   real(r8),    intent(in   ) :: x
   real(r8) :: intpoly
! local variables
   integer(i4) :: i, k
!
   intpoly = 0.0_r8
!
   if (n.eq.nl+1) then
     do i = 1,n
       intpoly = intpoly+arg_i(i)*(x**dble(i))/dble(i)
     enddo
   elseif (n.eq.nl) then
     do i = 1,n
       intpoly = intpoly+arg_m(i)*(x**dble(i))/dble(i)
     enddo
   endif
!
   end function intpoly
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module polynomial
