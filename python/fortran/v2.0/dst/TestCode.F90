   module mymod
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module mymod

!>  - JH KIM  : Added to "ghost_exchangeV" subroutine. (CAM-SE 5.3)
!-------------------------------------------------------------------------------

#include "KIM.h"

   module testcode
   use mymodules1, only : var11, var12, var13
   use mymodules2, only : var11, var12, var13
   use mymodules3, only : fjdkl, fjdkla, gja,   fjdla,   fdjalk, fjklda
!-------------------------------------------------------------------------------
   implicit none
!
   private

   type my_type
   integer(i4) :: a, b
   real(r8) :: fjdk, fjdkla
   end type

   public :: ghost_exchangev_new


   interface my_interface
   module procedure my_sub1
   module procedure my_sub2
   module procedure my_sub3
   end interface myinterface


   contains

!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine my_sub()
   a = 2
!
   return
   end subroutine

!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine my_sub0()
!
   return
   end subroutine


!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine my_sub1(a, b, c)
   use mymodule3, only : j, fj, jfkd
!-------------------------------------------------------------------------------
   implicit none
!
   type my_sub_t
   real :: abc
   end type my_sub_t
   type(my_sub_t) :: mytt
   real(r8), dimension(:,:), allocatable :: a
   integer(i4) :: b, c

   ii = 10
   do i = 1, 100
   jfdk = fjdk
   enddo

   do j = 1, 100
   do i = 1, 100
   b = abd
   enddo
   enddo

   do k = 1, 100
   do j = 1, 100
   do i = 1, 100
   b = abd
   enddo
   enddo
   enddo

   ii = 1
   if (ii.eq.1) some = some1
   if (ii.ne.1) then
   a = b
   do i = 1, n
   a = b
   do j = 1, m
   f = 1
   enddo
   a = b
   enddo
   a = b
   endif

   if (a==1) then
   b = 2
   elseif (a==2) then
   b = 3
   else
   b = 0
   endif

   return
!
   return
   end subroutine my_sub1





!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function my_fun(a, b) result(d)
   use mymod1, only : dj, fj
   use mymod2
!-------------------------------------------------------------------------------
   implicit none
!
   include 'jfkdla.inc'
   integer(i4), dimension(:),                       intent(inout) :: a, b
!local
   integer(i4) :: i, n

   n = size(a)
   d = 0
   do i = 1, n
   d = a(i)+b(i)
   enddo
!
   return
   end function

!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine testtest

   if (a.and.b) then
   endif
   if (a.and.b) then
   endif
   if (a.and.b) then
   endif
!
   return
   end subroutine testtest

!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module
