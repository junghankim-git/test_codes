!-------------------------------------------------------------------------------
:01:01:   module mymod
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:-1:00:   end module mymod
!-------------------------------------------------------------------------------
:-3:00:!>  - JH KIM  : Added to "ghost_exchangeV" subroutine. (CAM-SE 5.3)
:-3:00:!-------------------------------------------------------------------------------
:-3:00:#include "KIM.h"
!-------------------------------------------------------------------------------
:01:01:   module testcode
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  your   name    code comment
!
!  structure:
!
!-------------------------------------------------------------------------------
:03:03:   use mymodules1, only: var11, var12, var13
:03:03:   use mymodules2, only: var11, var12, var13
:03:03:   use mymodules3, only: fjdkl, fjdkla, gja, &
:-2:03:                                                                        fjdla, &
:-2:03:                                                                fdjalk, fjklda
:03:03:   implicit none
!
:04:04:   private
!
:05:05:   type my_type
:04:04:     integer(i4) :: a, b
:04:04:     integer(i4) :: c, d
:04:04:     integer(i4) :: e, f
:04:04:     integer(i4) :: g, h
:04:04:     real(r8) :: fjdk, fjdkla
:-1:04:   end type
:04:04:   public :: ghost_exchangev_new
!
:05:05:   interface my_interface
:04:04:     module procedure my_sub1
:04:04:     module procedure my_sub2
:04:04:     module procedure my_sub3
:-1:04:   end interface myinterface
!
:-1:01:   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   subroutine my_sub()
!
:06:06:   a = 2
!
:-1:01:   end subroutine
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   subroutine my_sub0()
!
:-1:01:   end subroutine
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   subroutine my_sub1(a, b, c)
!-------------------------------------------------------------------------------
:03:03:   use mymodule3, only: j, fj, jfkd
:03:03:   implicit none
!
:05:05:   type my_sub_t
:04:04:     real :: abc
:-1:03:   end type my_sub_t
!
:04:04:   type(my_sub_t) :: mytt
:04:04:   real(r8), dimension(:,:), allocatable :: a
:04:04:   integer(i4) :: b, c
!
:06:06:   a = 1
:05:05:   if (xxx) then
:-1:06:   endif
:06:06:   a = 1
:05:05:   if (xxx) then
:06:06:     a = 2
:-1:06:   endif
:06:06:   a = 1
:05:05:   if (xxx) then
:05:05:     if (xxx) then
:-1:05:     endif
:-1:06:   endif
:06:06:   a = 1
:05:05:   if (xxx) then
:05:05:     if (xxx) then
:06:06:       a = 2
:-1:05:     endif
:-1:06:   endif
:06:06:   a = 1
:05:05:   if (xxx) then
:05:05:     if (xxx) then
:05:05:       if (xxx) then
:-1:05:       endif
:-1:05:     endif
:-1:06:   endif
:06:06:   a = 1
:05:05:   if (xxx) then
:05:05:     if (xxx) then
:05:05:       if (xxx) then
:06:06:         a = 2
:-1:05:       endif
:-1:05:     endif
:-1:06:   endif
:05:05:   if (xxx) then
:05:05:     if (xxx) then
:06:06:       a = 2
:-1:05:     endif
:05:05:     if (xxx) then
:06:06:       a = 3
:-1:05:     endif
:-1:06:   endif
:05:05:   do i = 1,n
:05:05:     do j = 1, m
:05:05:       if (a.eq.b) then
:06:06:         a = 1
:05:05:       else
:06:06:         b = 1
:-1:05:       endif
:06:06:       call mysub()
:06:06:       a = b
:-1:05:     enddo
:-1:06:   enddo
:06:06:   return
!
:-1:01:   end subroutine my_sub1
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function my_fun(a, b) result(d)
!-------------------------------------------------------------------------------
:03:03:   use mymod1, only: dj, fj
:03:03:   use mymod2
:03:03:   implicit none
:03:03:   include 'jfkdla.inc'
!
:04:04:   integer(i4), dimension(:), intent(inout) :: a, b
:-3:04:!local
:04:04:   integer(i4) :: i, n
!
:06:06:   n = size(a)
:06:06:   d = 0
:05:05:   do i = 1,n
:06:06:     d = a(i)+b(i)
:-1:06:   enddo
!
:-1:01:   end function
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   subroutine testtest
!
:05:05:   if (xxx) then
:05:05:   elseif (xxx) then
:-1:02:   endif
!
:05:05:   if (xxx) then
:06:06:     a = 1
:05:05:   elseif (xxx) then
:06:06:     a = 2
:-1:02:   endif
!
:05:05:   if (xxx) then
:06:06:     a = 1
:06:06:     b = 2
:05:05:   elseif (xxx) then
:06:06:     a = 2
:06:06:     b = 3
:-1:02:   endif
!
:05:05:   if (xxx) then
:05:05:   elseif (xxx) then
:05:05:   elseif (xxx) then
:-1:02:   endif
!
:05:05:   if (xxx) then
:06:06:     a = 1
:05:05:   elseif (xxx) then
:06:06:     a = 2
:05:05:   elseif (xxx) then
:06:06:     a = 3
:-1:02:   endif
!
:05:05:   if (xxx) then
:06:06:     a = 1
:06:06:     b = 2
:05:05:   elseif (xxx) then
:06:06:     a = 2
:06:06:     b = 3
:05:05:   elseif (xxx) then
:06:06:     a = 3
:06:06:     b = 4
:-1:02:   endif
!
:05:05:   if (xxx) then
:05:05:     if (xxx) then
:05:05:     else
:-1:05:     endif
:05:05:     if (xxx) then
:05:05:     else
:-1:05:     endif
:05:05:   elseif (xxx) then
:05:05:     if (xxx) then
:05:05:     else
:-1:05:     endif
:05:05:     if (xxx) then
:05:05:     else
:-1:05:     endif
:-1:02:   endif
!
:05:05:   if (xxx) then
:05:05:     if (xxx) then
:06:06:       a = 1
:05:05:     else
:06:06:       a = 1
:-1:05:     endif
:05:05:     if (xxx) then
:06:06:       a = 1
:05:05:     else
:06:06:       a = 1
:-1:05:     endif
:05:05:   elseif (xxx) then
:05:05:     if (xxx) then
:06:06:       a = 1
:05:05:     else
:06:06:       a = 1
:-1:05:     endif
:05:05:     if (xxx) then
:06:06:       a = 1
:05:05:     else
:06:06:       a = 1
:-1:05:     endif
:-1:02:   endif
!
:-1:01:   end subroutine testtest
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   subroutine my_fun(a, b) result(d)
!-------------------------------------------------------------------------------
:03:03:   use mymod1, only: dj, fj
:03:03:   use mymod2
:03:03:   implicit none
:03:03:   include 'jfkdla.inc'
!
:04:04:   integer(i4), dimension(:), intent(inout) :: a, b
:-3:04:!if (... &
:-3:04:!   ...) then
:-3:04:! if(aa=bb)then
:-3:04:! else
:-3:04:!  endif
:-3:04:!endif
:-3:04:!print *, ' fda fda &
:-3:04:!  &  fdafafdasfda &
:-3:04:!  &  fdafafdasfda &
:-3:04:!  &  fdafafdasfda &
:-3:04:!  &  fdafafdasfda'
!
:-1:01:   end subroutine
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:-1:00:   end module
!-------------------------------------------------------------------------------
