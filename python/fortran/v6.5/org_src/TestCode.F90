module mymod
end module mymod

!>  - JH KIM  : Added to "ghost_exchangeV" subroutine. (CAM-SE 5.3)
!-------------------------------------------------------------------------------

#include "KIM.h"

module testCOde
    USE MyModules1, only    : var11, var12, var13
    use MyModules2, only : var11, var12, var13
USE mymodules3, only: fjdkl, fjdkla, gja, &
fjdla, &
                fdjalk, fjklda
  implicit none
  private

  type my_type
 integer (kind=i4) :: a, b
 integer (kind = i4) :: c, d
 integer (kind= i4) :: e,f
 integer (kind =i4) :: g,h
  real(r8) :: fjdk,fjdkla
end type

  public :: ghost_exchangeV_new


   interface my_interface
    module procedure my_sub1
    module procedure my_sub2
    module procedure my_sub3
 end interface myinterface


contains 

subroutine my_sub()
 a = 2
end subroutine

subroutine my_sub0()
end subroutine


subroutine my_sub1(a, b, c)
  use mymodule3, only: j, fj, jfkd
  implicit None
type my_sub_t
  real :: abc
end type my_sub_t
  type(my_sub_t) :: mytt
 real(r8), dimension(:, :), allocatable ::a
 integer(i4) :: b, c


a=1
if(xxx)then
endif
a=1
if(xxx)then
a=2
endif
a=1
if(xxx)then
if(xxx)then
endif
endif
a=1
if(xxx)then
if(xxx)then
a=2
endif
endif
a=1
if(xxx)then
if(xxx)then
if(xxx)then
endif
endif
endif
a=1
if(xxx)then
if(xxx)then
if(xxx)then
a=2
endif
endif
endif


 do i = 1,n
  do j = 1,m
    if (a.eq.b) then
      a = 1
    else
      b = 1
    endif
    call mysub()
    a = b
  enddo
 enddo

 return
end subroutine my_sub1





function my_fun(a, b) result(d)
 use mymod1, only: dj, fj
use mymod2
implicit none
 include 'jfkdla.inc'
integer(i4), dimension(:), INTENT(INOUT) :: a, b
!local
integer(i4) :: i, n

n = size(a)
d = 0
do i = 1, n
d = a(i) + b(i)
end do
end function

subroutine testtest


if(xxx) then
elseif (xxx) then
endif

if(xxx) then
a=1
elseif (xxx) then
a=2
endif

if(xxx) then
a=1
b=2
elseif (xxx) then
a=2
b=3
endif

if(xxx) then
elseif (xxx) then
elseif (xxx) then
endif

if(xxx) then
a=1
elseif (xxx) then
a=2
elseif (xxx) then
a=3
endif

if(xxx) then
a=1
b=2
elseif (xxx) then
a=2
b=3
elseif (xxx) then
a=3
b=4
endif

if(xxx) then
  if (xxx) then
  else
  endif
  if (xxx) then
  else
  endif
elseif (xxx) then
  if (xxx) then
  else
  endif
  if (xxx) then
  else
  endif
endif

if(xxx) then
  if (xxx) then
   a= 1
  else
   a= 1
  endif
  if (xxx) then
   a= 1
  else
   a= 1
  endif
elseif (xxx) then
  if (xxx) then
   a= 1
  else
   a= 1
  endif
  if (xxx) then
   a= 1
  else
   a= 1
  endif
endif

end subroutine testtest

   end module 
