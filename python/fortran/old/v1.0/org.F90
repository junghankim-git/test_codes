


:01:01:module mymod

:-1:00:end module mymod
:-3:00:!>  - JH KIM  : Added to "ghost_exchangeV" subroutine. (CAM-SE 5.3)
:-3:00:!-------------------------------------------------------------------------------
:-3:00:#include "KIM.h"



:01:01:module testCOde
:03:03: use MyModules1, only : var11, var12, var13
:03:03: use MyModules2, only : var11, var12, var13
:03:03: use mymodules3, only : fjdkl, fjdkla, gja,:-2:03:fjdla,:-2:03:fdjalk, fjklda
:03:03: implicit none
:04:04: private
:05:05: type my_type
:04:04:   integer(i4) :: a, b
:04:04:   real(r8) :: fjdk, fjdkla
:-1:04: end type
:04:04: public :: ghost_exchangeV_new
:05:05: interface my_interface
:04:04:   module procedure my_sub1
:04:04:   module procedure my_sub2
:04:04:   module procedure my_sub3
:-1:04: end interface myinterface


:-1:01:contains



:02:02: subroutine my_sub()
:06:06:   a = 2
:-1:01: end subroutine



:02:02: subroutine my_sub0()

:-1:01: end subroutine



:02:02: subroutine my_sub1(a, b, c)
:03:03:  use mymodule3, only : j, fj, jfkd
:03:03:  implicit none
:05:05:  type my_sub_t
:04:04:    real :: abc
:-1:03:  end type my_sub_t
:04:04:  type(my_sub_t) :: mytt
:04:04:  real(r8), dimension(:, :), allocatable :: a
:04:04:  integer(i4) :: b, c

:06:06:    ii = 10
:05:05:    do i = 1, 100
:06:06:      jfdk = fjdk
:-1:06:    end do
:05:05:    do j = 1, 100
:05:05:      do i = 1, 100
:06:06:        b = abd
:-1:05:      end do
:-1:06:    end do
:05:05:    do k = 1, 100
:05:05:      do j = 1, 100
:05:05:        do i = 1, 100
:06:06:          b = abd
:-1:05:        end do
:-1:05:        end do
:-1:06:    end do
:06:06:    ii = 1
:06:06:    if(ii .eq. 1) some = some1
:05:05:    if(ii .ne. 1) then
:06:06:      a = b
:05:05:      do i = 1, n
:06:06:        a = b
:05:05:        do j = 1, m
:06:06:          f = 1
:-1:06:        enddo
:06:06:        a = b
:-1:06:      enddo
:06:06:      a = b
:-1:06:    endif
:05:05:    if(a == 1) then
:06:06:      b = 2
:05:05:    else if(a == 2) then
:06:06:      b = 3
:05:05:    else
:06:06:      b = 0
:-1:06:    endif
:06:06:    return
:-1:01: end subroutine my_sub1



:02:02: function my_fun(a, b) result(d)
:03:03:  use mymod1, only : dj, fj
:03:03:  use mymod2
:03:03:  implicit none
:03:03:  include 'jfkdla.inc'
:04:04:  integer(i4), dimension(:), intent(inout) :: a, b
:-3:04:!local
:04:04:  integer(i4) :: i, n

:06:06:    n = size(a)
:06:06:    d = 0
:05:05:    do i = 1, n
:06:06:      d = a(i) + b(i)
:-1:06:    end do
:-1:01: end function

:-1:00:end module
