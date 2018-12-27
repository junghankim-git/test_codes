


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
:06:06: type my_type
:04:04:   integer(i4) :: a, b
:04:04:   real(r8) :: fjdk, fjdkla
:-1:04: end type
:04:04: public :: ghost_exchangeV_new
:06:06: interface my_interface
:04:04:   module procedure my_sub1
:04:04:   module procedure my_sub2
:04:04:   module procedure my_sub3
:-1:04: end interface myinterface


:-1:01:contains



:02:02: subroutine my_sub()
:05:05:   a = 2
:-1:01: end subroutine



:02:02: subroutine my_sub0()

:-1:01: end subroutine



:02:02: subroutine my_sub1(a, b, c)
:03:03:  use mymodule3, only : j, fj, jfkd
:03:03:  implicit none
:06:06:  type my_sub_t
:04:04:    real :: abc
:-1:03:  end type my_sub_t
:04:04:  type(my_sub_t) :: mytt
:04:04:  real(r8), dimension(:, :), allocatable :: a
:04:04:  integer(i4) :: b, c

:05:05:    ii = 10
:06:06:    do i = 1, 100
:05:05:      jfdk = fjdk
:-1:05:    end do
:06:06:    do j = 1, 100
:06:06:      do i = 1, 100
:05:05:        b = abd
:-1:06:      end do
:-1:05:    end do
:06:06:    do k = 1, 100
:06:06:      do j = 1, 100
:06:06:        do i = 1, 100
:05:05:          b = abd
:-1:06:        end do
:-1:06:        end do
:-1:05:    end do
:05:05:    ii = 1
:05:05:    if(ii .eq. 1) some = some1
:06:06:    if(ii .ne. 1) then
:05:05:      a = b
:06:06:      do i = 1, n
:05:05:        a = b
:06:06:        do j = 1, m
:05:05:          f = 1
:-1:05:        enddo
:05:05:        a = b
:-1:05:      enddo
:05:05:      a = b
:-1:05:    endif
:06:06:    if(a == 1) then
:05:05:      b = 2
:06:06:    else if(a == 2) then
:05:05:      b = 3
:06:06:    else
:05:05:      b = 0
:-1:05:    endif
:05:05:    return
:-1:01: end subroutine my_sub1



:02:02: function my_fun(a, b) result(d)
:03:03:  use mymod1, only : dj, fj
:03:03:  use mymod2
:03:03:  implicit none
:03:03:  include 'jfkdla.inc'
:04:04:  integer(i4), dimension(:), intent(inout) :: a, b
:-3:04:!local
:04:04:  integer(i4) :: i, n

:05:05:    n = size(a)
:05:05:    d = 0
:06:06:    do i = 1, n
:05:05:      d = a(i) + b(i)
:-1:05:    end do
:-1:01: end function

:-1:00:end module
