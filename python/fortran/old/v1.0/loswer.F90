:00:00:module mymod
:-1:-1:end module mymod





:-3:-1:!>  - JH KIM  : Added to "ghost_exchangeV" subroutine. (CAM-SE 5.3)
:-3:-1:!-------------------------------------------------------------------------------
:-3:-1:#include "KIM.h"
:00:00:module testCOde
:02:02: use MyModules1, only : var11, var12, var13
:02:02: use MyModules2, only : var11, var12, var13
:02:02: use mymodules3, only : fjdkl, fjdkla, gja,:-2:02:fjdla,:-2:02:fdjalk, fjklda
:02:02: implicit none
:03:03: private
:04:04: type my_type
:03:03:  integer(i4) :: a, b
:03:03:  real(r8) :: fjdk, fjdkla
:-1:03: end type
:03:03: public :: ghost_exchangeV_new
:04:04: interface my_interface
:03:03:  module procedure my_sub1
:03:03:  module procedure my_sub2
:03:03:  module procedure my_sub3
:-1:03: end interface myinterface
:-1:-1:contains





:01:01: subroutine my_sub()
:05:05:   a = 2
:-1:-1:end subroutine





:01:01: subroutine my_sub0()
:-1:-1:end subroutine





:01:01: subroutine my_sub1(a, b, c)
:02:02:  use mymodule3, only : j, fj, jfkd
:02:02:  implicit none
:04:04:  type my_sub_t
:03:03:   real :: abc
:-1:03:  end type my_sub_t
:03:03:  type(my_sub_t) :: mytt
:03:03:  real(r8), dimension(:, :), allocatable :: a
:03:03:  integer(i4) :: b, c
:05:05:    ii = 10
:04:04:    do i = 1, 100
:05:05:      jfdk = fjdk
:-1:05:    end do
:04:04:    do j = 1, 100
:04:04:      do i = 1, 100
:05:05:        b = abd
:-1:04:      end do
:-1:05:    end do
:04:04:    do k = 1, 100
:04:04:      do j = 1, 100
:04:04:        do i = 1, 100
:05:05:          b = abd
:-1:04:        end do
:-1:04:        end do
:-1:05:    end do
:05:05:    ii = 1
:05:05:    if(ii .eq. 1) some = some1
:04:04:    if(ii .ne. 1) then
:05:05:      a = b
:04:04:      do i = 1, n
:05:05:        a = b
:04:04:        do j = 1, m
:05:05:          f = 1
:-1:05:        enddo
:05:05:        a = b
:-1:05:      enddo
:05:05:      a = b
:-1:05:    endif
:04:04:    if(a == 1) then
:05:05:      b = 2
:04:04:    else if(a == 2) then
:05:05:      b = 3
:04:04:    else
:05:05:      b = 0
:-1:05:    endif
:05:05:    return
:-1:-1:end subroutine my_sub1





:01:01: function my_fun(a, b) result(d)
:02:02:  use mymod1, only : dj, fj
:02:02:  use mymod2
:02:02:  implicit none
:02:02:  include 'jfkdla.inc'
:03:03:  integer(i4), dimension(:), intent(inout) :: a, b
:-3:03:!local
:03:03:  integer(i4) :: i, n
:05:05:    n = size(a)
:05:05:    d = 0
:04:04:    do i = 1, n
:05:05:      d = a(i) + b(i)
:-1:05:    end do
:-1:-1:end function





