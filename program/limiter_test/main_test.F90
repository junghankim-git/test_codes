!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
   implicit none
!
   integer, parameter :: i4 = 4, r8 = 8
!
   integer(i4), parameter :: n = 2
   real(r8), dimension(n, n, 1) :: qq
!
   print*, 'Test 1'
   qq(1, 1, 1) = 1.0_r8
   qq(2, 1, 1) = 1.0_r8
   qq(1, 2, 1) = 1.0_r8
   qq(2, 2, 1) = 1.0_r8
   print*, qq
   call limiter2d_zero(qq, n, 1)
   print*, qq
   print*, ' '
!
   print*, 'Test 2'
   qq(1, 1, 1) =-1.0_r8
   qq(2, 1, 1) =-1.0_r8
   qq(1, 2, 1) =-1.0_r8
   qq(2, 2, 1) =-1.0_r8
   print*, qq
   call limiter2d_zero(qq, n, 1)
   print*, qq
   print*, ' '
!
   print*, 'Test 3'
   qq(1, 1, 1) = 1.0_r8
   qq(2, 1, 1) = 1.0_r8
   qq(1, 2, 1) = 1.0_r8
   qq(2, 2, 1) =-1.0_r8
   print*, qq
   call limiter2d_zero(qq, n, 1)
   print*, qq
   print*, ' '
!
   print*, 'Test 4'
   qq(1, 1, 1) = 1.0_r8
   qq(2, 1, 1) =-1.0_r8
   qq(1, 2, 1) =-1.0_r8
   qq(2, 2, 1) =-1.0_r8
   print*, qq
   call limiter2d_zero(qq, n, 1)
   print*, qq
   print*, ' '
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine limiter2d_zero(q, np, nlev)
!------------------------------------------
! mass conserving zero limiter (2D only).  to be called just before DSS
!
! this routine is called inside a DSS loop, and so Q had already
! been multiplied by the mass matrix.  Thus dont include the mass
! matrix when computing the mass = integral of Q over the element
!
! 20150717, SJ.CHOI
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),                       intent(in   ) :: np, nlev
   real(r8), dimension(np, np, nlev), intent(inout) :: q
! local variables
   integer :: i, j, k, kk
   real(r8) :: mass, mass_new, ml
!
   do k = 1, nlev
     mass = 0
     do j = 1, np
       do i = 1, np
         ml = q(i, j, k)
         mass = mass+ml
       enddo
     enddo
  
     ! negative mass.  so reduce all postive values to zero
     ! then increase negative values as much as possible
     if (mass<0) q(:,:, k) =-q(:,:, k)
     mass_new = 0
     do j = 1, np
     do i = 1, np
       if (q(i, j, k)<0) then
         q(i, j, k) = 0
       else
         ml = q(i, j, k)
         mass_new = mass_new+ml
       endif
     enddo
     enddo
  
     ! now scale the all positive values to restore mass
     if (mass_new>0) q(:,:, k) = q(:,:, k)*abs(mass)/mass_new
     if (mass<0) q(:,:, k) =-q(:,:, k)
   enddo

!
   return
   end subroutine limiter2d_zero
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program test
!-------------------------------------------------------------------------------
