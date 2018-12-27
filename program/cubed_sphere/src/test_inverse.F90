!-------------------------------------------------------------------------------
   program test
!
   use kinds,          only: i4, l4, r4, r8
   use matrix_control, only: matrix_inverse, matrix_chk_identity
!
   implicit none
!
   real(r8), parameter :: pi = acos(-1.0_r8)
   integer(i4), parameter :: n = 2
   real(r8), dimension(2,2) :: a, inva
   integer(i4) :: i
!
   a(1,1) = 2.3_r8
   a(2,1) = 1.3_r8
   a(1,2) = 6.3_r8
   a(2,2) =-4.3_r8
!
   call matrix_inverse(n,a,inva)
   print*,matrix_chk_identity(n,a,inva)
!
   end program test
!-------------------------------------------------------------------------------
