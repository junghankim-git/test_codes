!-------------------------------------------------------------------------------
   module bilinear
!-------------------------------------------------------------------------------
!
!  abstract :  bi-linear interpolation module
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only : i4, l4, r4, r8
   use matrix_control, only : matrix_inverse
!
   private
!
   type bilinear_t
     real(r8), dimension(4) :: coef
   end type bilinear_t
!
   public :: bilinear_t, bilinear_set, bilinear_interpolation
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine bilinear_set(bi_t, psi, x, y)
!-------------------------------------------------------------------------------
   implicit none
!
   type(bilinear_t),                 intent(inout) :: bi_t
   real(r8), dimension(4),           intent(in   ) :: psi ! target
   real(r8), dimension(4), optional, intent(in   ) :: x, y ! grid
! local variables
   integer(i4) :: i, j
   real(r8), dimension(:,:), allocatable :: mat, invmat
!
   allocate(invmat(4, 4))
!
! arbitrary domain (x1,y1),(x2,y2),(x3,y3),(x4,y4)
   if (present(x).and.present(y)) then !
!
     allocate(mat(4,4))

     do j = 1,4
       mat(1,j) = 1.0_r8; mat(2,j) = x(j) ; mat(3,j) = y(j) ; mat(4,j) = x(j)*y(j)
     enddo
!
     call matrix_inverse(4,mat,invmat)
!
     deallocate(mat)
!
   ! reference domain [-1,1]x[-1,1]
   else
!
     invmat(1,1) = 0.25_r8; invmat(2,1) = 0.25_r8; invmat(3,1) = 0.25_r8; invmat(4,1) = 0.25_r8
     invmat(1,2) =-0.25_r8; invmat(2,2) = 0.25_r8; invmat(3,2) = 0.25_r8; invmat(4,2) =-0.25_r8
     invmat(1,3) =-0.25_r8; invmat(2,3) =-0.25_r8; invmat(3,3) = 0.25_r8; invmat(4,3) = 0.25_r8
     invmat(1,4) = 0.25_r8; invmat(2,4) =-0.25_r8; invmat(3,4) = 0.25_r8; invmat(4,4) =-0.25_r8
!     bi_t%coef(1) = 0.25_r8*( psi(1)+psi(2)+psi(3)+psi(4))
!     bi_t%coef(2) = 0.25_r8*(-psi(1)+psi(2)+psi(3)-psi(4))
!     bi_t%coef(3) = 0.25_r8*(-psi(1)-psi(2)+psi(3)+psi(4))
!     bi_t%coef(4) = 0.25_r8*( psi(1)-psi(2)+psi(3)-psi(4))
   endif
!
   do j = 1,4
     bi_t%coef(j) = 0.0_r8
     do i = 1,4
       bi_t%coef(j) = bi_t%coef(j)+invmat(i,j)*psi(i)
     enddo
   enddo
!
   deallocate(invmat)
!
   return
   end subroutine bilinear_set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine bilinear_interpolation(bi_t, psi, x, y)
!-------------------------------------------------------------------------------
   implicit none
!
   type(bilinear_t), intent(in   ) :: bi_t
   real(r8),         intent(  out) :: psi
   real(r8),         intent(in   ) :: x, y
! local variables
   integer(i4) :: i, j
!
   psi = bi_t%coef(1)+bi_t%coef(2)*x+bi_t%coef(3)*y+bi_t%coef(4)*x*y
!
   return
   end subroutine bilinear_interpolation
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module bilinear
