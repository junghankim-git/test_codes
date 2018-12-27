!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
   use kinds,          only: i4, r8
   use pseudo_inverse, only: psedo_inv_initialize, psedo_inv_finalize, psedo_inv_get
   use pseudo_inverse, only: psedo_inv_print_svd, psedo_inv_check_svd
   use pseudo_inverse, only: psedo_inv_print, psedo_inv_check
!-------------------------------------------------------------------------------
   implicit none
! 
   integer(i4), parameter :: m = 8, n = 8
   real(r8),dimension(m, n) :: a
   real(r8),dimension(n, m) :: ainv
!
   call psedo_inv_initialize(m,n)
   call mymatrix_initialize(1,m,n,a)
   call psedo_inv_get(m,n,a,ainv)
   call psedo_inv_print_svd()
   call psedo_inv_check_svd()
   call psedo_inv_print()
   call psedo_inv_check()
   call psedo_inv_finalize()
!
   end program test
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine mymatrix_initialize(opt, m_in, n_in, a_in)
   use kinds, only : i4, r8
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                     intent(in   ) :: opt, m_in, n_in
   real(r8), dimension(m_in, n_in), intent(inout) :: a_in
! local variables
   integer(i4) :: i, j
!
   if (opt==0.and.m_in==4.and.n_in==5) then
     a_in(1,1) = 1.0_r8; a_in(2,1) = 0.0_r8; a_in(3,1) = 0.0_r8; a_in(4,1) = 0.0_r8
     a_in(1,2) = 0.0_r8; a_in(2,2) = 0.0_r8; a_in(3,2) = 0.0_r8; a_in(4,2) = 4.0_r8
     a_in(1,3) = 0.0_r8; a_in(2,3) = 3.0_r8; a_in(3,3) = 0.0_r8; a_in(4,3) = 0.0_r8
     a_in(1,4) = 0.0_r8; a_in(2,4) = 0.0_r8; a_in(3,4) = 0.0_r8; a_in(4,4) = 0.0_r8
     a_in(1,5) = 2.0_r8; a_in(2,5) = 0.0_r8; a_in(3,5) = 0.0_r8; a_in(4,5) = 0.0_r8
   elseif (opt==1.and.m_in==8.and.n_in==8) then
     open(unit=21,file='in.txt',form='formatted',status='old')
     !
     do j = 1,n_in
     do i = 1,m_in
       read(21,*) a_in(i,j)
     enddo
     enddo
     !
     close(21)
   endif
!
   return
   end subroutine mymatrix_initialize
!-------------------------------------------------------------------------------
