!-------------------------------------------------------------------------------
   module pseudo_inverse
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only : i4, r4, r8
   use matrix_control, only : print_matrix
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), public :: m, n
   real(r8), dimension(:,:), allocatable, public :: a
   ! SVD
   real(r8), dimension(:,:), allocatable, public :: u
   real(r8), dimension(:,:), allocatable, public :: vt
   real(r8), dimension(:), allocatable, public :: s
!
   real(r8), dimension(:,:), allocatable, public :: sigma
   real(r8), dimension(:,:), allocatable, public :: sigmat
   ! result
   real(r8), dimension(:,:), allocatable, public :: ainv
!
   public :: psedo_inv_initialize, psedo_inv_finalize, psedo_inv_get
   public :: psedo_inv_print_svd, psedo_inv_check_svd
   public :: psedo_inv_print, psedo_inv_check
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine psedo_inv_initialize(m_in, n_in)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: m_in, n_in
!
   m = m_in
   n = n_in
   allocate(a(m_in,n_in))
   allocate(u(m_in,m_in))
   allocate(vt(n_in,n_in))
   allocate(s(min(n_in,m_in)))
   allocate(sigma(m_in,n_in))
   allocate(sigmat(n_in,m_in))
   allocate(ainv(n_in,m_in))
!
   return
   end subroutine psedo_inv_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine psedo_inv_finalize()
!-------------------------------------------------------------------------------
   implicit none
!
   deallocate(a)
   deallocate(u)
   deallocate(vt)
   deallocate(s)
   deallocate(sigma)
   deallocate(sigmat)
   deallocate(ainv)
!
   return
   end subroutine psedo_inv_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine psedo_inv_print_svd()
!-------------------------------------------------------------------------------
   implicit none
!
   call print_matrix('# check matrix(a) = ',a,m,n)
!
   call print_matrix('# singular values(sigma) = ',s,1,min(m,n))
!
   call print_matrix('# left singular vector(u) = ',u,m,m)
!
   call print_matrix('# right singular vector(vt) = ',vt,n,n)
!
   return
   end subroutine psedo_inv_print_svd
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine psedo_inv_check_svd()
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(m, n) :: res, tmp
   real(r8), dimension(m, m) :: idenm
   real(r8), dimension(n, n) :: idenn
! local variables
   integer(i4) :: i, j, k
!
   idenm = 0.0_r8
   do j = 1,m
   do i = 1,m
   do k = 1,m
     idenm(i,j) = idenm(i,j)+u(i,k)*u(j,k)
   enddo
   enddo
   enddo
!
   call print_matrix('# check identity(u*ut) = ',idenm,m,m)
!
   idenn = 0.0_r8
   do j = 1,n
   do i = 1,n
   do k = 1,n
     idenn(i,j) = idenn(i,j)+vt(i,k)*vt(j,k)
   enddo
   enddo
   enddo
!
   call print_matrix('# check identity(vt*v) = ',idenn,n,n)
!
   tmp = 0.0_r8
   do j = 1,n
   do i = 1,m
   do k = 1,m
     tmp(i,j) = tmp(i,j)+u(i,k)*sigma(k,j)
   enddo
   enddo
   enddo
!
   res = 0.0_r8
   do j = 1,n
   do i = 1,m
   do k = 1,n
     res(i,j) = res(i,j)+tmp(i,k)*vt(k,j)
   enddo
   enddo
   enddo
   call print_matrix('# check matrix(u*s*vt) = ',res,m,n)
   call print_matrix('# check matrix(sigma)  = ',sigma,m,n)
   call print_matrix('# check matrix(sigmat) = ',sigmat,n,m)
!
   return
   end subroutine psedo_inv_check_svd
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine psedo_inv_print()
!-------------------------------------------------------------------------------
   implicit none
!
   call print_matrix('# check matrix(a)    = ',a,m,n)
!
   call print_matrix('# check matrix(a^-1) = ',ainv,n,m)
!
   return
   end subroutine psedo_inv_print
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine psedo_inv_check()
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(m, m) :: iden
   integer(i4) :: i, j, k
!
   iden = 0.0_r8
   do j = 1,m
   do i = 1,m
   do k = 1,n
     iden(i,j) = iden(i,j)+a(i,k)*ainv(k,j)
   enddo
   enddo
   enddo
   call print_matrix('# check identity(a*a^-1) = ',iden,m,m)
!
   iden = 0.0_r8
   do j = 1,m
   do i = 1,m
   do k = 1,n
     iden(i,j) = iden(i,j)+ainv(i,k)*a(k,j)
   enddo
   enddo
   enddo
   call print_matrix('# check identity(a^-1*a) = ',iden,m,m)
!
   return
   end subroutine psedo_inv_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine psedo_inv_get(m_in, n_in, a_in, ainv_in)
   use lapack, only : gesvd
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                     intent(in   ) :: m_in, n_in
   real(r8), dimension(m_in, n_in), intent(inout) :: a_in
   real(r8), dimension(n_in, m_in), intent(inout) :: ainv_in
! local variables
   integer(i4) :: lwork, info
   real(r8), dimension(5*m_in*n_in) :: work
   real(r8), dimension(m_in, n_in) :: a_tmp
   real(r8), dimension(n_in, m_in) :: tmp
   integer(i4) :: i, j, k
!
   a     = a_in
   a_tmp = a_in
   lwork = 5*m_in*n_in
!
   call dgesvd('A','A',m_in,n_in,a_tmp,m_in,s,u,m_in,vt,n_in,work,lwork,info)
   if (info.ne.0) then
     print*,'DGESVD : check arguments.... program will be stoped.'
     stop
   endif
!
   sigma = 0.0_r8
   sigmat = 0.0_r8
   do j = 1,n
   do i = 1,m
     !if(i==j.and.s(i).ne.0.0_r8) then
     if (i==j) then
       if ((s(i).le.-1.0d-14).or.(s(i).ge.1.0d-14)) then
         sigma(i,j) = s(i)
         sigmat(j,i) = 1.0_r8/s(i)
       endif
     endif
   enddo
   enddo
!
   tmp = 0.0_r8
   do j = 1,m
   do i = 1,n
   do k = 1,n
     tmp(i,j) = tmp(i,j)+vt(k,i)*sigmat(k,j)
   enddo
   enddo
   enddo
!
   ainv = 0.0_r8
   do j = 1,m
   do i = 1,n
   do k = 1,m
     ainv(i,j) = ainv(i,j)+tmp(i,k)*u(j,k)
   enddo
   enddo
   enddo
!
   ainv_in = ainv
!
   return
   end subroutine psedo_inv_get
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module pseudo_inverse
!-------------------------------------------------------------------------------
