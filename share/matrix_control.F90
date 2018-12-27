!-------------------------------------------------------------------------------
   module matrix_control
!-------------------------------------------------------------------------------
!
!  abstract : matrix control module 
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds
   use statistics, only: get_lerror
!
   implicit none
!
   public :: print_matrix, matrix_inverse, matrix_multiplicity, matrix_chk_identity
   public :: get_diagonal, tridiagonal_inverse, fast_tridiagonal_inverse, fast_pentadiagonal_inverse
   public :: diagonal_inverse_interface, diagonal_solve_interface
   public :: fast_pentadiagonal_dummy
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_matrix(string, m, n, mat)
!-------------------------------------------------------------------------------
   implicit none

   character(*),             intent(in   ) :: string
   real(r8), dimension(m,n), intent(in   ) :: mat
   integer(r4),              intent(in   ) :: m, n
! local variables
   integer(4)    :: i, j
   character(20) :: frmt
!
   write(frmt,'(a,i4.4,a)') '(',m,'(f6.2,x))'
!
   print*, string
!   do i = 1, m
!     write(*,trim(frmt)) (mat(i,j),j=1,n)
!   enddo
   do j = 1, n
     write(*,trim(frmt)) mat(:,j)
   enddo
   print*, ' '
!
   return
   end subroutine print_matrix
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine matrix_inverse(n, a, inva)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n,n), intent(in   ) :: a
   real(r8), dimension(n,n), intent(inout) :: inva
! local variables
   integer(i4) :: i, j
   real(r8) :: det, inv_det
! - for lapack
   real(r8), dimension(n) :: ipivot
   real(r8), dimension(n) :: work
   integer(i4) :: info
!
   inva(:,:) = a(:,:)
!
   call dgetrf(n,n,inva,n,ipivot,info)
   if (info.ne.0) then
     print*,'lapack error : dgetrf'
     stop
   endif
   call dgetri(n,inva,n,ipivot,work,n,info)
   if (info.ne.0) then
     print*, 'lapack error : dgetri'
     stop
   endif
!
   return
   end subroutine matrix_inverse
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function matrix_multiplicity(n, a, b) result(c)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n,n), intent(in   ) :: a, b
   real(r8), dimension(n,n)                :: c
! local variables
   integer(i4) :: i, j, k
!
   c(:,:) = 0.0
   do j = 1,n
   do i = 1,n
     do k = 1,n
       c(i,j) = c(i,j)+a(i,k)*b(k,j)
     enddo
   enddo
   enddo
!
   end function matrix_multiplicity
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function matrix_chk_identity(n, a, b) result(lerr)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n,n), intent(in   ) :: a, b
   real(r8)                                :: lerr
! local variables
   integer(i4) :: i, j, k
   real(r8), dimension(n,n) :: ii, jj
!
   ii(:,:) = 0.0_r8
   jj(:,:) = 0.0_r8
   do i = 1,n
     jj(i,i) = 1.0_r8
   enddo
   ii = matrix_multiplicity(n,a,b)
!
   lerr = get_lerror(n,n,ii,jj,'L2')
!
   end function matrix_chk_identity
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function matrix_chk_identity_old(n, a, b) result(isid)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n,n), intent(in   ) :: a, b
   logical(l4) :: isid
! local variables
   integer(i4) :: i, j, k
   real(r8), dimension(n, n) :: ii
!
   do j = 1,n
   do i = 1,n
     ii(i,j) = 0.0_r8
     do k = 1,n
       ii(i,j) = ii(i,j)+a(i,k)*b(k,j)
     enddo
   enddo
   enddo
!
   isid = .true.
   do j = 1,n
   do i = 1,n
     if (i.eq.j) then
       if (abs(ii(i,j)-1.0_r8).gt.1.0d-15) then
         isid = .false.
       endif
     else
       if (abs(ii(i,j)).gt.1.0d-15) then
         isid = .false.
       endif
     endif
   enddo
   enddo
!
   end function matrix_chk_identity_old
!-------------------------------------------------------------------------------
!
!   --> j
! | a1, c1, e1, 00, 00
! v b1, a2, c2, e2, 00
! i d1, b2, a3, c3, e3
!   00, d2, b3, a4, c4
!   00, 00, d3, b4, a5
!
!-------------------------------------------------------------------------------
   subroutine get_diagonal(n, a, b, c, mat, d, e)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                      intent(in   ) :: n
   real(r8), dimension(n),           intent(in   ) :: a, b, c
   real(r8), dimension(n,n),         intent(inout) :: mat
   real(r8), dimension(n), optional, intent(in   ) :: d, e
! local variables
   integer(i4) :: i, j
   logical(l4) :: penta
!
   if (present(d).and.present(e)) then
     penta = .true.
   else
     penta = .false.
   endif
!
   mat = 0.0_r8
   if (.not.penta) then
     do i=1,n
       if     (i.eq.1) then
         mat(  i,i) = a(i)
         mat(i+1,i) = b(i)
       elseif (i.eq.n) then
         mat(i-1,i) = c(i-1)
         mat(  i,i) = a(i)
       else
         mat(i-1,i) = c(i-1)
         mat(  i,i) = a(i)
         mat(i+1,i) = b(i)
       endif
     enddo
   else
     do i = 1,n
       if     (i.eq.1) then
         mat(  i,i) = a(i)
         mat(i+1,i) = b(i)
         mat(i+2,i) = d(i)
       elseif (i.eq.2) then
         mat(i-1,i) = c(i-1)
         mat(i  ,i) = a(i)
         mat(i+1,i) = b(i)
         mat(i+2,i) = d(i)
       elseif (i.eq.n) then
         mat(i-2,i) = e(i-2)
         mat(i-1,i) = c(i-1)
         mat(  i,i) = a(i)
       elseif (i.eq.n-1) then
         mat(i-2,i) = e(i-2)
         mat(i-1,i) = c(i-1)
         mat(  i,i) = a(i)
         mat(i+1,i) = b(i)
       else
         mat(i-2,i) = e(i-2)
         mat(i-1,i) = c(i-1)
         mat(  i,i) = a(i)
         mat(i+1,i) = b(i)
         mat(i+2,i) = d(i)
       endif
     enddo
   endif
!
   return
   end subroutine get_diagonal
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine diagonal_inverse_interface(im, ndia, n, a, b, c, imat, mat, d, e)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: im  ! 0(my), 1(fast), 2(LAPACK:general), 3(LAPACK:tri)
   integer(i4),              intent(in   ) :: ndia
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(inout) :: a, b, c
   real(r8), dimension(n,n), intent(inout) :: imat 
   real(r8), dimension(n),   optional, intent(inout) :: d, e
   real(r8), dimension(n,n), optional, intent(in   ) :: mat
! local variables
   integer(i4) :: i, err
   integer(i4), dimension(:)  , allocatable :: ipivot
   real(r8)   , dimension(:)  , allocatable :: work, a1, b1, c1, c2
   real(r8)   , dimension(:,:), allocatable :: ab
!
   if (im==2.and..not.present(mat)) then
     print*,'check option: mat'
     stop
   endif
   if (ndia/=3.and.ndia/=5) then
     print*,'check option: ndia'
     stop
   endif
   if (ndia==5.and.(.not.present(d).or..not.present(e))) then
     print*,'check option: ndia and d, e'
     stop
   endif
!
   if     (im==0) then
     !
     if (ndia==3) then
       call tridiagonal_inverse(n,a,b,c,imat)
     else
       print *, 'check'
       stop
     endif
     !
   elseif (im==1) then
     !
     if (ndia==3) then
       call fast_tridiagonal_inverse(n,a,b,c,imat)
     else
       call fast_pentadiagonal_inverse(n,a,b,c,d,e,imat)
     endif
     !
   elseif (im==2) then
     !
     allocate(ipivot(n),work(n))
     imat(:,:) = mat(:,:)
     call dgetrf(n,n,imat,n,ipivot,err)
#if 0
do i=1,n
print '(10(f6.2))',imat(i,:)
enddo
#endif
     call dgetri(n,imat,n,ipivot,work,n,err)
     deallocate(ipivot,work)
     !
   elseif (im==3) then
     !
     if (ndia==3) then
       !
       allocate(ipivot(n),work(n),a1(n),b1(n-1),c1(n-1),c2(1:n-2))
       a1 = a; b1 = b(1:n-1); c1 = c(1:n-1)
       call dgttrf(n,b1,a1,c1,c2,ipivot,err)
       imat = 0.0d0
       do i = 1,n
         if     (i==1) then
           imat(i  ,i) = a1(i)
           imat(i+1,i) = b1(i)
         elseif (i==2) then
           imat(i-1,i) = c1(i-1)
           imat(i  ,i) = a1(i)
           imat(i+1,i) = b1(i)
         elseif (i==n) then
           imat(i-2,i) = c2(i-2)
           imat(i-1,i) = c1(i-1)
           imat(i  ,i) = a1(i)
         else
           imat(i-2,i) = c2(i-2)
           imat(i-1,i) = c1(i-1)
           imat(i  ,i) = a1(i)
           imat(i+1,i) = b1(i)
         endif
       enddo
#if 0
do i=1,n
print '(10(f6.2))',imat(i,:)
enddo
#endif
       call dgetri(n,imat,n,ipivot,work,n,err)
       deallocate(ipivot,work,a1,b1,c1,c2)
       !
     else ! if ndia==5
       !
       allocate(ipivot(n),work(n),ab(4,n))
       ab(1,:) = 0.0_r8
       ab(2,:) = c(:)
       ab(3,:) = a(:)
       ab(4,:) = b(:)
       call dgbtrf(n,n,1,1,ab,n,ipivot,err)
       imat = 0.0d0
       do i = 1,n
         if     (i==1) then
           imat(i  ,i) = ab(3,i)
           imat(i+1,i) = ab(4,i)
         elseif (i==2) then
           imat(i-1,i) = ab(2,i-1)
           imat(i  ,i) = ab(3,i)
           imat(i+1,i) = ab(4,i)
         elseif (i==n) then
           imat(i-2,i) = ab(1,i-2)
           imat(i-1,i) = ab(2,i-1)
           imat(i  ,i) = ab(3,i)
         else
           imat(i-2,i) = ab(1,i-2)
           imat(i-1,i) = ab(2,i-1)
           imat(i  ,i) = ab(3,i)
           imat(i+1,i) = ab(4,i)
         endif
       enddo
#if 0
do i=1,n
print '(10(f6.2))',imat(i,:)
enddo
#endif
       call dgetri(n,imat,n,ipivot,work,n,err)
       deallocate(ipivot,work,ab)
       !
     endif ! if ndia == 5
     !
   else ! if im == 3
     print*,'check option'
     stop
   endif
!
   end subroutine diagonal_inverse_interface
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine diagonal_solve_interface(im, ndia, n, a, b, c, x, d, e)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: im  ! 0(my), 1(fast), 2(LAPACK:general), 3(LAPACK:tri)
   integer(i4),              intent(in   ) :: ndia
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(inout) :: a, b, c
   real(r8), dimension(n),   intent(inout) :: x
   real(r8), dimension(n),   optional, intent(inout) :: d, e
!   real(r8), dimension(n,n), optional, intent(in   ) :: mat
! local variables
   integer(i4) :: i, j, err, kl, ku, ldab
   integer(i4), dimension(:)  , allocatable :: ipivot
   real(r8)   , dimension(:)  , allocatable :: a1, b1, c1, y
   real(r8)   , dimension(:,:), allocatable :: mat, imat, ab
!
!   if (im==2.and..not.present(mat)) then
!     print*,'check option: mat'
!     stop
!   endif
   if (ndia/=3.and.ndia/=5) then
     print*,'check option: ndia'
     stop
   endif
   if (ndia==5.and.(.not.present(d).or..not.present(e))) then
     print*,'check option: ndia and d, e'
     stop
   endif
!
   if     (im==0) then
     !
     if (ndia==3) then
       allocate(imat(n,n))
       call tridiagonal_inverse(n,a,b,c,imat)
       deallocate(imat)
     else
       print *, 'check'
       stop
     endif
     !
   elseif (im==1) then
     !
     allocate(imat(n,n),y(n))
     if (ndia==3) then
       call fast_tridiagonal_inverse(n,a,b,c,imat)
     else
       call fast_pentadiagonal_inverse(n,a,b,c,d,e,imat)
     endif
     y(:) = 0.0_r8
     do j = 1,n
       do i = 1,n
         y(i) = y(i)+imat(i,j)*x(j)
       enddo
     enddo
     x(:) = y(:)
     deallocate(imat,y)
     !
   elseif (im==2) then
     !
     allocate(imat(n,n),ipivot(n))
     if (ndia==3) then
       call get_diagonal(n,a,b,c,imat)
     else
       call get_diagonal(n,a,b,c,imat,d,e)
     endif
     call dgesv(n,1,imat,n,ipivot,x,n,err)
     deallocate(imat,ipivot)
     !
   elseif (im==3) then
     !
     if (ndia==3) then
       !
       allocate(a1(n),b1(n-1),c1(n-1))
       a1 = a; b1 = b(1:n-1); c1 = c(1:n-1)
       call dgtsv(n,1,b1,a1,c1,x,n,err)
       deallocate(a1,b1,c1)
       !
     else ! if ndia==5
       !
       kl = 2; ku = 2
       ldab = 2*kl+ku+1
       allocate(ipivot(n),ab(ldab,n))
       ! ref: LAPACK manual
       do j=1,n
         do i=max(1,j-ku),min(n,j+kl)
           if     (i==j-2) then
             ab(kl+ku+1+i-j,j) = e(i)
           elseif (i==j-1) then
             ab(kl+ku+1+i-j,j) = c(i)
           elseif (i==j) then
             ab(kl+ku+1+i-j,j) = a(j)
           elseif (i==j+1) then
             ab(kl+ku+1+i-j,j) = b(j)
           elseif (i==j+2) then
             ab(kl+ku+1+i-j,j) = d(j)
           endif
         enddo
       enddo
       ! equivalent to
       !ab(1,:)     = 0.0_r8
       !ab(2,:)     = 0.0_r8
       !ab(3,3:n)   = e(1:n-2)
       !ab(4,2:n)   = c(1:n-1)
       !ab(5,:)     = a(:)
       !ab(6,1:n-1) = b(1:n-1)
       !ab(7,1:n-2) = d(1:n-2)
       !
       call dgbsv(n,kl,ku,1,ab,ldab,ipivot,x,n,err)
       !if (err/=0) then
       !  print *,'check dgbsv',err
       !  stop
       !endif
       deallocate(ipivot,ab)
       !
     endif ! if ndia == 5
     !
   else ! if im == 3
     print*,'check option'
     stop
   endif
!
   end subroutine diagonal_solve_interface
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine tridiagonal_inverse(n, a, b, c, inv)
!-------------------------------------------------------------------------------
   ! ref: Moawwad El-Mikkawy and Abdelrahman Karawia, 2005
   !      Inversion of general tridiagonal matrices
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(in   ) :: a, b, c
   real(r8), dimension(n,n), intent(inout) :: inv
! local variables
   integer(i4) :: i, j, k
   real(r8)    :: tmp, al(0:n), be(1:n+1)
!
   ! alpha
   al(0) = 1.0_r8
   al(1) = a(1)
   do i = 2,n
     al(i) = a(i)*al(i-1)-b(i)*c(i-1)*al(i-2)
   enddo
   ! beta
   be(n+1) = 1
   be(n  ) = a(n)
   do i = n-1,1,-1
     be(i) = a(i)*be(i+1)-b(i+1)*c(i)*be(i+2)
   enddo
   ! inverse matrix
   inv(1,1) = 1.0/( a(1)-(b(2)*c(1)*be(3)/be(2)) )
   inv(n,n) = 1.0/( a(n)-(b(n)*c(n-1)*al(n-2)/al(n-1)) )
   do i = 2,n-1
     inv(i,i) = 1.0/( a(i) -(b(i)*c(i-1)*al(i-2)/al(i-1)) -(b(i+1)*c(i)*be(i+2)/be(i+1)) )
   enddo
   do j = 1,n
   do i = 1,n
     if     (i.lt.j) then
       tmp = 1.0
       do k = 1,j-i
         tmp = tmp*c(j-k)
       enddo
       inv(i,j) = (-1.0)**dble(j-i)*tmp*(al(i-1)/al(j-1))*inv(j,j)
     elseif (i.gt.j) then
       tmp = 1.0
       do k = 1,i-j
         tmp = tmp*b(j+k)
       enddo
       inv(i,j) = (-1.0)**dble(i-j)*tmp*(be(i+1)/be(j+1))*inv(j,j)
     endif
   enddo
   enddo
   inv = transpose(inv)
!
   return
   end subroutine tridiagonal_inverse
!-------------------------------------------------------------------------------
!
! algorithem matrix
!   --> j
! | a1, c1, 00, 00
! v b1, a2, c2, 00
! i 00, b2, a3, c3
!   00, 00, b3, a4
!-------------------------------------------------------------------------------
   subroutine fast_tridiagonal_inverse(n, a, b, c, CC)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(in   ) :: a, b
   real(r8), dimension(n),   intent(inout) :: c
   real(r8), dimension(n,n), intent(inout) :: CC
! local variables
   integer(i4) :: i, j, k
   !real(r8)    :: CC(n,n), EE(n,n), AA(0:n)
   real(r8)    :: EE(n,n), AA(0:n)
!
   c(n) = 1.0_r8
! calculate CC(:,1:n)
   ! cal AA(0:n)
   AA(0) = 1.0_r8
   AA(1) = -a(1)*AA(0)/c(1)
   do i = 1,n-1
     !AA(i+1) = -(b(i+1)*AA(i-1)+a(i+1)*AA(i))/c(i+1)
     AA(i+1) = -(b(i)*AA(i-1)+a(i+1)*AA(i))/c(i+1)
   enddo
   do i = 1,n
     CC(i,n) = -AA(i-1)/AA(n)
   enddo
! calculate CC(:,1:n-1)
   EE(:,:) = 0.0_r8
   do i = 1,n
     EE(i,i) = 1.0_r8
   enddo
   CC(:,n-1) = 1.0_r8/c(n-1)*(EE(:,n)-a(n)*CC(:,n))
   do j = n-1,2,-1
     !CC(:,j-1) = 1.0/c(j-1)*(EE(:,j)-a(j)*CC(:,j)-b(j+1)*CC(:,j+1))
     CC(:,j-1) = 1.0_r8/c(j-1)*(EE(:,j)-a(j)*CC(:,j)-b(j)*CC(:,j+1))
   enddo
!
   return
   end subroutine fast_tridiagonal_inverse
!-------------------------------------------------------------------------------
!
! algorithem matrix
!   --> j
! | a1, c1, e1, 00, 00
! v b1, a2, c2, e2, 00
! i d1, b2, a3, c3, e3
!   00, d2, b3, a4, c4
!   00, 00, d3, b4, a5
!-------------------------------------------------------------------------------
   subroutine fast_pentadiagonal_inverse(n, a, b, c, d, e, CC)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(in   ) :: a, b, d
   real(r8), dimension(n),   intent(inout) :: c, e
   real(r8), dimension(n,n), intent(inout) :: CC
! local variables
   integer(i4) :: i, j, k
   !real(r8)    :: CC(n,n), EE(n,n), AA(0:n+1), BB(0:n+1), XX(0:n+1), YY(0:n+1)
   !real(r8)    :: EE(n,n), AA(0:n+1), BB(0:n+1), XX(0:n+1), YY(0:n+1)
   real(r8), dimension(:,:), allocatable :: EE
   real(r8), dimension(:),   allocatable :: AA, BB, XX, YY
!
   allocate(EE(n,n),AA(0:n+1), BB(0:n+1), XX(0:n+1), YY(0:n+1))
!
   e(n-1) = 1.0_r8
   e(n)   = 1.0_r8
   c(n)   = 0.0_r8
! calculate CC(:,1:n)
   ! cal AA(0:n)
   AA(0) = 0.0_r8
   AA(1) = 1.0_r8
   AA(2) = -( a(1)*AA(0)+c(1)*AA(1) )/e(1)
   !AA(3) = -( b(2)*AA(0)+a(2)*AA(1)+c(2)*AA(2) )/e(2)
   AA(3) = -( b(1)*AA(0)+a(2)*AA(1)+c(2)*AA(2) )/e(2)
   do i = 3,n
     !AA(i+1) = -(d(i)*AA(i-3)+b(i)*AA(i-2)+a(i)*AA(i-1)+c(i)*AA(i))/e(i)
     AA(i+1) = -(d(i-2)*AA(i-3)+b(i-1)*AA(i-2)+a(i)*AA(i-1)+c(i)*AA(i))/e(i)
   enddo
   BB(0) = 1.0_r8
   BB(1) = 0.0_r8
   BB(2) = -( a(1)*BB(0)+c(1)*BB(1) )/e(1)
   !BB(3) = -( b(2)*BB(0)+a(2)*BB(1)+c(2)*BB(2) )/e(2)
   BB(3) = -( b(1)*BB(0)+a(2)*BB(1)+c(2)*BB(2) )/e(2)
   do i = 3,n
     !BB(i+1) = -(d(i)*BB(i-3)+b(i)*BB(i-2)+a(i)*BB(i-1)+c(i)*BB(i))/e(i)
     BB(i+1) = -(d(i-2)*BB(i-3)+b(i-1)*BB(i-2)+a(i)*BB(i-1)+c(i)*BB(i))/e(i)
   enddo
   do i = 0,n+1
     XX(i) = AA(n+1)*BB(i)-(AA(i)*BB(n+1))
     YY(i) = AA(n)*BB(i)-(AA(i)*BB(n))
   enddo
!
   do i = 1,n
     CC(i,n)   = -YY(i-1)/YY(n+1)
     CC(i,n-1) = -XX(i-1)/XX(n)
   enddo
! calculate CC(:,1:n-1)
   EE(:,:) = 0.0_r8
   do i = 1,n
     EE(i,i) = 1.0_r8
   enddo
!
   CC(:,n-2) = 1.0_r8/e(n-2)*(EE(:,n)-a(n)*CC(:,n)-c(n-1)*CC(:,n-1))
   CC(:,n-3) = 1.0_r8/e(n-3)*(EE(:,n-1)-c(n-2)*CC(:,n-2)-a(n-1)*CC(:,n-1)-b(n-1)*CC(:,n))
   do j = n-2,3,-1
     !CC(:,j-2) = 1.0_r8/e(j-2)*(EE(:,j)-c(j-1)*CC(:,j-1)-a(j)*CC(:,j)-b(j+1)*CC(:,j+1)-d(j+2)*CC(:,j+2))
     CC(:,j-2) = 1.0_r8/e(j-2)*(EE(:,j)-c(j-1)*CC(:,j-1)-a(j)*CC(:,j)-b(j)*CC(:,j+1)-d(j)*CC(:,j+2))
   enddo
   deallocate(EE,AA,BB,XX,YY)
!
   return
   end subroutine fast_pentadiagonal_inverse
!-------------------------------------------------------------------------------
!
! algorithem matrix
!   --> j
! | a1, c1, e1, 00, 00
! v b1, a2, c2, e2, 00
! i d1, b2, a3, c3, e3
!   00, d2, b3, a4, c4
!   00, 00, d3, b4, a5
!-------------------------------------------------------------------------------
   subroutine fast_pentadiagonal_dummy(n, a, b, c, d, e, CC)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(in   ) :: a, b, d
   real(r8), dimension(n),   intent(inout) :: c, e
   real(r8), dimension(n,n), intent(inout) :: CC
! local variables
   integer(i4) :: i, j, k
   !real(r8)    :: CC(n,n), EE(n,n), AA(0:n+1), BB(0:n+1), XX(0:n+1), YY(0:n+1)
   !real(r8)    :: EE(n,n), AA(0:n+1), BB(0:n+1), XX(0:n+1), YY(0:n+1)
!   real(r8), dimension(:,:), allocatable :: EE
   real(r8), dimension(:),   allocatable :: AA, BB, XX, YY, EE
!
   allocate(EE(1:n),AA(0:n+1),BB(0:n+1),XX(0:n+1),YY(0:n+1))
!
   e(n-1) = 1.0_r8
   e(n)   = 1.0_r8
   c(n)   = 0.0_r8
! calculate CC(:,1:n)
   ! cal AA(0:n)
   AA(0) = 0.0_r8
   AA(1) = 1.0_r8
   AA(2) = -( a(1)*AA(0)+c(1)*AA(1) )/e(1)
   !AA(3) = -( b(2)*AA(0)+a(2)*AA(1)+c(2)*AA(2) )/e(2)
   AA(3) = -( b(1)*AA(0)+a(2)*AA(1)+c(2)*AA(2) )/e(2)
   do i = 3,n
     !AA(i+1) = -(d(i)*AA(i-3)+b(i)*AA(i-2)+a(i)*AA(i-1)+c(i)*AA(i))/e(i)
     AA(i+1) = -(d(i-2)*AA(i-3)+b(i-1)*AA(i-2)+a(i)*AA(i-1)+c(i)*AA(i))/e(i)
   enddo
   BB(0) = 1.0_r8
   BB(1) = 0.0_r8
   BB(2) = -( a(1)*BB(0)+c(1)*BB(1) )/e(1)
   !BB(3) = -( b(2)*BB(0)+a(2)*BB(1)+c(2)*BB(2) )/e(2)
   BB(3) = -( b(1)*BB(0)+a(2)*BB(1)+c(2)*BB(2) )/e(2)
   do i = 3,n
     !BB(i+1) = -(d(i)*BB(i-3)+b(i)*BB(i-2)+a(i)*BB(i-1)+c(i)*BB(i))/e(i)
     BB(i+1) = -(d(i-2)*BB(i-3)+b(i-1)*BB(i-2)+a(i)*BB(i-1)+c(i)*BB(i))/e(i)
   enddo
   do i = 0,n+1
     XX(i) = AA(n+1)*BB(i)-(AA(i)*BB(n+1))
     YY(i) = AA(n)*BB(i)-(AA(i)*BB(n))
   enddo
!
   do i = 1,n
     CC(i,n)   = -YY(i-1)/YY(n+1)
     CC(i,n-1) = -XX(i-1)/XX(n)
   enddo
! calculate CC(:,1:n-1)
   EE(:) = 0.0_r8
   do i = n,1,-1
     EE(i) = 1.0_r8
   enddo
!
   CC(:,n-2) = 1.0_r8/e(n-2)*(EE(n-1)-a(n)*CC(:,n)-c(n-1)*CC(:,n-1))
   CC(:,n-3) = 1.0_r8/e(n-3)*(EE(n)-c(n-2)*CC(:,n-2)-a(n-1)*CC(:,n-1)-b(n-1)*CC(:,n))
   do j = n-2,3,-1
     !CC(:,j-2) = 1.0_r8/e(j-2)*(EE(:,j)-c(j-1)*CC(:,j-1)-a(j)*CC(:,j)-b(j+1)*CC(:,j+1)-d(j+2)*CC(:,j+2))
     CC(:,j-2) = 1.0_r8/e(j-2)*(EE(n-j+1)-c(j-1)*CC(:,j-1)-a(j)*CC(:,j)-b(j)*CC(:,j+1)-d(j)*CC(:,j+2))
   enddo
   deallocate(EE,AA,BB,XX,YY)
!
   return
   end subroutine fast_pentadiagonal_dummy
#if 0
!-------------------------------------------------------------------------------
!
! algorithem matrix
!   --> j
! | a1, c1, e1, 00, 00
! v b2, a2, c2, e2, 00
! i d3, b3, a3, c3, e3
!   00, d4, b4, a4, c4
!   00, 00, d5, b5, a5
!-------------------------------------------------------------------------------
   subroutine fast_pentadiagonal_inverse(n, a, b, c, d, e, CC)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(in   ) :: a, b, d
   real(r8), dimension(n),   intent(inout) :: c, e
   real(r8), dimension(n,n), intent(inout) :: CC
! local variables
   integer(i4) :: i, j, k
   !real(r8)    :: CC(n,n), EE(n,n), AA(0:n+1), BB(0:n+1), XX(0:n+1), YY(0:n+1)
   real(r8)    :: EE(n,n), AA(0:n+1), BB(0:n+1), XX(0:n+1), YY(0:n+1)
!
   e(n-1) = 1.0
   e(n)   = 1.0
   c(n)   = 0.0
! calculate CC(:,1:n)
   ! cal AA(0:n)
   AA(0) = 0.0
   AA(1) = 1.0
   AA(2) = -( a(1)*AA(0)+c(1)*AA(1) )/e(1)
   AA(3) = -( b(2)*AA(0)+a(2)*AA(1)+c(2)*AA(2) )/e(2)
   do i = 3,n
     AA(i+1) = -(d(i)*AA(i-3)+b(i)*AA(i-2)+a(i)*AA(i-1)+c(i)*AA(i))/e(i)
   enddo
   BB(0) = 1.0
   BB(1) = 0.0
   BB(2) = -( a(1)*BB(0)+c(1)*BB(1) )/e(1)
   BB(3) = -( b(2)*BB(0)+a(2)*BB(1)+c(2)*BB(2) )/e(2)
   do i = 3,n
     BB(i+1) = -(d(i)*BB(i-3)+b(i)*BB(i-2)+a(i)*BB(i-1)+c(i)*BB(i))/e(i)
   enddo
   do i = 0,n+1
     XX(i) = AA(n+1)*BB(i)-(AA(i)*BB(n+1))
     YY(i) = AA(n)*BB(i)-(AA(i)*BB(n))
   enddo
!
   do i = 1,n
     CC(i,n)   = -YY(i-1)/YY(n+1)
     CC(i,n-1) = -XX(i-1)/XX(n)
   enddo
! calculate CC(:,1:n-1)
   EE(:,:) = 0.0
   do i = 1,n
     EE(i,i) = 1.0
   enddo
!
   CC(:,n-2) = 1.0/e(n-2)*(EE(:,n)-a(n)*CC(:,n)-c(n-1)*CC(:,n-1))
   CC(:,n-3) = 1.0/e(n-3)*(EE(:,n-1)-c(n-2)*CC(:,n-2)-a(n-1)*CC(:,n-1)-b(n)*CC(:,n))
   do j = n-2,3,-1
     CC(:,j-2) = 1.0/e(j-2)*(EE(:,j)-c(j-1)*CC(:,j-1)-a(j)*CC(:,j)-b(j+1)*CC(:,j+1)-d(j+2)*CC(:,j+2))
   enddo
!
   return
   end subroutine fast_pentadiagonal_inverse
#endif
#if 0
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function determinent(n, a) result(det)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n,n), intent(in   ) :: a
   real(r8) :: det
! local variables
   real(r8) :: plus, minus
   integer(i4) :: i, j
!
   end function determinent
#endif
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module matrix_control
