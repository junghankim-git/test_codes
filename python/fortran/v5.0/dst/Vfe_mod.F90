#include "KIM.h"
!
!-----------------------------------------------------------------------------------------!
! VFE_MOD: main routine for VFE
!
!-----------------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
   module vfe_mod
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  your   name    code comment
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kiapsbase, only: iulog=>kim_iu_log, r8=>kim_real8_kind, i4=>kim_int_kind
   use dimensions, only: plev=>nlev, plevp=>nlevp
   use kiapsparallel, only: kim_par
   implicit none
!
   private
   logical, public :: debug_vfe = .false.
! Variables for the De Boor algorithm
   integer(i4), parameter, public :: basis_degree = 1 ! b-spline degree:1 = linear, 3 = cubic
   integer(i4), parameter, public :: deri_order = 1 ! first-order derivatives
   integer(i4), parameter, public :: nquad = 5 ! gaussian quadrature order
   integer(i4), parameter, public :: basis_order = basis_degree+1
   integer(i4), parameter, public :: nlevbc = 50+2 ! nmber of eta levels(= nlev+2 for linear)
   integer(i4), parameter, public :: nknot = nlevbc+2*basis_degree ! number of knot components
   integer(i4), parameter, public :: nbasis = nknot-basis_order ! number of basis functions = nknot-basis_order
   real(r8), public, dimension(:,:), allocatable, save :: vopderi
   real(r8), public, dimension(:,:), allocatable, save :: vopinte
   public :: ini_vfe, fin_vfe, sub_verint, sub_verder
   public :: vfe_mass_correction
   public :: minv, ludcmp, lubksb
!
   contains
!
!-----------------------------------------------------------------------------------------!
! SUBROUTINE INI_VFE
!
!-----------------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ini_vfe(n)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: n
!
   allocate(vopderi(n,n),vopinte(n+1,n))
!
   end subroutine ini_vfe
!
!-----------------------------------------------------------------------------------------!
! SUBROUTINE FIN_VFE
!
!-----------------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine fin_vfe
!-------------------------------------------------------------------------------
   implicit none
!
   deallocate(vopderi)
   deallocate(vopinte)
!
   end subroutine fin_vfe
!
!-----------------------------------------------------------------------------------------!
! SUBROUTINE VERIN
! Vertical Integration using matrix multiplication
!-----------------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine sub_verint(n, vopinte, zin, zout)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: n
   real(r8), intent(in   ) :: vopinte(n+1, n), zin(n)
   real(r8), intent(  out) :: zout(n+1)
   integer(i4) :: i, j
   real(r8) :: total
!
   zout(:) = 0.0_r8
   do j = 1,n
     do i = 1, n+1
       zout(i) = zout(i)+vopinte(i,j)*zin(j)
     enddo
   enddo
!
   end subroutine sub_verint
!
!-----------------------------------------------------------------------------------------!
! SUBROUTINE SUB_VERDER
! Vertical Derivative using matrix multiplication
!-----------------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine sub_verder(n, vopderi, zin, zout)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: n
   real(r8), intent(in   ) :: vopderi(n, n), zin(n)
   real(r8), intent(  out) :: zout(n)
   integer :: i, j
!
   zout(:) = 0.0_r8
   do j = 1,n
     do i = 1, n
       zout(i) = zout(i)+vopderi(i,j)*zin(j)
     enddo
   enddo
!
   end subroutine sub_verder
!
!-----------------------------------------------------------------------------------------!
! SUBROUTINE VFE_MASS_CORRECTION
! mass(pressure) correction in VFE scheme
!-----------------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine vfe_mass_correction(etai, hyai, hybi, hyam, hybm, dhyam, dhybm)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(plev+1), intent(in   ) :: etai, hyai, hybi
   real(r8), dimension(plev), intent(  out) :: hyam, hybm, dhyam, dhybm
   real(r8) :: zfac
   real(r8), dimension(plev) :: del_etaf, del_af, del_bf, dvaf, dvbf
   real(r8), dimension(plev+1) :: vafc, vbfc, vaf, vbf
   integer :: iref, iter, i, j, k
!
   if (kim_par%ismasterproc) write(*,*)
   if (kim_par%ismasterproc) write(*,*) '***mass correction for vfe***'
   vafc(:) = 0.0_r8
   vbfc(:) = 0.0_r8
   dvaf(:) = 0.0_r8
   dvbf(:) = 0.0_r8
   vaf(:) = 0.0_r8
   vbf(:) = 0.0_r8
! 1. Define dB/dn in the usual FD way at middle levels
   do k = 1,plev
     del_etaf(k) = etai(k+1)-etai(k) ! dn
     del_af(k) = hyai(k+1)-hyai(k) ! da
     del_bf(k) = hybi(k+1)-hybi(k) ! db
     dvbf(k) = del_bf(k)/del_etaf(k) ! db/dn
   enddo
! 2. Integrate dB/dn through the whole column,VBFC(plev+1) has total integral
   if (kim_par%ismasterproc) write(*,*) 'Vcoord coefficients b:'
   call sub_verint(plev,vopinte,dvbf,vbfc)
   if (kim_par%ismasterproc) write(*,'(a,es20.10) ') " before-correction:total_int b = ",vbfc(plev+1)
! 3. To obtain a corrected vector,B' and dB'/dn
   do k = 1,plev
     dvbf(k) = dvbf(k)/vbfc(plev+1) ! db'/dn
   enddo
   call sub_verint(plev,vopinte,dvbf,vbf) ! b'
   do k = 1,plev
     if (vbf(k)<0.0_r8) vbf(k) = 0.0_r8 ! to correct negative value
   enddo
   if (kim_par%ismasterproc) write(*,'(a,es20.10) ') " after-correction:total_int b = ",vbf(plev+1)
! 4. To obtain a corrected vector,A' and dA'/dn: Use an iterative process to correct the original vector
   if (kim_par%ismasterproc) write(*,*) 'Vcoord coefficients a:'
   iter = 0
   do
     iter = iter+1
     do k = 1,1
       dvaf(k) = del_af(k)/del_etaf(k)
       if (del_af(k)<0.0_r8) iref = k
     enddo
     do k = 2,plev
       dvaf(k) = del_af(k)/del_etaf(k)
       if (del_af(k)<0.0_r8.and.del_af(k-1)>= 0.0_r8) iref = k
     enddo
     call sub_verint(plev,vopinte,dvaf,vafc) ! a'
     zfac =-1.0_r8*vafc(iref)/(vafc(plev+1)-vafc(iref))
     do k = 1,iref
       del_af(k) = del_af(k)/zfac
       dvaf(k) = del_af(k)/del_etaf(k) ! da'/dn
     enddo
     if (dabs(vafc(plev+1))<1.d-14) then
       goto 111
     else
       if (kim_par%ismasterproc) write(*,'(a,2i5,es20.10) ') " before-correction:iref,iter,total_int a = ",iref,iter,vafc(plev+1)
     endif
     if (iter>20) then
       if (kim_par%ismasterproc) write(*,*) ' stop:the number of iteration exceeds iter = 20:hybridmassvfe.f90 '
       stop
     endif
   enddo
   111 call sub_verint(plev,vopinte,dvaf,vaf) ! a'
   do k = 1,plev
     if (vaf(k)<0.0_r8) vaf(k) = 0.0_r8 ! to correct negative value
   enddo
   if (kim_par%ismasterproc) write(*,'(a,2i5,es20.10) ') " after-correction:iref,iter,total_int a = ",iref,iter,vaf(plev+1)
   if (kim_par%ismasterproc) write(*,*)
   hyam(:) = vaf(1:plev)
   hybm(:) = vbf(1:plev)
   dhyam(:) = dvaf(:)
   dhybm(:) = dvbf(:)
!
   end subroutine vfe_mass_correction
!
!---------------------------------------------------------------------!
!     Subroutines for matrix inversion using LU Decomposition
!---------------------------------------------------------------------!
!===================================================================!
!  !---- Further Comments ------------------------------------------!
!
!     A*X= B
!
!   CALL LUDCMP(A,N,M,INDX,D)
!   CALL LUBKSB(A,N,M,INDX,B) ! The answer X is returned in B.
!                              ! Original Matrix A is destroyed.
!                              ! For further calculation with
!                              ! different B : CALL LUBKSB once more.
!
!---- Further Comments ----------------------------------------------!
!====================================================================!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine minv(a, n, m, d, y)
!--------------------------------------------------------------------!
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: n, m
   real(r8), intent(in   ) :: a(m, m)
   real(r8), intent(  out) :: y(m, m)
   real(r8), intent(  out) :: d
   real(r8) :: aa(m, m)
   integer :: indx(m), i, j
!--------------------------------------------------------------------!
!
   aa(:,:) = a(:,:) ! to save a matrix
   do i = 1,n
     do j = 1, n
       y(i,j) = 0.
     enddo
       y(i,i) = 1.d0
   enddo
   call ludcmp(aa,n,m,indx,d)
   do j = 1,n
     d = d*aa(j,j) ! determinant of the original matrix a
   enddo
   do j = 1,n
     call lubksb(aa,n,m,indx,y(1,j)) ! y = inverse matrix of a
   enddo ! a is destroyed.
!--------------------------------------------------------------------!
!
   end subroutine minv
!====================================================================!
!
!====================================================================!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ludcmp(a, n, m, indx, d)
!--------------------------------------------------------------------!
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: n, m
   real(r8), intent(inout) :: a(m, m)
   real(r8), intent(inout) :: d
   real(r8) :: aamax, sum, dum, tiny, vv(999)
   integer :: indx(m), i, j, k, imax
!
   data tiny/1.0d-20/
!--------------------------------------------------------------------!
!
   if (n.gt.999) stop ' too small work-array size ! '
   d = 1.d0
   do i = 1,n
     aamax = 0.
     do j = 1,n
       if (dabs(a(i,j)).gt.aamax) aamax = dabs(a(i,j))
     enddo
     if (aamax.eq.0) stop 'Singular matrix in ludcmp'
     vv(i) = 1.d0/aamax
   enddo
!
   do j = 1,n ! 19
     do i = 1, j-1
       sum = a(i,j)
       do k = 1,i-1
         sum = sum-a(i,k)*a(k,j)
       enddo
       a(i,j) = sum
     enddo
       aamax = 0.
       do i = j,n ! 16
         sum = a(i,j)
         do k = 1,j-1
           sum = sum-a(i,k)*a(k,j)
         enddo
         a(i,j) = sum
         dum = vv(i)*dabs(sum)
         if (dum.ge.aamax) then
           imax = i
           aamax = dum
         endif
       enddo ! 16
!
       if (j.ne.imax) then
         do k = 1, n
           dum = a(imax,k)
           a(imax,k) = a(j,k)
           a(j,k) = dum
         enddo
           d =-d
           vv(imax) = vv(j)
       endif
       indx(j) = imax
       if (a(j,j).eq.0) a(j,j) = tiny
       if (j.ne.n) then
         dum = 1.d0/a(j,j)
         do i = j+1,n
           a(i,j) = a(i,j)*dum
         enddo
       endif
   enddo !---19
!--------------------------------------------------------------------!
!
   end subroutine ludcmp
!====================================================================!
!
!
!====================================================================!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine lubksb(a, n, m, indx, b)
!--------------------------------------------------------------------!
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: n, m
   real(r8), intent(inout) :: a(m, m)
   real(r8), intent(inout) :: b(n)
   real(r8) :: sum, dum, tiny
   integer :: indx(m), i, j, k, ii, ll
!--------------------------------------------------------------------!
!
   ii = 0
   do i = 1,n
     ll = indx(i)
     sum = b(ll)
     b(ll) = b(i)
     if (ii.ne.0) then
       do j = ii, i-1
         sum = sum-a(i,j)*b(j)
       enddo
         elseif (sum.ne.0) then
         ii = i
     endif
     b(i) = sum
   enddo
   do i = n,1,-1
     sum = b(i)
     do j = i+1,n
       sum = sum-a(i,j)*b(j)
     enddo
     b(i) = sum/a(i,i)
   enddo
!--------------------------------------------------------------------!
!
   end subroutine lubksb
!====================================================================!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module vfe_mod
!-------------------------------------------------------------------------------
