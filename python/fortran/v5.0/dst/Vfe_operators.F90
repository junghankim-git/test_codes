!-------------------------------------------------------------------------------
   module vfe_operators
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
   use vfe_mod, only: debug_vfe, basis_degree, deri_order, nquad, basis_order, nlevbc, nknot, nbasis, minv
   use vfe_basis_funs, only: compute_basis_funs_eta, compute_basis_funs_quad, define_quadrature_set
   use kiapsparallel, only: kim_par
   implicit none
!
   public :: sub_deri_operator, sub_inte_operator
!
   contains
!
!---------------------------------------------------------------------------------!
! SUBROUTINE SUB_DERI_OPERATOR: Vertical derivative operator (VOPDERI)
!---------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine sub_deri_operator(debug_vfe, etaf, knots, vopderi)
!-------------------------------------------------------------------------------
   implicit none
! Variables for hybrid coordinates
!
   logical, intent(in   ) :: debug_vfe
   real(r8), intent(in   ), dimension(0:plev+1) :: etaf
   real(r8), intent(in   ), dimension(nknot) :: knots
   real(r8), intent(  out), dimension(plev, plev) :: vopderi
! Variables for b-spline basis function
   real(r8), dimension(nlevbc, nbasis) :: bspline_funs
   real(r8), dimension(nquad*(plev+1), nbasis) :: ibspline_funs_quad
   real(r8), dimension(nquad*(plev+1), nbasis) :: kbspline_funs_quad
   real(r8), dimension(nquad*(plev+1), nbasis) :: jbspline_deri_quad
! Variables for VFE
   real(r8) :: zdet
   real(r8), dimension(nlevbc, nbasis) :: zp
   real(r8), dimension(nlevbc, nbasis) :: zpinv
   real(r8), dimension(nlevbc, nbasis) :: zs
   real(r8), dimension(nlevbc, nbasis) :: za
   real(r8), dimension(nlevbc, nbasis) :: zainv
   real(r8), dimension(nlevbc, nbasis) :: zb
   real(r8), dimension(nlevbc, nbasis) :: zderi
   real(r8), dimension(nlevbc, nbasis) :: ztemp
   integer(i4) :: i, j, k, ilev, icount, iquad
   real(r8) :: deta, zw1, zw2
   real(r8), dimension(nquad) :: quadx, quadw
!
!--------------------------------------------------------------------------------------------------------!
! STEP 1: Matrix P^-1
! ZP(i,j) = BASISFUNC_j(eta_i)) from FE space of the function to be differenciated into physical space
! i: row->eta points, j: column->basis function
!--------------------------------------------------------------------------------------------------------!
!
!
   call compute_basis_funs_eta(etaf,knots,bspline_funs)
   zp(:,:) = bspline_funs(:,:)
! Invert SF: ZSFINV projects the function to be integrated from physical space onto FE space
   zdet = 0.0_r8
   call minv(zp,nlevbc,nbasis,zdet,zpinv)
!
!--------------------------------------------------------------------------------------------------------!
! STEP 2: Matrix S
! ZS(i,j)=BASISDER_j(eta_i) to project the derivative from FE space to physical space
!--------------------------------------------------------------------------------------------------------!
!
   zs(:,:) = zp(:,:)
   if (debug_vfe) then
     write(*,'(a) ') ''
     write(*,'(a) ') ' matrix s'
     write(*,'(a) ') '----------------------------------------------------------------'
     do i = 1,nlevbc
       write(*,'(50f6.1) ')(100.*zs(i,j),j = 1,nbasis)
     enddo
   endif
!
! Compute basis functions E_i(eta),E_k(eta),and dE_j(eta)/deta at eta quadrature points (quadx)
!
   call define_quadrature_set(nquad,quadw,quadx)
   call compute_basis_funs_quad(etaf,knots,ibspline_funs_quad,jbspline_deri_quad)
   kbspline_funs_quad(:,:) = ibspline_funs_quad(:,:)
!
!--------------------------------------------------------------------------------------------------------!
! STEP 3: Mass matrix A
! Compute A_ik(eta): integrate the product of basis functions at all quadrature points
!--------------------------------------------------------------------------------------------------------!
!
   do i = 1,nbasis
     do k = 1, nbasis
       icount = 1
       za(i,k) = 0.0_r8
       do ilev = 1,nlevbc-1
         deta = etaf(ilev)-etaf(ilev-1)
         do iquad = 1,nquad
           za(i,k) = za(i,k)+0.50_r8*deta*quadw(iquad)*ibspline_funs_quad(icount,i)*kbspline_funs_quad(icount,k)
           icount = icount+1
         enddo
       enddo
     enddo
   enddo
!
! Invert mass matrix: A^-1
   zdet = 0.0_r8
   call minv(za,nbasis,nbasis,zdet,zainv)
   if (debug_vfe) then
     write(*,'(a) ') ''
     write(*,'(a) ') ' mass matrix:a'
     write(*,'(a) ') '----------------------------------------------------------------'
     do i = 1,nbasis
       write(*,'(50f6.1) ')(100.*za(i,j),j = 1,nbasis)
     enddo
   endif
!
!--------------------------------------------------------------------------------------------------------!
! STEP 4: Operator matrix B
! Compute B_kj(eta): integrate the product of basis function and its derivative at all quadrature points
!--------------------------------------------------------------------------------------------------------!
!
   do j = 1,nbasis
     do k = 1, nbasis
       icount = 1
       zb(k,j) = 0.0_r8
       do ilev = 1,nlevbc-1
         deta = etaf(ilev)-etaf(ilev-1)
         do iquad = 1,nquad
           zb(k,j) = zb(k,j)+0.50_r8*deta*quadw(iquad)*jbspline_deri_quad(icount,j)*kbspline_funs_quad(icount,k)
           icount = icount+1
         enddo
       enddo
     enddo
   enddo
   if (debug_vfe) then
     write(*,'(a) ') ''
     write(*,'(a) ') ' operator matrix:b'
     write(*,'(a) ') '----------------------------------------------------------------'
     do k = 1,nbasis
       write(*,'(50f6.1) ')(100.*zb(k,j),j = 1,nbasis)
     enddo
   endif
!
!--------------------------------------------------------------------------------------------------------!
! STEP 5: ZDERI
! Compute the product of S*A^-1*B*P^-1
!--------------------------------------------------------------------------------------------------------!
!
! B*P^-1
   zderi(:,:) = 0.0_r8
   do j = 1,nbasis
     do i = 1, nbasis
       do k = 1, nbasis
         zderi(i,j) = zderi(i,j)+zb(i,k)*zpinv(k,j)
       enddo
       enddo
   enddo
!
! A^-1*B*P^-1
   ztemp(:,:) = 0.0_r8
   do j = 1,nbasis
     do i = 1, nbasis
       do k = 1, nbasis
         ztemp(i,j) = ztemp(i,j)+zainv(i,k)*zderi(k,j)
       enddo
       enddo
   enddo
!
! ZDERI=S*A^-1*B*P^-1
   zderi(:,:) = 0.0_r8
   do j = 1,nbasis
     do i = 1, nbasis
       do k = 1, nbasis
         zderi(i,j) = zderi(i,j)+zs(i,k)*ztemp(k,j)
       enddo
       enddo
   enddo
   if (debug_vfe) then
     write(*,'(a) ') ''
     write(*,'(a) ') ' product matrix:zderi'
     write(*,'(a) ') '----------------------------------------------------------------'
     do k = 1,nbasis
       write(*,'(50f6.1) ')(100.*zderi(k,j),j = 1,nbasis)
     enddo
   endif
!
!--------------------------------------------------------------------------------------------------------!
! STEP 6: VOPDERI
! Reconstruction: ZDERI(N+2,N+2) -> VOPDERI(N,N)
!
! Boundary conditions: linear extrapolation to nodes 0 and plev+1 derivative at nodes 1 and plev
!--------------------------------------------------------------------------------------------------------!
!
   zw1 = 1.0_r8
   zw2 = 0.0_r8
   do i = 1,plev
     do j = 2, plev-1
       vopderi(i,j) = zderi(i+1,j+1)
     enddo
       vopderi(i,1) = zderi(i+1,2)+zw1*zderi(i+1,1)
       vopderi(i,2) = zderi(i+1,3)+zw2*zderi(i+1,1)
       vopderi(i,plev-1) = zderi(i+1,plev)+zw2*zderi(i+1,plev+2)
       vopderi(i,plev) = zderi(i+1,plev+1)+zw1*zderi(i+1,plev+2)
   enddo
   if (debug_vfe) then
     write(*,'(a) ') ''
     write(*,'(a) ') ' derivative operator matrix:vopderi(n,n) '
     write(*,'(a) ') '----------------------------------------------------------------'
     do i = 1,plev
       write(*,'(50f6.2) ')(vopderi(i,j),j = 1,plev)
     enddo
   endif
!
   end subroutine sub_deri_operator
!
!---------------------------------------------------------------------------------!
! SUBROUTINE SUB_INTE_OPERATOR: Vertical integral operator (VOPINTE)
!---------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine sub_inte_operator(debug_vfe, etaf, knots, vopinte)
!-------------------------------------------------------------------------------
   implicit none
! Variables for hybrid coordinates
!
   logical, intent(in   ) :: debug_vfe
   real(r8), intent(in   ), dimension(0:plev+1) :: etaf
   real(r8), intent(in   ), dimension(nknot) :: knots
   real(r8), intent(  out), dimension(plev+1, plev) :: vopinte
! Variables for b-spline basis function
   real(r8), dimension(nlevbc, nbasis) :: bspline_funs
   real(r8), dimension(nquad*(plev+1), nbasis) :: ebspline_funs_quad
   real(r8), dimension(nquad*(plev+1), nbasis) :: dbspline_funs_quad
   real(r8), dimension(nquad*(plev+1), nbasis) :: jbspline_deri_quad
! Variables for VFE
   real(r8) :: zdet
   real(r8), dimension(nlevbc-1, nbasis-1) :: zp
   real(r8), dimension(nlevbc-1, nbasis-1) :: zs
   real(r8), dimension(nlevbc-1, nbasis-1) :: zsinv
   real(r8), dimension(nbasis-1, nbasis-1) :: za
   real(r8), dimension(nbasis-1, nbasis-1) :: zainv
   real(r8), dimension(nbasis-1, nbasis-1) :: zb
   real(r8), dimension(nlevbc-1, nbasis-1) :: zinte
   real(r8), dimension(nlevbc-1, nbasis-1) :: ztemp
   integer(i4) :: i, j, k, ilev, icount, iquad
   real(r8) :: deta, sum1, sum2
   real(r8), dimension(nquad) :: quadx, quadw
   real(r8), dimension(0:plev+1) :: zd
!
!--------------------------------------------------------------------------------------------------------!
! STEP 1: Matrix P: ZP
! ZP(i,j) = BASISFUNC_j(eta_i)) from FE space of the function to be differenciated into physical space
! i: row->eta points, j: column->basis function
!--------------------------------------------------------------------------------------------------------!
!
!
   call compute_basis_funs_eta(etaf,knots,bspline_funs)
   zp(:,:) = 0.0_r8
   do i = 2,nlevbc
     do j = 2, nbasis
       zp(i-1,j-1) = bspline_funs(i,j)
     enddo
   enddo
!
!--------------------------------------------------------------------------------------------------------!
! STEP 2: Matrix S^-1: ZSINV
! ZS(i,j)=BASISDER_j(eta_i) to project the derivative from FE space to physical space
!--------------------------------------------------------------------------------------------------------!
!
   zs(:,:) = zp(:,:)
   zdet = 0.0_r8
   call minv(zs,nlevbc-1,nbasis-1,zdet,zsinv)
   if (debug_vfe) then
     write(*,'(a) ') ''
     write(*,'(a) ') ' matrix zp'
     write(10,'(a) ') ' matrix zp'
     write(*,'(a) ') '----------------------------------------------------------------'
     do i = 1,nlevbc-1
       write(*,'(50f6.1) ')(100.*zp(i,j),j = 1,nbasis-1)
     enddo
   endif
!
!--------------------------------------------------------------------------------------------------------!
! STEP 3: Mass matrix A^-1
! Compute A_ik(eta): integrate the product of basis functions at all quadrature points
!--------------------------------------------------------------------------------------------------------!
!
! Compute basis functions E_i(eta),E_k(eta),and dE_j(eta)/deta at eta quadrature points (quadx)
!
   za(:,:) = 0.0_r8
   call define_quadrature_set(nquad,quadw,quadx)
   call compute_basis_funs_quad(etaf,knots,ebspline_funs_quad,jbspline_deri_quad)
   dbspline_funs_quad(:,:) = ebspline_funs_quad(:,:)
   dbspline_funs_quad(:,1) = 0.0_r8
   do i = 2,nbasis
     do k = 2, nbasis
       icount = 1
       za(i,k) = 0.0_r8
       do ilev = 1,nlevbc-1
         deta = etaf(ilev)-etaf(ilev-1)
         do iquad = 1,nquad
           za(i-1,k-1) = za(i-1,k-1)+0.5d0*deta*quadw(iquad)*dbspline_funs_quad(icount,i)*dbspline_funs_quad(icount,k)
           icount = icount+1
         enddo
       enddo
     enddo
   enddo
!
! Invert mass matrix: A^-1
!
   zdet = 0.0_r8
   call minv(za,nbasis-1,nbasis-1,zdet,zainv)
!ZAINV(:,:) = ZAINV(:,:)*NLEVEL
   if (debug_vfe) then
     write(*,*) ' determinant of a*nlev = ',zdet
     write(*,'(a) ') ''
     write(*,'(a) ') ' mass matrix:a'
     write(*,'(a) ') '----------------------------------------------------------------'
     do i = 1,nbasis-1
!  WRITE (*, '(50F20.10)') (ZA(i,j), j = 1, nbasis-1)
       write(*,'(50f6.1) ')(100.*za(i,j),j = 1,nbasis-1)
     enddo
     write(*,'(a) ') ''
     write(*,'(a) ') ' mass matrix:ainv'
     write(*,'(a) ') '----------------------------------------------------------------'
     do i = 1,nbasis-1
       write(*,'(50f6.0) ')(zainv(i,j),j = 1,nbasis-1)
     enddo
   endif
!
!--------------------------------------------------------------------------------------------------------!
! STEP 4: Integral matrix B
! Compute B_ik(eta): integrate the product of basis functions at all quadrature points
!--------------------------------------------------------------------------------------------------------!
!
   zb(:,:) = 0.0_r8
   do j = 2,nbasis
     do i = 1, nbasis-1
       icount = 1
       sum1 = 0.0_r8
       do ilev = 1,j
         deta = etaf(ilev)-etaf(ilev-1)
         do iquad = 1,nquad
           sum1 = sum1+0.5d0*deta*quadw(iquad)*ebspline_funs_quad(icount,i)
           icount = icount+1
         enddo
       enddo
       icount = 1
       sum2 = 0.0_r8
       do ilev = 1,nlevbc-1
         deta = etaf(ilev)-etaf(ilev-1)
         do iquad = 1,nquad
           sum2 = sum2+0.5d0*deta*quadw(iquad)*dbspline_funs_quad(icount,j)*sum1
           icount = icount+1
         enddo
       enddo
       zb(j-1,i) = sum2
     enddo
   enddo
!
!
! Temporarily adopted a part of ZB from verfe1i.f90 for comparison ---------------------------------------------!
!
   zd(:) = 0.0_r8
   do i = 1,plev+1
     if (i==1) then
       zd(i) = etaf(1)
     elseif (i==plev+1) then
       zd(i) = 1.0_r8-etaf(plev)
     else
       zd(i) = etaf(i)-etaf(i-1)
     endif
   enddo
   do i = 1,plev+1
     do j = 1, plev+1
       if (i==j+2) then
         zb(j,i) = zd(j+1)**2/24.0_r8
       elseif (i==j+1) then
         zb(j,i) = zd(j)**2/8.0_r8+zd(j)*zd(j+1)/4.0_r8+zd(j+1)**2/8.0_r8
       elseif (i==j) then
         zb(j,i) = zd(j-1)*zd(j)/4.0_r8+5.0_r8*zd(j)**2/24.0_r8+zd(j+1)*(zd(j-1)+zd(j))/4.0_r8
       endif
       enddo
   enddo
   zb(plev+1,1) = zd(1)*zd(plev+1)/4.0_r8
   zb(plev,plev+1) = zd(plev)**2/8.0_r8+zd(plev)*zd(plev+1)/4.0_r8+zd(plev+1)**2/6.0_r8
   zb(plev+1,plev+1) = zd(plev)*zd(plev+1)/4.0_r8+zd(plev+1)**2/3.0_r8
!
! Temporarily adopted a part of ZA from verfe1i.f90 for comparison ---------------------------------------------!
!
!
   if (debug_vfe) then
     write(*,'(a) ') ''
     write(*,'(a) ') ' integral matrix:zb*1e6'
     write(*,'(a) ') '----------------------------------------------------------------'
     do i = 1,nbasis-1
       write(*,'(50f6.1) ')(1.e6*zb(i,j),j = 1,nbasis-1)
     enddo
   endif
!
!--------------------------------------------------------------------------------------------------------!
! STEP 5:
! Compute integral operator ZINTE = P*A^-1*B*S^-1
!--------------------------------------------------------------------------------------------------------!
!
! B*S^-1
   zinte(:,:) = 0.0_r8
   do j = 1,nbasis-1
     do i = 1, nbasis-1
       do k = 1, nbasis-1
         zinte(i,j) = zinte(i,j)+zb(i,k)*zsinv(k,j)
       enddo
       enddo
   enddo
!
! A^-1*B*P^-1
   ztemp(:,:) = 0.0_r8
   do j = 1,nbasis-1
     do i = 1, nbasis-1
       do k = 1, nbasis-1
         ztemp(i,j) = ztemp(i,j)+zainv(i,k)*zinte(k,j)
       enddo
       enddo
   enddo
!
! ZDERI=P*A^-1*B*S^-1
   zinte(:,:) = 0.0d0
   do j = 1,nbasis-1
     do i = 1, nbasis-1
       do k = 1, nbasis-1
         zinte(i,j) = zinte(i,j)+zp(i,k)*ztemp(k,j)
       enddo
       enddo
   enddo
!
!--------------------------------------------------------------------------------------------------------!
! STEP 6:
! Reconstruction: ZINTE(N+1,N+1) -> VOPINTE(N+1,N)
!--------------------------------------------------------------------------------------------------------!
!
   vopinte(:,:) = 0.0_r8
   do i = 1,plev
     do j = 1, plev+1
       if (i==1) then
         vopinte(j,i) = zinte(j,i)+zinte(j,i+1)
       else
         vopinte(j,i) = zinte(j,i+1)
       endif
       enddo
   enddo
   if (debug_vfe) then
     write(*,'(a) ') ''
     write(*,'(a) ') ' integral operator matrix:vopinte(n+1,n) '
     write(*,'(a) ') '----------------------------------------------------------------'
     do j = 1,plev+1
       write(*,'(50f6.1) ')(vopinte(j,i),i = 1,plev)
     enddo
     zdet = 0.0
     call minv(zinte,nbasis-1,nbasis-1,zdet,ztemp)
     write(*,*) ' determintant of zinte = ',zdet
     write(*,'(a) ') ' inverse of matrix for integrals:'
     write(*,'(a) ') '----------------------------------------------------------------'
     do j = 1,plev+1
       write(*,'(50f6.1) ')(ztemp(j,i),i = 1,plev+1)
     enddo
     write(*,*) ' derivatives of function f(x) = 1'
     do j = 1,plev
       sum1 = 0.0d0
       do i = 1,plev+1
         sum1 = sum1+ztemp(j+1,i)
       enddo
       write(*,*) ' level ',j,' derivative = ',sum1
     enddo
   endif
!
   end subroutine sub_inte_operator
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module vfe_operators
!-------------------------------------------------------------------------------
