!-------------------------------------------------------------------------------
   module vfe_basis_funs
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
   use vfe_mod, only: debug_vfe, basis_degree, deri_order, nquad, basis_order, nlevbc, nknot, nbasis
   implicit none
!
   public :: compute_basis_funs_eta, &
                                                      compute_basis_funs_quad, &
                                                                   basis_funs, &
                                                                  basis_deris, &
                                                               find_knot_span, &
                                                          define_quadrature_set
!
   contains
!
!---------------------------------------------------------------------------------!
! SUBROUTINE COMPUTE_BASIS_FUNS_ETA: To compute basis functions at eta full levels
!
!---------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine compute_basis_funs_eta(etaf, knots, basis_fun)
!   IMPLICIT NONE
!
   real(r8), intent(in   ), dimension(0:plev+1) :: etaf
   real(r8), intent(in   ), dimension(nknot) :: knots
   real(r8), intent(  out), dimension(nlevbc, nbasis) :: basis_fun
   integer(i4) :: ip
   integer(i4) :: ileft
   real(r8) :: eta_point
   real(r8), dimension(nbasis) :: basisf
!
!----------------------------------------------------------------------------------------------!
! B-spline basis functions at eta full level points
!----------------------------------------------------------------------------------------------!
!
!
   if (debug_vfe) then
     write(*,'(a) ') ''
     write(*,'(a) ') ' b-spline basis functions at eta full level points'
     write(*,'(a) ') ' i left eta b1(x) b2(x) b3(x),.....'
     write(*,'(a) ') '----------------------------------------------------------------'
   endif
!
   basisf(:) = 0.0_r8
   basis_fun(:,:) = 0.0_r8
   do ip = 1,nlevbc
     eta_point = etaf(ip-1)
! Find left-most points with x-point from the knot
     call find_knot_span(nbasis,basis_order,knots,eta_point,ileft)
! Get B(I,K)(X) in VALUES(1:N): B(LEFT-K+1,K)(X),...,B(LEFT,K)(X). All the others are zero.
     call basis_funs(basis_order,ileft,knots,eta_point,basisf(ileft+1-basis_order))
     basis_fun(ip,:) = basisf(:)
     if (debug_vfe) write(*,'(2i4,f15.10,100f15.10) ') ip,ileft,eta_point,basisf(1:nbasis)
! Zero out the values just computed in preparation for the next evaluation point
     basisf(:) = 0.0_r8
   enddo
!
   end subroutine compute_basis_funs_eta
!
!---------------------------------------------------------------------------------!
! SUBROUTINE COMPUTE_BASIS_FUNS_QUAD: To compute basis functions at eta quadrature points
!
!---------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine compute_basis_funs_quad(etaf, knots, basis_fun_quad, basis_der_quad)
!   IMPLICIT NONE
!
   real(r8), intent(in   ), dimension(0:plev+1) :: etaf
   real(r8), intent(in   ), dimension(nknot) :: knots
   real(r8), intent(  out), dimension(nquad*(plev+1), nbasis) :: basis_fun_quad
   real(r8), intent(  out), dimension(nquad*(plev+1), nbasis) :: basis_der_quad
   integer(i4) :: ip, icount, iquad, ileft
   real(r8) :: eta_point
   real(r8), dimension(nquad) :: quadx
   real(r8), dimension(nquad) :: quadw
   real(r8), dimension(nbasis) :: basisf
   real(r8), dimension(nbasis) :: basisd
!
!----------------------------------------------------------------------------------------------!
! B-spline basis functions at eta full level points
!----------------------------------------------------------------------------------------------!
!
!
   call define_quadrature_set(nquad,quadw,quadx)
   if (debug_vfe) then
     write(*,'(a) ') ''
     write(*,'(a) ') ' b-spline basis functions at gauss quadrature points'
     write(*,'(a) ') ' i left eta b1(x) b2(x) b3(x),.....'
     write(*,'(a) ') '----------------------------------------------------------------'
   endif
   icount = 1
   basisf(:) = 0.0_r8
   basisd(:) = 0.0_r8
   do ip = 1,nlevbc-1
     do iquad = 1, nquad
       eta_point =((1.0_r8-quadx(iquad))*etaf(ip-1)+(quadx(iquad)-(-1.0_r8))*etaf(ip))/(1.0_r8-(-1.0_r8))
! Find left-most points with x-point from the knot
       call find_knot_span(nbasis,basis_order,knots,eta_point,ileft)
       call basis_deris(basis_order,deri_order,ileft,knots,eta_point,basisf(ileft+1-basis_order),basisd(ileft+1-basis_order))
       if (debug_vfe) write(*,'(2i4,f10.5,100f12.7) ') icount,ileft,eta_point,basisf(1:nbasis)
       basis_fun_quad(icount,:) = basisf(:)
       basis_der_quad(icount,:) = basisd(:)
       icount = icount+1
       basisf(:) = 0.0_r8
       basisd(:) = 0.0_r8
     enddo
   enddo
!
   end subroutine compute_basis_funs_quad
!
!---------------------------------------------------------------------------------!
! SUBROUTINE BASIS_FUNS
!
!---------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine basis_funs(basis_order, ileft, knots, pts, basisf)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: basis_order
   integer(i4), intent(in   ) :: ileft
   real(r8), intent(in   ) :: pts
   real(r8), intent(in   ), dimension(ileft+basis_order) :: knots
   real(r8), intent(  out), dimension(basis_order) :: basisf
   integer(i4) :: i, j, k
   integer(i4) :: basis_degree
   real(r8) :: saved, temp
   real(r8), dimension(basis_order-1) :: delta_left
   real(r8), dimension(basis_order-1) :: delta_rght
   real(r8), dimension(0:basis_order-1, 0:basis_order-1) :: ndu
!
   basis_degree = basis_order-1
!
!---------------------------------------------------------------------------------!
! To compute B-spline bassis function
!---------------------------------------------------------------------------------!
!
   ndu(0,0) = 1.0_r8
   do j = 1,basis_degree
     delta_left(j) = pts-knots(ileft+1-j)
     delta_rght(j) = knots(ileft+j)-pts
     saved = 0.0_r8
     do k = 0,j-1
       ndu(j,k) = delta_rght(k+1)+delta_left(j-k)
       temp = ndu(k,j-1)/ndu(j,k)
       ndu(k,j) = saved+delta_rght(k+1)*temp
       saved = delta_left(j-k)*temp
     enddo
     ndu(j,j) = saved
   enddo
! B-spline basis function
   do i = 1,basis_order
     basisf(i) = ndu(i-1,basis_degree)
   enddo
!
   end subroutine basis_funs
!
!---------------------------------------------------------------------------------!
! SUBROUTINE BASIS_DERIS
!
!---------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine basis_deris(basis_order, deri_order, ileft, knots, pts, basisf, basisd)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: basis_order
   integer(i4), intent(in   ) :: deri_order
   integer(i4), intent(in   ) :: ileft
   real(r8), intent(in   ) :: pts
   real(r8), intent(in   ), dimension(ileft+basis_order) :: knots
   real(r8), intent(  out), dimension(basis_order) :: basisf
   real(r8), intent(  out), dimension(basis_order) :: basisd
   integer(i4) :: i, j, k, l
   integer(i4) :: row1, row2, diff_lk, diff_pk, low, high
   integer(i4) :: basis_degree
   real(r8) :: saved, temp
   real(r8), dimension(basis_order-1) :: delta_left
   real(r8), dimension(basis_order-1) :: delta_rght
   real(r8), dimension(0:deri_order, 0:basis_order) :: deris
   real(r8), dimension(0:1, 0:basis_order) :: a
   real(r8), dimension(0:basis_order-1, 0:basis_order-1) :: ndu
!
!
   basis_degree = basis_order-1
!
!---------------------------------------------------------------------------------!
! To compute B-spline bassis function
!---------------------------------------------------------------------------------!
!
   ndu(0,0) = 1.0_r8
   do j = 1,basis_order-1
     delta_left(j) = pts-knots(ileft+1-j)
     delta_rght(j) = knots(ileft+j)-pts
     saved = 0.0_r8
     do k = 0,j-1
       ndu(j,k) = delta_rght(k+1)+delta_left(j-k)
       temp = ndu(k,j-1)/ndu(j,k)
       ndu(k,j) = saved+delta_rght(k+1)*temp
       saved = delta_left(j-k)*temp
     enddo
     ndu(j,j) = saved
   enddo
!
!---------------------------------------------------------------------------------!
! To compute the k-th order derivative of B-spline bassis function
!---------------------------------------------------------------------------------!
!
   do j = 0,basis_degree
     deris(0,j) = ndu(j,basis_degree)
   enddo
   do l = 0,basis_degree
     row1 = 0
     row2 = 1
     a(0,0) = 1.0_r8
     do k = 1,deri_order
       temp = 0.0_r8
       diff_lk = l-k
       diff_pk = basis_degree-k
       if (diff_lk>= 0) then
         a(row2,0) = a(row1,0)/ndu(diff_pk+1,diff_lk)
         temp = a(row2,0)*ndu(diff_lk,diff_pk)
       endif
       if (diff_lk>=-1) then
         low = 1
       else
         low =-diff_lk
       endif
       if (l-1<=diff_pk) then
         high = k-1
       else
         high = basis_degree-l
       endif
       do j = low,high
         a(row2,j) =(a(row1,j)-a(row1,j-1))/ndu(diff_pk+1,diff_lk+j)
         temp = temp+a(row2,j)*ndu(diff_lk+j,diff_pk)
       enddo
       if (l<=diff_pk) then
         a(row2,k) =-a(row1,k-1)/ndu(diff_pk+1,l)
         temp = temp+a(row2,k)*ndu(l,diff_pk)
       endif
       deris(k,l) = temp
       j = row1
       row1 = row2
       row2 = j
     enddo
   enddo
   l = basis_degree
   do k = 1,deri_order
     do j = 0, basis_degree
       deris(k,j) = deris(k,j)*l
     enddo
       l = l*(basis_degree-k)
   enddo
   do i = 1,basis_order
     basisf(i) = deris(0,i-1) ! fisrt derivative of b-spline basis
     basisd(i) = deris(1,i-1) ! fisrt derivative of b-spline basis
   enddo
!
   end subroutine basis_deris
!
!---------------------------------------------------------------------------------!
! SUBROUTINE FIND_KNOT_SPAN
!
!---------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine find_knot_span(nbasis, basis_order, knots, pts, ileft)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), parameter :: tolerance = 5.d-15
   integer(i4), intent(in   ) :: nbasis
   integer(i4), intent(in   ) :: basis_order
   real(r8), intent(in   ) :: pts
   real(r8), intent(in   ), dimension(nbasis+basis_order) :: knots
   integer(i4) :: ileft
   integer(i4) :: basis_degree, nknot, low, high, mid
!
   nknot = nbasis+basis_order
   basis_degree = basis_order-1
   if ((pts<knots(1)-tolerance).or.(pts>knots(nknot)+tolerance)) then
     write(*,'(3x,a,f10.5) ') 'STOP! the point is out of range:',pts
     stop
     ileft =-1
   elseif (dabs(pts-knots(nknot-basis_degree))<tolerance) then
     ileft = nknot-basis_order
   else
     low = basis_degree
     high = nknot-basis_degree
     mid =(low+high)/2
     do while(pts<knots(mid).or.pts>= knots(mid+1))
       if (pts<knots(mid)) then
         high = mid
       else
         low = mid
       endif
         mid =(low+high)/2
     enddo
     ileft = mid
   endif
!
   end subroutine find_knot_span
!
!---------------------------------------------------------------------------------!
! SUBROUTINE DEFINE_QUADRATURE_SET
!
!---------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine define_quadrature_set(nquad, quadw, quadx)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: nquad
   real(r8), dimension(nquad), intent(  out) :: quadw, quadx
!
   if (nquad==1) then
     quadx(1) = 0.0_r8
     quadw(1) = 2.0_r8
   elseif (nquad==2) then
     quadx(1) =-0.5773502691896257645091487805020_r8
     quadx(2) = 0.5773502691896257645091487805020_r8
     quadw(1) = 1.0_r8
     quadw(2) = 1.0_r8
   elseif (nquad==3) then
     quadx(1) =-0.7745966692414833770358530799560_r8
     quadx(2) = 0.0_r8
     quadx(3) = 0.7745966692414833770358530799560_r8
     quadw(1) = 5.0_r8/9.0_r8
     quadw(2) = 8.0_r8/9.0_r8
     quadw(3) = 5.0_r8/9.0_r8
   elseif (nquad==4) then
     quadx(1) =-0.8611363115940525752239464888930_r8
     quadx(2) =-0.3399810435848562648026657591030_r8
     quadx(3) = 0.3399810435848562648026657591030_r8
     quadx(4) = 0.8611363115940525752239464888930_r8
     quadw(1) = 0.3478548451374538573730639492220_r8
     quadw(2) = 0.6521451548625461426269360507780_r8
     quadw(3) = 0.6521451548625461426269360507780_r8
     quadw(4) = 0.3478548451374538573730639492220_r8
   elseif (nquad==5) then
     quadx(1) =-0.9061798459386639927976268782990_r8
     quadx(2) =-0.5384693101056830910363144207000_r8
     quadx(3) = 0.0_r8
     quadx(4) = 0.5384693101056830910363144207000_r8
     quadx(5) = 0.9061798459386639927976268782990_r8
     quadw(1) = 0.2369268850561890875142640407200_r8
     quadw(2) = 0.4786286704993664680412915148360_r8
     quadw(3) = 0.5688888888888888888888888888890_r8
     quadw(4) = 0.4786286704993664680412915148360_r8
     quadw(5) = 0.2369268850561890875142640407200_r8
   elseif (nquad==6) then
     quadx(1) =-0.9324695142031520278123015544940_r8
     quadx(2) =-0.6612093864662645136613995950200_r8
     quadx(3) =-0.2386191860831969086305017216810_r8
     quadx(4) = 0.2386191860831969086305017216810_r8
     quadx(5) = 0.6612093864662645136613995950200_r8
     quadx(6) = 0.9324695142031520278123015544940_r8
     quadw(1) = 0.1713244923791703450402961421730_r8
     quadw(2) = 0.3607615730481386075698335138380_r8
     quadw(3) = 0.4679139345726910473898703439900_r8
     quadw(4) = 0.4679139345726910473898703439900_r8
     quadw(5) = 0.3607615730481386075698335138380_r8
     quadw(6) = 0.1713244923791703450402961421730_r8
   endif
!
   end subroutine define_quadrature_set
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module vfe_basis_funs
!-------------------------------------------------------------------------------
