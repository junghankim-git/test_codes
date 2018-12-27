#include "KIM.h"
!-------------------------------------------------------------------------------
   module vertical
! ------------------------------
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
   use kiapsbase, only: real_kind=>kim_real8_kind
! ------------------------------
   use dimensions, only: nlev
! ------------------------------
   implicit none
!
   private
!
   type, public :: vcoord_t
     sequence
     real(real_kind) :: a(nlev+1) ! hybrid pressure coordinate coefficients on interfaces
     real(real_kind) :: b(nlev+1) ! hybrid sigma coordinate(dp/deta) coefficient
     real(real_kind) :: amid(nlev) ! hybrid pressure coordinate coefficients on mid levels
     real(real_kind) :: bmid(nlev) ! hybrid sigma coordinate(dp/deta) coefficient on mid levels
     real(real_kind) :: db(nlev) ! delta b on mid levels
     real(real_kind) :: hmat(nlev, nlev) ! hydrostatic integral matrix(sigma coordinates)
     real(real_kind) :: cmat(nlev, nlev) ! energy convergence integral matrix(sigma coordinates)
   end type
   public :: hmat_init
   public :: cmat_init
   public :: vcoord_init
   public :: vcoord_print
   public :: hmat_print
   private :: ecmwf_hmat_init
   private :: ecmwf_cmat_init
   private :: ccm_hmat_init
   private :: ccm_cmat_init
!
   contains
! ====================================================
! subroutine ecmwf_Bmid encapsulates the
! ECMWF method for computing sigma mid levels.
!
! reference: Ritchie,et.al. and Burridge and Simmons.
! ====================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ecmwf_bmid(sigma, dsigma, sigmam)
!
   real(real_kind) :: sigma(nlev+1) ! interface sigma level
   real(real_kind) :: dsigma(nlev) ! delta sigma level
   real(real_kind) :: sigmam(nlev) ! mid point sigma level
! Local variables
   integer k
   real(real_kind) :: c
!
   c = 1.0d0
   sigmam(1) = 0.5d0*dsigma(1)
   do k = 2,nlev
     sigmam(k) = exp((1.0d0/dsigma(k))*(sigma(k+1)*log(sigma(k+1)) &
                                  -sigma(k)*log(sigma(k)))-c)
   enddo
!
   end subroutine ecmwf_bmid
! ====================================================
! subroutine ccm1_sigma_mid encapsulates the
! CCM1 method for computing sigma mid levels.
!
! reference: ???
! ====================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ccm1_bmid(sigma, dsigma, sigmam)
!
   real(real_kind) :: sigma(nlev+1) ! interface sigma level
   real(real_kind) :: dsigma(nlev) ! delta sigma level
   real(real_kind) :: sigmam(nlev) ! mid point sigma level
   integer k
!
   do k = 1, nlev
     sigmam(k) = 0.5d0*(sigma(k+1)+sigma(k))
   enddo
!
   end subroutine ccm1_bmid
!  ================================================
!  cmat_init:
!
!  ================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function cmat_init(p, k, l, npts) result(cmat)
!
   integer, intent(in   ) :: npts
   real(real_kind), intent(in   ) :: p(npts, npts, nlev)
   integer, intent(in   ) :: k, l
   real(real_kind) :: cmat(npts, npts)
   integer :: i, j
!
   if (l<k) then
     cmat(:,:) = 1.0d0/p(:,:,k)
   elseif (l==k) then
     cmat(:,:) = 0.5d0/p(:,:,k)
   else
     cmat(:,:) = 0.0d0
   endif
!
   end function cmat_init
!  ================================================
!  hmat_init:
!
!  Compute the Hydrostatic Matrix Element H(k,l)
!  using the relationships in eq. (3.a.109) of the
!  CCM-2 Description (June 1993), NCAR/TN-382+STR
!
!  ================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function hmat_init(logps, logp) result(hmat)
!-------------------------------------------------------------------------------
   use dimensions, only: np
!
   real(real_kind), intent(in   ) :: logps(np, np)
   real(real_kind), intent(in   ) :: logp(np, np, nlev)
   integer :: k, m
   real(real_kind) :: hmat(np, np, nlev, nlev)
   integer :: i, j
!
   do k = 2, nlev
!      if(m<k)
     do m = 1, k-1
       hmat(:,:,m,k) = 0.0_real_kind
     enddo
   enddo
!
   do k = 1, nlev-1
!      if(m==k)
     do j = 1, np
       do i = 1, np
         hmat(i,j,k,k) = 0.50d0*(logp(i,j,k+1)-logp(i,j,k))
       enddo
       enddo
!      m>k<nlev
         do m = k+1, nlev-1
           do j = 1, np
             do i = 1, np
               hmat(i,j,m,k) = 0.50d0*(logp(i,j,m+1)-logp(i,j,m-1))
             enddo
             enddo
             enddo
!      m==nlev!=k
               do j = 1, np
                 do i = 1, np
                   hmat(i,j,nlev,k) = logps(i,j)-0.50d0*(logp(i,j,nlev)+logp(i,j,nlev-1))
                 enddo
                 enddo
   enddo
!   m==k==nlev
!
   do j = 1, np
     do i = 1, np
       hmat(i,j,nlev,nlev) = logps(i,j)-logp(i,j,nlev)
     enddo
   enddo
!
   end function hmat_init
! ======================================================
! vcoord_print:
!
! print hybrid coordinate values (A and B coefficients)
! ======================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine vcoord_print(vert, formulation)
!
   type(vcoord_t) :: vert
   character(len=*) formulation ! formulation choice(ecmwf, ccm1)
   integer k
!
   write(6,10) formulation
   10 format(a6," hybrid vertical coordinate formulation")
   write(6,*) "---------------------------------------"
   write(6,*) "interfaces"
   do k = 1,nlev+1
     write(6,11) k-1,vert%a(k),k-1,vert%b(k)
   11 format("a(",i4,"+1/2) = ",e22.15," b(",i4,"+1/2) = ",e22.15)
   enddo
   write(6,*) ""
   write(6,*) "---------------------------------------"
   write(6,*) "mid levels:"
   do k = 1,nlev
     write(6,12) k,vert%amid(k),k-1,vert%bmid(k)
   12 format("a(",i4,") = ",e22.15," b(",i4,") = ",e22.15)
   enddo
   write(6,*) ""
!
   end subroutine vcoord_print
! =============================================================
! vcoord_init: Initialize vertical coordinate system...
! returns: 0 for success
! =============================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function vcoord_init(vert, formulation, vfile_mid, vfile_int) result(ierr)
!
   integer :: ierr
   type(vcoord_t) :: vert
   character(len=*) formulation ! formulation choice(ecmwf, ccm)
   character(len=*) vfile_mid ! mid levels vertical coordinate system
   character(len=*) vfile_int ! interface levels vertical coordinate system
! Local variables
   real(real_kind) :: s(nlev+1)
   integer ierr11, ierr12
   integer k
!
   ierr = 0
!$OMP CRITICAL
   open(unit = 11,file = vfile_int,form = 'unformatted',status = 'old',iostat = ierr11)
   if (ierr11/= 0) then
     ierr = ierr11
   endif
   open(unit = 12,file = vfile_mid,form = 'unformatted',status = 'old',iostat = ierr12)
   if (ierr12/= 0) then
     ierr = ierr12
   endif
   do k = 1,nlev+1
     read(11) vert%a(k)
     read(11) vert%b(k)
   enddo
   do k = 1,nlev
     read(12) vert%amid(k)
     read(12) vert%bmid(k)
   enddo
   close(11)
   close(12)
!$OMP END CRITICAL
   print*,"successfully read vertical coordinates"
   do k = 1,nlev
     vert%db(k) = vert%b(k+1)-vert%b(k)
   enddo
!
   end function vcoord_init
! ===============================================
! hmat_print prints the entries of the
! CCM1 Hydrostatic matrix (Hmat) in the case of
! sigma coordinates.
! ===============================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine hmat_print(vert)
!
   type(vcoord_t) :: vert
! Local variables
   integer k, l
!
   print*
   print*,"ccm hydrostatic matrix"
   do k = 1,nlev
     do l = 1, nlev
       write(6,10) k,l,vert%hmat(k,l)
   10 format("b(",i4,",",i4,") = ",e22.15)
     enddo
       print*
   enddo
!
   end subroutine hmat_print
! ===============================================
! ccm_hmat_init fills in entries of CCM1 H matrix.
! This applies only to the sigma coordinate
! implementation of CCM-1.
!
! see Description of CCM1 (NCAR/TN-285+STR),
! page 40-41 for details.
! ===============================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ccm_hmat_init(vert)
!
   type(vcoord_t) :: vert
! Local variables
   integer k, l
! ===================================
! eq (3.b.14) of Description of CCM1
!                (NCAR/TN-285+STR)
! ===================================
!
   vert%hmat(nlev,nlev) =-log(vert%bmid(nlev))
! ===================================
! eq (3.b.15),ibid.
! ===================================
   do k = 1,nlev-1
     vert%hmat(nlev,k) = 0.0d0
   enddo
! ===================================
! eq (3.b.16),ibid.
! ===================================
   do k = 1,nlev-1
     vert%hmat(k,k) = 0.50d0*log(vert%bmid(k+1)/vert%bmid(k))
   enddo
! ===================================
! eq (3.b.17),ibid.
! ===================================
   do k = 1,nlev-1
     vert%hmat(k,k+1) = vert%hmat(k+1,k+1)+0.50d0*log(vert%bmid(k+1)/vert%bmid(k))
   enddo
! ===================================
! eq (3.b.18),ibid.
! ===================================
   do l = 3,nlev
     do k = l-2, 1,-1
       vert%hmat(k,l) = vert%hmat(k+1,l)
     enddo
   enddo
! ===================================
! eq (3.b.19),ibid.
! ===================================
   do k = 1,nlev
     do l = 1, k-1
       vert%hmat(k,l) = 0.0d0
     enddo
   enddo
!
   end subroutine ccm_hmat_init
! ===============================================
! ccm1_cmat_init fills in entries of CCM1 C matrix.
! see Description of CCM1 (NCAR/TN-285+STR),
! page 41, eq. (3.b.21) for details.
!
! ===============================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ccm_cmat_init(vert)
!
   type(vcoord_t) :: vert
! Local variables
   integer j, k
!
   do k = 1, nlev
     do j = 1, nlev
       vert%cmat(k,j) = vert%hmat(j,k)*(vert%db(j)/vert%db(k))
     enddo
   enddo
!
   end subroutine ccm_cmat_init
! ===============================================
! ecmwf_hmat_init fills in entries of ECMWF hydrostatic
! integral or B matrix.
!  (see RDL nodes deriving this from
!   eqs. (2.21) and (2.22) of Ritchie, et. al.)
! ===============================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ecmwf_hmat_init(vert)
!
   type(vcoord_t) :: vert
! Local variables
   integer j, k
!
   do k = 1, nlev
     do j = 1, k-1
       vert%hmat(k,j) = 0.0d0
     enddo
       vert%hmat(k,k) = 1.0d0-(vert%b(k)/vert%db(k))*log(vert%b(k+1)/vert%b(k))
       do j = k+1,nlev
         vert%hmat(k,j) = log(vert%b(j+1)/vert%b(j))
       enddo
   enddo
! ======================================================
! Special case (see Ritchie discussion of alpha(1))
! ======================================================
!
   vert%hmat(1,1) = log(2.0d0)
!
   end subroutine ecmwf_hmat_init
! ===============================================
! ecmwf_cmat_init fills in entries of ECMWF energy
! conversion integral, or C matrix.
!  (see RDL nodes deriving this from
!   eqs. (2.25) of Ritchie, et. al.)
! ===============================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ecmwf_cmat_init(vert)
!
   type(vcoord_t) :: vert
! Local variables
   integer j, k
!
   do k = 1, nlev
     do j = 1, k-1
       vert%cmat(k,j) =(vert%db(j)/vert%db(k))*log(vert%b(k+1)/vert%b(k))
     enddo
       vert%cmat(k,k) = 1.0d0-(vert%b(k)/vert%db(k))*log(vert%b(k+1)/vert%b(k))
       do j = k+1,nlev
         vert%cmat(k,j) = 0.0d0
       enddo
   enddo
! ======================================================
! Special case (see Ritchie discussion of alpha(1))
! ======================================================
!
   vert%cmat(1,1) = log(2.0d0)
!
   end subroutine ecmwf_cmat_init
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module vertical
!-------------------------------------------------------------------------------
