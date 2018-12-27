!-------------------------------------------------------------------------------
!
! MODULE Control
!
!> @brief
!
!> @date 09MAR2015
!>  - Ja-Rin Park : Implementation of VFE scheme (Vertical coordinate reconstructed)
!> @date 01OCT2015
!>  - Junghan Kim : Add the p_top
!
!-------------------------------------------------------------------------------
#include "KIM.h"
!-------------------------------------------------------------------------------
   module hybvcoord
!```````````````````````````````````````````````````````
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
   use kiapsbase, only: iulog=>kim_iu_log, r8=>kim_real8_kind
   use dimensions, only: nlev, plev=>nlev, plevp=>nlevp
   use physicalconstants, only: p1000mb=>kim_p0, p0=>kim_p0
   use control, only: on_vfe, vfebasis_order
   use vfe_mod, only: debug_vfe, basis_order, basis_degree, nlevbc, nbasis, nknot, &
vopinte, vopderi, ini_vfe, fin_vfe, sub_verint, sub_verder, vfe_mass_correction
   use vfe_operators, only: sub_deri_operator, sub_inte_operator
!```````````````````````````````````````````````````````
   implicit none
!-------------------------------------------------------
!
   private
   real(r8), public :: p_top
#if defined(CORE_SW)
!-----------------------------------------------------------------------
!
! Suk-Jin CHOI
!
! Sequence....Bottom(1) to TOP(nlev) [B2T]
!
! Prototype definition :
!  pd(k) = hyai(k)*(ps0-ptop)+hybi(k)*(pds-ptop)+ptop
!  dp/deta(k) = hyai_e(k)*(ps0-ptop)+hybi_e(k)*(ps-ptop)
!
! There is several options for vertical coordinate
!
! option(1-0) (Default: #define vcoord_sigma)
!     : Pure Sigma coordinate with non-zero ptop
!     : eta <- hyai(k)+hybi(k), and then hyai = 0.
!     : p(k) = eta(k)*(ps0-ptop)+ptop
!     : dp/deta(k) = ps0-ptop
!
! option(2-0) (#define vcoord_hybridsigma, #undef hysig1)
!     : Hybrid-Sigma coordinate
!     : dp/deta(k) at mid-level = hyai_e(k)*(ps0-ptop)+hybi_e(k)*(ps-ptop)
!     : dp/deta(k) at int-level are interpolated using dp/deta(k) at mid-level
!
! option(2-1) (#define vcoord_hybridsigma, #define hysig1)
!     : Hybrid-Sigma coordinate
!     : dp/deta(k) at mid-level = hyai_e(k)*(ps0-ptop)+hybi_e(k)*(ps-ptop)
!     : dp/deta(k) at int-level = hyam_e(k)*(ps0-ptop)+hybm_e(k)*(ps-ptop)
!     : hyam_e(k) = (hyam(k) - hyam(k-1))/deta
!
! ------------------------------------------------------
!
!         TOP
!    ---------------  nlev+1/2 (nlev+1)
!    ===============  nlev     (nlev)
!    ---------------  nlev-1/2 (nlev)
!
!    ===============  nlev-1   (nlev-1)
!
!    ---------------  nlev-3/2 (nlev-1)
!          .
!          .
!          .
!          .
!    ===============   2+1/2 (3)
!
!    ---------------   2     (2)
!
!    ===============   1+1/2 (2)
!    ---------------   1     (1)
!    ===============    +1/2 (1)
!         SFC
!-----------------------------------------------------------------------
!
   type, public :: hvcoord_t
     real(r8) ps0 ! base state sfc pressure for level definitions
     real(r8) hyai(nlev+1) ! ps0 component of hybrid coordinate-interfaces
     real(r8) hyam(nlev) ! ps0 component of hybrid coordinate-midpoints
     real(r8) hybi(nlev+1) ! ps component of hybrid coordinate-interfaces
     real(r8) hybm(nlev) ! ps component of hybrid coordinate-midpoints
     real(r8) znw(nlev+1) ! etai = hyai+hybi
     real(r8) znu(nlev) ! etam = hyam+hybm
     real(r8) dnw(nlev)
     real(r8) rdnw(nlev)
     real(r8) dn(nlev)
     real(r8) rdn(nlev)
     real(r8) fnp(nlev)
     real(r8) fnm(nlev)
     real(r8) hyai_e(nlev)
     real(r8) hybi_e(nlev)
     real(r8) hyam_e(nlev) ! level of bottom(1) not in used..
     real(r8) hybm_e(nlev)
     real(r8) p_top
     real(r8) cf1
     real(r8) cf2
     real(r8) cf3
     real(r8) cfn
     real(r8) cfn1
   end type
#elif defined(CORE_SH)
!-----------------------------------------------------------------------
!
! Purpose: Hybrid level definitions: p = a*p0 + b*ps
!          interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
!          midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
!
!-----------------------------------------------------------------------
!
   type, public :: hvcoord_t
     real(r8) ps0 ! base state sfc pressure for level definitions
     real(r8) hyai(plevp) ! ps0 component of hybrid coordinate-interfaces
     real(r8) hyam(plev) ! ps0 component of hybrid coordinate-midpoints
     real(r8) hybi(plevp) ! ps component of hybrid coordinate-interfaces
     real(r8) hybm(plev) ! ps component of hybrid coordinate-midpoints
     real(r8) hybd(plev) ! difference in b(hybi) across layers
     real(r8) dhyam(plev) ! ps0 differential component of hybrid coordinate-midpoints
     real(r8) dhybm(plev) ! ps differential component of hybrid coordinate-midpoints
     real(r8) prsfac ! log pressure extrapolation factor(time, space independent)
     real(r8) etam(plev) ! eta. stored for conviencience
     real(r8) etai(plevp) !
     integer nprlev ! number of pure pressure levels at top
     integer pad
   end type
#endif
! ---- for physics -----------------------------------
   real(r8), public :: hyai1d(plevp), hybi1d(plevp), si1d(plevp), sl1d(plev)
! ----------------------------------------------------
   public :: hvcoord_init
!-------------------------------------------------------
!
   contains
!==========================================================
! =============================================================
! hvcoord_init: Initialize hybrid vertical coordinate system...
! returns: hv coordinate structure, ierr=0 means success
!
! =============================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function hvcoord_init(hvfile_mid, hvfile_int, lprint, masterproc, ierr) result(hvcoord)
!--------------------
!
   character(len=*), intent(in   ) :: hvfile_mid ! file containing mid levels vertical coordinate system
   character(len=*), intent(in   ) :: hvfile_int ! file containing interface levels vertical coordinate system
   logical, intent(in   ) :: lprint
   logical, intent(in   ) :: masterproc
   integer, intent(  out) :: ierr
   type(hvcoord_t) :: hvcoord
!local
   integer :: k
   integer :: ierr11, ierr12, plevp_in, plev_in, ln
   real(r8) :: hyai_loc(nlev+1)
   real(r8) :: hybi_loc(nlev+1)
   real(r8) :: hyam_loc(nlev)
   real(r8) :: hybm_loc(nlev)
   real(r8) :: cof1, cof2
   real(r8) :: ci(nlev+1)
   real(r8) :: del(nlev)
#if defined(CORE_SH)
   real(r8) amean, bmean, atest, btest, eps
#endif
   real(r8), dimension(0:nlev+1) :: vetaf
   real(r8), dimension(nknot) :: knots
!--------------------
!
   ierr = 0
   hvcoord%ps0 = p1000mb ! base state surface pressure [pa]
!sjchoi  ---------------------------------------------------------------------
#if defined(CORE_SW)
   hvcoord%p_top = p_top ! [pa] default
#if defined(IGW)  /* idealized test */
   hvcoord%p_top = 27358.273_r8 ! [pa] for igw ~10km top
#elif defined(ASP)
   hvcoord%p_top = 0._r8 ! [pa] for asp 2008
#elif defined(DCMIP_case1X)
   hvcoord%p_top = 25494.4_r8 ! [pa] for dcmip 11,12,13
#elif defined(DCMIP_case20)
   hvcoord%p_top = 20544.8_r8 ! [pa] for dcmip 20
#elif defined(DCMIP_case21) || defined(DCMIP_case22)
   hvcoord%p_top = 3281.8_r8 ! [pa] for dcmip 21,22
#elif defined(DCMIP_case31)
   hvcoord%p_top = 27391.9_r8 ! [pa] for dcmip 31
#endif /* idealized test */
#endif
!-----------------------------------------------------------------------------
   ln = len(trim(hvfile_int))
! check if file name ends with .ascii
   if (hvfile_int(ln-4:ln)=='ascii') then
     open(unit = 11,file = hvfile_int,form = 'formatted',status = 'old',iostat = ierr11)
     if (ierr11/= 0) then
       write(iulog,*) 'open() error:',__file__,__line__,hvfile_int
       ierr = ierr11
     endif
     if (ierr==0) then
       read(11,*) plevp_in
       if (plevp_in.ne.plevp) then
         write(iulog,*) 'Error:hyai input file and homme nlev+1 do not match',plevp,plevp_in
         ierr = 1
       endif
       read(11,*) hyai_loc(1:plevp)
       read(11,*) plevp_in
       if (plevp_in.ne.nlev+1) then
         write(iulog,*) 'Error:hybi input file and homme nlev+1 do not match',plevp,plevp_in
         ierr = 1
       endif
       read(11,*) hybi_loc(1:plevp)
       close(11)
       read(12,*) plev_in
       if (plev_in.ne.plev) then
         write(iulog,*) 'Error:hyam input file and homme plev do not match',plev,plev_in
         ierr = 1
       endif
       read(12,*) hyam_loc(1:plev)
       read(12,*) plev_in
       if (plev_in.ne.plev) then
         write(iulog,*) 'Error:hybi input file and homme plev do not match',plev,plev_in
         ierr = 1
       endif
       read(12,*) hybm_loc(1:plev)
       close(12)
#if defined(CORE_SH)
       do k = 1,nlev+1 ! top to bottom
         hvcoord%hyai(k) = hyai_loc(k)
         hvcoord%hybi(k) = hybi_loc(k)
! -------- for physics --------------------------
         hyai1d(k) = hyai_loc(nlev+2-k) ! bottom to top in physics
         hybi1d(k) = hybi_loc(nlev+2-k) ! bottom to top in physics
! ------------------------------------------------------
       enddo
       do k = 1,nlev ! top to bottom
         hvcoord%hyam(k) = hyam_loc(k)
         hvcoord%hybm(k) = hybm_loc(k)
       enddo
#endif
#if defined(CORE_SW)
#if defined(vcoord_sigma)
       if (masterproc) then
         write(iulog,'(a) ') '================================================'
         write(iulog,'(a) ') 'YOU are using sigma vertical coord.'
         write(iulog,'(a) ') '================================================'
       endif
       do k = 1,nlev+1 ! bottom to top
         hvcoord%hyai(k) = 0._r8
         hvcoord%hybi(k) = hyai_loc(nlev+1-k+1)+hybi_loc(nlev+1-k+1)
       enddo
       do k = 1,nlev ! bottom to top
         hvcoord%hyam(k) = 0._r8
         hvcoord%hybm(k) = .5_r8*(hvcoord%hyai(k+1)+hvcoord%hyai(k))+.5_r8*(hvcoord%hybi(k+1)+hvcoord%hybi(k))
       enddo
#elif defined(vcoord_hybridsigma)
       if (masterproc) then
         write(iulog,'(a) ') '================================================'
         write(iulog,'(a) ') 'YOU are using hybrid-sigma vertical coord.'
         write(iulog,'(a) ') '================================================'
       endif
       do k = 1,nlev+1 ! bottom to top
         hvcoord%hyai(k) = hyai_loc(nlev+1-k+1)
         hvcoord%hybi(k) = hybi_loc(nlev+1-k+1)
       enddo
       do k = 1,nlev ! bottom to top
         hvcoord%hyam(k) = .5_r8*(hvcoord%hyai(k+1)+hvcoord%hyai(k))
         hvcoord%hybm(k) = .5_r8*(hvcoord%hybi(k+1)+hvcoord%hybi(k))
       enddo
#endif /*  -- vcoord_sigma */
#endif
     endif
   else
     open(unit = 11,file = hvfile_int,form = 'unformatted',status = 'old',iostat = ierr11)
     if (ierr11/= 0) then
       write(iulog,*) 'open() error:',__file__,__line__,hvfile_int
       ierr = ierr11
     endif
     open(unit = 12,file = hvfile_mid,form = 'unformatted',access = 'sequential',status = 'old',iostat = ierr12)
     if (ierr12/= 0) then
       write(iulog,*) 'open() error:',__file__,__line__,hvfile_mid
       ierr = ierr12
     endif
     if (ierr==0) then
       do k = 1, nlev+1
         read(11) hyai_loc(k)
         read(11) hybi_loc(k)
       enddo
         close(11)
         do k = 1,nlev
           read(12) hyam_loc(k)
           read(12) hybm_loc(k)
         enddo
         close(12)
#if defined(CORE_SH)
         do k = 1,nlev+1 ! top to bottom
           hvcoord%hyai(k) = hyai_loc(k)
           hvcoord%hybi(k) = hybi_loc(k)
         enddo
         do k = 1,nlev ! top to bottom
           hvcoord%hyam(k) = hyam_loc(k)
           hvcoord%hybm(k) = hybm_loc(k)
         enddo
#endif
#if defined(CORE_SW)
#if defined(vcoord_sigma)
         if (masterproc) then
           write(iulog,'(a) ') '================================================'
           write(iulog,'(a) ') 'YOU are using sigma vertical coord.'
           write(iulog,'(a) ') '================================================'
         endif
         do k = 1,nlev+1 ! bottom to top
           hvcoord%hyai(k) = 0._r8
           hvcoord%hybi(k) = hyai_loc(nlev+1-k+1)+hybi_loc(nlev+1-k+1)
         enddo
         do k = 1,nlev ! bottom to top
           hvcoord%hyam(k) = 0._r8
           hvcoord%hybm(k) = .5_r8*(hvcoord%hyai(k+1)+hvcoord%hyai(k))+.5_r8*(hvcoord%hybi(k+1)+hvcoord%hybi(k))
         enddo
#elif defined(vcoord_hybridsigma)
         if (masterproc) then
           write(iulog,'(a) ') '================================================'
           write(iulog,'(a) ') 'YOU are using hybrid-sigma vertical coord.'
           write(iulog,'(a) ') '================================================'
         endif
         do k = 1,nlev+1 ! bottom to top
           hvcoord%hyai(k) = hyai_loc(nlev+1-k+1)
           hvcoord%hybi(k) = hybi_loc(nlev+1-k+1)
         enddo
         do k = 1,nlev ! bottom to top
           hvcoord%hyam(k) = .5_r8*(hvcoord%hyai(k+1)+hvcoord%hyai(k))
           hvcoord%hybm(k) = .5_r8*(hvcoord%hybi(k+1)+hvcoord%hybi(k))
         enddo
#endif /*  -- vcoord_sigma */
#endif
     endif
   endif
   if (ierr>0) return
#if defined(CORE_SH)
   eps = 1.d-05
   hvcoord%nprlev = 0
   do k = 1,plev
     hvcoord%etam(k) = hvcoord%hyam(k)+hvcoord%hybm(k)
! ------ for physics --------------------
     sl1d(k) = hvcoord%hyam(plev+1-k)+hvcoord%hybm(plev+1-k) ! bottom to top in physics
! ----------------------------------------
   enddo
   do k = 1,plevp
     hvcoord%etai(k) = hvcoord%hyai(k)+hvcoord%hybi(k)
! ------ for physics --------------------
     si1d(k) = hvcoord%hyai(plevp+1-k)+hvcoord%hybi(plevp+1-k) ! bottom to top in physics
! ----------------------------------------
   enddo
! ====================================
! Set vertical operator matrix for VFE
! ====================================
   if (on_vfe) then
     if (masterproc) then
       write(iulog,'(a) ') '================================================'
       write(iulog,'(a) ') 'Set vfe coordinate'
       write(iulog,'(a) ') '================================================'
     endif
       call ini_vfe(plev)
       vetaf(0) = 0.0_r8
       do k = 1,plev
         vetaf(k) = hvcoord%etam(k)
       enddo
       vetaf(plev+1) = 1.0_r8
       if (masterproc) write(iulog,*) 'vfe operator type:vfea(deboor algorithm) '
!Define the knots: 0,0,0,{0,.......,1},1,1,1 (in case of order=4)
       do k = 1,basis_degree
         knots(k) = 0.0_r8
       enddo
       do k = basis_degree+1,nlevbc+basis_degree
         knots(k) = vetaf(k-basis_degree-1)
       enddo
       do k = nlevbc+1+basis_degree,nlevbc+2*basis_degree
         knots(k) = 1.0_r8
       enddo
       if (masterproc) then
         write(iulog,'(a,i4) ') ' b-spline order = ',basis_order
         write(iulog,'(a,i4) ') ' number of basis functions = ',nbasis
! WRITE(iulog,'(a)'   ) '  Knot sequence:'
! DO k = 1,nlevbc+2*basis_degree
!    WRITE(iulog,'(I8,ES20.10)') k,knots(k)
! END DO
       endif
       vopinte(:,:) = 0.0
       vopderi(:,:) = 0.0
       if (vfebasis_order/= 1) stop '(stop) vfe linear basis is only working now... '
       call sub_deri_operator(debug_vfe,vetaf(0:plev+1),knots,vopderi)
       call sub_inte_operator(debug_vfe,vetaf(0:plev+1),knots,vopinte)
! Assign corrected vectors
       call vfe_mass_correction(hvcoord%etai(:),hvcoord%hyai(:),hvcoord%hybi(:),&
hvcoord%hyam(:),hvcoord%hybm(:),hvcoord%dhyam(:),hvcoord%dhybm(:))
! IF (masterproc) THEN
!    WRITE(*,*)  " VFE-hvcoord%(A,B,dA,dB)="
!    DO k = 1,plev
!       WRITE(*,'(I6,4ES20.10)') k,hvcoord%hyam(k), hvcoord%hybm(k),&
!                                   hvcoord%dhyam(k),hvcoord%dhybm(k)
!    END DO
! END IF
   endif !(on_vfe)
! ==================================
! Set layer locations
! ==================================
! ===============================================================
! Interfaces. Set nprlev to the interface above,the first time a
! nonzero surface pressure contribution is found. "nprlev"
! identifies the lowest pure pressure interface.
! ===============================================================
   do k = 1,plev
     if (hvcoord%nprlev==0.and.hvcoord%hybi(k).ne.0.0) hvcoord%nprlev = k-1
   enddo
!
! Set nprlev if no nonzero b's have been found. All interfaces are
! pure pressure. A pure pressure model requires other changes as well.
!
   if (hvcoord%nprlev==0) hvcoord%nprlev = plev+2
! =======================================
! Set delta sigma part of layer thickness
! pressures
! =======================================
   do k = 1,plev
     hvcoord%hybd(k) = hvcoord%hybi(k+1)-hvcoord%hybi(k)
   enddo
! ====================================================
! Calculate the log pressure extrapolation factor
! ====================================================
! JPE just prevents a compiler warning when plev==1
#if (PLEV>1)
   hvcoord%prsfac = log(hvcoord%hyam(plev)+hvcoord%hybm(plev))/&
log((hvcoord%hyam(plev)+hvcoord%hybm(plev))/(hvcoord%hyam(plev-1)+hvcoord%hybm(plev-1)))
#endif
! ======================================================================
! Test that A's and B's at full levels are arithmetic means of A's and
! B's at interfaces
! ======================================================================
   do k = 1,plev
     amean =(hvcoord%hyai(k+1)+hvcoord%hyai(k))*0.5d0
     bmean =(hvcoord%hybi(k+1)+hvcoord%hybi(k))*0.5d0
     if (amean==0..and.hvcoord%hyam(k)==0.) then
       atest = 0.
     else
       atest = abs(amean-hvcoord%hyam(k))/(0.5d0*(abs(amean+hvcoord%hyam(k))))
     endif
     if (bmean==0..and.hvcoord%hybm(k)==0.) then
       btest = 0.
     else
       btest = abs(bmean-hvcoord%hybm(k))/(0.5d0*(abs(bmean+hvcoord%hybm(k))))
     endif
     if (atest>eps) then
       if (masterproc) then
         write(iulog,9850)
         write(iulog,*) 'k,atest,eps = ',k,atest,eps
       endif
!        call endrun
     endif
     if (btest>eps) then
       if (masterproc) then
         write(iulog,9850)
         write(iulog,*) 'k,btest,eps = ',k,btest,eps
       endif
!        call endrun
     endif
   enddo
#endif
#if defined(CORE_SW)
!-- Vertical coordinate setting -----------------
   do k = 1,nlev+1 ! bottom to top
     hvcoord%znw(k) = hvcoord%hyai(k)+hvcoord%hybi(k)
!for physics
     si1d(k) = hvcoord%znw(k)
   enddo
   do k = 1,nlev ! bottom to top
     hvcoord%dnw(k) = hvcoord%znw(k+1)-hvcoord%znw(k)
     hvcoord%rdnw(k) = 1.0_r8/hvcoord%dnw(k)
     hvcoord%znu(k) = hvcoord%hyam(k)+hvcoord%hybm(k)
!for physics
     sl1d(k) = hvcoord%znu(k)
   enddo
   do k = 2,nlev
     hvcoord%dn(k) = 0.5_r8*(hvcoord%dnw(k)+hvcoord%dnw(k-1))
     hvcoord%rdn(k) = 1.0_r8/hvcoord%dn(k)
     hvcoord%fnp(k) = 0.5_r8*hvcoord%dnw(k)/hvcoord%dn(k)
     hvcoord%fnm(k) = 0.5_r8*hvcoord%dnw(k-1)/hvcoord%dn(k)
   enddo
#if defined(vcoord_sigma)
! sjchoi. 20150806
! this flag is for removing computational error.
! Ultimately it should be the same as that of the below section.
   do k = 1,nlev
     hvcoord%hyai_e(k) = 0._r8
     hvcoord%hybi_e(k) = 1._r8
   enddo
   do k = 2,nlev !!caution!! from 2
     hvcoord%hyam_e(k) = 0._r8
     hvcoord%hybm_e(k) = 1._r8
   enddo
#else
   do k = 1,nlev
     hvcoord%hyai_e(k) = hvcoord%rdnw(k)*(hvcoord%hyai(k+1)-hvcoord%hyai(k))
     hvcoord%hybi_e(k) = hvcoord%rdnw(k)*(hvcoord%hybi(k+1)-hvcoord%hybi(k))
   enddo
   do k = 2,nlev !!caution!! from 2
     hvcoord%hyam_e(k) = hvcoord%rdn(k)*(hvcoord%hyam(k)-hvcoord%hyam(k-1))
     hvcoord%hybm_e(k) = hvcoord%rdn(k)*(hvcoord%hybm(k)-hvcoord%hybm(k-1))
   enddo
#endif
!------------------------------------------------
!-- Extrapolation coefficient setting -----------
   cof1 =(2._r8*hvcoord%dn(2)+hvcoord%dn(3)) &
                                     /(hvcoord%dn(2)+hvcoord%dn(3)) &
                                       *hvcoord%dnw(1)/hvcoord%dn(2)
   cof2 = hvcoord%dn(2) &
                                     /(hvcoord%dn(2)+hvcoord%dn(3)) &
                                       *hvcoord%dnw(1)/hvcoord%dn(3)
   hvcoord%cf1 = hvcoord%fnp(2)+cof1
   hvcoord%cf2 = hvcoord%fnm(2)-cof1-cof2
   hvcoord%cf3 = cof2
! 외분점 구함....
   hvcoord%cfn =(0.5_r8*hvcoord%dnw(nlev)+hvcoord%dn(nlev))/hvcoord%dn(nlev)
   hvcoord%cfn1 =-0.5_r8*hvcoord%dnw(nlev)/hvcoord%dn(nlev)
!------------------------------------------------
#endif
   if (masterproc) then
     if (lprint) then
#if defined(CORE_SW)
       write(iulog,'(a) ') '================================================'
       write(iulog,'(a) ') 'vertical coord(bottom to top):hyai,hybi,etai*1000'
       write(iulog,'(a) ') ':hyam,hybm,etam*1000'
       write(iulog,'(a) ') ':hyai_e,hybi_e,etai_e*1000'
       write(iulog,'(a) ') ':hyam_e,hybm_e,etam_e*1000'
       write(iulog,'(a) ') '================================================'
       write(iulog,9800) 1,hvcoord%hyai(1),hvcoord%hybi(1),hvcoord%znw(1)
       write(iulog,9810) hvcoord%hyam(1),hvcoord%hybm(1),hvcoord%znu(1)
       write(iulog,9810) hvcoord%hyai_e(1),hvcoord%hybi_e(1),hvcoord%hyai_e(1)+hvcoord%hybi_e(1)
       do k = 2,nlev
         write(iulog,9800) k,hvcoord%hyai(k),hvcoord%hybi(k),hvcoord%znw(k)
         write(iulog,9810) hvcoord%hyam(k),hvcoord%hybm(k),hvcoord%znu(k)
         write(iulog,9810) hvcoord%hyai_e(k),hvcoord%hybi_e(k),hvcoord%hyai_e(k)+hvcoord%hybi_e(k)
         write(iulog,9810) hvcoord%hyam_e(k),hvcoord%hybm_e(k),hvcoord%hyam_e(k)+hvcoord%hybm_e(k)
       enddo
       write(iulog,9800) nlev+1,hvcoord%hyai(nlev+1),hvcoord%hybi(nlev+1),hvcoord%znw(nlev+1)
       open(unit = 5,file = "std.txt",form = "formatted",status = 'unknown')
       do k = 1,nlev
         write(5,*) k,-1.*hvcoord%dnw(k)
       enddo
       close(5)
       write(iulog,'(a) ') '================================================'
       write(iulog,'(a) ') 'p_top = '
       write(iulog,*) hvcoord%p_top
       write(iulog,'(a) ') '================================================'
#elif defined(CORE_SH)
       write(iulog,'(a) ') '================================================'
       write(iulog,'(a) ') 'vertical coord(bottom to top):hyai,hybi,eta*1000'
       write(iulog,'(a) ') '================================================'
       do k = 1,plev
         write(iulog,9800) k,hvcoord%hyai(k),hvcoord%hybi(k),hvcoord%hyai(k)+hvcoord%hybi(k)
         write(iulog,9810) hvcoord%hyam(k),hvcoord%hybm(k),hvcoord%hyam(k)+hvcoord%hybm(k)
       enddo
       write(iulog,9800) plevp,hvcoord%hyai(plevp),hvcoord%hybi(plevp),hvcoord%hyai(plevp)+hvcoord%hybi(plevp)
#endif
     else
! sanity check for endian problem with file:
       write(iulog,*) 'min/max hybm() coordinates:',minval(hvcoord%hybm(1:nlev)),maxval(hvcoord%hybm(1:nlev))
     endif
   endif
   9800 format(1x,i3,3p,3(f10.4,10x))
   9810 format(1x,3x,3p,4(10x,f10.4))
   9820 format(1x,i3,1x,f10.4)
   9850 format('HYCOEF:a and/or b vertical level coefficients at full',/,&
                 ' levels are not the arithmetic mean of half-level values')
!
   end function hvcoord_init
!==========================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module hybvcoord
!-------------------------------------------------------------------------------
