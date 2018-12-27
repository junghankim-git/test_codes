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

MODULE HybVCoord

  !```````````````````````````````````````````````````````
  USE KiapsBase,  ONLY : iulog => KIM_IU_LOG, r8 => KIM_REAL8_KIND
  USE Dimensions, ONLY : nlev, plev => nlev, plevp => nlevp
  USE PhysicalConstants, ONLY : p1000mb => KIM_P0, p0 => KIM_P0
  USE Control, ONLY : on_VFE, VFEbasis_order 
  USE VFE_MOD, ONLY : debug_vfe, basis_order, basis_degree, nlevbc, nbasis, nknot, &
                      VOPINTE, VOPDERI, INI_VFE, FIN_VFE, SUB_VERINT, SUB_VERDER, VFE_MASS_CORRECTION
  USE VFE_OPERATORS, ONLY : SUB_DERI_OPERATOR, SUB_INTE_OPERATOR
  !```````````````````````````````````````````````````````

  IMPLICIT NONE
  !-------------------------------------------------------

  PRIVATE


  REAL(r8), PUBLIC  :: p_top
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
  TYPE, PUBLIC :: hvcoord_t
      REAL(r8) ps0           ! base state sfc pressure for level definitions
      REAL(r8) hyai(nlev+1)  ! ps0 component of hybrid coordinate - interfaces
      REAL(r8) hyam(nlev)    ! ps0 component of hybrid coordinate - midpoints
      REAL(r8) hybi(nlev+1)  ! ps component of hybrid coordinate - interfaces
      REAL(r8) hybm(nlev)    ! ps component of hybrid coordinate - midpoints

      REAL(r8) znw(nlev+1)   ! etai = hyai + hybi
      REAL(r8) znu(nlev)     ! etam = hyam + hybm

      REAL(r8) dnw(nlev)
      REAL(r8) rdnw(nlev)

      REAL(r8) dn(nlev) 
      REAL(r8) rdn(nlev)

      REAL(r8) fnp(nlev)
      REAL(r8) fnm(nlev)

      REAL(r8) hyai_e(nlev) 
      REAL(r8) hybi_e(nlev) 
      REAL(r8) hyam_e(nlev) ! Level of Bottom(1) not in used.. 
      REAL(r8) hybm_e(nlev)

      REAL(r8) p_top 

      REAL(r8) cf1 
      REAL(r8) cf2 
      REAL(r8) cf3 
      REAL(r8) cfn 
      REAL(r8) cfn1
  END TYPE
#elif defined(CORE_SH)
!----------------------------------------------------------------------- 
! 
! Purpose: Hybrid level definitions: p = a*p0 + b*ps
!          interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
!          midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
! 
!-----------------------------------------------------------------------
  type, public :: hvcoord_t
      real(r8) ps0          ! base state sfc pressure for level definitions
      real(r8) hyai(plevp)  ! ps0 component of hybrid coordinate - interfaces
      real(r8) hyam(plev)   ! ps0 component of hybrid coordinate - midpoints
      real(r8) hybi(plevp)  ! ps component of hybrid coordinate - interfaces
      real(r8) hybm(plev)   ! ps component of hybrid coordinate - midpoints
      real(r8) hybd(plev)   ! difference  in b (hybi) across layers
      real(r8) dhyam(plev)  ! ps0 differential component of hybrid coordinate - midpoints
      real(r8) dhybm(plev)  ! ps differential component of hybrid coordinate - midpoints
      real(r8) prsfac       ! log pressure extrapolation factor (time, space independent)
      real(r8) etam(plev)   ! eta.  stored for conviencience
      real(r8) etai(plevp)  ! 
      integer  nprlev       ! number of pure pressure levels at top  
      integer  pad
  end type

#endif
  ! ---- for physics -----------------------------------
  REAL(KIND=r8), PUBLIC :: hyai1d(plevp), hybi1d(plevp), si1d(plevp), sl1d(plev) 
  ! ----------------------------------------------------

  PUBLIC :: hvcoord_init
  !-------------------------------------------------------

  
CONTAINS
!==========================================================

  ! =============================================================
  ! hvcoord_init: Initialize hybrid vertical coordinate system...
  ! returns: hv coordinate structure, ierr=0 means success
  ! 
  ! =============================================================

  FUNCTION hvcoord_init(hvfile_mid, hvfile_int, lprint, masterproc, ierr) result(hvcoord)

    !--------------------
    CHARACTER(LEN=*), INTENT(IN   ) :: hvfile_mid     ! file containing mid levels vertical coordinate system 
    CHARACTER(LEN=*), INTENT(IN   ) :: hvfile_int     ! file containing interface levels vertical coordinate system 
    LOGICAL         , INTENT(IN   ) :: lprint
    LOGICAL         , INTENT(IN   ) :: masterproc
    INTEGER         , INTENT(  OUT) :: ierr
    TYPE (hvcoord_t)                :: hvcoord

    !local
    INTEGER   :: k
    INTEGER   :: ierr11,ierr12, plevp_in, plev_in,ln
    REAL(r8)  :: hyai_loc(nlev+1)
    REAL(r8)  :: hybi_loc(nlev+1)
    REAL(r8)  :: hyam_loc(nlev)
    REAL(r8)  :: hybm_loc(nlev)
    REAL(r8)  :: cof1, cof2

    REAL(r8)  :: ci(nlev+1)
    REAL(r8)  :: del(nlev)
#if defined(CORE_SH)
    real(r8) amean,bmean,atest,btest,eps
#endif

    REAL(r8),DIMENSION(0:nlev+1) :: VETAF
    REAL(r8),DIMENSION(nknot)    :: knots
    !--------------------

    ierr=0

    hvcoord%ps0    = p1000mb       ! Base state surface pressure [Pa]

!sjchoi  ---------------------------------------------------------------------
#if defined(CORE_SW)
    hvcoord%p_top  = p_top          ! [Pa] Default

#if defined(IGW)  /* idealized test */
    hvcoord%p_top  = 27358.273_r8   ! [Pa] For IGW ~10km top
#elif defined(ASP)
    hvcoord%p_top  = 0._r8          ! [Pa] For ASP 2008 
#elif defined(DCMIP_case1X) 
    hvcoord%p_top  = 25494.4_r8     ! [Pa] For DCMIP 11,12,13
#elif defined(DCMIP_case20)
    hvcoord%p_top  = 20544.8_r8     ! [Pa] For DCMIP 20
#elif defined(DCMIP_case21) || defined(DCMIP_case22)
    hvcoord%p_top  = 3281.8_r8      ! [Pa] For DCMIP 21, 22
#elif defined(DCMIP_case31)
    hvcoord%p_top  = 27391.9_r8     ! [Pa] For DCMIP 31
#endif /* idealized test */

#endif
!-----------------------------------------------------------------------------

    ln = LEN(TRIM(hvfile_int))

    ! check if file name ends with .ascii
    IF ( hvfile_int(ln-4:ln) == 'ascii' ) THEN
       OPEN(UNIT=11,FILE=hvfile_int, form='formatted',status='old',iostat=ierr11)
       IF (ierr11 /= 0) THEN
          WRITE(iulog,*)'open() error:', __FILE__,__LINE__,hvfile_int
          ierr=ierr11
       END IF
       
       IF(ierr==0) THEN
          READ(11,*) plevp_in
          IF (plevp_in .ne. plevp) THEN
             WRITE(iulog,*) 'Error: hyai input file and HOMME nlev+1 do not match',plevp,plevp_in
             ierr=1
          END IF
          READ(11,*) hyai_loc(1:plevp)
          READ(11,*) plevp_in
          IF (plevp_in .ne. nlev+1) THEN
             WRITE(iulog,*) 'Error: hybi input file and HOMME nlev+1 do not match',plevp,plevp_in
             ierr=1
          END IF
          READ(11,*) hybi_loc(1:plevp)
          CLOSE(11)

          READ(12,*) plev_in
          IF (plev_in .ne. plev) THEN
             WRITE(iulog,*) 'Error: hyam input file and HOMME plev do not match',plev,plev_in
             ierr=1
          END IF
          READ(12,*) hyam_loc(1:plev)
          READ(12,*) plev_in
          IF (plev_in .ne. plev) THEN
             WRITE(iulog,*) 'Error: hybi input file and HOMME plev do not match',plev,plev_in
             ierr=1
          END IF
          READ(12,*) hybm_loc(1:plev)
          CLOSE(12)

#if defined(CORE_SH)
          DO k=1, nlev+1  ! TOP TO BOTTOM
             hvcoord%hyai(k) = hyai_loc(k)
             hvcoord%hybi(k) = hybi_loc(k)
             ! -------- for physics --------------------------
             hyai1d(k) = hyai_loc(nlev+2-k) ! BOTTOM TO TOP in PHYSICS
             hybi1d(k) = hybi_loc(nlev+2-k) ! BOTTOM TO TOP in PHYSICS
             ! ------------------------------------------------------
          END DO

          DO k=1, nlev    ! TOP TO BOTTOM
             hvcoord%hyam(k) = hyam_loc(k)
             hvcoord%hybm(k) = hybm_loc(k)
          END DO
#endif
#if defined(CORE_SW)
#if defined(vcoord_sigma)
          IF (masterproc) THEN
            WRITE(iulog,'(a)')'================================================'
            WRITE(iulog,'(a)')'YOU are USING SIGMA vertical coord.'
            WRITE(iulog,'(a)')'================================================'
          END IF
          DO k=1, nlev+1  ! BOTTOM TO TOP
             hvcoord%hyai(k) = 0._r8
             hvcoord%hybi(k) = hyai_loc(nlev+1-k+1) + hybi_loc(nlev+1-k+1)
          END DO
    
          DO k=1, nlev    ! BOTTOM TO TOP
             hvcoord%hyam(k) = 0._r8
             hvcoord%hybm(k) = .5_r8*(hvcoord%hyai(k+1)+hvcoord%hyai(k)) + .5_r8*(hvcoord%hybi(k+1)+hvcoord%hybi(k))
          END DO
#elif defined(vcoord_hybridsigma)
          IF (masterproc) THEN
            WRITE(iulog,'(a)')'================================================'
            WRITE(iulog,'(a)')'YOU are USING HYBRID-sigma vertical coord.'
            WRITE(iulog,'(a)')'================================================'
          END IF
          DO k=1, nlev+1  ! BOTTOM TO TOP
             hvcoord%hyai(k) = hyai_loc(nlev+1-k+1)
             hvcoord%hybi(k) = hybi_loc(nlev+1-k+1)
          END DO
    
          DO k=1, nlev    ! BOTTOM TO TOP
             hvcoord%hyam(k) = .5_r8*(hvcoord%hyai(k+1)+hvcoord%hyai(k))
             hvcoord%hybm(k) = .5_r8*(hvcoord%hybi(k+1)+hvcoord%hybi(k))
          END DO
#endif /*  -- vcoord_sigma */
#endif
       END IF

    ELSE
       OPEN(UNIT=11,FILE=hvfile_int, form='unformatted',status='old',iostat=ierr11)
       IF (ierr11 /= 0) THEN
          WRITE(iulog,*)'open() error:', __FILE__,__LINE__,hvfile_int
          ierr=ierr11
       END IF

       OPEN(UNIT=12,FILE=hvfile_mid, form='unformatted',access='sequential',status='old',iostat=ierr12)
       IF (ierr12 /= 0) THEN
          WRITE(iulog,*)'open() error:', __FILE__,__LINE__,hvfile_mid
          ierr=ierr12
       END IF


       IF(ierr==0) THEN	
          DO k=1, nlev+1
             READ(11) hyai_loc(k)
             READ(11) hybi_loc(k)
          END DO
          CLOSE(11)
          
          DO k=1, nlev
             READ(12) hyam_loc(k)
             READ(12) hybm_loc(k)
          END DO
          CLOSE(12)
          
#if defined(CORE_SH)
          DO k=1, nlev+1  ! TOP TO BOTTOM
             hvcoord%hyai(k) = hyai_loc(k)
             hvcoord%hybi(k) = hybi_loc(k)
          END DO

          DO k=1, nlev    ! TOP TO BOTTOM
             hvcoord%hyam(k) = hyam_loc(k)
             hvcoord%hybm(k) = hybm_loc(k)
          END DO
#endif
#if defined(CORE_SW)
#if defined(vcoord_sigma)
          IF (masterproc) THEN
            WRITE(iulog,'(a)')'================================================'
            WRITE(iulog,'(a)')'YOU are USING SIGMA vertical coord.'
            WRITE(iulog,'(a)')'================================================'
          END IF
          DO k=1, nlev+1  ! BOTTOM TO TOP
             hvcoord%hyai(k) = 0._r8
             hvcoord%hybi(k) = hyai_loc(nlev+1-k+1) + hybi_loc(nlev+1-k+1)
          END DO

          DO k=1, nlev    ! BOTTOM TO TOP
             hvcoord%hyam(k) = 0._r8
             hvcoord%hybm(k) = .5_r8*(hvcoord%hyai(k+1)+hvcoord%hyai(k)) + .5_r8*(hvcoord%hybi(k+1)+hvcoord%hybi(k))
          END DO
#elif defined(vcoord_hybridsigma)
          IF (masterproc) THEN
            WRITE(iulog,'(a)')'================================================'
            WRITE(iulog,'(a)')'YOU are USING HYBRID-sigma vertical coord.'
            WRITE(iulog,'(a)')'================================================'
          END IF
          DO k=1, nlev+1  ! BOTTOM TO TOP
             hvcoord%hyai(k) = hyai_loc(nlev+1-k+1)
             hvcoord%hybi(k) = hybi_loc(nlev+1-k+1)
          END DO

          DO k=1, nlev    ! BOTTOM TO TOP
             hvcoord%hyam(k) = .5_r8*(hvcoord%hyai(k+1)+hvcoord%hyai(k))
             hvcoord%hybm(k) = .5_r8*(hvcoord%hybi(k+1)+hvcoord%hybi(k))
          END DO
#endif /*  -- vcoord_sigma */
#endif
       END IF
    ENDIF

    IF(ierr>0) RETURN 

#if defined(CORE_SH)
  eps            = 1.D-05
  hvcoord%nprlev = 0

  do k=1,plev
     hvcoord%etam(k)=hvcoord%hyam(k)+hvcoord%hybm(k)
     ! ------ for physics --------------------
     sl1d(k)=hvcoord%hyam(plev+1-k)+hvcoord%hybm(plev+1-k)   ! BOTTOM TO TOP in PHYSICS
     ! ----------------------------------------
  end do
  do k=1,plevp
     hvcoord%etai(k)=hvcoord%hyai(k)+hvcoord%hybi(k)
     ! ------ for physics --------------------
     si1d(k)=hvcoord%hyai(plevp+1-k)+hvcoord%hybi(plevp+1-k) ! BOTTOM TO TOP in PHYSICS
     ! ----------------------------------------
  enddo

  ! ====================================
  ! Set vertical operator matrix for VFE
  ! ====================================
  IF (on_VFE) THEN

     IF (masterproc) THEN
        WRITE(iulog,'(a)')'================================================'
        WRITE(iulog,'(a)')'Set VFE coordinate' 
        WRITE(iulog,'(a)')'================================================'
     END IF

     CALL INI_VFE(plev)

     VETAF(0) = 0.0_r8
     DO k = 1, plev
        VETAF(k) = hvcoord%etam(k)
     END DO
     VETAF(plev+1) = 1.0_r8

     IF (masterproc) WRITE(iulog,*) 'vfe operator type : VFEA (Deboor algorithm)'

     !Define the knots: 0, 0, 0, {0, ......., 1}, 1, 1, 1 (in case of order=4)
     DO k = 1, basis_degree
        knots(k) = 0.0_r8
     END DO
     DO k = basis_degree+1, nlevbc+basis_degree
        knots(k) = VETAF(k-basis_degree-1)
     END DO
     DO k = nlevbc+1+basis_degree, nlevbc+2*basis_degree
        knots(k) = 1.0_r8
     END DO

     IF (masterproc) THEN
        WRITE(iulog, '(a,i4)') '  B-spline order = ', basis_order
        WRITE(iulog, '(a,i4)') '  Number of basis functions = ', nbasis
      ! WRITE(iulog, '(a)'   ) '  Knot sequence:'
      ! DO k = 1, nlevbc+2*basis_degree
      !    WRITE(iulog, '(I8,ES20.10)') k, knots(k)
      ! END DO
     END IF

     VOPINTE(:,:)=0.0
     VOPDERI(:,:)=0.0

     IF(VFEbasis_order/=1) STOP ' (stop) VFE linear basis is only working now... '
     CALL SUB_DERI_OPERATOR(debug_vfe, VETAF(0:plev+1), knots, VOPDERI)
     CALL SUB_INTE_OPERATOR(debug_vfe, VETAF(0:plev+1), knots, VOPINTE)

     ! Assign corrected vectors
     CALL VFE_MASS_CORRECTION( hvcoord%etai(:),hvcoord%hyai(:),hvcoord%hybi(:), &
                               hvcoord%hyam(:),hvcoord%hybm(:),hvcoord%dhyam(:),hvcoord%dhybm(:) )
 
   ! IF (masterproc) THEN 
   !    WRITE(*,*)  " VFE-hvcoord%(A,B,dA,dB)="
   !    DO k = 1, plev
   !       WRITE(*,'(I6,4ES20.10)') k, hvcoord%hyam(k),  hvcoord%hybm(k), &
   !                                   hvcoord%dhyam(k), hvcoord%dhybm(k)
   !    END DO
   ! END IF

  END IF !(on_VFE)


  ! ==================================
  ! Set layer locations
  ! ==================================

  ! ===============================================================
  ! Interfaces. Set nprlev to the interface above, the first time a 
  ! nonzero surface pressure contribution is found. "nprlev" 
  ! identifies the lowest pure pressure interface.
  ! ===============================================================

  do k=1,plev
     if (hvcoord%nprlev==0 .and. hvcoord%hybi(k).ne.0.0) hvcoord%nprlev = k - 1
  end do
!
! Set nprlev if no nonzero b's have been found. All interfaces are 
! pure pressure. A pure pressure model requires other changes as well. 
!
  if (hvcoord%nprlev==0) hvcoord%nprlev = plev + 2

  ! =======================================
  ! Set delta sigma part of layer thickness
  ! pressures
  ! =======================================

  do k=1,plev
     hvcoord%hybd(k) = hvcoord%hybi(k+1) - hvcoord%hybi(k)
  end do

  ! ====================================================
  ! Calculate the log pressure extrapolation factor
  ! ====================================================
! JPE just prevents a compiler warning when plev==1
#if (PLEV>1)
  hvcoord%prsfac = log( hvcoord%hyam(plev  ) + hvcoord%hybm(plev)) / &
                   log((hvcoord%hyam(plev  ) + hvcoord%hybm(plev)) / (hvcoord%hyam(plev-1) + hvcoord%hybm(plev-1)))
#endif
  ! ======================================================================
  ! Test that A's and B's at full levels are arithmetic means of A's and
  ! B's at interfaces
  ! ======================================================================

  do k = 1,plev
     amean = ( hvcoord%hyai(k+1) + hvcoord%hyai(k) )*0.5D0
     bmean = ( hvcoord%hybi(k+1) + hvcoord%hybi(k) )*0.5D0
     if(amean == 0. .and. hvcoord%hyam(k) == 0.) then
        atest = 0.
     else
        atest = abs( amean - hvcoord%hyam(k) )/ ( 0.5D0*( abs(amean + hvcoord%hyam(k)) ) )
     endif
     if(bmean == 0. .and. hvcoord%hybm(k) == 0.) then
        btest = 0.
     else
        btest = abs( bmean - hvcoord%hybm(k) )/ ( 0.5D0*( abs(bmean + hvcoord%hybm(k)) ) )
     endif
     if (atest > eps) then
        if (masterproc) then
           write(iulog,9850)
           write(iulog,*)'k,atest,eps=',k,atest,eps
        end if
!        call endrun
     endif

     if (btest > eps) then
        if (masterproc) then
           write(iulog,9850)
           write(iulog,*)'k,btest,eps=',k,btest,eps
        end if
!        call endrun
     endif
  end do
#endif
#if defined(CORE_SW)
    !-- Vertical coordinate setting -----------------
    DO k=1, nlev+1      ! BOTTOM TO TOP
       hvcoord%znw(k) = hvcoord%hyai(k)+hvcoord%hybi(k)
       !for physics
       si1d(k)        = hvcoord%znw(k)
    END DO
    DO k=1, nlev        ! BOTTOM TO TOP
       hvcoord%dnw(k)  = hvcoord%znw(k+1) - hvcoord%znw(k)
       hvcoord%rdnw(k) = 1.0_r8/hvcoord%dnw(k)
       hvcoord%znu(k)  = hvcoord%hyam(k) + hvcoord%hybm(k)
       !for physics
       sl1d(k)        = hvcoord%znu(k)
    END DO
    DO k=2, nlev
       hvcoord%dn(k)  = 0.5_r8*(hvcoord%dnw(k)+hvcoord%dnw(k-1))
       hvcoord%rdn(k) = 1.0_r8/hvcoord%dn(k)
       hvcoord%fnp(k) = 0.5_r8*hvcoord%dnw(k  )/hvcoord%dn(k)
       hvcoord%fnm(k) = 0.5_r8*hvcoord%dnw(k-1)/hvcoord%dn(k)
    END DO

#if defined(vcoord_sigma) 
    ! sjchoi. 20150806
    ! this flag is for removing computational error.
    ! Ultimately it should be the same as that of the below section.
    DO k=1, nlev
       hvcoord%hyai_e(k) = 0._r8
       hvcoord%hybi_e(k) = 1._r8
    END DO

    DO k=2, nlev !!CAUTION!! from 2 
       hvcoord%hyam_e(k) = 0._r8
       hvcoord%hybm_e(k) = 1._r8
    END DO
#else
    DO k=1, nlev
       hvcoord%hyai_e(k) = hvcoord%rdnw(k)*(hvcoord%hyai(k+1) - hvcoord%hyai(k))
       hvcoord%hybi_e(k) = hvcoord%rdnw(k)*(hvcoord%hybi(k+1) - hvcoord%hybi(k))
    END DO

    DO k=2, nlev !!CAUTION!! from 2 
       hvcoord%hyam_e(k) = hvcoord%rdn(k)*(hvcoord%hyam(k) - hvcoord%hyam(k-1))
       hvcoord%hybm_e(k) = hvcoord%rdn(k)*(hvcoord%hybm(k) - hvcoord%hybm(k-1))
    END DO
#endif
    !------------------------------------------------

    !-- Extrapolation coefficient setting -----------
    cof1 = (2._r8*hvcoord%dn(2)+hvcoord%dn(3)) &
           /(hvcoord%dn(2)+hvcoord%dn(3)) &
           *hvcoord%dnw(1)/hvcoord%dn(2)
    cof2 =     hvcoord%dn(2) &
           /(hvcoord%dn(2)+hvcoord%dn(3)) &
           *hvcoord%dnw(1)/hvcoord%dn(3)

    hvcoord%cf1  = hvcoord%fnp(2) + cof1
    hvcoord%cf2  = hvcoord%fnm(2) - cof1 - cof2
    hvcoord%cf3  = cof2

    ! 외분점 구함....
    hvcoord%cfn  = (0.5_r8*hvcoord%dnw(nlev )+hvcoord%dn(nlev ))/hvcoord%dn(nlev )
    hvcoord%cfn1 = -0.5_r8*hvcoord%dnw(nlev )/hvcoord%dn(nlev )
    !------------------------------------------------
#endif

    IF (masterproc) THEN
       IF (lprint) THEN
#if defined(CORE_SW)
          WRITE(iulog,'(a)')'================================================'
          WRITE(iulog,'(a)')'vertical coord (Bottom to Top) : hyai, hybi, etai *1000'
          WRITE(iulog,'(a)')'                               : hyam, hybm, etam *1000'
          WRITE(iulog,'(a)')'                               : hyai_e, hybi_e, etai_e *1000'
          WRITE(iulog,'(a)')'                               : hyam_e, hybm_e, etam_e *1000'
          WRITE(iulog,'(a)')'================================================'

             WRITE(iulog,9800)1,hvcoord%hyai(1),hvcoord%hybi(1),hvcoord%znw(1)
             WRITE(iulog,9810)  hvcoord%hyam(1),hvcoord%hybm(1),hvcoord%znu(1)
             WRITE(iulog,9810)  hvcoord%hyai_e(1), hvcoord%hybi_e(1), hvcoord%hyai_e(1)+hvcoord%hybi_e(1)
          DO k=2, nlev
             WRITE(iulog,9800)k,hvcoord%hyai(k),hvcoord%hybi(k),hvcoord%znw(k)
             WRITE(iulog,9810)  hvcoord%hyam(k),hvcoord%hybm(k),hvcoord%znu(k)
             WRITE(iulog,9810)  hvcoord%hyai_e(k), hvcoord%hybi_e(k), hvcoord%hyai_e(k)+hvcoord%hybi_e(k)
             WRITE(iulog,9810)  hvcoord%hyam_e(k), hvcoord%hybm_e(k), hvcoord%hyam_e(k)+hvcoord%hybm_e(k)
          END DO
          WRITE(iulog,9800)nlev+1,hvcoord%hyai(nlev+1),hvcoord%hybi(nlev+1),hvcoord%znw(nlev+1)

          OPEN (unit=5, file="STD.txt", form="formatted", status='unknown')
          DO k=1, nlev
             WRITE(5,*) k, -1.*hvcoord%dnw(k)
          END DO
          CLOSE(5)

          WRITE(iulog,'(a)')'================================================'
          WRITE(iulog,'(a)')'p_top = '
          WRITE(iulog,*) hvcoord%p_top
          WRITE(iulog,'(a)')'================================================'
#elif defined(CORE_SH)
          WRITE(iulog,'(a)')'================================================'
          WRITE(iulog,'(a)')'vertical coord (Bottom to Top) : hyai, hybi, eta *1000'
          WRITE(iulog,'(a)')'================================================'
          do k=1,plev
             write(iulog,9800)k,hvcoord%hyai(k),hvcoord%hybi(k),hvcoord%hyai(k)+hvcoord%hybi(k)
             write(iulog,9810) hvcoord%hyam(k), hvcoord%hybm(k), hvcoord%hyam(k)+hvcoord%hybm(k)
          end do
          write(iulog,9800)plevp,hvcoord%hyai(plevp),hvcoord%hybi(plevp),hvcoord%hyai(plevp)+hvcoord%hybi(plevp)
#endif
       ELSE
          ! sanity check for endian problem with file:
          WRITE(iulog,*)'min/max hybm() coordinates: ',minval(hvcoord%hybm(1:nlev)),maxval(hvcoord%hybm(1:nlev))
       END IF
    END IF
  
9800 format( 1x, i3, 3p, 3(f10.4,10x) )
9810 format( 1x, 3x, 3p, 4(10x,f10.4) )
9820 format( 1x, i3, 1x, f10.4 )
9850 format('HYCOEF: A and/or B vertical level coefficients at full',/, &
            ' levels are not the arithmetic mean of half-level values')

  END FUNCTION hvcoord_init
!==========================================================

END MODULE HybVCoord
