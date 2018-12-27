!===============================================================================
!  MODULE CostFunMinMod
!
!> @brief
!> - Modules for minimization of cost functions
!
!> @date 27MAY2013
!> - HJ SONG: Design and code
!> @date 12JUL2013
!> - HJ SONG: Insert parameter transformation modeuls and clean up
!> @date 30MAY2016
!> - HJ SONG: Hybrid 4D Ensemble Variational Code Written referring to
!>            Lorenc (2003), Lorenc et al. (2015), and Kuhl et al. (2013)
!===============================================================================
MODULE CostFunMinMod

!===============================================================================
! A. Common modules and variables
!===============================================================================
! 1. Common modules
  USE KiapsBase,      ONLY : int_kind  => KIM_int_kind,   &
                             real_kind => KIM_real8_kind, &
                             iulog     => KIM_iu_log
  USE Dimensions,     ONLY : np, nlev,                 &
                             nelemd, nets, nete,       &
                             nvar3d, nvar2d,           &
                             nlevar,                   &
                             zwn2,                     &
                             zwn2_a,                   & ! for alpha
                             neig,                     &
                             neig_chi,                 &
                             neig_Tps,                 &
                             neig_q,                   &
                             zwn_nproc, zwn_iproc,     &
                             zwn_nproc_a, zwn_iproc_a, & ! for alpha
                             ntime,                    & ! # of time bins
                             nt_itime                    ! itime for nature
  USE KiapsParallel,  ONLY : par    => KIM_par,    &
                             hybrid => KIM_hybrid, &
                             global_shared_buf,    &
                             global_shared_sum,    &
                             par_double_precision, &
                             par_character,        &
                             par_sum,              &
                             par_max,par_min
  USE KiapsParallel,  ONLY : BarrierPar
  USE KiapsParMPI,    ONLY : Par_Reduce
  USE KiapsWave,      ONLY : zwn_blength,    &
                             zwn_bstart,     &
                             zwn_blength_a,  & ! for alpha
                             zwn_bstart_a      ! for alpha
  USE GlobalNorms,    ONLY : Wrap_Repro_Sum
  USE Timing,         ONLY : TimingStart, TimingStop
  USE BackgroundMod,  ONLY : TempToM,               &
                             bg_md, nt_org,         &
                             bg_uv_rot,             &
                             !up
                             !errvar,                &
                             errvar_up,             &
                             bg_eig,                &
                             bcov_sp,               & ! eigenvectors of B for psi
                             eigv,                  & ! eigenvalues of B for psi
                             bcov_sp_chi,           & ! eigenvectors of B for chi
                             eigv_chi,              & ! eigenvalues of B for chi
                             bcov_sp_Tps,           & ! eigenvectors of B for (T,ps)
                             eigv_Tps,              & ! eigenvalues of B for (T,ps)
                             bcov_sp_q,             & ! eigenvectors of B for q
                             eigv_q,                & ! eigenvalues of B for q
                             bcov_sp_a,             & ! eigenvectors of B for alpha
                             eigv_a,                & ! eigenvalues of B for alpha
                             nsmpl,                 & 
                             cv_opt_hum,            &
                             is_moist_mr,           &
                             exst_nt,               &
                             bg_smpl,               &
                             b_clim_v,              &
                             b_ens_v
  USE ObservationMod, ONLY : obsmet,          &
                             sdmiss,          &
                             obcheck,         &
                             obcheck_outer,   &
                             crit_omb_amsua,  &
                             crit_omb_gpsro,  &
                             crit_omb_iasi,   &
                             crit_omb_cris,   &
                             crit_omb_atms,   &
                             crit_omb_atmswv, &
                             crit_omb_mhs,    &
                             crit_omb_csr,    &
                             singleobs,       &
                             da_sonde,        &
                             da_surface,      &
                             da_bogus,        &
                             da_aircraft,     &
                             da_amv,          &
                             da_scatwind,     &
                             da_amsua,        &
                             da_gpsro,        &
                             da_iasi,         &
                             da_cris,         &
                             da_atms,         &
                             da_atmswv,       &
                             da_mhs,          &
                             da_csr,          &
                             RTTOVused,       &
                             ROPPused
  USE ObsAllocateMod, ONLY : obsbg_type, &
                             obs_type
  USE CalRTTOVMod,    ONLY : CalRTTOV
  USE CalROPPMod,     ONLY : CalROPP
  USE ReadObsMod,     ONLY : obs
  USE CtrlToSpecMod,  ONLY : RedDof
  USE SpecToCtrlMod,  ONLY : SpecToCtrl, &
                             AdjSpecToCtrl_opt
  USE CtrlToModelMod, ONLY : TlmCtrlToModel, AdjCtrlToModel, &
                             PhiToTemp
  USE ModelToObsMod,     ONLY : ModelToObs
  USE TlmModelToObsMod,  ONLY : TlmModelToObs
  USE AdjModelToObsMod,  ONLY : AdjModelToObs
  USE MakebMod,       ONLY : MakebSingle,                  &
                             MakebSonde,    MakebSurface,  &
                             MakebAircraft, MakebAmv,      &
                             MakebScatwind,                &
                             MakebAmsua,    MakebGpsro,    &
                             MakebIasi,     MakebCris,     &
                             MakebAtms,     MakebMhs,      &
                             MakebCsr,      MakebBogus,    &
                             MakebAtmswv
  USE MultiplyRMod,   ONLY : MultiplyRSingle,                      &
                             MultiplyRSonde,    MultiplyRSurface,  &
                             MultiplyRAircraft, MultiplyRAmv,      &
                             MultiplyRScatwind,                    &
                             MultiplyRAmsua,    MultiplyRGpsro,    &
                             MultiplyRIasi,     MultiplyRCris,     &
                             MultiplyRAtms,     MultiplyRMhs,      &
                             MultiplyRCsr,      MultiplyRBogus,    &
                             MultiplyRAtmswv
!  USE ModelToEigMod,  ONLY : TlmModelToEig
  USE KiapsGrid,     ONLY : nUPs_l,           &
                            ConvertUP2EP_hyb, &
                            ConvertEP2UP_hyb, &
                            ConvertEP2EP0_hyb

! 2. Common variables
  IMPLICIT NONE
  PRIVATE
!===============================================================================

!===============================================================================
! B. Public modules and variables
!===============================================================================
! 1. Public modules
  PUBLIC :: MinConGraSpec
!!  PUBLIC :: EigToModel

! 2. Public variables
  REAL(KIND=real_kind), DIMENSION(:,:,:,:,:), ALLOCATABLE, &
     PUBLIC :: an_md                                       ! Analysis
  REAL(KIND=real_kind), &
     PUBLIC :: cgtol1, cgtol2
  INTEGER(KIND=int_kind), &
     PUBLIC :: noutl,     &
               niter,     &
               optmethod
  TYPE(obsbg_type), DIMENSION(:), ALLOCATABLE, &  ! for time dimension
     PUBLIC :: gs_obs, ad_obs, tl_obs, bg_obs, an_obs
  
  TYPE(obs_type), DIMENSION(:), ALLOCATABLE, &  ! for time dimension
     PUBLIC :: obs_org  ! for observation check every outer-loop

  REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE, &
     PUBLIC :: x_eig_acc,                            &
               x_eig_acc_chi,                        &
               x_eig_acc_Tps,                        &
               x_eig_acc_q
  REAL(KIND=real_kind), DIMENSION(:,:,:), ALLOCATABLE, &
     PUBLIC :: x_eig_acc_a  ! for alpha-variable

  INTEGER(KIND=int_kind) :: ierr

!===============================================================================

CONTAINS

!===============================================================================
!  SUBROUTINE MinConGraSpec
!
!> @brief
!> - Minimization of cost function using conjugate gradient
!> - in spectral space
!
!> @date 12JUL2013
!> - HJ SONG: Design and code
!> @date 30MAY2016
!> - HJ SONG: Hybrid 4D Ensemble Variational Code Written
!
!> @param[out] an_md
!===============================================================================
  SUBROUTINE MinConGraSpec

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE Edge,          ONLY : edgebuffer_t,   &
                              InitEdgeBuffer, &
                              EdgeVpack,      &
                              EdgeVunpack
    USE Bndry,         ONLY : Bndry_ExchangeV
 
  ! 2. Variables

  ! 2-2. Local
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime) &
       :: gs_md,                                               &
          x_md
    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc) &
       :: x_eig
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc) &
       :: x_eig_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc) &
       :: x_eig_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc) &
       :: x_eig_q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl) &
       :: x_eig_a  ! for alpha-variable

    INTEGER(KIND=int_kind)      &
       :: i, j, k, kk, ie,      &
          wn, vr,               &
          iter, ioutl,          &
          iproc, ierr, err,     &
          itime 
  !=============================================================================

! Extended to 4D
!! ALLOCATE Analysis
    ALLOCATE(an_md(np,np,nlevar,nelemd,ntime))
    ALLOCATE(gs_obs(ntime), &
             tl_obs(ntime), &
             ad_obs(ntime), &
             bg_obs(ntime), &
             an_obs(ntime), &
             obs_org(ntime))

  DO itime = 1, ntime
!! ALLOCATE observation information
    IF( singleobs )THEN
        ALLOCATE(gs_obs(itime)%single(obsmet(itime)%nsg_e_max,nelemd))
        ALLOCATE(tl_obs(itime)%single(obsmet(itime)%nsg_e_max,nelemd))
        ALLOCATE(ad_obs(itime)%single(obsmet(itime)%nsg_e_max,nelemd))
        ALLOCATE(bg_obs(itime)%single(obsmet(itime)%nsg_e_max,nelemd))
        ALLOCATE(an_obs(itime)%single(obsmet(itime)%nsg_e_max,nelemd))
        ALLOCATE(obs_org(itime)%single(obsmet(itime)%nsd_e_max,nelemd))
        obs_org(itime)%single(:,:) = obs(itime)%single(:,:)
    END IF
    IF (da_sonde) THEN
        ALLOCATE(gs_obs(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
        ALLOCATE(tl_obs(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
        ALLOCATE(ad_obs(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
        ALLOCATE(bg_obs(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
        ALLOCATE(an_obs(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
        ALLOCATE(obs_org(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
        obs_org(itime)%sonde(:,:) = obs(itime)%sonde(:,:)
    END IF
    IF (da_surface) THEN
        ALLOCATE(gs_obs(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
        ALLOCATE(tl_obs(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
        ALLOCATE(ad_obs(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
        ALLOCATE(bg_obs(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
        ALLOCATE(an_obs(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
        ALLOCATE(obs_org(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
        obs_org(itime)%surface(:,:) = obs(itime)%surface(:,:)
    END IF
    IF (da_bogus) THEN
        ALLOCATE(gs_obs(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
        ALLOCATE(tl_obs(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
        ALLOCATE(ad_obs(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
        ALLOCATE(bg_obs(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
        ALLOCATE(an_obs(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
        ALLOCATE(obs_org(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
        obs_org(itime)%bogus(:,:) = obs(itime)%bogus(:,:)
    END IF
    IF (da_aircraft) THEN
        ALLOCATE(gs_obs(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
        ALLOCATE(tl_obs(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
        ALLOCATE(ad_obs(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
        ALLOCATE(bg_obs(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
        ALLOCATE(an_obs(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
        ALLOCATE(obs_org(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
        obs_org(itime)%aircraft(:,:) = obs(itime)%aircraft(:,:)
    END IF
    IF (da_amv) THEN
        ALLOCATE(gs_obs(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
        ALLOCATE(tl_obs(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
        ALLOCATE(ad_obs(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
        ALLOCATE(bg_obs(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
        ALLOCATE(an_obs(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
        ALLOCATE(obs_org(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
        obs_org(itime)%amv(:,:) = obs(itime)%amv(:,:)
    END IF
    IF (da_scatwind) THEN
        ALLOCATE(gs_obs(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
        ALLOCATE(tl_obs(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
        ALLOCATE(ad_obs(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
        ALLOCATE(bg_obs(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
        ALLOCATE(an_obs(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
        ALLOCATE(obs_org(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
        obs_org(itime)%scatwind(:,:) = obs(itime)%scatwind(:,:)
    END IF
    IF (da_amsua) THEN
        ALLOCATE(gs_obs(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
        ALLOCATE(tl_obs(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
        ALLOCATE(ad_obs(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
        ALLOCATE(bg_obs(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
        ALLOCATE(an_obs(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
        ALLOCATE(obs_org(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
        obs_org(itime)%amsua(:,:,:) = obs(itime)%amsua(:,:,:)
    END IF
    IF (da_gpsro) THEN
        ALLOCATE(gs_obs(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
        ALLOCATE(tl_obs(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
        ALLOCATE(ad_obs(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
        ALLOCATE(bg_obs(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
        ALLOCATE(an_obs(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
        ALLOCATE(obs_org(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
        obs_org(itime)%gpsro(:,:) = obs(itime)%gpsro(:,:)
    END IF
    IF (da_iasi) THEN
        ALLOCATE(gs_obs(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
        ALLOCATE(tl_obs(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
        ALLOCATE(ad_obs(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
        ALLOCATE(bg_obs(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
        ALLOCATE(an_obs(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
        ALLOCATE(obs_org(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
        obs_org(itime)%iasi(:,:,:) = obs(itime)%iasi(:,:,:)
    END IF
    IF (da_cris) THEN
        ALLOCATE(gs_obs(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
        ALLOCATE(tl_obs(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
        ALLOCATE(ad_obs(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
        ALLOCATE(bg_obs(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
        ALLOCATE(an_obs(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
        ALLOCATE(obs_org(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
        obs_org(itime)%cris(:,:,:) = obs(itime)%cris(:,:,:)
    END IF
    IF (da_atms) THEN
        ALLOCATE(gs_obs(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
        ALLOCATE(tl_obs(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
        ALLOCATE(ad_obs(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
        ALLOCATE(bg_obs(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
        ALLOCATE(an_obs(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
        ALLOCATE(obs_org(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
        obs_org(itime)%atms(:,:,:) = obs(itime)%atms(:,:,:)
    END IF
    IF (da_atmswv) THEN
        ALLOCATE(gs_obs(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
        ALLOCATE(tl_obs(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
        ALLOCATE(ad_obs(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
        ALLOCATE(bg_obs(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
        ALLOCATE(an_obs(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
        ALLOCATE(obs_org(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
        obs_org(itime)%atmswv(:,:,:) = obs(itime)%atmswv(:,:,:)
    END IF
    IF (da_mhs) THEN
        ALLOCATE(gs_obs(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
        ALLOCATE(tl_obs(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
        ALLOCATE(ad_obs(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
        ALLOCATE(bg_obs(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
        ALLOCATE(an_obs(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
        ALLOCATE(obs_org(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
        obs_org(itime)%mhs(:,:,:) = obs(itime)%mhs(:,:,:)
    END IF
    IF (da_csr) THEN
        ALLOCATE(gs_obs(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
        ALLOCATE(tl_obs(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
        ALLOCATE(ad_obs(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
        ALLOCATE(bg_obs(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
        ALLOCATE(an_obs(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
        ALLOCATE(obs_org(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
        obs_org(itime)%csr(:,:,:) = obs(itime)%csr(:,:,:)
    END IF
  END DO ! itime

  !=============================================================================
  ! OUTER LOOP
  !=============================================================================
    ALLOCATE(x_eig_acc(neig,zwn_nproc))
    ALLOCATE(x_eig_acc_chi(neig_chi,zwn_nproc))
    ALLOCATE(x_eig_acc_Tps(neig_Tps,zwn_nproc))
    ALLOCATE(x_eig_acc_q(neig_q,zwn_nproc))
    ALLOCATE(x_eig_acc_a(nlev+1,zwn_nproc_a,nsmpl))

    x_eig_acc = 0.D0
    x_eig_acc_chi = 0.D0
    x_eig_acc_Tps = 0.D0
    x_eig_acc_q = 0.D0
    x_eig_acc_a = 0.D0  ! alpha
    x_eig = 0.D0

    gs_md = bg_md

    x_eig_chi = 0.D0
    x_eig_Tps = 0.D0
    x_eig_q = 0.D0
    x_eig_a = 0.D0  ! alpha
    DO ioutl = 1, noutl
       IF (par%ismasterproc) print*,'ioutl:',ioutl
       CALL Mpi_Barrier(par%comm, ierr)
    !==========================================================================
  
    !==========================================================================
    ! A. RCG
    !==========================================================================
       IF(optmethod .eq. 1) CALL RCGmethod(ioutl,       &
                                           gs_md,       &
                                           x_eig,       &
                                           x_eig_chi,   &
                                           x_eig_Tps,   &
                                           x_eig_q,     &
                                           x_eig_a)   ! alpha

    !==========================================================================
    
    !==========================================================================
    ! B. PreConditioned PRCG
    !==========================================================================
!       IF(optmethod .eq. 2) CALL PRCGmethod(ioutl, gs_md, x_eig)
  
    !==========================================================================
    ! C. L-BFGS
    !==========================================================================
!       IF (optmethod .eq. 3) CALL LBFGSmethod(ioutl, gs_md, x_eig)


    !==========================================================================

       x_md = 0.D0
       ierr = TimingStart('TlmEigToModel_outer')   !--- profiling
       CALL TlmEigToModel(gs_md(:,:,:,:,:),       &
                          x_eig,     &
                          x_eig_chi, &
                          x_eig_Tps, &
                          x_eig_q,   &
                          x_eig_a,   &
                          x_md)
       ierr = TimingStop('TlmEigToModel_outer')   !--- profiling

       gs_md(:,:,:,:,:) = gs_md(:,:,:,:,:) + x_md(:,:,:,:,:)
     !==========================================================================

    END DO ! for outer loop

    DEALLOCATE(x_eig_acc)
    DEALLOCATE(x_eig_acc_chi)
    DEALLOCATE(x_eig_acc_Tps)
    DEALLOCATE(x_eig_acc_q)
    DEALLOCATE(x_eig_acc_a)
  !=============================================================================

  !=============================================================================
  ! Transform eigen space to model space
  !=============================================================================
    an_md = gs_md
  !=============================================================================

    IF( exst_nt )THEN

#define Check_RMSE_on
#ifdef Check_RMSE_on     
       ierr = TimingStart('Check_RMSE')   !--- profiling
         CALL Check_RMSE
       ierr = TimingStop('Check_RMSE')   !--- profiling
#endif
    ENDIF

#define Check_OB_on
#ifdef Check_OB_on     
       ierr = TimingStart('Check_OB')   !--- profiling
       DO itime = 1, ntime
         CALL Check_OB(itime)
       END DO
       ierr = TimingStop('Check_OB')   !--- profiling
#endif

    DEALLOCATE(gs_obs, &
               tl_obs, &
               ad_obs, &
               bg_obs, &
               an_obs, &
               obs_org)

  !=============================================================================
  ! Convert specific humidity to mixing ratio
  !=============================================================================
    IF (is_moist_mr) THEN
       DO itime = 1, ntime
       DO ie = 1, nelemd
       DO k  = 1, nlev
       DO j  = 1, np
       DO i  = 1, np
          bg_md(i,j,k+3*nlev,ie,itime) = bg_md(i,j,k+3*nlev,ie,itime)/   &
                                         (1.D0-bg_md(i,j,k+3*nlev,ie,itime))  ! to mixing ratio
          an_md(i,j,k+3*nlev,ie,itime) = an_md(i,j,k+3*nlev,ie,itime)/   &
                                         (1.D0-an_md(i,j,k+3*nlev,ie,itime))  ! to mixing ratio
       END DO
       END DO
       END DO
       END DO
       END DO
    END IF

    IF(par%ismasterproc) THEN
       print *, '======= Successfully Finished 3dVar Analysis ======='
    END IF

    RETURN

  END SUBROUTINE MinConGraSpec
!===============================================================================
  SUBROUTINE RCGmethod(ioutl,     &
                       gs_md,     &
                       x_eig,     &
                       x_eig_chi, &
                       x_eig_Tps, &
                       x_eig_q,   &
                       x_eig_a)  ! for alpha (HJS)

  ! 2. Variables
  ! 2-1. Output
   INTEGER(KIND=int_kind), INTENT(IN) :: ioutl

   REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime), &
       INTENT(IN)    :: gs_md

   REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc), &
       INTENT(INOUT) :: x_eig
   REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc), &
       INTENT(INOUT) :: x_eig_chi
   REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc), &
       INTENT(INOUT) :: x_eig_Tps
   REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc), &
       INTENT(INOUT) :: x_eig_q
   REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl), &
       INTENT(INOUT) :: x_eig_a
 
  ! 2-2. Local
    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc) &           ! for psi
       :: b_eig, r_eig, p_eig, Ap_eig
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc) &       ! for chi
       :: b_eig_chi, r_eig_chi, p_eig_chi, Ap_eig_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc) &       ! for (T,ps)
       :: b_eig_Tps, r_eig_Tps, p_eig_Tps, Ap_eig_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc) &         ! for q
       :: b_eig_q, r_eig_q, p_eig_q, Ap_eig_q
    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl) &  !-- alpha 
       :: b_eig_a, r_eig_a, p_eig_a, Ap_eig_a

    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc,niter) & 
       :: r_eig_history                                         ! for psi
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc,niter) &
       :: r_eig_chi_history                                     ! for chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc,niter) &
       :: r_eig_Tps_history                                     ! for (T,ps)
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc,niter) &
       :: r_eig_q_history                                       ! for q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl,niter) &
       :: r_eig_q_history_a                                       !-- alpha

    REAL(KIND=real_kind) &
       :: rsold, rsnew,  &
          rs0, pAP,      &                 ! rs at initial
          alpha,beta,    &
          j_b, j_o, tmp
    INTEGER(KIND=int_kind)  &
       :: i, j, k, kk, ie,  &
          wn, vr, iter,     &
          iproc, ierr, err, &
          itime,            &
          ismpl
    REAL(KIND=real_kind) &
         :: rjrk
    REAL(KIND=real_kind)  &
       :: tmpmax, tmpmin, &
          maxv, minv

    INTEGER(KIND=int_kind) &
         :: info, nlanczos, nl, nw
    REAL(KIND=real_kind), DIMENSION(niter) &
         :: Lanczos_diag! diagonal part of niter-dimensional Lanczos tridiagonal matrix 
    REAL(KIND=real_kind), DIMENSION(niter-1) &
         :: Lanczos_subd! sub-diagonal part of niter-dimensional Lanczos tridiagonal matrix 
    REAL(KIND=real_kind) &
         :: alphaold, betaold, gamma, eps, ortho
    REAL(KIND=real_kind), DIMENSION(niter, niter) &
         :: Lanczos_orth! sub-diagonal part of niter-dimensional Lanczos tridiagonal matrix 
    REAL(KIND=real_kind), DIMENSION(2*niter-2) &
         :: work
    REAL(KIND=real_kind), allocatable, DIMENSION(:,:,:) &
         :: h_eig! the matrix of orthonormal vectors that approximate the Hessian eigenvectors
    INTEGER(KIND=int_kind), allocatable, DIMENSION(:) &
         :: eig, eig1
    REAL(KIND=real_kind), allocatable, DIMENSION(:) &
         :: eigenval
    logical :: eig_loc(niter)


    real(kind=real_kind)   :: criteria                 ! Convergence criterion
    REAL(kind=real_kind)   :: cgtol
  !=============================================================================
   !! convergence criteria
    IF (ioutl .eq. 1) cgtol = cgtol1
    IF (ioutl .gt. 1) THEN
      IF (ioutl .eq. 2) cgtol = cgtol2
    END IF

    IF (par%ismasterproc) print*,'Reorthogonalized CG'
    IF (par%ismasterproc) print*,'Convergence criterion:',cgtol

     !==========================================================================
     ! A. Check Maximum ABS Errors of Hx - Obs for V, U, T, Q, PS 
     !==========================================================================
#define Check_HxmObs_on
#ifdef Check_HxmObs_on     
       IF (par%ismasterproc) WRITE(iulog,*) 'Check Hx_bg - Obs'
       ierr = TimingStart('Check_HxmObs')  !--- profiling
       DO itime = 1, ntime
          IF (obcheck) THEN
             CALL rej_obs_omb(itime)  ! rejection of obs. based on omb check (HJS)
          END IF
          CALL Check_HxmObs(itime)  ! extended to 4D
       END DO
       IF (.NOT. obcheck_outer) THEN
          obcheck=.false.
       END IF
       ierr = TimingStop('Check_HxmObs')
#endif
        tmpmax= maxval(gs_md)
        tmpmin= minval(gs_md)
        CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
        CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
        IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [gs_md]:',minv, maxv

     !==========================================================================
     ! RTTOV, ROPP calculation for outer-loop
     !==========================================================================

        IF (ioutl .gt. 1) THEN
        !! RTTOV used
         IF (RTTOVused) THEN
           DO itime = 1, ntime
              CALL CalRTTOV(itime,gs_md)
           END DO
         END IF
         !! ROPP used
         IF (ROPPused) THEN
           DO itime = 1, ntime
              CALL CalROPP(itime,gs_md)
           END DO
         END IF
        END IF

     !==========================================================================
     ! B. Make b
     !==========================================================================
       IF (par%ismasterproc) WRITE(iulog,*) 'Make_RHS_b'
       ierr = TimingStart('Make_RHS_b')  !--- profiling
         x_eig_acc = x_eig_acc + x_eig
         x_eig_acc_chi = x_eig_acc_chi + x_eig_chi
         x_eig_acc_Tps = x_eig_acc_Tps + x_eig_Tps
         x_eig_acc_q = x_eig_acc_q + x_eig_q
         x_eig_acc_a = x_eig_acc_a + x_eig_a
         x_eig = 0.D0
         x_eig_chi = 0.D0
         x_eig_Tps = 0.D0
         x_eig_q = 0.D0
         x_eig_a = 0.D0
         CALL Make_RHS_b(gs_md,       &
                         b_eig,       &
                         b_eig_chi,   &
                         b_eig_Tps,   &
                         b_eig_q,     &
                         b_eig_a    &    ! modified for alpha-variable
                         )
       ierr = TimingStop('Make_RHS_b')
     !==========================================================================
   
     !==========================================================================
     ! C. Evaluation of adjoint operators
     !==========================================================================
#define Check_Adj_on
#ifdef Check_Adj_on     
       IF (par%ismasterproc) WRITE(iulog,*) 'Check_Adj'
       ierr = TimingStart('Check_Adj')  !--- profiling
         CALL Check_Adj(gs_md,       &
                        b_eig,       &
                        b_eig_chi,   &
                        b_eig_Tps,   &
                        b_eig_q,     &
                        b_eig_a)
       ierr = TimingStop('Check_Adj')
#endif
     !==========================================================================
   
     !==========================================================================
     ! D. Make Ax
     !==========================================================================
       ierr = TimingStart('Make_LHS_Ax0')  !--- profiling
         CALL Make_LHS_Ax(gs_md,       &
                          x_eig,       &
                          x_eig_chi,   &
                          x_eig_Tps,   &
                          x_eig_q,     &
                          x_eig_a,     &
                          Ap_eig,      &
                          Ap_eig_chi,  &
                          Ap_eig_Tps,  &
                          Ap_eig_q,    &
                          Ap_eig_a)
       ierr = TimingStop('Make_LHS_Ax0')
     !==========================================================================
   
     !==========================================================================
     ! E. Initialize r, p, and rsold
     !==========================================================================
     ! 1. Initialize r and p
       DO wn = 1, zwn_nproc
          r_eig(:,wn) = b_eig(:,wn) - Ap_eig(:,wn)
          r_eig_chi(:,wn) = b_eig_chi(:,wn) - Ap_eig_chi(:,wn)
          r_eig_Tps(:,wn) = b_eig_Tps(:,wn) - Ap_eig_Tps(:,wn)
          r_eig_q(:,wn) = b_eig_q(:,wn) - Ap_eig_q(:,wn)
          p_eig(:,wn) = r_eig(:,wn)
          p_eig_chi(:,wn) = r_eig_chi(:,wn)
          p_eig_Tps(:,wn) = r_eig_Tps(:,wn)
          p_eig_q(:,wn) = r_eig_q(:,wn)
       END DO
   
       DO ismpl = 1, nsmpl  !-- for alpha-variable
       DO wn = 1, zwn_nproc_a
          r_eig_a(:,wn,ismpl) = b_eig_a(:,wn,ismpl) - Ap_eig_a(:,wn,ismpl)
          p_eig_a(:,wn,ismpl) = r_eig_a(:,wn,ismpl)
       END DO
       END DO ! ismpl
   
     ! 2. Initialize rsold
       ierr = TimingStart('InnerProduct_r0')   !--- profiling
         CALL InnerProduct(neig,        &
                           neig_chi,    &
                           neig_Tps,    &
                           neig_q,      &
                           zwn_nproc,   &
                           zwn_nproc_a, &
                           r_eig,       &
                           r_eig_chi,   &
                           r_eig_Tps,   &
                           r_eig_q,     &
                           r_eig_a,     &   ! modified for alpha-variable
                           r_eig,       &
                           r_eig_chi,   &
                           r_eig_Tps,   &
                           r_eig_q,     &
                           r_eig_a,     & 
                           rsold)
       ierr = TimingStop('InnerProduct_r0')   !--- profiling
       IF (par%ismasterproc) WRITE(iulog,*) 'rsold,',rsold

       rs0   = rsold

     ! 3. Initialize r*d_eig_history : Save the residual history
       r_eig_history(:,:,:) = 0.d0
       r_eig_chi_history(:,:,:) = 0.d0
       r_eig_Tps_history(:,:,:) = 0.d0
       r_eig_q_history(:,:,:) = 0.d0

       r_eig_q_history_a(:,:,:,:) = 0.d0

       r_eig_history(:,:,1) = r_eig / DSQRT(rsold)
       r_eig_chi_history(:,:,1) = r_eig_chi / DSQRT(rsold)
       r_eig_Tps_history(:,:,1) = r_eig_Tps / DSQRT(rsold)
       r_eig_q_history(:,:,1) = r_eig_q / DSQRT(rsold)

       r_eig_q_history_a(:,:,:,1) = r_eig_a / DSQRT(rsold)
    !==========================================================================
   
    !==========================================================================
    ! F. Main iteration in conjugate gradient
    !==========================================================================
       ierr = TimingStart('Innerloop')   !--- profiling
       DO iter = 1, niter

!       IF (par%ismasterproc) WRITE(iulog,*) 'x_eig(:,1)'
!       IF (par%ismasterproc) WRITE(iulog,*)  x_eig(:,1)
        tmpmax= maxval(x_eig)
        tmpmin= minval(x_eig)
        CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
        CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
        IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [x_eig]:',minv, maxv
         
        tmpmax= maxval(x_eig_chi)
        tmpmin= minval(x_eig_chi)
        CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
        CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
        IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [x_eig_chi]:',minv, maxv

        tmpmax= maxval(x_eig_Tps)
        tmpmin= minval(x_eig_Tps)
        CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
        CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
        IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [x_eig_Tps]:',minv, maxv

        tmpmax= maxval(x_eig_q)
        tmpmin= minval(x_eig_q)
        CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
        CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
        IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [x_eig_q]:',minv, maxv
 
        tmpmax= maxval(x_eig_a)
        tmpmin= minval(x_eig_a)
        CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
        CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
        IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [x_eig_a]:',minv, maxv
 
!        IF( iter == 1 .or. mod(iter,10).EQ. 0)THEN
        ! 1. Calculate Cost Function
          ierr = TimingStart('Cal_CostFun')   !--- profiling
             CALL Cal_Cost_Func(gs_md,       &
                                x_eig,       &
                                x_eig_chi,   &
                                x_eig_Tps,   &
                                x_eig_q,     &
                                x_eig_a,     &
                                j_o,         &
                                j_b)
          ierr = TimingStop('Cal_CostFun')
!        ENDIF
        ! 2. Calculate Ap ------------------------------------------------------
        ierr = TimingStart('Make_LHS_Ax')  !--- profiling
           CALL Make_LHS_Ax(gs_md,        &
                            p_eig,        &
                            p_eig_chi,    &
                            p_eig_Tps,    &
                            p_eig_q,      &
                            p_eig_a,      &
                            Ap_eig,       &
                            Ap_eig_chi,   &
                            Ap_eig_Tps,   &
                            Ap_eig_q,     &
                            Ap_eig_a)
        ierr = TimingStop('Make_LHS_Ax')

!       IF (par%ismasterproc) WRITE(iulog,*) 'Ap_eig(:,1)'
!       IF (par%ismasterproc) WRITE(iulog,*)  Ap_eig(:,1)
!       tmpmax= maxval(Ap_eig)
!       tmpmin= minval(Ap_eig)
!       CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
!       CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
!       IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [Ap_eig]:',minv, maxv
   
        ! 3-1. Calculate alpha, x, r, and rsnew
        ierr = TimingStart('InnerProduct_pAp')   !--- profiling
           CALL InnerProduct(neig,        &
                             neig_chi,    &
                             neig_Tps,    &
                             neig_q,      &
                             zwn_nproc,   &
                             zwn_nproc_a, &
                             p_eig,       &
                             p_eig_chi,   &
                             p_eig_Tps,   &
                             p_eig_q,     &
                             p_eig_a,     &
                             Ap_eig,      &
                             Ap_eig_chi,  &
                             Ap_eig_Tps,  &
                             Ap_eig_q,    &
                             Ap_eig_a,    &
                             pAp)
        ierr = TimingStop('InnerProduct_pAp')   !--- profiling

!       IF (par%ismasterproc) WRITE(iulog,*) 'pAp,',pAp
!       tmp = 0.D0
!       DO wn = 1, zwn_nproc
!          IF( wn+zwn_iproc .LE. zwn2 )THEN
!             tmp = tmp + SUM(p_eig(:,wn)*Ap_eig(:,wn))
!          END IF
!       END DO
!       global_shared_buf(nets:nete,1) = tmp/DBLE(nete-nets+1)
!       CALL Wrap_Repro_Sum(nvars=1, comm=hybrid%par%comm)
!       pAp = global_shared_sum(1)
!       IF (par%ismasterproc) WRITE(iulog,*) 'pAp,',pAp

        alpha = rsold/pAp

         ! 3-2. Update the iterate and the residual        
        DO wn = 1, zwn_nproc
           x_eig(:,wn) = x_eig(:,wn) + alpha*p_eig(:,wn)
           x_eig_chi(:,wn) = x_eig_chi(:,wn) + alpha*p_eig_chi(:,wn)
           x_eig_Tps(:,wn) = x_eig_Tps(:,wn) + alpha*p_eig_Tps(:,wn)
           x_eig_q(:,wn) = x_eig_q(:,wn) + alpha*p_eig_q(:,wn)
           r_eig(:,wn) = r_eig(:,wn) - alpha*Ap_eig(:,wn)
           r_eig_chi(:,wn) = r_eig_chi(:,wn) - alpha*Ap_eig_chi(:,wn)
           r_eig_Tps(:,wn) = r_eig_Tps(:,wn) - alpha*Ap_eig_Tps(:,wn)
           r_eig_q(:,wn) = r_eig_q(:,wn) - alpha*Ap_eig_q(:,wn)
        END DO

        DO ismpl = 1, nsmpl
        DO wn = 1, zwn_nproc_a
           x_eig_a(:,wn,ismpl) = x_eig_a(:,wn,ismpl) + alpha*p_eig_a(:,wn,ismpl)
           r_eig_a(:,wn,ismpl) = r_eig_a(:,wn,ismpl) - alpha*Ap_eig_a(:,wn,ismpl)
        END DO
        END DO ! ismpl

!       IF (par%ismasterproc) WRITE(iulog,*) 'p_eig(:,1)'
!       IF (par%ismasterproc) WRITE(iulog,*)  p_eig(:,1)
!       tmpmax= maxval(p_eig)
!       tmpmin= minval(p_eig)
!       CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
!       CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
!       IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [p_eig]:',minv, maxv

!       IF (par%ismasterproc) WRITE(iulog,*) 'x_eig(:,1)'
!       IF (par%ismasterproc) WRITE(iulog,*)  x_eig(:,1)
!       tmpmax= maxval(x_eig)
!       tmpmin= minval(x_eig)
!       CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
!       CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
!       IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [x_eig]:',minv, maxv

!       IF (par%ismasterproc) WRITE(iulog,*) 'Ap_eig(:,1)'
!       IF (par%ismasterproc) WRITE(iulog,*)  Ap_eig(:,1)
!       tmpmax= maxval(Ap_eig)
!       tmpmin= minval(Ap_eig)
!       CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
!       CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
!       IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [Ap_eig]:',minv, maxv
 
         ! 3-3. Re-orthogonalization step 1 
         ! Remove component in previous residual direction
          DO k = 1, iter

             ! calculate the dot_product of r_k and r_j
             ! r_k : previous residual
             ! r_j : current  residual
             ierr = TimingStart('InnerProduct_rirj')   !--- profiling
               CALL InnerProduct(neig,                        &
                                 neig_chi,                    &
                                 neig_Tps,                    &
                                 neig_q,                      &
                                 zwn_nproc,                   &
                                 zwn_nproc_a,                 &
                                 r_eig,                       &
                                 r_eig_chi,                   &
                                 r_eig_Tps,                   &
                                 r_eig_q,                     &
                                 r_eig_a,                     &
                                 r_eig_history(:,:,k),        &
                                 r_eig_chi_history(:,:,k),    &
                                 r_eig_Tps_history(:,:,k),    &
                                 r_eig_q_history(:,:,k),      &
                                 r_eig_q_history_a(:,:,:,k),  &
                                 rjrk)
             ierr = TimingStop('InnerProduct_rirj')   !--- profiling
             DO wn = 1, zwn_nproc
                r_eig(:,wn)   = r_eig(:,wn) - rjrk*r_eig_history(:,wn,k) 
                r_eig_chi(:,wn) = r_eig_chi(:,wn) - rjrk*r_eig_chi_history(:,wn,k)
                r_eig_Tps(:,wn) = r_eig_Tps(:,wn) - rjrk*r_eig_Tps_history(:,wn,k)
                r_eig_q(:,wn) = r_eig_q(:,wn) - rjrk*r_eig_q_history(:,wn,k)
             END DO

             DO ismpl = 1, nsmpl
             DO wn = 1, zwn_nproc_a
                r_eig_a(:,wn,ismpl) = r_eig_a(:,wn,ismpl) - rjrk*r_eig_q_history_a(:,wn,ismpl,k)
             END DO
             END DO ! ismpl
          END DO
  
          ierr = TimingStart('InnerProduct_rsnew')   !--- profiling
            CALL InnerProduct(neig,        &
                              neig_chi,    &
                              neig_Tps,    &
                              neig_q,      &
                              zwn_nproc,   &
                              zwn_nproc_a, &
                              r_eig,       &
                              r_eig_chi,   &
                              r_eig_Tps,   &
                              r_eig_q,     &
                              r_eig_a,     &
                              r_eig,       &
                              r_eig_chi,   &
                              r_eig_Tps,   &
                              r_eig_q,     &
                              r_eig_a,     &
                              rsnew)
          ierr = TimingStop('InnerProduct_rsnew')   !--- profiling
          beta = rsnew/rsold
          criteria = rsnew/rs0
          rsold = rsnew
!             IF( iter == 1 )  rs0 = rsnew

        ! 4. Check rsnew
          IF (par%ismasterproc) WRITE(iulog,*) 'Check residual'
          IF (par%ismasterproc)                                           &
             WRITE(iulog,*) 'ioutl = ', ioutl,                            &
                            'iter = ', iter, 'rsnew = ', rsnew, criteria 
          IF (par%ismasterproc) WRITE(21,*) ioutl,iter,criteria,j_b+j_o,j_b,j_o
 
        ! Convergence criterion of criteri as normalized rs
          IF (criteria .LT. cgtol) THEN
             IF (par%ismasterproc)                               &
                WRITE(iulog,*) 'ioutl = ', ioutl,                &
                               'total iteration number = ', iter
             EXIT
          END IF
   
        ! 5. Calculate p and rsold
          DO wn = 1, zwn_nproc
             p_eig(:,wn) = r_eig(:,wn) + beta*p_eig(:,wn)
             p_eig_chi(:,wn) = r_eig_chi(:,wn) + beta*p_eig_chi(:,wn)
             p_eig_Tps(:,wn) = r_eig_Tps(:,wn) + beta*p_eig_Tps(:,wn)
             p_eig_q(:,wn) = r_eig_q(:,wn) + beta*p_eig_q(:,wn)
          END DO

          DO ismpl = 1, nsmpl
          DO wn = 1, zwn_nproc_a
             p_eig_a(:,wn,ismpl) = r_eig_a(:,wn,ismpl) + beta*p_eig_a(:,wn,ismpl)
          END DO
          END DO ! ismpl

        ! 6. Re-orthogonalization step 2 : 
        ! normalize residual vector and save it
          if( iter < niter ) then
             r_eig_history(:,:,iter+1)   = r_eig / DSQRT(rsnew)
             r_eig_chi_history(:,:,iter+1)   = r_eig_chi / DSQRT(rsnew)
             r_eig_Tps_history(:,:,iter+1)   = r_eig_Tps / DSQRT(rsnew)
             r_eig_q_history(:,:,iter+1)   = r_eig_q / DSQRT(rsnew)

             r_eig_q_history_a(:,:,:,iter+1)   = r_eig_a / DSQRT(rsnew)
          end if

!          IF (iter .EQ. 1 ) THEN
!              nlanczos = 1
!              Lanczos_diag = 0.d0
!              Lanczos_subd = 0.d0
!              Lanczos_diag(iter) = 1.d0 / alpha
!          ELSE
!              nlanczos = nlanczos + 1
!              Lanczos_diag(iter)   = 1.d0 / alpha + betaold / alphaold
!              Lanczos_subd(iter-1) = DSQRT(betaold) / alphaold
!          END IF
!          alphaold = alpha
!          betaold  = beta

       END DO ! for iter
       ierr = TimingStop('Innerloop')   !--- profiling
 
  END SUBROUTINE RCGmethod 
!===============================================================================
!===============================================================================
  SUBROUTINE Check_HxmObs(itime)
    INTEGER(KIND=int_kind), INTENT(IN) :: itime
  ! 2-2. Local
    REAL(KIND=real_kind), DIMENSION(17) :: tmpmax, tmpmin
    REAL(KIND=real_kind)                 &
         :: maxsdu,maxsdv,maxsdt,maxsdq, &
            maxsfu,maxsfv,maxsft,maxsfq, maxsfp, &
            minsdu,minsdv,minsdt,minsdq, &
            minsfu,minsfv,minsft,minsfq, minsfp, &
            maxaru,maxarv,maxart,        &
            minaru,minarv,minart,        &
            maxamvu, maxamvv,            &
            minamvu, minamvv,            &
            maxscatwindu, maxscatwindv,  &
            minscatwindu, minscatwindv,  &
            maxbgp, minbgp


    INTEGER(KIND=int_kind) &
       :: i, j, k, ie

! For verification -------------------------------------------------------------- ihKwon

    CALL ModelToObs(itime,bg_md(:,:,:,:,itime), gs_obs(itime))

    IF (par%ismasterproc) THEN
       print*,'Variable,     Obs,                  Background,                    O-B,                  Observation error'
       print*,'itime= ', itime
    ENDIF

    tmpmax = 0.
    tmpmin = 0.

    ! single observation test
    IF( singleobs )THEN
    DO ie = nets, nete
    DO i  = 1, obsmet(itime)%nsg_e(ie)
       IF ((gs_obs(itime)%single(i,ie)%u%val .NE. sdmiss) .AND. (obs(itime)%single(i,ie)%u%val .NE. sdmiss))THEN
          tmpmax(1)= max( obs(itime)%single(i,ie)%u%val-gs_obs(itime)%single(i,ie)%u%val,tmpmax(1) )
          tmpmin(1)= min( obs(itime)%single(i,ie)%u%val-gs_obs(itime)%single(i,ie)%u%val,tmpmin(1) )
!      IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Single U: ',   obs(itime)%single(i,ie)%u%val,gs_obs(itime)%single(i,ie)%u%val, &
                           obs(itime)%single(i,ie)%u%val-gs_obs(itime)%single(i,ie)%u%val,       &
                           sqrt(obs(itime)%single(i,ie)%u%error)
!      END IF
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsg_e(ie)
       IF ((gs_obs(itime)%single(i,ie)%v%val .NE. sdmiss) .AND. (obs(itime)%single(i,ie)%v%val .NE. sdmiss))THEN
          tmpmax(2)= max( obs(itime)%single(i,ie)%v%val-gs_obs(itime)%single(i,ie)%v%val,tmpmax(1) )
          tmpmin(2)= min( obs(itime)%single(i,ie)%v%val-gs_obs(itime)%single(i,ie)%v%val,tmpmin(1) )
!      IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Single V: ',   obs(itime)%single(i,ie)%v%val,gs_obs(itime)%single(i,ie)%v%val, &
                           obs(itime)%single(i,ie)%v%val-gs_obs(itime)%single(i,ie)%v%val,       &
                           sqrt(obs(itime)%single(i,ie)%v%error)
!      END IF
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsg_e(ie)
       IF ((gs_obs(itime)%single(i,ie)%t%val .NE. sdmiss) .AND. (obs(itime)%single(i,ie)%t%val .NE. sdmiss))THEN
          tmpmax(3)= max( obs(itime)%single(i,ie)%t%val-gs_obs(itime)%single(i,ie)%t%val,tmpmax(1) )
          tmpmin(4)= min( obs(itime)%single(i,ie)%t%val-gs_obs(itime)%single(i,ie)%t%val,tmpmin(1) )
!      IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Single T: ',   obs(itime)%single(i,ie)%t%val,gs_obs(itime)%single(i,ie)%t%val, &
                           obs(itime)%single(i,ie)%t%val-gs_obs(itime)%single(i,ie)%t%val,       &
                           sqrt(obs(itime)%single(i,ie)%t%error)
!      END IF
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsg_e(ie)
       IF ((gs_obs(itime)%single(i,ie)%q%val .NE. sdmiss) .AND. (obs(itime)%single(i,ie)%q%val .NE. sdmiss))THEN
          tmpmax(4)= max( obs(itime)%single(i,ie)%q%val-gs_obs(itime)%single(i,ie)%q%val,tmpmax(1) )
          tmpmin(4)= min( obs(itime)%single(i,ie)%q%val-gs_obs(itime)%single(i,ie)%q%val,tmpmin(1) )
!      IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Single Q: ',   obs(itime)%single(i,ie)%q%val,gs_obs(itime)%single(i,ie)%q%val, &
                           obs(itime)%single(i,ie)%q%val-gs_obs(itime)%single(i,ie)%q%val,       &
                           sqrt(obs(itime)%single(i,ie)%q%error)
!      END IF
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsg_e(ie)
       IF ((gs_obs(itime)%single(i,ie)%ps%val .NE. sdmiss) .AND. (obs(itime)%single(i,ie)%ps%val .NE. sdmiss))THEN
          tmpmax(5)= max( obs(itime)%single(i,ie)%ps%val-gs_obs(itime)%single(i,ie)%ps%val,tmpmax(1) )
          tmpmin(5)= min( obs(itime)%single(i,ie)%ps%val-gs_obs(itime)%single(i,ie)%ps%val,tmpmin(1) )
!      IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Single Ps: ',  obs(itime)%single(i,ie)%ps%val,gs_obs(itime)%single(i,ie)%ps%val, &
                           obs(itime)%single(i,ie)%ps%val-gs_obs(itime)%single(i,ie)%ps%val,       &
                           sqrt(obs(itime)%single(i,ie)%ps%error)
!      END IF
       END IF
    END DO
    END DO
    END IF

    ! sonde observation
    IF (da_sonde) THEN
    DO ie = nets, nete
    DO i  = 1, obsmet(itime)%nsd_e(1,ie)
       IF ((gs_obs(itime)%sonde(i,ie)%u%val .NE. sdmiss) .AND. (obs(itime)%sonde(i,ie)%u%val .NE. sdmiss))THEN
          tmpmax(1)= max( obs(itime)%sonde(i,ie)%u%val-gs_obs(itime)%sonde(i,ie)%u%val,tmpmax(1) )
          tmpmin(1)= min( obs(itime)%sonde(i,ie)%u%val-gs_obs(itime)%sonde(i,ie)%u%val,tmpmin(1) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Sonde U: ',   obs(itime)%sonde(i,ie)%u%val,gs_obs(itime)%sonde(i,ie)%u%val, &
                          obs(itime)%sonde(i,ie)%u%val-gs_obs(itime)%sonde(i,ie)%u%val,       &
                          sqrt(obs(itime)%sonde(i,ie)%u%error)
       END IF
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(2,ie)
       IF ((gs_obs(itime)%sonde(i,ie)%v%val .NE. sdmiss) .AND. (obs(itime)%sonde(i,ie)%v%val .NE. sdmiss))THEN
          tmpmax(2)= max( obs(itime)%sonde(i,ie)%v%val-gs_obs(itime)%sonde(i,ie)%v%val,tmpmax(2) )
          tmpmin(2)= min( obs(itime)%sonde(i,ie)%v%val-gs_obs(itime)%sonde(i,ie)%v%val,tmpmin(2) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Sonde V: ',   obs(itime)%sonde(i,ie)%v%val,gs_obs(itime)%sonde(i,ie)%v%val, &
                          obs(itime)%sonde(i,ie)%v%val-gs_obs(itime)%sonde(i,ie)%v%val,       &
                          sqrt(obs(itime)%sonde(i,ie)%v%error)
       END IF
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(3,ie)
       IF ((gs_obs(itime)%sonde(i,ie)%t%val .NE. sdmiss) .AND. (obs(itime)%sonde(i,ie)%t%val .NE. sdmiss))THEN
          tmpmax(3)= max( obs(itime)%sonde(i,ie)%t%val-gs_obs(itime)%sonde(i,ie)%t%val,tmpmax(3) )
          tmpmin(3)= min( obs(itime)%sonde(i,ie)%t%val-gs_obs(itime)%sonde(i,ie)%t%val,tmpmin(3) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Sonde T: ',   obs(itime)%sonde(i,ie)%t%val,gs_obs(itime)%sonde(i,ie)%t%val, &
                                obs(itime)%sonde(i,ie)%t%val-gs_obs(itime)%sonde(i,ie)%t%val, &
                                sqrt(obs(itime)%sonde(i,ie)%t%error)
       END IF
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(4,ie)
       IF ((gs_obs(itime)%sonde(i,ie)%q%val .NE. sdmiss) .AND. (obs(itime)%sonde(i,ie)%q%val .NE. sdmiss))THEN
          tmpmax(4)= max( obs(itime)%sonde(i,ie)%q%val-gs_obs(itime)%sonde(i,ie)%q%val,tmpmax(4) )
          tmpmin(4)= min( obs(itime)%sonde(i,ie)%q%val-gs_obs(itime)%sonde(i,ie)%q%val,tmpmin(4) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Sonde Q: ',   obs(itime)%sonde(i,ie)%q%val,gs_obs(itime)%sonde(i,ie)%q%val, &
                                obs(itime)%sonde(i,ie)%q%val-gs_obs(itime)%sonde(i,ie)%q%val, &
                                sqrt(obs(itime)%sonde(i,ie)%q%error)
       END IF
       END IF
    END DO
    END DO
    END IF

    ! surface observation
    IF (da_surface) THEN
    DO ie = nets, nete
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((gs_obs(itime)%surface(i,ie)%u%val .NE. sdmiss) .AND. (obs(itime)%surface(i,ie)%u%val .NE. sdmiss))THEN
          tmpmax(5)= max( obs(itime)%surface(i,ie)%u%val-gs_obs(itime)%surface(i,ie)%u%val,tmpmax(5) )
          tmpmin(5)= min( obs(itime)%surface(i,ie)%u%val-gs_obs(itime)%surface(i,ie)%u%val,tmpmin(5) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Surface U: ',   obs(itime)%surface(i,ie)%u%val,gs_obs(itime)%surface(i,ie)%u%val, &
                                  obs(itime)%surface(i,ie)%u%val-gs_obs(itime)%surface(i,ie)%u%val, &
                                  sqrt(obs(itime)%surface(i,ie)%u%error)
       END IF
       END IF
       IF ((gs_obs(itime)%surface(i,ie)%v%val .NE. sdmiss) .AND. (obs(itime)%surface(i,ie)%v%val .NE. sdmiss))THEN
          tmpmax(6)= max( obs(itime)%surface(i,ie)%v%val-gs_obs(itime)%surface(i,ie)%v%val,tmpmax(6) )
          tmpmin(6)= min( obs(itime)%surface(i,ie)%v%val-gs_obs(itime)%surface(i,ie)%v%val,tmpmin(6) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Surface V: ',   obs(itime)%surface(i,ie)%v%val,gs_obs(itime)%surface(i,ie)%v%val, &
                                  obs(itime)%surface(i,ie)%v%val-gs_obs(itime)%surface(i,ie)%v%val, &
                                  sqrt(obs(itime)%surface(i,ie)%v%error)
       END IF
       END IF
       IF ((gs_obs(itime)%surface(i,ie)%t%val .NE. sdmiss) .AND. (obs(itime)%surface(i,ie)%t%val .NE. sdmiss))THEN
          tmpmax(7)= max( obs(itime)%surface(i,ie)%t%val-gs_obs(itime)%surface(i,ie)%t%val,tmpmax(7) )
          tmpmin(7)= min( obs(itime)%surface(i,ie)%t%val-gs_obs(itime)%surface(i,ie)%t%val,tmpmin(7) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Surface T: ',   obs(itime)%surface(i,ie)%t%val,gs_obs(itime)%surface(i,ie)%t%val, &
                                  obs(itime)%surface(i,ie)%t%val-gs_obs(itime)%surface(i,ie)%t%val, &
                                  sqrt(obs(itime)%surface(i,ie)%t%error)
       END IF
       END IF
       IF ((gs_obs(itime)%surface(i,ie)%q%val .NE. sdmiss) .AND. (obs(itime)%surface(i,ie)%q%val .NE. sdmiss))THEN
          tmpmax(8)= max( obs(itime)%surface(i,ie)%q%val-gs_obs(itime)%surface(i,ie)%q%val,tmpmax(8) )
          tmpmin(8)= min( obs(itime)%surface(i,ie)%q%val-gs_obs(itime)%surface(i,ie)%q%val,tmpmin(8) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Surface Q: ',   obs(itime)%surface(i,ie)%q%val,gs_obs(itime)%surface(i,ie)%q%val, &
                                  obs(itime)%surface(i,ie)%q%val-gs_obs(itime)%surface(i,ie)%q%val, &
                                  sqrt(obs(itime)%surface(i,ie)%q%error)
       END IF
       END IF
       IF ((gs_obs(itime)%surface(i,ie)%p%val .NE. sdmiss) .AND. (obs(itime)%surface(i,ie)%p%val .NE. sdmiss))THEN
          tmpmax(9)= max( obs(itime)%surface(i,ie)%p%val-gs_obs(itime)%surface(i,ie)%p%val,tmpmax(9) )
          tmpmin(9)= min( obs(itime)%surface(i,ie)%p%val-gs_obs(itime)%surface(i,ie)%p%val,tmpmin(9) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Surface PS: ',   obs(itime)%surface(i,ie)%p%val,gs_obs(itime)%surface(i,ie)%p%val, &
                                   obs(itime)%surface(i,ie)%p%val-gs_obs(itime)%surface(i,ie)%p%val, &
                                   sqrt(obs(itime)%surface(i,ie)%p%error)
       END IF
       END IF
    END DO
    END DO
    END IF

    ! bogus observation
    IF (da_bogus) THEN
    DO ie = nets, nete
    DO i  = 1, obsmet(itime)%nbg_e(ie)
       IF ((gs_obs(itime)%bogus(i,ie)%p%val .NE. sdmiss) .AND. (obs(itime)%bogus(i,ie)%p%val .NE. sdmiss))THEN
          tmpmax(10)= max( obs(itime)%bogus(i,ie)%p%val-gs_obs(itime)%bogus(i,ie)%p%val,tmpmax(10) )
          tmpmin(10)= min( obs(itime)%bogus(i,ie)%p%val-gs_obs(itime)%bogus(i,ie)%p%val,tmpmin(10) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Surface_bg PS: ',obs(itime)%bogus(i,ie)%p%val,gs_obs(itime)%bogus(i,ie)%p%val, &
                                   obs(itime)%bogus(i,ie)%p%val-gs_obs(itime)%bogus(i,ie)%p%val, &
                                   sqrt(obs(itime)%bogus(i,ie)%p%error)
       END IF
       END IF
    END DO
    END DO
    END IF

    ! aircraft observation
    IF (da_aircraft) THEN
    DO ie = nets, nete
    DO i  = 1, obsmet(itime)%nar_e(1,ie)
       IF ((gs_obs(itime)%aircraft(i,ie)%u%val .NE. sdmiss) .AND. (obs(itime)%aircraft(i,ie)%u%val .NE. sdmiss))THEN
          tmpmax(11)= max( obs(itime)%aircraft(i,ie)%u%val-gs_obs(itime)%aircraft(i,ie)%u%val,tmpmax(11) )
          tmpmin(11)= min( obs(itime)%aircraft(i,ie)%u%val-gs_obs(itime)%aircraft(i,ie)%u%val,tmpmin(11) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Aircraft U: ',  obs(itime)%aircraft(i,ie)%u%val,gs_obs(itime)%aircraft(i,ie)%u%val, &
                                  obs(itime)%aircraft(i,ie)%u%val-gs_obs(itime)%aircraft(i,ie)%u%val, &
                                  sqrt(obs(itime)%aircraft(i,ie)%u%error)
       END IF
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nar_e(2,ie)
       IF ((gs_obs(itime)%aircraft(i,ie)%v%val .NE. sdmiss) .AND. (obs(itime)%aircraft(i,ie)%v%val .NE. sdmiss))THEN
          tmpmax(12)= max( obs(itime)%aircraft(i,ie)%v%val-gs_obs(itime)%aircraft(i,ie)%v%val,tmpmax(12) )
          tmpmin(12)= min( obs(itime)%aircraft(i,ie)%v%val-gs_obs(itime)%aircraft(i,ie)%v%val,tmpmin(12) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Aircraft V: ',  obs(itime)%aircraft(i,ie)%v%val,gs_obs(itime)%aircraft(i,ie)%v%val, &
                                  obs(itime)%aircraft(i,ie)%v%val-gs_obs(itime)%aircraft(i,ie)%v%val, &
                                  sqrt(obs(itime)%aircraft(i,ie)%v%error)
       END IF
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nar_e(3,ie)
       IF ((gs_obs(itime)%aircraft(i,ie)%t%val .NE. sdmiss) .AND. (obs(itime)%aircraft(i,ie)%t%val .NE. sdmiss))THEN
          tmpmax(13)= max( obs(itime)%aircraft(i,ie)%t%val-gs_obs(itime)%aircraft(i,ie)%t%val,tmpmax(13) )
          tmpmin(13)= min( obs(itime)%aircraft(i,ie)%t%val-gs_obs(itime)%aircraft(i,ie)%t%val,tmpmin(13) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Aircraft T: ',  obs(itime)%aircraft(i,ie)%t%val,gs_obs(itime)%aircraft(i,ie)%t%val, &
                                  obs(itime)%aircraft(i,ie)%t%val-gs_obs(itime)%aircraft(i,ie)%t%val, &
                                  sqrt(obs(itime)%aircraft(i,ie)%t%error)
       END IF
       END IF
    END DO
    END DO
    END IF

    ! amv observation
    IF (da_amv) THEN
    DO ie = nets, nete
    DO i  = 1, obsmet(itime)%namv_e(ie)
       IF ((gs_obs(itime)%amv(i,ie)%u%val .NE. sdmiss) .AND. (obs(itime)%amv(i,ie)%u%val .NE. sdmiss))THEN
          tmpmax(14)= max( obs(itime)%amv(i,ie)%u%val-gs_obs(itime)%amv(i,ie)%u%val,tmpmax(14) )
          tmpmin(14)= min( obs(itime)%amv(i,ie)%u%val-gs_obs(itime)%amv(i,ie)%u%val,tmpmin(14) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Amv U: ',  obs(itime)%amv(i,ie)%u%val,gs_obs(itime)%amv(i,ie)%u%val, &
                             obs(itime)%amv(i,ie)%u%val-gs_obs(itime)%amv(i,ie)%u%val, &
                             sqrt(obs(itime)%amv(i,ie)%u%error)
       END IF
       END IF
       IF ((gs_obs(itime)%amv(i,ie)%v%val .NE. sdmiss) .AND. (obs(itime)%amv(i,ie)%v%val .NE. sdmiss))THEN
          tmpmax(15)= max( obs(itime)%amv(i,ie)%v%val-gs_obs(itime)%amv(i,ie)%v%val,tmpmax(15) )
          tmpmin(15)= min( obs(itime)%amv(i,ie)%v%val-gs_obs(itime)%amv(i,ie)%v%val,tmpmin(15) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Amv V: ',  obs(itime)%amv(i,ie)%v%val,gs_obs(itime)%amv(i,ie)%v%val, &
                             obs(itime)%amv(i,ie)%v%val-gs_obs(itime)%amv(i,ie)%v%val, &
                             sqrt(obs(itime)%amv(i,ie)%v%error)
       END IF
       END IF
    END DO
    END DO
    END IF

    ! scatwind observation
    IF (da_scatwind) THEN
    DO ie = nets, nete
    DO i  = 1, obsmet(itime)%nscatwind_e(ie)
       IF ((gs_obs(itime)%scatwind(i,ie)%u%val .NE. sdmiss) .AND. (obs(itime)%scatwind(i,ie)%u%val .NE. sdmiss))THEN
          tmpmax(16)= max( obs(itime)%scatwind(i,ie)%u%val-gs_obs(itime)%scatwind(i,ie)%u%val,tmpmax(16) )
          tmpmin(16)= min( obs(itime)%scatwind(i,ie)%u%val-gs_obs(itime)%scatwind(i,ie)%u%val,tmpmin(16) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Scatwind U: ',  obs(itime)%scatwind(i,ie)%u%val,gs_obs(itime)%scatwind(i,ie)%u%val, &
                                  obs(itime)%scatwind(i,ie)%u%val-gs_obs(itime)%scatwind(i,ie)%u%val, &
                                  sqrt(obs(itime)%scatwind(i,ie)%u%error)
       END IF
       END IF
       IF ((gs_obs(itime)%scatwind(i,ie)%v%val .NE. sdmiss) .AND. (obs(itime)%scatwind(i,ie)%v%val .NE. sdmiss))THEN
          tmpmax(17)= max( obs(itime)%scatwind(i,ie)%v%val-gs_obs(itime)%scatwind(i,ie)%v%val,tmpmax(17) )
          tmpmin(17)= min( obs(itime)%scatwind(i,ie)%v%val-gs_obs(itime)%scatwind(i,ie)%v%val,tmpmin(17) )
       IF ( par%ismasterproc .and. i == 1 ) THEN
          print*,'Scatwind V: ',  obs(itime)%scatwind(i,ie)%v%val,gs_obs(itime)%scatwind(i,ie)%v%val, &
                                  obs(itime)%scatwind(i,ie)%v%val-gs_obs(itime)%scatwind(i,ie)%v%val, &
                                  sqrt(obs(itime)%scatwind(i,ie)%v%error)
       END IF
       END IF
    END DO
    END DO
    END IF

    IF( singleobs )THEN
        CALL Par_Reduce(tmpmax(1),maxsdu,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(2),maxsdv,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(3),maxsdt,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(4),maxsdq,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(5),maxsfp,1,PAR_MAX)
        CALL Par_Reduce(tmpmin(1),minsdu,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(2),minsdv,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(3),minsdt,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(4),minsdq,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(5),minsfp,1,PAR_MIN)
    END IF
    IF (da_sonde) THEN
        CALL Par_Reduce(tmpmax(1),maxsdu,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(2),maxsdv,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(3),maxsdt,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(4),maxsdq,1,PAR_MAX)
        CALL Par_Reduce(tmpmin(1),minsdu,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(2),minsdv,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(3),minsdt,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(4),minsdq,1,PAR_MIN)
    END IF
    IF (da_surface) THEN
        CALL Par_Reduce(tmpmax(5),maxsfu,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(6),maxsfv,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(7),maxsft,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(8),maxsfq,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(9),maxsfp,1,PAR_MAX)
        CALL Par_Reduce(tmpmin(5),minsfu,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(6),minsfv,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(7),minsft,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(8),minsfq,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(9),minsfp,1,PAR_MIN)
    END IF
    IF (da_bogus) THEN
        CALL Par_Reduce(tmpmax(10),maxbgp,1,PAR_MAX)
        CALL Par_Reduce(tmpmin(10),minbgp,1,PAR_MIN)
    END IF
    IF (da_aircraft) THEN
        CALL Par_Reduce(tmpmax(11), maxaru,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(12),maxarv,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(13),maxart,1,PAR_MAX)
        CALL Par_Reduce(tmpmin(11), minaru,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(12),minarv,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(13),minart,1,PAR_MIN)
    END IF
    IF (da_amv) THEN
        CALL Par_Reduce(tmpmax(14), maxamvu,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(15), maxamvv,1,PAR_MAX)
        CALL Par_Reduce(tmpmin(14), minamvu,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(15), minamvv,1,PAR_MIN)
    END IF
    IF (da_scatwind) THEN
        CALL Par_Reduce(tmpmax(16), maxscatwindu,1,PAR_MAX)
        CALL Par_Reduce(tmpmax(17), maxscatwindv,1,PAR_MAX)
        CALL Par_Reduce(tmpmin(16), minscatwindu,1,PAR_MIN)
        CALL Par_Reduce(tmpmin(17), minscatwindv,1,PAR_MIN)
    END IF

    IF (par%ismasterproc) THEN
       PRINT*,'Minimum and Maximum of Hx_bg - Obs'
       IF( singleobs )THEN
           PRINT*, 'Single U_bg - U_obs: min=',minsdu,' max=',maxsdu
           PRINT*, 'Single V_bg - V_obs: min=',minsdv,' max=',maxsdv
           PRINT*, 'Single T_bg - T_obs: min=',minsdt,' max=',maxsdt
           PRINT*, 'Single Q_bg - Q_obs: min=',minsdq,' max=',maxsdq
           PRINT*, 'Single Ps_bg - Ps_obs: min=',minsfp,' max=',maxsfp
       END IF
       IF (da_sonde) THEN
           PRINT*, 'Sonde U_bg - U_obs: min=',minsdu,' max=',maxsdu
           PRINT*, 'Sonde V_bg - V_obs: min=',minsdv,' max=',maxsdv
           PRINT*, 'Sonde T_bg - T_obs: min=',minsdt,' max=',maxsdt
           PRINT*, 'Sonde Q_bg - Q_obs: min=',minsdq,' max=',maxsdq
       END IF
       IF (da_surface) THEN
           PRINT*, 'Surface U_bg - U_obs: min=',minsfu,' max=',maxsfu
           PRINT*, 'Surface V_bg - V_obs: min=',minsfv,' max=',maxsfv
           PRINT*, 'Surface T_bg - T_obs: min=',minsft,' max=',maxsft
           PRINT*, 'Surface Q_bg - Q_obs: min=',minsfq,' max=',maxsfq
           PRINT*, 'Surface P_bg - P_obs: min=',minsfp,' max=',maxsfp
       END IF
       IF (da_bogus) THEN
           PRINT*, 'Surface Bogus P_bg - P_obs: min=',minbgp,' max=',maxbgp
       END IF
       IF (da_aircraft) THEN
           PRINT*, 'Aircraft U_bg - U_obs: min=',minaru,' max=',maxaru
           PRINT*, 'Aircraft V_bg - V_obs: min=',minarv,' max=',maxarv
           PRINT*, 'Aircraft T_bg - T_obs: min=',minart,' max=',maxart
       END IF
       IF (da_amv) THEN
           PRINT*, 'Amv U_bg - U_obs: min=',minamvu,' max=',maxamvu
           PRINT*, 'Amv V_bg - V_obs: min=',minamvv,' max=',maxamvv
       END IF
       IF (da_scatwind) THEN
           PRINT*, 'Scatwind U_bg - U_obs: min=',minscatwindu,' max=',maxscatwindu
           PRINT*, 'Scatwind V_bg - V_obs: min=',minscatwindv,' max=',maxscatwindv
       END IF
    END IF

    RETURN

  END SUBROUTINE Check_HxmObs
!===============================================================================
!===============================================================================
  SUBROUTINE rej_obs_omb(itime)

  USE check_omb_mod, ONLY : SondeOBCheck,    &
                            SurfaceOBCheck,  &
                            AircraftOBCheck, &       
                            AmvOBCheck,      &
                            ScatwindOBCheck, &
                            amsuaOBCheck,    &
                            gpsroOBCheck,    &
                            iasiOBCheck,     &
                            crisOBCheck,     &
                            atmsOBCheck,     &
                            atmswvOBCheck,   &
                            mhsOBCheck,      &
                            csrOBCheck

    INTEGER(KIND=int_kind), INTENT(IN) :: itime
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie

!=================================================================   
    CALL ModelToObs(itime,bg_md(:,:,:,:,itime), gs_obs(itime))

    IF (par%ismasterproc) THEN
       print*,'itime= ', itime
    ENDIF

    ! sonde observation
    IF (da_sonde) THEN
       obs(itime)%sonde(:,:) = obs_org(itime)%sonde(:,:)
       CALL SondeOBCheck(itime,"u",obs(itime)%sonde%u%val,gs_obs(itime)%sonde%u%val, 20.D0)
       CALL SondeOBCheck(itime,"v",obs(itime)%sonde%v%val,gs_obs(itime)%sonde%v%val, 20.D0)
       CALL SondeOBCheck(itime,"t",obs(itime)%sonde%t%val,gs_obs(itime)%sonde%t%val, 20.D0)
       IF (cv_opt_hum == 0) THEN
          CALL SondeOBCheck(itime,"q",obs(itime)%sonde%q%val,gs_obs(itime)%sonde%q%val, 1.D-2)
       ELSE IF(cv_opt_hum == 1) THEN
          CALL SondeOBCheck(itime,"q",obs(itime)%sonde%q%val,gs_obs(itime)%sonde%q%val, 2.D0)
       END IF
    END IF

    ! surface observation
    IF (da_surface) THEN
       obs(itime)%surface(:,:) = obs_org(itime)%surface(:,:)
       CALL SurfaceOBCheck(itime,obs(itime)%surface%u%val,gs_obs(itime)%surface%u%val, 20.D0)
       CALL SurfaceOBCheck(itime,obs(itime)%surface%v%val,gs_obs(itime)%surface%v%val, 20.D0)
       CALL SurfaceOBCheck(itime,obs(itime)%surface%t%val,gs_obs(itime)%surface%t%val, 20.D0)
       CALL SurfaceOBCheck(itime,obs(itime)%surface%q%val,gs_obs(itime)%surface%q%val, 1.D-2)
       CALL SurfaceOBCheck(itime,obs(itime)%surface%p%val,gs_obs(itime)%surface%p%val, 1000.D0)
    END IF

    ! aircraft observation
    IF (da_aircraft) THEN
       obs(itime)%aircraft(:,:) = obs_org(itime)%aircraft(:,:)
       CALL AircraftOBCheck(itime,"u",obs(itime)%aircraft%u%val,gs_obs(itime)%aircraft%u%val, 20.D0)
       CALL AircraftOBCheck(itime,"v",obs(itime)%aircraft%v%val,gs_obs(itime)%aircraft%v%val, 20.D0)
       CALL AircraftOBCheck(itime,"t",obs(itime)%aircraft%t%val,gs_obs(itime)%aircraft%t%val, 20.D0)
    END IF

    ! amv observation
    IF (da_amv) THEN
       obs(itime)%amv(:,:) = obs_org(itime)%amv(:,:)
       CALL AmvOBCheck(itime,obs(itime)%amv%u%val,gs_obs(itime)%amv%u%val, 15.D0)
       CALL AmvOBCheck(itime,obs(itime)%amv%v%val,gs_obs(itime)%amv%v%val, 15.D0)
    END IF

    ! scatwind observation
    IF (da_scatwind) THEN
       obs(itime)%scatwind(:,:) = obs_org(itime)%scatwind(:,:)
       CALL ScatwindOBCheck(itime,obs(itime)%scatwind%u%val,gs_obs(itime)%scatwind%u%val, 5.D0)
       CALL ScatwindOBCheck(itime,obs(itime)%scatwind%v%val,gs_obs(itime)%scatwind%v%val, 5.D0)
    END IF

    ! amsua
    IF (da_amsua) THEN
       obs(itime)%amsua(:,:,:) = obs_org(itime)%amsua(:,:,:)
       CALL amsuaOBCheck(itime,                                 &
                         obs(itime)%amsua%tb%val,               &
                         DSQRT(obs(itime)%amsua%tb%error),      &
                         obsmet(itime)%amsuabg%amsuahdx%tb%val, &
                         crit_omb_amsua)
    END IF

    ! gpsroOBCheck
    IF (da_gpsro) THEN
       obs(itime)%gpsro(:,:) = obs_org(itime)%gpsro(:,:)
       CALL gpsroOBCheck(itime,                                &
                         obs(itime)%gpsro%ba%val,              &
                         DSQRT(obs(itime)%gpsro%ba%error),     &
                         obsmet(itime)%gpsrobg%gpsrohdx%ba%val,&
                         crit_omb_gpsro)
    END IF

    ! iasiOBCheck
    IF (da_iasi) THEN
       obs(itime)%iasi(:,:,:) = obs_org(itime)%iasi(:,:,:)
       CALL iasiOBCheck(itime,                               &
                        obs(itime)%iasi%tb%val,              &
                        DSQRT(obs(itime)%iasi%tb%error),     &
                        obsmet(itime)%iasibg%iasihdx%tb%val, &
                        crit_omb_iasi)
    END IF

    ! crisOBCheck
    IF (da_cris) THEN
       obs(itime)%cris(:,:,:) = obs_org(itime)%cris(:,:,:)
       CALL crisOBCheck(itime,                               &
                        obs(itime)%cris%tb%val,              &
                        DSQRT(obs(itime)%cris%tb%error),     &
                        obsmet(itime)%crisbg%crishdx%tb%val, &
                        crit_omb_cris)
    END IF

    ! atmsOBCheck
    IF (da_atms) THEN
       obs(itime)%atms(:,:,:) = obs_org(itime)%atms(:,:,:)
       CALL atmsOBCheck(itime,                               &
                        obs(itime)%atms%tb%val,              &
                        DSQRT(obs(itime)%atms%tb%error),     &
                        obsmet(itime)%atmsbg%atmshdx%tb%val, &
                        crit_omb_atms)
    END IF

    ! atmswvOBCheck
    IF (da_atmswv) THEN
       obs(itime)%atmswv(:,:,:) = obs_org(itime)%atmswv(:,:,:)
       CALL atmswvOBCheck(itime,                                &
                        obs(itime)%atmswv%tb%val,               &
                        DSQRT(obs(itime)%atmswv%tb%error),      &
                        obsmet(itime)%atmswvbg%atmswvhdx%tb%val,&
                        crit_omb_atmswv)
    END IF

    ! mhsOBCheck
    IF (da_mhs) THEN
       obs(itime)%mhs(:,:,:) = obs_org(itime)%mhs(:,:,:)
       CALL mhsOBCheck(itime,                             &
                       obs(itime)%mhs%tb%val,             &
                       obsmet(itime)%mhsbg%mhshdx%tb%val, &
                       crit_omb_mhs)
    END IF

    ! csrOBCheck
    IF (da_csr) THEN
       obs(itime)%csr(:,:,:) = obs_org(itime)%csr(:,:,:)
       CALL csrOBCheck(itime,                             &
                       obs(itime)%csr%tb%val,             &
                       DSQRT(obs(itime)%csr%tb%error),    &
                       obsmet(itime)%csrbg%csrhdx%tb%val, &
                       crit_omb_csr)
    END IF

    RETURN

  END SUBROUTINE rej_obs_omb
!===============================================================================
!===============================================================================
! SUBROUTINE Check_RMSE(an_md,bg_md)
  SUBROUTINE Check_RMSE
  ! 2. Variables
  ! 2-1. Input
!   REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd), &
!      INTENT(IN) :: an_md
!      INTENT(IN) :: bg_md

  ! 2-2. Local
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie,     &
          wn, vr,          &
          n3d,n2d,         &
          nsonde,          &
          nsurfc,          &
          naircra
    real(kind=real_kind) :: rmse_ua, rmse_va, rmse_ta, rmse_qa, rmse_psa, &
                            rmse_ub, rmse_vb, rmse_tb, rmse_qb, rmse_psb, &
                            rmse_sduo, rmse_sdvo, rmse_sdto, rmse_sdqo
 
    rmse_ub   = 0.D0; rmse_vb   = 0.D0; rmse_tb   = 0.D0; rmse_qb   = 0.D0; rmse_psb  = 0.D0 
    rmse_ua   = 0.D0; rmse_va   = 0.D0; rmse_ta   = 0.D0; rmse_qa   = 0.D0; rmse_psa  = 0.D0
    rmse_sduo = 0.D0; rmse_sdvo = 0.D0; rmse_sdto = 0.D0; rmse_sdqo = 0.D0

    !    ierr = TimingStart('ModelToObs_Make_RMSE')  !--- profiling
    !CALL ModelToObs(an_md, gs_obs)
    !    ierr = TimingStop('ModelToObs_Make_RMSE')

    global_shared_buf(:,1:18) = 0.D0

!.. Calculate sonde observation RMSE error
!   IF (da_sonde) THEN
!   DO ie = nets, nete
!   DO i  = 1, nsd_e(1,ie)
!      IF ((obs%sonde(i,ie)%u%val .NE. sdmiss) .AND. (gs_obs%sonde(i,ie)%u%val .NE. sdmiss)) THEN
!         global_shared_buf(ie,2) = global_shared_buf(ie,2)               &
!                                 +(gs_obs%sonde(i,ie)%u%val - obs%sonde(i,ie)%u%val)**2.D0
!         global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
!      END IF
!   END DO
!   DO i  = 1, nsd_e(2,ie)
!      IF ((obs%sonde(i,ie)%v%val .NE. sdmiss) .AND. (gs_obs%sonde(i,ie)%v%val .NE. sdmiss)) THEN
!         global_shared_buf(ie,3) = global_shared_buf(ie,3)               &
!                                 +(gs_obs%sonde(i,ie)%v%val - obs%sonde(i,ie)%v%val)**2.D0
!         global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
!      END IF
!   END DO
!   DO i  = 1, nsd_e(3,ie)
!      IF ((obs%sonde(i,ie)%t%val .NE. sdmiss) .AND. (gs_obs%sonde(i,ie)%t%val .NE. sdmiss)) THEN
!         global_shared_buf(ie,4) = global_shared_buf(ie,4)               &
!                                 +(gs_obs%sonde(i,ie)%t%val - obs%sonde(i,ie)%t%val)**2.D0
!         global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
!      END IF
!   END DO
!   DO i  = 1, nsd_e(4,ie)
!      IF ((obs%sonde(i,ie)%q%val .NE. sdmiss) .AND. (gs_obs%sonde(i,ie)%q%val .NE. sdmiss)) THEN
!         global_shared_buf(ie,5) = global_shared_buf(ie,5)               &
!                                 +(gs_obs%sonde(i,ie)%q%val - obs%sonde(i,ie)%q%val)**2.D0
!         global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
!      END IF
!   END DO
!   END DO
!   END IF

   DO ie = 1, nelemd
   DO  k = 11, nlev
   DO  j = 1, np
   DO  i = 1, np
       global_shared_buf(ie,7) = global_shared_buf(ie,7) + 1.D0
!.. Calculate Background error (U, V, T, and Q)
       global_shared_buf(ie,8)= global_shared_buf(ie,8)                               &
                               + (bg_md(i,j,k,ie,nt_itime)       - nt_org(i,j,k,ie))**2.D0
       global_shared_buf(ie,9)= global_shared_buf(ie,9)                                   &
                               + (bg_md(i,j,nlev+k,ie,nt_itime)  -nt_org(i,j,nlev+k,ie))**2.D0
       global_shared_buf(ie,10)= global_shared_buf(ie,10)                                   &
                               + (bg_md(i,j,2*nlev+k,ie,nt_itime)-nt_org(i,j,2*nlev+k,ie))**2.D0
       global_shared_buf(ie,11)= global_shared_buf(ie,11)                                   &
                               + (bg_md(i,j,3*nlev+k,ie,nt_itime)-nt_org(i,j,3*nlev+k,ie))**2.D0
!.. Calculate Analysis error (U, V, T, and Q)
       global_shared_buf(ie,12)= global_shared_buf(ie,12)                         &
                               + (an_md(i,j,k,ie,nt_itime)        -nt_org(i,j,k,ie))**2.D0
       global_shared_buf(ie,13)= global_shared_buf(ie,13)                              &
                               + (an_md(i,j,nlev+k,ie,nt_itime)   -nt_org(i,j,nlev+k,ie))**2.D0
       global_shared_buf(ie,14)= global_shared_buf(ie,14)                               &
                              + (an_md(i,j,2*nlev+k,ie,nt_itime) -nt_org(i,j,2*nlev+k,ie))**2.D0
       global_shared_buf(ie,15)= global_shared_buf(ie,15)                                &
                               +  (an_md(i,j,3*nlev+k,ie,nt_itime)-nt_org(i,j,3*nlev+k,ie))**2.D0
   END DO
   END DO
   END DO
   END DO

   DO ie = 1, nelemd
   DO  j = 1, np
   DO  i = 1, np
       global_shared_buf(ie,16) = global_shared_buf(ie,16) + 1.D0
!.. Calculate Background error (PS)
       global_shared_buf(ie,17)= global_shared_buf(ie,17)                                       &
                               + (bg_md(i,j,4*nlev+1,ie,nt_itime)   - nt_org(i,j,4*nlev+1,ie))**2.D0
!.. Calculate Analysis error (PS)
       global_shared_buf(ie,18)= global_shared_buf(ie,18)                                       &
                               + (an_md(i,j,4*nlev+1,ie,nt_itime)       - nt_org(i,j,4*nlev+1,ie))**2.D0
   END DO
   END DO
   END DO

    ierr = TimingStart('Repro_Sum_Make_RMSE')  !--- profiling
    CALL Wrap_Repro_Sum(nvars=18, comm=hybrid%par%comm)
    ierr = TimingStop('Repro_Sum_Make_RMSE')  !--- profiling

!sonde observation
!   rmse_sduo = global_shared_sum(2)
!   rmse_sdvo = global_shared_sum(3)
!   rmse_sdto = global_shared_sum(4)
!   rmse_sdqo = global_shared_sum(5)
!   nsonde    = global_shared_sum(6)

    n3d       = global_shared_sum(7)
    rmse_ub   = global_shared_sum(8)
    rmse_vb   = global_shared_sum(9)
    rmse_tb   = global_shared_sum(10)
    rmse_qb   = global_shared_sum(11)
    rmse_ua   = global_shared_sum(12)
    rmse_va   = global_shared_sum(13)
    rmse_ta   = global_shared_sum(14)
    rmse_qa   = global_shared_sum(15)
    n2d       = global_shared_sum(16)
    rmse_psb  = global_shared_sum(17)
    rmse_psa  = global_shared_sum(18)

!   IF (da_sonde) THEN
!   IF (par%isMasterProc) WRITE(iulog,*) "RMSE against Sonde observation"
!   IF (par%isMasterProc) WRITE(iulog,*) "U  =", DSQRT(rmse_sduo/nsonde)
!   IF (par%isMasterProc) WRITE(iulog,*) "V  =", DSQRT(rmse_sdvo/nsonde)
!   IF (par%isMasterProc) WRITE(iulog,*) "T  =", DSQRT(rmse_sdto/nsonde)
!   IF (par%isMasterProc) WRITE(iulog,*) "Q  =", DSQRT(rmse_sdqo/nsonde)
!   END IF

    IF (par%isMasterProc) WRITE(iulog,*) "RMSE of U (from level 11 to bottom)"
    IF (par%isMasterProc) WRITE(iulog,*) "Background =", DSQRT(rmse_ub/n3d)
    IF (par%isMasterProc) WRITE(iulog,*) "Analysis   =", DSQRT(rmse_ua/n3d)
    IF (par%isMasterProc) WRITE(iulog,*) "RMSE of V (from level 11 to bottom)"
    IF (par%isMasterProc) WRITE(iulog,*) "Background =", DSQRT(rmse_vb/n3d)
    IF (par%isMasterProc) WRITE(iulog,*) "Analysis   =", DSQRT(rmse_va/n3d)
    IF (par%isMasterProc) WRITE(iulog,*) "RMSE of T (from level 11 to bottom)"
    IF (par%isMasterProc) WRITE(iulog,*) "Background =", DSQRT(rmse_tb/n3d)
    IF (par%isMasterProc) WRITE(iulog,*) "Analysis   =", DSQRT(rmse_ta/n3d)
    IF (par%isMasterProc) WRITE(iulog,*) "RMSE of Q (from level 11 to bottom)"
    IF (par%isMasterProc) WRITE(iulog,*) "Background =", DSQRT(rmse_qb/n3d)
    IF (par%isMasterProc) WRITE(iulog,*) "Analysis   =", DSQRT(rmse_qa/n3d)
    IF (par%isMasterProc) WRITE(iulog,*) "RMSE of PS"
    IF (par%isMasterProc) WRITE(iulog,*) "Background =", DSQRT(rmse_psb/n2d)
    IF (par%isMasterProc) WRITE(iulog,*) "Analysis   =", DSQRT(rmse_psa/n2d)

    RETURN   

  END SUBROUTINE Check_RMSE
!===============================================================================
!===============================================================================
! SUBROUTINE Check_OB(an_md, bg_md)
  SUBROUTINE Check_OB(itime)

    implicit none

    INTEGER(KIND=int_kind), INTENT(IN) :: itime
!   REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd), &
!      INTENT(IN) :: an_md, bg_md

  !  Local
    real(kind=real_kind) :: obmsdu, obmsdv, obmsdt, obmsdq          ! O-B, Average, Sonde
    real(kind=real_kind) :: obmsfu, obmsfv, obmsft, obmsfq, obmsfps ! O-B, Average, Surface
    real(kind=real_kind) :: obdsdu, obdsdv, obdsdt, obdsdq          ! O-B, SDV,     Sonde
    real(kind=real_kind) :: obdsfu, obdsfv, obdsft, obdsfq, obdsfps ! O-B, SDV,     Surface
    real(kind=real_kind) :: obmbgps                                 ! O-B, Average, Bogus
    real(kind=real_kind) :: obdbgps                                 ! O-B, SDV,     Bogus
    real(kind=real_kind) :: oamsdu, oamsdv, oamsdt, oamsdq          ! O-A, Average, Sonde
    real(kind=real_kind) :: oamsfu, oamsfv, oamsft, oamsfq, oamsfps ! O-A, Average, Surface
    real(kind=real_kind) :: oadsdu, oadsdv, oadsdt, oadsdq          ! O-A, SDV,     Sonde
    real(kind=real_kind) :: oadsfu, oadsfv, oadsft, oadsfq, oadsfps ! O-A, SDV,     Surface
    real(kind=real_kind) :: oambgps                                 ! O-A, Average, Bogus
    real(kind=real_kind) :: oadbgps                                 ! O-A, SDV,     Bogus
    real(kind=real_kind) :: nobsdu, nobsdv, nobsdt, nobsdq          ! number of O-B, Sonde
    real(kind=real_kind) :: nobsfu, nobsfv, nobsft, nobsfq, nobsfps ! number of O-B, Surface
    real(kind=real_kind) :: noasdu, noasdv, noasdt, noasdq          ! number of O-A, Sonde
    real(kind=real_kind) :: noasfu, noasfv, noasft, noasfq, noasfps ! number of O-A, Surface
    real(kind=real_kind) :: nobbgps                                 ! number of O-B, Bogus
    real(kind=real_kind) :: noabgps                                 ! number of O-A, Bogus

    integer(kind=int_kind) :: i, k, ie, ierr

    ! background
    CALL ModelToObs(itime,bg_md(:,:,:,:,itime), bg_obs(itime))

    ! analysis
    CALL ModelToObs(itime,an_md(:,:,:,:,itime), an_obs(itime))

!.. Mean of O-B and O-A for Single Observatino
    IF( singleobs )THEN

    DO ie = nets, nete
    DO i  = 1, obsmet(itime)%nsg_e(ie)
       print*,'i,ie:',i,ie
       print*,'bg_obs(itime)%single(i,ie)%u%val,an_obs(itime)%single(i,ie)%u%val', &
               bg_obs(itime)%single(i,ie)%u%val,an_obs(itime)%single(i,ie)%u%val
       print*,'bg_obs(itime)%single(i,ie)%v%val,an_obs(itime)%single(i,ie)%v%val', &
               bg_obs(itime)%single(i,ie)%v%val,an_obs(itime)%single(i,ie)%v%val
       print*,'bg_obs(itime)%single(i,ie)%t%val,an_obs(itime)%single(i,ie)%t%val', &
               bg_obs(itime)%single(i,ie)%t%val,an_obs(itime)%single(i,ie)%t%val
       print*,'bg_obs(itime)%single(i,ie)%q%val,an_obs(itime)%single(i,ie)%q%val', &
               bg_obs(itime)%single(i,ie)%q%val,an_obs(itime)%single(i,ie)%q%val
       print*,'bg_obs(itime)%single(i,ie)%ps%val,an_obs(itime)%single(i,ie)%ps%val', &
               bg_obs(itime)%single(i,ie)%ps%val,an_obs(itime)%single(i,ie)%ps%val
    END DO
    END DO

    END IF

!.. Mean of O-B and O-A for Sonde
    IF (da_sonde) THEN

    obmsdu= 0.D0; obmsdv= 0.D0; obmsdt= 0.D0; obmsdq= 0.D0
    obdsdu= 0.D0; obdsdv= 0.D0; obdsdt= 0.D0; obdsdq= 0.D0
    oamsdu= 0.D0; oamsdv= 0.D0; oamsdt= 0.D0; oamsdq= 0.D0
    oadsdu= 0.D0; oadsdv= 0.D0; oadsdt= 0.D0; oadsdq= 0.D0
    nobsdu= 0   ; nobsdv= 0   ; nobsdt= 0   ; nobsdq= 0
    noasdu= 0   ; noasdv= 0   ; noasdt= 0   ; noasdq= 0

    global_shared_buf(:,1:18) = 0.D0

    DO ie = nets, nete
   ! background
    DO i  = 1, obsmet(itime)%nsd_e(1,ie)
       IF ((obs(itime)%sonde(i,ie)%u%val .NE. sdmiss) .AND. (bg_obs(itime)%sonde(i,ie)%u%val .NE. sdmiss)) THEN
          global_shared_buf(ie,1) = global_shared_buf(ie,1)                           &
                                  +(obs(itime)%sonde(i,ie)%u%val - bg_obs(itime)%sonde(i,ie)%u%val)
          global_shared_buf(ie,11)= global_shared_buf(ie,11) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(2,ie)
       IF ((obs(itime)%sonde(i,ie)%v%val .NE. sdmiss) .AND. (bg_obs(itime)%sonde(i,ie)%v%val .NE. sdmiss)) THEN
          global_shared_buf(ie,2) = global_shared_buf(ie,2)                           &
                                  +(obs(itime)%sonde(i,ie)%v%val - bg_obs(itime)%sonde(i,ie)%v%val)
          global_shared_buf(ie,12)= global_shared_buf(ie,12) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(3,ie)
       IF ((obs(itime)%sonde(i,ie)%t%val .NE. sdmiss) .AND. (bg_obs(itime)%sonde(i,ie)%t%val .NE. sdmiss)) THEN
          global_shared_buf(ie,3) = global_shared_buf(ie,3)                           &
                                  +(obs(itime)%sonde(i,ie)%t%val - bg_obs(itime)%sonde(i,ie)%t%val)
          global_shared_buf(ie,13)= global_shared_buf(ie,13) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(4,ie)
       IF ((obs(itime)%sonde(i,ie)%q%val .NE. sdmiss) .AND. (bg_obs(itime)%sonde(i,ie)%q%val .NE. sdmiss)) THEN
          global_shared_buf(ie,4) = global_shared_buf(ie,4)                           &
                                  +(obs(itime)%sonde(i,ie)%q%val - bg_obs(itime)%sonde(i,ie)%q%val)
          global_shared_buf(ie,14)= global_shared_buf(ie,14) + 1.D0
       END IF
    END DO

   ! analysis
    DO i  = 1, obsmet(itime)%nsd_e(1,ie)
       IF ((obs(itime)%sonde(i,ie)%u%val .NE. sdmiss) .AND. (an_obs(itime)%sonde(i,ie)%u%val .NE. sdmiss)) THEN
          global_shared_buf(ie,5) = global_shared_buf(ie,5)                           &
                                  +(obs(itime)%sonde(i,ie)%u%val - an_obs(itime)%sonde(i,ie)%u%val)
          global_shared_buf(ie,15)= global_shared_buf(ie,15) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(2,ie)
       IF ((obs(itime)%sonde(i,ie)%v%val .NE. sdmiss) .AND. (an_obs(itime)%sonde(i,ie)%v%val .NE. sdmiss)) THEN
          global_shared_buf(ie,6) = global_shared_buf(ie,6)                           &
                                  +(obs(itime)%sonde(i,ie)%v%val - an_obs(itime)%sonde(i,ie)%v%val)
          global_shared_buf(ie,16)= global_shared_buf(ie,16) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(3,ie)
       IF ((obs(itime)%sonde(i,ie)%t%val .NE. sdmiss) .AND. (an_obs(itime)%sonde(i,ie)%t%val .NE. sdmiss)) THEN
          global_shared_buf(ie,7) = global_shared_buf(ie,7)                           &
                                  +(obs(itime)%sonde(i,ie)%t%val - an_obs(itime)%sonde(i,ie)%t%val)
          global_shared_buf(ie,17)= global_shared_buf(ie,17) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(4,ie)
       IF ((obs(itime)%sonde(i,ie)%q%val .NE. sdmiss) .AND. (an_obs(itime)%sonde(i,ie)%q%val .NE. sdmiss)) THEN
          global_shared_buf(ie,8) = global_shared_buf(ie,8)                           &
                                  +(obs(itime)%sonde(i,ie)%q%val - an_obs(itime)%sonde(i,ie)%q%val)
          global_shared_buf(ie,18)= global_shared_buf(ie,18) + 1.D0
       END IF
    END DO
    END DO  ! ie

    CALL Wrap_Repro_Sum(nvars=18, comm=hybrid%par%comm)

    ! O-B
    nobsdu   = global_shared_sum(11)
    nobsdv   = global_shared_sum(12)
    nobsdt   = global_shared_sum(13)
    nobsdq   = global_shared_sum(14)
    obmsdu   = global_shared_sum(1) /nobsdu
    obmsdv   = global_shared_sum(2) /nobsdv
    obmsdt   = global_shared_sum(3) /nobsdt
    obmsdq   = global_shared_sum(4) /nobsdq
    ! O-A
    noasdu   = global_shared_sum(15)
    noasdv   = global_shared_sum(16)
    noasdt   = global_shared_sum(17)
    noasdq   = global_shared_sum(18)
    oamsdu   = global_shared_sum(5) /noasdu
    oamsdv   = global_shared_sum(6) /noasdv
    oamsdt   = global_shared_sum(7) /noasdt
    oamsdq   = global_shared_sum(8) /noasdq

    END IF  ! da_sonde

!.. Mean of O-B and O-A for Surface observation
    IF (da_surface) THEN

    obmsfu= 0.D0; obmsfv= 0.D0; obmsft= 0.D0; obmsfq= 0.D0; obmsfps= 0.D0
    obdsfu= 0.D0; obdsfv= 0.D0; obdsft= 0.D0; obdsfq= 0.D0; obdsfps= 0.D0
    oamsfu= 0.D0; oamsfv= 0.D0; oamsft= 0.D0; oamsfq= 0.D0; oamsfps= 0.D0
    oadsfu= 0.D0; oadsfv= 0.D0; oadsft= 0.D0; oadsfq= 0.D0; oadsfps= 0.D0
    nobsfu= 0   ; nobsfv= 0   ; nobsft= 0   ; nobsfq= 0   ; nobsfps= 0
    noasfu= 0   ; noasfv= 0   ; noasft= 0   ; noasfq= 0   ; noasfps= 0

    global_shared_buf(:,1:20) = 0.D0

    DO ie = nets, nete
    ! background
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%u%val .NE. sdmiss) .AND. (bg_obs(itime)%surface(i,ie)%u%val .NE. sdmiss)) THEN
          global_shared_buf(ie,1) = global_shared_buf(ie,1)                               &
                                  +(obs(itime)%surface(i,ie)%u%val - bg_obs(itime)%surface(i,ie)%u%val)
          global_shared_buf(ie,11)= global_shared_buf(ie,11) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%v%val .NE. sdmiss) .AND. (bg_obs(itime)%surface(i,ie)%v%val .NE. sdmiss)) THEN
          global_shared_buf(ie,2) = global_shared_buf(ie,2)                               &
                                  +(obs(itime)%surface(i,ie)%v%val - bg_obs(itime)%surface(i,ie)%v%val)
          global_shared_buf(ie,12)= global_shared_buf(ie,12) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%t%val .NE. sdmiss) .AND. (bg_obs(itime)%surface(i,ie)%t%val .NE. sdmiss)) THEN
          global_shared_buf(ie,3) = global_shared_buf(ie,3)                               &
                                  +(obs(itime)%surface(i,ie)%t%val - bg_obs(itime)%surface(i,ie)%t%val)
          global_shared_buf(ie,13)= global_shared_buf(ie,13) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%q%val .NE. sdmiss) .AND. (bg_obs(itime)%surface(i,ie)%q%val .NE. sdmiss)) THEN
          global_shared_buf(ie,4) = global_shared_buf(ie,4)                               &
                                  +(obs(itime)%surface(i,ie)%q%val - bg_obs(itime)%surface(i,ie)%q%val)
          global_shared_buf(ie,14)= global_shared_buf(ie,14) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%p%val .NE. sdmiss) .AND. (bg_obs(itime)%surface(i,ie)%p%val .NE. sdmiss)) THEN
          global_shared_buf(ie,5) = global_shared_buf(ie,5)                               &
                                  +(obs(itime)%surface(i,ie)%p%val - bg_obs(itime)%surface(i,ie)%p%val)
          global_shared_buf(ie,15)= global_shared_buf(ie,15) + 1.D0
       END IF
    END DO

   ! analysis
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%u%val .NE. sdmiss) .AND. (an_obs(itime)%surface(i,ie)%u%val .NE. sdmiss)) THEN
          global_shared_buf(ie,6) = global_shared_buf(ie,6)                               &
                                  +(obs(itime)%surface(i,ie)%u%val - an_obs(itime)%surface(i,ie)%u%val)
          global_shared_buf(ie,16)= global_shared_buf(ie,16) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%v%val .NE. sdmiss) .AND. (an_obs(itime)%surface(i,ie)%v%val .NE. sdmiss)) THEN
          global_shared_buf(ie,7) = global_shared_buf(ie,7)                               &
                                  +(obs(itime)%surface(i,ie)%v%val - an_obs(itime)%surface(i,ie)%v%val)
          global_shared_buf(ie,17)= global_shared_buf(ie,17) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%t%val .NE. sdmiss) .AND. (an_obs(itime)%surface(i,ie)%t%val .NE. sdmiss)) THEN
          global_shared_buf(ie,8) = global_shared_buf(ie,8)                               &
                                  +(obs(itime)%surface(i,ie)%t%val - an_obs(itime)%surface(i,ie)%t%val)
          global_shared_buf(ie,18)= global_shared_buf(ie,18) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%q%val .NE. sdmiss) .AND. (an_obs(itime)%surface(i,ie)%q%val .NE. sdmiss)) THEN
          global_shared_buf(ie,9) = global_shared_buf(ie,9)                               &
                                  +(obs(itime)%surface(i,ie)%q%val - an_obs(itime)%surface(i,ie)%q%val)
          global_shared_buf(ie,19)= global_shared_buf(ie,19) + 1.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%p%val .NE. sdmiss) .AND. (an_obs(itime)%surface(i,ie)%p%val .NE. sdmiss)) THEN
          global_shared_buf(ie,10) = global_shared_buf(ie,10)                               &
                                  +(obs(itime)%surface(i,ie)%p%val - an_obs(itime)%surface(i,ie)%p%val)
          global_shared_buf(ie,20)= global_shared_buf(ie,20) + 1.D0
       END IF
    END DO
    END DO  ! ie

    CALL Wrap_Repro_Sum(nvars=20, comm=hybrid%par%comm)

    ! O-B
    nobsfu   = global_shared_sum(11)
    nobsfv   = global_shared_sum(12)
    nobsft   = global_shared_sum(13)
    nobsfq   = global_shared_sum(14)
    nobsfps  = global_shared_sum(15)
    obmsfu   = global_shared_sum(1) /nobsfu
    obmsfv   = global_shared_sum(2) /nobsfv
    obmsft   = global_shared_sum(3) /nobsft
    obmsfq   = global_shared_sum(4) /nobsfq
    obmsfps  = global_shared_sum(5) /nobsfps
    ! O-A
    noasfu   = global_shared_sum(16)
    noasfv   = global_shared_sum(17)
    noasft   = global_shared_sum(18)
    noasfq   = global_shared_sum(19)
    noasfps  = global_shared_sum(20)
    oamsfu   = global_shared_sum(6) /noasfu
    oamsfv   = global_shared_sum(7) /noasfv
    oamsft   = global_shared_sum(8) /noasft
    oamsfq   = global_shared_sum(9) /noasfq
    oamsfps  = global_shared_sum(10) /noasfps

    END IF  ! da_surface

    IF (da_bogus) THEN

    obmbgps= 0.D0
    obdbgps= 0.D0
    oambgps= 0.D0
    oadbgps= 0.D0
    nobbgps= 0
    noabgps= 0

    global_shared_buf(:,1:4) = 0.D0

    DO ie = nets, nete
    ! background
    DO i  = 1, obsmet(itime)%nbg_e(ie)
       IF ((obs(itime)%bogus(i,ie)%p%val .NE. sdmiss) .AND. (bg_obs(itime)%bogus(i,ie)%p%val .NE. sdmiss)) THEN
          global_shared_buf(ie,1) = global_shared_buf(ie,1)                               &
                                  +(obs(itime)%bogus(i,ie)%p%val - bg_obs(itime)%bogus(i,ie)%p%val)
          global_shared_buf(ie,3)= global_shared_buf(ie,3) + 1.D0
       END IF
    END DO

   ! analysis
    DO i  = 1, obsmet(itime)%nbg_e(ie)
       IF ((obs(itime)%bogus(i,ie)%p%val .NE. sdmiss) .AND. (an_obs(itime)%bogus(i,ie)%p%val .NE. sdmiss)) THEN
          global_shared_buf(ie,2) = global_shared_buf(ie,2)                               &
                                  +(obs(itime)%bogus(i,ie)%p%val - an_obs(itime)%bogus(i,ie)%p%val)
          global_shared_buf(ie,4)= global_shared_buf(ie,4) + 1.D0
       END IF
    END DO
    END DO

    CALL Wrap_Repro_Sum(nvars=4, comm=hybrid%par%comm)

    ! O-B
    nobbgps  = global_shared_sum(3)
    obmbgps  = global_shared_sum(1) /nobbgps
    ! O-A
    noabgps  = global_shared_sum(4)
    oambgps  = global_shared_sum(2) /noabgps

    END IF ! da_bogus

!.. Print out the mean of O-B and O-A
    IF (da_sonde) THEN
    IF (par%isMasterProc) WRITE(iulog,*) "Mean of O-B and O-A"
    IF (par%isMasterProc) WRITE(iulog,*) "itime= ", itime
    IF (par%isMasterProc) WRITE(iulog,*) "Against Sonde Observation"
    IF (par%isMasterProc) WRITE(iulog,*) "O-B U =", obmsdu
    IF (par%isMasterProc) WRITE(iulog,*) "O-A U =", oamsdu
    IF (par%isMasterProc) WRITE(iulog,*) "O-B V =", obmsdv
    IF (par%isMasterProc) WRITE(iulog,*) "O-A V =", oamsdv
    IF (par%isMasterProc) WRITE(iulog,*) "O-B T =", obmsdt
    IF (par%isMasterProc) WRITE(iulog,*) "O-A T =", oamsdt
    IF (par%isMasterProc) WRITE(iulog,*) "O-B Q =", obmsdq
    IF (par%isMasterProc) WRITE(iulog,*) "O-A Q =", oamsdq
    END IF
    IF (da_surface) THEN
    IF (par%isMasterProc) WRITE(iulog,*) "Agaist Surface Observation"
    IF (par%isMasterProc) WRITE(iulog,*) "O-B U  =", obmsfu
    IF (par%isMasterProc) WRITE(iulog,*) "O-A U  =", oamsfu
    IF (par%isMasterProc) WRITE(iulog,*) "O-B V  =", obmsfv
    IF (par%isMasterProc) WRITE(iulog,*) "O-A V  =", oamsfv
    IF (par%isMasterProc) WRITE(iulog,*) "O-B T  =", obmsft
    IF (par%isMasterProc) WRITE(iulog,*) "O-A T  =", oamsft
    IF (par%isMasterProc) WRITE(iulog,*) "O-B Q  =", obmsfq
    IF (par%isMasterProc) WRITE(iulog,*) "O-A Q  =", oamsfq
    IF (par%isMasterProc) WRITE(iulog,*) "O-B PS =", obmsfps
    IF (par%isMasterProc) WRITE(iulog,*) "O-A PS =", oamsfps
    END IF
    IF (da_bogus) THEN
    IF (par%isMasterProc) WRITE(iulog,*) "Agaist Surface Bogus Observation"
    IF (par%isMasterProc) WRITE(iulog,*) "O-B PS =", obmbgps
    IF (par%isMasterProc) WRITE(iulog,*) "O-A PS =", oambgps
    END IF

!.. Standard deviation of O-B and O-A for Sonde
    IF (da_sonde) THEN

    global_shared_buf(:,1:8) = 0.D0

    DO ie = nets, nete
   ! background
    DO i  = 1, obsmet(itime)%nsd_e(1,ie)
       IF ((obs(itime)%sonde(i,ie)%u%val .NE. sdmiss) .AND. (bg_obs(itime)%sonde(i,ie)%u%val .NE. sdmiss)) THEN
          global_shared_buf(ie,1) = global_shared_buf(ie,1)                                         &
                                  +((obs(itime)%sonde(i,ie)%u%val -bg_obs(itime)%sonde(i,ie)%u%val)-obmsdu)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(2,ie)
       IF ((obs(itime)%sonde(i,ie)%v%val .NE. sdmiss) .AND. (bg_obs(itime)%sonde(i,ie)%v%val .NE. sdmiss)) THEN
          global_shared_buf(ie,2) = global_shared_buf(ie,2)                                         &
                                  +((obs(itime)%sonde(i,ie)%v%val -bg_obs(itime)%sonde(i,ie)%v%val)-obmsdv)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(3,ie)
       IF ((obs(itime)%sonde(i,ie)%t%val .NE. sdmiss) .AND. (bg_obs(itime)%sonde(i,ie)%t%val .NE. sdmiss)) THEN
          global_shared_buf(ie,3) = global_shared_buf(ie,3)                                         &
                                  +((obs(itime)%sonde(i,ie)%t%val -bg_obs(itime)%sonde(i,ie)%t%val)-obmsdt)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(4,ie)
       IF ((obs(itime)%sonde(i,ie)%q%val .NE. sdmiss) .AND. (bg_obs(itime)%sonde(i,ie)%q%val .NE. sdmiss)) THEN
          global_shared_buf(ie,4) = global_shared_buf(ie,4)                                         &
                                  +((obs(itime)%sonde(i,ie)%q%val -bg_obs(itime)%sonde(i,ie)%q%val)-obmsdq)**2.D0
       END IF
    END DO
   ! analysis
    DO i  = 1, obsmet(itime)%nsd_e(1,ie)
       IF ((obs(itime)%sonde(i,ie)%u%val .NE. sdmiss) .AND. (an_obs(itime)%sonde(i,ie)%u%val .NE. sdmiss)) THEN
          global_shared_buf(ie,5) = global_shared_buf(ie,5)                                         &
                                  +((obs(itime)%sonde(i,ie)%u%val -an_obs(itime)%sonde(i,ie)%u%val)-oamsdu)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(2,ie)
       IF ((obs(itime)%sonde(i,ie)%v%val .NE. sdmiss) .AND. (an_obs(itime)%sonde(i,ie)%v%val .NE. sdmiss)) THEN
          global_shared_buf(ie,6) = global_shared_buf(ie,6)                                         &
                                  +((obs(itime)%sonde(i,ie)%v%val -an_obs(itime)%sonde(i,ie)%v%val)-oamsdv)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(3,ie)
       IF ((obs(itime)%sonde(i,ie)%t%val .NE. sdmiss) .AND. (an_obs(itime)%sonde(i,ie)%t%val .NE. sdmiss)) THEN
          global_shared_buf(ie,7) = global_shared_buf(ie,7)                                         &
                                  +((obs(itime)%sonde(i,ie)%t%val -an_obs(itime)%sonde(i,ie)%t%val)-oamsdt)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsd_e(4,ie)
       IF ((obs(itime)%sonde(i,ie)%q%val .NE. sdmiss) .AND. (an_obs(itime)%sonde(i,ie)%q%val .NE. sdmiss)) THEN
          global_shared_buf(ie,8) = global_shared_buf(ie,8)                                         &
                                  +((obs(itime)%sonde(i,ie)%q%val -an_obs(itime)%sonde(i,ie)%q%val)-oamsdq)**2.D0
       END IF
    END DO
    END DO    ! ie

    CALL Wrap_Repro_Sum(nvars=8, comm=hybrid%par%comm)

    ! O-B
    obdsdu   = sqrt(global_shared_sum(1) /nobsdu)
    obdsdv   = sqrt(global_shared_sum(2) /nobsdv)
    obdsdt   = sqrt(global_shared_sum(3) /nobsdt)
    obdsdq   = sqrt(global_shared_sum(4) /nobsdq)
    ! O-A
    oadsdu   = sqrt(global_shared_sum(5) /noasdu)
    oadsdv   = sqrt(global_shared_sum(6) /noasdv)
    oadsdt   = sqrt(global_shared_sum(7) /noasdt)
    oadsdq   = sqrt(global_shared_sum(8) /noasdq)

    END IF    ! da_sonde

!.. Standard deviation of O-B and O-A for Surface
    IF (da_surface) THEN

    global_shared_buf(:,1:10) = 0.D0

    DO ie = nets, nete
    ! background
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%u%val .NE. sdmiss) .AND. (bg_obs(itime)%surface(i,ie)%u%val .NE. sdmiss)) THEN
          global_shared_buf(ie,1) = global_shared_buf(ie,1)                                             &
                                  +((obs(itime)%surface(i,ie)%u%val -bg_obs(itime)%surface(i,ie)%u%val)-obmsfu)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%v%val .NE. sdmiss) .AND. (bg_obs(itime)%surface(i,ie)%v%val .NE. sdmiss)) THEN
          global_shared_buf(ie,2) = global_shared_buf(ie,2)                                             &
                                  +((obs(itime)%surface(i,ie)%v%val -bg_obs(itime)%surface(i,ie)%v%val)-obmsfv)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%t%val .NE. sdmiss) .AND. (bg_obs(itime)%surface(i,ie)%t%val .NE. sdmiss)) THEN
          global_shared_buf(ie,3) = global_shared_buf(ie,3)                                             &
                                  +((obs(itime)%surface(i,ie)%t%val -bg_obs(itime)%surface(i,ie)%t%val)-obmsft)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%q%val .NE. sdmiss) .AND. (bg_obs(itime)%surface(i,ie)%q%val .NE. sdmiss)) THEN
          global_shared_buf(ie,4) = global_shared_buf(ie,4)                                             &
                                  +((obs(itime)%surface(i,ie)%q%val -bg_obs(itime)%surface(i,ie)%q%val)-obmsfq)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%p%val .NE. sdmiss) .AND. (bg_obs(itime)%surface(i,ie)%p%val .NE. sdmiss)) THEN
          global_shared_buf(ie,5) = global_shared_buf(ie,5)                                             &
                                  +((obs(itime)%surface(i,ie)%p%val -bg_obs(itime)%surface(i,ie)%p%val)-obmsfps)**2.D0
       END IF
    END DO
    ! analysis
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%u%val .NE. sdmiss) .AND. (an_obs(itime)%surface(i,ie)%u%val .NE. sdmiss)) THEN
          global_shared_buf(ie,6) = global_shared_buf(ie,6)                                             &
                                  +((obs(itime)%surface(i,ie)%u%val -an_obs(itime)%surface(i,ie)%u%val)-oamsfu)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%v%val .NE. sdmiss) .AND. (an_obs(itime)%surface(i,ie)%v%val .NE. sdmiss)) THEN
          global_shared_buf(ie,7) = global_shared_buf(ie,7)                                             &
                                  +((obs(itime)%surface(i,ie)%v%val -an_obs(itime)%surface(i,ie)%v%val)-oamsfv)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%t%val .NE. sdmiss) .AND. (an_obs(itime)%surface(i,ie)%t%val .NE. sdmiss)) THEN
          global_shared_buf(ie,8) = global_shared_buf(ie,8)                                             &
                                  +((obs(itime)%surface(i,ie)%t%val -an_obs(itime)%surface(i,ie)%t%val)-oamsft)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%q%val .NE. sdmiss) .AND. (an_obs(itime)%surface(i,ie)%q%val .NE. sdmiss)) THEN
          global_shared_buf(ie,9) = global_shared_buf(ie,9)                                             &
                                  +((obs(itime)%surface(i,ie)%q%val -an_obs(itime)%surface(i,ie)%q%val)-oamsfq)**2.D0
       END IF
    END DO
    DO i  = 1, obsmet(itime)%nsf_e(ie)
       IF ((obs(itime)%surface(i,ie)%p%val .NE. sdmiss) .AND. (an_obs(itime)%surface(i,ie)%p%val .NE. sdmiss)) THEN
          global_shared_buf(ie,10) = global_shared_buf(ie,10)                                              &
                                  +((obs(itime)%surface(i,ie)%p%val -an_obs(itime)%surface(i,ie)%p%val)-oamsfps)**2.D0
       END IF
    END DO
    END DO  ! ie

    CALL Wrap_Repro_Sum(nvars=10, comm=hybrid%par%comm)

    ! O-B
    obdsfu   = sqrt(global_shared_sum(1) /nobsfu)
    obdsfv   = sqrt(global_shared_sum(2) /nobsfv)
    obdsft   = sqrt(global_shared_sum(3) /nobsft)
    obdsfq   = sqrt(global_shared_sum(4) /nobsfq)
    obdsfps  = sqrt(global_shared_sum(5) /nobsfps)
    ! O-A
    oadsfu   = sqrt(global_shared_sum(6) /noasfu)
    oadsfv   = sqrt(global_shared_sum(7) /noasfv)
    oadsft   = sqrt(global_shared_sum(8) /noasft)
    oadsfq   = sqrt(global_shared_sum(9) /noasfq)
    oadsfps  = sqrt(global_shared_sum(10) /noasfps)

    END IF   ! da_surface

!.. Standard deviation of O-B and O-A for Bogus
    IF (da_bogus) THEN

    global_shared_buf(:,1:2) = 0.D0

    DO ie = nets, nete
    ! background
    DO i  = 1, obsmet(itime)%nbg_e(ie)
       IF ((obs(itime)%bogus(i,ie)%p%val .NE. sdmiss) .AND. (bg_obs(itime)%bogus(i,ie)%p%val .NE. sdmiss)) THEN
          global_shared_buf(ie,1) = global_shared_buf(ie,1)                                             &
                                  +((obs(itime)%bogus(i,ie)%p%val -bg_obs(itime)%bogus(i,ie)%p%val)-obmbgps)**2.D0
       END IF
    END DO
    ! analysis
    DO i  = 1, obsmet(itime)%nbg_e(ie)
       IF ((obs(itime)%bogus(i,ie)%p%val .NE. sdmiss) .AND. (an_obs(itime)%bogus(i,ie)%p%val .NE. sdmiss)) THEN
          global_shared_buf(ie,2) = global_shared_buf(ie,2)                                              &
                                  +((obs(itime)%bogus(i,ie)%p%val -an_obs(itime)%bogus(i,ie)%p%val)-oambgps)**2.D0
       END IF
    END DO  ! da_bogus
    END DO  ! ie

    CALL Wrap_Repro_Sum(nvars=2, comm=hybrid%par%comm)

    ! O-B
    obdbgps   = sqrt(global_shared_sum(1) /nobbgps)
    ! O-A
    oadbgps   = sqrt(global_shared_sum(2) /noabgps)    ! O-A

    END IF   ! da_bogus


!.. Print out the standard deviation of O-B and O-A

    IF (da_sonde) THEN
    IF (par%isMasterProc) WRITE(iulog,*) "Standard deviation of O-B and O-A"
    IF (par%isMasterProc) WRITE(iulog,*) "Against Sonde Observation"
    IF (par%isMasterProc) WRITE(iulog,*) "O-B U =", obdsdu, "#obs=",nobsdu
    IF (par%isMasterProc) WRITE(iulog,*) "O-A U =", oadsdu, "#obs=",noasdu
    IF (par%isMasterProc) WRITE(iulog,*) "O-B V =", obdsdv, "#obs=",nobsdv
    IF (par%isMasterProc) WRITE(iulog,*) "O-A V =", oadsdv, "#obs=",noasdv
    IF (par%isMasterProc) WRITE(iulog,*) "O-B T =", obdsdt, "#obs=",nobsdt
    IF (par%isMasterProc) WRITE(iulog,*) "O-A T =", oadsdt, "#obs=",noasdt
    IF (par%isMasterProc) WRITE(iulog,*) "O-B Q =", obdsdq, "#obs=",nobsdq
    IF (par%isMasterProc) WRITE(iulog,*) "O-A Q =", oadsdq, "#obs=",noasdq
    END IF
    IF (da_surface) THEN
    IF (par%isMasterProc) WRITE(iulog,*) "Against Surface Observation"
    IF (par%isMasterProc) WRITE(iulog,*) "O-B U  =", obdsfu, "#obs=",nobsfu
    IF (par%isMasterProc) WRITE(iulog,*) "O-A U  =", oadsfu, "#obs=",noasfu
    IF (par%isMasterProc) WRITE(iulog,*) "O-B V  =", obdsfv, "#obs=",nobsfv
    IF (par%isMasterProc) WRITE(iulog,*) "O-A V  =", oadsfv, "#obs=",noasfv
    IF (par%isMasterProc) WRITE(iulog,*) "O-B T  =", obdsft, "#obs=",nobsft
    IF (par%isMasterProc) WRITE(iulog,*) "O-A T  =", oadsft, "#obs=",noasft
    IF (par%isMasterProc) WRITE(iulog,*) "O-B Q  =", obdsfq, "#obs=",nobsfq
    IF (par%isMasterProc) WRITE(iulog,*) "O-A Q  =", oadsfq, "#obs=",noasfq
    IF (par%isMasterProc) WRITE(iulog,*) "O-B PS =", obdsfps, "#obs=",nobsfps
    IF (par%isMasterProc) WRITE(iulog,*) "O-A PS =", oadsfps, "#obs=",noasfps
    END IF
    IF (da_bogus) THEN
    IF (par%isMasterProc) WRITE(iulog,*) "Against Surface Bogus Observation"
    IF (par%isMasterProc) WRITE(iulog,*) "O-B PS =", obdbgps, "#obs=",nobbgps
    IF (par%isMasterProc) WRITE(iulog,*) "O-A PS =", oadbgps, "#obs=",noabgps
    END IF

!    DO k = 0, par%nprocs-1
!    IF( par%iproc == k )THEN
!    DO ie= nets, nete
!    DO i = 1, nsf_e(ie)
!       IF( obs%surface(i,ie)%u%val .NE. sdmiss )THEN
!          WRITE(iulog,'(3I5,6F10.4)') par%iproc,i,ie,   &          
!                         obsloc%surface%all(ie)%lat(i), &
!                         obsloc%surface%all(ie)%lon(i), &
!                         obsloc%surface%all(ie)%hgt(i), &
!                         obs%surface(i,ie)%u%val,       &
!                         obs%surface(i,ie)%u%val        &
!                         -bg_obs%surface(i,ie)%u%val,   &
!                         obs%surface(i,ie)%u%val        &
!                         -an_obs%surface(i,ie)%u%val
!       ENDIF
!    ENDDO
!    ENDDO
!    ENDIF
!       CALL Mpi_Barrier(par%comm, ierr)
!    ENDDO
!
!    DO k = 0, par%nprocs-1
!    IF( par%iproc == k )THEN
!    DO ie= nets, nete
!    DO i = 1, nsf_e(ie)
!       IF( obs%surface(i,ie)%t%val .NE. sdmiss )THEN
!          WRITE(iulog,'(3I5,6F10.4)') par%iproc,i,ie,   &
!                         obsloc%surface%all(ie)%lat(i), &
!                         obsloc%surface%all(ie)%lon(i), &
!                         obsloc%surface%all(ie)%hgt(i), &
!                         obs%surface(i,ie)%t%val,       &
!                         obs%surface(i,ie)%t%val        &
!                         -bg_obs%surface(i,ie)%t%val,   &
!                         obs%surface(i,ie)%t%val        &
!                         -an_obs%surface(i,ie)%t%val
!       ENDIF
!    ENDDO
!    ENDDO
!    ENDIF
!       CALL Mpi_Barrier(par%comm, ierr)
!    ENDDO
!
!    DO k = 0, par%nprocs-1
!    IF( par%iproc == k )THEN
!    DO ie= nets, nete
!    DO i = 1, nsf_e(ie)
!       IF( obs%surface(i,ie)%p%val .NE. sdmiss )THEN
!          WRITE(iulog,'(3I5,3F10.4,3F10.2)') par%iproc,i,ie, &
!                         obsloc%surface%all(ie)%lat(i),      &
!                         obsloc%surface%all(ie)%lon(i),      &
!                         obsloc%surface%all(ie)%hgt(i),      &
!                         obs%surface(i,ie)%p%val,            &
!                         obs%surface(i,ie)%p%val             &
!                         -bg_obs%surface(i,ie)%p%val,        &
!                         obs%surface(i,ie)%p%val             &
!                         -an_obs%surface(i,ie)%p%val
!       ENDIF
!    ENDDO
!    ENDDO
!    ENDIF
!       CALL Mpi_Barrier(par%comm, ierr)
!    ENDDO
!       CALL Mpi_Barrier(par%comm, ierr)
!
    RETURN

  END SUBROUTINE Check_OB
!===============================================================================
!===============================================================================
  SUBROUTINE Check_Adj(gs_md,       &
                       b_eig,       &
                       b_eig_chi,   &
                       b_eig_Tps,   &
                       b_eig_q,     &
                       b_eig_a)
  ! 2. Variables
  ! 2-1. In/Output
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime), &
       INTENT(IN)  :: gs_md   

    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc), &
       INTENT(IN) :: b_eig
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc), &
       INTENT(IN) :: b_eig_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc), &
       INTENT(IN) :: b_eig_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc), &
       INTENT(IN) :: b_eig_q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl), &
       INTENT(IN) :: b_eig_a

  ! 2-2. Local
!   REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
!      :: tl_md
    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc) &
       :: ad_eig
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc) &
       :: ad_eig_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc) &
       :: ad_eig_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc) &
       :: ad_eig_q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl) &
       :: ad_eig_a

    REAL(KIND=real_kind) &
       :: tmp, yTy
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie,     &
          wn, vr,          &
          itime

  ! 1. Eigen space -> observation space
    ierr = TimingStart('TlmEigToObs_Check_Adj')   !--- profiling
    CALL TlmEigToObs(gs_md(:,:,:,:,:), &
                     b_eig,            &
                     b_eig_chi,        &
                     b_eig_Tps,        &
                     b_eig_q,          &
                     b_eig_a,          &
                     tl_obs(:))
    ierr = TimingStop('TlmEigToObs_Check_Adj')   !--- profiling

  DO itime = 1, ntime
    ! single observation
    IF( singleobs )THEN
           CALL MultiplyRSingle(itime,tl_obs(itime)%single%u%val, tl_obs(itime)%single%v%val, &
                                      tl_obs(itime)%single%t%val, tl_obs(itime)%single%q%val, &
                                      tl_obs(itime)%single%ps%val)
    END IF
    ! sonde observation
    IF (da_sonde) THEN
           CALL MultiplyRSonde(itime,tl_obs(itime)%sonde%u%val, tl_obs(itime)%sonde%v%val, &
                               tl_obs(itime)%sonde%t%val, tl_obs(itime)%sonde%q%val)
    END IF
    ! surface observation  
    IF (da_surface) THEN
           CALL MultiplyRSurface(itime,tl_obs(itime)%surface%u%val, tl_obs(itime)%surface%v%val, &
                                 tl_obs(itime)%surface%t%val, tl_obs(itime)%surface%q%val, tl_obs(itime)%surface%p%val)
    END IF
    ! bogus observation  
    IF (da_bogus) THEN
           CALL MultiplyRBogus(itime,tl_obs(itime)%bogus%p%val)
    END IF
    ! aircraft observation  
    IF (da_aircraft) THEN
           CALL MultiplyRAircraft(itime,tl_obs(itime)%aircraft%u%val, tl_obs(itime)%aircraft%v%val, &
                                  tl_obs(itime)%aircraft%t%val)
    END IF
      ! amv observation  
    IF (da_amv) THEN
           CALL MultiplyRAmv(itime,tl_obs(itime)%amv%u%val, tl_obs(itime)%amv%v%val)
    END IF
    ! scatwind observation  
    IF (da_scatwind) THEN
           CALL MultiplyRScatwind(itime,tl_obs(itime)%scatwind%u%val, tl_obs(itime)%scatwind%v%val)
    END IF
    ! amsua observation  
    IF (da_amsua) THEN
           CALL MultiplyRAmsua(itime,tl_obs(itime)%amsua%t%val)
    END IF
    ! gpsro observation  
    IF (da_gpsro) THEN
            CALL MultiplyRGpsro(itime,tl_obs(itime)%gpsro%t%val, tl_obs(itime)%gpsro%q%val)
    END IF
    ! iasi observation  
    IF (da_iasi) THEN
           CALL MultiplyRIasi(itime,tl_obs(itime)%iasi%t%val, tl_obs(itime)%iasi%q%val)
    END IF
    ! cris observation  
    IF (da_cris) THEN
           CALL MultiplyRCris(itime,tl_obs(itime)%cris%t%val,tl_obs(itime)%cris%q%val)
    END IF
    ! atms observation  
    IF (da_atms) THEN
           CALL MultiplyRAtms(itime,tl_obs(itime)%atms%t%val)
    END IF
    ! atmswv observation  
    IF (da_atmswv) THEN
           CALL MultiplyRAtmswv(itime,tl_obs(itime)%atmswv%t%val,tl_obs(itime)%atmswv%q%val)
    END IF
    ! mhs observation  
    IF (da_mhs) THEN
           CALL MultiplyRMhs(itime,tl_obs(itime)%mhs%t%val, tl_obs(itime)%mhs%q%val)
    END IF
    ! csr observation  
    IF (da_csr) THEN
           CALL MultiplyRCsr(itime,tl_obs(itime)%csr%t%val, tl_obs(itime)%csr%q%val)
    END IF


  ! 2. Adjoint of eigen space -> observation space
    ! Single observation
    IF( singleobs )THEN
       ad_obs(itime)%single  = tl_obs(itime)%single
           CALL MultiplyRSingle(itime,ad_obs(itime)%single%u%val, ad_obs(itime)%single%v%val, &
                                      ad_obs(itime)%single%t%val, ad_obs(itime)%single%q%val, &
                                      ad_obs(itime)%single%ps%val                            )
    END IF
    ! sonde observation
    IF (da_sonde) THEN
       ad_obs(itime)%sonde   = tl_obs(itime)%sonde
           CALL MultiplyRSonde(itime,ad_obs(itime)%sonde%u%val, ad_obs(itime)%sonde%v%val, &
                               ad_obs(itime)%sonde%t%val, ad_obs(itime)%sonde%q%val)
    END IF
    ! surface observation
    IF (da_surface) THEN
       ad_obs(itime)%surface   = tl_obs(itime)%surface
           CALL MultiplyRSurface(itime,ad_obs(itime)%surface%u%val, ad_obs(itime)%surface%v%val, &
                                 ad_obs(itime)%surface%t%val, ad_obs(itime)%surface%q%val, ad_obs(itime)%surface%p%val)
    END IF
    ! bogus observation
    IF (da_bogus) THEN
       ad_obs(itime)%bogus   = tl_obs(itime)%bogus
           CALL MultiplyRBogus(itime,ad_obs(itime)%bogus%p%val)
    END IF
    ! aircraft observation
    IF (da_aircraft) THEN
       ad_obs(itime)%aircraft  = tl_obs(itime)%aircraft
           CALL MultiplyRAircraft(itime,ad_obs(itime)%aircraft%u%val, ad_obs(itime)%aircraft%v%val, &
                                  ad_obs(itime)%aircraft%t%val )
    END IF
    ! amv observation
    IF (da_amv) THEN
       ad_obs(itime)%amv  = tl_obs(itime)%amv
            CALL MultiplyRAmv(itime,ad_obs(itime)%amv%u%val, ad_obs(itime)%amv%v%val)
    END IF
    ! scatwind observation
    IF (da_scatwind) THEN
       ad_obs(itime)%scatwind   = tl_obs(itime)%scatwind
           CALL MultiplyRScatwind(itime,ad_obs(itime)%scatwind%u%val, ad_obs(itime)%scatwind%v%val)
    END IF
    ! amsua observation  
    IF (da_amsua) THEN
       ad_obs(itime)%amsua  = tl_obs(itime)%amsua
           CALL MultiplyRAmsua(itime,ad_obs(itime)%amsua%t%val)
    END IF
    ! gpsro observation  
    IF (da_gpsro) THEN
       ad_obs(itime)%gpsro  = tl_obs(itime)%gpsro
            CALL MultiplyRGpsro(itime,ad_obs(itime)%gpsro%t%val, ad_obs(itime)%gpsro%q%val)
    END IF
    ! iasi observation  
    IF (da_iasi) THEN
       ad_obs(itime)%iasi  = tl_obs(itime)%iasi
           CALL MultiplyRIasi(itime,ad_obs(itime)%iasi%t%val, ad_obs(itime)%iasi%q%val)
    END IF
    ! cris observation  
    IF (da_cris) THEN
       ad_obs(itime)%cris  = tl_obs(itime)%cris
           CALL MultiplyRCris(itime,ad_obs(itime)%cris%t%val,ad_obs(itime)%cris%q%val)
    END IF
    ! atms observation  
    IF (da_atms) THEN
       ad_obs(itime)%atms  = tl_obs(itime)%atms
           CALL MultiplyRAtms(itime,ad_obs(itime)%atms%t%val)
    END IF
    ! atmswv observation  
    IF (da_atmswv) THEN
       ad_obs(itime)%atmswv  = tl_obs(itime)%atmswv
           CALL MultiplyRAtmswv(itime,ad_obs(itime)%atmswv%t%val,ad_obs(itime)%atmswv%q%val)
    END IF
    ! mhs observation  
    IF (da_mhs) THEN
       ad_obs(itime)%mhs  = tl_obs(itime)%mhs
           CALL MultiplyRMhs(itime,ad_obs(itime)%mhs%t%val, ad_obs(itime)%mhs%q%val)
    END IF
    ! csr observation  
    IF (da_csr) THEN
       ad_obs(itime)%csr  = tl_obs(itime)%csr
           CALL MultiplyRCsr(itime,ad_obs(itime)%csr%t%val, ad_obs(itime)%csr%q%val)
    END IF

  END DO ! itime   

       ierr = TimingStart('AdjEigToObs_Check_Adj')   !--- profiling
       ad_eig = 0.D0
       ad_eig_chi = 0.D0
       ad_eig_Tps = 0.D0
       ad_eig_q = 0.D0
       ad_eig_a = 0.D0
       CALL AdjEigToObs(gs_md(:,:,:,:,:), &
                        ad_obs(:),    &
                        ad_eig,       &
                        ad_eig_chi,   &
                        ad_eig_Tps,   &
                        ad_eig_q,     &
                        ad_eig_a)
       ierr = TimingStop('AdjEigToObs_Check_Adj')   !--- profiling

  ! 3. Evaluation
    global_shared_buf(:,1) = 0.D0

  DO itime = 1, ntime 
    ! single observation test
    IF( singleobs )THEN
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%nsg_e(ie)
             IF (tl_obs(itime)%single(i,ie)%u%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%single(i,ie)%u%val**2.D0
             END IF
             IF (tl_obs(itime)%single(i,ie)%v%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%single(i,ie)%v%val**2.D0
             END IF
             IF (tl_obs(itime)%single(i,ie)%t%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%single(i,ie)%t%val**2.D0
             END IF
             IF (tl_obs(itime)%single(i,ie)%q%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%single(i,ie)%q%val**2.D0
             END IF
             IF (tl_obs(itime)%single(i,ie)%ps%val .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%single(i,ie)%ps%val**2.D0
             END IF
       END DO
       END DO
    END IF
    ! sonde observation
    IF (da_sonde) THEN
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%nsd_e(1,ie)
             IF (tl_obs(itime)%sonde(i,ie)%u%val .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%sonde(i,ie)%u%val**2.D0
             END IF
       END DO
       DO i  = 1, obsmet(itime)%nsd_e(2,ie)
             IF (tl_obs(itime)%sonde(i,ie)%v%val .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%sonde(i,ie)%v%val**2.D0
             END IF
       END DO
       DO i  = 1, obsmet(itime)%nsd_e(3,ie)
             IF (tl_obs(itime)%sonde(i,ie)%t%val .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%sonde(i,ie)%t%val**2.D0
             END IF
       END DO
       DO i  = 1, obsmet(itime)%nsd_e(4,ie)
             IF (tl_obs(itime)%sonde(i,ie)%q%val .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%sonde(i,ie)%q%val**2.D0
             END IF
       END DO
       END DO
    END IF
    ! surface observation
    IF (da_surface) THEN
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%nsf_e(ie)
             IF (tl_obs(itime)%surface(i,ie)%u%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%surface(i,ie)%u%val**2.D0
             END IF
             IF (tl_obs(itime)%surface(i,ie)%v%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%surface(i,ie)%v%val**2.D0
             END IF
             IF (tl_obs(itime)%surface(i,ie)%t%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%surface(i,ie)%t%val**2.D0
             END IF
             IF (tl_obs(itime)%surface(i,ie)%q%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%surface(i,ie)%q%val**2.D0
             END IF
             IF (tl_obs(itime)%surface(i,ie)%p%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%surface(i,ie)%p%val**2.D0
             END IF
       END DO
       END DO
    END IF
    ! bogus observation
    IF (da_bogus) THEN
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%nbg_e(ie)
             IF (tl_obs(itime)%bogus(i,ie)%p%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%bogus(i,ie)%p%val**2.D0
             END IF
       END DO
       END DO
    END IF
    ! aircraft observation
    IF (da_aircraft) THEN
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%nar_e(1,ie)
             IF (tl_obs(itime)%aircraft(i,ie)%u%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%aircraft(i,ie)%u%val**2.D0
             END IF
       END DO
       DO i  = 1, obsmet(itime)%nar_e(2,ie)
             IF (tl_obs(itime)%aircraft(i,ie)%v%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%aircraft(i,ie)%v%val**2.D0
             END IF
       END DO
       DO i  = 1, obsmet(itime)%nar_e(3,ie)
             IF (tl_obs(itime)%aircraft(i,ie)%t%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%aircraft(i,ie)%t%val**2.D0
             END IF
       END DO
       END DO
    END IF
    ! amv observation
    IF (da_amv) THEN
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%namv_e(ie)
             IF (tl_obs(itime)%amv(i,ie)%u%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%amv(i,ie)%u%val**2.D0
             END IF
             IF (tl_obs(itime)%amv(i,ie)%v%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%amv(i,ie)%v%val**2.D0
             END IF
       END DO
       END DO
    END IF
    ! scatwind observation
    IF (da_scatwind) THEN
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%nscatwind_e(ie)
             IF (tl_obs(itime)%scatwind(i,ie)%u%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%scatwind(i,ie)%u%val**2.D0
             END IF
             IF (tl_obs(itime)%scatwind(i,ie)%v%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + tl_obs(itime)%scatwind(i,ie)%v%val**2.D0
             END IF
       END DO
       END DO
    END IF
    ! amsua observation
    IF (da_amsua) THEN
       DO j  = 1, obsmet(itime)%namsuach
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%namsua_e(ie)
             IF (tl_obs(itime)%amsua(i,j,ie)%t%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%amsua(i,j,ie)%t%val**2.D0
             END IF
       END DO
       END DO
       END DO
    END IF
    ! gpsro observation
    IF (da_gpsro) THEN
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%ngpsro_e(ie)
             IF (tl_obs(itime)%gpsro(i,ie)%t%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%gpsro(i,ie)%t%val**2.D0
             END IF
             IF (tl_obs(itime)%gpsro(i,ie)%q%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%gpsro(i,ie)%q%val**2.D0
             END IF
       END DO
       END DO
     END IF
    ! iasi observation
    IF (da_iasi) THEN
       DO j = 1, obsmet(itime)%niasich
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%niasi_e(ie)
             IF (tl_obs(itime)%iasi(i,j,ie)%t%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%iasi(i,j,ie)%t%val**2.D0
             END IF
             IF (tl_obs(itime)%iasi(i,j,ie)%q%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%iasi(i,j,ie)%q%val**2.D0
             END IF
       END DO
       END DO
       END DO
     END IF
    ! cris observation
    IF (da_cris) THEN
       DO j = 1, obsmet(itime)%ncrisch
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%ncris_e(ie)
             IF (tl_obs(itime)%cris(i,j,ie)%t%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%cris(i,j,ie)%t%val**2.D0
             END IF
             IF (tl_obs(itime)%cris(i,j,ie)%q%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%cris(i,j,ie)%q%val**2.D0
             END IF
       END DO
       END DO
       END DO
     END IF
    ! atms observation
    IF (da_atms) THEN
       DO j  = 1, obsmet(itime)%natmsch
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%natms_e(ie)
             IF (tl_obs(itime)%atms(i,j,ie)%t%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%atms(i,j,ie)%t%val**2.D0
             END IF
       END DO
       END DO
       END DO
    END IF
    ! atmswv observation
    IF (da_atmswv) THEN
       DO j  = 1, obsmet(itime)%natmswvch
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%natmswv_e(ie)
             IF (tl_obs(itime)%atmswv(i,j,ie)%t%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%atmswv(i,j,ie)%t%val**2.D0
             END IF
             IF (tl_obs(itime)%atmswv(i,j,ie)%q%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%atmswv(i,j,ie)%q%val**2.D0
             END IF
       END DO
       END DO
       END DO
    END IF
    ! mhs observation
    IF (da_mhs) THEN
       DO j  = 1, obsmet(itime)%nmhsch
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%nmhs_e(ie)
             IF (tl_obs(itime)%mhs(i,j,ie)%t%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%mhs(i,j,ie)%t%val**2.D0
             END IF
             IF (tl_obs(itime)%mhs(i,j,ie)%q%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%mhs(i,j,ie)%q%val**2.D0
             END IF
       END DO
       END DO
       END DO
    END IF
    ! csr observation
    IF (da_csr) THEN
       DO j  = 1, obsmet(itime)%ncsrch
       DO ie = nets, nete
       DO i  = 1, obsmet(itime)%ncsr_e(ie)
             IF (tl_obs(itime)%csr(i,j,ie)%t%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%csr(i,j,ie)%t%val**2.D0
             END IF
             IF (tl_obs(itime)%csr(i,j,ie)%q%val  .NE. sdmiss) THEN
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%csr(i,j,ie)%q%val**2.D0
             END IF
       END DO
       END DO
       END DO
    END IF

  END DO ! itime

    CALL Wrap_Repro_Sum(nvars=1, comm=hybrid%par%comm)
    IF (par%ismasterproc) WRITE(iulog,*) ' '
    IF (par%ismasterproc) WRITE(iulog,*) 'Adjoint test '
    IF (par%ismasterproc) THEN
       yTy = global_shared_sum(1)
       WRITE(iulog,*) '<Hx,Hx>     = ', yTy
    ENDIF

    CALL InnerProduct(neig,        &
                      neig_chi,    &
                      neig_Tps,    &
                      neig_q,      &
                      zwn_nproc,   &
                      zwn_nproc_a, &
                      b_eig,       &
                      b_eig_chi,   &
                      b_eig_Tps,   &
                      b_eig_q,     &
                      b_eig_a,   &
                      ad_eig,      &
                      ad_eig_chi,  &
                      ad_eig_Tps,  &
                      ad_eig_q,    &
                      ad_eig_a,  &
                      tmp)
    IF (par%ismasterproc) THEN
    WRITE(iulog,*) '<H^T Hx, x> = ', tmp 
    WRITE(iulog,*) '<Hx,Hx> - <H^T Hx, x> = ', abs( yTy - tmp) 
    ENDIF

    RETURN

  END SUBROUTINE Check_Adj
!===============================================================================
!===============================================================================
!  SUBROUTINE InnerProduct
!
!> @brief
!> - Calculate Inner Product
!
!> @date 15Sep2014
!> - S. KIM: Initial code
!> @date 30MAY2016
!> - HJ SONG: Modified for alpha-variable
!
!> @param[in]  vect1, vect2
!> @param[out] ddot_prodcut
!===============================================================================
  SUBROUTINE InnerProduct(neig,        &
                          neig_chi,    &
                          neig_Tps,    &
                          neig_q,      &
                          zwn_nproc,   &
                          zwn_nproc_a, &
                          vect1,       &
                          vect1_chi,   &
                          vect1_Tps,   &
                          vect1_q,     &
                          vect1_a,     &
                          vect2,       &
                          vect2_chi,   &
                          vect2_Tps,   &
                          vect2_q,     &
                          vect2_a,     &
                          ddot_product)
  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules

  ! 2. Variables
  ! 2-1. Input
    INTEGER(KIND=int_kind),          &
       INTENT(IN) :: neig,           &
                     neig_chi,       &
                     neig_Tps,       &
                     neig_q,         &
                     zwn_nproc,      &
                     zwn_nproc_a
    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc), &
       INTENT(IN) :: vect1,vect2
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc), &
       INTENT(IN) :: vect1_chi,vect2_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc), &
       INTENT(IN) :: vect1_Tps,vect2_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc), &
       INTENT(IN) :: vect1_q,vect2_q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl), &
       INTENT(IN) :: vect1_a,vect2_a  !--- for alpha (HJS)

  ! 2-2. Output
    REAL(KIND=real_kind), INTENT(OUT) :: ddot_product 

  ! 2-3. Local
    REAL(KIND=real_kind) :: tmp,tmp1
    INTEGER(KIND=int_kind)     &
       :: i, j, k, ie, wn,jump, &
          ismpl ! 4DEnVAR
  !=============================================================================

  !=============================================================================
  ! B. Main body
  !=============================================================================
    tmp1= 0.D0
    ddot_product = 0.D0
    jump=5
    k=mod(zwn_nproc,jump)
    j=k+1
    IF (k .NE. 0) Then
       DO wn = 1, k
          tmp1 =tmp1 + dot_product(vect1(:,wn),vect2(:,wn)) +         &
                       dot_product(vect1_chi(:,wn),vect2_chi(:,wn)) + &
                       dot_product(vect1_Tps(:,wn),vect2_Tps(:,wn)) + &
                       dot_product(vect1_q(:,wn),vect2_q(:,wn)) 
       END DO
       !IF (par%ismasterproc) WRITE(iulog,*) 'tmp 1= ', tmp1,k,zwn_nproc,jump
    END IF
 
    DO wn = j,zwn_nproc,jump
       tmp1=tmp1                                              &
            +dot_product(vect1(:,wn),vect2(:,wn))             &
            +dot_product(vect1(:,wn+1),vect2(:,wn+1))         &
            +dot_product(vect1(:,wn+2),vect2(:,wn+2))         &
            +dot_product(vect1(:,wn+3),vect2(:,wn+3))         &
            +dot_product(vect1(:,wn+4),vect2(:,wn+4))         &
            +dot_product(vect1_chi(:,wn),vect2_chi(:,wn))     &
            +dot_product(vect1_chi(:,wn+1),vect2_chi(:,wn+1)) &
            +dot_product(vect1_chi(:,wn+2),vect2_chi(:,wn+2)) &
            +dot_product(vect1_chi(:,wn+3),vect2_chi(:,wn+3)) &
            +dot_product(vect1_chi(:,wn+4),vect2_chi(:,wn+4)) &
            +dot_product(vect1_Tps(:,wn),vect2_Tps(:,wn))     &
            +dot_product(vect1_Tps(:,wn+1),vect2_Tps(:,wn+1)) &
            +dot_product(vect1_Tps(:,wn+2),vect2_Tps(:,wn+2)) &
            +dot_product(vect1_Tps(:,wn+3),vect2_Tps(:,wn+3)) &
            +dot_product(vect1_Tps(:,wn+4),vect2_Tps(:,wn+4)) &
            +dot_product(vect1_q(:,wn),vect2_q(:,wn))         &
            +dot_product(vect1_q(:,wn+1),vect2_q(:,wn+1))     &
            +dot_product(vect1_q(:,wn+2),vect2_q(:,wn+2))     &
            +dot_product(vect1_q(:,wn+3),vect2_q(:,wn+3))     &
            +dot_product(vect1_q(:,wn+4),vect2_q(:,wn+4))
    END DO            

    k=mod(zwn_nproc_a,jump)
    j=k+1
  DO ismpl = 1, nsmpl 
    IF (k .NE. 0) Then
       DO wn = 1, k
          tmp1 =tmp1 + dot_product(vect1_a(:,wn,ismpl),vect2_a(:,wn,ismpl))
       END DO
       !IF (par%ismasterproc) WRITE(iulog,*) 'tmp 1= ', tmp1,k,zwn_nproc,jump
    END IF

    DO wn = j,zwn_nproc_a,jump
       tmp1=tmp1                                                          &
            +dot_product(vect1_a(:,wn,ismpl),vect2_a(:,wn,ismpl))         &
            +dot_product(vect1_a(:,wn+1,ismpl),vect2_a(:,wn+1,ismpl))     &
            +dot_product(vect1_a(:,wn+2,ismpl),vect2_a(:,wn+2,ismpl))     &
            +dot_product(vect1_a(:,wn+3,ismpl),vect2_a(:,wn+3,ismpl))     &
            +dot_product(vect1_a(:,wn+4,ismpl),vect2_a(:,wn+4,ismpl))
    END DO            
  END DO ! ismpl

    global_shared_buf(:,1) = 0.D0
    global_shared_buf(nets:nete,1) = tmp1/DBLE(nete-nets+1)
    CALL Wrap_Repro_Sum(nvars=1, comm=hybrid%par%comm)
    ddot_product = global_shared_sum(1)
    !IF (par%ismasterproc) WRITE(iulog,*) 'tmp 1= ', tmp1,global_shared_sum(1)

    RETURN
  END SUBROUTINE InnerProduct
!===============================================================================
!===============================================================================
  SUBROUTINE Make_LHS_Ax(gs_md,       &
                         x_eig,       &
                         x_eig_chi,   &
                         x_eig_Tps,   &
                         x_eig_q,     &
                         x_eig_a,     &
                         Ap_eig,      &
                         Ap_eig_chi,  &
                         Ap_eig_Tps,  &
                         Ap_eig_q,    &
                         Ap_eig_a)
  ! In/Output
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime), &
       INTENT(IN) :: gs_md

    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc), &
       INTENT(IN) :: x_eig
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc), &
       INTENT(IN) :: x_eig_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc), &
       INTENT(IN) :: x_eig_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc), &
       INTENT(IN) :: x_eig_q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl), &
       INTENT(IN) :: x_eig_a

    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc), &
       INTENT(OUT) :: Ap_eig
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc), &
       INTENT(OUT) :: Ap_eig_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc), &
       INTENT(OUT) :: Ap_eig_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc), &
       INTENT(OUT) :: Ap_eig_q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl), &
       INTENT(OUT) :: Ap_eig_a
 
  ! Local
    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc) &
       ::ad_eig
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc) &
       ::ad_eig_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc) &
       ::ad_eig_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc) &
       ::ad_eig_q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl) &
       ::ad_eig_a

    REAL(KIND=real_kind) &
       :: tmp
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie, wn, &
          itime, ismpl

   ! 1. Eigen space -> observation space
     ierr = TimingStart('TlmEigToObs_Make_Ax')   !--- profiling
     CALL TlmEigToObs(gs_md,       &
                      x_eig,       &
                      x_eig_chi,   &
                      x_eig_Tps,   &
                      x_eig_q,     &
                      x_eig_a,     &
                      tl_obs(:))
     ierr = TimingStop('TlmEigToObs_Make_Ax')   !--- profiling

     DO itime = 1, ntime
     ! 2. Multiply with R^-1
     ! 2.0. single observation test
       IF( singleobs )THEN
           CALL MultiplyRSingle(itime,tl_obs(itime)%single%u%val, tl_obs(itime)%single%v%val, &
                                      tl_obs(itime)%single%t%val, tl_obs(itime)%single%q%val, &
                                      tl_obs(itime)%single%ps%val                    )
       END IF
     ! 2.1. sonde observation
       IF (da_sonde) THEN
           CALL MultiplyRSonde(itime,tl_obs(itime)%sonde%u%val, tl_obs(itime)%sonde%v%val, &
                               tl_obs(itime)%sonde%t%val, tl_obs(itime)%sonde%q%val)
       END IF
     ! 2.2. surface observation
       IF (da_surface) THEN
           CALL MultiplyRSurface(itime,tl_obs(itime)%surface%u%val, tl_obs(itime)%surface%v%val, &
                                 tl_obs(itime)%surface%t%val, tl_obs(itime)%surface%q%val, tl_obs(itime)%surface%p%val)
       END IF
     ! 2.3. bogus observation
       IF (da_bogus) THEN
           CALL MultiplyRBogus(itime,tl_obs(itime)%bogus%p%val)
       END IF
     ! 2.3. aircraft observation
       IF (da_aircraft) THEN
           CALL MultiplyRAircraft(itime,tl_obs(itime)%aircraft%u%val, tl_obs(itime)%aircraft%v%val, &
                                  tl_obs(itime)%aircraft%t%val)
       END IF
     ! 2.4. amv observation
       IF (da_amv) THEN
           CALL MultiplyRAmv(itime,tl_obs(itime)%amv%u%val, tl_obs(itime)%amv%v%val)
       END IF
     ! 2.5. scatwind observation
       IF (da_scatwind) THEN
           CALL MultiplyRScatwind(itime,tl_obs(itime)%scatwind%u%val, tl_obs(itime)%scatwind%v%val)
       END IF
     ! 2.6. amsua observation
       IF (da_amsua) THEN
           CALL MultiplyRAmsua(itime,tl_obs(itime)%amsua%t%val)
       END IF
     ! 2.7. gpsro observation
       IF (da_gpsro) THEN
           CALL MultiplyRGpsro(itime,tl_obs(itime)%gpsro%t%val, tl_obs(itime)%gpsro%q%val)
       END IF
     ! 2.8. iasi observation
       IF (da_iasi) THEN
           CALL MultiplyRIasi(itime,tl_obs(itime)%iasi%t%val, tl_obs(itime)%iasi%q%val)
       END IF
     ! 2.9. cris observation
       IF (da_cris) THEN
           CALL MultiplyRCris(itime,tl_obs(itime)%cris%t%val,tl_obs(itime)%cris%q%val)
       END IF
     ! 2.10. atms observation
       IF (da_atms) THEN
           CALL MultiplyRAtms(itime,tl_obs(itime)%atms%t%val)
       END IF
     ! 2.11. atmswv observation
       IF (da_atmswv) THEN
           CALL MultiplyRAtmswv(itime,tl_obs(itime)%atmswv%t%val,tl_obs(itime)%atmswv%q%val)
       END IF
     ! 2.12. mhs observation
       IF (da_mhs) THEN
           CALL MultiplyRMhs(itime,tl_obs(itime)%mhs%t%val, tl_obs(itime)%mhs%q%val)
       END IF
     ! 2.13. csr observation
       IF (da_csr) THEN
           CALL MultiplyRCsr(itime,tl_obs(itime)%csr%t%val, tl_obs(itime)%csr%q%val)
       END IF

     END DO ! itime

  ! 3. Adjoint of spectral space -> observation space
     ierr = TimingStart('AdjEigToObs_Make_Ax')   !--- profiling
     ad_eig = 0.D0
     ad_eig_chi = 0.D0
     ad_eig_Tps = 0.D0
     ad_eig_q = 0.D0
     ad_eig_a = 0.D0
     CALL AdjEigToObs(gs_md,        &
                      tl_obs(:),    &
                      ad_eig,       &
                      ad_eig_chi,   &
                      ad_eig_Tps,   &
                      ad_eig_q,     &
                      ad_eig_a      &
                     )
     ierr = TimingStop('AdjEigToObs_Make_Ax')   !--- profiling
 
  ! 4. Finalize calculation of Ax
     DO wn = 1, zwn_nproc
        Ap_eig(:,wn) = x_eig(:,wn) + ad_eig(:,wn)
        Ap_eig_chi(:,wn) = x_eig_chi(:,wn) + ad_eig_chi(:,wn)
        Ap_eig_Tps(:,wn) = x_eig_Tps(:,wn) + ad_eig_Tps(:,wn)
        Ap_eig_q(:,wn) = x_eig_q(:,wn) + ad_eig_q(:,wn)
     END DO

     DO ismpl = 1, nsmpl
     DO wn = 1, zwn_nproc_a
        Ap_eig_a(:,wn,ismpl) = x_eig_a(:,wn,ismpl) + ad_eig_a(:,wn,ismpl)
     END DO
     END DO ! ismpl

     RETURN 

  END SUBROUTINE Make_LHS_Ax
!===============================================================================
!===============================================================================
  SUBROUTINE  Make_RHS_b(gs_md,       &
                         b_eig,       &
                         b_eig_chi,   &
                         b_eig_Tps,   &
                         b_eig_q,     &
                         b_eig_a    &
                        )
  ! 2. Variables
  ! 2-1. In/Output
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime), &
       INTENT(IN)  :: gs_md   

    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc),&
       INTENT(OUT) :: b_eig
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc),&
       INTENT(OUT) :: b_eig_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc),&
       INTENT(OUT) :: b_eig_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc),&
       INTENT(OUT) :: b_eig_q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl),&
       INTENT(OUT) :: b_eig_a

  ! 2-2. Local
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
       :: tl_md
    REAL(KIND=real_kind) &
       :: tmp
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie,     &
          wn, vr,          &
          itime

    DO itime = 1, ntime
!     IF (par%ismasterproc) WRITE(iulog,*) ' In1 ModelToObs_Make_RHS_b'
      ierr = TimingStart('ModelToObs_Make_RHS_b')  !--- profiling
         CALL ModelToObs(itime,                &
                         gs_md(:,:,:,:,itime), &
                         gs_obs(itime))
      ierr = TimingStop('ModelToObs_Make_RHS_b')
!    IF (par%ismasterproc) WRITE(iulog,*) 'Out1 ModelToObs_Make_RHS_b'
 
      tl_md = 0.D0
      DO ie = 1, nelemd
      DO  k = 1, nlevar
      DO  j = 1, np
      DO  i = 1, np
          tl_md(i,j,k,ie) = bg_md(i,j,k,ie,itime) - gs_md(i,j,k,ie,itime)
      END DO
      END DO
      END DO
      END DO

      ierr = TimingStart('TlmModelToObs_Make_RHS_b')  !--- profiling
         CALL TlmModelToObs(itime,                &
                            gs_md(:,:,:,:,itime), &
                            tl_md(:,:,:,:),       &
                            tl_obs(itime))
      ierr = TimingStop('TlmModelToObs_Make_RHS_b')

       ! single observation test
       IF( singleobs )THEN
           CALL MakebSingle(itime,gs_obs(itime)%single%u%val, gs_obs(itime)%single%v%val, &
                                  gs_obs(itime)%single%t%val, gs_obs(itime)%single%q%val, &
                                  gs_obs(itime)%single%ps%val,                            &
                                  tl_obs(itime)%single%u%val, tl_obs(itime)%single%v%val, &
                                  tl_obs(itime)%single%t%val, tl_obs(itime)%single%q%val, &
                                  tl_obs(itime)%single%ps%val                    )
       END IF
       ! sonde observation
       IF (da_sonde) THEN
           CALL MakebSonde(itime,gs_obs(itime)%sonde%u%val, gs_obs(itime)%sonde%v%val, &
                           gs_obs(itime)%sonde%t%val, gs_obs(itime)%sonde%q%val, &
                           tl_obs(itime)%sonde%u%val, tl_obs(itime)%sonde%v%val, &
                           tl_obs(itime)%sonde%t%val, tl_obs(itime)%sonde%q%val)
       END IF
       ! surface observation
       IF (da_surface) THEN
           CALL MakebSurface(itime,gs_obs(itime)%surface%u%val, gs_obs(itime)%surface%v%val, &
                             gs_obs(itime)%surface%t%val, gs_obs(itime)%surface%q%val, gs_obs(itime)%surface%p%val, &
                             tl_obs(itime)%surface%u%val, tl_obs(itime)%surface%v%val, &
                             tl_obs(itime)%surface%t%val, tl_obs(itime)%surface%q%val, tl_obs(itime)%surface%p%val)
       END IF
       ! bogus observation
       IF (da_bogus) THEN
           CALL MakebBogus(itime,gs_obs(itime)%bogus%p%val, tl_obs(itime)%bogus%p%val)
       END IF
       ! aircraft observation
       IF (da_aircraft) THEN
           CALL MakebAircraft(itime,gs_obs(itime)%aircraft%u%val, gs_obs(itime)%aircraft%v%val, &
                              gs_obs(itime)%aircraft%t%val,                                    &
                              tl_obs(itime)%aircraft%u%val, &
                              tl_obs(itime)%aircraft%v%val, tl_obs(itime)%aircraft%t%val)
       END IF
       ! amv observation
       IF (da_amv) THEN
           CALL MakebAmv(itime,gs_obs(itime)%amv%u%val, gs_obs(itime)%amv%v%val, &
                         tl_obs(itime)%amv%u%val, tl_obs(itime)%amv%v%val)
       END IF
       ! scatwind observation
       IF (da_scatwind) THEN
           CALL MakebScatwind(itime,gs_obs(itime)%scatwind%u%val, gs_obs(itime)%scatwind%v%val, &
                                    tl_obs(itime)%scatwind%u%val, tl_obs(itime)%scatwind%v%val)
       END IF
       ! amsua observation
       IF (da_amsua) THEN
           CALL MakebAmsua(itime,gs_obs(itime)%amsua%t%val, &
                           tl_obs(itime)%amsua%t%val)
       END IF
       ! gpsro observation
       IF (da_gpsro) THEN
           CALL MakebGpsro(itime,gs_obs(itime)%gpsro%t%val, gs_obs(itime)%gpsro%q%val, &
                           tl_obs(itime)%gpsro%t%val, tl_obs(itime)%gpsro%q%val)
       END IF
       ! iasi observation
       IF (da_iasi) THEN
           CALL MakebIasi(itime,gs_obs(itime)%iasi%t%val, gs_obs(itime)%iasi%q%val, &
                                tl_obs(itime)%iasi%t%val, tl_obs(itime)%iasi%q%val)
       END IF
       ! cris observation
       IF (da_cris) THEN
           CALL MakebCris(itime,gs_obs(itime)%cris%t%val, gs_obs(itime)%cris%q%val, &
                                tl_obs(itime)%cris%t%val, tl_obs(itime)%cris%q%val)
       END IF
       ! atms observation
       IF (da_atms) THEN
           CALL MakebAtms(itime,gs_obs(itime)%atms%t%val, &
                          tl_obs(itime)%atms%t%val)
       END IF
       ! atmswv observation
       IF (da_atmswv) THEN
           CALL MakebAtmswv(itime,gs_obs(itime)%atmswv%t%val, gs_obs(itime)%atmswv%q%val,&
                                  tl_obs(itime)%atmswv%t%val, tl_obs(itime)%atmswv%q%val)
       END IF
       ! mhs observation
       IF (da_mhs) THEN
           CALL MakebMhs(itime,gs_obs(itime)%mhs%t%val, gs_obs(itime)%mhs%q%val, &
                               tl_obs(itime)%mhs%t%val, tl_obs(itime)%mhs%q%val)
       END IF
       ! csr observation
       IF (da_csr) THEN
           CALL MakebCsr(itime,gs_obs(itime)%csr%t%val, gs_obs(itime)%csr%q%val, &
                               tl_obs(itime)%csr%t%val, tl_obs(itime)%csr%q%val)
       END IF

     END DO ! itime

     ierr = TimingStart('AdjEigToObs_Make_RHS_b')   !--- profiling
     b_eig = 0.D0
     b_eig_chi = 0.D0
     b_eig_Tps = 0.D0
     b_eig_q = 0.D0
     b_eig_a = 0.D0
     CALL AdjEigToObs(gs_md(:,:,:,:,:), &
                      tl_obs(:),   &
                      b_eig,       &
                      b_eig_chi,   &
                      b_eig_Tps,   &
                      b_eig_q,     &
                      b_eig_a    &
                      )
     ierr = TimingStop('AdjEigToObs_Make_RHS_b')   !--- profiling

     b_eig = b_eig - x_eig_acc
     b_eig_chi = b_eig_chi - x_eig_acc_chi
     b_eig_Tps = b_eig_Tps - x_eig_acc_Tps
     b_eig_q = b_eig_q - x_eig_acc_q
     b_eig_a = b_eig_a - x_eig_acc_a

     RETURN

  END SUBROUTINE Make_RHS_b
!===============================================================================
!===============================================================================
  SUBROUTINE Cal_Cost_Func(gs_md,       &
                           x_eig,       &
                           x_eig_chi,   &
                           x_eig_Tps,   &
                           x_eig_q,     &
                           x_eig_a,     &
                           j_o,         &
                           j_b)
  ! 2. Variables
  ! 2-1. In/Output
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime), &
       INTENT(IN) :: gs_md

    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc), &
       INTENT(IN) :: x_eig   
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc), &
       INTENT(IN) :: x_eig_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc), &
       INTENT(IN) :: x_eig_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc), &
       INTENT(IN) :: x_eig_q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl), &
       INTENT(IN) :: x_eig_a

    REAL(KIND=real_kind),      &
       INTENT(OUT) :: j_b, j_o

  ! 2-2. Local
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
       :: tl_md
    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc) &
       :: tl_eig, tlx_eig
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc) &
       :: tl_eig_chi, tlx_eig_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc) &
       :: tl_eig_Tps, tlx_eig_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc) &
       :: tl_eig_q, tlx_eig_q
    REAL(KIND=real_kind)  &
       :: tmpmax, tmpmin, &
          maxv, minv
    REAL(KIND=real_kind)            &
       :: rmse_an_u, rmse_an_v,     &
          rmse_an_t, rmse_an_q,     &
          rmse_an_ps,tmp,           &
          jo_single,                &
          jo_sonde, jo_surface,     &
          jo_bogus,                 &
          jo_aircraft, jo_amv,      &
          jo_scatwind,              &
          jo_amsua, jo_gpsro,       &
          jo_iasi, jo_cris,         &
          jo_atms, jo_atmswv,       &
          jo_mhs, jo_csr
    INTEGER(KIND=int_kind)      &
       :: i, j, k, ie, wn,ierr, &
          itime, ismpl

    ! 1-1. Jo
    DO itime = 1, ntime
       ierr = TimingStart('ModelToObs_Cal_CostFun')   !--- profiling
          CALL ModelToObs(itime,gs_md(:,:,:,:,itime), gs_obs(itime))
       ierr = TimingStop('ModelToObs_Cal_CostFun')   !--- profiling
    END DO ! itime

    ierr = TimingStart('TlmEigToObs_Cal_CostFun')   !--- profiling
       CALL TlmEigToObs(gs_md(:,:,:,:,:), &
                        x_eig,            &
                        x_eig_chi,        &
                        x_eig_Tps,        &
                        x_eig_q,          &
                        x_eig_a,          &
                        tl_obs(:))
    ierr = TimingStop('TlmEigToObs_Cal_CostFun')   !--- profiling

    global_shared_buf(:,2:46) = 0.D0
    DO itime = 1, ntime
       ! Single observation test
       IF( singleobs )THEN
          DO ie = nets, nete
          DO i  = 1, obsmet(itime)%nsg_e(ie)
             IF ((obs(itime)%single(i,ie)%u%val .NE. sdmiss) .AND.   (gs_obs(itime)%single(i,ie)%u%val .NE. sdmiss) .AND. (tl_obs(itime)%single(i,ie)%u%val .NE. sdmiss)) THEN
                  global_shared_buf(ie,2) = global_shared_buf(ie,2)                 &
                                           +( gs_obs(itime)%single(i,ie)%u%val      &
                                             +tl_obs(itime)%single(i,ie)%u%val      &
                                             -obs(itime)%single(i,ie)%u%val )**2.D0 &
                                           /( obs(itime)%single(i,ie)%u%error )
             END IF
             IF ((obs(itime)%single(i,ie)%v%val .NE. sdmiss) .AND.   (gs_obs(itime)%single(i,ie)%v%val .NE. sdmiss) .AND. (tl_obs(itime)%single(i,ie)%v%val .NE. sdmiss)) THEN
                  global_shared_buf(ie,3) = global_shared_buf(ie,3)                 &
                                           +( gs_obs(itime)%single(i,ie)%v%val      &
                                             +tl_obs(itime)%single(i,ie)%v%val      &
                                             -obs(itime)%single(i,ie)%v%val )**2.D0 &
                                           /( obs(itime)%single(i,ie)%v%error )
             END IF
             IF ((obs(itime)%single(i,ie)%t%val .NE. sdmiss) .AND.    (gs_obs(itime)%single(i,ie)%t%val .NE. sdmiss) .AND. (tl_obs(itime)%single(i,ie)%t%val .NE. sdmiss)) THEN
                  global_shared_buf(ie,4) = global_shared_buf(ie,4)                 &
                                           +( gs_obs(itime)%single(i,ie)%t%val      &
                                             +tl_obs(itime)%single(i,ie)%t%val      &
                                             -obs(itime)%single(i,ie)%t%val )**2.D0 &
                                           /( obs(itime)%single(i,ie)%t%error )
             END IF
             IF ((obs(itime)%single(i,ie)%q%val .NE. sdmiss) .AND.    (gs_obs(itime)%single(i,ie)%q%val .NE. sdmiss) .AND. (tl_obs(itime)%single(i,ie)%q%val .NE. sdmiss)) THEN
                  global_shared_buf(ie,5) = global_shared_buf(ie,5)                 &
                                           +( gs_obs(itime)%single(i,ie)%q%val      &
                                             +tl_obs(itime)%single(i,ie)%q%val      &
                                             -obs(itime)%single(i,ie)%q%val )**2.D0 &
                                           /( obs(itime)%single(i,ie)%q%error )
             END IF
             IF ((obs(itime)%single(i,ie)%ps%val .NE. sdmiss) .AND.    (gs_obs(itime)%single(i,ie)%ps%val .NE. sdmiss) .AND. (tl_obs(itime)%single(i,ie)%ps%val .NE. sdmiss)) THEN
                  global_shared_buf(ie,6) = global_shared_buf(ie,6)                  &
                                           +( gs_obs(itime)%single(i,ie)%ps%val      &
                                             +tl_obs(itime)%single(i,ie)%ps%val      &
                                             -obs(itime)%single(i,ie)%ps%val )**2.D0 &
                                           /( obs(itime)%single(i,ie)%ps%error )
             END IF
          END DO
          END DO
       END IF

      ! sonde observation
      IF (da_sonde) THEN
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%nsd_e(1,ie)
            IF ((obs(itime)%sonde(i,ie)%u%val .NE. sdmiss) .AND. (gs_obs(itime)%sonde(i,ie)%u%val .NE. sdmiss) .AND. (tl_obs(itime)%sonde(i,ie)%u%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,2) = global_shared_buf(ie,2) +      &
                                           (gs_obs(itime)%sonde(i,ie)%u%val +    &
                                            tl_obs(itime)%sonde(i,ie)%u%val -    &
                                            obs(itime)%sonde(i,ie)%u%val)**2.D0/ &
                                           (obs(itime)%sonde(i,ie)%u%error)
                global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
             END IF
         END DO
         DO i  = 1, obsmet(itime)%nsd_e(2,ie)
             IF ((obs(itime)%sonde(i,ie)%v%val .NE. sdmiss) .AND.   (gs_obs(itime)%sonde(i,ie)%v%val .NE. sdmiss) .AND. (tl_obs(itime)%sonde(i,ie)%v%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,3) = global_shared_buf(ie,3) +      &
                                           (gs_obs(itime)%sonde(i,ie)%v%val +    &
                                            tl_obs(itime)%sonde(i,ie)%v%val -    &
                                            obs(itime)%sonde(i,ie)%v%val)**2.D0/ &
                                           (obs(itime)%sonde(i,ie)%v%error)
                global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
             END IF
          END DO
         DO i  = 1, obsmet(itime)%nsd_e(3,ie)
             IF ((obs(itime)%sonde(i,ie)%t%val .NE. sdmiss) .AND.   (gs_obs(itime)%sonde(i,ie)%t%val .NE. sdmiss) .AND. (tl_obs(itime)%sonde(i,ie)%t%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,4) = global_shared_buf(ie,4) +      &
                                           (gs_obs(itime)%sonde(i,ie)%t%val +    &
                                            tl_obs(itime)%sonde(i,ie)%t%val -    &
                                            obs(itime)%sonde(i,ie)%t%val)**2.D0/ &
                                           (obs(itime)%sonde(i,ie)%t%error)
                global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
             END IF
         END DO
         DO i  = 1, obsmet(itime)%nsd_e(4,ie)
             IF ((obs(itime)%sonde(i,ie)%q%val .NE. sdmiss) .AND.   (gs_obs(itime)%sonde(i,ie)%q%val .NE. sdmiss) .AND. (tl_obs(itime)%sonde(i,ie)%q%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,5) = global_shared_buf(ie,5) +      &
                                           (gs_obs(itime)%sonde(i,ie)%q%val +    &
                                            tl_obs(itime)%sonde(i,ie)%q%val -    &
                                            obs(itime)%sonde(i,ie)%q%val)**2.D0/ &
                                           (obs(itime)%sonde(i,ie)%q%error)
                global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
             END IF
         END DO
         END DO
      END IF

      ! surface observation
      IF (da_surface) THEN
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%nsf_e(ie)
            IF ((obs(itime)%surface(i,ie)%u%val    .NE. sdmiss) .AND. (gs_obs(itime)%surface(i,ie)%u%val .NE. sdmiss) .AND. (tl_obs(itime)%surface(i,ie)%u%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,7) = global_shared_buf(ie,7) +        &
                                           (gs_obs(itime)%surface(i,ie)%u%val +    &
                                            tl_obs(itime)%surface(i,ie)%u%val -    &
                                            obs(itime)%surface(i,ie)%u%val)**2.D0/ &
                                           (obs(itime)%surface(i,ie)%u%error)
                global_shared_buf(ie,12) = global_shared_buf(ie,12) + 1.D0
             END IF
             IF ((obs(itime)%surface(i,ie)%v%val   .NE. sdmiss) .AND. (gs_obs(itime)%surface(i,ie)%v%val .NE. sdmiss) .AND. (tl_obs(itime)%surface(i,ie)%v%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,8) = global_shared_buf(ie,8) +        &
                                           (gs_obs(itime)%surface(i,ie)%v%val +    &
                                            tl_obs(itime)%surface(i,ie)%v%val -    &
                                            obs(itime)%surface(i,ie)%v%val)**2.D0/ &
                                           (obs(itime)%surface(i,ie)%v%error)
                global_shared_buf(ie,12) = global_shared_buf(ie,12) + 1.D0
             END IF
             IF ((obs(itime)%surface(i,ie)%t%val   .NE. sdmiss) .AND. (gs_obs(itime)%surface(i,ie)%t%val .NE. sdmiss) .AND. (tl_obs(itime)%surface(i,ie)%t%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,9) = global_shared_buf(ie,9) +        &
                                           (gs_obs(itime)%surface(i,ie)%t%val +    &
                                            tl_obs(itime)%surface(i,ie)%t%val -    &
                                            obs(itime)%surface(i,ie)%t%val)**2.D0/ &
                                           (obs(itime)%surface(i,ie)%t%error)
                global_shared_buf(ie,12) = global_shared_buf(ie,12) + 1.D0
             END IF
             IF ((obs(itime)%surface(i,ie)%q%val   .NE. sdmiss) .AND. (gs_obs(itime)%surface(i,ie)%q%val .NE. sdmiss) .AND. (tl_obs(itime)%surface(i,ie)%q%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,10) = global_shared_buf(ie,10) +        &
                                           (gs_obs(itime)%surface(i,ie)%q%val +    &
                                            tl_obs(itime)%surface(i,ie)%q%val -    &
                                            obs(itime)%surface(i,ie)%q%val)**2.D0/ &
                                           (obs(itime)%surface(i,ie)%q%error)
                global_shared_buf(ie,12) = global_shared_buf(ie,12) + 1.D0
             END IF
             IF ((obs(itime)%surface(i,ie)%p%val .NE. sdmiss) .AND.   (gs_obs(itime)%surface(i,ie)%p%val .NE. sdmiss) .AND. (tl_obs(itime)%surface(i,ie)%p%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,11) = global_shared_buf(ie,11) +      &
                                           (gs_obs(itime)%surface(i,ie)%p%val +    &
                                            tl_obs(itime)%surface(i,ie)%p%val -    &
                                            obs(itime)%surface(i,ie)%p%val)**2.D0/ &
                                           (obs(itime)%surface(i,ie)%p%error)
                global_shared_buf(ie,12) = global_shared_buf(ie,12) + 1.D0
             END IF
         END DO
         END DO
      END IF

      ! bogus observation
      IF (da_bogus) THEN
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%nbg_e(ie)
             IF ((obs(itime)%bogus(i,ie)%p%val .NE. sdmiss) .AND.   (gs_obs(itime)%bogus(i,ie)%p%val .NE. sdmiss) .AND. (tl_obs(itime)%bogus(i,ie)%p%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,13) = global_shared_buf(ie,13) +      &
                                           (gs_obs(itime)%bogus(i,ie)%p%val +    &
                                            tl_obs(itime)%bogus(i,ie)%p%val -    &
                                            obs(itime)%bogus(i,ie)%p%val)**2.D0/ &
                                           (obs(itime)%bogus(i,ie)%p%error)
                global_shared_buf(ie,14) = global_shared_buf(ie,14) + 1.D0
             END IF
         END DO
         END DO
      END IF


      ! aircraft observation
      IF (da_aircraft) THEN
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%nar_e(1,ie)
            IF ((obs(itime)%aircraft(i,ie)%u%val    .NE. sdmiss) .AND. (gs_obs(itime)%aircraft(i,ie)%u%val .NE. sdmiss) .AND. (tl_obs(itime)%aircraft(i,ie)%u%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,15) = global_shared_buf(ie,15) +       &
                                           (gs_obs(itime)%aircraft(i,ie)%u%val +    &
                                            tl_obs(itime)%aircraft(i,ie)%u%val -    &
                                            obs(itime)%aircraft(i,ie)%u%val)**2.D0/ &
                                           (obs(itime)%aircraft(i,ie)%u%error)
                 global_shared_buf(ie,18) = global_shared_buf(ie,18) + 1.D0
             END IF
         END DO
         DO i  = 1, obsmet(itime)%nar_e(2,ie)
            IF ((obs(itime)%aircraft(i,ie)%v%val    .NE. sdmiss) .AND. (gs_obs(itime)%aircraft(i,ie)%v%val .NE. sdmiss) .AND. (tl_obs(itime)%aircraft(i,ie)%v%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,16) = global_shared_buf(ie,16) +       &
                                           (gs_obs(itime)%aircraft(i,ie)%v%val +    &
                                            tl_obs(itime)%aircraft(i,ie)%v%val -    &
                                            obs(itime)%aircraft(i,ie)%v%val)**2.D0/ &
                                           (obs(itime)%aircraft(i,ie)%v%error)
                 global_shared_buf(ie,18) = global_shared_buf(ie,18) + 1.D0
             END IF
         END DO
         DO i  = 1, obsmet(itime)%nar_e(3,ie)
            IF ((obs(itime)%aircraft(i,ie)%t%val    .NE. sdmiss) .AND. (gs_obs(itime)%aircraft(i,ie)%t%val .NE. sdmiss) .AND. (tl_obs(itime)%aircraft(i,ie)%t%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,17) = global_shared_buf(ie,17) +       &
                                           (gs_obs(itime)%aircraft(i,ie)%t%val +    &
                                            tl_obs(itime)%aircraft(i,ie)%t%val -    &
                                            obs(itime)%aircraft(i,ie)%t%val)**2.D0/ &
                                           (obs(itime)%aircraft(i,ie)%t%error)
                 global_shared_buf(ie,18) = global_shared_buf(ie,18) + 1.D0
             END IF
         END DO
         END DO
      END IF

         ! amv observation
         IF (da_amv) THEN
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%namv_e(ie)
            IF ((obs(itime)%amv(i,ie)%u%val    .NE. sdmiss) .AND. (gs_obs(itime)%amv(i,ie)%u%val .NE. sdmiss) .AND. (tl_obs(itime)%amv(i,ie)%u%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,19) = global_shared_buf(ie,19) +       &
                                           (gs_obs(itime)%amv(i,ie)%u%val +    &
                                            tl_obs(itime)%amv(i,ie)%u%val -    &
                                            obs(itime)%amv(i,ie)%u%val)**2.D0/ &
                                           (obs(itime)%amv(i,ie)%u%error)
                 global_shared_buf(ie,21) = global_shared_buf(ie,21) + 1.D0
             END IF
            IF ((obs(itime)%amv(i,ie)%v%val    .NE. sdmiss) .AND. (gs_obs(itime)%amv(i,ie)%v%val .NE. sdmiss) .AND. (tl_obs(itime)%amv(i,ie)%v%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,20) = global_shared_buf(ie,20) +       &
                                           (gs_obs(itime)%amv(i,ie)%v%val +    &
                                            tl_obs(itime)%amv(i,ie)%v%val -    &
                                            obs(itime)%amv(i,ie)%v%val)**2.D0/ &
                                           (obs(itime)%amv(i,ie)%v%error)
                 global_shared_buf(ie,21) = global_shared_buf(ie,21) + 1.D0
             END IF
         END DO
         END DO
         END IF

      ! scatwind observation
      IF (da_scatwind) THEN
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%nscatwind_e(ie)
            IF ((obs(itime)%scatwind(i,ie)%u%val    .NE. sdmiss) .AND. (gs_obs(itime)%scatwind(i,ie)%u%val .NE. sdmiss) .AND. (tl_obs(itime)%scatwind(i,ie)%u%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,22) = global_shared_buf(ie,22) +        &
                                           (gs_obs(itime)%scatwind(i,ie)%u%val +    &
                                            tl_obs(itime)%scatwind(i,ie)%u%val -    &
                                            obs(itime)%scatwind(i,ie)%u%val)**2.D0/ &
                                           (obs(itime)%scatwind(i,ie)%u%error)
                global_shared_buf(ie,24) = global_shared_buf(ie,24) + 1.D0
             END IF
             IF ((obs(itime)%scatwind(i,ie)%v%val   .NE. sdmiss) .AND. (gs_obs(itime)%scatwind(i,ie)%v%val .NE. sdmiss) .AND. (tl_obs(itime)%scatwind(i,ie)%v%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,23) = global_shared_buf(ie,23) +        &
                                           (gs_obs(itime)%scatwind(i,ie)%v%val +    &
                                            tl_obs(itime)%scatwind(i,ie)%v%val -    &
                                            obs(itime)%scatwind(i,ie)%v%val)**2.D0/ &
                                           (obs(itime)%scatwind(i,ie)%v%error)
                global_shared_buf(ie,24) = global_shared_buf(ie,24) + 1.D0
             END IF
         END DO
         END DO
      END IF

      ! amsua observation
      IF (da_amsua) THEN
         DO j  = 1, obsmet(itime)%namsuach
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%namsua_e(ie)
            IF ((obs(itime)%amsua(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%amsua(i,j,ie)%t%val .NE. sdmiss) .AND. (tl_obs(itime)%amsua(i,j,ie)%t%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,25) = global_shared_buf(ie,25) +      &
                                           (gs_obs(itime)%amsua(i,j,ie)%t%val +    &
                                            tl_obs(itime)%amsua(i,j,ie)%t%val -    &
                                            obs(itime)%amsua(i,j,ie)%tb%val)**2.D0/ &
                                           (obs(itime)%amsua(i,j,ie)%tb%error)
                 global_shared_buf(ie,26) = global_shared_buf(ie,26) + 1.D0
            END IF
         END DO
         END DO
         END DO
      END IF

      ! gpsro observation
      IF (da_gpsro) THEN
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%ngpsro_e(ie)
            IF ((obs(itime)%gpsro(i,ie)%ba%val    .NE. sdmiss) .AND. (gs_obs(itime)%gpsro(i,ie)%t%val .NE. sdmiss) .AND.  (tl_obs(itime)%gpsro(i,ie)%t%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,27) = global_shared_buf(ie,27) + &
                                           (gs_obs(itime)%gpsro(i,ie)%t%val +        &
                                            tl_obs(itime)%gpsro(i,ie)%t%val -        &
                                            obs(itime)%gpsro(i,ie)%ba%val)**2.D0/    &
                                           (obs(itime)%gpsro(i,ie)%ba%error)
                 global_shared_buf(ie,29) = global_shared_buf(ie,29) + 1.D0
            END IF
            IF ((obs(itime)%gpsro(i,ie)%ba%val    .NE. sdmiss) .AND. (gs_obs(itime)%gpsro(i,ie)%q%val .NE. sdmiss) .AND.  (tl_obs(itime)%gpsro(i,ie)%q%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,28) = global_shared_buf(ie,28) + &
                                           (gs_obs(itime)%gpsro(i,ie)%q%val +        &
                                            tl_obs(itime)%gpsro(i,ie)%q%val -        &
                                            obs(itime)%gpsro(i,ie)%ba%val)**2.D0/    &
                                           (obs(itime)%gpsro(i,ie)%ba%error)
                 global_shared_buf(ie,29) = global_shared_buf(ie,29) + 1.D0
            END IF
         END DO
         END DO
      END IF

      ! iasi observation
      IF (da_iasi) THEN
         DO j = 1, obsmet(itime)%niasich
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%niasi_e(ie)
            IF ((obs(itime)%iasi(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%iasi(i,j,ie)%t%val .NE. sdmiss) .AND.  (tl_obs(itime)%iasi(i,j,ie)%t%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,30) = global_shared_buf(ie,30) + &
                                           (gs_obs(itime)%iasi(i,j,ie)%t%val +        &
                                            tl_obs(itime)%iasi(i,j,ie)%t%val -        &
                                            obs(itime)%iasi(i,j,ie)%tb%val)**2.D0/    &
                                           (obs(itime)%iasi(i,j,ie)%tb%error)
                 global_shared_buf(ie,32) = global_shared_buf(ie,32) + 1.D0
            END IF
            IF ((obs(itime)%iasi(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%iasi(i,j,ie)%q%val .NE. sdmiss) .AND.  (tl_obs(itime)%iasi(i,j,ie)%q%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,31) = global_shared_buf(ie,31) + &
                                           (gs_obs(itime)%iasi(i,j,ie)%q%val +        &
                                            tl_obs(itime)%iasi(i,j,ie)%q%val -        &
                                            obs(itime)%iasi(i,j,ie)%tb%val)**2.D0/    &
                                           (obs(itime)%iasi(i,j,ie)%tb%error)
                 global_shared_buf(ie,32) = global_shared_buf(ie,32) + 1.D0
            END IF
         END DO
         END DO
         END DO
      END IF

     ! cris observation
      IF (da_cris) THEN
         DO j = 1, obsmet(itime)%ncrisch
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%ncris_e(ie)
            IF ((obs(itime)%cris(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%cris(i,j,ie)%t%val .NE. sdmiss) .AND.  (tl_obs(itime)%cris(i,j,ie)%t%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,33) = global_shared_buf(ie,33) + &
                                           (gs_obs(itime)%cris(i,j,ie)%t%val +        &
                                            tl_obs(itime)%cris(i,j,ie)%t%val -        &
                                            obs(itime)%cris(i,j,ie)%tb%val)**2.D0/    &
                                           (obs(itime)%cris(i,j,ie)%tb%error)
                 global_shared_buf(ie,35) = global_shared_buf(ie,35) + 1.D0
            END IF
            IF ((obs(itime)%cris(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%cris(i,j,ie)%q%val .NE. sdmiss) .AND.  (tl_obs(itime)%cris(i,j,ie)%q%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,34) = global_shared_buf(ie,34) + &
                                           (gs_obs(itime)%cris(i,j,ie)%q%val +        &
                                            tl_obs(itime)%cris(i,j,ie)%q%val -        &
                                            obs(itime)%cris(i,j,ie)%tb%val)**2.D0/    &
                                           (obs(itime)%cris(i,j,ie)%tb%error)
                 global_shared_buf(ie,35) = global_shared_buf(ie,35) + 1.D0
            END IF
         END DO
         END DO
         END DO
      END IF

      ! atms observation
      IF (da_atms) THEN
         DO j  = 1, obsmet(itime)%natmsch
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%natms_e(ie)
            IF ((obs(itime)%atms(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%atms(i,j,ie)%t%val .NE. sdmiss) .AND. (tl_obs(itime)%atms(i,j,ie)%t%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,36) = global_shared_buf(ie,36) +      &
                                           (gs_obs(itime)%atms(i,j,ie)%t%val +    &
                                            tl_obs(itime)%atms(i,j,ie)%t%val -    &
                                            obs(itime)%atms(i,j,ie)%tb%val)**2.D0/ &
                                           (obs(itime)%atms(i,j,ie)%tb%error)
                 global_shared_buf(ie,37) = global_shared_buf(ie,37) + 1.D0
            END IF
         END DO
         END DO
         END DO
      END IF

      ! atmswv observation
      IF (da_atmswv) THEN
         DO j  = 1, obsmet(itime)%natmswvch
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%natmswv_e(ie)
            IF ((obs(itime)%atmswv(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%atmswv(i,j,ie)%t%val .NE. sdmiss) .AND. (tl_obs(itime)%atmswv(i,j,ie)%t%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,38) = global_shared_buf(ie,38) +      &
                                           (gs_obs(itime)%atmswv(i,j,ie)%t%val +    &
                                            tl_obs(itime)%atmswv(i,j,ie)%t%val -    &
                                            obs(itime)%atmswv(i,j,ie)%tb%val)**2.D0/ &
                                           (obs(itime)%atmswv(i,j,ie)%tb%error)
                 global_shared_buf(ie,40) = global_shared_buf(ie,40) + 1.D0
            END IF
            IF ((obs(itime)%atmswv(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%atmswv(i,j,ie)%q%val .NE. sdmiss) .AND. (tl_obs(itime)%atmswv(i,j,ie)%q%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,39) = global_shared_buf(ie,39) +      &
                                           (gs_obs(itime)%atmswv(i,j,ie)%q%val +    &
                                            tl_obs(itime)%atmswv(i,j,ie)%q%val -    &
                                            obs(itime)%atmswv(i,j,ie)%tb%val)**2.D0/ &
                                           (obs(itime)%atmswv(i,j,ie)%tb%error)
                 global_shared_buf(ie,40) = global_shared_buf(ie,40) + 1.D0
            END IF

         END DO
         END DO
         END DO
      END IF

      ! mhs observation
      IF (da_mhs) THEN
         DO j  = 1, obsmet(itime)%nmhsch
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%nmhs_e(ie)
            IF ((obs(itime)%mhs(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%mhs(i,j,ie)%t%val .NE. sdmiss) .AND. (tl_obs(itime)%mhs(i,j,ie)%t%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,41) = global_shared_buf(ie,41) +      &
                                           (gs_obs(itime)%mhs(i,j,ie)%t%val +    &
                                            tl_obs(itime)%mhs(i,j,ie)%t%val -    &
                                            obs(itime)%mhs(i,j,ie)%tb%val)**2.D0/ &
                                           (obs(itime)%mhs(i,j,ie)%tb%error)
                 global_shared_buf(ie,43) = global_shared_buf(ie,43) + 1.D0
            END IF
            IF ((obs(itime)%mhs(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%mhs(i,j,ie)%q%val .NE. sdmiss) .AND. (tl_obs(itime)%mhs(i,j,ie)%q%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,42) = global_shared_buf(ie,42) +      &
                                           (gs_obs(itime)%mhs(i,j,ie)%q%val +    &
                                            tl_obs(itime)%mhs(i,j,ie)%q%val -    &
                                            obs(itime)%mhs(i,j,ie)%tb%val)**2.D0/ &
                                           (obs(itime)%mhs(i,j,ie)%tb%error)
                 global_shared_buf(ie,43) = global_shared_buf(ie,43) + 1.D0
            END IF
         END DO
         END DO
         END DO
      END IF

      ! csr observation
      IF (da_csr) THEN
         DO j  = 1, obsmet(itime)%ncsrch
         DO ie = nets, nete
         DO i  = 1, obsmet(itime)%ncsr_e(ie)
            IF ((obs(itime)%csr(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%csr(i,j,ie)%t%val .NE. sdmiss) .AND. (tl_obs(itime)%csr(i,j,ie)%t%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,44) = global_shared_buf(ie,44) +      &
                                           (gs_obs(itime)%csr(i,j,ie)%t%val +    &
                                            tl_obs(itime)%csr(i,j,ie)%t%val -    &
                                            obs(itime)%csr(i,j,ie)%tb%val)**2.D0/ &
                                           (obs(itime)%csr(i,j,ie)%tb%error)
                 global_shared_buf(ie,46) = global_shared_buf(ie,46) + 1.D0
            END IF
            IF ((obs(itime)%csr(i,j,ie)%tb%val    .NE. sdmiss) .AND. (gs_obs(itime)%csr(i,j,ie)%q%val .NE. sdmiss) .AND. (tl_obs(itime)%csr(i,j,ie)%q%val .NE. sdmiss)) THEN
                 global_shared_buf(ie,45) = global_shared_buf(ie,45) +      &
                                           (gs_obs(itime)%csr(i,j,ie)%q%val +    &
                                            tl_obs(itime)%csr(i,j,ie)%q%val -    &
                                            obs(itime)%csr(i,j,ie)%tb%val)**2.D0/ &
                                           (obs(itime)%csr(i,j,ie)%tb%error)
                 global_shared_buf(ie,46) = global_shared_buf(ie,46) + 1.D0
            END IF
         END DO
         END DO
         END DO
      END IF


    END DO ! itime

     ! 1-2. Jb
       tmp = 0.D0
       DO wn = 1, zwn_nproc
          tmp = tmp + SUM((x_eig(:,wn)-x_eig_acc(:,wn))**2.D0) +   &
                      SUM((x_eig_chi(:,wn)-x_eig_acc_chi(:,wn))**2.D0) + &
                      SUM((x_eig_Tps(:,wn)-x_eig_acc_Tps(:,wn))**2.D0) + &
                      SUM((x_eig_q(:,wn)-x_eig_acc_q(:,wn))**2.D0)
       END DO

       DO ismpl = 1, nsmpl
       DO wn = 1, zwn_nproc_a
          tmp = tmp + SUM((x_eig_a(:,wn,ismpl)-x_eig_a(:,wn,ismpl))**2.D0)
       END DO
       END DO ! ismpl

       global_shared_buf(:,1) = 0.D0
       global_shared_buf(nets:nete,1) = tmp/DBLE(nete-nets+1)

      ! 1-3. Communication and evaluation
     ierr = TimingStart('Wrap_Repro_Sum_rmse_Cal_CostFun')   !--- profiling
       CALL Wrap_Repro_Sum(nvars=46, comm=hybrid%par%comm)
     ierr = TimingStop('Wrap_Repro_Sum_rmse_Cal_CostFun')   !--- profiling

       j_b = 0.5D0*global_shared_sum(1)
       j_o = 0.5D0*(global_shared_sum(2) + global_shared_sum(3) +  &
                    global_shared_sum(4) + global_shared_sum(5) +  & ! sonde
                    global_shared_sum(7) + global_shared_sum(8) +  &
                    global_shared_sum(9) + global_shared_sum(10)+  &
                    global_shared_sum(11)+                         & ! surface
                    global_shared_sum(13)+                         & ! bogus
                    global_shared_sum(15)+ global_shared_sum(16)+  &
                    global_shared_sum(17)+                         & ! aircraft
                    global_shared_sum(19)+ global_shared_sum(20)+  & ! amv
                    global_shared_sum(22)+ global_shared_sum(23)+  & ! scatwind
                    global_shared_sum(25)+                         & ! amsua
                    global_shared_sum(27)+ global_shared_sum(28)+  & ! gpsro
                    global_shared_sum(30)+ global_shared_sum(31)+  & ! iasi
                    global_shared_sum(33)+ global_shared_sum(34)+  & ! cris
                    global_shared_sum(36)+                         & ! atms
                    global_shared_sum(38)+ global_shared_sum(39)+  & ! atmswv
                    global_shared_sum(41)+ global_shared_sum(42)+  & ! mhs
                    global_shared_sum(44)+ global_shared_sum(45)   ) ! csr

       IF( singleobs )THEN
          jo_single   = 0.5D0*(global_shared_sum(2) + global_shared_sum(3) + &
                               global_shared_sum(4) + global_shared_sum(5) + &
                               global_shared_sum(6)                       )
       END IF
       IF (da_sonde) THEN
          jo_sonde    = 0.5D0*(global_shared_sum(2) + global_shared_sum(3) + &
                               global_shared_sum(4) + global_shared_sum(5))
       END IF
       IF (da_surface) THEN
          jo_surface  = 0.5D0*(global_shared_sum(7) + global_shared_sum(8) + &
                               global_shared_sum(9) + global_shared_sum(10)+ global_shared_sum(11))
       END IF
       IF (da_bogus) THEN
          jo_bogus    = 0.5D0*(global_shared_sum(13))
       END IF
       IF (da_aircraft) THEN
          jo_aircraft = 0.5D0*(global_shared_sum(15) + global_shared_sum(16) + &
                               global_shared_sum(17))
       END IF
       IF (da_amv) THEN
          jo_amv      = 0.5D0*(global_shared_sum(19) + global_shared_sum(20))
       END IF
       IF (da_scatwind) THEN
          jo_scatwind = 0.5D0*(global_shared_sum(22) + global_shared_sum(23))
       END IF
       IF (da_amsua) THEN
          jo_amsua    = 0.5D0*(global_shared_sum(25))
       END IF
       IF (da_gpsro) THEN
          jo_gpsro    = 0.5D0*(global_shared_sum(27)+global_shared_sum(28))
       END IF
       IF (da_iasi) THEN
          jo_iasi     = 0.5D0*(global_shared_sum(30)+global_shared_sum(31))
       END IF
       IF (da_cris) THEN
          jo_cris     = 0.5D0*(global_shared_sum(33)+global_shared_sum(34))
       END IF
       IF (da_atms) THEN
          jo_atms     = 0.5D0*(global_shared_sum(36))
       END IF
       IF (da_atmswv) THEN
          jo_atmswv   = 0.5D0*(global_shared_sum(38)+global_shared_sum(39))
       END IF
       IF (da_mhs) THEN
          jo_mhs     = 0.5D0*(global_shared_sum(41)+global_shared_sum(42))
       END IF
       IF (da_csr) THEN
          jo_csr     = 0.5D0*(global_shared_sum(44)+global_shared_sum(45))
       END IF

       IF (par%ismasterproc) WRITE(iulog,*) ' '
       IF (par%ismasterproc) WRITE(iulog,*) 'Check J'
       IF (par%isMasterProc) WRITE(iulog,*) 'Jb = ', j_b
       IF (par%isMasterProc) WRITE(iulog,*) 'Jo = ', j_o
       IF (par%isMasterProc) WRITE(iulog,*) 'J  = ', j_b + j_o
       IF (par%ismasterproc) WRITE(iulog,*) ' '
       IF( singleobs )THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_single    = ', jo_single
       END IF
       IF (da_sonde) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_sonde    = ', jo_sonde
       END IF
       IF (da_surface) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_surface  = ', jo_surface
       END IF
       IF (da_bogus) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_bogus    = ', jo_bogus
       END IF
       IF (da_aircraft) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_aircraft = ', jo_aircraft
       END IF
       IF (da_amv) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_amv      = ', jo_amv
       END IF
      IF (da_scatwind) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_scatwind = ', jo_scatwind
       END IF
       IF (da_amsua) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_amsua    = ', jo_amsua
       END IF
       IF (da_gpsro) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_gpsro    = ', jo_gpsro
       END IF
       IF (da_iasi) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_iasi     = ', jo_iasi
       END IF
       IF (da_cris) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_cris     = ', jo_cris
       END IF
       IF (da_atms) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_atms     = ', jo_atms
       END IF
       IF (da_atmswv) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_atmswv   = ', jo_atmswv
       END IF
       IF (da_mhs) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_mhs      = ', jo_mhs
       END IF
       IF (da_csr) THEN
       IF (par%isMasterProc) WRITE(iulog,*) 'jo_csr      = ', jo_csr
       END IF
       IF (da_sonde) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of sonde observations    = ', DINT(global_shared_sum(6))
       END IF
       IF (da_surface) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of surface observations  = ', DINT(global_shared_sum(12))
       END IF
       IF (da_bogus) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of bogus observations    = ', DINT(global_shared_sum(14))
       END IF
       IF (da_aircraft) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of aircraft observations = ', DINT(global_shared_sum(18))
       END IF
       IF (da_amv) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of amv observations      = ', DINT(global_shared_sum(21))
       END IF
       IF (da_scatwind) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of scatwind observations = ', DINT(global_shared_sum(24))
       END IF
       IF (da_amsua) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of amsua observations    = ', DINT(global_shared_sum(26))
       END IF
       IF (da_gpsro) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of gpsro observations    = ', DINT(global_shared_sum(29))
       END IF
       IF (da_iasi) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of iasi observations     = ', DINT(global_shared_sum(32))
       END IF
       IF (da_cris) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of cris observations     = ', DINT(global_shared_sum(35))
       END IF
       IF (da_atms) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of atms observations     = ', DINT(global_shared_sum(37))
       END IF
       IF (da_atmswv) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of atmswv observations   = ', DINT(global_shared_sum(40))
       END IF
       IF (da_mhs) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of mhs observations      = ', DINT(global_shared_sum(43))
       END IF
       IF (da_csr) THEN
       IF (par%isMasterProc) WRITE(iulog,*) '# of csr observations      = ', DINT(global_shared_sum(46))
       END IF

  RETURN
  END SUBROUTINE Cal_Cost_Func
!===============================================================================

  !--------- modified for alpha (HJS)
  SUBROUTINE EigToSpec(l_zwn_nproc,   &
                       l_zwn2,        &
                       l_zwn_blength, &
                       l_zwn_bstart,  &
                       eigen,         &
                       eigen_chi,     &
                       eigen_Tps,     &
                       eigen_q,       &
                       sp)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 2. Variables
  ! 2-1. Input
    INTEGER(KIND=int_kind),       &
       INTENT(IN) :: l_zwn_nproc, &
                     l_zwn2
    INTEGER(KIND=int_kind), DIMENSION(:), ALLOCATABLE, &
       INTENT(IN) :: l_zwn_blength,                    &
                     l_zwn_bstart
    REAL(KIND=real_kind), DIMENSION(neig,l_zwn_nproc), &
       INTENT(IN) :: eigen
    REAL(KIND=real_kind), DIMENSION(neig_chi,l_zwn_nproc), &
       INTENT(IN) :: eigen_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,l_zwn_nproc), &
       INTENT(IN) :: eigen_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,l_zwn_nproc), &
       INTENT(IN) :: eigen_q

  ! 2-2. Output
    REAL(KIND=real_kind), DIMENSION(nlevar,l_zwn2), &
       INTENT(OUT) :: sp

  ! 2-3. Local
    REAL(KIND=real_kind), DIMENSION(nlevar,l_zwn_nproc) &
       :: sp_tmp
    INTEGER(KIND=int_kind) &
       :: k,               &
          wn,              &
          iproc,           &
          ierr
  !=============================================================================

  !=============================================================================
  ! B. Main body
  !=============================================================================
    ierr = TimingStart('Run BeforeBcast_fromEigToSpec')

    sp_tmp = 0.D0

    DO wn = 1, l_zwn_nproc
       DO k  = 1, neig
          sp_tmp(1:nlev,wn) =    &
             sp_tmp(1:nlev,wn) + &
             bcov_sp(1:nlev,k,wn)*     &
             DSQRT(eigv(k,wn))*        &
             eigen(k,wn)
       END DO
       DO k  = 1, neig_chi
          sp_tmp(nlev+1:2*nlev,wn) =    &
             sp_tmp(nlev+1:2*nlev,wn) + &
             bcov_sp_chi(1:nlev,k,wn)*        &
             DSQRT(eigv_chi(k,wn))*           &
             eigen_chi(k,wn)
       END DO
       DO k  = 1, neig_Tps
          sp_tmp(2*nlev+1:3*nlev,wn) =    &
             sp_tmp(2*nlev+1:3*nlev,wn) + &
            bcov_sp_Tps(1:nlev,k,wn)*     &
             DSQRT(eigv_Tps(k,wn))*       &
             eigen_Tps(k,wn)
          sp_tmp(4*nlev+1,wn) =    &
             sp_tmp(4*nlev+1,wn) + &
             bcov_sp_Tps(nlev+1,k,wn)*   &
             DSQRT(eigv_Tps(k,wn))*      &
             eigen_Tps(k,wn)
       END DO
       DO k  = 1, neig_q
          sp_tmp(3*nlev+1:4*nlev,wn) =    &
             sp_tmp(3*nlev+1:4*nlev,wn) + &
             bcov_sp_q(:,k,wn)*                 &
             DSQRT(eigv_q(k,wn))*               &
             eigen_q(k,wn)
       END DO
    END DO
    ierr = TimingStop('Run BeforeBcast_fromEigToSpec')

    ierr = TimingStart('Run Bcast4SpecToCtrl_fromEigToSpec')

    CALL Mpi_AllGatherV(sp_tmp,               &
                        nlevar*l_zwn_nproc,   &
                        par_double_precision, &
                        sp,                   &
                        nlevar*l_zwn_blength, &
                        nlevar*l_zwn_bstart,  &
                        par_double_precision, &
                        par%comm,             &
                        ierr)

    ierr = TimingStop('Run Bcast4SpecToCtrl_fromEigToSpec')

    RETURN

  END SUBROUTINE EigToSpec


  ! Desinged for alpha-variable
  SUBROUTINE EigToSpecEns(l_zwn_nproc,   &
                          l_zwn2,        &
                          l_zwn_blength, &
                          l_zwn_bstart,  &
                          eigen_a,       &
                          sp)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 2. Variables
  ! 2-1. Input
    INTEGER(KIND=int_kind),       &
       INTENT(IN) :: l_zwn_nproc, &
                     l_zwn2
    INTEGER(KIND=int_kind), DIMENSION(:), ALLOCATABLE, &
       INTENT(IN) :: l_zwn_blength,                    &
                     l_zwn_bstart
    REAL(KIND=real_kind), DIMENSION(nlev+1,l_zwn_nproc), &
       INTENT(IN) :: eigen_a

  ! 2-2. Output
    REAL(KIND=real_kind), DIMENSION(nlevar,l_zwn2), &
       INTENT(OUT) :: sp

  ! 2-3. Local
    REAL(KIND=real_kind), DIMENSION(nlevar,l_zwn_nproc) &
       :: sp_tmp
    INTEGER(KIND=int_kind) &
       :: k,               &
          wn,              &
          iproc,           &
          ierr
  !=============================================================================

  !=============================================================================
  ! B. Main body
  !=============================================================================
    ierr = TimingStart('Run BeforeBcast_fromEigToSpec')

    sp_tmp = 0.D0

    DO wn = 1, l_zwn_nproc
!       DO k  = 1, nlev
!          sp_tmp(3*nlev+1:4*nlev,wn) =    &
!             sp_tmp(3*nlev+1:4*nlev,wn) + &
!             bcov_sp_a(:,k)*              &
!             DSQRT(eigv_a(k,wn))*         &
!             eigen_a(k,wn)
!       END DO
       DO k  = 1, nlev+1
          sp_tmp(2*nlev+1:3*nlev,wn) =    &
             sp_tmp(2*nlev+1:3*nlev,wn) + &
            bcov_sp_a(1:nlev,k)*          &
             DSQRT(eigv_a(k,wn))*         &
             eigen_a(k,wn)
          sp_tmp(4*nlev+1,wn) =           &
             sp_tmp(4*nlev+1,wn) +        &
             bcov_sp_a(nlev+1,k)*         &
             DSQRT(eigv_a(k,wn))*         &
             eigen_a(k,wn)
       END DO
!       sp_tmp(1:nlev,wn)=sp_tmp(3*nlev+1:4*nlev,wn)
!       sp_tmp(nlev+1:2*nlev,wn)=sp_tmp(3*nlev+1:4*nlev,wn)
!       sp_tmp(2*nlev+1:3*nlev,wn)=sp_tmp(3*nlev+1:4*nlev,wn)
!       sp_tmp(4*nlev+1,wn)=sp_tmp(4*nlev,wn)
       sp_tmp(1:nlev,wn)=sp_tmp(2*nlev+1:3*nlev,wn)
       sp_tmp(nlev+1:2*nlev,wn)=sp_tmp(2*nlev+1:3*nlev,wn)
       sp_tmp(3*nlev+1:4*nlev,wn)=sp_tmp(2*nlev+1:3*nlev,wn)
    END DO
    ierr = TimingStop('Run BeforeBcast_fromEigToSpec')

    ierr = TimingStart('Run Bcast4SpecToCtrl_fromEigToSpec')
    CALL Mpi_AllGatherV(sp_tmp,               &
                        nlevar*l_zwn_nproc,   &
                        par_double_precision, &
                        sp,                   &
                        nlevar*l_zwn_blength, &
                        nlevar*l_zwn_bstart,  &
                        par_double_precision, &
                        par%comm,             &
                        ierr)
    ierr = TimingStop('Run Bcast4SpecToCtrl_fromEigToSpec')

    RETURN

  END SUBROUTINE EigToSpecEns

!===============================================================================
!  SUBROUTINE TlmEigToObs
!
!> @brief
!> - Tangent linear operator of eigen space -> observation space
!
!> @date 23JUL2013
!> - HJ SONG: Design and code
!> @date 30MAY2016
!> - HJ SONG: Modify to Hybrid 4DEnVAR
!
!> @param[in] eigen3d, eigen2d, tl_eigen3d, tl_eigen2d
!> @param[out] tl_sd, tl_sf
!===============================================================================
  SUBROUTINE TlmEigToObs(l_gs_md,        &
                         tl_eigen,       &
                         tl_eigen_chi,   &
                         tl_eigen_Tps,   &
                         tl_eigen_q,     &
                         tl_eigen_a,     &
                         l_tl_obs)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules

  ! 2. Variables
  ! 2-1. Input
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime), &
       INTENT(IN) :: l_gs_md
    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc), &
       INTENT(IN) :: tl_eigen
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc), &
       INTENT(IN) :: tl_eigen_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc), &
       INTENT(IN) :: tl_eigen_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc), &
       INTENT(IN) :: tl_eigen_q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl), &
       INTENT(IN) :: tl_eigen_a

  ! 2-2. Output
    TYPE(obsbg_type), DIMENSION(ntime), &
       INTENT(OUT) :: l_tl_obs

  ! 2-3. Local
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
       :: tl_ct, tl_ctt,                                 &
          tl_ct_a
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime) &
       :: tl_md
    REAL(KIND=real_kind), DIMENSION(nlevar,zwn2) &
       :: tl_sp
    REAL(KIND=real_kind), DIMENSION(nlevar,zwn2_a) &
       :: tl_sp_a
    INTEGER(KIND=int_kind) &
       :: i, j, k, l, ie,  &
          wn, vr,          & 
          iproc, ierr,     &
          itime
  !=============================================================================

  !=============================================================================
  ! B. Main body
  !=============================================================================
    CALL TlmEigToModel(l_gs_md(:,:,:,:,:), &
                       tl_eigen,           &
                       tl_eigen_chi,       &
                       tl_eigen_Tps,       &
                       tl_eigen_q,         &
                       tl_eigen_a,         &
                       tl_md(:,:,:,:,:))
  
    DO itime = 1, ntime
      ierr = TimingStart('Run TlmModelToObs_fromTlmEigToObs')
      CALL TlmModelToObs(itime,                 &
                         l_gs_md(:,:,:,:,itime),&     ! in
                         tl_md(:,:,:,:,itime),  &     ! inout
                         l_tl_obs(itime))             ! out
      ierr = TimingStop('Run TlmModelToObs_fromTlmEigToObs')
    END DO ! itime

    RETURN

  END SUBROUTINE TlmEigToObs

!===============================================================================
!  SUBROUTINE TlmEigToModel
!
!> @brief
!> - Tangent linear operator of eigen space -> observation space
!
!> @date 23JUL2013
!> - HJ SONG: Design and code
!> @date 30MAY2015
!> - HJ SONG: Modify to Hybrid 4DEnVAR
!
!> @param[in] eigen3d, eigen2d, tl_eigen3d, tl_eigen2d
!> @param[out] tl_sd, tl_sf
!===============================================================================
  SUBROUTINE TlmEigToModel(l_gs_md,        &
                           tl_eigen,       &
                           tl_eigen_chi,   &
                           tl_eigen_Tps,   &
                           tl_eigen_q,     &
                           tl_eigen_a,     &
                           tl_md)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules

  ! 2. Variables
  ! 2-1. Input
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime), &
       INTENT(IN) :: l_gs_md
    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc), &
       INTENT(IN) :: tl_eigen
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc), &
       INTENT(IN) :: tl_eigen_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc), &
       INTENT(IN) :: tl_eigen_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc), &
       INTENT(IN) :: tl_eigen_q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl), &
       INTENT(IN) :: tl_eigen_a

  ! 2-2. Output
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime), &
       INTENT(OUT) :: tl_md

  ! 2-3. Local
  ! up
!    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
    REAL(KIND=real_kind), DIMENSION(nUPs_l,nlevar) &
       :: tl_ct, tl_ctt,                           &
          tl_ct_a
  ! up
    !REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime) &
    REAL(KIND=real_kind), DIMENSION(nUPs_l,nlevar,ntime) &
       :: tl_ct_hyb
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
       :: tl_ct_hyb_ep
    REAL(KIND=real_kind), DIMENSION(nlevar,zwn2) &
       :: tl_sp
    REAL(KIND=real_kind), DIMENSION(nlevar,zwn2_a) &
       :: tl_sp_a
    INTEGER(KIND=int_kind) &
       :: i, j, k, l, ie,  &
          wn, vr,          & 
          iproc, ierr,     &
          itime, ismpl
  !=============================================================================

  !=============================================================================
  ! B. Main body
  !=============================================================================

!---------------------------------------------------------
! Climatological Part
    tl_sp  = 0.D0
    tl_ct  = 0.D0
    tl_ctt = 0.D0

    ierr = TimingStart('Run EigToSpec_fromTlmEigToObs')
    CALL EigToSpec(zwn_nproc,    &
                   zwn2,         &
                   zwn_blength,  &
                   zwn_bstart,   &
                   tl_eigen,     &
                   tl_eigen_chi, &
                   tl_eigen_Tps, &
                   tl_eigen_q,   &
                   tl_sp)
    ierr = TimingStop('Run EigToSpec_fromTlmEigToObs')

    ierr = TimingStart('Run SpecToCtrl_fromTlmEigToModel')
    CALL SpecToCtrl(zwn2,  &
                    tl_sp, &
                    tl_ct)
    ierr = TimingStop('Run SpecToCtrl_fromTlmEigToModel')

    !up
    !tl_ctt = tl_ct *DSQRT(errvar)
    tl_ctt = tl_ct *DSQRT(errvar_up)

    DO itime = 1, ntime
       !tl_ct_hyb(:,:,:,:,itime) = DSQRT(b_clim)*tl_ctt(:,:,:,:)
    DO k = 1, nlevar
       ! up
       !tl_ct_hyb(:,:,k,:,itime) = DSQRT(b_clim_v(k))*tl_ctt(:,:,k,:)
       tl_ct_hyb(:,k,itime) = DSQRT(b_clim_v(k))*tl_ctt(:,k)
    END DO
    END DO
!-----------------------------------------------------------

!---------------------------------------------------------------
! Add Ensemble Part
    DO ismpl = 1, nsmpl
       ierr = TimingStart('Run EigToSpecEns_fromTlmEigToModel')
       CALL EigToSpecEns(zwn_nproc_a,               &
                         zwn2_a,                    &
                         zwn_blength_a,             &
                         zwn_bstart_a,              &
                         tl_eigen_a(:,:,ismpl),     &
                         tl_sp_a(:,:))
       ierr = TimingStop('Run EigToSpecEns_fromTlmEigToModel')

       ierr = TimingStart('Run SpecToCtrl_fromTlmEigToModel')
       CALL SpecToCtrl(zwn2_a,  &
                       tl_sp_a, &
                       tl_ct_a)
       ierr = TimingStop('Run SpecToCtrl_fromTlmEigToModel')

    DO itime = 1, ntime
    DO k = 1, nlevar
       tl_ct_hyb(:,k,itime) = tl_ct_hyb(:,k,itime) + &
                              DSQRT(b_ens_v(k))*bg_smpl(:,k,itime,ismpl)*tl_ct_a(:,k)
    END DO ! k
    END DO ! itime
    END DO ! ismpl
!------------------------------------------------------------------

    DO itime = 1, ntime
!---------------Parameter Transform-------------------
      ierr = TimingStart('Run TlmCtrlToModel_fromTlmEigToModel')

      ierr = TimingStart('UP2EP_fromTlmEigToModel')
      CALL ConvertUP2EP_hyb(nlevar,tl_ct_hyb(:,:,itime),tl_ct_hyb_ep)
      ierr = TimingStop('UP2EP_fromTlmEigToModel')

      CALL TlmCtrlToModel(l_gs_md(:,:,:,:,itime),     &    ! in
                          bg_uv_rot(:,:,:,:,:,itime), &    ! in
                          tl_ct_hyb_ep(:,:,:,:),      &    ! in
                          tl_md(:,:,:,:,itime))            ! out
      ierr = TimingStop('Run TlmCtrlToModel_fromTlmEigToModel')
!-----------------------------------------------------
    END DO ! itime
  
    RETURN

  END SUBROUTINE TlmEigToModel

!===============================================================================
!  SUBROUTINE AdjEigToObs
!
!> @brief
!> - Adjoint operator of eigen space -> observation space
!
!> @date 23JUL2013
!> - HJ SONG: Design and code
!> @date 30MAY2016
!> - HJ SONG: A version for Hybrid 4DEnVAR
!
!> @param[in] eigen3d, eigen2d
!> @param[inout] ad_sd, ad_eigen3d, ad_eigen2d
!===============================================================================
  SUBROUTINE AdjEigToObs(l_gs_md,        &
                         l_ad_obs,       &
                         ad_eigen,       & ! for psi
                         ad_eigen_chi,   & ! for chi
                         ad_eigen_Tps,   & ! for (T,ps)
                         ad_eigen_q,     &
                         ad_eigen_a      &
                         )     

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules

  ! 2. Variables
  ! 2-1. Input
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime), &
       INTENT(IN) :: l_gs_md

  ! 2-2. Input/Output
    TYPE(obsbg_type), DIMENSION(ntime), &
       INTENT(INOUT) :: l_ad_obs
    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc), &
       INTENT(INOUT) :: ad_eigen                           ! adjoint variable for psi
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc), &
       INTENT(INOUT) :: ad_eigen_chi                       ! adjoint variable for chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc), &
       INTENT(INOUT) :: ad_eigen_Tps                       ! adjoint variable for (T,ps)
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc), &
       INTENT(INOUT) :: ad_eigen_q                         ! adjoint variable for q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl), &
       INTENT(INOUT) :: ad_eigen_a                         ! adjoint variable for alpha

  ! 2-3. Local
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
       :: ad_ct, ad_ctt
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime) &
       :: ad_md
    REAL(KIND=real_kind), DIMENSION(nlevar,zwn_nproc) &
       :: ad_sp
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie,     &
          wn, vr,          & 
          iproc, ierr,     &
          itime
  !=============================================================================

  !=============================================================================
  ! B. Prepare a basic state
  !=============================================================================
!   md = gs_md
!    ad_ct = 0.D0
  !=============================================================================

   ad_md   = 0.D0
   DO itime = ntime, 1, -1
    !=============================================================================
    ! C. Main body of adjoint process
    !=============================================================================
      ierr = TimingStart('Run AdjModelToObs_fromAdjEigToObs')
      CALL AdjModelToObs(itime,                  &
                         l_gs_md(:,:,:,:,itime), &    ! in
                         l_ad_obs(itime),        &    ! inout
                         ad_md(:,:,:,:,itime))       ! inout
      ierr = TimingStop('Run AdjModelToObs_fromAdjEigToObs')
    END DO ! itime

    CALL AdjEigToModel(l_gs_md(:,:,:,:,:), &
                       ad_md(:,:,:,:,:),   &
                       ad_eigen,           &
                       ad_eigen_chi,       &
                       ad_eigen_Tps,       &
                       ad_eigen_q,         &
                       ad_eigen_a          &
                       )
   
    RETURN

  END SUBROUTINE AdjEigToObs



!===============================================================================
!  SUBROUTINE AdjEigToModel
!
!> @brief
!> - Adjoint operator of eigen space -> observation space
!
!> @date 23JUL2013
!> - HJ SONG: Design and code
!> @date 30MAY2016
!> - HJ SONG: A version for Hybrid 4DEnVAR
!
!> @param[in] eigen3d, eigen2d
!> @param[inout] ad_sd, ad_eigen3d, ad_eigen2d
!===============================================================================
  SUBROUTINE AdjEigToModel(l_gs_md,        &
                           ad_md,          &
                           ad_eigen,       & ! for psi
                           ad_eigen_chi,   & ! for chi
                           ad_eigen_Tps,   & ! for (T,ps)
                           ad_eigen_q,     &
                           ad_eigen_a      &
                           )    

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules

  ! 2. Variables
  ! 2-1. Input
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime), &
       INTENT(IN) :: l_gs_md

  ! 2-2. Input/Output
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime), &
       INTENT(INOUT) :: ad_md
    REAL(KIND=real_kind), DIMENSION(neig,zwn_nproc), &
       INTENT(INOUT) :: ad_eigen                           ! adjoint variable for psi
    REAL(KIND=real_kind), DIMENSION(neig_chi,zwn_nproc), &
       INTENT(INOUT) :: ad_eigen_chi                       ! adjoint variable for chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,zwn_nproc), &
       INTENT(INOUT) :: ad_eigen_Tps                       ! adjoint variable for (T,ps)
    REAL(KIND=real_kind), DIMENSION(neig_q,zwn_nproc), &
       INTENT(INOUT) :: ad_eigen_q                         ! adjoint variable for q

    REAL(KIND=real_kind), DIMENSION(nlev+1,zwn_nproc_a,nsmpl), &
       INTENT(INOUT) :: ad_eigen_a                         ! adjoint variable for alpha

  ! 2-3. Local
    !up
    !REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
    REAL(KIND=real_kind), DIMENSION(nUPs_l,nlevar) &
       :: ad_ct, ad_ctt,                                 &
          ad_ct_a
    !up
    !REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime) &
    REAL(KIND=real_kind), DIMENSION(nUPs_l,nlevar,ntime) &
       :: ad_ct_hyb
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
       :: ad_ct_hyb_ep
    REAL(KIND=real_kind), DIMENSION(nlevar,zwn_nproc) &
       :: ad_sp,                                      &
          ad_sp_a
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie,     &
          wn, vr,          & 
          iproc, ierr,     &
          itime, ismpl
  !=============================================================================

  !=============================================================================
  ! B. Prepare a basic state
  !=============================================================================
    DO itime = ntime, 1, -1 
      !up
      !ad_ct_hyb(:,:,:,:,itime)  = 0.D0
      ad_ct_hyb_ep  = 0.D0
      ierr = TimingStart('Run AdjCtrlToModel_fromAdjEigToModel')
      CALL AdjCtrlToModel(l_gs_md(:,:,:,:,itime),     &
                          bg_uv_rot(:,:,:,:,:,itime), &   ! in
                          ad_md(:,:,:,:,itime),       &   ! inout
                          ad_ct_hyb_ep(:,:,:,:))          ! inout
      !up
      ierr = TimingStart('EP2UP_fromAdjEigToModel')
      CALL ConvertEP2EP0_hyb(nlevar,ad_ct_hyb_ep)
      CALL ConvertEP2UP_hyb(nlevar,ad_ct_hyb_ep,ad_ct_hyb(:,:,itime))
      ierr = TimingStop('EP2UP_fromAdjEigToModel')

      ierr = TimingStop('Run AdjCtrlToModel_fromAdjEigToModel')
    END DO ! itime

!--------------------------------------------------
! Adjoint of adding ensemble part
    DO ismpl = nsmpl, 1, -1
       ad_ct_a = 0.D0
    DO itime = ntime, 1, -1
    DO k = nlevar, 1, -1
       ad_ct_a(:,k) = ad_ct_a(:,k) + &
                 DSQRT(b_ens_v(k))*bg_smpl(:,k,itime,ismpl)*ad_ct_hyb(:,k,itime)
    END DO
    END DO

       ad_sp_a = 0.D0
       ierr = TimingStart('Run AdjSpecToCtrl_fromAdjEigToModel')
       CALL AdjSpecToCtrl_opt(zwn2_a,        &
                              zwn_nproc_a,   &
                              zwn_iproc_a,   &
                              zwn_bstart_a,  &
                              zwn_blength_a, &
                              ad_ct_a,       &
                              ad_sp_a)
       ierr = TimingStop('Run AdjSpecToCtrl_fromAdjEigToModel')
   
       ierr = TimingStart('Run AfterAdjSpecToCtrl_fromAdjEigToModel')
       CALL AdjEigToSpecEns(zwn_nproc_a,               &
                            ad_sp_a,                   &
                            ad_eigen_a(:,:,ismpl))
       ierr = TimingStop('Run AfterAdjSpecToCtrl_fromAdjEigToModel')
    END DO ! ismpl
!--------------------------------------------------------------------------

!---------------------------------------------------------------------
! Adjoint of climatological part
    ad_ctt = 0.D0
    DO itime = ntime, 1, -1
    DO k = nlevar, 1, -1
       ad_ctt(:,k) = ad_ctt(:,k) + DSQRT(b_clim_v(k))*ad_ct_hyb(:,k,itime)
       ad_ct_hyb(:,k,itime) = 0.D0
    END DO ! k
    END DO ! itime

    ad_ct = ad_ctt*DSQRT(errvar_up)

    ad_sp = 0.D0
    ierr = TimingStart('Run AdjSpecToCtrl_fromAdjEigToModel')
    CALL AdjSpecToCtrl_opt(zwn2,        &
                           zwn_nproc,   &
                           zwn_iproc,   &
                           zwn_bstart,  &
                           zwn_blength, &
                           ad_ct,       &
                           ad_sp)
    ierr = TimingStop('Run AdjSpecToCtrl_fromAdjEigToModel')

    ierr = TimingStart('Run AfterAdjSpecToCtrl_fromAdjEigToModel')
    CALL AdjEigToSpec(zwn_nproc,    &
                      ad_sp,        &
                      ad_eigen,     &
                      ad_eigen_chi, &
                      ad_eigen_Tps, &
                      ad_eigen_q)
    ierr = TimingStop('Run AfterAdjSpecToCtrl_fromAdjEigToModel')
!------------------------------------------------------------------------

    RETURN

  END SUBROUTINE AdjEigToModel




!===============================================================================
  SUBROUTINE AdjEigToSpec(l_zwn_nproc,  &
                          ad_sp,        &
                          ad_eigen,     &
                          ad_eigen_chi, &
                          ad_eigen_Tps, &
                          ad_eigen_q)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules

  ! 2. Variables
  ! 2-1. Input
    INTEGER(KIND=int_kind), &
       INTENT(IN) :: l_zwn_nproc
    REAL(KIND=real_kind), DIMENSION(nlevar,l_zwn_nproc), &
       INTENT(IN) :: ad_sp

  ! 2-2. Input/Output
    REAL(KIND=real_kind), DIMENSION(neig,l_zwn_nproc), &
       INTENT(INOUT) :: ad_eigen
    REAL(KIND=real_kind), DIMENSION(neig_chi,l_zwn_nproc), &
       INTENT(INOUT) :: ad_eigen_chi
    REAL(KIND=real_kind), DIMENSION(neig_Tps,l_zwn_nproc), &
       INTENT(INOUT) :: ad_eigen_Tps
    REAL(KIND=real_kind), DIMENSION(neig_q,l_zwn_nproc), &
       INTENT(INOUT) :: ad_eigen_q

  ! 2-3. Local
    INTEGER(KIND=int_kind) &
       :: k,               &
          wn
  !=============================================================================

  !=============================================================================
  ! C. Main body of adjoint process
  !=============================================================================
    DO wn = 1, l_zwn_nproc
       DO k  = 1, neig
          ad_eigen(k,wn) =                    &
             DSQRT(eigv(k,wn))*               &
             SUM(bcov_sp(1:nlev,k,wn)*        &
                 ad_sp(1:nlev,wn))
          ! Check
          !IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn) 
       END DO
       DO k  = 1, neig_chi
          ad_eigen_chi(k,wn) =        &
             DSQRT(eigv_chi(k,wn))*   &
             SUM(bcov_sp_chi(:,k,wn)* &
                 ad_sp(nlev+1:2*nlev,wn))
          ! Check
          !IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn) 
       END DO
       DO k  = 1, neig_Tps
          ad_eigen_Tps(k,wn) =                          &
             DSQRT(eigv_Tps(k,wn))*                     &
             (SUM(bcov_sp_Tps(1:nlev,k,wn)*             &
                  ad_sp(2*nlev+1:3*nlev,wn))+ &
              bcov_sp_Tps(nlev+1,k,wn)*                 &
              ad_sp(4*nlev+1,wn))
          ! Check
          !IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn) 
       END DO
       DO k  = 1, neig_q
          ad_eigen_q(k,wn) =        &
             DSQRT(eigv_q(k,wn))*   &
             SUM(bcov_sp_q(:,k,wn)* &
                 ad_sp(3*nlev+1:4*nlev,wn))
          ! Check
          !IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn) 
       END DO
    END DO

    RETURN

  END SUBROUTINE AdjEigToSpec


!===============================================================================
! for alpha-variable (HJS)
  SUBROUTINE AdjEigToSpecEns(l_zwn_nproc,  &
                             ad_sp,        &
                             ad_eigen_a)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules

  ! 2. Variables
  ! 2-1. Input
    INTEGER(KIND=int_kind), &
       INTENT(IN) :: l_zwn_nproc
    REAL(KIND=real_kind), DIMENSION(nlevar,l_zwn_nproc), &
       INTENT(INOUT) :: ad_sp

  ! 2-2. Input/Output
    REAL(KIND=real_kind), DIMENSION(nlev+1,l_zwn_nproc), &
       INTENT(INOUT) :: ad_eigen_a

  ! 2-3. Local
    INTEGER(KIND=int_kind) &
       :: k,               &
          wn
  !=============================================================================

  !=============================================================================
  ! C. Main body of adjoint process
  !=============================================================================
    DO wn = l_zwn_nproc, 1, -1
       ad_sp(2*nlev+1:3*nlev,wn)=ad_sp(2*nlev+1:3*nlev,wn) + &
                                 ad_sp(1:nlev,wn)
       ad_sp(1:nlev,wn)=0.D0
       ad_sp(2*nlev+1:3*nlev,wn)=ad_sp(2*nlev+1:3*nlev,wn) + &
                                 ad_sp(nlev+1:2*nlev,wn)
       ad_sp(nlev+1:2*nlev,wn)=0.D0
       ad_sp(2*nlev+1:3*nlev,wn)=ad_sp(2*nlev+1:3*nlev,wn) + &
                                 ad_sp(3*nlev+1:4*nlev,wn)
       ad_sp(3*nlev+1:4*nlev,wn)=0.D0

!       DO k  = 1, nlev
!          ad_eigen_a(k,wn) =      &
!             DSQRT(eigv_a(k,wn))* &
!             SUM(bcov_sp_a(:,k)*  &
!                 ad_sp(3*nlev+1:4*nlev,wn))
!          ! Check
!          !IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn) 
!       END DO

       DO k  = 1, nlev+1
          ad_eigen_a(k,wn) =                  &
             DSQRT(eigv_a(k,wn))*             &
             (SUM(bcov_sp_a(1:nlev,k)*        &
                  ad_sp(2*nlev+1:3*nlev,wn))+ &
              bcov_sp_a(nlev+1,k)*            &
              ad_sp(4*nlev+1,wn))
!          ! Check
!          !IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn) 
       END DO

    END DO

    RETURN

  END SUBROUTINE AdjEigToSpecEns

END MODULE CostFunMinMod
