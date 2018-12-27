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
!-------------------------------------------------------------------------------
   module costfunminmod
!===============================================================================
! A. Common modules and variables
!===============================================================================
! 1. Common modules
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
   use kiapsbase, only: int_kind=>kim_int_kind, real_kind=>kim_real8_kind, iulog=>kim_iu_log
   use dimensions, only: np, nlev, nelemd, nets, nete, nvar3d, nvar2d, nlevar, zwn2, zwn2_a, ! for alpha
!
   neig,neig_chi,neig_tps,neig_q,zwn_nproc,zwn_iproc,zwn_nproc_a,zwn_iproc_a,! for alpha
   ntime,! # of time bins
   nt_itime ! itime for nature
   use kiapsparallel, only: par=>kim_par, hybrid=>kim_hybrid, global_shared_buf, global_shared_sum, par_double_precision, par_character, par_sum, par_max, par_min
   use kiapsparallel, only: barrierpar
   use kiapsparmpi, only: par_reduce
   use kiapswave, only: zwn_blength, zwn_bstart, zwn_blength_a, ! for alpha
!
   zwn_bstart_a ! for alpha
   use globalnorms, only: wrap_repro_sum
   use timing, only: timingstart, timingstop
   use backgroundmod, only: temptom, bg_md, nt_org, bg_uv_rot, !up
!errvar,                &
!
   errvar_up,bg_eig,bcov_sp,! eigenvectors of b for psi
   eigv,! eigenvalues of b for psi
   bcov_sp_chi,! eigenvectors of b for chi
   eigv_chi,! eigenvalues of b for chi
   bcov_sp_tps,! eigenvectors of b for(t,ps)
   eigv_tps,! eigenvalues of b for(t,ps)
   bcov_sp_q,! eigenvectors of b for q
   eigv_q,! eigenvalues of b for q
   bcov_sp_a,! eigenvectors of b for alpha
   eigv_a,! eigenvalues of b for alpha
   nsmpl,cv_opt_hum,is_moist_mr,exst_nt,bg_smpl,b_clim_v,b_ens_v
   use observationmod, only: obsmet, sdmiss, obcheck, obcheck_outer, crit_omb_amsua, crit_omb_gpsro, crit_omb_iasi, crit_omb_cris, crit_omb_atms, crit_omb_atmswv, crit_omb_mhs, crit_omb_csr, singleobs, da_sonde, da_surface, da_bogus, da_aircraft, da_amv, da_scatwind, da_amsua, da_gpsro, da_iasi, da_cris, da_atms, da_atmswv, da_mhs, da_csr, rttovused, roppused
   use obsallocatemod, only: obsbg_type, obs_type
   use calrttovmod, only: calrttov
   use calroppmod, only: calropp
   use readobsmod, only: obs
   use ctrltospecmod, only: reddof
   use spectoctrlmod, only: spectoctrl, adjspectoctrl_opt
   use ctrltomodelmod, only: tlmctrltomodel, adjctrltomodel, phitotemp
   use modeltoobsmod, only: modeltoobs
   use tlmmodeltoobsmod, only: tlmmodeltoobs
   use adjmodeltoobsmod, only: adjmodeltoobs
   use makebmod, only: makebsingle, makebsonde, makebsurface, makebaircraft, makebamv, makebscatwind, makebamsua, makebgpsro, makebiasi, makebcris, makebatms, makebmhs, makebcsr, makebbogus, makebatmswv
   use multiplyrmod, only: multiplyrsingle, multiplyrsonde, multiplyrsurface, multiplyraircraft, multiplyramv, multiplyrscatwind, multiplyramsua, multiplyrgpsro, multiplyriasi, multiplyrcris, multiplyratms, multiplyrmhs, multiplyrcsr, multiplyrbogus, multiplyratmswv
!  USE ModelToEigMod,  ONLY : TlmModelToEig
   use kiapsgrid, only: nups_l, convertup2ep_hyb, convertep2up_hyb, convertep2ep0_hyb
! 2. Common variables
   implicit none
!
   private
!===============================================================================
!===============================================================================
! B. Public modules and variables
!===============================================================================
! 1. Public modules
   public :: mincongraspec
!!  PUBLIC :: EigToModel
! 2. Public variables
   real(real_kind), dimension(:,:,:,:,:), allocatable, public :: an_md ! analysis
   real(real_kind), public :: cgtol1, cgtol2
   integer(int_kind), public :: noutl, niter, optmethod
   type(obsbg_type), dimension(:), allocatable, ! for time dimension
   public :: gs_obs, ad_obs, tl_obs, bg_obs, an_obs
   type(obs_type), dimension(:), allocatable, ! for time dimension
   public :: obs_org ! for observation check every outer-loop
   real(real_kind), dimension(:,:), allocatable, public :: x_eig_acc, x_eig_acc_chi, x_eig_acc_tps, x_eig_acc_q
   real(real_kind), dimension(:,:,:), allocatable, public :: x_eig_acc_a ! for alpha-variable
   integer(int_kind) :: ierr
!===============================================================================
!
   contains
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
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine mincongraspec
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use edge, only: edgebuffer_t, initedgebuffer, edgevpack, edgevunpack
   use bndry, only: bndry_exchangev
! 2. Variables
! 2-2. Local
!
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime) :: gs_md, x_md
   real(real_kind), dimension(neig, zwn_nproc) :: x_eig
   real(real_kind), dimension(neig_chi, zwn_nproc) :: x_eig_chi
   real(real_kind), dimension(neig_tps, zwn_nproc) :: x_eig_tps
   real(real_kind), dimension(neig_q, zwn_nproc) :: x_eig_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl) :: x_eig_a ! for alpha-variable
   integer(int_kind) :: i, j, k, kk, ie, wn, vr, iter, ioutl, iproc, ierr, err, itime
!=============================================================================
! Extended to 4D
!! ALLOCATE Analysis
!
   allocate(an_md(np,np,nlevar,nelemd,ntime))
   allocate(gs_obs(ntime),tl_obs(ntime),ad_obs(ntime),bg_obs(ntime),an_obs(ntime),obs_org(ntime))
   do itime = 1,ntime
!! ALLOCATE observation information
     if (singleobs) then
       allocate(gs_obs(itime)%single(obsmet(itime)%nsg_e_max,nelemd))
       allocate(tl_obs(itime)%single(obsmet(itime)%nsg_e_max,nelemd))
       allocate(ad_obs(itime)%single(obsmet(itime)%nsg_e_max,nelemd))
       allocate(bg_obs(itime)%single(obsmet(itime)%nsg_e_max,nelemd))
       allocate(an_obs(itime)%single(obsmet(itime)%nsg_e_max,nelemd))
       allocate(obs_org(itime)%single(obsmet(itime)%nsd_e_max,nelemd))
       obs_org(itime)%single(:,:) = obs(itime)%single(:,:)
     endif
     if (da_sonde) then
       allocate(gs_obs(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
       allocate(tl_obs(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
       allocate(ad_obs(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
       allocate(bg_obs(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
       allocate(an_obs(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
       allocate(obs_org(itime)%sonde(obsmet(itime)%nsd_e_max,nelemd))
       obs_org(itime)%sonde(:,:) = obs(itime)%sonde(:,:)
     endif
     if (da_surface) then
       allocate(gs_obs(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
       allocate(tl_obs(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
       allocate(ad_obs(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
       allocate(bg_obs(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
       allocate(an_obs(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
       allocate(obs_org(itime)%surface(obsmet(itime)%nsf_e_max,nelemd))
       obs_org(itime)%surface(:,:) = obs(itime)%surface(:,:)
     endif
     if (da_bogus) then
       allocate(gs_obs(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
       allocate(tl_obs(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
       allocate(ad_obs(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
       allocate(bg_obs(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
       allocate(an_obs(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
       allocate(obs_org(itime)%bogus(obsmet(itime)%nbg_e_max,nelemd))
       obs_org(itime)%bogus(:,:) = obs(itime)%bogus(:,:)
     endif
     if (da_aircraft) then
       allocate(gs_obs(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
       allocate(tl_obs(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
       allocate(ad_obs(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
       allocate(bg_obs(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
       allocate(an_obs(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
       allocate(obs_org(itime)%aircraft(obsmet(itime)%nar_e_max,nelemd))
       obs_org(itime)%aircraft(:,:) = obs(itime)%aircraft(:,:)
     endif
     if (da_amv) then
       allocate(gs_obs(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
       allocate(tl_obs(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
       allocate(ad_obs(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
       allocate(bg_obs(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
       allocate(an_obs(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
       allocate(obs_org(itime)%amv(obsmet(itime)%namv_e_max,nelemd))
       obs_org(itime)%amv(:,:) = obs(itime)%amv(:,:)
     endif
     if (da_scatwind) then
       allocate(gs_obs(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
       allocate(tl_obs(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
       allocate(ad_obs(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
       allocate(bg_obs(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
       allocate(an_obs(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
       allocate(obs_org(itime)%scatwind(obsmet(itime)%nscatwind_e_max,nelemd))
       obs_org(itime)%scatwind(:,:) = obs(itime)%scatwind(:,:)
     endif
     if (da_amsua) then
       allocate(gs_obs(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
       allocate(tl_obs(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
       allocate(ad_obs(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
       allocate(bg_obs(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
       allocate(an_obs(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
       allocate(obs_org(itime)%amsua(obsmet(itime)%namsua_e_max,obsmet(itime)%namsuach,nelemd))
       obs_org(itime)%amsua(:,:,:) = obs(itime)%amsua(:,:,:)
     endif
     if (da_gpsro) then
       allocate(gs_obs(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
       allocate(tl_obs(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
       allocate(ad_obs(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
       allocate(bg_obs(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
       allocate(an_obs(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
       allocate(obs_org(itime)%gpsro(obsmet(itime)%ngpsro_e_max,nelemd))
       obs_org(itime)%gpsro(:,:) = obs(itime)%gpsro(:,:)
     endif
     if (da_iasi) then
       allocate(gs_obs(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
       allocate(tl_obs(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
       allocate(ad_obs(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
       allocate(bg_obs(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
       allocate(an_obs(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
       allocate(obs_org(itime)%iasi(obsmet(itime)%niasi_e_max,obsmet(itime)%niasich,nelemd))
       obs_org(itime)%iasi(:,:,:) = obs(itime)%iasi(:,:,:)
     endif
     if (da_cris) then
       allocate(gs_obs(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
       allocate(tl_obs(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
       allocate(ad_obs(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
       allocate(bg_obs(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
       allocate(an_obs(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
       allocate(obs_org(itime)%cris(obsmet(itime)%ncris_e_max,obsmet(itime)%ncrisch,nelemd))
       obs_org(itime)%cris(:,:,:) = obs(itime)%cris(:,:,:)
     endif
     if (da_atms) then
       allocate(gs_obs(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
       allocate(tl_obs(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
       allocate(ad_obs(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
       allocate(bg_obs(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
       allocate(an_obs(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
       allocate(obs_org(itime)%atms(obsmet(itime)%natms_e_max,obsmet(itime)%natmsch,nelemd))
       obs_org(itime)%atms(:,:,:) = obs(itime)%atms(:,:,:)
     endif
     if (da_atmswv) then
       allocate(gs_obs(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
       allocate(tl_obs(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
       allocate(ad_obs(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
       allocate(bg_obs(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
       allocate(an_obs(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
       allocate(obs_org(itime)%atmswv(obsmet(itime)%natmswv_e_max,obsmet(itime)%natmswvch,nelemd))
       obs_org(itime)%atmswv(:,:,:) = obs(itime)%atmswv(:,:,:)
     endif
     if (da_mhs) then
       allocate(gs_obs(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
       allocate(tl_obs(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
       allocate(ad_obs(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
       allocate(bg_obs(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
       allocate(an_obs(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
       allocate(obs_org(itime)%mhs(obsmet(itime)%nmhs_e_max,obsmet(itime)%nmhsch,nelemd))
       obs_org(itime)%mhs(:,:,:) = obs(itime)%mhs(:,:,:)
     endif
     if (da_csr) then
       allocate(gs_obs(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
       allocate(tl_obs(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
       allocate(ad_obs(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
       allocate(bg_obs(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
       allocate(an_obs(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
       allocate(obs_org(itime)%csr(obsmet(itime)%ncsr_e_max,obsmet(itime)%ncsrch,nelemd))
       obs_org(itime)%csr(:,:,:) = obs(itime)%csr(:,:,:)
     endif
   enddo ! itime
!=============================================================================
! OUTER LOOP
!=============================================================================
   allocate(x_eig_acc(neig,zwn_nproc))
   allocate(x_eig_acc_chi(neig_chi,zwn_nproc))
   allocate(x_eig_acc_tps(neig_tps,zwn_nproc))
   allocate(x_eig_acc_q(neig_q,zwn_nproc))
   allocate(x_eig_acc_a(nlev+1,zwn_nproc_a,nsmpl))
   x_eig_acc = 0.d0
   x_eig_acc_chi = 0.d0
   x_eig_acc_tps = 0.d0
   x_eig_acc_q = 0.d0
   x_eig_acc_a = 0.d0 ! alpha
   x_eig = 0.d0
   gs_md = bg_md
   x_eig_chi = 0.d0
   x_eig_tps = 0.d0
   x_eig_q = 0.d0
   x_eig_a = 0.d0 ! alpha
   do ioutl = 1,noutl
     if (par%ismasterproc) print*,'ioutl:',ioutl
     call mpi_barrier(par%comm,ierr)
!==========================================================================
!==========================================================================
! A. RCG
!==========================================================================
     if (optmethod.eq.1) call rcgmethod(ioutl,gs_md,x_eig,x_eig_chi,x_eig_tps,x_eig_q,x_eig_a) ! alpha
!==========================================================================
!==========================================================================
! B. PreConditioned PRCG
!==========================================================================
!       IF(optmethod .eq. 2) CALL PRCGmethod(ioutl,gs_md,x_eig)
!==========================================================================
! C. L-BFGS
!==========================================================================
!       IF (optmethod .eq. 3) CALL LBFGSmethod(ioutl,gs_md,x_eig)
!==========================================================================
     x_md = 0.d0
     ierr = timingstart('TlmEigToModel_outer') !---profiling
     call tlmeigtomodel(gs_md(:,:,:,:,:),x_eig,x_eig_chi,x_eig_tps,x_eig_q,x_eig_a,x_md)
     ierr = timingstop('TlmEigToModel_outer') !---profiling
     gs_md(:,:,:,:,:) = gs_md(:,:,:,:,:)+x_md(:,:,:,:,:)
!==========================================================================
   enddo ! for outer loop
   deallocate(x_eig_acc)
   deallocate(x_eig_acc_chi)
   deallocate(x_eig_acc_tps)
   deallocate(x_eig_acc_q)
   deallocate(x_eig_acc_a)
!=============================================================================
!=============================================================================
! Transform eigen space to model space
!=============================================================================
   an_md = gs_md
!=============================================================================
   if (exst_nt) then
#define Check_RMSE_on
#ifdef Check_RMSE_on
     ierr = timingstart('Check_RMSE') !---profiling
     call check_rmse
     ierr = timingstop('Check_RMSE') !---profiling
#endif
   endif
#define Check_OB_on
#ifdef Check_OB_on
   ierr = timingstart('Check_OB') !---profiling
   do itime = 1,ntime
     call check_ob(itime)
   enddo
   ierr = timingstop('Check_OB') !---profiling
#endif
   deallocate(gs_obs,tl_obs,ad_obs,bg_obs,an_obs,obs_org)
!=============================================================================
! Convert specific humidity to mixing ratio
!=============================================================================
   if (is_moist_mr) then
     do itime = 1, ntime
       do ie = 1, nelemd
         do k = 1, nlev
           do j = 1, np
             do i = 1, np
               bg_md(i,j,k+3*nlev,ie,itime) = bg_md(i,j,k+3*nlev,ie,itime)/(1.d0-bg_md(i,j,k+3*nlev,ie,itime)) ! to mixing ratio
               an_md(i,j,k+3*nlev,ie,itime) = an_md(i,j,k+3*nlev,ie,itime)/(1.d0-an_md(i,j,k+3*nlev,ie,itime)) ! to mixing ratio
             enddo
           enddo
         enddo
       enddo
     enddo
   endif
   if (par%ismasterproc) then
     print*,'=======successfully finished 3dvar analysis======='
   endif
   return
!
   end subroutine mincongraspec
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine rcgmethod(ioutl, gs_md, x_eig, x_eig_chi, x_eig_tps, x_eig_q, x_eig_a) ! for alpha(hjs)
! 2. Variables
! 2-1. Output
!
   integer(int_kind), intent(in   ) :: ioutl
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(in   ) :: gs_md
   real(real_kind), dimension(neig, zwn_nproc), intent(inout) :: x_eig
   real(real_kind), dimension(neig_chi, zwn_nproc), intent(inout) :: x_eig_chi
   real(real_kind), dimension(neig_tps, zwn_nproc), intent(inout) :: x_eig_tps
   real(real_kind), dimension(neig_q, zwn_nproc), intent(inout) :: x_eig_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl), intent(inout) :: x_eig_a
! 2-2. Local
   real(real_kind), dimension(neig, zwn_nproc) ! for psi
   :: b_eig, r_eig, p_eig, ap_eig
   real(real_kind), dimension(neig_chi, zwn_nproc) ! for chi
   :: b_eig_chi, r_eig_chi, p_eig_chi, ap_eig_chi
   real(real_kind), dimension(neig_tps, zwn_nproc) ! for(t, ps)
   :: b_eig_tps, r_eig_tps, p_eig_tps, ap_eig_tps
   real(real_kind), dimension(neig_q, zwn_nproc) ! for q
   :: b_eig_q, r_eig_q, p_eig_q, ap_eig_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl) !--alpha
   :: b_eig_a, r_eig_a, p_eig_a, ap_eig_a
   real(real_kind), dimension(neig, zwn_nproc, niter) :: r_eig_history ! for psi
   real(real_kind), dimension(neig_chi, zwn_nproc, niter) :: r_eig_chi_history ! for chi
   real(real_kind), dimension(neig_tps, zwn_nproc, niter) :: r_eig_tps_history ! for(t, ps)
   real(real_kind), dimension(neig_q, zwn_nproc, niter) :: r_eig_q_history ! for q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl, niter) :: r_eig_q_history_a !--alpha
   real(real_kind) :: rsold, rsnew, rs0, pap, ! rs at initial
!
   alpha,beta,j_b,j_o,tmp
   integer(int_kind) :: i, j, k, kk, ie, wn, vr, iter, iproc, ierr, err, itime, ismpl
   real(real_kind) :: rjrk
   real(real_kind) :: tmpmax, tmpmin, maxv, minv
   integer(int_kind) :: info, nlanczos, nl, nw
   real(real_kind), dimension(niter) :: lanczos_diag! diagonal part of niter-dimensional lanczos tridiagonal matrix
   real(real_kind), dimension(niter-1) :: lanczos_subd! sub-diagonal part of niter-dimensional lanczos tridiagonal matrix
   real(real_kind) :: alphaold, betaold, gamma, eps, ortho
   real(real_kind), dimension(niter, niter) :: lanczos_orth! sub-diagonal part of niter-dimensional lanczos tridiagonal matrix
   real(real_kind), dimension(2*niter-2) :: work
   real(real_kind), allocatable, dimension(:,:,:) :: h_eig! the matrix of orthonormal vectors that approximate the hessian eigenvectors
   integer(int_kind), allocatable, dimension(:) :: eig, eig1
   real(real_kind), allocatable, dimension(:) :: eigenval
   logical :: eig_loc(niter)
   real(real_kind) :: criteria ! convergence criterion
   real(real_kind) :: cgtol
!=============================================================================
!! convergence criteria
!
   if (ioutl.eq.1) cgtol = cgtol1
   if (ioutl.gt.1) then
     if (ioutl.eq.2) cgtol = cgtol2
   endif
   if (par%ismasterproc) print*,'Reorthogonalized cg'
   if (par%ismasterproc) print*,'Convergence criterion:',cgtol
!==========================================================================
! A. Check Maximum ABS Errors of Hx - Obs for V,U,T,Q,PS
!==========================================================================
#define Check_HxmObs_on
#ifdef Check_HxmObs_on
   if (par%ismasterproc) write(iulog,*) 'Check hx_bg-obs'
   ierr = timingstart('Check_HxmObs') !---profiling
   do itime = 1,ntime
     if (obcheck) then
       call rej_obs_omb(itime) ! rejection of obs. based on omb check(hjs)
     endif
     call check_hxmobs(itime) ! extended to 4d
   enddo
   if (.not.obcheck_outer) then
     obcheck = .false.
   endif
   ierr = timingstop('Check_HxmObs')
#endif
   tmpmax = maxval(gs_md)
   tmpmin = minval(gs_md)
   call par_reduce(tmpmax,maxv,1,par_max)
   call par_reduce(tmpmin,minv,1,par_min)
   if (par%ismasterproc) write(iulog,*) 'Min max [gs_md]:',minv,maxv
!==========================================================================
! RTTOV,ROPP calculation for outer-loop
!==========================================================================
   if (ioutl.gt.1) then
!! RTTOV used
     if (rttovused) then
       do itime = 1, ntime
         call calrttov(itime,gs_md)
       enddo
     endif
!! ROPP used
     if (roppused) then
       do itime = 1, ntime
         call calropp(itime,gs_md)
       enddo
     endif
   endif
!==========================================================================
! B. Make b
!==========================================================================
   if (par%ismasterproc) write(iulog,*) 'Make_RHS_b'
   ierr = timingstart('Make_RHS_b') !---profiling
   x_eig_acc = x_eig_acc+x_eig
   x_eig_acc_chi = x_eig_acc_chi+x_eig_chi
   x_eig_acc_tps = x_eig_acc_tps+x_eig_tps
   x_eig_acc_q = x_eig_acc_q+x_eig_q
   x_eig_acc_a = x_eig_acc_a+x_eig_a
   x_eig = 0.d0
   x_eig_chi = 0.d0
   x_eig_tps = 0.d0
   x_eig_q = 0.d0
   x_eig_a = 0.d0
   call make_rhs_b(gs_md,b_eig,b_eig_chi,b_eig_tps,b_eig_q,b_eig_a ! modified for alpha-variable
  )
   ierr = timingstop('Make_RHS_b')
!==========================================================================
!==========================================================================
! C. Evaluation of adjoint operators
!==========================================================================
#define Check_Adj_on
#ifdef Check_Adj_on
   if (par%ismasterproc) write(iulog,*) 'Check_Adj'
   ierr = timingstart('Check_Adj') !---profiling
   call check_adj(gs_md,b_eig,b_eig_chi,b_eig_tps,b_eig_q,b_eig_a)
   ierr = timingstop('Check_Adj')
#endif
!==========================================================================
!==========================================================================
! D. Make Ax
!==========================================================================
   ierr = timingstart('Make_LHS_Ax0') !---profiling
   call make_lhs_ax(gs_md,x_eig,x_eig_chi,x_eig_tps,x_eig_q,x_eig_a,ap_eig,ap_eig_chi,ap_eig_tps,ap_eig_q,ap_eig_a)
   ierr = timingstop('Make_LHS_Ax0')
!==========================================================================
!==========================================================================
! E. Initialize r,p,and rsold
!==========================================================================
! 1. Initialize r and p
   do wn = 1,zwn_nproc
     r_eig(:,wn) = b_eig(:,wn)-ap_eig(:,wn)
     r_eig_chi(:,wn) = b_eig_chi(:,wn)-ap_eig_chi(:,wn)
     r_eig_tps(:,wn) = b_eig_tps(:,wn)-ap_eig_tps(:,wn)
     r_eig_q(:,wn) = b_eig_q(:,wn)-ap_eig_q(:,wn)
     p_eig(:,wn) = r_eig(:,wn)
     p_eig_chi(:,wn) = r_eig_chi(:,wn)
     p_eig_tps(:,wn) = r_eig_tps(:,wn)
     p_eig_q(:,wn) = r_eig_q(:,wn)
   enddo
   do ismpl = 1,nsmpl !--for alpha-variable
     do wn = 1, zwn_nproc_a
       r_eig_a(:,wn,ismpl) = b_eig_a(:,wn,ismpl)-ap_eig_a(:,wn,ismpl)
       p_eig_a(:,wn,ismpl) = r_eig_a(:,wn,ismpl)
     enddo
   enddo ! ismpl
! 2. Initialize rsold
   ierr = timingstart('InnerProduct_r0') !---profiling
   call innerproduct(neig,neig_chi,neig_tps,neig_q,zwn_nproc,zwn_nproc_a,r_eig,r_eig_chi,r_eig_tps,r_eig_q,r_eig_a,! modified for alpha-variable
   r_eig,r_eig_chi,r_eig_tps,r_eig_q,r_eig_a,rsold)
   ierr = timingstop('InnerProduct_r0') !---profiling
   if (par%ismasterproc) write(iulog,*) 'rsold,',rsold
   rs0 = rsold
! 3. Initialize r*d_eig_history : Save the residual history
   r_eig_history(:,:,:) = 0.d0
   r_eig_chi_history(:,:,:) = 0.d0
   r_eig_tps_history(:,:,:) = 0.d0
   r_eig_q_history(:,:,:) = 0.d0
   r_eig_q_history_a(:,:,:,:) = 0.d0
   r_eig_history(:,:,1) = r_eig/dsqrt(rsold)
   r_eig_chi_history(:,:,1) = r_eig_chi/dsqrt(rsold)
   r_eig_tps_history(:,:,1) = r_eig_tps/dsqrt(rsold)
   r_eig_q_history(:,:,1) = r_eig_q/dsqrt(rsold)
   r_eig_q_history_a(:,:,:,1) = r_eig_a/dsqrt(rsold)
!==========================================================================
!==========================================================================
! F. Main iteration in conjugate gradient
!==========================================================================
   ierr = timingstart('Innerloop') !---profiling
   do iter = 1,niter
!       IF (par%ismasterproc) WRITE(iulog,*) 'x_eig(:,1)'
!       IF (par%ismasterproc) WRITE(iulog,*)  x_eig(:,1)
     tmpmax = maxval(x_eig)
     tmpmin = minval(x_eig)
     call par_reduce(tmpmax,maxv,1,par_max)
     call par_reduce(tmpmin,minv,1,par_min)
     if (par%ismasterproc) write(iulog,*) 'Min max [x_eig]:',minv,maxv
     tmpmax = maxval(x_eig_chi)
     tmpmin = minval(x_eig_chi)
     call par_reduce(tmpmax,maxv,1,par_max)
     call par_reduce(tmpmin,minv,1,par_min)
     if (par%ismasterproc) write(iulog,*) 'Min max [x_eig_chi]:',minv,maxv
     tmpmax = maxval(x_eig_tps)
     tmpmin = minval(x_eig_tps)
     call par_reduce(tmpmax,maxv,1,par_max)
     call par_reduce(tmpmin,minv,1,par_min)
     if (par%ismasterproc) write(iulog,*) 'Min max [x_eig_tps]:',minv,maxv
     tmpmax = maxval(x_eig_q)
     tmpmin = minval(x_eig_q)
     call par_reduce(tmpmax,maxv,1,par_max)
     call par_reduce(tmpmin,minv,1,par_min)
     if (par%ismasterproc) write(iulog,*) 'Min max [x_eig_q]:',minv,maxv
     tmpmax = maxval(x_eig_a)
     tmpmin = minval(x_eig_a)
     call par_reduce(tmpmax,maxv,1,par_max)
     call par_reduce(tmpmin,minv,1,par_min)
     if (par%ismasterproc) write(iulog,*) 'Min max [x_eig_a]:',minv,maxv
!        IF( iter == 1 .or. mod(iter,10).EQ. 0)THEN
! 1. Calculate Cost Function
     ierr = timingstart('Cal_CostFun') !---profiling
     call cal_cost_func(gs_md,x_eig,x_eig_chi,x_eig_tps,x_eig_q,x_eig_a,j_o,j_b)
     ierr = timingstop('Cal_CostFun')
!        ENDIF
! 2. Calculate Ap ------------------------------------------------------
     ierr = timingstart('Make_LHS_Ax') !---profiling
     call make_lhs_ax(gs_md,p_eig,p_eig_chi,p_eig_tps,p_eig_q,p_eig_a,ap_eig,ap_eig_chi,ap_eig_tps,ap_eig_q,ap_eig_a)
     ierr = timingstop('Make_LHS_Ax')
!       IF (par%ismasterproc) WRITE(iulog,*) 'Ap_eig(:,1)'
!       IF (par%ismasterproc) WRITE(iulog,*)  Ap_eig(:,1)
!       tmpmax= maxval(Ap_eig)
!       tmpmin= minval(Ap_eig)
!       CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
!       CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
!       IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [Ap_eig]:',minv,maxv
! 3-1. Calculate alpha,x,r,and rsnew
     ierr = timingstart('InnerProduct_pAp') !---profiling
     call innerproduct(neig,neig_chi,neig_tps,neig_q,zwn_nproc,zwn_nproc_a,p_eig,p_eig_chi,p_eig_tps,p_eig_q,p_eig_a,ap_eig,ap_eig_chi,ap_eig_tps,ap_eig_q,ap_eig_a,pap)
     ierr = timingstop('InnerProduct_pAp') !---profiling
!       IF (par%ismasterproc) WRITE(iulog,*) 'pAp,',pAp
!       tmp = 0.D0
!       DO wn = 1,zwn_nproc
!          IF( wn+zwn_iproc .LE. zwn2 )THEN
!             tmp = tmp + SUM(p_eig(:,wn)*Ap_eig(:,wn))
!          END IF
!       END DO
!       global_shared_buf(nets:nete,1) = tmp/DBLE(nete-nets+1)
!       CALL Wrap_Repro_Sum(nvars=1,comm=hybrid%par%comm)
!       pAp = global_shared_sum(1)
!       IF (par%ismasterproc) WRITE(iulog,*) 'pAp,',pAp
     alpha = rsold/pap
! 3-2. Update the iterate and the residual
     do wn = 1,zwn_nproc
       x_eig(:,wn) = x_eig(:,wn)+alpha*p_eig(:,wn)
       x_eig_chi(:,wn) = x_eig_chi(:,wn)+alpha*p_eig_chi(:,wn)
       x_eig_tps(:,wn) = x_eig_tps(:,wn)+alpha*p_eig_tps(:,wn)
       x_eig_q(:,wn) = x_eig_q(:,wn)+alpha*p_eig_q(:,wn)
       r_eig(:,wn) = r_eig(:,wn)-alpha*ap_eig(:,wn)
       r_eig_chi(:,wn) = r_eig_chi(:,wn)-alpha*ap_eig_chi(:,wn)
       r_eig_tps(:,wn) = r_eig_tps(:,wn)-alpha*ap_eig_tps(:,wn)
       r_eig_q(:,wn) = r_eig_q(:,wn)-alpha*ap_eig_q(:,wn)
     enddo
     do ismpl = 1,nsmpl
       do wn = 1, zwn_nproc_a
         x_eig_a(:,wn,ismpl) = x_eig_a(:,wn,ismpl)+alpha*p_eig_a(:,wn,ismpl)
         r_eig_a(:,wn,ismpl) = r_eig_a(:,wn,ismpl)-alpha*ap_eig_a(:,wn,ismpl)
       enddo
     enddo ! ismpl
!       IF (par%ismasterproc) WRITE(iulog,*) 'p_eig(:,1)'
!       IF (par%ismasterproc) WRITE(iulog,*)  p_eig(:,1)
!       tmpmax= maxval(p_eig)
!       tmpmin= minval(p_eig)
!       CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
!       CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
!       IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [p_eig]:',minv,maxv
!       IF (par%ismasterproc) WRITE(iulog,*) 'x_eig(:,1)'
!       IF (par%ismasterproc) WRITE(iulog,*)  x_eig(:,1)
!       tmpmax= maxval(x_eig)
!       tmpmin= minval(x_eig)
!       CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
!       CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
!       IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [x_eig]:',minv,maxv
!       IF (par%ismasterproc) WRITE(iulog,*) 'Ap_eig(:,1)'
!       IF (par%ismasterproc) WRITE(iulog,*)  Ap_eig(:,1)
!       tmpmax= maxval(Ap_eig)
!       tmpmin= minval(Ap_eig)
!       CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
!       CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
!       IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [Ap_eig]:',minv,maxv
! 3-3. Re-orthogonalization step 1
! Remove component in previous residual direction
     do k = 1,iter
! calculate the dot_product of r_k and r_j
! r_k : previous residual
! r_j : current  residual
       ierr = timingstart('InnerProduct_rirj') !---profiling
       call innerproduct(neig,neig_chi,neig_tps,neig_q,zwn_nproc,zwn_nproc_a,r_eig,r_eig_chi,r_eig_tps,r_eig_q,r_eig_a,r_eig_history(:,:,k),r_eig_chi_history(:,:,k),r_eig_tps_history(:,:,k),r_eig_q_history(:,:,k),r_eig_q_history_a(:,:,:,k),rjrk)
       ierr = timingstop('InnerProduct_rirj') !---profiling
       do wn = 1,zwn_nproc
         r_eig(:,wn) = r_eig(:,wn)-rjrk*r_eig_history(:,wn,k)
         r_eig_chi(:,wn) = r_eig_chi(:,wn)-rjrk*r_eig_chi_history(:,wn,k)
         r_eig_tps(:,wn) = r_eig_tps(:,wn)-rjrk*r_eig_tps_history(:,wn,k)
         r_eig_q(:,wn) = r_eig_q(:,wn)-rjrk*r_eig_q_history(:,wn,k)
       enddo
       do ismpl = 1,nsmpl
         do wn = 1, zwn_nproc_a
           r_eig_a(:,wn,ismpl) = r_eig_a(:,wn,ismpl)-rjrk*r_eig_q_history_a(:,wn,ismpl,k)
         enddo
       enddo ! ismpl
     enddo
     ierr = timingstart('InnerProduct_rsnew') !---profiling
     call innerproduct(neig,neig_chi,neig_tps,neig_q,zwn_nproc,zwn_nproc_a,r_eig,r_eig_chi,r_eig_tps,r_eig_q,r_eig_a,r_eig,r_eig_chi,r_eig_tps,r_eig_q,r_eig_a,rsnew)
     ierr = timingstop('InnerProduct_rsnew') !---profiling
     beta = rsnew/rsold
     criteria = rsnew/rs0
     rsold = rsnew
!             IF( iter == 1 )  rs0 = rsnew
! 4. Check rsnew
     if (par%ismasterproc) write(iulog,*) 'Check residual'
     if (par%ismasterproc) write(iulog,*) 'ioutl = ',ioutl,'iter = ',iter,'rsnew = ',rsnew,criteria
     if (par%ismasterproc) write(21,*) ioutl,iter,criteria,j_b+j_o,j_b,j_o
! Convergence criterion of criteri as normalized rs
     if (criteria.lt.cgtol) then
       if (par%ismasterproc) write(iulog,*) 'ioutl = ',ioutl,'total iteration number = ',iter
       exit
     endif
! 5. Calculate p and rsold
     do wn = 1,zwn_nproc
       p_eig(:,wn) = r_eig(:,wn)+beta*p_eig(:,wn)
       p_eig_chi(:,wn) = r_eig_chi(:,wn)+beta*p_eig_chi(:,wn)
       p_eig_tps(:,wn) = r_eig_tps(:,wn)+beta*p_eig_tps(:,wn)
       p_eig_q(:,wn) = r_eig_q(:,wn)+beta*p_eig_q(:,wn)
     enddo
     do ismpl = 1,nsmpl
       do wn = 1, zwn_nproc_a
         p_eig_a(:,wn,ismpl) = r_eig_a(:,wn,ismpl)+beta*p_eig_a(:,wn,ismpl)
       enddo
     enddo ! ismpl
! 6. Re-orthogonalization step 2 :
! normalize residual vector and save it
     if (iter<niter) then
       r_eig_history(:,:,iter+1) = r_eig/dsqrt(rsnew)
       r_eig_chi_history(:,:,iter+1) = r_eig_chi/dsqrt(rsnew)
       r_eig_tps_history(:,:,iter+1) = r_eig_tps/dsqrt(rsnew)
       r_eig_q_history(:,:,iter+1) = r_eig_q/dsqrt(rsnew)
       r_eig_q_history_a(:,:,:,iter+1) = r_eig_a/dsqrt(rsnew)
     endif
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
   enddo ! for iter
   ierr = timingstop('Innerloop') !---profiling
!
   end subroutine rcgmethod
!===============================================================================
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine check_hxmobs(itime)
!
   integer(int_kind), intent(in   ) :: itime
! 2-2. Local
   real(real_kind), dimension(17) :: tmpmax, tmpmin
   real(real_kind) :: maxsdu, maxsdv, maxsdt, maxsdq, maxsfu, maxsfv, maxsft, maxsfq, maxsfp, minsdu, minsdv, minsdt, minsdq, minsfu, minsfv, minsft, minsfq, minsfp, maxaru, maxarv, maxart, minaru, minarv, minart, maxamvu, maxamvv, minamvu, minamvv, maxscatwindu, maxscatwindv, minscatwindu, minscatwindv, maxbgp, minbgp
   integer(int_kind) :: i, j, k, ie
! For verification -------------------------------------------------------------- ihKwon
!
   call modeltoobs(itime,bg_md(:,:,:,:,itime),gs_obs(itime))
   if (par%ismasterproc) then
     print*,'Variable,obs,background,o-b,observation error'
     print*,'itime = ',itime
   endif
   tmpmax = 0.
   tmpmin = 0.
! single observation test
   if (singleobs) then
     do ie = nets, nete
       do i = 1, obsmet(itime)%nsg_e(ie)
         if ((gs_obs(itime)%single(i, ie)%u%val.ne.sdmiss).and.(obs(itime)%single(i, ie)%u%val.ne.sdmiss)) then
           tmpmax(1) = max(obs(itime)%single(i,ie)%u%val-gs_obs(itime)%single(i,ie)%u%val,tmpmax(1))
           tmpmin(1) = min(obs(itime)%single(i,ie)%u%val-gs_obs(itime)%single(i,ie)%u%val,tmpmin(1))
!      IF ( par%ismasterproc .and. i == 1 ) THEN
           print*,'Single u:',obs(itime)%single(i,ie)%u%val,gs_obs(itime)%single(i,ie)%u%val,obs(itime)%single(i,ie)%u%val-gs_obs(itime)%single(i,ie)%u%val,sqrt(obs(itime)%single(i,ie)%u%error)
!      END IF
         endif
       enddo
       do i = 1, obsmet(itime)%nsg_e(ie)
         if ((gs_obs(itime)%single(i, ie)%v%val.ne.sdmiss).and.(obs(itime)%single(i, ie)%v%val.ne.sdmiss)) then
           tmpmax(2) = max(obs(itime)%single(i,ie)%v%val-gs_obs(itime)%single(i,ie)%v%val,tmpmax(1))
           tmpmin(2) = min(obs(itime)%single(i,ie)%v%val-gs_obs(itime)%single(i,ie)%v%val,tmpmin(1))
!      IF ( par%ismasterproc .and. i == 1 ) THEN
           print*,'Single v:',obs(itime)%single(i,ie)%v%val,gs_obs(itime)%single(i,ie)%v%val,obs(itime)%single(i,ie)%v%val-gs_obs(itime)%single(i,ie)%v%val,sqrt(obs(itime)%single(i,ie)%v%error)
!      END IF
         endif
       enddo
       do i = 1, obsmet(itime)%nsg_e(ie)
         if ((gs_obs(itime)%single(i, ie)%t%val.ne.sdmiss).and.(obs(itime)%single(i, ie)%t%val.ne.sdmiss)) then
           tmpmax(3) = max(obs(itime)%single(i,ie)%t%val-gs_obs(itime)%single(i,ie)%t%val,tmpmax(1))
           tmpmin(4) = min(obs(itime)%single(i,ie)%t%val-gs_obs(itime)%single(i,ie)%t%val,tmpmin(1))
!      IF ( par%ismasterproc .and. i == 1 ) THEN
           print*,'Single t:',obs(itime)%single(i,ie)%t%val,gs_obs(itime)%single(i,ie)%t%val,obs(itime)%single(i,ie)%t%val-gs_obs(itime)%single(i,ie)%t%val,sqrt(obs(itime)%single(i,ie)%t%error)
!      END IF
         endif
       enddo
       do i = 1, obsmet(itime)%nsg_e(ie)
         if ((gs_obs(itime)%single(i, ie)%q%val.ne.sdmiss).and.(obs(itime)%single(i, ie)%q%val.ne.sdmiss)) then
           tmpmax(4) = max(obs(itime)%single(i,ie)%q%val-gs_obs(itime)%single(i,ie)%q%val,tmpmax(1))
           tmpmin(4) = min(obs(itime)%single(i,ie)%q%val-gs_obs(itime)%single(i,ie)%q%val,tmpmin(1))
!      IF ( par%ismasterproc .and. i == 1 ) THEN
           print*,'Single q:',obs(itime)%single(i,ie)%q%val,gs_obs(itime)%single(i,ie)%q%val,obs(itime)%single(i,ie)%q%val-gs_obs(itime)%single(i,ie)%q%val,sqrt(obs(itime)%single(i,ie)%q%error)
!      END IF
         endif
       enddo
       do i = 1, obsmet(itime)%nsg_e(ie)
         if ((gs_obs(itime)%single(i, ie)%ps%val.ne.sdmiss).and.(obs(itime)%single(i, ie)%ps%val.ne.sdmiss)) then
           tmpmax(5) = max(obs(itime)%single(i,ie)%ps%val-gs_obs(itime)%single(i,ie)%ps%val,tmpmax(1))
           tmpmin(5) = min(obs(itime)%single(i,ie)%ps%val-gs_obs(itime)%single(i,ie)%ps%val,tmpmin(1))
!      IF ( par%ismasterproc .and. i == 1 ) THEN
           print*,'Single ps:',obs(itime)%single(i,ie)%ps%val,gs_obs(itime)%single(i,ie)%ps%val,obs(itime)%single(i,ie)%ps%val-gs_obs(itime)%single(i,ie)%ps%val,sqrt(obs(itime)%single(i,ie)%ps%error)
!      END IF
         endif
       enddo
     enddo
   endif
! sonde observation
   if (da_sonde) then
     do ie = nets, nete
       do i = 1, obsmet(itime)%nsd_e(1, ie)
         if ((gs_obs(itime)%sonde(i, ie)%u%val.ne.sdmiss).and.(obs(itime)%sonde(i, ie)%u%val.ne.sdmiss)) then
           tmpmax(1) = max(obs(itime)%sonde(i,ie)%u%val-gs_obs(itime)%sonde(i,ie)%u%val,tmpmax(1))
           tmpmin(1) = min(obs(itime)%sonde(i,ie)%u%val-gs_obs(itime)%sonde(i,ie)%u%val,tmpmin(1))
           if (par%ismasterproc.and.i==1) then
             print*,'Sonde u:',obs(itime)%sonde(i,ie)%u%val,gs_obs(itime)%sonde(i,ie)%u%val,obs(itime)%sonde(i,ie)%u%val-gs_obs(itime)%sonde(i,ie)%u%val,sqrt(obs(itime)%sonde(i,ie)%u%error)
           endif
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(2, ie)
         if ((gs_obs(itime)%sonde(i, ie)%v%val.ne.sdmiss).and.(obs(itime)%sonde(i, ie)%v%val.ne.sdmiss)) then
           tmpmax(2) = max(obs(itime)%sonde(i,ie)%v%val-gs_obs(itime)%sonde(i,ie)%v%val,tmpmax(2))
           tmpmin(2) = min(obs(itime)%sonde(i,ie)%v%val-gs_obs(itime)%sonde(i,ie)%v%val,tmpmin(2))
           if (par%ismasterproc.and.i==1) then
             print*,'Sonde v:',obs(itime)%sonde(i,ie)%v%val,gs_obs(itime)%sonde(i,ie)%v%val,obs(itime)%sonde(i,ie)%v%val-gs_obs(itime)%sonde(i,ie)%v%val,sqrt(obs(itime)%sonde(i,ie)%v%error)
           endif
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(3, ie)
         if ((gs_obs(itime)%sonde(i, ie)%t%val.ne.sdmiss).and.(obs(itime)%sonde(i, ie)%t%val.ne.sdmiss)) then
           tmpmax(3) = max(obs(itime)%sonde(i,ie)%t%val-gs_obs(itime)%sonde(i,ie)%t%val,tmpmax(3))
           tmpmin(3) = min(obs(itime)%sonde(i,ie)%t%val-gs_obs(itime)%sonde(i,ie)%t%val,tmpmin(3))
           if (par%ismasterproc.and.i==1) then
             print*,'Sonde t:',obs(itime)%sonde(i,ie)%t%val,gs_obs(itime)%sonde(i,ie)%t%val,obs(itime)%sonde(i,ie)%t%val-gs_obs(itime)%sonde(i,ie)%t%val,sqrt(obs(itime)%sonde(i,ie)%t%error)
           endif
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(4, ie)
         if ((gs_obs(itime)%sonde(i, ie)%q%val.ne.sdmiss).and.(obs(itime)%sonde(i, ie)%q%val.ne.sdmiss)) then
           tmpmax(4) = max(obs(itime)%sonde(i,ie)%q%val-gs_obs(itime)%sonde(i,ie)%q%val,tmpmax(4))
           tmpmin(4) = min(obs(itime)%sonde(i,ie)%q%val-gs_obs(itime)%sonde(i,ie)%q%val,tmpmin(4))
           if (par%ismasterproc.and.i==1) then
             print*,'Sonde q:',obs(itime)%sonde(i,ie)%q%val,gs_obs(itime)%sonde(i,ie)%q%val,obs(itime)%sonde(i,ie)%q%val-gs_obs(itime)%sonde(i,ie)%q%val,sqrt(obs(itime)%sonde(i,ie)%q%error)
           endif
         endif
       enddo
     enddo
   endif
! surface observation
   if (da_surface) then
     do ie = nets, nete
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((gs_obs(itime)%surface(i, ie)%u%val.ne.sdmiss).and.(obs(itime)%surface(i, ie)%u%val.ne.sdmiss)) then
           tmpmax(5) = max(obs(itime)%surface(i,ie)%u%val-gs_obs(itime)%surface(i,ie)%u%val,tmpmax(5))
           tmpmin(5) = min(obs(itime)%surface(i,ie)%u%val-gs_obs(itime)%surface(i,ie)%u%val,tmpmin(5))
           if (par%ismasterproc.and.i==1) then
             print*,'Surface u:',obs(itime)%surface(i,ie)%u%val,gs_obs(itime)%surface(i,ie)%u%val,obs(itime)%surface(i,ie)%u%val-gs_obs(itime)%surface(i,ie)%u%val,sqrt(obs(itime)%surface(i,ie)%u%error)
           endif
         endif
         if ((gs_obs(itime)%surface(i, ie)%v%val.ne.sdmiss).and.(obs(itime)%surface(i, ie)%v%val.ne.sdmiss)) then
           tmpmax(6) = max(obs(itime)%surface(i,ie)%v%val-gs_obs(itime)%surface(i,ie)%v%val,tmpmax(6))
           tmpmin(6) = min(obs(itime)%surface(i,ie)%v%val-gs_obs(itime)%surface(i,ie)%v%val,tmpmin(6))
           if (par%ismasterproc.and.i==1) then
             print*,'Surface v:',obs(itime)%surface(i,ie)%v%val,gs_obs(itime)%surface(i,ie)%v%val,obs(itime)%surface(i,ie)%v%val-gs_obs(itime)%surface(i,ie)%v%val,sqrt(obs(itime)%surface(i,ie)%v%error)
           endif
         endif
         if ((gs_obs(itime)%surface(i, ie)%t%val.ne.sdmiss).and.(obs(itime)%surface(i, ie)%t%val.ne.sdmiss)) then
           tmpmax(7) = max(obs(itime)%surface(i,ie)%t%val-gs_obs(itime)%surface(i,ie)%t%val,tmpmax(7))
           tmpmin(7) = min(obs(itime)%surface(i,ie)%t%val-gs_obs(itime)%surface(i,ie)%t%val,tmpmin(7))
           if (par%ismasterproc.and.i==1) then
             print*,'Surface t:',obs(itime)%surface(i,ie)%t%val,gs_obs(itime)%surface(i,ie)%t%val,obs(itime)%surface(i,ie)%t%val-gs_obs(itime)%surface(i,ie)%t%val,sqrt(obs(itime)%surface(i,ie)%t%error)
           endif
         endif
         if ((gs_obs(itime)%surface(i, ie)%q%val.ne.sdmiss).and.(obs(itime)%surface(i, ie)%q%val.ne.sdmiss)) then
           tmpmax(8) = max(obs(itime)%surface(i,ie)%q%val-gs_obs(itime)%surface(i,ie)%q%val,tmpmax(8))
           tmpmin(8) = min(obs(itime)%surface(i,ie)%q%val-gs_obs(itime)%surface(i,ie)%q%val,tmpmin(8))
           if (par%ismasterproc.and.i==1) then
             print*,'Surface q:',obs(itime)%surface(i,ie)%q%val,gs_obs(itime)%surface(i,ie)%q%val,obs(itime)%surface(i,ie)%q%val-gs_obs(itime)%surface(i,ie)%q%val,sqrt(obs(itime)%surface(i,ie)%q%error)
           endif
         endif
         if ((gs_obs(itime)%surface(i, ie)%p%val.ne.sdmiss).and.(obs(itime)%surface(i, ie)%p%val.ne.sdmiss)) then
           tmpmax(9) = max(obs(itime)%surface(i,ie)%p%val-gs_obs(itime)%surface(i,ie)%p%val,tmpmax(9))
           tmpmin(9) = min(obs(itime)%surface(i,ie)%p%val-gs_obs(itime)%surface(i,ie)%p%val,tmpmin(9))
           if (par%ismasterproc.and.i==1) then
             print*,'Surface ps:',obs(itime)%surface(i,ie)%p%val,gs_obs(itime)%surface(i,ie)%p%val,obs(itime)%surface(i,ie)%p%val-gs_obs(itime)%surface(i,ie)%p%val,sqrt(obs(itime)%surface(i,ie)%p%error)
           endif
         endif
       enddo
     enddo
   endif
! bogus observation
   if (da_bogus) then
     do ie = nets, nete
       do i = 1, obsmet(itime)%nbg_e(ie)
         if ((gs_obs(itime)%bogus(i, ie)%p%val.ne.sdmiss).and.(obs(itime)%bogus(i, ie)%p%val.ne.sdmiss)) then
           tmpmax(10) = max(obs(itime)%bogus(i,ie)%p%val-gs_obs(itime)%bogus(i,ie)%p%val,tmpmax(10))
           tmpmin(10) = min(obs(itime)%bogus(i,ie)%p%val-gs_obs(itime)%bogus(i,ie)%p%val,tmpmin(10))
           if (par%ismasterproc.and.i==1) then
             print*,'Surface_bg ps:',obs(itime)%bogus(i,ie)%p%val,gs_obs(itime)%bogus(i,ie)%p%val,obs(itime)%bogus(i,ie)%p%val-gs_obs(itime)%bogus(i,ie)%p%val,sqrt(obs(itime)%bogus(i,ie)%p%error)
           endif
         endif
       enddo
     enddo
   endif
! aircraft observation
   if (da_aircraft) then
     do ie = nets, nete
       do i = 1, obsmet(itime)%nar_e(1, ie)
         if ((gs_obs(itime)%aircraft(i, ie)%u%val.ne.sdmiss).and.(obs(itime)%aircraft(i, ie)%u%val.ne.sdmiss)) then
           tmpmax(11) = max(obs(itime)%aircraft(i,ie)%u%val-gs_obs(itime)%aircraft(i,ie)%u%val,tmpmax(11))
           tmpmin(11) = min(obs(itime)%aircraft(i,ie)%u%val-gs_obs(itime)%aircraft(i,ie)%u%val,tmpmin(11))
           if (par%ismasterproc.and.i==1) then
             print*,'Aircraft u:',obs(itime)%aircraft(i,ie)%u%val,gs_obs(itime)%aircraft(i,ie)%u%val,obs(itime)%aircraft(i,ie)%u%val-gs_obs(itime)%aircraft(i,ie)%u%val,sqrt(obs(itime)%aircraft(i,ie)%u%error)
           endif
         endif
       enddo
       do i = 1, obsmet(itime)%nar_e(2, ie)
         if ((gs_obs(itime)%aircraft(i, ie)%v%val.ne.sdmiss).and.(obs(itime)%aircraft(i, ie)%v%val.ne.sdmiss)) then
           tmpmax(12) = max(obs(itime)%aircraft(i,ie)%v%val-gs_obs(itime)%aircraft(i,ie)%v%val,tmpmax(12))
           tmpmin(12) = min(obs(itime)%aircraft(i,ie)%v%val-gs_obs(itime)%aircraft(i,ie)%v%val,tmpmin(12))
           if (par%ismasterproc.and.i==1) then
             print*,'Aircraft v:',obs(itime)%aircraft(i,ie)%v%val,gs_obs(itime)%aircraft(i,ie)%v%val,obs(itime)%aircraft(i,ie)%v%val-gs_obs(itime)%aircraft(i,ie)%v%val,sqrt(obs(itime)%aircraft(i,ie)%v%error)
           endif
         endif
       enddo
       do i = 1, obsmet(itime)%nar_e(3, ie)
         if ((gs_obs(itime)%aircraft(i, ie)%t%val.ne.sdmiss).and.(obs(itime)%aircraft(i, ie)%t%val.ne.sdmiss)) then
           tmpmax(13) = max(obs(itime)%aircraft(i,ie)%t%val-gs_obs(itime)%aircraft(i,ie)%t%val,tmpmax(13))
           tmpmin(13) = min(obs(itime)%aircraft(i,ie)%t%val-gs_obs(itime)%aircraft(i,ie)%t%val,tmpmin(13))
           if (par%ismasterproc.and.i==1) then
             print*,'Aircraft t:',obs(itime)%aircraft(i,ie)%t%val,gs_obs(itime)%aircraft(i,ie)%t%val,obs(itime)%aircraft(i,ie)%t%val-gs_obs(itime)%aircraft(i,ie)%t%val,sqrt(obs(itime)%aircraft(i,ie)%t%error)
           endif
         endif
       enddo
     enddo
   endif
! amv observation
   if (da_amv) then
     do ie = nets, nete
       do i = 1, obsmet(itime)%namv_e(ie)
         if ((gs_obs(itime)%amv(i, ie)%u%val.ne.sdmiss).and.(obs(itime)%amv(i, ie)%u%val.ne.sdmiss)) then
           tmpmax(14) = max(obs(itime)%amv(i,ie)%u%val-gs_obs(itime)%amv(i,ie)%u%val,tmpmax(14))
           tmpmin(14) = min(obs(itime)%amv(i,ie)%u%val-gs_obs(itime)%amv(i,ie)%u%val,tmpmin(14))
           if (par%ismasterproc.and.i==1) then
             print*,'Amv u:',obs(itime)%amv(i,ie)%u%val,gs_obs(itime)%amv(i,ie)%u%val,obs(itime)%amv(i,ie)%u%val-gs_obs(itime)%amv(i,ie)%u%val,sqrt(obs(itime)%amv(i,ie)%u%error)
           endif
         endif
         if ((gs_obs(itime)%amv(i, ie)%v%val.ne.sdmiss).and.(obs(itime)%amv(i, ie)%v%val.ne.sdmiss)) then
           tmpmax(15) = max(obs(itime)%amv(i,ie)%v%val-gs_obs(itime)%amv(i,ie)%v%val,tmpmax(15))
           tmpmin(15) = min(obs(itime)%amv(i,ie)%v%val-gs_obs(itime)%amv(i,ie)%v%val,tmpmin(15))
           if (par%ismasterproc.and.i==1) then
             print*,'Amv v:',obs(itime)%amv(i,ie)%v%val,gs_obs(itime)%amv(i,ie)%v%val,obs(itime)%amv(i,ie)%v%val-gs_obs(itime)%amv(i,ie)%v%val,sqrt(obs(itime)%amv(i,ie)%v%error)
           endif
         endif
       enddo
     enddo
   endif
! scatwind observation
   if (da_scatwind) then
     do ie = nets, nete
       do i = 1, obsmet(itime)%nscatwind_e(ie)
         if ((gs_obs(itime)%scatwind(i, ie)%u%val.ne.sdmiss).and.(obs(itime)%scatwind(i, ie)%u%val.ne.sdmiss)) then
           tmpmax(16) = max(obs(itime)%scatwind(i,ie)%u%val-gs_obs(itime)%scatwind(i,ie)%u%val,tmpmax(16))
           tmpmin(16) = min(obs(itime)%scatwind(i,ie)%u%val-gs_obs(itime)%scatwind(i,ie)%u%val,tmpmin(16))
           if (par%ismasterproc.and.i==1) then
             print*,'Scatwind u:',obs(itime)%scatwind(i,ie)%u%val,gs_obs(itime)%scatwind(i,ie)%u%val,obs(itime)%scatwind(i,ie)%u%val-gs_obs(itime)%scatwind(i,ie)%u%val,sqrt(obs(itime)%scatwind(i,ie)%u%error)
           endif
         endif
         if ((gs_obs(itime)%scatwind(i, ie)%v%val.ne.sdmiss).and.(obs(itime)%scatwind(i, ie)%v%val.ne.sdmiss)) then
           tmpmax(17) = max(obs(itime)%scatwind(i,ie)%v%val-gs_obs(itime)%scatwind(i,ie)%v%val,tmpmax(17))
           tmpmin(17) = min(obs(itime)%scatwind(i,ie)%v%val-gs_obs(itime)%scatwind(i,ie)%v%val,tmpmin(17))
           if (par%ismasterproc.and.i==1) then
             print*,'Scatwind v:',obs(itime)%scatwind(i,ie)%v%val,gs_obs(itime)%scatwind(i,ie)%v%val,obs(itime)%scatwind(i,ie)%v%val-gs_obs(itime)%scatwind(i,ie)%v%val,sqrt(obs(itime)%scatwind(i,ie)%v%error)
           endif
         endif
       enddo
     enddo
   endif
   if (singleobs) then
     call par_reduce(tmpmax(1),maxsdu,1,par_max)
     call par_reduce(tmpmax(2),maxsdv,1,par_max)
     call par_reduce(tmpmax(3),maxsdt,1,par_max)
     call par_reduce(tmpmax(4),maxsdq,1,par_max)
     call par_reduce(tmpmax(5),maxsfp,1,par_max)
     call par_reduce(tmpmin(1),minsdu,1,par_min)
     call par_reduce(tmpmin(2),minsdv,1,par_min)
     call par_reduce(tmpmin(3),minsdt,1,par_min)
     call par_reduce(tmpmin(4),minsdq,1,par_min)
     call par_reduce(tmpmin(5),minsfp,1,par_min)
   endif
   if (da_sonde) then
     call par_reduce(tmpmax(1),maxsdu,1,par_max)
     call par_reduce(tmpmax(2),maxsdv,1,par_max)
     call par_reduce(tmpmax(3),maxsdt,1,par_max)
     call par_reduce(tmpmax(4),maxsdq,1,par_max)
     call par_reduce(tmpmin(1),minsdu,1,par_min)
     call par_reduce(tmpmin(2),minsdv,1,par_min)
     call par_reduce(tmpmin(3),minsdt,1,par_min)
     call par_reduce(tmpmin(4),minsdq,1,par_min)
   endif
   if (da_surface) then
     call par_reduce(tmpmax(5),maxsfu,1,par_max)
     call par_reduce(tmpmax(6),maxsfv,1,par_max)
     call par_reduce(tmpmax(7),maxsft,1,par_max)
     call par_reduce(tmpmax(8),maxsfq,1,par_max)
     call par_reduce(tmpmax(9),maxsfp,1,par_max)
     call par_reduce(tmpmin(5),minsfu,1,par_min)
     call par_reduce(tmpmin(6),minsfv,1,par_min)
     call par_reduce(tmpmin(7),minsft,1,par_min)
     call par_reduce(tmpmin(8),minsfq,1,par_min)
     call par_reduce(tmpmin(9),minsfp,1,par_min)
   endif
   if (da_bogus) then
     call par_reduce(tmpmax(10),maxbgp,1,par_max)
     call par_reduce(tmpmin(10),minbgp,1,par_min)
   endif
   if (da_aircraft) then
     call par_reduce(tmpmax(11),maxaru,1,par_max)
     call par_reduce(tmpmax(12),maxarv,1,par_max)
     call par_reduce(tmpmax(13),maxart,1,par_max)
     call par_reduce(tmpmin(11),minaru,1,par_min)
     call par_reduce(tmpmin(12),minarv,1,par_min)
     call par_reduce(tmpmin(13),minart,1,par_min)
   endif
   if (da_amv) then
     call par_reduce(tmpmax(14),maxamvu,1,par_max)
     call par_reduce(tmpmax(15),maxamvv,1,par_max)
     call par_reduce(tmpmin(14),minamvu,1,par_min)
     call par_reduce(tmpmin(15),minamvv,1,par_min)
   endif
   if (da_scatwind) then
     call par_reduce(tmpmax(16),maxscatwindu,1,par_max)
     call par_reduce(tmpmax(17),maxscatwindv,1,par_max)
     call par_reduce(tmpmin(16),minscatwindu,1,par_min)
     call par_reduce(tmpmin(17),minscatwindv,1,par_min)
   endif
   if (par%ismasterproc) then
     print*,'Minimum and maximum of hx_bg-obs'
     if (singleobs) then
       print*,'Single u_bg-u_obs:min = ',minsdu,' max = ',maxsdu
       print*,'Single v_bg-v_obs:min = ',minsdv,' max = ',maxsdv
       print*,'Single t_bg-t_obs:min = ',minsdt,' max = ',maxsdt
       print*,'Single q_bg-q_obs:min = ',minsdq,' max = ',maxsdq
       print*,'Single ps_bg-ps_obs:min = ',minsfp,' max = ',maxsfp
     endif
     if (da_sonde) then
       print*,'Sonde u_bg-u_obs:min = ',minsdu,' max = ',maxsdu
       print*,'Sonde v_bg-v_obs:min = ',minsdv,' max = ',maxsdv
       print*,'Sonde t_bg-t_obs:min = ',minsdt,' max = ',maxsdt
       print*,'Sonde q_bg-q_obs:min = ',minsdq,' max = ',maxsdq
     endif
     if (da_surface) then
       print*,'Surface u_bg-u_obs:min = ',minsfu,' max = ',maxsfu
       print*,'Surface v_bg-v_obs:min = ',minsfv,' max = ',maxsfv
       print*,'Surface t_bg-t_obs:min = ',minsft,' max = ',maxsft
       print*,'Surface q_bg-q_obs:min = ',minsfq,' max = ',maxsfq
       print*,'Surface p_bg-p_obs:min = ',minsfp,' max = ',maxsfp
     endif
     if (da_bogus) then
       print*,'Surface bogus p_bg-p_obs:min = ',minbgp,' max = ',maxbgp
     endif
     if (da_aircraft) then
       print*,'Aircraft u_bg-u_obs:min = ',minaru,' max = ',maxaru
       print*,'Aircraft v_bg-v_obs:min = ',minarv,' max = ',maxarv
       print*,'Aircraft t_bg-t_obs:min = ',minart,' max = ',maxart
     endif
     if (da_amv) then
       print*,'Amv u_bg-u_obs:min = ',minamvu,' max = ',maxamvu
       print*,'Amv v_bg-v_obs:min = ',minamvv,' max = ',maxamvv
     endif
     if (da_scatwind) then
       print*,'Scatwind u_bg-u_obs:min = ',minscatwindu,' max = ',maxscatwindu
       print*,'Scatwind v_bg-v_obs:min = ',minscatwindv,' max = ',maxscatwindv
     endif
   endif
   return
!
   end subroutine check_hxmobs
!===============================================================================
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine rej_obs_omb(itime)
!-------------------------------------------------------------------------------
   use check_omb_mod, only: sondeobcheck, surfaceobcheck, aircraftobcheck, amvobcheck, scatwindobcheck, amsuaobcheck, gpsroobcheck, iasiobcheck, crisobcheck, atmsobcheck, atmswvobcheck, mhsobcheck, csrobcheck
!
   integer(int_kind), intent(in   ) :: itime
   integer(int_kind) :: i, j, k, ie
!=================================================================
!
   call modeltoobs(itime,bg_md(:,:,:,:,itime),gs_obs(itime))
   if (par%ismasterproc) then
     print*,'itime = ',itime
   endif
! sonde observation
   if (da_sonde) then
     obs(itime)%sonde(:,:) = obs_org(itime)%sonde(:,:)
     call sondeobcheck(itime,"u",obs(itime)%sonde%u%val,gs_obs(itime)%sonde%u%val,20.d0)
     call sondeobcheck(itime,"v",obs(itime)%sonde%v%val,gs_obs(itime)%sonde%v%val,20.d0)
     call sondeobcheck(itime,"t",obs(itime)%sonde%t%val,gs_obs(itime)%sonde%t%val,20.d0)
     if (cv_opt_hum==0) then
       call sondeobcheck(itime,"q",obs(itime)%sonde%q%val,gs_obs(itime)%sonde%q%val,1.d-2)
     elseif (cv_opt_hum==1) then
       call sondeobcheck(itime,"q",obs(itime)%sonde%q%val,gs_obs(itime)%sonde%q%val,2.d0)
     endif
   endif
! surface observation
   if (da_surface) then
     obs(itime)%surface(:,:) = obs_org(itime)%surface(:,:)
     call surfaceobcheck(itime,obs(itime)%surface%u%val,gs_obs(itime)%surface%u%val,20.d0)
     call surfaceobcheck(itime,obs(itime)%surface%v%val,gs_obs(itime)%surface%v%val,20.d0)
     call surfaceobcheck(itime,obs(itime)%surface%t%val,gs_obs(itime)%surface%t%val,20.d0)
     call surfaceobcheck(itime,obs(itime)%surface%q%val,gs_obs(itime)%surface%q%val,1.d-2)
     call surfaceobcheck(itime,obs(itime)%surface%p%val,gs_obs(itime)%surface%p%val,1000.d0)
   endif
! aircraft observation
   if (da_aircraft) then
     obs(itime)%aircraft(:,:) = obs_org(itime)%aircraft(:,:)
     call aircraftobcheck(itime,"u",obs(itime)%aircraft%u%val,gs_obs(itime)%aircraft%u%val,20.d0)
     call aircraftobcheck(itime,"v",obs(itime)%aircraft%v%val,gs_obs(itime)%aircraft%v%val,20.d0)
     call aircraftobcheck(itime,"t",obs(itime)%aircraft%t%val,gs_obs(itime)%aircraft%t%val,20.d0)
   endif
! amv observation
   if (da_amv) then
     obs(itime)%amv(:,:) = obs_org(itime)%amv(:,:)
     call amvobcheck(itime,obs(itime)%amv%u%val,gs_obs(itime)%amv%u%val,15.d0)
     call amvobcheck(itime,obs(itime)%amv%v%val,gs_obs(itime)%amv%v%val,15.d0)
   endif
! scatwind observation
   if (da_scatwind) then
     obs(itime)%scatwind(:,:) = obs_org(itime)%scatwind(:,:)
     call scatwindobcheck(itime,obs(itime)%scatwind%u%val,gs_obs(itime)%scatwind%u%val,5.d0)
     call scatwindobcheck(itime,obs(itime)%scatwind%v%val,gs_obs(itime)%scatwind%v%val,5.d0)
   endif
! amsua
   if (da_amsua) then
     obs(itime)%amsua(:,:,:) = obs_org(itime)%amsua(:,:,:)
     call amsuaobcheck(itime,obs(itime)%amsua%tb%val,dsqrt(obs(itime)%amsua%tb%error),obsmet(itime)%amsuabg%amsuahdx%tb%val,crit_omb_amsua)
   endif
! gpsroOBCheck
   if (da_gpsro) then
     obs(itime)%gpsro(:,:) = obs_org(itime)%gpsro(:,:)
     call gpsroobcheck(itime,obs(itime)%gpsro%ba%val,dsqrt(obs(itime)%gpsro%ba%error),obsmet(itime)%gpsrobg%gpsrohdx%ba%val,crit_omb_gpsro)
   endif
! iasiOBCheck
   if (da_iasi) then
     obs(itime)%iasi(:,:,:) = obs_org(itime)%iasi(:,:,:)
     call iasiobcheck(itime,obs(itime)%iasi%tb%val,dsqrt(obs(itime)%iasi%tb%error),obsmet(itime)%iasibg%iasihdx%tb%val,crit_omb_iasi)
   endif
! crisOBCheck
   if (da_cris) then
     obs(itime)%cris(:,:,:) = obs_org(itime)%cris(:,:,:)
     call crisobcheck(itime,obs(itime)%cris%tb%val,dsqrt(obs(itime)%cris%tb%error),obsmet(itime)%crisbg%crishdx%tb%val,crit_omb_cris)
   endif
! atmsOBCheck
   if (da_atms) then
     obs(itime)%atms(:,:,:) = obs_org(itime)%atms(:,:,:)
     call atmsobcheck(itime,obs(itime)%atms%tb%val,dsqrt(obs(itime)%atms%tb%error),obsmet(itime)%atmsbg%atmshdx%tb%val,crit_omb_atms)
   endif
! atmswvOBCheck
   if (da_atmswv) then
     obs(itime)%atmswv(:,:,:) = obs_org(itime)%atmswv(:,:,:)
     call atmswvobcheck(itime,obs(itime)%atmswv%tb%val,dsqrt(obs(itime)%atmswv%tb%error),obsmet(itime)%atmswvbg%atmswvhdx%tb%val,crit_omb_atmswv)
   endif
! mhsOBCheck
   if (da_mhs) then
     obs(itime)%mhs(:,:,:) = obs_org(itime)%mhs(:,:,:)
     call mhsobcheck(itime,obs(itime)%mhs%tb%val,obsmet(itime)%mhsbg%mhshdx%tb%val,crit_omb_mhs)
   endif
! csrOBCheck
   if (da_csr) then
     obs(itime)%csr(:,:,:) = obs_org(itime)%csr(:,:,:)
     call csrobcheck(itime,obs(itime)%csr%tb%val,dsqrt(obs(itime)%csr%tb%error),obsmet(itime)%csrbg%csrhdx%tb%val,crit_omb_csr)
   endif
   return
!
   end subroutine rej_obs_omb
!===============================================================================
!===============================================================================
! SUBROUTINE Check_RMSE(an_md,bg_md)
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine check_rmse
! 2. Variables
! 2-1. Input
!   REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd), &
!      INTENT(IN) :: an_md
!      INTENT(IN) :: bg_md
! 2-2. Local
!
   integer(int_kind) :: i, j, k, ie, wn, vr, n3d, n2d, nsonde, nsurfc, naircra
   real(real_kind) :: rmse_ua, rmse_va, rmse_ta, rmse_qa, rmse_psa, rmse_ub, rmse_vb, rmse_tb, rmse_qb, rmse_psb, rmse_sduo, rmse_sdvo, rmse_sdto, rmse_sdqo
!
   rmse_ub = 0.d0; rmse_vb = 0.d0; rmse_tb = 0.d0; rmse_qb = 0.d0; rmse_psb = 0.d0
   rmse_ua = 0.d0; rmse_va = 0.d0; rmse_ta = 0.d0; rmse_qa = 0.d0; rmse_psa = 0.d0
   rmse_sduo = 0.d0; rmse_sdvo = 0.d0; rmse_sdto = 0.d0; rmse_sdqo = 0.d0
!    ierr = TimingStart('ModelToObs_Make_RMSE')  !--- profiling
!CALL ModelToObs(an_md,gs_obs)
!    ierr = TimingStop('ModelToObs_Make_RMSE')
   global_shared_buf(:,1:18) = 0.d0
!.. Calculate sonde observation RMSE error
!   IF (da_sonde) THEN
!   DO ie = nets,nete
!   DO i  = 1,nsd_e(1,ie)
!      IF ((obs%sonde(i,ie)%u%val .NE. sdmiss) .AND. (gs_obs%sonde(i,ie)%u%val .NE. sdmiss)) THEN
!         global_shared_buf(ie,2) = global_shared_buf(ie,2)               &
!                                 +(gs_obs%sonde(i,ie)%u%val - obs%sonde(i,ie)%u%val)**2.D0
!         global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
!      END IF
!   END DO
!   DO i  = 1,nsd_e(2,ie)
!      IF ((obs%sonde(i,ie)%v%val .NE. sdmiss) .AND. (gs_obs%sonde(i,ie)%v%val .NE. sdmiss)) THEN
!         global_shared_buf(ie,3) = global_shared_buf(ie,3)               &
!                                 +(gs_obs%sonde(i,ie)%v%val - obs%sonde(i,ie)%v%val)**2.D0
!         global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
!      END IF
!   END DO
!   DO i  = 1,nsd_e(3,ie)
!      IF ((obs%sonde(i,ie)%t%val .NE. sdmiss) .AND. (gs_obs%sonde(i,ie)%t%val .NE. sdmiss)) THEN
!         global_shared_buf(ie,4) = global_shared_buf(ie,4)               &
!                                 +(gs_obs%sonde(i,ie)%t%val - obs%sonde(i,ie)%t%val)**2.D0
!         global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
!      END IF
!   END DO
!   DO i  = 1,nsd_e(4,ie)
!      IF ((obs%sonde(i,ie)%q%val .NE. sdmiss) .AND. (gs_obs%sonde(i,ie)%q%val .NE. sdmiss)) THEN
!         global_shared_buf(ie,5) = global_shared_buf(ie,5)               &
!                                 +(gs_obs%sonde(i,ie)%q%val - obs%sonde(i,ie)%q%val)**2.D0
!         global_shared_buf(ie,6) = global_shared_buf(ie,6) + 1.D0
!      END IF
!   END DO
!   END DO
!   END IF
   do ie = 1,nelemd
     do k = 11, nlev
       do j = 1, np
         do i = 1, np
           global_shared_buf(ie,7) = global_shared_buf(ie,7)+1.d0
!.. Calculate Background error (U,V,T,and Q)
           global_shared_buf(ie,8) = global_shared_buf(ie,8)+(bg_md(i,j,k,ie,nt_itime)-nt_org(i,j,k,ie))**2.d0
           global_shared_buf(ie,9) = global_shared_buf(ie,9)+(bg_md(i,j,nlev+k,ie,nt_itime)-nt_org(i,j,nlev+k,ie))**2.d0
           global_shared_buf(ie,10) = global_shared_buf(ie,10)+(bg_md(i,j,2*nlev+k,ie,nt_itime)-nt_org(i,j,2*nlev+k,ie))**2.d0
           global_shared_buf(ie,11) = global_shared_buf(ie,11)+(bg_md(i,j,3*nlev+k,ie,nt_itime)-nt_org(i,j,3*nlev+k,ie))**2.d0
!.. Calculate Analysis error (U,V,T,and Q)
           global_shared_buf(ie,12) = global_shared_buf(ie,12)+(an_md(i,j,k,ie,nt_itime)-nt_org(i,j,k,ie))**2.d0
           global_shared_buf(ie,13) = global_shared_buf(ie,13)+(an_md(i,j,nlev+k,ie,nt_itime)-nt_org(i,j,nlev+k,ie))**2.d0
           global_shared_buf(ie,14) = global_shared_buf(ie,14)+(an_md(i,j,2*nlev+k,ie,nt_itime)-nt_org(i,j,2*nlev+k,ie))**2.d0
           global_shared_buf(ie,15) = global_shared_buf(ie,15)+(an_md(i,j,3*nlev+k,ie,nt_itime)-nt_org(i,j,3*nlev+k,ie))**2.d0
         enddo
       enddo
     enddo
   enddo
   do ie = 1,nelemd
     do j = 1, np
       do i = 1, np
         global_shared_buf(ie,16) = global_shared_buf(ie,16)+1.d0
!.. Calculate Background error (PS)
         global_shared_buf(ie,17) = global_shared_buf(ie,17)+(bg_md(i,j,4*nlev+1,ie,nt_itime)-nt_org(i,j,4*nlev+1,ie))**2.d0
!.. Calculate Analysis error (PS)
         global_shared_buf(ie,18) = global_shared_buf(ie,18)+(an_md(i,j,4*nlev+1,ie,nt_itime)-nt_org(i,j,4*nlev+1,ie))**2.d0
       enddo
     enddo
   enddo
   ierr = timingstart('Repro_Sum_Make_RMSE') !---profiling
   call wrap_repro_sum(nvars = 18,comm = hybrid%par%comm)
   ierr = timingstop('Repro_Sum_Make_RMSE') !---profiling
!sonde observation
!   rmse_sduo = global_shared_sum(2)
!   rmse_sdvo = global_shared_sum(3)
!   rmse_sdto = global_shared_sum(4)
!   rmse_sdqo = global_shared_sum(5)
!   nsonde    = global_shared_sum(6)
   n3d = global_shared_sum(7)
   rmse_ub = global_shared_sum(8)
   rmse_vb = global_shared_sum(9)
   rmse_tb = global_shared_sum(10)
   rmse_qb = global_shared_sum(11)
   rmse_ua = global_shared_sum(12)
   rmse_va = global_shared_sum(13)
   rmse_ta = global_shared_sum(14)
   rmse_qa = global_shared_sum(15)
   n2d = global_shared_sum(16)
   rmse_psb = global_shared_sum(17)
   rmse_psa = global_shared_sum(18)
!   IF (da_sonde) THEN
!   IF (par%isMasterProc) WRITE(iulog,*) "RMSE against Sonde observation"
!   IF (par%isMasterProc) WRITE(iulog,*) "U  =",DSQRT(rmse_sduo/nsonde)
!   IF (par%isMasterProc) WRITE(iulog,*) "V  =",DSQRT(rmse_sdvo/nsonde)
!   IF (par%isMasterProc) WRITE(iulog,*) "T  =",DSQRT(rmse_sdto/nsonde)
!   IF (par%isMasterProc) WRITE(iulog,*) "Q  =",DSQRT(rmse_sdqo/nsonde)
!   END IF
   if (par%ismasterproc) write(iulog,*) "rmse of u(from level 11 to bottom) "
   if (par%ismasterproc) write(iulog,*) "background = ",dsqrt(rmse_ub/n3d)
   if (par%ismasterproc) write(iulog,*) "analysis = ",dsqrt(rmse_ua/n3d)
   if (par%ismasterproc) write(iulog,*) "rmse of v(from level 11 to bottom) "
   if (par%ismasterproc) write(iulog,*) "background = ",dsqrt(rmse_vb/n3d)
   if (par%ismasterproc) write(iulog,*) "analysis = ",dsqrt(rmse_va/n3d)
   if (par%ismasterproc) write(iulog,*) "rmse of t(from level 11 to bottom) "
   if (par%ismasterproc) write(iulog,*) "background = ",dsqrt(rmse_tb/n3d)
   if (par%ismasterproc) write(iulog,*) "analysis = ",dsqrt(rmse_ta/n3d)
   if (par%ismasterproc) write(iulog,*) "rmse of q(from level 11 to bottom) "
   if (par%ismasterproc) write(iulog,*) "background = ",dsqrt(rmse_qb/n3d)
   if (par%ismasterproc) write(iulog,*) "analysis = ",dsqrt(rmse_qa/n3d)
   if (par%ismasterproc) write(iulog,*) "rmse of ps"
   if (par%ismasterproc) write(iulog,*) "background = ",dsqrt(rmse_psb/n2d)
   if (par%ismasterproc) write(iulog,*) "analysis = ",dsqrt(rmse_psa/n2d)
   return
!
   end subroutine check_rmse
!===============================================================================
!===============================================================================
! SUBROUTINE Check_OB(an_md, bg_md)
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine check_ob(itime)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(int_kind), intent(in   ) :: itime
!   REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd), &
!      INTENT(IN) :: an_md, bg_md
!  Local
   real(real_kind) :: obmsdu, obmsdv, obmsdt, obmsdq ! o-b, average, sonde
   real(real_kind) :: obmsfu, obmsfv, obmsft, obmsfq, obmsfps ! o-b, average, surface
   real(real_kind) :: obdsdu, obdsdv, obdsdt, obdsdq ! o-b, sdv, sonde
   real(real_kind) :: obdsfu, obdsfv, obdsft, obdsfq, obdsfps ! o-b, sdv, surface
   real(real_kind) :: obmbgps ! o-b, average, bogus
   real(real_kind) :: obdbgps ! o-b, sdv, bogus
   real(real_kind) :: oamsdu, oamsdv, oamsdt, oamsdq ! o-a, average, sonde
   real(real_kind) :: oamsfu, oamsfv, oamsft, oamsfq, oamsfps ! o-a, average, surface
   real(real_kind) :: oadsdu, oadsdv, oadsdt, oadsdq ! o-a, sdv, sonde
   real(real_kind) :: oadsfu, oadsfv, oadsft, oadsfq, oadsfps ! o-a, sdv, surface
   real(real_kind) :: oambgps ! o-a, average, bogus
   real(real_kind) :: oadbgps ! o-a, sdv, bogus
   real(real_kind) :: nobsdu, nobsdv, nobsdt, nobsdq ! number of o-b, sonde
   real(real_kind) :: nobsfu, nobsfv, nobsft, nobsfq, nobsfps ! number of o-b, surface
   real(real_kind) :: noasdu, noasdv, noasdt, noasdq ! number of o-a, sonde
   real(real_kind) :: noasfu, noasfv, noasft, noasfq, noasfps ! number of o-a, surface
   real(real_kind) :: nobbgps ! number of o-b, bogus
   real(real_kind) :: noabgps ! number of o-a, bogus
   integer(int_kind) :: i, k, ie, ierr
! background
!
   call modeltoobs(itime,bg_md(:,:,:,:,itime),bg_obs(itime))
! analysis
   call modeltoobs(itime,an_md(:,:,:,:,itime),an_obs(itime))
!.. Mean of O-B and O-A for Single Observatino
   if (singleobs) then
     do ie = nets, nete
       do i = 1, obsmet(itime)%nsg_e(ie)
         print*,'i,ie:',i,ie
         print*,'bg_obs(itime)%single(i,ie)%u%val,an_obs(itime)%single(i,ie)%u%val',bg_obs(itime)%single(i,ie)%u%val,an_obs(itime)%single(i,ie)%u%val
         print*,'bg_obs(itime)%single(i,ie)%v%val,an_obs(itime)%single(i,ie)%v%val',bg_obs(itime)%single(i,ie)%v%val,an_obs(itime)%single(i,ie)%v%val
         print*,'bg_obs(itime)%single(i,ie)%t%val,an_obs(itime)%single(i,ie)%t%val',bg_obs(itime)%single(i,ie)%t%val,an_obs(itime)%single(i,ie)%t%val
         print*,'bg_obs(itime)%single(i,ie)%q%val,an_obs(itime)%single(i,ie)%q%val',bg_obs(itime)%single(i,ie)%q%val,an_obs(itime)%single(i,ie)%q%val
         print*,'bg_obs(itime)%single(i,ie)%ps%val,an_obs(itime)%single(i,ie)%ps%val',bg_obs(itime)%single(i,ie)%ps%val,an_obs(itime)%single(i,ie)%ps%val
       enddo
     enddo
   endif
!.. Mean of O-B and O-A for Sonde
   if (da_sonde) then
     obmsdu = 0.d0; obmsdv = 0.d0; obmsdt = 0.d0; obmsdq = 0.d0
     obdsdu = 0.d0; obdsdv = 0.d0; obdsdt = 0.d0; obdsdq = 0.d0
     oamsdu = 0.d0; oamsdv = 0.d0; oamsdt = 0.d0; oamsdq = 0.d0
     oadsdu = 0.d0; oadsdv = 0.d0; oadsdt = 0.d0; oadsdq = 0.d0
     nobsdu = 0 ; nobsdv = 0 ; nobsdt = 0 ; nobsdq = 0
     noasdu = 0 ; noasdv = 0 ; noasdt = 0 ; noasdq = 0
     global_shared_buf(:,1:18) = 0.d0
     do ie = nets,nete
! background
       do i = 1, obsmet(itime)%nsd_e(1, ie)
         if ((obs(itime)%sonde(i, ie)%u%val.ne.sdmiss).and.(bg_obs(itime)%sonde(i, ie)%u%val.ne.sdmiss)) then
           global_shared_buf(ie,1) = global_shared_buf(ie,1)+(obs(itime)%sonde(i,ie)%u%val-bg_obs(itime)%sonde(i,ie)%u%val)
           global_shared_buf(ie,11) = global_shared_buf(ie,11)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(2, ie)
         if ((obs(itime)%sonde(i, ie)%v%val.ne.sdmiss).and.(bg_obs(itime)%sonde(i, ie)%v%val.ne.sdmiss)) then
           global_shared_buf(ie,2) = global_shared_buf(ie,2)+(obs(itime)%sonde(i,ie)%v%val-bg_obs(itime)%sonde(i,ie)%v%val)
           global_shared_buf(ie,12) = global_shared_buf(ie,12)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(3, ie)
         if ((obs(itime)%sonde(i, ie)%t%val.ne.sdmiss).and.(bg_obs(itime)%sonde(i, ie)%t%val.ne.sdmiss)) then
           global_shared_buf(ie,3) = global_shared_buf(ie,3)+(obs(itime)%sonde(i,ie)%t%val-bg_obs(itime)%sonde(i,ie)%t%val)
           global_shared_buf(ie,13) = global_shared_buf(ie,13)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(4, ie)
         if ((obs(itime)%sonde(i, ie)%q%val.ne.sdmiss).and.(bg_obs(itime)%sonde(i, ie)%q%val.ne.sdmiss)) then
           global_shared_buf(ie,4) = global_shared_buf(ie,4)+(obs(itime)%sonde(i,ie)%q%val-bg_obs(itime)%sonde(i,ie)%q%val)
           global_shared_buf(ie,14) = global_shared_buf(ie,14)+1.d0
         endif
       enddo
! analysis
       do i = 1, obsmet(itime)%nsd_e(1, ie)
         if ((obs(itime)%sonde(i, ie)%u%val.ne.sdmiss).and.(an_obs(itime)%sonde(i, ie)%u%val.ne.sdmiss)) then
           global_shared_buf(ie,5) = global_shared_buf(ie,5)+(obs(itime)%sonde(i,ie)%u%val-an_obs(itime)%sonde(i,ie)%u%val)
           global_shared_buf(ie,15) = global_shared_buf(ie,15)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(2, ie)
         if ((obs(itime)%sonde(i, ie)%v%val.ne.sdmiss).and.(an_obs(itime)%sonde(i, ie)%v%val.ne.sdmiss)) then
           global_shared_buf(ie,6) = global_shared_buf(ie,6)+(obs(itime)%sonde(i,ie)%v%val-an_obs(itime)%sonde(i,ie)%v%val)
           global_shared_buf(ie,16) = global_shared_buf(ie,16)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(3, ie)
         if ((obs(itime)%sonde(i, ie)%t%val.ne.sdmiss).and.(an_obs(itime)%sonde(i, ie)%t%val.ne.sdmiss)) then
           global_shared_buf(ie,7) = global_shared_buf(ie,7)+(obs(itime)%sonde(i,ie)%t%val-an_obs(itime)%sonde(i,ie)%t%val)
           global_shared_buf(ie,17) = global_shared_buf(ie,17)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(4, ie)
         if ((obs(itime)%sonde(i, ie)%q%val.ne.sdmiss).and.(an_obs(itime)%sonde(i, ie)%q%val.ne.sdmiss)) then
           global_shared_buf(ie,8) = global_shared_buf(ie,8)+(obs(itime)%sonde(i,ie)%q%val-an_obs(itime)%sonde(i,ie)%q%val)
           global_shared_buf(ie,18) = global_shared_buf(ie,18)+1.d0
         endif
       enddo
     enddo ! ie
     call wrap_repro_sum(nvars = 18,comm = hybrid%par%comm)
! O-B
     nobsdu = global_shared_sum(11)
     nobsdv = global_shared_sum(12)
     nobsdt = global_shared_sum(13)
     nobsdq = global_shared_sum(14)
     obmsdu = global_shared_sum(1)/nobsdu
     obmsdv = global_shared_sum(2)/nobsdv
     obmsdt = global_shared_sum(3)/nobsdt
     obmsdq = global_shared_sum(4)/nobsdq
! O-A
     noasdu = global_shared_sum(15)
     noasdv = global_shared_sum(16)
     noasdt = global_shared_sum(17)
     noasdq = global_shared_sum(18)
     oamsdu = global_shared_sum(5)/noasdu
     oamsdv = global_shared_sum(6)/noasdv
     oamsdt = global_shared_sum(7)/noasdt
     oamsdq = global_shared_sum(8)/noasdq
   endif ! da_sonde
!.. Mean of O-B and O-A for Surface observation
   if (da_surface) then
     obmsfu = 0.d0; obmsfv = 0.d0; obmsft = 0.d0; obmsfq = 0.d0; obmsfps = 0.d0
     obdsfu = 0.d0; obdsfv = 0.d0; obdsft = 0.d0; obdsfq = 0.d0; obdsfps = 0.d0
     oamsfu = 0.d0; oamsfv = 0.d0; oamsft = 0.d0; oamsfq = 0.d0; oamsfps = 0.d0
     oadsfu = 0.d0; oadsfv = 0.d0; oadsft = 0.d0; oadsfq = 0.d0; oadsfps = 0.d0
     nobsfu = 0 ; nobsfv = 0 ; nobsft = 0 ; nobsfq = 0 ; nobsfps = 0
     noasfu = 0 ; noasfv = 0 ; noasft = 0 ; noasfq = 0 ; noasfps = 0
     global_shared_buf(:,1:20) = 0.d0
     do ie = nets,nete
! background
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%u%val.ne.sdmiss).and.(bg_obs(itime)%surface(i, ie)%u%val.ne.sdmiss)) then
           global_shared_buf(ie,1) = global_shared_buf(ie,1)+(obs(itime)%surface(i,ie)%u%val-bg_obs(itime)%surface(i,ie)%u%val)
           global_shared_buf(ie,11) = global_shared_buf(ie,11)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%v%val.ne.sdmiss).and.(bg_obs(itime)%surface(i, ie)%v%val.ne.sdmiss)) then
           global_shared_buf(ie,2) = global_shared_buf(ie,2)+(obs(itime)%surface(i,ie)%v%val-bg_obs(itime)%surface(i,ie)%v%val)
           global_shared_buf(ie,12) = global_shared_buf(ie,12)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%t%val.ne.sdmiss).and.(bg_obs(itime)%surface(i, ie)%t%val.ne.sdmiss)) then
           global_shared_buf(ie,3) = global_shared_buf(ie,3)+(obs(itime)%surface(i,ie)%t%val-bg_obs(itime)%surface(i,ie)%t%val)
           global_shared_buf(ie,13) = global_shared_buf(ie,13)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%q%val.ne.sdmiss).and.(bg_obs(itime)%surface(i, ie)%q%val.ne.sdmiss)) then
           global_shared_buf(ie,4) = global_shared_buf(ie,4)+(obs(itime)%surface(i,ie)%q%val-bg_obs(itime)%surface(i,ie)%q%val)
           global_shared_buf(ie,14) = global_shared_buf(ie,14)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%p%val.ne.sdmiss).and.(bg_obs(itime)%surface(i, ie)%p%val.ne.sdmiss)) then
           global_shared_buf(ie,5) = global_shared_buf(ie,5)+(obs(itime)%surface(i,ie)%p%val-bg_obs(itime)%surface(i,ie)%p%val)
           global_shared_buf(ie,15) = global_shared_buf(ie,15)+1.d0
         endif
       enddo
! analysis
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%u%val.ne.sdmiss).and.(an_obs(itime)%surface(i, ie)%u%val.ne.sdmiss)) then
           global_shared_buf(ie,6) = global_shared_buf(ie,6)+(obs(itime)%surface(i,ie)%u%val-an_obs(itime)%surface(i,ie)%u%val)
           global_shared_buf(ie,16) = global_shared_buf(ie,16)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%v%val.ne.sdmiss).and.(an_obs(itime)%surface(i, ie)%v%val.ne.sdmiss)) then
           global_shared_buf(ie,7) = global_shared_buf(ie,7)+(obs(itime)%surface(i,ie)%v%val-an_obs(itime)%surface(i,ie)%v%val)
           global_shared_buf(ie,17) = global_shared_buf(ie,17)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%t%val.ne.sdmiss).and.(an_obs(itime)%surface(i, ie)%t%val.ne.sdmiss)) then
           global_shared_buf(ie,8) = global_shared_buf(ie,8)+(obs(itime)%surface(i,ie)%t%val-an_obs(itime)%surface(i,ie)%t%val)
           global_shared_buf(ie,18) = global_shared_buf(ie,18)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%q%val.ne.sdmiss).and.(an_obs(itime)%surface(i, ie)%q%val.ne.sdmiss)) then
           global_shared_buf(ie,9) = global_shared_buf(ie,9)+(obs(itime)%surface(i,ie)%q%val-an_obs(itime)%surface(i,ie)%q%val)
           global_shared_buf(ie,19) = global_shared_buf(ie,19)+1.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%p%val.ne.sdmiss).and.(an_obs(itime)%surface(i, ie)%p%val.ne.sdmiss)) then
           global_shared_buf(ie,10) = global_shared_buf(ie,10)+(obs(itime)%surface(i,ie)%p%val-an_obs(itime)%surface(i,ie)%p%val)
           global_shared_buf(ie,20) = global_shared_buf(ie,20)+1.d0
         endif
       enddo
     enddo ! ie
     call wrap_repro_sum(nvars = 20,comm = hybrid%par%comm)
! O-B
     nobsfu = global_shared_sum(11)
     nobsfv = global_shared_sum(12)
     nobsft = global_shared_sum(13)
     nobsfq = global_shared_sum(14)
     nobsfps = global_shared_sum(15)
     obmsfu = global_shared_sum(1)/nobsfu
     obmsfv = global_shared_sum(2)/nobsfv
     obmsft = global_shared_sum(3)/nobsft
     obmsfq = global_shared_sum(4)/nobsfq
     obmsfps = global_shared_sum(5)/nobsfps
! O-A
     noasfu = global_shared_sum(16)
     noasfv = global_shared_sum(17)
     noasft = global_shared_sum(18)
     noasfq = global_shared_sum(19)
     noasfps = global_shared_sum(20)
     oamsfu = global_shared_sum(6)/noasfu
     oamsfv = global_shared_sum(7)/noasfv
     oamsft = global_shared_sum(8)/noasft
     oamsfq = global_shared_sum(9)/noasfq
     oamsfps = global_shared_sum(10)/noasfps
   endif ! da_surface
   if (da_bogus) then
     obmbgps = 0.d0
     obdbgps = 0.d0
     oambgps = 0.d0
     oadbgps = 0.d0
     nobbgps = 0
     noabgps = 0
     global_shared_buf(:,1:4) = 0.d0
     do ie = nets,nete
! background
       do i = 1, obsmet(itime)%nbg_e(ie)
         if ((obs(itime)%bogus(i, ie)%p%val.ne.sdmiss).and.(bg_obs(itime)%bogus(i, ie)%p%val.ne.sdmiss)) then
           global_shared_buf(ie,1) = global_shared_buf(ie,1)+(obs(itime)%bogus(i,ie)%p%val-bg_obs(itime)%bogus(i,ie)%p%val)
           global_shared_buf(ie,3) = global_shared_buf(ie,3)+1.d0
         endif
       enddo
! analysis
       do i = 1, obsmet(itime)%nbg_e(ie)
         if ((obs(itime)%bogus(i, ie)%p%val.ne.sdmiss).and.(an_obs(itime)%bogus(i, ie)%p%val.ne.sdmiss)) then
           global_shared_buf(ie,2) = global_shared_buf(ie,2)+(obs(itime)%bogus(i,ie)%p%val-an_obs(itime)%bogus(i,ie)%p%val)
           global_shared_buf(ie,4) = global_shared_buf(ie,4)+1.d0
         endif
       enddo
     enddo
     call wrap_repro_sum(nvars = 4,comm = hybrid%par%comm)
! O-B
     nobbgps = global_shared_sum(3)
     obmbgps = global_shared_sum(1)/nobbgps
! O-A
     noabgps = global_shared_sum(4)
     oambgps = global_shared_sum(2)/noabgps
   endif ! da_bogus
!.. Print out the mean of O-B and O-A
   if (da_sonde) then
     if (par%ismasterproc) write(iulog,*) "mean of o-b and o-a"
     if (par%ismasterproc) write(iulog,*) "itime = ",itime
     if (par%ismasterproc) write(iulog,*) "against sonde observation"
     if (par%ismasterproc) write(iulog,*) "o-b u = ",obmsdu
     if (par%ismasterproc) write(iulog,*) "o-a u = ",oamsdu
     if (par%ismasterproc) write(iulog,*) "o-b v = ",obmsdv
     if (par%ismasterproc) write(iulog,*) "o-a v = ",oamsdv
     if (par%ismasterproc) write(iulog,*) "o-b t = ",obmsdt
     if (par%ismasterproc) write(iulog,*) "o-a t = ",oamsdt
     if (par%ismasterproc) write(iulog,*) "o-b q = ",obmsdq
     if (par%ismasterproc) write(iulog,*) "o-a q = ",oamsdq
   endif
   if (da_surface) then
     if (par%ismasterproc) write(iulog,*) "agaist surface observation"
     if (par%ismasterproc) write(iulog,*) "o-b u = ",obmsfu
     if (par%ismasterproc) write(iulog,*) "o-a u = ",oamsfu
     if (par%ismasterproc) write(iulog,*) "o-b v = ",obmsfv
     if (par%ismasterproc) write(iulog,*) "o-a v = ",oamsfv
     if (par%ismasterproc) write(iulog,*) "o-b t = ",obmsft
     if (par%ismasterproc) write(iulog,*) "o-a t = ",oamsft
     if (par%ismasterproc) write(iulog,*) "o-b q = ",obmsfq
     if (par%ismasterproc) write(iulog,*) "o-a q = ",oamsfq
     if (par%ismasterproc) write(iulog,*) "o-b ps = ",obmsfps
     if (par%ismasterproc) write(iulog,*) "o-a ps = ",oamsfps
   endif
   if (da_bogus) then
     if (par%ismasterproc) write(iulog,*) "agaist surface bogus observation"
     if (par%ismasterproc) write(iulog,*) "o-b ps = ",obmbgps
     if (par%ismasterproc) write(iulog,*) "o-a ps = ",oambgps
   endif
!.. Standard deviation of O-B and O-A for Sonde
   if (da_sonde) then
     global_shared_buf(:,1:8) = 0.d0
     do ie = nets,nete
! background
       do i = 1, obsmet(itime)%nsd_e(1, ie)
         if ((obs(itime)%sonde(i, ie)%u%val.ne.sdmiss).and.(bg_obs(itime)%sonde(i, ie)%u%val.ne.sdmiss)) then
           global_shared_buf(ie,1) = global_shared_buf(ie,1)+((obs(itime)%sonde(i,ie)%u%val-bg_obs(itime)%sonde(i,ie)%u%val)-obmsdu)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(2, ie)
         if ((obs(itime)%sonde(i, ie)%v%val.ne.sdmiss).and.(bg_obs(itime)%sonde(i, ie)%v%val.ne.sdmiss)) then
           global_shared_buf(ie,2) = global_shared_buf(ie,2)+((obs(itime)%sonde(i,ie)%v%val-bg_obs(itime)%sonde(i,ie)%v%val)-obmsdv)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(3, ie)
         if ((obs(itime)%sonde(i, ie)%t%val.ne.sdmiss).and.(bg_obs(itime)%sonde(i, ie)%t%val.ne.sdmiss)) then
           global_shared_buf(ie,3) = global_shared_buf(ie,3)+((obs(itime)%sonde(i,ie)%t%val-bg_obs(itime)%sonde(i,ie)%t%val)-obmsdt)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(4, ie)
         if ((obs(itime)%sonde(i, ie)%q%val.ne.sdmiss).and.(bg_obs(itime)%sonde(i, ie)%q%val.ne.sdmiss)) then
           global_shared_buf(ie,4) = global_shared_buf(ie,4)+((obs(itime)%sonde(i,ie)%q%val-bg_obs(itime)%sonde(i,ie)%q%val)-obmsdq)**2.d0
         endif
       enddo
! analysis
       do i = 1, obsmet(itime)%nsd_e(1, ie)
         if ((obs(itime)%sonde(i, ie)%u%val.ne.sdmiss).and.(an_obs(itime)%sonde(i, ie)%u%val.ne.sdmiss)) then
           global_shared_buf(ie,5) = global_shared_buf(ie,5)+((obs(itime)%sonde(i,ie)%u%val-an_obs(itime)%sonde(i,ie)%u%val)-oamsdu)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(2, ie)
         if ((obs(itime)%sonde(i, ie)%v%val.ne.sdmiss).and.(an_obs(itime)%sonde(i, ie)%v%val.ne.sdmiss)) then
           global_shared_buf(ie,6) = global_shared_buf(ie,6)+((obs(itime)%sonde(i,ie)%v%val-an_obs(itime)%sonde(i,ie)%v%val)-oamsdv)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(3, ie)
         if ((obs(itime)%sonde(i, ie)%t%val.ne.sdmiss).and.(an_obs(itime)%sonde(i, ie)%t%val.ne.sdmiss)) then
           global_shared_buf(ie,7) = global_shared_buf(ie,7)+((obs(itime)%sonde(i,ie)%t%val-an_obs(itime)%sonde(i,ie)%t%val)-oamsdt)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsd_e(4, ie)
         if ((obs(itime)%sonde(i, ie)%q%val.ne.sdmiss).and.(an_obs(itime)%sonde(i, ie)%q%val.ne.sdmiss)) then
           global_shared_buf(ie,8) = global_shared_buf(ie,8)+((obs(itime)%sonde(i,ie)%q%val-an_obs(itime)%sonde(i,ie)%q%val)-oamsdq)**2.d0
         endif
       enddo
     enddo ! ie
     call wrap_repro_sum(nvars = 8,comm = hybrid%par%comm)
! O-B
     obdsdu = sqrt(global_shared_sum(1)/nobsdu)
     obdsdv = sqrt(global_shared_sum(2)/nobsdv)
     obdsdt = sqrt(global_shared_sum(3)/nobsdt)
     obdsdq = sqrt(global_shared_sum(4)/nobsdq)
! O-A
     oadsdu = sqrt(global_shared_sum(5)/noasdu)
     oadsdv = sqrt(global_shared_sum(6)/noasdv)
     oadsdt = sqrt(global_shared_sum(7)/noasdt)
     oadsdq = sqrt(global_shared_sum(8)/noasdq)
   endif ! da_sonde
!.. Standard deviation of O-B and O-A for Surface
   if (da_surface) then
     global_shared_buf(:,1:10) = 0.d0
     do ie = nets,nete
! background
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%u%val.ne.sdmiss).and.(bg_obs(itime)%surface(i, ie)%u%val.ne.sdmiss)) then
           global_shared_buf(ie,1) = global_shared_buf(ie,1)+((obs(itime)%surface(i,ie)%u%val-bg_obs(itime)%surface(i,ie)%u%val)-obmsfu)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%v%val.ne.sdmiss).and.(bg_obs(itime)%surface(i, ie)%v%val.ne.sdmiss)) then
           global_shared_buf(ie,2) = global_shared_buf(ie,2)+((obs(itime)%surface(i,ie)%v%val-bg_obs(itime)%surface(i,ie)%v%val)-obmsfv)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%t%val.ne.sdmiss).and.(bg_obs(itime)%surface(i, ie)%t%val.ne.sdmiss)) then
           global_shared_buf(ie,3) = global_shared_buf(ie,3)+((obs(itime)%surface(i,ie)%t%val-bg_obs(itime)%surface(i,ie)%t%val)-obmsft)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%q%val.ne.sdmiss).and.(bg_obs(itime)%surface(i, ie)%q%val.ne.sdmiss)) then
           global_shared_buf(ie,4) = global_shared_buf(ie,4)+((obs(itime)%surface(i,ie)%q%val-bg_obs(itime)%surface(i,ie)%q%val)-obmsfq)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%p%val.ne.sdmiss).and.(bg_obs(itime)%surface(i, ie)%p%val.ne.sdmiss)) then
           global_shared_buf(ie,5) = global_shared_buf(ie,5)+((obs(itime)%surface(i,ie)%p%val-bg_obs(itime)%surface(i,ie)%p%val)-obmsfps)**2.d0
         endif
       enddo
! analysis
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%u%val.ne.sdmiss).and.(an_obs(itime)%surface(i, ie)%u%val.ne.sdmiss)) then
           global_shared_buf(ie,6) = global_shared_buf(ie,6)+((obs(itime)%surface(i,ie)%u%val-an_obs(itime)%surface(i,ie)%u%val)-oamsfu)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%v%val.ne.sdmiss).and.(an_obs(itime)%surface(i, ie)%v%val.ne.sdmiss)) then
           global_shared_buf(ie,7) = global_shared_buf(ie,7)+((obs(itime)%surface(i,ie)%v%val-an_obs(itime)%surface(i,ie)%v%val)-oamsfv)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%t%val.ne.sdmiss).and.(an_obs(itime)%surface(i, ie)%t%val.ne.sdmiss)) then
           global_shared_buf(ie,8) = global_shared_buf(ie,8)+((obs(itime)%surface(i,ie)%t%val-an_obs(itime)%surface(i,ie)%t%val)-oamsft)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%q%val.ne.sdmiss).and.(an_obs(itime)%surface(i, ie)%q%val.ne.sdmiss)) then
           global_shared_buf(ie,9) = global_shared_buf(ie,9)+((obs(itime)%surface(i,ie)%q%val-an_obs(itime)%surface(i,ie)%q%val)-oamsfq)**2.d0
         endif
       enddo
       do i = 1, obsmet(itime)%nsf_e(ie)
         if ((obs(itime)%surface(i, ie)%p%val.ne.sdmiss).and.(an_obs(itime)%surface(i, ie)%p%val.ne.sdmiss)) then
           global_shared_buf(ie,10) = global_shared_buf(ie,10)+((obs(itime)%surface(i,ie)%p%val-an_obs(itime)%surface(i,ie)%p%val)-oamsfps)**2.d0
         endif
       enddo
     enddo ! ie
     call wrap_repro_sum(nvars = 10,comm = hybrid%par%comm)
! O-B
     obdsfu = sqrt(global_shared_sum(1)/nobsfu)
     obdsfv = sqrt(global_shared_sum(2)/nobsfv)
     obdsft = sqrt(global_shared_sum(3)/nobsft)
     obdsfq = sqrt(global_shared_sum(4)/nobsfq)
     obdsfps = sqrt(global_shared_sum(5)/nobsfps)
! O-A
     oadsfu = sqrt(global_shared_sum(6)/noasfu)
     oadsfv = sqrt(global_shared_sum(7)/noasfv)
     oadsft = sqrt(global_shared_sum(8)/noasft)
     oadsfq = sqrt(global_shared_sum(9)/noasfq)
     oadsfps = sqrt(global_shared_sum(10)/noasfps)
   endif ! da_surface
!.. Standard deviation of O-B and O-A for Bogus
   if (da_bogus) then
     global_shared_buf(:,1:2) = 0.d0
     do ie = nets,nete
! background
       do i = 1, obsmet(itime)%nbg_e(ie)
         if ((obs(itime)%bogus(i, ie)%p%val.ne.sdmiss).and.(bg_obs(itime)%bogus(i, ie)%p%val.ne.sdmiss)) then
           global_shared_buf(ie,1) = global_shared_buf(ie,1)+((obs(itime)%bogus(i,ie)%p%val-bg_obs(itime)%bogus(i,ie)%p%val)-obmbgps)**2.d0
         endif
       enddo
! analysis
       do i = 1, obsmet(itime)%nbg_e(ie)
         if ((obs(itime)%bogus(i, ie)%p%val.ne.sdmiss).and.(an_obs(itime)%bogus(i, ie)%p%val.ne.sdmiss)) then
           global_shared_buf(ie,2) = global_shared_buf(ie,2)+((obs(itime)%bogus(i,ie)%p%val-an_obs(itime)%bogus(i,ie)%p%val)-oambgps)**2.d0
         endif
       enddo ! da_bogus
     enddo ! ie
     call wrap_repro_sum(nvars = 2,comm = hybrid%par%comm)
! O-B
     obdbgps = sqrt(global_shared_sum(1)/nobbgps)
! O-A
     oadbgps = sqrt(global_shared_sum(2)/noabgps) ! o-a
   endif ! da_bogus
!.. Print out the standard deviation of O-B and O-A
   if (da_sonde) then
     if (par%ismasterproc) write(iulog,*) "standard deviation of o-b and o-a"
     if (par%ismasterproc) write(iulog,*) "against sonde observation"
     if (par%ismasterproc) write(iulog,*) "o-b u = ",obdsdu,"#obs = ",nobsdu
     if (par%ismasterproc) write(iulog,*) "o-a u = ",oadsdu,"#obs = ",noasdu
     if (par%ismasterproc) write(iulog,*) "o-b v = ",obdsdv,"#obs = ",nobsdv
     if (par%ismasterproc) write(iulog,*) "o-a v = ",oadsdv,"#obs = ",noasdv
     if (par%ismasterproc) write(iulog,*) "o-b t = ",obdsdt,"#obs = ",nobsdt
     if (par%ismasterproc) write(iulog,*) "o-a t = ",oadsdt,"#obs = ",noasdt
     if (par%ismasterproc) write(iulog,*) "o-b q = ",obdsdq,"#obs = ",nobsdq
     if (par%ismasterproc) write(iulog,*) "o-a q = ",oadsdq,"#obs = ",noasdq
   endif
   if (da_surface) then
     if (par%ismasterproc) write(iulog,*) "against surface observation"
     if (par%ismasterproc) write(iulog,*) "o-b u = ",obdsfu,"#obs = ",nobsfu
     if (par%ismasterproc) write(iulog,*) "o-a u = ",oadsfu,"#obs = ",noasfu
     if (par%ismasterproc) write(iulog,*) "o-b v = ",obdsfv,"#obs = ",nobsfv
     if (par%ismasterproc) write(iulog,*) "o-a v = ",oadsfv,"#obs = ",noasfv
     if (par%ismasterproc) write(iulog,*) "o-b t = ",obdsft,"#obs = ",nobsft
     if (par%ismasterproc) write(iulog,*) "o-a t = ",oadsft,"#obs = ",noasft
     if (par%ismasterproc) write(iulog,*) "o-b q = ",obdsfq,"#obs = ",nobsfq
     if (par%ismasterproc) write(iulog,*) "o-a q = ",oadsfq,"#obs = ",noasfq
     if (par%ismasterproc) write(iulog,*) "o-b ps = ",obdsfps,"#obs = ",nobsfps
     if (par%ismasterproc) write(iulog,*) "o-a ps = ",oadsfps,"#obs = ",noasfps
   endif
   if (da_bogus) then
     if (par%ismasterproc) write(iulog,*) "against surface bogus observation"
     if (par%ismasterproc) write(iulog,*) "o-b ps = ",obdbgps,"#obs = ",nobbgps
     if (par%ismasterproc) write(iulog,*) "o-a ps = ",oadbgps,"#obs = ",noabgps
   endif
!    DO k = 0,par%nprocs-1
!    IF( par%iproc == k )THEN
!    DO ie= nets,nete
!    DO i = 1,nsf_e(ie)
!       IF( obs%surface(i,ie)%u%val .NE. sdmiss )THEN
!          WRITE(iulog,'(3I5,6F10.4)') par%iproc,i,ie,  &
!                         obsloc%surface%all(ie)%lat(i),&
!                         obsloc%surface%all(ie)%lon(i),&
!                         obsloc%surface%all(ie)%hgt(i),&
!                         obs%surface(i,ie)%u%val,      &
!                         obs%surface(i,ie)%u%val        &
!                         -bg_obs%surface(i,ie)%u%val,  &
!                         obs%surface(i,ie)%u%val        &
!                         -an_obs%surface(i,ie)%u%val
!       ENDIF
!    ENDDO
!    ENDDO
!    ENDIF
!       CALL Mpi_Barrier(par%comm,ierr)
!    ENDDO
!
!    DO k = 0,par%nprocs-1
!    IF( par%iproc == k )THEN
!    DO ie= nets,nete
!    DO i = 1,nsf_e(ie)
!       IF( obs%surface(i,ie)%t%val .NE. sdmiss )THEN
!          WRITE(iulog,'(3I5,6F10.4)') par%iproc,i,ie,  &
!                         obsloc%surface%all(ie)%lat(i),&
!                         obsloc%surface%all(ie)%lon(i),&
!                         obsloc%surface%all(ie)%hgt(i),&
!                         obs%surface(i,ie)%t%val,      &
!                         obs%surface(i,ie)%t%val        &
!                         -bg_obs%surface(i,ie)%t%val,  &
!                         obs%surface(i,ie)%t%val        &
!                         -an_obs%surface(i,ie)%t%val
!       ENDIF
!    ENDDO
!    ENDDO
!    ENDIF
!       CALL Mpi_Barrier(par%comm,ierr)
!    ENDDO
!
!    DO k = 0,par%nprocs-1
!    IF( par%iproc == k )THEN
!    DO ie= nets,nete
!    DO i = 1,nsf_e(ie)
!       IF( obs%surface(i,ie)%p%val .NE. sdmiss )THEN
!          WRITE(iulog,'(3I5,3F10.4,3F10.2)') par%iproc,i,ie,&
!                         obsloc%surface%all(ie)%lat(i),     &
!                         obsloc%surface%all(ie)%lon(i),     &
!                         obsloc%surface%all(ie)%hgt(i),     &
!                         obs%surface(i,ie)%p%val,           &
!                         obs%surface(i,ie)%p%val             &
!                         -bg_obs%surface(i,ie)%p%val,       &
!                         obs%surface(i,ie)%p%val             &
!                         -an_obs%surface(i,ie)%p%val
!       ENDIF
!    ENDDO
!    ENDDO
!    ENDIF
!       CALL Mpi_Barrier(par%comm,ierr)
!    ENDDO
!       CALL Mpi_Barrier(par%comm,ierr)
!
   return
!
   end subroutine check_ob
!===============================================================================
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine check_adj(gs_md, b_eig, b_eig_chi, b_eig_tps, b_eig_q, b_eig_a)
! 2. Variables
! 2-1. In/Output
!
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(in   ) :: gs_md
   real(real_kind), dimension(neig, zwn_nproc), intent(in   ) :: b_eig
   real(real_kind), dimension(neig_chi, zwn_nproc), intent(in   ) :: b_eig_chi
   real(real_kind), dimension(neig_tps, zwn_nproc), intent(in   ) :: b_eig_tps
   real(real_kind), dimension(neig_q, zwn_nproc), intent(in   ) :: b_eig_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl), intent(in   ) :: b_eig_a
! 2-2. Local
!   REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
!      :: tl_md
   real(real_kind), dimension(neig, zwn_nproc) :: ad_eig
   real(real_kind), dimension(neig_chi, zwn_nproc) :: ad_eig_chi
   real(real_kind), dimension(neig_tps, zwn_nproc) :: ad_eig_tps
   real(real_kind), dimension(neig_q, zwn_nproc) :: ad_eig_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl) :: ad_eig_a
   real(real_kind) :: tmp, yty
   integer(int_kind) :: i, j, k, ie, wn, vr, itime
! 1. Eigen space -> observation space
!
   ierr = timingstart('TlmEigToObs_Check_Adj') !---profiling
   call tlmeigtoobs(gs_md(:,:,:,:,:),b_eig,b_eig_chi,b_eig_tps,b_eig_q,b_eig_a,tl_obs(:))
   ierr = timingstop('TlmEigToObs_Check_Adj') !---profiling
   do itime = 1,ntime
! single observation
     if (singleobs) then
       call multiplyrsingle(itime,tl_obs(itime)%single%u%val,tl_obs(itime)%single%v%val,tl_obs(itime)%single%t%val,tl_obs(itime)%single%q%val,tl_obs(itime)%single%ps%val)
     endif
! sonde observation
     if (da_sonde) then
       call multiplyrsonde(itime,tl_obs(itime)%sonde%u%val,tl_obs(itime)%sonde%v%val,tl_obs(itime)%sonde%t%val,tl_obs(itime)%sonde%q%val)
     endif
! surface observation
     if (da_surface) then
       call multiplyrsurface(itime,tl_obs(itime)%surface%u%val,tl_obs(itime)%surface%v%val,tl_obs(itime)%surface%t%val,tl_obs(itime)%surface%q%val,tl_obs(itime)%surface%p%val)
     endif
! bogus observation
     if (da_bogus) then
       call multiplyrbogus(itime,tl_obs(itime)%bogus%p%val)
     endif
! aircraft observation
     if (da_aircraft) then
       call multiplyraircraft(itime,tl_obs(itime)%aircraft%u%val,tl_obs(itime)%aircraft%v%val,tl_obs(itime)%aircraft%t%val)
     endif
! amv observation
     if (da_amv) then
       call multiplyramv(itime,tl_obs(itime)%amv%u%val,tl_obs(itime)%amv%v%val)
     endif
! scatwind observation
     if (da_scatwind) then
       call multiplyrscatwind(itime,tl_obs(itime)%scatwind%u%val,tl_obs(itime)%scatwind%v%val)
     endif
! amsua observation
     if (da_amsua) then
       call multiplyramsua(itime,tl_obs(itime)%amsua%t%val)
     endif
! gpsro observation
     if (da_gpsro) then
       call multiplyrgpsro(itime,tl_obs(itime)%gpsro%t%val,tl_obs(itime)%gpsro%q%val)
     endif
! iasi observation
     if (da_iasi) then
       call multiplyriasi(itime,tl_obs(itime)%iasi%t%val,tl_obs(itime)%iasi%q%val)
     endif
! cris observation
     if (da_cris) then
       call multiplyrcris(itime,tl_obs(itime)%cris%t%val,tl_obs(itime)%cris%q%val)
     endif
! atms observation
     if (da_atms) then
       call multiplyratms(itime,tl_obs(itime)%atms%t%val)
     endif
! atmswv observation
     if (da_atmswv) then
       call multiplyratmswv(itime,tl_obs(itime)%atmswv%t%val,tl_obs(itime)%atmswv%q%val)
     endif
! mhs observation
     if (da_mhs) then
       call multiplyrmhs(itime,tl_obs(itime)%mhs%t%val,tl_obs(itime)%mhs%q%val)
     endif
! csr observation
     if (da_csr) then
       call multiplyrcsr(itime,tl_obs(itime)%csr%t%val,tl_obs(itime)%csr%q%val)
     endif
! 2. Adjoint of eigen space -> observation space
! Single observation
     if (singleobs) then
       ad_obs(itime)%single = tl_obs(itime)%single
       call multiplyrsingle(itime,ad_obs(itime)%single%u%val,ad_obs(itime)%single%v%val,ad_obs(itime)%single%t%val,ad_obs(itime)%single%q%val,ad_obs(itime)%single%ps%val)
     endif
! sonde observation
     if (da_sonde) then
       ad_obs(itime)%sonde = tl_obs(itime)%sonde
       call multiplyrsonde(itime,ad_obs(itime)%sonde%u%val,ad_obs(itime)%sonde%v%val,ad_obs(itime)%sonde%t%val,ad_obs(itime)%sonde%q%val)
     endif
! surface observation
     if (da_surface) then
       ad_obs(itime)%surface = tl_obs(itime)%surface
       call multiplyrsurface(itime,ad_obs(itime)%surface%u%val,ad_obs(itime)%surface%v%val,ad_obs(itime)%surface%t%val,ad_obs(itime)%surface%q%val,ad_obs(itime)%surface%p%val)
     endif
! bogus observation
     if (da_bogus) then
       ad_obs(itime)%bogus = tl_obs(itime)%bogus
       call multiplyrbogus(itime,ad_obs(itime)%bogus%p%val)
     endif
! aircraft observation
     if (da_aircraft) then
       ad_obs(itime)%aircraft = tl_obs(itime)%aircraft
       call multiplyraircraft(itime,ad_obs(itime)%aircraft%u%val,ad_obs(itime)%aircraft%v%val,ad_obs(itime)%aircraft%t%val)
     endif
! amv observation
     if (da_amv) then
       ad_obs(itime)%amv = tl_obs(itime)%amv
       call multiplyramv(itime,ad_obs(itime)%amv%u%val,ad_obs(itime)%amv%v%val)
     endif
! scatwind observation
     if (da_scatwind) then
       ad_obs(itime)%scatwind = tl_obs(itime)%scatwind
       call multiplyrscatwind(itime,ad_obs(itime)%scatwind%u%val,ad_obs(itime)%scatwind%v%val)
     endif
! amsua observation
     if (da_amsua) then
       ad_obs(itime)%amsua = tl_obs(itime)%amsua
       call multiplyramsua(itime,ad_obs(itime)%amsua%t%val)
     endif
! gpsro observation
     if (da_gpsro) then
       ad_obs(itime)%gpsro = tl_obs(itime)%gpsro
       call multiplyrgpsro(itime,ad_obs(itime)%gpsro%t%val,ad_obs(itime)%gpsro%q%val)
     endif
! iasi observation
     if (da_iasi) then
       ad_obs(itime)%iasi = tl_obs(itime)%iasi
       call multiplyriasi(itime,ad_obs(itime)%iasi%t%val,ad_obs(itime)%iasi%q%val)
     endif
! cris observation
     if (da_cris) then
       ad_obs(itime)%cris = tl_obs(itime)%cris
       call multiplyrcris(itime,ad_obs(itime)%cris%t%val,ad_obs(itime)%cris%q%val)
     endif
! atms observation
     if (da_atms) then
       ad_obs(itime)%atms = tl_obs(itime)%atms
       call multiplyratms(itime,ad_obs(itime)%atms%t%val)
     endif
! atmswv observation
     if (da_atmswv) then
       ad_obs(itime)%atmswv = tl_obs(itime)%atmswv
       call multiplyratmswv(itime,ad_obs(itime)%atmswv%t%val,ad_obs(itime)%atmswv%q%val)
     endif
! mhs observation
     if (da_mhs) then
       ad_obs(itime)%mhs = tl_obs(itime)%mhs
       call multiplyrmhs(itime,ad_obs(itime)%mhs%t%val,ad_obs(itime)%mhs%q%val)
     endif
! csr observation
     if (da_csr) then
       ad_obs(itime)%csr = tl_obs(itime)%csr
       call multiplyrcsr(itime,ad_obs(itime)%csr%t%val,ad_obs(itime)%csr%q%val)
     endif
   enddo ! itime
   ierr = timingstart('AdjEigToObs_Check_Adj') !---profiling
   ad_eig = 0.d0
   ad_eig_chi = 0.d0
   ad_eig_tps = 0.d0
   ad_eig_q = 0.d0
   ad_eig_a = 0.d0
   call adjeigtoobs(gs_md(:,:,:,:,:),ad_obs(:),ad_eig,ad_eig_chi,ad_eig_tps,ad_eig_q,ad_eig_a)
   ierr = timingstop('AdjEigToObs_Check_Adj') !---profiling
! 3. Evaluation
   global_shared_buf(:,1) = 0.d0
   do itime = 1,ntime
! single observation test
     if (singleobs) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nsg_e(ie)
           if (tl_obs(itime)%single(i, ie)%u%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%single(i,ie)%u%val**2.d0
           endif
           if (tl_obs(itime)%single(i, ie)%v%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%single(i,ie)%v%val**2.d0
           endif
           if (tl_obs(itime)%single(i, ie)%t%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%single(i,ie)%t%val**2.d0
           endif
           if (tl_obs(itime)%single(i, ie)%q%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%single(i,ie)%q%val**2.d0
           endif
           if (tl_obs(itime)%single(i, ie)%ps%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%single(i,ie)%ps%val**2.d0
           endif
         enddo
       enddo
     endif
! sonde observation
     if (da_sonde) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nsd_e(1, ie)
           if (tl_obs(itime)%sonde(i, ie)%u%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%sonde(i,ie)%u%val**2.d0
           endif
         enddo
         do i = 1, obsmet(itime)%nsd_e(2, ie)
           if (tl_obs(itime)%sonde(i, ie)%v%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%sonde(i,ie)%v%val**2.d0
           endif
         enddo
         do i = 1, obsmet(itime)%nsd_e(3, ie)
           if (tl_obs(itime)%sonde(i, ie)%t%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%sonde(i,ie)%t%val**2.d0
           endif
         enddo
         do i = 1, obsmet(itime)%nsd_e(4, ie)
           if (tl_obs(itime)%sonde(i, ie)%q%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%sonde(i,ie)%q%val**2.d0
           endif
         enddo
       enddo
     endif
! surface observation
     if (da_surface) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nsf_e(ie)
           if (tl_obs(itime)%surface(i, ie)%u%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%surface(i,ie)%u%val**2.d0
           endif
           if (tl_obs(itime)%surface(i, ie)%v%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%surface(i,ie)%v%val**2.d0
           endif
           if (tl_obs(itime)%surface(i, ie)%t%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%surface(i,ie)%t%val**2.d0
           endif
           if (tl_obs(itime)%surface(i, ie)%q%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%surface(i,ie)%q%val**2.d0
           endif
           if (tl_obs(itime)%surface(i, ie)%p%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%surface(i,ie)%p%val**2.d0
           endif
         enddo
       enddo
     endif
! bogus observation
     if (da_bogus) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nbg_e(ie)
           if (tl_obs(itime)%bogus(i, ie)%p%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%bogus(i,ie)%p%val**2.d0
           endif
         enddo
       enddo
     endif
! aircraft observation
     if (da_aircraft) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nar_e(1, ie)
           if (tl_obs(itime)%aircraft(i, ie)%u%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%aircraft(i,ie)%u%val**2.d0
           endif
         enddo
         do i = 1, obsmet(itime)%nar_e(2, ie)
           if (tl_obs(itime)%aircraft(i, ie)%v%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%aircraft(i,ie)%v%val**2.d0
           endif
         enddo
         do i = 1, obsmet(itime)%nar_e(3, ie)
           if (tl_obs(itime)%aircraft(i, ie)%t%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%aircraft(i,ie)%t%val**2.d0
           endif
         enddo
       enddo
     endif
! amv observation
     if (da_amv) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%namv_e(ie)
           if (tl_obs(itime)%amv(i, ie)%u%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%amv(i,ie)%u%val**2.d0
           endif
           if (tl_obs(itime)%amv(i, ie)%v%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%amv(i,ie)%v%val**2.d0
           endif
         enddo
       enddo
     endif
! scatwind observation
     if (da_scatwind) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nscatwind_e(ie)
           if (tl_obs(itime)%scatwind(i, ie)%u%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%scatwind(i,ie)%u%val**2.d0
           endif
           if (tl_obs(itime)%scatwind(i, ie)%v%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%scatwind(i,ie)%v%val**2.d0
           endif
         enddo
       enddo
     endif
! amsua observation
     if (da_amsua) then
       do j = 1, obsmet(itime)%namsuach
         do ie = nets, nete
           do i = 1, obsmet(itime)%namsua_e(ie)
             if (tl_obs(itime)%amsua(i, j, ie)%t%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%amsua(i,j,ie)%t%val**2.d0
             endif
           enddo
         enddo
       enddo
     endif
! gpsro observation
     if (da_gpsro) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%ngpsro_e(ie)
           if (tl_obs(itime)%gpsro(i, ie)%t%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%gpsro(i,ie)%t%val**2.d0
           endif
           if (tl_obs(itime)%gpsro(i, ie)%q%val.ne.sdmiss) then
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%gpsro(i,ie)%q%val**2.d0
           endif
         enddo
       enddo
     endif
! iasi observation
     if (da_iasi) then
       do j = 1, obsmet(itime)%niasich
         do ie = nets, nete
           do i = 1, obsmet(itime)%niasi_e(ie)
             if (tl_obs(itime)%iasi(i, j, ie)%t%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%iasi(i,j,ie)%t%val**2.d0
             endif
             if (tl_obs(itime)%iasi(i, j, ie)%q%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%iasi(i,j,ie)%q%val**2.d0
             endif
           enddo
         enddo
       enddo
     endif
! cris observation
     if (da_cris) then
       do j = 1, obsmet(itime)%ncrisch
         do ie = nets, nete
           do i = 1, obsmet(itime)%ncris_e(ie)
             if (tl_obs(itime)%cris(i, j, ie)%t%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%cris(i,j,ie)%t%val**2.d0
             endif
             if (tl_obs(itime)%cris(i, j, ie)%q%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%cris(i,j,ie)%q%val**2.d0
             endif
           enddo
         enddo
       enddo
     endif
! atms observation
     if (da_atms) then
       do j = 1, obsmet(itime)%natmsch
         do ie = nets, nete
           do i = 1, obsmet(itime)%natms_e(ie)
             if (tl_obs(itime)%atms(i, j, ie)%t%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%atms(i,j,ie)%t%val**2.d0
             endif
           enddo
         enddo
       enddo
     endif
! atmswv observation
     if (da_atmswv) then
       do j = 1, obsmet(itime)%natmswvch
         do ie = nets, nete
           do i = 1, obsmet(itime)%natmswv_e(ie)
             if (tl_obs(itime)%atmswv(i, j, ie)%t%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%atmswv(i,j,ie)%t%val**2.d0
             endif
             if (tl_obs(itime)%atmswv(i, j, ie)%q%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%atmswv(i,j,ie)%q%val**2.d0
             endif
           enddo
         enddo
       enddo
     endif
! mhs observation
     if (da_mhs) then
       do j = 1, obsmet(itime)%nmhsch
         do ie = nets, nete
           do i = 1, obsmet(itime)%nmhs_e(ie)
             if (tl_obs(itime)%mhs(i, j, ie)%t%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%mhs(i,j,ie)%t%val**2.d0
             endif
             if (tl_obs(itime)%mhs(i, j, ie)%q%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%mhs(i,j,ie)%q%val**2.d0
             endif
           enddo
         enddo
       enddo
     endif
! csr observation
     if (da_csr) then
       do j = 1, obsmet(itime)%ncsrch
         do ie = nets, nete
           do i = 1, obsmet(itime)%ncsr_e(ie)
             if (tl_obs(itime)%csr(i, j, ie)%t%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%csr(i,j,ie)%t%val**2.d0
             endif
             if (tl_obs(itime)%csr(i, j, ie)%q%val.ne.sdmiss) then
               global_shared_buf(ie,1) = global_shared_buf(ie,1)+tl_obs(itime)%csr(i,j,ie)%q%val**2.d0
             endif
           enddo
         enddo
       enddo
     endif
   enddo ! itime
   call wrap_repro_sum(nvars = 1,comm = hybrid%par%comm)
   if (par%ismasterproc) write(iulog,*) ' '
   if (par%ismasterproc) write(iulog,*) 'Adjoint test '
   if (par%ismasterproc) then
     yty = global_shared_sum(1)
     write(iulog,*) '<hx,hx>= ',yty
   endif
   call innerproduct(neig,neig_chi,neig_tps,neig_q,zwn_nproc,zwn_nproc_a,b_eig,b_eig_chi,b_eig_tps,b_eig_q,b_eig_a,ad_eig,ad_eig_chi,ad_eig_tps,ad_eig_q,ad_eig_a,tmp)
   if (par%ismasterproc) then
     write(iulog,*) '<h^t hx,x>= ',tmp
     write(iulog,*) '<hx,hx>-<h^t hx,x>= ',abs(yty-tmp)
   endif
   return
!
   end subroutine check_adj
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
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine innerproduct(neig, neig_chi, neig_tps, neig_q, zwn_nproc, zwn_nproc_a, vect1, vect1_chi, vect1_tps, vect1_q, vect1_a, vect2, vect2_chi, vect2_tps, vect2_q, vect2_a, ddot_product)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
! 2. Variables
! 2-1. Input
!
   integer(int_kind), intent(in   ) :: neig, neig_chi, neig_tps, neig_q, zwn_nproc, zwn_nproc_a
   real(real_kind), dimension(neig, zwn_nproc), intent(in   ) :: vect1, vect2
   real(real_kind), dimension(neig_chi, zwn_nproc), intent(in   ) :: vect1_chi, vect2_chi
   real(real_kind), dimension(neig_tps, zwn_nproc), intent(in   ) :: vect1_tps, vect2_tps
   real(real_kind), dimension(neig_q, zwn_nproc), intent(in   ) :: vect1_q, vect2_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl), intent(in   ) :: vect1_a, vect2_a !---for alpha(hjs)
! 2-2. Output
   real(real_kind), intent(  out) :: ddot_product
! 2-3. Local
   real(real_kind) :: tmp, tmp1
   integer(int_kind) :: i, j, k, ie, wn, jump, ismpl ! 4denvar
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!
   tmp1 = 0.d0
   ddot_product = 0.d0
   jump = 5
   k = mod(zwn_nproc,jump)
   j = k+1
   if (k.ne.0) then
     do wn = 1, k
       tmp1 = tmp1+dot_product(vect1(:,wn),vect2(:,wn))+dot_product(vect1_chi(:,wn),vect2_chi(:,wn))+dot_product(vect1_tps(:,wn),vect2_tps(:,wn))+dot_product(vect1_q(:,wn),vect2_q(:,wn))
     enddo
!IF (par%ismasterproc) WRITE(iulog,*) 'tmp 1= ', tmp1,k,zwn_nproc,jump
   endif
   do wn = j,zwn_nproc,jump
     tmp1 = tmp1+dot_product(vect1(:,wn),vect2(:,wn))+dot_product(vect1(:,wn+1),vect2(:,wn+1))+dot_product(vect1(:,wn+2),vect2(:,wn+2))+dot_product(vect1(:,wn+3),vect2(:,wn+3))+dot_product(vect1(:,wn+4),vect2(:,wn+4))+dot_product(vect1_chi(:,wn),vect2_chi(:,wn))+dot_product(vect1_chi(:,wn+1),vect2_chi(:,wn+1))+dot_product(vect1_chi(:,wn+2),vect2_chi(:,wn+2))+dot_product(vect1_chi(:,wn+3),vect2_chi(:,wn+3))+dot_product(vect1_chi(:,wn+4),vect2_chi(:,wn+4))+dot_product(vect1_tps(:,wn),vect2_tps(:,wn))+dot_product(vect1_tps(:,wn+1),vect2_tps(:,wn+1))+dot_product(vect1_tps(:,wn+2),vect2_tps(:,wn+2))+dot_product(vect1_tps(:,wn+3),vect2_tps(:,wn+3))+dot_product(vect1_tps(:,wn+4),vect2_tps(:,wn+4))+dot_product(vect1_q(:,wn),vect2_q(:,wn))+dot_product(vect1_q(:,wn+1),vect2_q(:,wn+1))+dot_product(vect1_q(:,wn+2),vect2_q(:,wn+2))+dot_product(vect1_q(:,wn+3),vect2_q(:,wn+3))+dot_product(vect1_q(:,wn+4),vect2_q(:,wn+4))
   enddo
   k = mod(zwn_nproc_a,jump)
   j = k+1
   do ismpl = 1,nsmpl
     if (k.ne.0) then
       do wn = 1, k
         tmp1 = tmp1+dot_product(vect1_a(:,wn,ismpl),vect2_a(:,wn,ismpl))
       enddo
!IF (par%ismasterproc) WRITE(iulog,*) 'tmp 1= ', tmp1,k,zwn_nproc,jump
     endif
     do wn = j, zwn_nproc_a, jump
       tmp1 = tmp1+dot_product(vect1_a(:,wn,ismpl),vect2_a(:,wn,ismpl))+dot_product(vect1_a(:,wn+1,ismpl),vect2_a(:,wn+1,ismpl))+dot_product(vect1_a(:,wn+2,ismpl),vect2_a(:,wn+2,ismpl))+dot_product(vect1_a(:,wn+3,ismpl),vect2_a(:,wn+3,ismpl))+dot_product(vect1_a(:,wn+4,ismpl),vect2_a(:,wn+4,ismpl))
     enddo
   enddo ! ismpl
   global_shared_buf(:,1) = 0.d0
   global_shared_buf(nets:nete,1) = tmp1/dble(nete-nets+1)
   call wrap_repro_sum(nvars = 1,comm = hybrid%par%comm)
   ddot_product = global_shared_sum(1)
!IF (par%ismasterproc) WRITE(iulog,*) 'tmp 1= ',tmp1,global_shared_sum(1)
   return
!
   end subroutine innerproduct
!===============================================================================
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine make_lhs_ax(gs_md, x_eig, x_eig_chi, x_eig_tps, x_eig_q, x_eig_a, ap_eig, ap_eig_chi, ap_eig_tps, ap_eig_q, ap_eig_a)
! In/Output
!
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(in   ) :: gs_md
   real(real_kind), dimension(neig, zwn_nproc), intent(in   ) :: x_eig
   real(real_kind), dimension(neig_chi, zwn_nproc), intent(in   ) :: x_eig_chi
   real(real_kind), dimension(neig_tps, zwn_nproc), intent(in   ) :: x_eig_tps
   real(real_kind), dimension(neig_q, zwn_nproc), intent(in   ) :: x_eig_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl), intent(in   ) :: x_eig_a
   real(real_kind), dimension(neig, zwn_nproc), intent(  out) :: ap_eig
   real(real_kind), dimension(neig_chi, zwn_nproc), intent(  out) :: ap_eig_chi
   real(real_kind), dimension(neig_tps, zwn_nproc), intent(  out) :: ap_eig_tps
   real(real_kind), dimension(neig_q, zwn_nproc), intent(  out) :: ap_eig_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl), intent(  out) :: ap_eig_a
! Local
   real(real_kind), dimension(neig, zwn_nproc) :: ad_eig
   real(real_kind), dimension(neig_chi, zwn_nproc) :: ad_eig_chi
   real(real_kind), dimension(neig_tps, zwn_nproc) :: ad_eig_tps
   real(real_kind), dimension(neig_q, zwn_nproc) :: ad_eig_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl) :: ad_eig_a
   real(real_kind) :: tmp
   integer(int_kind) :: i, j, k, ie, wn, itime, ismpl
! 1. Eigen space -> observation space
!
   ierr = timingstart('TlmEigToObs_Make_Ax') !---profiling
   call tlmeigtoobs(gs_md,x_eig,x_eig_chi,x_eig_tps,x_eig_q,x_eig_a,tl_obs(:))
   ierr = timingstop('TlmEigToObs_Make_Ax') !---profiling
   do itime = 1,ntime
! 2. Multiply with R^-1
! 2.0. single observation test
     if (singleobs) then
       call multiplyrsingle(itime,tl_obs(itime)%single%u%val,tl_obs(itime)%single%v%val,tl_obs(itime)%single%t%val,tl_obs(itime)%single%q%val,tl_obs(itime)%single%ps%val)
     endif
! 2.1. sonde observation
     if (da_sonde) then
       call multiplyrsonde(itime,tl_obs(itime)%sonde%u%val,tl_obs(itime)%sonde%v%val,tl_obs(itime)%sonde%t%val,tl_obs(itime)%sonde%q%val)
     endif
! 2.2. surface observation
     if (da_surface) then
       call multiplyrsurface(itime,tl_obs(itime)%surface%u%val,tl_obs(itime)%surface%v%val,tl_obs(itime)%surface%t%val,tl_obs(itime)%surface%q%val,tl_obs(itime)%surface%p%val)
     endif
! 2.3. bogus observation
     if (da_bogus) then
       call multiplyrbogus(itime,tl_obs(itime)%bogus%p%val)
     endif
! 2.3. aircraft observation
     if (da_aircraft) then
       call multiplyraircraft(itime,tl_obs(itime)%aircraft%u%val,tl_obs(itime)%aircraft%v%val,tl_obs(itime)%aircraft%t%val)
     endif
! 2.4. amv observation
     if (da_amv) then
       call multiplyramv(itime,tl_obs(itime)%amv%u%val,tl_obs(itime)%amv%v%val)
     endif
! 2.5. scatwind observation
     if (da_scatwind) then
       call multiplyrscatwind(itime,tl_obs(itime)%scatwind%u%val,tl_obs(itime)%scatwind%v%val)
     endif
! 2.6. amsua observation
     if (da_amsua) then
       call multiplyramsua(itime,tl_obs(itime)%amsua%t%val)
     endif
! 2.7. gpsro observation
     if (da_gpsro) then
       call multiplyrgpsro(itime,tl_obs(itime)%gpsro%t%val,tl_obs(itime)%gpsro%q%val)
     endif
! 2.8. iasi observation
     if (da_iasi) then
       call multiplyriasi(itime,tl_obs(itime)%iasi%t%val,tl_obs(itime)%iasi%q%val)
     endif
! 2.9. cris observation
     if (da_cris) then
       call multiplyrcris(itime,tl_obs(itime)%cris%t%val,tl_obs(itime)%cris%q%val)
     endif
! 2.10. atms observation
     if (da_atms) then
       call multiplyratms(itime,tl_obs(itime)%atms%t%val)
     endif
! 2.11. atmswv observation
     if (da_atmswv) then
       call multiplyratmswv(itime,tl_obs(itime)%atmswv%t%val,tl_obs(itime)%atmswv%q%val)
     endif
! 2.12. mhs observation
     if (da_mhs) then
       call multiplyrmhs(itime,tl_obs(itime)%mhs%t%val,tl_obs(itime)%mhs%q%val)
     endif
! 2.13. csr observation
     if (da_csr) then
       call multiplyrcsr(itime,tl_obs(itime)%csr%t%val,tl_obs(itime)%csr%q%val)
     endif
   enddo ! itime
! 3. Adjoint of spectral space -> observation space
   ierr = timingstart('AdjEigToObs_Make_Ax') !---profiling
   ad_eig = 0.d0
   ad_eig_chi = 0.d0
   ad_eig_tps = 0.d0
   ad_eig_q = 0.d0
   ad_eig_a = 0.d0
   call adjeigtoobs(gs_md,tl_obs(:),ad_eig,ad_eig_chi,ad_eig_tps,ad_eig_q,ad_eig_a)
   ierr = timingstop('AdjEigToObs_Make_Ax') !---profiling
! 4. Finalize calculation of Ax
   do wn = 1,zwn_nproc
     ap_eig(:,wn) = x_eig(:,wn)+ad_eig(:,wn)
     ap_eig_chi(:,wn) = x_eig_chi(:,wn)+ad_eig_chi(:,wn)
     ap_eig_tps(:,wn) = x_eig_tps(:,wn)+ad_eig_tps(:,wn)
     ap_eig_q(:,wn) = x_eig_q(:,wn)+ad_eig_q(:,wn)
   enddo
   do ismpl = 1,nsmpl
     do wn = 1, zwn_nproc_a
       ap_eig_a(:,wn,ismpl) = x_eig_a(:,wn,ismpl)+ad_eig_a(:,wn,ismpl)
     enddo
   enddo ! ismpl
   return
!
   end subroutine make_lhs_ax
!===============================================================================
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine make_rhs_b(gs_md, b_eig, b_eig_chi, b_eig_tps, b_eig_q, b_eig_a)
! 2. Variables
! 2-1. In/Output
!
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(in   ) :: gs_md
   real(real_kind), dimension(neig, zwn_nproc), intent(  out) :: b_eig
   real(real_kind), dimension(neig_chi, zwn_nproc), intent(  out) :: b_eig_chi
   real(real_kind), dimension(neig_tps, zwn_nproc), intent(  out) :: b_eig_tps
   real(real_kind), dimension(neig_q, zwn_nproc), intent(  out) :: b_eig_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl), intent(  out) :: b_eig_a
! 2-2. Local
   real(real_kind), dimension(np, np, nlevar, nelemd) :: tl_md
   real(real_kind) :: tmp
   integer(int_kind) :: i, j, k, ie, wn, vr, itime
!
   do itime = 1, ntime
!     IF (par%ismasterproc) WRITE(iulog,*) ' In1 ModelToObs_Make_RHS_b'
     ierr = timingstart('ModelToObs_Make_RHS_b') !---profiling
     call modeltoobs(itime,gs_md(:,:,:,:,itime),gs_obs(itime))
     ierr = timingstop('ModelToObs_Make_RHS_b')
!    IF (par%ismasterproc) WRITE(iulog,*) 'Out1 ModelToObs_Make_RHS_b'
     tl_md = 0.d0
     do ie = 1,nelemd
       do k = 1, nlevar
         do j = 1, np
           do i = 1, np
             tl_md(i,j,k,ie) = bg_md(i,j,k,ie,itime)-gs_md(i,j,k,ie,itime)
           enddo
         enddo
       enddo
     enddo
     ierr = timingstart('TlmModelToObs_Make_RHS_b') !---profiling
     call tlmmodeltoobs(itime,gs_md(:,:,:,:,itime),tl_md(:,:,:,:),tl_obs(itime))
     ierr = timingstop('TlmModelToObs_Make_RHS_b')
! single observation test
     if (singleobs) then
       call makebsingle(itime,gs_obs(itime)%single%u%val,gs_obs(itime)%single%v%val,gs_obs(itime)%single%t%val,gs_obs(itime)%single%q%val,gs_obs(itime)%single%ps%val,tl_obs(itime)%single%u%val,tl_obs(itime)%single%v%val,tl_obs(itime)%single%t%val,tl_obs(itime)%single%q%val,tl_obs(itime)%single%ps%val)
     endif
! sonde observation
     if (da_sonde) then
       call makebsonde(itime,gs_obs(itime)%sonde%u%val,gs_obs(itime)%sonde%v%val,gs_obs(itime)%sonde%t%val,gs_obs(itime)%sonde%q%val,tl_obs(itime)%sonde%u%val,tl_obs(itime)%sonde%v%val,tl_obs(itime)%sonde%t%val,tl_obs(itime)%sonde%q%val)
     endif
! surface observation
     if (da_surface) then
       call makebsurface(itime,gs_obs(itime)%surface%u%val,gs_obs(itime)%surface%v%val,gs_obs(itime)%surface%t%val,gs_obs(itime)%surface%q%val,gs_obs(itime)%surface%p%val,tl_obs(itime)%surface%u%val,tl_obs(itime)%surface%v%val,tl_obs(itime)%surface%t%val,tl_obs(itime)%surface%q%val,tl_obs(itime)%surface%p%val)
     endif
! bogus observation
     if (da_bogus) then
       call makebbogus(itime,gs_obs(itime)%bogus%p%val,tl_obs(itime)%bogus%p%val)
     endif
! aircraft observation
     if (da_aircraft) then
       call makebaircraft(itime,gs_obs(itime)%aircraft%u%val,gs_obs(itime)%aircraft%v%val,gs_obs(itime)%aircraft%t%val,tl_obs(itime)%aircraft%u%val,tl_obs(itime)%aircraft%v%val,tl_obs(itime)%aircraft%t%val)
     endif
! amv observation
     if (da_amv) then
       call makebamv(itime,gs_obs(itime)%amv%u%val,gs_obs(itime)%amv%v%val,tl_obs(itime)%amv%u%val,tl_obs(itime)%amv%v%val)
     endif
! scatwind observation
     if (da_scatwind) then
       call makebscatwind(itime,gs_obs(itime)%scatwind%u%val,gs_obs(itime)%scatwind%v%val,tl_obs(itime)%scatwind%u%val,tl_obs(itime)%scatwind%v%val)
     endif
! amsua observation
     if (da_amsua) then
       call makebamsua(itime,gs_obs(itime)%amsua%t%val,tl_obs(itime)%amsua%t%val)
     endif
! gpsro observation
     if (da_gpsro) then
       call makebgpsro(itime,gs_obs(itime)%gpsro%t%val,gs_obs(itime)%gpsro%q%val,tl_obs(itime)%gpsro%t%val,tl_obs(itime)%gpsro%q%val)
     endif
! iasi observation
     if (da_iasi) then
       call makebiasi(itime,gs_obs(itime)%iasi%t%val,gs_obs(itime)%iasi%q%val,tl_obs(itime)%iasi%t%val,tl_obs(itime)%iasi%q%val)
     endif
! cris observation
     if (da_cris) then
       call makebcris(itime,gs_obs(itime)%cris%t%val,gs_obs(itime)%cris%q%val,tl_obs(itime)%cris%t%val,tl_obs(itime)%cris%q%val)
     endif
! atms observation
     if (da_atms) then
       call makebatms(itime,gs_obs(itime)%atms%t%val,tl_obs(itime)%atms%t%val)
     endif
! atmswv observation
     if (da_atmswv) then
       call makebatmswv(itime,gs_obs(itime)%atmswv%t%val,gs_obs(itime)%atmswv%q%val,tl_obs(itime)%atmswv%t%val,tl_obs(itime)%atmswv%q%val)
     endif
! mhs observation
     if (da_mhs) then
       call makebmhs(itime,gs_obs(itime)%mhs%t%val,gs_obs(itime)%mhs%q%val,tl_obs(itime)%mhs%t%val,tl_obs(itime)%mhs%q%val)
     endif
! csr observation
     if (da_csr) then
       call makebcsr(itime,gs_obs(itime)%csr%t%val,gs_obs(itime)%csr%q%val,tl_obs(itime)%csr%t%val,tl_obs(itime)%csr%q%val)
     endif
   enddo ! itime
!
   ierr = timingstart('AdjEigToObs_Make_RHS_b') !---profiling
   b_eig = 0.d0
   b_eig_chi = 0.d0
   b_eig_tps = 0.d0
   b_eig_q = 0.d0
   b_eig_a = 0.d0
   call adjeigtoobs(gs_md(:,:,:,:,:),tl_obs(:),b_eig,b_eig_chi,b_eig_tps,b_eig_q,b_eig_a)
   ierr = timingstop('AdjEigToObs_Make_RHS_b') !---profiling
   b_eig = b_eig-x_eig_acc
   b_eig_chi = b_eig_chi-x_eig_acc_chi
   b_eig_tps = b_eig_tps-x_eig_acc_tps
   b_eig_q = b_eig_q-x_eig_acc_q
   b_eig_a = b_eig_a-x_eig_acc_a
   return
!
   end subroutine make_rhs_b
!===============================================================================
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine cal_cost_func(gs_md, x_eig, x_eig_chi, x_eig_tps, x_eig_q, x_eig_a, j_o, j_b)
! 2. Variables
! 2-1. In/Output
!
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(in   ) :: gs_md
   real(real_kind), dimension(neig, zwn_nproc), intent(in   ) :: x_eig
   real(real_kind), dimension(neig_chi, zwn_nproc), intent(in   ) :: x_eig_chi
   real(real_kind), dimension(neig_tps, zwn_nproc), intent(in   ) :: x_eig_tps
   real(real_kind), dimension(neig_q, zwn_nproc), intent(in   ) :: x_eig_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl), intent(in   ) :: x_eig_a
   real(real_kind), intent(  out) :: j_b, j_o
! 2-2. Local
   real(real_kind), dimension(np, np, nlevar, nelemd) :: tl_md
   real(real_kind), dimension(neig, zwn_nproc) :: tl_eig, tlx_eig
   real(real_kind), dimension(neig_chi, zwn_nproc) :: tl_eig_chi, tlx_eig_chi
   real(real_kind), dimension(neig_tps, zwn_nproc) :: tl_eig_tps, tlx_eig_tps
   real(real_kind), dimension(neig_q, zwn_nproc) :: tl_eig_q, tlx_eig_q
   real(real_kind) :: tmpmax, tmpmin, maxv, minv
   real(real_kind) :: rmse_an_u, rmse_an_v, rmse_an_t, rmse_an_q, rmse_an_ps, tmp, jo_single, jo_sonde, jo_surface, jo_bogus, jo_aircraft, jo_amv, jo_scatwind, jo_amsua, jo_gpsro, jo_iasi, jo_cris, jo_atms, jo_atmswv, jo_mhs, jo_csr
   integer(int_kind) :: i, j, k, ie, wn, ierr, itime, ismpl
! 1-1. Jo
!
   do itime = 1, ntime
     ierr = timingstart('ModelToObs_Cal_CostFun') !---profiling
     call modeltoobs(itime,gs_md(:,:,:,:,itime),gs_obs(itime))
     ierr = timingstop('ModelToObs_Cal_CostFun') !---profiling
   enddo ! itime
!
   ierr = timingstart('TlmEigToObs_Cal_CostFun') !---profiling
   call tlmeigtoobs(gs_md(:,:,:,:,:),x_eig,x_eig_chi,x_eig_tps,x_eig_q,x_eig_a,tl_obs(:))
   ierr = timingstop('TlmEigToObs_Cal_CostFun') !---profiling
   global_shared_buf(:,2:46) = 0.d0
   do itime = 1,ntime
! Single observation test
     if (singleobs) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nsg_e(ie)
           if ((obs(itime)%single(i, ie)%u%val.ne.sdmiss).and.(gs_obs(itime)%single(i, ie)%u%val.ne.sdmiss).and.(tl_obs(itime)%single(i, ie)%u%val.ne.sdmiss)) then
             global_shared_buf(ie,2) = global_shared_buf(ie,2)+(gs_obs(itime)%single(i,ie)%u%val+tl_obs(itime)%single(i,ie)%u%val-obs(itime)%single(i,ie)%u%val)**2.d0/(obs(itime)%single(i,ie)%u%error)
           endif
           if ((obs(itime)%single(i, ie)%v%val.ne.sdmiss).and.(gs_obs(itime)%single(i, ie)%v%val.ne.sdmiss).and.(tl_obs(itime)%single(i, ie)%v%val.ne.sdmiss)) then
             global_shared_buf(ie,3) = global_shared_buf(ie,3)+(gs_obs(itime)%single(i,ie)%v%val+tl_obs(itime)%single(i,ie)%v%val-obs(itime)%single(i,ie)%v%val)**2.d0/(obs(itime)%single(i,ie)%v%error)
           endif
           if ((obs(itime)%single(i, ie)%t%val.ne.sdmiss).and.(gs_obs(itime)%single(i, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%single(i, ie)%t%val.ne.sdmiss)) then
             global_shared_buf(ie,4) = global_shared_buf(ie,4)+(gs_obs(itime)%single(i,ie)%t%val+tl_obs(itime)%single(i,ie)%t%val-obs(itime)%single(i,ie)%t%val)**2.d0/(obs(itime)%single(i,ie)%t%error)
           endif
           if ((obs(itime)%single(i, ie)%q%val.ne.sdmiss).and.(gs_obs(itime)%single(i, ie)%q%val.ne.sdmiss).and.(tl_obs(itime)%single(i, ie)%q%val.ne.sdmiss)) then
             global_shared_buf(ie,5) = global_shared_buf(ie,5)+(gs_obs(itime)%single(i,ie)%q%val+tl_obs(itime)%single(i,ie)%q%val-obs(itime)%single(i,ie)%q%val)**2.d0/(obs(itime)%single(i,ie)%q%error)
           endif
           if ((obs(itime)%single(i, ie)%ps%val.ne.sdmiss).and.(gs_obs(itime)%single(i, ie)%ps%val.ne.sdmiss).and.(tl_obs(itime)%single(i, ie)%ps%val.ne.sdmiss)) then
             global_shared_buf(ie,6) = global_shared_buf(ie,6)+(gs_obs(itime)%single(i,ie)%ps%val+tl_obs(itime)%single(i,ie)%ps%val-obs(itime)%single(i,ie)%ps%val)**2.d0/(obs(itime)%single(i,ie)%ps%error)
           endif
         enddo
       enddo
     endif
! sonde observation
     if (da_sonde) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nsd_e(1, ie)
           if ((obs(itime)%sonde(i, ie)%u%val.ne.sdmiss).and.(gs_obs(itime)%sonde(i, ie)%u%val.ne.sdmiss).and.(tl_obs(itime)%sonde(i, ie)%u%val.ne.sdmiss)) then
             global_shared_buf(ie,2) = global_shared_buf(ie,2)+(gs_obs(itime)%sonde(i,ie)%u%val+tl_obs(itime)%sonde(i,ie)%u%val-obs(itime)%sonde(i,ie)%u%val)**2.d0/(obs(itime)%sonde(i,ie)%u%error)
             global_shared_buf(ie,6) = global_shared_buf(ie,6)+1.d0
           endif
         enddo
         do i = 1, obsmet(itime)%nsd_e(2, ie)
           if ((obs(itime)%sonde(i, ie)%v%val.ne.sdmiss).and.(gs_obs(itime)%sonde(i, ie)%v%val.ne.sdmiss).and.(tl_obs(itime)%sonde(i, ie)%v%val.ne.sdmiss)) then
             global_shared_buf(ie,3) = global_shared_buf(ie,3)+(gs_obs(itime)%sonde(i,ie)%v%val+tl_obs(itime)%sonde(i,ie)%v%val-obs(itime)%sonde(i,ie)%v%val)**2.d0/(obs(itime)%sonde(i,ie)%v%error)
             global_shared_buf(ie,6) = global_shared_buf(ie,6)+1.d0
           endif
         enddo
         do i = 1, obsmet(itime)%nsd_e(3, ie)
           if ((obs(itime)%sonde(i, ie)%t%val.ne.sdmiss).and.(gs_obs(itime)%sonde(i, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%sonde(i, ie)%t%val.ne.sdmiss)) then
             global_shared_buf(ie,4) = global_shared_buf(ie,4)+(gs_obs(itime)%sonde(i,ie)%t%val+tl_obs(itime)%sonde(i,ie)%t%val-obs(itime)%sonde(i,ie)%t%val)**2.d0/(obs(itime)%sonde(i,ie)%t%error)
             global_shared_buf(ie,6) = global_shared_buf(ie,6)+1.d0
           endif
         enddo
         do i = 1, obsmet(itime)%nsd_e(4, ie)
           if ((obs(itime)%sonde(i, ie)%q%val.ne.sdmiss).and.(gs_obs(itime)%sonde(i, ie)%q%val.ne.sdmiss).and.(tl_obs(itime)%sonde(i, ie)%q%val.ne.sdmiss)) then
             global_shared_buf(ie,5) = global_shared_buf(ie,5)+(gs_obs(itime)%sonde(i,ie)%q%val+tl_obs(itime)%sonde(i,ie)%q%val-obs(itime)%sonde(i,ie)%q%val)**2.d0/(obs(itime)%sonde(i,ie)%q%error)
             global_shared_buf(ie,6) = global_shared_buf(ie,6)+1.d0
           endif
         enddo
       enddo
     endif
! surface observation
     if (da_surface) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nsf_e(ie)
           if ((obs(itime)%surface(i, ie)%u%val.ne.sdmiss).and.(gs_obs(itime)%surface(i, ie)%u%val.ne.sdmiss).and.(tl_obs(itime)%surface(i, ie)%u%val.ne.sdmiss)) then
             global_shared_buf(ie,7) = global_shared_buf(ie,7)+(gs_obs(itime)%surface(i,ie)%u%val+tl_obs(itime)%surface(i,ie)%u%val-obs(itime)%surface(i,ie)%u%val)**2.d0/(obs(itime)%surface(i,ie)%u%error)
             global_shared_buf(ie,12) = global_shared_buf(ie,12)+1.d0
           endif
           if ((obs(itime)%surface(i, ie)%v%val.ne.sdmiss).and.(gs_obs(itime)%surface(i, ie)%v%val.ne.sdmiss).and.(tl_obs(itime)%surface(i, ie)%v%val.ne.sdmiss)) then
             global_shared_buf(ie,8) = global_shared_buf(ie,8)+(gs_obs(itime)%surface(i,ie)%v%val+tl_obs(itime)%surface(i,ie)%v%val-obs(itime)%surface(i,ie)%v%val)**2.d0/(obs(itime)%surface(i,ie)%v%error)
             global_shared_buf(ie,12) = global_shared_buf(ie,12)+1.d0
           endif
           if ((obs(itime)%surface(i, ie)%t%val.ne.sdmiss).and.(gs_obs(itime)%surface(i, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%surface(i, ie)%t%val.ne.sdmiss)) then
             global_shared_buf(ie,9) = global_shared_buf(ie,9)+(gs_obs(itime)%surface(i,ie)%t%val+tl_obs(itime)%surface(i,ie)%t%val-obs(itime)%surface(i,ie)%t%val)**2.d0/(obs(itime)%surface(i,ie)%t%error)
             global_shared_buf(ie,12) = global_shared_buf(ie,12)+1.d0
           endif
           if ((obs(itime)%surface(i, ie)%q%val.ne.sdmiss).and.(gs_obs(itime)%surface(i, ie)%q%val.ne.sdmiss).and.(tl_obs(itime)%surface(i, ie)%q%val.ne.sdmiss)) then
             global_shared_buf(ie,10) = global_shared_buf(ie,10)+(gs_obs(itime)%surface(i,ie)%q%val+tl_obs(itime)%surface(i,ie)%q%val-obs(itime)%surface(i,ie)%q%val)**2.d0/(obs(itime)%surface(i,ie)%q%error)
             global_shared_buf(ie,12) = global_shared_buf(ie,12)+1.d0
           endif
           if ((obs(itime)%surface(i, ie)%p%val.ne.sdmiss).and.(gs_obs(itime)%surface(i, ie)%p%val.ne.sdmiss).and.(tl_obs(itime)%surface(i, ie)%p%val.ne.sdmiss)) then
             global_shared_buf(ie,11) = global_shared_buf(ie,11)+(gs_obs(itime)%surface(i,ie)%p%val+tl_obs(itime)%surface(i,ie)%p%val-obs(itime)%surface(i,ie)%p%val)**2.d0/(obs(itime)%surface(i,ie)%p%error)
             global_shared_buf(ie,12) = global_shared_buf(ie,12)+1.d0
           endif
         enddo
       enddo
     endif
! bogus observation
     if (da_bogus) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nbg_e(ie)
           if ((obs(itime)%bogus(i, ie)%p%val.ne.sdmiss).and.(gs_obs(itime)%bogus(i, ie)%p%val.ne.sdmiss).and.(tl_obs(itime)%bogus(i, ie)%p%val.ne.sdmiss)) then
             global_shared_buf(ie,13) = global_shared_buf(ie,13)+(gs_obs(itime)%bogus(i,ie)%p%val+tl_obs(itime)%bogus(i,ie)%p%val-obs(itime)%bogus(i,ie)%p%val)**2.d0/(obs(itime)%bogus(i,ie)%p%error)
             global_shared_buf(ie,14) = global_shared_buf(ie,14)+1.d0
           endif
         enddo
       enddo
     endif
! aircraft observation
     if (da_aircraft) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nar_e(1, ie)
           if ((obs(itime)%aircraft(i, ie)%u%val.ne.sdmiss).and.(gs_obs(itime)%aircraft(i, ie)%u%val.ne.sdmiss).and.(tl_obs(itime)%aircraft(i, ie)%u%val.ne.sdmiss)) then
             global_shared_buf(ie,15) = global_shared_buf(ie,15)+(gs_obs(itime)%aircraft(i,ie)%u%val+tl_obs(itime)%aircraft(i,ie)%u%val-obs(itime)%aircraft(i,ie)%u%val)**2.d0/(obs(itime)%aircraft(i,ie)%u%error)
             global_shared_buf(ie,18) = global_shared_buf(ie,18)+1.d0
           endif
         enddo
         do i = 1, obsmet(itime)%nar_e(2, ie)
           if ((obs(itime)%aircraft(i, ie)%v%val.ne.sdmiss).and.(gs_obs(itime)%aircraft(i, ie)%v%val.ne.sdmiss).and.(tl_obs(itime)%aircraft(i, ie)%v%val.ne.sdmiss)) then
             global_shared_buf(ie,16) = global_shared_buf(ie,16)+(gs_obs(itime)%aircraft(i,ie)%v%val+tl_obs(itime)%aircraft(i,ie)%v%val-obs(itime)%aircraft(i,ie)%v%val)**2.d0/(obs(itime)%aircraft(i,ie)%v%error)
             global_shared_buf(ie,18) = global_shared_buf(ie,18)+1.d0
           endif
         enddo
         do i = 1, obsmet(itime)%nar_e(3, ie)
           if ((obs(itime)%aircraft(i, ie)%t%val.ne.sdmiss).and.(gs_obs(itime)%aircraft(i, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%aircraft(i, ie)%t%val.ne.sdmiss)) then
             global_shared_buf(ie,17) = global_shared_buf(ie,17)+(gs_obs(itime)%aircraft(i,ie)%t%val+tl_obs(itime)%aircraft(i,ie)%t%val-obs(itime)%aircraft(i,ie)%t%val)**2.d0/(obs(itime)%aircraft(i,ie)%t%error)
             global_shared_buf(ie,18) = global_shared_buf(ie,18)+1.d0
           endif
         enddo
       enddo
     endif
! amv observation
     if (da_amv) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%namv_e(ie)
           if ((obs(itime)%amv(i, ie)%u%val.ne.sdmiss).and.(gs_obs(itime)%amv(i, ie)%u%val.ne.sdmiss).and.(tl_obs(itime)%amv(i, ie)%u%val.ne.sdmiss)) then
             global_shared_buf(ie,19) = global_shared_buf(ie,19)+(gs_obs(itime)%amv(i,ie)%u%val+tl_obs(itime)%amv(i,ie)%u%val-obs(itime)%amv(i,ie)%u%val)**2.d0/(obs(itime)%amv(i,ie)%u%error)
             global_shared_buf(ie,21) = global_shared_buf(ie,21)+1.d0
           endif
           if ((obs(itime)%amv(i, ie)%v%val.ne.sdmiss).and.(gs_obs(itime)%amv(i, ie)%v%val.ne.sdmiss).and.(tl_obs(itime)%amv(i, ie)%v%val.ne.sdmiss)) then
             global_shared_buf(ie,20) = global_shared_buf(ie,20)+(gs_obs(itime)%amv(i,ie)%v%val+tl_obs(itime)%amv(i,ie)%v%val-obs(itime)%amv(i,ie)%v%val)**2.d0/(obs(itime)%amv(i,ie)%v%error)
             global_shared_buf(ie,21) = global_shared_buf(ie,21)+1.d0
           endif
         enddo
       enddo
     endif
! scatwind observation
     if (da_scatwind) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%nscatwind_e(ie)
           if ((obs(itime)%scatwind(i, ie)%u%val.ne.sdmiss).and.(gs_obs(itime)%scatwind(i, ie)%u%val.ne.sdmiss).and.(tl_obs(itime)%scatwind(i, ie)%u%val.ne.sdmiss)) then
             global_shared_buf(ie,22) = global_shared_buf(ie,22)+(gs_obs(itime)%scatwind(i,ie)%u%val+tl_obs(itime)%scatwind(i,ie)%u%val-obs(itime)%scatwind(i,ie)%u%val)**2.d0/(obs(itime)%scatwind(i,ie)%u%error)
             global_shared_buf(ie,24) = global_shared_buf(ie,24)+1.d0
           endif
           if ((obs(itime)%scatwind(i, ie)%v%val.ne.sdmiss).and.(gs_obs(itime)%scatwind(i, ie)%v%val.ne.sdmiss).and.(tl_obs(itime)%scatwind(i, ie)%v%val.ne.sdmiss)) then
             global_shared_buf(ie,23) = global_shared_buf(ie,23)+(gs_obs(itime)%scatwind(i,ie)%v%val+tl_obs(itime)%scatwind(i,ie)%v%val-obs(itime)%scatwind(i,ie)%v%val)**2.d0/(obs(itime)%scatwind(i,ie)%v%error)
             global_shared_buf(ie,24) = global_shared_buf(ie,24)+1.d0
           endif
         enddo
       enddo
     endif
! amsua observation
     if (da_amsua) then
       do j = 1, obsmet(itime)%namsuach
         do ie = nets, nete
           do i = 1, obsmet(itime)%namsua_e(ie)
             if ((obs(itime)%amsua(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%amsua(i, j, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%amsua(i, j, ie)%t%val.ne.sdmiss)) then
               global_shared_buf(ie,25) = global_shared_buf(ie,25)+(gs_obs(itime)%amsua(i,j,ie)%t%val+tl_obs(itime)%amsua(i,j,ie)%t%val-obs(itime)%amsua(i,j,ie)%tb%val)**2.d0/(obs(itime)%amsua(i,j,ie)%tb%error)
               global_shared_buf(ie,26) = global_shared_buf(ie,26)+1.d0
             endif
           enddo
         enddo
       enddo
     endif
! gpsro observation
     if (da_gpsro) then
       do ie = nets, nete
         do i = 1, obsmet(itime)%ngpsro_e(ie)
           if ((obs(itime)%gpsro(i, ie)%ba%val.ne.sdmiss).and.(gs_obs(itime)%gpsro(i, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%gpsro(i, ie)%t%val.ne.sdmiss)) then
             global_shared_buf(ie,27) = global_shared_buf(ie,27)+(gs_obs(itime)%gpsro(i,ie)%t%val+tl_obs(itime)%gpsro(i,ie)%t%val-obs(itime)%gpsro(i,ie)%ba%val)**2.d0/(obs(itime)%gpsro(i,ie)%ba%error)
             global_shared_buf(ie,29) = global_shared_buf(ie,29)+1.d0
           endif
           if ((obs(itime)%gpsro(i, ie)%ba%val.ne.sdmiss).and.(gs_obs(itime)%gpsro(i, ie)%q%val.ne.sdmiss).and.(tl_obs(itime)%gpsro(i, ie)%q%val.ne.sdmiss)) then
             global_shared_buf(ie,28) = global_shared_buf(ie,28)+(gs_obs(itime)%gpsro(i,ie)%q%val+tl_obs(itime)%gpsro(i,ie)%q%val-obs(itime)%gpsro(i,ie)%ba%val)**2.d0/(obs(itime)%gpsro(i,ie)%ba%error)
             global_shared_buf(ie,29) = global_shared_buf(ie,29)+1.d0
           endif
         enddo
       enddo
     endif
! iasi observation
     if (da_iasi) then
       do j = 1, obsmet(itime)%niasich
         do ie = nets, nete
           do i = 1, obsmet(itime)%niasi_e(ie)
             if ((obs(itime)%iasi(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%iasi(i, j, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%iasi(i, j, ie)%t%val.ne.sdmiss)) then
               global_shared_buf(ie,30) = global_shared_buf(ie,30)+(gs_obs(itime)%iasi(i,j,ie)%t%val+tl_obs(itime)%iasi(i,j,ie)%t%val-obs(itime)%iasi(i,j,ie)%tb%val)**2.d0/(obs(itime)%iasi(i,j,ie)%tb%error)
               global_shared_buf(ie,32) = global_shared_buf(ie,32)+1.d0
             endif
             if ((obs(itime)%iasi(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%iasi(i, j, ie)%q%val.ne.sdmiss).and.(tl_obs(itime)%iasi(i, j, ie)%q%val.ne.sdmiss)) then
               global_shared_buf(ie,31) = global_shared_buf(ie,31)+(gs_obs(itime)%iasi(i,j,ie)%q%val+tl_obs(itime)%iasi(i,j,ie)%q%val-obs(itime)%iasi(i,j,ie)%tb%val)**2.d0/(obs(itime)%iasi(i,j,ie)%tb%error)
               global_shared_buf(ie,32) = global_shared_buf(ie,32)+1.d0
             endif
           enddo
         enddo
       enddo
     endif
! cris observation
     if (da_cris) then
       do j = 1, obsmet(itime)%ncrisch
         do ie = nets, nete
           do i = 1, obsmet(itime)%ncris_e(ie)
             if ((obs(itime)%cris(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%cris(i, j, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%cris(i, j, ie)%t%val.ne.sdmiss)) then
               global_shared_buf(ie,33) = global_shared_buf(ie,33)+(gs_obs(itime)%cris(i,j,ie)%t%val+tl_obs(itime)%cris(i,j,ie)%t%val-obs(itime)%cris(i,j,ie)%tb%val)**2.d0/(obs(itime)%cris(i,j,ie)%tb%error)
               global_shared_buf(ie,35) = global_shared_buf(ie,35)+1.d0
             endif
             if ((obs(itime)%cris(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%cris(i, j, ie)%q%val.ne.sdmiss).and.(tl_obs(itime)%cris(i, j, ie)%q%val.ne.sdmiss)) then
               global_shared_buf(ie,34) = global_shared_buf(ie,34)+(gs_obs(itime)%cris(i,j,ie)%q%val+tl_obs(itime)%cris(i,j,ie)%q%val-obs(itime)%cris(i,j,ie)%tb%val)**2.d0/(obs(itime)%cris(i,j,ie)%tb%error)
               global_shared_buf(ie,35) = global_shared_buf(ie,35)+1.d0
             endif
           enddo
         enddo
       enddo
     endif
! atms observation
     if (da_atms) then
       do j = 1, obsmet(itime)%natmsch
         do ie = nets, nete
           do i = 1, obsmet(itime)%natms_e(ie)
             if ((obs(itime)%atms(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%atms(i, j, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%atms(i, j, ie)%t%val.ne.sdmiss)) then
               global_shared_buf(ie,36) = global_shared_buf(ie,36)+(gs_obs(itime)%atms(i,j,ie)%t%val+tl_obs(itime)%atms(i,j,ie)%t%val-obs(itime)%atms(i,j,ie)%tb%val)**2.d0/(obs(itime)%atms(i,j,ie)%tb%error)
               global_shared_buf(ie,37) = global_shared_buf(ie,37)+1.d0
             endif
           enddo
         enddo
       enddo
     endif
! atmswv observation
     if (da_atmswv) then
       do j = 1, obsmet(itime)%natmswvch
         do ie = nets, nete
           do i = 1, obsmet(itime)%natmswv_e(ie)
             if ((obs(itime)%atmswv(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%atmswv(i, j, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%atmswv(i, j, ie)%t%val.ne.sdmiss)) then
               global_shared_buf(ie,38) = global_shared_buf(ie,38)+(gs_obs(itime)%atmswv(i,j,ie)%t%val+tl_obs(itime)%atmswv(i,j,ie)%t%val-obs(itime)%atmswv(i,j,ie)%tb%val)**2.d0/(obs(itime)%atmswv(i,j,ie)%tb%error)
               global_shared_buf(ie,40) = global_shared_buf(ie,40)+1.d0
             endif
             if ((obs(itime)%atmswv(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%atmswv(i, j, ie)%q%val.ne.sdmiss).and.(tl_obs(itime)%atmswv(i, j, ie)%q%val.ne.sdmiss)) then
               global_shared_buf(ie,39) = global_shared_buf(ie,39)+(gs_obs(itime)%atmswv(i,j,ie)%q%val+tl_obs(itime)%atmswv(i,j,ie)%q%val-obs(itime)%atmswv(i,j,ie)%tb%val)**2.d0/(obs(itime)%atmswv(i,j,ie)%tb%error)
               global_shared_buf(ie,40) = global_shared_buf(ie,40)+1.d0
             endif
           enddo
         enddo
       enddo
     endif
! mhs observation
     if (da_mhs) then
       do j = 1, obsmet(itime)%nmhsch
         do ie = nets, nete
           do i = 1, obsmet(itime)%nmhs_e(ie)
             if ((obs(itime)%mhs(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%mhs(i, j, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%mhs(i, j, ie)%t%val.ne.sdmiss)) then
               global_shared_buf(ie,41) = global_shared_buf(ie,41)+(gs_obs(itime)%mhs(i,j,ie)%t%val+tl_obs(itime)%mhs(i,j,ie)%t%val-obs(itime)%mhs(i,j,ie)%tb%val)**2.d0/(obs(itime)%mhs(i,j,ie)%tb%error)
               global_shared_buf(ie,43) = global_shared_buf(ie,43)+1.d0
             endif
             if ((obs(itime)%mhs(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%mhs(i, j, ie)%q%val.ne.sdmiss).and.(tl_obs(itime)%mhs(i, j, ie)%q%val.ne.sdmiss)) then
               global_shared_buf(ie,42) = global_shared_buf(ie,42)+(gs_obs(itime)%mhs(i,j,ie)%q%val+tl_obs(itime)%mhs(i,j,ie)%q%val-obs(itime)%mhs(i,j,ie)%tb%val)**2.d0/(obs(itime)%mhs(i,j,ie)%tb%error)
               global_shared_buf(ie,43) = global_shared_buf(ie,43)+1.d0
             endif
           enddo
         enddo
       enddo
     endif
! csr observation
     if (da_csr) then
       do j = 1, obsmet(itime)%ncsrch
         do ie = nets, nete
           do i = 1, obsmet(itime)%ncsr_e(ie)
             if ((obs(itime)%csr(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%csr(i, j, ie)%t%val.ne.sdmiss).and.(tl_obs(itime)%csr(i, j, ie)%t%val.ne.sdmiss)) then
               global_shared_buf(ie,44) = global_shared_buf(ie,44)+(gs_obs(itime)%csr(i,j,ie)%t%val+tl_obs(itime)%csr(i,j,ie)%t%val-obs(itime)%csr(i,j,ie)%tb%val)**2.d0/(obs(itime)%csr(i,j,ie)%tb%error)
               global_shared_buf(ie,46) = global_shared_buf(ie,46)+1.d0
             endif
             if ((obs(itime)%csr(i, j, ie)%tb%val.ne.sdmiss).and.(gs_obs(itime)%csr(i, j, ie)%q%val.ne.sdmiss).and.(tl_obs(itime)%csr(i, j, ie)%q%val.ne.sdmiss)) then
               global_shared_buf(ie,45) = global_shared_buf(ie,45)+(gs_obs(itime)%csr(i,j,ie)%q%val+tl_obs(itime)%csr(i,j,ie)%q%val-obs(itime)%csr(i,j,ie)%tb%val)**2.d0/(obs(itime)%csr(i,j,ie)%tb%error)
               global_shared_buf(ie,46) = global_shared_buf(ie,46)+1.d0
             endif
           enddo
         enddo
       enddo
     endif
   enddo ! itime
! 1-2. Jb
   tmp = 0.d0
   do wn = 1,zwn_nproc
     tmp = tmp+sum((x_eig(:,wn)-x_eig_acc(:,wn))**2.d0)+sum((x_eig_chi(:,wn)-x_eig_acc_chi(:,wn))**2.d0)+sum((x_eig_tps(:,wn)-x_eig_acc_tps(:,wn))**2.d0)+sum((x_eig_q(:,wn)-x_eig_acc_q(:,wn))**2.d0)
   enddo
   do ismpl = 1,nsmpl
     do wn = 1, zwn_nproc_a
       tmp = tmp+sum((x_eig_a(:,wn,ismpl)-x_eig_a(:,wn,ismpl))**2.d0)
     enddo
   enddo ! ismpl
   global_shared_buf(:,1) = 0.d0
   global_shared_buf(nets:nete,1) = tmp/dble(nete-nets+1)
! 1-3. Communication and evaluation
   ierr = timingstart('Wrap_Repro_Sum_rmse_Cal_CostFun') !---profiling
   call wrap_repro_sum(nvars = 46,comm = hybrid%par%comm)
   ierr = timingstop('Wrap_Repro_Sum_rmse_Cal_CostFun') !---profiling
   j_b = 0.5d0*global_shared_sum(1)
   j_o = 0.5d0*(global_shared_sum(2)+global_shared_sum(3)+global_shared_sum(4)+global_shared_sum(5)+! sonde
   global_shared_sum(7)+global_shared_sum(8)+global_shared_sum(9)+global_shared_sum(10)+global_shared_sum(11)+! surface
   global_shared_sum(13)+! bogus
   global_shared_sum(15)+global_shared_sum(16)+global_shared_sum(17)+! aircraft
   global_shared_sum(19)+global_shared_sum(20)+! amv
   global_shared_sum(22)+global_shared_sum(23)+! scatwind
   global_shared_sum(25)+! amsua
   global_shared_sum(27)+global_shared_sum(28)+! gpsro
   global_shared_sum(30)+global_shared_sum(31)+! iasi
   global_shared_sum(33)+global_shared_sum(34)+! cris
   global_shared_sum(36)+! atms
   global_shared_sum(38)+global_shared_sum(39)+! atmswv
   global_shared_sum(41)+global_shared_sum(42)+! mhs
   global_shared_sum(44)+global_shared_sum(45)) ! csr
   if (singleobs) then
     jo_single = 0.5d0*(global_shared_sum(2)+global_shared_sum(3)+global_shared_sum(4)+global_shared_sum(5)+global_shared_sum(6))
   endif
   if (da_sonde) then
     jo_sonde = 0.5d0*(global_shared_sum(2)+global_shared_sum(3)+global_shared_sum(4)+global_shared_sum(5))
   endif
   if (da_surface) then
     jo_surface = 0.5d0*(global_shared_sum(7)+global_shared_sum(8)+global_shared_sum(9)+global_shared_sum(10)+global_shared_sum(11))
   endif
   if (da_bogus) then
     jo_bogus = 0.5d0*(global_shared_sum(13))
   endif
   if (da_aircraft) then
     jo_aircraft = 0.5d0*(global_shared_sum(15)+global_shared_sum(16)+global_shared_sum(17))
   endif
   if (da_amv) then
     jo_amv = 0.5d0*(global_shared_sum(19)+global_shared_sum(20))
   endif
   if (da_scatwind) then
     jo_scatwind = 0.5d0*(global_shared_sum(22)+global_shared_sum(23))
   endif
   if (da_amsua) then
     jo_amsua = 0.5d0*(global_shared_sum(25))
   endif
   if (da_gpsro) then
     jo_gpsro = 0.5d0*(global_shared_sum(27)+global_shared_sum(28))
   endif
   if (da_iasi) then
     jo_iasi = 0.5d0*(global_shared_sum(30)+global_shared_sum(31))
   endif
   if (da_cris) then
     jo_cris = 0.5d0*(global_shared_sum(33)+global_shared_sum(34))
   endif
   if (da_atms) then
     jo_atms = 0.5d0*(global_shared_sum(36))
   endif
   if (da_atmswv) then
     jo_atmswv = 0.5d0*(global_shared_sum(38)+global_shared_sum(39))
   endif
   if (da_mhs) then
     jo_mhs = 0.5d0*(global_shared_sum(41)+global_shared_sum(42))
   endif
   if (da_csr) then
     jo_csr = 0.5d0*(global_shared_sum(44)+global_shared_sum(45))
   endif
   if (par%ismasterproc) write(iulog,*) ' '
   if (par%ismasterproc) write(iulog,*) 'Check j'
   if (par%ismasterproc) write(iulog,*) 'Jb = ',j_b
   if (par%ismasterproc) write(iulog,*) 'Jo = ',j_o
   if (par%ismasterproc) write(iulog,*) 'J = ',j_b+j_o
   if (par%ismasterproc) write(iulog,*) ' '
   if (singleobs) then
     if (par%ismasterproc) write(iulog,*) 'jo_single = ',jo_single
   endif
   if (da_sonde) then
     if (par%ismasterproc) write(iulog,*) 'jo_sonde = ',jo_sonde
   endif
   if (da_surface) then
     if (par%ismasterproc) write(iulog,*) 'jo_surface = ',jo_surface
   endif
   if (da_bogus) then
     if (par%ismasterproc) write(iulog,*) 'jo_bogus = ',jo_bogus
   endif
   if (da_aircraft) then
     if (par%ismasterproc) write(iulog,*) 'jo_aircraft = ',jo_aircraft
   endif
   if (da_amv) then
     if (par%ismasterproc) write(iulog,*) 'jo_amv = ',jo_amv
   endif
   if (da_scatwind) then
     if (par%ismasterproc) write(iulog,*) 'jo_scatwind = ',jo_scatwind
   endif
   if (da_amsua) then
     if (par%ismasterproc) write(iulog,*) 'jo_amsua = ',jo_amsua
   endif
   if (da_gpsro) then
     if (par%ismasterproc) write(iulog,*) 'jo_gpsro = ',jo_gpsro
   endif
   if (da_iasi) then
     if (par%ismasterproc) write(iulog,*) 'jo_iasi = ',jo_iasi
   endif
   if (da_cris) then
     if (par%ismasterproc) write(iulog,*) 'jo_cris = ',jo_cris
   endif
   if (da_atms) then
     if (par%ismasterproc) write(iulog,*) 'jo_atms = ',jo_atms
   endif
   if (da_atmswv) then
     if (par%ismasterproc) write(iulog,*) 'jo_atmswv = ',jo_atmswv
   endif
   if (da_mhs) then
     if (par%ismasterproc) write(iulog,*) 'jo_mhs = ',jo_mhs
   endif
   if (da_csr) then
     if (par%ismasterproc) write(iulog,*) 'jo_csr = ',jo_csr
   endif
   if (da_sonde) then
     if (par%ismasterproc) write(iulog,*) '# of sonde observations = ',dint(global_shared_sum(6))
   endif
   if (da_surface) then
     if (par%ismasterproc) write(iulog,*) '# of surface observations = ',dint(global_shared_sum(12))
   endif
   if (da_bogus) then
     if (par%ismasterproc) write(iulog,*) '# of bogus observations = ',dint(global_shared_sum(14))
   endif
   if (da_aircraft) then
     if (par%ismasterproc) write(iulog,*) '# of aircraft observations = ',dint(global_shared_sum(18))
   endif
   if (da_amv) then
     if (par%ismasterproc) write(iulog,*) '# of amv observations = ',dint(global_shared_sum(21))
   endif
   if (da_scatwind) then
     if (par%ismasterproc) write(iulog,*) '# of scatwind observations = ',dint(global_shared_sum(24))
   endif
   if (da_amsua) then
     if (par%ismasterproc) write(iulog,*) '# of amsua observations = ',dint(global_shared_sum(26))
   endif
   if (da_gpsro) then
     if (par%ismasterproc) write(iulog,*) '# of gpsro observations = ',dint(global_shared_sum(29))
   endif
   if (da_iasi) then
     if (par%ismasterproc) write(iulog,*) '# of iasi observations = ',dint(global_shared_sum(32))
   endif
   if (da_cris) then
     if (par%ismasterproc) write(iulog,*) '# of cris observations = ',dint(global_shared_sum(35))
   endif
   if (da_atms) then
     if (par%ismasterproc) write(iulog,*) '# of atms observations = ',dint(global_shared_sum(37))
   endif
   if (da_atmswv) then
     if (par%ismasterproc) write(iulog,*) '# of atmswv observations = ',dint(global_shared_sum(40))
   endif
   if (da_mhs) then
     if (par%ismasterproc) write(iulog,*) '# of mhs observations = ',dint(global_shared_sum(43))
   endif
   if (da_csr) then
     if (par%ismasterproc) write(iulog,*) '# of csr observations = ',dint(global_shared_sum(46))
   endif
   return
!
   end subroutine cal_cost_func
!===============================================================================
!--------- modified for alpha (HJS)
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine eigtospec(l_zwn_nproc, l_zwn2, l_zwn_blength, l_zwn_bstart, eigen, eigen_chi, eigen_tps, eigen_q, sp)
!=============================================================================
! A. Declaration
!=============================================================================
! 2. Variables
! 2-1. Input
!
   integer(int_kind), intent(in   ) :: l_zwn_nproc, l_zwn2
   integer(int_kind), dimension(:), allocatable, intent(in   ) :: l_zwn_blength, l_zwn_bstart
   real(real_kind), dimension(neig, l_zwn_nproc), intent(in   ) :: eigen
   real(real_kind), dimension(neig_chi, l_zwn_nproc), intent(in   ) :: eigen_chi
   real(real_kind), dimension(neig_tps, l_zwn_nproc), intent(in   ) :: eigen_tps
   real(real_kind), dimension(neig_q, l_zwn_nproc), intent(in   ) :: eigen_q
! 2-2. Output
   real(real_kind), dimension(nlevar, l_zwn2), intent(  out) :: sp
! 2-3. Local
   real(real_kind), dimension(nlevar, l_zwn_nproc) :: sp_tmp
   integer(int_kind) :: k, wn, iproc, ierr
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!
   ierr = timingstart('Run beforebcast_fromeigtospec')
   sp_tmp = 0.d0
   do wn = 1,l_zwn_nproc
     do k = 1, neig
       sp_tmp(1:nlev,wn) = sp_tmp(1:nlev,wn)+bcov_sp(1:nlev,k,wn)*dsqrt(eigv(k,wn))*eigen(k,wn)
     enddo
     do k = 1, neig_chi
       sp_tmp(nlev+1:2*nlev,wn) = sp_tmp(nlev+1:2*nlev,wn)+bcov_sp_chi(1:nlev,k,wn)*dsqrt(eigv_chi(k,wn))*eigen_chi(k,wn)
     enddo
     do k = 1, neig_tps
       sp_tmp(2*nlev+1:3*nlev,wn) = sp_tmp(2*nlev+1:3*nlev,wn)+bcov_sp_tps(1:nlev,k,wn)*dsqrt(eigv_tps(k,wn))*eigen_tps(k,wn)
       sp_tmp(4*nlev+1,wn) = sp_tmp(4*nlev+1,wn)+bcov_sp_tps(nlev+1,k,wn)*dsqrt(eigv_tps(k,wn))*eigen_tps(k,wn)
     enddo
     do k = 1, neig_q
       sp_tmp(3*nlev+1:4*nlev,wn) = sp_tmp(3*nlev+1:4*nlev,wn)+bcov_sp_q(:,k,wn)*dsqrt(eigv_q(k,wn))*eigen_q(k,wn)
     enddo
   enddo
   ierr = timingstop('Run beforebcast_fromeigtospec')
   ierr = timingstart('Run bcast4spectoctrl_fromeigtospec')
   call mpi_allgatherv(sp_tmp,nlevar*l_zwn_nproc,par_double_precision,sp,nlevar*l_zwn_blength,nlevar*l_zwn_bstart,par_double_precision,par%comm,ierr)
   ierr = timingstop('Run bcast4spectoctrl_fromeigtospec')
   return
!
   end subroutine eigtospec
! Desinged for alpha-variable
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine eigtospecens(l_zwn_nproc, l_zwn2, l_zwn_blength, l_zwn_bstart, eigen_a, sp)
!=============================================================================
! A. Declaration
!=============================================================================
! 2. Variables
! 2-1. Input
!
   integer(int_kind), intent(in   ) :: l_zwn_nproc, l_zwn2
   integer(int_kind), dimension(:), allocatable, intent(in   ) :: l_zwn_blength, l_zwn_bstart
   real(real_kind), dimension(nlev+1, l_zwn_nproc), intent(in   ) :: eigen_a
! 2-2. Output
   real(real_kind), dimension(nlevar, l_zwn2), intent(  out) :: sp
! 2-3. Local
   real(real_kind), dimension(nlevar, l_zwn_nproc) :: sp_tmp
   integer(int_kind) :: k, wn, iproc, ierr
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!
   ierr = timingstart('Run beforebcast_fromeigtospec')
   sp_tmp = 0.d0
   do wn = 1,l_zwn_nproc
!       DO k  = 1, nlev
!          sp_tmp(3*nlev+1:4*nlev,wn) =    &
!             sp_tmp(3*nlev+1:4*nlev,wn) + &
!             bcov_sp_a(:,k)*              &
!             DSQRT(eigv_a(k,wn))*         &
!             eigen_a(k,wn)
!       END DO
     do k = 1, nlev+1
       sp_tmp(2*nlev+1:3*nlev,wn) = sp_tmp(2*nlev+1:3*nlev,wn)+bcov_sp_a(1:nlev,k)*dsqrt(eigv_a(k,wn))*eigen_a(k,wn)
       sp_tmp(4*nlev+1,wn) = sp_tmp(4*nlev+1,wn)+bcov_sp_a(nlev+1,k)*dsqrt(eigv_a(k,wn))*eigen_a(k,wn)
     enddo
!       sp_tmp(1:nlev,wn)=sp_tmp(3*nlev+1:4*nlev,wn)
!       sp_tmp(nlev+1:2*nlev,wn)=sp_tmp(3*nlev+1:4*nlev,wn)
!       sp_tmp(2*nlev+1:3*nlev,wn)=sp_tmp(3*nlev+1:4*nlev,wn)
!       sp_tmp(4*nlev+1,wn)=sp_tmp(4*nlev,wn)
     sp_tmp(1:nlev,wn) = sp_tmp(2*nlev+1:3*nlev,wn)
     sp_tmp(nlev+1:2*nlev,wn) = sp_tmp(2*nlev+1:3*nlev,wn)
     sp_tmp(3*nlev+1:4*nlev,wn) = sp_tmp(2*nlev+1:3*nlev,wn)
   enddo
   ierr = timingstop('Run beforebcast_fromeigtospec')
   ierr = timingstart('Run bcast4spectoctrl_fromeigtospec')
   call mpi_allgatherv(sp_tmp,nlevar*l_zwn_nproc,par_double_precision,sp,nlevar*l_zwn_blength,nlevar*l_zwn_bstart,par_double_precision,par%comm,ierr)
   ierr = timingstop('Run bcast4spectoctrl_fromeigtospec')
   return
!
   end subroutine eigtospecens
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
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine tlmeigtoobs(l_gs_md, tl_eigen, tl_eigen_chi, tl_eigen_tps, tl_eigen_q, tl_eigen_a, l_tl_obs)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
! 2. Variables
! 2-1. Input
!
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(in   ) :: l_gs_md
   real(real_kind), dimension(neig, zwn_nproc), intent(in   ) :: tl_eigen
   real(real_kind), dimension(neig_chi, zwn_nproc), intent(in   ) :: tl_eigen_chi
   real(real_kind), dimension(neig_tps, zwn_nproc), intent(in   ) :: tl_eigen_tps
   real(real_kind), dimension(neig_q, zwn_nproc), intent(in   ) :: tl_eigen_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl), intent(in   ) :: tl_eigen_a
! 2-2. Output
   type(obsbg_type), dimension(ntime), intent(  out) :: l_tl_obs
! 2-3. Local
   real(real_kind), dimension(np, np, nlevar, nelemd) :: tl_ct, tl_ctt, tl_ct_a
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime) :: tl_md
   real(real_kind), dimension(nlevar, zwn2) :: tl_sp
   real(real_kind), dimension(nlevar, zwn2_a) :: tl_sp_a
   integer(int_kind) :: i, j, k, l, ie, wn, vr, iproc, ierr, itime
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!
   call tlmeigtomodel(l_gs_md(:,:,:,:,:),tl_eigen,tl_eigen_chi,tl_eigen_tps,tl_eigen_q,tl_eigen_a,tl_md(:,:,:,:,:))
   do itime = 1,ntime
     ierr = timingstart('Run tlmmodeltoobs_fromtlmeigtoobs')
     call tlmmodeltoobs(itime,l_gs_md(:,:,:,:,itime),! in
     tl_md(:,:,:,:,itime),! inout
     l_tl_obs(itime)) ! out
     ierr = timingstop('Run tlmmodeltoobs_fromtlmeigtoobs')
   enddo ! itime
   return
!
   end subroutine tlmeigtoobs
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
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine tlmeigtomodel(l_gs_md, tl_eigen, tl_eigen_chi, tl_eigen_tps, tl_eigen_q, tl_eigen_a, tl_md)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
! 2. Variables
! 2-1. Input
!
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(in   ) :: l_gs_md
   real(real_kind), dimension(neig, zwn_nproc), intent(in   ) :: tl_eigen
   real(real_kind), dimension(neig_chi, zwn_nproc), intent(in   ) :: tl_eigen_chi
   real(real_kind), dimension(neig_tps, zwn_nproc), intent(in   ) :: tl_eigen_tps
   real(real_kind), dimension(neig_q, zwn_nproc), intent(in   ) :: tl_eigen_q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl), intent(in   ) :: tl_eigen_a
! 2-2. Output
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(  out) :: tl_md
! 2-3. Local
! up
!    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
   real(real_kind), dimension(nups_l, nlevar) :: tl_ct, tl_ctt, tl_ct_a
! up
!REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime) &
   real(real_kind), dimension(nups_l, nlevar, ntime) :: tl_ct_hyb
   real(real_kind), dimension(np, np, nlevar, nelemd) :: tl_ct_hyb_ep
   real(real_kind), dimension(nlevar, zwn2) :: tl_sp
   real(real_kind), dimension(nlevar, zwn2_a) :: tl_sp_a
   integer(int_kind) :: i, j, k, l, ie, wn, vr, iproc, ierr, itime, ismpl
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!---------------------------------------------------------
! Climatological Part
!
   tl_sp = 0.d0
   tl_ct = 0.d0
   tl_ctt = 0.d0
   ierr = timingstart('Run eigtospec_fromtlmeigtoobs')
   call eigtospec(zwn_nproc,zwn2,zwn_blength,zwn_bstart,tl_eigen,tl_eigen_chi,tl_eigen_tps,tl_eigen_q,tl_sp)
   ierr = timingstop('Run eigtospec_fromtlmeigtoobs')
   ierr = timingstart('Run spectoctrl_fromtlmeigtomodel')
   call spectoctrl(zwn2,tl_sp,tl_ct)
   ierr = timingstop('Run spectoctrl_fromtlmeigtomodel')
!up
!tl_ctt = tl_ct *DSQRT(errvar)
   tl_ctt = tl_ct*dsqrt(errvar_up)
   do itime = 1,ntime
!tl_ct_hyb(:,:,:,:,itime) = DSQRT(b_clim)*tl_ctt(:,:,:,:)
     do k = 1, nlevar
! up
!tl_ct_hyb(:,:,k,:,itime) = DSQRT(b_clim_v(k))*tl_ctt(:,:,k,:)
       tl_ct_hyb(:,k,itime) = dsqrt(b_clim_v(k))*tl_ctt(:,k)
     enddo
   enddo
!-----------------------------------------------------------
!---------------------------------------------------------------
! Add Ensemble Part
   do ismpl = 1,nsmpl
     ierr = timingstart('Run eigtospecens_fromtlmeigtomodel')
     call eigtospecens(zwn_nproc_a,zwn2_a,zwn_blength_a,zwn_bstart_a,tl_eigen_a(:,:,ismpl),tl_sp_a(:,:))
     ierr = timingstop('Run eigtospecens_fromtlmeigtomodel')
     ierr = timingstart('Run spectoctrl_fromtlmeigtomodel')
     call spectoctrl(zwn2_a,tl_sp_a,tl_ct_a)
     ierr = timingstop('Run spectoctrl_fromtlmeigtomodel')
     do itime = 1,ntime
       do k = 1, nlevar
         tl_ct_hyb(:,k,itime) = tl_ct_hyb(:,k,itime)+dsqrt(b_ens_v(k))*bg_smpl(:,k,itime,ismpl)*tl_ct_a(:,k)
       enddo ! k
     enddo ! itime
   enddo ! ismpl
!------------------------------------------------------------------
   do itime = 1,ntime
!---------------Parameter Transform-------------------
     ierr = timingstart('Run tlmctrltomodel_fromtlmeigtomodel')
     ierr = timingstart('UP2EP_fromTlmEigToModel')
     call convertup2ep_hyb(nlevar,tl_ct_hyb(:,:,itime),tl_ct_hyb_ep)
     ierr = timingstop('UP2EP_fromTlmEigToModel')
     call tlmctrltomodel(l_gs_md(:,:,:,:,itime),! in
     bg_uv_rot(:,:,:,:,:,itime),! in
     tl_ct_hyb_ep(:,:,:,:),! in
     tl_md(:,:,:,:,itime)) ! out
     ierr = timingstop('Run tlmctrltomodel_fromtlmeigtomodel')
!-----------------------------------------------------
   enddo ! itime
   return
!
   end subroutine tlmeigtomodel
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
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine adjeigtoobs(l_gs_md, l_ad_obs, ad_eigen, ! for psi
!
   ad_eigen_chi,! for chi
   ad_eigen_tps,! for(t,ps)
   ad_eigen_q,ad_eigen_a)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
! 2. Variables
! 2-1. Input
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(in   ) :: l_gs_md
! 2-2. Input/Output
   type(obsbg_type), dimension(ntime), intent(inout) :: l_ad_obs
   real(real_kind), dimension(neig, zwn_nproc), intent(inout) :: ad_eigen ! adjoint variable for psi
   real(real_kind), dimension(neig_chi, zwn_nproc), intent(inout) :: ad_eigen_chi ! adjoint variable for chi
   real(real_kind), dimension(neig_tps, zwn_nproc), intent(inout) :: ad_eigen_tps ! adjoint variable for(t, ps)
   real(real_kind), dimension(neig_q, zwn_nproc), intent(inout) :: ad_eigen_q ! adjoint variable for q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl), intent(inout) :: ad_eigen_a ! adjoint variable for alpha
! 2-3. Local
   real(real_kind), dimension(np, np, nlevar, nelemd) :: ad_ct, ad_ctt
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime) :: ad_md
   real(real_kind), dimension(nlevar, zwn_nproc) :: ad_sp
   integer(int_kind) :: i, j, k, ie, wn, vr, iproc, ierr, itime
!=============================================================================
!=============================================================================
! B. Prepare a basic state
!=============================================================================
!   md = gs_md
!    ad_ct = 0.D0
!=============================================================================
!
   ad_md = 0.d0
   do itime = ntime,1,-1
!=============================================================================
! C. Main body of adjoint process
!=============================================================================
     ierr = timingstart('Run adjmodeltoobs_fromadjeigtoobs')
     call adjmodeltoobs(itime,l_gs_md(:,:,:,:,itime),! in
     l_ad_obs(itime),! inout
     ad_md(:,:,:,:,itime)) ! inout
     ierr = timingstop('Run adjmodeltoobs_fromadjeigtoobs')
   enddo ! itime
   call adjeigtomodel(l_gs_md(:,:,:,:,:),ad_md(:,:,:,:,:),ad_eigen,ad_eigen_chi,ad_eigen_tps,ad_eigen_q,ad_eigen_a)
   return
!
   end subroutine adjeigtoobs
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
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine adjeigtomodel(l_gs_md, ad_md, ad_eigen, ! for psi
!
   ad_eigen_chi,! for chi
   ad_eigen_tps,! for(t,ps)
   ad_eigen_q,ad_eigen_a)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
! 2. Variables
! 2-1. Input
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(in   ) :: l_gs_md
! 2-2. Input/Output
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(inout) :: ad_md
   real(real_kind), dimension(neig, zwn_nproc), intent(inout) :: ad_eigen ! adjoint variable for psi
   real(real_kind), dimension(neig_chi, zwn_nproc), intent(inout) :: ad_eigen_chi ! adjoint variable for chi
   real(real_kind), dimension(neig_tps, zwn_nproc), intent(inout) :: ad_eigen_tps ! adjoint variable for(t, ps)
   real(real_kind), dimension(neig_q, zwn_nproc), intent(inout) :: ad_eigen_q ! adjoint variable for q
   real(real_kind), dimension(nlev+1, zwn_nproc_a, nsmpl), intent(inout) :: ad_eigen_a ! adjoint variable for alpha
! 2-3. Local
!up
!REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
   real(real_kind), dimension(nups_l, nlevar) :: ad_ct, ad_ctt, ad_ct_a
!up
!REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd,ntime) &
   real(real_kind), dimension(nups_l, nlevar, ntime) :: ad_ct_hyb
   real(real_kind), dimension(np, np, nlevar, nelemd) :: ad_ct_hyb_ep
   real(real_kind), dimension(nlevar, zwn_nproc) :: ad_sp, ad_sp_a
   integer(int_kind) :: i, j, k, ie, wn, vr, iproc, ierr, itime, ismpl
!=============================================================================
!=============================================================================
! B. Prepare a basic state
!=============================================================================
!
   do itime = ntime, 1,-1
!up
!ad_ct_hyb(:,:,:,:,itime)  = 0.D0
     ad_ct_hyb_ep = 0.d0
     ierr = timingstart('Run adjctrltomodel_fromadjeigtomodel')
     call adjctrltomodel(l_gs_md(:,:,:,:,itime),bg_uv_rot(:,:,:,:,:,itime),! in
     ad_md(:,:,:,:,itime),! inout
     ad_ct_hyb_ep(:,:,:,:)) ! inout
!up
     ierr = timingstart('EP2UP_fromAdjEigToModel')
     call convertep2ep0_hyb(nlevar,ad_ct_hyb_ep)
     call convertep2up_hyb(nlevar,ad_ct_hyb_ep,ad_ct_hyb(:,:,itime))
     ierr = timingstop('EP2UP_fromAdjEigToModel')
     ierr = timingstop('Run adjctrltomodel_fromadjeigtomodel')
   enddo ! itime
!--------------------------------------------------
! Adjoint of adding ensemble part
!
   do ismpl = nsmpl, 1,-1
     ad_ct_a = 0.d0
     do itime = ntime,1,-1
       do k = nlevar, 1,-1
         ad_ct_a(:,k) = ad_ct_a(:,k)+dsqrt(b_ens_v(k))*bg_smpl(:,k,itime,ismpl)*ad_ct_hyb(:,k,itime)
       enddo
     enddo
     ad_sp_a = 0.d0
     ierr = timingstart('Run adjspectoctrl_fromadjeigtomodel')
     call adjspectoctrl_opt(zwn2_a,zwn_nproc_a,zwn_iproc_a,zwn_bstart_a,zwn_blength_a,ad_ct_a,ad_sp_a)
     ierr = timingstop('Run adjspectoctrl_fromadjeigtomodel')
     ierr = timingstart('Run afteradjspectoctrl_fromadjeigtomodel')
     call adjeigtospecens(zwn_nproc_a,ad_sp_a,ad_eigen_a(:,:,ismpl))
     ierr = timingstop('Run afteradjspectoctrl_fromadjeigtomodel')
   enddo ! ismpl
!--------------------------------------------------------------------------
!---------------------------------------------------------------------
! Adjoint of climatological part
!
   ad_ctt = 0.d0
   do itime = ntime,1,-1
     do k = nlevar, 1,-1
       ad_ctt(:,k) = ad_ctt(:,k)+dsqrt(b_clim_v(k))*ad_ct_hyb(:,k,itime)
       ad_ct_hyb(:,k,itime) = 0.d0
     enddo ! k
   enddo ! itime
   ad_ct = ad_ctt*dsqrt(errvar_up)
   ad_sp = 0.d0
   ierr = timingstart('Run adjspectoctrl_fromadjeigtomodel')
   call adjspectoctrl_opt(zwn2,zwn_nproc,zwn_iproc,zwn_bstart,zwn_blength,ad_ct,ad_sp)
   ierr = timingstop('Run adjspectoctrl_fromadjeigtomodel')
   ierr = timingstart('Run afteradjspectoctrl_fromadjeigtomodel')
   call adjeigtospec(zwn_nproc,ad_sp,ad_eigen,ad_eigen_chi,ad_eigen_tps,ad_eigen_q)
   ierr = timingstop('Run afteradjspectoctrl_fromadjeigtomodel')
!------------------------------------------------------------------------
   return
!
   end subroutine adjeigtomodel
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine adjeigtospec(l_zwn_nproc, ad_sp, ad_eigen, ad_eigen_chi, ad_eigen_tps, ad_eigen_q)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
! 2. Variables
! 2-1. Input
!
   integer(int_kind), intent(in   ) :: l_zwn_nproc
   real(real_kind), dimension(nlevar, l_zwn_nproc), intent(in   ) :: ad_sp
! 2-2. Input/Output
   real(real_kind), dimension(neig, l_zwn_nproc), intent(inout) :: ad_eigen
   real(real_kind), dimension(neig_chi, l_zwn_nproc), intent(inout) :: ad_eigen_chi
   real(real_kind), dimension(neig_tps, l_zwn_nproc), intent(inout) :: ad_eigen_tps
   real(real_kind), dimension(neig_q, l_zwn_nproc), intent(inout) :: ad_eigen_q
! 2-3. Local
   integer(int_kind) :: k, wn
!=============================================================================
!=============================================================================
! C. Main body of adjoint process
!=============================================================================
!
   do wn = 1, l_zwn_nproc
     do k = 1, neig
       ad_eigen(k,wn) = dsqrt(eigv(k,wn))*sum(bcov_sp(1:nlev,k,wn)*ad_sp(1:nlev,wn))
! Check
!IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn)
     enddo
     do k = 1, neig_chi
       ad_eigen_chi(k,wn) = dsqrt(eigv_chi(k,wn))*sum(bcov_sp_chi(:,k,wn)*ad_sp(nlev+1:2*nlev,wn))
! Check
!IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn)
     enddo
     do k = 1, neig_tps
       ad_eigen_tps(k,wn) = dsqrt(eigv_tps(k,wn))*(sum(bcov_sp_tps(1:nlev,k,wn)*ad_sp(2*nlev+1:3*nlev,wn))+bcov_sp_tps(nlev+1,k,wn)*ad_sp(4*nlev+1,wn))
! Check
!IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn)
     enddo
     do k = 1, neig_q
       ad_eigen_q(k,wn) = dsqrt(eigv_q(k,wn))*sum(bcov_sp_q(:,k,wn)*ad_sp(3*nlev+1:4*nlev,wn))
! Check
!IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn)
     enddo
   enddo
!
   return
!
   end subroutine adjeigtospec
!===============================================================================
! for alpha-variable (HJS)
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine adjeigtospecens(l_zwn_nproc, ad_sp, ad_eigen_a)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
! 2. Variables
! 2-1. Input
!
   integer(int_kind), intent(in   ) :: l_zwn_nproc
   real(real_kind), dimension(nlevar, l_zwn_nproc), intent(inout) :: ad_sp
! 2-2. Input/Output
   real(real_kind), dimension(nlev+1, l_zwn_nproc), intent(inout) :: ad_eigen_a
! 2-3. Local
   integer(int_kind) :: k, wn
!=============================================================================
!=============================================================================
! C. Main body of adjoint process
!=============================================================================
!
   do wn = l_zwn_nproc, 1,-1
     ad_sp(2*nlev+1:3*nlev,wn) = ad_sp(2*nlev+1:3*nlev,wn)+ad_sp(1:nlev,wn)
     ad_sp(1:nlev,wn) = 0.d0
     ad_sp(2*nlev+1:3*nlev,wn) = ad_sp(2*nlev+1:3*nlev,wn)+ad_sp(nlev+1:2*nlev,wn)
     ad_sp(nlev+1:2*nlev,wn) = 0.d0
     ad_sp(2*nlev+1:3*nlev,wn) = ad_sp(2*nlev+1:3*nlev,wn)+ad_sp(3*nlev+1:4*nlev,wn)
     ad_sp(3*nlev+1:4*nlev,wn) = 0.d0
!       DO k  = 1,nlev
!          ad_eigen_a(k,wn) =      &
!             DSQRT(eigv_a(k,wn))* &
!             SUM(bcov_sp_a(:,k)*  &
!                 ad_sp(3*nlev+1:4*nlev,wn))
!          ! Check
!          !IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn)
!       END DO
     do k = 1,nlev+1
       ad_eigen_a(k,wn) = dsqrt(eigv_a(k,wn))*(sum(bcov_sp_a(1:nlev,k)*ad_sp(2*nlev+1:3*nlev,wn))+bcov_sp_a(nlev+1,k)*ad_sp(4*nlev+1,wn))
!          ! Check
!          !IF(par%isMasterProc) WRITE(iulog,*) 'ad_eigen = ',ad_eigen(k,wn)
     enddo
   enddo
!
   return
!
   end subroutine adjeigtospecens
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module costfunminmod
!-------------------------------------------------------------------------------
