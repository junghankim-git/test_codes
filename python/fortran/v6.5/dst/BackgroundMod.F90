!===============================================================================
!  MODULE: BackgroundMod
!
!> @brief
!> - Read and handle background ensemble
!>
!> @date 10JUN2013
!> - HJ SONG: Devise and code reading subroutines
!> @date 21JUL2013
!> - HJ SONG & JH KWUN: Devise and code transform routines
!> @date 03JUL2014
!> - IH KWON: Code revision for the B method
!> @date 12JAN2015
!> - JH KWUN: Code revision for regression coefficient and ModelToCtrl
!> @date 04FEB2016
!> - HJ SONG: Code revision for zonal average of error variance and regression
!>            coefficients and zonal-wavenumber-average of error covarinace
!> @date 30MAY2016
!> - HJ SONG: Modification for Hybrid 4DEnVAR
!> @date 09AUG2016
!> - JH KIM : I/O Refactoring with the KIO
!===============================================================================
!-------------------------------------------------------------------------------
   module backgroundmod
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
   use kiapsbase, only: int_kind=>kim_int_kind, real_kind=>kim_real8_kind, iulog=>kim_iu_log, len_path=>kim_string_len_path, len_arg=>kim_string_len_arg
   use dimensions, only: np, nlev, ne, nvar3d, nvar2d, nlevar, zwn, zwn2, zwn_nproc, zwn_iproc, zwn_a, zwn2_a, zwn_nproc_a, zwn_iproc_a, reg_nproc, reg_iproc, psi_nproc, psi_iproc, chi_nproc, chi_iproc, tps_nproc, tps_iproc, hum_nproc, hum_iproc, neig, neig_chi, neig_tps, neig_q, nelemd, nelem, nets, nete, nk, ntime, nt_itime ! itime of nature
   use kiapsgrid, only: elem, deriv, nups_l, nups_g, neps_l, neps_g, dofs, dofs_ep, convertep2up_hyb, convertup2ep_hyb, convertup2ep_hyb_l
   use kiapswave, only: zwn_blength, zwn_bstart, psi_blength, psi_bstart, chi_blength, chi_bstart, tps_blength, tps_bstart, hum_blength, hum_bstart, map_zwn, map_reg, map_psi, map_chi, map_tps, map_hum, genmap
   use kiapsreg, only: reg_blength, reg_bstart
   use timing, only: timingstart, timingstop
   use ctrltospecmod, only: ctrltospec_opt, lat_up, norm_fac, d_loc_h
   use kiapsparmpi, only: par_reduce
   use kiapsparallel, only: par=>kim_par, hybrid=>kim_hybrid, par_double_precision, par_character, par_integer, global_shared_buf, global_shared_sum, barrierpar, par_max, par_min, abortpar, kim_par
   use kio, only: kio_t, iniiosystem, finiosystem, openfile, closefile, getdimension, getvariable
   use globalnorms, only: wrap_repro_sum
   use physicalconstants, only: rvap=>kim_rvap, rgas=>kim_rair, ps0=>kim_p0, gravi=>kim_g, dd_pi, kim_rearth
!  USE BECOutputMod,      ONLY : InitBECOutput,  &
!                                WriteBECOutput, &
!                                FinBECOutput
   use bectuoutputmod, only: initbectuoutput, writebectuoutput, finbectuoutput
   use edge, only: edgebuffer_t, initedgebuffer, freeedgebuffer, edgevpack, edgevunpack
   use bndry, only: bndry_exchangev
! 2. Common variables
   implicit none
!
   private
   character(len=len_path), dimension(:), allocatable :: listfile ! sample file list
!CHARACTER(LEN=len_path),           &
!   PUBLIC :: natureFile,           & ! Nature file name
!             bgFile,               & ! Background file name
!             becFile                 ! Background error covariance file name
   character(len=len_path), dimension(:), allocatable, public :: bgfile ! background file name
   character(len=len_path), public :: naturefile, becfile ! background error covariance file name
   character(len=len_arg), public :: model_name, bmethod ! method for background error covariance
   real(real_kind), parameter :: tr = 270.d0
   integer(int_kind), public :: cv_opt_wind, cv_opt_hum, cv_opt_tu, cv_opt_psu, ens_opt_center
   real(real_kind), public :: beta1, infl_ens, distr_infl
   integer(int_kind), public :: levtrans
!===============================================================================
!===============================================================================
! B. Public modules and variables
!===============================================================================
! 1. Public modules
   public :: setreadback
!PUBLIC :: ReadBackEns
!PUBLIC :: ReadBackHyb
!PUBLIC :: ReadBackNMC
   public :: readbackh4dev
   public :: calphib
   public :: invlapcg
   public :: temptom
   public :: calnlbe
   public :: calpsb
   public :: calpsbnlbe
   public :: modeltoctrl
   public :: tlmtemptophi
!PUBLIC :: TlmTempToM
   public :: psitorotwind
   public :: adjustp
! 2. Public variables
!REAL(KIND=real_kind), DIMENSION(:,:,:,:), ALLOCATABLE, &
!   PUBLIC :: bg_md,                                    & ! Background
!             nt_org,                                   & ! Nature
!             errvar,                                   & ! variance of backgrouond error
!             z_org,                                    & ! Geopotential height for each layer
!             bg_md_q,                                  & ! Specific humidity for background
!             pmid,                                     & ! Middle level pressure
!             pint                                        ! Interface level pressure
!up
!REAL(KIND=real_kind), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, & ! extended to 4D
   real(real_kind), dimension(:,:,:,:), allocatable, public :: bg_smpl
   real(real_kind), dimension(:,:,:,:,:,:), allocatable, public :: bg_uv_rot ! rotational wind in background field
   real(real_kind), dimension(:,:,:,:,:), allocatable, public :: bg_md, pmid
   real(real_kind), dimension(:,:,:,:), allocatable, public :: nt_org, errvar, z_org, bg_md_q, pint ! interface level pressure
   real(real_kind), dimension(:,:,:), allocatable, public :: bcov_sp, bcov_sp_chi, bcov_sp_tps, bcov_sp_q, phis, zsfc, rl_org
   real(real_kind), dimension(:,:,:,:), allocatable, public :: u10m, v10m, t2m, q2m, tsfc, sfctype
   real(real_kind), dimension(:,:), allocatable, public :: eigv, eigv_chi, eigv_tps, eigv_q, eigv_a, bcov_sp_a, !up
!
   bg_eig,errvar_up
   real(real_kind), dimension(:), allocatable, public :: hyam, hybm, hyai, hybi, znu, znw, b_clim_v, b_ens_v
   real(real_kind), public :: cgtol_il, rescalepsi, rescalechi, rescalet, rescaleq, rescaleps
   real(real_kind), public :: b_clim, b_ens
   integer(int_kind), public :: nsmpl, lwork, niter_il
   character(len=len_path), public :: list_file, bgfilelist
   logical, public :: exst_bg, exst_ens, exst_bec, exst_nt, bottom2top
   logical, public :: vlocal
   logical, public :: is_moist_mr
   real(real_kind), public :: vlenscale
   real(real_kind), dimension(:,:,:), allocatable, public :: regress_up
   real(real_kind), dimension(:,:,:), allocatable, public :: pmid_up, pmid_ens_mean
   real(real_kind), dimension(:,:), allocatable, public :: t_ens_mean, q_ens_mean
   real(real_kind), dimension(:), allocatable, public :: zsfc_up, zsfc_ens
   integer(int_kind) :: ierr
!===============================================================================
!
   contains
!===============================================================================
!  SUBROUTINE: SetReadBack
!
!> @brief
!> - Set reading the names of background ensemble files
!>
!> @date 11JUN2013
!> - HJ SONG: Devise and code
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine setreadback
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
! 2. Local variables
!
   type(kio_t) :: io_t
   integer(int_kind) :: i, j, k, ie, ismpl, err, ierr, wn, itime
!=============================================================================
!=============================================================================
! B. Main Body (set vertical coordinates)
!=============================================================================
!.. Inquiry of file existence
!
   exst_bg = .false.
   if (len_trim(bgfilelist)/= 0) inquire(file = trim(bgfilelist),exist = exst_bg)
   exst_ens = .false.
   if (len_trim(list_file)/= 0) inquire(file = trim(list_file),exist = exst_ens)
   exst_bec = .false.
   if (len_trim(becfile)/= 0) inquire(file = trim(becfile),exist = exst_bec)
   call mpi_barrier(par%comm,ierr)
   if (par%ismasterproc) then
     print*,trim(model_name)," is selected for background"
     print*,trim(bmethod)," is selected for backgroun error estimateion"
     print*,trim(list_file)," is the list file"
     print*,'nsmpl:',nsmpl
   endif
   allocate(bgfile(ntime))
   if (par%ismasterproc) then
     open(unit = 11,file = trim(bgfilelist),form = 'FORMATTED',status = 'OLD')
     do itime = 1,ntime
       read(11,'(a) ') bgfile(itime)
     enddo
     close(11)
   endif
   call iniiosystem(io_t,par%comm,par%nprocs,par%iproc)
   if (bmethod=='h4dev'.or.bmethod=='3dvar'.or.bmethod=='4dclv') then !--for hybrid or nmc
     if (exst_bg) then
!.. Open a netcdf file of background
       call openfile(io_t,trim(bgfile(nt_itime)))
       if (par%ismasterproc) print*,'Open ',trim(bgfile(nt_itime)),' to set vertical coordinates'
     else
       allocate(listfile(nsmpl))
       open(unit = 10,file = trim(list_file),form = 'FORMATTED',status = 'OLD')
       do ismpl = 1,nsmpl
         read(10,'(a) ') listfile(ismpl)
       enddo
       close(10)
!.. Open a netcdf file of one of the samples
       call openfile(io_t,trim(listfile(1)))
       if (par%ismasterproc) print*,'Open ',trim(listfile(1)),' to set vertical coordinates'
     endif
   else !-----------for ens
     allocate(listfile(nsmpl))
     open(unit = 10,file = trim(list_file),form = 'FORMATTED',status = 'OLD')
     do ismpl = 1,nsmpl
       read(10,'(a) ') listfile(ismpl)
     enddo
     close(10)
!.. Open a netcdf file of Ensemblem sample
     call openfile(io_t,trim(listfile(1)))
     if (par%ismasterproc) print*,'Open ',trim(listfile(1)),' to set vertical coordinates'
   endif
!.. Read hyam and hybm
   allocate(hyam(nlev),hybm(nlev))
   call getvariable(io_t,'hyam',hyam)
   call getvariable(io_t,'hybm',hybm)
   bottom2top = .false.
   if (hybm(1)>hybm(nlev)) then
     bottom2top = .true.
   endif
   if (bottom2top) then
     hyam = hyam(nlev:1:-1)
     hybm = hybm(nlev:1:-1)
   endif
!.. Read hyai and hybi
   allocate(hyai(nlev+1),hybi(nlev+1))
   call getvariable(io_t,'hyai',hyai)
   call getvariable(io_t,'hybi',hybi)
   if (bottom2top) then
     hyai = hyai(nlev+1:1:-1)
     hybi = hybi(nlev+1:1:-1)
   endif
   if (model_name=='kim-sw') then
     allocate(znu(nlev),znw(nlev+1))
     do k = 1,nlev
       znu(k) = hyam(k)+hybm(k)
     enddo
     do k = 1,nlev+1
       znw(k) = hyai(k)+hybi(k)
     enddo
!      IF( par%ismasterproc )THEN
!         print*,"Vertical coord (Top to Bottom): hyam,hybm,eta *1000"
!         DO k = 1,nlev
!            print*,k,hyam(k),hybm(k),znu(k)*1000.
!         ENDDO
!      ENDIF
!   ELSE
!      IF( par%ismasterproc )THEN
!         print*,"Vertical coord (Top to Bottom): hyam,hybm,pressure (Pa)"
!         DO k = 1,nlev
!            print*,k,hyam(k),hybm(k),hyam(k)*ps0 + hybm(k)*1000.D2
!         ENDDO
!       ENDIF
   endif ! model_name
   call mpi_barrier(par%comm,ierr)
!.. Check zwn_nproc and zwn_iproc
   if (par%ismasterproc) print*,'Number of wave coefficients',zwn2
   print*,'proc,zwn_nproc,begin,end:',par%iproc,zwn_nproc,zwn_iproc+1,zwn_iproc+zwn_nproc
   if (par%ismasterproc) print*,'Number of wave coefficients for alpha',zwn2_a
   print*,'proc,zwn_nproc_a,begin,end:',par%iproc,zwn_nproc_a,zwn_iproc_a+1,zwn_iproc_a+zwn_nproc_a
   if (par%ismasterproc) print*,'End of setreadback'
   call closefile(io_t)
   call finiosystem(io_t)
!=============================================================================
   return
!
   end subroutine setreadback
!===============================================================================
!  SUBROUTINE ReadBackH4DEV
!
!> @brief
!> - Read a background ensemble, transform background states,
!> - and make error covariance for Hybrid 4D Ensemble Variational Code.
!>
!> @date 30MAY2016
!> - HJ SONG: First written
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine readbackh4dev
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
! 2. Local variables
!
   real(real_kind), dimension(:,:,:,:,:), allocatable :: bg_smpl_q, bg_smpl_p, bg_smpl_chi, bg_smpl_pint
   real(real_kind), dimension(:,:,:,:,:), allocatable :: !bg_smpl_mean, bg_smpl_q_mean, bg_smpl_p_mean, bg_smpl_pint_mean
   real(real_kind), dimension(:,:,:,:), allocatable :: errvar_md, bg_smpl_ep
   real(real_kind), dimension(:,:,:,:), allocatable :: bg_ct
   real(real_kind), dimension(:,:,:), allocatable :: bg_sp_smpl, reg_avg, bg_smpl_mean
   real(real_kind), dimension(:,:), allocatable :: bg_sp, bg_md_up
   real(real_kind), dimension(:), allocatable :: l_reg_avg, l_eigf_psi, l_eigf_chi, l_eigf_tps, l_eigf_hum
   real(real_kind) :: psfact, distr
   real(real_kind) :: tmpmax, tmpmin, maxv, minv
   integer(int_kind) :: i, j, k, l, ie, ig, ks, ismpl, wn, vr, iproc, vertnum, icnt, ierr, itime
!=============================================================================
!
   if (par%ismasterproc) print*,'Background error samples for nmc method'
!=============================================================================
!=============================================================================
! A. Read background
!=============================================================================
!.. Allocate  variables
   allocate(bg_md(np,np,nlevar,nelemd,ntime))
   allocate(bg_uv_rot(np,np,nlev,2,nelemd,ntime))
   allocate(bg_md_q(np,np,nlev,nelemd))
   allocate(pint(np,np,nlev+1,nelemd))
   allocate(pmid(np,np,nlev,nelemd,ntime))
   allocate(z_org(np,np,nlev,nelemd))
   allocate(zsfc(np,np,nelemd))
   allocate(phis(np,np,nelemd))
   allocate(rl_org(np,np,nelemd))
   allocate(u10m(np,np,nelemd,ntime))
   allocate(v10m(np,np,nelemd,ntime))
   allocate(t2m(np,np,nelemd,ntime))
   allocate(q2m(np,np,nelemd,ntime))
   allocate(tsfc(np,np,nelemd,ntime))
   allocate(sfctype(np,np,nelemd,ntime))
   bg_md = 0.d0
   bg_md_q = 0.d0
   if (par%ismasterproc) write(iulog,*) 'model_name:',model_name
   if (model_name=='kim-old') then ! read background for kim old version
     if (par%ismasterproc) write(iulog,*) 'Call readback_kimold'
   elseif (model_name=='kim-sw'.or.model_name=='kim-sh') then ! read background for kim-sw or-sh
     if (par%ismasterproc) write(iulog,*) 'Call readback_kim'
     do itime = 1,ntime
       call readback_kim(itime)
     enddo
   else
     if (par%ismasterproc) then
       stop "unrecognized model name"
     endif
   endif
   if (par%ismasterproc) write(iulog,*) 'End of reading background file'
   if (par%ismasterproc) write(iulog,*) 'bg_md(1,1,:,1,1) ',bg_md(1,1,:,1,1)
!=============================================================================
!=============================================================================
! B. Read background error covariance or generate it
!=============================================================================
! 1. Allocate
   allocate(errvar_md(np,np,nlevar,nelemd))
   allocate(errvar(np,np,nlevar,nelemd))
   allocate(errvar_up(nups_l,nlevar))
   allocate(eigv(neig,zwn_nproc))
   allocate(bcov_sp(nlev,neig,zwn_nproc))
   allocate(eigv_chi(neig_chi,zwn_nproc))
   allocate(bcov_sp_chi(nlev,neig_chi,zwn_nproc))
   allocate(eigv_tps(neig_tps,zwn_nproc))
   allocate(bcov_sp_tps(nlev+1,neig_tps,zwn_nproc))
   allocate(eigv_q(neig_q,zwn_nproc))
   allocate(bcov_sp_q(nlev,neig_q,zwn_nproc))
   errvar_md = 0.d0
   errvar = 0.d0
   errvar_up = 0.d0
   eigv = 0.d0
   eigv_chi = 0.d0
   eigv_tps = 0.d0
   eigv_q = 0.d0
   bcov_sp = 0.d0
   bcov_sp_chi = 0.d0
   bcov_sp_tps = 0.d0
   bcov_sp_q = 0.d0
   if (cv_opt_tu==1.or.cv_opt_tu==2.or.cv_opt_tu==3) then
     allocate(regress_up(nlev,nups_l,2*nlev+1))
     allocate(l_reg_avg(reg_nproc))
     allocate(l_eigf_psi(psi_nproc),l_eigf_chi(chi_nproc),l_eigf_tps(tps_nproc),l_eigf_hum(hum_nproc))
   endif
!=============================================================================
! B-1  Read background error covariance from file
   if (exst_bec) then
!=============================================================================
!   Read background error covariance
     call readbackerrcov(l_reg_avg,l_eigf_psi,l_eigf_chi,l_eigf_tps,l_eigf_hum,errvar_md)
! Additional Output: errvar,eigv,bcov_sp,Regress,errvar_md
!up
     call convertep2up_hyb(nlevar,errvar,errvar_up)
!.. Pop the locally distributed variables
     call regresspop(l_reg_avg,regress_up)
     call eigfpop(l_eigf_psi,l_eigf_chi,l_eigf_tps,l_eigf_hum)
! bcov_sp,    &
! bcov_sp_chi,&
! bcov_sp_Tps,&
! bcov_sp_hum)
     if (.not.allocated(listfile)) allocate(listfile(nsmpl*ntime))
     if (par%ismasterproc) then
       open(unit = 11,file = trim(list_file),form = 'FORMATTED',status = 'OLD')
       do ismpl = 1,nsmpl*ntime
         read(11,'(a) ') listfile(ismpl)
       enddo
       close(11)
     endif
     call mpi_barrier(par%comm,ierr)
     do ismpl = 1,nsmpl*ntime
       call mpi_bcast(listfile(ismpl),256,par_character,par%root,par%comm,ierr)
     enddo
!up
!ALLOCATE( bg_smpl(np,np,nlevar,nelemd,ntime,nsmpl) )
     allocate(bg_smpl(nups_l,nlevar,ntime,nsmpl))
     bg_smpl = 0.d0
     if (bmethod=='h4dev') then
! up
!ALLOCATE( bg_smpl_mean(np,np,nlevar,nelemd,ntime) )
       allocate(bg_smpl_mean(nups_l,nlevar,ntime))
! 1. Read samples (PIO input) and calculate mean field -----------
       call readsamples(nsmpl,listfile,bg_smpl) ! out
       deallocate(listfile)
! 2. Perturbation matrix ---------------------------------------------
       bg_smpl_mean = 0.d0
       if (par%ismasterproc) print*,'ens_opt_center = ',ens_opt_center
       if (ens_opt_center.eq.0) then
         do itime = 1, ntime
!up
!DO ie = 1, nelemd
           do k = 1, nlevar
!DO j  = 1, np
!DO i  = 1, np
             do i = 1, nups_l
!.. Calculate Ensemble mean of u, v, T, humidity (q or rh), and ps
               do ismpl = 1, nsmpl
!bg_smpl_mean(i,j,k,ie,itime) = bg_smpl_mean(i,j,k,ie,itime)  &
!                             + bg_smpl(i,j,k,ie,itime,ismpl)
                 bg_smpl_mean(i,k,itime) = bg_smpl_mean(i,k,itime)+bg_smpl(i,k,itime,ismpl)
               enddo
!bg_smpl_mean(i,j,k,ie,itime) = bg_smpl_mean(i,j,k,ie,itime) /DBLE(nsmpl)
               bg_smpl_mean(i,k,itime) = bg_smpl_mean(i,k,itime)/dble(nsmpl)
!END DO
!END DO
             enddo
           enddo
         enddo ! itime
       elseif (ens_opt_center.eq.1) then
!up
         allocate(bg_md_up(nups_l,nlevar))
         do itime = 1,ntime
           call convertep2up_hyb(nlevar,bg_md(:,:,:,:,itime),bg_md_up)
           bg_smpl_mean(:,:,itime) = bg_md_up
         enddo
       else
         if (par%ismasterproc) print*,'Wrong option for the center of ensemble samples'
       endif
!.. Calculate perturbations,delta_X (bg_smpl)
       do itime = 1,ntime
!DO ie = 1, nelemd
         do k = 1, nlevar
!DO j  = 1, np
!DO i  = 1, np
           do i = 1, nups_l
             do ismpl = 1, nsmpl
!bg_smpl(i,j,k,ie,itime,ismpl) = bg_smpl(i,j,k,ie,itime,ismpl) &
!                              - bg_smpl_mean(i,j,k,ie,itime)
               bg_smpl(i,k,itime,ismpl) = bg_smpl(i,k,itime,ismpl)-bg_smpl_mean(i,k,itime)
             enddo
!END DO
!END DO
           enddo
         enddo
       enddo ! itime
       if (par%ismasterproc) write(iulog,*) "inflation factor for ensemble:",infl_ens
!up
!IF (par%ismasterproc) WRITE(iulog,*)  "bg_smpl_mean(1,1,:,1,:):",bg_smpl_mean(1,1,:,1,:)
       if (par%ismasterproc) write(iulog,*) "bg_smpl_mean(1,:,:):",bg_smpl_mean(1,:,:)
       deallocate(bg_smpl_mean)
!.. Model space (U,V,T,Q) -> controal varible space
       allocate(bg_smpl_ep(np,np,nlevar,nelemd))
       allocate(bg_ct(np,np,nlevar,nelemd))
       bg_ct = 0.d0
       do ismpl = 1,nsmpl
         do itime = 1, ntime
           if (par%ismasterproc) write(iulog,'(a,i4,a) ') 'To control variable(ismpl = ',ismpl,') '
!CALL ModelToCtrl(bg_smpl(:,:,:,:,itime,ismpl),    &
           call convertup2ep_hyb(nlevar,bg_smpl(:,:,itime,ismpl),bg_smpl_ep)
           call modeltoctrl(bg_smpl_ep(:,:,:,:),bg_md(:,:,:,:,itime),bg_uv_rot(:,:,:,:,:,itime),bg_ct)
!up
!bg_smpl(:,:,:,:,itime,ismpl) = bg_ct(:,:,:,:)
           call convertep2up_hyb(nlevar,bg_ct(:,:,:,:),bg_smpl(:,:,itime,ismpl))
         enddo ! itime
       enddo ! ismpl
       deallocate(bg_ct)
       deallocate(bg_smpl_ep)
       bg_smpl = bg_smpl/dsqrt(dble(nsmpl-1))
       bg_smpl = infl_ens*bg_smpl
     endif ! bmethod=='3dvar'
     call eigalpha()
   else
     if (par%ismasterproc) write(iulog,*) "you need to generate climatological b."
   endif !! exst_bec
!=============================================================================
   call mpi_barrier(par%comm,ierr)
!   print*,par%iproc,'bg_md',bg_md(1,1,:,1)
!=============================================================================
! C. Rescaling variance of background error
!=============================================================================
!.. psi
   if (rescalepsi/= 1.d0) then
     do ie = 1, nelemd
       errvar(:,:,1:nlev*1,ie) = errvar(:,:,1:nlev*1,ie)*rescalepsi
     enddo
     if (errvar_md(1, 1, 1, 1)/= 0.d0) then
       do ie = 1, nelemd
         errvar_md(:,:,1:nlev*1,ie) = errvar_md(:,:,1:nlev*1,ie)*rescalepsi
       enddo
     endif
     if (par%ismasterproc) write(iulog,*) 'Psi variance was rescaled,ratio = ',rescalepsi
   endif
!.. Unbalanced chi
   if (rescalechi/= 1.d0) then
     do ie = 1, nelemd
       errvar(:,:,nlev+1:nlev*2,ie) = errvar(:,:,nlev+1:nlev*2,ie)*rescalechi
     enddo
     if (errvar_md(1, 1, 1, 1)/= 0.d0) then
       do ie = 1, nelemd
         errvar_md(:,:,nlev+1:nlev*2,ie) = errvar_md(:,:,nlev+1:nlev*2,ie)*rescalechi
       enddo
     endif
     if (par%ismasterproc) write(iulog,*) 'Chi variance was rescaled,ratio = ',rescalechi
   endif
!.. Unbalanced T
   if (rescalet/= 1.d0) then
     do ie = 1, nelemd
       errvar(:,:,nlev*2+1:nlev*3,ie) = errvar(:,:,nlev*2+1:nlev*3,ie)*rescalet
     enddo
     if (errvar_md(1, 1, 1, 1)/= 0.d0) then
       do ie = 1, nelemd
         errvar_md(:,:,nlev*2+1:nlev*3,ie) = errvar_md(:,:,nlev*2+1:nlev*3,ie)*rescalet
       enddo
     endif
     if (par%ismasterproc) write(iulog,*) 'T variance was rescaled,ratio = ',rescalet
   endif
!.. Unbalanced Q
   if (rescaleq/= 1.d0) then
     do ie = 1, nelemd
       errvar(:,:,nlev*3+1:nlev*4,ie) = errvar(:,:,nlev*3+1:nlev*4,ie)*rescaleq
     enddo
     if (errvar_md(1, 1, 1, 1)/= 0.d0) then
       do ie = 1, nelemd
         errvar_md(:,:,nlev*3+1:nlev*4,ie) = errvar_md(:,:,nlev*3+1:nlev*4,ie)*rescaleq
       enddo
     endif
     if (par%ismasterproc) write(iulog,*) 'Q variance was rescaled,ratio = ',rescaleq
   endif
!.. Unbalanced Ps
   if (rescaleps/= 1.d0) then
     do ie = 1, nelemd
       errvar(:,:,nlev*4+1,ie) = errvar(:,:,nlev*4+1,ie)*rescaleps
     enddo
     if (errvar_md(1, 1, 1, 1)/= 0.d0) then
       do ie = 1, nelemd
         errvar_md(:,:,nlev*4+1,ie) = errvar_md(:,:,nlev*4+1,ie)*rescaleps
       enddo
     endif
     if (par%ismasterproc) write(iulog,*) 'Ps variance was rescaled,ratio = ',rescaleps
   endif
!.. Calculate vertically-varying weight for climatological and ensemble B
!   levTrans: default = 25
   allocate(b_ens_v(nlevar))
   allocate(b_clim_v(nlevar))
   do k = 1,levtrans
     distr = dble(levtrans+1-k)/dble(levtrans)*distr_infl
     b_ens_v(k) =(1.d0/exp(distr**2.d0))*b_ens
   enddo
   do k = levtrans+1,nlev
     b_ens_v(k) = b_ens
   enddo
   b_ens_v(nlev+1:2*nlev) = b_ens_v(1:nlev)
   b_ens_v(2*nlev+1:3*nlev) = b_ens_v(1:nlev)
   b_ens_v(3*nlev+1:4*nlev) = b_ens_v(1:nlev)
   b_ens_v(4*nlev+1) = b_ens
   b_clim_v =(b_clim+b_ens)-b_ens_v
   if (par%ismasterproc) write(iulog,*) 'b_clim:',b_clim
   if (par%ismasterproc) write(iulog,*) 'b_ens_v:',b_ens_v
   if (par%ismasterproc) write(iulog,*) 'b_clim_v:',b_clim_v
!=============================================================================
! D. Read Nature file (for verification of OSSE test)
!=============================================================================
   exst_nt = .false.
   if (len_trim(naturefile)/= 0) inquire(file = trim(naturefile),exist = exst_nt)
   if (exst_nt) then
     allocate(nt_org(np,np,nlevar,nelemd))
     nt_org = 0.d0
     call readnature
   endif
!=============================================================================
! E. Check background error
!=============================================================================
   if (exst_nt.and.errvar_md(1,1,1,1)/= 0.d0) then
     global_shared_buf(:,2:11) = 0.d0
     do ie = 1,nelemd
!   Calculate variance between background and nature (From level 11 to bottom)
       global_shared_buf(ie,2) = global_shared_buf(ie,2)+sum((bg_md(:,:,11:nlev,ie,nt_itime)-nt_org(:,:,11:nlev,ie))**2.d0)
       global_shared_buf(ie,3) = global_shared_buf(ie,3)+sum((bg_md(:,:,11+nlev*1:nlev*2,ie,nt_itime)-nt_org(:,:,11+nlev*1:nlev*2,ie))**2.d0)
       global_shared_buf(ie,4) = global_shared_buf(ie,4)+sum((bg_md(:,:,11+nlev*2:nlev*3,ie,nt_itime)-nt_org(:,:,11+nlev*2:nlev*3,ie))**2.d0)
       global_shared_buf(ie,5) = global_shared_buf(ie,5)+sum((bg_md(:,:,11+nlev*3:nlev*4,ie,nt_itime)-nt_org(:,:,11+nlev*3:nlev*4,ie))**2.d0)
       global_shared_buf(ie,6) = global_shared_buf(ie,6)+sum((bg_md(:,:,1+nlev*4,ie,nt_itime)-nt_org(:,:,1+nlev*4,ie))**2.d0)
!   Calculate Variance of Background Error covariance (From level 11 to bottom)
       global_shared_buf(ie,7) = global_shared_buf(ie,7)+sum(errvar_md(:,:,11:nlev,ie))
       global_shared_buf(ie,8) = global_shared_buf(ie,8)+sum(errvar_md(:,:,11+nlev*1:nlev*2,ie))
       global_shared_buf(ie,9) = global_shared_buf(ie,9)+sum(errvar_md(:,:,11+nlev*2:nlev*3,ie))
       global_shared_buf(ie,10) = global_shared_buf(ie,10)+sum(errvar_md(:,:,11+nlev*3:nlev*4,ie))
       global_shared_buf(ie,11) = global_shared_buf(ie,11)+sum(errvar_md(:,:,1+nlev*4,ie))
     enddo
     call mpi_barrier(par%comm,ierr)
     call wrap_repro_sum(nvars = 11,comm = hybrid%par%comm)
     if (par%ismasterproc) write(iulog,*) 'Variance between background and nature'
     if (par%ismasterproc) write(iulog,*) 'vari_bg_u = ',global_shared_sum(2)/dble(np*np*nelem*nlev)
     if (par%ismasterproc) write(iulog,*) 'vari_bg_v = ',global_shared_sum(3)/dble(np*np*nelem*nlev)
     if (par%ismasterproc) write(iulog,*) 'vari_bg_T = ',global_shared_sum(4)/dble(np*np*nelem*nlev)
     if (par%ismasterproc) write(iulog,*) 'vari_bg_q = ',global_shared_sum(5)/dble(np*np*nelem*nlev)
     if (par%ismasterproc) write(iulog,*) 'vari_bg_ps = ',global_shared_sum(6)/dble(np*np*nelem*1)
     if (par%ismasterproc) write(iulog,*) 'Variance of background error covariance'
     if (par%ismasterproc) write(iulog,*) 'vari_B_u = ',global_shared_sum(7)/dble(np*np*nelem*nlev)
     if (par%ismasterproc) write(iulog,*) 'vari_B_v = ',global_shared_sum(8)/dble(np*np*nelem*nlev)
     if (par%ismasterproc) write(iulog,*) 'vari_B_T = ',global_shared_sum(9)/dble(np*np*nelem*nlev)
     if (par%ismasterproc) write(iulog,*) 'vari_B_q = ',global_shared_sum(10)/dble(np*np*nelem*nlev)
     if (par%ismasterproc) write(iulog,*) 'vari_B_ps = ',global_shared_sum(11)/dble(np*np*nelem*1)
   endif ! exst_nt
!=============================================================================
   if (par%ismasterproc) write(iulog,*) 'exst_nt = ',exst_nt
   deallocate(errvar_md)
   deallocate(errvar)
   deallocate(l_reg_avg)
   deallocate(l_eigf_psi,l_eigf_chi,l_eigf_tps,l_eigf_hum)
   if (par%ismasterproc) write(iulog,*) 'ReadBackH4DEV end'
!
   end subroutine readbackh4dev
!===============================================================================
!===============================================================================
!
! SUBROUTINE EigAlpha
!
!> @brief
!> - make localization matix for alpha using
!> - Gaussian filter for horizon (Kuhl et al. 2013) and
!> - Gaspari and Cohn's function for vertical
!>
!> @date 30MAY2016
!> - HJ SONG: Initial version
!> @date 05SEP2016
!> - HJ SONG: Modification to link ps and column varables
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine eigalpha()
! Local variables
!
   real(real_kind), dimension(nlev+1) :: pmid_avg, tmp
   real(real_kind), dimension(nlev+1, nlev+1) :: vloc
   real(real_kind), dimension(:), allocatable :: work
   real(real_kind) :: dist, dist_zero2, ratio
   integer(int_kind) :: i, j, ie, l, k, wn, info
   integer(int_kind) :: ierr
!
   do k = 1, nlev
     global_shared_buf(:,1) = 0.d0
     do ie = 1,nelemd
       global_shared_buf(ie,1) = global_shared_buf(ie,1)+sum(pmid(:,:,k,ie,1))
     enddo
     call mpi_barrier(par%comm,ierr)
     call wrap_repro_sum(nvars = 1,comm = hybrid%par%comm)
     pmid_avg(k) = global_shared_sum(1)/dble(np*np*nelem)
   enddo
!
   global_shared_buf(:,1) = 0.d0
   do ie = 1,nelemd
     global_shared_buf(ie,1) = global_shared_buf(ie,1)+sum(bg_md(:,:,4*nlev+1,ie,1))
   enddo
   call mpi_barrier(par%comm,ierr)
   call wrap_repro_sum(nvars = 1,comm = hybrid%par%comm)
   pmid_avg(nlev+1) = global_shared_sum(1)/dble(np*np*nelem)
   if (par%ismasterproc) write(iulog,*) 'averaged pmid = ',pmid_avg(:)
!   Make localization matrix
   vloc = 0.d0
!   vlenScale= 0.2D0 (default)
   dist_zero2 = dsqrt(10.d0/3.d0)*vlenscale
   do l = 1,nlev+1
     do k = 1, nlev+1
       dist = dabs(dlog(pmid_avg(k))-dlog(pmid_avg(l)))
       ratio = dist/dist_zero2
       if (dist.le.dist_zero2) then
         vloc(k,l) =-1.d0/4.d0*(ratio**5.d0)+1.d0/2.d0*(ratio**4.d0)+5.d0/8.d0*(ratio**3.d0)-5.d0/3.d0*(ratio**2.d0)+1.d0
       else
         if (dist.le.2.d0*dist_zero2) then
           vloc(k,l) = 1.d0/12.d0*(ratio**5.d0)-1.d0/2.d0*(ratio**4.d0)+5.d0/8.d0*(ratio**3.d0)+5.d0/3.d0*(ratio**2.d0)-5.d0*ratio+4.d0-2.d0/3.d0/ratio
         else
           vloc(k,l) = 0.d0
         endif
       endif
     enddo
   enddo
   if (par%ismasterproc) then
     print*,'vloc'
     do k = 1,nlev+1
       write(111,*),(vloc(k,l),l = 1,nlev+1)
     enddo
     print*,' '
   endif
   call mpi_barrier(par%comm,ierr)
!   Eigendecomposition and make eigen-pairs for alpha
   allocate(bcov_sp_a(nlev+1,nlev+1))
   allocate(eigv_a(nlev+1,zwn_nproc_a))
   allocate(work(lwork))
   call dsyev('V','L',nlev+1,vloc,nlev+1,tmp,work,lwork,info)
   do k = 1,nlev+1
     bcov_sp_a(:,k) = vloc(:,nlev+1-k+1)
     do wn = 1,zwn_nproc_a
       eigv_a(k,wn) = tmp(nlev+1-k+1)*norm_fac*dexp(-dint(dsqrt(dble(wn+zwn_iproc_a-1)))**2.d0/d_loc_h**2.d0)
     enddo
   enddo
   if (par%ismasterproc) then
     do k = 1, nlev+1
       write(iulog,*) 'eigNum,eigv_a(eignum,1):',k,eigv_a(k,1)
     enddo
   endif
   call mpi_barrier(par%comm,ierr)
   return
!
   end subroutine eigalpha
!========================================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine zonavg(errvar) ! inout
!-------------------------------------------------------------------------------
   include 'mpif.h'
!
   real(real_kind), dimension(np, np, nlevar, nelemd), intent(inout) :: errvar
   integer(int_kind), dimension(:), allocatable :: n_sum, n_avg
   real(real_kind), dimension(:,:), allocatable :: errvar_sum, errvar_avg, errvar_up
   integer(int_kind) :: i, j, ie, ilat, k, ierr, iup
!
   allocate(errvar_up(nups_l,nlevar))
   ierr = timingstart('EP2UP_fromZonAvg')
   call convertep2up_hyb(nlevar,errvar,errvar_up)
   ierr = timingstop('EP2UP_fromZonAvg')
   allocate(errvar_sum(nlevar,180),n_sum(180))
   errvar_sum = 0.d0
   n_sum = 0
! Conduct zonal-averaging of variance
   do iup = 1,nups_l
     do ilat = 1, 180
       if (lat_up(iup, 1).ge.dble(ilat-91)/180.d0*dd_pi.and.lat_up(iup, 1).lt.dble(ilat-90)/180.d0*dd_pi) then
         errvar_sum(:,ilat) = errvar_sum(:,ilat)+errvar_up(iup,:)
         n_sum(ilat) = n_sum(ilat)+1
       endif
     enddo
     if (lat_up(iup, 1).eq.dble(180-90)/180.d0*dd_pi) then
       errvar_sum(:,180) = errvar_sum(:,180)+errvar_up(iup,:)
       n_sum(180) = n_sum(180)+1
     endif
   enddo
   call mpi_barrier(par%comm,ierr)
   allocate(errvar_avg(nlevar,180),n_avg(180))
   errvar_avg = 0.d0
   n_avg = 0
   call mpi_allreduce(errvar_sum,errvar_avg,nlevar*180,par_double_precision,mpi_sum,par%comm,ierr)
   call mpi_allreduce(n_sum,n_avg,180,par_integer,mpi_sum,par%comm,ierr)
   do ilat = 1,180
     if (par%ismasterproc) then
       print*,'iLat = ',ilat
       print*,'errvar_sum(1,ilat) = ',errvar_sum(1,ilat)
       print*,'n_sum(ilat) = ',n_sum(ilat)
       print*,'errvar_avg(1,ilat) = ',errvar_avg(1,ilat)
       print*,'n_avg(ilat) = ',n_avg(ilat)
     endif
     if (n_avg(ilat).gt.0) then
       errvar_avg(:,ilat) = errvar_avg(:,ilat)/dble(n_avg(ilat))
     endif
     if (par%ismasterproc) then
       print*,'errvar_avg(1,ilat) ',errvar_avg(1,ilat)
     endif
   enddo
   call mpi_barrier(par%comm,ierr)
   deallocate(errvar_sum,n_sum)
   deallocate(n_avg)
   errvar_up = 0.d0
   do iup = 1,nups_l
     do ilat = 1, 180
       if (lat_up(iup, 1).ge.dble(ilat-91)/180.d0*dd_pi.and.lat_up(iup, 1).lt.dble(ilat-90)/180.d0*dd_pi) then
         errvar_up(iup,:) = errvar_avg(:,ilat)
       endif
     enddo
     if (lat_up(iup, 1).eq.dble(180-90)/180.d0*dd_pi) then
       errvar_up(iup,:) = errvar_avg(:,180)
     endif
     do k = 1, nlevar
       if (errvar_up(iup, k).eq.0.d0) then
         print*,'errvar zero at iup,k = ',iup,k
       endif
     enddo
   enddo
   call mpi_barrier(par%comm,ierr)
   deallocate(errvar_avg)
   ierr = timingstart('UP2EP_fromZonAvg')
   call convertup2ep_hyb(nlevar,errvar_up,errvar)
   ierr = timingstop('UP2EP_fromZonAvg')
   deallocate(errvar_up)
!
   end subroutine zonavg
!==================================================================================
!====================================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine eigfpop(l_eigf_psi, l_eigf_chi, l_eigf_tps, l_eigf_hum) ! in
!bcov_sp,     & ! out
!bcov_sp_chi, & ! out
!bcov_sp_Tps, & ! out
!bcov_sp_q    & ! out
!-------------------------------------------------------------------------------
   include 'mpif.h'
! 1. Input and Output
!
   real(real_kind), dimension(psi_nproc), intent(in   ) :: l_eigf_psi
   real(real_kind), dimension(chi_nproc), intent(in   ) :: l_eigf_chi
   real(real_kind), dimension(tps_nproc), intent(in   ) :: l_eigf_tps
   real(real_kind), dimension(hum_nproc), intent(in   ) :: l_eigf_hum
! 2. Local variables
   real(real_kind), dimension(:,:,:), allocatable :: errcov_sp_psi, errcov_sp_chi, errcov_sp_tps, errcov_sp_hum
   real(real_kind), dimension(:,:,:), allocatable :: eigf_swn_psi, eigf_swn_chi, eigf_swn_tps, eigf_swn_hum
   real(real_kind), dimension(:), allocatable :: eigf_psi1d, eigf_chi1d, eigf_tps1d, eigf_hum1d
   integer(int_kind) :: i, k, l, ismpl, wn, swn, iswn, iproc, icnt, info, ierr
! -> Global variable
!
   allocate(eigf_psi1d(nlev*neig*zwn))
   allocate(eigf_chi1d(nlev*neig_chi*zwn))
   allocate(eigf_tps1d((nlev+1)*neig_tps*zwn))
   allocate(eigf_hum1d(nlev*neig_q*zwn))
   call mpi_allgatherv(l_eigf_psi,psi_nproc,par_double_precision,eigf_psi1d,psi_blength,psi_bstart,par_double_precision,par%comm,ierr)
   call mpi_allgatherv(l_eigf_chi,chi_nproc,par_double_precision,eigf_chi1d,chi_blength,chi_bstart,par_double_precision,par%comm,ierr)
   call mpi_allgatherv(l_eigf_tps,tps_nproc,par_double_precision,eigf_tps1d,tps_blength,tps_bstart,par_double_precision,par%comm,ierr)
   call mpi_allgatherv(l_eigf_hum,hum_nproc,par_double_precision,eigf_hum1d,hum_blength,hum_bstart,par_double_precision,par%comm,ierr)
! -> Spherical wavenumber structure
   allocate(eigf_swn_psi(nlev,neig,zwn),eigf_swn_chi(nlev,neig_chi,zwn),eigf_swn_tps(nlev+1,neig_tps,zwn),eigf_swn_hum(nlev,neig_q,zwn))
   eigf_swn_psi = 0.d0
   eigf_swn_chi = 0.d0
   eigf_swn_tps = 0.d0
   eigf_swn_hum = 0.d0
   i = 0
   do iswn = 1,zwn
     do l = 1, neig
       do k = 1, nlev
         i = i+1
         eigf_swn_psi(k,l,iswn) = eigf_psi1d(i)
       enddo
     enddo
   enddo
   i = 0
   do iswn = 1,zwn
     do l = 1, neig_chi
       do k = 1, nlev
         i = i+1
         eigf_swn_chi(k,l,iswn) = eigf_chi1d(i)
       enddo
     enddo
   enddo
   i = 0
   do iswn = 1,zwn
     do l = 1, neig_tps
       do k = 1, nlev+1
         i = i+1
         eigf_swn_tps(k,l,iswn) = eigf_tps1d(i)
       enddo
     enddo
   enddo
   i = 0
   do iswn = 1,zwn
     do l = 1, neig_q
       do k = 1, nlev
         i = i+1
         eigf_swn_hum(k,l,iswn) = eigf_hum1d(i)
       enddo
     enddo
   enddo
   deallocate(eigf_psi1d)
   deallocate(eigf_chi1d)
   deallocate(eigf_tps1d)
   deallocate(eigf_hum1d)
! -> Distributed wavenumber space
   do wn = 1,zwn_nproc
     swn = dint(dsqrt(dble(wn+zwn_iproc-1))) ! spherical wavenumber
     bcov_sp(:,:,wn) = eigf_swn_psi(:,:,swn+1)
     bcov_sp_chi(:,:,wn) = eigf_swn_chi(:,:,swn+1)
     bcov_sp_tps(:,:,wn) = eigf_swn_tps(:,:,swn+1)
     bcov_sp_q(:,:,wn) = eigf_swn_hum(:,:,swn+1)
   enddo
   deallocate(eigf_swn_psi)
   deallocate(eigf_swn_chi)
   deallocate(eigf_swn_tps)
   deallocate(eigf_swn_hum)
   return
!
   end subroutine eigfpop
!========================================================================================
!========================================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine regressloc(reg_avg, l_reg_avg)
!
   real(real_kind), dimension(nlev, 2*nlev+1, 180), intent(in   ) :: reg_avg
   real(real_kind), dimension(reg_nproc), intent(  out) :: l_reg_avg
   real(real_kind), dimension(:), allocatable :: reg_tmp
   integer(int_kind) :: i, ilat, ireg, k
! -> local variable
!
   allocate(reg_tmp(nlev*(2*nlev+1)*180))
   i = 0
   do ilat = 1,180
     do ireg = 1, 2*nlev+1
       do k = 1, nlev
         i = i+1
         reg_tmp(i) = reg_avg(k,ireg,ilat)
       enddo
     enddo
   enddo
   l_reg_avg(1:reg_nproc) = reg_tmp(reg_iproc+1:reg_iproc+reg_nproc)
   deallocate(reg_tmp)
!
   end subroutine regressloc
!==================================================================================
!==================================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine regresspop(l_reg_avg, regress_up)
!-------------------------------------------------------------------------------
   include 'mpif.h'
!
   real(real_kind), dimension(reg_nproc), intent(in   ) :: l_reg_avg
   real(real_kind), dimension(nlev, nups_l, nk), intent(  out) :: regress_up
   real(real_kind), dimension(:), allocatable :: reg_tmp
   real(real_kind), dimension(:,:,:), allocatable :: lat, reg_avg
   integer(int_kind) :: i, j, ie, ilat, ireg, k, ierr, iup
! -> global variable
!
   allocate(reg_tmp(nlev*nk*180))
   call mpi_allgatherv(l_reg_avg,reg_nproc,par_double_precision,reg_tmp,reg_blength,reg_bstart,par_double_precision,par%comm,ierr)
   allocate(reg_avg(nlev,nk,180))
   i = 0
   do ilat = 1,180
     do ireg = 1, nk
       do k = 1, nlev
         i = i+1
         reg_avg(k,ireg,ilat) = reg_tmp(i)
       enddo
     enddo
   enddo
   deallocate(reg_tmp)
   regress_up = 0.d0
   do iup = 1,nups_l
     do ilat = 1, 180
       if (lat_up(iup, 1).ge.dble(ilat-91)/180.d0*dd_pi.and.lat_up(iup, 1).lt.dble(ilat-90)/180.d0*dd_pi) then
         regress_up(:,iup,:) = reg_avg(:,:,ilat)
       endif
     enddo
     if (lat_up(iup, 1).eq.dble(180-90)/180.d0*dd_pi) then
       regress_up(:,iup,:) = reg_avg(:,:,180)
     endif
   enddo
   call mpi_barrier(par%comm,ierr)
!
   end subroutine regresspop
!==========================================================================================
!==========================================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine reddof_miss(x, x_reduced)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Variables
! 1-1. Input
!
   real(real_kind), dimension(np, np, nelemd), intent(in   ) :: x
! 1-2. Output
   real(real_kind), dimension(np, np, nelemd), intent(  out) :: x_reduced
! 1-3. Local
   integer(int_kind) :: l_ie, ie, ist, ien, jst, jen, i, j
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
! 1. Zero out duplicated points
!
   x_reduced =-999.d0
   do l_ie = nets,nete
     ie = elem(l_ie)%globalid
     jst = 1
     jen = np-1
     ist = 1
     ien = np-1
     if (((ne*ne*2+1).le.ie).and.(ie.le.(ne*ne*4))) then
       jst = 2
       jen = np
     endif
     if ((ne*ne*5+1).le.ie) then
       ist = 2
       ien = np
     endif
     x_reduced(ist:ien,jst:jen,l_ie) = x(ist:ien,jst:jen,l_ie)
     if (ie.eq.(ne*(ne-1)+1)) x_reduced(1,np,l_ie) = x(1,np,l_ie)
     if (ie.eq.ne*(ne+1)) x_reduced(np,1,l_ie) = x(np,1,l_ie)
   enddo
!=============================================================================
!
   end subroutine reddof_miss
!=======================================================================================
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine readsamples(nsmpl, infilelist, bg_smpl)
!===============================================================================
!-------------------------------------------------------------------------------
   use netcdf
   use vertintmod, only: vertint
   use obscorrmod, only: sfcprscorr
!
   integer(int_kind), intent(in   ) :: nsmpl
   character(len=len_path), dimension(nsmpl*ntime), intent(in   ) :: infilelist
!up
!REAL(KIND=real_kind),DIMENSION(np,np,nlevar,nelemd,ntime,nsmpl), &
   real(real_kind), dimension(nups_l, nlevar, ntime, nsmpl), intent(  out) :: bg_smpl
! Local variables
!up
!REAL(KIND=real_kind), DIMENSION(:,:,:), ALLOCATABLE &
   real(real_kind), dimension(:), allocatable :: ps_tmp
!REAL(KIND=real_kind), DIMENSION(:,:,:,:), ALLOCATABLE &
   real(real_kind), dimension(:,:), allocatable :: pmid_ens, bg_smpl_tmp
   real(real_kind), dimension(:), allocatable :: x3d, x3di, x2d
   real(real_kind), dimension(:,:,:,:), allocatable :: zsfc_ep
!   REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE &
!       :: zsfc_up, zsfc_ens
   real(real_kind), dimension(:,:), allocatable ! :: pmid_up, bg_md_up
   :: bg_md_up
   real(real_kind), dimension(:,:,:,:), allocatable :: tq_smpl_tmp
   real(real_kind) :: psfact
   integer(int_kind) :: i, j, k, ie, ic, icnt, itime
   integer(int_kind) :: ierr, tmp
   integer(int_kind) :: vr
   integer(int_kind) :: vertnum
   integer(int_kind) :: status, ncid, varid
   logical :: kimswsamples
   type(kio_t) :: io_t
   real(real_kind), dimension(:,:), allocatable :: ups_1d, eps_1d
   real(real_kind), dimension(:,:), allocatable :: ups_2d, eps_2d
!=============================================================================
! Read Samples
!=============================================================================
!    ALLOCATE( UPs_1D(np*np*nelemd) )
!    ALLOCATE( UPs_2D(np*np*nelemd,nlev) )
!
   allocate(ups_1d(nups_l,1))
   allocate(ups_2d(nups_l,nlev))
!up
!ALLOCATE( EPs_1D(nEPs_l,1) )
!ALLOCATE( EPs_2D(nEPs_l,nlev) )
!    CALL IniInputAPIs
!.. Check for existance of model variable's variance
   status = nf90_open(trim(infilelist(1)),nf90_nowrite,ncid)
   status = nf90_inq_varid(ncid,'pint',varid)
   if (status==0) then
     if (par%ismasterproc) write(iulog,*) trim(becfile),"samples from kim-sw"
     kimswsamples = .true.
   else
     kimswsamples = .false.
   endif
   status = nf90_close(ncid)
!    ALLOCATE( x3d(np*np*nelemd*nlev) )
!    ALLOCATE( x3di(np*np*nelemd*(nlev+1)) )
!    ALLOCATE( x2d(np*np*nelemd) )
   bg_smpl = 0.d0
!up
!ALLOCATE( ps_tmp(np,np,nelemd) )
!ALLOCATE( pmid_ens(np,np,nlev,nelemd) )
!ALLOCATE( bg_smpl_tmp(np,np,nlev,nelemd) )
   allocate(ps_tmp(nups_l))
   allocate(pmid_ens(nups_l,nlev))
   allocate(bg_smpl_tmp(nups_l,nlev))
   allocate(tq_smpl_tmp(nups_l,2,ntime,nsmpl))
   allocate(zsfc_ens(nups_l))
   allocate(zsfc_ep(np,np,1,nelemd))
   allocate(zsfc_up(nups_l))
   zsfc_ep(:,:,1,:) = zsfc
   call convertep2up_hyb(1,zsfc_ep,zsfc_up)
!    ALLOCATE(pmid_up(nUPs_l,nlev))
   allocate(pmid_up(nups_l,nlev,ntime))
   allocate(bg_md_up(nups_l,4*nlev+1))
   allocate(pmid_ens_mean(nups_l,nlev,ntime))
   allocate(t_ens_mean(nups_l,ntime))
   allocate(q_ens_mean(nups_l,ntime))
   pmid_ens_mean = 0.d0
   t_ens_mean = 0.d0
   q_ens_mean = 0.d0
   do itime = 1,ntime
!       CALL ConvertEP2UP_hyb(nlev,pmid(:,:,:,:,itime),pmid_up)
     call convertep2up_hyb(nlev,pmid(:,:,:,:,itime),pmid_up(:,:,itime))
     call convertep2up_hyb(4*nlev+1,bg_md(:,:,:,:,itime),bg_md_up)
     do ic = 1,nsmpl
!..    Open a pio input file
       call iniiosystem(io_t,par%comm,par%nprocs,par%iproc)
       call openfile(io_t,trim(infilelist(ic+nsmpl*(itime-1))))
       call getdimension(io_t,'ncol',tmp)
!IF (tmp .LE. 0 .OR. tmp .EQ. nUPs_g) THEN
       if (tmp.le.0) then
         if (par%ismasterproc) print*,'Check dimension(ncol) or ep file:',trim(infilelist(ic+nsmpl*(itime-1)))
         call abortpar()
       endif
!CALL GetDimension(io_t,'ncol',tmp,DOFs_EP)
       call getdimension(io_t,'ncol',tmp,dofs)
       if (par%ismasterproc) print*,'Read samples from ',trim(infilelist(ic+nsmpl*(itime-1)))
!CALL Read2DDoubleInput('p',nlev,l_2dDoubleValues,1,ierr)
       call getvariable(io_t,'p',ups_2d)
!up
!CALL ConvertUP2EP_hyb_l(nlev,UPs_2D,EPs_2D)
       if (par%ismasterproc) write(iulog,*) 'Read pressure at middle level'
       icnt = 0
!DO ie= 1,nelemd
!DO j = 1,np
!DO i = 1,np
       do i = 1,nups_l
         icnt = icnt+1
         do k = 1,nlev
!up
!pmid_ens(i,j,nlev+1-k,ie) = EPs_2D(icnt,k)
           pmid_ens(i,nlev+1-k) = ups_2d(icnt,k)
!pmid(i,j,nlev+1-k,ie,itime) = x3d(icnt)
         enddo
!END DO
!END DO
       enddo
       pmid_ens_mean(:,:,itime) = pmid_ens_mean(:,:,itime)+pmid_ens/dble(nsmpl)
       do vr = 1,nvar3d
         select case(vr)
         case(1)
!..    Read U wind (m/s) and calculate the mean
           call getvariable(io_t,'u',ups_2d)
         case(2)
!..    Read V wind (m/s) and calculate the mean
           call getvariable(io_t,'v',ups_2d)
         case(3)
!..    Read Temperature (K) and calculate the mean
           call getvariable(io_t,'T',ups_2d)
         case(4)
!..    Read Specific humidity (kg/kg) and calculate the mean
           call getvariable(io_t,'q',ups_2d)
           if (is_moist_mr) then
             do k = 1, nlev
               do i = 1, nups_l
                 ups_2d(i,k) = ups_2d(i,k)/(1.d0+ups_2d(i,k)) ! to specific humidity
               enddo
             enddo
           endif
         end select
!up
!CALL ConvertUP2EP_hyb_l(nlev, UPs_2D, EPs_2D)
         if (bottom2top) then
           icnt = 0
!up
!DO ie = 1,nelemd
!DO j  = 1,np
!DO i  = 1,np
           do i = 1,nups_l
             icnt = icnt+1
             do k = 1,nlev
!bg_smpl(i,j,nlev+1-k+nlev*(vr-1),ie,itime,ic) = &
!                          EPs_2D(icnt,k)
!up
!bg_smpl_tmp(i,j,nlev+1-k,ie) = &
!                          EPs_2D(icnt,k)
               bg_smpl_tmp(i,nlev+1-k) = ups_2d(icnt,k)
!bg_smpl(i,j,nlev+1-k+nlev*(vr-1),ie,itime,ic) = x3d(icnt)
             enddo
!END DO
!END DO
           enddo
         else
           icnt = 0
!up
!DO ie = 1,nelemd
!DO j  = 1,np
!DO i  = 1,np
           do i = 1,nups_l
             icnt = icnt+1
             do k = 1,nlev
!bg_smpl(i,j,k+nlev*(vr-1),ie,itime,ic) = EPs_2D(icnt,k)
!up
!bg_smpl_tmp(i,j,k,ie) = EPs_2D(icnt,k)
               bg_smpl_tmp(i,k) = ups_2d(icnt,k)
!bg_smpl(i,j,k+nlev*(vr-1),ie,itime,ic) = x3d(icnt)
             enddo
!END DO
!END DO
           enddo
         endif
!--------------------
! Adjust pressure of LETKF samples to H-4DEV background
!up
!DO ie = 1, nelemd
!DO j  = 1, np
!DO i  = 1, np
         do i = 1, nups_l
           do k = 1, nlev
!CALL VertInt(p_bgh_v, p_ob, nlev, obmiss, &
!             bgh,                         &
!             bgv)
!up
!CALL VertInt(pmid_ens(i,j,:,ie),                &
!             pmid(i,j,k,ie,itime),              &
             call vertint(pmid_ens(i,:),! pmid_up(i,k),pmid_up(i,k,itime),nlev,!bg_md(i,j,k+nlev*(vr-1),ie,itime),!up
!bg_smpl_tmp(i,j,k,ie),            &
!bg_smpl_tmp(i,j,:,ie),            &
!bg_smpl(i,j,k+nlev*(vr-1),ie,itime,ic))
             bg_smpl_tmp(i,k),bg_smpl_tmp(i,:),bg_smpl(i,k+nlev*(vr-1),itime,ic))
           enddo
         enddo
!END DO
!END DO
         if (vr.eq.3) then
           tq_smpl_tmp(:,1,itime,ic) = bg_smpl_tmp(:,nlev)
         elseif (vr.eq.4) then
           tq_smpl_tmp(:,2,itime,ic) = bg_smpl_tmp(:,nlev)
         endif
!          IF (itime.EQ.1 .AND. ic.EQ.1 .AND. par%isMasterProc)  THEN
!             PRINT *, 'vr', 'ie', 'k', 'j', 'i', &
!                      'pmid_ens', 'pmid',        &
!                      'bg_smpl_tmp', 'bg_smpl'
!             DO ie = 1, nelemd
!             DO j  = 1, np
!             DO i  = 1, np
!             DO k  = 1, nlev
!                   PRINT *, vr, ie, j, i, k,                          &
!                            pmid_ens(i,j,k,ie), pmid(i,j,k,ie,itime), &
!                            bg_smpl_tmp(i,j,k,ie), bg_smpl(i,j,k+nlev*(vr-1),ie,itime,ic)
!             END DO
!             END DO
!             END DO
!             END DO
!          END IF
!-----------------------------------------------------------------
       enddo ! vr
       if (par%ismasterproc) print*,'Read u,v,t,q'
       do vr = 1,nvar2d
         select case(vr)
         case(1)
!..    Read surface pressure (hPa -> Pa) and calculate the mean
!CALL GetVariable(io_t, 'ps', UPs_1D)
           call getvariable(io_t,'ps',ups_1d(:,1))
!up
!CALL ConvertUP2EP_hyb_l(1,UPs_1D,EPs_1D)
           psfact = 1.d0
           if (par%ismasterproc) print*,'Read ps,psfact = ',psfact
         end select
         icnt = 0
!up
!DO ie = 1,nelemd
!DO j  = 1,np
!DO i  = 1,np
         do i = 1,nups_l
           icnt = icnt+1
!bg_smpl(i,j,nlev*nvar3d+vr,ie,itime,ic) = UPs_1D(icnt)*psfact
!ps_tmp(i,j,ie) = UPs_1D(icnt)*psfact
!up
!ps_tmp(i,j,ie) = EPs_1D(icnt,1)*psfact
           ps_tmp(i) = ups_1d(icnt,1)*psfact
!bg_smpl(i,j,nlev*nvar3d+vr,ie,itime,ic) = x2d(icnt)*psfact
!END DO
!END DO
         enddo
       enddo
!------------------
! For ps adjustment according to topography difference
!-----------------
! Read topography of an ensemble sample
       if (itime.eq.1.and.ic.eq.1) then
         if (par%ismasterproc) print*,'Read topo'
!CALL GetVariable(io_t,'topo',UPs_1D)
         call getvariable(io_t,'topo',ups_1d(:,1))
!up
!CALL ConvertUP2EP_hyb_l(1,UPs_1D,EPs_1D)
         icnt = 0
!DO ie = 1,nelemd
!DO j  = 1,np
!DO i  = 1,np
         do i = 1,nups_l
           icnt = icnt+1
!zsfc_ens(i,j,ie) = EPs_1D(icnt,1)
           zsfc_ens(i) = ups_1d(icnt,1)
!END DO
!END DO
         enddo
       endif
! Adjusting ps
!up
!DO ie = 1,nelemd
!DO j  = 1,np
!DO i  = 1,np
       do i = 1,nups_l
!CALL SfcPrsCorr (bg_md(i,j,3*nlev,ie,itime),      &
!                 bg_md(i,j,4*nlev,ie,itime),      &
!                 zsfc(i,j,ie),                    &
!                 bg_smpl(i,j,3*nlev,ie,itime,ic), &
!                 bg_smpl(i,j,4*nlev,ie,itime,ic), &
!                 zsfc_ens(i,j,ie),                &
!                 ps_tmp(i,j,ie),                  &
!                 bg_smpl(i,j,4*nlev+1,ie,itime,ic))
         call sfcprscorr(bg_md_up(i,3*nlev),bg_md_up(i,4*nlev),zsfc_up(i),!bg_smpl(i,3*nlev,itime,ic),!bg_smpl(i,4*nlev,itime,ic),tq_smpl_tmp(i,1,itime,ic),tq_smpl_tmp(i,2,itime,ic),zsfc_ens(i),ps_tmp(i),bg_smpl(i,4*nlev+1,itime,ic))
!END DO
!END DO
       enddo
       t_ens_mean(:,itime) = t_ens_mean(:,itime)+tq_smpl_tmp(:,1,itime,ic)/dble(nsmpl)
       q_ens_mean(:,itime) = q_ens_mean(:,itime)+tq_smpl_tmp(:,2,itime,ic)/dble(nsmpl)
! Checking adjusted ps
       if (itime.eq.1.and.ic.eq.1) then
         if (par%ismasterproc) print*,'zsfc','zsfc_ens','bg_md_ps','ps_ens','bg_smpl_ps'
!up
!DO ie = 1,nelemd
!DO j  = 1,np
!DO i  = 1,np
         do i = 1,nups_l
           if (dabs(zsfc_up(i)-zsfc_ens(i)).gt.500.d0) then
             print*,zsfc_up(i),zsfc_ens(i),bg_md_up(i,4*nlev+1),ps_tmp(i),bg_smpl(i,4*nlev+1,itime,ic)
           endif
!END DO
!END DO
         enddo
       endif
!--------------------------------------------------------------
!..    Close the pio file
       call closefile(io_t)
       call finiosystem(io_t)
     enddo !! ic
   enddo !! itime
   deallocate(ups_1d)
   deallocate(ups_2d)
   deallocate(zsfc_ep)
!    DEALLOCATE( zsfc_up )
!    DEALLOCATE( zsfc_ens )
!    DEALLOCATE( pmid_up )
   deallocate(bg_smpl_tmp)
   deallocate(tq_smpl_tmp)
   deallocate(bg_md_up)
!DEALLOCATE( EPs_1D )
!DEALLOCATE( EPs_2D )
!CALL FinInputAPIs
!
   end subroutine readsamples
!===============================================================================
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine adjustp(an_ep)
!===============================================================================
!-------------------------------------------------------------------------------
   use vertintmod, only: vertint
   use obscorrmod, only: sfcprscorr
!
   real(real_kind), dimension(np, np, nlevar, nelemd, ntime), intent(inout) :: an_ep
! Local variables
   real(real_kind), dimension(:,:), allocatable :: an_up, an_up_adj
   integer(int_kind) :: i, k, vr, itime
!
   allocate(an_up(nups_l,4*nlev+1))
   allocate(an_up_adj(nups_l,4*nlev+1))
   do itime = 1,ntime
     call convertep2up_hyb(4*nlev+1,an_ep(:,:,:,:,itime),an_up)
     do vr = 1,nvar3d
       do i = 1, nups_l
         do k = 1, nlev
           call vertint(pmid_up(i,:,itime),pmid_ens_mean(i,k,itime),nlev,an_up(i,k+nlev*(vr-1)),an_up(i,1+nlev*(vr-1):nlev+nlev*(vr-1)),an_up_adj(i,k+nlev*(vr-1)))
         enddo
       enddo
     enddo
     do i = 1,nups_l
       call sfcprscorr(t_ens_mean(i,itime),q_ens_mean(i,itime),zsfc_ens(i),an_up(i,3*nlev),an_up(i,4*nlev),zsfc_up(i),an_up(i,4*nlev+1),an_up_adj(i,4*nlev+1))
     enddo
     call convertup2ep_hyb(4*nlev+1,an_up_adj,an_ep(:,:,:,:,itime))
   enddo
   deallocate(pmid_up,pmid_ens_mean)
   deallocate(an_up,an_up_adj)
   deallocate(t_ens_mean,q_ens_mean)
   deallocate(zsfc_ens,zsfc_up)
   return
!
   end subroutine adjustp
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine readback_kim(itime)
!===============================================================================
! Local variables
!
   integer(int_kind), intent(in   ) :: itime
   real(real_kind), dimension(:), allocatable :: x3d, x3di, x2d
   type(kio_t) :: io_t
   real(real_kind), dimension(:,:), allocatable :: ups_1d, eps_1d
   real(real_kind), dimension(:,:), allocatable :: ups_2d, ups_2d_p1, eps_2d, eps_2d_p1
   real(real_kind), dimension(:,:,:,:), allocatable :: bg_ct
   real(real_kind) :: psfact
   integer(int_kind) :: i, j, k, ie, icnt
   integer(int_kind) :: ierr, tmp
   integer(int_kind) :: vr
   integer(int_kind) :: vertnum
!=============================================================================
! Prepare ReadPio
!=============================================================================
!=============================================================================
! Read Background
!=============================================================================
!ALLOCATE( UPs_1D(np*np*nelemd) )
!ALLOCATE( UPs_2D(np*np*nelemd,nlev) )
!ALLOCATE( UPs_2D_p1(np*np*nelemd,nlev+1) )
!
   allocate(ups_1d(nups_l,1))
   allocate(ups_2d(nups_l,nlev))
   allocate(ups_2d_p1(nups_l,nlev+1))
   allocate(eps_1d(neps_l,1))
   allocate(eps_2d(neps_l,nlev))
   allocate(eps_2d_p1(neps_l,nlev+1))
   call iniiosystem(io_t,par%comm,par%nprocs,par%iproc)
   call openfile(io_t,trim(bgfile(itime)))
   call getdimension(io_t,'ncol',tmp)
!IF (tmp .LE. 0 .OR. tmp .EQ. nUPs_g) THEN
   if (tmp.le.0) then
     if (par%ismasterproc) print*,'Check dimension or ep file:',trim(bgfile(itime))
     call abortpar()
   endif
   call getdimension(io_t,'ncol',tmp,dofs)
   if (par%ismasterproc) write(iulog,*) 'nlev = ',nlev
   if (par%ismasterproc) write(iulog,*) 'Bacground from kim model'
   if (par%ismasterproc) write(iulog,*) 'File of background = ',bgfile(itime)
   do vr = 1,nvar3d
     select case(vr)
     case(1)
       call getvariable(io_t,'u',ups_2d)
       if (par%ismasterproc) write(iulog,*) 'Read u'
     case(2)
       call getvariable(io_t,'v',ups_2d)
       if (par%ismasterproc) write(iulog,*) 'Read v'
     case(3)
       call getvariable(io_t,'T',ups_2d)
       if (par%ismasterproc) write(iulog,*) 'Read t'
     case(4)
       call getvariable(io_t,'q',ups_2d)
       if (par%ismasterproc) write(iulog,*) 'Read q'
       if (is_moist_mr) then
         do k = 1, nlev
           do i = 1, nups_l
             ups_2d(i,k) = ups_2d(i,k)/(1.d0+ups_2d(i,k)) ! to specific humidity
           enddo
         enddo
       endif
     end select
     call convertup2ep_hyb_l(nlev,ups_2d,eps_2d)
     do k = 1,nlev
       icnt = 0
       do ie = 1,nelemd
         do j = 1, np
           do i = 1, np
             icnt = icnt+1
             bg_md(i,j,nlev+1-k+nlev*(vr-1),ie,itime) = eps_2d(icnt,k)
!bg_md(i,j,nlev+1-k+nlev*(vr-1),ie,itime) = UPs_2D(icnt,k)
!bg_md(i,j,nlev+1-k+nlev*(vr-1),ie,itime) = x3d(icnt,k)
           enddo
         enddo
       enddo
     enddo
!       IF (par%ismasterproc) WRITE(iulog,*) 'vr    = ',vr
!       IF (par%ismasterproc) WRITE(iulog,*) 'input = ',UPs_2D(1,:)
!       IF (par%ismasterproc) WRITE(iulog,*) 'save  = ',bg_md(1,1,nlev*(vr-1)+1:nlev*(vr-1)+nlev,1,itime)
   enddo
   bg_md_q(:,:,:,:) = bg_md(:,:,nlev*3+1:nlev*4,:,itime)
   if (cv_opt_hum==1) then
     call getvariable(io_t,'rh',ups_2d)
     call convertup2ep_hyb_l(nlev,ups_2d,eps_2d)
     if (par%ismasterproc) write(iulog,*) 'Read rh'
     vr = 4
     do k = 1,nlev
       icnt = 0
       do ie = 1,nelemd
         do j = 1, np
           do i = 1, np
             icnt = icnt+1
             bg_md(i,j,nlev+1-k+nlev*(vr-1),ie,itime) = eps_2d(icnt,k)
!bg_md(i,j,nlev+1-k+nlev*(vr-1),ie,itime) = x3d(icnt)
           enddo
         enddo
       enddo
     enddo
   endif
   do vr = 1,nvar2d
     select case(vr)
     case(1)
       call getvariable(io_t,'ps',ups_1d(:,1))
       call convertup2ep_hyb_l(1,ups_1d,eps_1d)
       if (par%ismasterproc) write(iulog,*) 'Read ps'
     end select
     psfact = 1.d0
     if (eps_1d(1,1)<2000.) psfact = 1.d+2 ! if unit of ps is hpa,convert into pa
!IF( x2d(1) < 2000. ) psfact=1.D+2  ! If unit of PS is hPa,convert into Pa
     icnt = 0
     do ie = 1,nelemd
       do j = 1, np
         do i = 1, np
           icnt = icnt+1
           bg_md(i,j,nlevar,ie,itime) = eps_1d(icnt,1)*psfact
!bg_md(i,j,nlevar,ie,itime) = x2d(icnt)*psfact
         enddo
       enddo
     enddo
   enddo
!! u10m
   call getvariable(io_t,'u10m',ups_1d(:,1))
   call convertup2ep_hyb_l(1,ups_1d,eps_1d)
   if (par%ismasterproc) write(iulog,*) 'Read u10m'
   icnt = 0
   do ie = 1,nelemd
     do j = 1, np
       do i = 1, np
         icnt = icnt+1
         u10m(i,j,ie,itime) = eps_1d(icnt,1)
       enddo
     enddo
   enddo
!! v10m
   call getvariable(io_t,'v10m',ups_1d(:,1))
   call convertup2ep_hyb_l(1,ups_1d,eps_1d)
   if (par%ismasterproc) write(iulog,*) 'Read v10m'
   icnt = 0
   do ie = 1,nelemd
     do j = 1, np
       do i = 1, np
         icnt = icnt+1
         v10m(i,j,ie,itime) = eps_1d(icnt,1)
       enddo
     enddo
   enddo
!! t2m
   call getvariable(io_t,'t2m',ups_1d(:,1))
   call convertup2ep_hyb_l(1,ups_1d,eps_1d)
   if (par%ismasterproc) write(iulog,*) 'Read t2m'
   icnt = 0
   do ie = 1,nelemd
     do j = 1, np
       do i = 1, np
         icnt = icnt+1
         t2m(i,j,ie,itime) = eps_1d(icnt,1)
       enddo
     enddo
   enddo
!! q2m
   call getvariable(io_t,'q2m',ups_1d(:,1))
   if (is_moist_mr) then
     do i = 1, nups_l
       ups_1d(i,1) = ups_1d(i,1)/(1.d0+ups_1d(i,1)) ! to specific humidity
     enddo
   endif
   call convertup2ep_hyb_l(1,ups_1d,eps_1d)
   if (par%ismasterproc) write(iulog,*) 'Read q2m'
   icnt = 0
   do ie = 1,nelemd
     do j = 1, np
       do i = 1, np
         icnt = icnt+1
         q2m(i,j,ie,itime) = eps_1d(icnt,1)
       enddo
     enddo
   enddo
!! tsfc
   call getvariable(io_t,'tsfc',ups_1d(:,1))
   call convertup2ep_hyb_l(1,ups_1d,eps_1d)
   if (par%ismasterproc) write(iulog,*) 'Read tsfc'
   icnt = 0
   do ie = 1,nelemd
     do j = 1, np
       do i = 1, np
         icnt = icnt+1
         tsfc(i,j,ie,itime) = eps_1d(icnt,1)
       enddo
     enddo
   enddo
!! sfctype
   call getvariable(io_t,'slmsk',ups_1d(:,1))
   call convertup2ep_hyb_l(1,ups_1d,eps_1d)
   if (par%ismasterproc) write(iulog,*) 'Read sfctype'
   icnt = 0
   do ie = 1,nelemd
     do j = 1, np
       do i = 1, np
         icnt = icnt+1
         sfctype(i,j,ie,itime) = eps_1d(icnt,1)
       enddo
     enddo
   enddo
!! znt
   call getvariable(io_t,'znt',ups_1d(:,1))
   call convertup2ep_hyb_l(1,ups_1d,eps_1d)
   if (par%ismasterproc) write(iulog,*) 'Read znt'
   icnt = 0
   do ie = 1,nelemd
     do j = 1, np
       do i = 1, np
         icnt = icnt+1
         rl_org(i,j,ie) = eps_1d(icnt,1)
!rl_org(i,j,ie) = x2d(icnt)
       enddo
     enddo
   enddo
!! topo
   call getvariable(io_t,'topo',ups_1d(:,1))
   call convertup2ep_hyb_l(1,ups_1d,eps_1d)
   if (par%ismasterproc) write(iulog,*) 'Read topo'
   icnt = 0
   do ie = 1,nelemd
     do j = 1, np
       do i = 1, np
         icnt = icnt+1
         zsfc(i,j,ie) = eps_1d(icnt,1)
!zsfc(i,j,ie) = x2d(icnt)
       enddo
     enddo
   enddo
   if (model_name=='kim-sw') then ! read p and pint for kim-sw
     call getvariable(io_t,'p',ups_2d)
     call convertup2ep_hyb_l(nlev,ups_2d,eps_2d)
     if (par%ismasterproc) write(iulog,*) 'Read pressure at middle level'
     do k = 1,nlev
       icnt = 0
       do ie = 1,nelemd
         do j = 1, np
           do i = 1, np
             icnt = icnt+1
             pmid(i,j,nlev+1-k,ie,itime) = eps_2d(icnt,k)
!pmid(i,j,nlev+1-k,ie,itime) = x3d(icnt)
           enddo
         enddo
       enddo
     enddo
     call getvariable(io_t,'pint',ups_2d_p1)
     call convertup2ep_hyb_l(nlev+1,ups_2d_p1,eps_2d_p1)
     if (par%ismasterproc) write(iulog,*) 'Read pressure at interface level'
     do k = 1,nlev+1
       icnt = 0
       do ie = 1,nelemd
         do j = 1, np
           do i = 1, np
             icnt = icnt+1
             pint(i,j,nlev+2-k,ie) = eps_2d_p1(icnt,k)
!pint(i,j,nlev+2-k,ie) = x3di(icnt)
           enddo
         enddo
       enddo
     enddo
   else ! calculate p and pint for kim-sh
     do k = 1, nlev
       do ie = 1, nelemd
         do j = 1, np
           do i = 1, np
             pmid(i,j,k,ie,itime) = hyam(k)*ps0+hybm(k)*bg_md(i,j,nlevar,ie,itime)
             pint(i,j,k,ie) = hyai(k)*ps0+hybi(k)*bg_md(i,j,nlevar,ie,itime)
           enddo
         enddo
       enddo
     enddo
     k = nlev+1
     do ie = 1,nelemd
       do j = 1, np
         do i = 1, np
           pint(i,j,k,ie) = hyai(k)*ps0+hybi(k)*bg_md(i,j,nlevar,ie,itime)
         enddo
       enddo
     enddo
   endif ! model_name
   if (model_name=='kim-sw') then ! read hgt for kim-sw
     call getvariable(io_t,'hgt',ups_2d)
     call convertup2ep_hyb_l(nlev,ups_2d,eps_2d)
     if (par%ismasterproc) write(iulog,*) 'Read hgt'
     do k = 1,nlev
       icnt = 0
       do ie = 1,nelemd
         do j = 1, np
           do i = 1, np
             icnt = icnt+1
             z_org(i,j,nlev+1-k,ie) = eps_2d(icnt,k)
!z_org(i,j,nlev+1-k,ie) = x3d(icnt)
           enddo
         enddo
       enddo
     enddo
   else ! calculate hgt for kim-sh
     if (par%ismasterproc) write(iulog,*) 'Calculate hgt'
   endif ! model_name
   if (cv_opt_wind.eq.1) then
     allocate(bg_ct(np,np,nlev,nelemd))
     bg_uv_rot(:,:,:,1,:,itime) = bg_md(:,:,1:nlev,:,itime)
     bg_uv_rot(:,:,:,2,:,itime) = bg_md(:,:,nlev+1:2*nlev,:,itime)
     call calpsi(bg_uv_rot(:,:,:,:,:,itime),bg_ct)
     call psitorotwind(bg_ct,bg_uv_rot(:,:,:,:,:,itime))
     deallocate(bg_ct)
   endif
   phis = zsfc*gravi
   if (par%ismasterproc) write(iulog,*) 'pmid(1,1,:,1,itime):',pmid(1,1,:,1,itime)
   if (par%ismasterproc) write(iulog,*) 'zsfc(1,1,1):',zsfc(1,1,1)
   if (par%ismasterproc) write(iulog,*) 'z_org(1,1,:,1):',z_org(1,1,:,1)
   deallocate(ups_1d)
   deallocate(ups_2d)
   deallocate(ups_2d_p1)
   deallocate(eps_1d)
   deallocate(eps_2d)
   deallocate(eps_2d_p1)
   call closefile(io_t)
   call finiosystem(io_t)
!
   end subroutine readback_kim
!===============================================================================
!===============================================================================
! Modified by HJ Song for less-weighted background error covariance I/O
! On 29JAN2016
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine readbackerrcov(l_reg_avg, l_eigf_psi, l_eigf_chi, l_eigf_tps, l_eigf_hum, errvar_md)
!-------------------------------------------------------------------------------
   use netcdf
!
   type(kio_t) :: io_t
   real(real_kind), dimension(reg_nproc), intent(  out) :: l_reg_avg
   real(real_kind), dimension(psi_nproc), intent(  out) :: l_eigf_psi
   real(real_kind), dimension(chi_nproc), intent(  out) :: l_eigf_chi
   real(real_kind), dimension(tps_nproc), intent(  out) :: l_eigf_tps
   real(real_kind), dimension(hum_nproc), intent(  out) :: l_eigf_hum
   real(real_kind), dimension(np, np, nlevar, nelemd), intent(  out), optional :: errvar_md
! Local variables
   real(real_kind), dimension(:), allocatable :: w3d, w3d_chi, w3d_tps, w3d_q, x2d, c3d, reg_tmp
   real(real_kind), dimension(:,:), allocatable :: t3d, w2d, w2d_chi, w2d_tps, w2d_q, x3d
   real(real_kind) :: psfact
   real(real_kind) :: tmpmin, tmpmax
   real(real_kind) :: minv, maxv
   integer(int_kind) :: i, j, k, l
   integer(int_kind) :: ie, icnt, ig, wn, swn
   integer(int_kind) :: ierr, tmp
   integer(int_kind) :: vr
   integer(int_kind) :: vertnum
   integer(int_kind) :: status, ncid, varid
   logical :: exst_varmd
!=============================================================================
! Read Backgroun Error Covariance
!=============================================================================
!.. Check for existance of model variable's variance
!
   status = nf90_open(trim(becfile),nf90_nowrite,ncid)
   status = nf90_inq_varid(ncid,'varU',varid)
   if (status/= 0) then
     if (par%ismasterproc) write(iulog,*) trim(becfile),"no variance of model variables in backerrcov"
     exst_varmd = .false.
   else
     exst_varmd = .true.
   endif
   status = nf90_close(ncid)
   allocate(t3d(np*np*nelemd,nlevar))
!ALLOCATE( w2d(zwn2*neig) )
   allocate(w2d(zwn_nproc,neig))
!ALLOCATE( w3d(zwn2*neig*nlev) )
   allocate(w3d(psi_nproc))
!ALLOCATE( w2d_chi(zwn2*neig_chi) )
   allocate(w2d_chi(zwn_nproc,neig_chi))
!ALLOCATE( w3d_chi(zwn2*neig_chi*nlev) )
   allocate(w3d_chi(chi_nproc))
!ALLOCATE( w2d_Tps(zwn2*neig_Tps) )
   allocate(w2d_tps(zwn_nproc,neig_tps))
!ALLOCATE( w3d_Tps(zwn2*neig_Tps*(nlev+1)) )
   allocate(w3d_tps(tps_nproc))
!ALLOCATE( w2d_q(zwn2*neig_q) )
   allocate(w2d_q(zwn_nproc,neig_q))
!ALLOCATE( w3d_q(zwn2*neig_q*nlev) )
   allocate(w3d_q(hum_nproc))
!.. Open PIO for Backgroun Error Covariance file
   call iniiosystem(io_t,par%comm,par%nprocs,par%iproc)
   call openfile(io_t,trim(becfile))
   call getdimension(io_t,'ncol',tmp)
   if (tmp.le.0.or.tmp.eq.nups_g) then
     if (par%ismasterproc) print*,'Check dimension or ep file:',trim(becfile)
     call abortpar()
   endif
   call getdimension(io_t,'ncol',tmp,dofs_ep)
   call getdimension(io_t,'nwav',tmp,map_zwn)
   call genmap(reg_nproc,reg_iproc,map_reg)
   call getdimension(io_t,'nreg',tmp,map_reg)
   call getdimension(io_t,'npsi',tmp,map_psi)
   call getdimension(io_t,'nchi',tmp,map_chi)
   call getdimension(io_t,'nTps',tmp,map_tps)
   call getdimension(io_t,'nhum',tmp,map_hum)
!PRINT *,'map zwn 1 = ',map_zwn(1)
!PRINT *,'map reg 1 = ',map_reg(1)
   if (par%ismasterproc) write(iulog,*) 'File of background error covariance:',becfile
   call getvariable(io_t,'vari',t3d)
   if (par%ismasterproc) write(iulog,*) 'Read variance of control variables'
   do k = 1,nlevar
     icnt = 0
     do ie = 1,nelemd
       do j = 1, np
         do i = 1, np
           icnt = icnt+1
           errvar(i,j,k,ie) = t3d(icnt,k)
         enddo
       enddo
     enddo
   enddo
   tmpmax = maxval(errvar)
   tmpmin = minval(errvar)
   call par_reduce(tmpmax,maxv,1,par_max)
   call par_reduce(tmpmin,minv,1,par_min)
! maxv = GlobalMAX(errvar) ! USE ParFunctions
! minv = GlobalMIN(errvar) !
   if (par%ismasterproc) write(iulog,*) 'Min !.. read eigenmodes
!.. Read eigenvalues
!.. Read eigenvalues of psi
   call getvariable(io_t,'eigv',w2d)
   if (par%ismasterproc) write(iulog,*) 'Read eigenvalue'
!!!ALLOCATE( w2d(zwn_nproc,neig) )
   do ig = 1,neig
     do wn = 1, zwn_nproc
       eigv(ig,wn) = w2d(wn,ig)
     enddo
   enddo
!PRINT *,SIZE(eigv,DIM=1),SIZE(eigv,DIM=2),SIZE(w2d,DIM=1),SIZE(w2d,DIM=2)
!eigv(:,:) = w2d(:,:)
   tmpmax = maxval(eigv)
   tmpmin = minval(eigv)
   call par_reduce(tmpmax,maxv,1,par_max)
   call par_reduce(tmpmin,minv,1,par_min)
   if (par%ismasterproc) write(iulog,*) 'Min !.. read eigenvalues of chi
   call getvariable(io_t,'egvc',w2d_chi)
   if (par%ismasterproc) write(iulog,*) 'Read eigenvalue_chi'
!icnt = 0
   do ig = 1,neig_chi
     do wn = 1, zwn_nproc
       eigv_chi(ig,wn) = w2d_chi(wn,ig)
     enddo
   enddo
   tmpmax = maxval(eigv_chi)
   tmpmin = minval(eigv_chi)
   call par_reduce(tmpmax,maxv,1,par_max)
   call par_reduce(tmpmin,minv,1,par_min)
   if (par%ismasterproc) write(iulog,*) 'Min !.. read eigenvalues of(t,ps)
   call getvariable(io_t,'egvT',w2d_tps)
   if (par%ismasterproc) write(iulog,*) 'Read eigenvalue_tps'
   do ig = 1,neig_tps
     do wn = 1, zwn_nproc
       eigv_tps(ig,wn) = w2d_tps(wn,ig)
     enddo
   enddo
   tmpmax = maxval(eigv_tps)
   tmpmin = minval(eigv_tps)
   call par_reduce(tmpmax,maxv,1,par_max)
   call par_reduce(tmpmin,minv,1,par_min)
   if (par%ismasterproc) write(iulog,*) 'Min !.. read eigenvalues of q
   call getvariable(io_t,'egvq',w2d_q)
   if (par%ismasterproc) write(iulog,*) 'Read eigenvalue_q'
   do ig = 1,neig_q
     do wn = 1, zwn_nproc
       eigv_q(ig,wn) = w2d_q(wn,ig)
     enddo
   enddo
   tmpmax = maxval(eigv_q)
   tmpmin = minval(eigv_q)
   call par_reduce(tmpmax,maxv,1,par_max)
   call par_reduce(tmpmin,minv,1,par_min)
   if (par%ismasterproc) write(iulog,*) 'Min ! read eigenvectors
! Read eigenvectors of psi
   call getvariable(io_t,'eigf',w3d)
   l_eigf_psi = w3d
   if (par%ismasterproc) write(iulog,*) 'Read eigenfunction'
!    icnt = 0
!    DO k =1,nlev
!       DO ig=1,neig
!       DO wn=1,zwn_nproc
!          icnt= icnt + 1
!          bcov_sp(k,ig,wn)= w3d(icnt)
!       END DO
!       END DO
!    END DO
!.. Read eigenvectors of chi
   call getvariable(io_t,'egfc',w3d_chi)
   l_eigf_chi = w3d_chi
   if (par%ismasterproc) write(iulog,*) 'Read eigenfunction_chi'
!    icnt = 0
!    DO k =1,nlev
!       DO ig=1,neig_chi
!       DO wn=1,zwn_nproc
!          icnt= icnt + 1
!          bcov_sp_chi(k,ig,wn)= w3d_chi(icnt)
!       END DO
!       END DO
!    END DO
!.. Read eigenvectors of Tps
   call getvariable(io_t,'egfT',w3d_tps)
   l_eigf_tps = w3d_tps
   if (par%ismasterproc) write(iulog,*) 'Read eigenfunction_tps'
!    icnt = 0
!    DO k =1,nlev+1
!       DO ig=1,neig_Tps
!       DO wn=1,zwn_nproc
!          icnt= icnt + 1
!          bcov_sp_Tps(k,ig,wn)= w3d_Tps(icnt)
!       END DO
!       END DO
!    END DO
!.. Read eigenvectors of q
   call getvariable(io_t,'egfq',w3d_q)
   l_eigf_hum = w3d_q
   if (par%ismasterproc) write(iulog,*) 'Read eigenfunction_q'
!    icnt = 0
!    DO k =1,nlev
!       DO ig=1,neig_q
!       DO wn=1,zwn_nproc
!          icnt= icnt + 1
!          bcov_sp_q(k,ig,wn)= w3d_q(icnt)
!       END DO
!       END DO
!    END DO
!.. Check eigenfunctions
   tmpmax = maxval(l_eigf_psi)
   tmpmin = minval(l_eigf_psi)
   call par_reduce(tmpmax,maxv,1,par_max)
   call par_reduce(tmpmin,minv,1,par_min)
   if (par%ismasterproc) write(iulog,*) 'Min tmpmax = maxval(l_eigf_chi)
   tmpmin = minval(l_eigf_chi)
   call par_reduce(tmpmax,maxv,1,par_max)
   call par_reduce(tmpmin,minv,1,par_min)
   if (par%ismasterproc) write(iulog,*) 'Min tmpmax = maxval(l_eigf_tps)
   tmpmin = minval(l_eigf_tps)
   call par_reduce(tmpmax,maxv,1,par_max)
   call par_reduce(tmpmin,minv,1,par_min)
   if (par%ismasterproc) write(iulog,*) 'Min tmpmax = maxval(l_eigf_hum)
   tmpmin = minval(l_eigf_hum)
   call par_reduce(tmpmax,maxv,1,par_max)
   call par_reduce(tmpmin,minv,1,par_min)
   if (par%ismasterproc) write(iulog,*) 'Min ! read variance
   if (exst_varmd) then
     allocate(x2d(np*np*nelemd))
     allocate(x3d(np*np*nelemd,nlev))
     do vr = 1,nvar3d
       select case(vr)
       case(1)
         call getvariable(io_t,'varU',x3d)
         if (par%ismasterproc) write(iulog,*) 'Read variance of u'
       case(2)
         call getvariable(io_t,'varV',x3d)
         if (par%ismasterproc) write(iulog,*) 'Read variance of v'
       case(3)
         call getvariable(io_t,'varT',x3d)
         if (par%ismasterproc) write(iulog,*) 'Read variance of t'
       case(4)
         call getvariable(io_t,'varQ',x3d)
         if (par%ismasterproc) write(iulog,*) 'Read variance of q'
       end select
       do k = 1, nlev
         icnt = 0
         do ie = 1,nelemd
           do j = 1, np
             do i = 1, np
               icnt = icnt+1
               errvar_md(i,j,k+nlev*(vr-1),ie) = x3d(icnt,k)
             enddo
           enddo
         enddo
       enddo
     enddo
     call getvariable(io_t,'varP',x2d)
     if (par%ismasterproc) write(iulog,*) 'Read variance of ps'
     icnt = 0
     do ie = 1,nelemd
       do j = 1, np
         do i = 1, np
           icnt = icnt+1
           errvar_md(i,j,nlevar,ie) = x2d(icnt)
         enddo
       enddo
     enddo
     deallocate(x2d,x3d)
   endif ! model_name
   if (cv_opt_tu==1.or.cv_opt_tu==2.or.cv_opt_tu==3) then
!ALLOCATE( c3d(np*np*nelem*nlev*nk) )
     allocate(c3d(reg_nproc))
     call getvariable(io_t,'regn',c3d)
     l_reg_avg = c3d
     if (par%ismasterproc) write(iulog,*) 'Read regression coefficient'
     tmpmax = maxval(l_reg_avg)
     tmpmin = minval(l_reg_avg)
     call par_reduce(tmpmax,maxv,1,par_max)
     call par_reduce(tmpmin,minv,1,par_min)
     if (par%ismasterproc) write(iulog,*) 'Min deallocate(c3d)
!DEALLOCATE( dimSizesRegCoef )
   endif ! cv_opt_tu
   call closefile(io_t)
   call finiosystem(io_t)
   deallocate(t3d,w2d,w3d)
   deallocate(w2d_chi,w3d_chi)
   deallocate(w2d_tps,w3d_tps)
   deallocate(w2d_q,w3d_q)
   return
!
   end subroutine readbackerrcov
!===============================================================================
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine readnature
!===============================================================================
! Local variables
!
   type(kio_t) :: io_t
   real(real_kind), dimension(:,:), allocatable :: x3d
   real(real_kind), dimension(:), allocatable :: x2d
   real(real_kind) :: psfact
   integer(int_kind) :: i, j, k, ie, icnt
   integer(int_kind) :: ierr, tmp
   integer(int_kind) :: vr
   integer(int_kind) :: vertnum
!=============================================================================
! Read Nature
!=============================================================================
!
   allocate(x3d(np*np*nelemd,nlev))
   allocate(x2d(np*np*nelemd))
   call iniiosystem(io_t,par%comm,par%nprocs,par%iproc)
   call openfile(io_t,trim(naturefile))
   call getdimension(io_t,'ncol',tmp)
   if (tmp.le.0.or.tmp.eq.nups_g) then
     if (par%ismasterproc) print*,'Check dimension or ep file:',trim(naturefile)
     call abortpar()
   endif
   call getdimension(io_t,'ncol',tmp,dofs_ep)
   if (par%ismasterproc) write(iulog,*) 'File of nature = ',naturefile
   do vr = 1,nvar3d
     select case(vr)
     case(1)
       call getvariable(io_t,'u',x3d)
     case(2)
       call getvariable(io_t,'v',x3d)
     case(3)
       call getvariable(io_t,'T',x3d)
     case(4)
       if (cv_opt_hum==0) then
         call getvariable(io_t,'q',x3d)
       elseif (cv_opt_hum==1) then
         call getvariable(io_t,'rh',x3d)
       endif
     end select
     if (model_name=='kim-old') then
       do k = 1, nlev
         icnt = 0
         do ie = 1,nelemd
           do j = 1, np
             do i = 1, np
               icnt = icnt+1
               nt_org(i,j,k+nlev*(vr-1),ie) = x3d(icnt,k)
             enddo
           enddo
         enddo
       enddo
     else
       do k = 1, nlev
         icnt = 0
         do ie = 1,nelemd
           do j = 1, np
             do i = 1, np
               icnt = icnt+1
               nt_org(i,j,nlev+1-k+nlev*(vr-1),ie) = x3d(icnt,k)
             enddo
           enddo
         enddo
       enddo
     endif
   enddo
   do vr = 1,nvar2d
     select case(vr)
     case(1)
       call getvariable(io_t,'ps',x2d)
     end select
     psfact = 1.d0
     if (x2d(1)<10000.) psfact = 1.d+2 ! if unit of ps is hpa,
! it will be converted into Pa
     icnt = 0
     do ie = 1,nelemd
       do j = 1, np
         do i = 1, np
           icnt = icnt+1
           nt_org(i,j,nlev*nvar3d+vr,ie) = x2d(icnt)*psfact
         enddo
       enddo
     enddo
   enddo
   call closefile(io_t)
   call finiosystem(io_t)
   deallocate(x3d,x2d)
!
   end subroutine readnature
!===============================================================================
!===============================================================================
!  SUBROUTINE TempToM
!
!> @brief
!> - Transformation of temperature -> geopotential
!
!> @date 21JUN2013
!> - JH KWUN: Write the initial version based on HOMME code
!> @date 16JUL2013
!> - HJ SONG: Modify the argument from p and dp to ps, and plug in HybDA
!
!> @param[in] ps
!> @param[inout] tphi
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine temptom(t, q, m)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use physicalconstants, only: rgas=>kim_rair, ps0=>kim_p0
! 2. Variables
! 2-1. Input variables
!
   real(real_kind), dimension(np, np, nlev, nelemd), intent(in   ) :: q, t
! 2-2. Input/output variables
   real(real_kind), dimension(np, np, nlev, nelemd), intent(  out) :: m
! 2-3. Local variables
   real(real_kind), dimension(np, np, nlev, nelemd) :: p, dp, phi, phii, tphi
   real(real_kind) :: hkk, hkl
   integer(int_kind) :: i, j, k, ie
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!
   call barrierpar()
! 0. T -> Tv
   print*,'T->tv',par%iproc
   tphi = 0.d0
   m = 0.d0
   do ie = nets,nete
     do k = 1, nlev
       do j = 1, np
         do i = 1, np
           tphi(i,j,k,ie) = t(i,j,k,ie)*(1.d0+(rvap/rgas-1.d0)*q(i,j,k,ie))
         enddo
       enddo
     enddo
   enddo
   call barrierpar()
! 1. pint -> dp
   print*,'pint->dp',par%iproc
   dp = 0.d0
   do ie = nets,nete
     do j = 1, np
       do i = 1, np
         do k = 1, nlev
           dp(i,j,k,ie) = pint(i,j,k+1,ie)-pint(i,j,k,ie)
         enddo
       enddo
     enddo
   enddo
   call barrierpar()
! 2. Temperature(Tv) -> Phi
   print*,'t->phi',par%iproc
   phi = 0.d0
   phii = 0.d0
   hkk = 0.d0
   hkl = 0.d0
   do ie = nets,nete
     do j = 1, np
       do i = 1, np
         hkk = dp(i,j,nlev,ie)*0.5d0/p(i,j,nlev,ie)
         hkl = 2.d0*hkk
         phii(i,j,nlev,ie) = rgas*tphi(i,j,nlev,ie)*hkl
         phi(i,j,nlev,ie) = phis(i,j,ie)+rgas*tphi(i,j,nlev,ie)*hkk
       enddo
       do k = nlev-1, 2,-1
         do i = 1, np
           hkk = dp(i,j,k,ie)*0.5d0/p(i,j,k,ie)
           hkl = 2.d0*hkk
           phii(i,j,k,ie) = phii(i,j,k+1,ie)+rgas*tphi(i,j,k,ie)*hkl
           phi(i,j,k,ie) = phis(i,j,ie)+phii(i,j,k+1,ie)+rgas*tphi(i,j,k,ie)*hkk
         enddo
       enddo
       do i = 1, np
         hkk = 0.5d0*dp(i,j,1,ie)/p(i,j,1,ie)
         phi(i,j,1,ie) = phis(i,j,ie)+phii(i,j,2,ie)+rgas*tphi(i,j,1,ie)*hkk
       enddo
     enddo
   enddo
   call barrierpar()
! M = phi + RTv*ln(p)
! Tr = Tvr : virtual temperature at reference level
   print*,'phi->m',par%iproc
   do ie = nets,nete
     do k = 1, nlev
       do j = 1, np
         do i = 1, np
           m(i,j,k,ie) = phi(i,j,k,ie)+rgas*tr*dlog(p(i,j,k,ie))
         enddo
       enddo
     enddo
   enddo
!=============================================================================
   return
!
   end subroutine temptom
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine modeltoctrl_nochi(tl_bg_md, tl_bg_chi, gs_md, tl_bg_ct)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Variables
! 1-1. Input variables
!
   real(real_kind), dimension(np, np, nlevar, nelemd), intent(in   ) :: tl_bg_md, gs_md
   real(real_kind), dimension(np, np, nlev, nelemd), intent(in   ) :: tl_bg_chi
! 1-2. Output variables
   real(real_kind), dimension(np, np, nlevar, nelemd), intent(  out) :: tl_bg_ct
! 1-3. Local variables
   real(real_kind), dimension(np, np, nlevar, nelemd) :: md
   real(real_kind), dimension(np, np, nlev, 2, nelemd) :: tl_uv, bg_uv
   real(real_kind), dimension(np, np, nlev, nelemd) :: tl_t, tl_q, tl_m, bg_t, tl_rh, bg_rh
   real(real_kind), dimension(np, np, nelemd) :: tl_ps, bg_ps
   real(real_kind), dimension(np, np, nlev, nelemd) :: tl_ct
   real(real_kind), dimension(np, np, 2*nlev+1, nelemd) :: tl_var
   real(real_kind), dimension(nups_l, nlev) :: tl_ct_up
   real(real_kind), dimension(nups_l, 2*nlev+1) :: tl_var_up
   real(real_kind), dimension(2*nlev+1) :: tmp
   integer(int_kind) :: ie, i, j, k, iup
!=============================================================================
!
   md = gs_md
   tl_bg_ct = 0.d0
!=============================================================================
! Main body
!=============================================================================
   tl_uv(:,:,:,1,:) = tl_bg_md(:,:,1:nlev,:)
   tl_uv(:,:,:,2,:) = tl_bg_md(:,:,nlev+1:2*nlev,:)
   tl_t(:,:,:,:) = tl_bg_md(:,:,2*nlev+1:3*nlev,:)
   tl_rh(:,:,:,:) = tl_bg_md(:,:,3*nlev+1:4*nlev,:)
   tl_ps(:,:,:) = tl_bg_md(:,:,nvar3d*nlev+1,:)
!    bg_uv(:,:,:,1,:) = md(:,:,1:nlev,:)
!    bg_uv(:,:,:,2,:) = md(:,:,nlev+1:2*nlev,:)
   bg_uv(:,:,:,1,:) = md(:,:,1:nlev,:)
   bg_uv(:,:,:,2,:) = md(:,:,nlev+1:2*nlev,:)
   bg_t(:,:,:,:) = md(:,:,2*nlev+1:3*nlev,:)
   bg_rh(:,:,:,:) = md(:,:,3*nlev+1:4*nlev,:)
   bg_ps(:,:,:) = md(:,:,nvar3d*nlev+1,:)
!-----------------------------------------------------------------------------
   tl_bg_ct = tl_bg_md
   if (cv_opt_wind==1) then
     call calpsi(tl_uv,tl_ct)
     tl_bg_ct(:,:,1:nlev,:) = tl_ct
   endif
!    IF(cv_opt_psu == 1) THEN
!       !tl_ct = 0.D0
!       CALL CalPsb(tl_uv,&
!                   tl_ct)
!    ELSE IF(cv_opt_psu == 2) THEN
!       CALL CalPsbNlbe(tl_uv,bg_uv,bg_ps,&
!                       tl_ct)
!    END IF
!
!    tl_ct(:,:,nlev,:) = tl_ct(:,:,nlev,:)*bg_ps(:,:,:)
!
!    ! psu = ps - psb
!    tl_bg_ct(:,:,nvar3d*nlev+1,:) = tl_bg_md(:,:,nvar3d*nlev+1,:) - &
!                                    tl_ct(:,:,nlev,:)
   if (cv_opt_tu.le.2) then
     call psitorotwind(tl_bg_ct(:,:,1:nlev,:),tl_uv)
     if (cv_opt_tu==1) then
       call calphib(tl_uv,tl_ct)
     elseif (cv_opt_tu==2) then
       call calnlbe(tl_uv,bg_uv,tl_ct)
     endif
   endif !(cv_opt_tu.le.2)
   ierr = timingstart('EP2UP_fromModelToCtrl_nochi')
   call convertep2up_hyb(nlev,tl_ct,tl_ct_up)
   ierr = timingstop('EP2UP_fromModelToCtrl_nochi')
! Mb -> Tb,PSb
   tl_var_up = 0.d0
   do iup = 1,nups_l
     tmp = 0.d0
     tmp = matmul(tl_ct_up(iup,:),regress_up(:,iup,:)) ! 2*nlev+1
     do k = 1,2*nlev+1
       tl_var_up(iup,k) = tmp(k)
     enddo
   enddo
   ierr = timingstart('UP2EP_fromModelToCtrl_nochi')
   call convertup2ep_hyb(2*nlev+1,tl_var_up,tl_var)
   ierr = timingstop('UP2EP_fromModelToCtrl_nochi')
! Tu' = T' - Tb'
   tl_bg_ct(:,:,2*nlev+1:3*nlev,:) = tl_bg_md(:,:,2*nlev+1:3*nlev,:)-tl_var(:,:,1:nlev,:)
! psu = ps - psb
   tl_bg_ct(:,:,nvar3d*nlev+1,:) = tl_bg_md(:,:,nvar3d*nlev+1,:)-!tl_ct(:,:,nlev,:)
   tl_var(:,:,nlev+1,:)
!.. Make unbalanced velocity potential
!    tl_uv(:,:,:,1,:) = tl_bg_md(:,:,1:nlev,:)
!    tl_uv(:,:,:,2,:) = tl_bg_md(:,:,nlev+1:2*nlev,:)
!    CALL CalChi(tl_uv,&
!                tl_ct)
   tl_bg_ct(:,:,nlev+1:2*nlev,:) = tl_bg_chi(:,:,:,:)-tl_var(:,:,nlev+2:2*nlev+1,:)
!=============================================================================
   return
!
   end subroutine modeltoctrl_nochi
!===============================================================================
!  SUBROUTINE ModelToCtrl
!
!> @brief
!> - Transformation of model variables -> control variables
!
!> @date 21JUN2013
!> - JH KWUN: Write an initial version
!> @date 16JUL2013
!> - HJ SONG: Rewrite and plug in HybDA
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine modeltoctrl(tl_bg_md, gs_md, bg_uv, tl_bg_ct)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Variables
! 1-1. Input variables
!
   real(real_kind), dimension(np, np, nlevar, nelemd), intent(in   ) :: tl_bg_md, gs_md
   real(real_kind), dimension(np, np, nlev, 2, nelemd), intent(in   ) :: bg_uv
! 1-2. Output variables
   real(real_kind), dimension(np, np, nlevar, nelemd), intent(  out) :: tl_bg_ct
! 1-3. Local variables
   real(real_kind), dimension(np, np, nlevar, nelemd) :: md
   real(real_kind), dimension(np, np, nlev, 2, nelemd) :: tl_uv !, bg_uv
   real(real_kind), dimension(np, np, nlev, nelemd) :: tl_t, tl_q, tl_m, bg_t, tl_rh, bg_rh
   real(real_kind), dimension(np, np, nelemd) :: tl_ps, bg_ps
   real(real_kind), dimension(np, np, nlev, nelemd) :: tl_ct
   real(real_kind), dimension(np, np, 2*nlev+1, nelemd) :: tl_var
   real(real_kind), dimension(2*nlev+1) :: tmp
   real(real_kind), dimension(nups_l, nlev) :: tl_ct_up
   real(real_kind), dimension(nups_l, 2*nlev+1) :: tl_var_up
   integer(int_kind) :: ie, i, j, k, iup
!=============================================================================
!
   md = gs_md
   tl_bg_ct = 0.d0
!=============================================================================
! Main body
!=============================================================================
   tl_uv(:,:,:,1,:) = tl_bg_md(:,:,1:nlev,:)
   tl_uv(:,:,:,2,:) = tl_bg_md(:,:,nlev+1:2*nlev,:)
   tl_t(:,:,:,:) = tl_bg_md(:,:,2*nlev+1:3*nlev,:)
   tl_rh(:,:,:,:) = tl_bg_md(:,:,3*nlev+1:4*nlev,:)
   tl_ps(:,:,:) = tl_bg_md(:,:,nvar3d*nlev+1,:)
!    bg_uv(:,:,:,1,:) = md(:,:,1:nlev,:)
!    bg_uv(:,:,:,2,:) = md(:,:,nlev+1:2*nlev,:)
   bg_t(:,:,:,:) = md(:,:,2*nlev+1:3*nlev,:)
   bg_rh(:,:,:,:) = md(:,:,3*nlev+1:4*nlev,:)
   bg_ps(:,:,:) = md(:,:,nvar3d*nlev+1,:)
!-----------------------------------------------------------------------------
   tl_bg_ct = tl_bg_md
   if (cv_opt_wind==1) then
     call calpsi(tl_uv,tl_ct)
     tl_bg_ct(:,:,1:nlev,:) = tl_ct
   endif
!    IF(cv_opt_psu == 1) THEN
!       !tl_ct = 0.D0
!       CALL CalPsb(tl_uv,&
!                   tl_ct)
!
!    ELSE IF(cv_opt_psu == 2) THEN
!       CALL CalPsbNlbe(tl_uv,bg_uv,bg_ps,&
!                       tl_ct)
!    END IF
!    tl_ct(:,:,nlev,:) = tl_ct(:,:,nlev,:)*bg_ps(:,:,:)
!
!    ! psu = ps - psb
!    tl_bg_ct(:,:,nvar3d*nlev+1,:) = tl_bg_md(:,:,nvar3d*nlev+1,:) - &
!                                    tl_ct(:,:,nlev,:)
   if (cv_opt_tu.le.2) then
     call psitorotwind(tl_bg_ct(:,:,1:nlev,:),tl_uv)
     if (cv_opt_tu==1) then
       call calphib(tl_uv,tl_ct)
     elseif (cv_opt_tu==2) then
       call calnlbe(tl_uv,bg_uv,tl_ct)
     endif
   endif !(cv_opt_tu.le.2)
   ierr = timingstart('EP2UP_fromModelToCtrl')
   call convertep2up_hyb(nlev,tl_ct,tl_ct_up)
   ierr = timingstop('EP2UP_fromModelToCtrl')
! Mb -> Tb,PSb
   tl_var_up = 0.d0
   do iup = 1,nups_l
     tmp = 0.d0
     tmp = matmul(tl_ct_up(iup,:),regress_up(:,iup,:)) ! 2*nlev+1
     do k = 1,2*nlev+1
       tl_var_up(iup,k) = tmp(k)
     enddo
   enddo
   ierr = timingstart('UP2EP_fromModelToCtrl')
   call convertup2ep_hyb(2*nlev+1,tl_var_up,tl_var)
   ierr = timingstop('UP2EP_fromModelToCtrl')
! Tu' = T' - Tb'
   tl_bg_ct(:,:,2*nlev+1:3*nlev,:) = tl_bg_md(:,:,2*nlev+1:3*nlev,:)-tl_var(:,:,1:nlev,:)
! psu = ps - psb
   tl_bg_ct(:,:,nvar3d*nlev+1,:) = tl_bg_md(:,:,nvar3d*nlev+1,:)-!tl_ct(:,:,nlev,:)
   tl_var(:,:,nlev+1,:)
!.. Make unbalanced velocity potential
   tl_uv(:,:,:,1,:) = tl_bg_md(:,:,1:nlev,:)
   tl_uv(:,:,:,2,:) = tl_bg_md(:,:,nlev+1:2*nlev,:)
   call calchi(tl_uv,tl_ct)
   tl_bg_ct(:,:,nlev+1:2*nlev,:) = tl_ct-tl_var(:,:,nlev+2:2*nlev+1,:)
!=============================================================================
   return
!
   end subroutine modeltoctrl
!===============================================================================
!  SUBROUTINE PsiToRotWind
!
!> @brief
!> - Transformation of psi -> rotational winds
!
!> @date 21JUN2013
!> - HJ SONG & JH KWUN: Devise and code
!> @date 16JUL2013
!> - HJ SONG: Clean up and plug in to HybDA
!> @date 02DEC2015
!> - HJ SONG: Implement
!
!> @param[in] psi
!> @param[out] urot
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine psitorotwind(psi, urot)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use derivative, only: gradient_sphere
   use edge, only: edgebuffer_t, initedgebuffer, freeedgebuffer, edgevpack, edgevunpack
   use bndry, only: bndry_exchangev
! 2. Variables
! 2-1. Input variables
!
   real(real_kind), dimension(np, np, nlev, nelemd), intent(in   ) :: psi
! 2-2. Output variables
   real(real_kind), dimension(np, np, nlev, 2, nelemd), intent(  out) :: urot
! 2-3. Local variables
   type(edgebuffer_t) :: edge2
   real(real_kind), dimension(np, np, 2) :: utmp
   integer(int_kind) :: i, j, k, ie
!=============================================================================
!=============================================================================
! B. Psi -> rotational winds
!=============================================================================
!
   call initedgebuffer(edge2,2*nlev)
   do k = 1,nlev
     utmp = 0.d0
     do ie = nets,nete
       utmp = gradient_sphere(psi(:,:,k,ie),deriv(hybrid%ithr),elem(ie)%dinv)
       urot(:,:,k,1,ie) =-utmp(:,:,2)
       urot(:,:,k,2,ie) = utmp(:,:,1)
     enddo
   enddo
   do ie = nets,nete
     do k = 1, nlev
       urot(:,:,k,1,ie) = elem(ie)%spheremp(:,:)*urot(:,:,k,1,ie)
       urot(:,:,k,2,ie) = elem(ie)%spheremp(:,:)*urot(:,:,k,2,ie)
     enddo
     call edgevpack(edge2,urot(:,:,:,1,ie),nlev,0,elem(ie)%desc)
     call edgevpack(edge2,urot(:,:,:,2,ie),nlev,nlev,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edge2)
   do ie = nets,nete
     call edgevunpack(edge2,urot(:,:,:,1,ie),nlev,0,elem(ie)%desc)
     call edgevunpack(edge2,urot(:,:,:,2,ie),nlev,nlev,elem(ie)%desc)
     do k = 1,nlev
       urot(:,:,k,1,ie) = elem(ie)%rspheremp(:,:)*urot(:,:,k,1,ie)
       urot(:,:,k,2,ie) = elem(ie)%rspheremp(:,:)*urot(:,:,k,2,ie)
     enddo
   enddo
   call freeedgebuffer(edge2)
!=============================================================================
   return
!
   end subroutine psitorotwind
!===============================================================================
!  SUBROUTINE CalPsi
!
!> @brief
!> - Transformation of horzontal winds -> stream function
!
!> @date 21JUN2013
!> - HJ SONG & JH KWUN: Devise and code
!> @date 16JUL2013
!> - HJ SONG: Clean up and plug in HybDA
!
!> @param[in] uv
!> @param[out] psi
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine calpsi(uv, psi)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use derivative, only: vorticity_sphere_wk
! 2. Variables
! 2-1. Input variables
!
   real(real_kind), dimension(np, np, nlev, 2, nelemd), intent(in   ) :: uv
! 2-2. Output variables
   real(real_kind), dimension(np, np, nlev, nelemd), intent(  out) :: psi
! 2-3. Local variables
   type(edgebuffer_t) :: edge2
   real(real_kind), dimension(np, np, nlev, nelemd) :: b
   real(real_kind), dimension(np, np, nlev, 2) :: uv_ll
   integer(int_kind) :: i, j, k, ie
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!
   call initedgebuffer(edge2,nlev)
   do ie = nets,nete
     do k = 1, nlev
       do j = 1, np
         do i = 1, np
           uv_ll(i,j,k,1) = uv(i,j,k,1,ie)
           uv_ll(i,j,k,2) = uv(i,j,k,2,ie)
         enddo
       enddo
       b(:,:,k,ie) = vorticity_sphere_wk(uv_ll(:,:,k,:),deriv(hybrid%ithr),elem(ie))
     enddo
     call edgevpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edge2)
   do ie = nets,nete
     call edgevunpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
     do k = 1,nlev
       b(:,:,k,ie) = elem(ie)%rspheremp(:,:)*b(:,:,k,ie)
     enddo
   enddo
!    psi = 0.D0
!    CALL InvLapCg(b,  &
!                  psi)
   call invlapspec(b,nlev,psi) ! out:solution
   call freeedgebuffer(edge2)
!=============================================================================
   return
!
   end subroutine calpsi
!===============================================================================
!  SUBROUTINE CalChi
!
!> @brief
!> - Transformation of horzontal winds -> velocity potential
!
!> @date 21JUN2013
!> - HJ SONG & JH KWUN: Devise and code
!> @date 16JUL2013
!> - HJ SONG: Clean up and plug in HybDA
!
!> @param[in] uv
!> @param[out] chi
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine calchi(uv, chi)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use derivative, only: divergence_sphere_wk
! 2. Variables
! 2-1. Input variables
!
   real(real_kind), dimension(np, np, nlev, 2, nelemd), intent(in   ) :: uv
! 2-2. Output variables
   real(real_kind), dimension(np, np, nlev, nelemd), intent(  out) :: chi
! 2-3. Local variables
   type(edgebuffer_t) :: edge2
   real(real_kind), dimension(np, np, nlev, nelemd) :: b
   real(real_kind), dimension(np, np, nlev, 2) :: uv_ll
   integer(int_kind) :: i, j, k, ie
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!
   call initedgebuffer(edge2,nlev)
   do ie = nets,nete
     do k = 1, nlev
       do j = 1, np
         do i = 1, np
           uv_ll(i,j,k,1) = uv(i,j,k,1,ie)
           uv_ll(i,j,k,2) = uv(i,j,k,2,ie)
         enddo
       enddo
       b(:,:,k,ie) = divergence_sphere_wk(uv_ll(:,:,k,:),deriv(hybrid%ithr),elem(ie))
     enddo
     call edgevpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edge2)
   do ie = nets,nete
     call edgevunpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
     do k = 1,nlev
       b(:,:,k,ie) = elem(ie)%rspheremp(:,:)*b(:,:,k,ie)
     enddo
   enddo
!    chi = 0.D0
!    CALL InvLapCg(b,&
!                  chi)
   call invlapspec(b,nlev,chi)
   call freeedgebuffer(edge2)
!=============================================================================
   return
!
   end subroutine calchi
!===============================================================================
!  SUBROUTINE CalPhib
!
!> @brief
!> - Transformation of horzontal winds -> balanced geopotential
!
!> @date 21JUN2013
!> - HJ SONG & JH KWUN: Devise and code
!> @date 16JUL2013
!> - HJ SONG: Clean up and plug in HybDA
!
!> @param[in] uv
!> @param[out] phib
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine calphib(uv, phib)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use derivative, only: divergence_sphere_wk
! 2. Variables
! 2-1. Input variables
!
   real(real_kind), dimension(np, np, nlev, 2, nelemd), intent(in   ) :: uv
! 2-2. Output variables
   real(real_kind), dimension(np, np, nlev, nelemd), intent(  out) :: phib
! 2-3. Local variables
   type(edgebuffer_t) :: edgep2
   real(real_kind), dimension(np, np, nlev, nelemd) :: b
   real(real_kind), dimension(np, np, nlev, 2) :: uv_ll, ucor
   integer(int_kind) :: i, j, k, ie
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!
   call initedgebuffer(edgep2,nlev)
   do ie = nets,nete
     do k = 1, nlev
       do j = 1, np
         do i = 1, np
           uv_ll(i,j,k,1) = uv(i,j,k,1,ie)
           uv_ll(i,j,k,2) = uv(i,j,k,2,ie)
         enddo
       enddo
       do j = 1, np
         do i = 1, np
           ucor(i,j,k,1) = elem(ie)%fcor(i,j)*uv_ll(i,j,k,2)
           ucor(i,j,k,2) = elem(ie)%fcor(i,j)*(-uv_ll(i,j,k,1))
         enddo
       enddo
       b(:,:,k,ie) = divergence_sphere_wk(ucor(:,:,k,:),deriv(hybrid%ithr),elem(ie))
     enddo
     call edgevpack(edgep2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edgep2)
   do ie = nets,nete
     call edgevunpack(edgep2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
     do k = 1,nlev
       b(:,:,k,ie) = elem(ie)%rspheremp(:,:)*b(:,:,k,ie)
     enddo
   enddo
!    phib = 0.D0
!    CALL InvLapCg(b,   &
!                  phib)
   call invlapspec(b,nlev,phib)
   call freeedgebuffer(edgep2)
!=============================================================================
   return
!
   end subroutine calphib
!===============================================================================
!  SUBROUTINE CalPsb
!
!> @brief
!> - Transformation of horzontal winds -> balanced surface pressure
!
!> @date 21JUN2013
!> - HJ SONG & JH KWUN: Devise and code
!> @date 16JUL2013
!> - HJ SONG: Clean up and plug in HybDA
!> @date  9SEP2013
!> - JH KWUN: Code
!
!> @param[in] uv
!> @param[out] psb
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine calpsb(uv, psb)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use derivative, only: divergence_sphere_wk
! 2. Variables
! 2-1. Input variables
!
   real(real_kind), dimension(np, np, nlev, 2, nelemd), intent(in   ) :: uv
   real(real_kind), dimension(np, np, nelemd) :: ps
! 2-2. Output variables
   real(real_kind), dimension(np, np, nlev, nelemd), intent(  out) :: psb
! 2-3. Local variables
   type(edgebuffer_t) :: edge2
   real(real_kind), dimension(np, np, nlev, nelemd) :: b
   real(real_kind), dimension(np, np, nlev, 2) :: uv_ll, ucor
   integer(int_kind) :: i, j, k, ie
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!
   call initedgebuffer(edge2,1)
   b(:,:,:,:) = 0.d0
   do ie = nets,nete
     do j = 1, np
       do i = 1, np
         uv_ll(i,j,nlev,1) = uv(i,j,nlev,1,ie)
         uv_ll(i,j,nlev,2) = uv(i,j,nlev,2,ie)
       enddo
     enddo
     do j = 1, np
       do i = 1, np
         ucor(i,j,nlev,1) = elem(ie)%fcor(i,j)*uv_ll(i,j,nlev,2)
         ucor(i,j,nlev,2) = elem(ie)%fcor(i,j)*(-uv_ll(i,j,nlev,1))
       enddo
     enddo
     b(:,:,nlev,ie) = divergence_sphere_wk(ucor(:,:,nlev,:),deriv(hybrid%ithr),elem(ie))
     call edgevpack(edge2,b(:,:,nlev,ie),1,0,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edge2)
   do ie = nets,nete
     call edgevunpack(edge2,b(:,:,nlev,ie),1,0,elem(ie)%desc)
     b(:,:,nlev,ie) = elem(ie)%rspheremp(:,:)*b(:,:,nlev,ie)
   enddo
   b(:,:,nlev,:) = b(:,:,nlev,:)/(rgas*tr)
!    psb(:,:,:,:) = 0.D0
!    CALL InvLapCg(b,  &
!                  psb)
!.. Spectral inverse
   call invlapspec(b(:,:,nlev,:),1,psb(:,:,nlev,:))
   call freeedgebuffer(edge2)
!=============================================================================
   return
!
   end subroutine calpsb
!===============================================================================
!  SUBROUTINE CalPsbNlbe
!
!> @brief
!> - Transformation of horzontal winds -> balanced surface pressure
!
!> @date 21JUN2013
!> - HJ SONG & JH KWUN: Devise and code
!> @date 16JUL2013
!> - HJ SONG: Clean up and plug in HybDA
!> @date  9SEP2013
!> - JH KWUN: Code
!
!> @param[in] uv
!> @param[out] psb
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine calpsbnlbe(uv, uvb, ps, psb)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use derivative, only: divergence_sphere_wk, gradient_sphere
! 2. Variables
! 2-1. Input variables
!
   real(real_kind), dimension(np, np, nlev, 2, nelemd), intent(in   ) :: uv, uvb
   real(real_kind), dimension(np, np, nelemd), intent(in   ) :: ps
! 2-2. Output variables
   real(real_kind), dimension(np, np, nlev, nelemd), intent(  out) :: psb
! 2-3. Local variables
   type(edgebuffer_t) :: edgep2
   real(real_kind), dimension(np, np, nlev, 2, nelemd) :: uadv
   real(real_kind), dimension(np, np, nlev, nelemd) :: b
   real(real_kind), dimension(np, np, nlev, 2) :: uv_ll, ucor, adv
   real(real_kind), dimension(np, np, 2) :: gradv, gradvb, ugradv, ugradvb
   integer(int_kind) :: i, j, k, ie, comp
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!
   call initedgebuffer(edgep2,2)
   b(:,:,:,:) = 0.d0
! 2. advection
! - (uvb*grad uv + uv*grad uvb)
   uadv = 0.d0
   do ie = nets,nete
     gradv = 0.d0
     gradvb = 0.d0
     do k = nlev,nlev
       do comp = 1, 2
         gradv = gradient_sphere(uv(:,:,k,comp,ie),deriv(hybrid%ithr),elem(ie)%dinv)
         gradvb = gradient_sphere(uvb(:,:,k,comp,ie),deriv(hybrid%ithr),elem(ie)%dinv)
         ugradv(:,:,comp) = uvb(:,:,k,1,ie)*gradv(:,:,1)+uvb(:,:,k,2,ie)*gradv(:,:,2)
         ugradvb(:,:,comp) = uv(:,:,k,1,ie)*gradvb(:,:,1)+uv(:,:,k,2,ie)*gradvb(:,:,2)
       enddo
       uadv(:,:,k,1,ie) = ugradv(:,:,1)+ugradvb(:,:,1)
       uadv(:,:,k,2,ie) = ugradv(:,:,2)+ugradvb(:,:,2)
     enddo
   enddo
! DSS
   do ie = nets,nete
     do k = nlev, nlev
       uadv(:,:,k,1,ie) = elem(ie)%spheremp(:,:)*uadv(:,:,k,1,ie)
       uadv(:,:,k,2,ie) = elem(ie)%spheremp(:,:)*uadv(:,:,k,2,ie)
     enddo
     call edgevpack(edgep2,uadv(:,:,nlev,1,ie),1,0,elem(ie)%desc)
     call edgevpack(edgep2,uadv(:,:,nlev,2,ie),1,1,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edgep2)
   do ie = nets,nete
     call edgevunpack(edgep2,uadv(:,:,nlev,1,ie),1,0,elem(ie)%desc)
     call edgevunpack(edgep2,uadv(:,:,nlev,2,ie),1,1,elem(ie)%desc)
     do k = nlev,nlev
       uadv(:,:,k,1,ie) = elem(ie)%rspheremp(:,:)*uadv(:,:,k,1,ie)
       uadv(:,:,k,2,ie) = elem(ie)%rspheremp(:,:)*uadv(:,:,k,2,ie)
     enddo
   enddo
! 1. linear balance equation
   do ie = nets,nete
     ucor = 0.d0
     do j = 1,np
       do i = 1, np
         ucor(i,j,nlev,1) = elem(ie)%fcor(i,j)*uv(i,j,nlev,2,ie)
         ucor(i,j,nlev,2) = elem(ie)%fcor(i,j)*(-uv(i,j,nlev,1,ie))
         adv(i,j,nlev,1) =-uadv(i,j,nlev,1,ie)
         adv(i,j,nlev,2) =-uadv(i,j,nlev,2,ie)
       enddo
     enddo
     ucor(:,:,nlev,1) = ucor(:,:,nlev,1)+adv(:,:,nlev,1)
     ucor(:,:,nlev,2) = ucor(:,:,nlev,2)+adv(:,:,nlev,2)
     b(:,:,nlev,ie) = divergence_sphere_wk(ucor(:,:,nlev,:),deriv(hybrid%ithr),elem(ie))
     call edgevpack(edgep2,b(:,:,nlev,ie),1,0,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edgep2)
   do ie = nets,nete
     call edgevunpack(edgep2,b(:,:,nlev,ie),1,0,elem(ie)%desc)
     b(:,:,nlev,ie) = elem(ie)%rspheremp(:,:)*b(:,:,nlev,ie)
   enddo
   b(:,:,nlev,:) = b(:,:,nlev,:)/(rgas*tr)
!    psb(:,:,:,:) = 0.D0
!    CALL InvLapCg(b,  &
!                  psb)
!.. Spectral inverse
   call invlapspec(b(:,:,nlev,:),1,psb(:,:,nlev,:))
   call freeedgebuffer(edgep2)
!=============================================================================
   return
!
   end subroutine calpsbnlbe
!===============================================================================
!  SUBROUTINE CalNLBE
!
!> @brief
!> - Transformation of horzontal winds -> balanced geopotential + advection term
!
!> @date 24SEP2013
!> - HJ SONG & JH KWUN: Devise and code
!
!> @param[in] uv
!> @param[out] phib
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine calnlbe(uv, uvb, phib)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use derivative, only: divergence_sphere_wk, gradient_sphere
! 2. Variables
! 2-1. Input variables
!
   real(real_kind), dimension(np, np, nlev, 2, nelemd), intent(in   ) :: uv, uvb
! 2-2. Output variables
   real(real_kind), dimension(np, np, nlev, nelemd), intent(  out) :: phib
! 2-3. Local variables
   type(edgebuffer_t) :: edge2
   real(real_kind), dimension(np, np, nlev, 2, nelemd) :: uadv
   real(real_kind), dimension(np, np, nlev, nelemd) :: b
   real(real_kind), dimension(np, np, nlev, 2) :: uv_ll, ucor, adv
   real(real_kind), dimension(np, np, 2) :: gradv, gradvb, ugradv, ugradvb
   integer(int_kind) :: i, j, k, ie, comp
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
!
   call initedgebuffer(edge2,2*nlev)
   b = 0.d0
! 1. advection
! - (uvb*grad uv + uv*grad uvb)
   uadv = 0.d0
   do ie = nets,nete
     gradv = 0.d0
     gradvb = 0.d0
     do k = 1,nlev
       do comp = 1, 2
         gradv = gradient_sphere(uv(:,:,k,comp,ie),deriv(hybrid%ithr),elem(ie)%dinv)
         gradvb = gradient_sphere(uvb(:,:,k,comp,ie),deriv(hybrid%ithr),elem(ie)%dinv)
         ugradv(:,:,comp) = uvb(:,:,k,1,ie)*gradv(:,:,1)+uvb(:,:,k,2,ie)*gradv(:,:,2)
         ugradvb(:,:,comp) = uv(:,:,k,1,ie)*gradvb(:,:,1)+uv(:,:,k,2,ie)*gradvb(:,:,2)
       enddo
       uadv(:,:,k,1,ie) = ugradv(:,:,1)+ugradvb(:,:,1)
       uadv(:,:,k,2,ie) = ugradv(:,:,2)+ugradvb(:,:,2)
     enddo
   enddo
! DSS
   do ie = nets,nete
     do k = 1, nlev
       uadv(:,:,k,1,ie) = elem(ie)%spheremp(:,:)*uadv(:,:,k,1,ie)
       uadv(:,:,k,2,ie) = elem(ie)%spheremp(:,:)*uadv(:,:,k,2,ie)
     enddo
     call edgevpack(edge2,uadv(:,:,:,1,ie),nlev,0,elem(ie)%desc)
     call edgevpack(edge2,uadv(:,:,:,2,ie),nlev,nlev,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edge2)
   do ie = nets,nete
     call edgevunpack(edge2,uadv(:,:,:,1,ie),nlev,0,elem(ie)%desc)
     call edgevunpack(edge2,uadv(:,:,:,2,ie),nlev,nlev,elem(ie)%desc)
     do k = 1,nlev
       uadv(:,:,k,1,ie) = elem(ie)%rspheremp(:,:)*uadv(:,:,k,1,ie)
       uadv(:,:,k,2,ie) = elem(ie)%rspheremp(:,:)*uadv(:,:,k,2,ie)
     enddo
   enddo
! 1. linear balance equation
   do ie = nets,nete
     do k = 1, nlev
       ucor = 0.d0
       do j = 1,np
         do i = 1, np
           ucor(i,j,k,1) = elem(ie)%fcor(i,j)*uv(i,j,k,2,ie)
           ucor(i,j,k,2) = elem(ie)%fcor(i,j)*(-uv(i,j,k,1,ie))
           adv(i,j,k,1) =-uadv(i,j,k,1,ie)
           adv(i,j,k,2) =-uadv(i,j,k,2,ie)
         enddo
       enddo
       ucor(:,:,k,1) = ucor(:,:,k,1)+adv(:,:,k,1)
       ucor(:,:,k,2) = ucor(:,:,k,2)+adv(:,:,k,2)
       b(:,:,k,ie) = divergence_sphere_wk(ucor(:,:,k,:),deriv(hybrid%ithr),elem(ie))
     enddo
     call edgevpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edge2)
   do ie = nets,nete
     call edgevunpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
     do k = 1,nlev
       b(:,:,k,ie) = elem(ie)%rspheremp(:,:)*b(:,:,k,ie)
     enddo
   enddo
!    phib = 0.D0
!    CALL InvLapCg(b,   &
!                  phib)
!.. Spectral inverse
   call invlapspec(b,nlev,phib)
   call freeedgebuffer(edge2)
!=============================================================================
   return
!
   end subroutine calnlbe
!===============================================================================
!  SUBROUTINE TlmTempToPhi
!
!> @brief
!> - Transformation of temperature -> geopotential
!
!> @date 21JUN2013
!> - JH KWUN: Write the initial version based on HOMME code
!> @date 16JUL2013
!> - HJ SONG: Modify the argument from p and dp to ps, and plug in HybDA
!
!> @param[in] ps
!> @param[inout] tphi
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine tlmtemptophi(tl_ps, tl_q, tl_t, ps, q, t, tl_phi)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use physicalconstants, only: rvap=>kim_rvap
! 2. Variables
! 2-1. Input variables
!
   real(real_kind), dimension(np, np, nelemd), intent(in   ) :: ps, tl_ps
   real(real_kind), dimension(np, np, nlev, nelemd), intent(in   ) :: q, tl_q, t, tl_t
   real(real_kind), dimension(np, np, nlev, nelemd) :: tphi
! 2-2. Input/output variables
   real(real_kind), dimension(np, np, nlev, nelemd), intent(  out) :: tl_phi
! 2-3. Local variables
   real(real_kind), dimension(np, np, nlev+1, nelemd) :: ph, tl_ph
   real(real_kind), dimension(np, np, nlev, nelemd) :: p, dp, phi, phii, tl_p, tl_dp, tl_phii, tl_tphi
   real(real_kind) :: hkk, hkl, tl_hkk, tl_hkl
   integer(int_kind) :: i, j, k, ie
!=============================================================================
!=============================================================================
! B. Main body
!=============================================================================
! 0. T -> Tv
!
   tphi = 0
   tl_tphi = 0
   do ie = nets,nete
     do k = 1, nlev
       do j = 1, np
         do i = 1, np
           tphi(i,j,k,ie) = t(i,j,k,ie)*(1.d0+(rvap/rgas-1.d0)*q(i,j,k,ie))
           tl_tphi(i,j,k,ie) = tl_t(i,j,k,ie)+tl_t(i,j,k,ie)*q(i,j,k,ie)*(rvap/rgas-1.d0)+t(i,j,k,ie)*tl_q(i,j,k,ie)*(rvap/rgas-1.d0)
         enddo
       enddo
     enddo
   enddo
! 1. ps -> p and dp
   do ie = nets,nete
     do k = 1, nlev+1
       do j = 1, np
         do i = 1, np
           ph(i,j,k,ie) = hyai(k)*ps0+hybi(k)*ps(i,j,ie)
         enddo
       enddo
     enddo
   enddo
   do ie = nets,nete
     do k = 1, nlev
       do j = 1, np
         do i = 1, np
           p(i,j,k,ie) = hyam(k)*ps0+hybm(k)*ps(i,j,ie)
           dp(i,j,k,ie) = ph(i,j,k+1,ie)-ph(i,j,k,ie)
         enddo
       enddo
     enddo
   enddo
   do ie = nets,nete
     do k = 1, nlev+1
       do j = 1, np
         do i = 1, np
           tl_ph(i,j,k,ie) = hybi(k)*tl_ps(i,j,ie)
         enddo
       enddo
     enddo
   enddo
   do ie = nets,nete
     do k = 1, nlev
       do j = 1, np
         do i = 1, np
           tl_p(i,j,k,ie) = hybm(k)*tl_ps(i,j,ie)
           tl_dp(i,j,k,ie) = tl_ph(i,j,k+1,ie)-tl_ph(i,j,k,ie)
         enddo
       enddo
     enddo
   enddo
! 2. Temperature (Tv) -> Phi
   phi = 0.d0
   phii = 0.d0
   hkk = 0.d0
   tl_phi = 0.d0
   tl_phii = 0.d0
   tl_hkk = 0.d0
   do ie = nets,nete
     do j = 1, np
       do i = 1, np
         hkk = dp(i,j,nlev,ie)*0.5d0/p(i,j,nlev,ie)
         hkl = 2.d0*hkk
         tl_hkk = tl_dp(i,j,nlev,ie)*0.5d0/p(i,j,nlev,ie)-dp(i,j,nlev,ie)*0.5d0*tl_p(i,j,nlev,ie)/(p(i,j,nlev,ie)**2.d0)
         tl_hkl = 2.d0*tl_hkk
         tl_phii(i,j,nlev,ie) = rgas*tl_tphi(i,j,nlev,ie)*hkl+rgas*tphi(i,j,nlev,ie)*tl_hkl
         tl_phi(i,j,nlev,ie) = rgas*tl_tphi(i,j,nlev,ie)*hkk+rgas*tphi(i,j,nlev,ie)*tl_hkk
       enddo
       do k = nlev-1, 2,-1
         do i = 1, np
           hkk = dp(i,j,k,ie)*0.5d0/p(i,j,k,ie)
           hkl = 2.d0*hkk
           tl_hkk = tl_dp(i,j,k,ie)*0.5d0/p(i,j,k,ie)-dp(i,j,k,ie)*0.5d0*tl_p(i,j,k,ie)/(p(i,j,k,ie)**2.d0) ! tl_hkl = 2.d0*tl_hkk
           tl_phii(i,j,k,ie) = tl_phii(i,j,k+1,ie)+rgas*tl_tphi(i,j,k,ie)*hkl+rgas*tphi(i,j,k,ie)*tl_hkl
           tl_phi(i,j,k,ie) = tl_phii(i,j,k+1,ie)+rgas*tl_tphi(i,j,k,ie)*hkk+rgas*tphi(i,j,k,ie)*tl_hkk
         enddo
       enddo
       do i = 1, np
         hkk = 0.5d0*dp(i,j,1,ie)/p(i,j,1,ie)
         tl_hkk = tl_dp(i,j,1,ie)*0.5d0/p(i,j,1,ie)-dp(i,j,1,ie)*0.5d0*tl_p(i,j,1,ie)/(p(i,j,1,ie)**2.d0) ! tl_phi(i,j,1,ie) = tl_phii(i,j,2,ie)+rgas*tl_tphi(i,j,1,ie)*hkk+rgas*tphi(i,j,1,ie)*tl_hkk
       enddo
     enddo
   enddo
! M = phi + RTv*ln(p)
! Tr = Tvr : virtual temperature at reference level
   do ie = nets,nete
     do k = 1, nlev
       do j = 1, np
         do i = 1, np
!phi(i,j,k,ie) = phi(i,j,k,ie) + rgas*Tr*DLOG(p(i,j,k,ie))
           tl_phi(i,j,k,ie) = tl_phi(i,j,k,ie)+rgas*tr*tl_p(i,j,k,ie)/p(i,j,k,ie)
         enddo
       enddo
     enddo
   enddo
!tphi = phi
!tl_tphi = tl_phi
!=============================================================================
   return
!
   end subroutine tlmtemptophi
!===============================================================================
! SUBROUTINE InvLapCg
!
!> @brief
!> - Inverse a laplace operator defined on cubed-sphere grids
!> - using conjugate gradient
!
!> @date 27MAR2013
!> - HJ SONG: Devise and code
!> @date 15MAY2013
!> - JH KWUN: Extract the code for inverse of Laplace operator
!> @date 07JUN2013
!> - HJ SONG: Organize the code and extend it into 3-dimension
!> @date 12JUN2013
!> - JH KWUN: Debug the usage of boundary exchange modules
!> @date 21JUL2013
!> - HJ SONG: Normalize b and A, and plug in HybDA
!
!> @param[in] b_org
!> @param[out] x
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine invlapcg(b_org, x)
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use derivative, only: laplace_sphere_wk_hyb
   use kiapsparallel, only: global_shared_buf, global_shared_sum
   use globalnorms, only: wrap_repro_sum
! 2. Variables
! 2-1. Input
!
   real(real_kind), dimension(np, np, nlev, nelemd), intent(in   ) :: b_org
! 2-2. Input/Output
   real(real_kind), dimension(np, np, nlev, nelemd), intent(inout) :: x
! 2-3. Local
   type(edgebuffer_t) :: edge2
   real(real_kind), dimension(np, np, nlev, nelemd) :: ap, r, p, b
   real(real_kind), dimension(:,:), pointer :: viscosity=>null()
   real(real_kind) :: rsold, rsnew, alpha, b_norm
   integer(int_kind) :: i, j, k, ie, iter
!=============================================================================
!=============================================================================
! B. Prepare inverse of Laplace operator
!=============================================================================
! 1. Initialize edge buffer
!
   call initedgebuffer(edge2,2*nlev)
! 2. Make b
   do ie = nets,nete
     do k = 1, nlev
       b(:,:,k,ie) = elem(ie)%rspheremp(:,:)*b_org(:,:,k,ie)
     enddo
     call edgevpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edge2)
   do ie = nets,nete
     call edgevunpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
     do k = 1,nlev
       b(:,:,k,ie) = elem(ie)%spheremp(:,:)*b(:,:,k,ie)
     enddo
   enddo
   do ie = nets,nete
     global_shared_buf(ie,1) = sum(b(:,:,:,ie)**2.d0)
   enddo
   call wrap_repro_sum(nvars = 1,comm = hybrid%par%comm)
   b_norm = dsqrt(global_shared_sum(1))
   b = b/b_norm
! 3. Make Ax
   do ie = nets,nete
     do k = 1, nlev
       x(:,:,k,ie) = elem(ie)%spheremp(:,:)*x(:,:,k,ie)
     enddo
     call edgevpack(edge2,x(:,:,:,ie),nlev,0,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edge2)
   do ie = nets,nete
     call edgevunpack(edge2,x(:,:,:,ie),nlev,0,elem(ie)%desc)
     do k = 1,nlev
       x(:,:,k,ie) = elem(ie)%rspheremp(:,:)*x(:,:,k,ie)
       ap(:,:,k,ie) = laplace_sphere_wk_hyb(x(:,:,k,ie),deriv(hybrid%ithr),elem(ie),viscosity)
     enddo
   enddo
   do ie = nets,nete
     do k = 1, nlev
       ap(:,:,k,ie) = elem(ie)%rspheremp(:,:)*ap(:,:,k,ie)
     enddo
     call edgevpack(edge2,ap(:,:,:,ie),nlev,0,elem(ie)%desc)
   enddo
   call bndry_exchangev(hybrid,edge2)
   do ie = nets,nete
     call edgevunpack(edge2,ap(:,:,:,ie),nlev,0,elem(ie)%desc)
     do k = 1,nlev
       ap(:,:,k,ie) = elem(ie)%spheremp(:,:)*ap(:,:,k,ie)
     enddo
   enddo
   ap = ap/b_norm
! 4. Initialize r,p,and rsold
   r = b-ap
   p = r
   do ie = nets,nete
     global_shared_buf(ie,1) = sum(r(:,:,:,ie)**2.d0)
   enddo
   call wrap_repro_sum(nvars = 1,comm = hybrid%par%comm)
   rsold = global_shared_sum(1)
!=============================================================================
!=============================================================================
! C. Iteration of conjugate gradient
!=============================================================================
   do iter = 1,niter_il
! 1. Make Ap
     do ie = nets, nete
       do k = 1, nlev
         p(:,:,k,ie) = elem(ie)%spheremp(:,:)*p(:,:,k,ie)
       enddo
       call edgevpack(edge2,p(:,:,:,ie),nlev,0,elem(ie)%desc)
     enddo
     call bndry_exchangev(hybrid,edge2)
     do ie = nets,nete
       call edgevunpack(edge2,p(:,:,:,ie),nlev,0,elem(ie)%desc)
       do k = 1,nlev
         p(:,:,k,ie) = elem(ie)%rspheremp(:,:)*p(:,:,k,ie)
         ap(:,:,k,ie) = laplace_sphere_wk_hyb(p(:,:,k,ie),deriv(hybrid%ithr),elem(ie),viscosity)
       enddo
     enddo
     do ie = nets,nete
       do k = 1, nlev
         ap(:,:,k,ie) = elem(ie)%rspheremp(:,:)*ap(:,:,k,ie)
       enddo
       call edgevpack(edge2,ap(:,:,:,ie),nlev,0,elem(ie)%desc)
     enddo
     call bndry_exchangev(hybrid,edge2)
     do ie = nets,nete
       call edgevunpack(edge2,ap(:,:,:,ie),nlev,0,elem(ie)%desc)
       do k = 1,nlev
         ap(:,:,k,ie) = elem(ie)%spheremp(:,:)*ap(:,:,k,ie)
       enddo
     enddo
     ap = ap/b_norm
! 2. Calculate alpha,new x,new r,and rsnew
     do ie = nets,nete
       global_shared_buf(ie,1) = sum(p(:,:,:,ie)*ap(:,:,:,ie))
     enddo
     call wrap_repro_sum(nvars = 1,comm = hybrid%par%comm)
     if (global_shared_sum(1).eq.0.d0) exit
     alpha = rsold/global_shared_sum(1)
     x = x+alpha*p
     r = r-alpha*ap
     do ie = nets,nete
       global_shared_buf(ie,1) = sum(r(:,:,:,ie)**2.d0)
     enddo
     call wrap_repro_sum(nvars = 1,comm = hybrid%par%comm)
     rsnew = global_shared_sum(1)
! 3. Check rsnew
!       IF (par%ismasterproc) &
!          WRITE(iulog,*) 'iter = ',iter,'rsnew = ',rsnew
     if (rsnew.lt.cgtol_il) then
!          IF (par%ismasterproc) &
!             WRITE(iulog,*) '# of total iteration in InvLapCg = ', iter
       exit
     endif
! 4. Calculate new p and replace rsold
     p = r+(rsnew/rsold)*p
     rsold = rsnew
   enddo
   call freeedgebuffer(edge2)
!=============================================================================
   return
!
   end subroutine invlapcg
!====================================================================
! Spectral inverse module for solving Poisson equation
! By HJ Song
! On 10/08/2015
!====================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine invlapspec(b, vert_lev, b_inv) ! out:solution
!=============================================================================
! A. Declaration
!=============================================================================
! 1. Local modules
!-------------------------------------------------------------------------------
   use ctrltospecmod, only: ctrltospec_inv
   use spectoctrlmod, only: spectoctrl_inv
! 2. Variables
! 2-1. Input variables
!
   integer(int_kind), intent(in   ) :: vert_lev
   real(real_kind), dimension(np, np, vert_lev, nelemd), intent(in   ) :: b
   real(real_kind), dimension(np, np, vert_lev, nelemd), intent(  out) :: b_inv
   real(real_kind), dimension(nups_l, vert_lev) :: b_up, b_inv_up
   real(real_kind), dimension(vert_lev, zwn2) :: b_sp
   integer(int_kind) :: i, j, l, k, m, ie, iproc, ierr
!=============================================================================
! 1. Calculate spectral coefficients
!
   ierr = timingstart('EP2UP_fromInvLapSpec')
   call convertep2up_hyb(vert_lev,b,b_up)
   ierr = timingstop('EP2UP_fromInvLapSpec')
   ierr = timingstart('CtrlToSpec_inv_fromInvLapSpec')
   call ctrltospec_inv(b_up,vert_lev,b_sp)
   ierr = timingstop('CtrlToSpec_inv_fromInvLapSpec')
! 2. Inverse Laplacian in the spectral space (new implementation for efficiency)
   ierr = timingstart('MultiplayInvCoeff_fromInvLapSpec')
   b_sp(:,1) = 0.d0
   do i = 2,zwn2
     l = int(dsqrt(dble(i-1)))
     b_sp(:,i) = b_sp(:,i)*(-kim_rearth**2.d0/(dble(l)*dble(l+1)))
   enddo
   ierr = timingstop('MultiplayInvCoeff_fromInvLapSpec')
! 3. To grid space
   ierr = timingstart('SpecToCtrl_inv_fromInvLapSpec')
   call spectoctrl_inv(b_sp,vert_lev,b_inv_up)
   ierr = timingstop('SpecToCtrl_inv_fromInvLapSpec')
   ierr = timingstart('UP2EP_fromInvLapSpec')
   call convertup2ep_hyb(vert_lev,b_inv_up,b_inv)
   ierr = timingstop('UP2EP_fromInvLapSpec')
!
   end subroutine invlapspec
!===============================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module backgroundmod
!-------------------------------------------------------------------------------
