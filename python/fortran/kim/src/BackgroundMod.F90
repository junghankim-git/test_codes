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
MODULE BackgroundMod

!===============================================================================
! A. Common modules and variables
!===============================================================================
! 1. Common modules
  USE KiapsBase,         ONLY : int_kind  => KIM_int_kind,        &
                                real_kind => KIM_real8_kind,      &
                                iulog     => KIM_iu_log,          &
                                len_path  => KIM_string_len_path, &
                                len_arg   => KIM_string_len_arg 
  USE Dimensions,        ONLY : np, nlev, ne,              &
                                nvar3d, nvar2d,            &
                                nlevar,                    &
                                zwn, zwn2,                 &
                                zwn_nproc,zwn_iproc,       &
                                zwn_a, zwn2_a,             & ! wavenumber of alpha
                                zwn_nproc_a,zwn_iproc_a,   & ! for alpha
                                reg_nproc,reg_iproc,       &
                                psi_nproc,psi_iproc,       &
                                chi_nproc,chi_iproc,       &
                                Tps_nproc,Tps_iproc,       &
                                hum_nproc,hum_iproc,       &
                                neig,                      & ! # of eigenmodes for psi
                                neig_chi,                  & ! # of eigenmodes for chi
                                neig_Tps,                  & ! # of eigenmodes for T,ps
                                neig_q,                    & ! # of eigenmodes for q
                                nelemd, nelem, nets, nete, &
                                nk,                        &
                                ntime,                     & ! # of time bins in 4D 
                                nt_itime                     ! itime of nature
  USE KiapsGrid,         ONLY : elem, deriv,  &
                                nUPs_l, nUPs_g,   &
                                nEPs_l, nEPs_g,   &
                                DOFs, DOFs_EP,    &
                                ConvertEP2UP_hyb, & 
                                ConvertUP2EP_hyb, &
                                ConvertUP2EP_hyb_l
  USE KiapsWave,         ONLY : zwn_blength, &
                                zwn_bstart,  &
                                psi_blength, &
                                psi_bstart,  &
                                chi_blength, &
                                chi_bstart,  &
                                Tps_blength, &
                                Tps_bstart,  &
                                hum_blength, &
                                hum_bstart,  &
                                map_zwn,     &
                                map_reg,     &
                                map_psi,     &
                                map_chi,     &
                                map_Tps,     &
                                map_hum,     &
                                GenMap
  USE KiapsReg,          ONLY : reg_blength, &
                                reg_bstart
  USE Timing,            ONLY : TimingStart, TimingStop
  USE CtrlToSpecMod,     ONLY : CtrlToSpec_opt, &
                                lat_up,         &
                                norm_fac,       &
                                d_loc_h
  USE KiapsParMPI,       ONLY : Par_Reduce
  USE KiapsParallel,     ONLY : par    => KIM_par,    &
                                hybrid => KIM_hybrid, &
                                par_double_precision, &
                                par_character,        &       
                                par_integer,          &       
                                global_shared_buf,    &
                                global_shared_sum,    &
                                BarrierPar,           &
                                par_max,par_min,      &
                                AbortPar,             &
                                KIM_Par
  USE KIO,               ONLY : KIO_t,         &
                                IniIOSystem,   &
                                FinIOSystem,   &
                                OpenFile,      &
                                CloseFile,     &
                                GetDimension,  &
                                GetVariable
  USE GlobalNorms,       ONLY : Wrap_Repro_Sum
  USE PhysicalConstants, ONLY : rvap => KIM_rvap, &
                                rgas => KIM_rair, &
                                ps0 => KIM_p0,    &
                                gravi=> KIM_g,    &
                                DD_PI,            &
                                KIM_REARTH
!  USE BECOutputMod,      ONLY : InitBECOutput,  &
!                                WriteBECOutput, &
!                                FinBECOutput
  USE BECTuOutputMod,    ONLY : InitBECTuOutput,  &
                                WriteBECTuOutput, &
                                FinBECTuOutput
  USE Edge,              ONLY : edgebuffer_t,   &
                                InitEdgeBuffer, &
                                FreeEdgeBuffer, &
                                EdgeVpack,      &
                                EdgeVunpack
  USE Bndry,             ONLY : Bndry_ExchangeV


! 2. Common variables
  IMPLICIT NONE
  PRIVATE
  CHARACTER(LEN=len_path), DIMENSION(:), ALLOCATABLE &
     :: listFile                                       ! Sample file list
  !CHARACTER(LEN=len_path),           &
  !   PUBLIC :: natureFile,           & ! Nature file name
  !             bgFile,               & ! Background file name
  !             becFile                 ! Background error covariance file name
  CHARACTER(LEN=len_path), DIMENSION(:), ALLOCATABLE, &
     PUBLIC :: bgFile                 ! Background file name
  CHARACTER(LEN=len_path),           &
     PUBLIC :: natureFile,           & ! Nature file name
               becFile                 ! Background error covariance file name
  CHARACTER(LEN=len_arg  ),          &
     PUBLIC :: model_name,           & ! Model which generates the background
               bmethod                 ! Method for Background error covariance
  REAL(KIND=real_kind),PARAMETER &
     :: Tr = 270.D0
  INTEGER(KIND=int_kind), PUBLIC &
     :: cv_opt_wind,             &
        cv_opt_hum,              &
        cv_opt_Tu,               &
        cv_opt_psu,              &
        ens_opt_center         
  REAL(KIND=real_kind), PUBLIC &
     :: beta1,                 &
        infl_ens,              &
        distr_infl
  INTEGER(KIND=int_kind), PUBLIC &
     :: levTrans
  
!===============================================================================

!===============================================================================
! B. Public modules and variables
!===============================================================================
! 1. Public modules
  PUBLIC :: SetReadBack
  !PUBLIC :: ReadBackEns
  !PUBLIC :: ReadBackHyb
  !PUBLIC :: ReadBackNMC
  PUBLIC :: ReadBackH4DEV
  PUBLIC :: CalPhib
  PUBLIC :: InvLapCg
  PUBLIC :: TempToM
  PUBLIC :: CalNLBE
  PUBLIC :: CalPsb
  PUBLIC :: CalPsbNlbe
  PUBLIC :: ModelToCtrl
  PUBLIC :: TlmTempToPhi
  !PUBLIC :: TlmTempToM
  PUBLIC :: PsiToRotWind
  PUBLIC :: AdjustP

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
  REAL(KIND=real_kind), DIMENSION(:,:,:,:), ALLOCATABLE, & ! extended to 4D
     PUBLIC :: bg_smpl
  REAL(KIND=real_kind), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, &
     PUBLIC :: bg_uv_rot                                     ! Rotational wind in background field
  REAL(KIND=real_kind), DIMENSION(:,:,:,:,:), ALLOCATABLE, &
     PUBLIC :: bg_md,                                    & ! Background
               pmid
  REAL(KIND=real_kind), DIMENSION(:,:,:,:), ALLOCATABLE, &
     PUBLIC :: nt_org,                                   & ! Nature
               errvar,                                   & ! variance of backgrouond error
               z_org,                                    & ! Geopotential height for each layer
               bg_md_q,                                  & ! Specific humidity for background
               pint                                        ! Interface level pressure
  REAL(KIND=real_kind), DIMENSION(:,:,:), ALLOCATABLE,   &
     PUBLIC :: bcov_sp,                                  & ! Eigenfunction of backErrCor
               bcov_sp_chi,                              &
               bcov_sp_Tps,                              &
               bcov_sp_q,                                &
               phis,                                     &
               zsfc,                                     &
               rl_org                              
  REAL(KIND=real_kind), DIMENSION(:,:,:,:), ALLOCATABLE, &
     PUBLIC :: u10m,                                     & ! 10m zonal wind
               v10m,                                     & ! 10m meridional wind
               t2m,                                      & ! 2m temperature
               q2m,                                      & ! 2m specific humidity
               tsfc,                                     & ! surface temperature
               sfctype
  REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE, &
     PUBLIC :: eigv,                                 &   ! Eigenvalue of backErrCor
               eigv_chi,                             &
               eigv_Tps,                             &
               eigv_q,                               &
               eigv_a,                               &   ! for alpha-variable
               bcov_sp_a,                            &   ! for alpha-variable
               !up
               bg_eig,                               &
               errvar_up
  REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE, &
     PUBLIC :: hyam, hybm,                         &
               hyai, hybi,                         &
               znu , znw,                          &
               b_clim_v, b_ens_v 
  REAL(KIND=real_kind),    &
     PUBLIC :: cgtol_il,   &
               rescalePsi, &
               rescaleChi, &
               rescaleT,   &
               rescaleQ,   &
               rescalePs
  REAL(KIND=real_kind),    &
     PUBLIC :: b_clim,     &
               b_ens
  INTEGER(KIND=int_kind), &
     PUBLIC :: nsmpl,     &
               lwork,     &
               niter_il
  CHARACTER(LEN=len_path), &
     PUBLIC :: list_file,  &
               bgFileList
  LOGICAL,                &
     PUBLIC :: exst_bg,   &
               exst_ens,  &
               exst_bec,  &
               exst_nt,   &
               bottom2top
  LOGICAL,            &
     PUBLIC :: vlocal
  LOGICAL,            &
     PUBLIC :: is_moist_mr
  REAL(KIND=real_kind),  &
     PUBLIC :: vlenScale

  REAL(KIND=real_kind), DIMENSION(:,:,:), ALLOCATABLE, &
     PUBLIC :: Regress_up

  REAL(KIND=real_kind), DIMENSION(:,:,:), ALLOCATABLE, &
     PUBLIC :: pmid_up,                                &
               pmid_ens_mean
  REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE, &
     PUBLIC :: T_ens_mean,                           &
               q_ens_mean
   REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE, &
     PUBLIC :: zsfc_up,                             &
               zsfc_ens

  INTEGER(KIND=int_kind) :: ierr

!===============================================================================

CONTAINS

!===============================================================================
!  SUBROUTINE: SetReadBack
!
!> @brief
!> - Set reading the names of background ensemble files
!>
!> @date 11JUN2013
!> - HJ SONG: Devise and code
!===============================================================================
  SUBROUTINE SetReadBack

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules

  ! 2. Local variables
    TYPE(KIO_t) &
       :: io_t
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie,     &
          ismpl,           &
          err, ierr, wn,   &
          itime
  !=============================================================================

  !=============================================================================
  ! B. Main Body (set vertical coordinates)
  !=============================================================================

!.. Inquiry of file existence
    exst_bg = .FALSE. 
    IF( LEN_TRIM(bgFileList)        /= 0 )                           &
                 INQUIRE(FILE=trim(bgFileList)       ,EXIST=exst_bg)

    exst_ens= .FALSE.
    IF( LEN_TRIM(list_file) /= 0 )                           &
                 INQUIRE(FILE=trim(list_file),EXIST=exst_ens)

    exst_bec= .FALSE.
    IF( LEN_TRIM(becFile)       /= 0 )                           &
                 INQUIRE(FILE=trim(becFile      ),EXIST=exst_bec)

    CALL Mpi_Barrier(par%comm, ierr)

    IF( par%ismasterproc )THEN
       print*, TRIM(model_name), " is selected for Background"   
       print*, TRIM(bmethod), " is selected for Backgroun error estimateion"   
       print*, TRIM(list_file), " is the list file"
       print*, 'nsmpl:',nsmpl
    ENDIF
    
    ALLOCATE(bgFile(ntime))    
    IF (par%ismasterproc) THEN
      OPEN(UNIT=11, FILE=TRIM(bgFileList), FORM='FORMATTED', STATUS='OLD')
         DO itime = 1, ntime
            READ(11,'(A)') bgFile(itime)
         END DO
      CLOSE(11)
    END IF

    CALL IniIOSystem(io_t, par%comm, par%nprocs, par%iproc)

    IF( bmethod == 'h4dev' .OR. bmethod == '3dvar' .OR. bmethod == '4dclv')THEN   !-- For hybrid or nmc

       IF( exst_bg )THEN
!.. Open a netcdf file of background      
          CALL OpenFile(io_t, TRIM(bgFile(nt_itime)))
          IF (par%ismasterproc) print*,                             &
             'Open ',TRIM(bgFile(nt_itime)),' to set vertical coordinates'
       ELSE
          ALLOCATE(listFile(nsmpl))
          OPEN(UNIT=10, FILE=TRIM(list_file), FORM='FORMATTED', STATUS='OLD')
          DO ismpl = 1, nsmpl
             READ(10,'(A)') listFile(ismpl)
          END DO
          CLOSE(10)
!.. Open a netcdf file of one of the samples
          CALL OpenFile(io_t, TRIM(listFile(1)))
          IF (par%ismasterproc) print*,                             &
             'Open ',TRIM(listFile(1)),' to set vertical coordinates'
       ENDIF

    ELSE                                                  !----------- For ens

       ALLOCATE(listFile(nsmpl))
       OPEN(UNIT=10, FILE=TRIM(list_file), FORM='FORMATTED', STATUS='OLD')
       DO ismpl = 1, nsmpl
          READ(10,'(A)') listFile(ismpl)
       END DO
       CLOSE(10)
!.. Open a netcdf file of Ensemblem sample
       CALL OpenFile(io_t, TRIM(listFile(1)))
       IF (par%ismasterproc) print*,                             &
          'Open ',TRIM(listFile(1)),' to set vertical coordinates'
    ENDIF

!.. Read hyam and hybm
    ALLOCATE( hyam(nlev),hybm(nlev) )
    CALL GetVariable(io_t, 'hyam', hyam)
    CALL GetVariable(io_t, 'hybm', hybm)

       bottom2top = .false.
    IF( hybm(1) > hybm(nlev) )THEN
       bottom2top = .true.
    ENDIF

    IF( bottom2top )THEN 
       hyam= hyam(nlev:1:-1)
       hybm= hybm(nlev:1:-1)
    ENDIF

!.. Read hyai and hybi
    ALLOCATE(hyai(nlev+1), hybi(nlev+1))
    CALL GetVariable(io_t, 'hyai', hyai)
    CALL GetVariable(io_t, 'hybi', hybi)
    IF( bottom2top )THEN 
       hyai= hyai(nlev+1:1:-1)
       hybi= hybi(nlev+1:1:-1)
    ENDIF

    IF( model_name == 'kim-sw' )THEN

       ALLOCATE( znu(nlev), znw(nlev+1) )

       DO k= 1, nlev
          znu(k)= hyam(k) + hybm(k)
       ENDDO
       DO k= 1, nlev+1
          znw(k)= hyai(k) + hybi(k)
       ENDDO
!      IF( par%ismasterproc )THEN
!         print*,"Vertical coord (Top to Bottom): hyam, hybm, eta *1000"
!         DO k = 1, nlev
!            print*, k, hyam(k), hybm(k), znu(k)*1000.
!         ENDDO
!      ENDIF
!   ELSE
!      IF( par%ismasterproc )THEN
!         print*,"Vertical coord (Top to Bottom): hyam, hybm, pressure (Pa)"
!         DO k = 1, nlev
!            print*, k, hyam(k), hybm(k), hyam(k)*ps0 + hybm(k)*1000.D2
!         ENDDO
!       ENDIF

    ENDIF  ! model_name

    CALL Mpi_Barrier(par%comm, ierr)

!.. Check zwn_nproc and zwn_iproc
    IF (par%ismasterproc) print*,'Number of wave coefficients',zwn2
    print*,'proc,zwn_nproc,begin,end:',                                 &
                    par%iproc,zwn_nproc,zwn_iproc+1,zwn_iproc+zwn_nproc
    IF (par%ismasterproc) print*,'Number of wave coefficients for alpha',zwn2_a
    print*,'proc,zwn_nproc_a,begin,end:',                                 &
                    par%iproc,zwn_nproc_a,zwn_iproc_a+1,zwn_iproc_a+zwn_nproc_a

    IF (par%ismasterproc) print*,'End of SetReadBack'

    CALL CloseFile(io_t)
    CALL FinIOSystem(io_t)
  !=============================================================================

    RETURN

  END SUBROUTINE SetReadBack

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
  SUBROUTINE ReadBackH4DEV

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules

  ! 2. Local variables
    REAL(KIND=real_kind), DIMENSION(:,:,:,:,:), ALLOCATABLE &
       :: bg_smpl_q,                                        &
          bg_smpl_p,                                        &
          bg_smpl_chi,                                      &
          bg_smpl_pint
    REAL(KIND=real_kind), DIMENSION(:,:,:,:,:), ALLOCATABLE &
       :: & !bg_smpl_mean,                                   &
          bg_smpl_q_mean,                                 &
          bg_smpl_p_mean,                                 &
          bg_smpl_pint_mean
    REAL(KIND=real_kind), DIMENSION(:,:,:,:), ALLOCATABLE &
       :: errvar_md, bg_smpl_ep
    REAL(KIND=real_kind), DIMENSION(:,:,:,:), ALLOCATABLE &
       :: bg_ct
    REAL(KIND=real_kind), DIMENSION(:,:,:), ALLOCATABLE &
       :: bg_sp_smpl, Reg_avg, bg_smpl_mean
    REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE &
       :: bg_sp, bg_md_up
    REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE &
       :: l_Reg_avg,                                &
          l_eigf_psi,                               &
          l_eigf_chi,                               &
          l_eigf_Tps,                               &
          l_eigf_hum                                   
    REAL(KIND=real_kind)  &
       :: psfact, distr
    REAL(KIND=real_kind)  &
       :: tmpmax, tmpmin, &
          maxv, minv
    INTEGER(KIND=int_kind)        &
       :: i, j, k, l, ie, ig, ks, &
          ismpl, wn, vr,          &
          iproc, vertnum,         &
          icnt, ierr,             &
          itime
  !=============================================================================
    IF (par%ismasterproc) print*,'Background error samples for NMC method'
  !=============================================================================

  !=============================================================================
  ! A. Read background
  !=============================================================================

!.. Allocate  variables
    ALLOCATE( bg_md(np,np,nlevar,nelemd,ntime) )
    ALLOCATE( bg_uv_rot(np,np,nlev,2,nelemd,ntime) )
    ALLOCATE( bg_md_q(np,np,nlev,nelemd))
    ALLOCATE( pint(np,np,nlev+1,nelemd) )
    ALLOCATE( pmid(np,np,nlev,nelemd,ntime) )
    ALLOCATE( z_org(np,np,nlev,nelemd) )
    ALLOCATE( zsfc(np,np,nelemd) )
    ALLOCATE( phis(np,np,nelemd) )
    ALLOCATE( rl_org(np,np,nelemd) )
    ALLOCATE( u10m(np,np,nelemd,ntime) )
    ALLOCATE( v10m(np,np,nelemd,ntime) )
    ALLOCATE( t2m(np,np,nelemd,ntime) )
    ALLOCATE( q2m(np,np,nelemd,ntime) )
    ALLOCATE( tsfc(np,np,nelemd,ntime) )
    ALLOCATE( sfctype(np,np,nelemd,ntime) )

    bg_md   = 0.D0
    bg_md_q = 0.D0

    IF (par%ismasterproc) WRITE(iulog,*) 'model_name:',model_name

    IF( model_name == 'kim-old' )THEN       ! Read Background for KIM old version
    IF (par%ismasterproc) WRITE(iulog,*) 'Call ReadBack_KIMold'

    ELSE IF( model_name == 'kim-sw'       &
        .OR. model_name == 'kim-sh' )THEN   ! Read Background for KIM-SW or -SH
    IF (par%ismasterproc) WRITE(iulog,*) 'Call ReadBack_KIM'
       DO itime = 1, ntime
          CALL ReadBack_KIM(itime)
       END DO
    ELSE
       IF( par%ismasterproc )THEN
          STOP "Unrecognized model name"
       ENDIF
    END IF

    IF (par%ismasterproc) WRITE(iulog,*) 'End of Reading Background File'
    IF (par%ismasterproc) WRITE(iulog,*) 'bg_md(1,1,:,1,1)',bg_md(1,1,:,1,1)

  !=============================================================================

  !=============================================================================
  ! B. Read background error covariance or generate it
  !=============================================================================
  ! 1. Allocate
    ALLOCATE( errvar_md(np,np,nlevar,nelemd) )
    ALLOCATE( errvar(np,np,nlevar,nelemd) )
    ALLOCATE( errvar_up(nUPs_l,nlevar) )
    ALLOCATE( eigv(neig,zwn_nproc) )
    ALLOCATE( bcov_sp(nlev,neig,zwn_nproc) )
    ALLOCATE( eigv_chi(neig_chi,zwn_nproc) )
    ALLOCATE( bcov_sp_chi(nlev,neig_chi,zwn_nproc) )
    ALLOCATE( eigv_Tps(neig_Tps,zwn_nproc) )
    ALLOCATE( bcov_sp_Tps(nlev+1,neig_Tps,zwn_nproc) )
    ALLOCATE( eigv_q(neig_q,zwn_nproc) )
    ALLOCATE( bcov_sp_q(nlev,neig_q,zwn_nproc) )

    errvar_md  = 0.D0
    errvar     = 0.D0
    errvar_up  = 0.D0
    eigv       = 0.D0
    eigv_chi   = 0.D0
    eigv_Tps   = 0.D0
    eigv_q     = 0.D0
    bcov_sp    = 0.D0
    bcov_sp_chi= 0.D0
    bcov_sp_Tps= 0.D0
    bcov_sp_q  = 0.D0

    IF (cv_opt_Tu == 1 .OR. cv_opt_Tu == 2 .OR. cv_opt_Tu ==3) THEN
       ALLOCATE( Regress_up(nlev,nUPs_l,2*nlev+1))
       ALLOCATE( l_Reg_avg(reg_nproc) )
       ALLOCATE( l_eigf_psi(psi_nproc), &
                 l_eigf_chi(chi_nproc), &
                 l_eigf_Tps(Tps_nproc), &
                 l_eigf_hum(hum_nproc) )
    END IF

  !=============================================================================
  ! B-1  Read background error covariance from file
    IF( exst_bec ) THEN
  !=============================================================================

!   Read background error covariance
    CALL ReadBackErrCov( l_Reg_avg,  &
                         l_eigf_psi, & 
                         l_eigf_chi, &
                         l_eigf_Tps, &
                         l_eigf_hum, &
                         errvar_md )
        ! Additional Output: errvar, eigv, bcov_sp, Regress, errvar_md

    !up
    CALL ConvertEP2UP_hyb(nlevar,errvar,errvar_up)

!.. Pop the locally distributed variables
    CALL RegressPop(l_Reg_avg, &
                    Regress_up)

    CALL eigfPop(l_eigf_psi,  &
                 l_eigf_chi,  &
                 l_eigf_Tps,  &
                 l_eigf_hum)
                ! bcov_sp,     &
                ! bcov_sp_chi, &
                ! bcov_sp_Tps, &
                ! bcov_sp_hum)

    IF( .not. ALLOCATED(listFile) ) ALLOCATE(listFile(nsmpl*ntime))

    IF( par%ismasterproc )THEN
       OPEN(UNIT=11, FILE=TRIM(list_file), FORM='FORMATTED', STATUS='OLD')
       DO ismpl = 1, nsmpl*ntime
          READ(11,'(A)') listFile(ismpl)
       END DO
       CLOSE(11)
    END IF
    CALL Mpi_Barrier(par%comm, ierr)
    DO ismpl = 1, nsmpl*ntime
       CALL MPI_Bcast(listFile(ismpl), 256, par_character, &
                      par%root, par%comm, ierr)
    END DO

    !up
    !ALLOCATE( bg_smpl(np,np,nlevar,nelemd,ntime,nsmpl) )
    ALLOCATE( bg_smpl(nUPs_l,nlevar,ntime,nsmpl) )
    bg_smpl           = 0.D0

  IF (bmethod == 'h4dev') THEN
    ! up
    !ALLOCATE( bg_smpl_mean(np,np,nlevar,nelemd,ntime) )
    ALLOCATE( bg_smpl_mean(nUPs_l,nlevar,ntime) )


! 1. Read samples (PIO input) and calculate mean field -----------

    CALL ReadSamples( nsmpl,         &  ! In
                      listFile,      &  ! In
                      bg_smpl)          ! Out

    DEALLOCATE( listFile )

! 2. Perturbation matrix ---------------------------------------------

    bg_smpl_mean      = 0.D0
 IF (par%isMasterProc) PRINT *, 'ens_opt_center=', ens_opt_center
 IF (ens_opt_center .EQ. 0) THEN
    DO itime = 1, ntime
    !up
    !DO ie = 1, nelemd
    DO k  = 1, nlevar
    !DO j  = 1, np
    !DO i  = 1, np
    DO i  = 1, nUPs_l
!.. Calculate Ensemble mean of u, v, T, humidity (q or rh), and ps
       DO ismpl = 1, nsmpl
          !bg_smpl_mean(i,j,k,ie,itime) = bg_smpl_mean(i,j,k,ie,itime)  &
          !                             + bg_smpl(i,j,k,ie,itime,ismpl)
          bg_smpl_mean(i,k,itime) = bg_smpl_mean(i,k,itime)  &
                                       + bg_smpl(i,k,itime,ismpl)
       END DO
       !bg_smpl_mean(i,j,k,ie,itime) = bg_smpl_mean(i,j,k,ie,itime) /DBLE(nsmpl)
       bg_smpl_mean(i,k,itime) = bg_smpl_mean(i,k,itime) /DBLE(nsmpl)
    !END DO
    !END DO
    END DO
    END DO
    END DO ! itime
 ELSE IF (ens_opt_center .EQ. 1) THEN
    !up
    ALLOCATE(bg_md_up(nUPs_l,nlevar))
    DO itime = 1, ntime
       CALL ConvertEP2UP_hyb(nlevar,bg_md(:,:,:,:,itime),bg_md_up)
       bg_smpl_mean(:,:,itime) = bg_md_up
    END DO
 ELSE
    IF (par%isMasterProc) PRINT *, 'Wrong option for the center of ensemble samples'
 END IF
   
!.. Calculate perturbations, delta_X (bg_smpl)
    DO itime = 1, ntime
    !DO ie = 1, nelemd
    DO k  = 1, nlevar
    !DO j  = 1, np
    !DO i  = 1, np
    DO i  = 1, nUPs_l
       DO ismpl = 1, nsmpl
          !bg_smpl(i,j,k,ie,itime,ismpl) = bg_smpl(i,j,k,ie,itime,ismpl) &
          !                              - bg_smpl_mean(i,j,k,ie,itime)
          bg_smpl(i,k,itime,ismpl) = bg_smpl(i,k,itime,ismpl) &
                                        - bg_smpl_mean(i,k,itime)
       END DO
    !END DO
    !END DO
    END DO
    END DO
    END DO ! itime

    IF (par%ismasterproc) WRITE(iulog,*)  "inflation factor for ensemble:",infl_ens
    !up
    !IF (par%ismasterproc) WRITE(iulog,*)  "bg_smpl_mean(1,1,:,1,:):",bg_smpl_mean(1,1,:,1,:)
    IF (par%ismasterproc) WRITE(iulog,*)  "bg_smpl_mean(1,:,:):",bg_smpl_mean(1,:,:)

    DEALLOCATE( bg_smpl_mean )

!.. Model space (U,V,T,Q) -> controal varible space
    ALLOCATE(bg_smpl_ep(np,np,nlevar,nelemd))
    ALLOCATE(bg_ct(np,np,nlevar,nelemd))
    bg_ct     = 0.D0

    DO ismpl = 1, nsmpl
    DO itime = 1, ntime
       IF (par%ismasterproc) WRITE(iulog,'(a,i4,a)')               &
                          'To control variable (ismpl = ', ismpl,')'
       !CALL ModelToCtrl(bg_smpl(:,:,:,:,itime,ismpl),     &
       CALL ConvertUP2EP_hyb(nlevar,bg_smpl(:,:,itime,ismpl),bg_smpl_ep) 
       CALL ModelToCtrl(bg_smpl_ep(:,:,:,:),        &
                        bg_md(:,:,:,:,itime),       &
                        bg_uv_rot(:,:,:,:,:,itime), &
                        bg_ct)
       !up
       !bg_smpl(:,:,:,:,itime,ismpl) = bg_ct(:,:,:,:)
       CALL ConvertEP2UP_hyb(nlevar,bg_ct(:,:,:,:),bg_smpl(:,:,itime,ismpl))
    END DO ! itime
    END DO ! ismpl
    DEALLOCATE(bg_ct)
    DEALLOCATE(bg_smpl_ep)

    bg_smpl = bg_smpl/DSQRT(DBLE(nsmpl-1))
    bg_smpl = infl_ens*bg_smpl
  END IF ! bmethod == '3dvar'

    CALL EigAlpha()

    ELSE

    IF (par%ismasterproc) WRITE(iulog,*) "You need to generate climatological B."

    ENDIF !! exst_bec
  !=============================================================================

    CALL Mpi_Barrier(par%comm, ierr)

!   print*,par%iproc,'bg_md',bg_md(1,1,:,1)

  !=============================================================================
  ! C. Rescaling variance of background error
  !=============================================================================

!.. psi
    IF (rescalePsi /= 1.D0) THEN
    DO ie = 1, nelemd
       errvar(:,:,1:nlev*1,ie)                      &
              = errvar(:,:,1:nlev*1,ie) *rescalePsi
    END DO
       IF( errvar_md(1,1,1,1) /= 0.D0 )THEN
       DO ie = 1, nelemd
          errvar_md(:,:,1:nlev*1,ie)                      &
                 = errvar_md(:,:,1:nlev*1,ie) *rescalePsi
       END DO
       END IF
    IF (par%ismasterproc) WRITE(iulog,*) 'Psi Variance was rescaled, ratio= ',rescalePsi
    END IF

!.. Unbalanced chi
    IF (rescaleChi /= 1.D0) THEN
    DO ie = 1, nelemd
       errvar(:,:,nlev+1:nlev*2,ie)                      &
              = errvar(:,:,nlev+1:nlev*2,ie) *rescaleChi
    END DO
       IF( errvar_md(1,1,1,1) /= 0.D0 )THEN
       DO ie = 1, nelemd
          errvar_md(:,:,nlev+1:nlev*2,ie)                      &
                 = errvar_md(:,:,nlev+1:nlev*2,ie) *rescaleChi
       END DO
       END IF
    IF (par%ismasterproc) WRITE(iulog,*) 'Chi Variance was rescaled, ratio= ',rescaleChi
    END IF

!.. Unbalanced T
    IF (rescaleT /= 1.D0) THEN
    DO ie = 1, nelemd
       errvar(:,:,nlev*2+1:nlev*3,ie)                    &
              = errvar(:,:,nlev*2+1:nlev*3,ie) *rescaleT
    END DO
       IF( errvar_md(1,1,1,1) /= 0.D0 )THEN
       DO ie = 1, nelemd
          errvar_md(:,:,nlev*2+1:nlev*3,ie)                    &
                 = errvar_md(:,:,nlev*2+1:nlev*3,ie) *rescaleT
       END DO
       END IF
    IF (par%ismasterproc) WRITE(iulog,*) 'T Variance was rescaled, ratio= ',rescaleT
    END IF

!.. Unbalanced Q
    IF (rescaleQ /= 1.D0) THEN
    DO ie = 1, nelemd
       errvar(:,:,nlev*3+1:nlev*4,ie)                    &
              = errvar(:,:,nlev*3+1:nlev*4,ie) *rescaleQ
    END DO
       IF( errvar_md(1,1,1,1) /= 0.D0 )THEN
       DO ie = 1, nelemd
          errvar_md(:,:,nlev*3+1:nlev*4,ie)                    &
                 = errvar_md(:,:,nlev*3+1:nlev*4,ie) *rescaleQ
       END DO
       END IF
    IF (par%ismasterproc) WRITE(iulog,*) 'Q Variance was rescaled, ratio= ',rescaleQ
    END IF

!.. Unbalanced Ps
    IF (rescalePs /= 1.D0) THEN
    DO ie = 1, nelemd
       errvar(:,:,nlev*4+1,ie)                    &
              = errvar(:,:,nlev*4+1,ie) *rescalePs
    END DO
       IF( errvar_md(1,1,1,1) /= 0.D0 )THEN
       DO ie = 1, nelemd 
          errvar_md(:,:,nlev*4+1,ie)                    &
                 = errvar_md(:,:,nlev*4+1,ie) *rescalePs
       END DO
       END IF
    IF (par%ismasterproc) WRITE(iulog,*) 'Ps Variance was rescaled, ratio= ',rescalePs
    END IF

!.. Calculate vertically-varying weight for climatological and ensemble B
!   levTrans: default = 25
    ALLOCATE(b_ens_v(nlevar))
    ALLOCATE(b_clim_v(nlevar))
    DO k = 1, levTrans
       distr = DBLE(levTrans+1-k)/DBLE(levTrans)*distr_infl
       b_ens_v(k) = (1.D0/EXP(distr**2.D0))*b_ens
    END DO
    DO k = levTrans+1, nlev
       b_ens_v(k) = b_ens
    END DO

    b_ens_v(nlev+1:2*nlev)=b_ens_v(1:nlev)
    b_ens_v(2*nlev+1:3*nlev)=b_ens_v(1:nlev)
    b_ens_v(3*nlev+1:4*nlev)=b_ens_v(1:nlev)
    b_ens_v(4*nlev+1)=b_ens

    b_clim_v = (b_clim+b_ens) - b_ens_v

    IF( par%ismasterproc ) WRITE(iulog,*) 'b_clim:',b_clim
    IF( par%ismasterproc ) WRITE(iulog,*) 'b_ens_v:',b_ens_v
    IF( par%ismasterproc ) WRITE(iulog,*) 'b_clim_v:',b_clim_v

  !=============================================================================
  ! D. Read Nature file (for verification of OSSE test)
  !=============================================================================

    exst_nt = .FALSE.
    IF( LEN_TRIM(natureFile) /= 0 )                            &
                INQUIRE(FILE=trim(natureFile) ,EXIST=exst_nt)

    IF( exst_nt )THEN

       ALLOCATE(nt_org(np,np,nlevar,nelemd))
       nt_org = 0.D0

       CALL ReadNature

    ENDIF

  !=============================================================================
  ! E. Check background error
  !=============================================================================

    IF( exst_nt .and. errvar_md(1,1,1,1) /= 0.D0 )THEN

    global_shared_buf(:,2:11) = 0.D0

    DO ie = 1, nelemd

!   Calculate variance between background and nature (From level 11 to bottom)
       global_shared_buf(ie,2) = global_shared_buf(ie,2)                       &
                               + SUM( (bg_md(:,:,11:nlev,ie,nt_itime)                   &
                                     -nt_org(:,:,11:nlev         ,ie))**2.D0 )
       global_shared_buf(ie,3) = global_shared_buf(ie,3)                       &
                               + SUM( (bg_md(:,:,11+nlev*1:nlev*2,ie,nt_itime)          &
                                     -nt_org(:,:,11+nlev*1:nlev*2,ie))**2.D0 )
       global_shared_buf(ie,4) = global_shared_buf(ie,4)                       &
                               + SUM( (bg_md(:,:,11+nlev*2:nlev*3,ie,nt_itime)          &
                                     -nt_org(:,:,11+nlev*2:nlev*3,ie))**2.D0 )
       global_shared_buf(ie,5) = global_shared_buf(ie,5)                       &
                               + SUM( (bg_md(:,:,11+nlev*3:nlev*4,ie,nt_itime)          &
                                     -nt_org(:,:,11+nlev*3:nlev*4,ie))**2.D0 )
       global_shared_buf(ie,6) = global_shared_buf(ie,6)                       &
                               + SUM( (bg_md(:,:,1+nlev*4,ie,nt_itime)                  &
                                     -nt_org(:,:,1+nlev*4        ,ie))**2.D0 )
!   Calculate Variance of Background Error covariance (From level 11 to bottom)
       global_shared_buf(ie,7) = global_shared_buf(ie,7)                 &
                               + SUM(errvar_md(:,:,11:nlev         ,ie))
       global_shared_buf(ie,8) = global_shared_buf(ie,8)                 &
                               + SUM(errvar_md(:,:,11+nlev*1:nlev*2,ie))
       global_shared_buf(ie,9) = global_shared_buf(ie,9)                 &
                               + SUM(errvar_md(:,:,11+nlev*2:nlev*3,ie))
       global_shared_buf(ie,10)= global_shared_buf(ie,10)                &
                               + SUM(errvar_md(:,:,11+nlev*3:nlev*4,ie))
       global_shared_buf(ie,11)= global_shared_buf(ie,11)                &
                               + SUM(errvar_md(:,:,1+nlev*4        ,ie))

    ENDDO

    CALL Mpi_Barrier(par%comm, ierr)

    CALL Wrap_Repro_Sum(nvars=11, comm=hybrid%par%comm)

    IF (par%ismasterproc) WRITE(iulog,*) 'Variance between Background and Nature'
    IF (par%isMasterProc) WRITE(iulog,*) 'vari_bg_u = ', global_shared_sum(2) /DBLE(np*np*nelem*nlev)
    IF (par%isMasterProc) WRITE(iulog,*) 'vari_bg_v = ', global_shared_sum(3) /DBLE(np*np*nelem*nlev)
    IF (par%isMasterProc) WRITE(iulog,*) 'vari_bg_T = ', global_shared_sum(4) /DBLE(np*np*nelem*nlev)
    IF (par%isMasterProc) WRITE(iulog,*) 'vari_bg_q = ', global_shared_sum(5) /DBLE(np*np*nelem*nlev)
    IF (par%isMasterProc) WRITE(iulog,*) 'vari_bg_ps = ',global_shared_sum(6) /DBLE(np*np*nelem*1 )

    IF (par%ismasterproc) WRITE(iulog,*) 'Variance of Background Error Covariance'
    IF (par%isMasterProc) WRITE(iulog,*) 'vari_B_u  = ',global_shared_sum(7) /DBLE(np*np*nelem*nlev)
    IF (par%isMasterProc) WRITE(iulog,*) 'vari_B_v  = ',global_shared_sum(8) /DBLE(np*np*nelem*nlev)
    IF (par%isMasterProc) WRITE(iulog,*) 'vari_B_T  = ',global_shared_sum(9) /DBLE(np*np*nelem*nlev)
    IF (par%isMasterProc) WRITE(iulog,*) 'vari_B_q  = ',global_shared_sum(10)/DBLE(np*np*nelem*nlev)
    IF (par%isMasterProc) WRITE(iulog,*) 'vari_B_ps = ',global_shared_sum(11)/DBLE(np*np*nelem*1 )

    ENDIF  ! exst_nt

  !=============================================================================

    IF (par%ismasterproc) WRITE(iulog,*) 'exst_nt=',exst_nt

    DEALLOCATE( errvar_md )
    DEALLOCATE( errvar )
    DEALLOCATE(l_Reg_avg)
    DEALLOCATE(l_eigf_psi, &
               l_eigf_chi, &
               l_eigf_Tps, &
               l_eigf_hum) 

    IF (par%isMasterProc) WRITE(iulog,*) 'ReadBackH4DEV End'

    END SUBROUTINE ReadBackH4DEV
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

  SUBROUTINE EigAlpha()

  ! Local variables
    REAL(KIND=real_kind), DIMENSION(nlev+1) &
       :: pmid_avg,                         &
          tmp
    REAL(KIND=real_kind), DIMENSION(nlev+1,nlev+1) &
       :: vloc
    REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE &
       :: work
    REAL(KIND=real_kind)  &
       :: dist,           &
          dist_zero2,     &
          ratio
    INTEGER(KIND=int_kind) &
       :: i, j, ie, l, k,  &
          wn,              &
          info
    INTEGER(KIND=int_kind) &
       :: ierr

    DO k = 1, nlev
       global_shared_buf(:,1) = 0.D0

       DO ie = 1, nelemd
          global_shared_buf(ie,1) = global_shared_buf(ie,1) &
                                  + SUM(pmid(:,:,k,ie,1))
       ENDDO
       CALL Mpi_Barrier(par%comm, ierr)
       CALL Wrap_Repro_Sum(nvars=1, comm=hybrid%par%comm)
       pmid_avg(k)=global_shared_sum(1)/DBLE(np*np*nelem)
    ENDDO

    global_shared_buf(:,1) = 0.D0

    DO ie = 1, nelemd
       global_shared_buf(ie,1) = global_shared_buf(ie,1) &
                               + SUM(bg_md(:,:,4*nlev+1,ie,1))
    ENDDO
    CALL Mpi_Barrier(par%comm, ierr)
    CALL Wrap_Repro_Sum(nvars=1, comm=hybrid%par%comm)
    pmid_avg(nlev+1)=global_shared_sum(1)/DBLE(np*np*nelem)

    IF (par%isMasterProc) WRITE(iulog,*) 'averaged pmid = ', pmid_avg(:)

!   Make localization matrix
    vloc=0.D0

!   vlenScale= 0.2D0 (default)
    dist_zero2=DSQRT(10.D0/3.D0)*vlenScale

    DO l = 1, nlev+1
    DO k = 1, nlev+1

       dist=DABS(DLOG(pmid_avg(k))-DLOG(pmid_avg(l)))
       ratio=dist/dist_zero2

       IF (dist .le. dist_zero2) THEN
          vloc(k,l) = -1.D0/4.D0*(ratio**5.D0) + &
                       1.D0/2.D0*(ratio**4.D0) + &
                       5.D0/8.D0*(ratio**3.D0) - &
                       5.D0/3.D0*(ratio**2.D0) + &
                       1.D0
       ELSE
          IF (dist .le. 2.D0*dist_zero2) THEN
             vloc(k,l) = 1.D0/12.D0*(ratio**5.D0) - &
                         1.D0/2.D0*(ratio**4.D0) +  &
                         5.D0/8.D0*(ratio**3.D0) +  &
                         5.D0/3.D0*(ratio**2.D0) -  &
                         5.D0*ratio +               &
                         4.D0 -                     &
                         2.D0/3.D0/ratio
          ELSE
             vloc(k,l) = 0.D0
          END IF
       END IF

    END DO
    END DO

    IF (par%isMasterProc) THEN
        PRINT *, 'vloc'
        DO k = 1, nlev+1
           WRITE(111,*), (vloc(k,l),l=1,nlev+1)
        END DO
        PRINT *, ' '
    END IF
    CALL Mpi_Barrier(par%comm, ierr)

!   Eigendecomposition and make eigen-pairs for alpha
    ALLOCATE(bcov_sp_a(nlev+1,nlev+1))
    ALLOCATE(eigv_a(nlev+1,zwn_nproc_a))
    ALLOCATE(work(lwork))

    CALL DsyEv('V', 'L', nlev+1, vloc, nlev+1, &
               tmp, work, lwork, info)

    DO k = 1, nlev+1
       bcov_sp_a(:,k) = vloc(:,nlev+1-k+1)
    DO wn = 1, zwn_nproc_a
       eigv_a(k,wn) = tmp(nlev+1-k+1)*                                 &
                      norm_fac*                                        &
                      DEXP(-DINT(DSQRT(DBLE(wn+zwn_iproc_a-1)))**2.D0/ &
                            d_loc_h**2.D0)
    END DO
    END DO

    IF (par%ismasterproc) THEN
    DO k = 1, nlev+1
       WRITE(iulog,*) 'eigNum, eigv_a(eigNum,1):',k,eigv_a(k,1)
    END DO
    END IF
    CALL Mpi_Barrier(par%comm, ierr)

    RETURN

  END SUBROUTINE EigAlpha


!========================================================================================
  SUBROUTINE ZonAvg(errvar) ! inout

    INCLUDE 'mpif.h'

    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd), &
      INTENT(INOUT) :: errvar

    INTEGER(KIND=int_kind), DIMENSION(:), ALLOCATABLE &
      :: n_sum,                                       &
         n_avg
    REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE &
      :: errvar_sum,                                  &
         errvar_avg,                                  &
         errvar_up

    INTEGER(KIND=int_kind) &
      :: i, j, ie,         &
         ilat, k,          &
         ierr, iup

    ALLOCATE(errvar_up(nUPs_l,nlevar))
    ierr = TimingStart('EP2UP_fromZonAvg')
    CALL ConvertEP2UP_hyb(nlevar,errvar,errvar_up)
    ierr = TimingStop('EP2UP_fromZonAvg')

    ALLOCATE(errvar_sum(nlevar,180),n_sum(180))

    errvar_sum = 0.D0
    n_sum = 0

    ! Conduct zonal-averaging of variance
    DO  iup = 1, nUPs_l
       DO ilat = 1, 180
       IF (lat_up(iup,1).GE.DBLE(ilat-91)/180.D0*dd_pi .AND. lat_up(iup,1).LT.DBLE(ilat-90)/180.D0*dd_pi) THEN
          errvar_sum(:,ilat) = errvar_sum(:,ilat) + errvar_up(iup,:)
          n_sum(ilat) = n_sum(ilat) + 1
       END IF
       END DO
       IF (lat_up(iup,1).EQ.DBLE(180-90)/180.D0*dd_pi) THEN
          errvar_sum(:,180) = errvar_sum(:,180) + errvar_up(iup,:)
          n_sum(180) = n_sum(180) + 1
       END IF
    END DO
    CALL Mpi_Barrier(par%comm,ierr)

    ALLOCATE(errvar_avg(nlevar,180),n_avg(180))

    errvar_avg = 0.D0
    n_avg = 0

    CALL Mpi_AllReduce(errvar_sum,           &
                       errvar_avg,           &
                       nlevar*180,           &
                       par_double_precision, &
                       Mpi_Sum,              &
                       par%comm,             &
                       ierr)

    CALL Mpi_AllReduce(n_sum,       &
                       n_avg,       &
                       180,         &
                       par_integer, &
                       Mpi_Sum,     &
                       par%comm,    &
                       ierr)

    DO ilat = 1, 180
       IF (par%isMasterProc) THEN
          print *, 'iLat= ', iLat
          print *, 'errvar_sum(1,ilat) =', errvar_sum(1,ilat)
          print *, 'n_sum(ilat) =', n_sum(ilat)
          print *, 'errvar_avg(1,ilat) =', errvar_avg(1,ilat)
          print *, 'n_avg(ilat) =', n_avg(ilat)
       END IF
       IF (n_avg(ilat).GT.0) THEN
          errvar_avg(:,ilat) = errvar_avg(:,ilat)/DBLE(n_avg(ilat))
       END IF
       IF (par%isMasterProc) THEN
          print *, 'errvar_avg(1,ilat)', errvar_avg(1,ilat)
       END IF
    END DO
    CALL Mpi_Barrier(par%comm, ierr)

    DEALLOCATE(errvar_sum,n_sum)
    DEALLOCATE(n_avg)

    errvar_up = 0.D0
    DO  iup = 1, nUPs_l

       DO ilat = 1, 180
       IF (lat_up(iup,1).GE.DBLE(iLat-91)/180.D0*dd_pi .AND. lat_up(iup,1).LT.DBLE(iLat-90)/180.D0*dd_pi) THEN
          errvar_up(iup,:) = errvar_avg(:,ilat)
       END IF
       END DO
       IF (lat_up(iup,1).EQ.DBLE(180-90)/180.D0*dd_pi) THEN
          errvar_up(iup,:) = errvar_avg(:,180)
       END IF
      
       DO k = 1, nlevar
          IF (errvar_up(iup,k).eq.0.D0) THEN
             PRINT *, 'errvar zero at iup, k =', iup, k
          END IF
       END DO
    END DO
    CALL Mpi_Barrier(par%comm, ierr)

    DEALLOCATE(errvar_avg)

    ierr = TimingStart('UP2EP_fromZonAvg')
    CALL ConvertUP2EP_hyb(nlevar,errvar_up,errvar)
    ierr = TimingStop('UP2EP_fromZonAvg')

    DEALLOCATE(errvar_up)

  END SUBROUTINE ZonAvg
!==================================================================================

!====================================================================================
  SUBROUTINE eigfPop(l_eigf_psi,  & ! in
                     l_eigf_chi,  & ! in 
                     l_eigf_Tps,  & ! in
                     l_eigf_hum)    ! in
                     !bcov_sp,     & ! out
                     !bcov_sp_chi, & ! out
                     !bcov_sp_Tps, & ! out
                     !bcov_sp_q    & ! out

    INCLUDE 'mpif.h'

  ! 1. Input and Output
    REAL(KIND=real_kind), DIMENSION(psi_nproc), &
      INTENT(IN) :: l_eigf_psi
    REAL(KIND=real_kind), DIMENSION(chi_nproc), &
      INTENT(IN) :: l_eigf_chi
    REAL(KIND=real_kind), DIMENSION(Tps_nproc), &
      INTENT(IN) :: l_eigf_Tps
    REAL(KIND=real_kind), DIMENSION(hum_nproc), &
      INTENT(IN) :: l_eigf_hum

  ! 2. Local variables
    REAL(KIND=real_kind), DIMENSION(:,:,:), ALLOCATABLE &
       :: errcov_sp_psi,                                &
          errcov_sp_chi,                                &
          errcov_sp_Tps,                                &
          errcov_sp_hum
    REAL(KIND=real_kind), DIMENSION(:,:,:), ALLOCATABLE &
       :: eigf_swn_psi,                                 &
          eigf_swn_chi,                                 &
          eigf_swn_Tps,                                 &
          eigf_swn_hum 
    REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE &
       :: eigf_psi1d,                               &
          eigf_chi1d,                               &
          eigf_Tps1d,                               &
          eigf_hum1d
    INTEGER(KIND=int_kind) &
       :: i, k, l,         &
          ismpl,           &
          wn,              &
          swn, iswn,       &
          iproc,           &
          icnt,            &
          info, ierr

! -> Global variable
    ALLOCATE(eigf_psi1d(nlev*neig*zwn))
    ALLOCATE(eigf_chi1d(nlev*neig_chi*zwn))
    ALLOCATE(eigf_Tps1d((nlev+1)*neig_Tps*zwn))
    ALLOCATE(eigf_hum1d(nlev*neig_q*zwn))

    CALL Mpi_AllGatherV(l_eigf_psi,           &
                        psi_nproc,            &
                        par_double_precision, &
                        eigf_psi1d,           &
                        psi_blength,          &
                        psi_bstart,           &
                        par_double_precision, &
                        par%comm,             &
                        ierr)
    CALL Mpi_AllGatherV(l_eigf_chi,           &
                        chi_nproc,            &
                        par_double_precision, &
                        eigf_chi1d,           &
                        chi_blength,          &
                        chi_bstart,           &
                        par_double_precision, &
                        par%comm,             &
                        ierr)
    CALL Mpi_AllGatherV(l_eigf_Tps,           &
                        Tps_nproc,            &
                        par_double_precision, &
                        eigf_Tps1d,           &
                        Tps_blength,          &
                        Tps_bstart,           &
                        par_double_precision, &
                        par%comm,             &
                        ierr)
    CALL Mpi_AllGatherV(l_eigf_hum,           &
                        hum_nproc,            &
                        par_double_precision, &
                        eigf_hum1d,           &
                        hum_blength,          &
                        hum_bstart,           &
                        par_double_precision, &
                        par%comm,             &
                        ierr)

! -> Spherical wavenumber structure
    ALLOCATE(eigf_swn_psi(nlev,neig,zwn),       &
             eigf_swn_chi(nlev,neig_chi,zwn),   &
             eigf_swn_Tps(nlev+1,neig_Tps,zwn), &
             eigf_swn_hum(nlev,neig_q,zwn))

    eigf_swn_psi=0.D0
    eigf_swn_chi=0.D0
    eigf_swn_Tps=0.D0
    eigf_swn_hum=0.D0

    i=0
    DO iswn = 1, zwn
    DO l    = 1, neig
    DO k    = 1, nlev
       i = i + 1
       eigf_swn_psi(k,l,iswn)=eigf_psi1d(i)
    END DO
    END DO
    END DO

    i=0
    DO iswn = 1, zwn
    DO l    = 1, neig_chi
    DO k    = 1, nlev
       i = i + 1
       eigf_swn_chi(k,l,iswn)=eigf_chi1d(i)
    END DO
    END DO
    END DO

    i=0
    DO iswn = 1, zwn
    DO l    = 1, neig_Tps
    DO k    = 1, nlev+1
       i = i + 1
       eigf_swn_Tps(k,l,iswn)=eigf_Tps1d(i)
    END DO
    END DO
    END DO

    i=0
    DO iswn = 1, zwn
    DO l    = 1, neig_q
    DO k    = 1, nlev
       i = i + 1
       eigf_swn_hum(k,l,iswn)=eigf_hum1d(i)
    END DO
    END DO
    END DO

    DEALLOCATE(eigf_psi1d)
    DEALLOCATE(eigf_chi1d)
    DEALLOCATE(eigf_Tps1d)
    DEALLOCATE(eigf_hum1d)

! -> Distributed wavenumber space
    DO wn = 1, zwn_nproc
       swn = DINT(DSQRT(DBLE(wn+zwn_iproc-1))) ! spherical wavenumber

       bcov_sp(:,:,wn)    =eigf_swn_psi(:,:,swn+1)
       bcov_sp_chi(:,:,wn)=eigf_swn_chi(:,:,swn+1)
       bcov_sp_Tps(:,:,wn)=eigf_swn_Tps(:,:,swn+1)
       bcov_sp_q(:,:,wn)  =eigf_swn_hum(:,:,swn+1)
    END DO

    DEALLOCATE(eigf_swn_psi)
    DEALLOCATE(eigf_swn_chi)
    DEALLOCATE(eigf_swn_Tps)
    DEALLOCATE(eigf_swn_hum)

    RETURN
  END SUBROUTINE eigfPop
!========================================================================================

!========================================================================================
  SUBROUTINE RegressLoc(Reg_avg, &
                        l_Reg_avg)

    REAL(KIND=real_kind), DIMENSION(nlev,2*nlev+1,180), &
      INTENT(IN) :: Reg_avg

    REAL(KIND=real_kind), DIMENSION(reg_nproc), &
      INTENT(OUT) :: l_Reg_avg

    REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE &
      :: Reg_tmp

    INTEGER(KIND=int_kind) &
      :: i,                &
         ilat, ireg, k

! -> local variable
    ALLOCATE(Reg_tmp(nlev*(2*nlev+1)*180))

    i = 0
    DO ilat = 1, 180
    DO ireg = 1, 2*nlev+1
    DO k    = 1, nlev
       i = i + 1
       Reg_tmp(i)=Reg_avg(k,ireg,ilat)
    END DO
    END DO
    END DO

    l_Reg_avg(1:reg_nproc) = Reg_tmp(reg_iproc+1:reg_iproc+reg_nproc)

    DEALLOCATE(Reg_tmp)


  END SUBROUTINE RegressLoc
!==================================================================================

!==================================================================================
  SUBROUTINE RegressPop(l_Reg_avg, &
                        Regress_up)

    INCLUDE 'mpif.h'

    REAL(KIND=real_kind), DIMENSION(reg_nproc), &
      INTENT(IN) :: l_Reg_avg

    REAL(KIND=real_kind), DIMENSION(nlev,nUPs_l,nk), &
      INTENT(OUT) :: Regress_up

    REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE &
      :: Reg_tmp
    REAL(KIND=real_kind), DIMENSION(:,:,:), ALLOCATABLE &
      :: lat,                                           &
         Reg_avg
    INTEGER(KIND=int_kind) &
      :: i, j, ie,         &
         ilat, ireg, k,    &
         ierr, iup

! -> global variable
    ALLOCATE(Reg_tmp(nlev*nk*180))

    CALL Mpi_AllGatherV(l_Reg_avg,            &
                        reg_nproc,            &
                        par_double_precision, &
                        Reg_tmp,              &
                        reg_blength,          &
                        reg_bstart,           &
                        par_double_precision, &
                        par%comm,             &
                        ierr)

    ALLOCATE(Reg_avg(nlev,nk,180))

    i = 0
    DO ilat = 1, 180
    DO ireg = 1, nk
    DO k    = 1, nlev
       i = i + 1
       Reg_avg(k,ireg,ilat) = Reg_tmp(i)
    END DO
    END DO
    END DO

    DEALLOCATE(Reg_tmp)

    Regress_up = 0.D0
    DO  iup = 1, nUPs_l
       DO ilat = 1, 180
       IF (lat_up(iup,1).GE.DBLE(ilat-91)/180.D0*dd_pi .AND. lat_up(iup,1).LT.DBLE(ilat-90)/180.D0*dd_pi) THEN
          Regress_up(:,iup,:) = Reg_avg(:,:,ilat)
       END IF
       END DO
       IF (lat_up(iup,1).EQ.DBLE(180-90)/180.D0*dd_pi) THEN
          Regress_up(:,iup,:) = Reg_avg(:,:,180)
       END IF
    END DO
    CALL Mpi_Barrier(par%comm, ierr)

  END SUBROUTINE RegressPop
!==========================================================================================

!==========================================================================================
  SUBROUTINE RedDof_miss(x, &
                         x_reduced)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Variables
  ! 1-1. Input
    REAL(KIND=real_kind), DIMENSION(np,np,nelemd), &
       INTENT(IN) :: x

  ! 1-2. Output
    REAL(KIND=real_kind), DIMENSION(np,np,nelemd), &
       INTENT(OUT) :: x_reduced

  ! 1-3. Local
    INTEGER(KIND=int_kind)    &
       :: l_ie, ie,           &
          ist, ien, jst, jen, &
          i, j
  !=============================================================================

  !=============================================================================
  ! B. Main body
  !=============================================================================
  ! 1. Zero out duplicated points
    x_reduced = -999.D0
    DO l_ie = nets, nete
       ie = elem(l_ie)%GlobalId

       jst = 1
       jen = np-1
       ist = 1
       ien = np-1
       IF (((ne*ne*2 + 1) .LE. ie) .AND. (ie .LE. (ne*ne*4))) THEN
          jst = 2
          jen = np
       END IF
       IF ((ne*ne*5 + 1) .LE. ie) THEN
          ist = 2
          ien = np
       END IF

       x_reduced(ist:ien,jst:jen,l_ie) = x(ist:ien,jst:jen,l_ie)
       IF (ie .EQ. (ne*(ne-1)+1)) x_reduced(1,np,l_ie) = x(1,np,l_ie)
       IF (ie .EQ. ne*(ne+1)) x_reduced(np,1,l_ie) = x(np,1,l_ie)
    END DO
  !=============================================================================

  END SUBROUTINE RedDof_miss
!=======================================================================================

!===============================================================================
  SUBROUTINE ReadSamples( nsmpl,         &
                          infilelist,    &
                          bg_smpl)
!===============================================================================

   USE NETCDF
 
   USE VertIntMod,      ONLY : VertInt 
   USE ObsCorrMod,      ONLY : SfcPrsCorr                                 

   INTEGER(KIND=int_kind), &
        INTENT(IN)   :: nsmpl
   CHARACTER(LEN=len_path),DIMENSION(nsmpl*ntime), &
        INTENT(IN)   :: infilelist
   !up
   !REAL(KIND=real_kind),DIMENSION(np,np,nlevar,nelemd,ntime,nsmpl), &
   REAL(KIND=real_kind),DIMENSION(nUPs_l,nlevar,ntime,nsmpl), &
        INTENT(OUT)  :: bg_smpl

! Local variables

   !up
   !REAL(KIND=real_kind), DIMENSION(:,:,:), ALLOCATABLE &
   REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE &
       :: ps_tmp
   !REAL(KIND=real_kind), DIMENSION(:,:,:,:), ALLOCATABLE &
   REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE &
       :: pmid_ens, bg_smpl_tmp
   REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE &
       :: x3d, x3di, x2d
   REAL(KIND=real_kind), DIMENSION(:,:,:,:), ALLOCATABLE &
       :: zsfc_ep
!   REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE &
!       :: zsfc_up, zsfc_ens
   REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE &
!       :: pmid_up, bg_md_up
       :: bg_md_up
   REAL(KIND=real_kind), DIMENSION(:,:,:,:), ALLOCATABLE &
       :: tq_smpl_tmp

   REAL(KIND=real_kind)   :: psfact
   INTEGER(KIND=int_kind) :: i, j, k, ie, ic, icnt, &
                             itime
   INTEGER(KIND=int_kind) :: ierr, tmp
   INTEGER(KIND=int_kind) :: vr
   INTEGER(KIND=int_kind) :: vertnum
   INTEGER(KIND=int_kind) :: status,ncid,varid
   LOGICAL                :: kimswSamples
   TYPE(KIO_t)            :: io_t

   REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE :: UPs_1D, EPs_1D
   REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE :: UPs_2D, EPs_2D

  !=============================================================================
  ! Read Samples
  !=============================================================================
!    ALLOCATE( UPs_1D(np*np*nelemd) )
!    ALLOCATE( UPs_2D(np*np*nelemd,nlev) )
    ALLOCATE( UPs_1D(nUPs_l,1) )
    ALLOCATE( UPs_2D(nUPs_l,nlev) )
    !up
    !ALLOCATE( EPs_1D(nEPs_l,1) )
    !ALLOCATE( EPs_2D(nEPs_l,nlev) )

!    CALL IniInputAPIs

!.. Check for existance of model variable's variance
    status = nf90_open( trim(infilelist(1)),NF90_NOWRITE,ncid )
    status = nf90_inq_varid( ncid,'pint',varid )
    IF( status == 0 )THEN
       IF( par%ismasterproc ) WRITE(iulog,*)                      &
       trim(becFile),"Samples from KIM-SW"
       kimswSamples = .TRUE.
    ELSE
       kimswSamples = .FALSE.
    END IF
    status = nf90_close(ncid)

!    ALLOCATE( x3d(np*np*nelemd*nlev) )
!    ALLOCATE( x3di(np*np*nelemd*(nlev+1)) )
!    ALLOCATE( x2d(np*np*nelemd) )

    bg_smpl        = 0.D0

    !up
    !ALLOCATE( ps_tmp(np,np,nelemd) )
    !ALLOCATE( pmid_ens(np,np,nlev,nelemd) )
    !ALLOCATE( bg_smpl_tmp(np,np,nlev,nelemd) )
    ALLOCATE( ps_tmp(nUPs_l) )
    ALLOCATE( pmid_ens(nUPs_l,nlev) )
    ALLOCATE( bg_smpl_tmp(nUPs_l,nlev) )
    ALLOCATE( tq_smpl_tmp(nUPs_l,2,ntime,nsmpl) )

    ALLOCATE( zsfc_ens(nUPs_l) )
    ALLOCATE(zsfc_ep(np,np,1,nelemd))
    ALLOCATE(zsfc_up(nUPs_l))
    zsfc_ep(:,:,1,:)=zsfc
    CALL ConvertEP2UP_hyb(1,zsfc_ep,zsfc_up)

!    ALLOCATE(pmid_up(nUPs_l,nlev))
    ALLOCATE(pmid_up(nUPs_l,nlev,ntime))
    ALLOCATE(bg_md_up(nUPs_l,4*nlev+1))

    ALLOCATE(pmid_ens_mean(nUPs_l,nlev,ntime))
    ALLOCATE(T_ens_mean(nUPs_l,ntime))
    ALLOCATE(q_ens_mean(nUPs_l,ntime))

    pmid_ens_mean = 0.D0
    T_ens_mean = 0.D0
    q_ens_mean = 0.D0

    DO itime = 1, ntime
!       CALL ConvertEP2UP_hyb(nlev,pmid(:,:,:,:,itime),pmid_up)
       CALL ConvertEP2UP_hyb(nlev,pmid(:,:,:,:,itime),pmid_up(:,:,itime))
       CALL ConvertEP2UP_hyb(4*nlev+1,bg_md(:,:,:,:,itime),bg_md_up)
    DO ic = 1, nsmpl
!..    Open a pio input file
       CALL IniIOSystem(io_t, par%comm, par%nprocs, par%iproc)
       CALL OpenFile(io_t, TRIM(infilelist(ic+nsmpl*(itime-1))))
       CALL GetDimension(io_t, 'ncol', tmp)
       !IF (tmp .LE. 0 .OR. tmp .EQ. nUPs_g) THEN
       IF (tmp .LE. 0) THEN
         IF (par%ismasterproc) print*,'Check dimension(ncol) or EP file: ',TRIM(infilelist(ic+nsmpl*(itime-1)))
         CALL AbortPar()
       END IF
       !CALL GetDimension(io_t, 'ncol', tmp, DOFs_EP)
       CALL GetDimension(io_t, 'ncol', tmp, DOFs)

       IF (par%ismasterproc) print*,'Read samples from ',TRIM(infilelist(ic+nsmpl*(itime-1)))

       !CALL Read2DDoubleInput('p', nlev, l_2dDoubleValues, 1, ierr)
       CALL GetVariable(io_t, 'p', UPs_2D)
       !up
       !CALL ConvertUP2EP_hyb_l(nlev, UPs_2D, EPs_2D)
       IF (par%ismasterproc) WRITE(iulog,*) 'Read pressure at middle level'
       icnt = 0
       !DO ie= 1, nelemd
       !DO j = 1, np
       !DO i = 1, np
       DO i = 1, nUPs_l
          icnt = icnt + 1
          DO k = 1, nlev
             !up
             !pmid_ens(i,j,nlev+1-k,ie) = EPs_2D(icnt,k)
             pmid_ens(i,nlev+1-k) = UPs_2D(icnt,k)
             !pmid(i,j,nlev+1-k,ie,itime) = x3d(icnt)
          END DO
       !END DO
       !END DO
       END DO

       pmid_ens_mean(:,:,itime)=pmid_ens_mean(:,:,itime)+pmid_ens/DBLE(nsmpl)

       DO vr = 1, nvar3d
          SELECT CASE (vr)
          CASE (1)
!..    Read U wind (m/s) and calculate the mean
             CALL GetVariable(io_t, 'u', UPs_2D)
          CASE (2)
!..    Read V wind (m/s) and calculate the mean
             CALL GetVariable(io_t, 'v', UPs_2D)
          CASE (3)
!..    Read Temperature (K) and calculate the mean
             CALL GetVariable(io_t, 'T', UPs_2D)
          CASE (4)
!..    Read Specific humidity (kg/kg) and calculate the mean
             CALL GetVariable(io_t, 'q', UPs_2D)
             IF (is_moist_mr) THEN
                DO k = 1, nlev
                DO i = 1, nUPs_l
                   UPs_2D(i,k) = UPs_2D(i,k)/(1.D0+UPs_2D(i,k))  ! to specific humidity
                END DO
                END DO
             END IF
          END SELECT
          !up
          !CALL ConvertUP2EP_hyb_l(nlev, UPs_2D, EPs_2D)
          IF( bottom2top )THEN
             icnt = 0
             !up
             !DO ie = 1, nelemd
             !DO j  = 1, np
             !DO i  = 1, np
             DO i  = 1, nUPs_l
                icnt = icnt + 1
                DO k  = 1, nlev
                   !bg_smpl(i,j,nlev+1-k+nlev*(vr-1),ie,itime,ic) = &
                   !                          EPs_2D(icnt,k)
                   !up
                   !bg_smpl_tmp(i,j,nlev+1-k,ie) = &
                   !                          EPs_2D(icnt,k)
                   bg_smpl_tmp(i,nlev+1-k) = &
                                             UPs_2D(icnt,k)
                !bg_smpl(i,j,nlev+1-k+nlev*(vr-1),ie,itime,ic) = x3d(icnt)
                END DO
             !END DO
             !END DO
             END DO
          ELSE
             icnt = 0
             !up
             !DO ie = 1, nelemd
             !DO j  = 1, np
             !DO i  = 1, np
             DO i  = 1, nUPs_l
                icnt = icnt + 1
                DO k  = 1, nlev
                  !bg_smpl(i,j,k+nlev*(vr-1),ie,itime,ic) = EPs_2D(icnt,k)
                  !up
                  !bg_smpl_tmp(i,j,k,ie) = EPs_2D(icnt,k)
                  bg_smpl_tmp(i,k) = UPs_2D(icnt,k)
                !bg_smpl(i,j,k+nlev*(vr-1),ie,itime,ic) = x3d(icnt)
                END DO
             !END DO
             !END DO
             END DO
          ENDIF

!--------------------
! Adjust pressure of LETKF samples to H-4DEV background
          !up
          !DO ie = 1, nelemd
          !DO j  = 1, np
          !DO i  = 1, np
          DO i  = 1, nUPs_l
          DO k  = 1, nlev
             !CALL VertInt(p_bgh_v, p_ob, nlev, obmiss, &
             !             bgh,                         &
             !             bgv)
             !up
             !CALL VertInt(pmid_ens(i,j,:,ie),                &
             !             pmid(i,j,k,ie,itime),              &
             CALL VertInt(pmid_ens(i,:),                &
!                          pmid_up(i,k),                 &
                          pmid_up(i,k,itime),           &
                          nlev,                         &
                          !bg_md(i,j,k+nlev*(vr-1),ie,itime), &
                          !up
                          !bg_smpl_tmp(i,j,k,ie),             &
                          !bg_smpl_tmp(i,j,:,ie),             &
                          !bg_smpl(i,j,k+nlev*(vr-1),ie,itime,ic))
                          bg_smpl_tmp(i,k),             &
                          bg_smpl_tmp(i,:),             &
                          bg_smpl(i,k+nlev*(vr-1),itime,ic))
          END DO
          END DO
          !END DO
          !END DO

          IF (vr.EQ.3) THEN
             tq_smpl_tmp(:,1,itime,ic)=bg_smpl_tmp(:,nlev)
          ELSE IF (vr.EQ.4) THEN
             tq_smpl_tmp(:,2,itime,ic)=bg_smpl_tmp(:,nlev)
          END IF

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

       END DO   ! vr

       IF (par%ismasterproc) print*,'Read u, v, T, q'

       DO vr = 1, nvar2d
          SELECT CASE (vr)
          CASE (1)
!..    Read surface pressure (hPa -> Pa) and calculate the mean
          !CALL GetVariable(io_t, 'ps', UPs_1D)
          CALL GetVariable(io_t, 'ps', UPs_1D(:,1))
          !up
          !CALL ConvertUP2EP_hyb_l(1, UPs_1D, EPs_1D)
          psfact=1.D0
          IF (par%ismasterproc) print*,'Read ps, psfact=',psfact
          END SELECT

          icnt = 0
          !up
          !DO ie = 1, nelemd
          !DO j  = 1, np
          !DO i  = 1, np
          DO i  = 1, nUPs_l
             icnt = icnt + 1
             !bg_smpl(i,j,nlev*nvar3d+vr,ie,itime,ic) = UPs_1D(icnt)*psfact
             !ps_tmp(i,j,ie) = UPs_1D(icnt)*psfact
             !up
             !ps_tmp(i,j,ie) = EPs_1D(icnt,1)*psfact
             ps_tmp(i) = UPs_1D(icnt,1)*psfact
             !bg_smpl(i,j,nlev*nvar3d+vr,ie,itime,ic) = x2d(icnt)*psfact
          !END DO
          !END DO
          END DO
       END DO

!------------------
! For ps adjustment according to topography difference
!-----------------
! Read topography of an ensemble sample
       IF (itime.EQ.1 .AND. ic.EQ.1)  THEN
          IF (par%isMasterProc) PRINT *, 'Read topo'
          !CALL GetVariable(io_t, 'topo', UPs_1D)
          CALL GetVariable(io_t, 'topo', UPs_1D(:,1))
          !up
          !CALL ConvertUP2EP_hyb_l(1, UPs_1D, EPs_1D)

          icnt = 0
          !DO ie = 1, nelemd
          !DO j  = 1, np
          !DO i  = 1, np
          DO i  = 1, nUPs_l
             icnt = icnt + 1
             !zsfc_ens(i,j,ie) = EPs_1D(icnt,1)
             zsfc_ens(i) = UPs_1D(icnt,1)
          !END DO
          !END DO
          END DO
       END IF

! Adjusting ps
       !up
       !DO ie = 1, nelemd
       !DO j  = 1, np
       !DO i  = 1, np
       DO i  = 1, nUPs_l
         !CALL SfcPrsCorr (bg_md(i,j,3*nlev,ie,itime),      &
         !                 bg_md(i,j,4*nlev,ie,itime),      &
         !                 zsfc(i,j,ie),                    &
         !                 bg_smpl(i,j,3*nlev,ie,itime,ic), &
         !                 bg_smpl(i,j,4*nlev,ie,itime,ic), &
         !                 zsfc_ens(i,j,ie),                &
         !                 ps_tmp(i,j,ie),                  &
         !                 bg_smpl(i,j,4*nlev+1,ie,itime,ic))
         CALL SfcPrsCorr (bg_md_up(i,3*nlev),        &
                          bg_md_up(i,4*nlev),        &
                          zsfc_up(i),                &
                          !bg_smpl(i,3*nlev,itime,ic), &
                          !bg_smpl(i,4*nlev,itime,ic), &
                          tq_smpl_tmp(i,1,itime,ic), &
                          tq_smpl_tmp(i,2,itime,ic), &
                          zsfc_ens(i),               &
                          ps_tmp(i),                 &
                          bg_smpl(i,4*nlev+1,itime,ic))
       !END DO
       !END DO
       END DO

       T_ens_mean(:,itime)=T_ens_mean(:,itime)+tq_smpl_tmp(:,1,itime,ic)/DBLE(nsmpl)
       q_ens_mean(:,itime)=q_ens_mean(:,itime)+tq_smpl_tmp(:,2,itime,ic)/DBLE(nsmpl)

! Checking adjusted ps
       IF (itime.EQ.1 .AND. ic.EQ.1)  THEN
          IF (par%isMasterProc) PRINT *, 'zsfc', 'zsfc_ens', &
                                         'bg_md_ps', 'ps_ens', 'bg_smpl_ps'
          !up
          !DO ie = 1, nelemd
          !DO j  = 1, np
          !DO i  = 1, np
          DO i  = 1, nUPs_l
             IF (DABS(zsfc_up(i)-zsfc_ens(i)).GT.500.D0) THEN
                PRINT *, zsfc_up(i), zsfc_ens(i), &
                         bg_md_up(i,4*nlev+1), ps_tmp(i), bg_smpl(i,4*nlev+1,itime,ic)
             END IF
          !END DO
          !END DO
          END DO
       END IF
!--------------------------------------------------------------


!..    Close the pio file
       CALL CloseFile(io_t)
       CALL FinIOSystem(io_t)

    END DO    !! ic
    END DO    !! itime

    DEALLOCATE( UPs_1D )
    DEALLOCATE( UPs_2D )
    DEALLOCATE( zsfc_ep )
!    DEALLOCATE( zsfc_up )
!    DEALLOCATE( zsfc_ens )
!    DEALLOCATE( pmid_up )
    DEALLOCATE( bg_smpl_tmp )
    DEALLOCATE( tq_smpl_tmp )
    DEALLOCATE( bg_md_up )
    !DEALLOCATE( EPs_1D )
    !DEALLOCATE( EPs_2D )
    !CALL FinInputAPIs

  END SUBROUTINE ReadSamples

!===============================================================================

!===============================================================================
  SUBROUTINE AdjustP(an_ep         &
                     )
!===============================================================================

   USE VertIntMod,      ONLY : VertInt
   USE ObsCorrMod,      ONLY : SfcPrsCorr

   REAL(KIND=real_kind),DIMENSION(np,np,nlevar,nelemd,ntime), &
        INTENT(INOUT)  :: an_ep

! Local variables
   REAL(KIND=real_kind),DIMENSION(:,:), ALLOCATABLE &
        :: an_up,                                   &
           an_up_adj
   INTEGER(KIND=int_kind) :: i, k, vr, &
                             itime

   ALLOCATE(an_up(nUPs_l,4*nlev+1))
   ALLOCATE(an_up_adj(nUPs_l,4*nlev+1))

    DO itime = 1, ntime
       CALL ConvertEP2UP_hyb(4*nlev+1,an_ep(:,:,:,:,itime),an_up)

       DO vr = 1, nvar3d
          DO i  = 1, nUPs_l
          DO k  = 1, nlev
             CALL VertInt(                                         &
                          pmid_up(i,:,itime),                      &
                          pmid_ens_mean(i,k,itime),                &
                          nlev,                                    &
                          an_up(i,k+nlev*(vr-1)),                  &
                          an_up(i,1+nlev*(vr-1):nlev+nlev*(vr-1)), &
                          an_up_adj(i,k+nlev*(vr-1)))
          END DO
          END DO
       END DO

       DO i  = 1, nUPs_l
         CALL SfcPrsCorr (T_ens_mean(i,itime),   &
                          q_ens_mean(i,itime),   &
                          zsfc_ens(i),           &
                          an_up(i,3*nlev),       &
                          an_up(i,4*nlev),       &
                          zsfc_up(i),            &
                          an_up(i,4*nlev+1),     &
                          an_up_adj(i,4*nlev+1))
       END DO

       CALL ConvertUP2EP_hyb(4*nlev+1,an_up_adj,an_ep(:,:,:,:,itime))
    END DO

    DEALLOCATE(pmid_up, pmid_ens_mean)
    DEALLOCATE(an_up, an_up_adj)
    DEALLOCATE(T_ens_mean, q_ens_mean)
    DEALLOCATE(zsfc_ens, zsfc_up)

    RETURN

  END SUBROUTINE AdjustP


!===============================================================================
  SUBROUTINE ReadBack_KIM(itime)
!===============================================================================


! Local variables

   INTEGER(KIND=int_kind), INTENT(IN) :: itime
   REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE &
       :: x3d, x3di, x2d

   TYPE(KIO_t)            :: io_t
   REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE :: UPs_1D, EPs_1D
   REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE :: UPs_2D, UPs_2D_p1, &
                                                        EPs_2D, EPs_2D_p1

   REAL(KIND=real_kind), DIMENSION(:,:,:,:), ALLOCATABLE &
       :: bg_ct

   REAL(KIND=real_kind)   :: psfact
   INTEGER(KIND=int_kind) :: i, j, k, ie, icnt
   INTEGER(KIND=int_kind) :: ierr, tmp
   INTEGER(KIND=int_kind) :: vr
   INTEGER(KIND=int_kind) :: vertnum

  !=============================================================================
  ! Prepare ReadPio
  !=============================================================================
  !=============================================================================
  ! Read Background
  !=============================================================================

    !ALLOCATE( UPs_1D(np*np*nelemd) )
    !ALLOCATE( UPs_2D(np*np*nelemd,nlev) )
    !ALLOCATE( UPs_2D_p1(np*np*nelemd,nlev+1) )
    ALLOCATE( UPs_1D(nUPs_l,1) )
    ALLOCATE( UPs_2D(nUPs_l,nlev) )
    ALLOCATE( UPs_2D_p1(nUPs_l,nlev+1) )
    ALLOCATE( EPs_1D(nEPs_l,1) )
    ALLOCATE( EPs_2D(nEPs_l,nlev) )
    ALLOCATE( EPs_2D_p1(nEPs_l,nlev+1) )

    CALL IniIOSystem(io_t, par%comm, par%nprocs, par%iproc)
    CALL OpenFile(io_t, TRIM(bgFile(itime)))
    CALL GetDimension(io_t, 'ncol', tmp)
    !IF (tmp .LE. 0 .OR. tmp .EQ. nUPs_g) THEN
    IF (tmp .LE. 0) THEN
      IF (par%ismasterproc) print*,'Check dimension or EP file: ',TRIM(bgFile(itime))
      CALL AbortPar()
    END IF
    CALL GetDimension(io_t, 'ncol', tmp, DOFs)

    IF (par%ismasterproc) WRITE(iulog,*) 'nlev = ', nlev
    IF (par%ismasterproc) WRITE(iulog,*) 'Bacground from KIM model'
    IF (par%ismasterproc) WRITE(iulog,*) 'File of Background=', bgFile(itime)

    Do vr = 1, nvar3d
       SELECT CASE (vr)
       CASE (1)
          CALL GetVariable(io_t, 'u', UPs_2D)
          IF (par%ismasterproc) WRITE(iulog,*) 'Read u'
       CASE (2)
          CALL GetVariable(io_t, 'v', UPs_2D)
          IF (par%ismasterproc) WRITE(iulog,*) 'Read v'
       CASE (3)
          CALL GetVariable(io_t, 'T', UPs_2D)
          IF (par%ismasterproc) WRITE(iulog,*) 'Read T'
       CASE (4)
          CALL GetVariable(io_t, 'q', UPs_2D)
          IF (par%ismasterproc) WRITE(iulog,*) 'Read q'
          IF (is_moist_mr) THEN
             DO k = 1, nlev
             DO i = 1, nUPs_l
                UPs_2D(i,k) = UPs_2D(i,k)/(1.D0+UPs_2D(i,k))  ! to specific humidity
             END DO
             END DO
          END IF
       END SELECT
       CALL ConvertUP2EP_hyb_l(nlev, UPs_2D, EPs_2D)
       DO k = 1, nlev
       icnt = 0
       DO ie= 1, nelemd
       DO j = 1, np
       DO i = 1, np
          icnt = icnt + 1
          bg_md(i,j,nlev+1-k+nlev*(vr-1),ie,itime) = EPs_2D(icnt,k)
          !bg_md(i,j,nlev+1-k+nlev*(vr-1),ie,itime) = UPs_2D(icnt,k)
          !bg_md(i,j,nlev+1-k+nlev*(vr-1),ie,itime) = x3d(icnt,k)
       END DO
       END DO
       END DO
       END DO
!       IF (par%ismasterproc) WRITE(iulog,*) 'vr    = ', vr
!       IF (par%ismasterproc) WRITE(iulog,*) 'input = ', UPs_2D(1,:)
!       IF (par%ismasterproc) WRITE(iulog,*) 'save  = ', bg_md(1,1,nlev*(vr-1)+1:nlev*(vr-1)+nlev,1,itime)
    END DO
    
    bg_md_q(:,:,:,:) = bg_md(:,:,nlev*3+1:nlev*4,:,itime)

    IF(cv_opt_hum == 1) THEN
       CALL GetVariable(io_t, 'rh', UPs_2D)
       CALL ConvertUP2EP_hyb_l(nlev, UPs_2D, EPs_2D)
       IF (par%ismasterproc) WRITE(iulog,*) 'Read rh'
       vr = 4
       DO k = 1, nlev
       icnt = 0
       DO ie= 1, nelemd
       DO j = 1, np
       DO i = 1, np
         icnt = icnt + 1
         bg_md(i,j,nlev+1-k+nlev*(vr-1),ie,itime) = EPs_2D(icnt,k)
         !bg_md(i,j,nlev+1-k+nlev*(vr-1),ie,itime) = x3d(icnt)
       END DO
       END DO
       END DO
       END DO
    END IF

    DO vr = 1, nvar2d
       SELECT CASE (vr)
       CASE (1)
          CALL GetVariable(io_t, 'ps', UPs_1D(:,1))
          CALL ConvertUP2EP_hyb_l(1, UPs_1D, EPs_1D)
          IF (par%ismasterproc) WRITE(iulog,*) 'Read ps'
       END SELECT
       psfact=1.D0
       IF( EPs_1D(1,1) < 2000. ) psfact=1.D+2  ! If unit of PS is hPa, convert into Pa
       !IF( x2d(1) < 2000. ) psfact=1.D+2  ! If unit of PS is hPa, convert into Pa
       icnt = 0
       DO ie= 1, nelemd
       DO j = 1, np
       DO i = 1, np
          icnt = icnt + 1
          bg_md(i,j,nlevar,ie,itime) = EPs_1D(icnt,1)*psfact
          !bg_md(i,j,nlevar,ie,itime) = x2d(icnt)*psfact
       END DO
       END DO
       END DO
    END DO

!! u10m
    CALL GetVariable(io_t, 'u10m', UPs_1D(:,1))
    CALL ConvertUP2EP_hyb_l(1, UPs_1D, EPs_1D)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read u10m'
    icnt = 0
    DO ie= 1, nelemd
    DO j = 1, np
    DO i = 1, np
       icnt = icnt + 1
       u10m(i,j,ie,itime) = EPs_1D(icnt,1)
    END DO
    END DO
    END DO

!! v10m
    CALL GetVariable(io_t, 'v10m', UPs_1D(:,1))
    CALL ConvertUP2EP_hyb_l(1, UPs_1D, EPs_1D)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read v10m'
    icnt = 0
    DO ie= 1, nelemd
    DO j = 1, np
    DO i = 1, np
       icnt = icnt + 1
       v10m(i,j,ie,itime) = EPs_1D(icnt,1)
    END DO
    END DO
    END DO

!! t2m
    CALL GetVariable(io_t, 't2m', UPs_1D(:,1))
    CALL ConvertUP2EP_hyb_l(1, UPs_1D, EPs_1D)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read t2m'
    icnt = 0
    DO ie= 1, nelemd
    DO j = 1, np
    DO i = 1, np
       icnt = icnt + 1
       t2m(i,j,ie,itime) = EPs_1D(icnt,1)
    END DO
    END DO
    END DO

!! q2m
    CALL GetVariable(io_t, 'q2m', UPs_1D(:,1))
    IF (is_moist_mr) THEN
       DO i = 1, nUPs_l 
          UPs_1D(i,1) = UPs_1D(i,1)/(1.D0+UPs_1D(i,1))  ! to specific humidity
       END DO
    END IF
    CALL ConvertUP2EP_hyb_l(1, UPs_1D, EPs_1D)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read q2m'
    icnt = 0
    DO ie= 1, nelemd
    DO j = 1, np
    DO i = 1, np
       icnt = icnt + 1
       q2m(i,j,ie,itime) = EPs_1D(icnt,1)
    END DO
    END DO
    END DO

!! tsfc
    CALL GetVariable(io_t, 'tsfc', UPs_1D(:,1))
    CALL ConvertUP2EP_hyb_l(1, UPs_1D, EPs_1D)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read tsfc'
    icnt = 0
    DO ie= 1, nelemd
    DO j = 1, np
    DO i = 1, np
       icnt = icnt + 1
       tsfc(i,j,ie,itime) = EPs_1D(icnt,1)
    END DO
    END DO
    END DO

!! sfctype
    CALL GetVariable(io_t, 'slmsk', UPs_1D(:,1))
    CALL ConvertUP2EP_hyb_l(1, UPs_1D, EPs_1D)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read sfctype'
    icnt = 0
    DO ie= 1, nelemd
    DO j = 1, np
    DO i = 1, np
       icnt = icnt + 1
       sfctype(i,j,ie,itime) = EPs_1D(icnt,1)
    END DO
    END DO
    END DO

!! znt
    CALL GetVariable(io_t, 'znt', UPs_1D(:,1))
    CALL ConvertUP2EP_hyb_l(1, UPs_1D, EPs_1D)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read znt'
    icnt = 0
    DO ie= 1, nelemd
    DO j = 1, np
    DO i = 1, np
       icnt = icnt + 1
       rl_org(i,j,ie) = EPs_1D(icnt,1)
       !rl_org(i,j,ie) = x2d(icnt)
    END DO
    END DO
    END DO

!! topo
    CALL GetVariable(io_t, 'topo', UPs_1D(:,1))
    CALL ConvertUP2EP_hyb_l(1, UPs_1D, EPs_1D)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read topo'
    icnt = 0
    DO ie= 1, nelemd
    DO j = 1, np
    DO i = 1, np
       icnt = icnt + 1
       zsfc(i,j,ie) = EPs_1D(icnt,1)
       !zsfc(i,j,ie) = x2d(icnt)
    END DO
    END DO
    END DO

    IF( model_name == 'kim-sw' )THEN          ! Read p and pint for KIM-SW
    
       CALL GetVariable(io_t, 'p', UPs_2D)
       CALL ConvertUP2EP_hyb_l(nlev, UPs_2D, EPs_2D)
       IF (par%ismasterproc) WRITE(iulog,*) 'Read pressure at middle level'
       DO k = 1, nlev
       icnt = 0
       DO ie= 1, nelemd
       DO j = 1, np
       DO i = 1, np
          icnt = icnt + 1
          pmid(i,j,nlev+1-k,ie,itime) = EPs_2D(icnt,k)
          !pmid(i,j,nlev+1-k,ie,itime) = x3d(icnt)
       END DO
       END DO
       END DO
       END DO
   
       CALL GetVariable(io_t, 'pint', UPs_2D_p1)
       CALL ConvertUP2EP_hyb_l(nlev+1, UPs_2D_p1, EPs_2D_p1)
       IF (par%ismasterproc) WRITE(iulog,*) 'Read pressure at interface level'
       DO k = 1, nlev+1            
       icnt = 0                    
       DO ie= 1, nelemd         
       DO j = 1, np      
       DO i = 1, np
          icnt = icnt + 1
          pint(i,j,nlev+2-k,ie) = EPs_2D_p1(icnt,k)
          !pint(i,j,nlev+2-k,ie) = x3di(icnt)
       END DO
       END DO
       END DO
       END DO

    ELSE                                      ! Calculate p and pint for KIM-SH

       DO k = 1, nlev
          DO ie= 1, nelemd
          DO j = 1, np
          DO i = 1, np
             pmid(i,j,k,ie,itime) = hyam(k)*ps0 + hybm(k)*bg_md(i,j,nlevar,ie,itime)
             pint(i,j,k,ie) = hyai(k)*ps0 + hybi(k)*bg_md(i,j,nlevar,ie,itime)
          END DO
          END DO
          END DO
       END DO
          k = nlev+1
          DO ie= 1, nelemd
          DO j = 1, np
          DO i = 1, np
             pint(i,j,k,ie) = hyai(k)*ps0 + hybi(k)*bg_md(i,j,nlevar,ie,itime)
          END DO
          END DO
          END DO

    ENDIF   ! model_name

    IF( model_name == 'kim-sw' )THEN           ! Read hgt for KIM-SW

       CALL GetVariable(io_t, 'hgt', UPs_2D)
       CALL ConvertUP2EP_hyb_l(nlev, UPs_2D, EPs_2D)
       IF (par%ismasterproc) WRITE(iulog,*) 'Read hgt'
       DO k = 1, nlev
       icnt = 0
       DO ie= 1, nelemd
       DO j = 1, np
       DO i = 1, np
          icnt = icnt + 1
          z_org(i,j,nlev+1-k,ie) = EPs_2D(icnt,k)
          !z_org(i,j,nlev+1-k,ie) = x3d(icnt)
       END DO
       END DO
       END DO
       END DO

    ELSE                                      ! Calculate hgt for KIM-SH

       IF (par%ismasterproc) WRITE(iulog,*) 'Calculate hgt'

    ENDIF   ! model_name

    IF (cv_opt_wind .eq. 1) THEN
       ALLOCATE(bg_ct(np,np,nlev,nelemd))
       bg_uv_rot(:,:,:,1,:,itime)=bg_md(:,:,1:nlev,:,itime)
       bg_uv_rot(:,:,:,2,:,itime)=bg_md(:,:,nlev+1:2*nlev,:,itime)
       CALL CalPsi(bg_uv_rot(:,:,:,:,:,itime), &
                   bg_ct)
       CALL PsiToRotWind(bg_ct, &
                         bg_uv_rot(:,:,:,:,:,itime))
       DEALLOCATE(bg_ct)
    END IF

    phis = zsfc *gravi

    IF (par%ismasterproc) WRITE(iulog,*) 'pmid(1,1,:,1,itime):',pmid(1,1,:,1,itime)
    IF (par%ismasterproc) WRITE(iulog,*) 'zsfc(1,1,1):',zsfc(1,1,1)
    IF (par%ismasterproc) WRITE(iulog,*) 'z_org(1,1,:,1):',z_org(1,1,:,1)

    DEALLOCATE( UPs_1D )
    DEALLOCATE( UPs_2D )
    DEALLOCATE( UPs_2D_p1 )
    DEALLOCATE( EPs_1D )
    DEALLOCATE( EPs_2D )
    DEALLOCATE( EPs_2D_p1 )
    CALL CloseFile(io_t)
    CALL FinIOSystem(io_t)

  END SUBROUTINE ReadBack_KIM
!===============================================================================

!===============================================================================
! Modified by HJ Song for less-weighted background error covariance I/O
! On 29JAN2016
!===============================================================================
  SUBROUTINE ReadBackErrCov( l_Reg_avg,  &
                             l_eigf_psi, &
                             l_eigf_chi, &
                             l_eigf_Tps, &
                             l_eigf_hum, &
                             errvar_md )
   USE NETCDF

   TYPE(KIO_t) &
       :: io_t
   REAL(KIND=real_kind), DIMENSION(reg_nproc), INTENT(OUT) &
       :: l_Reg_avg
   REAL(KIND=real_kind), DIMENSION(psi_nproc), INTENT(OUT) &
       :: l_eigf_psi
   REAL(KIND=real_kind), DIMENSION(chi_nproc), INTENT(OUT) &
       :: l_eigf_chi
   REAL(KIND=real_kind), DIMENSION(Tps_nproc), INTENT(OUT) &
       :: l_eigf_Tps
   REAL(KIND=real_kind), DIMENSION(hum_nproc), INTENT(OUT) &
       :: l_eigf_hum

   REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd), INTENT(OUT), OPTIONAL &
       :: errvar_md

! Local variables

   REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE   &
       :: w3d,                                       &
          w3d_chi,                                   &
          w3d_Tps,                                   &
          w3d_q,                                     &
          x2d,                                       &
          c3d,                                       &
          Reg_tmp
   REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE &
       :: t3d, w2d,                                  &
          w2d_chi,                                   &
          w2d_Tps,                                   &
          w2d_q,                                     &
          x3d      

   REAL(KIND=real_kind)   :: psfact
   REAL(KIND=real_kind)   :: tmpmin, tmpmax
   REAL(KIND=real_kind)   :: minv  , maxv
   INTEGER(KIND=int_kind) :: i, j, k, l
   INTEGER(KIND=int_kind) :: ie, icnt, ig, wn, swn
   INTEGER(KIND=int_kind) :: ierr, tmp
   INTEGER(KIND=int_kind) :: vr
   INTEGER(KIND=int_kind) :: vertnum
   INTEGER(KIND=int_kind) :: status,ncid,varid
   LOGICAL                :: exst_varmd

  !=============================================================================
  ! Read Backgroun Error Covariance
  !=============================================================================

!.. Check for existance of model variable's variance
    status = nf90_open( trim(becFile),NF90_NOWRITE,ncid )
    status = nf90_inq_varid( ncid,'varU',varid )
    IF( status /= 0 )THEN
       IF( par%ismasterproc ) WRITE(iulog,*)                      &
       trim(becFile),"No variance of model variables in BackErrCov"  
       exst_varmd = .FALSE.
    ELSE
       exst_varmd = .TRUE.
    END IF
    status = nf90_close(ncid)

    ALLOCATE( t3d(np*np*nelemd,nlevar) )
    !ALLOCATE( w2d(zwn2*neig) )
    ALLOCATE( w2d(zwn_nproc,neig) )
    !ALLOCATE( w3d(zwn2*neig*nlev) )
    ALLOCATE( w3d(psi_nproc) )
    !ALLOCATE( w2d_chi(zwn2*neig_chi) )
    ALLOCATE( w2d_chi(zwn_nproc,neig_chi) )
    !ALLOCATE( w3d_chi(zwn2*neig_chi*nlev) )
    ALLOCATE( w3d_chi(chi_nproc) )
    !ALLOCATE( w2d_Tps(zwn2*neig_Tps) )
    ALLOCATE( w2d_Tps(zwn_nproc,neig_Tps) )
    !ALLOCATE( w3d_Tps(zwn2*neig_Tps*(nlev+1)) )
    ALLOCATE( w3d_Tps(Tps_nproc) )
    !ALLOCATE( w2d_q(zwn2*neig_q) )
    ALLOCATE( w2d_q(zwn_nproc,neig_q) )
    !ALLOCATE( w3d_q(zwn2*neig_q*nlev) )
    ALLOCATE( w3d_q(hum_nproc) )

!.. Open PIO for Backgroun Error Covariance file
    CALL IniIOSystem(io_t, par%comm, par%nprocs, par%iproc)
    CALL OpenFile(io_t, TRIM(becFile))
    CALL GetDimension(io_t, 'ncol', tmp)
    IF (tmp .LE. 0 .OR. tmp .EQ. nUPs_g) THEN
      IF (par%ismasterproc) print*,'Check dimension or EP file: ',TRIM(becFile)
      CALL AbortPar()
    END IF
    CALL GetDimension(io_t, 'ncol', tmp, DOFs_EP)
    CALL GetDimension(io_t, 'nwav', tmp, map_zwn)
    CALL GenMap(reg_nproc,  reg_iproc, map_reg)
    CALL GetDimension(io_t, 'nreg', tmp, map_reg)
    CALL GetDimension(io_t, 'npsi', tmp, map_psi)
    CALL GetDimension(io_t, 'nchi', tmp, map_chi)
    CALL GetDimension(io_t, 'nTps', tmp, map_Tps)
    CALL GetDimension(io_t, 'nhum', tmp, map_hum)

!PRINT *, 'map zwn 1 = ', map_zwn(1)
!PRINT *, 'map reg 1 = ', map_reg(1)

    IF (par%ismasterproc) WRITE(iulog,*) 'File of background error covariance:',becFile

    CALL GetVariable(io_t, 'vari', t3d)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read variance of control variables'

    DO k  = 1, nlevar
       icnt = 0
       DO ie = 1, nelemd
       DO j  = 1, np
       DO i  = 1, np
          icnt = icnt + 1
          errvar(i,j,k,ie) = t3d(icnt,k)
       END DO
       END DO
       END DO
    END DO

    tmpmax= maxval(errvar)
    tmpmin= minval(errvar)
    CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
    CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
    ! maxv = GlobalMAX(errvar) ! USE ParFunctions
    ! minv = GlobalMIN(errvar) !
    IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [errvar]:',minv, maxv

!.. Read eigenmodes

!.. Read eigenvalues
!.. Read eigenvalues of psi
    CALL GetVariable(io_t, 'eigv', w2d)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read eigenvalue'

    !!!ALLOCATE( w2d(zwn_nproc,neig) )
    DO ig=1,neig
    DO wn=1,zwn_nproc
       eigv(ig,wn)= w2d(wn,ig)
    END DO
    END DO
!PRINT *, SIZE(eigv,DIM=1), SIZE(eigv,DIM=2), SIZE(w2d,DIM=1), SIZE(w2d,DIM=2)
    !eigv(:,:) = w2d(:,:)

    tmpmax= maxval(eigv)
    tmpmin= minval(eigv)
    CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
    CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
    IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [eigenValue]:',minv, maxv

!.. Read eigenvalues of chi
    CALL GetVariable(io_t, 'egvc', w2d_chi)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read eigenvalue_chi'

    !icnt = 0
    DO ig=1,neig_chi
    DO wn=1,zwn_nproc
       eigv_chi(ig,wn)= w2d_chi(wn,ig)
    END DO
    END DO

    tmpmax= maxval(eigv_chi)
    tmpmin= minval(eigv_chi)
    CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
    CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
    IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [eigenValue_chi]:',minv, maxv

!.. Read eigenvalues of (T,ps)
    CALL GetVariable(io_t, 'egvT', w2d_Tps)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read eigenvalue_Tps'

    DO ig=1,neig_Tps
    DO wn=1,zwn_nproc
       eigv_Tps(ig,wn)= w2d_Tps(wn,ig)
    END DO
    END DO

    tmpmax= maxval(eigv_Tps)
    tmpmin= minval(eigv_Tps)
    CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
    CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
    IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [eigenValue_Tps]:',minv, maxv

!.. Read eigenvalues of q
    CALL GetVariable(io_t, 'egvq', w2d_q)
    IF (par%ismasterproc) WRITE(iulog,*) 'Read eigenvalue_q'

    DO ig=1,neig_q
    DO wn=1,zwn_nproc
       eigv_q(ig,wn)= w2d_q(wn,ig)
    END DO
    END DO

    tmpmax= maxval(eigv_q)
    tmpmin= minval(eigv_q)
    CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
    CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
    IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [eigenValue_q]:',minv, maxv

! Read eigenvectors
! Read eigenvectors of psi
    CALL GetVariable(io_t, 'eigf', w3d)
    l_eigf_psi=w3d
    IF (par%ismasterproc) WRITE(iulog,*) 'Read eigenfunction'

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
    CALL GetVariable(io_t, 'egfc', w3d_chi)
    l_eigf_chi=w3d_chi
    IF (par%ismasterproc) WRITE(iulog,*) 'Read eigenfunction_chi'

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
    CALL GetVariable(io_t, 'egfT', w3d_Tps)
    l_eigf_Tps=w3d_Tps
    IF (par%ismasterproc) WRITE(iulog,*) 'Read eigenfunction_Tps'

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
    CALL GetVariable(io_t, 'egfq', w3d_q)
    l_eigf_hum=w3d_q
    IF (par%ismasterproc) WRITE(iulog,*) 'Read eigenfunction_q'

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
    tmpmax= maxval(l_eigf_psi)
    tmpmin= minval(l_eigf_psi)
    CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
    CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
    IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [eigenFunction_psi]:',minv, maxv

    tmpmax= maxval(l_eigf_chi)
    tmpmin= minval(l_eigf_chi)
    CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
    CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
    IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [eigenFunction_chi]:',minv, maxv

    tmpmax= maxval(l_eigf_Tps)
    tmpmin= minval(l_eigf_Tps)
    CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
    CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
    IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [eigenFunction_Tps]:',minv, maxv

    tmpmax= maxval(l_eigf_hum)
    tmpmin= minval(l_eigf_hum)
    CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
    CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
    IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [eigenFunction_q]:',minv, maxv

! Read variance 
    IF( exst_varmd )THEN

       ALLOCATE( x2d(np*np*nelemd) )
       ALLOCATE( x3d(np*np*nelemd,nlev) )

       Do vr = 1, nvar3d
          SELECT CASE (vr)
          CASE (1)
             CALL GetVariable(io_t, 'varU', x3d)
             IF (par%ismasterproc) WRITE(iulog,*) 'Read variance of u'
          CASE (2)
             CALL GetVariable(io_t, 'varV', x3d)
             IF (par%ismasterproc) WRITE(iulog,*) 'Read variance of v'
          CASE (3)
             CALL GetVariable(io_t, 'varT', x3d)
             IF (par%ismasterproc) WRITE(iulog,*) 'Read variance of T'
          CASE (4)
             CALL GetVariable(io_t, 'varQ', x3d)
             IF (par%ismasterproc) WRITE(iulog,*) 'Read variance of q'
          END SELECT
          DO k = 1, nlev
             icnt = 0
             DO ie= 1, nelemd
             DO j = 1, np
             DO i = 1, np
                icnt = icnt + 1
                errvar_md(i,j,k+nlev*(vr-1),ie) = x3d(icnt,k)
             END DO
             END DO
             END DO
          END DO
       END DO

       CALL GetVariable(io_t, 'varP', x2d)
       IF (par%ismasterproc) WRITE(iulog,*) 'Read variance of ps'

       icnt = 0
       DO ie= 1, nelemd
       DO j = 1, np
       DO i = 1, np
          icnt = icnt + 1
          errvar_md(i,j,nlevar,ie) = x2d(icnt)
       END DO
       END DO
       END DO

       DEALLOCATE( x2d, x3d )

    END IF  ! model_name

    IF (cv_opt_Tu == 1 .OR. cv_opt_Tu == 2 .OR. cv_opt_Tu ==3) THEN

       !ALLOCATE( c3d(np*np*nelem*nlev*nk) )
       ALLOCATE(c3d(reg_nproc))

       CALL GetVariable(io_t, 'regn', c3d)
       l_Reg_avg=c3d
       IF (par%ismasterproc) WRITE(iulog,*) 'Read Regression Coefficient'

       tmpmax= maxval(l_Reg_avg)
       tmpmin= minval(l_Reg_avg)
       CALL Par_Reduce(tmpmax,maxv,1,PAR_MAX)
       CALL Par_Reduce(tmpmin,minv,1,PAR_MIN)
       IF( par%ismasterproc ) WRITE(iulog,*) 'Min & Max [RegressCoef]:',minv, maxv

       DEALLOCATE( c3d )
       !DEALLOCATE( dimSizesRegCoef )

    END IF  ! cv_opt_Tu

    CALL CloseFile(io_t)
    CALL FinIOSystem(io_t)

    DEALLOCATE( t3d, w2d, w3d )
    DEALLOCATE( w2d_chi, w3d_chi )
    DEALLOCATE( w2d_Tps, w3d_Tps )
    DEALLOCATE( w2d_q, w3d_q )

    RETURN

  END SUBROUTINE ReadBackErrCov
!===============================================================================

!===============================================================================
  SUBROUTINE ReadNature
!===============================================================================

! Local variables

   TYPE(KIO_t) &
       :: io_t
   REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE &
       :: x3d
   REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE &
       :: x2d

   REAL(KIND=real_kind)   :: psfact
   INTEGER(KIND=int_kind) :: i, j, k, ie, icnt
   INTEGER(KIND=int_kind) :: ierr, tmp
   INTEGER(KIND=int_kind) :: vr
   INTEGER(KIND=int_kind) :: vertnum

  !=============================================================================
  ! Read Nature
  !=============================================================================

    ALLOCATE( x3d(np*np*nelemd,nlev) )
    ALLOCATE( x2d(np*np*nelemd) )

    CALL IniIOSystem(io_t, par%comm, par%nprocs, par%iproc)
    CALL OpenFile(io_t, TRIM(natureFile))
    CALL GetDimension(io_t, 'ncol', tmp)
    IF (tmp .LE. 0 .OR. tmp .EQ. nUPs_g) THEN
      IF (par%ismasterproc) print*,'Check dimension or EP file: ',TRIM(natureFile)
      CALL AbortPar()
    END IF
    CALL GetDimension(io_t, 'ncol', tmp, DOFs_EP)

    IF (par%ismasterproc) WRITE(iulog,*) 'File of nature=', natureFile
    DO vr = 1, nvar3d
          SELECT CASE (vr)
          CASE (1)
             CALL GetVariable(io_t, 'u', x3d)
          CASE (2)
             CALL GetVariable(io_t, 'v', x3d)
          CASE (3)
             CALL GetVariable(io_t, 'T', x3d)
          CASE (4)
          IF (cv_opt_hum == 0) THEN
          CALL GetVariable(io_t, 'q', x3d)
          ELSE IF(cv_opt_hum == 1) THEN
          CALL GetVariable(io_t, 'rh', x3d)
          END IF
       END SELECT
       IF( model_name == 'kim-old' )THEN
          DO k  = 1, nlev
          icnt = 0
          DO ie = 1, nelemd
          DO j  = 1, np
          DO i  = 1, np
             icnt = icnt + 1
             nt_org(i,j,k+nlev*(vr-1),ie) = x3d(icnt,k)
          END DO
          END DO
          END DO
          END DO
       ELSE
          DO k  = 1, nlev
          icnt = 0
          DO ie = 1, nelemd
          DO j  = 1, np
          DO i  = 1, np
             icnt = icnt + 1
             nt_org(i,j,nlev+1-k+nlev*(vr-1),ie) = x3d(icnt,k)
          END DO
          END DO
          END DO
          END DO
       ENDIF
    END DO

    DO vr = 1, nvar2d
       SELECT CASE (vr)
          CASE (1)
             CALL GetVariable(io_t, 'ps', x2d)
       END SELECT
       psfact=1.D0
       IF( x2d(1) < 10000. ) psfact=1.D+2  ! If unit of PS is hPa,
                                           ! it will be converted into Pa
       icnt = 0
       DO ie = 1, nelemd
       DO j  = 1, np
       DO i  = 1, np
          icnt = icnt + 1
          nt_org(i,j,nlev*nvar3d+vr,ie) = x2d(icnt) *psfact
       END DO
       END DO
       END DO
    END DO

    CALL CloseFile(io_t)
    CALL FinIOSystem(io_t)

    DEALLOCATE( x3d, x2d )

  END SUBROUTINE ReadNature
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
  SUBROUTINE TempToM( t, &
                      q, &
                      m)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE PhysicalConstants, ONLY : rgas => KIM_rair, ps0 => KIM_p0

  ! 2. Variables
  ! 2-1. Input variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(IN) :: q, t

  ! 2-2. Input/output variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(OUT) :: m

  ! 2-3. Local variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: p, dp, phi, phii, tphi
    REAL(KIND=real_kind) &
       :: hkk, hkl
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie
  !=============================================================================

  !=============================================================================
  ! B. Main body
  !=============================================================================
  
    CALL BarrierPar()
  
  ! 0. T -> Tv
    PRINT *, 'T->Tv', par%iProc
    tphi = 0.D0
    m = 0.D0
    DO ie = nets, nete
    DO k  = 1, nlev
    DO j  = 1, np
    DO i  = 1, np
       tphi(i,j,k,ie) = t(i,j,k,ie)*(1.D0 + (rvap/rgas - 1.D0)* &
                        q(i,j,k,ie))
    END DO
    END DO
    END DO
    END DO

    CALL BarrierPar()

  ! 1. pint -> dp
    PRINT *, 'pint -> dp', par%iProc
    dp = 0.D0
    DO ie = nets, nete
    DO j = 1, np
    DO i = 1, np
    DO k = 1, nlev
       dp(i,j,k,ie) = pint(i,j,k+1,ie) - pint(i,j,k,ie)
    END DO
    END DO
    END DO
    END DO

    CALL BarrierPar()

  ! 2. Temperature(Tv) -> Phi
    PRINT *, 't->phi', par%iProc
    phi  = 0.D0
    phii = 0.D0
    hkk  = 0.D0
    hkl  = 0.D0
    DO ie = nets, nete
    DO j  = 1, np
       DO i = 1, np
          hkk = dp(i,j,nlev,ie)*0.5D0/p(i,j,nlev,ie)
          hkl = 2.D0*hkk
          phii(i,j,nlev,ie) = rgas*tphi(i,j,nlev,ie)*hkl
          phi(i,j,nlev,ie)  = phis(i,j,ie) + rgas*tphi(i,j,nlev,ie)*hkk
       END DO

       DO k = nlev-1,2,-1
          DO i = 1, np
             hkk = dp(i,j,k,ie)*0.5D0/p(i,j,k,ie)
             hkl = 2.D0*hkk
             phii(i,j,k,ie) = phii(i,j,k+1,ie) + rgas*tphi(i,j,k,ie)*hkl
             phi(i,j,k,ie)  = phis(i,j,ie) + phii(i,j,k+1,ie) + &
                              rgas*tphi(i,j,k,ie)*hkk
          END DO
       END DO

       DO i = 1, np
          hkk = 0.5D0*dp(i,j,1,ie)/p(i,j,1,ie)
          phi(i,j,1,ie) = phis(i,j,ie) + phii(i,j,2,ie) + &
                          rgas*tphi(i,j,1,ie)*hkk
       END DO

    END DO
    END DO

    CALL BarrierPar()

    ! M = phi + RTv*ln(p)
    ! Tr = Tvr : virtual temperature at reference level 
    PRINT *, 'phi->m', par%iProc
    DO ie = nets, nete
    DO k = 1, nlev
    DO j = 1, np
    DO i = 1, np
       m(i,j,k,ie) = phi(i,j,k,ie) + rgas*Tr*DLOG(p(i,j,k,ie))
    END DO
    END DO
    END DO
    END DO
  !=============================================================================

    RETURN

  END SUBROUTINE TempToM



  SUBROUTINE ModelToCtrl_nochi(tl_bg_md,  &
                               tl_bg_chi, &
                               gs_md,     &
                               tl_bg_ct)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Variables
  ! 1-1. Input variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd), &
       INTENT(IN) :: tl_bg_md, gs_md
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(IN) :: tl_bg_chi

  ! 1-2. Output variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd), &
       INTENT(OUT) :: tl_bg_ct

  ! 1-3. Local variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
       :: md
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd) &
       :: tl_uv, bg_uv
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: tl_t, tl_q, tl_m, bg_t, tl_rh, bg_rh
    REAL(KIND=real_kind), DIMENSION(np,np,nelemd) &
       :: tl_ps, bg_ps
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: tl_ct
    REAL(KIND=real_kind), DIMENSION(np,np,2*nlev+1,nelemd) &
       :: tl_var
    REAL(KIND=real_kind), DIMENSION(nUPs_l,nlev) &
       :: tl_ct_up
    REAL(KIND=real_kind), DIMENSION(nUPs_l,2*nlev+1) &
       :: tl_var_up
    REAL(KIND=real_kind), DIMENSION(2*nlev+1) &
       :: tmp
    INTEGER(KIND=int_kind) &
       :: ie, i, j, k, iup
       
  !=============================================================================

    md = gs_md
    tl_bg_ct = 0.D0

  !=============================================================================
  ! Main body
  !=============================================================================

    tl_uv(:,:,:,1,:) = tl_bg_md(:,:,1:nlev,:)
    tl_uv(:,:,:,2,:) = tl_bg_md(:,:,nlev+1:2*nlev,:)
    tl_t(:,:,:,:)    = tl_bg_md(:,:,2*nlev+1:3*nlev,:)
    tl_rh(:,:,:,:)   = tl_bg_md(:,:,3*nlev+1:4*nlev,:)
    tl_ps(:,:,:)     = tl_bg_md(:,:,nvar3d*nlev+1,:)

!    bg_uv(:,:,:,1,:) = md(:,:,1:nlev,:)
!    bg_uv(:,:,:,2,:) = md(:,:,nlev+1:2*nlev,:)
    bg_uv(:,:,:,1,:) = md(:,:,1:nlev,:)
    bg_uv(:,:,:,2,:) = md(:,:,nlev+1:2*nlev,:)
    bg_t(:,:,:,:)    = md(:,:,2*nlev+1:3*nlev,:)
    bg_rh(:,:,:,:)   = md(:,:,3*nlev+1:4*nlev,:)
    bg_ps(:,:,:)     = md(:,:,nvar3d*nlev+1,:)

!-----------------------------------------------------------------------------
    tl_bg_ct = tl_bg_md

    IF (cv_opt_wind == 1) THEN
       CALL CalPsi(tl_uv, &
                   tl_ct)
       tl_bg_ct(:,:,1:nlev,:) = tl_ct
    END IF
     
!    IF(cv_opt_psu == 1) THEN
!       !tl_ct = 0.D0
!       CALL CalPsb(tl_uv, &
!                   tl_ct)
!    ELSE IF(cv_opt_psu == 2) THEN
!       CALL CalPsbNlbe(tl_uv, bg_uv, bg_ps, &
!                       tl_ct)
!    END IF
!   
!    tl_ct(:,:,nlev,:) = tl_ct(:,:,nlev,:)*bg_ps(:,:,:)
!
!    ! psu = ps - psb
!    tl_bg_ct(:,:,nvar3d*nlev+1,:) = tl_bg_md(:,:,nvar3d*nlev+1,:) - &
!                                    tl_ct(:,:,nlev,:)

    IF (cv_opt_Tu .LE. 2) THEN

       CALL PsiToRotWind(tl_bg_ct(:,:,1:nlev,:), &
                         tl_uv)
 
       IF(cv_opt_Tu == 1) THEN
          CALL CalPhib(tl_uv, &
                       tl_ct)

       ELSE IF(cv_opt_Tu == 2) THEN
          CALL CalNLBE(tl_uv, bg_uv, &
                       tl_ct)
       END IF

    END IF ! (cv_opt_Tu .LE. 2)

    ierr = TimingStart('EP2UP_fromModelToCtrl_nochi')
    CALL ConvertEP2UP_hyb(nlev,tl_ct,tl_ct_up)
    ierr = TimingStop('EP2UP_fromModelToCtrl_nochi')

    ! Mb -> Tb, PSb
    tl_var_up = 0.D0
    DO iup  = 1, nUPs_l
       tmp = 0.D0
       tmp = MATMUL(tl_ct_up(iup,:),Regress_up(:,iup,:))    ! 2*nlev+1
    DO k = 1, 2*nlev+1
       tl_var_up(iup,k) = tmp(k)
    END DO
    END DO
    ierr = TimingStart('UP2EP_fromModelToCtrl_nochi')
    CALL ConvertUP2EP_hyb(2*nlev+1,tl_var_up,tl_var)
    ierr = TimingStop('UP2EP_fromModelToCtrl_nochi')

    ! Tu' = T' - Tb'
    tl_bg_ct(:,:,2*nlev+1:3*nlev,:) = tl_bg_md(:,:,2*nlev+1:3*nlev,:) - &
                                      tl_var(:,:,1:nlev,:)    

    ! psu = ps - psb
    tl_bg_ct(:,:,nvar3d*nlev+1,:) = tl_bg_md(:,:,nvar3d*nlev+1,:) - &
                                    !tl_ct(:,:,nlev,:)
                                    tl_var(:,:,nlev+1,:)

!.. Make unbalanced velocity potential
!    tl_uv(:,:,:,1,:) = tl_bg_md(:,:,1:nlev,:)
!    tl_uv(:,:,:,2,:) = tl_bg_md(:,:,nlev+1:2*nlev,:)
!    CALL CalChi(tl_uv, &
!                tl_ct)

    tl_bg_ct(:,:,nlev+1:2*nlev,:) = tl_bg_chi(:,:,:,:) - &
                                    tl_var(:,:,nlev+2:2*nlev+1,:)

!=============================================================================

    RETURN

  END SUBROUTINE ModelToCtrl_nochi

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
  SUBROUTINE ModelToCtrl(tl_bg_md, &
                         gs_md,    &
                         bg_uv,    &
                         tl_bg_ct)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Variables
  ! 1-1. Input variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd), &
       INTENT(IN) :: tl_bg_md, gs_md
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd), &
       INTENT(IN) :: bg_uv

  ! 1-2. Output variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd), &
       INTENT(OUT) :: tl_bg_ct

  ! 1-3. Local variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlevar,nelemd) &
       :: md
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd) &
       :: tl_uv          !, bg_uv
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: tl_t, tl_q, tl_m, bg_t, tl_rh, bg_rh
    REAL(KIND=real_kind), DIMENSION(np,np,nelemd) &
       :: tl_ps, bg_ps
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: tl_ct
    REAL(KIND=real_kind), DIMENSION(np,np,2*nlev+1,nelemd) &
       :: tl_var
    REAL(KIND=real_kind), DIMENSION(2*nlev+1) &
       :: tmp
    REAL(KIND=real_kind), DIMENSION(nUPs_l,nlev) &
       :: tl_ct_up
    REAL(KIND=real_kind), DIMENSION(nUPs_l,2*nlev+1) &
       :: tl_var_up
    INTEGER(KIND=int_kind) &
       :: ie, i, j, k, iup
  !=============================================================================

    md = gs_md
    tl_bg_ct = 0.D0

  !=============================================================================
  ! Main body
  !=============================================================================

    tl_uv(:,:,:,1,:) = tl_bg_md(:,:,1:nlev,:)
    tl_uv(:,:,:,2,:) = tl_bg_md(:,:,nlev+1:2*nlev,:)
    tl_t(:,:,:,:)    = tl_bg_md(:,:,2*nlev+1:3*nlev,:)
    tl_rh(:,:,:,:)   = tl_bg_md(:,:,3*nlev+1:4*nlev,:)
    tl_ps(:,:,:)     = tl_bg_md(:,:,nvar3d*nlev+1,:)

!    bg_uv(:,:,:,1,:) = md(:,:,1:nlev,:)
!    bg_uv(:,:,:,2,:) = md(:,:,nlev+1:2*nlev,:)
    bg_t(:,:,:,:)    = md(:,:,2*nlev+1:3*nlev,:)
    bg_rh(:,:,:,:)   = md(:,:,3*nlev+1:4*nlev,:)
    bg_ps(:,:,:)     = md(:,:,nvar3d*nlev+1,:)

!-----------------------------------------------------------------------------
    tl_bg_ct = tl_bg_md

    IF (cv_opt_wind == 1) THEN
       CALL CalPsi(tl_uv, &
                   tl_ct)
       tl_bg_ct(:,:,1:nlev,:) = tl_ct
    END IF
     
!    IF(cv_opt_psu == 1) THEN
!       !tl_ct = 0.D0
!       CALL CalPsb(tl_uv, &
!                   tl_ct)
!
!    ELSE IF(cv_opt_psu == 2) THEN
!       CALL CalPsbNlbe(tl_uv, bg_uv, bg_ps, &
!                       tl_ct)
!    END IF
!    tl_ct(:,:,nlev,:) = tl_ct(:,:,nlev,:)*bg_ps(:,:,:)
!   
!    ! psu = ps - psb
!    tl_bg_ct(:,:,nvar3d*nlev+1,:) = tl_bg_md(:,:,nvar3d*nlev+1,:) - &
!                                    tl_ct(:,:,nlev,:)

    IF (cv_opt_Tu .LE. 2) THEN
   
       CALL PsiToRotWind(tl_bg_ct(:,:,1:nlev,:), &
                         tl_uv)
 
       IF(cv_opt_Tu == 1) THEN
          CALL CalPhib(tl_uv, &
                       tl_ct)

       ELSE IF(cv_opt_Tu == 2) THEN
          CALL CalNLBE(tl_uv, bg_uv, &
                       tl_ct)
  
       END IF

    END IF ! (cv_opt_Tu .LE. 2)

    ierr = TimingStart('EP2UP_fromModelToCtrl')
    CALL ConvertEP2UP_hyb(nlev,tl_ct,tl_ct_up)
    ierr = TimingStop('EP2UP_fromModelToCtrl')

    ! Mb -> Tb, PSb
    tl_var_up = 0.D0
    DO iup  = 1, nUPs_l
       tmp = 0.D0
       tmp = MATMUL(tl_ct_up(iup,:),Regress_up(:,iup,:))    ! 2*nlev+1
    DO k = 1, 2*nlev+1
       tl_var_up(iup,k) = tmp(k)
    END DO
    END DO
    ierr = TimingStart('UP2EP_fromModelToCtrl')
    CALL ConvertUP2EP_hyb(2*nlev+1,tl_var_up,tl_var)
    ierr = TimingStop('UP2EP_fromModelToCtrl')

    ! Tu' = T' - Tb'
    tl_bg_ct(:,:,2*nlev+1:3*nlev,:) = tl_bg_md(:,:,2*nlev+1:3*nlev,:) - &
                                      tl_var(:,:,1:nlev,:)

    ! psu = ps - psb
    tl_bg_ct(:,:,nvar3d*nlev+1,:) = tl_bg_md(:,:,nvar3d*nlev+1,:) - &
                                    !tl_ct(:,:,nlev,:)
                                    tl_var(:,:,nlev+1,:)

!.. Make unbalanced velocity potential
    tl_uv(:,:,:,1,:) = tl_bg_md(:,:,1:nlev,:)
    tl_uv(:,:,:,2,:) = tl_bg_md(:,:,nlev+1:2*nlev,:)
    CALL CalChi(tl_uv, &
                tl_ct)

    tl_bg_ct(:,:,nlev+1:2*nlev,:) = tl_ct - &
                                    tl_var(:,:,nlev+2:2*nlev+1,:)

!=============================================================================

    RETURN

  END SUBROUTINE ModelToCtrl

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
  SUBROUTINE PsiToRotWind(psi, &
                          urot)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE Derivative, ONLY : Gradient_Sphere
    USE Edge,       ONLY : edgebuffer_t,   &
                           InitEdgeBuffer, &
                           FreeEdgeBuffer, &
                           EdgeVpack, EdgeVunpack
    USE Bndry,      ONLY : Bndry_ExchangeV

  ! 2. Variables
  ! 2-1. Input variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(IN) :: psi

  ! 2-2. Output variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd), &
       INTENT(OUT) :: urot

  ! 2-3. Local variables
    TYPE(edgebuffer_t) &
       :: edge2
    REAL(KIND=real_kind), DIMENSION(np,np,2) &
       :: utmp
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie
  !=============================================================================

  !=============================================================================
  ! B. Psi -> rotational winds
  !=============================================================================
    CALL InitEdgeBuffer(edge2, 2*nlev)
    DO k  = 1, nlev
       utmp = 0.D0
       DO ie = nets, nete
          utmp = Gradient_Sphere(psi(:,:,k,ie),deriv(hybrid%ithr), &
                                 elem(ie)%dinv)
          urot(:,:,k,1,ie) = -utmp(:,:,2)
          urot(:,:,k,2,ie) =  utmp(:,:,1)
       END DO
    END DO

    DO ie = nets, nete
       DO k = 1, nlev
          urot(:,:,k,1,ie) = elem(ie)%spheremp(:,:)*urot(:,:,k,1,ie)
          urot(:,:,k,2,ie) = elem(ie)%spheremp(:,:)*urot(:,:,k,2,ie)
       END DO
       CALL EdgeVpack(edge2, urot(:,:,:,1,ie), nlev, &
                      0, elem(ie)%desc)
       CALL EdgeVpack(edge2, urot(:,:,:,2,ie), nlev, &
                      nlev, elem(ie)%desc)
    END DO

    CALL Bndry_ExchangeV(hybrid, &
                         edge2)

    DO ie = nets, nete
       CALL EdgeVunpack(edge2, urot(:,:,:,1,ie), nlev, &
                        0, elem(ie)%desc)
       CALL EdgeVunpack(edge2, urot(:,:,:,2,ie), nlev, &
                        nlev, elem(ie)%desc)
       DO k = 1, nlev
          urot(:,:,k,1,ie) = elem(ie)%rspheremp(:,:)*urot(:,:,k,1,ie)
          urot(:,:,k,2,ie) = elem(ie)%rspheremp(:,:)*urot(:,:,k,2,ie)
       END DO
    END DO
    CALL FreeEdgeBuffer(edge2)
  !=============================================================================

    RETURN

  END SUBROUTINE PsiToRotWind


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
  SUBROUTINE CalPsi(uv,  &
                    psi)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE Derivative, ONLY : vorticity_sphere_wk

  ! 2. Variables
  ! 2-1. Input variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd), &
       INTENT(IN) :: uv

  ! 2-2. Output variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(OUT) :: psi

  ! 2-3. Local variables
    TYPE(edgebuffer_t) &
       :: edge2
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: b
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2) &
       :: uv_ll
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie
  !=============================================================================

  !=============================================================================
  ! B. Main body
  !=============================================================================
    CALL InitEdgeBuffer(edge2, nlev)

    DO ie = nets, nete
    DO k  = 1, nlev
       DO j  = 1, np
       DO i  = 1, np
          uv_ll(i,j,k,1) = uv(i,j,k,1,ie)
          uv_ll(i,j,k,2) = uv(i,j,k,2,ie)
       END DO
       END DO
       b(:,:,k,ie) = vorticity_sphere_wk(uv_ll(:,:,k,:),               &
                                         deriv(hybrid%ithr), elem(ie))
    END DO
       CALL EdgeVpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
    END DO

    CALL Bndry_ExchangeV(hybrid, edge2)

    DO ie = nets, nete
       CALL EdgeVunpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
       DO k = 1, nlev
          b(:,:,k,ie)=elem(ie)%rspheremp(:,:)*b(:,:,k,ie)
       END DO
    END DO

!    psi = 0.D0
!    CALL InvLapCg(b,   &
!                  psi)

    CALL InvLapSpec(b,    & ! In: forcing for Poisson equation
                    nlev, & ! In: vertical level
                    psi)    ! Out: solution
    CALL FreeEdgeBuffer(edge2)
  !=============================================================================

    RETURN

  END SUBROUTINE CalPsi

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
  SUBROUTINE CalChi(uv,  &
                    chi)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE Derivative, ONLY : divergence_sphere_wk

  ! 2. Variables
  ! 2-1. Input variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd), &
       INTENT(IN) :: uv

  ! 2-2. Output variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(OUT) :: chi

  ! 2-3. Local variables
    TYPE(edgebuffer_t) &
       :: edge2
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: b
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2) &
       :: uv_ll
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie
  !=============================================================================

  !=============================================================================
  ! B. Main body 
  !=============================================================================
    CALL InitEdgeBuffer(edge2, nlev)

    DO ie = nets, nete
    DO k  = 1, nlev
       DO j  = 1, np
       DO i  = 1, np
          uv_ll(i,j,k,1) = uv(i,j,k,1,ie)
          uv_ll(i,j,k,2) = uv(i,j,k,2,ie)
       END DO
       END DO
       b(:,:,k,ie) = divergence_sphere_wk(uv_ll(:,:,k,:),               &
                                          deriv(hybrid%ithr), elem(ie))
    END DO
       CALL EdgeVpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
    END DO

    CALL Bndry_ExchangeV(hybrid, edge2)

    DO ie = nets, nete
       CALL EdgeVunpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
       DO k = 1, nlev
          b(:,:,k,ie)=elem(ie)%rspheremp(:,:)*b(:,:,k,ie)
       END DO
    END DO

!    chi = 0.D0
!    CALL InvLapCg(b, &
!                  chi)

    CALL InvLapSpec(b,    &
                    nlev, &
                    chi)
    CALL FreeEdgeBuffer(edge2)
  !=============================================================================

    RETURN

  END SUBROUTINE CalChi


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
  SUBROUTINE CalPhib(uv,   &
                     phib)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE Derivative, ONLY : divergence_sphere_wk

  ! 2. Variables
  ! 2-1. Input variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd), &
       INTENT(IN) :: uv

  ! 2-2. Output variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(OUT) :: phib

  ! 2-3. Local variables
    TYPE(edgebuffer_t) &
       :: edgep2
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: b
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2) &
       :: uv_ll, ucor
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie
  !=============================================================================

  !=============================================================================
  ! B. Main body
  !=============================================================================
    CALL InitEdgeBuffer(edgep2, nlev)

    DO ie = nets, nete
    DO k  = 1, nlev
      DO j  = 1, np
      DO i  = 1, np
         uv_ll(i,j,k,1) = uv(i,j,k,1,ie)
         uv_ll(i,j,k,2) = uv(i,j,k,2,ie)
      END DO
      END DO
      DO j = 1, np
      DO i = 1, np
         ucor(i,j,k,1) = elem(ie)%fcor(i,j)*  uv_ll(i,j,k,2)
         ucor(i,j,k,2) = elem(ie)%fcor(i,j)*(-uv_ll(i,j,k,1))
      END DO
      END DO
      b(:,:,k,ie) = divergence_sphere_wk(ucor(:,:,k,:),                &
                                         deriv(hybrid%ithr), elem(ie))
    END DO
       CALL EdgeVpack(edgep2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
    END DO

    CALL Bndry_ExchangeV(hybrid, edgep2)

    DO ie = nets, nete
       CALL EdgeVunpack(edgep2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
       DO k = 1, nlev
          b(:,:,k,ie)=elem(ie)%rspheremp(:,:)*b(:,:,k,ie)
       END DO
    END DO

!    phib = 0.D0
!    CALL InvLapCg(b,    &
!                  phib)

    CALL InvLapSpec(b,    &
                    nlev, &
                    phib)
    CALL FreeEdgeBuffer(edgep2)
  !=============================================================================

    RETURN

  END SUBROUTINE CalPhib


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
  SUBROUTINE CalPsb(uv,  &
                    psb)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE Derivative, ONLY : divergence_sphere_wk

  ! 2. Variables
  ! 2-1. Input variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd), &
       INTENT(IN) :: uv
    REAL(KIND=real_kind), DIMENSION(np,np,nelemd) &
       :: ps

  ! 2-2. Output variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(OUT) :: psb

  ! 2-3. Local variables
    TYPE(edgebuffer_t) &
       :: edge2
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: b
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2) &
       :: uv_ll, ucor
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie
  !=============================================================================

  !=============================================================================
  ! B. Main body
  !=============================================================================
    CALL InitEdgeBuffer(edge2, 1)

    b(:,:,:,:) = 0.D0
    DO ie = nets, nete
       DO j  = 1, np
       DO i  = 1, np
          uv_ll(i,j,nlev,1) = uv(i,j,nlev,1,ie)
          uv_ll(i,j,nlev,2) = uv(i,j,nlev,2,ie)
       END DO
       END DO
       DO j = 1, np
       DO i = 1, np
          ucor(i,j,nlev,1) = elem(ie)%fcor(i,j)*  uv_ll(i,j,nlev,2)
          ucor(i,j,nlev,2) = elem(ie)%fcor(i,j)*(-uv_ll(i,j,nlev,1))
       END DO
       END DO

       b(:,:,nlev,ie) = divergence_sphere_wk(ucor(:,:,nlev,:),   &
                                             deriv(hybrid%ithr), &
                                             elem(ie))

       CALL EdgeVpack(edge2,b(:,:,nlev,ie),1,0,elem(ie)%desc)
    END DO

    CALL Bndry_ExchangeV(hybrid, edge2)

    DO ie = nets, nete
       CALL EdgeVunpack(edge2,b(:,:,nlev,ie),1,0,elem(ie)%desc)
       b(:,:,nlev,ie)=elem(ie)%rspheremp(:,:)*b(:,:,nlev,ie)
    END DO

    b(:,:,nlev,:) = b(:,:,nlev,:)/(rgas*Tr)

!    psb(:,:,:,:) = 0.D0
!    CALL InvLapCg(b,   &
!                  psb)

!.. Spectral inverse
    CALL InvLapSpec(b(:,:,nlev,:), &
                    1,             &
                    psb(:,:,nlev,:))
    CALL FreeEdgeBuffer(edge2)
  !=============================================================================

    RETURN

  END SUBROUTINE CalPsb


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
  SUBROUTINE CalPsbNlbe(uv, uvb, ps, &
                        psb)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE Derivative, ONLY : divergence_sphere_wk, &
                           Gradient_Sphere

  ! 2. Variables
  ! 2-1. Input variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd), &
       INTENT(IN) :: uv, uvb
    REAL(KIND=real_kind), DIMENSION(np,np,nelemd), &
       INTENT(IN) :: ps

  ! 2-2. Output variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(OUT) :: psb

  ! 2-3. Local variables
    TYPE(edgebuffer_t) &
       :: edgep2
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd) &
       :: uadv
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: b
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2) &
       :: uv_ll, ucor, adv
    REAL(KIND=real_kind), DIMENSION(np,np,2) &
       :: gradv, gradvb, ugradv, ugradvb
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie, comp
  !=============================================================================

 !=============================================================================
  ! B. Main body
  !=============================================================================
    CALL InitEdgeBuffer(edgep2, 2)
    b(:,:,:,:) = 0.D0

    ! 2. advection 
    ! - (uvb*grad uv + uv*grad uvb)
    uadv = 0.D0
    DO ie = nets, nete
       gradv  = 0.D0
       gradvb = 0.D0
       DO k  = nlev, nlev
          DO comp = 1, 2
          gradv   = Gradient_Sphere(uv(:,:,k,comp,ie), deriv(hybrid%ithr), &
                                  elem(ie)%dinv)
          gradvb  = Gradient_Sphere(uvb(:,:,k,comp,ie), deriv(hybrid%ithr),&
                                  elem(ie)%dinv)
          ugradv(:,:,comp)  = uvb(:,:,k,1,ie)*gradv(:,:,1) + uvb(:,:,k,2,ie)*gradv(:,:,2)
          ugradvb(:,:,comp) = uv(:,:,k,1,ie)*gradvb(:,:,1) + uv(:,:,k,2,ie)*gradvb(:,:,2)
          END DO
          uadv(:,:,k,1,ie) = ugradv(:,:,1) + ugradvb(:,:,1)
          uadv(:,:,k,2,ie) = ugradv(:,:,2) + ugradvb(:,:,2)
       END DO
    END DO

   ! DSS
    DO ie = nets, nete
       DO k = nlev, nlev
          uadv(:,:,k,1,ie) = elem(ie)%spheremp(:,:)*uadv(:,:,k,1,ie)
          uadv(:,:,k,2,ie) = elem(ie)%spheremp(:,:)*uadv(:,:,k,2,ie)
       END DO
       CALL EdgeVpack(edgep2, uadv(:,:,nlev,1,ie), 1, &
                      0, elem(ie)%desc)
       CALL EdgeVpack(edgep2, uadv(:,:,nlev,2,ie), 1, &
                      1, elem(ie)%desc)
    END DO

    CALL Bndry_ExchangeV(hybrid, &
                         edgep2)

    DO ie = nets, nete
       CALL EdgeVunpack(edgep2, uadv(:,:,nlev,1,ie), 1, &
                        0, elem(ie)%desc)
       CALL EdgeVunpack(edgep2, uadv(:,:,nlev,2,ie), 1, &
                        1, elem(ie)%desc)
       DO k = nlev, nlev
          uadv(:,:,k,1,ie) = elem(ie)%rspheremp(:,:)*uadv(:,:,k,1,ie)
          uadv(:,:,k,2,ie) = elem(ie)%rspheremp(:,:)*uadv(:,:,k,2,ie)
       END DO
    END DO

 ! 1. linear balance equation
    DO ie = nets, nete
         ucor = 0.D0
       DO j  = 1, np
       DO i  = 1, np
          ucor(i,j,nlev,1) = elem(ie)%fcor(i,j)*  uv(i,j,nlev,2,ie)
          ucor(i,j,nlev,2) = elem(ie)%fcor(i,j)*(-uv(i,j,nlev,1,ie))
          adv(i,j,nlev,1)  = -uadv(i,j,nlev,1,ie)
          adv(i,j,nlev,2)  = -uadv(i,j,nlev,2,ie)
       END DO
       END DO

       ucor(:,:,nlev,1)  =  ucor(:,:,nlev,1) + adv(:,:,nlev,1)
       ucor(:,:,nlev,2)  =  ucor(:,:,nlev,2) + adv(:,:,nlev,2)
       b(:,:,nlev,ie)    =  divergence_sphere_wk(ucor(:,:,nlev,:),             &
                                                deriv(hybrid%ithr), elem(ie))
       CALL EdgeVpack(edgep2,b(:,:,nlev,ie),1,0,elem(ie)%desc)
    END DO

    CALL Bndry_ExchangeV(hybrid, edgep2)

    DO ie = nets, nete
       CALL EdgeVunpack(edgep2,b(:,:,nlev,ie),1,0,elem(ie)%desc)
       b(:,:,nlev,ie)=elem(ie)%rspheremp(:,:)*b(:,:,nlev,ie)
    END DO

    b(:,:,nlev,:) = b(:,:,nlev,:)/(rgas*Tr)

!    psb(:,:,:,:) = 0.D0
!    CALL InvLapCg(b,   &
!                  psb)

!.. Spectral inverse
    CALL InvLapSpec(b(:,:,nlev,:), &
                    1,             &
                    psb(:,:,nlev,:))
    CALL FreeEdgeBuffer(edgep2)
  !=============================================================================

    RETURN

  END SUBROUTINE CalPsbNlbe


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
  SUBROUTINE CalNLBE(uv, uvb, &
                     phib)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE Derivative, ONLY : divergence_sphere_wk, &
                           Gradient_Sphere

  ! 2. Variables
  ! 2-1. Input variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd), &
       INTENT(IN) :: uv, uvb

  ! 2-2. Output variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(OUT) :: phib

  ! 2-3. Local variables
    TYPE(edgebuffer_t) &
       :: edge2
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2,nelemd) &
       :: uadv
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: b
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,2) &
       :: uv_ll, ucor, adv
    REAL(KIND=real_kind), DIMENSION(np,np,2) &
       :: gradv, gradvb, ugradv, ugradvb
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie, comp
  !=============================================================================

!=============================================================================
  ! B. Main body
  !=============================================================================
    CALL InitEdgeBuffer(edge2, 2*nlev)
    b = 0.D0

    ! 1. advection 
    ! - (uvb*grad uv + uv*grad uvb)
    uadv = 0.D0
    DO ie = nets, nete
       gradv  = 0.D0
       gradvb = 0.D0
       DO k  = 1, nlev
          DO comp = 1, 2
          gradv   = Gradient_Sphere(uv(:,:,k,comp,ie), deriv(hybrid%ithr), &
                                  elem(ie)%dinv)
          gradvb  = Gradient_Sphere(uvb(:,:,k,comp,ie), deriv(hybrid%ithr),&
                                  elem(ie)%dinv)
          ugradv(:,:,comp)  = uvb(:,:,k,1,ie)*gradv(:,:,1) + uvb(:,:,k,2,ie)*gradv(:,:,2)
          ugradvb(:,:,comp) = uv(:,:,k,1,ie)*gradvb(:,:,1) + uv(:,:,k,2,ie)*gradvb(:,:,2)
          END DO
          uadv(:,:,k,1,ie) = ugradv(:,:,1) + ugradvb(:,:,1)
          uadv(:,:,k,2,ie) = ugradv(:,:,2) + ugradvb(:,:,2)
       END DO
    END DO

    ! DSS
    DO ie = nets, nete
       DO k = 1, nlev
          uadv(:,:,k,1,ie) = elem(ie)%spheremp(:,:)*uadv(:,:,k,1,ie)
          uadv(:,:,k,2,ie) = elem(ie)%spheremp(:,:)*uadv(:,:,k,2,ie)
       END DO
       CALL EdgeVpack(edge2, uadv(:,:,:,1,ie), nlev, &
                      0, elem(ie)%desc)
       CALL EdgeVpack(edge2, uadv(:,:,:,2,ie), nlev, &
                      nlev, elem(ie)%desc)
    END DO

    CALL Bndry_ExchangeV(hybrid, &
                         edge2)

    DO ie = nets, nete
       CALL EdgeVunpack(edge2, uadv(:,:,:,1,ie), nlev, &
                        0, elem(ie)%desc)
       CALL EdgeVunpack(edge2, uadv(:,:,:,2,ie), nlev, &
                        nlev, elem(ie)%desc)
       DO k = 1, nlev
          uadv(:,:,k,1,ie) = elem(ie)%rspheremp(:,:)*uadv(:,:,k,1,ie)
          uadv(:,:,k,2,ie) = elem(ie)%rspheremp(:,:)*uadv(:,:,k,2,ie)
       END DO
    END DO

! 1. linear balance equation
    DO ie = nets, nete
    DO k  = 1, nlev
       ucor = 0.D0
       DO j  = 1, np
       DO i  = 1, np
          ucor(i,j,k,1) =  elem(ie)%fcor(i,j)*  uv(i,j,k,2,ie)
          ucor(i,j,k,2) =  elem(ie)%fcor(i,j)*(-uv(i,j,k,1,ie))
          adv(i,j,k,1)  =  -uadv(i,j,k,1,ie)
          adv(i,j,k,2)  =  -uadv(i,j,k,2,ie)
       END DO
       END DO
 
       ucor(:,:,k,1)  =  ucor(:,:,k,1) + adv(:,:,k,1)
       ucor(:,:,k,2)  =  ucor(:,:,k,2) + adv(:,:,k,2)
       b(:,:,k,ie)    =  divergence_sphere_wk(ucor(:,:,k,:),                &
                                              deriv(hybrid%ithr), elem(ie))
    END DO
       CALL EdgeVpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
    END DO

    CALL Bndry_ExchangeV(hybrid, edge2)

    DO ie = nets, nete
       CALL EdgeVunpack(edge2,b(:,:,:,ie),nlev,0,elem(ie)%desc)
       DO k = 1, nlev
          b(:,:,k,ie)=elem(ie)%rspheremp(:,:)*b(:,:,k,ie)
       END DO
    END DO

!    phib = 0.D0
!    CALL InvLapCg(b,    &
!                  phib)

!.. Spectral inverse
    CALL InvLapSpec(b,    &
                    nlev, &
                    phib)
    CALL FreeEdgeBuffer(edge2)
  !=============================================================================

    RETURN

  END SUBROUTINE CalNLBE


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
  SUBROUTINE TlmTempToPhi(tl_ps,  &
                          tl_q,   &
                          tl_t,   &
                          ps,     &
                          q,      &
                          t,      &
                          tl_phi)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE PhysicalConstants, ONLY : rvap => KIM_rvap

  ! 2. Variables
  ! 2-1. Input variables
    REAL(KIND=real_kind), DIMENSION(np,np,nelemd), &
       INTENT(IN) :: ps, tl_ps
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(IN) :: q, tl_q, t, tl_t
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: tphi

  ! 2-2. Input/output variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(OUT) :: tl_phi

  ! 2-3. Local variables
    REAL(KIND=real_kind), DIMENSION(np,np,nlev+1,nelemd) &
       :: ph, tl_ph
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: p, dp, phi, phii,                            &
          tl_p, tl_dp, tl_phii, tl_tphi
    REAL(KIND=real_kind)           &
       :: hkk, hkl, tl_hkk, tl_hkl
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie
  !=============================================================================

  !=============================================================================
  ! B. Main body
  !=============================================================================
    ! 0. T -> Tv
    tphi = 0
    tl_tphi = 0
    DO ie = nets, nete
    DO k  = 1, nlev
    DO j  = 1, np
    DO i  = 1, np
       tphi(i,j,k,ie) = t(i,j,k,ie)*(1.D0 + (rvap/rgas - 1.D0)* &
                        q(i,j,k,ie))
       tl_tphi(i,j,k,ie) = tl_t(i,j,k,ie) + &
                           tl_t(i,j,k,ie)*q(i,j,k,ie)*(rvap/rgas - 1.D0) + &
                           t(i,j,k,ie)*tl_q(i,j,k,ie)*(rvap/rgas - 1.D0)
    END DO
    END DO
    END DO
    END DO

    ! 1. ps -> p and dp
    DO ie = nets, nete
    DO k = 1, nlev+1
    DO j = 1, np
    DO i = 1, np
       ph(i,j,k,ie) = hyai(k)*ps0 + hybi(k)*ps(i,j,ie)
    END DO
    END DO
    END DO
    END DO

    DO ie = nets, nete
    DO k = 1, nlev
    DO j = 1, np
    DO i = 1, np
       p(i,j,k,ie) = hyam(k)*ps0 + hybm(k)*ps(i,j,ie)
       dp(i,j,k,ie) = ph(i,j,k+1,ie) - ph(i,j,k,ie)
    END DO
    END DO
    END DO
    END DO

    DO ie = nets, nete
    DO k = 1, nlev+1
    DO j = 1, np
    DO i = 1, np
       tl_ph(i,j,k,ie) = hybi(k)*tl_ps(i,j,ie)
    END DO
    END DO
    END DO
    END DO

    DO ie = nets, nete
    DO k = 1, nlev
    DO j = 1, np
    DO i = 1, np
       tl_p(i,j,k,ie) = hybm(k)*tl_ps(i,j,ie)
       tl_dp(i,j,k,ie) = tl_ph(i,j,k+1,ie) - tl_ph(i,j,k,ie)
    END DO
    END DO
    END DO
    END DO

    ! 2. Temperature (Tv) -> Phi
    phi  = 0.D0
    phii = 0.D0
    hkk  = 0.D0
    tl_phi  = 0.D0
    tl_phii = 0.D0
    tl_hkk  = 0.D0

    DO ie = nets, nete
    DO j  = 1, np
       DO i = 1, np
          hkk = dp(i,j,nlev,ie)*0.5D0/p(i,j,nlev,ie)
          hkl = 2.D0*hkk
          tl_hkk = tl_dp(i,j,nlev,ie)*0.5D0/p(i,j,nlev,ie) -     &
                   dp(i,j,nlev,ie)*0.5D0*tl_p(i,j,nlev,ie)/      &
                   (p(i,j,nlev,ie)**2.D0)
          tl_hkl = 2.D0*tl_hkk
          tl_phii(i,j,nlev,ie) = rgas*tl_tphi(i,j,nlev,ie)*hkl + &
                                 rgas*tphi(i,j,nlev,ie)*tl_hkl
          tl_phi(i,j,nlev,ie) = rgas*tl_tphi(i,j,nlev,ie)*hkk + &
                                rgas*tphi(i,j,nlev,ie)*tl_hkk
       END DO

       DO k = nlev-1,2,-1
          DO i = 1, np
             hkk = dp(i,j,k,ie)*0.5D0/p(i,j,k,ie)
             hkl = 2.D0*hkk
             tl_hkk = tl_dp(i,j,k,ie)*0.5D0/p(i,j,k,ie) -     &
                      dp(i,j,k,ie)*0.5D0*tl_p(i,j,k,ie)/(p(i,j,k,ie)**2.D0) ! &
             tl_hkl = 2.D0*tl_hkk
             tl_phii(i,j,k,ie) = tl_phii(i,j,k+1,ie) + &
                                 rgas*tl_tphi(i,j,k,ie)*hkl + &
                                 rgas*tphi(i,j,k,ie)*tl_hkl
             tl_phi(i,j,k,ie)  = tl_phii(i,j,k+1,ie) + &
                                 rgas*tl_tphi(i,j,k,ie)*hkk + &
                                 rgas*tphi(i,j,k,ie)*tl_hkk
          END DO
       END DO

       DO i = 1, np
          hkk = 0.5D0*dp(i,j,1,ie)/p(i,j,1,ie)
          tl_hkk = tl_dp(i,j,1,ie)*0.5D0/p(i,j,1,ie) -     &
                   dp(i,j,1,ie)*0.5D0*tl_p(i,j,1,ie)/(p(i,j,1,ie)**2.D0) ! &
          tl_phi(i,j,1,ie) = tl_phii(i,j,2,ie) + &
                          rgas*tl_tphi(i,j,1,ie)*hkk + &
                          rgas*tphi(i,j,1,ie)*tl_hkk
       END DO
    END DO
    END DO

   ! M = phi + RTv*ln(p)
   ! Tr = Tvr : virtual temperature at reference level 
    DO ie = nets, nete
    DO k = 1, nlev
    DO j = 1, np
    DO i = 1, np
        !phi(i,j,k,ie) = phi(i,j,k,ie) + rgas*Tr*DLOG(p(i,j,k,ie))
        tl_phi(i,j,k,ie) = tl_phi(i,j,k,ie) + rgas*Tr*tl_p(i,j,k,ie)/p(i,j,k,ie)
    END DO
    END DO
    END DO
    END DO

    !tphi = phi
    !tl_tphi = tl_phi
  !=============================================================================

    RETURN

  END SUBROUTINE TlmTempToPhi


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
  SUBROUTINE InvLapCg(b_org, &
                      x)

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE Derivative,    ONLY : laplace_sphere_wk_hyb
    USE KiapsParallel, ONLY : global_shared_buf,    &
                              global_shared_sum
    USE GlobalNorms,   ONLY : Wrap_Repro_Sum

  ! 2. Variables
  ! 2-1. Input
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(IN) :: b_org

  ! 2-2. Input/Output
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd), &
       INTENT(INOUT) :: x

  ! 2-3. Local
    TYPE(EdgeBuffer_t) &
       :: edge2
    REAL(KIND=real_kind), DIMENSION(np,np,nlev,nelemd) &
       :: Ap, r, p, b
    REAL(KIND=real_kind), DIMENSION(:,:), POINTER &
       :: viscosity => NULL()
    REAL(KIND=real_kind) &
       :: rsold, rsnew, alpha, b_norm
    INTEGER(KIND=int_kind) &
       :: i, j, k, ie,     &
          iter
  !=============================================================================

  !=============================================================================
  ! B. Prepare inverse of Laplace operator
  !=============================================================================
  ! 1. Initialize edge buffer
    CALL InitEdgeBuffer(edge2, 2*nlev)

  ! 2. Make b
    DO ie = nets, nete
       DO k = 1, nlev
          b(:,:,k,ie)  = elem(ie)%rspheremp(:,:)*b_org(:,:,k,ie)
       END DO
       CALL EdgeVpack(edge2, b(:,:,:,ie), nlev, &
                      0, elem(ie)%desc)
    END DO
    CALL Bndry_ExchangeV(hybrid, &
                         edge2)
    DO ie = nets, nete
       CALL EdgeVunpack(edge2, b(:,:,:,ie), nlev, &
                        0, elem(ie)%desc)
       DO k = 1, nlev
          b(:,:,k,ie)  = elem(ie)%spheremp(:,:)*b(:,:,k,ie)
       END DO
    END DO
    DO ie = nets, nete
       global_shared_buf(ie,1) = SUM(b(:,:,:,ie)**2.D0)
    END DO
    CALL Wrap_Repro_Sum(nvars=1, comm=hybrid%par%comm)
    b_norm = DSQRT(global_shared_sum(1))
    b = b/b_norm

  ! 3. Make Ax
    DO ie = nets, nete
       DO k = 1, nlev
          x(:,:,k,ie) = elem(ie)%spheremp(:,:)*x(:,:,k,ie)
       END DO
       CALL EdgeVpack(edge2, x(:,:,:,ie), nlev, &
                      0, elem(ie)%desc)
    END DO
    CALL Bndry_ExchangeV(hybrid, &
                         edge2)
    DO ie = nets, nete
       CALL EdgeVunpack(edge2,x(:,:,:,ie), nlev, &
                        0, elem(ie)%desc)
       DO k = 1, nlev
          x(:,:,k,ie) = elem(ie)%rspheremp(:,:)*x(:,:,k,ie)
          Ap(:,:,k,ie) = laplace_sphere_wk_hyb(x(:,:,k,ie), deriv(hybrid%ithr), &
                                           elem(ie), viscosity)
       END DO
    END DO
    DO ie = nets, nete
       DO k = 1, nlev
          Ap(:,:,k,ie) = elem(ie)%rspheremp(:,:)*Ap(:,:,k,ie)
       END DO
       CALL EdgeVpack(edge2, Ap(:,:,:,ie), nlev, &
                      0, elem(ie)%desc)
    END DO
    CALL Bndry_ExchangeV(hybrid, &
                         edge2)
    DO ie = nets, nete
       CALL EdgeVunpack(edge2, Ap(:,:,:,ie), nlev, &
                        0, elem(ie)%desc)
       DO k = 1, nlev
          Ap(:,:,k,ie) = elem(ie)%spheremp(:,:)*Ap(:,:,k,ie)
       END DO
    END DO
    Ap = Ap/b_norm

  ! 4. Initialize r, p, and rsold
    r = b - Ap
    p = r
    DO ie = nets, nete
       global_shared_buf(ie,1) = SUM(r(:,:,:,ie)**2.D0)
    END DO
    CALL Wrap_Repro_Sum(nvars=1, comm=hybrid%par%comm)
    rsold = global_shared_sum(1)
  !=============================================================================

  !=============================================================================
  ! C. Iteration of conjugate gradient
  !=============================================================================
    DO iter = 1, niter_il

  ! 1. Make Ap
       DO ie = nets, nete
          DO k = 1, nlev
             p(:,:,k,ie) = elem(ie)%spheremp(:,:)*p(:,:,k,ie)
          END DO
          CALL EdgeVpack(edge2, p(:,:,:,ie), nlev, &
                         0, elem(ie)%desc)
       END DO
       CALL Bndry_ExchangeV(hybrid, &
                            edge2)
       DO ie = nets, nete
          CALL EdgeVunpack(edge2,p(:,:,:,ie), nlev, &
                           0, elem(ie)%desc)
          DO k = 1, nlev
             p(:,:,k,ie) = elem(ie)%rspheremp(:,:)*p(:,:,k,ie)
             Ap(:,:,k,ie) = laplace_sphere_wk_hyb(p(:,:,k,ie), deriv(hybrid%ithr), &
                                               elem(ie), viscosity)
          END DO
       END DO
       DO ie = nets, nete
          DO k = 1, nlev
             Ap(:,:,k,ie) = elem(ie)%rspheremp(:,:)*Ap(:,:,k,ie)
          END DO
          CALL edgeVpack(edge2, Ap(:,:,:,ie), nlev, &
                         0, elem(ie)%desc)
       END DO
       CALL Bndry_ExchangeV(hybrid, &
                            edge2)
       DO ie = nets, nete
          CALL EdgeVunpack(edge2, Ap(:,:,:,ie), nlev, &
                           0, elem(ie)%desc)
          DO k = 1, nlev
             Ap(:,:,k,ie) = elem(ie)%spheremp(:,:)*Ap(:,:,k,ie)
          END DO
       END DO
       Ap = Ap/b_norm

  ! 2. Calculate alpha, new x, new r, and rsnew
       DO ie = nets, nete
          global_shared_buf(ie,1) = SUM(p(:,:,:,ie)*Ap(:,:,:,ie))
       END DO
       CALL Wrap_Repro_Sum(nvars=1, comm=hybrid%par%comm)
       IF (global_shared_sum(1) .EQ. 0.D0) EXIT
       alpha = rsold/global_shared_sum(1)
       x = x + alpha*p
       r = r - alpha*Ap
       DO ie = nets, nete
          global_shared_buf(ie,1) = SUM(r(:,:,:,ie)**2.D0)
       END DO
       CALL Wrap_Repro_Sum(nvars=1, comm=hybrid%par%comm)
       rsnew = global_shared_sum(1)

  ! 3. Check rsnew
!       IF (par%ismasterproc) &
!          WRITE(iulog,*) 'iter = ', iter, 'rsnew = ', rsnew
       IF (rsnew .LT. cgtol_il) THEN
!          IF (par%ismasterproc) &
!             WRITE(iulog,*) '# of total iteration in InvLapCg = ', iter
          EXIT
       END IF

  ! 4. Calculate new p and replace rsold
       p = r + (rsnew/rsold)*p
       rsold = rsnew

    END DO
    CALL FreeEdgeBuffer(edge2)
  !=============================================================================

    RETURN

  END SUBROUTINE InvLapCg

  !====================================================================
  ! Spectral inverse module for solving Poisson equation
  ! By HJ Song
  ! On 10/08/2015
  !====================================================================
  SUBROUTINE InvLapSpec(b,        & ! In: forcing for Poisson equation
                        vert_lev, & ! In: vertical level of forcing
                        b_inv)      ! Out: solution

  !=============================================================================
  ! A. Declaration
  !=============================================================================
  ! 1. Local modules
    USE CtrlToSpecMod,     ONLY : CtrlToSpec_inv
    USE SpecToCtrlMod,     ONLY : SpecToCtrl_inv

  ! 2. Variables
  ! 2-1. Input variables
    INTEGER(KIND=int_kind), &
      INTENT(IN) :: vert_lev
    REAL(KIND=real_kind), DIMENSION(np,np,vert_lev,nelemd), &
      INTENT(IN) :: b
    REAL(KIND=real_kind), DIMENSION(np,np,vert_lev,nelemd), &
      INTENT(OUT) :: b_inv

    REAL(KIND=real_kind), DIMENSION(nUPs_l,vert_lev) &
      :: b_up,                                       &
         b_inv_up

    REAL(KIND=real_kind), DIMENSION(vert_lev,zwn2) :: b_sp
    INTEGER(KIND=int_kind) &
       :: i, j, l, k, m, ie, iproc, ierr
  !=============================================================================
  ! 1. Calculate spectral coefficients
    ierr = TimingStart('EP2UP_fromInvLapSpec')
    CALL ConvertEP2UP_hyb(vert_lev,b,b_up)
    ierr = TimingStop('EP2UP_fromInvLapSpec')

    ierr = TimingStart('CtrlToSpec_inv_fromInvLapSpec')
    CALL CtrlToSpec_inv(b_up,     & ! _inv version for using different vertical levels
                        vert_lev, &
                        b_sp)
    ierr = TimingStop('CtrlToSpec_inv_fromInvLapSpec')

  ! 2. Inverse Laplacian in the spectral space (new implementation for efficiency)
    ierr = TimingStart('MultiplayInvCoeff_fromInvLapSpec')

    b_sp(:,1)=0.D0
    DO i=2,zwn2
       l=INT(DSQRT(DBLE(i-1)))
       b_sp(:,i)=b_sp(:,i)*(-KIM_REARTH**2.D0/(DBLE(l)*DBLE(l+1)))
    END DO

    ierr = TimingStop('MultiplayInvCoeff_fromInvLapSpec')

  ! 3. To grid space
    ierr = TimingStart('SpecToCtrl_inv_fromInvLapSpec')
    CALL SpecToCtrl_inv(b_sp,     &
                        vert_lev, &
                        b_inv_up)
    ierr = TimingStop('SpecToCtrl_inv_fromInvLapSpec')

    ierr = TimingStart('UP2EP_fromInvLapSpec')
    CALL ConvertUP2EP_hyb(vert_lev,b_inv_up,b_inv)
    ierr = TimingStop('UP2EP_fromInvLapSpec')

  END SUBROUTINE InvLapSpec
!===============================================================================

END MODULE BackgroundMod
