!-------------------------------------------------------------------------------
!
! MODULE KiapsGrid
!
!> @brief
!>  - Module to control model variables
!>
!> @date 17OCT2012
!>  - Junghan Kim : First written
!> @date 27NOV2012
!>  - Junghan Kim : KIAPSGM + CoreCarEGM
!> @date 13MAR2013
!>  - Junghan Kim : First written new KiapsGrid module
!> @date 12MAR2014
!>  - Sukjin Choi : Revise + add CoreNon
!> @date 30JAN2015
!>  - Junghan Kim : Some functions for FVM tracer. (just temporary)
!> @date 25FEB2015
!>  - In-Sun Song : Comments added and codes cleaned up
!> @date 03JUL2015
!>  - Junghan KIM : nUPs_l, nEPs_l, nUPs_g, nEPs_g are added, code clean-up
!> @date 08SEP2015
!>  - Junghan KIM : ConvertEP2UP, ConvertUP2EP
!
!-------------------------------------------------------------------------------
!
#include <KIM.h>
!
!-------------------------------------------------------------------------------
   module kiapsgrid
!
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
   use kiapsbase, only: i4=>kim_int_kind, &
                                                          l4=>kim_log_kind, &
                                                        r4=>kim_real4_kind, &
                                                        r8=>kim_real8_kind, &
                                                    r16=>kim_longdouble_kind
   use kiapsparallel, only: kim_par
   use dimensions, only: nc
   use element, only: element_t
   use derivative, only: derivative_t
   use hybvcoord, only: hvcoord_t, &
                                                                   hvcoord_init
   use logger, only: debugmessage, &
                                                                   logmessage, &
                                                                  infomessage, &
                                                                  warnmessage, &
                                                                 errormessage, &
                                                                 fatalmessage, &
                                                                       tostring
#ifdef FVM_TRAC
   use fvm_control_volume_mod, only: fvm_struct
   use fvm_mod, only: fvm_init1, fvm_init2, fvm_init3
#endif
!
   implicit none
!
!
   private
!
   logical(l4), public :: linited = .false.
   logical(l4), public :: lsetted = .false.
   integer(i4), public :: nups_l, neps_l, nups_g, neps_g
   type(element_t), dimension(:), pointer, public :: elem
   type(derivative_t), dimension(:), allocatable, public :: deriv
   type(hvcoord_t), public :: hvcoord
#ifdef FVM_TRAC
   type(fvm_struct), dimension(:), pointer, public :: fvm
   real(r16) :: fvm_corners(nc+1)
   real(r16) :: fvm_points(nc)
#endif
!
   public :: inigrid, fingrid, convertep2up, convertup2ep
!
!
   interface convertep2up
     module procedure convertep2up_integer
     module procedure convertep2up_real
     module procedure convertep2up_double
   end interface convertep2up
!
   interface convertup2ep
!    MODULE PROCEDURE ConvertUP2EP_Integer
!    MODULE PROCEDURE ConvertUP2EP_Real
     module procedure convertup2ep_double
   end interface convertup2ep
!
   contains
!
!-------------------------------------------------------------------------------
!
!  SUBROUTINE IniGrid
!
!> @brief
!>  - Initialize the KiapsGrid module
!>  - Domain decomposition
!>  - Get the resource values from resourde DB of KiapsFramework module
!>
!> @date 17OCT2012
!>  - Junghan Kim : First written
!> @date 14MAR2013
!>  - Junghan Kim :
!>
!> @param[in] none
!> @param[out] ierr
!> @return none
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine inigrid(ierr)
!
!-------------------------------------------------------------------------------
   use physicalconstants, only: dd_pi
   use kiapsbase, only: kim_iu_log
   use dimensions, only: np, nlev, nelem, nelemd, nelemdmax, nets, nete, &
                                          globaluniquecols, ntrac, qsize, nc
   use kiapsparallel, only: iam, parallel_t, barrierpar, abortpar, &
                                global_shared_buf, nrepro_vars, par_integer, &
                              par_double_precision, par_max, par_sum, finpar
   use kiapsparallel, only: nthreads, omp_get_thread_num, omp_set_num_threads, &
                                                                     kim_hybrid
   use kiapsparmpi, only: par_allreduce
   use errorhandle, only: errorhandler, kim_success, kim_error_warning
   use control, only: integration, &
                                  topology, partmethod, vfile_mid, vfile_int
   use timestep, only: tl, nmax, timelevel_t, time_at, timelevel_init
   use massmatrix, only: mass_matrix
   use cube, only: cubeedgecount, cubeelemcount, cubetopology, &
                                      cube_init_atomic, rotation_init_atomic, &
                           set_corner_coordinates, assign_node_numbers_to_elem
   use metagraph, only: metavertex_t, metaedge_t, &
                                                 localelemcount, initmetagraph
   use gridgraph, only: gridvertex_t, gridedge_t
   use scheduler, only: schedule, genedgesched
   use quadrature, only: quadrature_t, test_gauss, test_gausslobatto, gausslobatto
   use spacecurve, only: genspacepart
   use params, only: sfcurve
   use dof, only: global_dof, createuniqueindex, setelemoffset
   use reprosummod, only: reprosum, reprosumdefaultopts, reprosumsetopts
   use derivative, only: derivinit
#ifdef FVM_TRAC
   use fvm_mod, only: fvm_init1, fvm_init2, fvm_init3
#endif
!
!   USE KiapsParallel,  ONLY : AbortPar
!   USE Cube,           ONLY : CubeElemCount
!   USE ReproSumMod,    ONLY : ReproSumDefaultOpts, ReproSumSetOpts
!
   implicit none
!
!
   integer(i4), intent(  out), optional :: ierr
!
   type(quadrature_t) :: gp
   type(gridvertex_t), dimension(:), allocatable, target :: gridvertex
   type(gridedge_t), dimension(:), allocatable, target :: gridedge
   type(metavertex_t), dimension(:), allocatable, target :: metavertex
   type(metaedge_t), dimension(:), allocatable, target :: metaedge
   real(r8), dimension(:,:), allocatable :: aratio
   real(r8), dimension(1) :: area
   real(r16), dimension(nc+1) :: cslam_corners
   real(r16), dimension(nc) :: cslam_points(nc)
   real(r8) :: xtmp
!
   character(len=80) rot_type ! cube edge rotation type
!
   integer, dimension(:), allocatable :: tailpartition(:)
   integer, dimension(:), allocatable :: headpartition(:)
!
   integer :: total_nelem
   real(r8) :: approx_elements_per_task
   logical :: repro_sum_use_ddpdd, repro_sum_recompute
   real(r8) :: repro_sum_rel_diff_max
!
   integer :: nelem_edge, nedge
   integer :: i, ii, ie, ith, l_ierr
!
! Begins
! ------
!
!
   call infomessage('Start:grid initialize..')
   linited = .true.
!
   total_nelem = cubeelemcount()
   approx_elements_per_task = dble(total_nelem)/dble(kim_par%nprocs)
   if (approx_elements_per_task<1.0d0) then
     call fatalmessage('Number of elements = '//tostring(total_nelem))
     call fatalmessage('Number of processes = '//tostring(kim_par%nprocs))
     call abortpar(message = 'There is not enough parallelism in '//&
                                  'the job,that is,there is less than '//&
                                                     'one elements per task.')
   endif
!
! Initializes reproducible sum module
! -----------------------------------
!
   call reprosumdefaultopts(&
                                repro_sum_use_ddpdd_out = repro_sum_use_ddpdd,&
                          repro_sum_rel_diff_max_out = repro_sum_rel_diff_max,&
                                repro_sum_recompute_out = repro_sum_recompute)
   call reprosumsetopts(&
                                 repro_sum_use_ddpdd_in = repro_sum_use_ddpdd,&
                           repro_sum_rel_diff_max_in = repro_sum_rel_diff_max,&
                                 repro_sum_recompute_in = repro_sum_recompute)
!
   call logmessage('Initialized repro_sum')
!
   call logmessage(' ')
   call logmessage('Total simulation time = '//tostring(time_at(nmax)))
   call logmessage(' ')
!
!   IF (KIM_par%isMasterProc) then
!      ! =============================================
!      ! Compute total simulated time...
!      ! =============================================
!      write(KIM_IU_LOG,*)" "
!      write(KIM_IU_LOG,*)" total simulated time = ",Time_at(nmax)
!      write(KIM_IU_LOG,*)" "
!      ! =============================================
!      ! Perform Gauss/Gauss Lobatto tests...
!      ! =============================================
!      CALL test_gauss(np)
!      CALL test_gausslobatto(np)
!   END IF
!
! ===============================================================
! Allocate and initialize the graph (array of GridVertex_t types)
! ===============================================================
   if (topology=="cube") then
!if (KIM_par%isMasterProc) then
!   write(KIM_IU_LOG,*)"creating cube topology..."
!END IF
     call logmessage("creating cube topology:start...")
     nelem = cubeelemcount()
     nelem_edge = cubeedgecount()
     allocate(gridvertex(nelem))
     allocate(gridedge(nelem_edge))
     call cubetopology(gridedge,gridvertex)
!if(KIM_par%isMasterProc)       write(KIM_IU_LOG,*)"...done."
     call logmessage("creating cube topology:done...")
   endif
!if(KIM_par%isMasterProc) write(KIM_IU_LOG,*)"partitioning graph..."
   call logmessage("partitioning graph:start...")
   if (partmethod.eq.sfcurve) then
     call genspacepart(gridedge,gridvertex)
   else
!CALL FinPar('not support Metis')
     call fatalmessage("not support metis...")
   endif
! ===========================================================
! given partition,count number of local element descriptors
! ===========================================================
   allocate(metavertex(1))
   allocate(schedule(1))
   nelem_edge = size(gridedge)
   allocate(tailpartition(nelem_edge))
   allocate(headpartition(nelem_edge))
   do i = 1,nelem_edge
     tailpartition(i) = gridedge(i)%tail%processor_number
     headpartition(i) = gridedge(i)%head%processor_number
   enddo
! ====================================================
!  Generate the communication graph
! ====================================================
   call initmetagraph(iam,metavertex(1),gridvertex,gridedge)
   nelemd = localelemcount(metavertex(1))
   if (nelemd.le.0) then
     call fatalmessage("not yet ready to handle nelemd = 0 yet")
!       CALL AbortPar(message='Not yet ready to handle nelemd = 0 yet' )
!       stop
   endif
#ifdef _MPI
   call mpi_allreduce(nelemd,nelemdmax,1,&
                                par_integer,par_max,kim_par%comm,l_ierr)
#else
   nelemdmax = nelemd
#endif
   allocate(elem(nelemd))
   nets = 1
   nete = nelemd
#ifdef FVM_TRAC
   if (ntrac>0) then
     allocate(fvm(nelemd))
   else
     allocate(fvm(0))
   endif
   do i = 1,nc+1
     fvm_corners(i) = 2*(i-1)/dble(nc)-1 ! [-1,1] including end points
   enddo
   do i = 1,nc
     fvm_points(i) =(fvm_corners(i)+fvm_corners(i+1))/2
   enddo
#endif
! ====================================================
!  Generate the communication schedule
! ====================================================
   call genedgesched(elem,iam,schedule(1),metavertex(1))
   allocate(global_shared_buf(nelemd,nrepro_vars))
   call barrierpar(par = kim_par)
   gp = gausslobatto(np) ! gll points
! CSLAM nodes are equally spaced in alpha/beta
! HOMME with equ-angular gnomonic projection maps alpha/beta space
! to the reference element via simple scale + translation
! thus,CSLAM nodes in reference element [-1,1] are a tensor product of
! array 'cslam_corners(:)' computed below:
   xtmp = nc
   do i = 1,nc+1
     cslam_corners(i) = 2*(i-1)/xtmp-1 ! [-1,1] including end points
   enddo
   do i = 1,nc
     cslam_points(i) =(cslam_corners(i)+cslam_corners(i+1))/2
   enddo
   if (topology=="cube") then
!IF(KIM_par%isMasterProc) write(KIM_IU_LOG,*) "initializing cube elements..."
     call logmessage("initializing cube elements:start...")
     do ie = 1,nelemd
       call set_corner_coordinates(elem(ie))
     enddo
     call assign_node_numbers_to_elem(elem,gridvertex)
     do ie = 1,nelemd
       call cube_init_atomic(elem(ie),gp%points)
     enddo
   endif
! =================================================================
! Initialize mass_matrix
! =================================================================
!IF(KIM_par%isMasterProc) write(KIM_IU_LOG,*) 'running mass_matrix'
   call logmessage("initializing mass_matrix:start...")
   call mass_matrix(kim_par,elem)
   allocate(aratio(nelemd,1))
!
   rot_type = "contravariant"
!
   if (topology=="cube") then
     area = 0
     do ie = 1,nelemd
       aratio(ie,1) = sum(elem(ie)%mp(:,:)*elem(ie)%metdet(:,:))
     enddo
     call reprosum(aratio,area,nelemd,nelemd,1,commid = kim_par%comm)
     area(1) = 4*dd_pi/area(1) ! ratio correction
     deallocate(aratio)
!       IF (KIM_par%isMasterProc) &
!            write(KIM_IU_LOG,'(a,f20.17)') " re-initializing cube elements: area correction=",area(1)
     call logmessage("re-initializing cube elements:area correction = "//tostring(area(1)))
     do ie = 1,nelemd
       call cube_init_atomic(elem(ie),gp%points,area(1))
       call rotation_init_atomic(elem(ie),rot_type)
     enddo
   endif
!    IF(KIM_par%isMasterProc) write(KIM_IU_LOG,*) 're-running mass_matrix'
   call logmessage("re-running mass_matrix")
   call mass_matrix(kim_par,elem)
! =================================================================
! Determine the global degree of freedome for each gridpoint
! =================================================================
!    IF(KIM_Par%isMasterProc) write(KIM_IU_LOG,*) 'running global_dof'
   call logmessage("running global_dof")
   call global_dof(kim_par,elem)
! =================================================================
! Create Unique Indices
! =================================================================
!IF(KIM_par%isMasterProc) write(KIM_IU_LOG,*) 'create unique indices'
   call logmessage("create unique indices")
   do ie = 1,nelemd
     call createuniqueindex(elem(ie)%globalid,elem(ie)%gdofp,elem(ie)%idxp)
   enddo
!IF(KIM_par%isMasterProc) write(KIM_IU_LOG,*) 'apply offset'
   call logmessage("apply offset")
   call setelemoffset(kim_par,elem,globaluniquecols)
   do ie = 1,nelemd
     elem(ie)%idxv=>elem(ie)%idxp
   enddo
   deallocate(gridedge)
   deallocate(gridvertex)
#ifdef FVM_TRAC
   if (ntrac>0) then
     call fvm_init1(kim_par)
   endif
#endif
! ==================================
! Initialize derivative structure
! ==================================
   allocate(deriv(0:nthreads-1))
#ifdef FVM_TRAC
   call derivinit(deriv(kim_hybrid%ithr),fvm_corners,fvm_points)
   if (ntrac>0) then
     call fvm_init2(elem,fvm,kim_hybrid,1,nelemd,tl)
   endif
#else
   call derivinit(deriv(kim_hybrid%ithr),cslam_corners)
#endif
! ==================================
! Initialize the vertical coordinate  (cam initializes hvcoord externally)
! ==================================
   hvcoord = hvcoord_init(vfile_mid,vfile_int,.true.,kim_hybrid%ismasterthr,l_ierr)
   if (l_ierr/= 0) then
!       CALL FinPar(message="error in hvcoord_init")
     call fatalmessage("error in hvcoord_init")
   endif
   if (present(ierr)) ierr = 0
#if defined(CORE_SW)
!sjchoi------------------- 140930 start
   do ie = 1,nelemd
     elem(ie)%state%rn2d = 0._r8
     elem(ie)%state%rc2d = 0._r8
     elem(ie)%state%sn2d = 0._r8
   enddo
!sjchoi------------------- 140930 end
#endif
! ==================================
! Calculate number of unique and entire points
! ==================================
   neps_l = np*np*nelemd
   neps_g = np*np*nelem
   nups_l = 0
   do ie = 1,nelemd
     nups_l = nups_l+elem(ie)%idxp%numuniquepts
   enddo
   call par_allreduce(nups_l,nups_g,1,par_sum)
   call logmessage('nEPs_l = '//tostring(neps_l))
   call logmessage('nEPs_g = '//tostring(neps_g))
   call logmessage('nUPs_l = '//tostring(nups_l))
   call logmessage('nUPs_g = '//tostring(nups_g))
! ==================================
! Boundary Exchange for lat,lon
! ==================================
   call synclatlon()
   call infomessage('End:grid initialize..')
   return
!
   end subroutine inigrid
!==========================================================
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Fin
!
!> @brief
!>  - Finalize the KiapsGrid module
!>
!> @date 17OCT2012
!>  - JH KIM : First written
!> @date 27NOV2012
!>  - JH KIM : add FinGrid for CoreCarEGM
!> @date 14MAR2013
!>  - JH KIM :
!>
!> @param[in] none
!> @param[out] ierr
!> @return none
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine fingrid(ierr)
!-------------------------------------------------------------------------------
   implicit none
!--------------------
!
   integer(i4), optional, intent(  out) :: ierr
!--------------------
!
   call infomessage('Start:grid finalize..')
   deallocate(elem)
   deallocate(deriv)
   if (present(ierr)) ierr = 0
   call infomessage('End:grid finalize..')
   return
!
   end subroutine fingrid
!==========================================================
!-------------------------------------------------------------------------------
!
!   SUBROUTINE SyncLatLon
!
!> @brief
!>  -
!>
!> @date 21JUL2015
!>  - JH KIM : First written
!>
!> @param[in] none
!> @param[out] ierr
!> @return none
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine synclatlon()
!-------------------------------------------------------------------------------
   use dimensions, only: np, nlev, nelem, nelemd, nets, nete
   use edge, only: edgebuffer_t, initedgebuffer, freeedgebuffer, edgevpack, edgevunpack, edgevunpackmin
   use bndry, only: bndry_exchangev
   implicit none
!
   integer(i4) :: ie, ip, jp, ilev, kptr
   type(edgebuffer_t) :: l_edge2d
!
   call initedgebuffer(l_edge2d,2)
   do ie = nets,nete
     kptr = 0
     call edgevpack(l_edge2d,elem(ie)%spherep(:,:)%lat,1,kptr,elem(ie)%desc)
     kptr = 1
     call edgevpack(l_edge2d,elem(ie)%spherep(:,:)%lon,1,kptr,elem(ie)%desc)
   enddo
   call bndry_exchangev(kim_par,l_edge2d)
   do ie = nets,nete
     kptr = 0
     call edgevunpackmin(l_edge2d,elem(ie)%spherep(:,:)%lat,1,kptr,elem(ie)%desc)
     kptr = 1
     call edgevunpackmin(l_edge2d,elem(ie)%spherep(:,:)%lon,1,kptr,elem(ie)%desc)
   enddo
   call freeedgebuffer(l_edge2d)
!
   end subroutine synclatlon
!==========================================================
!-------------------------------------------------------------------------------
!
!   SUBROUTINE ConvertEP2UP
!
!> @brief
!>  -
!>
!> @date 07SEP2015
!>  - JH KIM : First written
!>
!> @param[in] none
!> @param[out] ierr
!> @return none
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine convertep2up_integer(nlevs, eps, ups)
!-------------------------------------------------------------------------------
   use dimensions, only: np, nlev, nelem, nelemd
   implicit none
!
   integer(i4), intent(in   ) :: nlevs
   integer(i4), dimension(np, np, nelemd, nlevs), intent(in   ) :: eps
   integer(i4), dimension(nups_l, nlevs), intent(  out) :: ups
!INTEGER(i4), DIMENSION(nUPs_l, nlevs), INTENT(OUT) :: UPs
!local
   integer(i4) :: ip, jp, ie, ilev, iup_e, iep, iup
!
   do ilev = 1, nlevs
     iup = 0
     do ie = 1,nelemd
       do iup_e = 1, elem(ie)%idxp%numuniquepts
         ip = elem(ie)%idxp%ia(iup_e)
         jp = elem(ie)%idxp%ja(iup_e)
!iep = ip + (jp-1)*np + (ie-1)*np*np
         iup = iup+1
         ups(iup,ilev) = eps(ip,jp,ie,ilev)
       enddo
     enddo
   enddo
!
   end subroutine convertep2up_integer
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine convertep2up_real(nlevs, eps, ups)
!-------------------------------------------------------------------------------
   use dimensions, only: np, nlev, nelem, nelemd
   implicit none
!
   integer(i4), intent(in   ) :: nlevs
   real(r4), dimension(np, np, nelemd, nlevs), intent(in   ) :: eps
   real(r4), dimension(nups_l, nlevs), intent(  out) :: ups
!local
   integer(i4) :: ip, jp, ie, ilev, iup_e, iep, iup
!
   do ilev = 1, nlevs
     iup = 0
     do ie = 1,nelemd
       do iup_e = 1, elem(ie)%idxp%numuniquepts
         ip = elem(ie)%idxp%ia(iup_e)
         jp = elem(ie)%idxp%ja(iup_e)
!iep = ip + (jp-1)*np + (ie-1)*np*np
         iup = iup+1
         ups(iup,ilev) = eps(ip,jp,ie,ilev)
       enddo
     enddo
   enddo
!
   end subroutine convertep2up_real
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine convertep2up_double(nlevs, eps, ups)
!-------------------------------------------------------------------------------
   use dimensions, only: np, nlev, nelem, nelemd
   implicit none
!
   integer(i4), intent(in   ) :: nlevs
   real(r8), dimension(np, np, nelemd, nlevs), intent(in   ) :: eps
   real(r8), dimension(nups_l, nlevs), intent(  out) :: ups
!local
   integer(i4) :: ip, jp, ie, ilev, iup_e, iep, iup
!
   do ilev = 1, nlevs
     iup = 0
     do ie = 1,nelemd
       do iup_e = 1, elem(ie)%idxp%numuniquepts
         ip = elem(ie)%idxp%ia(iup_e)
         jp = elem(ie)%idxp%ja(iup_e)
!iep = ip + (jp-1)*np + (ie-1)*np*np
         iup = iup+1
         ups(iup,ilev) = eps(ip,jp,ie,ilev)
       enddo
     enddo
   enddo
!
   end subroutine convertep2up_double
!==========================================================
!-------------------------------------------------------------------------------
!
!   SUBROUTINE ConvertUP2EP
!
!> @brief
!>  -
!>
!> @date 07SEP2015
!>  - JH KIM : First written
!>
!> @param[in] none
!> @param[out] ierr
!> @return none
!
!-------------------------------------------------------------------------------
!    SUBROUTINE ConvertUP2EP_Integer(nlevs, UPs, EPs)
!    USE Dimensions,     ONLY: np, nlev, nelem, nelemd
!    USE Edge,           ONLY: EdgeBuffer_t, InitEdgeBuffer, FreeEdgeBuffer, EdgeVpack, EdgeVunpack, EdgeVunpackMAX, EdgeVunpackCheck
!    USE Bndry,          ONLY: bndry_exchangeV
!    IMPLICIT NONE
!    INTEGER(i4), INTENT(IN)                                 :: nlevs
!    INTEGER(i4), DIMENSION(nUPs_l, nlevs), INTENT(IN)       :: UPs
!    INTEGER(i4), DIMENSION(np,np,nelemd,nlevs), INTENT(OUT) :: EPs
!
!    !local
!    INTEGER(i4)          :: ip, jp, ie, ilev, iup_e, iep, iup
!    TYPE(EdgeBuffer_t)   :: edge3d
!
!      CALL InitEdgeBuffer(edge3d, nlevs)
!
!      EPs(:,:,:,:) = -HUGE(1)
!
!      DO ilev = 1, nlevs
!        iup = 0
!        DO ie = 1, nelemd
!          DO iup_e = 1, elem(ie)%idxP%NumUniquePts
!            ip=elem(ie)%idxP%ia(iup_e)
!            jp=elem(ie)%idxP%ja(iup_e)
!            !iep = ip + (jp-1)*np + (ie-1)*np*np
!            iup = iup + 1
!            EPs(ip,jp,ie,ilev) = UPs(iup,ilev)
!          END DO
!        END DO
!      END DO
!
!
!      DO ie = 1, nelemd
!        ilev = 0
!        CALL EdgeVpack(edge3d, EPs(:,:,ie,:),nlev,ilev,elem(ie)%desc)
!      END DO
!
!      CALL Bndry_ExchangeV(KIM_Par,edge3d)
!
!      DO ie = 1, nelemd
!        ilev = 0
!        CALL EdgeVunpackMAX(edge3d, EPs(:,:,ie,:),nlev,ilev,elem(ie)%desc)
!      END DO
!
!      CALL FreeEdgeBuffer(edge3d)
!
!    END SUBROUTINE ConvertUP2EP_Integer
!
!
!    SUBROUTINE ConvertUP2EP_Real(nlevs, UPs, EPs)
!    USE Dimensions,     ONLY: np, nlev, nelem, nelemd
!    USE Edge,           ONLY: EdgeBuffer_t, InitEdgeBuffer, FreeEdgeBuffer, EdgeVpack, EdgeVunpack, EdgeVunpackMAX, EdgeVunpackCheck
!    USE Bndry,          ONLY: bndry_exchangeV
!    IMPLICIT NONE
!    INTEGER(i4), INTENT(IN)                              :: nlevs
!    REAL(r4), DIMENSION(nUPs_l, nlevs), INTENT(IN)       :: UPs
!    REAL(r4), DIMENSION(np,np,nelemd,nlevs), INTENT(OUT) :: EPs
!
!    !local
!    INTEGER(i4)          :: ip, jp, ie, ilev, iup_e, iep, iup
!    TYPE(EdgeBuffer_t)   :: edge3d
!
!      CALL InitEdgeBuffer(edge3d, nlevs)
!
!      EPs(:,:,:,:) = -HUGE(1.0)
!
!      DO ilev = 1, nlevs
!        iup = 0
!        DO ie = 1, nelemd
!          DO iup_e = 1, elem(ie)%idxP%NumUniquePts
!            ip=elem(ie)%idxP%ia(iup_e)
!            jp=elem(ie)%idxP%ja(iup_e)
!            !iep = ip + (jp-1)*np + (ie-1)*np*np
!            iup = iup + 1
!            EPs(ip,jp,ie,ilev) = UPs(iup,ilev)
!          END DO
!        END DO
!      END DO
!
!
!      DO ie = 1, nelemd
!        ilev = 0
!        CALL EdgeVpack(edge3d, EPs(:,:,ie,:),nlev,ilev,elem(ie)%desc)
!      END DO
!
!      CALL Bndry_ExchangeV(KIM_Par,edge3d)
!
!      DO ie = 1, nelemd
!        ilev = 0
!        CALL EdgeVunpackMAX(edge3d, EPs(:,:,ie,:),nlev,ilev,elem(ie)%desc)
!      END DO
!
!      CALL FreeEdgeBuffer(edge3d)
!
!    END SUBROUTINE ConvertUP2EP_Real
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine convertup2ep_double(nlevs, ups, eps)
!-------------------------------------------------------------------------------
   use dimensions, only: np, nlev, nelem, nelemd
   use edge, only: edgebuffer_t, initedgebuffer, freeedgebuffer, edgevpack, edgevunpack, edgevunpackmax, edgevunpackcheck
   use bndry, only: bndry_exchangev
   implicit none
!
   integer(i4), intent(in   ) :: nlevs
   real(r8), dimension(nups_l, nlevs), intent(in   ) :: ups
   real(r8), dimension(np, np, nelemd, nlevs), intent(  out) :: eps
!local
   integer(i4) :: ip, jp, ie, ilev, iup_e, iep, iup
   type(edgebuffer_t) :: edge3d
!
   call initedgebuffer(edge3d,nlevs)
   eps(:,:,:,:) =-huge(1.0_r8)
   do ilev = 1,nlevs
     iup = 0
     do ie = 1,nelemd
       do iup_e = 1, elem(ie)%idxp%numuniquepts
         ip = elem(ie)%idxp%ia(iup_e)
         jp = elem(ie)%idxp%ja(iup_e)
!iep = ip + (jp-1)*np + (ie-1)*np*np
         iup = iup+1
         eps(ip,jp,ie,ilev) = ups(iup,ilev)
       enddo
     enddo
   enddo
   do ie = 1,nelemd
     ilev = 0
     call edgevpack(edge3d,eps(:,:,ie,:),nlevs,ilev,elem(ie)%desc)
   enddo
   call bndry_exchangev(kim_par,edge3d)
   do ie = 1,nelemd
     ilev = 0
     call edgevunpackmax(edge3d,eps(:,:,ie,:),nlevs,ilev,elem(ie)%desc)
   enddo
   call freeedgebuffer(edge3d)
!
   end subroutine convertup2ep_double
!==========================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module kiapsgrid
!-------------------------------------------------------------------------------
