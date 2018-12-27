#include <KIM.h>
MODULE KiapsGrid
  USE KiapsBase,                 ONLY : i4         => KIM_INT_KIND,     &   ! i4 is ....
                                        l4         => KIM_LOG_KIND,     &   ! l4 is ....
                                        r4         => KIM_REAL4_KIND,   &    ! r4 is ....
                                        r8         => KIM_REAL8_KIND,   &   ! r8 is ....
                                        r16        => KIM_LONGDOUBLE_KIND     ! i4 is ....
  USE KiapsParallel,             ONLY : KIM_Par
  USE Dimensions,                ONLY : nc
  USE Element,                   ONLY : element_t
  USE Derivative,                ONLY : derivative_t
  USE HybVCoord,                 ONLY : hvcoord_t,  &
                                        hvcoord_init
  USE Logger,                    ONLY : DebugMessage,  &    ! i4 is ....
                                        LogMessage,    &    ! i4 is ....
                                        InfoMessage,   &    ! i4 is ....
                                        WarnMessage,   &    ! i4 is ....
                                        ErrorMessage,  &    ! i4 is ....
                                        FatalMessage,  &    ! i4 is ....
                                        ToString    ! i4 is ....
#ifdef FVM_TRAC
  USE fvm_control_volume_mod,    ONLY : fvm_struct
  USE fvm_mod,                   ONLY : fvm_init1, fvm_init2, fvm_init3
#endif
  IMPLICIT NONE
  PRIVATE
  LOGICAL(l4),                                    PUBLIC :: linited = .FALSE.
  LOGICAL(l4),                                    PUBLIC :: lsetted = .FALSE.
  INTEGER(i4),                                    PUBLIC :: nUPs_l, nEPs_l, nUPs_g, nEPs_g
  TYPE (element_t),    DIMENSION(:), POINTER,     PUBLIC :: elem
  TYPE (derivative_t), DIMENSION(:), ALLOCATABLE, PUBLIC :: deriv
  TYPE (hvcoord_t),                               PUBLIC :: hvcoord
#ifdef FVM_TRAC
  TYPE (fvm_struct), DIMENSION(:), POINTER, PUBLIC      :: fvm
  REAL(r16)  :: fvm_corners(nc+1)
  REAL(r16)  :: fvm_points(nc)
#endif
  PUBLIC :: IniGrid, FinGrid, ConvertEP2UP, ConvertUP2EP

  INTERFACE ConvertEP2UP
    MODULE PROCEDURE ConvertEP2UP_Integer
    MODULE PROCEDURE ConvertEP2UP_Real
    MODULE PROCEDURE ConvertEP2UP_Double
  END INTERFACE ConvertEP2UP

  INTERFACE ConvertUP2EP
    MODULE PROCEDURE ConvertUP2EP_Double
  END INTERFACE ConvertUP2EP

CONTAINS
  SUBROUTINE IniGrid(ierr)
    USE PhysicalConstants, ONLY : DD_PI
    USE KiapsBase,      ONLY : KIM_IU_LOG
    USE Dimensions,     ONLY : np, nlev, nelem, nelemd, nelemdmax, nets, nete,    &
                               GlobalUniqueCols, ntrac, qsize, nc
    USE KiapsParallel,  ONLY : iam, parallel_t, BarrierPar, AbortPar,             &
                               Global_Shared_Buf, nRePro_Vars, PAR_INTEGER,       &
                               PAR_DOUBLE_PRECISION, PAR_MAX, PAR_SUM, FinPar
    USE KiapsParallel,  ONLY : nThreads, OMP_GET_THREAD_NUM, OMP_SET_NUM_THREADS, &
                               KIM_Hybrid
    USE KiapsParMPI,    ONLY : Par_AllReduce
    USE ErrorHandle,    ONLY : ErrorHandler, KIM_SUCCESS, KIM_ERROR_WARNING
    USE Control,        ONLY : integration, &
                               topology, partmethod, vfile_mid, vfile_int
    USE TimeStep,       ONLY : tl, nmax, TimeLevel_t, Time_At, TimeLevel_Init
    USE MassMatrix,     ONLY : mass_matrix
    USE Cube,           ONLY : cubeedgecount , cubeelemcount, cubetopology,       &
                               cube_init_atomic, rotation_init_atomic,            &
                               set_corner_coordinates, assign_node_numbers_to_elem
    USE MetaGraph,      ONLY : MetaVertex_t, MetaEdge_t,  &
                               localelemcount, InitMetagraph
    USE Gridgraph,      ONLY : GridVertex_t, GridEdge_t
    USE Scheduler,      ONLY : schedule, genEdgeSched
    USE Quadrature,     ONLY : Quadrature_t, test_gauss, test_gausslobatto, gausslobatto
    USE SpaceCurve,     ONLY : genspacepart
    USE Params,         ONLY : SFCURVE
    USE DOF,            ONLY : global_dof, CreateUniqueIndex, SetElemOffset
    USE ReproSumMod,    ONLY : ReproSum, ReproSumDefaultOpts, ReproSumSetOpts
    USE Derivative,     ONLY : DerivInit
#ifdef FVM_TRAC
    USE fvm_mod,        ONLY : fvm_init1, fvm_init2, fvm_init3
#endif
    IMPLICIT NONE
    INTEGER(i4), INTENT(OUT), OPTIONAL :: ierr
    TYPE(Quadrature_t)  :: gp
    TYPE(GridVertex_t), DIMENSION(:),   ALLOCATABLE, TARGET :: GridVertex
    TYPE(GridEdge_t),   DIMENSION(:),   ALLOCATABLE, TARGET :: Gridedge
    TYPE(MetaVertex_t), DIMENSION(:),   ALLOCATABLE, TARGET :: MetaVertex
    TYPE(MetaEdge_t),   DIMENSION(:),   ALLOCATABLE, TARGET :: MetaEdge
    REAL(r8),           DIMENSION(:,:), ALLOCATABLE         :: aratio
    REAL(r8),           DIMENSION(1)                        :: area
    REAL(r16),          DIMENSION(nc+1) :: cslam_corners
    REAL(r16),          DIMENSION(nc)   :: cslam_points(nc)
    REAL(r8) :: xtmp
    CHARACTER(LEN=80) rot_type   ! cube edge rotation type
    INTEGER, DIMENSION(:), ALLOCATABLE :: TailPartition(:)
    INTEGER, DIMENSION(:), ALLOCATABLE :: HeadPartition(:)
    INTEGER  :: total_nelem
    REAL(r8) :: approx_elements_per_task
    LOGICAL  :: repro_sum_use_ddpdd, repro_sum_recompute
    REAL(r8) :: repro_sum_rel_diff_max
    INTEGER :: nelem_edge, nedge
    INTEGER :: i, ii, ie, ith, l_ierr
    CALL InfoMessage('Start : Grid initialize..')
    linited = .TRUE.
    total_nelem = CubeElemCount()
    approx_elements_per_task = DBLE(total_nelem)/DBLE(KIM_par%nProcs)
    IF (approx_elements_per_task < 1.0D0) THEN
      CALL FatalMessage('Number of elements = '//ToString(total_nelem))
      CALL FatalMessage('Number of processes = '//ToString(KIM_Par%nProcs))
      CALL AbortPar(message='There is not enough parallelism in '//  &
                            'the job, that is, there is less than '//  &
                            'one elements per task.')
    END IF
    CALL ReproSumDefaultOpts(                               &
         repro_sum_use_ddpdd_out=repro_sum_use_ddpdd,       &
         repro_sum_rel_diff_max_out=repro_sum_rel_diff_max, &
         repro_sum_recompute_out=repro_sum_recompute       )
    CALL ReproSumSetOpts(                                   &
         repro_sum_use_ddpdd_in=repro_sum_use_ddpdd,        &
         repro_sum_rel_diff_max_in=repro_sum_rel_diff_max,  &
         repro_sum_recompute_in=repro_sum_recompute        )
    CALL LogMessage('Initialized repro_sum')
    CALL LogMessage(' ')
    CALL LogMessage('Total simulation time = '//ToString(Time_at(nmax)))
    CALL LogMessage(' ')

    IF (topology=="cube") then
       CALL LogMessage("creating cube topology: start...")

       nelem      = CubeElemCount()
       nelem_edge = CubeEdgeCount()

       ALLOCATE(GridVertex(nelem))
       ALLOCATE(GridEdge(nelem_edge))

       CALL CubeTopology(GridEdge,GridVertex)

       CALL LogMessage("creating cube topology: done...")
    END IF


    CALL LogMessage("partitioning graph: start...")
    IF(partmethod .eq. SFCURVE) then
       CALL genspacepart(GridEdge,GridVertex)
    else
       CALL FatalMessage("not support Metis...")
    endif

    ALLOCATE(MetaVertex(1))
    ALLOCATE(Schedule(1))

    nelem_edge=SIZE(GridEdge)

    ALLOCATE(TailPartition(nelem_edge))
    ALLOCATE(HeadPartition(nelem_edge))
    DO i=1,nelem_edge
       TailPartition(i)=GridEdge(i)%tail%processor_number
       HeadPartition(i)=GridEdge(i)%head%processor_number
    END DO

    CALL initMetaGraph(iam,MetaVertex(1),GridVertex,GridEdge)

    nelemd = LocalElemCount(MetaVertex(1))

    IF(nelemd .le. 0) then
       CALL FatalMessage("Not yet ready to handle nelemd = 0 yet")
    endif


#ifdef _MPI 
    CALL mpi_allreduce(nelemd,nelemdmax,1, &
         PAR_INTEGER,PAR_MAX,KIM_par%comm,l_ierr)
#else
    nelemdmax=nelemd
#endif

    ALLOCATE(elem(nelemd))
    nets = 1
    nete = nelemd
#ifdef FVM_TRAC
    IF (ntrac > 0) THEN
      ALLOCATE(fvm(nelemd))
    ELSE
      ALLOCATE(fvm(0))
    END IF
    DO i=1,nc+1
       fvm_corners(i)= 2*(i-1)/DBLE(nc)- 1  ! [-1,1] including end points
    END DO
    DO i=1,nc
       fvm_points(i)= ( fvm_corners(i)+fvm_corners(i+1) ) /2
    END DO
#endif

    CALL genEdgeSched(elem,iam,Schedule(1),MetaVertex(1))


    ALLOCATE(global_shared_buf(nelemd,nrepro_vars))

    CALL BarrierPar(par=KIM_par)

    gp=gausslobatto(np)  ! GLL points
    xtmp=nc
    DO i=1,nc+1
       cslam_corners(i)= 2*(i-1)/xtmp - 1  ! [-1,1] including end points
    END DO
    DO i=1,nc
       cslam_points(i)= ( cslam_corners(i)+cslam_corners(i+1) ) /2
    END DO

    IF (topology=="cube") then
       CALL LogMessage("initializing cube elements: start...")
       DO ie=1,nelemd
          CALL set_corner_coordinates(elem(ie))
       END DO
          CALL assign_node_numbers_to_elem(elem, GridVertex)
       DO ie=1,nelemd
          CALL cube_init_atomic(elem(ie),gp%points)
       END DO
    END IF

    CALL LogMessage("initializing mass_matrix: start...")
    CALL mass_matrix(KIM_par,elem)
    ALLOCATE(aratio(nelemd,1))
    rot_type="contravariant"
    IF (topology=="cube") then
       area = 0
       DO ie=1,nelemd
          aratio(ie,1) = sum(elem(ie)%mp(:,:)*elem(ie)%metdet(:,:))
       END DO
       CALL ReproSum(aratio, area, nelemd, nelemd, 1, commid=KIM_par%comm)
       area(1) = 4*dd_pi/area(1)  ! ratio correction
       DEALLOCATE(aratio)
       CALL LogMessage("re-initializing cube elements: area correction="//toString(area(1)))

       DO ie=1,nelemd
          CALL cube_init_atomic(elem(ie),gp%points,area(1))
          CALL rotation_init_atomic(elem(ie),rot_type)
       END DO
    END IF
    CALL LogMessage("re-running mass_matrix")
    CALL mass_matrix(KIM_par,elem)

    CALL LogMessage("running global_dof")
    CALL global_dof(KIM_Par,elem)

    CALL LogMessage("create unique indices")
    DO ie=1,nelemd
       CALL CreateUniqueIndex(elem(ie)%GlobalId,elem(ie)%gdofP,elem(ie)%idxP)
    END DO

    CALL LogMessage("apply offset")
    CALL SetElemOffset(KIM_Par,elem, GlobalUniqueCols)

    DO ie=1,nelemd
       elem(ie)%idxV=>elem(ie)%idxP
    END DO

    DEALLOCATE(GridEdge)
    DEALLOCATE(GridVertex)

#ifdef FVM_TRAC
    IF (ntrac > 0) THEN
       CALL fvm_init1(KIM_par)
    END IF
#endif

    ALLOCATE(deriv(0:nThreads-1))
#ifdef FVM_TRAC
    CALL derivinit(deriv(KIM_Hybrid%ithr),fvm_corners,fvm_points)
    IF (ntrac > 0) THEN
       CALL fvm_init2(elem,fvm,KIM_Hybrid,1,nelemd,tl)
    END IF
#else
    CALL derivinit(deriv(KIM_Hybrid%ithr),cslam_corners)
#endif


    hvcoord = hvcoord_init(vfile_mid, vfile_int, .true., KIM_Hybrid%isMasterThr, l_ierr)
    IF (l_ierr /= 0) THEN
       CALL FatalMessage("error in hvcoord_init")
    END IF

    IF (PRESENT(ierr)) ierr = 0

#if defined(CORE_SW)
    DO ie=1,nelemd
       elem(ie)%state%rn2d = 0._r8
       elem(ie)%state%rc2d = 0._r8
       elem(ie)%state%sn2d = 0._r8
    END DO
#endif

    nEPs_l = np*np*nelemd
    nEPs_g = np*np*nelem

    nUPs_l = 0
    DO ie = 1, nelemd
       nUPs_l = nUPs_l + elem(ie)%idxP%NumUniquePts
    ENDDO
    CALL Par_AllReduce(nUPs_l, nUPs_g, 1, PAR_SUM)

    CALL LogMessage('nEPs_l = '//toString(nEPs_l))
    CALL LogMessage('nEPs_g = '//toString(nEPs_g))
    CALL LogMessage('nUPs_l = '//toString(nUPs_l))
    CALL LogMessage('nUPs_g = '//toString(nUPs_g))


    CALL SyncLatLon()




    CALL InfoMessage('End   : Grid initialize..')

    RETURN

    END SUBROUTINE IniGrid

    SUBROUTINE FinGrid(ierr)

    IMPLICIT NONE
    INTEGER(i4), OPTIONAL, INTENT(OUT)        :: ierr
    CALL InfoMessage('Start : Grid finalize..')

    DEALLOCATE(elem)
    DEALLOCATE(deriv)
    IF (PRESENT(ierr)) ierr = 0

    CALL InfoMessage('End   : Grid finalize..')
    RETURN

    END SUBROUTINE FinGrid

    SUBROUTINE SyncLatLon()
    USE Dimensions,  ONLY : np, nlev, nelem, nelemd, nets, nete
    USE Edge,        ONLY: EdgeBuffer_t, InitEdgeBuffer, FreeEdgeBuffer, EdgeVpack, EdgeVunpack, EdgeVunpackMIN
    USE Bndry,       ONLY: Bndry_ExchangeV
    IMPLICIT NONE
    INTEGER(i4)          :: ie, ip, jp, ilev, kptr
    TYPE(EdgeBuffer_t)   :: l_edge2d

      CALL InitEdgeBuffer(l_edge2d, 2)
 
      DO ie = nets, nete
        kptr = 0
        CALL EdgeVpack(l_edge2d, elem(ie)%spherep(:,:)%lat,1,kptr,elem(ie)%desc)
        kptr = 1
        CALL EdgeVpack(l_edge2d, elem(ie)%spherep(:,:)%lon,1,kptr,elem(ie)%desc)
      END DO

      CALL Bndry_ExchangeV(KIM_Par,l_edge2d)


      DO ie = nets, nete
        kptr = 0
        CALL EdgeVunpackMIN(l_edge2d, elem(ie)%spherep(:,:)%lat,1,kptr,elem(ie)%desc)
        kptr = 1
        CALL EdgeVunpackMIN(l_edge2d, elem(ie)%spherep(:,:)%lon,1,kptr,elem(ie)%desc)
      END DO

      CALL FreeEdgeBuffer(l_edge2d)

    END SUBROUTINE SyncLatLon





    SUBROUTINE ConvertEP2UP_Integer(nlevs, EPs, UPs)
    USE Dimensions,     ONLY : np, nlev, nelem, nelemd
    IMPLICIT NONE
    INTEGER(i4), INTENT(IN)                                :: nlevs
    INTEGER(i4), DIMENSION(np,np,nelemd,nlevs), INTENT(IN) :: EPs
    INTEGER(i4), DIMENSION(nUPs_l, nlevs), INTENT(OUT)     :: UPs
    INTEGER(i4)          :: ip, jp, ie, ilev, iup_e, iep, iup

      DO ilev = 1, nlevs
        iup = 0
        DO ie = 1, nelemd
          DO iup_e = 1, elem(ie)%idxP%NumUniquePts
            ip=elem(ie)%idxP%ia(iup_e)
            jp=elem(ie)%idxP%ja(iup_e)
            iup = iup + 1
            UPs(iup,ilev) = EPs(ip,jp,ie,ilev)
          END DO
        END DO
      END DO

    END SUBROUTINE ConvertEP2UP_Integer


    SUBROUTINE ConvertEP2UP_Real(nlevs, EPs, UPs)
    USE Dimensions,     ONLY : np, nlev, nelem, nelemd
    IMPLICIT NONE
    INTEGER(i4), INTENT(IN)                                :: nlevs
    REAL(r4), DIMENSION(np,np,nelemd,nlevs), INTENT(IN)    :: EPs
    REAL(r4), DIMENSION(nUPs_l, nlevs), INTENT(OUT)        :: UPs
    
    INTEGER(i4)          :: ip, jp, ie, ilev, iup_e, iep, iup

      DO ilev = 1, nlevs
        iup = 0
        DO ie = 1, nelemd
          DO iup_e = 1, elem(ie)%idxP%NumUniquePts
            ip=elem(ie)%idxP%ia(iup_e)
            jp=elem(ie)%idxP%ja(iup_e)
            iup = iup + 1
            UPs(iup,ilev) = EPs(ip,jp,ie,ilev)
          END DO
        END DO
      END DO

    END SUBROUTINE ConvertEP2UP_Real


    SUBROUTINE ConvertEP2UP_Double(nlevs, EPs, UPs)
    USE Dimensions,     ONLY : np, nlev, nelem, nelemd
    IMPLICIT NONE
    INTEGER(i4), INTENT(IN)                             :: nlevs
    REAL(r8), DIMENSION(np,np,nelemd,nlevs), INTENT(IN) :: EPs
    REAL(r8), DIMENSION(nUPs_l, nlevs), INTENT(OUT)     :: UPs
    
    INTEGER(i4)          :: ip, jp, ie, ilev, iup_e, iep, iup

      DO ilev = 1, nlevs
        iup = 0
        DO ie = 1, nelemd
          DO iup_e = 1, elem(ie)%idxP%NumUniquePts
            ip=elem(ie)%idxP%ia(iup_e)
            jp=elem(ie)%idxP%ja(iup_e)
            iup = iup + 1
            UPs(iup,ilev) = EPs(ip,jp,ie,ilev)
          END DO
        END DO
      END DO


    END SUBROUTINE ConvertEP2UP_Double






    SUBROUTINE ConvertUP2EP_Double(nlevs, UPs, EPs)
    USE Dimensions,     ONLY: np, nlev, nelem, nelemd
    USE Edge,           ONLY: EdgeBuffer_t, InitEdgeBuffer, FreeEdgeBuffer, EdgeVpack, EdgeVunpack, EdgeVunpackMAX, EdgeVunpackCheck
    USE Bndry,          ONLY: bndry_exchangeV
    IMPLICIT NONE
    INTEGER(i4), INTENT(IN)                              :: nlevs
    REAL(r8), DIMENSION(nUPs_l, nlevs), INTENT(IN)       :: UPs
    REAL(r8), DIMENSION(np,np,nelemd,nlevs), INTENT(OUT) :: EPs
    
    INTEGER(i4)          :: ip, jp, ie, ilev, iup_e, iep, iup
    TYPE(EdgeBuffer_t)   :: edge3d

      CALL InitEdgeBuffer(edge3d, nlevs)

      EPs(:,:,:,:) = -HUGE(1.0_r8)

      DO ilev = 1, nlevs
        iup = 0
        DO ie = 1, nelemd
          DO iup_e = 1, elem(ie)%idxP%NumUniquePts
            ip=elem(ie)%idxP%ia(iup_e)
            jp=elem(ie)%idxP%ja(iup_e)
            iup = iup + 1
            EPs(ip,jp,ie,ilev) = UPs(iup,ilev)
          END DO
        END DO
      END DO


      DO ie = 1, nelemd
        ilev = 0
        CALL EdgeVpack(edge3d, EPs(:,:,ie,:), nlevs, ilev, elem(ie)%desc)
      END DO

      CALL Bndry_ExchangeV(KIM_Par, edge3d)

      DO ie = 1, nelemd
        ilev = 0
        CALL EdgeVunpackMAX(edge3d, EPs(:,:,ie,:), nlevs, ilev, elem(ie)%desc)
      END DO

      CALL FreeEdgeBuffer(edge3d)

    END SUBROUTINE ConvertUP2EP_Double
END MODULE KiapsGrid
