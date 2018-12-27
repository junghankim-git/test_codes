#include <KIM.h>
!
   module kiapsgrid
!
   use kiapsbase, only: i4=>kim_int_kind, l4=>kim_log_kind, r4=>kim_real4_kind, r8=>kim_real8_kind, r16=>kim_longdouble_kind
   use kiapsparallel, only: kim_par
   use dimensions, only: nc
   use element, only: element_t
   use derivative, only: derivative_t
   use hybvcoord, only: hvcoord_t, hvcoord_init
   use logger, only: debugmessage, logmessage, infomessage, warnmessage, errormessage, fatalmessage, tostring
#ifdef FVM_TRAC
   use fvm_control_volume_mod, only: fvm_struct
   use fvm_mod, only: fvm_init1, fvm_init2, fvm_init3
#endif
   implicit none
!
   private
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
   public :: inigrid, fingrid, convertep2up, convertup2ep
!
   interface convertep2up
     module procedure convertep2up_integer
     module procedure convertep2up_real
     module procedure convertep2up_double
   end interface convertep2up
!
   interface convertup2ep
     module procedure convertup2ep_double
   end interface convertup2ep
!
   contains
!
!
!
!
   subroutine inigrid(ierr)
!
   use physicalconstants, only: dd_pi
   use kiapsbase, only: kim_iu_log
   use dimensions, only: np, nlev, nelem, nelemd, nelemdmax, nets, nete, globaluniquecols, ntrac, qsize, nc
   use kiapsparallel, only: iam, parallel_t, barrierpar, abortpar, global_shared_buf, nrepro_vars, par_integer, par_double_precision, par_max, par_sum, finpar
   use kiapsparallel, only: nthreads, omp_get_thread_num, omp_set_num_threads, kim_hybrid
   use kiapsparmpi, only: par_allreduce
   use errorhandle, only: errorhandler, kim_success, kim_error_warning
   use control, only: integration, topology, partmethod, vfile_mid, vfile_int
   use timestep, only: tl, nmax, timelevel_t, time_at, timelevel_init
   use massmatrix, only: mass_matrix
   use cube, only: cubeedgecount, cubeelemcount, cubetopology, cube_init_atomic, rotation_init_atomic, set_corner_coordinates, assign_node_numbers_to_elem
   use metagraph, only: metavertex_t, metaedge_t, localelemcount, initmetagraph
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
   implicit none
!
   integer(i4), intent(  out), optional :: ierr
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
   character(len=80) rot_type ! cube edge rotation type
   integer, dimension(:), allocatable :: tailpartition(:)
   integer, dimension(:), allocatable :: headpartition(:)
   integer :: total_nelem
   real(r8) :: approx_elements_per_task
   logical :: repro_sum_use_ddpdd, repro_sum_recompute
   real(r8) :: repro_sum_rel_diff_max
   integer :: nelem_edge, nedge
   integer :: i, ii, ie, ith, l_ierr
!
   call infomessage('Start:grid initialize..')
   linited = .true.
   total_nelem = cubeelemcount()
   approx_elements_per_task = dble(total_nelem)/dble(kim_par%nprocs)
   if (approx_elements_per_task<1.0d0) then
     call fatalmessage('Number of elements = '//tostring(total_nelem))
     call fatalmessage('Number of processes = '//tostring(kim_par%nprocs))
     call abortpar(message = 'There is not enough parallelism in '//'the job, that is, there is less than '//'one elements per task.')
   endif
   call reprosumdefaultopts(repro_sum_use_ddpdd_out = repro_sum_use_ddpdd, repro_sum_rel_diff_max_out = repro_sum_rel_diff_max, repro_sum_recompute_out = repro_sum_recompute)
   call reprosumsetopts(repro_sum_use_ddpdd_in = repro_sum_use_ddpdd, repro_sum_rel_diff_max_in = repro_sum_rel_diff_max, repro_sum_recompute_in = repro_sum_recompute)
   call logmessage('Initialized repro_sum')
   call logmessage(' ')
   call logmessage('Total simulation time = '//tostring(time_at(nmax)))
   call logmessage(' ')
   if (topology=="cube") then
     call logmessage("creating cube topology:start...")
     nelem = cubeelemcount()
     nelem_edge = cubeedgecount()
     allocate(gridvertex(nelem))
     allocate(gridedge(nelem_edge))
     call cubetopology(gridedge, gridvertex)
     call logmessage("creating cube topology:done...")
   endif
   call logmessage("partitioning graph:start...")
   if (partmethod.eq.sfcurve) then
     call genspacepart(gridedge, gridvertex)
   else
     call fatalmessage("not support metis...")
   endif
   allocate(metavertex(1))
   allocate(schedule(1))
   nelem_edge = size(gridedge)
   allocate(tailpartition(nelem_edge))
   allocate(headpartition(nelem_edge))
   do i = 1, nelem_edge
     tailpartition(i) = gridedge(i)%tail%processor_number
     headpartition(i) = gridedge(i)%head%processor_number
   enddo
   call initmetagraph(iam, metavertex(1), gridvertex, gridedge)
   nelemd = localelemcount(metavertex(1))
   if (nelemd.le.0) then
     call fatalmessage("not yet ready to handle nelemd = 0 yet")
   endif
#ifdef _MPI
   call mpi_allreduce(nelemd, nelemdmax, 1, par_integer, par_max, kim_par%comm, l_ierr)
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
   do i = 1, nc+1
     fvm_corners(i) = 2*(i-1)/dble(nc)-1 ! [-1, 1] including end points
   enddo
   do i = 1, nc
     fvm_points(i) =(fvm_corners(i)+fvm_corners(i+1))/2
   enddo
#endif
   call genedgesched(elem, iam, schedule(1), metavertex(1))
   allocate(global_shared_buf(nelemd, nrepro_vars))
   call barrierpar(par = kim_par)
   gp = gausslobatto(np) ! gll points
   xtmp = nc
   do i = 1, nc+1
     cslam_corners(i) = 2*(i-1)/xtmp-1 ! [-1, 1] including end points
   enddo
   do i = 1, nc
     cslam_points(i) =(cslam_corners(i)+cslam_corners(i+1))/2
   enddo
   if (topology=="cube") then
     call logmessage("initializing cube elements:start...")
     do ie = 1, nelemd
       call set_corner_coordinates(elem(ie))
     enddo
     call assign_node_numbers_to_elem(elem, gridvertex)
     do ie = 1, nelemd
       call cube_init_atomic(elem(ie), gp%points)
     enddo
   endif
   call logmessage("initializing mass_matrix:start...")
   call mass_matrix(kim_par, elem)
   allocate(aratio(nelemd, 1))
   rot_type = "contravariant"
   if (topology=="cube") then
     area = 0
     do ie = 1, nelemd
       aratio(ie, 1) = sum(elem(ie)%mp(:,:)*elem(ie)%metdet(:,:))
     enddo
     call reprosum(aratio, area, nelemd, nelemd, 1, commid = kim_par%comm)
     area(1) = 4*dd_pi/area(1) ! ratio correction
     deallocate(aratio)
     call logmessage("re-initializing cube elements:area correction = "//tostring(area(1)))
     do ie = 1, nelemd
       call cube_init_atomic(elem(ie), gp%points, area(1))
       call rotation_init_atomic(elem(ie), rot_type)
     enddo
   endif
   call logmessage("re-running mass_matrix")
   call mass_matrix(kim_par, elem)
   call logmessage("running global_dof")
   call global_dof(kim_par, elem)
   call logmessage("create unique indices")
   do ie = 1, nelemd
     call createuniqueindex(elem(ie)%globalid, elem(ie)%gdofp, elem(ie)%idxp)
   enddo
   call logmessage("apply offset")
   call setelemoffset(kim_par, elem, globaluniquecols)
   do ie = 1, nelemd
     elem(ie)%idxv=>elem(ie)%idxp
   enddo
   deallocate(gridedge)
   deallocate(gridvertex)
#ifdef FVM_TRAC
   if (ntrac>0) then
     call fvm_init1(kim_par)
   endif
#endif
   allocate(deriv(0:nthreads-1))
#ifdef FVM_TRAC
   call derivinit(deriv(kim_hybrid%ithr), fvm_corners, fvm_points)
   if (ntrac>0) then
     call fvm_init2(elem, fvm, kim_hybrid, 1, nelemd, tl)
   endif
#else
   call derivinit(deriv(kim_hybrid%ithr), cslam_corners)
#endif
   hvcoord = hvcoord_init(vfile_mid, vfile_int, .true., kim_hybrid%ismasterthr, l_ierr)
   if (l_ierr/= 0) then
     call fatalmessage("error in hvcoord_init")
   endif
   if (present(ierr)) ierr = 0
#if defined(CORE_SW)
   do ie = 1, nelemd
     elem(ie)%state%rn2d = 0._r8
     elem(ie)%state%rc2d = 0._r8
     elem(ie)%state%sn2d = 0._r8
   enddo
#endif
   neps_l = np*np*nelemd
   neps_g = np*np*nelem
   nups_l = 0
   do ie = 1, nelemd
     nups_l = nups_l+elem(ie)%idxp%numuniquepts
   enddo
   call par_allreduce(nups_l, nups_g, 1, par_sum)
   call logmessage('nEPs_l = '//tostring(neps_l))
   call logmessage('nEPs_g = '//tostring(neps_g))
   call logmessage('nUPs_l = '//tostring(nups_l))
   call logmessage('nUPs_g = '//tostring(nups_g))
   call synclatlon()
   call infomessage('End:grid initialize..')
   return
!
   end subroutine inigrid
!
!
!
!
   subroutine fingrid(ierr)
!
   implicit none
!
   integer(i4), optional, intent(  out) :: ierr
!
   call infomessage('Start:grid finalize..')
   deallocate(elem)
   deallocate(deriv)
   if (present(ierr)) ierr = 0
   call infomessage('End:grid finalize..')
   return
!
   end subroutine fingrid
!
!
!
!
   subroutine synclatlon()
!
   use dimensions, only: np, nlev, nelem, nelemd, nets, nete
   use edge, only: edgebuffer_t, initedgebuffer, freeedgebuffer, edgevpack, edgevunpack, edgevunpackmin
   use bndry, only: bndry_exchangev
   implicit none
!
   integer(i4) :: ie, ip, jp, ilev, kptr
   type(edgebuffer_t) :: l_edge2d
!
   call initedgebuffer(l_edge2d, 2)
   do ie = nets, nete
     kptr = 0
     call edgevpack(l_edge2d, elem(ie)%spherep(:,:)%lat, 1, kptr, elem(ie)%desc)
     kptr = 1
     call edgevpack(l_edge2d, elem(ie)%spherep(:,:)%lon, 1, kptr, elem(ie)%desc)
   enddo
   call bndry_exchangev(kim_par, l_edge2d)
   do ie = nets, nete
     kptr = 0
     call edgevunpackmin(l_edge2d, elem(ie)%spherep(:,:)%lat, 1, kptr, elem(ie)%desc)
     kptr = 1
     call edgevunpackmin(l_edge2d, elem(ie)%spherep(:,:)%lon, 1, kptr, elem(ie)%desc)
   enddo
   call freeedgebuffer(l_edge2d)
!
   end subroutine synclatlon
!
!
!
!
   subroutine convertep2up_integer(nlevs, eps, ups)
!
   use dimensions, only: np, nlev, nelem, nelemd
   implicit none
!
   integer(i4), intent(in   ) :: nlevs
   integer(i4), dimension(np, np, nelemd, nlevs), intent(in   ) :: eps
   integer(i4), dimension(nups_l, nlevs), intent(  out) :: ups
   integer(i4) :: ip, jp, ie, ilev, iup_e, iep, iup
!
   do ilev = 1, nlevs
     iup = 0
     do ie = 1, nelemd
       do iup_e = 1, elem(ie)%idxp%numuniquepts
         ip = elem(ie)%idxp%ia(iup_e)
         jp = elem(ie)%idxp%ja(iup_e)
         iup = iup+1
         ups(iup, ilev) = eps(ip, jp, ie, ilev)
       enddo
     enddo
   enddo
!
   end subroutine convertep2up_integer
!
!
!
!
   subroutine convertep2up_real(nlevs, eps, ups)
!
   use dimensions, only: np, nlev, nelem, nelemd
   implicit none
!
   integer(i4), intent(in   ) :: nlevs
   real(r4), dimension(np, np, nelemd, nlevs), intent(in   ) :: eps
   real(r4), dimension(nups_l, nlevs), intent(  out) :: ups
   integer(i4) :: ip, jp, ie, ilev, iup_e, iep, iup
!
   do ilev = 1, nlevs
     iup = 0
     do ie = 1, nelemd
       do iup_e = 1, elem(ie)%idxp%numuniquepts
         ip = elem(ie)%idxp%ia(iup_e)
         jp = elem(ie)%idxp%ja(iup_e)
         iup = iup+1
         ups(iup, ilev) = eps(ip, jp, ie, ilev)
       enddo
     enddo
   enddo
!
   end subroutine convertep2up_real
!
!
!
!
   subroutine convertep2up_double(nlevs, eps, ups)
!
   use dimensions, only: np, nlev, nelem, nelemd
   implicit none
!
   integer(i4), intent(in   ) :: nlevs
   real(r8), dimension(np, np, nelemd, nlevs), intent(in   ) :: eps
   real(r8), dimension(nups_l, nlevs), intent(  out) :: ups
   integer(i4) :: ip, jp, ie, ilev, iup_e, iep, iup
!
   do ilev = 1, nlevs
     iup = 0
     do ie = 1, nelemd
       do iup_e = 1, elem(ie)%idxp%numuniquepts
         ip = elem(ie)%idxp%ia(iup_e)
         jp = elem(ie)%idxp%ja(iup_e)
         iup = iup+1
         ups(iup, ilev) = eps(ip, jp, ie, ilev)
       enddo
     enddo
   enddo
!
   end subroutine convertep2up_double
!
!
!
!
   subroutine convertup2ep_double(nlevs, ups, eps)
!
   use dimensions, only: np, nlev, nelem, nelemd
   use edge, only: edgebuffer_t, initedgebuffer, freeedgebuffer, edgevpack, edgevunpack, edgevunpackmax, edgevunpackcheck
   use bndry, only: bndry_exchangev
   implicit none
!
   integer(i4), intent(in   ) :: nlevs
   real(r8), dimension(nups_l, nlevs), intent(in   ) :: ups
   real(r8), dimension(np, np, nelemd, nlevs), intent(  out) :: eps
   integer(i4) :: ip, jp, ie, ilev, iup_e, iep, iup
   type(edgebuffer_t) :: edge3d
!
   call initedgebuffer(edge3d, nlevs)
   eps(:,:,:,:) =-huge(1.0_r8)
   do ilev = 1, nlevs
     iup = 0
     do ie = 1, nelemd
       do iup_e = 1, elem(ie)%idxp%numuniquepts
         ip = elem(ie)%idxp%ia(iup_e)
         jp = elem(ie)%idxp%ja(iup_e)
         iup = iup+1
         eps(ip, jp, ie, ilev) = ups(iup, ilev)
       enddo
     enddo
   enddo
   do ie = 1, nelemd
     ilev = 0
     call edgevpack(edge3d, eps(:,:, ie, :), nlevs, ilev, elem(ie)%desc)
   enddo
   call bndry_exchangev(kim_par, edge3d)
   do ie = 1, nelemd
     ilev = 0
     call edgevunpackmax(edge3d, eps(:,:, ie, :), nlevs, ilev, elem(ie)%desc)
   enddo
   call freeedgebuffer(edge3d)
!
   end subroutine convertup2ep_double
!
!
!
!
   end module kiapsgrid
