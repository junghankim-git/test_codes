#include "KIM.h"
#undef DEBUGPART
!-------------------------------------------------------------------------------
   module scheduler
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
   use metagraph, only: metaedge_t
   use kiapsbase, only: int_kind=>kim_int_kind, log_kind=>kim_log_kind, real_kind=>kim_real8_kind, iulog=>kim_iu_log
   implicit none
!
   private
!
   type, public :: cycle_t
     sequence
     integer :: tag
     integer :: dest
     integer :: source
     integer :: lengthp
     integer :: lengthp_ghost
     integer :: type
     integer :: ptrp
     integer :: ptrp_ghost
     type(metaedge_t), pointer :: edge
   end type cycle_t
!
   type, public :: schedule_t
     sequence
     integer :: ncycles
     integer :: nelemd
     integer :: placeholder ! total integer count should be even
     integer :: nsendcycles
     integer :: nrecvcycles
     integer :: padding
     integer, pointer :: local2global(:)
     type(cycle_t), pointer :: cycle(:)
     type(cycle_t), pointer :: sendcycle(:)
     type(cycle_t), pointer :: recvcycle(:)
     type(cycle_t), pointer :: movecycle(:)
   end type schedule_t
!
   type, public :: graphstats_t
     sequence
     integer(int_kind) :: offnode
     integer(int_kind) :: onnode
     integer(int_kind) :: lb
     integer(int_kind) :: padding
   end type graphstats_t
   type(schedule_t), public, allocatable, target :: schedule(:)
   type(schedule_t), public, allocatable, target :: gschedule(:)
   type(schedule_t), public, allocatable, target :: sschedule(:)
   integer, public, parameter :: bndry_exchange_message = 10
   integer, private, allocatable, target :: global2local(:)
!KJH ADD for MPI
   integer, public, allocatable :: status(:,:)
   integer, public, allocatable :: rrequest(:)
   integer, public, allocatable :: srequest(:)
   integer :: minnelemd, maxnelemd
   public :: genedgesched ! setup the communication schedule for the edge based boundary exchange
   public :: printschedule, printcycle
   public :: checkschedule
   public :: findbufferslot
!  public :: MessageStats
!  public :: PrimMessageStats
!
   contains
!********************** GENCOMSCHED.F ******************************************
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine genedgesched(elem, partnumber, lschedule, metavertex)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use metagraph, only: metavertex_t
   use dimensions, only: nelem
!#ifdef _MPI
!    use parallel_mod, only : nComPoints, iam, PAR_STATUS_SIZE, rrequest, srequest, !	status, npackpoints
!#else
!    use parallel_mod, only : nComPoints, iam
!#endif
   use kiapsparallel, only: ncompoints, npackpoints, iam, par_status_size
   implicit none
!
   type(element_t), intent(inout) :: elem(:)
   integer, intent(in   ) :: partnumber
   type(schedule_t), intent(inout) :: lschedule
   type(metavertex_t), intent(inout) :: metavertex
   integer :: lengthp, total_length, lengthp_ghost
   integer :: i, j, is, ir, ncycle
   integer :: il, ie, ig
   integer :: nelemd0
   integer :: jmd
   integer :: inbr
   integer :: nsched
   integer, allocatable :: tmpp(:,:)
   integer, allocatable :: tmpp_ghost(:,:)
   integer :: nsend, nrecv, nedges
   integer :: icycle
   integer :: isched
   logical, parameter :: verboseprint = .false.
   logical, parameter :: debug = .false.
   integer :: ierr
!
   nsched = size(schedule)
! ================================================
! allocate some arrays for the call to MPI_gatherv
! ================================================
   minnelemd = nelem
   maxnelemd = 0
! =====================================================
! It looks like this is only used in this routine...
! so no need to put it in the schedule data-structure
! =====================================================
   allocate(global2local(nelem))
   if (debug) write(iulog,*) 'genEdgeSched:point #1'
   isched = partnumber
   nelemd0 = metavertex%nmembers
   if (verboseprint) then
     if (iam.eq.1) write(iulog,*) 'genEdgeSched:part # ',i,' has ',nelemd0,' elements '
   endif
   maxnelemd = amax0(maxnelemd,nelemd0)
   minnelemd = amin0(minnelemd,nelemd0)
   if (debug) write(iulog,*) 'genEdgeSched:point #2'
   if (debug) write(iulog,*) 'genEdgeSched:point #3'
   lschedule%ncycles = metavertex%nedges
   lschedule%nelemd = nelemd0
   if (debug) write(iulog,*) 'genEdgeSched:point #4'
!  Note the minus one is for the internal node
   nedges = metavertex%nedges
   if (2*(nedges/2).eq.nedges) then
     nedges = nedges/2
   else
     nedges =(nedges-1)/2
   endif
   lschedule%nsendcycles = nedges
   lschedule%nrecvcycles = nedges
   if (debug) write(iulog,*) 'genEdgeSched:point #5'
! Temporary array to calculate the Buffer Slot
   allocate(tmpp(2,nedges+1))
   allocate(tmpp_ghost(2,nedges+1))
   tmpp(1,:) =-1
   tmpp(2,:) = 0
   tmpp_ghost(1,:) =-1
   tmpp_ghost(2,:) = 0
!  Allocate all the cycle structures
   allocate(lschedule%sendcycle(nedges))
   allocate(lschedule%recvcycle(nedges))
   allocate(lschedule%movecycle(1))
   if (debug) write(iulog,*) 'genEdgeSched:point #6'
!==================================================================
!  Allocate and initalized the index translation arrays
   global2local =-1
   allocate(lschedule%local2global(nelemd0))
   if (debug) write(iulog,*) 'genEdgeSched:point #7'
   do il = 1,nelemd0
     ig = metavertex%members(il)%number
     global2local(ig) = il
     lschedule%local2global(il) = ig
#ifndef _PREDICT
     elem(il)%desc%putmapp =-1
     elem(il)%desc%getmapp =-1
     elem(il)%desc%putmapp_ghost =-1
     elem(il)%desc%getmapp_ghost =-1
     elem(il)%desc%reverse = .false.
#endif
   enddo
!==================================================================
   if (debug) write(iulog,*) 'genEdgeSched:point #8'
   total_length = 0
   ncycle = lschedule%ncycles
   is = 1
   ir = 1
   do j = 1,ncycle
     lengthp = metavertex%edges(j)%wgtp
     lengthp_ghost = metavertex%edges(j)%wgtp_ghost
     if ((metavertex%edges(j)%headvertex==partnumber).and.(metavertex%edges(j)%tailvertex==partnumber)) then
       inbr = partnumber
       if (debug) write(iulog,*) 'genEdgeSched:point #9',iam
       lschedule%movecycle%ptrp = findbufferslot(inbr,lengthp,tmpp)
       lschedule%movecycle%ptrp_ghost = findbufferslot(inbr,lengthp_ghost,tmpp_ghost)
       call setcycle(elem,lschedule,lschedule%movecycle(1),metavertex%edges(j))
       if (debug) write(iulog,*) 'genEdgeSched:point #10',iam
     elseif (metavertex%edges(j)%tailvertex==partnumber) then
       inbr = metavertex%edges(j)%headvertex
       if (debug) write(iulog,*) 'genEdgeSched:point #11',iam
       lschedule%sendcycle(is)%ptrp = findbufferslot(inbr,lengthp,tmpp)
       lschedule%sendcycle(is)%ptrp_ghost = findbufferslot(inbr,lengthp_ghost,tmpp_ghost)
       call setcycle(elem,lschedule,lschedule%sendcycle(is),metavertex%edges(j))
       if (debug) write(iulog,*) 'genEdgeSched:point #12',iam
       is = is+1
     elseif (metavertex%edges(j)%headvertex==partnumber) then
       inbr = metavertex%edges(j)%tailvertex
       if (debug) write(iulog,*) 'genEdgeSched:point #13',iam
       lschedule%recvcycle(ir)%ptrp = findbufferslot(inbr,lengthp,tmpp)
       lschedule%recvcycle(ir)%ptrp_ghost = findbufferslot(inbr,lengthp_ghost,tmpp_ghost)
       call setcycle(elem,lschedule,lschedule%recvcycle(ir),metavertex%edges(j))
       if (debug) write(iulog,*) 'genEdgeSched:point #14',iam
       ir = ir+1
     endif
   enddo
   deallocate(tmpp)
   deallocate(tmpp_ghost)
   do ie = 1,nelemd0
     elem(ie)%vertex = metavertex%members(ie)
     ig = metavertex%members(ie)%number
     elem(ie)%globalid = ig
     elem(ie)%localid = ie
#if 0
     call llinsertedge(eroot,ig,jmd)
!DBG write(iulog,*)'After call to LLInsertEdge in schedule: ie,ig ',ie,ig,jmd
#endif
   enddo
   deallocate(global2local)
!S-JMD call CheckSchedule()
#ifdef _MPI
!================================================================
!     Allocate a couple of structures for bndry_exchange
!        done here to remove it from the critical path
!================================================================
   ncompoints = 0
   nsend = nedges
   nrecv = nedges
   allocate(rrequest(nrecv))
   allocate(srequest(nsend))
   allocate(status(par_status_size,nrecv))
!===============================================================
!   Number of communication points ... to be used later to
!    setup the size of the communication buffer for MPI_Ibsend
!===============================================================
   do icycle = 1,nsend
     ncompoints = ncompoints+lschedule%sendcycle(icycle)%lengthp
   enddo
   npackpoints = ncompoints+lschedule%movecycle(1)%lengthp
!   nbuf = 4*(nv+1)*nelemd*8*4*nlev
!   write(iulog,*)'before call to allocate(combuffer) ',nbuf
!   allocate(combuffer(nbuf))
!   write(iulog,*)'IAM: ',iSched,'Before call to MPI_Buffer_Attach '
!   call MPI_Buffer_Attach(combuffer,nbuf,ierr)
!   write(iulog,*)'IAM: ',iSched,'After call to MPI_Buffer_Attach '
#endif
#ifdef DEBUGPART
   call haltmp("genedgesched:just testing the partitioning algorithms")
#endif
!
   end subroutine genedgesched
#ifdef DOTHIS
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine messagestats(nlyr)
!-------------------------------------------------------------------------------
   use dimensions, only: nmpi_per_node, nnodes, npart, nlev, np, nelem
   use perfmodel, only: perf_t, commtime_t, setptopnetworkparams, & ! _external
                              setserialparamsexplicit, setserialparamsimplicit
   use control, only: integration
!-----------------
   implicit none
!
   integer, intent(in   ) :: nlyr ! number of 2d layers in the communication
   integer :: icycle, ip
   real(real_kind) :: lb_nelemd, lb_volume
   integer :: length, nsend, nrecv
   integer :: i
   real(real_kind), allocatable :: time_total(:), time_calc(:), time_comm(:)
   real(real_kind), allocatable :: time_calc1(:), time_calc2(:)
   type(commtime_t), allocatable :: bndry1(:), bndry2(:), bndry3(:), bndryf(:)
   real(real_kind) :: time_comm1, time_commf, time_comm2, time_comm3
   integer(int_kind) :: bytes_per_point
   real(real_kind) :: time_per_elem
   real(real_kind) :: time_per_iter
   real(real_kind) :: latency_tmp, bandwidth_tmp
   real(real_kind) :: avg_cg_iters
   type(commtime_t) :: offnode, onnode
   type(perf_t) :: tcv
   type(graphstats_t), allocatable :: count(:)
   real(real_kind) :: time_serial, time_parallel, speedup
   integer(int_kind), allocatable :: offnode_count(:), onnode_count(:)
   integer(int_kind), allocatable :: lb_count(:)
   integer(int_kind) :: edgecut_offnode, edgecut_onnode, edgecut_total
   logical(log_kind) :: foundnetwork, foundmachine
   character(len=80) :: networkname, machinename
   integer(int_kind) :: indx_min, indx_comm(1), indx_calc(1)
   integer(8) :: imin_volume, imax_volume
   integer(int_kind), parameter :: configuration = 5
   integer :: node1, node2, nbr, in
   integer :: nelemd0
   logical, parameter :: debug = .false.
   logical, parameter :: predictperformance = .true.
!
   if (predictperformance) then
     networkname = "protobgl"
     call setptopnetworkparams(offnode,onnode,networkname,foundnetwork)
     write(iulog,*) 'MessageStats:after setserialparamsexplicit foundnetwork:',foundnetwork
     machinename = "protobgl"
!JMD integration = "explicit"
     write(iulog,*) 'MessageStats:integration:',integration
     if (integration=="explicit") then
       call setserialparamsexplicit(time_per_elem,np,nlev,machinename,foundmachine)
       write(iulog,*) 'MessageStats:after setserialparamsexplicit foundmachine:',foundmachine
     elseif (integration=="semi_imp") then
       call setserialparamsimplicit(time_per_elem,time_per_iter,avg_cg_iters,np,nelem,nlev,machinename,foundmachine)
     endif
   endif
!
   allocate(time_total(npart))
   allocate(time_calc(npart))
   allocate(time_comm(npart))
   allocate(time_calc1(npart))
   allocate(time_calc2(npart))
   allocate(bndry1(npart))
   allocate(bndry2(npart))
   allocate(bndry3(npart))
   allocate(bndryf(npart))
   write(iulog,*) 'MessageStats:nlev and npart nnodes nmpi_per_node are:',nlev,npart,nnodes,nmpi_per_node
   if (npart.gt.1) then
     allocate(count(nnodes))
     call foo(schedule,count)
!----------------------------------------------------
! Call fooCalc for each boundary exchange
!----------------------------------------------------
     write(iulog,*) 'MessageStats:'
     write(iulog,*) 'MessageStats:boundary exchange #1 '
     bytes_per_point =(3*nlev)*8
     call foocalc(schedule,bytes_per_point,offnode,onnode,count,bndry1)
     write(iulog,*) 'MessageStats:'
     write(iulog,*) 'MessageStats:filter boundary exchange '
     bytes_per_point =(2*nlev)*8
     call foocalc(schedule,bytes_per_point,offnode,onnode,count,bndryf)
     if (integration=="semi_imp") then
       write(iulog,*) 'MessageStats:'
       write(iulog,*) 'MessageStats:boundary exchange #2 '
       bytes_per_point =(2*nlev)*8
       call foocalc(schedule,bytes_per_point,offnode,onnode,count,bndry2)
       write(iulog,*) 'MessageStats:'
       write(iulog,*) 'MessageStats:helmholtz boundary exchange '
       bytes_per_point =(2*nlev)*8
       call foocalc(schedule,bytes_per_point,offnode,onnode,count,bndry3)
     endif
#ifdef _PREDICT
     if (predictperformance) then
!------------------------------------------------------
! Calculate information about execution time prediction
!------------------------------------------------------
       if (integration=="explicit") then
         do ip = 1, npart
           nelemd0 = schedule(ip)%nelemd
           time_calc(ip) = time_per_elem*nelemd0
           if (filter_freq>0) then
             time_commf =(bndryf(ip)%latency+bndryf(ip)%bandwidth)/real(filter_freq,real_kind)
           else
             time_commf = 0
           endif
           time_comm(ip) = time_commf+bndry1(ip)%latency+bndry1(ip)%bandwidth
         enddo
         elseif (integration=="semi_imp") then
           do ip = 1, npart
             nelemd0 = schedule(ip)%nelemd
             time_calc1(ip) = time_per_elem*nelemd0
             time_calc2(ip) = time_per_iter*nelemd0*avg_cg_iters
             time_calc(ip) = time_calc1(ip)+time_calc2(ip)
             if (filter_freq>0) then
               time_commf =(bndryf(ip)%latency+bndryf(ip)%bandwidth)/real(filter_freq,real_kind)
             else
               time_commf = 0
             endif
             time_comm(ip) = time_commf+bndry1(ip)%latency+bndry1(ip)%bandwidth+bndry2(ip)%latency+bndry2(ip)%bandwidth+avg_cg_iters*(bndry3(ip)%latency+bndry3(ip)%bandwidth)
           enddo
           endif
!-----------------------------------------------------------
! Print out information about the execution time prediction
!-----------------------------------------------------------
             indx_calc = maxloc(time_calc)
             indx_comm = maxloc(time_comm)
             time_parallel = time_calc(indx_calc(1))+time_comm(indx_comm(1))
             write(iulog,*) 'MessageStats:'
             write(iulog,*) 'MessageStats:execution time prediction '
             if (foundmachine) write(iulog,*) 'MessageStats:for machine:',trim(machinename)
             if (foundnetwork) write(iulog,*) 'MessageStats:with network:',trim(networkname)
             if (integration=="explicit") then
               time_serial = real(nelem,real_kind)*time_per_elem
               speedup = time_serial/time_parallel
               write(iulog,201) time_serial,1.0d-6*time_serial*(secpday/tstep)
!---------------------------------------
! First Boundary exchange in advance
!---------------------------------------
               if (foundnetwork) write(iulog,65) bndry1(indx_comm(1))%latency+bndry1(indx_comm(1))%bandwidth,bndry1(indx_comm(1))%latency,bndry1(indx_comm(1))%bandwidth
!---------------------------------------
! Filter Boundary exchange in advance
!---------------------------------------
               if (filter_freq>0) then
                 time_commf =(bndryf(indx_comm(1))%latency+bndryf(indx_comm(1))%bandwidth)/real(filter_freq,real_kind)
                 if (foundnetwork) write(iulog,68) time_commf,bndryf(indx_comm(1))%latency,bndryf(indx_comm(1))%bandwidth
               endif
!---------------------------------------
!  Total Communication Cost
!---------------------------------------
               if (foundmachine.and.foundnetwork) write(iulog,55) time_parallel,time_calc(indx_calc(1)),time_comm(indx_comm(1))
             elseif (integration=="semi_imp") then
               time_serial = real(nelem,real_kind)*(time_per_elem+avg_cg_iters*time_per_iter)
               speedup = time_serial/time_parallel
               write(iulog,202) avg_cg_iters
               write(iulog,201) time_serial,1.0d-6*time_serial*(secpday/tstep)
               write(iulog,35) maxval(time_calc1)
               write(iulog,45) maxval(time_calc2)
!---------------------------------------
! Filter Boundary exchange in advance_si
!---------------------------------------
               if (filter_freq>0) then
                 time_commf =(bndryf(indx_comm(1))%latency+bndryf(indx_comm(1))%bandwidth)/real(filter_freq,real_kind)
                 if (foundnetwork) write(iulog,68) time_commf,bndryf(indx_comm(1))%latency,bndryf(indx_comm(1))%bandwidth
               endif
!---------------------------------------
! First Boundary exchange in advance_si
!---------------------------------------
               time_comm1 = bndry1(indx_comm(1))%latency+bndry1(indx_comm(1))%bandwidth
               if (foundnetwork) write(iulog,65) time_comm1,bndry1(indx_comm(1))%latency,bndry1(indx_comm(1))%bandwidth
!---------------------------------------
! Second Boundary exchange in advance_si
!---------------------------------------
               time_comm2 = bndry2(indx_comm(1))%latency+bndry2(indx_comm(1))%bandwidth
               if (foundnetwork) write(iulog,66) time_comm2,bndry2(indx_comm(1))%latency,bndry2(indx_comm(1))%bandwidth
!---------------------------------------
! Boundary exchange in solver_mod
!---------------------------------------
               time_comm3 = bndry3(indx_comm(1))%latency+bndry3(indx_comm(1))%bandwidth
               if (foundnetwork) then
                 write(iulog,67) avg_cg_iters*time_comm3,avg_cg_iters*bndry3(indx_comm(1))%latency,avg_cg_iters*bndry3(indx_comm(1))%bandwidth
               endif
!---------------------------------------
!  Total Communication Cost
!---------------------------------------
               if (foundmachine.and.foundnetwork) write(iulog,55) time_parallel,time_calc(indx_calc(1)),time_comm1+time_commf+time_comm2+avg_cg_iters*time_comm3
             endif
             if (foundmachine.and.foundnetwork) write(iulog,95) speedup
     endif
     call haltmp("messagestats:just testing the partitioning algorithms")
#endif
   endif
   110 format(1x,a,i10,i10,f12.2,f12.2)
   100 format(1x,a,i8,i8,f10.2,f8.5)
   35 format(" messagestats:predicted si step time(usec/step) ",f12.1)
   45 format(" messagestats:predicted pcg solver time(usec/step) ",f12.1)
   55 format(" messagestats:predicted total time(usec/step) total:= ",f12.1," calc:= ",f12.1," comm:= ",f12.1)
   65 format(" messagestats:predicted bndry #1 time(usec/step) total:= ",f12.1," latency:= ",f12.1," band:= ",f12.1)
   165 format(" messagestats:predicted bndry #1 ",i5," time(usec/step) total:= ",f12.1," latency:= ",f12.1," band:= ",f12.1)
   66 format(" messagestats:predicted bndry #2 time(usec/step) total:= ",f12.1," latency:= ",f12.1," band:= ",f12.1)
   67 format(" messagestats:predicted helm bndry time(usec/step) total:= ",f12.1," latency:= ",f12.1," band:= ",f12.1)
   68 format(" messagestats:predicted fltr bndry time(usec/step) total:= ",f12.1," latency:= ",f12.1," band:= ",f12.1)
   75 format(" messagestats:edgecut total:",i10," onnode:",i10," offnode:",i10)
   95 format(" messagestats:speedup:",f10.1)
   201 format(" messagestats:serial execution time(usec/step) ",f12.1,"(sec/day) ",f12.1)
   202 format(" messagestats:average cg iterations:",f6.3)
   120 format(a,i8,i8,f8.5)
!
   end subroutine messagestats
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine primmessagestats(nlyr)
!-------------------------------------------------------------------------------
   use perfmodel, only: perf_t, commtime_t, setptopnetworkparams, & ! _external
                                                    setserialparamsprimexplicit
   use dimensions, only: nmpi_per_node, nnodes, npart, nlev, np, nelem, nelemd
   use control, only: integration
!-----------------
   implicit none
!
   integer, intent(in   ) :: nlyr ! number of 2d layers in the communication
   integer :: icycle, ip
   real(real_kind) :: lb_nelemd, lb_volume
   integer :: length, nsend, nrecv
   integer :: i
   real(real_kind), allocatable :: time_total(:), time_calc(:), time_comm(:)
   real(real_kind), allocatable :: time_calc1(:), time_calc2(:)
   type(commtime_t), allocatable :: bndry1(:), bndry2(:), bndry3(:), bndryf(:)
   real(real_kind) :: time_comm1, time_commf, time_comm2, time_comm3
   integer(int_kind) :: bytes_per_point
   real(real_kind) :: time_per_elem
   real(real_kind) :: time_per_iter
   real(real_kind) :: latency_tmp, bandwidth_tmp
   real(real_kind) :: avg_cg_iters
   type(commtime_t) :: offnode, onnode
   type(perf_t) :: tcv
   type(graphstats_t), allocatable :: count(:)
   real(real_kind) :: time_serial, time_parallel, speedup
   integer(int_kind), allocatable :: offnode_count(:), onnode_count(:)
   integer(int_kind), allocatable :: lb_count(:)
   integer(int_kind) :: edgecut_offnode, edgecut_onnode, edgecut_total
   logical(log_kind) :: foundnetwork, foundmachine
   character(len=80) :: networkname, machinename
   integer(int_kind) :: indx_min, indx_comm(1), indx_calc(1)
   integer(8) :: imin_volume, imax_volume
   integer(int_kind), parameter :: configuration = 5
   integer :: node1, node2, nbr, in
   integer :: nelemd0
   logical, parameter :: debug = .false.
   logical, parameter :: predictperformance = .true.
!
   if (predictperformance) then
     networkname = "protobgl"
     call setptopnetworkparams(offnode,onnode,networkname,foundnetwork)
     write(iulog,*) 'PrimMessageStats:after setserialparamsexplicit foundnetwork:',foundnetwork
     machinename = "protobgl"
     integration = "explicit"
     write(iulog,*) 'PrimMessageStats:integration:',integration
     if (integration=="explicit") then
       call setserialparamsprimexplicit(time_per_elem,np,nlev,machinename,foundmachine)
       write(iulog,*) 'PrimMessageStats:after setserialparamsprimexplicit foundmachine:',foundmachine
     endif
   endif
!
   allocate(time_total(npart))
   allocate(time_calc(npart))
   allocate(time_comm(npart))
   allocate(time_calc1(npart))
   allocate(time_calc2(npart))
   allocate(bndry1(npart))
   allocate(bndry2(npart))
   allocate(bndry3(npart))
   allocate(bndryf(npart))
   write(iulog,*) 'PrimMessageStats:nlev and npart nnodes nmpi_per_node are:',nlev,npart,nnodes,nmpi_per_node
   if (npart.gt.1) then
     allocate(count(nnodes))
     call foo(schedule,count)
!----------------------------------------------------
! Call fooCalc for each boundary exchange
!----------------------------------------------------
     write(iulog,*) 'PrimMessageStats:'
     write(iulog,*) 'PrimMessageStats:boundary exchange #1 '
     bytes_per_point = 3*8
     call foocalc(schedule,bytes_per_point,offnode,onnode,count,bndry1)
     write(iulog,*) 'PrimMessageStats:'
     write(iulog,*) 'PrimMessageStats:boundary exchange #2 '
     bytes_per_point =(3*nlev+1)*8
     call foocalc(schedule,bytes_per_point,offnode,onnode,count,bndry2)
     write(iulog,*) 'PrimMessageStats:'
     write(iulog,*) 'PrimMessageStats:filter boundary exchange '
     bytes_per_point =(3*nlev+1)*8
     call foocalc(schedule,bytes_per_point,offnode,onnode,count,bndryf)
#ifdef _PREDICT
     if (predictperformance) then
!------------------------------------------------------
! Calculate information about execution time prediction
!------------------------------------------------------
       if (integration=="explicit") then
         do ip = 1, npart
           nelemd0 = schedule(ip)%nelemd
           time_calc(ip) = time_per_elem*nelemd0
           if (filter_freq>0) then
             time_commf =(bndryf(ip)%latency+bndryf(ip)%bandwidth)/real(filter_freq,real_kind)
           else
             time_commf = 0
           endif
           time_comm(ip) = time_commf+bndry1(ip)%latency+bndry1(ip)%bandwidth+4.*(bndry2(ip)%latency+bndry2(ip)%bandwidth)
         enddo
         endif
!-----------------------------------------------------------
! Print out information about the execution time prediction
!-----------------------------------------------------------
           indx_calc = maxloc(time_calc)
           indx_comm = maxloc(time_comm)
           time_parallel = time_calc(indx_calc(1))+time_comm(indx_comm(1))
           write(iulog,*) 'PrimMessageStats:'
           write(iulog,*) 'PrimMessageStats:execution time prediction '
           if (foundmachine) write(iulog,*) 'PrimMessageStats:for machine:',trim(machinename)
           if (foundnetwork) write(iulog,*) 'PrimMessageStats:with network:',trim(networkname)
           if (integration=="explicit") then
             time_serial = real(nelem,real_kind)*time_per_elem
             speedup = time_serial/time_parallel
             write(iulog,201) time_serial,1.0d-6*time_serial*(secpday/tstep)
!---------------------------------------
! Small Boundary exchange in prim_advance
!---------------------------------------
             if (foundnetwork) write(iulog,65) bndry1(indx_comm(1))%latency+bndry1(indx_comm(1))%bandwidth,bndry1(indx_comm(1))%latency,bndry1(indx_comm(1))%bandwidth
!-------------------------------------------
! 4 Large Boundary exchanges in prim_advance
!-------------------------------------------
             if (foundnetwork) write(iulog,65) 4.*(bndry2(indx_comm(1))%latency+bndry2(indx_comm(1))%bandwidth),4.*bndry2(indx_comm(1))%latency,4.*bndry1(indx_comm(1))%bandwidth
!---------------------------------------
! Filter Boundary exchange in advance
!---------------------------------------
             if (filter_freq>0) then
               time_commf =(bndryf(indx_comm(1))%latency+bndryf(indx_comm(1))%bandwidth)/real(filter_freq,real_kind)
               if (foundnetwork) write(iulog,68) time_commf,bndryf(indx_comm(1))%latency,bndryf(indx_comm(1))%bandwidth
             endif
!---------------------------------------
!  Total Communication Cost
!---------------------------------------
             if (foundmachine.and.foundnetwork) write(iulog,55) time_parallel,time_calc(indx_calc(1)),time_comm(indx_comm(1))
           elseif (integration=="semi_imp") then
             write(iulog,*) 'PrimMessageStats:not yet support for semi-implicit time integration'
           endif
           if (foundmachine.and.foundnetwork) write(iulog,95) speedup
     endif
     call haltmp("primmessagestats:just testing the partitioning algorithms")
#endif
   endif
   110 format(1x,a,i10,i10,f12.2,f12.2)
   100 format(1x,a,i8,i8,f10.2,f8.5)
   35 format(" primmessagestats:predicted si step time(usec/step) ",f12.1)
   45 format(" primmessagestats:predicted pcg solver time(usec/step) ",f12.1)
   55 format(" primmessagestats:predicted total time(usec/step) total:= ",f12.1," calc:= ",f12.1," comm:= ",f12.1)
   65 format(" primmessagestats:predicted bndry #1 time(usec/step) total:= ",f12.1," latency:= ",f12.1," band:= ",f12.1)
   165 format(" primmessagestats:predicted bndry #1 ",i5," time(usec/step) total:= ",f12.1," latency:= ",f12.1," band:= ",f12.1)
   66 format(" primmessagestats:predicted bndry #2 time(usec/step) total:= ",f12.1," latency:= ",f12.1," band:= ",f12.1)
   67 format(" primmessagestats:predicted helm bndry time(usec/step) total:= ",f12.1," latency:= ",f12.1," band:= ",f12.1)
   68 format(" primmessagestats:predicted fltr bndry time(usec/step) total:= ",f12.1," latency:= ",f12.1," band:= ",f12.1)
   75 format(" primmessagestats:edgecut total:",i10," onnode:",i10," offnode:",i10)
   95 format(" primmessagestats:speedup:",f10.1)
   201 format(" primmessagestats:serial execution time(usec/step) ",f12.1,"(sec/day) ",f12.1)
   202 format(" primmessagestats:average cg iterations:",f6.3)
   120 format(a,i8,i8,f8.5)
!
   end subroutine primmessagestats
#endif
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine checkschedule()
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: i, nsched, nbufferwords_1, nbufferwords_2
   type(schedule_t), pointer :: pschedule
!
   nsched = size(schedule)
   do i = 1,nsched
     pschedule=>schedule(i)
     nbufferwords_1 = sum(pschedule%sendcycle(:)%lengthp)
     nbufferwords_2 = sum(pschedule%recvcycle(:)%lengthp)
     if (nbufferwords_1.ne.nbufferwords_2) then
       write(*,100) i,nbufferwords_1,nbufferwords_2
     endif
   enddo
   100 format('CheckSchedule:err iam:',i3,' sizeof(sendbuffer):',i10,' ! = sizeof(recvbuffer):',i10)
   110 format('CheckSchedule:err iam:',i3,' length(sendbuffer):',i10,' ! = length(recvbuffer):',i10)
!
   end subroutine checkschedule
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine printschedule(schedule)
!-------------------------------------------------------------------------------
   use gridgraph, only: printgridedge
   implicit none
!
   type(schedule_t), intent(in   ), target :: schedule(:)
   type(schedule_t), pointer :: pschedule
   type(cycle_t), pointer :: pcycle
   integer :: i, j, nsched
!
   nsched = size(schedule)
   write(*,*) '------new schedule format---------------------'
   do i = 1,nsched
     pschedule=>schedule(i)
     write(*,*)
     write(*,*) '----------------------------------------------'
     write(*,90) i,pschedule%ncycles
     write(*,*) '----------------------------------------------'
     write(*,*) '-----------send-------------------------------'
     do j = 1,pschedule%nsendcycles
       pcycle=>pschedule%sendcycle(j)
       call printcycle(pcycle)
       call printgridedge(pcycle%edge%members)
     enddo
     write(*,*) '-----------recv-------------------------------'
     do j = 1,pschedule%nrecvcycles
       pcycle=>pschedule%recvcycle(j)
       call printcycle(pcycle)
       call printgridedge(pcycle%edge%members)
     enddo
     write(*,*) '-----------move-------------------------------'
     pcycle=>pschedule%movecycle(1)
     call printcycle(pcycle)
     call printgridedge(pcycle%edge%members)
   enddo
   90 format('NODE # ',i2,2x,'NCYCLES ',i2)
   97 format(10x,'EDGE #',i2,2x,'TYPE ',i1,2x,'G.EDGES',i4,2x,'WORDS ',i5,2x,'SRC ',i3,2x,'DEST ',i3,2x,'PTR ',i4)
   100 format(15x,i4,5x,i3,1x,'(',i1,')--',i1,'-->',i3,1x,'(',i1,') ')
!
   end subroutine printschedule
#ifdef DOTHIS
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine foocalc(schedule, bytes_per_point, offnode, onnode, count, time)
!-------------------------------------------------------------------------------
   use perfmodel, only: commtime_t, perf_t ! _external
   use flops_mod, only: min_number, max_number, avg_number, & ! _external
min_volume, max_volume, avg_volume, min_message, max_message, avg_message, localcomvolume, totalcomvolume, tnsend
   use kiapsparallel, only: iam, ncompoints
   use dimensions, only: nmpi_per_node, npart
! results:: TnSend,{min,max}_{volume,number},Time_{total,calc,latency,bandwidth}
! inputs:: results + Schedule + bytes_per_point
!
   type(schedule_t), intent(in   ) :: schedule(:)
   integer(int_kind), intent(in   ) :: bytes_per_point
   type(commtime_t), intent(in   ) :: offnode, onnode
   type(graphstats_t), intent(in   ) :: count(:)
   type(commtime_t), intent(inout) :: time(:)
   real(real_kind) :: bandwidth_tmp, latency_tmp
   integer(int_kind) :: nbr, node1, node2, nsend, ip, icycle, nrecv
   integer(int_kind) :: length
   real(real_kind) :: lb_volume
   integer(int_kind) :: imax_volume, imin_volume
   type(perf_t) :: tcv
   integer :: nsched
   integer :: nelemd0
   logical, parameter :: debug = .false.
!
   nsched = size(schedule)
   min_message = 1000000
   min_volume = 1000000
   max_message = 0
   max_volume = 0
   totalcomvolume = 0
   tnsend = 0
   tcv%onnode = 0.0
   tcv%offnode = 0.0
   do ip = 1,nsched
     nsend = schedule(ip)%nsendcycles
     nrecv = schedule(ip)%nrecvcycles
     nelemd0 = schedule(ip)%nelemd
     ncompoints = 0
     time(ip)%latency = 0
     time(ip)%bandwidth = 0
     do icycle = 1,nsend
       length = schedule(ip)%sendcycle(icycle)%lengthv
       ncompoints = ncompoints+length
       min_message = min(min_message,bytes_per_point*length)
       max_message = max(max_message,bytes_per_point*length)
     enddo
     do icycle = 1,nrecv
       length = schedule(ip)%recvcycle(icycle)%lengthv
       nbr = schedule(ip)%recvcycle(icycle)%source
       node1 =((ip-1)/nmpi_per_node)+1
       node2 =((nbr-1)/nmpi_per_node)+1
       if (node1.eq.node2) then
! Message is set to an on-node processor
         tcv%onnode = tcv%onnode+length*bytes_per_point
         bandwidth_tmp = onnode%bandwidth
         latency_tmp = onnode%latency
       else
! Message is set to an off-node processor
         tcv%offnode = tcv%offnode+length*bytes_per_point
         bandwidth_tmp = real(count(node1)%offnode,real_kind)*offnode%bandwidth
         latency_tmp = offnode%latency
         if (debug) write(iulog,*) 'fooCalc:node # ',node1,' offnode count is ',real(count(node1)%offnode,real_kind)
         if (debug) write(iulog,*) 'fooCalc:node # ',node1,' offnode = {total,effective} ',offnode%bandwidth,bandwidth_tmp
       endif
       time(ip)%latency = time(ip)%latency+latency_tmp
       if (debug) then
         write(iulog,*) 'fooCalc:node # ',node1,' contribution in time ',length*bytes_per_point,length*bytes_per_point*bandwidth_tmp
       endif
       time(ip)%bandwidth = time(ip)%bandwidth+length*bytes_per_point*bandwidth_tmp
     enddo
!JMD  Note: I am removing the size of the communication buffer out of the computation
     localcomvolume = bytes_per_point*ncompoints
     min_volume = min(min_volume,localcomvolume)
     max_volume = max(max_volume,localcomvolume)
     tnsend = tnsend+nsend
     totalcomvolume = totalcomvolume+localcomvolume
   enddo
!------------------------------------------------------
! This Outputs information about each boundary exchange
!------------------------------------------------------
   avg_volume = dble(totalcomvolume)/dble(npart)
   lb_volume =(max_volume-avg_volume)/max_volume
   avg_message = dble(totalcomvolume)/dble(tnsend)
   imin_volume = nint(min_volume)
   imax_volume = nint(max_volume)
   if (iam.eq.1) then
     write(iulog,*) 'MessageStats:total message volume is(mbytes):',1.0e-6*totalcomvolume
#if 1
     write(iulog,*) 'MessageStats:imin_volume',imin_volume
     write(iulog,*) 'MessageStats:imax_volume',imax_volume
     write(iulog,*) 'MessageStats:avg_volume',avg_volume
     write(iulog,*) 'MessageStats:lb_volume',lb_volume
#else
!JMD For some reason on AIX this line breaks when the numbers are high
     write(iulog,110) 'MessageStats:single process volume(bytes) {min,max,avg,lb} =:',imin_volume,imax_volume,avg_volume,lb_volume
#endif
     write(iulog,*) 'MessageStats:single message(bytes) {min,max,avg} =:',min_message,max_message,avg_message
     write(iulog,85) 1.0e-3*totalcomvolume,1.0e-3*tcv%onnode,1.0e-3*tcv%offnode
   endif
   85 format(" messagestats:tcv total:",f10.1," onnode:",f10.1," offnode:",f10.1)
   110 format(1x,a,i10,i10,f12.2,f12.2)
!
   end subroutine foocalc
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine foo(schedule, count)
!-------------------------------------------------------------------------------
   use flops_mod, only: min_number, max_number, avg_number ! _external
   use dimensions, only: nmpi_per_node, npart, nelem
   use kiapsparallel, only: iam
!
   type(schedule_t), intent(in   ) :: schedule(:)
   type(graphstats_t), intent(inout) :: count(:)
   integer :: node2, nbr, ip, node1, icycle, nrecv, nsend
   logical :: offnode, found_offnode, found_onnode, found_lbnode
   integer(int_kind) :: edgecut_total
   real(real_kind) :: lb_nelemd, avg_nelemd
   integer(int_kind) :: edgecut_offnode, edgecut_onnode
   integer(int_kind) :: tnsend
   integer :: nsched
   integer :: nelemd0
   logical(log_kind), parameter :: debug = .false.
!
   nsched = size(schedule)
   count(:)%offnode = 0
   count(:)%onnode = 0
   count(:)%lb = 0
   edgecut_offnode = 0
   edgecut_onnode = 0
   min_number = 1000000
   max_number = 0
   tnsend = 0
   do ip = 1,nsched
     node1 =((ip-1)/nmpi_per_node)+1
     found_offnode = .false.
     found_onnode = .false.
     found_lbnode = .false.
     nrecv = schedule(ip)%nrecvcycles
     nsend = schedule(ip)%nsendcycles
     nelemd0 = schedule(ip)%nelemd
     do icycle = 1,nrecv
       nbr = schedule(ip)%recvcycle(icycle)%source
       node2 =((nbr-1)/nmpi_per_node)+1
       if (node2.ne.node1) then
! Message is set off-node
         if (nelemd0.eq.maxnelemd) found_lbnode = .true.
         found_offnode = .true.
         edgecut_offnode = edgecut_offnode+1
       else
! Message is set to an on-node processor
         edgecut_onnode = edgecut_onnode+1
       endif
       if (node2.eq.node1) found_onnode = .true.
     enddo
     if (found_offnode) count(node1)%offnode = count(node1)%offnode+1
     if (found_onnode) count(node1)%onnode = count(node1)%onnode+1
     if (found_lbnode) count(node1)%lb = count(node1)%lb+1
     min_number = min(min_number,nsend)
     max_number = max(max_number,nsend)
     tnsend = tnsend+nsend
   enddo
   if (debug) write(iulog,*) 'count(:)%offnode:',count(:)%offnode
   if (debug) write(iulog,*) 'count(:)%lb:',count(:)%lb
!-------------------------------------------------
! This Outputs general information about the Graph
!-------------------------------------------------
   edgecut_offnode = edgecut_offnode/2
   edgecut_onnode = edgecut_onnode/2
   edgecut_total = edgecut_offnode+edgecut_onnode
   avg_nelemd = dble(nelem)/dble(npart)
   lb_nelemd =(maxnelemd-avg_nelemd)/maxnelemd
   avg_number = dble(tnsend)/dble(npart)
   if (iam.eq.1) then
     write(iulog,*) 'MessageStats:number of mpi processes are:',npart
     write(iulog,100) 'MessageStats:nelemd is {min,max,avg,lb}:',minnelemd,maxnelemd,avg_nelemd,lb_nelemd
     write(iulog,*) 'MessageStats:number of neighbors {min,max,avg} ',min_number,max_number,avg_number
     write(iulog,75) edgecut_total,edgecut_onnode,edgecut_offnode
   endif
   100 format(1x,a,i8,i8,f10.2,f8.5)
   75 format(" messagestats:edgecut total:",i10," onnode:",i10," offnode:",i10)
!
   end subroutine foo
#endif
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine printcycle(cycle)
!-------------------------------------------------------------------------------
   implicit none
!
   type(cycle_t), intent(in   ), target :: cycle
!
   write(*,97) cycle%edge%number,cycle%type,cycle%edge%nmembers,cycle%lengthp,cycle%source,cycle%dest,cycle%ptrp
   97 format(5x,'METAEDGE #',i2,2x,'TYPE ',i1,2x,'G.EDGES',i4,2x,'WORDS ',i5,2x,'SRC ',i3,2x,'DEST ',i3,2x,'PTR ',i5)
!
   end subroutine printcycle
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine setcycle(elem, schedule, cycle, edge)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use dimensions, only: max_corner_elem, max_neigh_edges
   use kiapsparallel, only: abortpar, iam
   implicit none
!
   type(element_t), intent(inout) :: elem(:)
   type(schedule_t), intent(inout) :: schedule
   type(cycle_t), intent(inout) :: cycle
   type(metaedge_t), intent(in   ), target :: edge
   integer :: i, il, face
!  allocate and initalize the index into the buffer slots
!JMD    allocate(Cycle%sEdgeIndex(Edge%nmembers))
!JMD    allocate(Cycle%rEdgeIndex(Edge%nmembers))
#ifndef _PREDICT
!
   do i = 1, edge%nmembers
!   Setup send index
     il = global2local(edge%members(i)%tail%number)
     face = edge%members(i)%tail_face
     if (face.ge.5) then
       face = face+(face-5)*(max_corner_elem-1)+(edge%members(i)%tail_ind-1)
     endif
     if (il.gt.0) then
!          print*, "i = ", i, "face = ", face
!          print*, "putmapV = ", elem(il)%desc%putmapV(face)
!          print*, '-------------'
       elem(il)%desc%putmapp(face) = edge%edgeptrp(i)+cycle%ptrp-1 ! offset,so start at 0
       elem(il)%desc%putmapp_ghost(face) = edge%edgeptrp_ghost(i)+cycle%ptrp_ghost ! index,start at 1
       elem(il)%desc%reverse(face) = edge%members(i)%reverse
     endif
!   Setup receive index
     il = global2local(edge%members(i)%head%number)
     face = edge%members(i)%head_face
     if (face.ge.5) then
       face = face+(face-5)*(max_corner_elem-1)+(edge%members(i)%head_ind-1)
       if (face>max_neigh_edges) then
         print*,__file__,__line__,iam,face,i,max_corner_elem,edge%members(i)%head_ind
         call abortpar(message = 'Face value out of bounds')
       endif
     endif
     if (il.gt.0) then
       elem(il)%desc%getmapp(face) = edge%edgeptrp(i)+cycle%ptrp-1
       elem(il)%desc%getmapp_ghost(face) = edge%edgeptrp_ghost(i)+cycle%ptrp_ghost
     endif
   enddo
#endif
!
   cycle%edge=>edge
   cycle%type = edge%type
   cycle%dest = edge%headvertex
   cycle%source = edge%tailvertex
   cycle%tag = bndry_exchange_message
   cycle%lengthp = edge%wgtp
   cycle%lengthp_ghost = edge%wgtp_ghost
!
   end subroutine setcycle
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function findbufferslot(inbr, length, tmp) result(ptr)
!
   integer :: ptr
   integer, intent(in   ) :: inbr, length
   integer, intent(inout) :: tmp(:,:)
   integer :: i, n
!
   n = size(tmp,2)
   ptr = 0
   do i = 1,n
     if (tmp(1, i)==inbr) then
       ptr = tmp(2,i)
       return
     endif
       if (tmp(1, i)==-1) then
         tmp(1,i) = inbr
         if (i.eq.1) tmp(2,i) = 1
         ptr = tmp(2,i)
         if (i.ne.n) tmp(2,i+1) = ptr+length
         return
       endif
   enddo
!
   end function findbufferslot
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module scheduler
!-------------------------------------------------------------------------------
