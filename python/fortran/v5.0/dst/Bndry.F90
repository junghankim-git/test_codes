!-------------------------------------------------------------------------------
!>
!> @brief
!>  - bndry module.
!>
!> @date ?????2012
!>  - JH KIM  : First written from the HOMME and was modified for KIAPSGM framework.
!> @date 30JAN2015
!>  - JH KIM  : Added to "ghost_exchangeV" subroutine. (CAM-SE 5.3)
!-------------------------------------------------------------------------------
#include "KIM.h"
!-------------------------------------------------------------------------------
   module bndry
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
   use kiapsbase, only: int_kind=>kim_int_kind, real_kind=>kim_real8_kind, log_kind=>kim_log_kind
!USE KiapsParallel, ONLY : iam
   use kiapsparallel
   implicit none
!
   private
   public :: bndry_exchangev, ghost_exchangevfull, compute_ghost_corner_orientation
   public :: ghost_exchangev, ghost_exchangev3d
   public :: ghost_exchangev_new
!
   interface bndry_exchangev
     module procedure bndry_exchangev_nonth
     module procedure long_bndry_exchangev_nonth
     module procedure bndry_exchangev_thsave
   end interface
!
   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bndry_exchangev_nonth(par, buffer)
!-------------------------------------------------------------------------------
   use edge, only: edgebuffer_t
   use scheduler, only: schedule_t, cycle_t, schedule, status, srequest, rrequest, status, srequest, rrequest
   use kiapsparallel, only: parallel_t, abortpar, omp_in_parallel, par_success, par_status_size, par_double_precision, par_real8
   use kiapsparmpi, only: par_isend, par_irecv, par_allwait, par_wait, nreqs
!
   type(parallel_t) :: par
   type(edgebuffer_t) :: buffer
   type(schedule_t), pointer :: pschedule
   type(cycle_t), pointer :: pcycle
   integer :: dest, length, tag
   integer :: icycle, ierr
   integer :: iptr, source, nlyr
   integer :: nsendcycles, nrecvcycles
   integer :: errorcode, errorlen
   character*(80) errorstring
   logical(log_kind), parameter :: debug = .false.
   integer :: i
#ifdef _MPI
!
   if (omp_in_parallel()) then
     print*,'bndry_exchangeV:warning you are calling a non-thread safe'
     print*,' routine inside a threaded region.... '
     print*,' results are not predictable!! '
   endif
! Setup the pointer to proper Schedule
#ifdef _PREDICT
!
   pschedule=>schedule(iam)
#else
   pschedule=>schedule(1)
#endif
   nlyr = buffer%nlyr
   nsendcycles = pschedule%nsendcycles
   nrecvcycles = pschedule%nrecvcycles
!==================================================
!  Fire off the sends
!==================================================
   do icycle = 1,nsendcycles
     pcycle=>pschedule%sendcycle(icycle)
     dest = pcycle%dest-1
     length = nlyr*pcycle%lengthp
     tag = pcycle%tag
     iptr = pcycle%ptrp
!DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
     call par_isend(buffer%buf(1,iptr),length,dest,ierr = ierr)
!       call MPI_Isend(buffer%buf(1,iptr),length,PAR_REAL8,dest,tag,par%comm,Srequest(icycle),ierr)
     if (ierr.ne.par_success) then
       errorcode = ierr
       call mpi_error_string(errorcode,errorstring,errorlen,ierr)
       print*,'bndry_exchangeV:error after call to mpi_isend:',errorstring
     endif
   enddo ! icycle
!==================================================
!  Post the Receives
!==================================================
   do icycle = 1,nrecvcycles
     pcycle=>pschedule%recvcycle(icycle)
     source = pcycle%source-1
     length = nlyr*pcycle%lengthp
     tag = pcycle%tag
     iptr = pcycle%ptrp
!DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
     call par_irecv(buffer%receive(1,iptr),length,source,ierr = ierr)
!       call MPI_Irecv(buffer%receive(1,iptr),length,PAR_REAL8,&
!            source,tag,par%comm,Rrequest(icycle),ierr)
     if (ierr.ne.par_success) then
       errorcode = ierr
       call mpi_error_string(errorcode,errorstring,errorlen,ierr)
       print*,'bndry_exchangeV:error after call to mpi_irecv:',errorstring
     endif
   enddo ! icycle
!==================================================
!  Wait for all the receives to complete
!==================================================
   call par_allwait()
!    CALL Par_Wait()
!    call MPI_Waitall(nSendCycles,Srequest,status,ierr)
!    call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
   do icycle = 1,nrecvcycles
     pcycle=>pschedule%recvcycle(icycle)
     length = pcycle%lengthp
     iptr = pcycle%ptrp
     do i = 0,length-1
       buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
     enddo
   enddo ! icycle
#endif
!
   end subroutine bndry_exchangev_nonth
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine long_bndry_exchangev_nonth(par, buffer)
!-------------------------------------------------------------------------------
   use edge, only: longedgebuffer_t
   use scheduler, only: schedule_t, cycle_t, schedule, status, srequest, rrequest, status, srequest, rrequest
   use kiapsparallel, only: parallel_t, abortpar, omp_in_parallel, par_success, par_status_size, par_double_precision, par_integer
   use kiapsparmpi, only: par_isend, par_irecv, par_allwait
!
   type(parallel_t) :: par
   type(longedgebuffer_t) :: buffer
   type(schedule_t), pointer :: pschedule
   type(cycle_t), pointer :: pcycle
   integer :: dest, length, tag
   integer :: icycle, ierr
   integer :: iptr, source, nlyr
   integer :: nsendcycles, nrecvcycles
   integer :: errorcode, errorlen
   character*(80) errorstring
   logical(log_kind), parameter :: debug = .false.
   integer :: i
#ifdef _MPI
!
   if (omp_in_parallel()) then
     print*,'bndry_exchangeV:warning you are calling a non-thread safe'
     print*,' routine inside a threaded region.... '
     print*,' results are not predictable!! '
   endif
! Setup the pointer to proper Schedule
#ifdef _PREDICT
!
   pschedule=>schedule(iam)
#else
   pschedule=>schedule(1)
#endif
   nlyr = buffer%nlyr
   nsendcycles = pschedule%nsendcycles
   nrecvcycles = pschedule%nrecvcycles
!==================================================
!  Fire off the sends
!==================================================
   do icycle = 1,nsendcycles
     pcycle=>pschedule%sendcycle(icycle)
     dest = pcycle%dest-1
     length = nlyr*pcycle%lengthp
     tag = pcycle%tag
     iptr = pcycle%ptrp
!DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
     call par_isend(buffer%buf(1,iptr),length,dest,ierr = ierr)
     if (ierr.ne.par_success) then
       errorcode = ierr
       call mpi_error_string(errorcode,errorstring,errorlen,ierr)
       print*,'bndry_exchangeV:error after call to mpi_isend:',errorstring
     endif
   enddo ! icycle
!==================================================
!  Post the Receives
!==================================================
   do icycle = 1,nrecvcycles
     pcycle=>pschedule%recvcycle(icycle)
     source = pcycle%source-1
     length = nlyr*pcycle%lengthp
     tag = pcycle%tag
     iptr = pcycle%ptrp
!DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
     call par_irecv(buffer%receive(1,iptr),length,source,ierr = ierr)
     if (ierr.ne.par_success) then
       errorcode = ierr
       call mpi_error_string(errorcode,errorstring,errorlen,ierr)
       print*,'bndry_exchangeV:error after call to mpi_irecv:',errorstring
     endif
   enddo ! icycle
!==================================================
!  Wait for all the receives to complete
!==================================================
   call par_allwait()
   do icycle = 1,nrecvcycles
     pcycle=>pschedule%recvcycle(icycle)
     length = pcycle%lengthp
     iptr = pcycle%ptrp
     do i = 0,length-1
       buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
     enddo
   enddo ! icycle
#endif
!
   end subroutine long_bndry_exchangev_nonth
!********************************************************************************
!
!********************************************************************************
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine bndry_exchangev_thsave(hybrid, buffer)
!    use hybrid_mod, only : hybrid_t
!-------------------------------------------------------------------------------
   use edge, only: edgebuffer_t
   use scheduler, only: schedule_t, cycle_t, schedule, status, srequest, rrequest, status, srequest, rrequest
   use dimensions, only: nelemd, np
!+isong
!   use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
!-isong
   use timing, only: timingstart, timingstop
   use kiapsparallel, only: parallel_t, hybrid_t, abortpar, omp_in_parallel, par_success, par_status_size, par_double_precision
   use kiapsparmpi, only: par_isend, par_irecv, par_allwait
   implicit none
!
   type(hybrid_t) :: hybrid
   type(edgebuffer_t) :: buffer
   type(schedule_t), pointer :: pschedule
   type(cycle_t), pointer :: pcycle
   integer :: dest, length, tag
   integer :: icycle, ierr
   integer :: iptr, source, nlyr
   integer :: nsendcycles, nrecvcycles
   integer :: errorcode, errorlen
   character*(80) errorstring
   integer :: i, l_ierr
   logical(log_kind), parameter :: debug = .false.
!$OMP MASTER
!+isong
!   call t_startf('bndry_exchange')
!-isong
!    l_ierr = TimingStart('bndry_exchange')
!$OMP END MASTER
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
!
   if (hybrid%ithr==0) then
#ifdef _MPI
! Setup the pointer to proper Schedule
#ifdef _PREDICT
     pschedule=>schedule(iam)
#else
     pschedule=>schedule(1)
#endif
     nlyr = buffer%nlyr
     nsendcycles = pschedule%nsendcycles
     nrecvcycles = pschedule%nrecvcycles
!==================================================
!  Fire off the sends
!==================================================
     do icycle = 1,nsendcycles
       pcycle=>pschedule%sendcycle(icycle)
       dest = pcycle%dest-1
       length = nlyr*pcycle%lengthp
       tag = pcycle%tag
       iptr = pcycle%ptrp
!DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call par_isend(buffer%buf(1,iptr),length,dest,ierr = ierr)
       if (ierr.ne.par_success) then
         errorcode = ierr
         call mpi_error_string(errorcode,errorstring,errorlen,ierr)
         print*,'bndry_exchangeV:error after call to mpi_isend:',errorstring
       endif
     enddo ! icycle
!==================================================
!  Post the Receives
!==================================================
     do icycle = 1,nrecvcycles
       pcycle=>pschedule%recvcycle(icycle)
       source = pcycle%source-1
       length = nlyr*pcycle%lengthp
       tag = pcycle%tag
       iptr = pcycle%ptrp
!DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call par_irecv(buffer%receive(1,iptr),length,source,ierr = ierr)
       if (ierr.ne.par_success) then
         errorcode = ierr
         call mpi_error_string(errorcode,errorstring,errorlen,ierr)
         print*,'bndry_exchangeV:error after call to mpi_irecv:',errorstring
       endif
     enddo ! icycle
!==================================================
!  Wait for all the receives to complete
!==================================================
     call par_allwait()
     do icycle = 1,nrecvcycles
       pcycle=>pschedule%recvcycle(icycle)
       length = pcycle%lengthp
       iptr = pcycle%ptrp
       do i = 0,length-1
         buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
       enddo
     enddo ! icycle
#endif
   endif ! if (hybrid%ithr==0)
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
!$OMP MASTER
!+isong
!   call t_stopf('bndry_exchange')
!-isong
!    l_ierr = TimingStop('bndry_exchange')
!$OMP END MASTER
!
   end subroutine bndry_exchangev_thsave
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ghost_exchangevfull(hybrid, buffer, nc)
!
!   MT 2011:  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
!
!    use hybrid_mod, only : hybrid_t
!-------------------------------------------------------------------------------
   use edge, only: ghostbuffer_t
   use scheduler, only: schedule_t, cycle_t, schedule, status, srequest, rrequest, status, srequest, rrequest
   use dimensions, only: nelemd
   use kiapsparallel, only: parallel_t, hybrid_t, abortpar, omp_in_parallel, par_success, par_status_size, par_double_precision
   use kiapsparmpi, only: par_isend, par_irecv, par_allwait
   implicit none
!
   type(hybrid_t) :: hybrid
   type(ghostbuffer_t) :: buffer
   integer :: nc
   type(schedule_t), pointer :: pschedule
   type(cycle_t), pointer :: pcycle
   integer :: dest, length, tag
   integer :: icycle, ierr
   integer :: iptr, source, nlyr
   integer :: nsendcycles, nrecvcycles
   integer :: errorcode, errorlen
   character*(80) errorstring
   integer :: i, i1, i2
   logical(log_kind), parameter :: debug = .false.
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
!
   if (hybrid%ithr==0) then
#ifdef _MPI
! Setup the pointer to proper Schedule
#ifdef _PREDICT
     pschedule=>schedule(iam)
#else
     pschedule=>schedule(1)
#endif
     nlyr = buffer%nlyr
     nsendcycles = pschedule%nsendcycles
     nrecvcycles = pschedule%nrecvcycles
!==================================================
!  Fire off the sends
!==================================================
     do icycle = 1,nsendcycles
       pcycle=>pschedule%sendcycle(icycle)
       dest = pcycle%dest-1
       length = nlyr*pcycle%lengthp_ghost*nc*nc
       tag = pcycle%tag
       iptr = pcycle%ptrp_ghost
!print *,'ghost_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call par_isend(buffer%buf(1,1,1,iptr),length,dest,ierr = ierr)
       if (ierr.ne.par_success) then
         errorcode = ierr
         call mpi_error_string(errorcode,errorstring,errorlen,ierr)
         print*,'bndry_exchangeV:error after call to mpi_isend:',errorstring
       endif
     enddo ! icycle
!==================================================
!  Post the Receives
!==================================================
     do icycle = 1,nrecvcycles
       pcycle=>pschedule%recvcycle(icycle)
       source = pcycle%source-1
       length = nlyr*pcycle%lengthp_ghost*nc*nc
       tag = pcycle%tag
       iptr = pcycle%ptrp_ghost
!print *,'ghost_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call par_irecv(buffer%receive(1,1,1,iptr),length,source,ierr = ierr)
       if (ierr.ne.par_success) then
         errorcode = ierr
         call mpi_error_string(errorcode,errorstring,errorlen,ierr)
         print*,'bndry_exchangeV:error after call to mpi_irecv:',errorstring
       endif
     enddo ! icycle
!==================================================
!  Wait for all the receives to complete
!==================================================
     call par_allwait()
     do icycle = 1,nrecvcycles
       pcycle=>pschedule%recvcycle(icycle)
       length = pcycle%lengthp_ghost
       iptr = pcycle%ptrp_ghost
       do i = 0,length-1
         do i2 = 1, nc
           do i1 = 1, nc
             buffer%buf(i1,i2,1:nlyr,iptr+i) = buffer%receive(i1,i2,1:nlyr,iptr+i)
           enddo
           enddo
       enddo
     enddo ! icycle
#endif
   endif ! if (hybrid%ithr==0)
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
!
   end subroutine ghost_exchangevfull
! ===========================================
!  GHOST_EXCHANGEV:
!  Author: Christoph Erath
!  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ghost_exchangev(hybrid, buffer, nhc, npoints)
!
!   2011:  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
!
!    use hybrid_mod, only : hybrid_t
!-------------------------------------------------------------------------------
   use edge, only: ghostbuffertr_t
   use scheduler, only: schedule_t, cycle_t, schedule, status, srequest, rrequest, status, srequest, rrequest
   use dimensions, only: nelemd
   use kiapsparallel, only: parallel_t, hybrid_t, abortpar, omp_in_parallel, par_success, par_status_size, par_double_precision
   use kiapsparmpi, only: par_isend, par_irecv, par_allwait
   implicit none
!
   type(hybrid_t) :: hybrid
   type(ghostbuffertr_t) :: buffer
   integer :: nhc, npoints
   type(schedule_t), pointer :: pschedule
   type(cycle_t), pointer :: pcycle
   integer :: dest, length, tag
   integer :: icycle, ierr
   integer :: iptr, source, nlyr, ntrac
   integer :: nsendcycles, nrecvcycles
   integer :: errorcode, errorlen
   character*(80) errorstring
   integer :: i, i1, i2
   logical(log_kind), parameter :: debug = .false.
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
!
   if (hybrid%ithr==0) then
#ifdef _MPI
! Setup the pointer to proper Schedule
#ifdef _PREDICT
     pschedule=>schedule(iam)
#else
     pschedule=>schedule(1)
#endif
     nlyr = buffer%nlyr
     ntrac = buffer%ntrac
     nsendcycles = pschedule%nsendcycles
     nrecvcycles = pschedule%nrecvcycles
!==================================================
!  Fire off the sends
!==================================================
     do icycle = 1,nsendcycles
       pcycle=>pschedule%sendcycle(icycle)
       dest = pcycle%dest-1
       length = nlyr*ntrac*pcycle%lengthp_ghost*nhc*npoints
       tag = pcycle%tag
       iptr = pcycle%ptrp_ghost
!print *,'ghost_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call par_isend(buffer%buf(1,1,1,1,iptr),length,dest,ierr = ierr)
       if (ierr.ne.par_success) then
         errorcode = ierr
         call mpi_error_string(errorcode,errorstring,errorlen,ierr)
         print*,'bndry_exchangeV:error after call to mpi_isend:',errorstring
       endif
     enddo ! icycle
!==================================================
!  Post the Receives
!==================================================
     do icycle = 1,nrecvcycles
       pcycle=>pschedule%recvcycle(icycle)
       source = pcycle%source-1
       length = nlyr*ntrac*pcycle%lengthp_ghost*nhc*npoints
       tag = pcycle%tag
       iptr = pcycle%ptrp_ghost
!print *,'ghost_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call par_irecv(buffer%receive(1,1,1,1,iptr),length,source,ierr = ierr)
       if (ierr.ne.par_success) then
         errorcode = ierr
         call mpi_error_string(errorcode,errorstring,errorlen,ierr)
         print*,'bndry_exchangeV:error after call to mpi_irecv:',errorstring
       endif
     enddo ! icycle
!==================================================
!  Wait for all the receives to complete
!==================================================
     call par_allwait()
     do icycle = 1,nrecvcycles
       pcycle=>pschedule%recvcycle(icycle)
       length = pcycle%lengthp_ghost
       iptr = pcycle%ptrp_ghost
       do i = 0,length-1
         do i2 = 1, nhc
           do i1 = 1, npoints
             buffer%buf(i1,i2,1:nlyr,1:ntrac,iptr+i) = buffer%receive(i1,i2,1:nlyr,1:ntrac,iptr+i)
           enddo
           enddo
       enddo
     enddo ! icycle
#endif
   endif ! if (hybrid%ithr==0)
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
!
   end subroutine ghost_exchangev
! ===========================================
!  GHOST_EXCHANGEV:
!  Author: Christoph Erath
!  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ghost_exchangev_new(hybrid, buffer, nhc, npoints, ntrac)
!
!   2011:  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
!
!-------------------------------------------------------------------------------
   use kiapsparallel, only: hybrid_t, abortpar
   use kiapsbase, only: log_kind=>kim_log_kind
   use edge, only: ghostbuffertr_t
   use scheduler, only: schedule_t, cycle_t, schedule
   use dimensions, only: nelemd
#ifdef _MPI
!    use KiapsParallel, only : AbortPar, status, srequest, rrequest, &
!         mpireal_t, mpiinteger_t, mpi_success
   use kiapsparmpi, only: par_isend, par_irecv, par_allwait
#else
#endif
   implicit none
!
   type(hybrid_t) :: hybrid
   type(ghostbuffertr_t) :: buffer
   integer :: nhc, npoints, ntrac
   type(schedule_t), pointer :: pschedule
   type(cycle_t), pointer :: pcycle
   integer :: dest, length, tag
   integer :: icycle, ierr
   integer :: iptr, source, nlyr
   integer :: nsendcycles, nrecvcycles
   integer :: errorcode, errorlen
   character*(80) errorstring
   integer :: i, i1, i2
   logical(log_kind), parameter :: debug = .false.
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
!
   if (hybrid%ithr==0) then
#ifdef _MPI
! Setup the pointer to proper Schedule
#ifdef _PREDICT
     pschedule=>schedule(iam)
#else
     pschedule=>schedule(1)
#endif
     nlyr = buffer%nlyr
     nsendcycles = pschedule%nsendcycles
     nrecvcycles = pschedule%nrecvcycles
!==================================================
!  Fire off the sends
!==================================================
     do icycle = 1,nsendcycles
       pcycle=>pschedule%sendcycle(icycle)
       dest = pcycle%dest-1
       length = nlyr*ntrac*pcycle%lengthp_ghost*nhc*npoints
       tag = pcycle%tag
       iptr = pcycle%ptrp_ghost
!print *,'ghost_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
!call MPI_Isend(buffer%buf(1,1,1,1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
       call par_isend(buffer%buf(1,1,1,1,iptr),length,dest,ierr = ierr)
       if (ierr.ne.par_success) then
         errorcode = ierr
         call mpi_error_string(errorcode,errorstring,errorlen,ierr)
         print*,'bndry_exchangeV:error after call to mpi_isend:',errorstring
       endif
     enddo ! icycle
!==================================================
!  Post the Receives
!==================================================
     do icycle = 1,nrecvcycles
       pcycle=>pschedule%recvcycle(icycle)
       source = pcycle%source-1
       length = nlyr*ntrac*pcycle%lengthp_ghost*nhc*npoints
       tag = pcycle%tag
       iptr = pcycle%ptrp_ghost
!print *,'ghost_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
!call MPI_Irecv(buffer%receive(1,1,1,1,iptr),length,MPIreal_t,&
!     source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
       call par_irecv(buffer%receive(1,1,1,1,iptr),length,source,ierr = ierr)
       if (ierr.ne.par_success) then
         errorcode = ierr
         call mpi_error_string(errorcode,errorstring,errorlen,ierr)
         print*,'bndry_exchangeV:error after call to mpi_irecv:',errorstring
       endif
     enddo ! icycle
!==================================================
!  Wait for all the receives to complete
!==================================================
!call MPI_Waitall(nSendCycles,Srequest,status,ierr)
!call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
     call par_allwait()
     do icycle = 1,nrecvcycles
       pcycle=>pschedule%recvcycle(icycle)
       length = pcycle%lengthp_ghost
       iptr = pcycle%ptrp_ghost
       do i = 0,length-1
         do i2 = 1, nhc
           do i1 = 1, npoints
             buffer%buf(i1,i2,1:nlyr,1:ntrac,iptr+i) = buffer%receive(i1,i2,1:nlyr,1:ntrac,iptr+i)
           enddo
           enddo
       enddo
     enddo ! icycle
#endif
   endif ! if (hybrid%ithr==0)
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
!
   end subroutine ghost_exchangev_new
! ===========================================
!  GHOST_EXCHANGEV3d:
!  Author: James overflet
!  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine ghost_exchangev3d(hybrid, buffer)
!
!   2011:  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
!
!    use hybrid_mod,    only : hybrid_t
!-------------------------------------------------------------------------------
   use edge, only: ghostbuffer3d_t
   use scheduler, only: schedule_t, cycle_t, schedule, status, srequest, rrequest, status, srequest, rrequest
   use dimensions, only: nelemd
   use kiapsparallel, only: parallel_t, hybrid_t, abortpar, omp_in_parallel, par_double_precision, par_success, par_status_size, par_double_precision
   use kiapsparmpi, only: par_isend, par_irecv, par_allwait
   implicit none
!
   type(hybrid_t) :: hybrid
   type(ghostbuffer3d_t) :: buffer
   integer :: nhc, np
   type(schedule_t), pointer :: pschedule
   type(cycle_t), pointer :: pcycle
   integer :: dest, length, tag
   integer :: icycle, ierr
   integer :: iptr, source, nlyr
   integer :: nsendcycles, nrecvcycles
   integer :: errorcode, errorlen
   character*(80) errorstring
   integer :: i, i1, i2
   logical(log_kind), parameter :: debug = .false.
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
!
   if (hybrid%ithr==0) then
#ifdef _MPI
! Setup the pointer to proper Schedule
#ifdef _PREDICT
     pschedule=>schedule(iam)
#else
     pschedule=>schedule(1)
#endif
     nlyr = buffer%nlyr
     nhc = buffer%nhc
     np = buffer%np
     nsendcycles = pschedule%nsendcycles
     nrecvcycles = pschedule%nrecvcycles
!==================================================
!  Fire off the sends
!==================================================
     do icycle = 1,nsendcycles
       pcycle=>pschedule%sendcycle(icycle)
       dest = pcycle%dest-1
       length = nlyr*pcycle%lengthp_ghost*(nhc+1)*np
       tag = pcycle%tag
       iptr = pcycle%ptrp_ghost
!print *,'ghost_exchangeV3d: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call par_isend(buffer%buf(1,1,1,iptr),length,dest,ierr = ierr)
       if (ierr.ne.par_success) then
         errorcode = ierr
         call mpi_error_string(errorcode,errorstring,errorlen,ierr)
         print*,'ghost_exchangeV3d:error after call to mpi_isend:',errorstring
       endif
     enddo ! icycle
!==================================================
!  Post the Receives
!==================================================
     do icycle = 1,nrecvcycles
       pcycle=>pschedule%recvcycle(icycle)
       source = pcycle%source-1
       length = nlyr*pcycle%lengthp_ghost*(nhc+1)*np
       tag = pcycle%tag
       iptr = pcycle%ptrp_ghost
!print *,'ghost_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call par_irecv(buffer%receive(1,1,1,iptr),length,source,ierr = ierr)
       if (ierr.ne.par_success) then
         errorcode = ierr
         call mpi_error_string(errorcode,errorstring,errorlen,ierr)
         print*,'ghost_exchangeV3d:error after call to mpi_irecv:',errorstring
       endif
     enddo ! icycle
!==================================================
!  Wait for all the receives to complete
!==================================================
     call par_allwait()
     do icycle = 1,nrecvcycles
       pcycle=>pschedule%recvcycle(icycle)
       length = pcycle%lengthp_ghost
       iptr = pcycle%ptrp_ghost
       do i = 0,length-1
         buffer%buf(:,:,:,iptr+i) = buffer%receive(:,:,:,iptr+i)
!            do i2=1,nhc
!               do i1=1,np
!                  buffer%buf(i1,i2,1:nlyr,iptr+i) = buffer%receive(i1,i2,1:nlyr,iptr+i)
!               enddo
!            enddo
       enddo
     enddo ! icycle
#endif
   endif ! if (hybrid%ithr==0)
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
!
   end subroutine ghost_exchangev3d
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine compute_ghost_corner_orientation(hybrid, elem, nets, nete)
!-------------------------------------------------------------------------------
   use dimensions, only: nelemd, np
   use kiapsparallel, only: parallel_t, hybrid_t, abortpar, barrierpar, par_success, par_double_precision
   use element, only: element_t
   use edge, only: ghostbuffer_t, ghostvpackfull, ghostvunpackfull, initghostbufferfull, &
                                                                freeghostbuffer
   use control, only: north, south, east, west, neast, nwest, seast, swest
   implicit none
!
   type(hybrid_t), intent(in   ) :: hybrid
   type(element_t), intent(inout), target :: elem(:)
   integer :: nets, nete
   type(ghostbuffer_t) :: ghostbuf, ghostbuf_cv
   real(real_kind) :: cin(2, 2, 1, nets:nete) !ce:cslam tracer
   real(real_kind) :: cout(-1:4,-1:4, 1, nets:nete) !ce:cslam tracer
   integer :: i, j, ie, kptr, np1, np2, nc, nc1, nc2, k, nlev
   logical :: fail, fail1, fail2
   real(real_kind) :: tol = .1
!
   call barrierpar(hybrid%par)
!   if (hybrid%par%masterproc) print *,'computing ghost cell corner orientations'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first test on the Gauss Grid with same number of ghost cells:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nc = 2 ! test using gll interior points
   nc1 =-1
   nc2 = 4
   nlev = 1
   call initghostbufferfull(ghostbuf_cv,nlev,nc)
   do ie = nets,nete
     cin(1,1,1,ie) = elem(ie)%gdofp(1,1)
     cin(nc,nc,1,ie) = elem(ie)%gdofp(np,np)
     cin(1,nc,1,ie) = elem(ie)%gdofp(1,np)
     cin(nc,1,1,ie) = elem(ie)%gdofp(np,1)
   enddo
   cout = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  run ghost exchange on c array to get corner orientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do ie = nets,nete
     kptr = 0
     call ghostvpackfull(ghostbuf_cv,cin(:,:,:,ie),1,nc,nc,nlev,kptr,elem(ie)%desc)
   enddo
   call ghost_exchangevfull(hybrid,ghostbuf_cv,nc)
   do ie = nets,nete
     kptr = 0
     call ghostvunpackfull(ghostbuf_cv,cout(:,:,:,ie),nc1,nc2,nc,nlev,kptr,elem(ie)%desc)
   enddo
!       nc +--------+
!        ^ | nw  ne |
!     j  | |        |
!        1 | sw  se |
!          +--------+
!           1 --> nc
!              i
! check SW corner
   do ie = nets,nete
     fail1 = .false.
     fail2 = .false.
     if (elem(ie)%desc%putmapp_ghost(swest)/=-1) then
       if (abs(cout(nc1,1,1,ie)-cout(nc1,0,1,ie)).gt.tol) fail1 = .true.
       if (abs(cout(1,nc1,1,ie)-cout(0,nc1,1,ie)).gt.tol) fail2 = .true.
     endif
     if (fail1 .neqv. fail2) stop 'ghost exchange sw orientation failure'
     if (fail1) then
       elem(ie)%desc%reverse(swest) = .true.
!print *,'reversion sw orientation ie',ie
!print *,elem(ie)%desc%reverse(nwest),elem(ie)%desc%reverse(north),elem(ie)%desc%reverse(neast)
!print *,elem(ie)%desc%reverse(west),' ',elem(ie)%desc%reverse(east)
!print *,elem(ie)%desc%reverse(swest),elem(ie)%desc%reverse(south),elem(ie)%desc%reverse(seast)
     endif
   enddo
! check SE corner
   do ie = nets,nete
     fail1 = .false.
     fail2 = .false.
     if (elem(ie)%desc%putmapp_ghost(seast)/=-1) then
       if (abs(cout(nc2,1,1,ie)-cout(nc2,0,1,ie)).gt.tol) fail1 = .true.
       if (abs(cout(nc+1,nc1,1,ie)-cout(nc,nc1,1,ie)).gt.tol) fail2 = .true.
     endif
     if (fail1 .neqv. fail2) stop 'ghost exchange se orientation failure'
     if (fail1) then
       elem(ie)%desc%reverse(seast) = .true.
     endif
   enddo
! check NW corner
   do ie = nets,nete
     fail1 = .false.
     fail2 = .false.
     if (elem(ie)%desc%putmapp_ghost(nwest)/=-1) then
       if (abs(cout(nc1,nc+1,1,ie)-cout(nc1,nc,1,ie)).gt.tol) fail1 = .true.
       if (abs(cout(1,nc2,1,ie)-cout(0,nc2,1,ie)).gt.tol) fail2 = .true.
     endif
     if (fail1 .neqv. fail2) stop 'ghost exchange nw orientation failure'
     if (fail1) then
       elem(ie)%desc%reverse(nwest) = .true.
     endif
   enddo
! check NE corner
   do ie = nets,nete
     fail1 = .false.
     fail2 = .false.
     if (elem(ie)%desc%putmapp_ghost(neast)/=-1) then
       if (abs(cout(nc2,nc+1,1,ie)-cout(nc2,nc,1,ie)).gt.tol) fail1 = .true.
       if (abs(cout(nc+1,nc2,1,ie)-cout(nc,nc2,1,ie)).gt.tol) fail2 = .true.
     endif
     if (fail1 .neqv. fail2) stop 'ghost exchange ne orientation failure'
     if (fail1) then
       elem(ie)%desc%reverse(neast) = .true.
     endif
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  end ghost exchange corner orientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   end subroutine
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module bndry
!-------------------------------------------------------------------------------
