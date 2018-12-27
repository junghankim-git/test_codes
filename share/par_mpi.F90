!#include <KIM.h>
!-------------------------------------------------------------------------------
   module par_mpi
!-------------------------------------------------------------------------------
!
!  abstract :
!    Parallelization module for framework
!  history log :
!    2012-12-06   junghan kim    First written
!    2013-03-19   junghan kim    Change module structure 
!                                integrated scalar & vector APIs
!
!-------------------------------------------------------------------------------
   use parallel
   use parallel, only: i4=>i4_par,r4=>r4_par,r8=>r8_par,l4=>l4_par
   use par_env

   integer(i4), parameter, public :: max_reqs = 50
   integer(i4),            public :: stat(par_status_size)
   integer(i4),            public :: ireqs = 0
   integer(i4),            public :: nreqs = 0
   integer(i4),            public :: req(max_reqs)
   integer(i4),            public :: stats(par_status_size, max_reqs)

   interface par_send
     module procedure par_send_int
     module procedure par_send_real4
     module procedure par_send_real8
     module procedure par_send_char
     module procedure par_send_log
   end interface 

   interface par_recv
     module procedure par_recv_int
     module procedure par_recv_real4
     module procedure par_recv_real8
     module procedure par_recv_char
     module procedure par_recv_log
   end interface 

   interface par_isend
     module procedure par_isend_int
     module procedure par_isend_real4
     module procedure par_isend_real8
     module procedure par_isend_char
     module procedure par_isend_log
   end interface 

   interface par_irecv
     module procedure par_irecv_int
     module procedure par_irecv_real4
     module procedure par_irecv_real8
     module procedure par_irecv_char
     module procedure par_irecv_log
   end interface 

   interface par_bcast
     module procedure par_bcast_int
     module procedure par_bcast_real4
     module procedure par_bcast_real8
     module procedure par_bcast_char
     module procedure par_bcast_log
   end interface 

   interface par_gather
     module procedure par_gather_int
     module procedure par_gather_real4
     module procedure par_gather_real8
     module procedure par_gather_char
     module procedure par_gather_log
   end interface 

   interface par_gatherv
     module procedure par_gatherv_int
     module procedure par_gatherv_real4
     module procedure par_gatherv_real8
     module procedure par_gatherv_char
     module procedure par_gatherv_log
   end interface 

   interface par_allgather
     module procedure par_allgather_int
     module procedure par_allgather_real4
     module procedure par_allgather_real8
     module procedure par_allgather_char
     module procedure par_allgather_log
   end interface 

   interface par_allgatherv
     module procedure par_allgatherv_int
     module procedure par_allgatherv_real4
     module procedure par_allgatherv_real8
     module procedure par_allgatherv_char
     module procedure par_allgatherv_log
   end interface 

   interface par_reduce
     module procedure par_reduce_int
     module procedure par_reduce_real4
     module procedure par_reduce_real8
     module procedure par_reduce_char
     module procedure par_reduce_log
   end interface 

   interface par_allreduce
     module procedure par_allreduce_int
     module procedure par_allreduce_real4
     module procedure par_allreduce_real8
     module procedure par_allreduce_char
     module procedure par_allreduce_log
   end interface 

   interface par_scatter
     module procedure par_scatter_int
     module procedure par_scatter_real4
     module procedure par_scatter_real8
     module procedure par_scatter_char
     module procedure par_scatter_log
   end interface 

   interface par_scatterv
     module procedure par_scatterv_int
     module procedure par_scatterv_real4
     module procedure par_scatterv_real8
     module procedure par_scatterv_char
     module procedure par_scatterv_log
   end interface 
!
   private
!
   public :: par_send, par_recv, par_isend, par_irecv, par_wait, par_allwait
   public :: par_bcast, par_gather, par_gatherv, par_allgather, par_allgatherv
   public :: par_reduce, par_allreduce, par_scatter, par_scatterv
   public :: par_gather_mchar
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_send_int(sendbuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Send (vector, integer type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_send(sendbuf, ncounts, par_integer, rank, l_tag, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_send_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_send_real4(sendbuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Send (vector, real4 type)
! 
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),              intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_send(sendbuf, ncounts, par_real, rank, l_tag, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_send_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_send_real8(sendbuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Send (vector, real8 type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),              intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_send(sendbuf, ncounts, par_real8, rank, l_tag, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_send_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_send_char(sendbuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Send (character type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=ncounts), intent(in   ) :: sendbuf
   integer(i4),            intent(in   ) :: ncounts, rank
   integer(i4), optional,  intent(in   ) :: comm
   integer(i4), optional,  intent(  out) :: ierr
   integer(i4)                           :: l_tag, l_comm
   integer(i4)                           :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_send(sendbuf, ncounts, par_character, rank, l_tag, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_send_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_send_log(sendbuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Send (vector, logical type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none

   logical(l4),           intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_send(sendbuf, ncounts, par_logical, rank, l_tag, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_send_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_recv_int(recvbuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Recv (vector, integer type)
! 
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(inout) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_recv(recvbuf, ncounts, par_integer, rank, l_tag, l_comm, stat, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_recv_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_recv_real4(recvbuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Recv (vector, real4 type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),              intent(inout) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_recv(recvbuf, ncounts, par_real, rank, l_tag, l_comm, stat, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_recv_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_recv_real8(recvbuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Recv (vector, real8 type)
! 
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),              intent(inout) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_recv(recvbuf, ncounts, par_real8, rank, l_tag, l_comm, stat, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_recv_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_recv_char(recvbuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Recv (character type)
! 
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=ncounts), intent(inout) :: recvbuf
   integer(i4),            intent(in   ) :: ncounts, rank
   integer(i4), optional,  intent(in   ) :: comm
   integer(i4), optional,  intent(  out) :: ierr
   integer(i4)                           :: l_tag, l_comm
   integer(i4)                           :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_recv(recvbuf,ncounts,par_character,rank,l_tag,l_comm,stat,l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_recv_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_recv_log(recvbuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Recv (vector, logical type)
! 
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),           intent(inout) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_recv(recvbuf,ncounts,par_logical,rank,l_tag,l_comm,stat,l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_recv_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_isend_int(sendbuf, ncounts, rank, request, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Isend (vector, integer type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
    implicit none
!
    integer(i4),           intent(in   ) :: sendbuf
    integer(i4),           intent(in   ) :: ncounts, rank
    integer(i4), optional, intent(  out) :: request
    integer(i4), optional, intent(in   ) :: comm
    integer(i4), optional, intent(  out) :: ierr
    integer(i4)                          :: l_tag, l_comm
    integer(i4)                          :: l_req
    integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_isend(sendbuf,ncounts,par_integer,rank,l_tag,l_comm,l_req,l_ierr)
!
   if (present(request)) then
     request = l_req
   else
     nreqs = nreqs + 1
     req(nreqs) = l_req
   endif
!
   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_isend_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_isend_real4(sendbuf, ncounts, rank, request, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Isend (vector, real4 type)
! 
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),              intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(  out) :: request
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_req
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_isend(sendbuf,ncounts,par_real,rank,l_tag,l_comm,l_req,l_ierr)
!
   if (present(request)) then
     request = l_req
   else
     nreqs = nreqs + 1
     req(nreqs) = l_req
   endif
!
   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_isend_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_isend_real8(sendbuf, ncounts, rank, request, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Isend (vector, real8 type)
! 
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),              intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(  out) :: request
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_req
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_isend(sendbuf,ncounts,par_real8,rank,l_tag,l_comm,l_req,l_ierr)
!
   if (present(request)) then
     request = l_req
   else
     nreqs = nreqs + 1
     req(nreqs) = l_req
   endif
!
   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_isend_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_isend_char(sendbuf, ncounts, rank, request, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Isend (character type)
! 
!  history log :
!    2013-01-17   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=ncounts), intent(in   ) :: sendbuf
   integer(i4),            intent(in   ) :: ncounts, rank
   integer(i4), optional,  intent(  out) :: request
   integer(i4), optional,  intent(in   ) :: comm
   integer(i4), optional,  intent(  out) :: ierr
   integer(i4)                           :: l_tag, l_comm
   integer(i4)                           :: l_req
   integer(i4)                           :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_isend(sendbuf,ncounts,par_character,rank,l_tag,l_comm,l_req,l_ierr)
!
   if (present(request)) then
     request = l_req
   else
     nreqs = nreqs + 1
     req(nreqs) = l_req
   endif
!
   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_isend_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_isend_log(sendbuf, ncounts, rank, request, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Isend (vector, logical type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),           intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(  out) :: request
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_req
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_isend(sendbuf,ncounts,par_logical,rank,l_tag,l_comm,l_req,l_ierr)
!
   if (present(request)) then
     request = l_req
   else
     nreqs = nreqs + 1
     req(nreqs) = l_req
   endif
!
   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_isend_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_irecv_int(recvbuf, ncounts, rank, request, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Irecv (vector, integer type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none

   integer(i4),           intent(inout) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(  out) :: request
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_req
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
    if (present(comm)) then
      l_comm = comm
    else
      l_comm = kim_par%comm
    endif
!
    l_tag = 111
    call mpi_irecv(recvbuf,ncounts,par_integer,rank,l_tag,l_comm,l_req,l_ierr)
!
    if (present(request)) then
      request = l_req
    else
      nreqs = nreqs + 1
      req(nreqs) = l_req
    endif
!
    if (present(ierr)) ierr = l_ierr
    return

    end subroutine par_irecv_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_irecv_real4(recvbuf, ncounts, rank, request, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_Irecv (vector, real4 type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),              intent(inout) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(  out) :: request
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_req
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_irecv(recvbuf,ncounts,par_real,rank,l_tag,l_comm,l_req,l_ierr)
!
   if (present(request)) then
     request = l_req
   else
     nreqs = nreqs + 1
     req(nreqs) = l_req
   endif
!
   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_irecv_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_irecv_real8(recvbuf, ncounts, rank, request, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Irecv (vector, real8 type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),              intent(inout) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(  out) :: request
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_req
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_irecv(recvbuf,ncounts,par_real8,rank,l_tag,l_comm,l_req,l_ierr)
!
   if (present(request)) then
     request = l_req
   else
     nreqs = nreqs + 1
     req(nreqs) = l_req
   endif
!
   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_irecv_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_irecv_char(recvbuf, ncounts, rank, request, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Irecv (character type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=ncounts), intent(inout) :: recvbuf
   integer(i4),            intent(in   ) :: ncounts, rank
   integer(i4), optional,  intent(  out) :: request
   integer(i4), optional,  intent(in   ) :: comm
   integer(i4), optional,  intent(  out) :: ierr
   integer(i4)                           :: l_tag, l_comm
   integer(i4)                           :: l_req
   integer(i4)                           :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_irecv(recvbuf,ncounts,par_character,rank,l_tag,l_comm,l_req,l_ierr)
!
   if (present(request)) then
     request = l_req
   else
     nreqs = nreqs + 1
     req(nreqs) = l_req
   endif
!
   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_irecv_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_irecv_log(recvbuf, ncounts, rank, request, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Irecv (vector, logical type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),           intent(inout) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts, rank
   integer(i4), optional, intent(  out) :: request
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_tag, l_comm
   integer(i4)                          :: l_req
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   l_tag = 111
   call mpi_irecv(recvbuf,ncounts,par_logical,rank,l_tag,l_comm,l_req,l_ierr)
!
   if (present(request)) then
     request = l_req
   else
     nreqs = nreqs + 1
     req(nreqs) = l_req
   endif
!
   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_irecv_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_wait(request, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Wait
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), optional, intent(in   ) :: request
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_req, l_i
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(request)) then
     l_req = request
   else
     if (nreqs .eq. 0) return
     l_req = req(1)
     do l_i = 2, nreqs
       req(l_i-1) = req(l_i)
     enddo
     nreqs = nreqs - 1
   endif
!
   call mpi_wait(l_req, stat, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_wait
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allwait(ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Wait
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_req, l_i
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
   if (nreqs .eq. 0) return
!   do l_i = 1, nreqs
!     call mpi_wait(req(l_i), stats(:,l_i), l_ierr)
!   enddo
   call mpi_waitall(nreqs, req(:), stats(:,:), l_ierr)
   nreqs = 0

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allwait
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_bcast_int(combuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Bcast (vector, integer type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(inout) :: combuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_bcast(combuf, ncounts, par_integer, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_bcast_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_bcast_real4(combuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Bcast (vector, real4 type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),              intent(inout) :: combuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_bcast(combuf, ncounts, par_real, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_bcast_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_bcast_real8(combuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Bcast (vector, real8 type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
    implicit none
!
    real(r8),              intent(inout) :: combuf
    integer(i4),           intent(in   ) :: ncounts
    integer(i4), optional, intent(in   ) :: rank, comm
    integer(i4), optional, intent(  out) :: ierr
    integer(i4)                          :: l_rank, l_comm
    integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_bcast(combuf,ncounts,par_double_precision,l_rank,l_comm,l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_bcast_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_bcast_char(combuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Bcast (character type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),            intent(in   ) :: ncounts
   !character(len=ncounts), intent(inout) :: combuf
   character(len=1),       intent(inout) :: combuf
   integer(i4), optional,  intent(in   ) :: rank, comm
   integer(i4), optional,  intent(  out) :: ierr
   integer(i4)                           :: l_rank, l_comm
   integer(i4)                           :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_bcast(combuf,ncounts,par_character,l_rank,l_comm,l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_bcast_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_bcast_log(combuf, ncounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Bcast (vector, logical type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),           intent(inout) :: combuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_bcast(combuf,ncounts,par_logical,l_rank,l_comm,l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_bcast_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_gather_int(sendbuf, ncounts, recvbuf, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Gather (vector, integer type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4),           intent(  out) :: recvbuf
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_gather(sendbuf, ncounts, par_integer, recvbuf, ncounts,            &
                   par_integer, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_gather_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_gather_real4(sendbuf, ncounts, recvbuf, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Gather (vector, real4 type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),              intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts
   real(r4),              intent(  out) :: recvbuf
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_gather(sendbuf, ncounts, par_real, recvbuf, ncounts,               &
                   par_real, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_gather_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_gather_real8(sendbuf, ncounts, recvbuf, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Gather (vector, real8 type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),              intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts
   real(r8),              intent(  out) :: recvbuf
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_gather(sendbuf, ncounts, par_double_precision, recvbuf, ncounts,   &
                   par_double_precision, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_gather_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_gather_char(sendbuf, ncounts, recvbuf, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Gather (character type)
!
!  history log :
!    2013-01-18   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=ncounts), intent(in   ) :: sendbuf
   integer(i4),            intent(in   ) :: ncounts
   character(len=*),       intent(  out) :: recvbuf
   integer(i4), optional,  intent(in   ) :: rank, comm
   integer(i4), optional,  intent(  out) :: ierr
   integer(i4)                           :: l_rank, l_comm
   integer(i4)                           :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_gather(sendbuf, ncounts, par_character, recvbuf, ncounts,          &
                   par_character, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_gather_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_gather_mchar(sendbuf, ncounts, recvbuf, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Gather (multi dimension character type)
!
!  history log :
!    2013-01-11   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=ncounts),               intent(in   ) :: sendbuf
   integer(i4),                          intent(in   ) :: ncounts
   character(len=ncounts), dimension(:), intent(  out) :: recvbuf
   integer(i4), optional,                intent(in   ) :: rank, comm
   integer(i4), optional,                intent(  out) :: ierr
   integer(i4)                                         :: l_rank, l_comm
   integer(i4)                                         :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_gather(sendbuf, ncounts, par_character, recvbuf, ncounts,          &
                   par_character, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_gather_mchar
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_gather_log(sendbuf, ncounts, recvbuf, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Gather (vector, logical type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),           intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts
   logical(l4),           intent(  out) :: recvbuf
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_gather(sendbuf, ncounts, par_logical, recvbuf, ncounts,  &
                   par_logical, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_gather_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_gatherv_int(sendbuf, sendcounts, recvbuf, recvcounts, displs,&
                              rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_GatherV (vector, integer type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: sendbuf
   integer(i4),               intent(in   ) :: sendcounts
   integer(i4),               intent(  out) :: recvbuf
   integer(i4), dimension(:), intent(in   ) :: recvcounts
   integer(i4), dimension(:), intent(in   ) :: displs 
   integer(i4), optional,     intent(in   ) :: rank, comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_rank, l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_gatherv(sendbuf, sendcounts, par_integer, recvbuf, recvcounts,     &
                    displs, par_integer, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_gatherv_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_gatherv_real4(sendbuf, sendcounts, recvbuf, recvcounts,      &
                                displs, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_GatherV (vector, real4 type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),                  intent(in   ) :: sendbuf
   integer(i4),               intent(in   ) :: sendcounts
   real(r4),                  intent(  out) :: recvbuf
   integer(i4), dimension(:), intent(in   ) :: recvcounts
   integer(i4), dimension(:), intent(in   ) :: displs 
   integer(i4), optional,     intent(in   ) :: rank, comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_rank, l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_gatherv(sendbuf, sendcounts, par_real, recvbuf, recvcounts,        &
                    displs, par_real, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_gatherv_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_gatherv_real8(sendbuf, sendcounts, recvbuf, recvcounts,      &
                                displs, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_GatherV (vector, real8 type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),                  intent(in   ) :: sendbuf
   integer(i4),               intent(in   ) :: sendcounts
   real(r8),                  intent(  out) :: recvbuf
   integer(i4), dimension(:), intent(in   ) :: recvcounts
   integer(i4), dimension(:), intent(in   ) :: displs 
   integer(i4), optional,     intent(in   ) :: rank, comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_rank, l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_gatherv(sendbuf, sendcounts, par_double_precision, recvbuf,        &
                    recvcounts, displs, par_double_precision, l_rank, l_comm,  &
                    l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_gatherv_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_gatherv_char(sendbuf, sendcounts, recvbuf, recvcounts,       &
                               displs, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_GatherV (character type)
!
!  history log :
!    2013-01-17   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: sendcounts
   character(len=sendcounts), intent(in   ) :: sendbuf
   character(len=*),          intent(  out) :: recvbuf
   integer(i4), dimension(:), intent(in   ) :: recvcounts
   integer(i4), dimension(:), intent(in   ) :: displs 
   integer(i4), optional,     intent(in   ) :: rank, comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_rank, l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_gatherv(sendbuf, sendcounts, par_character, recvbuf, recvcounts,   &
                    displs, par_character, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_gatherv_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_gatherv_log(sendbuf, sendcounts, recvbuf, recvcounts,        &
                              displs, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_GatherV (vector, logical type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),               intent(in   ) :: sendbuf
   integer(i4),               intent(in   ) :: sendcounts
   logical(l4),               intent(  out) :: recvbuf
   integer(i4), dimension(:), intent(in   ) :: recvcounts
   integer(i4), dimension(:), intent(in   ) :: displs 
   integer(i4), optional,     intent(in   ) :: rank, comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_rank, l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_gatherv(sendbuf, sendcounts, par_logical, recvbuf, recvcounts,     &
                    displs, par_logical, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_gatherv_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allgather_int(sendbuf, sendcounts, recvbuf, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllGather (vector, integer type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: sendcounts
   integer(i4),           intent(  out) :: recvbuf
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allgather(sendbuf, sendcounts, par_integer, recvbuf, sendcounts,   &
                      par_integer, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allgather_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allgather_real4(sendbuf, sendcounts, recvbuf, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllGather (vector, real4 type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),              intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: sendcounts
   real(r4),              intent(  out) :: recvbuf
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allgather(sendbuf, sendcounts, par_real, recvbuf, sendcounts,      &
                      par_real, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allgather_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allgather_real8(sendbuf, sendcounts, recvbuf, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllGather (vector, real8 type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),              intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: sendcounts
   real(r8),              intent(  out) :: recvbuf
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allgather(sendbuf, sendcounts, par_double_precision, recvbuf,      &
                      sendcounts, par_double_precision, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allgather_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allgather_char(sendbuf, sendcounts, recvbuf, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllGather (character type)
!
!  history log :
!    2013-01-17   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: sendcounts
   character(len=sendcounts), intent(in   ) :: sendbuf
   character(len=*),          intent(  out) :: recvbuf
   integer(i4), optional,     intent(in   ) :: comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allgather(sendbuf, sendcounts, par_character, recvbuf, sendcounts, &
                      par_character, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allgather_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allgather_log(sendbuf, sendcounts, recvbuf, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllGather (vector, logical type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),           intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: sendcounts
   logical(l4),           intent(  out) :: recvbuf
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allgather(sendbuf, sendcounts, par_logical, recvbuf, sendcounts,   &
                      par_logical, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allgather_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allgatherv_int(sendbuf, sendcounts, recvbuf, recvcounts,     &
                                 displs, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllGatherV (vector, integer type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: sendbuf
   integer(i4),               intent(in   ) :: sendcounts
   integer(i4),               intent(  out) :: recvbuf
   integer(i4), dimension(:), intent(in   ) :: recvcounts
   integer(i4), dimension(:), intent(in   ) :: displs 
   integer(i4), optional,     intent(in   ) :: comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allgatherv(sendbuf, sendcounts, par_integer, recvbuf, recvcounts,  &
                       displs, par_integer, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allgatherv_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allgatherv_real4(sendbuf, sendcounts, recvbuf, recvcounts,   &
                                   displs, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllGatherV (vector, real4 type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),                  intent(in   ) :: sendbuf
   integer(i4),               intent(in   ) :: sendcounts
   real(r4),                  intent(  out) :: recvbuf
   integer(i4), dimension(:), intent(in   ) :: recvcounts
   integer(i4), dimension(:), intent(in   ) :: displs 
   integer(i4), optional,     intent(in   ) :: comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allgatherv(sendbuf, sendcounts, par_real, recvbuf, recvcounts, &
               displs, par_real, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allgatherv_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allgatherv_real8(sendbuf, sendcounts, recvbuf, recvcounts,   &
                                   displs, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllGatherV (vector, real8 type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),                  intent(in   ) :: sendbuf
   integer(i4),               intent(in   ) :: sendcounts
   real(r8),                  intent(  out) :: recvbuf
   integer(i4), dimension(:), intent(in   ) :: recvcounts
   integer(i4), dimension(:), intent(in   ) :: displs 
   integer(i4), optional,     intent(in   ) :: comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allgatherv(sendbuf, sendcounts, par_double_precision, recvbuf,     &
                       recvcounts, displs, par_double_precision, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allgatherv_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allgatherv_char(sendbuf, sendcounts, recvbuf, recvcounts,    &
                                  displs, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllGatherV (character type)
!
!  history log :
!    2013-01-17   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: sendcounts
   character(len=sendcounts), intent(in   ) :: sendbuf
   character(len=*),          intent(  out) :: recvbuf
   integer(i4), dimension(:), intent(in   ) :: recvcounts
   integer(i4), dimension(:), intent(in   ) :: displs 
   integer(i4), optional,     intent(in   ) :: comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allgatherv(sendbuf, sendcounts, par_character, recvbuf, recvcounts,&
                       displs, par_character, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allgatherv_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------

   subroutine par_allgatherv_log(sendbuf, sendcounts, recvbuf, recvcounts,     &
                                 displs, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllGatherV (vector, logical type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),               intent(in   ) :: sendbuf
   integer(i4),               intent(in   ) :: sendcounts
   logical(l4),               intent(  out) :: recvbuf
   integer(i4), dimension(:), intent(in   ) :: recvcounts
   integer(i4), dimension(:), intent(in   ) :: displs 
   integer(i4), optional,     intent(in   ) :: comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allgatherv(sendbuf, sendcounts, par_logical, recvbuf, recvcounts,  &
                       displs, par_logical, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allgatherv_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_reduce_int(sendbuf, recvbuf, ncounts, par_op, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Reduce (vector, integer type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: sendbuf
   integer(i4),           intent(  out) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4), optional, intent(in   ) :: par_op
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_par_op, l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(par_op)) then
     l_par_op = par_op
   else
     l_par_op = par_sum
   endif
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_reduce(sendbuf, recvbuf, ncounts, par_integer,                     &
                   l_par_op, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_reduce_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_reduce_real4(sendbuf, recvbuf, ncounts, par_op, rank, comm,  &
                               ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Reduce (vector, real4 type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),              intent(in   ) :: sendbuf
   real(r4),              intent(  out) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4), optional, intent(in   ) :: par_op
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_par_op, l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(par_op)) then
     l_par_op = par_op
   else
     l_par_op = par_sum
   endif
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_reduce(sendbuf, recvbuf, ncounts, par_real4,                       &
                   l_par_op, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_reduce_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_reduce_real8(sendbuf, recvbuf, ncounts, par_op, rank, comm,  &
                               ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Reduce (vector, real8 type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),              intent(in   ) :: sendbuf
   real(r8),              intent(  out) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4), optional, intent(in   ) :: par_op
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_par_op, l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(par_op)) then
     l_par_op = par_op
   else
     l_par_op = par_sum
   endif
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_reduce(sendbuf, recvbuf, ncounts, par_double_precision,            &
                   l_par_op, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_reduce_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_reduce_char(sendbuf, recvbuf, ncounts, par_op, rank, comm,   &
                              ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Reduce (character type)
!
!  history log :
!    2013-01-17   junghan kim    First written
!
!------------------------------------------------------------------------------
   implicit none
!
   character(len=ncounts), intent(in   ) :: sendbuf
   character(len=ncounts), intent(  out) :: recvbuf
   integer(i4),            intent(in   ) :: ncounts
   integer(i4), optional,  intent(in   ) :: par_op
   integer(i4), optional,  intent(in   ) :: rank, comm
   integer(i4), optional,  intent(  out) :: ierr
   integer(i4)                           :: l_par_op, l_rank, l_comm
   integer(i4)                           :: l_ierr
!------------------------------------------------------------------------------
!
   if (present(par_op)) then
     l_par_op = par_op
   else
     l_par_op = par_sum
   endif
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_reduce(sendbuf, recvbuf, ncounts, par_character,                   &
                   l_par_op, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_reduce_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_reduce_log(sendbuf, recvbuf, ncounts, par_op, rank, comm,    &
                             ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Reduce (vector, logical type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),           intent(in   ) :: sendbuf
   logical(l4),           intent(  out) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4), optional, intent(in   ) :: par_op
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_par_op, l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(par_op)) then
     l_par_op = par_op
   else
     l_par_op = par_sum
   endif
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_reduce(sendbuf, recvbuf, ncounts, par_logical,                     &
                   l_par_op, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_reduce_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allreduce_int(sendbuf, recvbuf, ncounts, par_op, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllReduce (vector, integer type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: sendbuf
   integer(i4),           intent(  out) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4), optional, intent(in   ) :: par_op
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_par_op, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(par_op)) then
     l_par_op = par_op
   else
     l_par_op = par_sum
   endif
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allreduce(sendbuf, recvbuf, ncounts, par_integer,                  &
                      l_par_op, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allreduce_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allreduce_real4(sendbuf, recvbuf, ncounts, par_op, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllReduce (vector, real4 type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),              intent(in   ) :: sendbuf
   real(r4),              intent(  out) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4), optional, intent(in   ) :: par_op
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_par_op, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(par_op)) then
     l_par_op = par_op
   else
     l_par_op = par_sum
   endif
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allreduce(sendbuf, recvbuf, ncounts, par_real4,                    &
                      l_par_op, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allreduce_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allreduce_real8(sendbuf, recvbuf, ncounts, par_op, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllReduce (vector, real8 type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),              intent(in   ) :: sendbuf
   real(r8),              intent(  out) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4), optional, intent(in   ) :: par_op
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_par_op, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(par_op)) then
     l_par_op = par_op
   else
     l_par_op = par_sum
   endif
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allreduce(sendbuf, recvbuf, ncounts, par_double_precision,         &
                      l_par_op, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allreduce_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allreduce_char(sendbuf, recvbuf, ncounts, par_op, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllReduce (vector, character type)
!
!  history log :
!    2013-01-17   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=ncounts), intent(in   ) :: sendbuf
   character(len=ncounts), intent(  out) :: recvbuf
   integer(i4),            intent(in   ) :: ncounts
   integer(i4), optional,  intent(in   ) :: par_op
   integer(i4), optional,  intent(in   ) :: comm
   integer(i4), optional,  intent(  out) :: ierr
   integer(i4)                           :: l_par_op, l_comm
   integer(i4)                           :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(par_op)) then
     l_par_op = par_op
   else
     l_par_op = par_sum
   endif
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allreduce(sendbuf, recvbuf, ncounts, par_character,                &
                      l_par_op, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allreduce_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_allreduce_log(sendbuf, recvbuf, ncounts, par_op, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_AllReduce (vector, logical type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),           intent(in   ) :: sendbuf
   logical(l4),           intent(  out) :: recvbuf
   integer(i4),           intent(in   ) :: ncounts
   integer(i4), optional, intent(in   ) :: par_op
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_par_op, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(par_op)) then
     l_par_op = par_op
   else
     l_par_op = par_sum
   endif
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   call mpi_allreduce(sendbuf, recvbuf, ncounts, par_logical,                  &
                      l_par_op, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_allreduce_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_scatter_int(sendbuf, ncounts, recvbuf, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Scatter (vector, integer type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   )                :: sendbuf
   integer(i4), intent(in   )                :: ncounts
   integer(i4), intent(  out)               :: recvbuf
   integer(i4), optional, intent(in   )         :: rank, comm
   integer(i4), optional, intent(  out)        :: ierr

   ! local variables
   integer(i4)                         :: l_rank, l_comm
   integer(i4)                         :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_scatter(sendbuf, ncounts, par_integer, recvbuf, ncounts,           &
                    par_integer, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_scatter_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_scatter_real4(sendbuf, ncounts, recvbuf, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Scatter (vector, real4 type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), intent(in   )                 :: sendbuf
   integer(i4), intent(in   )                :: ncounts
   real(r4),                intent(  out) :: recvbuf
   integer(i4), optional, intent(in   )         :: rank, comm
   integer(i4), optional, intent(  out)        :: ierr

   ! local variable
   integer(i4)                         :: l_rank, l_comm
   integer(i4)                         :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_scatter(sendbuf, ncounts, par_real, recvbuf, ncounts,              &
                    par_real, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_scatter_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_scatter_real8(sendbuf, ncounts, recvbuf, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Scatter (vector, real8 type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),              intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts
   real(r8),              intent(  out) :: recvbuf
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_scatter(sendbuf, ncounts, par_double_precision, recvbuf, ncounts,  &
                    par_double_precision, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_scatter_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_scatter_char(sendbuf, ncounts, recvbuf, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Scatter (character type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),       intent(in   ) :: sendbuf
   integer(i4),            intent(in   ) :: ncounts
   character(len=ncounts), intent(  out) :: recvbuf
   integer(i4), optional,  intent(in   ) :: rank, comm
   integer(i4), optional,  intent(  out) :: ierr
   integer(i4)                           :: l_rank, l_comm
   integer(i4)                           :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_scatter(sendbuf, ncounts, par_character, recvbuf, ncounts,         &
                    par_character, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_scatter_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_scatter_log(sendbuf, ncounts, recvbuf, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_Scatter (vector, logical type)
!
!  history log :
!    2012-12-06   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),           intent(in   ) :: sendbuf
   integer(i4),           intent(in   ) :: ncounts
   logical(l4),           intent(  out) :: recvbuf
   integer(i4), optional, intent(in   ) :: rank, comm
   integer(i4), optional, intent(  out) :: ierr
   integer(i4)                          :: l_rank, l_comm
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_scatter(sendbuf, ncounts, par_logical, recvbuf, ncounts,           &
                    par_logical, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_scatter_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_scatterv_int(sendbuf, sendcounts, displs, recvbuf,           &
                               recvcounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract 
!    Wrapped routine for MPI_ScatterV (vector, integer type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: sendbuf
   integer(i4), dimension(:), intent(in   ) :: sendcounts, displs 
   integer(i4),               intent(  out) :: recvbuf
   integer(i4),               intent(in   ) :: recvcounts
   integer(i4), optional,     intent(in   ) :: rank, comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_rank, l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_scatterv(sendbuf, sendcounts, displs, par_integer, recvbuf,        &
                     recvcounts, par_integer, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_scatterv_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_scatterv_real4(sendbuf, sendcounts, displs, recvbuf,         &
                                 recvcounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_ScatterV (vector, real4 type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),                  intent(in   ) :: sendbuf
   integer(i4), dimension(:), intent(in   ) :: sendcounts, displs 
   real(r4),                  intent(  out) :: recvbuf
   integer(i4),               intent(in   ) :: recvcounts
   integer(i4), optional,     intent(in   ) :: rank, comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_rank, l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_scatterv(sendbuf, sendcounts, displs, par_real, recvbuf,           &
                     recvcounts, par_real, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_scatterv_real4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_scatterv_real8(sendbuf, sendcounts, displs, recvbuf,         &
                                 recvcounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_ScatterV (vector, real8 type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),                  intent(in   ) :: sendbuf
   integer(i4), dimension(:), intent(in   ) :: sendcounts, displs 
   real(r8),                  intent(  out) :: recvbuf
   integer(i4),               intent(in   ) :: recvcounts
   integer(i4), optional,     intent(in   ) :: rank, comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_rank, l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_scatterv(sendbuf, sendcounts, displs, par_double_precision,        &
                     recvbuf, recvcounts,                                      &
                     par_double_precision, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_scatterv_real8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_scatterv_char(sendbuf, sendcounts, displs, recvbuf,          &
                                recvcounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_ScatterV (character type)
!
!  history log :
!    2013-01-17   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),          intent(in   ) :: sendbuf
   integer(i4), dimension(:), intent(in   ) :: sendcounts, displs 
   integer(i4),               intent(in   ) :: recvcounts
   character(len=recvcounts), intent(  out) :: recvbuf
   integer(i4), optional,     intent(in   ) :: rank, comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_rank, l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_scatterv(sendbuf, sendcounts, displs,                              &
                 par_character, recvbuf, recvcounts,                           &
                 par_character, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_scatterv_char
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_scatterv_log(sendbuf, sendcounts, displs, recvbuf,           &
                               recvcounts, rank, comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Wrapped routine for MPI_ScatterV (vector, logical type)
!
!  history log :
!    2013-01-14   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4),               intent(in   ) :: sendbuf
   integer(i4), dimension(:), intent(in   ) :: sendcounts, displs 
   logical(l4),               intent(  out) :: recvbuf
   integer(i4),               intent(in   ) :: recvcounts
   integer(i4), optional,     intent(in   ) :: rank, comm
   integer(i4), optional,     intent(  out) :: ierr
   integer(i4)                              :: l_rank, l_comm
   integer(i4)                              :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(comm)) then
     l_comm = comm
   else
     l_comm = kim_par%comm
   endif
!
   if (present(rank)) then
     l_rank = rank
   else
     l_rank = kim_par%root
   endif
!
   call mpi_scatterv(sendbuf, sendcounts, displs, par_logical, recvbuf,        &
                     recvcounts, par_logical, l_rank, l_comm, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_scatterv_log
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module par_mpi
!-------------------------------------------------------------------------------
