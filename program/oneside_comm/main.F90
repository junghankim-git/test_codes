!-------------------------------------------------------------------------------
   program main
!-------------------------------------------------------------------------------
   use kinds, only: i4, r4, r8, l4
   implicit none
!
   include 'mpif.h'
   ! mpi
   integer(i4) :: comm, nprocs, rank, root, err
   logical(l4) :: ismaster
   ! var
   integer(i4) :: w_task_id, w_task_win, cnt, var
   integer(mpi_address_kind) :: w_task_size, dsp
   logical(l4) :: myturn
!
! mpi
   call mpi_init(err)
   comm = mpi_comm_world
   call mpi_comm_size(comm,nprocs,err)
   call mpi_comm_rank(comm,  rank,err)
   root = 0
   ismaster = .false.
   if (rank.eq.root) ismaster = .true.
print *,' DEBUG step 0',rank
!
! code here
   ! initialize
!   if (ismaster) then
     w_task_id = -1
     call mpi_type_extent(mpi_integer,w_task_size,err)
     call check_mpi(err,comm,'mpi_type_extent')
     call mpi_win_create(w_task_id,w_task_size,w_task_size,mpi_info_null,comm,w_task_win,err)
     call check_mpi(err,comm,'mpi_win_create')
!   else
!     w_task_size = 0
!     call mpi_win_create(mpi_bottom,w_task_size,1,mpi_info_null,comm,w_task_win,err)
!     call check_mpi(err,comm,'mpi_win_create')
!   endif
!
   call mpi_barrier(comm,err)
   ! run code
print *,' DEBUG step 1',rank
!   if (.not.ismaster) then
     myturn = .false.
!     do while (.not.myturn)
       call mpi_win_fence(0,w_task_win,err)
       cnt = 1; dsp = 0
var = -99
       call mpi_get(var,cnt,mpi_integer,root,dsp,cnt,mpi_integer,w_task_win,err)
  print *, 'DEBUG var = ', var
       if (var.lt.0) then
         var = rank; cnt = 1; dsp = 0
         call mpi_put(var,cnt,mpi_integer,root,dsp,cnt,mpi_integer,w_task_win,err)
         myturn = .true.
       endif
       call mpi_win_fence(0,w_task_win,err)
       !call sleep(1)
!     enddo
print *,' DEBUG step 2',rank
     !
     print *,'process rank: ',rank
     !
!     call mpi_win_fence(0,w_task_win,err)
!     var = -1; cnt = 1; dsp = 0
!     call mpi_put(var,cnt,mpi_integer,root,dsp,cnt,mpi_integer,w_task_win,err)
!     call mpi_win_fence(0,w_task_win,err)
!   endif
!
print *,' DEBUG step 3',rank
   ! finalize
   call mpi_win_free(0,w_task_win,err)
   call check_mpi(err,comm,'mpi_win_free')
print *,' DEBUG step 4',rank
!
   call mpi_finalize(err)

!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_mpi(err,comm,message)
   implicit none
   integer(i4)     , intent(in   ) :: err, comm
   character(len=*), intent(in   ) :: message
! local variables
   integer(i4) :: ierr, errlen
   character(len=256) :: errstring
!
   if (err.ne.mpi_success) then
     call mpi_error_string(err,errstring,errlen,ierr)
     print*,trim(message)//': '//trim(errstring)
     call mpi_abort(comm,ierr)
   endif
!
   end subroutine check_mpi
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program main
!-------------------------------------------------------------------------------
