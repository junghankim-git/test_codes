!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  junghan kim    initial setup
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds,       only: i4, l4, r4, r8
   use task_control, only: task_t, task_initialize, task_driver, task_profile
   implicit none
   include 'mpif.h'
!
   integer(i4)  :: comm, nprocs, rank, root, ierr
   logical(l4)  :: ismaster
!
   integer(i4)  :: wtask
   logical(l4)  :: issync, isserial
!
   type(task_t) :: task
   integer(i4)  :: ntask = 100
!
! mpi setting
   call mpi_init(ierr)
   comm = mpi_comm_world
   root = 0
   call mpi_comm_size(comm,nprocs,ierr)
   call mpi_comm_rank(comm,  rank,ierr)
   ismaster = .false.
   if (rank.eq.root) ismaster = .true.
!
   wtask    = -1
   issync   = .false.
   isserial = .false.
   call task_initialize(task,ntask,comm,nprocs,rank,isserial=isserial,wtask=wtask)
   if (issync) then
     call task_driver_sync(task)
   else
     call task_driver(task)
   endif
   call task_profile(task)
!
   call mpi_finalize(ierr)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine test_sub(n)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: n
! local variables
   integer :: i
!
!
   return
   end subroutine test_sub
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program test
!-------------------------------------------------------------------------------
