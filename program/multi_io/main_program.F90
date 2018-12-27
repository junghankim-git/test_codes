!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    2017-05-18  junghan kim    initial setup
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds,           only: i4, l4, r4, r8
   use multi_io,        only: io_initialize, io_finalize, get_io_wtime
   use task_control_io, only: task_t, task_initialize, task_driver, task_profile
   implicit none
   include 'mpif.h'
!
   integer(i4)  :: comm, nprocs, rank, root, ierr
   logical(l4)  :: ismaster
!
   type(task_t) :: job
   integer(i4)  :: nfiles, nvars, ncols, nlevs
!
! mpi setting
   call mpi_init(ierr)
   comm = mpi_comm_world
   root = 0
   call mpi_comm_size(comm,nprocs,ierr)
   call mpi_comm_rank(comm,  rank,ierr)
   ismaster = .false.
   if (rank.eq.root) ismaster = .true.
! multi io
   !call io_initialize(50,6,200000,50,'/scratch/jhkim/data/multi_io') ! 458 MB x 50
!   nfiles = 350; nvars =   6; ncols = 200000; nlevs =    50 ! 50 x 458 MB = 23 GB
!   nfiles =  50; nvars =   6; ncols = 200000; nlevs =    50 ! 50 x 458 MB = 23 GB
!   nfiles =   1; nvars = 300; ncols = 200000; nlevs =    50
!   nfiles =   1; nvars =   1; ncols = 200000; nlevs = 15000
!   nfiles =   1; nvars =   6; ncols = 200000; nlevs =  2500
!
   nfiles = 100; nvars =   6; ncols = 200000; nlevs =    50 ! 50 x 458 MB = 23 GB
   call io_initialize(nfiles,nvars,ncols,nlevs,'/scratch/jhkim/data/multi_io')
!
!- START: job control
!
   !call task_initialize(job,nfiles,comm,nprocs,rank,.true.)
   call task_initialize(job,nfiles,comm,nprocs,rank)
!
   call task_driver(job)
!
!-  END : job control
!
!   print '(a,f6.2)', 'Total I/O wtime = ', get_io_wtime()
   call task_profile(job)
   call io_finalize()
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
