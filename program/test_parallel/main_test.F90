!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
   use parallel,  only : kim_par, par_initialize, par_finalize
   use par_func,  only : globalmin, globalmax, globalsum, decompose1d
   use repro_sum, only : ddpdd
!-------------------------------------------------------------------------------
   implicit none
!
   include 'mpif.h'
!
   integer, parameter :: i4 = 4, r8 = 8
   integer, parameter :: nsize = 100
   integer(i4) :: comm, nprocs, rank, ierr, stat(mpi_status_size)
!
   integer(i4) :: i, j, ista, iend
   integer(i4) :: lsize
!
   integer(i4), dimension(:), allocatable :: sizes, istas, iends
   real(r8), dimension(nsize) :: garray
   real(r8), dimension(:), allocatable :: larray
   real(r8) :: lsum, gsum, gmin, gmax, mmax, mmin
!
   call par_initialize()
   comm = kim_par%comm
   nprocs = kim_par%nprocs
   rank = kim_par%iproc
!
   allocate(sizes(0 : nprocs-1))
   allocate(istas(0 : nprocs-1))
   allocate(iends(0 : nprocs-1))
!
   lsize = decompose1d(1, nsize, nprocs, rank, ista, iend)
!
   allocate(larray(lsize))
   do i = 0, nprocs-1
     sizes(i) = decompose1d(1, nsize, nprocs, i, istas(i), iends(i))
   enddo
!
   if (rank==0) then
     do i = 1, nsize
       call random_number(garray(i))
     enddo
   endif
!
   if (rank==0) then
     larray(:) = garray(ista : iend)
     do i = 1, nprocs-1
       call mpi_send(garray(istas(i) : iends(i)), sizes(i), mpi_real8, i, 111, comm, ierr)
     enddo
   else
     call mpi_recv(larray(:), lsize, mpi_real8, 0, 111, comm, stat, ierr)
   endif
!
   lsum = sum(larray)
   call mpi_reduce(lsum, gsum, 1, mpi_real8, mpi_sum, 0, comm, ierr)
!
   if (rank==0) then
     mmin = minval(garray)
     mmax = maxval(garray)
   endif
   gmin = globalmin(larray)
   gmax = globalmax(larray)
!
   if (rank==0) then
     print*, 'min/max in master = ', mmin, mmax
     print*, 'min/max in global = ', gmin, gmax
     print*, 'sum in global = ', gsum
   endif
!
   gsum = globalsum(larray)
   if (rank==0) then
     print*, 'sum in global = ', gsum
   endif
!
   deallocate(larray)
   deallocate(istas)
   deallocate(iends)
   deallocate(sizes)
!
   call par_finalize()
!
   end program test
