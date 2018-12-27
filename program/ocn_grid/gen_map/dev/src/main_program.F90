!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
   use kinds,      only: i4, r8
   use ocn_domain, only: comm, nprocs, rank, ipr, jpr, ijpr, mproc, nproc, mnproc,&
                         mm, nn, i0, j0, ii, jj,                               &
                         initialize, finalize, gen_map, print_map
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
   implicit none
!
   include 'mpif.h'
   integer(i4) :: i, err
   integer(i4), dimension(:), allocatable :: map
!
   call initialize()
!
   call mpi_barrier(comm,err)
   map = gen_map(mm,nn,i0,ii,j0,jj)
!
   call print_map(map,1)
   call mpi_barrier(comm,err)
!
   deallocate(map)
   call finalize()
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
   integer(i4),           intent(in   ) :: n
! local variables
   integer(i4) :: i
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
