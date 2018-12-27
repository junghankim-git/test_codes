!-------------------------------------------------------------------------------
   module ocn_domain
!-------------------------------------------------------------------------------
   use kinds,    only: i4, r8
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
   include 'mpif.h'
!
   integer(i4), parameter :: mm = 16, nn = 12
   integer(i4) :: comm, nprocs, rank
   integer(i4) :: ipr, jpr, ijpr, mproc, nproc, mnproc
   integer(i4) :: i0, j0, ii, jj
   integer(i4), dimension(:,:), allocatable :: i0_pe, ii_pe, j0_pe, jj_pe
!
   public :: comm, nprocs, rank, ipr, jpr, ijpr, mproc, nproc, mnproc
   public :: i0, j0, ii, jj
   public :: initialize, finalize, gen_map, print_map
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine initialize()
!-------------------------------------------------------------------------------
   implicit none
!
!   integer(i4), intent(in   ) :: nprocs, rank
! local variables
   integer(i4) :: err, i, j, iend, jend
   integer(i4), dimension(:,:), allocatable :: tmp
!
! mpi
   call mpi_init(err)
   comm = mpi_comm_world
   call mpi_comm_size(comm,nprocs,err)
   call mpi_comm_rank(comm,  rank,err)
!
   if (mod(nprocs,4).ne.0) then
     if (rank.eq.0) then
       print*,'error in initialize'
       call flush()
     endif
     call mpi_finalize(err)
     stop
   endif
!
! proc. domain
   ijpr   = nprocs
   ipr    = 4
   jpr    = nprocs/4
!
   mnproc = rank+1
   mproc  = mod(mnproc,ipr)
   if (mproc.eq.0) mproc = 4
   nproc  = (mnproc-1)/ipr+1
!
   allocate(tmp(nprocs,6))
!
   call mpi_gather(mproc,1,mpi_integer,tmp(:,1),1,mpi_integer,0,comm,err)
   call mpi_gather(nproc,1,mpi_integer,tmp(:,2),1,mpi_integer,0,comm,err)
   if (rank.eq.0) then
     do i = 1,nprocs
       print '(i2,a,i3,x,i3)',i-1,': ',tmp(i,1),tmp(i,2)
     enddo
     print *, ' '
   endif
   call flush()
   call mpi_barrier(comm,err)
!
! grid domain
   ii = decompose_1d(1,mm,ipr,mproc-1,i0,iend)
   i0 = i0-1
   jj = decompose_1d(1,nn,jpr,nproc-1,j0,jend)
   j0 = j0-1
   call mpi_gather(  i0,1,mpi_integer,tmp(:,1),1,mpi_integer,0,comm,err)
   call mpi_gather(iend,1,mpi_integer,tmp(:,2),1,mpi_integer,0,comm,err)
   call mpi_gather(  ii,1,mpi_integer,tmp(:,3),1,mpi_integer,0,comm,err)
   call mpi_gather(  j0,1,mpi_integer,tmp(:,4),1,mpi_integer,0,comm,err)
   call mpi_gather(jend,1,mpi_integer,tmp(:,5),1,mpi_integer,0,comm,err)
   call mpi_gather(  jj,1,mpi_integer,tmp(:,6),1,mpi_integer,0,comm,err)
!
   if (rank.eq.0) then
     do i = 1,nprocs
       print '(i2,a,3(i3,x),a,3(i3,x),a)',i-1,': (',tmp(i,1),tmp(i,2),tmp(i,3),'), (',tmp(i,4),tmp(i,5),tmp(i,6),')'
     enddo
     print *, ' '
   endif
!
   deallocate(tmp)
!
   return
   end subroutine initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine finalize()
!-------------------------------------------------------------------------------
   implicit none
!
! local variables
   integer(i4) :: err
!
   call mpi_finalize(err)
!
   return
   end subroutine finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function gen_map(mmm,nnn,i00,iii,j00,jjj) result(map)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: mmm, nnn, i00, iii, j00, jjj
   integer(i4), dimension(:), allocatable :: map
! local variables
   integer(i4) :: err, i, j, v
!
   allocate(map(iii*jjj))
   v = 1
   do j =j00+1,j00+jjj
     do i = i00+1,i00+iii
       map(v) = i+mmm*(j-1)
       v = v+1
     enddo
   enddo
!
   return
   end function gen_map
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_map(map, vrank)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), dimension(:), intent(in   ) :: map
   integer(i4), optional,     intent(in   ) :: vrank
! local variables
   integer(i4) :: err, i, nmap, vproc
   character(len=512) :: string
!
   if (present(vrank)) then
     vproc = vrank
   else
     vproc = 0
   endif
!
   nmap = size(map)
   string = ''
   do i = 1,nmap
     if (i.eq.1) then
       write(string,'(a,i3)') ' ',map(i)
     else
       if (map(i)-map(i-1).eq.1) then
         write(string,'(a,i3)') trim(string)//', ',map(i)
       else
         write(string,'(a,i3)') trim(string)//new_line('')//'  ',map(i)
       endif
     endif
   enddo
   if (rank.eq.vproc) print*,trim(string)//new_line('')
!
   return
   end subroutine print_map
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function decompose_1d(n1, n2, nprocs, rank, ista, iend) result(res)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: n1, n2
   integer(i4),           intent(in   ) :: nprocs, rank
   integer(i4), optional, intent(  out) :: ista, iend
   integer(i4)                          :: res
!
! local variables
!
   integer(i4) :: l_sta, l_end, domain, extra
!
   domain = n2-n1+1
   res    = domain/nprocs
   extra  = mod(domain,nprocs)
   l_sta  = res*rank+n1+min(rank,extra)
   l_end  = l_sta+res-1
   if (rank<=extra-1) then
     l_end = l_end+1
   endif
   res = l_end-l_sta+1
!
   if (present(ista)) then
     ista = l_sta
   endif
!
   if (present(iend)) then
     iend = l_end
   endif
!
   return
!
   end function decompose_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module ocn_domain
!-------------------------------------------------------------------------------
