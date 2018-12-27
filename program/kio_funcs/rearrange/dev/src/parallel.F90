!-------------------------------------------------------------------------------
   module parallel
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    2017-08-08  junghan kim    initial setup
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds, only: i4, l4, r4, r8
   implicit none
   private
   include 'mpif.h'
!
   type par_t
     integer(i4) :: comm, nprocs, rank, root
     logical(l4) :: ismaster
   end type par_t
!
   public :: par_t, par_init, par_finalize, decompose_1d
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_init(par, comm)
!-------------------------------------------------------------------------------
   implicit none
   type(par_t),           intent(inout) :: par
   integer(i4), optional, intent(in   ) :: comm
! local variables
   integer :: err
!
   if (present(comm)) then
     par%comm = comm
   else
     call mpi_init(err)
     par%comm = mpi_comm_world
   endif
   call mpi_comm_size(par%comm,par%nprocs,err)
   call mpi_comm_rank(par%comm,par%rank  ,err)
   par%root = 0
   if (par%rank.eq.par%root) then
     par%ismaster = .true.
   else
     par%ismaster = .false.
   endif
!
   return
   end subroutine par_init
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_finalize(par)
!-------------------------------------------------------------------------------
   implicit none
   type(par_t), intent(inout) :: par
! local variables
   integer :: err
!
   call mpi_finalize(err)
!
   return
   end subroutine par_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function decompose_1d(n1, n2, nps, ip, ista, iend) result(res)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: n1, n2
   integer(i4),           intent(in   ) :: nps, ip
   integer(i4), optional, intent(  out) :: ista, iend
   integer(i4)                          :: res
! local variables
   integer(i4) :: l_sta, l_end, domain, extra
!
   domain = n2-n1+1
   res    = domain/nps
   extra  = mod(domain,nps)
   l_sta  = res*ip+n1+min(ip,extra)
   l_end  = l_sta+res-1
   if (ip<=extra-1) then
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
   end function decompose_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module parallel
!-------------------------------------------------------------------------------
