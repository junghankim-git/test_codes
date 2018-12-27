!------------------------------------------------------------------------------- 
module parallel
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    2017-03-09  junghan kim    initial setup
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds,    only: i4, l4, r8
!
   implicit none
   include 'mpif.h'
!
   private
!
   type par_t
     integer(i4) :: nprocs
     integer(i4) :: rank
     integer(i4) :: comm
     integer(i4) :: root
     logical(l4) :: ismaster
   end type par_t
   type(par_t) :: par
!
   public :: par, par_initialize, par_finalize, decompose_1d, onesided_comm
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_initialize()
!-------------------------------------------------------------------------------
   implicit none
!
! local variables
!
   integer(i4) :: ierr
!
   call mpi_init(ierr)
   par%comm = mpi_comm_world
   call mpi_comm_size(par%comm,par%nprocs,ierr)
   call mpi_comm_rank(par%comm,par%rank  ,ierr)
   par%root     = 0
   par%ismaster = .false.
   if (par%rank.eq.par%root) par%ismaster = .true.
!
   return
!
   end subroutine par_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_finalize()
!-------------------------------------------------------------------------------
   implicit none
!
! local variables
!
   integer(i4) :: ierr
!
   call mpi_finalize(par%comm,ierr)
!
   return
!
   end subroutine par_finalize
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
   function map2mpitype(nz, nmaps, map, max_map) result(dertype)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                   intent(in   ) :: nz, nmaps, max_map
   integer(i4), dimension(nmaps), intent(in   ) :: map
   integer(i4) :: dertype
!
! local variables
!
   integer(i4) :: i, j, ierr, type1, t_size
   integer(i4) :: ncnts
   integer(i4), dimension(:), allocatable :: bsize, disps, types
!
   if (minval(map)==1) then
     ncnts = nmaps
   else
     ncnts = nmaps+1
   endif
!
   if (maxval(map)==max_map) then
     ncnts = ncnts
   else
     ncnts = ncnts+1
   endif
!
   allocate(bsize(ncnts))
   allocate(disps(ncnts))
   allocate(types(ncnts))
!
   if (minval(map)==1) then
     do i = 1,ncnts-1
       bsize(i) = 1
       disps(i) = 8*(map(i)-1)
       types(i) = mpi_real8
     enddo
     if (maxval(map)==max_map) then
       bsize(ncnts) = 1
       disps(ncnts) = 8*(map(ncnts)-1)
       types(ncnts) = mpi_real8
     else
       bsize(ncnts) = 1
       disps(ncnts) = 8*(max_map) ! why?
       types(ncnts) = mpi_ub
     endif
   else
     bsize(1) = 1
     disps(1) = 0
     types(1) = mpi_lb
     do i = 2,ncnts-1
       bsize(i) = 1
       disps(i) = 8*(map(i-1)-1)
       types(i) = mpi_real8
     enddo
     if (maxval(map)==max_map) then
       bsize(ncnts) = 1
       disps(ncnts) = 8*(map(ncnts-1)-1)
       types(ncnts) = mpi_real8
     else
       bsize(ncnts) = 1
       disps(ncnts) = 8*(max_map) ! why?
       types(ncnts) = mpi_ub
     endif
   endif
!
   call mpi_type_struct(ncnts,bsize,disps,types,type1,ierr)
   call mpi_type_commit(type1,ierr)
   call mpi_type_contiguous(nz,type1,dertype,ierr)
   call mpi_type_commit(dertype,ierr)
!
   deallocate(bsize)
   deallocate(disps)
   deallocate(types)
!
   end function map2mpitype
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine onesided_comm(gnx, gny, nx, map, lvar, gvar)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),                     intent(in   ) :: gnx, gny, nx
   integer(i4), dimension( nx),     intent(in   ) :: map
   real(r8),    dimension( nx,gny), intent(in   ) :: lvar
   real(r8),    dimension(gnx,gny), intent(inout) :: gvar
!
! local variables
!
   integer(i4) :: i, j, ierr, gsize, lsize
   integer(i4) :: cnts, mpitype, var_mpitype, var_bytesize, gvar_win, win_bytesize
   integer(mpi_address_kind) :: disp, win_tsize, var1_lb, var1_bytesize
!
   gsize = gnx*gny
   lsize =  nx*gny
print *, 'DEBUG : step 1'
   var_mpitype = map2mpitype(gny,nx,map,gnx)
!
! local copy
   do j = 1,gny
     do i = 1,nx
       !gvar(gnx*(j-1)+map(i)) = lvar(i,j)
       gvar(map(i),j) = lvar(i,j)
     enddo
   enddo
!
! one-sided communication
print *, 'DEBUG : step 2'
   mpitype = mpi_real8
   call mpi_type_get_extent(mpitype,var1_lb,var1_bytesize,ierr)
   !if (par%ismaster) then
     win_tsize = var1_bytesize*gsize
     call mpi_win_create(gvar,win_tsize,var1_bytesize,mpi_info_null,par%comm,gvar_win,ierr)
   !else
   !  win_tsize = 0
   !  call mpi_win_create(mpi_bottom,win_tsize,1,mpi_info_null,par%comm,gvar_win,ierr)
   !endif
!
print *, 'DEBUG : step 3'
   call mpi_win_fence(0,gvar_win,ierr)
   !call mpi_win_fence(mpi_mode_noprecede,gvar_win,ierr)
!
print *, 'DEBUG : step 4'
   if (.not.par%ismaster) then
     disp = 0 ! careful
     cnts = 1
     !call mpi_win_lock(mpi_lock_exclusive,1,0,gvar_win,ierr)
     call mpi_put(lvar,lsize,mpitype,par%root,disp,cnts,var_mpitype,gvar_win,ierr)
     !call mpi_win_unlock(1,gvar_win,ierr)
   endif
!
print *, 'DEBUG : step 5', par%rank
   call mpi_win_fence(0,gvar_win,ierr)
   !call mpi_win_fence(mpi_mode_nostore,gvar_win,ierr)
!
print *, 'DEBUG : step 6', par%rank
   call mpi_win_free(gvar_win,ierr)
print *, 'DEBUG : step 7', par%rank
   call mpi_type_free(var_mpitype,ierr)
!
   return
!
   end subroutine onesided_comm
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module parallel
!-------------------------------------------------------------------------------
