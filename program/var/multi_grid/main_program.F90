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
   use mg_comm, only: parallel_t, decompose_1d, gen_send_info, gen_recv_info,  &
                      print_comm_info, convert_grid, convert_check, gather_var
   implicit none
   include 'mpif.h'
!
   integer,     parameter   :: i4 = 4, l4 = 4, r4 = 4, r8 = 8
   integer(i4), parameter   :: nz = 3, gn = 20, gm = 10
   integer(i4)              :: ln, lm, sn, en, sm, em    ! domain
   type(parallel_t)         :: par
   integer(i4)              :: i, j, ierr, comm, nprocs, rank
   integer(i4), dimension(:),   allocatable :: lns, sns, ens ! ln, sn, en for all procs.
   integer(i4), dimension(:),   allocatable :: lms, sms, ems ! ln, sn, en for all procs.
                                       ! n->m, procs, size,  start, end
   integer(i4), dimension(:),   allocatable :: nnm_s, snm_s, enm_s
   integer(i4), dimension(:),   allocatable :: nnm_r, snm_r, enm_r
                                       ! m->n, procs, size,  start, end
   integer(i4), dimension(:),   allocatable :: nmn_s, smn_s, emn_s
   integer(i4), dimension(:),   allocatable :: nmn_r, smn_r, emn_r
   ! variable
   real(r8),    dimension(:,:), allocatable :: var_grid1, var_ret
   real(r8),    dimension(:,:), allocatable :: var_grid2
   real(r8),    dimension(:,:), allocatable :: gvar
   ! profile
   real(r8) :: stime, etime

!
   call mpi_init(ierr)
   comm = mpi_comm_world
   call mpi_comm_size(comm,nprocs,ierr)
   call mpi_comm_rank(comm,  rank,ierr)
   par%comm     = comm
   par%nprocs   = nprocs
   par%rank     = rank
   par%root     = 0
   par%ismaster = .false.
   if (par%rank.eq.par%root) par%ismaster = .true.
!
! 0) domain information
   ln = decompose_1d(1,gn,par%nprocs,par%rank,sn,en)
   lm = decompose_1d(1,gm,par%nprocs,par%rank,sm,em)
   allocate(var_grid1(nz,ln),var_ret(nz,ln))
   allocate(var_grid2(nz,lm))
   allocate(gvar(nz,gn))
   call var_initialize(nz,ln,sn,en,var_grid1)
   var_grid2(:,:) = 0.0
   var_ret(:,:)   = 0.0
   gvar(:,:)      = 0.0
!
   allocate(lns(nprocs),sns(nprocs),ens(nprocs))
   allocate(lms(nprocs),sms(nprocs),ems(nprocs))
   do i = 0,nprocs-1
     lns(i+1) = decompose_1d(1,gn,nprocs,i,sns(i+1),ens(i+1))
     lms(i+1) = decompose_1d(1,gm,nprocs,i,sms(i+1),ems(i+1))
   enddo
!
   allocate(nnm_s(nprocs),snm_s(nprocs),enm_s(nprocs))
   allocate(nnm_r(nprocs),snm_r(nprocs),enm_r(nprocs))
   allocate(nmn_s(nprocs),smn_s(nprocs),emn_s(nprocs))
   allocate(nmn_r(nprocs),smn_r(nprocs),emn_r(nprocs))
!
! 1) communication information
   call gen_send_info(nprocs,sn,en,sms,ems,nnm_s,snm_s,enm_s)
   call gen_recv_info(nprocs,sm,em,sns,ens,nnm_r,snm_r,enm_r)
   call print_comm_info(par,nnm_s,snm_s,enm_s,nnm_r,snm_r,enm_r)
   call gen_send_info(nprocs,sm,em,sns,ens,nmn_s,smn_s,emn_s)
   call gen_recv_info(nprocs,sn,en,sms,ems,nmn_r,smn_r,emn_r)
   call print_comm_info(par,nmn_s,smn_s,emn_s,nmn_r,smn_r,emn_r)
!
! 2) convert grid
   stime = mpi_wtime(ierr)
   call convert_grid(par,nz,ln,lm,nnm_s,snm_s,enm_s,nnm_r,snm_r,enm_r,var_grid1,var_grid2)
   call convert_grid(par,nz,lm,ln,nmn_s,smn_s,emn_s,nmn_r,smn_r,emn_r,var_grid2,var_ret)
   !call convert_grid(par,nz,ln,lm,nnm_s,snm_s,enm_s,nnm_r,snm_r,enm_r,var_grid1,var_grid2,.true.)
   !call convert_grid(par,nz,lm,ln,nmn_s,smn_s,emn_s,nmn_r,smn_r,emn_r,var_grid2,var_ret,.true.)
   etime = mpi_wtime(ierr)
   print *, 'wall-clock time = ', etime-stime
!
! 3) checking
   if (rank.eq.0) then
     print *, ' # local  variable (grid1)'
   endif
   do i = 1,nprocs
     if (i-1.eq.rank) then
       print *, '  - rank = ', rank
       do j = 1,ln
         print '(a,3(f5.1,x))', '  ', var_grid1(:,j)
       enddo
     endif
     call mpi_barrier(comm,ierr)
   enddo
   if (rank.eq.0) then
     print *, ' '
     print *, ' # local  variable (grid2)'
   endif
   do i = 1,nprocs
     if (i-1.eq.rank) then
       print *, '  - rank = ', rank
       do j = 1,lm
         print '(a,3(f5.1,x))', '  ', var_grid2(:,j)
       enddo
     endif
     call mpi_barrier(comm,ierr)
   enddo
!
   call mpi_barrier(comm,ierr)
   call gather_var(par,nz,gn,ln,var_ret,gvar)
   if (rank==0) then
     print *, ' '
     print *, ' # global variable'
     do j = 1,gn
       print '(a,3(f5.1,x))', '  ', gvar(:,j)
     enddo
   endif
   call flush(6)
   call mpi_barrier(comm,ierr)
   call convert_check(par,nz,ln,gm,sn,en,var_grid1,var_ret)
!
   deallocate(nnm_s,snm_s,enm_s)
   deallocate(nnm_r,snm_r,enm_r)
   deallocate(nmn_s,smn_s,emn_s)
   deallocate(nmn_r,smn_r,emn_r)
   deallocate(lns,sns,ens)
   deallocate(lms,sms,ems)
   deallocate(var_grid1,var_ret)
   deallocate(var_grid2)
   deallocate(gvar)
   call mpi_finalize(ierr)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine var_initialize(nx,ny,sy,ey,var)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                intent(in   ) :: nx, ny, sy, ey
   real(r8), dimension(nx,ny), intent(inout) :: var
! local variables
   integer(i4) :: i, j
!
   do j = sy,ey
     do i = 1,nx
       var(i,j-sy+1) = real(j,r8)+real(i,r8)*100.0_r8
     enddo
   enddo
!
   return
   end subroutine var_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program test
!-------------------------------------------------------------------------------
