!-------------------------------------------------------------------------------
   module mg_comm
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
   private
!
   integer, parameter :: i4 = 4, l4 = 4, r4 = 4, r8 = 8
   type parallel_t
     integer(i4) :: comm, nprocs, rank, root
     logical(l4) :: ismaster
   end type parallel_t
!
   public :: parallel_t, decompose_1d, gen_send_info, gen_recv_info,           &
             print_comm_info, convert_grid, convert_check, gather_var
!
   contains
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
! ss, ee: src start end |  r_sss, r_ees: dst start end
!-------------------------------------------------------------------------------
   subroutine gen_send_info(nprocs,ss,ee,r_sss,r_ees,nns,sss,ees)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                    intent(in   ) :: nprocs, ss, ee
   integer(i4), dimension(nprocs), intent(in   ) :: r_sss, r_ees
   integer(i4), dimension(nprocs), intent(inout) :: nns, sss, ees
! local variables
   integer(i4) :: i, j
!
   nns(:) = 0; sss(:) = -1; ees(:) = -1
   do i = ss,ee
     do j = 1,nprocs
       if (i>=r_sss(j).and.i<=r_ees(j)) then
         nns(j) = nns(j)+1
         if (nns(j).eq.1) then
           sss(j) = i-ss+1
         else
           ees(j) = i-ss+1
         endif
       endif
     enddo
   enddo
!
   return
   end subroutine gen_send_info
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine gen_recv_info(nprocs,ss,ee,s_sss,s_ees,nns,sss,ees)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                    intent(in   ) :: nprocs, ss, ee
   integer(i4), dimension(nprocs), intent(in   ) :: s_sss, s_ees
   integer(i4), dimension(nprocs), intent(inout) :: nns, sss, ees
! local variables
   integer(i4) :: i, j
   logical(l4) :: isused(nprocs)
!
   nns(:) = 0; sss(:) = -1; ees(:) = -1
   isused(:) = .false.
   do i = ss,ee
     do j = 1,nprocs
       if (i>=s_sss(j).and.i<=s_ees(j)) then
         nns(j) = nns(j)+1
         if (nns(j).eq.1) then
           isused(j) = .true.
         endif
       endif
     enddo
   enddo
   do j = 1,nprocs
     if (nns(j).gt.0) then
       if (j.eq.1) then
         sss(j) = 1
         ees(j) = nns(j)
       else
         if (isused(j-1)) then
           sss(j) = ees(j-1)+1
           ees(j) = sss(j)+nns(j)-1
         else
           sss(j) = 1
           ees(j) = nns(j)
         endif
       endif
     endif
   enddo
!
   return
   end subroutine gen_recv_info
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_comm_info(par,sizes_s,starts_s,ends_s,sizes_r,starts_r,ends_r)
!-------------------------------------------------------------------------------
   implicit none
!
   type(parallel_t),                   intent(in   ) :: par
   integer(i4), dimension(par%nprocs), intent(in   ) :: sizes_s, starts_s, ends_s
   integer(i4), dimension(par%nprocs), intent(in   ) :: sizes_r, starts_r, ends_r
! local variables
   integer(i4) :: comm, nprocs, rank
   integer(i4) :: i, j, ierr
!
   comm   = par%comm
   nprocs = par%nprocs
   rank   = par%rank
!
   call mpi_barrier(comm,ierr)
!
   do i = 1,nprocs
     if (i-1.eq.rank) then
       do j = 1,nprocs
         print '(i2,a,i2,a,i2,a,i2,a,i2)', rank, ': to ', j-1, ': size = ',    &
                      sizes_s(j), ', start:end = ', starts_s(j), ':', ends_s(j)
       enddo
       print *, ' '
     endif
     call mpi_barrier(comm,ierr)
   enddo
   call flush(6)

   call mpi_barrier(comm,ierr)
   do i = 1,nprocs
     if (i-1.eq.rank) then
       do j = 1,nprocs
         print '(i2,a,i2,a,i2,a,i2,a,i2)', rank, ': from ', j-1, ': size = ',  &
                      sizes_r(j), ', start:end = ', starts_r(j), ':', ends_r(j)
       enddo
       print *, ' '
     endif
     call mpi_barrier(comm,ierr)
   enddo
   call flush(6)

   call mpi_barrier(comm,ierr)
!
   return
   end subroutine print_comm_info
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine convert_grid(par,nx,lny,lmy,nn_s,ss_s,ee_s,nn_r,ss_r,ee_r,src,dst,use_alltoall)
!-------------------------------------------------------------------------------
   implicit none
!
   type(parallel_t),                   intent(in   ) :: par
   integer(i4),                        intent(in   ) :: nx, lny, lmy
   integer(i4), dimension(par%nprocs), intent(in   ) :: nn_s, ss_s, ee_s, nn_r, ss_r, ee_r
   real(r8),    dimension(nx,lny),     intent(in   ) :: src
   real(r8),    dimension(nx,lmy),     intent(in   ) :: dst
   logical(l4), optional,              intent(in   ) :: use_alltoall
! local variables
   integer(i4) :: comm, nprocs, rank, i, ierr
   integer(i4), dimension(par%nprocs) :: reqs_s, reqs_r
   logical(l4) :: use_ata
!
   comm   = par%comm
   nprocs = par%nprocs
   rank   = par%rank
!
   if (present(use_alltoall)) then
     use_ata = use_alltoall
   else
     use_ata = .false.
   endif
!
   if (use_ata) then
     !
     call mpi_alltoallv(src,nx*nn_s,nx*(ss_s-1),mpi_real8,dst,nx*nn_r,nx*(ss_r-1),mpi_real8,comm,ierr)
     !
   else
     !
     ! send cycle
     do i = 1,nprocs
       if (nn_s(i).gt.0) then
         call mpi_isend(src(:,ss_s(i):ee_s(i)),nx*nn_s(i),mpi_real8,i-1,100+i,comm,reqs_s(i),ierr)
       endif
     enddo
     ! recv cycle
     do i = 1,nprocs
       if (nn_r(i).gt.0) then
         call mpi_irecv(dst(:,ss_r(i):ee_r(i)),nx*nn_r(i),mpi_real8,i-1,mpi_any_tag,comm,reqs_r(i),ierr)
       endif
     enddo
     !
     ! wait
     do i = 1,nprocs
       if (nn_r(i).gt.0) then
         call mpi_wait(reqs_r(i),mpi_status_ignore,ierr)
       endif
     enddo
     do i = 1,nprocs
       if (nn_s(i).gt.0) then
         call mpi_wait(reqs_s(i),mpi_status_ignore,ierr)
       endif
     enddo
     !
   endif
!
   return
   end subroutine convert_grid
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine gather_var(par,nx,ny,lny,lvar,gvar)
!-------------------------------------------------------------------------------
   implicit none
!
   type(parallel_t),               intent(in   ) :: par
   integer(i4),                    intent(in   ) :: nx, ny, lny
   real(r8),    dimension(nx,lny), intent(in   ) :: lvar
   real(r8),    dimension(nx,ny) , intent(inout) :: gvar
! local variables
   integer(i4) :: comm, nprocs, rank
   integer(i4) :: i, ierr
   integer(i4) :: lnys(par%nprocs), cnts(par%nprocs), dspl(par%nprocs)
!
   comm   = par%comm
   nprocs = par%nprocs
   rank   = par%rank
!
   call mpi_allgather(lny,1,mpi_integer,lnys,1,mpi_integer,comm,ierr)
!
   dspl(1) = 0
   do i=1,nprocs
     cnts(i) = nx*lnys(i)
     if (i.ne.1) then
       dspl(i) = dspl(i-1)+cnts(i-1)
     endif
   enddo
!
   call mpi_gatherv(lvar,nx*lny,mpi_real8,gvar,cnts,dspl,mpi_real8,0,comm,ierr)
!
   return
   end subroutine gather_var
!-------------------------------------------------------------------------------
!    
!  
!-------------------------------------------------------------------------------
   subroutine convert_check(par,nx,ny,cy,ss,ee,var1,var2)
!-------------------------------------------------------------------------------
   implicit none
!  
   type(parallel_t),              intent(in   ) :: par
   integer(i4),                   intent(in   ) :: nx, ny, cy, ss, ee
   real(r8),    dimension(nx,ny), intent(in   ) :: var1
   real(r8),    dimension(nx,ny), intent(in   ) :: var2
! local variables
   integer(i4) :: i, j, jj, ierr
   integer(i4) :: comm, nprocs, rank
   logical(l4) :: issame, issame_g
!  
   comm   = par%comm
   nprocs = par%nprocs
   rank   = par%rank
!
   issame = .true.
!  
   do j = ss,ee
     if (ee.le.cy) then
       jj = j-ss+1
       do i = 1,nx
         if (var1(i,jj).ne.var2(i,jj)) then
           !if(rank.eq.0) print *, i, jj, var1(i,jj), var2(i,jj)
           issame = .false.
           exit
         endif
       enddo
     endif
   enddo
!  
   call mpi_reduce(issame,issame_g,1,mpi_logical,mpi_land,0,comm,ierr)
!
   if (rank.eq.0) then
     print *, ' '
     if (issame_g) then
       print *, 'two variables are same!'
     else  
       print *, 'two variables are different!'
     endif 
   endif
!      
   return
   end subroutine convert_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module mg_comm
!-------------------------------------------------------------------------------
