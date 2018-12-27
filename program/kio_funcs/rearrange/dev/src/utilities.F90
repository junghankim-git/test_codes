!-------------------------------------------------------------------------------
   module utilities
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
   use kinds   , only: i4, l4, r4, r8
   use parallel, only: par_t, par_init, par_finalize, decompose_1d
   implicit none
   include 'mpif.h'
   private
!
   type comm_t ! all-to-all communication
     integer(i4)              :: nsend, nrecv
     integer(i4), allocatable :: sproc(:), ssize(:), sista(:) ! send (nrprocs)
     integer(i4), allocatable :: rproc(:), rsize(:), rtype(:) ! recv (nprocs)
   end type comm_t
!
   type meta_t
     ! mpi
     type(par_t) :: gp, rp
     integer(i4) :: nrprocs
     logical(l4) :: isrproc
     ! grid
     integer(i4) :: gsize, lsize
     integer(i4), dimension(:), allocatable :: map
     ! reduced grid
     integer(i4) :: rsize
     integer(i4), dimension(:), allocatable :: rmap
     ! com
     type(comm_t) :: cm
   end type meta_t
!
   public :: meta_t, initialize, finalize, processing
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine initialize(meta)
!-------------------------------------------------------------------------------
   implicit none
   type(meta_t), intent(inout) :: meta
! local variables
   integer(i4) :: comm, rank, nprocs, nrprocs, rcomm, err
   integer(i4), allocatable :: ranks_r(:), sizes_r(:), istas_r(:)
!
   ! parallel
   call par_init(meta%gp)
   if (meta%gp%rank.eq.0) write(*,'(a)') 'log) initialize (start)'
   comm   = meta%gp%comm
   nprocs = meta%gp%nprocs
   rank   = meta%gp%rank
!
   ! reduced process
   meta%nrprocs = nprocs/2
   if (meta%nrprocs.le.0) meta%nrprocs = 1
   nrprocs = meta%nrprocs
   allocate(ranks_r(nrprocs),sizes_r(nrprocs),istas_r(nrprocs))
   call make_rproc(comm,nprocs,rank,nrprocs,ranks_r,meta%isrproc,rcomm)
   if (meta%isrproc) then
     call par_init(meta%rp,rcomm)
   endif
!
   ! grid
   meta%gsize = 20
!
   ! local map
   call make_map(meta%gsize,nprocs,rank,meta%lsize,meta%map,.true.)
!
   ! reduced grid & map
   call make_rmap(meta%gsize,nrprocs,meta%rp%rank,ranks_r,                     &
                                          istas_r,sizes_r,meta%rsize,meta%rmap)
!
   call print_status(comm,nprocs,rank,meta%map,meta%isrproc,                   &
                            meta%rp%comm,meta%rp%nprocs,meta%rp%rank,meta%rmap)
   call mpi_barrier(comm,err)
!
   ! make communication address
   call make_comm_address(meta%map,meta%isrproc,comm,nprocs,nrprocs,           &
                                           ranks_r,      sizes_r,      istas_r,&
                                                   meta%cm%nsend,meta%cm%nrecv,&
                                     meta%cm%sproc,meta%cm%ssize,meta%cm%sista,&
                                     meta%cm%rproc,meta%cm%rsize,meta%cm%rtype)
   call print_address(comm,nprocs,rank,rcomm,meta%isrproc,                     &
                                                   meta%cm%nsend,meta%cm%nrecv,&
                                     meta%cm%sproc,meta%cm%sista,meta%cm%ssize,&
                                     meta%cm%rproc,meta%cm%rtype,meta%cm%rsize)
!
   deallocate(ranks_r,sizes_r,istas_r)
   if (meta%gp%rank.eq.0) write(*,'(a)') 'log) initialize (done)'//new_line('')
   call mpi_barrier(meta%gp%comm,err)
!
   return
   end subroutine initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine finalize(meta)
!-------------------------------------------------------------------------------
   implicit none
   type(meta_t), intent(inout) :: meta
! local variables
   integer :: err
!
   if (meta%gp%rank.eq.0) write(*,'(a)') 'log) finalize (start)'
   deallocate(meta%map)
   if (meta%isrproc) deallocate(meta%rmap)
   deallocate(meta%cm%sproc,meta%cm%ssize,meta%cm%sista)
   deallocate(meta%cm%rproc,meta%cm%rsize,meta%cm%rtype)
   call par_finalize(meta%gp)
   if (meta%gp%rank.eq.0) write(*,'(a)') 'log) finalize (done)'//new_line('')
!
   return
   end subroutine finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine processing(meta)
!-------------------------------------------------------------------------------
   implicit none
   type(meta_t), intent(inout) :: meta
! local variables
   integer :: err
!
   if (meta%gp%rank.eq.0) write(*,'(a)') 'log) processing (start)'
!
!
   if (meta%gp%rank.eq.0) write(*,'(a)') 'log) processing (done)'//new_line('')
   call mpi_barrier(meta%gp%comm,err)
!
   return
   end subroutine processing
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_rproc(comm,nprocs,rank,nrprocs,ranks_r,isrproc,gcomm)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4), intent(in   ) :: comm, nprocs, rank, nrprocs
   integer(i4), intent(  out) :: ranks_r(nrprocs)
   logical(l4), intent(  out) :: isrproc
   integer(i4), intent(  out) :: gcomm
! local variables
   integer(i4) :: err, i, tmp, cgrp, grp
!
   isrproc = .false.
   do i = 1,nrprocs
     tmp = decompose_1d(0,nprocs-1,nrprocs,i-1,ranks_r(i))
     if (rank.eq.ranks_r(i)) then
       isrproc = .true.
     endif
   enddo
!
   call mpi_comm_group(comm,cgrp,err)
   call mpi_group_incl(cgrp,nrprocs,ranks_r,grp,err)
   call mpi_comm_create(comm,grp,gcomm,err)
   call mpi_group_free(cgrp,err)
   call mpi_group_free(grp,err)
!
   return
   end subroutine make_rproc
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_map(gsize,nprocs,rank,lsize,map,block)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),              intent(in   ) :: gsize, nprocs, rank
   integer(i4),              intent(inout) :: lsize
   integer(i4), allocatable, intent(inout) :: map(:)
   logical(l4), optional   , intent(in   ) :: block
! local variables
   integer(i4) :: ip, ii, i, ista, iend
   logical(l4) :: isblock
!
   if (present(block)) then
     isblock = block
   else
     isblock = .false.
   endif
!
   lsize = decompose_1d(1,gsize,nprocs,rank,ista,iend)
   allocate(map(lsize))
!
   if (isblock) then
     do i = ista,iend
       map(i-ista+1) = i
     enddo
   else
     ip = 0
     ii = 0
     do i = 1,gsize
       if (ip.eq.rank) then
         ii = ii + 1
         map(ii) = i
       endif
       ip = ip+1
       if (ip.eq.nprocs) ip = 0
     enddo
   endif
!
   return
   end subroutine make_map
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_rmap(gsize,nrprocs,rrank,ranks_r,istas_r,sizes_r,rsize,rmap)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),              intent(in   ) :: gsize
   integer(i4),              intent(in   ) :: nrprocs, rrank
   integer(i4),              intent(in   ) :: ranks_r(nrprocs)
   integer(i4),              intent(inout) :: istas_r(nrprocs), sizes_r(nrprocs)
   integer(i4),              intent(inout) :: rsize
   integer(i4), allocatable, intent(inout) :: rmap(:)
! local variables
   integer(i4) :: ip,i
!
   do ip = 0,nrprocs-1
     sizes_r(ip+1) = decompose_1d(1,gsize,nrprocs,ip,istas_r(ip+1))
     if (ip.eq.rrank) then
       rsize = sizes_r(ip+1)
       allocate(rmap(rsize))
       do i = 1,rsize
         rmap(i) = istas_r(ip+1)+i-1
       enddo
     endif
   enddo
!
   return
   end subroutine make_rmap
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_status(comm,nprocs,rank,map,isrproc,rcomm,nrprocs,rrank,rmap)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4), intent(in   ) :: comm, nprocs, rank
   logical(l4), intent(in   ) :: isrproc
   integer(i4), intent(in   ) :: rcomm, nrprocs, rrank
   integer(i4), intent(in   ) :: map(:), rmap(:)
! local variables
!
   ! gproc
   call print_proc(comm,nprocs,rank,isrproc,nrprocs)
!
   ! print map
   call print_map(comm,nprocs,rank,map,'map')
   if (isrproc) then
     call print_map(rcomm,nrprocs,rrank,rmap,'rmap')
   endif
!
   return
   end subroutine print_status
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_proc(comm,nprocs,rank,isrproc,nrprocs)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4), intent(in   ) :: comm, nprocs, rank, nrprocs
   logical(l4), intent(in   ) :: isrproc
! local variables
   integer(i4) :: i, err
   logical(l4) :: isg(nprocs)
!
   ! gproc
   call mpi_gather(isrproc,1,mpi_logical,isg,1,mpi_logical,0,comm,err)
   if (rank.eq.0) then
     write(*,'(a)')           'log)  info.(process)'
     write(*,'(a,i3,a,i3,a)') 'log)   nprocs = ',nprocs,' (nrprocs =',nrprocs,')'
     do i = 1,nprocs
       write(*,'(a,i3,a,l3)') 'log)   rank   = ',i-1,': ',isg(i)
     enddo
   endif
!
   return
   end subroutine print_proc
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_map(comm,nprocs,rank,map,message)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4)     , intent(in   ) :: comm, nprocs, rank
   integer(i4)     , intent(in   ) :: map(:)
   character(len=*), intent(in   ) :: message
! local variables
   integer(i4) :: i, err, max_map
   integer(i4), allocatable :: nmaps(:), disps(:), maps(:)
   character         :: num
   character(len=32) :: frmt
!
   max_map = maxval(map)
   call mpi_reduce(max_map,max_map,1,mpi_integer,mpi_max,0,comm,err)
   if (rank.eq.0) then
     if (max_map.ge.1.and.max_map.le.9) then
       num = '1'
     elseif (max_map.ge.10.and.max_map.le.99) then
       num = '2'
     elseif (max_map.ge.100) then
       num = '3'
     endif
   endif
   ! print map
   call gather_map(comm,nprocs,rank,map,nmaps,disps,maps)
   if (rank.eq.0) then
     write(*,'(a)')           'log)  info.('//trim(message)//')'
     do i = 1,nprocs
       write(frmt,'(a,i3,a)') '(a,i'//num//',a,',nmaps(i),'(i'//num//',x))'
       write(*,trim(frmt))  'log)   rank   = ',i-1,': ',maps(disps(i)+1:disps(i)+nmaps(i))
     enddo
   endif
   deallocate(nmaps,disps,maps)
!
   return
   end subroutine print_map
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine gather_map(comm,nprocs,rank,map,nmaps,disps,maps)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),              intent(in   ) :: comm, nprocs, rank
   integer(i4),              intent(in   ) :: map(:)
   integer(i4), allocatable, intent(inout) :: nmaps(:),disps(:),maps(:)
! local variables
   integer(i4) :: i, err, nmaps_proc
!
   ! gather map
   allocate(nmaps(nprocs),disps(nprocs))
   call mpi_gather(size(map),1,mpi_integer,nmaps,1,mpi_integer,0,comm,err)
   if (rank.eq.0) then
     disps(1) = 0
     do i = 2,nprocs
       disps(i) = disps(i-1)+nmaps(i-1)
     enddo
     nmaps_proc = sum(nmaps)
   else
     nmaps_proc = 1
   endif
   allocate(maps(nmaps_proc))
   call mpi_gatherv(map,size(map),mpi_integer,maps,nmaps,disps,mpi_integer,0,comm,err)
!
   return
   end subroutine gather_map
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_comm_address(map,isrproc,comm,nprocs,nrprocs,ranks_r,sizes_r,istas_r,&
                               nsend,nrecv,sproc,ssize,sista,rproc,rsize,rtype)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4), dimension(:)             , intent(in   ) :: map
   logical(l4),                            intent(in   ) :: isrproc
   integer(i4),                            intent(in   ) :: comm, nprocs, nrprocs
   integer(i4), dimension(nrprocs)       , intent(in   ) :: ranks_r, sizes_r, istas_r
   integer(i4),                            intent(inout) :: nsend, nrecv
   integer(i4), dimension(:), allocatable, intent(inout) :: sproc, ssize, sista
   integer(i4), dimension(:), allocatable, intent(inout) :: rproc, rsize, rtype
! local variables
   integer(i4) :: nmap, i, ip, err, sreg(nrprocs), rreg(nprocs)
   integer(i4) :: sstat(mpi_status_size,nrprocs), rstat(mpi_status_size,nprocs)
   integer(i4), dimension(nrprocs) :: sproc1, ssize1, sista1
   integer(i4), dimension(nprocs)  :: rproc1, rsize1, rtype1
!
   ! make send address
   sproc1(:) = ranks_r(:)
   ssize1 = 0
   nmap = size(map)
   do i =1,nmap
     do ip = 0,nrprocs-1
       if (map(i).ge.istas_r(ip+1).and.map(i).le.istas_r(ip+1)+sizes_r(ip+1)-1) then
         ssize1(ip+1) = ssize1(ip+1)+1
         if (ssize1(ip+1).eq.1) then
           sista1(ip+1) = i
         endif
       endif
     enddo
   enddo
!
   ! make send address
   rtype1(:) = -1
   do ip = 1,nrprocs
     call mpi_isend(ssize1(ip),1,mpi_integer,sproc1(ip),111,comm,sreg(ip),err)
   enddo
   if (isrproc) then
     do ip = 1,nprocs
       rproc1(ip) = ip-1
       call mpi_irecv(rsize1(ip),1,mpi_integer,rproc1(ip),111,comm,rreg(ip),err)
     enddo
     call mpi_waitall(nprocs,rreg,rstat,err)
   endif
   call mpi_waitall(nrprocs,sreg,sstat,err)
!
   ! make reduced address
   nsend = count(ssize1.gt.0)
   nrecv = count(rsize1.gt.0)
   allocate(sproc(nsend),ssize(nsend),sista(nsend))
   allocate(rproc(nrecv),rsize(nrecv),rtype(nrecv))
   i = 1
   do ip = 1,nrprocs
     if (ssize1(ip).gt.0) then
       sproc(i) = sproc1(ip)
       ssize(i) = ssize1(ip)
       sista(i) = sista1(ip)
       i = i+1
     endif
   enddo
   if (isrproc) then
     i = 1
     do ip = 1,nprocs
       if (rsize1(ip).gt.0) then
         rproc(i) = rproc1(ip)
         rsize(i) = rsize1(ip)
         rtype(i) = rtype1(ip)
         i = i+1
       endif
     enddo
   endif
   ! map

!
   return
   end subroutine make_comm_address
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_address(scomm,nprocs,rank,rcomm,isrproc,nsend,nrecv,sproc,sista,ssize,rproc,rtype,rsize)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),                   intent(in   ) :: scomm, nprocs, rank, rcomm, nsend, nrecv
   logical(l4),                   intent(in   ) :: isrproc
   integer(i4), dimension(nsend), intent(in   ) :: sproc, sista, ssize
   integer(i4), dimension(nrecv), intent(in   ) :: rproc, rtype, rsize
! local variables
   integer(i4) :: i, ip, err
!
   ! print send cycle
   if (rank.eq.0) write(*,'(a)')         'log)  info.(send address)'
   do ip = 0,nprocs-1
     if (ip.eq.rank) then
       write(*,'(a,i3)') 'log)   rank   = ',rank
       do i = 1,nsend
         write(*,'(a,i3,i3,i3)') 'log)    p, s, i = ', sproc(i), ssize(i), sista(i)
       enddo
     endif
     call mpi_barrier(scomm,err)
   enddo
   call flush()
   call mpi_barrier(scomm,err)
!
   ! print recv cycle
   if (rank.eq.0) write(*,'(a)')         'log)  info.(recv address)'
   do ip = 0,nprocs-1
     if (ip.eq.rank) then
       if (isrproc) then
         write(*,'(a,i3)') 'log)   rank   = ',rank
         do i = 1,nrecv
           write(*,'(a,i3,i3,i3)') 'log)    p, s, i = ', rproc(i), rsize(i), rtype(i)
         enddo
       endif
     endif
     call mpi_barrier(scomm,err)
   enddo
   call flush()
   call mpi_barrier(scomm,err)
!
   return
   end subroutine print_address
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module utilities
!-------------------------------------------------------------------------------
