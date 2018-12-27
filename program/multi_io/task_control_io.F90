!-------------------------------------------------------------------------------
   module task_control_io  ! USER
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    2017-05-22  junghan kim    initial setup
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds,    only: i4, l4, r4, r8
   use multi_io, only: read_file, write_file, copy_file ! USER
!-------------------------------------------------------------------------------
   implicit none
   include 'mpif.h'
!
   private
!
   ! signal: -2(post), -1(exit), 0,1,2,3,...(user signal)
   integer(i4), parameter :: sig_exit_=-1
   integer(i4), parameter :: sig_post_=-2
   ! profile
   integer(i4), parameter :: nwtimes_ = 4   ! run, post, comm
   integer(i4), parameter :: w_all_ = 1, w_task_ = 2, w_post_= 3, w_comm_= 4
!
   type task_t
     integer(i4) :: ntasks
     integer(i4) :: comm, nprocs, rank, root
     logical(l4) :: ismaster, isserial
     ! profile
     integer(i4) :: nsend, nrecv, nsend_r, nrecv_r
     real(r8)    :: wtime(nwtimes_)
     ! task option
     integer(i4) :: wtask ! -2(random), -1(1 second), 0,1,2,3,...(user task)
   end type task_t
!
   public :: task_t, task_initialize, task_driver, task_driver_sync, task_profile
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine task_initialize(task, ntasks, comm, nprocs, rank, isserial, wtask)
!-------------------------------------------------------------------------------
   implicit none
!
   type(task_t),          intent(inout) :: task
   integer(i4),           intent(in   ) :: ntasks, comm, nprocs, rank
   logical(l4), optional, intent(in   ) :: isserial
   integer(i4), optional, intent(in   ) :: wtask
! local variables
!
   task%ntasks    = ntasks
   task%comm     = comm
   task%nprocs   = nprocs
   task%rank     = rank
   task%root     = 0
   task%ismaster = .false.
   if (task%rank.eq.task%root) task%ismaster = .true.
!
   if (present(isserial)) then
     task%isserial = isserial
   else
     task%isserial = .false.
   endif
!
   task%wtime(:) = 0.0_r8
   task%nsend    = 0
   task%nrecv    = 0
   task%nsend_r  = 0
   task%nrecv_r  = 0
!
   if (present(wtask)) then
     task%wtask = wtask
   else
     task%wtask = 0
   endif
!
   return
   end subroutine task_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine task_driver(task)
!-------------------------------------------------------------------------------
   implicit none
!
   type(task_t), intent(inout) :: task
! local variables
   integer(i4) :: itask, irank, sig, i, err
   real(r8)    :: st, et
!
   st = mpi_wtime(err)
!
   if (task%nprocs.gt.1.and..not.task%isserial) then
   !
     if (task%ismaster) then
     !
       do itask = 1,task%ntasks
         call recv_ready(task,irank)
         call send_signal(task,irank,itask)
       enddo
       !
       call recv_ready_bcast_signal(task,sig_exit_)
     !
     else
     !
       do while (.true.)
         call send_ready(task)
         call recv_signal(task,sig)
         if (sig.eq.sig_exit_) then
           exit
         else
           call run_task(task,sig)
         endif
       enddo
     !
     endif
   !
   else
   !
     if (task%ismaster) then
       do itask = 1,task%ntasks
         call run_task(task,itask)
       enddo
     endif
   !
   endif
!
   et = mpi_wtime(err)
   task%wtime(w_all_) = task%wtime(w_all_)+(et-st)
!
   return
   end subroutine task_driver
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine task_driver_sync(task)
!-------------------------------------------------------------------------------
   implicit none
!
   type(task_t), intent(inout) :: task
! local variables
   integer(i4) :: itask, irank, sig, i, err
   real(r8)    :: st, et
!
   st = mpi_wtime(err)
!
   if (task%nprocs.gt.1.and..not.task%isserial) then
   !
     if (task%ismaster) then
     !
       irank = 0
       do itask = 1,task%ntasks
         irank = irank+1
         call send_signal(task,irank,itask)
         if (irank.eq.task%nprocs-1) then ! full
           do i = 1,irank
             call bcast_signal(task,sig_post_)
             call post_task(task,i)
           enddo
           irank = 0
         endif
       enddo
       !
       do i = 1,irank
         call bcast_signal(task,sig_post_)
         call post_task(task,i)
       enddo
       !
       call bcast_signal(task,sig_exit_)
     !
     else
     !
       do while (.true.)
         call recv_signal(task,sig)
         if      (sig.eq.sig_exit_) then
           exit
         else if (sig.eq.sig_post_) then
           call post_task(task)
         else
           call run_task(task,sig)
         endif
       enddo
     !
     endif
   !
   else
   !
     if (task%ismaster) then
       do itask = 1,task%ntasks
         call run_task(task,itask)
       enddo
     endif
   !
   endif
!
   et = mpi_wtime(err)
   task%wtime(w_all_) = task%wtime(w_all_)+(et-st)
!
   return
   end subroutine task_driver_sync
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine send_signal(task, irank, signal)
!-------------------------------------------------------------------------------
   implicit none
!
   type(task_t), intent(inout) :: task
   integer(i4),  intent(in   ) :: irank, signal
! local variables
   real(r8)    :: st, et
   integer(i4) :: err
!
   st = mpi_wtime(err)
!
   print '(a,i4.4,a,i4.4)',' the signal(',signal,') was sent to rank: ',irank;call flush()
   call mpi_send(signal,1,mpi_integer,irank,101,task%comm,err)
!
   et = mpi_wtime(err)
   task%wtime(w_comm_) = task%wtime(w_comm_)+(et-st)
   task%nsend = task%nsend+1
!
   return
   end subroutine send_signal
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine recv_signal(task, signal)
!-------------------------------------------------------------------------------
   implicit none
!
   type(task_t), intent(inout) :: task
   integer(i4),  intent(inout) :: signal
! local variables
   real(r8)    :: st, et
   integer(i4) :: err, stat(mpi_status_size)
!
   st = mpi_wtime(err)
!
   !print '(a,i4.4,a,i4.4)',' rank ',task%rank,' received the signal ',signal;call flush()
   call mpi_recv(signal,1,mpi_integer,task%root,mpi_any_tag,task%comm,stat,err)
!
   et = mpi_wtime(err)
   task%wtime(w_comm_) = task%wtime(w_comm_)+(et-st)
   task%nrecv = task%nrecv+1
!
   return
   end subroutine recv_signal
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine send_ready(task)
!-------------------------------------------------------------------------------
   implicit none
!
   type(task_t), intent(inout) :: task
! local variables
   real(r8)    :: st, et
   integer(i4) :: err
!
   st = mpi_wtime(err)
!
   call mpi_send(task%rank,1,mpi_integer,task%root,101,task%comm,err)
!
   et = mpi_wtime(err)
   task%wtime(w_comm_) = task%wtime(w_comm_)+(et-st)
   task%nsend_r = task%nsend_r+1
!
   return
   end subroutine send_ready
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine recv_ready(task, irank)
!-------------------------------------------------------------------------------
   implicit none
!
   type(task_t), intent(inout) :: task
   integer(i4),  intent(inout) :: irank
! local variables
   real(r8)    :: st, et
   integer(i4) :: err, stat(mpi_status_size)
!
   st = mpi_wtime(err)
!
   call mpi_recv(irank,1,mpi_integer,mpi_any_source,mpi_any_tag,task%comm,stat,err)
!
   et = mpi_wtime(err)
   task%wtime(w_comm_) = task%wtime(w_comm_)+(et-st)
   task%nrecv_r = task%nrecv_r+1
!
   return
   end subroutine recv_ready
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine bcast_signal(task, signal)
!-------------------------------------------------------------------------------
   implicit none
!
   type(task_t), intent(inout) :: task
   integer(i4),  intent(in   ) :: signal
! local variables
   real(r8)    :: st, et
   integer(i4) :: irank, err
   integer(i4), dimension(:), allocatable :: regs
!
   st = mpi_wtime(err)
!
   allocate(regs(task%nprocs-1))
   do irank = 1,task%nprocs-1
     call mpi_isend(signal,1,mpi_integer,irank,101,task%comm,regs(irank),err)
     task%nsend = task%nsend+1
   enddo
   call mpi_waitall(task%nprocs-1,regs,mpi_status_ignore,err)
   deallocate(regs)
!
   et = mpi_wtime(err)
   task%wtime(w_comm_) = task%wtime(w_comm_)+(et-st)
!
   return
   end subroutine bcast_signal
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine recv_ready_bcast_signal(task, signal)
!-------------------------------------------------------------------------------
   implicit none
!
   type(task_t), intent(inout) :: task
   integer(i4),  intent(in   ) :: signal
! local variables
   real(r8)    :: st, et
   integer(i4) :: i, irank, err
   integer(i4), dimension(:), allocatable :: regs
!
   st = mpi_wtime(err)
!
print*,'step 1'; call flush()
   do i = 1,task%nprocs-1
     call recv_ready(task,irank)
   enddo
print*,'step 2'; call flush()
!
   allocate(regs(task%nprocs-1))
   do irank = 1,task%nprocs-1
     call mpi_isend(signal,1,mpi_integer,irank,101,task%comm,regs(irank),err)
     task%nsend = task%nsend+1
   enddo
   call mpi_waitall(task%nprocs-1,regs,mpi_status_ignore,err)
   deallocate(regs)
print*,'step 3'; call flush()
!
   et = mpi_wtime(err)
   task%wtime(w_comm_) = task%wtime(w_comm_)+(et-st)
!
   return
   end subroutine recv_ready_bcast_signal
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine post_task(task, irank)
!-------------------------------------------------------------------------------
   implicit none
!
   type(task_t),          intent(inout) :: task
   integer(i4), optional, intent(in   ) :: irank
! local variables
   real(r8)    :: st, et
   integer(i4) :: err, rank
!
   if (task%ismaster.and..not.present(irank)) then
     print *, 'some error in post_task'
   endif
!
   st = mpi_wtime(err)
!
   if (present(irank)) then
     rank = irank
   endif
   call mpi_bcast(rank,1,mpi_integer,task%root,task%comm,err)
!
   et = mpi_wtime(err)
   task%wtime(w_comm_) = task%wtime(w_comm_)+(et-st)
!
   st = mpi_wtime(err)
! user code here
!   call mpi_barrier(task%comm,err)
   et = mpi_wtime(err)
   task%wtime(w_post_) = task%wtime(w_post_)+(et-st)
!
   return
   end subroutine post_task
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine run_task(task, signal)
!-------------------------------------------------------------------------------
   implicit none
!
   type(task_t), intent(inout) :: task
   integer(i4),  intent(in   ) :: signal
! local variables
   real(r8)    :: st, et, wait
   integer(i4) :: err, nseed
   integer(i4), allocatable :: seed(:)
!
   st = mpi_wtime(err)
   !print '(a,i4.4,a,i4.4)',' rank ',task%rank,' starts the task: ',signal;call flush()
!
   if      (task%wtask.eq.-1) then
     call sleep(1)
   else if (task%wtask.eq.-2) then
     call random_seed(size=nseed)
     allocate(seed(nseed))
     seed(:) = task%rank
     call random_seed(put=seed)
     call random_number(wait)
     wait = wait+1.0_r8
     call sleep(int(2.0*wait))
     deallocate(seed)
   else
! user code here
     !call write_file(signal) ! USER
     !call read_file(signal) ! USER
     call copy_file(signal) ! USER
   endif
!
   et = mpi_wtime(err)
   task%wtime(w_task_) = task%wtime(w_task_)+(et-st)
!
   return
   end subroutine run_task
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine task_profile(task)
!-------------------------------------------------------------------------------
   implicit none
   type(task_t), intent(inout) :: task
! local vairiables
   integer(i4) :: i, err
   real(r8), dimension(:,:), allocatable :: wtime_all
   real(r8), dimension(:)  , allocatable :: wtime, wtime_min, wtime_max
!
   allocate(wtime_all(nwtimes_,task%nprocs))
   allocate(wtime(nwtimes_),wtime_min(nwtimes_),wtime_max(nwtimes_))
!
   if (task%nprocs.eq.1.or.task%isserial) then
     if (task%ismaster) then
       wtime(:)     = task%wtime(:)
       wtime_min(:) = task%wtime(:)
       wtime_max(:) = task%wtime(:)
     endif
   else
     call mpi_gather(task%wtime,nwtimes_,mpi_double,wtime_all,nwtimes_,         &
                                              mpi_double,task%root,task%comm,err)
     if (task%ismaster) then
       wtime(:) = 0.0_r8
       do i = 2,task%nprocs
         wtime(:) = wtime(:)+wtime_all(:,i)
       enddo
       wtime(:) = wtime(:)/real(task%nprocs-1,r8)
       do i = 1,nwtimes_
         wtime_min(i) = minval(wtime_all(i,2:))
         wtime_max(i) = maxval(wtime_all(i,2:))
       enddo
     endif
   endif
   if (task%ismaster) then
     print '(a)'                     , ' '
     print '(a)'                     , ' * profile :'
     print '(a,a7  ,a,a7  ,a,a7  ,a)', '           ',  'mean',' (',       'min',',',       'max',')'
     print '(a,f7.3,a,f7.3,a,f7.3,a)', '  -  all = ',wtime(w_all_) ,' (',wtime_min(w_all_) ,',',wtime_max(w_all_) ,')'
     print '(a,f7.3,a,f7.3,a,f7.3,a)', '  - task = ',wtime(w_task_),' (',wtime_min(w_task_),',',wtime_max(w_task_),')'
     print '(a,f7.3,a,f7.3,a,f7.3,a)', '  - post = ',wtime(w_post_),' (',wtime_min(w_post_),',',wtime_max(w_post_),')'
     print '(a,f7.3,a,f7.3,a,f7.3,a)', '  - comm = ',wtime(w_comm_),' (',wtime_min(w_comm_),',',wtime_max(w_comm_),')'
     print '(a)'                     , ' '
   endif
   call flush()
!
#if 0
   call mpi_barrier(task%comm,err)
   do i = 0,task%nprocs-1
     if (i.eq.task%rank) then
       print '(a,i3.3)'     , ' rank   : ',task%rank; call flush()
       print '(a,i3,a,i2,a)', ' ready  : send(',task%nsend_r,'), recv(',task%nrecv_r,')'; call flush()
       print '(a,i3,a,i2,a)', ' signal : send(',task%nsend,  '), recv(',task%nrecv,  ')'; call flush()
       call mpi_barrier(task%comm,err)
     endif
   enddo
#endif
!
   deallocate(wtime_all,wtime,wtime_min,wtime_max)
!
   return
   end subroutine task_profile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module task_control_io
!-------------------------------------------------------------------------------
