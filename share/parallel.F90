!#include <KIM.h>
!-------------------------------------------------------------------------------
   module parallel
!-------------------------------------------------------------------------------
!
!  abstract : base parallel module for framework
!
!  history log :
!    2013-12-10  junghan kim   first written
!    2013-03-08  junghan kim   add thread & hybrid(thread+process) 
!                              programming model
!    2013-03-08  insun song    comments corrected and codes cleaned up
!
!-------------------------------------------------------------------------------
   include 'mpif.h'
!
   private
!
! kind parameter
!
   integer(4), parameter, public        :: i4_par      = 4
   integer(4), parameter, public        :: i8_par      = 8
   integer(4), parameter, public        :: l4_par      = 4
   integer(4), parameter, public        :: r4_par      = 4
   integer(4), parameter, public        :: r8_par      = 8
   integer(4), parameter, public        :: r16_par     = 8
   integer(4), parameter, public        :: len_string  = 255

   integer(4), parameter :: i4          = i4_par
   integer(4), parameter :: i8          = i8_par
   integer(4), parameter :: l4          = l4_par
   integer(4), parameter :: r4          = r4_par
   integer(4), parameter :: r8          = r8_par
   integer(4), parameter :: r16         = r16_par
!
! master process rank
!
   integer(4), parameter, public        :: masterprocrank = 0
   integer(4), parameter, public        :: maxintercomm   = 10
!
! mpi parameters for communicator
!
   integer(i4), parameter, public :: par_comm_world  = mpi_comm_world
!
! MPI parameters for communication
!
   integer(i4), parameter, public :: par_proc_null   = mpi_proc_null
   integer(i4), parameter, public :: par_status_size = mpi_status_size
   integer(i4), parameter, public :: par_tag_ub      = mpi_tag_ub
   integer(i4), parameter, public :: par_any_source  = mpi_any_source
   integer(i4), parameter, public :: par_any_tag     = mpi_any_tag
   integer(i4), parameter, public :: par_tag         = mpi_tag
   integer(i4), parameter, public :: par_source      = mpi_source
   integer(i4), parameter, public :: par_sum         = mpi_sum
   integer(i4), parameter, public :: par_prod        = mpi_prod
   integer(i4), parameter, public :: par_min         = mpi_min
   integer(i4), parameter, public :: par_max         = mpi_max
   integer(i4), parameter, public :: par_minloc      = mpi_minloc
   integer(i4), parameter, public :: par_maxloc      = mpi_maxloc
   integer(i4), parameter, public :: par_land        = mpi_land
   integer(i4), parameter, public :: par_lor         = mpi_lor
   integer(i4), parameter, public :: par_lxor        = mpi_lxor
   integer(i4), parameter, public :: par_band        = mpi_band
   integer(i4), parameter, public :: par_bor         = mpi_bor
   integer(i4), parameter, public :: par_bxor        = mpi_bxor
!
! MPI parameters for parallel environment
!
   integer(i4), parameter, public :: par_address_kind      = mpi_address_kind
   integer(i4), parameter, public :: par_error             = mpi_error
   integer(i4), parameter, public :: par_success           = mpi_success
   integer(i4), parameter, public :: par_info_null         = mpi_info_null
   integer(i4), parameter, public :: par_request_null      = mpi_request_null
   integer(i4), parameter, public :: par_thread_multiple   = mpi_thread_multiple
   integer(i4), parameter, public :: par_thread_single     = mpi_thread_single
   integer(i4), parameter, public :: par_thread_funneled   = mpi_thread_funneled
   integer(i4), parameter, public :: par_thread_serialized =                    &
                                     mpi_thread_serialized
   integer(i4), parameter, public :: par_offset_kind       = mpi_offset_kind
   integer(i4), parameter, public :: par_order_fortran     = mpi_order_fortran
   integer(i4), parameter, public :: par_mode_wronly       = mpi_mode_wronly
   integer(i4), parameter, public :: par_mode_rdonly       = mpi_mode_rdonly
   integer(i4), parameter, public :: par_mode_create       = mpi_mode_create
!
! MPI parameters for data type
!
   integer(i4), parameter, public :: par_byte             = mpi_byte
   integer(i4), parameter, public :: par_character        = mpi_character
   integer(i4), parameter, public :: par_logical          = mpi_logical
   integer(i4), parameter, public :: par_integer          = mpi_integer
   integer(i4), parameter, public :: par_real             = mpi_real
   integer(i4), parameter, public :: par_double_precision = mpi_double_precision
   integer(i4), parameter, public :: par_real4            = mpi_real4
   integer(i4), parameter, public :: par_real8            = mpi_real8
!
! Topology
!
  integer(i4), parameter, public :: par_cart              = mpi_cart
  integer(i4), parameter, public :: par_graph             = mpi_graph
  integer(i4), parameter, public :: par_undefined         = mpi_undefined
!
! Type definition for MPI and hybrid parallelism
!
   type, public :: parallel_t
     sequence
     integer(i4)                           :: nprocs
     integer(i4)                           :: iproc
     integer(i4)                           :: root
     logical(l4)                           :: ismaster
     integer(i4)                           :: comm
     integer(i4)                           :: nintercomm
     integer(i4), dimension(maxintercomm)  :: intercomm
   end type parallel_t
!
   type, public :: hybrid_t
     sequence
     type(parallel_t)      :: par
     integer(i4) :: ithr
     integer(i4) :: nthrs
     logical(l4) :: ismasterthr
   end type hybrid_t
!
! Internal Object
! ---------------
!
  type(parallel_t), target, public :: kim_par
  type(hybrid_t),   target, public :: kim_hybrid
!
! Etc
! ---
!
  integer(i4), parameter,   public :: nrepro_vars=10
  integer(i4), parameter,   public :: ordered = 1
  integer(i4), parameter,   public :: fast = 2
!
  integer(i4),              public :: nthreads
  integer(i4),              public :: iam
  integer(i4),              public :: ncompoints, npackpoints
  integer(i4),              public :: nmpispernode
  real   (r8), allocatable, public :: global_shared_buf(:,:)
  real   (r8),              public :: global_shared_sum(nrepro_vars)
  logical(l4),              public :: partitionfornodes, partitionforframe
!
#ifdef _OPENMP
  interface omp_get_thread_num
    integer function omp_get_thread_num()
    end function omp_get_thread_num
  end interface
  interface omp_in_parallel
    logical function omp_in_parallel()
    end function omp_in_parallel
  end interface
  interface omp_set_num_threads
    subroutine omp_set_num_threads(nthreads_in)
    integer nthreads_in
    end subroutine omp_set_num_threads
  end interface
  interface omp_get_num_threads
    integer function omp_get_num_threads()
    end function omp_get_num_threads
  end interface
  interface omp_get_max_threads
    integer function omp_get_max_threads()
    end function omp_get_max_threads
  end interface
#endif
!
   public :: par_initialize, par_finalize, par_create, setintercomm,            &
             getintercommindex
   public :: par_abort, par_barrier
   public :: omp_get_thread_num, omp_in_parallel
   public :: omp_set_num_threads, omp_get_max_threads, omp_get_num_threads
   public :: hybrid_initialize, hybrid_create, hybrid_setup
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_initialize(comm, ierr)
!-------------------------------------------------------------------------------
!
!  abstract : Initialize the ParBase module & MPI execution environment
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), optional, intent(in   ) :: comm
   integer(i4), optional, intent(  out) :: ierr
!
! Local variables
!
   integer(i4)                                          :: l_ierr
   logical(l4)                                          :: l_initted
   character(len=mpi_max_processor_name)                :: myprocname
   character(len=mpi_max_processor_name), allocatable   :: allprocnames(:)
   integer(i4)                                          :: namelen, i, tmp
!-------------------------------------------------------------------------------
   call mpi_initialized(l_initted, l_ierr)
!
   if (.not. l_initted) then
     call mpi_init(l_ierr)
   endif
!
   if (present(comm)) then
     kim_par%comm = comm
   else
     kim_par%comm = par_comm_world
   endif
!
   call mpi_comm_size(kim_par%comm, kim_par%nprocs, l_ierr)
   call mpi_comm_rank(kim_par%comm, kim_par%iproc, l_ierr)
!
   kim_par%root = masterprocrank
!
   if (kim_par%iproc .eq. kim_par%root) then
     kim_par%ismaster = .true.
   else
     kim_par%ismaster = .false.
   endif
!
   kim_par%nintercomm = 0
!
   iam = kim_par%iproc + 1
!
! determines where the current mpi process is running.
! then uses the information to determine the number of
! mpi processes per node
!
   myprocname(:) = ''
   call mpi_get_processor_name(myprocname,namelen,l_ierr)
!
   allocate(allprocnames(kim_par%nprocs))
   do i = 1, kim_par%nprocs
     allprocnames(i)(:) =  ''
   enddo
!
! collects all the machine names
!
   call mpi_allgather(myprocname, mpi_max_processor_name, mpi_character,       &
   allprocnames, mpi_max_processor_name, mpi_character, kim_par%comm, l_ierr)
!
! calculates how many other mpi processes are in a node
!
   nmpispernode = 0
   do i = 1, kim_par%nprocs
     if (trim(adjustl(myprocname)) .eq. trim(adjustl(allprocnames(i)))) then
       nmpispernode = nmpispernode + 1
     endif
   enddo
!
! verifies every process agrees the nmpispernode
! if not, multi-level partitioning is disabled.
!
   call mpi_allreduce(nmpispernode, tmp, 1, par_integer, par_band,             &
                      kim_par%comm, l_ierr)
   if (tmp .ne. nmpispernode) then
!    write(kim_iu_log,*)'initmp:  disagrement accross nodes for nmpispernode'
     nmpispernode = 1
     partitionfornodes= .false.
   else
     partitionfornodes= .true.
   endif
!
   deallocate(allprocnames)
   if (present(ierr)) ierr = l_ierr
!
   return
   end subroutine par_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_finalize(message, ierr)
!-------------------------------------------------------------------------------
!
!  abstarct :
!    Finalize the ParBase module & MPI execution environment
!
!  history log :
!    2013-01-10  junghan kim   first written
!    2015-02-18  insun song    Minor bug for handling optional 
!                              variable message is fixed
!
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),         optional :: message
   integer(i4), intent(out), optional :: ierr
   integer(i4) :: l_ierr
!-------------------------------------------------------------------------------
   if (present(message)) then
     if (iam .eq. 1) then
       write(*,*) message
     endif
   endif
!
   call mpi_finalize(l_ierr)
!
   if (present(ierr)) ierr = l_ierr
!
   return
   end subroutine par_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function par_create(comm, ierr)
!
!  abstract :
!    Make new parallel_t object & Initialize that object
!
!  history log :
!    2013-01-10  junghan kim   first written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   )           :: comm
   integer(i4), intent(  out), optional :: ierr
   integer(i4)                          :: l_ierr
   type(parallel_t)                     :: par_create
!-------------------------------------------------------------------------------
   par_create%comm = comm
!
   call mpi_comm_size(comm, par_create%nprocs, l_ierr)
   call mpi_comm_rank(comm, par_create%iproc, l_ierr)
!
   par_create%root = masterprocrank
!
   if (par_create%iproc == par_create%root) then
     par_create%ismaster = .true.
   else
     par_create%ismaster = .false.
   endif
!
   par_create%nintercomm = 0
!
   if (present(ierr)) ierr = l_ierr
!
   return
   end function par_create
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine setintercomm(par, intercomm, ierr)
!-------------------------------------------------------------------------------
!
!  abstarct : 
!    Add a interCommunicator
!
!  history log :
!    2013-01-10  junghan kim   first written
!
!-------------------------------------------------------------------------------
   implicit none
!
   type(parallel_t), intent(inout), optional, target :: par
   integer(i4),      intent(in   )                   :: intercomm
   integer(i4),      intent(  out), optional         :: ierr
   type(parallel_t), pointer :: l_par
   integer(i4)               :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(par)) then
     l_par => par
   else
     l_par => kim_par
   endif
!
   if (l_par%nintercomm <= maxintercomm) then
     l_par%nInterComm = l_par%ninterComm + 1
     l_par%interComm = interComm
   else
     call par_abort(ierr = l_ierr)
   endif
!
   nullify(l_par)
!
   if (present(ierr)) ierr = l_ierr
!
   return
   end subroutine setintercomm
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_abort(message, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Abort the ParBase module 
!
!  history log :
!    2013-01-10  junghan kim   first written
!    2015-02-18  insun song    Minor bug fix for handling optional variable
!
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),         optional :: message
   integer(i4), intent(out), optional :: ierr
   integer(i4) :: l_info
   integer(i4) :: l_ierr
!-------------------------------------------------------------------------------
!
   if (present(message)) then
     write(*,*) ' aborting with error: ', message
   endif
!
   call mpi_abort(par_comm_world, l_info, l_ierr)
   call mpi_finalize(l_ierr)
!
   if (present(ierr)) ierr = l_ierr
!
   return
   end subroutine par_abort
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_barrier(par, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Blocks until all processes in the communicator have reached this routine
! 
!  history log :
!    2013-01-10  junghan kim   first written
!
!-------------------------------------------------------------------------------
   implicit none
!
   type(parallel_t), intent(in),  optional, target :: par
   integer(i4),      intent(out), optional         :: ierr
!
   type(parallel_t), pointer :: l_par
   integer(i4)               :: l_ierr
!-------------------------------------------------------------------------------
   nullify(l_par)
!
   if (present(par)) then
     l_par => par
   else
     l_par => kim_par
   endif
!
   call mpi_barrier(l_par%comm, l_ierr)
!
   nullify(l_par)
   if (present(ierr)) ierr = l_ierr
!
   return
   end subroutine par_barrier
!-------------------------------------------------------------------------------
!
!  abstract :
!    FUNCTIONS for single thread program
!
!  history log :
!    2013-01-10  junghan kim   first written
!
!-------------------------------------------------------------------------------
#ifndef _OPENMP
!
  function omp_get_thread_num() result(ithr)
    integer(i4)            ::  ithr
    ithr=0
  end function omp_get_thread_num
!
  function omp_get_num_threads() result(ithr)
    integer(i4)            ::  ithr
    ithr=1
  end function omp_get_num_threads
!
  function omp_in_parallel() result(ans)
    logical(l4)            :: ans
    ans=.false.
  end function omp_in_parallel
!
  subroutine omp_set_num_threads(nthrs)
    integer(i4), intent(in)::  nthrs
    nthreads=1
  end subroutine omp_set_num_threads
!
  integer function omp_get_max_threads()
    omp_get_max_threads=1
  end function omp_get_max_threads
!
#endif
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine hybrid_initialize(par, hybrid, nthr, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Make hybrid_t object & Initialize that object
!
!  history log :
!    2013-03-08  junghan kim   first written
!    2014-06-18  junghan kim   Bug fix
!
!-------------------------------------------------------------------------------
   implicit none
!
   type(parallel_t), intent(inout), optional, target :: par
   type(hybrid_t),   intent(inout), optional, target :: hybrid
   integer(i4),      intent(in   ), optional         :: nthr
   integer(i4),      intent(  out), optional         :: ierr

   type(parallel_t), pointer :: l_par
   type(hybrid_t),   pointer :: l_hybrid
   integer(i4)               :: l_ierr
!-------------------------------------------------------------------------------
   nullify(l_par)
   nullify(l_hybrid)
!
   if (present(par)) then
     l_par => par
   else
     l_par => kim_par
   endif
   if (present(hybrid)) then
     l_hybrid => hybrid
   else
     l_hybrid => kim_hybrid
   endif
!
   if (present(nthr)) then
     call omp_set_num_threads(nthr)
     nthreads = nthr
   else
!$OMP PARALLEL
#ifdef _OPENMP
   nthreads = omp_get_num_threads()
#else
   nthreads = 1
#endif
!$OMP END PARALLEL
   endif
!
   l_hybrid%par         = l_par
   l_hybrid%ithr        = omp_get_thread_num()
   l_hybrid%nthrs       = nthreads
   l_hybrid%ismasterthr = (l_par%ismaster.and.(l_hybrid%ithr.eq.0))
!
   nullify(l_par)
   nullify(l_hybrid)
!
   if (present(ierr)) ierr = l_ierr
!
   return
   end subroutine hybrid_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function hybrid_create(par, nthr, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Create hybrid_t object and initialize that object
! 
!  history log :
!    2014-06-13  junghan kim   first written
!
!-------------------------------------------------------------------------------
   implicit none
!
   type(parallel_t), intent(in   ), optional, target :: par
   integer(i4),      intent(in   ), optional         :: nthr
   integer(i4),      intent(  out), optional         :: ierr
   type(hybrid_t)                                    :: hybrid_create
!
   type(parallel_t), pointer :: l_par
   integer(i4)               :: l_nthr, l_ierr
!-------------------------------------------------------------------------------
   nullify(l_par)
!
   if (present(par)) then
     l_par => par
   else
     l_par => kim_par
   endif
!
   if (present(nthr)) then
     call omp_set_num_threads(nthr)
     l_nthr = nthr
   else
     l_nthr = omp_get_num_threads()
   endif
!
   hybrid_create%par         = l_par
   hybrid_create%ithr        = omp_get_thread_num()
   hybrid_create%nthrs       = l_nthr
   hybrid_create%ismasterthr = (l_par%ismaster.and.(hybrid_create%ithr.eq.0))
!
   nullify(l_par)
   if (present(ierr)) ierr = l_ierr
!
   return
   end function hybrid_create
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine hybrid_setup(hybrid, nthrs, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Make hybrid_t object & Initialize that object
!
!  history log :
!    2013-03-08  junghan kim   first written
!
!-------------------------------------------------------------------------------
   implicit none
!
   type(hybrid_t), intent(inout), optional, target :: hybrid
   integer(i4),    intent(in   ), optional, target :: nthrs
   integer(i4),    intent(  out), optional         :: ierr
!
   type(hybrid_t), pointer :: l_hybrid
   integer(i4)             :: l_nthrs
   integer(i4)             :: l_ierr
!-------------------------------------------------------------------------------
   nullify(l_hybrid)
!
   if (present(hybrid)) then
     l_hybrid => hybrid
   else
     l_hybrid => kim_hybrid
   endif
   if (present(nthrs)) then
     l_nthrs = nthrs
   else
     l_nthrs = nthreads
   endif
!
   l_hybrid%ithr        = omp_get_thread_num()
   l_hybrid%nthrs       = l_nthrs
   l_hybrid%ismasterthr = (l_hybrid%par%ismaster.and.(l_hybrid%ithr.eq.0))
!
   nullify(l_hybrid)
   if (present(ierr)) ierr = l_ierr
!
   return
   end subroutine hybrid_setup
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module parallel
!-------------------------------------------------------------------------------
