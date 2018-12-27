!-------------------------------------------------------------------------------
   module par_env
!-------------------------------------------------------------------------------
!
!  abstract :
!    Parallel(MPI derived type & topology) module for framework
!
!  history log :
!    2013-12-10   junghan kim    First written
!    2013-03-19   junghan kim    Adapt a new Parallel module
!
!-------------------------------------------------------------------------------
   use parallel
   use parallel, only: i4=>i4_par, r8=>r8_par, l4=>l4_par
!
   include 'mpif.h'
!
!   integer(i4), parameter, public :: par_cart_maxdims      = 5
   integer(i4), parameter, public :: par_cart_ndims        = 2
   integer(i4), parameter, public :: par_cart_neighbors    = 4
   integer(i4), parameter, public :: par_cart_up           = 1
   integer(i4), parameter, public :: par_cart_down         = 2
   integer(i4), parameter, public :: par_cart_left         = 3
   integer(i4), parameter, public :: par_cart_right        = 3
   integer(i4), parameter, public :: max_graph_indexeds    = 5

   type, public :: parallel_cart_t
     integer(i4)                               :: nprocs
     integer(i4)                               :: iproc
     integer(i4)                               :: comm
     integer(i4)                               :: ndims=par_cart_ndims
     integer(i4), dimension(par_cart_ndims)    :: dims
     logical(l4), dimension(par_cart_ndims)    :: period
     integer(i4), dimension(par_cart_ndims)    :: cartrank
     integer(i4), dimension(par_cart_neighbors):: neighbor
   end type

   type, public :: parallel_graph_t
     integer(i4)                               :: iproc
     integer(i4)                               :: comm
     integer(i4)                               :: nnodes
     integer(i4)                               :: nneighbors
     integer(i4), dimension(max_graph_indexeds):: neighbor
   end type
!
! internal objects
!
   type(parallel_cart_t),  target, public :: kim_par_cart
   type(parallel_graph_t), target, public :: kim_par_graph
!
   private
!
   public :: par_type_commit, par_type_free, par_type_contiguous, par_type_vector
   public :: par_type_indexed, par_type_struct, par_type_subarray
   public :: par_type_extent, par_type_size
   public :: par_cart_create, par_graph_create
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_type_commit(utype, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Commit the MPI derived type
!
!  history log :
!    2013-01-10   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(inout)           :: utype
   integer(i4), intent(  out), optional :: ierr
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
   call mpi_type_commit(utype, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_type_commit
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_type_free(utype, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Free the MPI derived type
!
!  history log :
!    2013-01-10   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   )           :: utype
   integer(i4), intent(  out), optional :: ierr
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
   call mpi_type_free(utype, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_type_free
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_type_contiguous(ncounts, otype, utype, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Creates a contiguous datatype
!
!  history log :
!    2013-01-10   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   )           :: ncounts, otype
   integer(i4), intent(  out)           :: utype
   integer(i4), intent(  out), optional :: ierr
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
   call mpi_type_contiguous(ncounts, otype, utype, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_type_contiguous 
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_type_vector(ncounts, blocklength, stride, otype, utype, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Creates a vector datatype
!
!  history log :
!    2013-01-10   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   )           :: ncounts, blocklength
   integer(i4), intent(in   )           :: stride, otype
   integer(i4), intent(  out)           :: utype
   integer(i4), intent(  out), optional :: ierr
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------

   call mpi_type_vector(ncounts, blocklength, stride, otype, utype, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_type_vector
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_type_indexed(ncounts, array_blocklength, array_disp, otype,  &
                               utype, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Creates a indexed datatype
!
!  history log :
!    2013-01-10   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   )            :: ncounts
   integer(i4), dimension(:), intent(in   )            :: array_blocklength
   integer(i4), dimension(:), intent(in   )            :: array_disp
   integer(i4),               intent(in   )            :: otype
   integer(i4),               intent(  out)            :: utype
   integer(i4),               intent(  out),  optional :: ierr
   integer(i4)                                         :: l_ierr
!-------------------------------------------------------------------------------
   call mpi_type_indexed(ncounts, array_blocklength, array_disp, otype,        &
                         utype, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_type_indexed
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_type_struct(ncounts, array_blocklength, array_disp,          &
                              array_otype, utype, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Creates a structed datatype
!
!  history log :
!    2013-01-10   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   )           :: ncounts
   integer(i4), dimension(:), intent(in   )           :: array_blocklength
   integer(i4), dimension(:), intent(in   )           :: array_disp
   integer(i4), dimension(:), intent(in   )           :: array_otype
   integer(i4),               intent(  out)           :: utype
   integer(i4),               intent(  out), optional :: ierr
   integer(i4)                                        :: l_ierr
!-------------------------------------------------------------------------------
   call mpi_type_struct(ncounts, array_blocklength, array_disp, array_otype,   &
                        utype, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_type_struct
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_type_subarray(ndims, array_size, array_subsize, array_start,&
                                otype, utype, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Creates a subarray datatype
!
!  history log :
!    2013-01-10   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   )           :: ndims
   integer(i4), dimension(:), intent(in   )           :: array_size
   integer(i4), dimension(:), intent(in   )           :: array_subsize
   integer(i4), dimension(:), intent(in   )           :: array_start
   integer(i4),               intent(in   )           :: otype
   integer(i4),               intent(  out)           :: utype
   integer(i4),               intent(  out), optional :: ierr
   integer(i4)                                        :: l_order
   integer(i4)                                        :: l_ierr
!-------------------------------------------------------------------------------
   l_order = mpi_order_fortran

   call mpi_type_create_subarray(ndims, array_size, array_subsize,array_start, &
                                 l_order, otype, utype, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_type_subarray
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_type_extent(utype, extent, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Compute the extent of derived type (not include spaces)
!
!  history log :
!    2013-01-10   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   )           :: utype
   integer(i4), intent(  out)           :: extent
   integer(i4), intent(  out), optional :: ierr
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------
   call mpi_type_extent(utype, extent, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_type_extent
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_type_size(utype, nsize, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Compute the extent of derived type (include spaces)
!
!  history log :
!    2013-01-10   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   )           :: utype
   integer(i4), intent(  out)           :: nsize
   integer(i4), intent(  out), optional :: ierr
   integer(i4)                          :: l_ierr
!-------------------------------------------------------------------------------

   call mpi_type_size(utype, nsize, l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_type_size
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine par_cart_create(dimsize, period, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Make a communicator for cartesian coordinator
!
!  history log :
!    2013-01-10   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   !integer(i4), intent(in)              :: comm, ndims
   integer(i4), dimension(2), intent(in   )            :: dimsize
   logical(l4), dimension(2), intent(in   )            :: period
   integer(i4),               intent(  out),  optional :: ierr
   logical(l4)                                         :: l_reorder
   integer(i4)                                         :: l_newcomm
   integer(i4)                                         :: l_ierr
!-------------------------------------------------------------------------------
   l_reorder = .false.

   call mpi_cart_create(kim_par%comm, kim_par_cart%ndims, dimsize, period,     &
                        l_reorder, l_newcomm, l_ierr)

   kim_par_cart%nprocs    = kim_par%nprocs
   kim_par_cart%iproc     = kim_par%iproc
   kim_par_cart%comm      = l_newcomm
   kim_par_cart%dims(:)   = dimsize(:)
   kim_par_cart%period(:) = period(:)

   call mpi_cart_coords(l_newcomm, kim_par_cart%iproc, kim_par_cart%ndims,     &
                        kim_par_cart%cartrank(:), l_ierr)
   call mpi_cart_shift(l_newcomm, 0, 1, kim_par_cart%neighbor(par_cart_left),  &
                       kim_par_cart%neighbor(par_cart_right), l_ierr)
   call mpi_cart_shift(l_newcomm, 1, 1, kim_par_cart%neighbor(par_cart_down),  &
                       kim_par_cart%neighbor(par_cart_up), l_ierr)

   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_cart_create
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
    SUBROUTINE Par_Graph_Create(indexes, edges, ierr)
!-------------------------------------------------------------------------------
!
!  abstract :
!    Make a communicator for graph topology
!
!  history log :
!    2013-01-10   junghan kim    First written
!
!-------------------------------------------------------------------------------
   implicit none
!
   !integer(i4), intent(in)              :: nodes
   integer(i4), dimension(:), intent(in   )           :: indexes
   integer(i4), dimension(:), intent(in   )           :: edges
   integer(i4),               intent(  out), optional :: ierr
   logical(l4)                                        :: l_reorder
   integer(i4)                                        :: l_newcomm
   integer(i4)                                        :: l_ierr, l_i, l_start
!-------------------------------------------------------------------------------
   l_reorder = .false.

   call mpi_graph_create(kim_par%comm, kim_par%nprocs, indexes, edges,         &
                         l_reorder, l_newcomm, l_ierr)
    
   kim_par_graph%iproc      = kim_par%iproc
   kim_par_graph%comm       = l_newcomm
   kim_par_graph%nnodes     = kim_par%nprocs
   kim_par_graph%nneighbors = indexes(kim_par_graph%iproc + 1)
!
   if (kim_par_graph%iproc .ne. 0) then
     kim_par_graph%nneighbors = kim_par_graph%nneighbors -                     &
                                indexes(kim_par_graph%iproc)
   endif
!
   kim_par_graph%neighbor(:) = par_proc_null
   l_start = indexes(kim_par_graph%iproc + 1) - kim_par_graph%nneighbors
!
   do l_i = 1, kim_par_graph%nneighbors
     kim_par_graph%neighbor(l_i) = edges(l_start + l_i)
   enddo
!
   if (present(ierr)) ierr = l_ierr
   return

   end subroutine par_graph_create
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module par_env
!-------------------------------------------------------------------------------
