!-------------------------------------------------------------------------------
!>
!> @brief
!>  - edge module.
!>
!> @date ?????2012
!>  - JH KIM  : First written from the HOMME and was modified for KIAPSGM framework.
!> @date 30JAN2015
!>  - JH KIM  : Added to "initEdgeBuffer", "ghostVpack_new", "ghostVunpack_new", "ghostVpack2D", "ghostVunpack2D" subroutines. (CAM-SE 5.3)
!-------------------------------------------------------------------------------
#include "KIM.h"
!-------------------------------------------------------------------------------
   module edge
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  your   name    code comment
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kiapsbase, only: int_kind=>kim_int_kind, log_kind=>kim_log_kind, real_kind=>kim_real8_kind
   use dimensions, only: max_neigh_edges
!+isong
! use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
!-isong
   use timing, only: timingstart, timingstop
   use kiapsparallel, only: finpar, abortpar, omp_get_num_threads, omp_get_thread_num
   implicit none
!
   private
!
   type, public :: rotation_t
     sequence
     integer :: nbr ! nbr direction:north south east west
     integer :: reverse ! 0 = do not reverse order
! 1 = reverse order
     real(real_kind), dimension(:,:,:), pointer :: r ! rotation matrix
   end type rotation_t
!
   type, public :: edgedescriptor_t
     sequence
     integer(int_kind) :: use_rotation
     integer(int_kind) :: padding
     integer(int_kind) :: putmapp(max_neigh_edges)
     integer(int_kind) :: getmapp(max_neigh_edges)
     integer(int_kind) :: putmapp_ghost(max_neigh_edges)
     integer(int_kind) :: getmapp_ghost(max_neigh_edges)
     logical(log_kind) :: reverse(max_neigh_edges)
     type(rotation_t), dimension(:), pointer :: rot=>null() ! identifies list of edges
!  that must be rotated, and how
   end type edgedescriptor_t
!
   type, public :: edgebuffer_t
     sequence
     real(real_kind), dimension(:,:), pointer :: buf
     real(real_kind), dimension(:,:), pointer :: receive
     integer :: nlyr ! number of layers
     integer :: nbuf ! size of the horizontal dimension of the buffers.
   end type edgebuffer_t
!
   type, public :: longedgebuffer_t
     sequence
     integer :: nlyr
     integer :: nbuf
     integer(int_kind), dimension(:,:), pointer :: buf
     integer(int_kind), dimension(:,:), pointer :: receive
   end type longedgebuffer_t
   public :: initedgebuffer, initlongedgebuffer
   public :: freeedgebuffer, freelongedgebuffer
   public :: edgevpack, edgedgvpack, longedgevpack
   public :: edgevunpack, edgedgvunpack, edgevunpackvert
   public :: edgevunpackmin, longedgevunpackmin
   public :: edgevunpackmax, edgevunpackcheck
   public :: edgerotate
   public :: buffermap
   logical, private :: threadsafe = .true.
!
   type, public :: ghostbuffer_t
     sequence
     real(real_kind), dimension(:,:,:,:), pointer :: buf
     real(real_kind), dimension(:,:,:,:), pointer :: receive
     integer :: nlyr ! number of layers
     integer :: ntrac ! number of tracers
     integer :: nbuf ! size of the horizontal dimension of the buffers.
   end type ghostbuffer_t
!
   type, public :: ghostbuffertr_t
     sequence
     real(real_kind), dimension(:,:,:,:,:), pointer :: buf
     real(real_kind), dimension(:,:,:,:,:), pointer :: receive
     integer :: nlyr ! number of layers
     integer :: ntrac ! number of tracers
     integer :: nbuf ! size of the horizontal dimension of the buffers.
   end type ghostbuffertr_t
!
   type, public :: ghostbuffer3d_t
     sequence
     real(real_kind), dimension(:,:,:,:), pointer :: buf
     real(real_kind), dimension(:,:,:,:), pointer :: receive
     integer :: nlyr ! number of layers
     integer :: nhc ! number of layers of ghost cells
     integer :: np ! number of points in a cell
     integer :: nbuf ! size of the horizontal dimension of the buffers.
   end type ghostbuffer3d_t
   public :: initghostbufferfull
   public :: freeghostbuffer
   public :: freeghostbuffertr
   public :: freeghostbuffer3d
   public :: ghostvpackfull
   public :: ghostvunpackfull
   public :: ghostvpack
   public :: ghostvunpack
   public :: ghostvpack2d
   public :: ghostvunpack2d
!for VLAG
   public :: ghostvpack_new
   public :: ghostvunpack_new
   public :: ghostvpack2d_level
   public :: ghostvunpack2d_level
   public :: initghostbuffer
   public :: ghostvpack3d
   public :: ghostvunpack3d
   public :: initghostbuffer3d
   integer*8 :: edgebuff_addr(0:1)
! Wrap pointer so we can make an array of them.
!
   type :: wrap_ptr
     real(real_kind), dimension(:,:), pointer :: ptr=>null()
   end type wrap_ptr
   type(wrap_ptr) :: edgebuff_ptrs(0:1)
   integer :: l_ierr
!
   contains
! =========================================
! initEdgeBuffer:
!
! create an Real based communication buffer
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine initedgebuffer(edge, nlyr, buf_ptr, receive_ptr)
!-------------------------------------------------------------------------------
   use dimensions, only: np, nelemd, max_corner_elem
   implicit none
!
   integer, intent(in   ) :: nlyr
   type(edgebuffer_t), intent(  out), target :: edge
   real(real_kind), optional, pointer :: buf_ptr(:), receive_ptr(:)
! Notes about the buf_ptr/receive_ptr options:
!
! You can pass in 1D pointers to this function. If they are not
! associated, they will be allocated and used as buffer space. If they
! are associated, their targets will be used as buffer space.
!
! The pointers must not be thread-private.
!
! If an EdgeBuffer_t object is initialized from pre-existing storage
! (i.e. buf_ptr is provided and not null), it must *not* be freed,
! and must not be used if the underlying storage has been deallocated.
!
! All these restrictions also applied to the old newbuf and newreceive
! options.
! Workaround for NAG bug.
! NAG 5.3.1 dies if you use pointer bounds remapping to set
! a pointer that is also a component. So remap to temporary,
! then use that to set component pointer.
   real(real_kind), pointer :: tmp_ptr(:,:)
! Local variables
   integer :: nbuf, ith
!
   nbuf = 4*(np+max_corner_elem)*nelemd
   edge%nlyr = nlyr
   edge%nbuf = nbuf
   if (nlyr==0) return ! tracer code might call initedgebuffer() with zero tracers
!   only master thread should allocate the buffer
#if (! defined ELEMENT_OPENMP)
!$OMP MASTER
#endif
   if (present(buf_ptr)) then
! If buffer is passed in but not allocated, allocate it.
     if (.not.associated(buf_ptr)) allocate(buf_ptr(nlyr*nbuf))
! Verify dimensions
     if (size(buf_ptr)<nlyr*nbuf) then
       print*,'size(buf_ptr),nlyr,nbuf = ',size(buf_ptr),nlyr,nbuf
       call abortpar('Error:user provided edge buffer is too small')
     endif
#ifdef HAVE_F2003_PTR_BND_REMAP
     tmp_ptr(1:nlyr,1:nbuf)=>buf_ptr
     edge%buf=>tmp_ptr
#else
! call F77 routine which will reshape array.
     call remap_2d_ptr_buf(edge,nlyr,nbuf,buf_ptr)
#endif
   else
     allocate(edge%buf(nlyr,nbuf))
   endif
   if (present(receive_ptr)) then
! If buffer is passed in but not allocated, allocate it.
     if (.not.associated(receive_ptr)) allocate(receive_ptr(nlyr*nbuf))
! Verify dimensions
     if (size(receive_ptr)<nlyr*nbuf) then
       print*,'size(receive_ptr),nlyr,nbuf = ',size(receive_ptr),nlyr,nbuf
       call abortpar('Error:user provided edge buffer is too small')
     endif
#ifdef HAVE_F2003_PTR_BND_REMAP
     tmp_ptr(1:nlyr,1:nbuf)=>receive_ptr
     edge%receive=>tmp_ptr
#else
! call F77 routine which will reshape array.
     call remap_2d_ptr_receive(edge,nlyr,nbuf,receive_ptr)
#endif
   else
     allocate(edge%receive(nlyr,nbuf))
   endif
   edge%buf(:,:) = 0.0d0
   edge%receive(:,:) = 0.0d0
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
#endif
!   make sure all threads wait until buffer is allocated
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
! sanity check on edge.  edge must NOT be thread-prive,but must be shared by all threads
! the calling program needs to instantiate 'edge' outside the threaded region.
! if 'edge' is thread private,it creates flaky openMP problems that are difficut to debug
! so lets try and detect it here:
   if (omp_get_num_threads()>1) then
     ith = omp_get_thread_num()
     if (ith<=1) then
       edgebuff_ptrs(ith)%ptr=>edge%buf
     endif
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
     if (.not.associated(edgebuff_ptrs(0)%ptr,edgebuff_ptrs(1)%ptr)) then
       call abortpar('ERROR:edge struct appears to be thread-private.')
     endif
#if (! defined ELEMENT_OPENMP)
!$OMP ENDMASTER
#endif
   endif
!
   end subroutine initedgebuffer
! =========================================
! initEdgeBuffer:
!
! create an Real based communication buffer
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine initedgebuffer_old(edge, nlyr)
!-------------------------------------------------------------------------------
   use dimensions, only: np, nelemd, max_corner_elem
   implicit none
!
   integer, intent(in   ) :: nlyr
   type(edgebuffer_t), intent(  out) :: edge
! Local variables
   integer :: nbuf, ith
! sanity check on edge.  edge must NOT be thread-prive, but must be shared by all threads
! the calling program needs to instantiate 'edge' outside the threaded region.
! if 'edge' is thread private, it creates flaky openMP problems that are difficut to debug
! so lets try and detect it here:
!
   if (omp_get_num_threads()>1) then
     ith = omp_get_thread_num()
     if (ith<=1) then
       edgebuff_addr(ith) = loc(edge)
     endif
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
     if (edgebuff_addr(0).ne.edgebuff_addr(1)) then
       call finpar(message = 'ERROR:edge struct appears to be thread-private.')
     endif
#if (! defined ELEMENT_OPENMP)
!$OMP ENDMASTER
#endif
   endif
!
   nbuf = 4*(np+max_corner_elem)*nelemd
!   only master thread should allocate the buffer
#if (! defined ELEMENT_OPENMP)
!$OMP MASTER
#endif
   edge%nlyr = nlyr
   edge%nbuf = nbuf
   allocate(edge%buf(nlyr,nbuf))
   allocate(edge%receive(nlyr,nbuf))
   edge%buf(:,:) = 0.0d0
   edge%receive(:,:) = 0.0d0
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
#endif
!   make sure all threads wait until buffer is allocated
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
!
   end subroutine initedgebuffer_old
! =========================================
! initLongEdgeBuffer:
!
! create an Integer based communication buffer
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine initlongedgebuffer(edge, nlyr)
!-------------------------------------------------------------------------------
   use dimensions, only: np, nelemd, max_corner_elem
   implicit none
!
   integer, intent(in   ) :: nlyr
   type(longedgebuffer_t), intent(  out) :: edge
! Local variables
   integer :: nbuf
! sanity check for threading
!
   if (omp_get_num_threads()>1) then
     call finpar(message = 'ERROR:initlongedgebuffer must be called before threaded reagion')
   endif
!
   nbuf = 4*(np+max_corner_elem)*nelemd
   edge%nlyr = nlyr
   edge%nbuf = nbuf
   allocate(edge%buf(nlyr,nbuf))
   edge%buf(:,:) = 0
   allocate(edge%receive(nlyr,nbuf))
   edge%receive(:,:) = 0
!
   end subroutine initlongedgebuffer
! =========================================
! edgeDGVpack:
!
! Pack edges of v into buf for DG stencil
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine edgedgvpack(edge, v, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: np
!
   type(edgebuffer_t) :: edge
   integer, intent(in   ) :: vlyr
   real(real_kind), intent(in   ) :: v(np, np, vlyr)
   integer, intent(in   ) :: kptr
   type(edgedescriptor_t) :: desc
! =========================================
! This code is just a wrapper call the
!   normal edgeVpack
! =========================================
!
   call edgevpack(edge,v,vlyr,kptr,desc)
!
   end subroutine edgedgvpack
! ===========================================
!  FreeEdgeBuffer:
!
!  Freed an edge communication buffer
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine freeedgebuffer(edge)
!-------------------------------------------------------------------------------
   implicit none
!
   type(edgebuffer_t), intent(inout) :: edge
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
!
   edge%nbuf = 0
   edge%nlyr = 0
   deallocate(edge%buf)
   deallocate(edge%receive)
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
#endif
!
   end subroutine freeedgebuffer
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine freeghostbuffer3d(buffer)
!-------------------------------------------------------------------------------
   implicit none
!
   type(ghostbuffer3d_t), intent(inout) :: buffer
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
!
   buffer%nbuf = 0
   buffer%nlyr = 0
   deallocate(buffer%buf)
   deallocate(buffer%receive)
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
#endif
!
   end subroutine freeghostbuffer3d
! ===========================================
!  FreeLongEdgeBuffer:
!
!  Freed an edge communication buffer
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine freelongedgebuffer(edge)
!-------------------------------------------------------------------------------
   implicit none
!
   type(longedgebuffer_t), intent(inout) :: edge
!
   edge%nbuf = 0
   edge%nlyr = 0
   deallocate(edge%buf)
   deallocate(edge%receive)
!
   end subroutine freelongedgebuffer
! =========================================
!
!> @brief Pack edges of v into an edge buffer for boundary exchange.
!
!> This subroutine packs for one or more vertical layers into an edge
!! buffer. If the buffer associated with edge is not large enough to
!! hold all vertical layers you intent to pack, the method will
!! halt the program with a call to parallel_mod::FinPar(message=).
!! @param[in] edge Edge Buffer into which the data will be packed.
!! This buffer must be previously allocated with initEdgeBuffer().
!! @param[in] v The data to be packed.
!! @param[in] vlyr Number of vertical level coming into the subroutine
!! for packing for input v.
!! @param[in] kptr Vertical pointer to the place in the edge buffer where
!! data will be located.
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine edgevpack(edge, v, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: np, max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(edgebuffer_t) :: edge
   integer, intent(in   ) :: vlyr
   real(real_kind), intent(in   ) :: v(np, np, vlyr)
   integer, intent(in   ) :: kptr
   type(edgedescriptor_t), intent(in   ) :: desc
! Local variables
   logical, parameter :: useunroll = .true.
   integer :: i, k, ir, l
   integer :: is, ie, in, iw
!$OMP MASTER
!+isong
!   call t_startf('edge_pack')
!-isong
!    l_ierr = TimingStart('edge_pack')
!$OMP END MASTER
!
   if (.not.threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
     threadsafe = .true.
   endif
!
   is = desc%putmapp(south)
   ie = desc%putmapp(east)
   in = desc%putmapp(north)
   iw = desc%putmapp(west)
   if (edge%nlyr<(kptr+vlyr)) then
     call finpar(message = 'edgeVpack:buffer overflow:size of the vertical dimension must be increased!')
   endif
   if (modulo(np,2)==0.and.useunroll) then
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i)
#endif
     do k = 1, vlyr
       do i = 1, np, 2
         edge%buf(kptr+k,is+i) = v(i,1,k)
         edge%buf(kptr+k,is+i+1) = v(i+1,1,k)
         edge%buf(kptr+k,ie+i) = v(np,i,k)
         edge%buf(kptr+k,ie+i+1) = v(np,i+1,k)
         edge%buf(kptr+k,in+i) = v(i,np,k)
         edge%buf(kptr+k,in+i+1) = v(i+1,np,k)
         edge%buf(kptr+k,iw+i) = v(1,i,k)
         edge%buf(kptr+k,iw+i+1) = v(1,i+1,k)
       enddo
       enddo
       else
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i)
#endif
         do k = 1, vlyr
           do i = 1, np
             edge%buf(kptr+k,is+i) = v(i,1,k)
             edge%buf(kptr+k,ie+i) = v(np,i,k)
             edge%buf(kptr+k,in+i) = v(i,np,k)
             edge%buf(kptr+k,iw+i) = v(1,i,k)
           enddo
           enddo
   endif
!  This is really kludgy way to setup the index reversals
!  But since it is so a rare event not real need to spend time optimizing
   if (desc%reverse(south)) then
     is = desc%putmapp(south)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,ir)
#endif
     do k = 1,vlyr
       do i = 1, np
         ir = np-i+1
         edge%buf(kptr+k,is+ir) = v(i,1,k)
       enddo
     enddo
   endif
   if (desc%reverse(east)) then
     ie = desc%putmapp(east)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,ir)
#endif
     do k = 1,vlyr
       do i = 1, np
         ir = np-i+1
         edge%buf(kptr+k,ie+ir) = v(np,i,k)
       enddo
     enddo
   endif
   if (desc%reverse(north)) then
   in = desc%putmapp(north)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,ir)
#endif
     do k = 1, vlyr
       do i = 1, np
         ir = np-i+1
         edge%buf(kptr+k,in+ir) = v(i,np,k)
       enddo
       enddo
   endif
   if (desc%reverse(west)) then
     iw = desc%putmapp(west)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,ir)
#endif
     do k = 1,vlyr
       do i = 1, np
         ir = np-i+1
         edge%buf(kptr+k,iw+ir) = v(1,i,k)
       enddo
     enddo
   endif
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (desc%putmapp(l)/=-1) then
       do k = 1, vlyr
         edge%buf(kptr+k,desc%putmapp(l)+1) = v(1,1,k)
       enddo
       endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (desc%putmapp(l)/=-1) then
       do k = 1, vlyr
         edge%buf(kptr+k,desc%putmapp(l)+1) = v(np,1,k)
       enddo
       endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (desc%putmapp(l)/=-1) then
       do k = 1, vlyr
         edge%buf(kptr+k,desc%putmapp(l)+1) = v(np,np,k)
       enddo
       endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (desc%putmapp(l)/=-1) then
       do k = 1, vlyr
         edge%buf(kptr+k,desc%putmapp(l)+1) = v(1,np,k)
       enddo
       endif
   enddo
!$OMP MASTER
!+isong
!   call t_stopf('edge_pack')
!-isong
!    l_ierr = TimingStop('edge_pack')
!$OMP END MASTER
!
   end subroutine edgevpack
! =========================================
! LongEdgeVpack:
!
! Pack edges of v into buf...
! =========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine longedgevpack(edge, v, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use control, only: north, south, east, west, neast, nwest, seast, swest
   use dimensions, only: np, max_corner_elem
!
   type(longedgebuffer_t) :: edge
   integer, intent(in   ) :: vlyr
   integer(int_kind), intent(in   ) :: v(np, np, vlyr)
   integer, intent(in   ) :: kptr
   type(edgedescriptor_t), intent(in   ) :: desc
! Local variables
   logical, parameter :: useunroll = .true.
   integer :: i, k, ir, l
   integer :: is, ie, in, iw
!
   if (.not.threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
     threadsafe = .true.
   endif
!
   is = desc%putmapp(south)
   ie = desc%putmapp(east)
   in = desc%putmapp(north)
   iw = desc%putmapp(west)
   if (modulo(np,2)==0.and.useunroll) then
     do k = 1, vlyr
       do i = 1, np, 2
         edge%buf(kptr+k,is+i) = v(i,1,k)
         edge%buf(kptr+k,is+i+1) = v(i+1,1,k)
         edge%buf(kptr+k,ie+i) = v(np,i,k)
         edge%buf(kptr+k,ie+i+1) = v(np,i+1,k)
         edge%buf(kptr+k,in+i) = v(i,np,k)
         edge%buf(kptr+k,in+i+1) = v(i+1,np,k)
         edge%buf(kptr+k,iw+i) = v(1,i,k)
         edge%buf(kptr+k,iw+i+1) = v(1,i+1,k)
       enddo
       enddo
       else
         do k = 1, vlyr
           do i = 1, np
             edge%buf(kptr+k,is+i) = v(i,1,k)
             edge%buf(kptr+k,ie+i) = v(np,i,k)
             edge%buf(kptr+k,in+i) = v(i,np,k)
             edge%buf(kptr+k,iw+i) = v(1,i,k)
           enddo
           enddo
   endif
!  This is really kludgy way to setup the index reversals
!  But since it is so a rare event not real need to spend time optimizing
   if (desc%reverse(south)) then
     is = desc%putmapp(south)
     do k = 1,vlyr
       do i = 1, np
         ir = np-i+1
         edge%buf(kptr+k,is+ir) = v(i,1,k)
       enddo
     enddo
   endif
   if (desc%reverse(east)) then
     ie = desc%putmapp(east)
     do k = 1,vlyr
       do i = 1, np
         ir = np-i+1
         edge%buf(kptr+k,ie+ir) = v(np,i,k)
       enddo
     enddo
   endif
   if (desc%reverse(north)) then
   in = desc%putmapp(north)
     do k = 1, vlyr
       do i = 1, np
         ir = np-i+1
         edge%buf(kptr+k,in+ir) = v(i,np,k)
       enddo
       enddo
   endif
   if (desc%reverse(west)) then
     iw = desc%putmapp(west)
     do k = 1,vlyr
       do i = 1, np
         ir = np-i+1
         edge%buf(kptr+k,iw+ir) = v(1,i,k)
       enddo
     enddo
   endif
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (desc%putmapp(l)/=-1) then
       do k = 1, vlyr
         edge%buf(kptr+k,desc%putmapp(l)+1) = v(1,1,k)
       enddo
       endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (desc%putmapp(l)/=-1) then
       do k = 1, vlyr
         edge%buf(kptr+k,desc%putmapp(l)+1) = v(np,1,k)
       enddo
       endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (desc%putmapp(l)/=-1) then
       do k = 1, vlyr
         edge%buf(kptr+k,desc%putmapp(l)+1) = v(np,np,k)
       enddo
       endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (desc%putmapp(l)/=-1) then
       do k = 1, vlyr
         edge%buf(kptr+k,desc%putmapp(l)+1) = v(1,np,k)
       enddo
       endif
   enddo
!
   end subroutine longedgevpack
! ========================================
! edgeVunpack:
!
! Unpack edges from edge buffer into v...
! ========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine edgevunpack(edge, v, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: np, max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(edgebuffer_t), intent(in   ) :: edge
   integer, intent(in   ) :: vlyr
   real(real_kind), intent(inout) :: v(np, np, vlyr)
   integer, intent(in   ) :: kptr
   type(edgedescriptor_t) :: desc
! Local
!logical, parameter :: UseUnroll = .TRUE.
   logical, parameter :: useunroll = .false. !'.true.':not apply to reproduced sum.
   real(real_kind) :: err_11(vlyr), err_41(vlyr), err_44(vlyr), err_14(vlyr)
   integer :: i, k, l
   integer :: is, ie, in, iw
!$OMP MASTER
!call t_startf('edge_unpack')
!$OMP END MASTER
!
   threadsafe = .false.
   is = desc%getmapp(south)
   ie = desc%getmapp(east)
   in = desc%getmapp(north)
   iw = desc%getmapp(west)
   err_11 = 0.0_real_kind
   err_41 = 0.0_real_kind
   err_44 = 0.0_real_kind
   err_14 = 0.0_real_kind
   if (modulo(np,2)==0.and.useunroll) then
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i)
#endif
     do k = 1, vlyr
       do i = 1, np, 2
         v(i,1,k) = v(i,1,k)+edge%buf(kptr+k,is+i)
         v(i+1,1,k) = v(i+1,1,k)+edge%buf(kptr+k,is+i+1)
         v(np,i,k) = v(np,i,k)+edge%buf(kptr+k,ie+i)
         v(np,i+1,k) = v(np,i+1,k)+edge%buf(kptr+k,ie+i+1)
         v(i,np,k) = v(i,np,k)+edge%buf(kptr+k,in+i)
         v(i+1,np,k) = v(i+1,np,k)+edge%buf(kptr+k,in+i+1)
         v(1,i,k) = v(1,i,k)+edge%buf(kptr+k,iw+i)
         v(1,i+1,k) = v(1,i+1,k)+edge%buf(kptr+k,iw+i+1)
       enddo
       enddo
       else
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i)
#endif
         do k = 1, vlyr
           do i = 2, np-1
             v(i,1,k) = v(i,1,k)+edge%buf(kptr+k,is+i)
             v(np,i,k) = v(np,i,k)+edge%buf(kptr+k,ie+i)
             v(i,np,k) = v(i,np,k)+edge%buf(kptr+k,in+i)
             v(1,i,k) = v(1,i,k)+edge%buf(kptr+k,iw+i)
           enddo
!v(1  ,1  ,k) = v(1  ,1  ,k)+edge%buf(kptr+k,is+1  )
!v(1  ,1  ,k) = v(1  ,1  ,k)+edge%buf(kptr+k,iw+1  )
             call ddpddlocal(v(1,1,k),err_11(k),edge%buf(kptr+k,is+1))
             call ddpddlocal(v(1,1,k),err_11(k),edge%buf(kptr+k,iw+1))
!v(np ,1  ,k) = v(np ,1  ,k)+edge%buf(kptr+k,ie+1  )
!v(np ,1  ,k) = v(np ,1  ,k)+edge%buf(kptr+k,is+np )
             call ddpddlocal(v(np,1,k),err_41(k),edge%buf(kptr+k,ie+1))
             call ddpddlocal(v(np,1,k),err_41(k),edge%buf(kptr+k,is+np))
!v(1  ,np ,k) = v(1  ,np ,k)+edge%buf(kptr+k,in+1  )
!v(1  ,np ,k) = v(1  ,np ,k)+edge%buf(kptr+k,iw+np )
             call ddpddlocal(v(1,np,k),err_14(k),edge%buf(kptr+k,in+1))
             call ddpddlocal(v(1,np,k),err_14(k),edge%buf(kptr+k,iw+np))
!v(np ,np ,k) = v(np ,np ,k)+edge%buf(kptr+k,ie+np )
!v(np ,np ,k) = v(np ,np ,k)+edge%buf(kptr+k,in+np )
             call ddpddlocal(v(np,np,k),err_44(k),edge%buf(kptr+k,ie+np))
             call ddpddlocal(v(np,np,k),err_44(k),edge%buf(kptr+k,in+np))
           enddo
   endif
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
!v(1  ,1 ,k)=v(1 ,1 ,k)+edge%buf(kptr+k,desc%getmapP(l)+1)
         call ddpddlocal(v(1,1,k),err_11(k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
!v(np ,1 ,k)=v(np,1 ,k)+edge%buf(kptr+k,desc%getmapP(l)+1)
         call ddpddlocal(v(np,1,k),err_41(k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
!v(np ,np,k)=v(np,np,k)+edge%buf(kptr+k,desc%getmapP(l)+1)
         call ddpddlocal(v(np,np,k),err_44(k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
!v(1  ,np,k)=v(1 ,np,k)+edge%buf(kptr+k,desc%getmapP(l)+1)
         call ddpddlocal(v(1,np,k),err_14(k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
!$OMP MASTER
!call t_stopf('edge_unpack')
!$OMP END MASTER
!
   end subroutine edgevunpack
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine edgevunpack_old(edge, v, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: np, max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(edgebuffer_t), intent(in   ) :: edge
   integer, intent(in   ) :: vlyr
   real(real_kind), intent(inout) :: v(np, np, vlyr)
   integer, intent(in   ) :: kptr
   type(edgedescriptor_t) :: desc
! Local
   logical, parameter :: useunroll = .true.
   integer :: i, k, l
   integer :: is, ie, in, iw
   integer :: l_ierr
!$OMP MASTER
!+isong
!   call t_startf('edge_unpack')
!-isong
!    l_ierr = TimingStart('edge_unpack')
!$OMP END MASTER
!
   threadsafe = .false.
   is = desc%getmapp(south)
   ie = desc%getmapp(east)
   in = desc%getmapp(north)
   iw = desc%getmapp(west)
   if (modulo(np,2)==0.and.useunroll) then
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i)
#endif
     do k = 1, vlyr
       do i = 1, np, 2
         v(i,1,k) = v(i,1,k)+edge%buf(kptr+k,is+i)
         v(i+1,1,k) = v(i+1,1,k)+edge%buf(kptr+k,is+i+1)
         v(np,i,k) = v(np,i,k)+edge%buf(kptr+k,ie+i)
         v(np,i+1,k) = v(np,i+1,k)+edge%buf(kptr+k,ie+i+1)
         v(i,np,k) = v(i,np,k)+edge%buf(kptr+k,in+i)
         v(i+1,np,k) = v(i+1,np,k)+edge%buf(kptr+k,in+i+1)
         v(1,i,k) = v(1,i,k)+edge%buf(kptr+k,iw+i)
         v(1,i+1,k) = v(1,i+1,k)+edge%buf(kptr+k,iw+i+1)
       enddo
       enddo
       else
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i)
#endif
         do k = 1, vlyr
           do i = 1, np
             v(i,1,k) = v(i,1,k)+edge%buf(kptr+k,is+i)
             v(np,i,k) = v(np,i,k)+edge%buf(kptr+k,ie+i)
             v(i,np,k) = v(i,np,k)+edge%buf(kptr+k,in+i)
             v(1,i,k) = v(1,i,k)+edge%buf(kptr+k,iw+i)
           enddo
           enddo
   endif
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(1,1,k) = v(1,1,k)+edge%buf(kptr+k,desc%getmapp(l)+1)
       enddo
       endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(np,1,k) = v(np,1,k)+edge%buf(kptr+k,desc%getmapp(l)+1)
       enddo
       endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(np,np,k) = v(np,np,k)+edge%buf(kptr+k,desc%getmapp(l)+1)
       enddo
       endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(1,np,k) = v(1,np,k)+edge%buf(kptr+k,desc%getmapp(l)+1)
       enddo
       endif
   enddo
!$OMP MASTER
!+isong
!   call t_stopf('edge_unpack')
!-isong
!    l_ierr = TimingStop('edge_unpack')
!$OMP END MASTER
!
   end subroutine edgevunpack_old
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine edgevunpackvert(edge, v, desc)
!-------------------------------------------------------------------------------
   use control, only: north, south, east, west, neast, nwest, seast, swest
   use dimensions, only: np, max_corner_elem, ne
   use coordinatesystems, only: cartesian3d_t
!
   type(edgebuffer_t), intent(inout) :: edge
   type(cartesian3d_t), intent(  out) :: v(:,:,:)
   type(edgedescriptor_t) :: desc
! Local
   logical, parameter :: useunroll = .true.
   integer :: i, k, l
   integer :: is, ie, in, iw
!
   threadsafe = .false.
   if (max_corner_elem.ne.1.and.ne==0) then
! MNL: this is used to construct the dual grid on the cube,
!      currently only supported for the uniform grid. If
!      this is desired on a refined grid, a little bit of
!      work will be required.
     call finpar(message = "edgevunpackvert should not be called with unstructured meshes")
   endif
   is = desc%getmapp(south)
   ie = desc%getmapp(east)
   in = desc%getmapp(north)
   iw = desc%getmapp(west)
! N+S
   do i = 1,np/2
! North
     v(3,i,np)%x = edge%buf(1,in+i)
     v(3,i,np)%y = edge%buf(2,in+i)
     v(3,i,np)%z = edge%buf(3,in+i)
! South
     v(2,i,1)%x = edge%buf(1,is+i)
     v(2,i,1)%y = edge%buf(2,is+i)
     v(2,i,1)%z = edge%buf(3,is+i)
   enddo
   do i = np/2+1,np
! North
     v(4,i,np)%x = edge%buf(1,in+i)
     v(4,i,np)%y = edge%buf(2,in+i)
     v(4,i,np)%z = edge%buf(3,in+i)
! South
     v(1,i,1)%x = edge%buf(1,is+i)
     v(1,i,1)%y = edge%buf(2,is+i)
     v(1,i,1)%z = edge%buf(3,is+i)
   enddo
   do i = 1,np/2
! East
     v(3,np,i)%x = edge%buf(1,ie+i)
     v(3,np,i)%y = edge%buf(2,ie+i)
     v(3,np,i)%z = edge%buf(3,ie+i)
! West
     v(4,1,i)%x = edge%buf(1,iw+i)
     v(4,1,i)%y = edge%buf(2,iw+i)
     v(4,1,i)%z = edge%buf(3,iw+i)
   enddo
   do i = np/2+1,np
! East
     v(2,np,i)%x = edge%buf(1,ie+i)
     v(2,np,i)%y = edge%buf(2,ie+i)
     v(2,np,i)%z = edge%buf(3,ie+i)
! West
     v(1,1,i)%x = edge%buf(1,iw+i)
     v(1,1,i)%y = edge%buf(2,iw+i)
     v(1,1,i)%z = edge%buf(3,iw+i)
   enddo
! SWEST
   do l = swest,swest+max_corner_elem-1
! find the one active corner, then exist
     if (desc%getmapp(l)/=-1) then
       v(1,1,1)%x = edge%buf(1,desc%getmapp(l)+1)
       v(1,1,1)%y = edge%buf(2,desc%getmapp(l)+1)
       v(1,1,1)%z = edge%buf(3,desc%getmapp(l)+1)
       exit
     else
       v(1,1,1)%x = 0_real_kind
       v(1,1,1)%y = 0_real_kind
       v(1,1,1)%z = 0_real_kind
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
! find the one active corner, then exist
     if (desc%getmapp(l)/=-1) then
       v(2,np,1)%x = edge%buf(1,desc%getmapp(l)+1)
       v(2,np,1)%y = edge%buf(2,desc%getmapp(l)+1)
       v(2,np,1)%z = edge%buf(3,desc%getmapp(l)+1)
       exit
     else
       v(2,np,1)%x = 0_real_kind
       v(2,np,1)%y = 0_real_kind
       v(2,np,1)%z = 0_real_kind
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
! find the one active corner, then exist
     if (desc%getmapp(l)/=-1) then
       v(3,np,np)%x = edge%buf(1,desc%getmapp(l)+1)
       v(3,np,np)%y = edge%buf(2,desc%getmapp(l)+1)
       v(3,np,np)%z = edge%buf(3,desc%getmapp(l)+1)
       exit
     else
       v(3,np,np)%x = 0_real_kind
       v(3,np,np)%y = 0_real_kind
       v(3,np,np)%z = 0_real_kind
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
! find the one active corner, then exist
     if (desc%getmapp(l)/=-1) then
       v(4,1,np)%x = edge%buf(1,desc%getmapp(l)+1)
       v(4,1,np)%y = edge%buf(2,desc%getmapp(l)+1)
       v(4,1,np)%z = edge%buf(3,desc%getmapp(l)+1)
       exit
     else
       v(4,1,np)%x = 0_real_kind
       v(4,1,np)%y = 0_real_kind
       v(4,1,np)%z = 0_real_kind
     endif
   enddo
! Fill the missing vertex info
   do i = 2,np/2
! North
     v(4,i,np)%x = v(3,i-1,np)%x
     v(4,i,np)%y = v(3,i-1,np)%y
     v(4,i,np)%z = v(3,i-1,np)%z
! South
     v(1,i,1)%x = v(2,i-1,1)%x
     v(1,i,1)%y = v(2,i-1,1)%y
     v(1,i,1)%z = v(2,i-1,1)%z
   enddo
   do i = np/2+1,np-1
! North
     v(3,i,np)%x = v(4,i+1,np)%x
     v(3,i,np)%y = v(4,i+1,np)%y
     v(3,i,np)%z = v(4,i+1,np)%z
! South
     v(2,i,1)%x = v(1,i+1,1)%x
     v(2,i,1)%y = v(1,i+1,1)%y
     v(2,i,1)%z = v(1,i+1,1)%z
   enddo
   do i = 2,np/2
! East
     v(2,np,i)%x = v(3,np,i-1)%x
     v(2,np,i)%y = v(3,np,i-1)%y
     v(2,np,i)%z = v(3,np,i-1)%z
! West
     v(1,1,i)%x = v(4,1,i-1)%x
     v(1,1,i)%y = v(4,1,i-1)%y
     v(1,1,i)%z = v(4,1,i-1)%z
   enddo
   do i = np/2+1,np-1
! East
     v(3,np,i)%x = v(2,np,i+1)%x
     v(3,np,i)%y = v(2,np,i+1)%y
     v(3,np,i)%z = v(2,np,i+1)%z
! West
     v(4,1,i)%x = v(1,1,i+1)%x
     v(4,1,i)%y = v(1,1,i+1)%y
     v(4,1,i)%z = v(1,1,i+1)%z
   enddo
!
   end subroutine edgevunpackvert
! ========================================
! edgeDGVunpack:
!
! Unpack edges from edge buffer into v...
! ========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine edgedgvunpack(edge, v, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: np
   use control, only: north, south, east, west
!
   type(edgebuffer_t), intent(in   ) :: edge
   integer, intent(in   ) :: vlyr
   real(real_kind), intent(  out) :: v(0:np+1, 0:np+1, vlyr)
   integer, intent(in   ) :: kptr
   type(edgedescriptor_t) :: desc
! Local
   integer :: i, k
   integer :: is, ie, in, iw
!
   threadsafe = .false.
   is = desc%getmapp(south)
   ie = desc%getmapp(east)
   in = desc%getmapp(north)
   iw = desc%getmapp(west)
   do k = 1,vlyr
     do i = 1, np
       v(i,0,k) = edge%buf(kptr+k,is+i)
       v(np+1,i,k) = edge%buf(kptr+k,ie+i)
       v(i,np+1,k) = edge%buf(kptr+k,in+i)
       v(0,i,k) = edge%buf(kptr+k,iw+i)
     enddo
   enddo
!
   end subroutine edgedgvunpack
! ========================================
! edgeVunpackMIN/MAX:
!
! Finds the Min/Max edges from edge buffer into v...
! ========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine edgevunpackmax(edge, v, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: np, max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(edgebuffer_t), intent(in   ) :: edge
   integer, intent(in   ) :: vlyr
   real(real_kind), intent(inout) :: v(np, np, vlyr)
   integer, intent(in   ) :: kptr
   type(edgedescriptor_t), intent(in   ) :: desc
! Local
   integer :: i, k, l
   integer :: is, ie, in, iw
!
   threadsafe = .false.
   is = desc%getmapp(south)
   ie = desc%getmapp(east)
   in = desc%getmapp(north)
   iw = desc%getmapp(west)
   do k = 1,vlyr
     do i = 1, np
       v(i,1,k) = max(v(i,1,k),edge%buf(kptr+k,is+i))
       v(np,i,k) = max(v(np,i,k),edge%buf(kptr+k,ie+i))
       v(i,np,k) = max(v(i,np,k),edge%buf(kptr+k,in+i))
       v(1,i,k) = max(v(1,i,k),edge%buf(kptr+k,iw+i))
     enddo
   enddo
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(1,1,k) = max(v(1,1,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(np,1,k) = max(v(np,1,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(np,np,k) = max(v(np,np,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(1,np,k) = max(v(1,np,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
!
   end subroutine edgevunpackmax
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine edgevunpackmin(edge, v, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: np, max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(edgebuffer_t), intent(in   ) :: edge
   integer, intent(in   ) :: vlyr
   real(real_kind), intent(inout) :: v(np, np, vlyr)
   integer, intent(in   ) :: kptr
   type(edgedescriptor_t), intent(in   ) :: desc
! Local
   integer :: i, k, l
   integer :: is, ie, in, iw
!
   threadsafe = .false.
   is = desc%getmapp(south)
   ie = desc%getmapp(east)
   in = desc%getmapp(north)
   iw = desc%getmapp(west)
   do k = 1,vlyr
     do i = 1, np
       v(i,1,k) = min(v(i,1,k),edge%buf(kptr+k,is+i))
       v(np,i,k) = min(v(np,i,k),edge%buf(kptr+k,ie+i))
       v(i,np,k) = min(v(i,np,k),edge%buf(kptr+k,in+i))
       v(1,i,k) = min(v(1,i,k),edge%buf(kptr+k,iw+i))
     enddo
   enddo
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(1,1,k) = min(v(1,1,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(np,1,k) = min(v(np,1,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(np,np,k) = min(v(np,np,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(1,np,k) = min(v(1,np,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
!
   end subroutine edgevunpackmin
! ========================================
! LongEdgeVunpackMIN:
!
! Finds the Min edges from edge buffer into v...
! ========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine longedgevunpackmin(edge, v, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use control, only: north, south, east, west, neast, nwest, seast, swest
   use dimensions, only: np, max_corner_elem
!
   type(longedgebuffer_t), intent(in   ) :: edge
   integer, intent(in   ) :: vlyr
   integer(int_kind), intent(inout) :: v(np, np, vlyr)
   integer, intent(in   ) :: kptr
   type(edgedescriptor_t), intent(in   ) :: desc
! Local
   integer :: i, k, l
   integer :: is, ie, in, iw
!
   threadsafe = .false.
   is = desc%getmapp(south)
   ie = desc%getmapp(east)
   in = desc%getmapp(north)
   iw = desc%getmapp(west)
   do k = 1,vlyr
     do i = 1, np
       v(i,1,k) = min(v(i,1,k),edge%buf(kptr+k,is+i))
       v(np,i,k) = min(v(np,i,k),edge%buf(kptr+k,ie+i))
       v(i,np,k) = min(v(i,np,k),edge%buf(kptr+k,in+i))
       v(1,i,k) = min(v(1,i,k),edge%buf(kptr+k,iw+i))
     enddo
   enddo
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(1,1,k) = min(v(1,1,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(np,1,k) = min(v(np,1,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(np,np,k) = min(v(np,np,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
         v(1,np,k) = min(v(1,np,k),edge%buf(kptr+k,desc%getmapp(l)+1))
       enddo
       endif
   enddo
!
   end subroutine longedgevunpackmin
! =============================
! edgerotate:
!
! Rotate edges in buffer...
! =============================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine edgerotate(edge, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: np
!
   type(edgebuffer_t) :: edge ! edge struct
   integer, intent(in   ) :: vlyr ! number of 2d vector fields to rotate
   integer, intent(in   ) :: kptr ! layer pointer into edge buffer
   type(edgedescriptor_t) :: desc
! Local variables
   integer :: i, k, k1, k2
   integer :: irot, ia, nbr
   real(real_kind), dimension(:,:,:), pointer :: r
   real(real_kind) :: tmp1, tmp2
#ifdef _USEASSOCIATED
!
   if (associated(rot)) then
#else
     if (desc%use_rotation==1) then
#endif
       do irot = 1, size(desc%rot)
         nbr = desc%rot(irot)%nbr
         r=>desc%rot(irot)%r
         ia = desc%putmapp(nbr)
! ========================================
! If nbr direction is (1-4) => is an edge
! ========================================
         if (nbr<=4) then
! ========================================================
!  Is an edge. Rotate it in place
! ========================================================
           do i = 1, np
             do k = 1, vlyr, 2
               k1 = kptr+k
               k2 = k1+1
               tmp1 = r(1,1,i)*edge%buf(k1,ia+i)+r(1,2,i)*edge%buf(k2,ia+i)
               tmp2 = r(2,1,i)*edge%buf(k1,ia+i)+r(2,2,i)*edge%buf(k2,ia+i)
               edge%buf(k1,ia+i) = tmp1
               edge%buf(k2,ia+i) = tmp2
             enddo
             enddo
             else
! ===================================================
! Is an element corner point, but not a cube corner
! point, just rotate it in place.
! ===================================================
               if (ia/=-1) then
                 do k = 1, vlyr, 2
                   k1 = kptr+k
                   k2 = k1+1
                   tmp1 = r(1,1,1)*edge%buf(k1,ia+1)+r(1,2,1)*edge%buf(k2,ia+1)
                   tmp2 = r(2,1,1)*edge%buf(k1,ia+1)+r(2,2,1)*edge%buf(k2,ia+1)
                   edge%buf(k1,ia+1) = tmp1
                   edge%buf(k2,ia+1) = tmp2
                 enddo
                 endif
         endif
       enddo
       endif
   end subroutine edgerotate
! =============================================
! buffermap:
!
! buffermap translates element number, inum and
! element edge/corner, facet, into an edge buffer
! memory location, loc.
! =============================================
   function buffermap(inum, facet) result(loc)
!-------------------------------------------------------------------------------
   use dimensions, only: np
!
   integer, intent(in   ) :: inum
   integer, intent(in   ) :: facet
   integer :: loc
!
   if (facet>4) then
     if (inum==-1) then
       loc = inum
     else
       loc =(inum-1)*(4*np+4)+4*np+(facet-5)
     endif
     else
       loc =(inum-1)*(4*np+4)+np*(facet-1)
   endif
   end function buffermap
! =========================================
! initghostbufferfull:
!
! create an Real based communication buffer
! =========================================
   subroutine initghostbufferfull(edge, nlyr, nc)
!-------------------------------------------------------------------------------
   use dimensions, only: nelemd, max_neigh_edges
   implicit none
!
   integer, intent(in   ) :: nlyr, nc
   type(ghostbuffer_t), intent(  out) :: edge
! Local variables
   integer :: nbuf
! sanity check for threading
!
   if (omp_get_num_threads()>1) then
     call finpar(message = 'ERROR:initghostbufferfull must be called before threaded region')
   endif
!
   nbuf = max_neigh_edges*nelemd
#if (! defined ELEMENT_OPENMP)
!$OMP MASTER
#endif
   edge%nlyr = nlyr
   edge%nbuf = nbuf
   allocate(edge%buf(nc,nc,nlyr,nbuf))
   edge%buf = 0
   allocate(edge%receive(nc,nc,nlyr,nbuf))
   edge%receive = 0
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif
   end subroutine initghostbufferfull
! ===========================================
!  FreeGhostBuffer:
!  Author: Christoph Erath, Mark Taylor
!  Freed an ghostpoints communication buffer
! =========================================
   subroutine freeghostbuffer(ghost)
!-------------------------------------------------------------------------------
   implicit none
!
   type(ghostbuffer_t), intent(inout) :: ghost
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
!
   ghost%nbuf = 0
   ghost%nlyr = 0
   deallocate(ghost%buf)
   deallocate(ghost%receive)
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
#endif
   end subroutine freeghostbuffer
! ===========================================
!  FreeGhostBuffer:
!  Author: Christoph Erath, Mark Taylor
!  Freed an ghostpoints communication buffer
! =========================================
   subroutine freeghostbuffertr(ghost)
!-------------------------------------------------------------------------------
   implicit none
!
   type(ghostbuffertr_t), intent(inout) :: ghost
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
!
   ghost%nbuf = 0
   ghost%nlyr = 0
   deallocate(ghost%buf)
   deallocate(ghost%receive)
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
#endif
   end subroutine freeghostbuffertr
! =========================================
! =========================================
!
!> @brief Pack edges of v into an edge buffer for boundary exchange.
!
!> This subroutine packs for one or more vertical layers into an edge
!! buffer. If the buffer associated with edge is not large enough to
!! hold all vertical layers you intent to pack, the method will
!! halt the program with a call to parallel_mod::FinPar(message=).
!! @param[in] edge Ghost Buffer into which the data will be packed.
!! This buffer must be previously allocated with initghostbufferfull().
!! @param[in] v The data to be packed.
!! @param[in] vlyr Number of vertical level coming into the subroutine
!! for packing for input v.
!! @param[in] kptr Vertical pointer to the place in the edge buffer where
!! data will be located.
! =========================================
   subroutine ghostvpackfull(edge, v, nc1, nc2, nc, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(ghostbuffer_t) :: edge
   integer, intent(in   ) :: vlyr
   real(real_kind), intent(in   ) :: v(nc1:nc2, nc1:nc2, vlyr)
   integer, intent(in   ) :: nc1, nc2, nc, kptr
   type(edgedescriptor_t), intent(in   ) :: desc
! Local variables
   integer :: i, k, ir, l, e
   integer :: is, ie, in, iw
!
   if (.not.threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
     threadsafe = .true.
   endif
! Example convenction for buffer to the north:
!  buf(:,,:,i,e)
!   each "edge" is a row of data (i=1,np) in the element
!     north most row of data goes into e=1
!     next row of data goes into e=2
!     ....
!     south most row of data goes into e=np
!   We need to pack this way to preserve the orientation
!   so the data can be unpacked correctly
! note: we think of buf as dimensioned buf(k,is,i,e)
! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
!
!
   is = desc%putmapp_ghost(south)
   ie = desc%putmapp_ghost(east)
   in = desc%putmapp_ghost(north)
   iw = desc%putmapp_ghost(west)
#if 1
   if (is>edge%nbuf) stop 'error is = '
   if (ie>edge%nbuf) stop 'error ie = '
   if (in>edge%nbuf) stop 'error in = '
   if (iw>edge%nbuf) stop 'error iw = '
   if (is<1) stop 'error is = 0'
   if (ie<1) stop 'error ie = 0'
   if (in<1) stop 'error in = 0'
   if (iw<1) stop 'error iw = 0'
#endif
!    print *,nc,is,ie,in,iw
   do k = 1,vlyr
     do e = 1, nc
       do i = 1, nc
         edge%buf(i,e,kptr+k,is) = v(i,e,k)
         edge%buf(i,e,kptr+k,ie) = v(nc-e+1,i,k)
         edge%buf(i,e,kptr+k,in) = v(i,nc-e+1,k)
         edge%buf(i,e,kptr+k,iw) = v(e,i,k)
       enddo
       enddo
   enddo
!  This is really kludgy way to setup the index reversals
!  But since it is so a rare event not real need to spend time optimizing
!  Check if the edge orientation of the recieving element is different
!  if it is,swap the order of data in the edge
   if (desc%reverse(south)) then
     is = desc%putmapp_ghost(south)
     do e = 1,nc
       do k = 1, vlyr
         do i = 1, nc
           ir = nc-i+1
           edge%buf(ir,e,kptr+k,is) = v(i,e,k)
         enddo
         enddo
     enddo
   endif
   if (desc%reverse(east)) then
     ie = desc%putmapp_ghost(east)
     do e = 1,nc
       do k = 1, vlyr
         do i = 1, nc
           ir = nc-i+1
           edge%buf(ir,e,kptr+k,ie) = v(nc-e+1,i,k)
         enddo
         enddo
     enddo
   endif
   if (desc%reverse(north)) then
   in = desc%putmapp_ghost(north)
     do e = 1, nc
       do k = 1, vlyr
         do i = 1, nc
           ir = nc-i+1
           edge%buf(ir,e,kptr+k,in) = v(i,nc-e+1,k)
         enddo
         enddo
         enddo
   endif
   if (desc%reverse(west)) then
     iw = desc%putmapp_ghost(west)
     do e = 1,nc
       do k = 1, vlyr
         do i = 1, nc
           ir = nc-i+1
           edge%buf(ir,e,kptr+k,iw) = v(e,i,k)
         enddo
         enddo
     enddo
   endif
! corners.  this is difficult because we dont know the orientaton
! of the corners,and this which (i,j) dimension maps to which dimension
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (l.ne.swest) stop 'ERROR1:swest ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
! edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(1,1 ,k)
         do e = 1, nc
           edge%buf(:,e,kptr+k,desc%putmapp_ghost(l)) = v(:,e,k)
         enddo
         enddo
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (l.ne.seast) stop 'ERROR1:seast ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
! edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(nc ,1 ,k)
         do e = 1, nc
           edge%buf(e,:,kptr+k,desc%putmapp_ghost(l)) = v(nc-e+1,:,k)
         enddo
         enddo
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (l.ne.neast) stop 'ERROR1:neast ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
!edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(nc ,nc,k)
         do e = 1, nc
           do i = 1, nc
             edge%buf(i,e,kptr+k,desc%putmapp_ghost(l)) = v(nc-i+1,nc-e+1,k)
           enddo
           enddo
           enddo
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (l.ne.nwest) stop 'ERROR1:nwest ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
!edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(1  ,nc,k)
         do e = 1, nc
           edge%buf(:,e,kptr+k,desc%putmapp_ghost(l)) = v(:,nc-e+1,k)
         enddo
         enddo
     endif
   enddo
   end subroutine ghostvpackfull
! ========================================
! edgeVunpack:
!
! Unpack edges from edge buffer into v...
! ========================================
   subroutine ghostvunpackfull(edge, v, nc1, nc2, nc, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(ghostbuffer_t), intent(in   ) :: edge
   integer, intent(in   ) :: vlyr
   real(real_kind), intent(inout) :: v(nc1:nc2, nc1:nc2, vlyr)
   integer, intent(in   ) :: kptr, nc1, nc2, nc
   type(edgedescriptor_t) :: desc
! Local
   logical, parameter :: useunroll = .true.
   integer :: i, k, l, e
   integer :: is, ie, in, iw, ic
   logical :: reverse
! make sure buffer is big enough:
!
   if ((nc2-nc1+1)<3*nc) then
     call finpar(message = "ghostvunpack:insufficient ghost cell region")
   endif
!
   threadsafe = .false.
   is = desc%getmapp_ghost(south)
   ie = desc%getmapp_ghost(east)
   in = desc%getmapp_ghost(north)
   iw = desc%getmapp_ghost(west)
! example for north buffer
! first row ('edge') goes in v(:,np+1,k)
! 2nd   row ('edge') goes in v(:,np+2,k)
! etc...
   do e = 1,nc
     do k = 1, vlyr
       do i = 1, nc
         v(i,1-e,k) = edge%buf(i,e,kptr+k,is)
         v(nc+e,i,k) = edge%buf(i,e,kptr+k,ie)
         v(i,nc+e,k) = edge%buf(i,e,kptr+k,in)
         v(1-e,i,k) = edge%buf(i,e,kptr+k,iw)
       enddo
       enddo
   enddo
! SWEST
   do l = swest,swest+max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do k = 1, vlyr
!v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
           do e = 1, nc
             do i = 1, nc
               v(1-e,1-i,k) = edge%buf(i,e,kptr+k,ic)
             enddo
             enddo
             enddo
             else
               do k = 1, vlyr
!v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                 do e = 1, nc
                   do i = 1, nc
                     v(1-e,1-i,k) = edge%buf(e,i,kptr+k,ic)
                   enddo
                   enddo
                   enddo
       endif
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do k = 1, vlyr
!v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
           do e = 1, nc
             do i = 1, nc
               v(nc+i,1-e,k) = edge%buf(e,i,kptr+k,ic)
             enddo
             enddo
             enddo
             else
               do k = 1, vlyr
!v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                 do e = 1, nc
                   do i = 1, nc
                     v(nc+i,1-e,k) = edge%buf(i,e,kptr+k,ic)
                   enddo
                   enddo
                   enddo
       endif
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do k = 1, vlyr
!v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
           do e = 1, nc
             do i = 1, nc
               v(nc+i,nc+e,k) = edge%buf(e,i,kptr+k,ic)
             enddo
             enddo
             enddo
             else
               do k = 1, vlyr
!v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                 do e = 1, nc
                   do i = 1, nc
                     v(nc+i,nc+e,k) = edge%buf(i,e,kptr+k,ic)
                   enddo
                   enddo
                   enddo
       endif
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do k = 1, vlyr
!v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
           do e = 1, nc
             do i = 1, nc
               v(1-i,nc+e,k) = edge%buf(e,i,kptr+k,ic)
             enddo
             enddo
             enddo
             else
               do k = 1, vlyr
!v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                 do e = 1, nc
                   do i = 1, nc
                     v(1-i,nc+e,k) = edge%buf(i,e,kptr+k,ic)
                   enddo
                   enddo
                   enddo
       endif
     endif
   enddo
   end subroutine ghostvunpackfull
! =========================================
! initGhostBuffer:
! Author: Christoph Erath
! create an Real based communication buffer
! npoints is the number of points on one side
! nhc is the deep of the ghost/halo zone
! =========================================
   subroutine initghostbuffer(ghost, nlyr, ntrac, nhc, npoints)
!-------------------------------------------------------------------------------
   use dimensions, only: nelemd, max_neigh_edges
   implicit none
!
   integer, intent(in   ) :: nlyr, ntrac, nhc, npoints
   type(ghostbuffertr_t), intent(  out) :: ghost
! Local variables
   integer :: nbuf
! make sure buffer is big enough:
!
   if (nhc>npoints) then
     call finpar(message = "intghostbuffer:halo region can not be larger then element size")
   endif
!
   if (ntrac<1) then
     call finpar(message = "intghostbuffer:you have to consider at least one tracer")
   endif
! sanity check for threading
!
   if (omp_get_num_threads()>1) then
     call finpar(message = 'ERROR:initghostbuffer must be called before threaded region')
   endif
!
   nbuf = max_neigh_edges*nelemd
#if (! defined ELEMENT_OPENMP)
!$OMP MASTER
#endif
   ghost%nlyr = nlyr
   ghost%ntrac = ntrac
   ghost%nbuf = nbuf
   allocate(ghost%buf(npoints,nhc,nlyr,ntrac,nbuf))
   ghost%buf = 0
   allocate(ghost%receive(npoints,nhc,nlyr,ntrac,nbuf))
   ghost%receive = 0
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif
   end subroutine initghostbuffer
! =========================================
! Christoph Erath
!> Packs the halo zone from v
! =========================================
   subroutine ghostvpack(edge, v, nhc, npoints, vlyr, ntrac, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
   use kiapsparallel, only: abortpar
!
   type(ghostbuffertr_t) :: edge
   integer, intent(in   ) :: vlyr
   integer, intent(in   ) :: ntrac
   real(real_kind), intent(in   ) :: v(1-nhc:npoints+nhc, 1-nhc:npoints+nhc, vlyr, ntrac)
   integer, intent(in   ) :: npoints, nhc, kptr
   type(edgedescriptor_t), intent(in   ) :: desc
! Local variables
   integer :: i, j, k, ir, l, itr
   integer :: is, ie, in, iw
!
   if (.not.threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
     threadsafe = .true.
   endif
! Example convenction for buffer to the north:
!  buf(:,,:,i,e)
!   each "edge" is a row of data (i=1,np) in the element
!     north most row of data goes into e=1
!     next row of data goes into e=2
!     ....
!     south most row of data goes into e=np
!   We need to pack this way to preserve the orientation
!   so the data can be unpacked correctly
! note: we think of buf as dimensioned buf(k,is,i,e)
! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
!
!
   is = desc%putmapp_ghost(south)
   ie = desc%putmapp_ghost(east)
   in = desc%putmapp_ghost(north)
   iw = desc%putmapp_ghost(west)
!    print *,nc,is,ie,in,iw
   do itr = 1,ntrac
     do k = 1, vlyr
       do i = 1, npoints
         do j = 1, nhc
           edge%buf(i,j,kptr+k,itr,is) = v(i,j,k,itr)
           edge%buf(i,j,kptr+k,itr,ie) = v(npoints-j+1,i,k,itr)
           edge%buf(i,j,kptr+k,itr,in) = v(i,npoints-j+1,k,itr)
           edge%buf(i,j,kptr+k,itr,iw) = v(j,i,k,itr)
         enddo
         enddo
         enddo
   enddo
!  This is really kludgy way to setup the index reversals
!  But since it is so a rare event not real need to spend time optimizing
!  Check if the edge orientation of the recieving element is different
!  if it is,swap the order of data in the edge
   if (desc%reverse(south)) then
!      is = desc%putmapP_ghost(south)
     do itr = 1, ntrac
       do k = 1, vlyr
         do i = 1, npoints
           do j = 1, nhc
             ir = npoints-i+1
             edge%buf(ir,j,kptr+k,itr,is) = v(i,j,k,itr)
           enddo
           enddo
           enddo
           enddo
   endif
   if (desc%reverse(east)) then
!      ie = desc%putmapP_ghost(east)
     do itr = 1, ntrac
       do k = 1, vlyr
         do i = 1, npoints
           do j = 1, nhc
             ir = npoints-i+1
             edge%buf(ir,j,kptr+k,itr,ie) = v(npoints-j+1,i,k,itr)
           enddo
           enddo
           enddo
           enddo
   endif
   if (desc%reverse(north)) then
!      in = desc%putmapP_ghost(north)
     do itr = 1, ntrac
       do k = 1, vlyr
         do i = 1, npoints
           do j = 1, nhc
             ir = npoints-i+1
             edge%buf(ir,j,kptr+k,itr,in) = v(i,npoints-j+1,k,itr)
           enddo
           enddo
           enddo
           enddo
   endif
   if (desc%reverse(west)) then
!      iw = desc%putmapP_ghost(west)
     do itr = 1, ntrac
       do k = 1, vlyr
         do i = 1, npoints
           do j = 1, nhc
             ir = npoints-i+1
             edge%buf(ir,j,kptr+k,itr,iw) = v(j,i,k,itr)
           enddo
           enddo
           enddo
           enddo
   endif
! corners.  this is difficult because we dont know the orientaton
! of the corners,and this which (i,j) dimension maps to which dimension
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (l.ne.swest) stop 'ERROR2:swest ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do itr = 1, ntrac
         do k = 1, vlyr
! edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(1,1 ,k)
           do i = 1, nhc
             do j = 1, nhc
               edge%buf(i,j,kptr+k,itr,desc%putmapp_ghost(l)) = v(i,j,k,itr)
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (l.ne.seast) stop 'ERROR2:seast ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do itr = 1, ntrac
         do k = 1, vlyr
! edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(nc ,1 ,k)
           do i = 1, nhc
             do j = 1, nhc
               edge%buf(i,j,kptr+k,itr,desc%putmapp_ghost(l)) = v(npoints-i+1,j,k,itr)
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (l.ne.neast) stop 'ERROR2:neast ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do itr = 1, ntrac
         do k = 1, vlyr
!edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(nc ,nc,k)
           do i = 1, nhc
             do j = 1, nhc
               edge%buf(i,j,kptr+k,itr,desc%putmapp_ghost(l)) = v(npoints-i+1,npoints-j+1,k,itr)
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (l.ne.nwest) stop 'ERROR2:nwest ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do itr = 1, ntrac
         do k = 1, vlyr
!edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(1  ,nc,k)
           do i = 1, nhc
             do j = 1, nhc
               edge%buf(i,j,kptr+k,itr,desc%putmapp_ghost(l)) = v(i,npoints-j+1,k,itr)
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
   end subroutine ghostvpack
! ========================================
! Christoph Erath
!
! Unpack the halo zone into v
! ========================================
   subroutine ghostvunpack(edge, v, nhc, npoints, vlyr, ntrac, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(ghostbuffertr_t), intent(in   ) :: edge
   integer, intent(in   ) :: vlyr
   integer, intent(in   ) :: ntrac
   real(real_kind), intent(inout) :: v(1-nhc:npoints+nhc, 1-nhc:npoints+nhc, vlyr, ntrac)
   integer, intent(in   ) :: kptr, nhc, npoints
   type(edgedescriptor_t) :: desc
! Local
   logical, parameter :: useunroll = .true.
   integer :: i, j, k, l, itr
   integer :: is, ie, in, iw, ic
   logical :: reverse
   real(real_kind) :: nan =-1.0
!
   nan = sqrt(nan)
   threadsafe = .false.
   is = desc%getmapp_ghost(south)
   ie = desc%getmapp_ghost(east)
   in = desc%getmapp_ghost(north)
   iw = desc%getmapp_ghost(west)
! example for north buffer
! first row ('edge') goes in v(:,np+1,k)
! 2nd   row ('edge') goes in v(:,np+2,k)
! etc...
   do itr = 1,ntrac
     do k = 1, vlyr
       do i = 1, npoints
         do j = 1, nhc
           v(i,1-j,k,itr) = edge%buf(i,j,kptr+k,itr,is)
           v(npoints+j,i,k,itr) = edge%buf(i,j,kptr+k,itr,ie)
           v(i,npoints+j,k,itr) = edge%buf(i,j,kptr+k,itr,in)
           v(1-j,i,k,itr) = edge%buf(i,j,kptr+k,itr,iw)
         enddo
         enddo
         enddo
   enddo
! SWEST
   do l = swest,swest+max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do itr = 1, ntrac
           do k = 1, vlyr
!v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
             do j = 1, nhc
               do i = 1, nhc
                 v(1-j,1-i,k,itr) = edge%buf(i,j,kptr+k,itr,ic)
               enddo
               enddo
               enddo
               enddo
               else
                 do itr = 1, ntrac
                   do k = 1, vlyr
!v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                     do j = 1, nhc
                       do i = 1, nhc
                         v(1-j,1-i,k,itr) = edge%buf(j,i,kptr+k,itr,ic)
                       enddo
                       enddo
                       enddo
                       enddo
       endif
     else
       do itr = 1, ntrac
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(1-i,1-j,k,itr) = nan
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do itr = 1, ntrac
           do k = 1, vlyr
!v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
             do j = 1, nhc
               do i = 1, nhc
                 v(npoints+i,1-j,k,itr) = edge%buf(j,i,kptr+k,itr,ic)
               enddo
               enddo
               enddo
               enddo
               else
                 do itr = 1, ntrac
                   do k = 1, vlyr
!v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                     do j = 1, nhc
                       do i = 1, nhc
                         v(npoints+i,1-j,k,itr) = edge%buf(i,j,kptr+k,itr,ic)
                       enddo
                       enddo
                       enddo
                       enddo
       endif
     else
       do itr = 1, ntrac
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(npoints+i,1-j,k,itr) = nan
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do itr = 1, ntrac
           do k = 1, vlyr
!v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
             do j = 1, nhc
               do i = 1, nhc
                 v(npoints+i,npoints+j,k,itr) = edge%buf(j,i,kptr+k,itr,ic)
               enddo
               enddo
               enddo
               enddo
               else
                 do itr = 1, ntrac
                   do k = 1, vlyr
!v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                     do j = 1, nhc
                       do i = 1, nhc
                         v(npoints+i,npoints+j,k,itr) = edge%buf(i,j,kptr+k,itr,ic)
                       enddo
                       enddo
                       enddo
                       enddo
       endif
     else
       do itr = 1, ntrac
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(npoints+i,npoints+j,k,itr) = nan
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do itr = 1, ntrac
           do k = 1, vlyr
!v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
             do j = 1, nhc
               do i = 1, nhc
                 v(1-i,npoints+j,k,itr) = edge%buf(j,i,kptr+k,itr,ic)
               enddo
               enddo
               enddo
               enddo
               else
                 do itr = 1, ntrac
                   do k = 1, vlyr
!v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                     do j = 1, nhc
                       do i = 1, nhc
                         v(1-i,npoints+j,k,itr) = edge%buf(i,j,kptr+k,itr,ic)
                       enddo
                       enddo
                       enddo
                       enddo
       endif
     else
       do itr = 1, ntrac
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(1-i,npoints+j,k,itr) = nan
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
   end subroutine ghostvunpack
!for VLAG
! =========================================
! Christoph Erath
!> Packs the halo zone from v
! =========================================
! NOTE: I have to give timelevels as argument, because element_mod is not compiled first
! and the array call has to be done in this way because of performance reasons!!!
   subroutine ghostvpack_new(edge, v, nhc, npoints, vlyr, ntrac, kptr, tn0, timelevels, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem, ntrac_d
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(ghostbuffertr_t) :: edge
   integer, intent(in   ) :: vlyr
   integer, intent(in   ) :: ntrac
   integer, intent(in   ) :: npoints, nhc, kptr, tn0, timelevels
   real(real_kind), intent(in   ) :: v(1-nhc:npoints+nhc, 1-nhc:npoints+nhc, vlyr, ntrac_d, timelevels)
   type(edgedescriptor_t), intent(in   ) :: desc
! Local variables
   integer :: i, j, k, ir, l, itr
   integer :: is, ie, in, iw
!
   if (.not.threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
     threadsafe = .true.
   endif
! Example convenction for buffer to the north:
!  buf(:,,:,i,e)
!   each "edge" is a row of data (i=1,np) in the element
!     north most row of data goes into e=1
!     next row of data goes into e=2
!     ....
!     south most row of data goes into e=np
!   We need to pack this way to preserve the orientation
!   so the data can be unpacked correctly
! note: we think of buf as dimensioned buf(k,is,i,e)
! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
!
!
   is = desc%putmapp_ghost(south)
   ie = desc%putmapp_ghost(east)
   in = desc%putmapp_ghost(north)
   iw = desc%putmapp_ghost(west)
!    print *,nc,is,ie,in,iw
   do itr = 1,ntrac
     do k = 1, vlyr
       do i = 1, npoints
         do j = 1, nhc
           edge%buf(i,j,kptr+k,itr,is) = v(i,j,k,itr,tn0)
           edge%buf(i,j,kptr+k,itr,ie) = v(npoints-j+1,i,k,itr,tn0)
           edge%buf(i,j,kptr+k,itr,in) = v(i,npoints-j+1,k,itr,tn0)
           edge%buf(i,j,kptr+k,itr,iw) = v(j,i,k,itr,tn0)
         enddo
         enddo
         enddo
   enddo
!  This is really kludgy way to setup the index reversals
!  But since it is so a rare event not real need to spend time optimizing
!  Check if the edge orientation of the recieving element is different
!  if it is,swap the order of data in the edge
   if (desc%reverse(south)) then
!      is = desc%putmapP_ghost(south)
     do itr = 1, ntrac
       do k = 1, vlyr
         do i = 1, npoints
           do j = 1, nhc
             ir = npoints-i+1
             edge%buf(ir,j,kptr+k,itr,is) = v(i,j,k,itr,tn0)
           enddo
           enddo
           enddo
           enddo
   endif
   if (desc%reverse(east)) then
!      ie = desc%putmapP_ghost(east)
     do itr = 1, ntrac
       do k = 1, vlyr
         do i = 1, npoints
           do j = 1, nhc
             ir = npoints-i+1
             edge%buf(ir,j,kptr+k,itr,ie) = v(npoints-j+1,i,k,itr,tn0)
           enddo
           enddo
           enddo
           enddo
   endif
   if (desc%reverse(north)) then
!      in = desc%putmapP_ghost(north)
     do itr = 1, ntrac
       do k = 1, vlyr
         do i = 1, npoints
           do j = 1, nhc
             ir = npoints-i+1
             edge%buf(ir,j,kptr+k,itr,in) = v(i,npoints-j+1,k,itr,tn0)
           enddo
           enddo
           enddo
           enddo
   endif
   if (desc%reverse(west)) then
!      iw = desc%putmapP_ghost(west)
     do itr = 1, ntrac
       do k = 1, vlyr
         do i = 1, npoints
           do j = 1, nhc
             ir = npoints-i+1
             edge%buf(ir,j,kptr+k,itr,iw) = v(j,i,k,itr,tn0)
           enddo
           enddo
           enddo
           enddo
   endif
! corners.  this is difficult because we dont know the orientaton
! of the corners,and this which (i,j) dimension maps to which dimension
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (l.ne.swest) call abortpar('ERROR2:swest ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapp_ghost(l)/=-1) then
       do itr = 1, ntrac
         do k = 1, vlyr
! edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(1,1 ,k)
           do i = 1, nhc
             do j = 1, nhc
               edge%buf(i,j,kptr+k,itr,desc%putmapp_ghost(l)) = v(i,j,k,itr,tn0)
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (l.ne.seast) call abortpar('ERROR2:seast ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapp_ghost(l)/=-1) then
       do itr = 1, ntrac
         do k = 1, vlyr
! edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(nc ,1 ,k)
           do i = 1, nhc
             do j = 1, nhc
               edge%buf(i,j,kptr+k,itr,desc%putmapp_ghost(l)) = v(npoints-i+1,j,k,itr,tn0)
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (l.ne.neast) call abortpar('ERROR2:neast ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapp_ghost(l)/=-1) then
       do itr = 1, ntrac
         do k = 1, vlyr
!edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(nc ,nc,k)
           do i = 1, nhc
             do j = 1, nhc
               edge%buf(i,j,kptr+k,itr,desc%putmapp_ghost(l)) = v(npoints-i+1,npoints-j+1,k,itr,tn0)
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (l.ne.nwest) call abortpar('ERROR2:nwest ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapp_ghost(l)/=-1) then
       do itr = 1, ntrac
         do k = 1, vlyr
!edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(1  ,nc,k)
           do i = 1, nhc
             do j = 1, nhc
               edge%buf(i,j,kptr+k,itr,desc%putmapp_ghost(l)) = v(i,npoints-j+1,k,itr,tn0)
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
   end subroutine ghostvpack_new
! ========================================
! Christoph Erath
!
! Unpack the halo zone into v
! ========================================
! NOTE: I have to give timelevels as argument, because element_mod is not compiled first
! and the array call has to be done in this way because of performance reasons!!!
   subroutine ghostvunpack_new(edge, v, nhc, npoints, vlyr, ntrac, kptr, tn0, timelevels, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem, ntrac_d
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(ghostbuffertr_t), intent(in   ) :: edge
   integer, intent(in   ) :: vlyr
   integer, intent(in   ) :: ntrac
   integer, intent(in   ) :: kptr, nhc, npoints, tn0, timelevels
   real(real_kind), intent(inout) :: v(1-nhc:npoints+nhc, 1-nhc:npoints+nhc, vlyr, ntrac_d, timelevels)
   type(edgedescriptor_t) :: desc
! Local
   logical, parameter :: useunroll = .true.
   integer :: i, j, k, l, itr
   integer :: is, ie, in, iw, ic
   logical :: reverse
   real(real_kind) :: nan =-1.0
!
   nan = sqrt(nan)
   threadsafe = .false.
   is = desc%getmapp_ghost(south)
   ie = desc%getmapp_ghost(east)
   in = desc%getmapp_ghost(north)
   iw = desc%getmapp_ghost(west)
! example for north buffer
! first row ('edge') goes in v(:,np+1,k)
! 2nd   row ('edge') goes in v(:,np+2,k)
! etc...
   do itr = 1,ntrac
     do k = 1, vlyr
       do i = 1, npoints
         do j = 1, nhc
           v(i,1-j,k,itr,tn0) = edge%buf(i,j,kptr+k,itr,is)
           v(npoints+j,i,k,itr,tn0) = edge%buf(i,j,kptr+k,itr,ie)
           v(i,npoints+j,k,itr,tn0) = edge%buf(i,j,kptr+k,itr,in)
           v(1-j,i,k,itr,tn0) = edge%buf(i,j,kptr+k,itr,iw)
         enddo
         enddo
         enddo
   enddo
! SWEST
   do l = swest,swest+max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do itr = 1, ntrac
           do k = 1, vlyr
!v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
             do j = 1, nhc
               do i = 1, nhc
                 v(1-j,1-i,k,itr,tn0) = edge%buf(i,j,kptr+k,itr,ic)
               enddo
               enddo
               enddo
               enddo
               else
                 do itr = 1, ntrac
                   do k = 1, vlyr
!v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                     do j = 1, nhc
                       do i = 1, nhc
                         v(1-j,1-i,k,itr,tn0) = edge%buf(j,i,kptr+k,itr,ic)
                       enddo
                       enddo
                       enddo
                       enddo
       endif
     else
       do itr = 1, ntrac
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(1-i,1-j,k,itr,tn0) = nan
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do itr = 1, ntrac
           do k = 1, vlyr
!v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
             do j = 1, nhc
               do i = 1, nhc
                 v(npoints+i,1-j,k,itr,tn0) = edge%buf(j,i,kptr+k,itr,ic)
               enddo
               enddo
               enddo
               enddo
               else
                 do itr = 1, ntrac
                   do k = 1, vlyr
!v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                     do j = 1, nhc
                       do i = 1, nhc
                         v(npoints+i,1-j,k,itr,tn0) = edge%buf(i,j,kptr+k,itr,ic)
                       enddo
                       enddo
                       enddo
                       enddo
       endif
     else
       do itr = 1, ntrac
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(npoints+i,1-j,k,itr,tn0) = nan
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do itr = 1, ntrac
           do k = 1, vlyr
!v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
             do j = 1, nhc
               do i = 1, nhc
                 v(npoints+i,npoints+j,k,itr,tn0) = edge%buf(j,i,kptr+k,itr,ic)
               enddo
               enddo
               enddo
               enddo
               else
                 do itr = 1, ntrac
                   do k = 1, vlyr
!v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                     do j = 1, nhc
                       do i = 1, nhc
                         v(npoints+i,npoints+j,k,itr,tn0) = edge%buf(i,j,kptr+k,itr,ic)
                       enddo
                       enddo
                       enddo
                       enddo
       endif
     else
       do itr = 1, ntrac
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(npoints+i,npoints+j,k,itr,tn0) = nan
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do itr = 1, ntrac
           do k = 1, vlyr
!v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
             do j = 1, nhc
               do i = 1, nhc
                 v(1-i,npoints+j,k,itr,tn0) = edge%buf(j,i,kptr+k,itr,ic)
               enddo
               enddo
               enddo
               enddo
               else
                 do itr = 1, ntrac
                   do k = 1, vlyr
!v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                     do j = 1, nhc
                       do i = 1, nhc
                         v(1-i,npoints+j,k,itr,tn0) = edge%buf(i,j,kptr+k,itr,ic)
                       enddo
                       enddo
                       enddo
                       enddo
       endif
     else
       do itr = 1, ntrac
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(1-i,npoints+j,k,itr,tn0) = nan
             enddo
             enddo
             enddo
             enddo
     endif
   enddo
   end subroutine ghostvunpack_new
! =================================================================================
! GHOSTVPACK2D
! AUTHOR: Christoph Erath (from a subroutine of Mark Taylor, ghostvpack)
! Pack edges of v into an ghost buffer for boundary exchange.
!
! This subroutine packs for one vertical layers into an ghost
! buffer. It is for cartesian points (v is only two dimensional).
! If the buffer associated with edge is not large enough to
! hold all vertical layers you intent to pack, the method will
! halt the program with a call to parallel_mod::FinPar(message=).
! INPUT:
! - ghost Buffer into which the data will be packed.
!   This buffer must be previously allocated with initGhostBuffer().
! - v The data to be packed.
! - nhc deep of ghost/halo zone
! - npoints number of points on on side
! - kptr Vertical pointer to the place in the edge buffer where
! data will be located.
! =================================================================================
   subroutine ghostvpack2d(ghost, v, nhc, npoints, kptr, desc, cellcenter)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(ghostbuffer_t) :: ghost
   real(real_kind), intent(in   ) :: v(1-nhc:npoints+nhc, 1-nhc:npoints+nhc)
   integer, intent(in   ) :: nhc, npoints, kptr
   type(edgedescriptor_t), intent(in   ) :: desc
   logical, optional, intent(in   ) :: cellcenter
! Local variables
   integer :: i, j, ir, l, e
   integer :: idxshift
   integer :: is, ie, in, iw
!
   if (present(cellcenter)) then
     idxshift = 1
   else
     idxshift = 0
   endif
!
   if (.not.threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
     threadsafe = .true.
   endif
! Example convenction for buffer to the north:
!  buf(:,,:,i,e)
!   each "edge" is a row of data (i=1,np) in the element
!     north most row of data goes into e=1
!     next row of data goes into e=2
!     ....
!     south most row of data goes into e=np
!   We need to pack this way to preserve the orientation
!   so the data can be unpacked correctly
! note: we think of buf as dimensioned buf(k,is,i,e)
! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
!
!
   is = desc%putmapp_ghost(south)
   ie = desc%putmapp_ghost(east)
   in = desc%putmapp_ghost(north)
   iw = desc%putmapp_ghost(west)
   do i = 1,npoints
     do j = 1, nhc
       ghost%buf(i,j,kptr,is) = v(i,j+1-idxshift)
       ghost%buf(i,j,kptr,ie) = v(npoints-j+idxshift,i)
       ghost%buf(i,j,kptr,in) = v(i,npoints-j+idxshift)
       ghost%buf(i,j,kptr,iw) = v(j+1-idxshift,i)
     enddo
   enddo
!  This is really kludgy way to setup the index reversals
!  But since it is so a rare event not real need to spend time optimizing
!  Check if the edge orientation of the recieving element is different
!  if it is,swap the order of data in the edge
   if (desc%reverse(south)) then
!        is = desc%putmapP_ghost(south)
     do i = 1, npoints
       do j = 1, nhc
         ir = npoints-i+1
         ghost%buf(ir,j,kptr,is) = v(i,j+1-idxshift)
       enddo
       enddo
   endif
   if (desc%reverse(east)) then
!        ie = desc%putmapP_ghost(east)
     do i = 1, npoints
       do j = 1, nhc
         ir = npoints-i+1
         ghost%buf(ir,j,kptr,ie) = v(npoints-j+idxshift,i)
       enddo
       enddo
   endif
   if (desc%reverse(north)) then
!        in = desc%putmapP_ghost(north)
     do i = 1, npoints
       do j = 1, nhc
         ir = npoints-i+1
         ghost%buf(ir,j,kptr,in) = v(i,npoints-j+idxshift)
       enddo
       enddo
   endif
   if (desc%reverse(west)) then
!        iw = desc%putmapP_ghost(west)
     do i = 1, npoints
       do j = 1, nhc
         ir = npoints-i+1
         ghost%buf(ir,j,kptr,iw) = v(j+1-idxshift,i)
       enddo
       enddo
   endif
! corners.  this is difficult because we dont know the orientaton
! of the corners,and this which (i,j) dimension maps to which dimension
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (l.ne.swest) stop 'ERROR3:swest ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do i = 1, nhc
         do j = 1, nhc
           ghost%buf(i,j,kptr,desc%putmapp_ghost(l)) = v(i+1-idxshift,j+1-idxshift)
         enddo
         enddo
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (l.ne.seast) stop 'ERROR3:seast ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do i = 1, nhc
         do j = 1, nhc
           ghost%buf(i,j,kptr,desc%putmapp_ghost(l)) = v(npoints-i+idxshift,j+1-idxshift)
         enddo
         enddo
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (l.ne.neast) stop 'ERROR3:neast ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do i = 1, nhc
         do j = 1, nhc
           ghost%buf(i,j,kptr,desc%putmapp_ghost(l)) = v(npoints-i+idxshift,npoints-j+idxshift)
         enddo
         enddo
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (l.ne.nwest) stop 'ERROR3:nwest ghost cell update requires ne>0 cubed-sphere mesh'
     if (desc%putmapp_ghost(l)/=-1) then
       do i = 1, nhc
         do j = 1, nhc
           ghost%buf(i,j,kptr,desc%putmapp_ghost(l)) = v(i+1-idxshift,npoints-j+idxshift)
         enddo
         enddo
     endif
   enddo
   end subroutine ghostvpack2d
! =================================================================================
! GHOSTVUNPACK2D
! AUTHOR: Christoph Erath (from a subroutine of Mark Taylor, ghostVunpack)
! Unpack ghost points from ghost buffer into v...
! It is for cartesian points (v is only two dimensional).
! INPUT SAME arguments as for GHOSTVPACK2d
! =================================================================================
   subroutine ghostvunpack2d(ghost, v, nhc, npoints, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(ghostbuffer_t), intent(in   ) :: ghost
   real(real_kind), intent(inout) :: v(1-nhc:npoints+nhc, 1-nhc:npoints+nhc)
   integer, intent(in   ) :: kptr, nhc, npoints
   type(edgedescriptor_t) :: desc
! Local
   logical, parameter :: useunroll = .true.
   integer :: i, j, l
   integer :: is, ie, in, iw, ic
   logical :: reverse
   real(real_kind) :: nan =-1.0
!
   nan = sqrt(nan)
   threadsafe = .false.
   is = desc%getmapp_ghost(south)
   ie = desc%getmapp_ghost(east)
   in = desc%getmapp_ghost(north)
   iw = desc%getmapp_ghost(west)
! example for north buffer
! first row ('edge') goes in v(:,np+1)
! 2nd   row ('edge') goes in v(:,np+2)
! etc...
   do i = 1,npoints
     do j = 1, nhc
       v(i,1-j) = ghost%buf(i,j,kptr,is)
       v(npoints+j,i) = ghost%buf(i,j,kptr,ie)
       v(i,npoints+j) = ghost%buf(i,j,kptr,in)
       v(1-j,i) = ghost%buf(i,j,kptr,iw)
     enddo
   enddo
! SWEST
   do l = swest,swest+max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do j = 1, nhc
           do i = 1, nhc
             v(1-i,1-j) = ghost%buf(j,i,kptr,ic)
           enddo
           enddo
           else
             do j = 1, nhc
               do i = 1, nhc
                 v(1-i,1-j) = ghost%buf(i,j,kptr,ic)
               enddo
               enddo
       endif
     else
       do j = 1, nhc
         do i = 1, nhc
           v(1-i,1-j) = nan
         enddo
         enddo
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do j = 1, nhc
           do i = 1, nhc
             v(npoints+i,1-j) = ghost%buf(j,i,kptr,ic)
           enddo
           enddo
           else
             do j = 1, nhc
               do i = 1, nhc
                 v(npoints+i,1-j) = ghost%buf(i,j,kptr,ic)
               enddo
               enddo
       endif
     else
       do j = 1, nhc
         do i = 1, nhc
           v(npoints+i,1-j) = nan
         enddo
         enddo
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do j = 1, nhc
           do i = 1, nhc
             v(npoints+i,npoints+j) = ghost%buf(j,i,kptr,ic)
           enddo
           enddo
           else
             do j = 1, nhc
               do i = 1, nhc
                 v(npoints+i,npoints+j) = ghost%buf(i,j,kptr,ic)
               enddo
               enddo
       endif
     else
       do j = 1, nhc
         do i = 1, nhc
           v(npoints+i,npoints+j) = nan
         enddo
         enddo
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do j = 1, nhc
           do i = 1, nhc
             v(1-i,npoints+j) = ghost%buf(j,i,kptr,ic)
           enddo
           enddo
           else
             do j = 1, nhc
               do i = 1, nhc
                 v(1-i,npoints+j) = ghost%buf(i,j,kptr,ic)
               enddo
               enddo
       endif
     else
       do j = 1, nhc
         do i = 1, nhc
           v(1-i,npoints+j) = nan
         enddo
         enddo
     endif
   enddo
   end subroutine ghostvunpack2d
! =========================================
! initGhostBuffer3d:
! Author: James Overfelt
! create an Real based communication buffer
! npoints is the number of points on one side
! nhc is the deep of the ghost/halo zone
! =========================================
   subroutine initghostbuffer3d(ghost, nlyr, np, nhc)
!-------------------------------------------------------------------------------
   use dimensions, only: nelemd, max_neigh_edges
   implicit none
!
   integer, intent(in   ) :: nlyr, nhc, np
   type(ghostbuffer3d_t), intent(  out) :: ghost
! Local variables
   integer :: nbuf
! sanity check for threading
!
   if (omp_get_num_threads()>1) then
     call finpar(message = 'ERROR:initghostbuffer must be called before threaded region')
   endif
!
   nbuf = max_neigh_edges*nelemd
#if (! defined ELEMENT_OPENMP)
!$OMP MASTER
#endif
   ghost%nlyr = nlyr
   ghost%nhc = nhc
   ghost%np = np
   ghost%nbuf = nbuf
   allocate(ghost%buf(np,(nhc+1),nlyr,nbuf))
   allocate(ghost%receive(np,(nhc+1),nlyr,nbuf))
   ghost%buf = 0
   ghost%receive = 0
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif
   end subroutine initghostbuffer3d
! =================================================================================
! GHOSTVPACK3D
! AUTHOR: James Overfelt (from a subroutine of Christoph Erath, ghostvpack2D)
! Pack edges of v into an ghost buffer for boundary exchange.
!
! This subroutine packs for many vertical layers into an ghost
! buffer.
! If the buffer associated with edge is not large enough to
! hold all vertical layers you intent to pack, the method will
! halt the program with a call to parallel_mod::FinPar(message=).
! INPUT:
! - ghost Buffer into which the data will be packed.
!   This buffer must be previously allocated with initGhostBuffer().
! - v The data to be packed.
! - nhc deep of ghost/halo zone
! - npoints number of points on on side
! - kptr Vertical pointer to the place in the edge buffer where
! data will be located.
! =================================================================================
   subroutine ghostvpack3d(ghost, v, vlyr, kptr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
   use kiapsparallel, only: abortpar
!
   type(ghostbuffer3d_t) :: ghost
   integer, intent(in   ) :: kptr, vlyr
   real(real_kind), intent(in   ) :: v(ghost%np, ghost%np, vlyr)
   type(edgedescriptor_t), intent(in   ) :: desc
   integer :: nhc, np
! Local variables
   integer :: i, j, k, ir, l, e
   integer :: is, ie, in, iw
!
   if (.not.threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
     threadsafe = .true.
   endif
! Example convenction for buffer to the north:
!  buf(:,,:,i,e)
!   each "edge" is a row of data (i=1,np) in the element
!     north most row of data goes into e=1
!     next row of data goes into e=2
!     ....
!     south most row of data goes into e=np
!   We need to pack this way to preserve the orientation
!   so the data can be unpacked correctly
! note: we think of buf as dimensioned buf(k,is,i,e)
! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
!
!
   nhc = ghost%nhc
   np = ghost%np
   is = desc%putmapp_ghost(south)
   ie = desc%putmapp_ghost(east)
   in = desc%putmapp_ghost(north)
   iw = desc%putmapp_ghost(west)
   do k = 1,vlyr
     do j = 1, nhc
       do i = 1, np
         ghost%buf(i,j,kptr+k,is) = v(i,j+1,k)
         ghost%buf(i,j,kptr+k,ie) = v(np-j,i,k)
         ghost%buf(i,j,kptr+k,in) = v(i,np-j,k)
         ghost%buf(i,j,kptr+k,iw) = v(j+1,i,k)
       enddo
       enddo
   enddo
!  This is really kludgy way to setup the index reversals
!  But since it is so a rare event not real need to spend time optimizing
!  Check if the edge orientation of the recieving element is different
!  if it is,swap the order of data in the edge
   if (desc%reverse(south)) then
     do k = 1, vlyr
       do j = 1, nhc
         do i = 1, np
           ir = np-i+1
           ghost%buf(ir,j,kptr+k,is) = v(i,j+1,k)
         enddo
         enddo
         enddo
   endif
   if (desc%reverse(east)) then
     do k = 1, vlyr
       do j = 1, nhc
         do i = 1, np
           ir = np-i+1
           ghost%buf(ir,j,kptr+k,ie) = v(np-j,i,k)
         enddo
         enddo
         enddo
   endif
   if (desc%reverse(north)) then
     do k = 1, vlyr
       do j = 1, nhc
         do i = 1, np
           ir = np-i+1
           ghost%buf(ir,j,kptr+k,in) = v(i,np-j,k)
         enddo
         enddo
         enddo
   endif
   if (desc%reverse(west)) then
     do k = 1, vlyr
       do j = 1, nhc
         do i = 1, np
           ir = np-i+1
           ghost%buf(ir,j,kptr+k,iw) = v(j+1,i,k)
         enddo
         enddo
         enddo
   endif
! corners.  this is difficult because we dont know the orientaton
! of the corners,and this which (i,j) dimension maps to which dimension
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
         do j = 1, nhc+1
           do i = 1, nhc+1
             ghost%buf(i,j,kptr+k,desc%putmapp_ghost(l)) = v(i,j,k)
           enddo
           enddo
           enddo
           endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
         do j = 1, nhc+1
           do i = 1, nhc+1
             ghost%buf(i,j,kptr+k,desc%putmapp_ghost(l)) = v(np-i+1,j,k)
           enddo
           enddo
           enddo
           endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
         do j = 1, nhc+1
           do i = 1, nhc+1
             ghost%buf(i,j,kptr+k,desc%putmapp_ghost(l)) = v(np-i+1,np-j+1,k)
           enddo
           enddo
           enddo
           endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
         do j = 1, nhc+1
           do i = 1, nhc+1
             ghost%buf(i,j,kptr+k,desc%putmapp_ghost(l)) = v(i,np-j+1,k)
           enddo
           enddo
           enddo
           endif
   enddo
   end subroutine ghostvpack3d
! =================================================================================
! GHOSTVUNPACK3D
! AUTHOR: James Overfelt (from a subroutine of Christoph Erath, ghostVunpack2d)
! Unpack ghost points from ghost buffer into v...
! It is for cartesian points (v is only two dimensional).
! INPUT SAME arguments as for GHOSTVPACK
! =================================================================================
   subroutine ghostvunpack3d(g, v, vlyr, kptr, desc, sw, se, nw, ne, mult)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(ghostbuffer3d_t), intent(in   ) :: g
   integer, intent(in   ) :: kptr, vlyr
   real(real_kind), intent(inout) :: v(1-g%nhc:g%np+g%nhc, 1-g%nhc:g%np+g%nhc, vlyr)
   integer, intent(  out) :: mult(5:8)
   real(real_kind), intent(  out) :: sw(1-g%nhc:1, 1-g%nhc:1, vlyr, max_corner_elem-1)
   real(real_kind), intent(  out) :: se(g%np:g%np+g%nhc, 1-g%nhc:1, vlyr, max_corner_elem-1)
   real(real_kind), intent(  out) :: ne(g%np:g%np+g%nhc, g%np:g%np+g%nhc, vlyr, max_corner_elem-1)
   real(real_kind), intent(  out) :: nw(1-g%nhc:1, g%np:g%np+g%nhc, vlyr, max_corner_elem-1)
   type(edgedescriptor_t) :: desc
   integer :: nhc, np
! Local
   logical, parameter :: useunroll = .true.
   integer :: i, j, k, l
   integer :: is, ie, in, iw, ic
   logical :: reverse
   real(real_kind) :: nan =-1.0
!
   nan = sqrt(nan)
   threadsafe = .false.
   nhc = g%nhc
   np = g%np
   is = desc%getmapp_ghost(south)
   ie = desc%getmapp_ghost(east)
   in = desc%getmapp_ghost(north)
   iw = desc%getmapp_ghost(west)
! fill in optional values with NaN
   do k = 1,vlyr
     do j = 1, nhc
       do i = 1, nhc
         v(1-i,1-j,k) = nan
         v(np+i,1-j,k) = nan
         v(np+i,np+j,k) = nan
         v(1-i,np+j,k) = nan
       enddo
       enddo
   enddo
! example for north buffer
! first row ('edge') goes in v(:,np+1)
! 2nd   row ('edge') goes in v(:,np+2)
! etc...
   do k = 1,vlyr
     do j = 1, nhc
       do i = 1, np
         v(i,1-j,k) = g%buf(i,j,kptr+k,is)
         v(np+j,i,k) = g%buf(i,j,kptr+k,ie)
         v(i,np+j,k) = g%buf(i,j,kptr+k,in)
         v(1-j,i,k) = g%buf(i,j,kptr+k,iw)
       enddo
       enddo
   enddo
! four sides are always just one
   mult(swest) = 0
   mult(seast) = 0
   mult(neast) = 0
   mult(nwest) = 0
! SWEST
   do l = swest,swest+max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (mult(swest).eq.0) then
         if (reverse) then
           do k = 1, vlyr
             do j = 1, nhc
               do i = 1, nhc
                 v(1-i,1-j,k) = g%buf(j+1,i+1,kptr+k,ic)
               enddo
               enddo
               enddo
               else
                 do k = 1, vlyr
                   do j = 1, nhc
                     do i = 1, nhc
                       v(1-i,1-j,k) = g%buf(i+1,j+1,kptr+k,ic)
                     enddo
                     enddo
                     enddo
                     endif
                     else
                       if (reverse) then
                         do k = 1, vlyr
                           do j = 0, nhc
                             do i = 0, nhc
                               sw(1-i,1-j,k,mult(swest)) = g%buf(j+1,i+1,kptr+k,ic)
                             enddo
                             enddo
                             enddo
                             else
                               do k = 1, vlyr
                                 do j = 0, nhc
                                   do i = 0, nhc
                                     sw(1-i,1-j,k,mult(swest)) = g%buf(i+1,j+1,kptr+k,ic)
                                   enddo
                                   enddo
                                   enddo
                                   endif
       endif
       mult(swest) = mult(swest)+1
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (mult(seast).eq.0) then
         if (reverse) then
           do k = 1, vlyr
             do j = 1, nhc
               do i = 1, nhc
                 v(np+i,1-j,k) = g%buf(j+1,i+1,kptr+k,ic)
               enddo
               enddo
               enddo
               else
                 do k = 1, vlyr
                   do j = 1, nhc
                     do i = 1, nhc
                       v(np+i,1-j,k) = g%buf(i+1,j+1,kptr+k,ic)
                     enddo
                     enddo
                     enddo
                     endif
                     else
                       if (reverse) then
                         do k = 1, vlyr
                           do j = 0, nhc
                             do i = 0, nhc
                               se(np+i,1-j,k,mult(seast)) = g%buf(j+1,i+1,kptr+k,ic)
                             enddo
                             enddo
                             enddo
                             else
                               do k = 1, vlyr
                                 do j = 0, nhc
                                   do i = 0, nhc
                                     se(np+i,1-j,k,mult(seast)) = g%buf(i+1,j+1,kptr+k,ic)
                                   enddo
                                   enddo
                                   enddo
                                   endif
       endif
       mult(seast) = mult(seast)+1
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (mult(neast).eq.0) then
         if (reverse) then
           do k = 1, vlyr
             do j = 1, nhc
               do i = 1, nhc
                 v(np+i,np+j,k) = g%buf(j+1,i+1,kptr+k,ic)
               enddo
               enddo
               enddo
               else
                 do k = 1, vlyr
                   do j = 1, nhc
                     do i = 1, nhc
                       v(np+i,np+j,k) = g%buf(i+1,j+1,kptr+k,ic)
                     enddo
                     enddo
                     enddo
                     endif
                     else
                       if (reverse) then
                         do k = 1, vlyr
                           do j = 0, nhc
                             do i = 0, nhc
                               ne(np+i,np+j,k,mult(neast)) = g%buf(j+1,i+1,kptr+k,ic)
                             enddo
                             enddo
                             enddo
                             else
                               do k = 1, vlyr
                                 do j = 0, nhc
                                   do i = 0, nhc
                                     ne(np+i,np+j,k,mult(neast)) = g%buf(i+1,j+1,kptr+k,ic)
                                   enddo
                                   enddo
                                   enddo
                                   endif
       endif
       mult(neast) = mult(neast)+1
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (mult(nwest).eq.0) then
         if (reverse) then
           do k = 1, vlyr
             do j = 1, nhc
               do i = 1, nhc
                 v(1-i,np+j,k) = g%buf(j+1,i+1,kptr+k,ic)
               enddo
               enddo
               enddo
               else
                 do k = 1, vlyr
                   do j = 1, nhc
                     do i = 1, nhc
                       v(1-i,np+j,k) = g%buf(i+1,j+1,kptr+k,ic)
                     enddo
                     enddo
                     enddo
                     endif
                     else
                       if (reverse) then
                         do k = 1, vlyr
                           do j = 0, nhc
                             do i = 0, nhc
                               nw(1-i,np+j,k,mult(nwest)) = g%buf(j+1,i+1,kptr+k,ic)
                             enddo
                             enddo
                             enddo
                             else
                               do k = 1, vlyr
                                 do j = 0, nhc
                                   do i = 0, nhc
                                     nw(1-i,np+j,k,mult(nwest)) = g%buf(i+1,j+1,kptr+k,ic)
                                   enddo
                                   enddo
                                   enddo
                                   endif
       endif
       mult(nwest) = mult(nwest)+1
     endif
   enddo
   end subroutine ghostvunpack3d
! =================================================================================
   subroutine ghostvpack2d_level(ghost, v, kptr, nhc, npoints, vlyr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem, ntrac_d
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(ghostbuffertr_t) :: ghost
   integer, intent(in   ) :: vlyr, kptr
   integer, intent(in   ) :: nhc, npoints
   real(real_kind), intent(in   ) :: v(1-nhc:npoints+nhc, 1-nhc:npoints+nhc, vlyr)
   type(edgedescriptor_t), intent(in   ) :: desc
! Local variables
   integer :: i, j, ir, l, e
   integer :: itr, k
   integer :: is, ie, in, iw
!
   if (.not.threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
     threadsafe = .true.
   endif
! Example convenction for buffer to the north:
!  buf(:,,:,i,e)
!   each "edge" is a row of data (i=1,np) in the element
!     north most row of data goes into e=1
!     next row of data goes into e=2
!     ....
!     south most row of data goes into e=np
!   We need to pack this way to preserve the orientation
!   so the data can be unpacked correctly
! note: we think of buf as dimensioned buf(k,is,i,e)
! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
!
!
   is = desc%putmapp_ghost(south)
   ie = desc%putmapp_ghost(east)
   in = desc%putmapp_ghost(north)
   iw = desc%putmapp_ghost(west)
   do k = 1,vlyr
     do i = 1, npoints
       do j = 1, nhc
         ghost%buf(i,j,k,kptr,is) = v(i,j+1,k)
         ghost%buf(i,j,k,kptr,ie) = v(npoints-j,i,k)
         ghost%buf(i,j,k,kptr,in) = v(i,npoints-j,k)
         ghost%buf(i,j,k,kptr,iw) = v(j+1,i,k)
       enddo
       enddo
   enddo
!  This is really kludgy way to setup the index reversals
!  But since it is so a rare event not real need to spend time optimizing
!  Check if the edge orientation of the recieving element is different
!  if it is,swap the order of data in the edge
   if (desc%reverse(south)) then
!        is = desc%putmapP_ghost(south)
     do k = 1, vlyr
       do i = 1, npoints
         do j = 1, nhc
           ir = npoints-i+1
           ghost%buf(ir,j,k,kptr,is) = v(i,j+1,k)
         enddo
         enddo
         enddo
   endif
   if (desc%reverse(east)) then
!        ie = desc%putmapP_ghost(east)
     do k = 1, vlyr
       do i = 1, npoints
         do j = 1, nhc
           ir = npoints-i+1
           ghost%buf(ir,j,k,kptr,ie) = v(npoints-j,i,k)
         enddo
         enddo
         enddo
   endif
   if (desc%reverse(north)) then
!        in = desc%putmapP_ghost(north)
     do k = 1, vlyr
       do i = 1, npoints
         do j = 1, nhc
           ir = npoints-i+1
           ghost%buf(ir,j,k,kptr,in) = v(i,npoints-j,k)
         enddo
         enddo
         enddo
   endif
   if (desc%reverse(west)) then
!        iw = desc%putmapP_ghost(west)
     do k = 1, vlyr
       do i = 1, npoints
         do j = 1, nhc
           ir = npoints-i+1
           ghost%buf(ir,j,k,kptr,iw) = v(j+1,i,k)
         enddo
         enddo
         enddo
   endif
! corners.  this is difficult because we dont know the orientaton
! of the corners,and this which (i,j) dimension maps to which dimension
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (l.ne.swest) call abortpar('ERROR3:swest ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
         do i = 1, nhc
           do j = 1, nhc
             ghost%buf(i,j,k,kptr,desc%putmapp_ghost(l)) = v(i+1,j+1,k)
           enddo
           enddo
           enddo
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (l.ne.seast) call abortpar('ERROR3:seast ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
         do i = 1, nhc
           do j = 1, nhc
             ghost%buf(i,j,k,kptr,desc%putmapp_ghost(l)) = v(npoints-i,j+1,k)
           enddo
           enddo
           enddo
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (l.ne.neast) call abortpar('ERROR3:neast ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
         do i = 1, nhc
           do j = 1, nhc
             ghost%buf(i,j,k,kptr,desc%putmapp_ghost(l)) = v(npoints-i,npoints-j,k)
           enddo
           enddo
           enddo
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (l.ne.nwest) call abortpar('ERROR3:nwest ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapp_ghost(l)/=-1) then
       do k = 1, vlyr
         do i = 1, nhc
           do j = 1, nhc
             ghost%buf(i,j,k,kptr,desc%putmapp_ghost(l)) = v(i+1,npoints-j,k)
           enddo
           enddo
           enddo
     endif
   enddo
   end subroutine ghostvpack2d_level
! =================================================================================
! GHOSTVUNPACK2D
! AUTHOR: Christoph Erath
! Unpack ghost points from ghost buffer into v...
! It is for cartesian points (v is only two dimensional).
! INPUT SAME arguments as for GHOSTVPACK2d
! =================================================================================
   subroutine ghostvunpack2d_level(ghost, v, kptr, nhc, npoints, vlyr, desc)
!-------------------------------------------------------------------------------
   use dimensions, only: max_corner_elem, ntrac_d
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   type(ghostbuffertr_t), intent(in   ) :: ghost
   integer, intent(in   ) :: nhc, npoints, vlyr, kptr
   real(real_kind), intent(inout) :: v(1-nhc:npoints+nhc, 1-nhc:npoints+nhc, vlyr)
   type(edgedescriptor_t) :: desc
! Local
   logical, parameter :: useunroll = .true.
   integer :: i, j, l, itr, k
   integer :: is, ie, in, iw, ic
   logical :: reverse
   real(real_kind) :: nan =-1.0
!
   nan = sqrt(nan)
   threadsafe = .false.
   is = desc%getmapp_ghost(south)
   ie = desc%getmapp_ghost(east)
   in = desc%getmapp_ghost(north)
   iw = desc%getmapp_ghost(west)
! example for north buffer
! first row ('edge') goes in v(:,np+1)
! 2nd   row ('edge') goes in v(:,np+2)
! etc...
   do k = 1,vlyr
     do i = 1, npoints
       do j = 1, nhc
         v(i,1-j,k) = ghost%buf(i,j,k,kptr,is)
         v(npoints+j,i,k) = ghost%buf(i,j,k,kptr,ie)
         v(i,npoints+j,k) = ghost%buf(i,j,k,kptr,in)
         v(1-j,i,k) = ghost%buf(i,j,k,kptr,iw)
       enddo
       enddo
   enddo
! SWEST
   do l = swest,swest+max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(1-i,1-j,k) = ghost%buf(j,i,k,kptr,ic)
             enddo
             enddo
             enddo
             else
               do k = 1, vlyr
                 do j = 1, nhc
                   do i = 1, nhc
                     v(1-i,1-j,k) = ghost%buf(i,j,k,kptr,ic)
                   enddo
                   enddo
                   enddo
       endif
     else
       do k = 1, vlyr
         do j = 1, nhc
           do i = 1, nhc
             v(1-i,1-j,k) = nan
           enddo
           enddo
           enddo
     endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(npoints+i,1-j,k) = ghost%buf(j,i,k,kptr,ic)
             enddo
             enddo
             enddo
             else
               do k = 1, vlyr
                 do j = 1, nhc
                   do i = 1, nhc
                     v(npoints+i,1-j,k) = ghost%buf(i,j,k,kptr,ic)
                   enddo
                   enddo
                   enddo
       endif
     else
       do k = 1, vlyr
         do j = 1, nhc
           do i = 1, nhc
             v(npoints+i,1-j,k) = nan
           enddo
           enddo
           enddo
     endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(npoints+i,npoints+j,k) = ghost%buf(j,i,k,kptr,ic)
             enddo
             enddo
             enddo
             else
               do k = 1, vlyr
                 do j = 1, nhc
                   do i = 1, nhc
                     v(npoints+i,npoints+j,k) = ghost%buf(i,j,k,kptr,ic)
                   enddo
                   enddo
                   enddo
       endif
     else
       do k = 1, vlyr
         do j = 1, nhc
           do i = 1, nhc
             v(npoints+i,npoints+j,k) = nan
           enddo
           enddo
           enddo
     endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     ic = desc%getmapp_ghost(l)
     if (ic/=-1) then
       reverse = desc%reverse(l)
       if (reverse) then
         do k = 1, vlyr
           do j = 1, nhc
             do i = 1, nhc
               v(1-i,npoints+j,k) = ghost%buf(j,i,k,kptr,ic)
             enddo
             enddo
             enddo
             else
               do k = 1, vlyr
                 do j = 1, nhc
                   do i = 1, nhc
                     v(1-i,npoints+j,k) = ghost%buf(i,j,k,kptr,ic)
                   enddo
                   enddo
                   enddo
       endif
     else
       do k = 1, vlyr
         do j = 1, nhc
           do i = 1, nhc
             v(1-i,npoints+j,k) = nan
           enddo
           enddo
           enddo
     endif
   enddo
   end subroutine ghostvunpack2d_level
   subroutine edgevunpackcheck(iee, edge, v, vlyr, kptr, desc, ncount, ncount_t, isprt)
!-------------------------------------------------------------------------------
   use dimensions, only: np, max_corner_elem
   use control, only: north, south, east, west, neast, nwest, seast, swest
!
   integer, intent(in   ) :: iee
   type(edgebuffer_t), intent(in   ) :: edge
   integer, intent(in   ) :: vlyr
   real(real_kind), intent(inout) :: v(np, np, vlyr)
   integer, intent(in   ) :: kptr
   type(edgedescriptor_t), intent(in   ) :: desc
   integer, intent(  out) :: ncount, ncount_t
   logical, intent(in   ) :: isprt
! Local
   integer :: i, k, l
   integer :: is, ie, in, iw, l_ncount, l_ncount_t
!
   l_ncount = 0
   l_ncount_t = 0
   threadsafe = .false.
   is = desc%getmapp(south)
   ie = desc%getmapp(east)
   in = desc%getmapp(north)
   iw = desc%getmapp(west)
   do k = 1,vlyr
     do i = 1, np
       call checkedge(iee,i,1,k,south,v(i,1,k),edge%buf(kptr+k,is+i),l_ncount,l_ncount_t,isprt)
       call checkedge(iee,np,i,k,east,v(np,i,k),edge%buf(kptr+k,ie+i),l_ncount,l_ncount_t,isprt)
       call checkedge(iee,i,np,k,north,v(i,np,k),edge%buf(kptr+k,in+i),l_ncount,l_ncount_t,isprt)
       call checkedge(iee,1,i,k,west,v(1,i,k),edge%buf(kptr+k,iw+i),l_ncount,l_ncount_t,isprt)
     enddo
   enddo
! SWEST
   do l = swest,swest+max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
!v(1  ,1 ,k)=MAX(v(1 ,1 ,k),edge%buf(kptr+k,desc%getmapP(l)+1))
         call checkedge(iee,1,1,k,swest,v(1,1,k),edge%buf(kptr+k,desc%getmapp(l)+1),l_ncount,l_ncount_t,isprt)
       enddo
       endif
   enddo
! SEAST
   do l = swest+max_corner_elem,swest+2*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
!v(np ,1 ,k)=MAX(v(np,1 ,k),edge%buf(kptr+k,desc%getmapP(l)+1))
         call checkedge(iee,np,1,k,seast,v(np,1,k),edge%buf(kptr+k,desc%getmapp(l)+1),l_ncount,l_ncount_t,isprt)
       enddo
       endif
   enddo
! NEAST
   do l = swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
!v(np ,np,k)=MAX(v(np,np,k),edge%buf(kptr+k,desc%getmapP(l)+1))
         call checkedge(iee,np,np,k,neast,v(np,np,k),edge%buf(kptr+k,desc%getmapp(l)+1),l_ncount,l_ncount_t,isprt)
       enddo
       endif
   enddo
! NWEST
   do l = swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (desc%getmapp(l)/=-1) then
       do k = 1, vlyr
!v(1  ,np,k)=MAX(v(1 ,np,k),edge%buf(kptr+k,desc%getmapP(l)+1))
         call checkedge(iee,1,np,k,nwest,v(1,np,k),edge%buf(kptr+k,desc%getmapp(l)+1),l_ncount,l_ncount_t,isprt)
       enddo
       endif
   enddo
   ncount = l_ncount
   ncount_t = l_ncount_t
   end subroutine edgevunpackcheck
   subroutine checkedge(ie, ip, jp, kp, dir, org, tgt, ncount, ncount_t, isprt)
!-------------------------------------------------------------------------------
   use kiapsbase, only: i4=>kim_int_kind, r8=>kim_real8_kind
   implicit none
!
   integer(i4), intent(in   ) :: ie, ip, jp, kp, dir
   real(r8), intent(in   ) :: org, tgt
   integer, intent(inout) :: ncount
   integer, intent(inout) :: ncount_t
   logical, intent(in   ) :: isprt
!
   ncount_t = ncount_t+1
   if (org.ne.tgt) then
!      PRINT '(A)', 'Edge diff: ie, i, j, k, direction, org, tgt:'
     ncount = ncount+1
     if ((ie==31.or.ie==32.or.ie==61.or.ie==62).and.kp==1) then
       if (isprt) print '(i4,x,i1,x,i1,x,i2,x,i2,x,e24.18,x,e24.18,x,e24.18) ',ie,ip,jp,kp,dir,org,tgt,(tgt-org)/org*100.0_r8
     endif
   endif
   end subroutine checkedge
   subroutine scs(add_a, err, add_b)
!-------------------------------------------------------------------------------
   use kiapsbase, only: r8=>kim_real8_kind
   implicit none
!
   real(r8), intent(inout) :: add_a, err
   real(r8), intent(in   ) :: add_b
!local
   real(r8) :: tsum, new_err, tmp
!
   tmp = err+add_b
   tsum = add_a+tmp
   new_err = tmp-(tsum-add_a)
   add_a = tsum
   err = new_err
   end subroutine scs
   subroutine ddpddlocal(add_a, err, add_b)
!-------------------------------------------------------------------------------
   use kiapsbase, only: r8=>kim_real8_kind
   implicit none
!
   real(r8), intent(inout) :: add_a, err
   real(r8), intent(in   ) :: add_b
!local
   real(r8) :: tsum, new_err
!
   tsum = add_a+add_b
!    new_err = (add_b-(tsum-add_a))+(add_a-(tsum-(tsum-add_a)))+err
!    CALL SCS(tsum,new_err,0.0_r8)
   new_err =(add_b-(tsum-add_a))+(add_a-(tsum-(tsum-add_a)))
   call scs(tsum,err,new_err)
   add_a = tsum
!    err     = new_err
   end subroutine ddpddlocal
!
   end module edge
!-------------------------------------------------------------------------------
#ifndef HAVE_F2003_PTR_BND_REMAP
!
! subroutine to allow sharing edge buffers
! this has to be outside a module to allow us to (F77 style) access the same chunk
! of memory with a different shape
!
! some compilers dont allow the 'target' attribute to be used in a F77 subroutine
! such as cray.  if that is the case, try compiling with -DHAVE_F2003_PTR_BND_REMAP
!
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine remap_2d_ptr_buf(edge, nlyr, nbuf, src_array)
!-------------------------------------------------------------------------------
   use kiapsbase, only: real_kind=>kim_real8_kind
   use edge, only: edgebuffer_t ! _external
! input
!
   type(edgebuffer_t) :: edge
   integer :: nlyr, nbuf
   real(real_kind), target :: src_array(nlyr, nbuf)
!
   edge%buf=>src_array
!
   end subroutine remap_2d_ptr_buf
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine remap_2d_ptr_receive(edge, nlyr, nbuf, src_array)
!-------------------------------------------------------------------------------
   use kiapsbase, only: real_kind=>kim_real8_kind
   use edge, only: edgebuffer_t ! _external
! input
!
   type(edgebuffer_t) :: edge
   integer :: nlyr, nbuf
   real(real_kind), target :: src_array(nlyr, nbuf)
!
   edge%receive=>src_array
!
   end subroutine remap_2d_ptr_receive
#endif
