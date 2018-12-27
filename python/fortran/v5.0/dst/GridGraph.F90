!-------------------------------------------------------------------------------
!>
!> @brief
!>  - grid graph module.
!>
!> @date ?????2012
!>  - JH KIM  : First written from the HOMME and was modified for KIAPSGM framework.
!> @date 30JAN2015
!>  - JH KIM  : Added to "nbrs_face", "nbrs_ptr" variables and some subroutines. (CAM-SE 5.3)
!-------------------------------------------------------------------------------
#include "KIM.h"
!-------------------------------------------------------------------------------
   module gridgraph
!-------------------------
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
   use kiapsbase, only: iulog=>kim_iu_log, real_kind=>kim_real8_kind
!-------------------------
   use control, only: north, south, east, west, neast, nwest, seast, swest
!-----
   implicit none
!
   private
   integer, public, parameter :: num_neighbors = 8
!
   type, public :: ptr
     sequence
     integer, dimension(:), pointer :: n=>null() ! number of neighbor element
     integer, dimension(:), pointer :: f=>null() ! cube face number of neighbor element
   end type ptr
!
   type, public :: gridvertex_t
     sequence
     type(ptr) :: nbrs(num_neighbors)
     integer :: wgtp(num_neighbors) ! the weights for edges defined by neighbors array
     integer :: wgtp_ghost(num_neighbors) ! the weights for edges defined by neighbors array
     integer, pointer :: nbrs_face(:)=>null() ! the cube face number of the neighbor element(nbrs array)
     integer :: nbrs_ptr(num_neighbors+1) ! for vertically lagrangian
     integer :: face_number ! which face of the cube this vertex is on
     integer :: number ! element number
     integer :: processor_number ! processor number
     integer :: spacecurve ! index in space-filling curve
   end type gridvertex_t
!
   type, public :: edgeindex_t
     sequence
     integer, pointer :: ixp(:)=>null()
     integer, pointer :: iyp(:)=>null()
   end type edgeindex_t
!
   type, public :: gridedge_t
     sequence
     integer :: head_face ! needed if head vertex has shape(i.e. square)
     integer :: tail_face ! needed if tail vertex has shape(i.e. square)
     integer :: head_ind !
     integer :: tail_ind !
     type(gridvertex_t), pointer :: head=>null() ! edge head vertex
     type(gridvertex_t), pointer :: tail=>null() ! edge tail vertex
     logical :: reverse
   end type gridedge_t
! ==========================================
! Public Interfaces
! ==========================================
   public :: set_gridvertex_number
   public :: printgridvertex
   public :: initgridedge
   public :: gridedge_search
   public :: gridedge_type
   public :: grid_edge_uses_vertex
   public :: printgridedge
   public :: checkgridneighbors
   public :: printchecksum
   public :: createsubgridgraph
   public :: freegraph
!
   interface assignment(=)
     module procedure copy_gridedge
     module procedure copy_edgeindex
     module procedure copy_gridvertex
   end interface
!
   contains
! =====================================
! copy edge:
! copy device for overloading = sign.
! =====================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   recursive subroutine copy_gridedge(edge2, edge1)
!
   type(gridedge_t), intent(  out) :: edge2
   type(gridedge_t), intent(in   ) :: edge1
!
   edge2%tail_face = edge1%tail_face
   edge2%head_face = edge1%head_face
   edge2%reverse = edge1%reverse
   if (associated(edge1%tail)) then
     edge2%tail=>edge1%tail
   endif
   if (associated(edge1%head)) then
     edge2%head=>edge1%head
   endif
#ifdef TESTGRID
   edge2%wgtv = edge1%wgtv
   edge2%wgtp = edge1%wgtp
   edge2%tailindex = edge1%tailindex
   edge2%headindex = edge1%headindex
#endif
!
   end subroutine copy_gridedge
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   recursive subroutine copy_gridvertex(vertex2, vertex1)
!-------------------------------------------------------------------------------
   implicit none
!
   type(gridvertex_t), intent(  out) :: vertex2
   type(gridvertex_t), intent(in   ) :: vertex1
   integer :: i, j, n
!JMD     vertex2%size      = vertex1%size
!
   n = size(vertex1%nbrs)
   do i = 1,n
     vertex2%wgtp(i) = vertex1%wgtp(i)
     vertex2%wgtp_ghost(i) = vertex1%wgtp_ghost(i)
     if (associated(vertex2%nbrs(i)%n)) &
                                           deallocate(vertex2%nbrs(i)%n)
     if (associated(vertex2%nbrs(i)%f)) &
                                           deallocate(vertex2%nbrs(i)%f)
     allocate(vertex2%nbrs(i)%n(size(vertex1%nbrs(i)%n)))
     allocate(vertex2%nbrs(i)%f(size(vertex1%nbrs(i)%n)))
     do j = 1,size(vertex1%nbrs(i)%n)
       vertex2%nbrs(i)%n(j) = vertex1%nbrs(i)%n(j)
       vertex2%nbrs(i)%f(j) = vertex1%nbrs(i)%f(j)
     enddo
   enddo
   do i = 1,num_neighbors+1
     vertex2%nbrs_ptr(i) = vertex1%nbrs_ptr(i)
   enddo
   vertex2%number = vertex1%number
   vertex2%processor_number = vertex1%processor_number
   vertex2%spacecurve = vertex1%spacecurve
!
   end subroutine copy_gridvertex
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   recursive subroutine copy_edgeindex(index2, index1)
!
   type(edgeindex_t), intent(  out) :: index2
   type(edgeindex_t), intent(in   ) :: index1
!
   if (associated(index1%iyp)) then
     index2%iyp=>index1%iyp
   endif
!
   if (associated(index1%ixp)) then
     index2%ixp=>index1%ixp
   endif
!
   end subroutine copy_edgeindex
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine freegraph(vertex)
!-------------------------------------------------------------------------------
   implicit none
!
   type(gridvertex_t) :: vertex(:)
   integer :: i, nelem
!
   nelem = size(vertex)
!JMD     do i=1,nelem
!JMD        deallocate(Vertex(i)%wgtV)
!JMD        deallocate(Vertex(i)%wgtG)
!JMD        deallocate(Vertex(i)%nbrs)
!JMD     enddo
!
   end subroutine freegraph
!===========================
! search edge list for match
!===========================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function gridedge_search(nvert1, nvert2, edge) result(number)
!
   integer, intent(in   ) :: nvert1
   integer, intent(in   ) :: nvert2
   type(gridedge_t), intent(in   ) :: edge(:)
   integer :: number
   integer :: tmp
   integer :: head
   integer :: tail
   integer :: nedge
   integer :: i
!
   nedge = size(edge)
   tail = nvert1
   head = nvert2
   if (tail>head) then
     tmp = tail
     tail = head
     head = tmp
   endif
   do i = 1,nedge
     if (edge(i)%tail%number==tail.and.edge(i)%head%number==head) then
       number = i
     endif
   enddo
!
   end function gridedge_search
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function gridedge_type(edge) result(type)
!-------------------------------------------------------------------------------
   use params, only: internal_edge, external_edge
!
   type(gridedge_t), intent(in   ) :: edge
   integer :: type
!
   if (edge%head%processor_number==edge%tail%processor_number) then
     type = internal_edge
     else
       type = external_edge
       endif
       end function gridedge_type
     function grid_edge_uses_vertex(vertex, edge) result(log)
!
   type(gridvertex_t), intent(in   ) :: vertex
   type(gridedge_t), intent(in   ) :: edge
   logical :: log
   integer :: number
!
   number = vertex%number
   if (number==edge%head%number.or.number==edge%tail%number) then
   log = .true.
   else
   log = .false.
   endif
     end function grid_edge_uses_vertex
   subroutine printchecksum(testpattern, checksum)
!-------------------------------------------------------------------------------
   use dimensions, only: nlev, nelemd, np
   implicit none
!
   real(real_kind), target, intent(in   ) :: testpattern(:,:,:,:)
   real(real_kind), target, intent(in   ) :: checksum(:,:,:,:)
   integer :: i, k, ix, iy
!
   print*
   write(iulog,*) 'checksums:'
   do i = 1,nelemd
!  Lets start out only looking at the first element
     write(iulog,*)
     do k = 1,nlev
       do iy = 1, np
         do ix = 1, np
           write(iulog,*) int(testpattern(ix,iy,k,i))," checksum = ",int(checksum(ix,iy,k,i))
         enddo
         enddo
     enddo
   enddo
   end subroutine printchecksum
   subroutine createsubgridgraph(vertex, svertex, local2global)
!-------------------------------------------------------------------------------
   implicit none
!
   type(gridvertex_t), intent(in   ) :: vertex(:)
   type(gridvertex_t), intent(inout) :: svertex(:)
   integer, intent(in   ) :: local2global(:)
   integer :: nelem, nelem_s, n, ncount
   integer :: inbr, i, ig, j, k
   integer, allocatable :: global2local(:)
   logical, parameter :: debug = .false.
!
   nelem = size(vertex)
   nelem_s = size(svertex)
   if (debug) write(iulog,*) 'CreateSubGridGraph:point #1'
   allocate(global2local(nelem))
   if (debug) write(iulog,*) 'CreateSubGridGraph:point #2'
   global2local(:) = 0
   do i = 1,nelem_s
     ig = local2global(i)
     global2local(ig) = i
   enddo
   if (debug) write(iulog,*) 'CreateSubGridGraph:point #3'
   do i = 1,nelem_s
     ig = local2global(i)
     if (debug) write(iulog,*) 'CreateSubGridGraph:point #4'
     call copy_gridvertex(svertex(i),vertex(ig))
     n = size(svertex(i)%nbrs(:))
! ==============================================
! Apply the correction to the neighbors list to
! reflect new subgraph numbers
! ==============================================
     ncount = 0
     if (debug) write(iulog,*) 'CreateSubGridGraph:point #5'
     do j = 1,n
       if (associated(svertex(i)%nbrs(j)%n)) then
         do k = 1, size(svertex(i)%nbrs(j)%n)
           if (debug) write(iulog,*) 'CreateSubGridGraph:point #5.1 size(global2local) svertex(i)%nbrs(j) ',&
                     size(global2local),svertex(i)%nbrs(j)%n(k)
           inbr = global2local(svertex(i)%nbrs(j)%n(k))
           if (debug) write(iulog,*) 'CreateSubGridGraph:point #5.2'
           if (inbr.gt.0) then
             svertex(i)%nbrs(j)%n(k) = inbr
             ncount = ncount+1
           else
             svertex(i)%wgtp(j) = 0
             svertex(i)%wgtp_ghost(j) = 0
             deallocate(svertex(i)%nbrs(j)%n)
             deallocate(svertex(i)%nbrs(j)%f)
           endif
         enddo
         endif
           if (debug) write(iulog,*) 'CreateSubGridGraph:point #5.3'
     enddo
     if (debug) write(iulog,*) 'CreateSubGridGraph:point #6'
     svertex(i)%number = i
   enddo
   if (debug) write(iulog,*) 'CreateSubGridGraph:point #7'
   deallocate(global2local)
   if (debug) write(iulog,*) 'CreateSubGridGraph:point #8'
   end subroutine createsubgridgraph
   subroutine printgridedge(edge)
!-------------------------------------------------------------------------------
   implicit none
!
   type(gridedge_t), intent(in   ) :: edge(:)
   integer :: i, nedge, ii, wgtp
!
   nedge = size(edge)
   write(iulog,95)
   do i = 1,nedge
     ii = edge(i)%tail_face
     wgtp = edge(i)%tail%wgtp(ii)
     write(iulog,100) i,&
                       edge(i)%tail%number,edge(i)%tail_face,wgtp,&
 edge(i)%head%number,edge(i)%head_face,gridedge_type(edge(i))
   enddo
   95 format(5x,'GRIDEDGE #',3x,'Tail(face) ',5x,'Head(face) ',3x,'Type')
   100 format(10x,i6,8x,i4,1x,'(',i1,')--',i2,'-->',i6,1x,'(',i1,') ',5x,'[',i1,']')
   end subroutine printgridedge
! ==========================================
! set_GridVertex_neighbors:
!
! Set global element number for element elem
! ==========================================
   subroutine set_gridvertex_number(elem, number)
!
   type(gridvertex_t) :: elem
   integer :: number
!
   elem%number = number
   end subroutine set_gridvertex_number
   subroutine printgridvertex(vertex)
!-------------------------------------------------------------------------------
   implicit none
!
   type(gridvertex_t), intent(in   ), target :: vertex(:)
   integer :: i, nvert
!
   nvert = size(vertex)
   write(iulog,98)
   do i = 1,nvert
     write(iulog,99) vertex(i)%number,vertex(i)%processor_number,&
           vertex(i)%nbrs(west)%n(:),vertex(i)%wgtp(west),&
           vertex(i)%nbrs(east)%n(:),vertex(i)%wgtp(east),&
         vertex(i)%nbrs(south)%n(:),vertex(i)%wgtp(south),&
         vertex(i)%nbrs(north)%n(:),vertex(i)%wgtp(north),&
         vertex(i)%nbrs(swest)%n(:),vertex(i)%wgtp(swest),&
         vertex(i)%nbrs(seast)%n(:),vertex(i)%wgtp(seast),&
         vertex(i)%nbrs(nwest)%n(:),vertex(i)%wgtp(nwest),&
           vertex(i)%nbrs(neast)%n(:),vertex(i)%wgtp(neast)
   enddo
   98 format(5x,'GRIDVERTEX #',2x,'PART',2x,'DEG',4x,'W',8x,'E',8x,&
               'S',8x,'N',7x,'SW',7x,'SE',7x,'NW',7x,'NE')
   99 format(10x,i3,8x,i1,2x,i1,2x,8(2x,i3,1x,'(',i2,') '))
   end subroutine printgridvertex
   subroutine checkgridneighbors(vertex)
!-------------------------------------------------------------------------------
   implicit none
!
   type(gridvertex_t), intent(in   ) :: vertex(:)
   integer :: i, j, k, l, m, nnbrs, inbrs, nvert
!
   nvert = size(vertex)
   do i = 1,nvert
     nnbrs = size(vertex(i)%nbrs)
     do j = 1,nnbrs
       do l = 1, size(vertex(i)%nbrs(j)%n)
         inbrs = vertex(i)%nbrs(j)%n(l)
         if (inbrs>0) then
           do k = 1, nnbrs
             do m = 1, size(vertex(i)%nbrs(k)%n)
               if (inbrs.eq.vertex(i)%nbrs(k)%n(m).and.(j/= k.or.l/= m)) &
write(iulog,*) 'CheckGridNeighbors:error identical neighbors detected for vertex ',i
             enddo
             enddo
         endif
       enddo
     enddo
   enddo
   end subroutine checkgridneighbors
   subroutine initgridedge(gridedge, gridvertex)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   implicit none
!
   type(gridedge_t), intent(inout) :: gridedge(:)
   type(gridvertex_t), intent(in   ), target :: gridvertex(:)
   integer :: i, j, k, m, l, iptr, wgtv, wgtp
   integer :: nelem, nelem_edge, inbr
   logical :: verbose = .false.
!
   nelem = size(gridvertex)
   nelem_edge = size(gridedge)
   gridedge(:)%reverse = .false.
   iptr = 1
   do j = 1,nelem
     do i = 1, num_neighbors
       if ((gridvertex(j)%wgtp(i).gt.0).and.(associated(gridvertex(j)%nbrs(i)%n))) then ! do this only if has a non-zero weight
         do l = 1, size(gridvertex(j)%nbrs(i)%n)
           if (nelem_edge<iptr) call abortpar(message = 'Error in initgridedge:number of edges greater than expected.')
           gridedge(iptr)%tail=>gridvertex(j)
           gridedge(iptr)%tail_face = i
           gridedge(iptr)%tail_ind = l
           inbr = gridvertex(j)%nbrs(i)%n(l)
           gridedge(iptr)%head=>gridvertex(inbr)
#ifdef TESTGRID
           wgtv = gridvertex(j)%wgtv(i)
           gridedge(iptr)%wgtv = wgtv
           wgtp = gridvertex(j)%wgtp(i)
           gridedge(iptr)%wgtp = wgtp
! allocate the indirect addressing arrays
           allocate(gridedge(iptr)%tailindex%ixv(wgtv))
           allocate(gridedge(iptr)%tailindex%ixp(wgtp))
           allocate(gridedge(iptr)%tailindex%iyv(wgtv))
           allocate(gridedge(iptr)%tailindex%iyp(wgtp))
           allocate(gridedge(iptr)%headindex%ixv(wgtv))
           allocate(gridedge(iptr)%headindex%ixp(wgtp))
           allocate(gridedge(iptr)%headindex%iyv(wgtv))
           allocate(gridedge(iptr)%headindex%iyp(wgtp))
#endif
! ===========================================
! Need this aweful piece of code to determine
! which "face" of the neighbor element the
! edge links (i.e. the "head_face")
! ===========================================
           do k = 1,num_neighbors
             if (associated(gridvertex(inbr)%nbrs(k)%n)) then
               do m = 1, size(gridvertex(inbr)%nbrs(k)%n)
                 if (gridvertex(inbr)%nbrs(k)%n(m)==gridvertex(j)%number) then
                   gridedge(iptr)%head_face = k
                   gridedge(iptr)%head_ind = m
                 endif
                 enddo
                 endif
           enddo
           iptr = iptr+1
         enddo
         endif
         enddo
   enddo
   if (nelem_edge+1/= iptr) call abortpar(message = 'Error in initgridedge:number of edges less than expected.')
   if (verbose) then
     print*
     write(iulog,*) "element edge tail,head list :(test) "
     do i = 1,nelem_edge
       write(iulog,*) gridedge(i)%tail%number,gridedge(i)%head%number
     enddo
     print*
     write(iulog,*) "element edge tail_face,head_face list :(test) "
     do i = 1,nelem_edge
       write(iulog,*) gridedge(i)%tail_face,gridedge(i)%head_face
     enddo
   endif
   end subroutine initgridedge
   end module gridgraph
!-------------------------------------------------------------------------------
