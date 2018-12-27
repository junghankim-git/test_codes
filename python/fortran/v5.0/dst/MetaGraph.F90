#include "KIM.h"
!************************MetaGraph.F****************************************
!-------------------------------------------------------------------------------
   module metagraph
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
   use kiapsbase, only: int_kind=>kim_int_kind, iulog=>kim_iu_log
   use gridgraph, only: gridvertex_t, gridedge_t
   use pio_types ! _external
   implicit none
!
   private
!
   type, public :: metaedge_t
     sequence
     type(gridedge_t), pointer :: members(:)
     integer, pointer :: edgeptrp(:)
     integer, pointer :: edgeptrp_ghost(:)
     integer :: number
     integer :: type
     integer :: wgtp ! sum of lengths of all messages to pack for edges
     integer :: wgtp_ghost ! sum of lengths of all messages to pack for ghost cells
     integer :: headvertex ! processor number to send to
     integer :: tailvertex ! processor number to send from
     integer :: nmembers ! number of messages to(un) pack(out) into this buffer
     integer :: padding ! just to quite compiler
   end type metaedge_t
!
   type, public :: metavertex_t ! one for each processor
     sequence
     integer :: number ! useless just the local processor number
     integer :: nmembers ! number of elements on this processor
     type(gridvertex_t), pointer :: members(:) ! array of elements on this processor
     type(metaedge_t), pointer :: edges(:) ! description of messages to send/receive
     integer :: nedges ! number of processors to communicate with(length of edges)
     integer :: padding ! just to quite compiler
   end type metavertex_t
! public :: findedge
   public :: edge_uses_vertex
   public :: printmetaedge, printmetavertex
   public :: localelemcount
!public :: MetaEdgeCount
   public :: initmetagraph
!
   interface assignment(=)
     module procedure copy_metaedge
   end interface
!
   contains
! =====================================
! copy vertex:
! copy device for overloading = sign.
! =====================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   recursive subroutine copy_metaedge(edge2, edge1)
!
   type(metaedge_t), intent(  out) :: edge2
   type(metaedge_t), intent(in   ) :: edge1
   integer i
!
   edge2%number = edge1%number
   edge2%type = edge1%type
   edge2%wgtp = edge1%wgtp
   edge2%wgtp_ghost = edge1%wgtp_ghost
   edge2%nmembers = edge1%nmembers
   if (associated(edge1%members)) then
     allocate(edge2%members(edge2%nmembers))
     do i = 1,edge2%nmembers
       edge2%members(i) = edge1%members(i)
     enddo
   endif
   if (associated(edge1%edgeptrp)) then
     allocate(edge2%edgeptrp(edge2%nmembers))
     allocate(edge2%edgeptrp_ghost(edge2%nmembers))
     do i = 1,edge2%nmembers
       edge2%edgeptrp(i) = edge1%edgeptrp(i)
       edge2%edgeptrp_ghost(i) = edge1%edgeptrp_ghost(i)
     enddo
   endif
   edge2%headvertex = edge1%headvertex
   edge2%tailvertex = edge1%tailvertex
!
   end subroutine copy_metaedge
! function findedge(mEdge,Edge) result(number)
!   type(MetaEdge_t), intent(inout) :: mEdge(:)
!   type(GridEdge_t), intent(in)    :: Edge
!   integer :: number
!   integer :: head,tail,exist
!   integer :: nedge
!   integer :: i
!   nedge=SIZE(mEdge)
!   number = 0
!   tail=Edge%tail%processor_number
!   head=Edge%head%processor_number
!   exist=0
!   do i=1,nedge
!      !       write(iulog,*)'mEdge(i)%number: ',mEdge(i)%number
!      if(mEdge(i)%number .ne. 0) then
!         if  ((mEdge(i)%TailVertex==tail .and. mEdge(i)%HeadVertex==head) ) then
!            number=i
!         end if
!         exist=exist+1
!      endif
!   end do
!   if(number == 0) number = exist + 1
! end function findedge
! function MetaEdgeCount(Edge) result(nMedge)
!   implicit none
!   type (GridEdge_t),intent(in)  :: Edge(:)
!   integer                       :: nMedge
!   integer                       :: nedge,i,j,maxedges
!   integer                       :: head_processor_number,tail_processor_number
!   integer, allocatable          :: tmp(:,:)
!   logical                       :: found
!   nedge = SIZE(Edge)
!   maxedges=nedge
!   allocate(tmp(2,maxedges))
!   tmp=0
!   nMedge=0
!   do i=1,nedge
!      head_processor_number = Edge(i)%head%processor_number
!      tail_processor_number = Edge(i)%tail%processor_number
!      found = .FALSE.
!      do j=1,nMedge
!         if((tmp(1,j) .eq. head_processor_number).and. !              (tmp(2,j) .eq. tail_processor_number)) found=.TRUE.
!      enddo
!      if(.NOT. found) then
!         nMedge=nMedge+1
!         tmp(1,nMedge) = head_processor_number
!         tmp(2,nMedge) = tail_processor_number
!      endif
!   enddo
!   !mem    write(iulog,*)'MetaEdgeCount: before call to deallocate(tmp)'
!   deallocate(tmp)
! end function MetaEdgeCount
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function localelemcount(vertex) result(nelemd)
!-------------------------------------------------------------------------------
   implicit none
!
   type(metavertex_t), intent(in   ) :: vertex
   integer :: nelemd
!
   nelemd = vertex%nmembers
!
   end function localelemcount
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function edge_uses_vertex(vertex, edge) result(log)
!-------------------------------------------------------------------------------
   implicit none
!
   type(metavertex_t), intent(in   ) :: vertex
   type(metaedge_t), intent(in   ) :: edge
   logical :: log
   integer :: number
!
   number = vertex%number
   if (number==edge%headvertex.or.number==edge%tailvertex) then
   log = .true.
   else
   log = .false.
   endif
!
   end function edge_uses_vertex
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine printmetaedge(edge)
!-------------------------------------------------------------------------------
   use gridgraph, only: printgridedge
   implicit none
!
   type(metaedge_t), intent(in   ) :: edge(:)
   integer :: i, nedge
!
   nedge = size(edge)
   do i = 1,nedge
     print*
     write(iulog,90) edge(i)%number,edge(i)%type,edge(i)%wgtp,edge(i)%nmembers,edge(i)%tailvertex,edge(i)%headvertex
     if (associated(edge(i)%members)) then
       call printgridedge(edge(i)%members)
     endif
   enddo
   90 format('METAEDGE #',i4,2x,'TYPE ',i1,2x,'WGT ',i4,2x,'NUM ',i6,2x,'Processors ',i4,'--->',i4)
!
   end subroutine printmetaedge
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine printmetavertex(vertex)
!-------------------------------------------------------------------------------
   use gridgraph, only: printgridvertex
   implicit none
!
   type(metavertex_t), intent(in   ), target :: vertex
   integer :: i, j, npart
!
   write(iulog,*)
   write(iulog,90) i
   write(iulog,95) vertex%nmembers
   call printgridvertex(vertex%members)
   write(iulog,96) vertex%nedges
   if (associated(vertex%edges)) then
     do j = 1, vertex%nedges
       write(iulog,97) vertex%edges(j)%number,vertex%edges(j)%type,vertex%edges(j)%wgtp,vertex%edges(j)%headvertex,vertex%edges(j)%tailvertex
     enddo
   endif
   90 format('METAVERTEX #',i2,2x)
   95 format(5x,i2,' member grid vertices')
   96 format(5x,i2,' incident meta edges ')
   97 format(10x,'METAEDGE #',i2,2x,'TYPE ',i1,2x,'WGT ',i4,2x,'Processors ',i2,'--->',i2)
   98 format(5x,'GRIDVERTEX #',5x,'West',5x,'East',5x,'South',5x,'North')
   99 format(10x,i2,7x,4(2x,i3,1x,'(',i1,') '))
!
   end subroutine printmetavertex
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine initmetagraph(thisprocessornumber, metavertex, gridvertex, gridedge)
!-------------------------------------------------------------------------------
   use ll_mod, only: root_t, llsetedgecount, llfree, llinsertedge, llgetedgecount, llfindedge
   use gridgraph, only: gridedge_type, printgridvertex
!------------------
!------------------
   implicit none
!
   integer(int_kind), intent(in   ) :: thisprocessornumber
   type(metavertex_t), intent(  out) :: metavertex
   type(gridvertex_t), intent(in   ), target :: gridvertex(:)
   type(gridedge_t), intent(in   ), target :: gridedge(:)
   type(metaedge_t), allocatable :: metaedge(:)
   integer :: nelem, nelem_edge, nedges
   integer, allocatable :: icount(:)
   integer :: ic, i, j, ii
   integer :: npart
   integer :: head_processor_number
   integer :: tail_processor_number
   integer :: nedge_active, enum
   logical :: found
   integer itail, ihead, wgtp
   type(root_t) :: medgelist ! root_t = c++std :: set<std :: pair<int, int>>
   logical :: verbose = .false.
   logical :: debug = .false.
!
   if (debug) write(iulog,*) 'initMetagraph:point #1'
!  Number of grid vertices
   nelem = size(gridvertex)
!  Number of grid edges
   nelem_edge = size(gridedge)
   medgelist%number = thisprocessornumber
   nullify(medgelist%first)
   call llsetedgecount(0)
   do i = 1,nelem_edge
     tail_processor_number = gridedge(i)%tail%processor_number
     head_processor_number = gridedge(i)%head%processor_number
     if (tail_processor_number.eq.thisprocessornumber.or.head_processor_number.eq.thisprocessornumber) then
       call llinsertedge(medgelist,tail_processor_number,head_processor_number,enum)
     endif
   enddo
   call llgetedgecount(nedges)
   allocate(metaedge(nedges))
! Initalize the Meta Vertices to zero... probably should be done
! in a separate routine
   metavertex%nmembers = 0
   metavertex%number = 0
   metavertex%nedges = 0
   if (debug) write(iulog,*) 'initMetagraph:point #2'
   nullify(metavertex%edges)
!  Give some identity to the Meta_vertex
   metavertex%number = thisprocessornumber
   if (debug) write(iulog,*) 'initMetagraph:point #3'
!  Look through all the small_vertices and determine the number of
!  member vertices
   if (debug) call printgridvertex(gridvertex)
   if (debug) write(iulog,*) 'initMetagraph:after call to printgridvertex point #3.1'
   if (debug) write(iulog,*) 'initMetaGraph:thisprocessornumber is ',thisprocessornumber
   do j = 1,nelem ! count number of elements on this processor
     if (gridvertex(j)%processor_number.eq.thisprocessornumber) then
       metavertex%nmembers = metavertex%nmembers+1
     endif
   enddo
   if (debug) write(iulog,*) 'initMetagraph:point #4 '
!  Allocate space for the members of the MetaVertices
   if (debug) write(iulog,*) 'initMetagraph:point #4.1 i,metavertex%nmembers',i,metavertex%nmembers
   allocate(metavertex%members(metavertex%nmembers))
   if (debug) write(iulog,*) 'initMetagraph:point #5'
!  Set the identity of the members of the MetaVertices
   ic = 1
   do j = 1,nelem
     if (gridvertex(j)%processor_number.eq.thisprocessornumber) then
       metavertex%members(ic) = gridvertex(j)
       ic = ic+1
     endif
   enddo
   nedges = size(metaedge)
   if (debug) write(iulog,*) 'initMetagraph:point #6 nedges',nedges
!  Zero out all the edge numbers ... this should probably be
!  move to some initalization routine
   metaedge%number = 0
   metaedge%nmembers = 0
   metaedge%wgtp = 0
   metaedge%wgtp_ghost = 0
   do i = 1,nedges
     nullify(metaedge(i)%members)
   enddo
   if (debug) write(iulog,*) 'initMetagraph:point #7'
! Insert all the grid edges into the Meta Edges
   do i = 1,nelem_edge
!  Which Meta Edge does this grid edge belong
     head_processor_number = gridedge(i)%head%processor_number
     tail_processor_number = gridedge(i)%tail%processor_number
     call llfindedge(medgelist,tail_processor_number,head_processor_number,j,found)
     if (found) then
!  Increment the number of grid edges contained in the grid edge
!  and setup the pointers
       if (debug) write(iulog,*) 'initMetagraph:point #8'
       ii = gridedge(i)%tail_face
       wgtp = gridedge(i)%tail%wgtp(ii)
       metaedge(j)%nmembers = metaedge(j)%nmembers+1
       metaedge(j)%wgtp = metaedge(j)%wgtp+wgtp
       metaedge(j)%wgtp_ghost = metaedge(j)%wgtp_ghost+gridedge(i)%tail%wgtp_ghost(ii)
       if (debug) write(iulog,*) 'initMetagraph:point #9'
!  If this the first grid edge to be inserted into the Meta Edge
!  do some more stuff
       if (metaedge(j)%nmembers.eq.1) then
         if (debug) write(iulog,*) 'initMetagraph:point #10'
         metaedge(j)%number = j ! its identity
         metaedge(j)%type = gridedge_type(gridedge(i)) ! type of grid edge
         if (debug) write(iulog,*) 'initMetagraph:point #11'
!  Setup the pointer to the head and tail of the Vertex
         metaedge(j)%headvertex = head_processor_number
         metaedge(j)%tailvertex = tail_processor_number
         if (debug) write(iulog,*) 'initMetagraph:point #12'
!  Determine the number of edges for the Meta_Vertex
!  This is the number of processors to communicate with
         metavertex%nedges = metavertex%nedges+1
         if (debug) write(iulog,*) 'initMetagraph:point #13'
       endif
     endif
   enddo
   do i = 1,nedges
!  Allocate space for the member edges and edge index
     allocate(metaedge(i)%members(metaedge(i)%nmembers))
     allocate(metaedge(i)%edgeptrp(metaedge(i)%nmembers))
     allocate(metaedge(i)%edgeptrp_ghost(metaedge(i)%nmembers))
     metaedge(i)%edgeptrp(:) = 0
     metaedge(i)%edgeptrp_ghost(:) = 0
   enddo
   if (debug) write(iulog,*) 'initMetagraph:point #14'
!  Insert the edges into the proper meta edges
   allocate(icount(nelem_edge))
   icount = 1
   do i = 1,nelem_edge
     head_processor_number = gridedge(i)%head%processor_number
     tail_processor_number = gridedge(i)%tail%processor_number
     call llfindedge(medgelist,tail_processor_number,head_processor_number,j,found)
     if (found) then
       metaedge(j)%members(icount(j)) = gridedge(i)
       if (icount(j)+1.le.metaedge(j)%nmembers) then
         ii = gridedge(i)%tail_face
         wgtp = gridedge(i)%tail%wgtp(ii)
         metaedge(j)%edgeptrp(icount(j)+1) = metaedge(j)%edgeptrp(icount(j))+wgtp
         wgtp = gridedge(i)%tail%wgtp_ghost(ii)
         metaedge(j)%edgeptrp_ghost(icount(j)+1) = metaedge(j)%edgeptrp_ghost(icount(j))+wgtp
       endif
       if (debug) write(iulog,*) 'initMetagraph:point #15'
       icount(j) = icount(j)+1
     endif
   enddo
   deallocate(icount)
   if (debug) write(iulog,*) 'initMetagraph:point #16'
!  Fill in the additional edge information for the Meta Vertex
!  Allocate the number of edges incident to each Meta Vertex
   allocate(metavertex%edges(metavertex%nedges))
   ic = 1
   do j = 1,nedges
     metavertex%edges(j) = metaedge(j)
   enddo
   if (verbose) then
     print*
     write(iulog,*) "edge bundle list :(initmetagraph) "
     call printmetaedge(metaedge)
     write(iulog,*) 'initmetagrap:before last call to printmetavertex'
     call printmetavertex(metavertex)
   endif
   deallocate(metaedge)
   call llfree(medgelist)
   90 format('EDGE #',i2,2x,'TYPE ',i1,2x,'Processor numbers ',i2,'--->',i2)
   100 format(10x,i2,1x,'(',i1,')--->',i2,1x,'(',i1,') ')
!
   end subroutine initmetagraph
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module metagraph
!-------------------------------------------------------------------------------
