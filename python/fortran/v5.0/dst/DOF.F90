#include "KIM.h"
!-------------------------------------------------------------------------------
   module dof
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
   use kiapsbase, only: int_kind=>kim_int_kind, long_kind=>kim_long_kind, &
         real_kind=>kim_real8_kind, longdouble_kind=>kim_longdouble_kind
   use dimensions, only: np, npsq, nelem, nelemd
   use quadrature, only: quadrature_t
   use element, only: element_t, index_t
   use kiapsparallel, only: parallel_t, par_integer
   use edge, only: longedgebuffer_t, initlongedgebuffer, freelongedgebuffer, &
                                             longedgevpack, longedgevunpackmin
   use bndry, only: bndry_exchangev
   implicit none
!
   private
! public data
! public subroutines
   public :: global_dof
   public :: genlocaldof
   public :: printdofp
   public :: uniquepoints
   public :: putuniquepoints
   public :: uniquencolsp
   public :: uniquecoords
   public :: createuniqueindex
   public :: setelemoffset
   public :: createmetadata
!
   interface uniquepoints
     module procedure uniquepoints2d
     module procedure uniquepoints3d
     module procedure uniquepoints4d
   end interface
!
   interface putuniquepoints
     module procedure putuniquepoints2d
     module procedure putuniquepoints3d
     module procedure putuniquepoints4d
   end interface
!
   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine genlocaldof(ig, npts, ldof)
!
   integer(int_kind), intent(in   ) :: ig
   integer(int_kind), intent(in   ) :: npts
   integer(int_kind), intent(inout) :: ldof(:,:)
   integer(int_kind) :: i, j, npts2
!
   npts2 = npts*npts
   do j = 1,npts
     do i = 1, npts
       ldof(i,j) =(ig-1)*npts2+(j-1)*npts+i
     enddo
   enddo
!
   end subroutine genlocaldof
! ===========================================
! global_dof
!
! Compute the global degree of freedom for each element...
! ===========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine global_dof(par, elem)
!
   type(parallel_t), intent(in   ) :: par
   type(element_t) :: elem(:)
   type(longedgebuffer_t) :: edge
   real(real_kind) da ! area element
   type(quadrature_t) :: gp
   integer(int_kind) :: ldofp(np, np, nelemd)
   integer ii
   integer i, j, ig, ie
   integer kptr
   integer iptr
! ===================
! begin code
! ===================
!
   call initlongedgebuffer(edge,1)
! =================================================
! mass matrix on the velocity grid
! =================================================
   do ie = 1,nelemd
     ig = elem(ie)%vertex%number
     call genlocaldof(ig,np,ldofp(:,:,ie))
     kptr = 0
     call longedgevpack(edge,ldofp(:,:,ie),1,kptr,elem(ie)%desc)
   enddo
! ==============================
! Insert boundary exchange here
! ==============================
   call bndry_exchangev(par,edge)
   do ie = 1,nelemd
! we should unpack directly into elem(ie)%gdofV, but we dont have
! a VunpackMIN that takes integer*8.  gdofV integer*8 means
! more than 2G grid points.
     kptr = 0
     call longedgevunpackmin(edge,ldofp(:,:,ie),1,kptr,elem(ie)%desc)
     elem(ie)%gdofp(:,:) = ldofp(:,:,ie)
   enddo
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
   call freelongedgebuffer(edge)
!
   end subroutine global_dof
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine uniquepoints2d(idxunique, src, dest)
!
   type(index_t) :: idxunique
   real(real_kind) :: src(:,:)
   real(real_kind) :: dest(:)
   integer(int_kind) :: i, j, ii
!
   do ii = 1, idxunique%numuniquepts
     i = idxunique%ia(ii)
     j = idxunique%ja(ii)
     dest(ii) = src(i,j)
   enddo
!
   end subroutine uniquepoints2d
! putUniquePoints first zeros out the destination array, then fills the unique points of the
! array with values from src.  A boundary communication should then be called to fill in the
! redundent points of the array
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine putuniquepoints2d(idxunique, src, dest)
!
   type(index_t) :: idxunique
   real(real_kind), intent(in   ) :: src(:)
   real(real_kind), intent(  out) :: dest(:,:)
   integer(int_kind) :: i, j, ii
!
   dest = 0.0d0
   do ii = 1,idxunique%numuniquepts
     i = idxunique%ia(ii)
     j = idxunique%ja(ii)
     dest(i,j) = src(ii)
   enddo
!
   end subroutine putuniquepoints2d
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine uniquencolsp(elem, idxunique, cid)
!-------------------------------------------------------------------------------
   use element, only: getcolumnidp, element_t
!
   type(element_t), intent(in   ) :: elem
   type(index_t), intent(in   ) :: idxunique
   integer, intent(  out) :: cid(:)
   integer(int_kind) :: i, j, ii
!
   do ii = 1, idxunique%numuniquepts
     i = idxunique%ia(ii)
     j = idxunique%ja(ii)
     cid(ii) = getcolumnidp(elem,i,j)
   enddo
!
   end subroutine uniquencolsp
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine uniquecoords(idxunique, src, lat, lon)
!-------------------------------------------------------------------------------
   use coordinatesystems, only: spherical_polar_t
!
   type(index_t), intent(in   ) :: idxunique
   type(spherical_polar_t) :: src(:,:)
   real(real_kind), intent(  out) :: lat(:)
   real(real_kind), intent(  out) :: lon(:)
   integer(int_kind) :: i, j, ii
!
   do ii = 1, idxunique%numuniquepts
     i = idxunique%ia(ii)
     j = idxunique%ja(ii)
     lat(ii) = src(i,j)%lat
     lon(ii) = src(i,j)%lon
   enddo
!
   end subroutine uniquecoords
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine uniquepoints3d(idxunique, nlyr, src, dest)
!
   type(index_t) :: idxunique
   integer(int_kind) :: nlyr
   real(real_kind) :: src(:,:,:)
   real(real_kind) :: dest(:,:)
   integer(int_kind) :: i, j, k, ii
!
   do ii = 1, idxunique%numuniquepts
     i = idxunique%ia(ii)
     j = idxunique%ja(ii)
     do k = 1,nlyr
       dest(ii,k) = src(i,j,k)
     enddo
   enddo
!
   end subroutine uniquepoints3d
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine uniquepoints4d(idxunique, d3, d4, src, dest)
!
   type(index_t) :: idxunique
   integer(int_kind) :: d3, d4
   real(real_kind) :: src(:,:,:,:)
   real(real_kind) :: dest(:,:,:)
   integer(int_kind) :: i, j, k, n, ii
!
   do n = 1, d4
     do k = 1, d3
       do ii = 1, idxunique%numuniquepts
         i = idxunique%ia(ii)
         j = idxunique%ja(ii)
         dest(ii,k,n) = src(i,j,k,n)
       enddo
       enddo
   enddo
!
   end subroutine uniquepoints4d
! putUniquePoints first zeros out the destination array, then fills the unique points of the
! array with values from src.  A boundary communication should then be called to fill in the
! redundent points of the array
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine putuniquepoints3d(idxunique, nlyr, src, dest)
!
   type(index_t) :: idxunique
   integer(int_kind) :: nlyr
   real(real_kind), intent(in   ) :: src(:,:)
   real(real_kind), intent(  out) :: dest(:,:,:)
   integer(int_kind) :: i, j, k, ii
!
   dest = 0.0d0
   do k = 1,nlyr
     do ii = 1, idxunique%numuniquepts
       i = idxunique%ia(ii)
       j = idxunique%ja(ii)
       dest(i,j,k) = src(ii,k)
     enddo
   enddo
!
   end subroutine putuniquepoints3d
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine putuniquepoints4d(idxunique, d3, d4, src, dest)
!
   type(index_t) :: idxunique
   integer(int_kind) :: d3, d4
   real(real_kind), intent(in   ) :: src(:,:,:)
   real(real_kind), intent(  out) :: dest(:,:,:,:)
   integer(int_kind) :: i, j, k, n, ii
!
   dest = 0.0d0
   do n = 1,d4
     do k = 1, d3
       do ii = 1, idxunique%numuniquepts
         i = idxunique%ia(ii)
         j = idxunique%ja(ii)
         dest(i,j,k,n) = src(ii,k,n)
       enddo
       enddo
   enddo
!
   end subroutine putuniquepoints4d
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine setelemoffset(par, elem, globaluniquecolsp)
#ifdef _MPI
!-------------------------------------------------------------------------------
   use kiapsparallel, only: par_sum
   use kiapsparmpi, only: par_allreduce
#endif
!
   type(parallel_t) :: par
   type(element_t) :: elem(:)
   integer, intent(  out) :: globaluniquecolsp
   integer(int_kind), allocatable :: numelemp(:), numelem2p(:)
   integer(int_kind), allocatable :: numelemv(:), numelem2v(:)
   integer(int_kind), allocatable :: goffset(:)
   integer(int_kind) :: ie, ig, nprocs, ierr
   logical, parameter :: debug = .false.
!
   nprocs = par%nprocs
   allocate(numelemp(nelem))
   allocate(numelem2p(nelem))
   allocate(goffset(nelem))
   numelemp = 0;numelem2p = 0;goffset = 0
   do ie = 1,nelemd
     ig = elem(ie)%globalid
     numelemp(ig) = elem(ie)%idxp%numuniquepts
   enddo
#ifdef _MPI
!     call MPI_Allreduce(numElemP,numElem2P,nelem,PAR_INTEGER,PAR_SUM,par%comm,ierr)
   call par_allreduce(numelemp(1),numelem2p(1),nelem,par_sum)
#else
   numelem2p = numelemp
#endif
   goffset(1) = 1
   do ig = 2,nelem
     goffset(ig) = goffset(ig-1)+numelem2p(ig-1)
   enddo
   do ie = 1,nelemd
     ig = elem(ie)%globalid
     elem(ie)%idxp%uniqueptoffset = goffset(ig)
   enddo
   globaluniquecolsp = goffset(nelem)+numelem2p(nelem)-1
   deallocate(numelemp)
   deallocate(numelem2p)
   deallocate(goffset)
!
   end subroutine setelemoffset
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine createuniqueindex(ig, gdof, idx)
!
   integer(int_kind) :: ig
   type(index_t) :: idx
   integer(long_kind) :: gdof(:,:)
   integer, allocatable :: ldof(:,:)
   integer :: i, j, ii, npts
!
   npts = size(gdof,dim=1)
   allocate(ldof(npts,npts))
! ====================
! Form the local DOF
! ====================
   call genlocaldof(ig,npts,ldof)
   ii = 1
   do j = 1,npts
     do i = 1, npts
! ==========================
! check for point ownership
! ==========================
       if (gdof(i, j).eq.ldof(i, j)) then
         idx%ia(ii) = i
         idx%ja(ii) = j
         ii = ii+1
       endif
       enddo
   enddo
   idx%numuniquepts = ii-1
   deallocate(ldof)
!
   end subroutine createuniqueindex
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine createmetadata(par, elem, subelement_corners, fdofp)
!
   type(parallel_t), intent(in   ) :: par
   type(element_t), target :: elem(:)
   integer, intent(  out), optional :: subelement_corners((np-1)*(np-1)*nelemd, 4)
   integer(int_kind), optional :: fdofp(np, np, nelemd)
   type(index_t), pointer :: idx
   type(longedgebuffer_t) :: edge
   integer :: i, j, ii, ie, base
   integer(long_kind), pointer :: gdof(:,:)
   integer :: fdofp_local(np, np, nelemd)
!
   call initlongedgebuffer(edge,1)
   fdofp_local = 0
   do ie = 1,nelemd
     idx=>elem(ie)%idxp
     do ii = 1,idx%numuniquepts
       i = idx%ia(ii)
       j = idx%ja(ii)
       fdofp_local(i,j,ie) =-(idx%uniqueptoffset+ii-1)
     enddo
     call longedgevpack(edge,fdofp_local(:,:,ie),1,0,elem(ie)%desc)
   enddo
   call bndry_exchangev(par,edge)
   do ie = 1,nelemd
     base =(ie-1)*(np-1)*(np-1)
     call longedgevunpackmin(edge,fdofp_local(:,:,ie),1,0,elem(ie)%desc)
     if (present(subelement_corners)) then
       ii = 0
       do j = 1,np-1
         do i = 1, np-1
           ii = ii+1
           subelement_corners(base+ii,1) =-fdofp_local(i,j,ie)
           subelement_corners(base+ii,2) =-fdofp_local(i,j+1,ie)
           subelement_corners(base+ii,3) =-fdofp_local(i+1,j+1,ie)
           subelement_corners(base+ii,4) =-fdofp_local(i+1,j,ie)
         enddo
       enddo
     endif
   enddo
   if (present(fdofp)) then
     fdofp =-fdofp_local
   endif
!
   end subroutine createmetadata
! ==========================================
!  PrintDofP
!
!   Prints the degree of freedom
! ==========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine printdofp(elem)
!-------------------------------------------------------------------------------
   implicit none
!
   type(element_t), intent(in   ) :: elem(:)
   integer :: ie, nse, i, j
!
   nse = size(elem)
   do ie = 1,nse
     print*,'Element # ',elem(ie)%vertex%number
     do j = np,1,-1
       write(6,*)(elem(ie)%gdofp(i,j),i = 1,np)
     enddo
   enddo
   10 format('I5')
!
   end subroutine printdofp
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module dof
!-------------------------------------------------------------------------------
