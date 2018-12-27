!-------------------------------------------------------------------------------
   module w3gdatmd
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
   use w3odatmd, only: iaproc, naproc
   implicit none
!
   integer :: nx, ny, nsea, nseal
   integer, dimension(:,:), allocatable :: mapsf, mapsta
   public :: nx, ny, nsea, nseal, mapsf, mapsta
   public :: initialize, finalize
!
   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine initialize()
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: ix, iy, idx, iproc
!
!
   nx = 10
   ny = 5
!
   allocate(mapsf(nx*ny,2))
   allocate(mapsta(ny,nx))
   mapsta(:,:) = 1
   mapsta(2,2) = 0
   mapsta(2,3) = 0
   mapsta(2,4) = 0
   mapsta(3,2) = 0
   mapsta(3,3) = 0
   mapsta(4,6) = 0
   mapsta(4,7) = 0
   mapsta(4,8) = 0
   mapsta(4,9) = 0
   mapsta(5,3) = 0
   mapsta(5,6) = 0
   mapsta(5,7) = 0
!
   nsea = 0
   nseal = 0
   iproc = 0
   do iy = 1,ny
     do ix = 1, nx
       idx = ix+(iy-1)*nx
       mapsf(idx,1) = ix
       mapsf(idx,2) = iy
       if (mapsta(iy,ix).eq.1) then
         nsea = nsea+1
         iproc = iproc+1
         if (iproc.eq.naproc+1) iproc = 1
         if (iproc.eq.iaproc) then
           nseal = nseal+1
         endif
       endif
     enddo
   enddo
!
   return
!
   end subroutine initialize
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine finalize()
!-------------------------------------------------------------------------------
   implicit none
!
   deallocate(mapsf,mapsta)
   return
!
   end subroutine finalize
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module w3gdatmd
!-------------------------------------------------------------------------------
