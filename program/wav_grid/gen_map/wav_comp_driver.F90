!-------------------------------------------------------------------------------
   module wav_comp_driver
!-------------------------------------------------------------------------------
!
!  abstract : wave model component
! 
!  history log : 
!    2016-03-18   junghan kim     initial setup
!    2016-04-15   junghan kim     the version for application to the KIM 
!                                 (KIM-CPL v1.0)
!    2017-06-12   junghan kim     component for wave-watch3
!
!-------------------------------------------------------------------------------
   use kinds, only: i4, r8
   use w3odatmd, only: iaproc, naproc
   use w3gdatmd, only: nx, ny, nsea, nseal, mapsf, mapsta
   implicit none
!
   private
!
! variables for cpl-ww3
   integer(i4)                            :: gvarsize, lvarsize
   integer(i4), dimension(:), allocatable :: map, map_mask
   real(r8)   , dimension(:), allocatable :: cpl_u10m, cpl_v10m, cpl_znt
!
   public :: make_model_map, convert_wav2cpl, convert_cpl2wav
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine make_model_map(gsize, lsize, umap, mask)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),                            intent(inout) :: gsize, lsize
   integer(i4), dimension(:), allocatable, intent(inout) :: umap, mask
! local variables
   integer(i4) :: gnx, gny, nland, nlandl, eny
   integer(i4) :: i, j, js, ix, iy, iix, iiy, iproc
   integer(i4), dimension(:), allocatable :: ocn_map, lnd_map, tmp_map
!
   gnx = nx
   if (ny.eq.5) then
     gny = 7
   else
     if (iaproc.eq.0) then
       print *, 'not support (nx, ny)'
       stop
     endif
   endif
   gsize = gnx*gny
   nland = gsize-nsea
   eny   = gny-ny
!
! step 1: ocean map
   allocate(ocn_map(nseal))
   do i = 1,nseal
     j  = iaproc+(i-1)*naproc
     ix = mapsf(j,1)
     iy = mapsf(j,2)+(eny/2)
     ocn_map(i) = ix+(iy-1)*gnx
   enddo
!
! step 2: land map
   allocate(tmp_map(nland))
!print *, 'DEBUG: gsize, nsea, nland, eny/2 = ', gsize, nsea, nland, eny/2
!print *, ' '
   ! ww3
   nlandl = 0
   iproc  = 0
   do ix = 1,gnx
     do iy = 1,ny
       iix = ix; iiy = eny/2+iy
       if (mapsta(iy,ix).eq.0) then  ! check
         iproc = iproc+1
         if (iproc.eq.naproc+1) iproc = 1
         if (iproc.eq.iaproc) then
           nlandl = nlandl+1
           tmp_map(nlandl) = iix+(iiy-1)*gnx
         endif
       endif
     enddo
   enddo
   ! poles (south)
   do iy = 1,eny/2
     do ix = 1,gnx
       iproc = iproc+1
       if (iproc.eq.naproc+1) iproc = 1
       if (iproc.eq.iaproc) then
         nlandl = nlandl+1
         tmp_map(nlandl) = ix+(iy-1)*gnx
       endif
     enddo
   enddo
   ! poles (north)
   do iy = ny+eny/2+1,gny
     do ix = 1,gnx
       iproc = iproc+1
       if (iproc.eq.naproc+1) iproc = 1
       if (iproc.eq.iaproc) then
         nlandl = nlandl+1
         tmp_map(nlandl) = ix+(iy-1)*gnx
       endif
     enddo
   enddo
!
   allocate(lnd_map(nlandl))
   do i = 1,nlandl
     lnd_map(i) = tmp_map(i)
   enddo
   deallocate(tmp_map)
!
! entire domain map
   lsize = nseal+nlandl
   allocate(umap(lsize))
   umap(1:nseal)       = ocn_map(:)
   umap(nseal+1:lsize) = lnd_map(:)
   call quick_sort(umap,1,lsize)
!print *, 'DEBUG: ocn_map = '
!print *, ocn_map
!print *, 'DEBUG: lnd_map = '
!print *, lnd_map
!
! ocn mask
   allocate(mask(lsize))
   mask(:) = -1
   js = 1
   do i = 1,lsize
     do j = js,nseal
       if (umap(i).eq.ocn_map(j)) then
         mask(i) = 1
         js      = j
         exit
       endif
     enddo
#if 0
     if (mask(i).eq.-1) then
!       if (iaproc.eq.0) then
         print *, 'check umap and ocn_map'
         stop
!       endif
     endif
#endif
   enddo
!
   deallocate(ocn_map)
   deallocate(lnd_map)
!
   return
   end subroutine make_model_map
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cyclic_decompose(nprocs, nsize, grid)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),                   intent(in   ) :: nprocs, nsize
   integer(i4), dimension(nsize), intent(inout) :: grid
! local variables
   integer :: i, iproc
!
   iproc = 1
   do i = 1,nsize
     grid(i) = iproc
     iproc   = iproc+1
     if (iproc>nprocs) iproc = 1
   enddo
!
   return
   end subroutine cyclic_decompose
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   recursive subroutine quick_sort(array, first, last)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), dimension(:), intent(inout) :: array
   integer(i4),               intent(in   ) :: first, last
!
   real(r8) :: x, t
   integer  :: i, j
!
   x = array((first+last)/2)
   i = first
   j = last
!
   do
     do while (array(i)<x)
       i = i+1
     end do
     do while (x<array(j))
       j = j-1
     end do
     if (i>=j) exit
     t        = array(i)
     array(i) = array(j)
     array(j) = t
     i = i+1
     j = j-1
   end do
!
   if (first<i-1) call quick_sort(array,first,i-1)
   if (j+1<last)  call quick_sort(array,j+1,last)
!
   return
   end subroutine quick_sort
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine convert_wav2cpl(mask, inbuf, oubuf)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4), dimension(:), intent(in   ) :: mask
   real(r8)   , dimension(:), intent(in   ) :: inbuf
   real(r8)   , dimension(:), intent(  out) :: oubuf
! local variables
   integer :: i, cnt, nmask
!
   nmask = size(mask)
   if (nmask.ne.size(oubuf)) then
     if (iaproc.eq.0) then
       print *, 'check size(mask) and size(oubuf)..'
       stop
     endif
   endif
!
   oubuf(:) = 0.0
   cnt = 0
   do i = 1,nmask
      if (mask(i).eq.1) then
         cnt = cnt+1
         oubuf(i) = inbuf(cnt)
      endif
   enddo
!
   return
   end subroutine convert_wav2cpl
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine convert_cpl2wav(mask, inbuf, oubuf)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4), dimension(:), intent(in   ) :: mask
   real(r8)   , dimension(:), intent(in   ) :: inbuf
   real(r8)   , dimension(:), intent(  out) :: oubuf
! local variables
   integer :: i, cnt, nmask
!
   nmask = size(mask)
   if (nmask.ne.size(inbuf)) then
     if (iaproc.eq.0) then
       print *, 'check size(mask) and size(inbuf)..'
       stop
     endif
   endif
!
   cnt = 0
   do i = 1,size(inbuf)
      if (mask(i).eq.1) then
          cnt = cnt+1
          oubuf(cnt) = inbuf(i)
      endif
   enddo
!
   return
   end subroutine convert_cpl2wav
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module wav_comp_driver
!-------------------------------------------------------------------------------
