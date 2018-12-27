!-------------------------------------------------------------------------------
   module space_filling_curve
!-------------------------------------------------------------------------------
!
!  abstract : 
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only : i4, l4, r4, r8
   use mesh_control, only : ww, ee, ss, nn, axis_x, axis_y, axis_xy, axis_mxy, axis_0
   use mesh_control, only : mesh_convert, mesh_dir_convert, mesh_print
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   type curve_t
     integer(i4), dimension(:,:), allocatable :: major  !(2, 2) or(3, 3) or(5, 5)
     integer(i4), dimension(:,:), allocatable :: joiner !(2, 2) or(3, 3) or(5, 5)
   end type curve_t
!
   integer(i4), parameter :: ncrvs  = 3 ! # of curves(2, 3, 5)
   integer(i4), parameter :: nprims = 4 ! # of primitives
!
   type sfc_t
     integer(i4) :: n
     ! factorization
     logical(i4) :: isfactorized
     integer(i4) :: nfacs
     integer(i4), dimension(:), allocatable :: factors
     ! curve and primitives
     integer(i4) :: imajor, ijoiner
     integer(i4) :: ncrvs  = ncrvs  !3 ! # of curves(2, 3, 5)
     integer(i4) :: nprims = nprims !4 ! # of primitives
     integer(i4), dimension(ncrvs) :: nums
     type(curve_t), dimension(nprims, ncrvs) :: prim !(nprims, ncrvs)
     ! result: major and joiner
     type(curve_t) :: curve ! result major and joiner
     ! mover (for indexing)
     integer(i4), dimension(2, 4) :: move
   end type sfc_t
!
   public :: sfc_t, sfc_initialize, sfc_finalize, sfc_gen_curve, sfc_get_sfc
   public :: do_factorize
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine sfc_initialize(sfc, nn)
!-------------------------------------------------------------------------------
   implicit none
!
   type(sfc_t), intent(inout) :: sfc
   integer(i4), intent(in   ) :: nn
! local variables
   integer(i4) :: iprim, icrv
!
   sfc%n = nn
   sfc%isfactorized = do_factorize(nn,sfc%nfacs,sfc%factors)
   if (.not.sfc%isfactorized) then
     print*,'n is not factorable...',sfc%n
     print*,'nfac, facs = ',sfc%nfacs,sfc%factors
     stop
   endif
!
   sfc%imajor = 2
   sfc%ijoiner = sfc%imajor
   sfc%nums =(/2,3,5/)
   do icrv = 1,sfc%ncrvs
   do iprim = 1,sfc%nprims
     allocate(sfc%prim(iprim,icrv)%major(sfc%nums(icrv),sfc%nums(icrv)))
     allocate(sfc%prim(iprim,icrv)%joiner(sfc%nums(icrv),sfc%nums(icrv)))
   enddo
   enddo
!
   sfc%move = reshape((/-1,0,1,0,0,-1,0,1/),shape =(/2,4/))
!
   allocate(sfc%curve%major(nn,nn))
   allocate(sfc%curve%joiner(nn,nn))
!
   call gen_primitive_curve(sfc)

!
   return
   end subroutine sfc_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine sfc_finalize(sfc)
!-------------------------------------------------------------------------------
   implicit none
!
   type(sfc_t), intent(  out) :: sfc
! local variables
   integer(i4) :: iprim, icrv
!
   if (allocated(sfc%factors)) deallocate(sfc%factors)
!
   do icrv = 1,sfc%ncrvs
   do iprim = 1,sfc%nprims
     if (allocated(sfc%prim(iprim,icrv)%major)) deallocate(sfc%prim(iprim,icrv)%major)
     if (allocated(sfc%prim(iprim,icrv)%joiner)) deallocate(sfc%prim(iprim,icrv)%joiner)
   enddo
   enddo
!
   if (allocated(sfc%curve%major)) deallocate(sfc%curve%major)
   if (allocated(sfc%curve%joiner)) deallocate(sfc%curve%joiner)
!
   return
   end subroutine sfc_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function do_factorize(nn, nfacs, factors) result(res)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                            intent(in   ) :: nn
   integer(i4),                            intent(  out) :: nfacs
   integer(i4), dimension(:), allocatable, intent(  out) :: factors
   logical(l4) :: res
! local variables
   integer(i4) :: tmp, tmpd, i
   logical(l4) :: found
   integer(i4), dimension(100) :: maxfactors
   !integer(i4), dimension(3) :: base =(/2,3,5/)
   integer(i4), dimension(3) :: base =(/5,3,2/)
!
   nfacs = 0
   do i = 1,3
     tmp = nn
     found = .true.
     do while(found)
       found = .false.
       tmpd = tmp/base(i)
       if (tmpd*base(i).eq.tmp) then
         nfacs = nfacs+1
         maxfactors(nfacs) = base(i)
         tmp = tmpd
         found = .true.
       endif
     enddo
   enddo
!
   if (nfacs.gt.0) then
     res = .true.
     if (allocated(factors)) deallocate(factors)
     allocate(factors(nfacs))
     factors(:) = maxfactors(1:nfacs)
   else
     res = .false.
   endif
!
   tmp = 1
   do i = 1,nfacs
     tmp = tmp*factors(i)
   enddo
   if (tmp.ne.nn) res = .false.
!
   end function do_factorize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine gen_primitive_curve(sfc)
!-------------------------------------------------------------------------------
   implicit none
!
   type(sfc_t), intent(inout) :: sfc
! local variables
   integer(i4) :: ll, rr, uu, dd, icrv, iprim
   integer(i4) :: axis(4)
!
   ll = ww
   rr = ee
   uu = nn
   dd = ss
   axis =(/-1,axis_0,axis_xy,axis_mxy/)
!
   ! 2: primitive 1
   sfc%prim(1,1)%major(:,:) = reshape((/ll,ll,uu,dd/),shape =(/2,2/))
   sfc%prim(1,1)%joiner(:,:) = reshape((/uu,ll,-1,dd/),shape =(/2,2/))
   ! 3: primitive 1
   sfc%prim(1,2)%major(:,:) = reshape((/ll,ll,ll,rr,uu,dd,ll,uu,dd/),shape =(/3,3/))
   sfc%prim(1,2)%joiner(:,:) = reshape((/uu,ll,ll,rr,uu,dd,-1,ll,dd/),shape =(/3,3/))
   ! 5: primitive 1
   sfc%prim(1,3)%major(:,:) = reshape((/ll,ll,ll,ll,ll,rr,uu,dd,uu,dd,ll,uu,dd,uu,dd,rr,uu,dd,rr,rr,ll,uu,dd,ll,ll/),shape =(/5,5/))
   sfc%prim(1,3)%joiner(:,:) = reshape((/uu,ll,ll,uu,ll,rr,uu,dd,ll,dd,uu,ll,rr,uu,dd,rr,uu,dd,rr,dd,-1,ll,dd,ll,ll/),shape =(/5,5/))
!
   do icrv = 1,sfc%ncrvs
   do iprim = 2,sfc%nprims
     call mesh_dir_convert(axis(iprim),sfc%prim(1,icrv)%major(:,:),sfc%prim(iprim,icrv)%major(:,:))
     call mesh_dir_convert(axis(iprim),sfc%prim(1,icrv)%joiner(:,:),sfc%prim(iprim,icrv)%joiner(:,:))
   enddo
   enddo
!
   return
   end subroutine gen_primitive_curve
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine sfc_gen_curve(sfc)
!-------------------------------------------------------------------------------
   implicit none
!
   type(sfc_t), intent(inout) :: sfc
!local
   integer(i4) :: nn, i, j, nfacs, imajor, ijoiner, idx, nlevs, ilev, ista, iend
   integer(i4), dimension(:), allocatable :: size_levs, size_prim_levs
   type(curve_t), dimension(:), allocatable :: curves
!
   imajor = sfc%imajor
   ijoiner = sfc%ijoiner
!
   idx = 1
   ilev = 0
   ista = 0
   iend = 0
!
   nfacs = sfc%nfacs
!
   nlevs = nfacs+1
!
   allocate(size_prim_levs(nlevs))
   allocate(size_levs(nlevs))
   size_prim_levs(1) = 1
   do ilev = 2,nlevs
     size_prim_levs(ilev) = sfc%factors(nfacs-ilev+2)
   enddo
!
   size_levs(:) = 1
   do ilev = 2,nlevs
     size_levs(ilev) = size_prim_levs(ilev)*size_levs(ilev-1)
   enddo
!
   allocate(curves(nlevs))
   do ilev = 1,nlevs
     nn = size_levs(ilev)
     allocate(curves(ilev)%major(nn,nn))
     allocate(curves(ilev)%joiner(nn,nn))
   enddo
!
  ! main algorithm
   curves(1)%major(1,1) = imajor
   curves(1)%joiner(1,1) = ijoiner
   do ilev = 2,nlevs
     nn = size_levs(ilev-1)
     do j = 1,nn
     do i = 1,nn
       imajor = curves(ilev-1)%major(i,j)
       ijoiner = curves(ilev-1)%joiner(i,j)
       call gen_next_major_joiner(sfc,ilev,size_prim_levs(ilev),i,j,imajor,ijoiner,curves(ilev)%major,curves(ilev)%joiner)
     enddo
     enddo
   enddo
!
   sfc%curve%major = curves(nlevs)%major
   sfc%curve%joiner = curves(nlevs)%joiner
!
   deallocate(size_prim_levs)
   deallocate(size_levs)
   do ilev = 1,nlevs
     deallocate(curves(ilev)%major)
     deallocate(curves(ilev)%joiner)
   enddo
   deallocate(curves)
!
   return
   end subroutine sfc_gen_curve
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine gen_next_major_joiner(sfc, ilev, nn, ipos, jpos, imajor, ijoiner, major, joiner)
!-------------------------------------------------------------------------------
   implicit none
!
   type(sfc_t),  intent(in   ) :: sfc
   integer(i4),  intent(in   ) :: ilev, nn, ipos, jpos, imajor, ijoiner
   integer(i4), dimension(:,:) :: major, joiner
! local variables
   integer(i4) :: is, js, ie, je, ibase, i, j, ii, jj
!
   call get_start_end(nn,imajor,is,js,ie,je)
!
   if (nn.eq.2) then
     ibase = 1
   elseif (nn.eq.3) then
     ibase = 2
   elseif (nn.eq.5) then
     ibase = 3
   endif
!
   do j = 1,nn
   do i = 1,nn
     ii = nn*(ipos-1)+i
     jj = nn*(jpos-1)+j
     major(ii,jj) = sfc%prim(imajor,ibase)%major(i,j)
     if (i.eq.ie.and.j.eq.je) then
       joiner(ii,jj) = ijoiner
     else
       joiner(ii,jj) = sfc%prim(imajor,ibase)%joiner(i,j)
     endif
   enddo
   enddo
!
   return
   end subroutine gen_next_major_joiner
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_start_end(n, imajor, is, js, ie, je)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: n, imajor
   integer(i4), intent(  out) :: is, js, ie, je
!
   if (imajor.eq.1) then
     is = n
     js = n
     ie = 1
     je = n
   elseif (imajor.eq.2) then
     is = 1
     js = 1
     ie = n
     je = 1
   elseif (imajor.eq.3) then
     is = n
     js = n
     ie = n
     je = 1
   elseif (imajor.eq.4) then
     is = 1
     js = 1
     ie = 1
     je = n
   else
     print*,'check imajor in get_start_end...'
     stop
   endif
!
   return
   end subroutine get_start_end
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine sfc_get_sfc(sfc, mesh)
!-------------------------------------------------------------------------------
   implicit none
!
   type(sfc_t),                              intent(in   ) :: sfc
   integer(i4), dimension(:,:), allocatable, intent(  out) :: mesh
! local variables
   integer(i4) :: i, n, is, js, ie, je, idx, ip, jp, iinc, jinc, nlevs, dir
!
   n = sfc%n
   nlevs = sfc%nfacs+1
   if (allocated(mesh)) deallocate(mesh)
   allocate(mesh(n,n))
   mesh(:,:) =-1
!
   call get_start_end(n,sfc%imajor,is,js,ie,je)
!
   idx = 1
   ip = is
   jp = js
   do i = 1,n*n
     mesh(ip,jp) = idx
     dir = sfc%curve%joiner(ip,jp)
     iinc = sfc%move(1,dir)
     jinc = sfc%move(2,dir)
     ip = ip+iinc
     jp = jp+jinc
     idx = idx+1
   enddo
!
   return
   end subroutine sfc_get_sfc
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module space_filling_curve
