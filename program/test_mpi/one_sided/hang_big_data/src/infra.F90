!-------------------------------------------------------------------------------
   module infra
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    2017-03-09  junghan kim    initial setup
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds,    only: i4, l4, r8
   use parallel, only: decompose_1d
!
   implicit none
!
   integer(i4), dimension(:),   allocatable :: map
   real(r8),    dimension(:,:), allocatable :: var
   real(r8),    dimension(:,:), allocatable :: gvar
!
   public :: map, var, gvar
   public :: infra_initialize, infra_finalize, check_var
!
   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function infra_initialize(nprocs, rank, gnx, gny, continuous) result(nx)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: nprocs, rank
   integer(i4), intent(in   ) :: gnx, gny
   logical(l4), intent(in   ) :: continuous
   integer(i4)                :: nx
!
! local variables
!
!
   nx = decompose_1d(1,gnx,nprocs,rank)
!
   call gen_map(nprocs,rank,gnx,continuous,map)
   call gen_var(gnx,gny,nx,map,var)
   allocate(gvar(gnx,gny))
   gvar = -1
!
   return
!
   end function infra_initialize
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine infra_finalize()
!-------------------------------------------------------------------------------
   implicit none
!
   deallocate(map)
   deallocate(var)
   deallocate(gvar)
!
   return
!
   end subroutine infra_finalize
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine gen_map(nprocs, rank, gnx, continuous, umap)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: nprocs, rank
   integer(i4), intent(in   ) :: gnx
   logical(l4), intent(in   ) :: continuous
   integer(i4), dimension(:), allocatable, intent(inout) :: umap
!
! local variables
!
   integer(i4) :: j
   integer(i4) :: lsizes, ista, iend
!
! set local size,ista,iend => size of maps
!
   lsizes = decompose_1d(1,gnx,nprocs,rank,ista,iend)
   if (allocated(umap)) deallocate(umap)
   allocate(umap(lsizes))
   if (continuous) then
! 0 0 0 1 1 2 2
     do j = 1, lsizes
       umap(j) = ista+j-1
     enddo
   else
! 0 1 2 0 1 2 0
     do j = 1, lsizes
       umap(j) = 1+rank+nprocs*(j-1)
     enddo
   endif
!
   return
!
   end subroutine gen_map
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine gen_var(gnx, gny, nx, umap, uvar)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                              intent(in   ) :: gnx, gny, nx
   integer(i4), dimension(nx),               intent(in   ) :: umap
   real(r8),    dimension(:,:), allocatable, intent(inout) :: uvar
!
! local variables
!
   integer(i4) :: i, j
!
   allocate(uvar(nx,gny))
   do j = 1,gny
   do i = 1,nx
     uvar(i,j) = umap(i)+1000*j
   enddo
   enddo
!
   return
!
   end subroutine gen_var
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function check_var(gnx, gny, uvar) result(issame)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: gnx, gny
   real(r8), dimension(gnx,gny), intent(in   ) :: uvar
   logical(l4)                                 :: issame
!
   integer(i4) :: i, j
!
   issame = .true.
   do j = 1,gny
   do i = 1,gnx
     if (uvar(i,j).ne.real(i+1000*j,r8)) then
       issame = .false.
       exit
     endif
   enddo
   enddo
!
   return
!
   end function check_var
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module infra
!-------------------------------------------------------------------------------
