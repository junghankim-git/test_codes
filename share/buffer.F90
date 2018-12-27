!-------------------------------------------------------------------------------
   module buffer
!-------------------------------------------------------------------------------
!
!  abstract :  bi-linear interpolation module
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!    2017-07-19  junghan kim    add min/max
!
!  structure : 
!
!-------------------------------------------------------------------------------
!   use netcdf
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   integer, parameter :: i4 = 4, r8 = 8, l4 = 4
!
   type buffer_t
   integer(i4) :: ndims =-1
   integer(i4), dimension(:),      allocatable :: dimsizes
   real(r8), dimension(:),         allocatable :: d1
   real(r8), dimension(:,:),       allocatable :: d2
   real(r8), dimension(:,:,:),     allocatable :: d3
   real(r8), dimension(:,:,:,:),   allocatable :: d4
   real(r8), dimension(:,:,:,:,:), allocatable :: d5
   end type buffer_t
!
   public :: buffer_t
   public :: buffer_initialize, buffer_finalize, buffer_min, buffer_max
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function buffer_min(buffer) result(minv)
!-------------------------------------------------------------------------------
   implicit none
!
   type(buffer_t), intent(in   ) :: buffer
   real(r8)                      :: minv
! local variables
   integer(i4) :: ndims
!
   ndims = buffer%ndims
!
   if     (ndims==1) then
     minv = minval(buffer%d1)
   elseif (ndims==2) then
     minv = minval(buffer%d2)
   elseif (ndims==3) then
     minv = minval(buffer%d3)
   elseif (ndims==4) then
     minv = minval(buffer%d4)
   elseif (ndims==5) then
     minv = minval(buffer%d5)
   else
     print*,'check ndims in buffer_initialize...',ndims
     stop
   endif
!
   return
   end function buffer_min
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function buffer_max(buffer) result(maxv)
!-------------------------------------------------------------------------------
   implicit none
!
   type(buffer_t), intent(in   ) :: buffer
   real(r8)                      :: maxv
! local variables
   integer(i4) :: ndims
!
   ndims = buffer%ndims
!
   if     (ndims==1) then
     maxv = maxval(buffer%d1)
   elseif (ndims==2) then
     maxv = maxval(buffer%d2)
   elseif (ndims==3) then
     maxv = maxval(buffer%d3)
   elseif (ndims==4) then
     maxv = maxval(buffer%d4)
   elseif (ndims==5) then
     maxv = maxval(buffer%d5)
   else
     print*,'check ndims in buffer_initialize...',ndims
     stop
   endif
!
   return
   end function buffer_max
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine buffer_initialize(buffer, ndims, dimsizes)
!-------------------------------------------------------------------------------
   implicit none
!
   type(buffer_t),                intent(  out) :: buffer
   integer(i4),                   intent(in   ) :: ndims
   integer(i4), dimension(ndims), intent(in   ) :: dimsizes
!
   buffer%ndims = ndims
   allocate(buffer%dimsizes(ndims))
   buffer%dimsizes = dimsizes
!
   if (ndims==1) then
     allocate(buffer%d1(dimsizes(1)))
   elseif (ndims==2) then
     allocate(buffer%d2(dimsizes(1),dimsizes(2)))
   elseif (ndims==3) then
     allocate(buffer%d3(dimsizes(1),dimsizes(2),dimsizes(3)))
   elseif (ndims==4) then
     allocate(buffer%d4(dimsizes(1),dimsizes(2),dimsizes(3),dimsizes(4)))
   elseif (ndims==5) then
     allocate(buffer%d5(dimsizes(1),dimsizes(2),dimsizes(3),dimsizes(4),dimsizes(5)))
   else
     print*,'check ndims in buffer_initialize...',ndims
     stop
   endif
!
   return
   end subroutine buffer_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine buffer_finalize(buffer)
!-------------------------------------------------------------------------------
   implicit none
!
   type(buffer_t), intent(inout) :: buffer
!
   buffer%ndims =-1
   if (allocated(buffer%dimsizes)) deallocate(buffer%dimsizes)
   if (allocated(buffer%d1)) deallocate(buffer%d1)
   if (allocated(buffer%d2)) deallocate(buffer%d2)
   if (allocated(buffer%d3)) deallocate(buffer%d3)
   if (allocated(buffer%d4)) deallocate(buffer%d4)
   if (allocated(buffer%d5)) deallocate(buffer%d5)
!
   return
   end subroutine buffer_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module buffer
