!-------------------------------------------------------------------------------
   module mesh_control
!-------------------------------------------------------------------------------
!
!  abstract : mesh control module
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only : i4, l4, r4, r8
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   integer(i4), parameter :: ww = 1
   integer(i4), parameter :: ee = 2
   integer(i4), parameter :: ss = 3
   integer(i4), parameter :: nn = 4
   integer(i4), parameter :: ws = 5
   integer(i4), parameter :: es = 6
   integer(i4), parameter :: wn = 7
   integer(i4), parameter :: en = 8
!
   integer(i4), parameter :: axis_x   = 1 ! x-axis
   integer(i4), parameter :: axis_y   = 2 ! y-axis
   integer(i4), parameter :: axis_xy  = 3 ! xy-axis
   integer(i4), parameter :: axis_mxy = 4 !-x, y x,-y
   integer(i4), parameter :: axis_0   = 5 !(0, 0)
   character(len=1), dimension(4) :: dstr =(/'w', 'e', 's', 'n'/)
!
   public :: ww, ee, ss, nn, ws, es, wn, en
   public :: axis_x, axis_y, axis_xy, axis_mxy, axis_0
!
   public :: mesh_print
   public :: mesh_convert, mesh_dir_convert
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine mesh_convert(axis, mesh, cmesh)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                 intent(in   ) :: axis
   integer(i4), dimension(:,:), intent(in   ) :: mesh
   integer(i4), dimension(:,:), intent(  out) :: cmesh
! local variables
   integer(i4) :: n, i, j
!
   n = size(mesh,dim=1)
!
   if (n.ne.size(mesh,dim=2).or.n.ne.size(cmesh,dim=1).or.n.ne.size(cmesh,dim=2)) then
     print*,'check dimensions of mesh in mesh_convert...'
     stop
   endif
!
   do j = 1,n
     do i = 1,n
       if (axis.eq.0) then
         cmesh(i,j) = mesh(i,j)
       elseif (axis.eq.axis_x) then
         cmesh(i,j) = mesh(i,n-j+1)
       elseif (axis.eq.axis_y) then
         cmesh(i,j) = mesh(n-i+1,j)
       elseif (axis.eq.axis_xy) then
         cmesh(i,j) = mesh(j,i)
       elseif (axis.eq.axis_mxy) then
         cmesh(i,j) = mesh(n-j+1,n-i+1)
       elseif (axis.eq.axis_0) then
         cmesh(i,j) = mesh(n-i+1,n-j+1)
       else
         print*,'check axis in mesh_convert...',axis
         stop
       endif
     enddo
   enddo
!
   return
   end subroutine mesh_convert
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine mesh_dir_convert(axis, mesh, cmesh)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                 intent(in   ) :: axis
   integer(i4), dimension(:,:), intent(in   ) :: mesh
   integer(i4), dimension(:,:), intent(  out) :: cmesh
! local variables
   integer(i4) :: n, i, j
   integer(i4), dimension(:,:), allocatable :: tmp
!
   n = size(mesh,dim=1)
!
   if (n.ne.size(mesh,dim=2).or.n.ne.size(cmesh,dim=1).or.n.ne.size(cmesh,dim=2)) then
     print*,'check dimensions of mesh in mesh_dir_convert...'
     stop
   endif
!
   allocate(tmp(n,n))
   call mesh_convert(axis,mesh,tmp)
   cmesh(:,:) = tmp(:,:)
   if (axis.eq.0) then
   elseif (axis.eq.axis_x) then
     do j = 1,n
     do i = 1,n
       if (tmp(i,j).eq.ss) cmesh(i,j) = nn
       if (tmp(i,j).eq.nn) cmesh(i,j) = ss
     enddo
     enddo
   elseif (axis.eq.axis_y) then
     do j = 1,n
     do i = 1,n
       if (tmp(i,j).eq.ww) cmesh(i,j) = ee
       if (tmp(i,j).eq.ee) cmesh(i,j) = ww
     enddo
     enddo
   elseif (axis.eq.axis_xy) then
     do j = 1,n
     do i = 1,n
       if (tmp(i,j).eq.ww) cmesh(i,j) = ss
       if (tmp(i,j).eq.ee) cmesh(i,j) = nn
       if (tmp(i,j).eq.ss) cmesh(i,j) = ww
       if (tmp(i,j).eq.nn) cmesh(i,j) = ee
     enddo
     enddo
   elseif (axis.eq.axis_mxy) then
     do j = 1,n
     do i = 1,n
       if (tmp(i,j).eq.ww) cmesh(i,j) = nn
       if (tmp(i,j).eq.ee) cmesh(i,j) = ss
       if (tmp(i,j).eq.ss) cmesh(i,j) = ee
       if (tmp(i,j).eq.nn) cmesh(i,j) = ww
     enddo
     enddo
   elseif (axis.eq.axis_0) then
     do j = 1,n
     do i = 1,n
       if (tmp(i,j).eq.ww) cmesh(i,j) = ee
       if (tmp(i,j).eq.ee) cmesh(i,j) = ww
       if (tmp(i,j).eq.ss) cmesh(i,j) = nn
       if (tmp(i,j).eq.nn) cmesh(i,j) = ss
     enddo
     enddo
   else
     print*,'check axis in mesh_dir_convert...',axis
     stop
   endif
   deallocate(tmp)
!
   return
   end subroutine mesh_dir_convert
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine mesh_print(mesh, isdir)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), dimension(:,:), intent(in   ) :: mesh
   logical(l4), optional :: isdir
! local variables
   integer(i4) :: n, i, j
   logical(l4) :: l_isdir
   character(len=128) :: var_str
!
   n = size(mesh,dim=1)
!
   if (n.ne.size(mesh,dim=2)) then
     print*,'check dimensions of mesh in mesh_print...'
     stop
   endif
!
   if (present(isdir)) then
     l_isdir = isdir
   else
     l_isdir = .false.
   endif
!
   do j = n,1,-1
     call value2string(mesh(:,j),var_str,l_isdir)
     print*,var_str
   enddo
   print*,' '
!
   return
   end subroutine mesh_print
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine value2string(val, string, isdir)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), dimension(:), intent(in   ) :: val
   character(len=*),          intent(  out) :: string
   logical(l4) :: isdir
! local variables
   character(len=16) :: fmts
   integer(i4) :: n, i
!
   n = size(val)
!
   if (isdir) then
     string = ''
     do i = 1,n
       if (val(i).eq.-1) then
         write(string,'(a,x,a1) ') trim(string),'-'
       else
         write(string,'(a,x,a1) ') trim(string),dstr(val(i))
       endif
     enddo
   else
     write(fmts,'(a,i3,a) ') '(',n,'(i4,x)) '
     write(string,trim(fmts)) val
   endif
!
   return
   end subroutine value2string
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module mesh_control
