!-------------------------------------------------------------------------------
   module semi_lagrangian
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    201?-09-30  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds,    only: i4, l4, r4, r8
   use logger,   only: printlog, tostring
   use controls, only: np, ne, ntimelevs
   use element,  only: element_t, nups
   use element,  only: elem2up_pos, elem2up_state
   use element,  only: up2elem_pos, up2elem_state
   use piecewise_remap, only: piecewise_t, piecewise_initialize, piecewise_finalize
   use piecewise_remap, only: piecewise_set_x, piecewise_set, piecewise_get
   use spline_remap, only: spline_t, spline_initialize, spline_finalize
   use spline_remap, only: spline_set_x, spline_set, spline_get, spline_cyclic_interface
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   integer(i4) :: elev = 0
   integer(i4) :: llev = 2
!
   ! 1d domain
   type, public :: semilagrangian_t
     real(r8), dimension(:), allocatable :: grid
     real(r8), dimension(:), allocatable :: cell
     real(r8), dimension(:), allocatable :: v, t
     real(r8) :: cmin, cmax
     real(r8) :: dt, mass
     integer(i4) :: method
     logical(l4) :: mono
     integer(i4) :: bndry
   end type semilagrangian_t
!
   type(piecewise_t) :: pw
   type(spline_t)    :: sr
!
   public :: sl_initialize, sl_finalize, sl_gen_ctl_volume, sl_run, sl_get_mass
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine sl_initialize(sl, dt, mono, bndry)
!-------------------------------------------------------------------------------
   implicit none
!
   type(semilagrangian_t), intent(inout) :: sl
   real(r8),               intent(in   ) :: dt
   logical(l4),            intent(in   ) :: mono
   integer(i4),            intent(in   ) :: bndry
!
   allocate(sl%grid(nups))
   allocate(sl%cell(0:nups))
   allocate(sl%v(nups))
   allocate(sl%t(nups))
!
   sl%method = 1
   sl%dt     = dt
   sl%mono   = mono
   sl%bndry  = bndry
!
   return
   end subroutine sl_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine sl_finalize(sl)
!-------------------------------------------------------------------------------
   implicit none
!
   type(semilagrangian_t), intent(inout) :: sl
!
   deallocate(sl%grid)
   deallocate(sl%cell)
   deallocate(sl%v)
   deallocate(sl%t)
!
   return
   end subroutine sl_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine sl_gen_ctl_volume(sl, elem)
!-------------------------------------------------------------------------------
   implicit none
!
   type(semilagrangian_t),        intent(inout) :: sl
   type(element_t), dimension(:), intent(inout) :: elem
! local variables
   integer(i4) :: ie,ip
!
   do ie = 1,ne
     call sl_gen_ctl_volume_elem(np,elem(ie)%pos%x(:),elem(ie)%pos%cx(:))
     !call check_grid_cell(np,elem(ie)%pos%x(:),elem(ie)%pos%cx(:))
   enddo
!
   call elem2up_pos(elem,sl%grid,sl%cell)
   sl%cmin = sl%cell(0)
   sl%cmax = sl%cell(nups)
!
   return
   end subroutine sl_gen_ctl_volume
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine sl_gen_ctl_volume_elem(n, grid, cell)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(in   ) :: grid
   real(r8), dimension(0:n), intent(  out) :: cell
! local variables
   integer(i4) :: i
   logical(l4) :: leven
   real(r8) :: dx0, dxc
!
   if (n<2) then
     print*,'must be greater then 2 ...'
   endif
!
   if ((n/2)*2==n) then
     leven = .true.
   else
     leven = .false.
   endif
!
   if (leven) then
!
     dxc = 0.5_r8*(grid(n/2+1)-grid(n/2))
     ! center
     cell(n/2) = grid(n/2)+dxc
     !
     do i = n/2-1,0,-1
       cell(i) = grid(i+1)-(cell(i+1)-grid(i+1))
       cell(n-i) = grid(n-i)-(cell(n-i-1)-grid(n-i))
     enddo
!
   else
!
     if (n==2) then
       dx0 = 0.5_r8*(grid(2)-grid(1))
     else
       dx0 = 0.5_r8*(grid(2)-grid(1))**2.0_r8/(grid(3)-grid(2))
     endif
     cell(0) = grid(1)-dx0
     cell(1) = grid(1)+dx0
     cell(n-1) = grid(n)-dx0
     cell(n) = grid(n)+dx0
     !
     do i = 2,n/2
       !2-3
       cell(i) = grid(i)+(grid(i)-cell(i-1))
       !n-2,n-3
       !n-2      = grid(n-2) + (grid(n-2) - cell(n-3))
       cell(n-i) = grid(n-i+1)-(cell(n-i+1)-grid(n-i+1))
     enddo
!
   endif
!
   return
   end subroutine sl_gen_ctl_volume_elem
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_grid_cell(n, grid, cell)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(in   ) :: grid
   real(r8), dimension(0:n), intent(in   ) :: cell
! local variables
   integer(i4) :: i
   integer(i4), parameter :: width = 9 ! width of decimal number(include plus, minus, point),(must be less than 10)
!   integer(i4), parameter :: bel_dec_width = 3      ! width of below decimal number (must be less than or equal to width-3)
   integer(i4), parameter :: up_dec_width = 4 ! width of upper decimal number(must be greater than or equal to 3)
   character(len=width) :: str_g, str_c
   character(len=width*2+1) :: two_nums
   character(len =(n+2)*(width*2+1)) :: str_pos, str_val
   ! formats
   character(len=4) :: grid_idx_fmt
   character(len=9) :: cell_idx_fmt
   character(len=6) :: each_val_fmt
   character(len=7) :: two_nums_fmt
!
   ! formats
   write(grid_idx_fmt,'(a,i1,a) ') '(i',width,') '
   write(cell_idx_fmt,'(a,i1,a) ') '(i',width-2,'.1,a2) '
   write(each_val_fmt,'(a,i1,a,i1,a) ') '(f',width,'.',width-up_dec_width,') '
   write(two_nums_fmt,'(a) ') '(a,x,a) '
#if 0
   print*,'grid_idx_fmt(eg. 1) = ',grid_idx_fmt
   print*,'cell_idx_fmt(eg. 3/2) = ',cell_idx_fmt
   print*,'each_val_fmt(eg. 3.12) = ',each_val_fmt
   print*,'two_nums_fmt(eg. 3.12 2.22) = ',two_nums_fmt
   print*,'grid_idx_fmt(eg. 1) = ',grid_idx_fmt
   print*,'cell_idx_fmt(eg. 3/2) = ',cell_idx_fmt
   print*,'each_val_fmt(eg. 3.12) = ',each_val_fmt
   print*,'two_nums_fmt(eg. 3.12 2.22) = ',two_nums_fmt
#endif
!
   ! position
   write(str_c,cell_idx_fmt)(0*2+1),'/2'
   write(two_nums,two_nums_fmt) ' ',str_c
   str_pos = two_nums
   do i = 1,n
     write(str_g,grid_idx_fmt) i
     write(str_c,cell_idx_fmt)(i*2+1),'/2'
     write(two_nums,two_nums_fmt)(str_g),(str_c)
     str_pos = trim(str_pos)//' '//two_nums
   enddo
   ! value
   write(str_c,each_val_fmt) cell(0)
   write(two_nums,two_nums_fmt) ' ',str_c
   str_val = two_nums
   do i = 1,n
     write(str_g,each_val_fmt) grid(i)
     write(str_c,each_val_fmt) cell(i)
     write(two_nums,two_nums_fmt)(str_g),(str_c)
     str_val = trim(str_val)//' '//two_nums
   enddo
   ! check grid,cell
   call printlog(' ',elev)
   call printlog('START:check position of points...',elev,llev)
   call printlog('Positions of points...',elev,llev+1)
   call printlog(str_pos,elev,llev+1)
   call printlog(str_val,elev,llev+1)
!
   do i = 2,n
     if (grid(i)<=grid(i-1)) then
       call printlog('check grids...,index:'//tostring(i),-2,llev+1)
     endif
   enddo
!
   do i = 1,n
     if (cell(i)<=cell(i-1)) then
       call printlog('check cells...,index:'//tostring(i),-2,llev+1)
     endif
   enddo
!
   do i = 1,n/2
     if (cell(i)-cell(i-1)/= cell(n+1-i)-cell(n+1-i-1)) then
       call printlog('check delta cells...,index:'//tostring(i),-1,llev+1)
       print*,cell(i)-cell(i-1),cell(n+1-i)-cell(n+1-i-1)
     endif
   enddo
!
   call printlog('END:check position of points...',elev,llev)
   call printlog(' ',elev)
!
   return
   end subroutine check_grid_cell
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine sl_run(sl, elem)
!-------------------------------------------------------------------------------
   implicit none
!
   type(semilagrangian_t),        intent(inout) :: sl
   type(element_t), dimension(:), intent(inout) :: elem
! local variables
   integer(i4) :: i, ns, nl, nr
   real(r8), dimension(:), allocatable :: scell, tcell
   real(r8), dimension(:), allocatable :: src
!
   allocate(tcell(0:nups))
!
   ! elem%v -> sl%v
   call printlog('elem%v->sl%v',elev,llev)
   call elem2up_state(elem,sl%v,1)
!
   ! cell' = -v*cell
   call printlog("cell' =-v*cell",elev,llev)
   call back_trajectory(sl,tcell)
!
   ! elem%T -> sl%T
   call printlog('elem%T->sl%t',elev,llev)
   call elem2up_state(elem,sl%t,2)
!
   ! expand elements,include allocations of 'scell','src'
   call printlog('expand elements:',elev,llev)
   call expand_source(sl,elem,tcell,ns,nl,nr,scell,src)
!
   ! call run
   call printlog('PPM remapping:',elev,llev)
   if (sl%method.eq.1) then
     call piecewise_initialize(pw,ns,2,sl%mono,sl%bndry,(np-1)*nl,(np-1)*nr)
     call piecewise_set_x(pw,scell)
     call piecewise_set(pw,src)
     call piecewise_get(pw,nups,tcell,sl%t)
     call piecewise_finalize(pw)
   else
#if 0
     call spline_initialize(sr,ns,2,sl%mono,sl%bndry,(np-1)*nl,(np-1)*nr)
     call spline_set_x(sr,scell)
     call spline_set(sr,src)
     call spline_get(sr,nups,tcell,sl%t)
     call spline_finalize(sr)
#else
     call spline_cyclic_interface(nups,sl%cell,sl%t,nups,tcell,sl%t,2,sl%mono,sl%bndry)
#endif
   endif
!
   ! sl%T to elem%T
   call printlog('restore state',elev,llev)
   call up2elem_state(sl%t,elem,2)
!
   ! get mass
   sl%mass = sl_get_mass(nups,sl%cell,sl%t)
!
   deallocate(tcell)
   deallocate(scell)
   deallocate(src)
!
   return
   end subroutine sl_run
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine back_trajectory(sl, tcell)
!-------------------------------------------------------------------------------
   implicit none
!
   type(semilagrangian_t),      intent(in   ) :: sl
   real(r8), dimension(0:nups), intent(  out) :: tcell
! local variables
   integer(i4) :: i
!
   ! temporary (maybe need interpolation of velocity)
   tcell(0) = sl%cell(0)-sl%v(1)*sl%dt
   do i = 1,nups
     tcell(i) = sl%cell(i)-sl%v(i)*sl%dt
   enddo
!
   return
   end subroutine back_trajectory
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine expand_source(sl, elem, tcell, ns, neleml, nelemr, scell, src)
!-------------------------------------------------------------------------------
   implicit none
!
   type(semilagrangian_t),              intent(in   ) :: sl
   type(element_t), dimension(:),       intent(in   ) :: elem
   real(r8), dimension(0:nups),         intent(in   ) :: tcell
   integer(i4),                         intent(inout) :: ns, neleml, nelemr
   real(r8), dimension(:), allocatable, intent(  out) :: scell
   real(r8), dimension(:), allocatable, intent(  out) :: src
! local variables
   integer(i4) :: i, ie, idx, nexpl, nexpr
   real(r8) :: dxl, dxr
!
   dxl = sl%cell(0)-tcell(0)
   dxr = tcell(nups)-sl%cell(nups)
!
   neleml = 0 ! # of left element
   nelemr = 0 ! # of right element
   nexpl = 0 ! # of left points
   nexpr = 0 ! # of right points
   ns = nups
!
   if (dxl>0) then
     neleml = int(dxl/elem(1)%dx)
     if (mod(dxl,elem(1)%dx).ne.0) neleml = neleml+1
     nexpl =(np-1)*neleml
     ns = ns+nexpl
   endif
   if (dxr>0) then
     nelemr = int(dxr/elem(1)%dx)
     if (mod(dxr,elem(1)%dx).ne.0) nelemr = nelemr+1
     nexpr =(np-1)*nelemr
     ns = ns+nexpr
   endif

   if (neleml>ne.or.nelemr>ne) then
     call printlog('nelemL or nelemr was greater than ne...',-2,llev+1)
   endif
   call printlog('nelemL = '//tostring(neleml),elev,llev+1)
   call printlog('nelemR = '//tostring(nelemr),elev,llev+1)
   allocate(scell(0:ns))
   allocate(src(ns))

   ! left
   if (neleml>0) then
!
     ! copy source (left)
     do i = nexpl,nexpl+nups
       scell(i) = sl%cell(i-nexpl)
     enddo
     do i = nexpl+1,nexpl+nups
       src(i) = sl%t(i-nexpl)
     enddo
     ! expand source (left)
     do i = 0,nexpl-1
       !scell(i) = sl%cmin-(sl%cmax-sl%cell(nups-nexpl+i))
       scell(i) = sl%cmin-(sl%cell(nups-1)-sl%cell(nups-1-nexpl+i))
     enddo
     do i = 1,nexpl
       src(i) = sl%t(nups-nexpl-1+i)
     enddo
!
   endif
!
   ! right
   if (nelemr>0) then
!
     ! copy source (right)
     do i = 0,nups
       scell(i) = sl%cell(i)
     enddo
     do i = 1,nups
       src(i) = sl%t(i)
     enddo
     ! expand source (right)
     do i = 0,nexpr-1
       !scell(nups+1+i) = sl%cmax+(sl%cell(1+i)-sl%cmin)
       scell(nups+1+i) = sl%cmax+(sl%cell(1+i+1)-sl%cell(1))
     enddo
     do i = 1,nexpr
       src(nups+i) = sl%t(1+i)
     enddo
!
   endif
!
   !call print_points(nups+1,nexpl,nexpr,sl%cell,scell,.false.)
   !call print_points(nups,  nexpl,nexpr,sl%t,   src,  .true. )
!
   return
   end subroutine expand_source
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_points(nn, nl, nr, pts, exp_pts, isgrid)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                   intent(in   ) :: nn, nl, nr
   real(r8), dimension(nn),       intent(in   ) :: pts
   real(r8), dimension(nn+nl+nr), intent(in   ) :: exp_pts
   logical,                       intent(in   ) :: isgrid
! local variables
   integer(i4) :: ts, ps
   character(len=5) :: head
   character(len=32) :: prtfmt
   character(len=2048) :: string
!
   ts = 7 ! total size of floating point
   ps = 3
!
   if (isgrid) then
     write(head,'(a,i2.2,a) ') '(',ts/2,'X,'
   else
     head = '('
   endif
!
   ! orginal
   if (nl>0.and.nr<=0) then
   write(prtfmt,'(a,i3.3,a,i3.3,a,i2.2,a,i2.2,a) ') trim(head),ts*nl,'X,',nn,'F',ts,'.',ps,') '
   elseif (nl<=0.and.nr>0) then
   write(prtfmt,'(a,i3.3,a,i2.2,a,i2.2,a,i3.3,a) ') trim(head),nn+1,'F',ts,'.',ps,',',ts*nr,'X) '
   elseif (nl>0.and.nr>0) then
   write(prtfmt,'(a,i3.3,a,i3.3,a,i2.2,a,i2.2,a,i3.3,a) ') trim(head),ts*nl,'X,',nn,'F',ts,'.',ps,',',ts*nr,'X) '
   else
   write(prtfmt,'(a,i3.3,a,i2.2,a,i2.2,a) ') trim(head),nn,'F',ts,'.',ps,') '
   endif
   write(string,trim(prtfmt)) pts
   call printlog(string,elev,llev+1)
!
   ! expand
   write(prtfmt,'(a,i3.3,a,i2.2,a,i2.2,a) ') trim(head),nl+nn+nr,'F',ts,'.',ps,') '
   write(string,prtfmt) exp_pts
   call printlog(string,elev,llev+1)
!
   return
   end subroutine print_points
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function sl_get_mass(n, cell, values) result(mass)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(0:n), intent(in   ) :: cell
   real(r8), dimension(n),   intent(in   ) :: values
   real(r8) :: mass
! local variables
   integer(i4) :: i
   real(r8) :: dx
!
   mass = 0.0_r8
   do i = 1,n-1
     dx = cell(i)-cell(i-1)
     mass = mass+values(i)*dx
   enddo
!
   end function sl_get_mass
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module semi_lagrangian
!-------------------------------------------------------------------------------
