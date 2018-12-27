!-------------------------------------------------------------------------------
   module element
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    2015-09-30  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds, only:i4, l4, r4, r8
   use logger, only:printlog, tostring
   use controls, only:np, ne, ntimelevs
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   integer(i4) :: elev = 0
   integer(i4) :: llev = 2
!
   integer(i4), parameter, public :: which_x = 0
   integer(i4), parameter, public :: which_v = 1
   integer(i4), parameter, public :: which_t = 2
! state type
   type, public :: elem_state_t
     real(r8), dimension(:,:,:), allocatable :: v ! velocity(np, 1, ntimelev)
     real(r8), dimension(:,:),   allocatable :: t ! temperature(np, ntimelev)
   end type
! position type
   type, public :: cartesian_t
     real(r8), dimension(:), allocatable :: x ! np
     real(r8), dimension(:), allocatable :: cx ! np+1
   end type
! element
   type, public :: element_t
     type(elem_state_t) :: state
     type(cartesian_t)  :: pos
     real(r8) :: dx
   end type
! 1d domain
   integer(i4), public :: nups, nups_l
!
   public :: elem_initialize, elem_finalize, elem_print, elem_difference, elem_set_position, elem_set_state
   public :: elem2up_pos, elem2up_state, up2elem_pos, up2elem_state
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine elem_initialize(elem)
!-------------------------------------------------------------------------------
   implicit none
!
   type(element_t), dimension(:), allocatable, intent(inout) :: elem
! local variables
   integer(i4) :: ie, ip
   character(len=80) :: string
!
   allocate(elem(ne))
   do ie = 1,ne
     allocate(elem(ie)%pos%x(np))
     allocate(elem(ie)%pos%cx(0:np))
     allocate(elem(ie)%state%v(np,1,ntimelevs))
     allocate(elem(ie)%state%t(np,ntimelevs))
   enddo
!
   nups = ne*(np-1)+1
   write(string,'(a,i02,a,i2,a,i02) ') 'Element was initialized:ne = ',ne,',np = ',np,',nups = ',nups
   call printlog(string,elev,llev)
!
   return
   end subroutine elem_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine elem_finalize(elem)
!-------------------------------------------------------------------------------
   implicit none
!
   type(element_t), dimension(:), allocatable, intent(inout) :: elem
! local variables
   integer(i4) :: ie
!
   do ie = 1,ne
     deallocate(elem(ie)%pos%x)
     deallocate(elem(ie)%pos%cx)
     deallocate(elem(ie)%state%v)
     deallocate(elem(ie)%state%t)
   enddo
   deallocate(elem)
!
   return
   end subroutine elem_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine elem_print(elem)
!-------------------------------------------------------------------------------
   implicit none
!
   type(element_t), dimension(:), allocatable, intent(in   ) :: elem
! local variables
   integer(i4) :: ie, ip
   character(len=50) :: string
!
   do ie = 1,ne
     call printlog('elem:'//trim(tostring(ie)),elev,llev)
     do ip = 1,np
       write(string,'(a,i2.2,a,i01,a,f7.2)') 'elem(',ie,') :: x(',ip,') = ',elem(ie)%pos%x(ip)
       call printlog(trim(string),elev,llev+1)
     enddo
     do ip = 1,np
       write(string,'(a,i2.2,a,i01,a,f7.2)') 'elem(',ie,') :: v(',ip,') = ',elem(ie)%state%v(ip,1,1)
       call printlog(trim(string),elev,llev+1)
     enddo
     do ip = 1,np
       write(string,'(a,i2.2,a,i01,a,f7.2)') 'elem(',ie,') :: t(',ip,') = ',elem(ie)%state%t(ip,1)
       call printlog(trim(string),elev,llev+1)
     enddo
   enddo
!
   return
   end subroutine elem_print
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine elem_difference(eleml, elemr)
!-------------------------------------------------------------------------------
   implicit none
!
   type(element_t), dimension(:), allocatable, intent(in   ) :: eleml, elemr
! local variables
   integer(i4) :: ie, ip
   character(len=50) :: string
!
   call printlog('Checking positions ...',elev,llev)
   do ie = 1,ne
   do ip = 1,np
     if (eleml(ie)%pos%x(ip).ne.elemr(ie)%pos%x(ip)) then
       write(string,'(a,i2.2,a,i01,a,x,f5.2,x,f5.2) ') 'diff ie = ',ie,',ip = ',ip,':',eleml(ie)%pos%x(ip),elemr(ie)%pos%x(ip)
       call printlog(trim(string),elev,llev+1)
     endif
   enddo
   enddo
!
   call printlog('Checking v ...',elev,llev)
   do ie = 1,ne
   do ip = 1,np
     if (eleml(ie)%state%v(ip,1,1).ne.elemr(ie)%state%v(ip,1,1)) then
     write(string,'(a,i2.2,a,i01,a,x,f5.2,x,f5.2) ') 'diff ie = ',ie,',ip = ',ip,':',eleml(ie)%state%v(ip,1,1),elemr(ie)%state%v(ip,1,1)
     call printlog(trim(string),elev,llev+1)
   endif
   enddo
   enddo
!
   call printlog('Checking t ...',elev,llev)
   do ie = 1,ne
   do ip = 1,np
     if (eleml(ie)%state%t(ip,1).ne.elemr(ie)%state%t(ip,1)) then
     write(string,'(a,i2.2,a,i01,a,x,f5.2,x,f5.2) ') 'diff ie = ',ie,',ip = ',ip,':',eleml(ie)%state%t(ip,1),elemr(ie)%state%t(ip,1)
     call printlog(trim(string),elev,llev+1)
   endif
   enddo
   enddo
!
   return
   end subroutine elem_difference
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine elem_set_position(elem, xmin, xmax)
   use quadrature, only:quadrature_t, gausslobatto
!-------------------------------------------------------------------------------
   implicit none
!
   type(element_t), dimension(:), intent(inout) :: elem
   real(r8),                      intent(in   ) :: xmin, xmax
! local variables
   type(quadrature_t) :: qp
   integer(i4) :: ie, ip
   real(r8) :: xlen, xlen_e, fracx
!
   qp = gausslobatto(np)
!
   xlen=xmax-xmin
   xlen_e = xlen/dble(ne)
   fracx = xlen_e/(qp%points(np)-qp%points(1))
   do ie = 1,ne
     do ip = 1,np
       elem(ie)%pos%x(ip) = fracx*qp%points(ip)+xlen_e*(dble(ie)-0.5_r8)+xmin
     enddo
     elem(ie)%dx = elem(ie)%pos%x(np)-elem(ie)%pos%x(1)
   enddo
!
   call bndryexchange(elem,0)
!
   return
   end subroutine elem_set_position
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine elem_set_state(elem, iexp)
   use controls, only:velocity, amp, mu, sigma
!-------------------------------------------------------------------------------
   implicit none
!
   type(element_t), dimension(:), intent(inout) :: elem
   integer(i4),                   intent(in   ) :: iexp
! local variables
   integer(i4) :: ie, ip
   real(r8) :: x
!
   if (iexp==1) then
!
     do ie = 1,ne
       elem(ie)%state%v(:,1,1) = velocity
       do ip = 1,np
         x = elem(ie)%pos%x(ip)
         if (x>mu-2.0_r8*sigma.and.x<mu+2.0_r8*sigma) then
           elem(ie)%state%t(ip,1) = amp
         else
           elem(ie)%state%t(ip,1) = 0.0_r8
         endif
       enddo
     enddo
!
   elseif (iexp==2) then
!
     do ie = 1,ne
       elem(ie)%state%v(:,1,1) = velocity
       do ip = 1,np
         x = elem(ie)%pos%x(ip)
         if (x>mu-5.0_r8*sigma.and.x<mu+5.0_r8*sigma) then
           elem(ie)%state%t(ip,1) = amp*dexp(-1.0_r8*(x-mu)*(x-mu)/(2.0_r8*sigma*sigma))
         else
           elem(ie)%state%t(ip,1) = 0.0_r8
         endif
       enddo
     enddo
!
   else


   endif
!
!   call bndryexchange(elem,1)
!   call bndryexchange(elem,2)
!
   return
   end subroutine elem_set_state
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine bndryexchange(elem, which)
   use controls, only:velocity
!-------------------------------------------------------------------------------
   implicit none
!
   type(element_t), dimension(:), intent(inout) :: elem
   integer(i4),                   intent(in   ) :: which
! local variables
   integer(i4) :: ie, tie
!
   if (which==which_x) then
!
     do ie = 1,ne-1
       tie = ie+1
       elem(ie)%pos%x(np) = 0.5_r8*(elem(ie)%pos%x(np)+elem(tie)%pos%x(1))
       elem(tie)%pos%x(1) = elem(ie)%pos%x(np)
     enddo
!     elem( 1)%pos%x( 1) = 0.5_r8*(elem( 1)%pos%x( 1)+elem(ne)%pos%x(np))
!     elem(ne)%pos%x(np) = elem( 1)%pos%x( 1)
!
   elseif (which==which_v) then
!
     do ie = 1,ne-1
       tie = ie+1
       elem(ie)%state%v(np,1,1) = 0.5_r8*(elem(ie)%state%v(np,1,1)+elem(tie)%state%v(1,1,1))
       elem(tie)%state%v(1,1,1) = elem(ie)%state%v(np,1,1)
     enddo
     elem(1)%state%v(1,1,1) = 0.5_r8*(elem(1)%state%v(1,1,1)+elem(ne)%state%v(np,1,1))
     elem(ne)%state%v(np,1,1) = elem(1)%state%v(1,1,1)
!
   elseif (which==which_t) then
!
     do ie = 1,ne-1
       tie = ie+1
       elem(ie)%state%t(np,1) = 0.5_r8*(elem(ie)%state%t(np,1)+elem(tie)%state%t(1,1))
       elem(tie)%state%t(1,1) = elem(ie)%state%t(np,1)
     enddo
     elem(1)%state%t(1,1) = 0.5_r8*(elem(1)%state%t(1,1)+elem(ne)%state%t(np,1))
     elem(ne)%state%t(np,1) = elem(1)%state%t(1,1)
!
   endif
!
   return
   end subroutine bndryexchange
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine elem2up_pos(elem, grid, cell)
!-------------------------------------------------------------------------------
   implicit none
!
   type(element_t), dimension(:), intent(in   ) :: elem
   real(r8), dimension(nups),     intent(  out) :: grid
   real(r8), dimension(0:nups),   intent(  out) :: cell
! local variables
   integer(i4) :: ie, ip, idx
!
   idx = 0
   do ie = 1,ne
   do ip = 1,np-1
     idx = idx+1
     grid(idx) = elem(ie)%pos%x(ip)
   enddo
   enddo
   grid(nups) = elem(ne)%pos%x(np)
!
   idx = 0
   do ie = 1,ne
   do ip = 1,np-1
     idx = idx+1
     cell(idx) = elem(ie)%pos%cx(ip)
   enddo
   enddo
   cell(0) = elem(1)%pos%cx(0)
   cell(nups) = elem(ne)%pos%cx(np)
!
   return
   end subroutine elem2up_pos
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine elem2up_state(elem, ups, which)
!-------------------------------------------------------------------------------
   implicit none
!
   type(element_t), dimension(:), intent(in   ) :: elem
   real(r8), dimension(:),        intent(  out) :: ups
   integer(i4),                   intent(in   ) :: which
! local variables
   integer(i4) :: ie, ip, idx
!
   if (which==which_v) then
!
     idx = 0
     do ie = 1,ne
     do ip = 1,np-1
       idx = idx+1
       ups(idx) = elem(ie)%state%v(ip,1,1)
     enddo
     enddo
     ups(nups) = elem(ne)%state%v(np,1,1)
!
   elseif (which==which_t) then
!
     idx = 0
     do ie = 1,ne
     do ip = 1,np-1
       idx = idx+1
       ups(idx) = elem(ie)%state%t(ip,1)
     enddo
     enddo
     ups(nups) = elem(ne)%state%t(np,1)
!
   endif
!
   return
   end subroutine elem2up_state
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine up2elem_pos(grid, cell, elem)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(nups),     intent(in   ) :: grid
   real(r8), dimension(0:nups), intent(in   ) :: cell
   type(element_t), dimension(:), intent(inout) :: elem
! local variables
   integer(i4) :: ie, ip, idx
!
   idx = 0
   do ie = 1,ne
     do ip = 1,np-1
       idx = idx+1
       elem(ie)%pos%x(ip) = grid(idx)
     enddo
     elem(ie)%pos%x(np) = grid(idx+1)
   enddo
!
   idx = 0
   do ie = 1,ne
     do ip = 1,np-1
       idx = idx+1
       elem(ie)%pos%cx(ip) = cell(idx)
     enddo
   enddo
   elem(1)%pos%cx(0) = cell(0)
   elem(ne)%pos%cx(np) = cell(nups)
!
   return
   end subroutine up2elem_pos
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine up2elem_state(ups, elem, which)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:),        intent(in   ) :: ups
   type(element_t), dimension(:), intent(inout) :: elem
   integer(i4),                   intent(in   ) :: which
! local variables
   integer(i4) :: ie, ip, idx
!
   if (which==which_v) then
!
   idx = 0
   do ie = 1,ne
     do ip = 1,np-1
       idx = idx+1
       elem(ie)%state%v(ip,1,1) = ups(idx)
     enddo
     elem(ie)%state%v(np,1,1) = ups(idx+1)
   enddo
!
   elseif (which==which_t) then
!
   idx = 0
   do ie = 1,ne
     do ip = 1,np-1
       idx = idx+1
       elem(ie)%state%t(ip,1) = ups(idx)
     enddo
     elem(ie)%state%t(np,1) = ups(idx+1)
   enddo
!
   endif
!
   return
   end subroutine up2elem_state
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module element
