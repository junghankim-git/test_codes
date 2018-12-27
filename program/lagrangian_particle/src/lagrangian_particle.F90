!-------------------------------------------------------------------------------
   module lagrangian_particle
!-------------------------------------------------------------------------------
!
!  abstract : 
!
!  history log :
!    2018-12-12  junghan kim    initial setup
!
!  structure :
!
!-------------------------------------------------------------------------------
   use kinds,   only: i4, i8, r4, r8, l4
   use netcdf
   private
!
   real(r8), parameter ::  g = -9.80665_r8
   real(r8), parameter :: pi = acos(-1.0_r8)
   logical(l4)         :: debug = .false.
!
   ! ivis  : 0(no vis.), 
   ! uptype: 0(euler), 1(leapfrog), 2(rk2), 3(rk3)
   ! bndry : 0(no bndry), 1(periodic), 2(reflection)
   type lagrangian_t
     ! meta
     integer(i4) :: ndims, npart, ivis, uptype
     real(r8)    :: dt 
     integer(i4) :: nstep
     integer(i4) :: bndry(3) ! 1:left&right, 2:top&down
     integer(i4) :: istep
     ! values
     real(r8), dimension(:),   allocatable :: xmin, xmax
     real(r8), dimension(:,:), allocatable :: x0, v0, a0
     real(r8), dimension(:,:), allocatable :: x, v
     real(r8), dimension(:,:), allocatable :: ax, av
     ! file
     integer(i4) :: fid, xid, vid
   end type lagrangian_t
!
   public :: lagrangian_t, initialize, run, finalize
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine initialize(lag, npart, ivis, uptype, dt, nstep, xmin, xmax, ymin, ymax, bndry)
   implicit none
   type(lagrangian_t), intent(inout) :: lag
   integer(i4)       , intent(in   ) :: npart, ivis, uptype, nstep, bndry(2)
   real(r8)          , intent(in   ) :: dt, xmin, xmax, ymin, ymax
! local variables
   integer(i4) :: ndims
   real(r8)    :: v0, vang
   integer(i4) :: ip, id
!
   ndims      = 2
   lag%ndims  = ndims
   lag%npart  = npart
   lag%ivis   = ivis
   lag%uptype = uptype
   lag%dt     = dt
   lag%nstep  = nstep
   lag%bndry(1:2) = bndry
!
   allocate(lag%xmin(ndims))
   allocate(lag%xmax(ndims))
   allocate(lag%x0(ndims,npart))
   allocate(lag%v0(ndims,npart))
   allocate(lag%a0(ndims,npart))
   allocate(lag%x(ndims,npart))
   allocate(lag%v(ndims,npart))
   allocate(lag%ax(nstep+1,ndims))
   allocate(lag%av(nstep+1,ndims))
!
   lag%xmin(1) = xmin
   lag%xmin(2) = ymin
   lag%xmax(1) = xmax
   lag%xmax(2) = ymax
!
   v0 = 10.0_r8
   if (npart.eq.1) then
     lag%x0(1,1) = -50.0_r8
     lag%x0(2,1) =  30.0_r8
     lag%v0(1,1) =       v0
     lag%v0(2,1) =   0.0_r8
     lag%a0(1,1) =   0.0_r8
     lag%a0(2,1) =        g
   elseif (npart.eq.2) then
     lag%x0(1,1) = -30.0_r8
     lag%x0(2,1) =  -1.0_r8
     lag%v0(1,1) =       v0
     lag%v0(2,1) =   0.0_r8
     lag%a0(1,1) =   0.0_r8
     lag%a0(2,1) =   0.0_r8
     !
     lag%x0(1,2) = +30.0_r8
     lag%x0(2,2) =  +1.0_r8
     lag%v0(1,2) =      -v0
     lag%v0(2,2) =   0.0_r8
     lag%a0(1,2) =   0.0_r8
     lag%a0(2,2) =   0.0_r8
   elseif (npart.eq.3) then
     lag%x0(1,1) = -30.0_r8
     lag%x0(2,1) =   9.0_r8
     lag%v0(1,1) =       v0
     lag%v0(2,1) =   0.0_r8
     lag%a0(1,1) =   0.0_r8
     lag%a0(2,1) =   0.0_r8
     !
     lag%x0(1,2) = +30.0_r8
     lag%x0(2,2) =  11.0_r8
     lag%v0(1,2) =      -v0
     lag%v0(2,2) =   0.0_r8
     lag%a0(1,2) =   0.0_r8
     lag%a0(2,2) =   0.0_r8
     !
     lag%x0(1,3) =   0.0_r8
     lag%x0(2,3) =  30.0_r8
     lag%v0(1,3) =   0.0_r8
     lag%v0(2,3) =      -v0
     lag%a0(1,3) =   0.0_r8
     lag%a0(2,3) =   0.0_r8
   else
     call random_number(lag%x0(:,:))
     do ip = 1,npart
       do id = 1,ndims
         lag%x0(id,ip) = (lag%xmax(id)-lag%xmin(id))*lag%x0(id,ip)+lag%xmin(id)
       enddo
       call random_number(vang)
       vang = 2.0_r8*pi*vang
       lag%v0(1,ip) = v0*cos(vang)
       lag%v0(2,ip) = v0*sin(vang)
       !
       lag%a0(1,ip) = 0.0_r8
       lag%a0(2,ip) = 0.0_r8
     enddo
   endif
   lag%x  = lag%x0
   lag%v  = lag%v0
   if (debug) call print_status(lag)
!
   lag%ax = 0.0_r8
   lag%av = 0.0_r8
   call analytic_solution(lag)
!
   lag%istep = 0
   lag%istep = lag%istep+1
!
   call create_file(lag)
!
   return
   end subroutine initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine finalize(lag)
   implicit none
   type(lagrangian_t), intent(inout) :: lag
! local variables
!
   call close_file(lag)
!
   deallocate(lag%xmin)
   deallocate(lag%xmax)
   deallocate(lag%x0)
   deallocate(lag%v0)
   deallocate(lag%a0)
   deallocate(lag%x)
   deallocate(lag%v)
   deallocate(lag%ax)
   deallocate(lag%av)
!
   return
   end subroutine finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_status(lag)
   implicit none
   type(lagrangian_t), intent(in   ) :: lag
! local variables
   integer(i4) :: ndims, npart, id, ip
!
   ndims = lag%ndims
   npart = lag%npart
   do ip = 1,npart
     print '(a,i2,x,f10.4,x,f10.4,a,f10.4,x,f10.4)',                           &
                'log: ',ip,lag%x(1,ip),lag%x(2,ip),', ',lag%v(1,ip),lag%v(2,ip)
     !print *,'log: ',ip,lag%x(1,ip),lag%x(2,ip),', ',lag%v(1,ip),lag%v(2,ip)
   enddo
   print *,' '
!
   return
   end subroutine print_status
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine run(lag)
   implicit none
   type(lagrangian_t), intent(inout) :: lag
! local variables
   integer(i4) :: ndims, npart, is
   real(r8), dimension(:,:), allocatable :: xrhs, vrhs
!
   ndims = lag%ndims
   npart = lag%npart
!
   allocate(xrhs(ndims,npart))
   allocate(vrhs(ndims,npart))
!
   do is = 1,lag%nstep
     print *, 'into istep: ',is
     call compute_rhs(lag%ivis,ndims,npart,lag%x,lag%v,lag%a0,xrhs,vrhs)
     call apply_rhs(ndims,npart,lag%dt,lag%x,lag%v,xrhs,vrhs,lag%x,lag%v)
     call apply_bndry(lag%bndry,lag%uptype,ndims,npart,lag%xmin,lag%xmax,lag%x,lag%v)
     lag%istep = lag%istep+1
     call write_xv_in_file(lag)
     if (debug) call print_status(lag)
   enddo
!
   deallocate(xrhs)
   deallocate(vrhs)
!
   return
   end subroutine run
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine compute_rhs(ivis, ndims, npart, x, v, a, xrhs, vrhs)
   implicit none
   integer(i4),                      intent(in   ) :: ivis, ndims, npart
   real(r8), dimension(ndims,npart), intent(in   ) :: x, v, a
   real(r8), dimension(ndims,npart), intent(  out) :: xrhs, vrhs
! local variables
   integer(i4) :: ip, id
   real(r8), dimension(:,:), allocatable :: vis, inter, aa
!
   allocate(aa(ndims,npart))
   aa(:,:) = a(:,:)
   !
   if (npart.gt.1) then
     allocate(inter(ndims,npart))
     inter = interaction(ndims,npart,x,v)
     aa(:,:) = aa(:,:)+inter(:,:)
   endif
   !
   if (ivis.gt.0) then
     allocate(vis(ndims,npart))
     vis = get_viscosity(ndims,npart,v)
     aa(:,:) = aa(:,:)+vis(:,:)
   endif
!
   do ip = 1,npart
     do id = 1,ndims
       xrhs(id,ip) = v(id,ip)
       vrhs(id,ip) = aa(id,ip)
     enddo
   enddo
!
   deallocate(aa)
   if (npart.gt.1) then
     deallocate(inter)
   endif
   if (ivis.gt.0) then
     deallocate(vis)
   endif
!
   return
   end subroutine compute_rhs
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine apply_rhs(ndims, npart, dt, x0, v0, xrhs, vrhs, x, v)
   implicit none
   integer(i4),                      intent(in   ) :: ndims, npart
   real(r8),                         intent(in   ) :: dt
   real(r8), dimension(ndims,npart), intent(in   ) :: x0, v0
   real(r8), dimension(ndims,npart), intent(in   ) :: xrhs, vrhs
   real(r8), dimension(ndims,npart), intent(  out) :: x, v
! local variables
   integer(i4) :: ip, id
!
   do ip = 1,npart
     do id = 1,ndims
       x(id,ip) = x0(id,ip)+dt*xrhs(id,ip)
       v(id,ip) = v0(id,ip)+dt*vrhs(id,ip)
     enddo
   enddo
!
   return
   end subroutine apply_rhs
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine apply_bndry(bndry, uptype, ndims, npart, xmin, xmax, x, v)
   implicit none
   integer(i4), dimension(ndims),    intent(in   ) :: bndry
   integer(i4),                      intent(in   ) :: uptype
   integer(i4),                      intent(in   ) :: ndims, npart
   real(r8), dimension(ndims),       intent(in   ) :: xmin, xmax
   real(r8), dimension(ndims,npart), intent(inout) :: x, v
! local variables
   integer(i4) :: ip, id
!
   do ip = 1,npart
     do id = 1,ndims
       if     (x(id,ip)<xmin(id)) then
         if     (bndry(id).eq.1) then
           x(id,ip) = xmax(id)-(xmin(id)-x(id,ip))
         elseif (bndry(id).eq.2) then
           x(id,ip) = 2.0*xmin(id)-x(id,ip)
           v(id,ip) = -v(id,ip)
         endif
       elseif (x(id,ip)>xmax(id)) then
         if     (bndry(id).eq.1) then
           x(id,ip) = xmin(id)+(x(id,ip)-xmax(id))
         elseif (bndry(id).eq.2) then
           x(id,ip) = 2.0*xmax(id)-x(id,ip)
           v(id,ip) = -v(id,ip)
         endif
       endif
     enddo
   enddo
!
   return
   end subroutine apply_bndry
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_magnitude(ndims, vec) result(mag)
   implicit none
   integer(i4),                intent(in   ) :: ndims
   real(r8), dimension(ndims), intent(in   ) :: vec
   real(r8)                                  :: mag
! local variables
   integer(i4) :: id
!
   mag = 0.0_r8
   do id = 1,ndims
     mag = mag+vec(id)*vec(id)
   enddo
   mag = mag**(1.0_r8/real(ndims,r8))
!
   return
   end function get_magnitude
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_viscosity(ndims, npart, v) result(vis)
   implicit none
   integer(i4),                      intent(in   ) :: ndims, npart
   real(r8), dimension(ndims,npart), intent(in   ) :: v
   real(r8), dimension(ndims,npart)                :: vis
! local variables
   integer(i4) :: ip, id
   real(r8)    :: mu, mag
!
   mu = 0.002
!
   do ip = 1,npart
     mag = get_magnitude(ndims,v)
     do id = 1,ndims
       vis(id,ip) = -mu*mag*v(id,ip)
     enddo
   enddo
!
   return
   end function get_viscosity
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function interaction(ndims, npart, x, v) result(inter)
   implicit none
   integer(i4),                      intent(in   ) :: ndims, npart
   real(r8), dimension(ndims,npart), intent(in   ) :: x, v
   real(r8), dimension(ndims,npart)                :: inter
! local variables
   integer(i4) :: ip, jp, id
   real(r8)    :: mu, dis, dir(ndims)
!
   if (npart.le.3) then
     mu = 50.0_r8
   else
     mu = 50.0_r8
   endif
!
   do ip = 1,npart
     inter(:,ip) = 0.0_r8
     do jp = 1,npart
       if (ip.ne.jp) then
         do id = 1,ndims
           dir(id) = x(id,ip)-x(id,jp)
         enddo
         dis = get_magnitude(ndims,dir)
         do id = 1,ndims
           inter(id,ip) = inter(id,ip)+mu*dir(id)/(dis**2.0_r8)
         enddo
       endif
     enddo
   enddo
!
   return
   end function interaction
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine const_acc(ndims, t, x0, v0, a0, x, v)
   implicit none
   integer(i4)               , intent(in   ) :: ndims
   real(r8)                  , intent(in   ) :: t
   real(r8), dimension(ndims), intent(in   ) :: x0, v0, a0
   real(r8), dimension(ndims), intent(  out) :: x, v
! local variables
   integer(i4) :: id
!
   do id = 1,ndims
     x(id) = x0(id)+t*v0(id)+0.5_r8*a0(id)*t*t
     v(id) = v0(id)+t*a0(id)
   enddo
!
   return
   end subroutine const_acc
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine analytic_solution(lag)
   implicit none
   type(lagrangian_t), intent(inout) :: lag
! local variables
   integer(i4) :: ndims, nstep, it, id
   real(r8)    :: tt, tb, xb, tr
   real(r8), dimension(:)  , allocatable :: xr, vr
   real(r8), dimension(:,:), allocatable :: ax, av
!
   ndims = lag%ndims
   nstep = lag%nstep
!
   tb = sqrt(-2.0_r8*(lag%x0(2,1)-lag%xmin(2))/g)
   xb = lag%xmin(1)+lag%v0(1,1)*tb
!
   allocate(xr(ndims))
   allocate(vr(ndims))
   allocate(ax(ndims,nstep+1))
   allocate(av(ndims,nstep+1))
!
   do it = 1,nstep+1
     tt = real(it-1)*lag%dt
     if (tt.le.tb) then
       call const_acc(ndims,tt,lag%x0(:,1),lag%v0(:,1),lag%a0(:,1),ax(:,it),av(:,it))
     else
       call const_acc(ndims,tb,lag%x0(:,1),lag%v0(:,1),lag%a0(:,1),ax(:,it),av(:,it))
       tr    = tt-tb
       xr(:) =  ax(:,it)
       vr(1) =  av(1,it)
       vr(2) = -av(2,it)
       call const_acc(ndims,tr,xr,vr,lag%a0(:,1),ax(:,it),av(:,it))
     endif
     lag%ax(it,:) = ax(:,it)
     lag%av(it,:) = av(:,it)
   enddo
!
   deallocate(xr)
   deallocate(vr)
   deallocate(ax)
   deallocate(av)
!
   return
   end subroutine analytic_solution
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine create_file(lag)
   implicit none
   type(lagrangian_t), intent(inout) :: lag
! local variables
   integer(i4) :: err, did1, did2, did3, didu, vid_xmin, vid_xmax, vid_ax, vid_av
   ! did1(ndims), did2(npart), did3(nstep+1), didu(time)
!
   err = nf90_create('result.nc',ior(nf90_clobber,nf90_64bit_offset),lag%fid)
   call nc_check(err,'nc: create')
   err = nf90_put_att(lag%fid,nf90_global,'ivis',lag%ivis)
   call nc_check(err,'nc: att(ivis)')
   err = nf90_put_att(lag%fid,nf90_global,'uptype',lag%uptype)
   call nc_check(err,'nc: att(uptype)')
   err = nf90_put_att(lag%fid,nf90_global,'bndry',lag%bndry)
   call nc_check(err,'nc: att(bndry)')
   err = nf90_put_att(lag%fid,nf90_global,'nstep',lag%nstep)
   call nc_check(err,'nc: att(nstep)')
   err = nf90_put_att(lag%fid,nf90_global,'dt',lag%dt)
   call nc_check(err,'nc: att(dt)')
   err = nf90_def_dim(lag%fid,'ndims',lag%ndims,did1)
   call nc_check(err,'nc: dim(ndims)')
   err = nf90_def_dim(lag%fid,'npart',lag%npart,did2)
   call nc_check(err,'nc: dim(npart)')
   err = nf90_def_dim(lag%fid,'nstep1',lag%nstep+1,did3)
   call nc_check(err,'nc: dim(nstep+1)')
   err = nf90_def_dim(lag%fid,'time',nf90_unlimited,didu)
   call nc_check(err,'nc: dim(time)')
   err = nf90_def_var(lag%fid,'xmin',nf90_real8,(/did1/),vid_xmin)
   call nc_check(err,'nc: var(xmin)')
   err = nf90_def_var(lag%fid,'xmax',nf90_real8,(/did1/),vid_xmax)
   call nc_check(err,'nc: var(xmax)')
   err = nf90_def_var(lag%fid,'ax',nf90_real8,(/did3,did1/),vid_ax)
   call nc_check(err,'nc: var(ax)')
   err = nf90_def_var(lag%fid,'av',nf90_real8,(/did3,did1/),vid_av)
   call nc_check(err,'nc: var(av)')
   err = nf90_def_var(lag%fid,'x',nf90_real8,(/did1,did2,didu/),lag%xid)
   call nc_check(err,'nc: var(x)')
   err = nf90_def_var(lag%fid,'v',nf90_real8,(/did1,did2,didu/),lag%vid)
   call nc_check(err,'nc: var(v)')
   err = nf90_enddef(lag%fid)
   call nc_check(err,'nc: enddef')
   err = nf90_put_var(lag%fid,vid_xmin,lag%xmin)
   call nc_check(err,'nc: put(xmin)')
   err = nf90_put_var(lag%fid,vid_xmax,lag%xmax)
   call nc_check(err,'nc: put(xmax)')
   err = nf90_put_var(lag%fid,vid_ax,lag%ax)
   call nc_check(err,'nc: put(ax)')
   err = nf90_put_var(lag%fid,vid_av,lag%av)
   call nc_check(err,'nc: put(av)')
   call write_xv_in_file(lag)
!
   return
   end subroutine create_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_xv_in_file(lag)
   implicit none
   type(lagrangian_t), intent(inout) :: lag
! local variables
   integer(i4) :: err
   integer(i4) :: sta(3), cnt(3)
!
   sta(1) = 1; cnt(1) = lag%ndims
   sta(2) = 1; cnt(2) = lag%npart
   sta(3) = lag%istep; cnt(3) = 1
!
   err = nf90_put_var(lag%fid,lag%xid,lag%x,start=sta,count=cnt)
   call nc_check(err,'nc: put(x)')
   err = nf90_put_var(lag%fid,lag%vid,lag%v,start=sta,count=cnt)
   call nc_check(err,'nc: put(v)')
!
   return
   end subroutine write_xv_in_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine close_file(lag)
   implicit none
   type(lagrangian_t), intent(inout) :: lag
! local variables
   integer(i4) :: err
!
   err = nf90_close(lag%fid)
   call nc_check(err,'nc: close')
!
   return
   end subroutine close_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_check(err,message)
   implicit none
   integer(i4),                intent(in   ) :: err
   character(len=*), optional, intent(in   ) :: message
! local variables
!
   if (err.ne.nf90_noerr) then
     if (present(message)) print*,trim(message)//': '//trim(nf90_strerror(err))
     stop
   endif
!
   end subroutine nc_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module lagrangian_particle
!-------------------------------------------------------------------------------
