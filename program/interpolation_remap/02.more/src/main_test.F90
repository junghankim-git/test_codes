!-------------------------------------------------------------------------------
   program testmain
!
   use kinds,           only: i4, l4, r8
   use statistics,      only: get_lerror
   use io_module,       only: create_ncfile, write_ncfile
   use gen_test_func,   only: get_func, get_func_area
   use piecewise_remap, only: piecewise_t,                                     &
                              piecewise_initialize, piecewise_finalize,        &
                              piecewise_set_x,      piecewise_set,             &
                              piecewise_get_fun, piecewise_get_fun_integral,   &
                              piecewise_get, piecewise_interface
   use spline_remap,    only: spline_t,                                        &
                              spline_initialize, spline_finalize,              &
                              spline_set_x,      spline_set,                   &
                              spline_get_fun, spline_get_fun_integral,         &
                              spline_get, spline_interface, spline_cyclic_interface
   use spline_remap,    only: remap_psm, remap_psm_cyclic, test_remap_psm
!-------------------------------------------------------------------------------
   implicit none
   include 'mpif.h'
!
   type(piecewise_t) :: pr
   type(spline_t)    :: sr
   ! source
   integer(i4) :: n
   real(r8), dimension(:), allocatable :: sx, scx, sy, siy
   ! target: grid
   integer(i4) :: m
   real(r8), dimension(:), allocatable :: tx, tcx, ty, ty_ana, tiy
   ! target: analytic
   integer(i4) :: l
   real(r8), dimension(:), allocatable :: x_ana, ty_fun, ty_int
   real(r8), dimension(:), allocatable :: y_ana, y1_ana, y2_ana, iy_ana
   ! boundary derivative
   logical(l4) :: use_spline
   integer(i4) :: unitnum, testcase
   integer(i4) :: order, bndry ! setting
   logical     :: ismono
   integer(i4) :: i, j
   real(r8)    :: stime, etime, elapse(5), mass1, mass2, l2_t, l2_a, st, et
   character(len=256) :: filename
   logical     :: iscam = .true.
!
   print*,'###################################'
   print*,'######## start : main test ########'
   print*,' '
!
   unitnum  = 21
!
   use_spline = .true.
   order      = 2   ! 0:pcm, 1:plm, 2:ppm
   bndry      = 2   ! 0:mirror, 1:mirror+derivative, 2:periodic(need base)
   ismono     = .false.
   testcase   = 1
   n          = 11
   m          = 10
   l          = 500
   print*,'(info) read namelist (start)'
   call read_namelist('./inputs.nl')
   print*,'(info) read namelist (done)'
   print*,' '
!
   print*,'(info) allocate variables (start)'
   allocate(sx(n))       ! source x
   allocate(scx(0:n))    ! source x (control volume)
   allocate(sy(n))       ! source y (samples)
   allocate(siy(0:n))    ! source x (control volume)
   allocate(tx(m))       ! target x
   allocate(tcx(0:m))    ! target x (control volume)
   allocate(ty(m))       ! target y
   allocate(ty_ana(m))   ! target y
   allocate(tiy(0:m))    ! target x (control volume)
   allocate(x_ana(l))    ! x for analytic function
   allocate(ty_fun(l))   ! y from remap function
   allocate(ty_int(l))   ! integral value from remap function
   allocate(y_ana(l))    ! y     (analytic)
   allocate(y1_ana(l))   ! y'    (analytic)
   allocate(y2_ana(l))   ! y''   (analytic)
   allocate(iy_ana(l))   ! int y (analytic)
   print*,'(info) allocate variables (done)'
   print*,' '
!
   elapse = 0.0
   stime = mpi_wtime()
   if (.not.use_spline) then
     if (.not.iscam) call piecewise_initialize(pr,n,order,ismono,bndry)
   else
     call spline_initialize(sr,n,order,ismono,bndry)
   endif
   etime = mpi_wtime()
   elapse(1) = elapse(1)+etime-stime
!
   call set_experiment(testcase,n,sx,scx,sy,m,tx,tcx,ty_ana,l,x_ana,y_ana,y1_ana,y2_ana,iy_ana)
!
! source
!
   stime = mpi_wtime()
   if (.not.use_spline) then
     if (.not.iscam) call piecewise_set_x(pr,scx)
   else
     call spline_set_x(sr,scx)
   endif
   etime = mpi_wtime()
   elapse(2) = elapse(2)+etime-stime
!
   stime = mpi_wtime()
   if (.not.use_spline) then
     if (.not.iscam) call piecewise_set(pr,sy)
   else
     call spline_set(sr,sy)
   endif
   etime = mpi_wtime()
   elapse(3) = elapse(3)+etime-stime
!
!
   ! spline value
   if (.not.use_spline) then
     if (.not.iscam) call piecewise_get_fun(pr,l,x_ana,ty_fun,0)
   else
     call spline_get_fun(sr,l,x_ana,ty_fun,0)
   endif
   ! integration: spline
   if (.not.use_spline) then
     if (.not.iscam) call piecewise_get_fun_integral(pr,l,x_ana,ty_int)
   else
     call spline_get_fun_integral(sr,l,x_ana,ty_int)
   endif
   do i = 2,l
     ty_int(i) = ty_int(i-1)+ty_int(i)
   enddo
!
! target & remap
   stime = mpi_wtime()
   if (.not.use_spline) then
     if (.not.iscam) call piecewise_get(pr,m,tcx,ty)
   else
     call spline_get(sr,m,tcx,ty)
   endif
   etime = mpi_wtime()
   elapse(4) = elapse(4)+etime-stime
!
   stime = mpi_wtime()
   if (.not.use_spline) then
     if (.not.iscam) then
       call piecewise_finalize(pr)
     else
       !call remap_psm(n,scx,sy,m,tcx,ty,ismono,bndry,l=l,ftx=x_ana,fty=ty_fun)
st = mpi_wtime()
       call remap_psm(n,scx,sy,m,tcx,ty,ismono,bndry)
et = mpi_wtime()
print'(a,f11.8)',' (cam) elapse time (total)   = ',et-st
       !call remap_psm_cyclic(n,scx,sy,m,tcx,ty,ismono,l=l,ftx=x_ana,fty=ty_fun)
     endif
   else
     call spline_finalize(sr)
   endif
   etime = mpi_wtime()
   elapse(5) = elapse(5)+etime-stime
!
! mass
   siy(0) = 0.0
   do i = 1,n
     siy(i) = siy(i-1)+sy(i)*(scx(i)-scx(i-1))
   enddo
   tiy(0) = 0.0
   do i = 1,m
     tiy(i) = tiy(i-1)+ty(i)*(tcx(i)-tcx(i-1))
   enddo
!
   call create_ncfile('result/result.nc',use_spline,order,bndry,ismono,testcase,n,m,l)
   call write_ncfile(sx,scx,sy,siy,tx,tcx,ty,tiy,ty_ana,x_ana,y_ana,ty_fun,ty_int)
   call system('cp result/result.nc '//trim(filename))
!
! get info
   mass1 = siy(n)
   mass2 = tiy(m)
   l2_t  = get_lerror(m,ty_ana,ty,'L2')
   l2_a  = get_lerror(l,y_ana,ty_fun,'L2')
!
   print*,'(info) deallocate variables (start)'
   deallocate(sx)
   deallocate(scx)
   deallocate(sy)
   deallocate(siy)
   deallocate(tx)
   deallocate(tcx)
   deallocate(ty)
   deallocate(ty_ana)
   deallocate(tiy)
   deallocate(x_ana)
   deallocate(ty_fun)
   deallocate(ty_int)
   deallocate(y_ana)
   deallocate(y1_ana)
   deallocate(y2_ana)
   deallocate(iy_ana)
   print*,'(info) deallocate variables (done)'
!
   print*,'-------------------------------------------------------'
   print*,'(info) mass (initial)        = ',mass1
   print*,'(info) mass (remap)          = ',mass2
   print*,'-------------------------------------------------------'
   !print*,'(info) L2 error (target)     = ',l2_t
   !print*,'(info) L2 error (function)   = ',l2_a
   print'(a,f11.8)',' (info) L2 error (target)     = ',l2_t
   print'(a,f11.8)',' (info) L2 error (function)   = ',l2_a
   print*,'-------------------------------------------------------'
   !print*,'(info) elapse time (total)   = ',sum(elapse)
   !print*,'(info) elapse time (initial) = ',elapse(1)
   !print*,'(info) elapse time (set)     = ',elapse(2)
   !print*,'(info) elapse time (get)     = ',elapse(3)
   !print*,'(info) elapse time (final)   = ',elapse(4)
   print'(a,f11.8)',' (info) elapse time (total)   = ',sum(elapse)
   print'(a,f11.8)',' (info) elapse time (initial) = ',elapse(1)
   print'(a,f11.8)',' (info) elapse time (set_x)   = ',elapse(2)
   print'(a,f11.8)',' (info) elapse time (set)     = ',elapse(3)
   print'(a,f11.8)',' (info) elapse time (get)     = ',elapse(4)
   print'(a,f11.8)',' (info) elapse time (final)   = ',elapse(5)
   print*,'-------------------------------------------------------'
!
   print*,' '
   print*,'######## end : main test ########'
   print*,'###################################'
   print*,' '
!   call test_remap_psm()
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_namelist(nlfilename)
!-------------------------------------------------------------------------------
   implicit none
   character(len=*), optional, intent(in   ) :: nlfilename
! local
   integer(i4) :: spline, mono
   namelist/inputs/use_spline, order, bndry, ismono, testcase, n, m, l
!
   if (present(nlfilename)) then
     open(21,file=trim(nlfilename),status='old')
     read(21,nml=inputs)
     close(21)
   else
     read(*,nml = inputs)
   endif
!
   spline = 0; mono = 0
   if (use_spline) spline = 1
   if (ismono) mono = 1
   write(filename,'(a,i1,a,i1,a,i1,a,i1,a,i1,a)')                              &
           'result/s',spline,'_o',order,'_b',bndry,'_m',mono,'_t',testcase,'.nc'
!
   print*,' use spline = ',use_spline
   print*,' order      = ',order
   print*,' bndry      = ',bndry
   print*,' ismono     = ',ismono
   print*,' testcase   = ',testcase
   print*,' n          = ',n
   print*,' m          = ',m
   print*,' l          = ',l
   print*,' file       = ',filename
!
   return
   end subroutine read_namelist
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_experiment(iexp,n,sx,scx,sy,m,tx,tcx,ty_ana,l,x_ana,y_ana,y1_ana,y2_ana,iy_ana)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: iexp, n, m, l
   real(r8), dimension(n),   intent(inout) :: sx, sy
   real(r8), dimension(0:n), intent(inout) :: scx
   real(r8), dimension(m),   intent(inout) :: tx, ty_ana
   real(r8), dimension(0:m), intent(inout) :: tcx
   real(r8), dimension(l),   intent(inout) :: x_ana, y_ana, y1_ana, y2_ana, iy_ana
   integer(i4) :: i
   real(r8) :: xmin, xmax, dsx, dx_ana
   real(r8) :: a, mu, tmp(2)
!
   xmin = 0.0_r8
   xmax = 1.0_r8
!   if (n==5) then
!     xmin = 0.5_r8
!     xmax = 5.5_r8
!   endif
!
   ! set x in source
   dsx = (xmax-xmin)/dble(n)
   !mu = 0.5_r8*(xmax+xmin)+dsx
   mu = 0.5_r8*(xmax+xmin)
   do i = 0,n
     scx(i) = dsx*dble(i)+xmin
   enddo
   do i = 1,n
     sx(i) = 0.5_r8*(scx(i-1)+scx(i))
   enddo
!
   ! set x in target
   dx_ana =(xmax-xmin)/dble(m)
   do i = 0,m
     tcx(i) = dx_ana*dble(i)+xmin
   enddo
   do i = 1,m
     tx(i) = 0.5*(tcx(i-1)+tcx(i))
   enddo
!
   ! set x in target (for checking functions)
   dx_ana =(xmax-xmin)/dble(l-1)
   do i = 1,l
     x_ana(i) = dx_ana*dble(i-1)+xmin
   enddo
   x_ana(l) = x_ana(l)-1.d-16
!
   call get_func_area(iexp,n,scx,sy,xmin,xmax,mu)
   call get_func_area(iexp,m,tcx,ty_ana,xmin,xmax,mu)
!   if (n==5) then
!     sy = (/1.0,2.0,1.0,3.0,2.0/)
!   endif
!
   call get_func(iexp,l,x_ana,y_ana,xmin,xmax,mu,0)
   call get_func(iexp,l,x_ana,y1_ana,xmin,xmax,mu,1)
   call get_func(iexp,l,x_ana,y2_ana,xmin,xmax,mu,2)
   call get_func(iexp,l,x_ana,iy_ana,xmin,xmax,mu,-1)
!
   return
   end subroutine set_experiment
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_mass(n, cx, y) result(mass)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(0:n), intent(in   ) :: cx
   real(r8), dimension(  n), intent(in   ) :: y
   real(r8)                                :: mass
! local
   integer(i4) :: i
!
   mass = 0.0
   do i = 1,n
     mass = mass+y(i)*(cx(i)-cx(i-1))
   enddo
!
   end function get_mass
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program
!-------------------------------------------------------------------------------
