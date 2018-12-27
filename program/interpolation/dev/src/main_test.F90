!-------------------------------------------------------------------------------
   program testmain
!-------------------------------------------------------------------------------
   use kinds,    only: i4, l4, r8
!   use io_ascii, only: open_file, close_file, write_ascii, write_ascii_lines
   use quadrature, only : quadrature_t, gausslobatto
   use gen_test_func, only: get_func,                                          &
                            get_poly1, get_poly2, get_poly3,                   &
                            get_harmonic, get_gaussian, get_random
   use lagrange_interp, only: lagrange_t, lagrange_initialize,                 &
                              lagrange_finalize, lagrange_set_x, lagrange_set, &
                              lagrange_get, lagrange_integral
   use newton_interp, only: newton_t, newton_initialize,                       &
                            newton_finalize, newton_set_x, newton_set,         &
                            newton_get
   use spline_interp, only: spline_t, spline_initialize, spline_finalize,      &
                            spline_set_x, spline_set, spline_get, spline_integral
   use io_module,     only: create_ncfile
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), parameter :: pi = 3.14159265358979323846264_r8
   integer(i4), parameter :: nfun = 4, norder = 4
   type(lagrange_t) :: lp
   type(newton_t)   :: np
   type(spline_t)   :: sp_l
   type(spline_t)   :: sp_q
   type(spline_t)   :: sp_c
! source
   integer(i4) :: n = 10
   real(r8), dimension(:), allocatable :: sx, sy
   real(r8) :: sdx
! target
   integer(i4) :: m = 500
   real(r8), dimension(:),   allocatable :: tx, ty_a, typ_a, typp_a, tint_a
   real(r8), dimension(:,:), allocatable :: ty, typ, typp, tint
!
   ! boundary derivative
   real(r8) :: yp1, ypp1, yppn
!
   integer(i4) :: unitnum, testcase
   logical     :: isgll
   integer(i4) :: i, j
!
   print*,'######## start interpolation ########'
   print*,' '
!
   unitnum = 21
   testcase = 1
   yp1 = 0.0_r8
   ypp1 = 0.0_r8
   yppn = 0.0_r8
   n = 10
   m = 500
   call read_namelist('./inputs.nl')
   allocate(sx(n))
   allocate(sy(n))
   allocate(tx(m))
   allocate(ty_a(m))
   allocate(typ_a(m))
   allocate(typp_a(m))
   allocate(tint_a(m))
   allocate(ty(m,nfun))
   allocate(typ(m,nfun))
   allocate(typp(m,nfun))
   allocate(tint(m,nfun))
!
   call set_experiment(testcase,n,sx,sy,m,tx,ty_a,typ_a,typp_a,tint_a)
!
   call lagrange_initialize(lp,n)
   call newton_initialize(np,n)
   call spline_initialize(sp_l,n,1)
   call spline_initialize(sp_q,n,2)
   call spline_initialize(sp_c,n,3)
!
   call lagrange_set_x(lp,sx)
   call newton_set_x(np,sx)
   call spline_set_x(sp_l,sx)
   call spline_set_x(sp_q,sx)
   call spline_set_x(sp_c,sx)
   ! set y in source
   call lagrange_set(lp,sy)
   call newton_set(np,sy)
   call spline_set(sp_l,sy,yp1,ypp1,yppn)
   call spline_set(sp_q,sy,yp1,ypp1,yppn)
   call spline_set(sp_c,sy,yp1,ypp1,yppn)
!
   ! get y in target: 0,1,2
   ty = 0.0_r8
   ! lagragian
   call lagrange_get(lp,m,tx,ty(:,1),0)
   call lagrange_get(lp,m,tx,typ(:,1),1)
   call lagrange_get(lp,m,tx,typp(:,1),2)
   call lagrange_get(lp,m,tx,tint(:,1),2)
! spline: linear
   call spline_get(sp_l,m,tx,ty(:,2),0)
   call spline_get(sp_l,m,tx,typ(:,2),1)
   call spline_get(sp_l,m,tx,typp(:,2),2)
   call spline_integral(sp_l,m,tx,tint(:,2),-1)
   do i = 2,m
     tint(i,2) = tint(i-1,2)+tint(i,2)
   enddo
! spline: quadric
   call spline_get(sp_q,m,tx,ty(:,3),0)
   call spline_get(sp_q,m,tx,typ(:,3),1)
   call spline_get(sp_q,m,tx,typp(:,3),2)
   call spline_integral(sp_q,m,tx,tint(:,3),-1)
   do i = 2,m
     tint(i,3) = tint(i-1,3)+tint(i,3)
   enddo
! spline: cubic
   call spline_get(sp_c,m,tx,ty(:,4),0)
   call spline_get(sp_c,m,tx,typ(:,4),1)
   call spline_get(sp_c,m,tx,typp(:,4),2)
   call spline_integral(sp_c,m,tx,tint(:,4),-1)
   do i = 2,m
     tint(i,4) = tint(i-1,4)+tint(i,4)
   enddo
!
   call lagrange_finalize(lp)
   call newton_finalize(np)
   call spline_finalize(sp_l)
   call spline_finalize(sp_q)
   call spline_finalize(sp_c)
!
! write
#if 0
   call open_file('Result.dat',unitnum)
   call write_ascii(unitnum,testcase)
   call write_ascii(unitnum,nfun)
   call write_ascii(unitnum,'Lagrangian')
   call write_ascii(unitnum,'Linear Spline')
   call write_ascii(unitnum,'Quadratic Spline')
   call write_ascii(unitnum,'Cubic Spline')
   call write_ascii(unitnum,n)
   call write_ascii_lines(unitnum,n,sx)
   call write_ascii_lines(unitnum,n,sy)
   call write_ascii(unitnum,m)
   call write_ascii_lines(unitnum,m,tx)
   call write_ascii_lines(unitnum,m,ty_a)
   call write_ascii_lines(unitnum,m,typ_a)
   call write_ascii_lines(unitnum,m,typp_a)
   call write_ascii_lines(unitnum,m,tint_a)
   call write_ascii_lines(unitnum,m,ty(:,1))
   call write_ascii_lines(unitnum,m,typ(:,1))
   call write_ascii_lines(unitnum,m,typp(:,1))
   call write_ascii_lines(unitnum,m,tint(:,1))
   call write_ascii_lines(unitnum,m,ty(:,2))
   call write_ascii_lines(unitnum,m,typ(:,2))
   call write_ascii_lines(unitnum,m,typp(:,2))
   call write_ascii_lines(unitnum,m,tint(:,2))
   call write_ascii_lines(unitnum,m,ty(:,3))
   call write_ascii_lines(unitnum,m,typ(:,3))
   call write_ascii_lines(unitnum,m,typp(:,3))
   call write_ascii_lines(unitnum,m,tint(:,3))
   call write_ascii_lines(unitnum,m,ty(:,4))
   call write_ascii_lines(unitnum,m,typ(:,4))
   call write_ascii_lines(unitnum,m,typp(:,4))
   call write_ascii_lines(unitnum,m,tint(:,4))
   call close_file(unitnum)
#endif
!
   call create_ncfile('result/result.nc',testcase,isgll,nfun,norder,n,m,       &
                                            sx,sy,tx,ty_a,typ_a,typp_a,tint_a, &
                                            ty,typ,typp,tint)
!
   deallocate(sx)
   deallocate(sy)
   deallocate(tx)
   deallocate(ty_a)
   deallocate(typ_a)
   deallocate(typp_a)
   deallocate(tint_a)
   deallocate(ty)
   deallocate(typ)
   deallocate(typp)
   deallocate(tint)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_namelist(nlfilename)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*), optional, intent(in   ) :: nlfilename
   namelist/inputs/n, m, testcase, isgll, yp1, ypp1, yppn
!
   if (present(nlfilename)) then
     open(21,file = trim(nlfilename),status = 'old')
     read(21,nml = inputs)
     close(21)
   else
     read(*,nml = inputs)
   endif
!
   print*,'n        = ',n
   print*,'m        = ',m
   print*,'testcase = ',testcase
   print*,'isGLL    = ',isgll
   print*,'yp1      = ',yp1
   print*,'ypp1     = ',ypp1
   print*,'yppn     = ',yppn
!
   return
   end subroutine read_namelist
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_experiment(iexp, n_s, x_s, y_s, m_t, x_t, y_a, yp_a, ypp_a, tint_a)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: iexp, n_s, m_t
   real(r8), dimension(n_s), intent(inout) :: x_s, y_s
   real(r8), dimension(m_t), intent(inout) :: x_t, y_a, yp_a, ypp_a, tint_a
! local variables
   type(quadrature_t) :: qp
   integer(i4) :: i
   real(r8) :: xmin, xmax, dx_s, dx_t
   real(r8) :: a, mu
!
   if (isgll) then
     qp = gausslobatto(n_s)
     do i = 1,n_s
       x_s(i) = qp%points(i)
     enddo
     xmin = x_s(1)
     xmax = x_s(n_s)
   else
     xmin = 0.0_r8
     xmax = 4.0_r8*pi
     dx_s =(xmax-xmin)/dble(n_s-1)
     ! set x in source
     do i = 1,n_s
       x_s(i) = dx_s*dble(i-1)+xmin
     enddo
   endif
   mu =(xmax-xmin)*5.0_r8/(4.0_r8*pi)+xmin
!
   call get_func(iexp,n_s,x_s,y_s,xmin,xmax,mu,0)
!
   ! set x in target
   dx_t =(xmax-xmin)/dble(m_t-1)
   do i = 1,m_t
     x_t(i) = dx_t*dble(i-1)+xmin
   enddo
!
   call get_func(iexp,m_t,x_t,y_a,xmin,xmax,mu,0)
   call get_func(iexp,m_t,x_t,yp_a,xmin,xmax,mu,1)
   call get_func(iexp,m_t,x_t,ypp_a,xmin,xmax,mu,2)
   call get_func(iexp,m_t,x_t,tint_a,xmin,xmax,mu,-1)
!
   return
   end subroutine set_experiment
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program testmain
!-------------------------------------------------------------------------------
