#define USE_LAPACK
!#undef  USE_LAPACK
#undef PROFILE
!#define PROFILE
!-------------------------------------------------------------------------------
   module spline_remap
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    2017-02-09  junghan kim    initial
!
!  structure:
!
!-------------------------------------------------------------------------------
   implicit none
#ifdef PROFILE
   include 'mpif.h'
#endif
!
   private
!
   integer , parameter :: i4 = 4, l4 = 4, r4 = 4, r8 = 8
   real(r8), parameter :: zero=0.0d0, one=1.0d0, vtiny=1d-12, qmax=1d50
   logical(l4), parameter :: islimit = .true.
#ifdef PROFILE
   real(r8)    :: time1, time2
#endif
!
   type, public :: spline_t
     ! order (input)
     integer(i4) :: order
     ! monotonic (input)
     logical(l4) :: mono
     ! bndry condition (input)
     integer(i4) :: bndry
     ! size
     integer(i4) :: n, nn, ha ! ha(halo)
     ! array order
     integer(i4) :: array_order
     ! values, derivative
     real(r8), dimension(:),   allocatable :: x, y, cx, dx, rdx
     real(r8), dimension(:),   allocatable :: cy
     real(r8), dimension(:,:), allocatable :: a
     real(r8), dimension(:),   allocatable :: inty
   end type
!
   public :: spline_initialize, spline_finalize
   public :: spline_set_x, spline_set
   public :: spline_get_fun, spline_get_fun_integral
   public :: spline_get, spline_interface, spline_cyclic_interface
   public :: unittest_main
   public :: remap_psm, remap_psm_cyclic, test_remap_psm
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_initialize(sr, n, order, mono, bndry)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),        intent(inout) :: sr
   integer(i4),           intent(in   ) :: n
   integer(i4), optional, intent(in   ) :: order
   logical(i4), optional, intent(in   ) :: mono
   integer(i4), optional, intent(in   ) :: bndry
! local variables
   integer(i4) :: i, k, ha, nn
   logical :: mono_l
!
   sr%n = n
!
   if (present(mono)) then
     sr%mono = mono
   else
     sr%mono = .false.
   endif
!
   if (present(order)) then
     sr%order = order
   else
     sr%order = 2 ! parabolic
   endif
!
   if (present(bndry)) then
     sr%bndry = bndry
   else
     sr%bndry = -1
   endif
!
   if (sr%bndry.lt.-1.or.sr%bndry.gt.2) then
     print *,'error: not support bndry ',sr%bndry
     stop
   endif
!
   if (order.ne.2.and.order.ne.4) then
     print *,'error: not support order ',order
     stop
   endif
!
   if (order==4.and.sr%bndry.ge.0) then
     ha = 2
   else
     ha = 0
   endif
   sr%ha = ha
   nn    = sr%n+2*ha
   sr%nn = nn
!
   allocate(sr%x(1:nn))     ! x at grid point
   allocate(sr%y(1:nn))     ! y at grid point
   allocate(sr%cx(0:nn))    ! x at cell point
   allocate(sr%cy(0:nn))    ! y at cell point
   allocate(sr%dx(1:nn))    ! dx between cells
   allocate(sr%rdx(1:nn))   ! 1.0/dx between cells
   allocate(sr%inty(0:nn))  ! integral value at cell point:dx_i*y_i
   if (sr%order==2) then
     allocate(sr%a(1:nn,0:2)) ! y = a0+a1x1+a2x2
   elseif (sr%order==4) then
     allocate(sr%a(1:nn,0:4)) ! y = a0+a1x1+a2x2+a3x3+a4x4
   endif
!
   return
   end subroutine spline_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_finalize(sr)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t), intent(inout) :: sr
!
   deallocate(sr%x)
   deallocate(sr%y)
   deallocate(sr%cx)
   deallocate(sr%cy)
   deallocate(sr%dx)
   deallocate(sr%rdx)
   deallocate(sr%a)
   deallocate(sr%inty)
!
   return
   end subroutine spline_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_set_x(sr, cx_s)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),              intent(inout) :: sr
   real(r8), dimension(0:sr%n), intent(in   ) :: cx_s
! local variables
   integer(i4) :: i, n, nn, ha
!
! check array order
   sr%array_order = check_array_order(cx_s)
   if (sr%array_order==0) then
     print *,'error: check array order (source)'
     stop
   endif
!
   n  = sr%n
   nn = sr%nn
   ha = sr%ha
!
! set cell x
   sr%cx(0+ha:nn-ha) = cx_s(:)
   ! halo
   if (sr%bndry==-1) then ! mirror
   elseif (sr%bndry==0.or.sr%bndry==1) then ! mirror
     do i = 0,ha-1
       sr%cx(i)    = sr%cx(ha)-(sr%cx(2*ha-i)-sr%cx(ha))
       sr%cx(nn-i) = sr%cx(nn-ha)+(sr%cx(nn-ha)-sr%cx(nn-2*ha+i))
     enddo
   elseif (sr%bndry==2) then ! periodic
     do i = 0,ha-1
       sr%cx(i)    = sr%cx(ha)-(sr%cx(nn-ha)-sr%cx(nn-2*ha+i))
       sr%cx(nn-i) = sr%cx(nn-ha)+(sr%cx(2*ha-i)-sr%cx(ha))
     enddo
   else
     if (sr%order==2) then
       print *,'error: not support bndry condition...'
       stop
     endif
   endif
!
! set grid x
   !do i = 1,nn
   !  sr%x(i) = 0.5_r8*(sr%cx(i-1)+sr%cx(i))
   !enddo
   sr%x(1:nn) = 0.5_r8*(sr%cx(0:nn-1)+sr%cx(1:nn))
!
! set dx
   !do i = 1,nn
   !  sr%dx(i) = sr%cx(i)-sr%cx(i-1)
   !enddo
   sr%dx(1:nn)  = sr%cx(1:nn)-sr%cx(0:nn-1)
   sr%rdx(1:nn) = 1.0_r8/sr%dx(1:nn)
!
   return
   end subroutine spline_set_x
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_set(sr, y_s)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),            intent(inout) :: sr
   real(r8), dimension(sr%n), intent(in   ) :: y_s
! local variables
   integer(i4) :: i, n, nn, ha, err
   real(r8)    :: tmp, ls, rs, ts, lls, rrs, tts, yl, yu
   real(r8), allocatable, dimension(:) :: a, b, c, d, e, h, m
   ! for qartic monotone
   real(r8)    :: xp, xm, w, sp, sm, s0, s1, mu, dy
   integer(i4) :: nu
   integer(i4), dimension(:), allocatable :: do_filter
!
#ifdef PROFILE
time1 = mpi_wtime()
#endif
   n  = sr%n
   nn = sr%nn
   ha = sr%ha
!
! set grid y and halo
   sr%y(1+ha:nn-ha) = y_s(:)
! why write? below this code two times
   sr%inty(0) = 0
   do i = 1,nn
     sr%inty(i) = sr%inty(i-1)+sr%y(i)*sr%dx(i)
   enddo
!
#ifdef PROFILE
time2 = mpi_wtime()
print '(a30,f11.8)','initial = ',time2-time1
#endif
!
   if (sr%order==0) then
     print *,'not yet developed...'
     stop
   elseif (sr%order==1) then
     print *,'not yet developed...'
     stop
   elseif (sr%order==2) then
     !
     ! boundary
#ifdef PROFILE
time1 = mpi_wtime()
#endif
     if     (sr%bndry==0) then
       sr%cy(0) = -1.0*(sr%y(2)-sr%y(1))/(sr%x(2)-sr%x(1))*(sr%cx(0)-sr%x(1))+sr%y(1)
       sr%cy(n) = -1.0*(sr%y(n)-sr%y(n-1))/(sr%x(n)-sr%x(n-1))*(sr%cx(n)-sr%x(n-1))+sr%y(n-1)
     elseif (sr%bndry==1) then
       sr%cy(0) = (sr%y(2)-sr%y(1))/(sr%x(2)-sr%x(1))*(sr%cx(0)-sr%x(1))+sr%y(1)
       sr%cy(n) = (sr%y(n)-sr%y(n-1))/(sr%x(n)-sr%x(n-1))*(sr%cx(n)-sr%x(n-1))+sr%y(n-1)
     elseif (sr%bndry==2) then
       tmp = sr%cx(n)+(sr%x(1)-sr%cx(0))
       sr%cy(n) = (sr%y(1)-sr%y(n))/(tmp-sr%x(n))*(tmp-sr%cx(n))+sr%y(n)
       sr%cy(0) = sr%cy(n)
     endif
#ifdef PROFILE
time2 = mpi_wtime()
print '(a30,f11.8)','gen boundary = ',time2-time1
#endif
     !
     allocate(a(1:n-1),b(1:n-1),c(1:n-1))
     b(1:n-2) = sr%rdx(2:n-1)
     a(1:n-1) = 2.0_r8*(sr%rdx(1:n-1)+sr%rdx(2:n))
     c(1:n-1) = sr%rdx(2:n)
     !psi (just save at cy for linear equation)
     sr%cy(1:n-1) = 3.0_r8*(sr%y(1:n-1)*sr%rdx(1:n-1)+sr%y(2:n)*sr%rdx(2:n))
     sr%cy(1)     = sr%cy(1)-b(1)*sr%cy(1-1)
     sr%cy(n-1)   = sr%cy(n-1)-c(n-1)*sr%cy(n)
#ifdef PROFILE
time1 = mpi_wtime()
print '(a30,f11.8)','gen matrix = ',time1-time2
#endif
     !call inverse_matrix_to_linear_equation(n-1,a,b,c,sr%cy(1:n-1))
     call linear_equation_solver(n-1,a,b,c,sr%cy(1:n-1))
     deallocate(a,b,c)
#ifdef PROFILE
time2 = mpi_wtime()
print '(a30,f11.8)','linear eq = ',time2-time1
#endif
     !
     ! monotonicity (grid scale)
     allocate(do_filter(n))
     if (sr%mono) then
       !call monotone_grid_scale(n,sr%y,sr%cy)
       call monotone_grid_scale(n,sr%y,sr%cy,do_filter)
     endif
#ifdef PROFILE
time1 = mpi_wtime()
print '(a30,f11.8)','monotone (grid) = ',time1-time2
#endif
     ! cal. arguments
     sr%a(1:n,0) =      sr%cy(0:n-1)
     sr%a(1:n,1) = -4.0*sr%cy(0:n-1)-2.0*sr%cy(1:n)+6.0*sr%y(1:n)
     sr%a(1:n,2) =  3.0*sr%cy(0:n-1)+3.0*sr%cy(1:n)-6.0*sr%y(1:n)
#ifdef PROFILE
time2 = mpi_wtime()
print '(a30,f11.8)','cal a0, a1, a2  = ',time2-time1
#endif
     ! monotonicity (subgrid scale)
     if (sr%mono) then
       !call monotone_subgrid_scale_psm(n,sr%y,sr%cy,sr%a)
       call monotone_subgrid_scale_psm(n,sr%y,sr%cy,sr%a,do_filter,sr%array_order)
     endif
     deallocate(do_filter)
#ifdef PROFILE
time1 = mpi_wtime()
print '(a30,f11.8)','monotone (subgrid) = ',time1-time2
#endif
   !
   elseif (sr%order==3) then
     print *,'not yet developed...'
     stop
   elseif (sr%order==4) then
     !
     ! boundary
     if (sr%bndry==0) then ! mirror
       do i = 1,ha
         sr%y(i)       = sr%y(2*ha+1-i)
         sr%y(nn-ha+i) = sr%y(nn-ha+1-i)
       enddo
     elseif (sr%bndry==1) then ! mirror+derivative
       do i = 1,ha
         sr%y(i)       = sr%y(ha+1)-(sr%y(2*ha+2-i)-sr%y(ha+1))
         sr%y(nn-ha+i) = sr%y(nn-ha)+(sr%y(nn-ha)-sr%y(nn-ha-i))
       enddo
     elseif (sr%bndry==2) then ! periodic
       do i = 1,ha
         sr%y(i)       = sr%y(nn-2*ha+i)
         sr%y(nn-ha+i) = sr%y(ha+i)
       enddo
     else
       print *,'error: not support bndry condition...'
       stop
     endif
     !
     allocate(a(1:nn-1),b(1:nn-1),c(1:nn-1),d(1:nn-1),e(1:nn-1))
     allocate(h(1:nn))
     allocate(m(0:nn))
     h(1:nn) = sr%dx(1:nn)
     m(0) = 0.0
     m(nn) = 0.0
     ! pennta matrix
     a(1) = 1.0/120.0*( 2.0*h(1)*(3.0*h(1)+2.0*h(2))                           &
                       +(h(2)+h(3))*(10.0*h(1)+7.0*h(2))                       &
                       +h(2)*h(3)*(2.0*h(2)+3.0*h(3))/(h(2)+h(3)) )
     c(1) = 1.0/120.0*( (3.0*h(2)+4.0*h(3))*h(3)                               &
                       +3.0*h(2)*(h(2)+h(3))                                   &
                       +h(1)*h(2)**2.0/(h(1)+h(2))  )
     e(1) = 1.0/120.0*h(3)**3.0/(h(2)+h(3))
     do i=2,nn-2
       ls     = h(i-1)+h(i)
       rs     = h(i)+h(i+1)
       rrs    = h(i+1)+h(i+2)
       ts     = h(i-1)+h(i)+h(i+1)
       tts    = h(i)+h(i+1)+h(i+2)
if (i.ne.2) d(i-2)   = h(i-1)**3.0/(ls*ts)
       b(i-1)   = (4.0*h(i-1)+3.0*h(i)+ls/h(i-1)*(h(i)*h(i+1)/rs+2.0*h(i)))*h(i-1)/ts &
               +h(i)**2.0*rrs/(rs*tts)
       a(i)   = (h(i-1)*h(i)/ls+2.0*h(i)+ls/h(i-1)*(3.0*h(i)+4.0*h(i+1)))*h(i-1)/ts &
               +(h(i+1)*h(i+2)/rrs+2.0*h(i+1)+rrs/h(i+2)*(4.0*h(i)+3.0*h(i+1)))*h(i+2)/tts
       c(i)   = ls*h(i+1)**2.0/(rs*ts)                                              &
               +(3.0*h(i+1)+4.0*h(i+2)+rrs/h(i+2)*(h(i)*h(i+1)/rs+2.0*h(i+1)))*h(i+2)/tts
       e(i)   = h(i+2)**3.0/(rrs*tts)
     enddo
     d(nn-3) = 1.0/120.0*h(nn-2)**3.0/(h(nn-2)+h(nn-1))
     b(nn-2) = 1.0/120.0*( h(nn-2)*(4.0*h(nn-2)+3.0*h(nn-1))                       &
                         +3.0*h(nn-1)*(h(nn-2)+h(nn-1))                           &
                         +h(nn-1)**2.0*h(nn)/(h(nn-1)+h(nn))  ) 
     a(nn-1) = 1.0/120.0*( (3.0*h(nn-2)+2.0*h(nn-1))*h(nn-2)*h(nn-1)/(h(nn-2)+h(nn-1))&
                         +2.0*h(nn)*(2.0*h(nn-1)+3.0*h(nn))                       &
                         +(h(nn-2)+h(nn-1))*(7.0*h(nn-1)+10.0*h(nn)) )
     !psi (just save at m for linear equation)
     m(1)   = (sr%y(3)-sr%y(2))/(h(2)+h(3))-(sr%y(2)-sr%y(1))/(h(1)+h(2))
     do i=2,nn-2
       m(i) = 120.0*(h(i-1)+h(i)+h(i+1)+h(i+2))*divided_difference(nn+1,sr%cx(0:nn),sr%y(1:nn),i-2,i+2)
     enddo
     m(nn-1) = (sr%y(nn-1)-sr%y(nn-2))/(h(nn-2)+h(nn-1))-(sr%y(nn)-sr%y(nn-1))/(h(nn-1)+h(nn))
     !call inverse_matrix_to_linear_equation(nn-1,a,b,c,m(1:nn-1),d,e)
     call linear_equation_solver(nn-1,a,b,c,m(1:nn-1),d,e)
     deallocate(a,b,c,d,e)
     ! cal y_i+(1/2)
     m(0) = 0.0
     m(nn) = 0.0
     sr%cy(0)   = (2.0-h(2)/(h(1)+h(2)))*sr%y(1)                               &
                 -h(1)/(h(1)+h(2))*sr%y(2)                                     &
                 +h(1)/120.0*(3.0*h(1)**2.0+2.0*h(2)*(3.0*h(1)+2*h(2)))*m(1)   &
                 +h(1)*h(2)**3.0/(120.0*(h(1)+h(2)))*m(2)
     sr%cy(1)   = h(2)/(h(1)+h(2))*sr%y(1)                                     &
                 +h(1)/(h(1)+h(2))*sr%y(2)                                     &
                 -h(1)*h(2)*(3.0*h(1)+2.0*h(2))/60.0*m(1)                      &
                 -h(1)*h(2)**3.0/(120.0*(h(1)+h(2)))*m(2)
     do i=2,nn-2
       ls  = h(i-1)+h(i)
       rs  = h(i)+h(i+1)
       rrs = h(i+1)+h(i+2)
       ts  = h(i-1)+h(i)+h(i+1)
       tts = h(i)+h(i+1)+h(i+2)
       tmp = ts/(ls*h(i+1))+tts/(h(i)*rrs)
       sr%cy(i) =-h(i)/ls**2.0*sr%y(i-1)                                       &
                 +(3.0-h(i-1)*(2.0*h(i-1)+3.0*h(i))/ls**2.0)/h(i)*sr%y(i)      &
                 +(3.0-h(i+2)*(3.0*h(i+1)+2.0*h(i+2))/rrs**2.0)/h(i+1)*sr%y(i+1)&
                 -h(i+1)/rrs**2.0*sr%y(i+2)                                    &
                 -1.0/120.0*h(i-1)**3.0*h(i)/ls**2.0*m(i-2)                    &
                 -h(i)/120.0*(h(i-1)*(4.0*h(i-1)+3.0*h(i))/ls+2.0*h(i))*m(i-1) &
                 +1.0/120.0*(h(i+1)**2.0*(3.0+h(i+2)*(2.0*h(i+1)+3.0*h(i+2))/rs**2.0)  &
                            -h(i)**2.0*(3.0+h(i-1)*(3.0*h(i-1)+2.0*h(i))/ls**2.0) )*m(i) &
                 +h(i+1)/120.0*(2.0*h(i+1)+h(i+2)*(3.0*h(i+1)+4.0*h(i+2))/rs)*m(i+1) &
                 +1.0/120.0*h(i+1)*h(i+2)**3.0/rrs**2.0*m(i+2)
       sr%cy(i) =sr%cy(i)/tmp
     enddo
     sr%cy(nn-1) = h(nn-1)/(h(nn-1)+h(nn))*sr%y(nn)                                 &
                 +h(nn)/(h(nn-1)+h(nn))*sr%y(nn-1)                                 &
                 +h(nn-1)**3.0*h(nn)/(120.0*(h(nn-1)+h(nn)))*m(nn-2)                &
                 +h(nn-1)*h(nn)*(2.0*h(nn-1)+3.0*h(nn))/60.0*m(nn-1)
     sr%cy(nn)   = (2.0-h(nn-1)/(h(nn-1)+h(nn)))*sr%y(nn)                           &
                 -h(nn)/(h(nn-1)+h(nn))*sr%y(nn-1)                                 &
                 -h(nn-1)**3.0*h(nn)/(120.0*(h(nn-1)+h(nn)))*m(nn-2)                &
                 -h(nn)/120.0*(2.0*h(nn-1)*(2.0*h(nn-1)+3.0*h(nn))+3.0*h(nn)**2.0)*m(nn-1)
     ! monotonicity (grid scale)
     if (sr%mono) then
       call monotone_grid_scale(nn,sr%y,sr%cy)
     endif
     ! cal. argumennts
     do i=1,nn
       sr%a(i,0) =      sr%cy(i-1)
       sr%a(i,1) = -4.0*sr%cy(i-1)-2.0*sr%cy(i)+h(i)**3*(m(i-1)/20.0+m(i)/30.0)+6.0*sr%y(i)
       sr%a(i,2) =  3.0*sr%cy(i-1)+3.0*sr%cy(i)-h(i)**3*(m(i-1)*7.0/40.0+m(i)*3.0/40.0)-6.0*sr%y(i)
       sr%a(i,3) =  h(i)**3/6.0*m(i-1)
       sr%a(i,4) =  h(i)**3/24.0*(m(i)-m(i-1))
     enddo
!
     ! monotonicity (subgrid scale)
     if (sr%mono) then
       do i=1,nn
         w  = sr%a(i,3)**2.0-8.0*sr%a(i,2)*sr%a(i,4)/3.0
         xm = (-sr%a(i,3)-sqrt(w))/(4.0*sr%a(i,4))
         xp = (-sr%a(i,3)+sqrt(w))/(4.0*sr%a(i,4))
         sm = sr%a(i,1)+2.0*sr%a(i,2)*xm+3.0*sr%a(i,3)*xm**2.0+4.0*sr%a(i,4)*xm**3.0
         sp = sr%a(i,1)+2.0*sr%a(i,2)*xp+3.0*sr%a(i,3)*xp**2.0+4.0*sr%a(i,4)*xp**3.0
         s0 = sr%a(i,1)
         s1 = sr%a(i,1)+2.0*sr%a(i,2)+3.0*sr%a(i,3)+4.0*sr%a(i,4)
         mu = s0*s1
         if     (mu<0.0) then
           if (w>=0.0.and.(xm*(1.0-xm)>0).and.(xp*(1.0-xp)>0).and.(s0*sm<0.0)) then
             nu = 3
           else
             nu = sign(1.0d0,s0)
           endif
         elseif (mu>0.0) then
           if (w>=0.0 .and. (((xm*(1.0-xm)>0).and.(s0*sm<0.0)) .or. &
                             ((xp*(1.0-xp)>0).and.(s0*sp<0.0)) .or. &
                             ((xm*(1.0-xm)>0).and.(xp*(1.0-xp)>0).and.(sm*sp<0.0)) &
                                                                                )) then
             nu = 2
           else
             nu = 0
           endif
         elseif (mu==0.0) then
           if (w>=0.0.and.(xm*(1.0-xm)>0).and.(xp*(1.0-xp)>0).and.(s0*sm+s1*sp<0.0)) then
             nu = 2
           elseif (w>=0.0.and.( (xm*(1.0-xm)>0.and.(s0+s1)*sm<0.0).or.(xp*(1.0-xp)>0.and.(s0+s1)*sp<0.0) )) then
             nu = sign(1.0d0,s0-s1)
           elseif (w>=0.0.and.(xm*(1.0-xm)>0).and.(xp*(1.0-xp)>0).and.(s0+s1==0.0)) then
             nu = sign(1.0d0,sm)
           else
             nu = 0
           endif
         endif
         if (nu.ne.0) then
           dy = sr%cy(i)-sr%cy(i-1)
           if     ( (sr%y(i)-sr%cy(i-1))*(sr%y(i)-(sr%cy(i-1)+dy/5.0))<=0.0 ) then
             sr%a(i,0) =  sr%cy(i-1)
             sr%a(i,1) =  0.0
             sr%a(i,2) =  0.0
             sr%a(i,3) =  0.0
             sr%a(i,4) = -5.0*sr%cy(i-1)+5.0*sr%y(i)
           elseif ( (sr%y(i)-(sr%cy(i-1)+dy*1.0/5.0))*(sr%y(i)-(sr%cy(i-1)+dy*2.0/5.0))<=0.0 ) then
             sr%a(i,0) =  sr%cy(i-1)
             sr%a(i,1) =  0.0
             sr%a(i,2) =  0.0
             sr%a(i,3) = -16.0*sr%cy(i-1)-4.0*sr%cy(i)+20.0*sr%y(i)
             sr%a(i,4) =  15.0*sr%cy(i-1)+5.0*sr%cy(i)-20.0*sr%y(i)
           elseif ( (sr%y(i)-(sr%cy(i-1)+dy*2.0/5.0))*(sr%y(i)-(sr%cy(i-1)+dy*3.0/5.0))<=0.0 ) then
             sr%a(i,0) =  sr%cy(i-1)
             sr%a(i,1) =  0.0
             sr%a(i,2) = -18.0*sr%cy(i-1)-12.0*sr%cy(i)+30.0*sr%y(i)
             sr%a(i,3) =  32.0*sr%cy(i-1)+28.0*sr%cy(i)-60.0*sr%y(i)
             sr%a(i,4) = -15.0*sr%cy(i-1)-15.0*sr%cy(i)+30.0*sr%y(i)
           elseif ( (sr%y(i)-(sr%cy(i-1)+dy*3.0/5.0))*(sr%y(i)-(sr%cy(i-1)+dy*4.0/5.0))<=0.0 ) then
             sr%a(i,0) =  sr%cy(i-1)
             sr%a(i,1) =  -8.0*sr%cy(i-1)-12.0*sr%cy(i)+20.0*sr%y(i)
             sr%a(i,2) =  18.0*sr%cy(i-1)+42.0*sr%cy(i)-60.0*sr%y(i)
             sr%a(i,3) = -16.0*sr%cy(i-1)-44.0*sr%cy(i)+60.0*sr%y(i)
             sr%a(i,4) =   5.0*sr%cy(i-1)+15.0*sr%cy(i)-20.0*sr%y(i)
           elseif ( (sr%y(i)-(sr%cy(i-1)+dy*4.0/5.0))*(sr%y(i)-sr%cy(i))<=0 ) then
             sr%a(i,0) =  -4.0*sr%cy(i)+ 5.0*sr%y(i)
             sr%a(i,1) =  20.0*sr%cy(i)-20.0*sr%y(i)
             sr%a(i,2) = -30.0*sr%cy(i)+30.0*sr%y(i)
             sr%a(i,3) =  20.0*sr%cy(i)-20.0*sr%y(i)
             sr%a(i,4) =  -5.0*sr%cy(i)+ 5.0*sr%y(i)
           else
             sr%a(i,0) =  sr%y(i)
             sr%a(i,1) =  0.0
             sr%a(i,2) =  0.0
             sr%a(i,3) =  0.0
             sr%a(i,4) =  0.0
           endif
         endif ! if (nu!=0)
       enddo ! do i = 1,n
     endif ! if (sr%mono)
!
! integral
   sr%inty(0) = 0
   do i = 1,nn
     sr%inty(i) = sr%inty(i-1)+sr%y(i)*sr%dx(i)
   enddo
!
     deallocate(h)
     deallocate(m)
   endif ! sr%order==4
!
   return
   end subroutine spline_set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine monotone_grid_scale(n, y, cy, do_filter)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),                 intent(in   ) :: n
   real(r8)   , dimension(n)  , intent(in   ) :: y
   real(r8)   , dimension(0:n), intent(inout) :: cy
   integer(i4), dimension(n)  , intent(inout), optional :: do_filter
! local variables
   integer(i4) :: i
   ! for cam
   integer(i4) :: im2, im1, ii, ip1, ip2, iym, iyp, t1, t2, t3
   real(r8), dimension(:), allocatable :: dy
!
   if (.not.present(do_filter)) then ! original
     !
     do i = 1,n-1
       if ((cy(i)-y(i))*(y(i+1)-cy(i)).lt.0.0) then
         if     (abs(cy(i)-y(i)).le.abs(cy(i)-y(i+1))) then
           cy(i) = y(i)
         elseif (abs(cy(i)-y(i)).gt.abs(cy(i)-y(i+1))) then
           cy(i) = y(i+1)
         endif
       endif
     enddo
     !
   else ! cam psm monotonicity
     ! 
     allocate(dy(n))
     do_filter = 0
     dy(1:n-1) = y(2:n)-y(1:n-1)
     dy(n)     = dy(n-1)
     dy = merge(zero,dy,abs(dy)<vtiny)
     ! 
     do i=0,n-1
       iym = max(1,i  )
       iyp = min(n,i+1)

       im2 = max(1,i-2)
       im1 = max(1,i-1)
       ii  = max(1,i  )
       ip1 = min(n,i+1)
       ip2 = min(n,i+2)
       ! t1: cy need monotone(0), t2: cy is local overshot(0), t3: left(0) or right(1)
       t1  = merge(1,0,(y(iyp)-cy(i))*(cy(i)-y(iym))>=0)
       t2  = merge(1,0,dy(im1)*(cy(i)-y(iym))>0.and.dy(im1)*dy(im2)>0.and.     &
                              dy(ip1)*dy(ip2)>0.and.dy(im1)*dy(ip1)<0)
       t3  = merge(1,0,abs(cy(i)-y(iym))>abs(cy(i)-y(iyp)))
       ! if do_filter==1: need monotone
       do_filter(i+1) = merge(0,1,t1+t2>0)
       cy(i) = (1-do_filter(i+1))*cy(i)+do_filter(i+1)*(t3*y(iyp)+(1-t3)*y(iym))
       do_filter(ii) = max(do_filter(ii),do_filter(i+1))
     enddo
     !
     if (islimit) then
       cy = merge(qmax,cy,cy>qmax)
       cy = merge(zero,cy,cy<zero)
     endif
     deallocate(dy)
     ! 
   endif
!
   end subroutine monotone_grid_scale
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine monotone_subgrid_scale_psm(n, y, cy, arg, do_filter, array_order)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),                   intent(in   ) :: n
   real(r8)   , dimension(n)    , intent(in   ) :: y
   real(r8)   , dimension(0:n)  , intent(in   ) :: cy
   real(r8)   , dimension(n,0:2), intent(inout) :: arg
   integer(i4), dimension(n)    , intent(inout), optional :: do_filter
   integer(i4),                   intent(in   ), optional :: array_order
! local variables
   integer(i4) :: i, k
   real(r8)    :: yl, yu
   ! for cam
   real(r8), dimension(n) :: dy,cyr
   real(r8)    :: f_xm,lev1,lev2,lev3,lev4,lev5,peaks_min,peaks_max,xm,xm_d
   integer(i4) :: peaks,im2,im1,ip1,ip2,lt1,lt2,lt3,t0,t1,t2,t3,t4,tm,tp
!
   if (.not.present(do_filter)) then ! original
     !
     do i=1,n
       !
       if (cy(i-1)<cy(i)) then
         yl = 2.0/3.0*cy(i-1)+1.0/3.0*cy(i)
         yu = 1.0/3.0*cy(i-1)+2.0/3.0*cy(i)
         if     (y(i).lt.yl) then
           arg(i,0) =  cy(i-1)
           arg(i,1) =  0.0
           arg(i,2) = -3.0*cy(i-1)+3.0*y(i)
         elseif (y(i).gt.yu) then
           arg(i,0) = -2.0*cy(i)+3.0*y(i)
           arg(i,1) =  6.0*cy(i)-6.0*y(i)
           arg(i,2) = -3.0*cy(i)+3.0*y(i)
         endif
       endif
       !
       if (cy(i-1)>cy(i)) then
         yu = 2.0/3.0*cy(i-1)+1.0/3.0*cy(i)
         yl = 1.0/3.0*cy(i-1)+2.0/3.0*cy(i)
         if     (y(i).lt.yl) then
           arg(i,0) = -2.0*cy(i)+3.0*y(i)
           arg(i,1) =  6.0*cy(i)-6.0*y(i)
           arg(i,2) = -3.0*cy(i)+3.0*y(i)
         elseif (y(i).gt.yu) then
           arg(i,0) =  cy(i-1)
           arg(i,1) =  0.0
           arg(i,2) = -3.0*cy(i-1)+3.0*y(i)
         endif
       endif
       !
       !if ((y(i)-cy(i-1))*(y(i)-cy(i)).gt.0.0) then
       if ((y(i)-cy(i-1))*(cy(i)-y(i)).le.0.0) then
         arg(i,0) = y(i)
         arg(i,1) = 0.0
         arg(i,2) = 0.0
       endif
       !
     enddo
     !
   else ! cam psm monotonicity
     !
     dy(1:n) = cy(1:n)-cy(0:n-1)
     dy = merge(zero,dy,abs(dy)<vtiny)
   
     cyr = cy(1:n)
   
     do k=1,n
       xm_d = merge(one,2*arg(k,2),abs(arg(k,2))<vtiny)
       if (array_order==1) then
         xm = merge(zero,-arg(k,1)/xm_d,abs(arg(k,2))<vtiny)
       else
         xm = merge(1.d0,-arg(k,1)/xm_d,abs(arg(k,2))<vtiny)
       endif
       f_xm = arg(k,0)+arg(k,1)*xm+arg(k,2)*xm**2
   
       t1 = merge(1,0,abs(arg(k,2))>vtiny)
       t2 = merge(1,0,xm<=zero.or.xm >=1)
       t3 = merge(1,0,arg(k,2)>zero)
       t4 = merge(1,0,arg(k,2)<zero)
       tm = merge(1,0,t1*((1-t2)+t3)==2)          ! local min
       tp = merge(1,0,t1*((1-t2)+(1-t3)+t4)==3)   ! local max
   
       ! peaks = -1(local min), 1(local max), 0(other)
       peaks = 0
       peaks = merge(-1,peaks,tm==1)
       peaks = merge(+1,peaks,tp==1)
       peaks_min = merge(f_xm,min(arg(k,0),arg(k,0)+arg(k,1)+arg(k,2)),tm==1)
       peaks_max = merge(f_xm,max(arg(k,0),arg(k,0)+arg(k,1)+arg(k,2)),tp==1)
   
       im2=max(1,k-2)
       im1=max(1,k-1)
       ip1=min(n,k+1)
       ip2=min(n,k+2)
   
       t1 = merge(abs(peaks),0,(dy(im2)*dy(im1)<=vtiny).or.                          &
                           (dy(ip1)*dy(ip2)<=vtiny).or.(dy(im1)*dy(ip1)>=vtiny).or.  &
                           (dy(im1)*real(peaks,r8)<=vtiny))
   
       if (islimit) then
         do_filter(k) = merge(1,t1+(1-t1)*do_filter(k),(cy(k-1)>=qmax).or.     &
                      (cy(k-1)<=zero).or.(peaks_max>qmax).or.(peaks_min<vtiny))
       else
         do_filter(k) = t1+(1-t1)*do_filter(k)
       endif
   
       if (do_filter(k)>0) then
         lev1 = cy(k-1)
         lev2 = (2*cy(k-1)+cyr(k))/3_r8
         lev3 = 0.5*(cy(k-1)+cyr(k))
         lev4 = (1/3_r8)*cy(k-1)+2*(1/3_r8)*cyr(k)
         lev5 = cyr(k)
   
         t1 = merge(1,0,cyr(k)>=cy(k-1))
         t2 = merge(1,0,y(k)<=lev1.or. y(k)>=lev5)
         t3 = merge(1,0,y(k)> lev1.and.y(k)< lev2)
         t4 = merge(1,0,y(k)> lev4.and.y(k)< lev5)
   
         lt1 = t1*t2
         lt2 = t1*(1-t2+t3)
         lt3 = t1*(1-t2+1-t3+t4)
   
         arg(k,0) = merge(y(k),arg(k,0),lt1==1)
         arg(k,1) = merge(zero,arg(k,1),lt1==1)
         arg(k,2) = merge(zero,arg(k,2),lt1==1)
   
         arg(k,0) = merge(cy(k-1),arg(k,0),lt2==2)
         arg(k,1) = merge(zero,arg(k,1),lt2==2)
         arg(k,2) = merge(3*(y(k)-cy(k-1)),arg(k,2),lt2==2)
   
         arg(k,0) = merge(-2*cyr(k)+3*y(k),arg(k,0),lt3==3)
         arg(k,1) = merge(+6*cyr(k)-6*y(k),arg(k,1),lt3==3)
         arg(k,2) = merge(-3*cyr(k)+3*y(k),arg(k,2),lt3==3)
   
         t2 = merge(1,0,y(k)>=lev1.or. y(k)<=lev5)
         t3 = merge(1,0,y(k)< lev1.and.y(k)> lev2)
         t4 = merge(1,0,y(k)< lev4.and.y(k)> lev5)
   
         lt1 = (1-t1)*t2
         lt2 = (1-t1)*(1-t2+t3)
         lt3 = (1-t1)*(1-t2+1-t3+t4)
   
         arg(k,0) = merge(y(k),arg(k,0),lt1==1)
         arg(k,1) = merge(zero,arg(k,1),lt1==1)
         arg(k,2) = merge(zero,arg(k,2),lt1==1)
   
         arg(k,0) = merge(cy(k-1),arg(k,0),lt2==2)
         arg(k,1) = merge(zero,arg(k,1),lt2==2)
         arg(k,2) = merge(3*(y(k)-cy(k-1)),arg(k,2),lt2==2)
   
         arg(k,0) = merge(-2*cyr(k)+3*y(k),arg(k,0),lt3==3)
         arg(k,1) = merge(+6*cyr(k)-6*y(k),arg(k,1),lt3==3)
         arg(k,2) = merge(-3*cyr(k)+3*y(k),arg(k,2),lt3==3)
       endif
     enddo
     !
   endif
!
   end subroutine monotone_subgrid_scale_psm
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   recursive function divided_difference(n, x, y, is, ie) result(res)
!  recursive algorithm of divided difference (used newton polynomial)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(0:n), intent(in   ) :: x
   real(r8), dimension(n),   intent(in   ) :: y
   integer(i4),              intent(in   ) :: is, ie
   real(r8)                                :: res
! local variables
   integer(i4) :: i, j
!
   if ((is.lt.0.or.is.gt.n).or.(ie.lt.0.or.ie.gt.n)) then
     print *,'error: check is or ie in divided_difference...'
     stop
   endif
!
   if (is+1==ie) then
     res = y(ie)
   else
     res = (divided_difference(n,x,y,is+1,ie)-divided_difference(n,x,y,is,ie-1))/(x(ie)-x(is))
   endif
!
   end function divided_difference
#if 0
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_print_status(sr)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t), intent(inout) :: sr
! local variables
   integer :: k
   character(len=20) :: form_idx, form_n, form_n_p1, form_n_halo, form_n_halo_p1, form_n_halo_m1
   real(r8) :: mass
!
   write(form_idx,'(a,i0.2,a)') '(a,x,',sr%n+2*sr%ha,'(i7.3)) '
   write(form_n,'(a,i0.2,a)') '(a,x,',sr%n,'(f7.3)) '
   write(form_n_p1,'(a,i0.2,a)') '(a,x,',sr%n+1,'(f7.3)) '
   write(form_n_halo_p1,'(a,i0.2,a)') '(a,x,',sr%n+2*sr%ha+1,'(f7.3)) '
   write(form_n_halo_m1,'(a,i0.2,a)') '(a,x,',sr%n+2*sr%ha-1,'(f7.3)) '
   write(form_n_halo,'(a,i0.2,a)') '(a,x,',sr%n+2*sr%ha,'(f7.3)) '
!
   print form_idx,'idx = ',(/-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21/)
   print form_n_halo,' x = ',sr%x(:)
   print form_n_halo_p1,'cx = ',sr%cx(:)
   print form_n_halo,'dx = ',sr%dx(:)
   print form_n_halo,' y = ',sr%y(:)
   print form_n_p1,'cy = ',sr%cy(:)
   print form_n_halo,'dy = ',sr%dy(:)
   print form_n_halo,'dely = ',sr%dely(:)
!
   return
   end subroutine spline_print_status
#endif
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_spline_value_in_piece_ptr(sr, kp, x, vorder) result(value)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t), intent(inout) :: sr
   integer(i4),    intent(in   ) :: kp
   real(r8),       intent(in   ) :: x
   integer(i4),    intent(in   ) :: vorder ! 0:function, 1:first derivative, 2:second derivative
   real(r8)                      :: value
! local variables
   integer(i4) :: n
   real(r8)    :: xi
!
   xi = (x-sr%cx(kp-1))*sr%rdx(kp)
!
   ! integration
   if (vorder==-1) then
   !
     if     (sr%order==0) then

     elseif (sr%order==1) then

     elseif (sr%order==2) then
       value = sr%dx(kp)*(sr%a(kp,0)*xi+sr%a(kp,1)*xi**2.0/2.0+sr%a(kp,2)*xi**3.0/3.0)
     elseif (sr%order==3) then

     elseif (sr%order==4) then
       value = sr%dx(kp)*(sr%a(kp,0)*xi+sr%a(kp,1)*xi**2.0/2.0+sr%a(kp,2)*xi**3.0/3.0 &
                         +sr%a(kp,3)*xi**4.0/4.0+sr%a(kp,4)*xi**5.0/5.0)
     endif
   !
   ! 0th value: interpolation value
   elseif (vorder==0) then
   !
     if     (sr%order==0) then
     elseif (sr%order==1) then
     elseif (sr%order==2) then
       value = sr%a(kp,0)+sr%a(kp,1)*xi+sr%a(kp,2)*xi**2.0
     elseif (sr%order==3) then
     elseif (sr%order==4) then
       value = sr%a(kp,0)+sr%a(kp,1)*xi+sr%a(kp,2)*xi**2.0+sr%a(kp,3)*xi**3.0+sr%a(kp,4)*xi**4.0
     endif
   !
   ! 1th value: y'
   elseif (vorder==1) then
   !
     if     (sr%order==0) then
     elseif (sr%order==1) then
     elseif (sr%order==2) then
       value = sr%a(kp,1)+2.0*sr%a(kp,2)*xi
     elseif (sr%order==3) then
     elseif (sr%order==4) then
       value = sr%a(kp,1)+2.0*sr%a(kp,2)*xi+3.0*sr%a(kp,3)*xi**2.0+4.0*sr%a(kp,4)*xi**3.0
     endif
   !
   ! 2th value: y''
   elseif (vorder==2) then
   !
     if     (sr%order==0) then
     elseif (sr%order==1) then
     elseif (sr%order==2) then
       value = 2.0*sr%a(kp,2)
     elseif (sr%order==3) then
     elseif (sr%order==4) then
       value = 2.0*sr%a(kp,2)+6.0*sr%a(kp,3)*xi+12.0*sr%a(kp,4)*xi**2.0
     endif
   !
   ! 3th value: y'''
   elseif (vorder==3) then
   !
     if     (sr%order==0) then
     elseif (sr%order==1) then
     elseif (sr%order==2) then
       value = 0.0
     elseif (sr%order==3) then
     elseif (sr%order==4) then
       value = 6.0*sr%a(kp,3)+24.0*sr%a(kp,4)*xi
     endif
!
   endif ! vorder
!
   end function get_spline_value_in_piece_ptr
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_spline_value_in_piece(sr, vorder, m, kp, x, value)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),            intent(inout) :: sr
   integer(i4),               intent(in   ) :: vorder ! 0:function, 1:first derivative, 2:second derivative
   integer(i4),               intent(in   ) :: m
   integer(i4), dimension(m), intent(in   ) :: kp
   real(r8)   , dimension(m), intent(in   ) :: x
   real(r8)   , dimension(m), intent(inout) :: value
! local variables
   integer(i4) :: nn
   real(r8)    :: xi(m)
!   real(r8), dimension(:), allocatable :: a, b, c
!
   nn = sr%nn
   !xi = (x-sr%cx(kp-1))*sr%rdx(kp)
   xi = (x-sr%cx(kp-1))*sr%rdx(kp)
!
#if 0
   allocate(a(nn),b(nn),c(nn))
   a(:) = sr%a(:,0)
   b(:) = sr%a(:,1)
   c(:) = sr%a(:,2)
time1 = mpi_wtime()
   value = sr%dx(kp)*(a(kp)*xi+b(kp)*xi**2.0/2.0+c(kp)*xi**3.0/3.0)
   deallocate(a,b,c)
time2 = mpi_wtime()
print '(a30,f11.8)','my_psm_get_idx(func) = ',time2-time1
#endif
!
   ! integration
   if (vorder==-1) then
   !
     if     (sr%order==0) then

     elseif (sr%order==1) then

     elseif (sr%order==2) then
       value = sr%dx(kp)*(sr%a(kp,0)*xi+sr%a(kp,1)*xi**2.0/2.0+sr%a(kp,2)*xi**3.0/3.0)
     elseif (sr%order==3) then

     elseif (sr%order==4) then
       value = sr%dx(kp)*(sr%a(kp,0)*xi+sr%a(kp,1)*xi**2.0/2.0+sr%a(kp,2)*xi**3.0/3.0 &
                         +sr%a(kp,3)*xi**4.0/4.0+sr%a(kp,4)*xi**5.0/5.0)
     endif
   !
   ! 0th value: interpolation value
   elseif (vorder==0) then
   !
     if     (sr%order==0) then
     elseif (sr%order==1) then
     elseif (sr%order==2) then
       value = sr%a(kp,0)+sr%a(kp,1)*xi+sr%a(kp,2)*xi**2.0
     elseif (sr%order==3) then
     elseif (sr%order==4) then
       value = sr%a(kp,0)+sr%a(kp,1)*xi+sr%a(kp,2)*xi**2.0+sr%a(kp,3)*xi**3.0+sr%a(kp,4)*xi**4.0
     endif
   !
   ! 1th value: y'
   elseif (vorder==1) then
   !
     if     (sr%order==0) then
     elseif (sr%order==1) then
     elseif (sr%order==2) then
       value = sr%a(kp,1)+2.0*sr%a(kp,2)*xi
     elseif (sr%order==3) then
     elseif (sr%order==4) then
       value = sr%a(kp,1)+2.0*sr%a(kp,2)*xi+3.0*sr%a(kp,3)*xi**2.0+4.0*sr%a(kp,4)*xi**3.0
     endif
   !
   ! 2th value: y''
   elseif (vorder==2) then
   !
     if     (sr%order==0) then
     elseif (sr%order==1) then
     elseif (sr%order==2) then
       value = 2.0*sr%a(kp,2)
     elseif (sr%order==3) then
     elseif (sr%order==4) then
       value = 2.0*sr%a(kp,2)+6.0*sr%a(kp,3)*xi+12.0*sr%a(kp,4)*xi**2.0
     endif
   !
   ! 3th value: y'''
   elseif (vorder==3) then
   !
     if     (sr%order==0) then
     elseif (sr%order==1) then
     elseif (sr%order==2) then
       value = 0.0
     elseif (sr%order==3) then
     elseif (sr%order==4) then
       value = 6.0*sr%a(kp,3)+24.0*sr%a(kp,4)*xi
     endif
!
   endif ! vorder
!
   end subroutine get_spline_value_in_piece
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_get_fun(sr, m, x, y, vorder)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),         intent(inout) :: sr
   integer(i4),            intent(in   ) :: m
   real(r8), dimension(m), intent(in   ) :: x
   real(r8), dimension(m), intent(  out) :: y
   integer(i4), optional,  intent(in   ) :: vorder ! 0:function, 1:first derivative, 2:second derivative
! local variables
   integer(i4) :: kp(m), k, i, nn, l_vorder
   logical(l4) :: in_src
   real(r8), dimension(m) :: xx
!
   if (present(vorder)) then
     l_vorder = vorder
   else
     l_vorder = 0
   endif
   if (l_vorder.ne.-1.and.l_vorder.ne.0.and.l_vorder.ne.1.and.l_vorder.ne.2.and.l_vorder.ne.3) then
     print *,'vorder must be -1,0,1,2,3 in spline_get_fun...'
     stop
   endif
!
   nn = sr%nn
!
! check boundary
   if (sr%bndry/=2) then
     in_src = check_array_boundary(sr%cx(0:nn),x,sr%array_order)
     if (.not.in_src) then
       print*,'error: check boundary.... in spline_get_fun'
       stop
     endif
   endif
!
   xx = x
   call get_idx_in_array(nn+1,sr%cx(0:nn),m,xx,sr%array_order,sr%bndry,kp)
   call get_spline_value_in_piece(sr,l_vorder,m,kp,xx,y)
!
   return
   end subroutine spline_get_fun
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_get_fun_integral(sr, m, x, y)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),         intent(inout) :: sr
   integer(i4),            intent(in   ) :: m
   real(r8), dimension(m), intent(in   ) :: x
   real(r8), dimension(m), intent(  out) :: y
! local variables
   integer(i4) :: kp(m)
   integer(i4) :: nk, k, i, nn, ncycle
   real(r8) :: area, tarea
   real(r8), dimension(m) :: xx, areap(m)
!
   nn = sr%nn
   xx = x  ! for cyclic boundary condition
!
!#define PROFILE_GET
#ifdef PROFILE_GET
time1 = mpi_wtime()
#endif
   call get_idx_in_array(nn+1,sr%cx(0:nn),m,xx,sr%array_order,sr%bndry,kp)
#ifdef PROFILE_GET
time2 = mpi_wtime()
print '(a30,f11.8)','my_psm_get_idx = ',time2-time1
time1 = mpi_wtime()
#endif
!
   call get_spline_value_in_piece(sr,-1,m,kp,xx,areap)
!
   tarea  = sr%inty(kp(1)-1)+areap(1)
   y(1)   = 0.0_r8
   ncycle = 0
   do i  = 1,m-1
     area = sr%inty(kp(i+1)-1)+areap(i+1)
     if (kp(i+1).lt.kp(i).or.ncycle>0) then
       if (kp(i+1).lt.kp(i)) ncycle  = ncycle+1
       area    = real(ncycle,r8)*sr%inty(nn)+area
     endif
     !
     y(i+1) = area-tarea
     tarea  = area
   enddo ! do i=1,m-1
#ifdef PROFILE_GET
time2 = mpi_wtime()
print '(a30,f11.8)','my_psm_get_int = ',time2-time1
#endif
!
   return
   end subroutine spline_get_fun_integral
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_get(sr, m, cx, y)
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t),           intent(inout) :: sr
   integer(i4),              intent(in   ) :: m
   real(r8), dimension(0:m), intent(in   ) :: cx
   real(r8), dimension(m),   intent(  out) :: y
! local variables
   real(r8), dimension(0:m) :: icy
   integer(i4) :: k, nn, array_order
   logical(l4) :: in_src
!
   nn = sr%nn
!
   ! check array order
   array_order = check_array_order(cx)
   if (sr%array_order.ne.array_order) then
     print *,'error: check array order (target)',sr%array_order,array_order
     stop
   endif
   ! check boundary
   if (sr%bndry/=2) then
     in_src = check_array_boundary(sr%cx(0:nn),cx,array_order)
     if (.not.in_src) then
       print*,'error: check boundary.... in spline_get'
       stop
     endif
   endif
!
! integration
   call spline_get_fun_integral(sr,m+1,cx,icy)
   y(1:m) = icy(1:m)/(cx(1:m)-cx(0:m-1))
!
   return
   end subroutine spline_get
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_interface(n, scx, sy, m, tcx, ty, method, ismono, bndry, l, ftx, fty)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(0:n), intent(in   ) :: scx
   real(r8), dimension(n)  , intent(in   ) :: sy
   integer(i4),              intent(in   ) :: m
   real(r8), dimension(0:m), intent(in   ) :: tcx
   real(r8), dimension(m)  , intent(  out) :: ty
   integer(i4), optional   , intent(in   ) :: l
   real(r8)   , optional   , intent(in   ) :: ftx(:)
   real(r8)   , optional   , intent(  out) :: fty(:)
   ! method = 2:parabolic
   integer(i4),              intent(in   ) :: method
   logical(i4),              intent(in   ) :: ismono
   integer(i4),              intent(in   ) :: bndry
! local variables
   type(spline_t) :: sr
!
   call spline_initialize(sr,n,method,ismono,bndry)
   call spline_set_x(sr,scx)
   call spline_set(sr,sy)
   call spline_get(sr,m,tcx,ty)
   if (present(l).and.present(ftx).and.present(fty)) then
     call spline_get_fun(sr,l,ftx,fty,0)
   endif
   call spline_finalize(sr)
!
   return
   end subroutine spline_interface
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine spline_cyclic_interface(n, scx, sy, m, tcx, ty, method, ismono, l, ftx, fty)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(0:n), intent(in   ) :: scx
   real(r8), dimension(n)  , intent(in   ) :: sy
   integer(i4),              intent(in   ) :: m
   real(r8), dimension(0:m), intent(in   ) :: tcx
   real(r8), dimension(m)  , intent(  out) :: ty
   integer(i4), optional   , intent(in   ) :: l
   real(r8)   , optional   , intent(in   ) :: ftx(:)
   real(r8)   , optional   , intent(  out) :: fty(:)
   ! method = 0:constant, 1:linear, 2:parabolic
   integer(i4),              intent(in   ) :: method
   logical(i4),              intent(in   ) :: ismono
! local variables
   type(spline_t) :: sr
   real(r8)       :: lbase, rbase
   integer(i4)    :: i, nn, array_order
   real(r8), allocatable, dimension(:) :: escx, esy
!
   nn = 3*n
   allocate(escx(0:nn),esy(nn))
!
! expand x (source)
   array_order = check_array_order(scx)
   if (array_order==0) then
     print*,'check array order of source in spline_cyclic_interface.'
     stop
   endif
   rbase           = scx(n)-scx(0)
   lbase           = scx(0)-scx(n)
   escx(n:2*n)     = scx(:)
   escx(2*n+1:3*n) = scx(1:n)+rbase   ! right
   escx(0:n-1)     = scx(0:n-1)+lbase ! left
!
! expand y (source)
   esy(n+1:2*n)   = sy(:)
   esy(1:n)       = sy(1:n)
   esy(2*n+1:3*n) = sy(1:n)
!
   call spline_initialize(sr,nn,method,ismono,2)
   call spline_set_x(sr,escx)
   call spline_set(sr,esy)
   call spline_get(sr,m,tcx,ty)
   if (present(l).and.present(ftx).and.present(fty)) then
     call spline_get_fun(sr,l,ftx,fty,0)
   endif
   call spline_finalize(sr)
!
   deallocate(escx,esy)
!
   return
   end subroutine spline_cyclic_interface
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function check_array_order(array) result(order) ! order=0(mix), order=1(ascend), order=-1(descend)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:), intent(in   ) :: array
   integer(i4)                           :: order
! local variables
   integer(i4) :: n, i
!
   n = size(array)
   order = 0
!
   !do i = 2,n
   do i = 2,5
     if (i==2) then
       if (array(i).ge.array(i-1)) then
         order = 1
       else
         order = -1
       endif
     else
       if ((array(i).gt.array(i-1).and.order==-1).or.(array(i).lt.array(i-1).and.order==1)) then
         order = 0
         exit
       endif
     endif
   enddo
!
   end function check_array_order
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function check_array_boundary(src,tgt,order) result(in_src)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(:), intent(in   ) :: src, tgt
   integer(i4),            intent(in   ) :: order
   logical(l4)                           :: in_src
! local variables
   integer(i4) :: n, m
!
   n = size(src)
   m = size(tgt)
   in_src = .true.
   if (order==1) then
     if (tgt(1).lt.src(1).or.tgt(m).gt.src(n)) then
       in_src = .false.
     endif
   elseif (order==-1) then
     if (tgt(1).gt.src(1).or.tgt(m).lt.src(n)) then
       in_src = .false.
     endif
   else
     print *,'error: check_array_boundary'
     stop
   endif
!
   end function check_array_boundary
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_idx_in_array(n, array, m, value, array_order, bndry, idx)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),               intent(in   ) :: n
   real(r8)   , dimension(n), intent(in   ) :: array
   integer(i4),               intent(in   ) :: m
   real(r8)   , dimension(m), intent(inout) :: value
   integer(i4),               intent(in   ) :: array_order, bndry
   integer(i4), dimension(m), intent(inout) :: idx
! local variables
   integer(i4) :: i, j, ip
   real(r8)    :: xmin, xmax, tmp
!
   if (bndry==2) then
     if (array_order==1) then
       xmin = array(1)
       xmax = array(n)
     else
       xmin = array(n)
       xmax = array(1)
     endif
   endif
!
   idx = -1
   ip  = 2
   global: do j = 1,m
     !value(j) = value(j)
     if (bndry==2) then
       if     (value(j).lt.xmin) then
         value(j) = value(j)+(xmax-xmin)
       elseif (value(j).gt.xmax) then
         value(j) = value(j)-(xmax-xmin)
       endif
     endif
     !
     if (value(j)==array(n)) then
         idx(j) = n-1
     else
       !
       if (array_order==1) then
         if (j>1) then
           if (value(j)<value(j-1)) ip = 2
         endif
         do i = ip,n
         !do i = 2,n
           if (value(j).ge.array(i-1).and.value(j).lt.array(i)) then
                                     !if (value(j).lt.array(i)) then ! o.k.
             idx(j) = i-1
             ip     = i
             cycle global
           endif
         enddo
       else    ! if (array_order==-1)
         if (j>1) then
           if (value(j)>value(j-1)) ip = 2
         endif
         do i = ip,n
         !do i = 2,n
           if (value(j).le.array(i-1).and.value(j).gt.array(i)) then
                                     !if (value(j).gt.array(i)) then ! o.k.
             idx(j) = i-1
             ip     = i
             cycle global
           endif
         enddo
       endif
       !
     endif
     if (idx(j)==-1) then
       print *,'error: check value range...',ip,n,value(j)
       print *,'                           ',idx
       print *,'     value(1), value(m)    ',value(1),value(m)
       print *,'     array(1), array(n-1:n)',array(1),array(n-1:n)
!       stop
     endif
   enddo global
!
   end subroutine get_idx_in_array
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_idx_in_array_ptr(n, array, val, array_order, bndry) result(idx)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),            intent(in   ) :: n
   real(r8), dimension(n), intent(in   ) :: array
   real(r8),               intent(inout) :: val
   integer(i4),            intent(in   ) :: array_order, bndry
   integer(i4)                           :: idx
! local variables
   integer(i4) :: i
   real(r8)    :: xmin, xmax
!
   if (bndry==2) then
     if (array_order==1) then
       xmin = array(1)
       xmax = array(n)
     else
       xmin = array(n)
       xmax = array(1)
     endif
     if     (val.lt.xmin) then
       val = val+(xmax-xmin)
     elseif (val.gt.xmax) then
       val = val-(xmax-xmin)
     endif
   endif
!
   idx = -1
   if (val==array(n)) then
       idx = n-1
   else
     if (array_order==1) then
       do i = 2,n
         if (val.ge.array(i-1).and.val.lt.array(i)) then
           idx = i-1
           exit
         endif
       enddo
     else    ! if (array_order==-1)
       do i = 2,n
         if (val.le.array(i-1).and.val.gt.array(i)) then
           idx = i-1
           exit
         endif
       enddo
     endif
   endif
   if (idx==-1) then
     print *,'error: check val range...'
   endif
!
   end function get_idx_in_array_ptr
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine linear_equation_solver(n, a, b, c, x, d, e)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(inout) :: a, b, c
   real(r8), dimension(n),   intent(inout) :: x
   real(r8), dimension(n),   optional, intent(inout) :: d, e
! local variables
   integer(i4) :: im    ! 0(my), 1(fast), 2(LAPACK:general), 3(LAPACK:tri)
   integer(i4) :: ndia  ! 3(tri), 5(penta)
   integer(i4) :: i, j, info, kl, ku, ldab
   logical(l4) :: quartic
   integer(i4), dimension(:)  , allocatable :: ipivot
   real(r8)   , dimension(:)  , allocatable :: a1, b1, c1, y
   real(r8)   , dimension(:,:), allocatable :: mat, imat, ab
!
#ifdef USE_LAPACK
   ! im = 0(nomal), 1(fast algorithm), 2(lapack:general), 3(tri and band)
   im = 3
#else
   im = 1
#endif
!
   if (present(d).and.present(e)) then
     ndia = 5
     quartic = .true.
   else
     ndia = 3
     quartic = .false.
   endif
!
   if     (im==0) then
     !
     if (.not.quartic) then
       allocate(imat(n,n))
       call tridiagonal_inverse(n,a,b,c,imat)
       deallocate(imat)
     else
       print *, 'check'
       stop
     endif
     !
   elseif (im==1) then
     !
     allocate(imat(n,n))
     allocate(y(n))
     if (.not.quartic) then
       call fast_tridiagonal_inverse(n,a,b,c,imat)
     else
       call fast_pentadiagonal_inverse(n,a,b,c,d,e,imat)
     endif
     y(:) = x(:)
     x(:) = 0.0_r8
     do j = 1,n
       do i = 1,n
         x(i) = x(i)+imat(i,j)*y(j)
       enddo
     enddo
     deallocate(imat)
     deallocate(y)
     !
   elseif (im==2) then
     !
     allocate(imat(n,n),ipivot(n))
     imat = 0.0
     if (.not.quartic) then
       call get_diagonal(n,a,b,c,imat)
     else
       call get_diagonal(n,a,b,c,imat,d,e)
     endif
     call dgesv(n,1,imat,n,ipivot,x,n,info)
     deallocate(imat,ipivot)
     !
   elseif (im==3) then
     !
     if (.not.quartic) then
       !
       allocate(a1(n),b1(n-1),c1(n-1))
       a1 = a; b1 = b(1:n-1); c1 = c(1:n-1)
       call dgtsv(n,1,b1,a1,c1,x,n,info)
       deallocate(a1,b1,c1)
       !
     else ! if ndia==5
       !
       kl = 2; ku = 2
       ldab = 2*kl+ku+1
       allocate(ipivot(n),ab(ldab,n))
       ! ref: LAPACK manual
       do j=1,n
         do i=max(1,j-ku),min(n,j+kl)
           if     (i==j-2) then
             ab(kl+ku+1+i-j,j) = e(i)
           elseif (i==j-1) then
             ab(kl+ku+1+i-j,j) = c(i)
           elseif (i==j) then
             ab(kl+ku+1+i-j,j) = a(j)
           elseif (i==j+1) then
             ab(kl+ku+1+i-j,j) = b(j)
           elseif (i==j+2) then
             ab(kl+ku+1+i-j,j) = d(j)
           endif
         enddo
       enddo
       !==ivalent to
       !ab(1,:)     = 0.0_r8
       !ab(2,:)     = 0.0_r8
       !ab(3,3:n)   = e(1:n-2)
       !ab(4,2:n)   = c(1:n-1)
       !ab(5,:)     = a(:)
       !ab(6,1:n-1) = b(1:n-1)
       !ab(7,1:n-2) = d(1:n-2)
       !
       call dgbsv(n,kl,ku,1,ab,ldab,ipivot,x,n,info)
       deallocate(ipivot,ab)
       !
     endif ! if ndia == 5
     !
   else ! if im == 3
     print*,'check option'
     stop
   endif
!
   end subroutine linear_equation_solver
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine inverse_matrix_to_linear_equation(n,a,b,c,x,d,e)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                      intent(in   ) :: n
   real(r8), dimension(n),           intent(in   ) :: a, b, c
   real(r8), dimension(n),           intent(inout) :: x
   real(r8), dimension(n), optional, intent(in   ) :: d, e
! local variables
   integer(i4) :: i, j
   logical(l4) :: quartic
   real(r8), allocatable, dimension(:,:) :: imat
   integer(i4) :: info
   real(r8), allocatable, dimension(:)   :: ipivot, work
!
   if (present(d).and.present(e)) then
     quartic = .true.
   else
     quartic = .false.
   endif
!
   allocate(imat(n,n))
   ! generate mat
   allocate(ipivot(n),work(n))
   imat = 0.0_r8
   if (.not.quartic) then
     call get_diagonal(n,a,b,c,imat)
   else
     call get_diagonal(n,a,b,c,imat,d,e)
   endif
#if 0
   ! inverse mat
   call dgetrf(n,n,imat,n,ipivot,info)
   if (info.ne.0) then
     print*,'lapack error: dgetrf'
     stop
   endif
   call dgetri(n,imat,n,ipivot,work,n,info)
   if (info.ne.0) then
     print*,'lapack error: dgetri'
     stop
   endif
   deallocate(ipivot,work)
   ! cal. result
   x = 0.0
   do j=1,n
     do i=1,n
       x(i) = x(i)+imat(i,j)*y(j)
     enddo
   enddo
   !
#else
   call dgesv(n,1,imat,n,ipivot,x,n,info)
#endif
   deallocate(imat)
!
   return
   end subroutine inverse_matrix_to_linear_equation
!-------------------------------------------------------------------------------
!
!   --> j
! | a1, c1, e1, 00, 00
! v b1, a2, c2, e2, 00
! i d1, b2, a3, c3, e3
!   00, d2, b3, a4, c4
!   00, 00, d3, b4, a5
!
!-------------------------------------------------------------------------------
   subroutine get_diagonal(n, a, b, c, mat, d, e)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                      intent(in   ) :: n
   real(r8), dimension(n),           intent(in   ) :: a, b, c
   real(r8), dimension(n,n),         intent(inout) :: mat
   real(r8), dimension(n), optional, intent(in   ) :: d, e
! local variables
   integer(i4) :: i, j
   logical(l4) :: penta
!
   if (present(d).and.present(e)) then
     penta = .true.
   else
     penta = .false.
   endif
!
   mat = 0.0_r8
   if (.not.penta) then
     do i=1,n
       if     (i==1) then
         mat(  i,i) = a(i)
         mat(i+1,i) = b(i)
       elseif (i==n) then
         mat(i-1,i) = c(i-1)
         mat(  i,i) = a(i)
       else
         mat(i-1,i) = c(i-1)
         mat(  i,i) = a(i)
         mat(i+1,i) = b(i)
       endif
     enddo
   else
     do i = 1,n
       if     (i==1) then
         mat(  i,i) = a(i)
         mat(i+1,i) = b(i)
         mat(i+2,i) = d(i)
       elseif (i==2) then
         mat(i-1,i) = c(i-1)
         mat(i  ,i) = a(i)
         mat(i+1,i) = b(i)
         mat(i+2,i) = d(i)
       elseif (i==n) then
         mat(i-2,i) = e(i-2)
         mat(i-1,i) = c(i-1)
         mat(  i,i) = a(i)
       elseif (i==n-1) then
         mat(i-2,i) = e(i-2)
         mat(i-1,i) = c(i-1)
         mat(  i,i) = a(i)
         mat(i+1,i) = b(i)
       else
         mat(i-2,i) = e(i-2)
         mat(i-1,i) = c(i-1)
         mat(  i,i) = a(i)
         mat(i+1,i) = b(i)
         mat(i+2,i) = d(i)
       endif
     enddo
   endif
!
   return
   end subroutine get_diagonal
!-------------------------------------------------------------------------------
!
! algorithem matrix
!   --> j
! | a1, c1, 00, 00
! v b1, a2, c2, 00
! i 00, b2, a3, c3
!   00, 00, b3, a4
!-------------------------------------------------------------------------------
   subroutine fast_tridiagonal_inverse(n, a, b, c, CC)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(in   ) :: a, b
   real(r8), dimension(n),   intent(inout) :: c
   real(r8), dimension(n,n), intent(inout) :: CC
! local variables
   integer(i4) :: i, j, k
   !real(r8)    :: CC(n,n), EE(n,n), AA(0:n)
   real(r8)    :: EE(n,n), AA(0:n)
!
   c(n) = 1.0_r8
! calculate CC(:,1:n)
   ! cal AA(0:n)
   AA(0) = 1.0_r8
   AA(1) = -a(1)*AA(0)/c(1)
   do i = 1,n-1
     !AA(i+1) = -(b(i+1)*AA(i-1)+a(i+1)*AA(i))/c(i+1)
     AA(i+1) = -(b(i)*AA(i-1)+a(i+1)*AA(i))/c(i+1)
   enddo
   do i = 1,n
     CC(i,n) = -AA(i-1)/AA(n)
   enddo
! calculate CC(:,1:n-1)
   EE(:,:) = 0.0_r8
   do i = 1,n
     EE(i,i) = 1.0_r8
   enddo
   CC(:,n-1) = 1.0_r8/c(n-1)*(EE(:,n)-a(n)*CC(:,n))
   do j = n-1,2,-1
     !CC(:,j-1) = 1.0/c(j-1)*(EE(:,j)-a(j)*CC(:,j)-b(j+1)*CC(:,j+1))
     CC(:,j-1) = 1.0_r8/c(j-1)*(EE(:,j)-a(j)*CC(:,j)-b(j)*CC(:,j+1))
   enddo
!
   return
   end subroutine fast_tridiagonal_inverse
!-------------------------------------------------------------------------------
!
! algorithem matrix
!   --> j
! | a1, c1, e1, 00, 00
! v b1, a2, c2, e2, 00
! i d1, b2, a3, c3, e3
!   00, d2, b3, a4, c4
!   00, 00, d3, b4, a5
!-------------------------------------------------------------------------------
   subroutine fast_pentadiagonal_inverse(n, a, b, c, d, e, CC)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n),   intent(in   ) :: a, b, d
   real(r8), dimension(n),   intent(inout) :: c, e
   real(r8), dimension(n,n), intent(inout) :: CC
! local variables
   integer(i4) :: i, j, k
   !real(r8)    :: CC(n,n), EE(n,n), AA(0:n+1), BB(0:n+1), XX(0:n+1), YY(0:n+1)
   !real(r8)    :: EE(n,n), AA(0:n+1), BB(0:n+1), XX(0:n+1), YY(0:n+1)
   real(r8), dimension(:,:), allocatable :: EE
   real(r8), dimension(:),   allocatable :: AA, BB, XX, YY
!
   allocate(EE(n,n),AA(0:n+1), BB(0:n+1), XX(0:n+1), YY(0:n+1))
!
   e(n-1) = 1.0_r8
   e(n)   = 1.0_r8
   c(n)   = 0.0_r8
! calculate CC(:,1:n)
   ! cal AA(0:n)
   AA(0) = 0.0_r8
   AA(1) = 1.0_r8
   AA(2) = -( a(1)*AA(0)+c(1)*AA(1) )/e(1)
   !AA(3) = -( b(2)*AA(0)+a(2)*AA(1)+c(2)*AA(2) )/e(2)
   AA(3) = -( b(1)*AA(0)+a(2)*AA(1)+c(2)*AA(2) )/e(2)
   do i = 3,n
     !AA(i+1) = -(d(i)*AA(i-3)+b(i)*AA(i-2)+a(i)*AA(i-1)+c(i)*AA(i))/e(i)
     AA(i+1) = -(d(i-2)*AA(i-3)+b(i-1)*AA(i-2)+a(i)*AA(i-1)+c(i)*AA(i))/e(i)
   enddo
   BB(0) = 1.0_r8
   BB(1) = 0.0_r8
   BB(2) = -( a(1)*BB(0)+c(1)*BB(1) )/e(1)
   !BB(3) = -( b(2)*BB(0)+a(2)*BB(1)+c(2)*BB(2) )/e(2)
   BB(3) = -( b(1)*BB(0)+a(2)*BB(1)+c(2)*BB(2) )/e(2)
   do i = 3,n
     !BB(i+1) = -(d(i)*BB(i-3)+b(i)*BB(i-2)+a(i)*BB(i-1)+c(i)*BB(i))/e(i)
     BB(i+1) = -(d(i-2)*BB(i-3)+b(i-1)*BB(i-2)+a(i)*BB(i-1)+c(i)*BB(i))/e(i)
   enddo
   do i = 0,n+1
     XX(i) = AA(n+1)*BB(i)-(AA(i)*BB(n+1))
     YY(i) = AA(n)*BB(i)-(AA(i)*BB(n))
   enddo
!
   do i = 1,n
     CC(i,n)   = -YY(i-1)/YY(n+1)
     CC(i,n-1) = -XX(i-1)/XX(n)
   enddo
! calculate CC(:,1:n-1)
   EE(:,:) = 0.0_r8
   do i = 1,n
     EE(i,i) = 1.0_r8
   enddo
!
   CC(:,n-2) = 1.0_r8/e(n-2)*(EE(:,n)-a(n)*CC(:,n)-c(n-1)*CC(:,n-1))
   CC(:,n-3) = 1.0_r8/e(n-3)*(EE(:,n-1)-c(n-2)*CC(:,n-2)-a(n-1)*CC(:,n-1)-b(n-1)*CC(:,n))
   do j = n-2,3,-1
     !CC(:,j-2) = 1.0_r8/e(j-2)*(EE(:,j)-c(j-1)*CC(:,j-1)-a(j)*CC(:,j)-b(j+1)*CC(:,j+1)-d(j+2)*CC(:,j+2))
     CC(:,j-2) = 1.0_r8/e(j-2)*(EE(:,j)-c(j-1)*CC(:,j-1)-a(j)*CC(:,j)-b(j)*CC(:,j+1)-d(j)*CC(:,j+2))
   enddo
   deallocate(EE,AA,BB,XX,YY)
!
   return
   end subroutine fast_pentadiagonal_inverse
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
   subroutine unittest_remap()
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4)       :: use_spline = .true., reverse = .false.
   type(spline_t)    :: sr
   ! source
   integer(i4) :: n
   real(r8), dimension(:), allocatable :: sx, scx, sdx, sy
   ! target
   integer(i4) :: m
   real(r8), dimension(:), allocatable :: tx, tcx, tdx, ty
   ! function
   integer(i4) :: l
   real(r8), dimension(:), allocatable :: ftx, fty
   ! boundary derivative
   integer(i4) :: order, bndry
   logical     :: ismono
   integer(i4) :: i, j
   real(r8)    :: mass1, mass2, dx
!
   order  = 2 ! 0:pcm, 1:plm, 2:ppm
   bndry  = 2 ! 0:mirror, 1:mirror+derivative, 2:periodic(need base)
   ismono = .false.
!
   n = 10
   m = 9
   l = 11
!
   allocate(sx(n))       ! source x
   allocate(scx(0:n))    ! source x (control volume)
   allocate(sdx(n))      ! source dx
   allocate(sy(n))       ! source y (samples)
   allocate(tx(m))       ! target x
   allocate(tcx(0:m))    ! target x (control volume)
   allocate(tdx(m))      ! target dx
   allocate(ty(m))       ! target y
   allocate(ftx(l))      ! spline function
   allocate(fty(l))      ! spline function
!
#define REVERSE
#ifndef REVERSE
   sdx    = 1.0
   scx(0) = 0.0
   !
   tdx    = (/1.0,1.0,1.0,1.0,5.0,0.25,0.25,0.25,0.25/)
   tcx(0) = 0.0
#else
   sdx    = -1.0
   scx(0) = 10.0
   !
   tdx    = (/-1.0,-1.0,-1.0,-1.0,-5.0,-0.25,-0.25,-0.25,-0.25/)
   tcx(0) = 10.0
#endif
   do i = 1,n
     scx(i) = scx(i-1)+sdx(i)
   enddo
   do i = 1,m
     tcx(i) = tcx(i-1)+tdx(i)
   enddo
   !
   dx = (scx(n)-scx(0))/dble(l-1)
   do i = 1,l
     ftx(i) = scx(0)+dx*dble(i-1)
   enddo
!
   do i = 1,n
     sy(i) = real(i,r8)-0.5
   enddo
!
   tcx = tcx-5.0
   call spline_interface(n,scx,sy,m,tcx,ty,order,ismono,bndry,l,ftx,fty)
   !call spline_cyclic_interface(n,scx,sy,m,tcx,ty,order,ismono,l,ftx,fty)
   !call remap_psm(n,scx,sy,m,tcx,ty,ismono,bndry,l,ftx,fty)
   !call remap_psm_cyclic(n,scx,sy,m,tcx,ty,ismono,l,ftx,fty)
!
! get mass
   print *, 'scx'
   print '(11(f6.2))', scx
   print '(a,10(f6.2))','   ',sy
   print *, 'tcx'
   print '(10(f6.2))', tcx
   print '(a,09(f6.2))','   ',ty
   mass1 = get_mass(n,scx,sy)
   mass2 = get_mass(n,tcx,ty)
   print *,'mass1 = ',mass1
   print *,'mass2 = ',mass2
   print '(11(f6.2))',ftx
   print '(11(f6.2))',fty
!
   deallocate(sx)
   deallocate(scx)
   deallocate(sdx)
   deallocate(sy)
   deallocate(tx)
   deallocate(tcx)
   deallocate(tdx)
   deallocate(ty)
   deallocate(ftx)
   deallocate(fty)
!
!-------------------------------------------------------------------------------
   end subroutine unittest_remap
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine unittest_main()
!-------------------------------------------------------------------------------
   implicit none
!
   type(spline_t) :: sr
   ! source
   integer(i4) :: n
   real(r8), dimension(:), allocatable :: scx, sdx, sy
   ! target
   integer(i4) :: m
   real(r8), dimension(:), allocatable :: tcx, tdx, ty
   ! condition (grid)
   integer(i4), parameter :: norder=1, nshift=4, nreverse=2, nbndry=3, nmono=2
   integer(i4) :: order(norder)
   integer(i4) :: shift(nshift)
   logical(l4) :: reverse(nreverse)
   ! condition (remap)
   integer(i4) :: bndry(nbndry)
   logical(l4) :: mono(nmono)
   integer(i4) :: io, is, ir, ib, im
   !
   ! mass
   real(r8)    :: mass1, mass2
   ! function
   integer(i4)                         :: l
   real(r8), dimension(:), allocatable :: ftx, fty
!
   order   = (/2/) ! 0:pcm, 1:plm, 2:ppm
   shift   = (/0,-1,1,99/)
   reverse = (/.false.,.true./)
   bndry   = (/0,1,2/) ! 0:mirror, 1:mirror+derivative, 2:periodic(need base)
   mono    = (/.false.,.true./)
!
   n = 10
   m = 9
   l = 11
!
   allocate(scx(0:n))    ! source x (control volume)
   allocate(sdx(n))      ! source dx
   allocate(sy(n))       ! source y (samples)
   allocate(tcx(0:m))    ! target x (control volume)
   allocate(tdx(m))      ! target dx
   allocate(ty(m))       ! source y (samples)
   allocate(ftx(l))      ! spline function
   allocate(fty(l))      ! spline function
!
   do io = 1,norder
    do is = 1,nshift
     do ir = 1,nreverse
      do ib = 1,nbndry
       do im = 1,nmono
         if (shift(is)/=0.and.bndry(ib)/=2) then
         else
           call unittest_core(reverse(ir),shift(is),order(io),bndry(ib),mono(im),&
                            n,scx,sy,m,tcx,ty,mass1,mass2)
         endif
       enddo ! im
      enddo ! ib
     enddo ! ir
    enddo ! is
   enddo ! io
!
#if 0
   !call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2,l,ftx,fty)
   !call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
   order = 2
! ppm,          grid(F, 0), bndry(0), mono(F)
   reverse=.false.;shift= 0;bndry=0;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(F, 0), bndry(0), mono(T)
   reverse=.false.;shift= 0;bndry=0;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(F, 0), bndry(1), mono(F)
   reverse=.false.;shift= 0;bndry=1;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(F, 0), bndry(1), mono(T)
   reverse=.false.;shift= 0;bndry=1;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(F, 0), bndry(2), mono(F)
   reverse=.false.;shift= 0;bndry=2;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(F, 0), bndry(2), mono(T)
   reverse=.false.;shift= 0;bndry=2;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T, 0), bndry(0), mono(F)
   reverse=.true. ;shift= 0;bndry=0;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T, 0), bndry(0), mono(T)
   reverse=.true. ;shift= 0;bndry=0;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T, 0), bndry(1), mono(F)
   reverse=.true. ;shift= 0;bndry=1;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T, 0), bndry(1), mono(T)
   reverse=.true. ;shift= 0;bndry=1;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T, 0), bndry(2), mono(F)
   reverse=.true. ;shift= 0;bndry=2;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T, 0), bndry(2), mono(T)
   reverse=.true. ;shift= 0;bndry=2;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
!
! ppm,          grid(F,-1), bndry(2), mono(F)
   reverse=.false.;shift=-1;bndry=2;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(F,-1), bndry(2), mono(T)
   reverse=.false.;shift=-1;bndry=2;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(F, 1), bndry(2), mono(F)
   reverse=.false.;shift= 1;bndry=2;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(F, 1), bndry(2), mono(T)
   reverse=.false.;shift= 1;bndry=2;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(F,99), bndry(2), mono(F)
   reverse=.false.;shift=99;bndry=2;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(F,99), bndry(2), mono(T)
   reverse=.false.;shift=99;bndry=2;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T,-1), bndry(2), mono(F)
   reverse=.true. ;shift=-1;bndry=2;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T,-1), bndry(2), mono(T)
   reverse=.true. ;shift=-1;bndry=2;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T, 1), bndry(2), mono(F)
   reverse=.true. ;shift= 1;bndry=2;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T, 1), bndry(2), mono(T)
   reverse=.true. ;shift= 1;bndry=2;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T,99), bndry(2), mono(F)
   reverse=.true. ;shift=99;bndry=2;mono=.false.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
! ppm,          grid(T,99), bndry(2), mono(T)
   reverse=.true. ;shift=99;bndry=2;mono=.true.
   call unittest_core(reverse,shift,order,bndry,mono,n,scx,sy,m,tcx,ty,mass1,mass2)
#endif
!
   deallocate(scx)
   deallocate(sdx)
   deallocate(sy)
   deallocate(tcx)
   deallocate(tdx)
   deallocate(ty)
   deallocate(ftx)
   deallocate(fty)
!
!-------------------------------------------------------------------------------
   end subroutine unittest_main
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine unittest_core(reverse, shift, order, bndry, ismono,              &
                            n, scx, sy, m, tcx, ty, mass1, mass2, l, ftx, fty)
!-------------------------------------------------------------------------------
   implicit none
   logical(l4),                 intent(in   ) :: reverse
   integer(i4),                 intent(in   ) :: shift
   integer(i4),                 intent(in   ) :: order, bndry
   logical(l4),                 intent(in   ) :: ismono
   integer(i4),                 intent(in   ) :: n
   real(r8)   , dimension(0:n), intent(inout) :: scx
   real(r8)   , dimension(n)  , intent(inout) :: sy
   integer(i4),                 intent(in   ) :: m
   real(r8)   , dimension(0:m), intent(inout) :: tcx
   real(r8)   , dimension(m),   intent(inout) :: ty
   real(r8)   ,                 intent(  out) :: mass1, mass2
   integer(i4),                 intent(in   ), optional :: l
   real(r8)   , dimension(:),   intent(inout), optional :: ftx, fty
!
   logical(l4)    :: is_fun, check
   type(spline_t) :: sr
   integer(i4)    :: i
   real(r8)       :: dx
!
! print
   print '(a)', '---------------------------------------------------------------------'
   print '(a,i2,a,l,a,i2,a,i2,a,l,a)','order = (',order,'), reverse = (',reverse, &
                    '), shift = (',shift,'), bndry = (',bndry,'), ismono = (',ismono,')'
!
   call unittest_gen_grid(reverse,shift,n,scx,sy,m,tcx)
!
   if (present(l).and.present(ftx).and.present(fty)) then
     is_fun = .true.
   else
     is_fun = .false.
   endif
!
   if (is_fun) then
     dx = (scx(n)-scx(0))/dble(l-1)
     do i = 1,l
       ftx(i) = scx(0)+dx*dble(i-1)
     enddo
   endif
!
   if (is_fun) then
     call spline_interface(n,scx,sy,m,tcx,ty,order,ismono,bndry,l,ftx,fty)
   else
     call spline_interface(n,scx,sy,m,tcx,ty,order,ismono,bndry)
   endif
!
! get mass
   mass1 = get_mass(n,scx,sy)
   mass2 = get_mass(n,tcx,ty)
!
   check = .true.
   if (shift==0.or.shift==-1.or.shift==1) then
     if(abs(mass2-mass1)>1.0d-14) check = .false.
   else
     if(abs(mass2-2.0_r8*mass1)>1.0d-14) check = .false.
   endif
! print
   if (.not.check) then
     print '(a)','scx'
     print '(11(f6.2))',scx
     print '(a,10(f6.2))','   ',sy
     print '(a)', 'tcx'
     print '(10(f6.2))',tcx
     print '(a,09(f6.2))','   ',ty
     print '(a,f15.11,x,f15.11)','mass = ',mass1,mass2
     if (is_fun) then
       print '(11(f6.2))',ftx
       print '(11(f6.2))',fty
     endif
   else
     print '(a)','o.k.'
   endif
   print '(a)', '---------------------------------------------------------------------'
   print *, ' '
!
!-------------------------------------------------------------------------------
   end subroutine unittest_core
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine unittest_gen_grid(reverse, shift, n, scx, sy, m, tcx)
!-------------------------------------------------------------------------------
   implicit none
   logical(l4),                 intent(in   ) :: reverse
   integer(i4),                 intent(in   ) :: shift
   integer(i4),                 intent(in   ) :: n
   real(r8)   , dimension(0:n), intent(inout) :: scx
   real(r8)   , dimension(n)  , intent(inout) :: sy
   integer(i4),                 intent(in   ) :: m
   real(r8)   , dimension(0:m), intent(inout) :: tcx
!
   real(r8), dimension(:), allocatable :: sdx
   real(r8), dimension(:), allocatable :: tdx
   integer(i4) :: i, j
!
   allocate(sdx(n))      ! source dx
   allocate(tdx(m))      ! target dx
!
   if (.not.reverse) then
     sdx(:) = 1.0
     scx(0) = 0.0
     tdx    = (/1.0,1.0,1.0,1.0,5.0,0.25,0.25,0.25,0.25/)
     tcx(0) = 0.0
   else
     sdx(:) = -1.0
     scx(0) = 10.0
     tdx    = (/-1.0,-1.0,-1.0,-1.0,-5.0,-0.25,-0.25,-0.25,-0.25/)
     tcx(0) = 10.0
   endif
   do i = 1,n
     scx(i) = scx(i-1)+sdx(i)
   enddo
   do i = 1,m
     tcx(i) = tcx(i-1)+tdx(i)
   enddo
   if     (shift== 0) then
     !tcx(:) = tcx(:)
   elseif (shift==-1) then
     tcx(:) = tcx(:)-0.25_r8
   elseif (shift== 1) then
     tcx(:) = tcx(:)+0.25_r8
   else
     if (.not.reverse) then
       tcx(0) = tcx(0)-5.0_r8
       tcx(m) = tcx(m)+5.0_r8
     else
       tcx(0) = tcx(0)+5.0_r8
       tcx(m) = tcx(m)-5.0_r8
     endif
   endif
!
   do i = 1,n
     sy(i) = real(i,r8)-0.5_r8
   enddo
   if (reverse) sy(:) = -sy(:)
!
   deallocate(sdx)      ! source dx
   deallocate(tdx)      ! target dx
!
!-------------------------------------------------------------------------------
   end subroutine unittest_gen_grid
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine remap_psm(n,scx,sy,m,tcx,ty,ismono,bndry,l,ftx,fty)
!-------------------------------------------------------------------------------
   implicit none
   integer(i4), intent(in   ) :: n
   real(r8)   , intent(in   ) :: scx(n+1)
   real(r8)   , intent(in   ) :: sy(n)
   integer(i4), intent(in   ) :: m
   real(r8)   , intent(in   ) :: tcx(m+1)
   real(r8)   , intent(  out) :: ty(m)
   logical(i4),           intent(in   ) :: ismono
   integer(i4),           intent(in   ) :: bndry
   integer(i4), optional, intent(in   ) :: l
   real(r8)   , optional, intent(in   ) :: ftx(:)
   real(r8)   , optional, intent(  out) :: fty(:)
! local variables
   real(r8), dimension(n+1)   :: scy,lower_diag,diag,upper_diag,q_diag,inty
   real(r8), dimension(n)     :: sx,dx,rdx,sydx,dy,scyr,za0,za1,za2
   real(r8), dimension(n,0:2) :: arg
   real(r8), dimension(m+1)   :: xi
   real(r8)    :: tmp_cal,inty1,inty2,tmp
   integer(i4) :: zkr(m+1),do_filter(n),i,ilev,jk,k,q,array_order
   logical(i4) :: abort = .false.
   ! for function 
   integer(i4), dimension(:), allocatable :: fzkr
   real(r8)   , dimension(:), allocatable :: fxi
!
!   if (tcx(1).lt.scx(1).or.tcx(m+1).gt.scx(n+1)) then
!     write(6,*) 'CHECK TARGET BOUNDARY'
!     abort=.true.
!   endif
! 
   array_order = check_array_order(scx)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! quadratic splies with UK met office monotonicity constraints  !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   dx   = scx(2:n+1)-scx(1:n)
   sx   = 0.5*(scx(2:n+1)+scx(1:n))
   rdx  = 1/dx
   sydx = sy*dx
!
   if (bndry==-1) then
     scy(1)   = sy(1)
     scy(n+1) = sy(n)
   else if (bndry==0) then
     scy(1)   = -1.0*(sy(2)-sy(1))/(sx(2)-sx(1))*(scx(1)-sx(1))+sy(1)
     scy(n+1) = -1.0*(sy(n)-sy(n-1))/(sx(n)-sx(n-1))*(scx(n+1)-sx(n-1))+sy(n-1)
   else if (bndry==1) then
     scy(1)   = (sy(2)-sy(1))/(sx(2)-sx(1))*(scx(1)-sx(1))+sy(1)
     scy(n+1) = (sy(n)-sy(n-1))/(sx(n)-sx(n-1))*(scx(n+1)-sx(n-1))+sy(n-1)
   else if (bndry==2) then
     tmp      = scx(n+1)+(sx(1)-scx(1))
     scy(n+1) = (sy(1)-sy(n))/(tmp-sx(n))*(tmp-scx(n+1))+sy(n)
     scy(1)   = scy(n+1)
   endif
   scy(1)   = 3.0*scy(1)
   scy(n+1) = 3.0*scy(n+1)
!
   scy(2:n) = 3*(sy(2:n)*rdx(2:n)+sy(1:n-1)*rdx(1:n-1))
! 
   lower_diag      = 0
   diag            = 0
   upper_diag      = 0
   lower_diag(1)   = 1
   lower_diag(2:n) = rdx(1:n-1)
   lower_diag(n+1) = 1
   diag(1)   = 2
   diag(2:n) = 2*(rdx(2:n)+rdx(1:n-1))
   diag(n+1) = 2
   upper_diag(1)   = 1
   upper_diag(2:n) = rdx(2:n)
   upper_diag(n+1) = 0
   q_diag(1)  = -upper_diag(1)/diag(1)
   scy(1)     = scy(1)/diag(1)
! 
   do k=2,n+1
     tmp_cal   =  1/(diag(k)+lower_diag(k)*q_diag(k-1))
     q_diag(k) = -upper_diag(k)*tmp_cal
     scy(k)    =  (scy(k)-lower_diag(k)*scy(k-1))*tmp_cal
   enddo
   do k=n,1,-1
     scy(k) = scy(k)+q_diag(k)*scy(k+1)
   enddo
! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  monotonicity modifications  !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! grid scale
   if (ismono) then
     call monotone_grid_scale(n,sy,scy,do_filter)
   endif
!
   za0 =    scy(1:n)
   za1 = -4*scy(1:n)-2*scy(2:n+1)+6*sy
   za2 =  3*scy(1:n)+3*scy(2:n+1)-6*sy
!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Compute the 3 quadratic spline coeffients {za0, za1, za2}                 !!
   !! knowing the quadratic spline parameters {rho_left,rho_right,sy}         !!
   !! Zerroukat et.al., Q.J.R. Meteorol. Soc., Vol. 128, pp. 2801-2820 (2002).   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! subgrid scale
   if (ismono) then
     !
     arg(:,0) = za0(:)
     arg(:,1) = za1(:)
     arg(:,2) = za2(:)
     call monotone_subgrid_scale_psm(n,sy,scy,arg,do_filter,array_order)
     za0(:) = arg(:,0)
     za1(:) = arg(:,1)
     za2(:) = arg(:,2)
     !
   endif ! if (ismono) then
! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! start iteration from top to bottom of atmosphere !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!#define PROFILE_CAM
#ifdef PROFILE_CAM
time1 = mpi_wtime()
#endif
   ! get index of source for target
   zkr      = -1
   ilev     = 2
   zkr(1)   = 1
   zkr(m+1) = n
   if (array_order==1) then
     kloop: do k = 1,m+1
       do jk = ilev,n+1
         if (scx(jk).ge.tcx(k)) then
           ilev   = jk
           zkr(k) = jk-1
           cycle kloop
         endif
       enddo
     enddo kloop
   else
     kloop1: do k = 1,m+1
       do jk = ilev,n+1
         if (scx(jk).le.tcx(k)) then
           ilev   = jk
           zkr(k) = jk-1
           cycle kloop1
         endif
       enddo
     enddo kloop1
   endif
#ifdef PROFILE_CAM
time2 = mpi_wtime()
print '(a30,f11.8)','cam_psm_get_idx = ',time2-time1
time1 = mpi_wtime()
#endif
! 
!#if 0
   ! int{f(x)} at source cells 
   inty(1) = 0
   do k = 1,n
     inty(k+1) = inty(k)+sydx(k)
   enddo
!
   ! xi
   xi = (tcx(1:m+1)-scx(zkr))/(scx(zkr+1)-scx(zkr))
! 
   !inty1 = 0
   inty1 = inty(zkr(1))+(za0(zkr(1))*xi(1)+(za1(zkr(1))/2)*(xi(1)**2)+ &
                                       (za2(zkr(1))/3)*(xi(1)**3))*dx(zkr(1))
   do k=1,m
     if (xi(k+1)>1.0_r8) then
       write(*,*) 'r not in [0:1]', xi(k+1)
       abort = .true.
     endif
     inty2 = inty(zkr(k+1))+(za0(zkr(k+1))*xi(k+1)+(za1(zkr(k+1))/2)*(xi(k+1)**2)+ &
                                       (za2(zkr(k+1))/3)*(xi(k+1)**3))*dx(zkr(k+1))
     ty(k) = (inty2-inty1)
     inty1 = inty2
   enddo
!#endif
!
#ifdef PROFILE_CAM
time2 = mpi_wtime()
print '(a30,f11.8)','cam_psm_get_int = ',time2-time1
#endif
!
   do k = 1,m
      ty(k) = ty(k)/(tcx(k+1)-tcx(k))
   enddo
!
   if (abort) then
     print *,'bad levels in remap1.  usually cfl violatioin'
     stop
   endif
!
! out remap function
!..............
   if (present(l).and.present(ftx).and.present(fty)) then
     allocate(fzkr(l),fxi(l))
     !
     fzkr    = 1000000
     ilev    = 2
     fzkr(1) = 1
     fzkr(l) = n
     if (array_order==1) then
       fkloop: do k = 1,l
         do jk = ilev,n+1
           if (scx(jk).ge.ftx(k)) then
             ilev    = jk
             fzkr(k) = jk-1
             cycle fkloop
           endif
         enddo
       enddo fkloop
     else
       fkloop1: do k = 1,l
         do jk = ilev,n+1
           if (scx(jk).le.ftx(k)) then
             ilev    = jk
             fzkr(k) = jk-1
             cycle fkloop1
           endif
         enddo
       enddo fkloop1
     endif
     ! 
     ! xi
     fxi = (ftx(1:l)-scx(fzkr))/(scx(fzkr+1)-scx(fzkr))
     !
     do k = 1,l
       fty(k) = za0(fzkr(k))+za1(fzkr(k))*fxi(k)+za2(fzkr(k))*fxi(k)**2.0
     enddo
!
     deallocate(fzkr,fxi)
   endif
!..............
!
!-------------------------------------------------------------------------------
   end subroutine remap_psm
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine remap_psm_cyclic(n, scx, sy, m, tcx, ty, ismono, l, ftx, fty)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(0:n), intent(in   ) :: scx
   real(r8), dimension(n),   intent(in   ) :: sy
   integer(i4),              intent(in   ) :: m
   real(r8), dimension(0:m), intent(in   ) :: tcx
   real(r8), dimension(m),   intent(  out) :: ty
   logical(i4),              intent(in   ) :: ismono
   integer(i4), optional   , intent(in   ) :: l
   real(r8)   , optional   , intent(in   ) :: ftx(:)
   real(r8)   , optional   , intent(  out) :: fty(:)
! local variables
   type(spline_t) :: sr
   real(r8)       :: vmin, lbase, rbase
   integer(i4)    :: i, nn, array_order
   real(r8), allocatable, dimension(:) :: escx, esy
   ! function
   integer(i4), dimension(:), allocatable :: fzkr
   real(r8)   , dimension(:), allocatable :: fxi
!
   nn = 3*n
   allocate(escx(0:nn),esy(nn))
!
! expand x (source)
   array_order = check_array_order(scx)
   if (array_order==0) then
     print*,'check array order of source in spline_cyclic_interface.'
     stop
   endif
   rbase           = scx(n)-scx(0)
   lbase           = scx(0)-scx(n)
   escx(n:2*n)     = scx(:)
   escx(2*n+1:3*n) = scx(1:n)+rbase   ! right
   escx(0:n-1)     = scx(0:n-1)+lbase ! left
!
! expand y (source)
   esy(n+1:2*n)   = sy(:)
   esy(1:n)       = sy(1:n)
   esy(2*n+1:3*n) = sy(1:n)
!
   if (present(l).and.present(ftx).and.present(fty)) then
     call remap_psm(nn,escx,esy,m,tcx,ty,ismono,2,l,ftx,fty)
   else
     call remap_psm(nn,escx,esy,m,tcx,ty,ismono,2)
   endif
!
   deallocate(escx,esy)
!
   return
   end subroutine remap_psm_cyclic
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine test_remap_psm()
!-------------------------------------------------------------------------------
   implicit none
! local variables
   integer(i4) :: n, m, l, i, j
   real(r8), dimension(:), allocatable :: sdx, scx, sy
   real(r8), dimension(:), allocatable :: tdx, tcx, ty
   real(r8), dimension(:), allocatable :: ftx, fty
!
   n = 10
   m = 9
   l = 10
   allocate(sdx(n),scx(n+1),sy(n))
   allocate(tdx(m),tcx(m+1),ty(m))
   allocate(ftx(l),fty(l))
!#define REVERSE
!
#ifndef REVERSE
   sdx    = 1.0
   scx(1) = 0.0
   !
   tdx    = (/1.0,1.0,1.0,1.0,5.0,0.25,0.25,0.25,0.25/)
   tcx(1) = 0.0
#else
   sdx    = -1.0
   scx(1) = 10.0
   !
   tdx    = (/-1.0,-1.0,-1.0,-1.0,-5.0,-0.25,-0.25,-0.25,-0.25/)
   tcx(1) = 10.0
#endif
   do i = 1,n
     scx(i+1) = scx(i)+sdx(i)
   enddo
   do i = 1,m
     tcx(i+1) = tcx(i)+tdx(i)
   enddo
!
   do i = 1,n
     sy(i)    = real(i,r8)-0.5
   enddo
!
   ! function
   do i = 1,l
     ftx(i) = 0.5*(scx(i)+scx(i+1))
   enddo
!
   print *, 'scx'
   print '(11(f6.2))', scx
!
   print '(a,10(f6.2))','   ',sy
   print *, 'tcx'
   print '(11(f6.2))', tcx
   !call remap_psm(n,scx,sy,m,tcx,ty,ismono=.false.,bndry=1,l=l,ftx=ftx,fty=fty) 
   call remap_psm(n,scx,sy,m,tcx,ty,ismono=.true.,bndry=-1,l=l,ftx=ftx,fty=fty) 
   print '(a,10(f6.2))','   ',ty
!   print '(a,10(f6.2))','   ',fty
!
   deallocate(sdx,scx,sy)
   deallocate(tdx,tcx,ty)
   deallocate(ftx,fty)
!
!-------------------------------------------------------------------------------
   end subroutine test_remap_psm
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module spline_remap
!-------------------------------------------------------------------------------
