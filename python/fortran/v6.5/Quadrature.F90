:-3:00:#include "KIM.h"
:-3:00:#undef _GAUSS_TABLE
:-3:00:#undef _QUAD_DBG
!-------------------------------------------------------------------------------
:01:01:   module quadrature
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  your   name    code comment
!
!  structure:
!
!-------------------------------------------------------------------------------
:03:03:   use kiapsbase, only: iulog=>kim_iu_log, real_kind=>kim_real8_kind, longdouble_kind=>kim_longdouble_kind
:03:03:   implicit none
!
:04:04:   private
!
:05:05:   type, public :: quadrature_t
:04:04:     real(longdouble_kind), dimension(:), pointer :: points
:04:04:     real(longdouble_kind), dimension(:), pointer :: weights
:-1:04:   end type quadrature_t
:04:04:   public :: gausslobatto
:04:04:   public :: test_gausslobatto
:04:04:   public :: gauss
:04:04:   public :: test_gauss
:04:04:   public :: legendre
:04:04:   public :: jacobi
:04:04:   public :: quad_norm
:04:04:   public :: trapezoid
:04:04:   private :: trapn
:04:04:   public :: simpsons
:04:04:   public :: gaussian_int
:04:04:   private :: gausslobatto_pts
:04:04:   private :: gausslobatto_wts
:04:04:   private :: gauss_pts
:04:04:   private :: gauss_wts
:04:04:   private :: jacobi_polynomials
:04:04:   private :: jacobi_derivatives
!
:-1:01:   contains
:-3:01:! ==============================================================
:-3:01:! gauss:
:-3:01:!
:-3:01:! Find the Gauss collocation points and the corresponding weights.
:-3:01:!
:-3:01:! ==============================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function gauss(npts) result(gs)
!
:04:04:   integer, intent(in   ) :: npts
:04:04:   type(quadrature_t) :: gs
!
:06:06:   allocate(gs%points(npts))
:06:06:   allocate(gs%weights(npts))
:06:06:   gs%points = gauss_pts(npts)
:06:06:   gs%weights = gauss_wts(npts,gs%points)
!
:-1:01:   end function gauss
:-3:01:#if defined(_GAUSS_TABLE)
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function gauss_pts(npts) result(pts)
!
:04:04:   integer, intent(in   ) :: npts
:04:04:   real(longdouble_kind) :: pts(npts)
!
:06:06:   pts(1) =-0.93246951420315202781d0
:06:06:   pts(2) =-0.66120938646626451366d0
:06:06:   pts(3) =-0.23861918608319690863d0
:06:06:   pts(4) =-pts(3)
:06:06:   pts(5) =-pts(2)
:06:06:   pts(6) =-pts(1)
!
:-1:01:   end function gauss_pts
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function gauss_wts(npts, pts) result(wts)
!
:04:04:   integer, intent(in   ) :: npts
:04:04:   real(longdouble_kind) :: pts(npts)
:04:04:   real(longdouble_kind) :: wts(npts)
!
:06:06:   wts(1) = 0.17132449237917034504d0
:06:06:   wts(2) = 0.36076157304813860756d0
:06:06:   wts(3) = 0.46791393457269104738d0
:06:06:   wts(4) = wts(3)
:06:06:   wts(5) = wts(2)
:06:06:   wts(6) = wts(1)
!
:-1:01:   end function gauss_wts
:-3:01:#else
:-3:01:! ==============================================================
:-3:01:! gauss_pts:
:-3:01:!
:-3:01:! Compute the Gauss Collocation points
:-3:01:! for Jacobi Polynomials
:-3:01:!
:-3:01:! ==============================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function gauss_pts(np1) result(pts)
!-------------------------------------------------------------------------------
:03:03:   use physicalconstants, only: qq_pi
!
:04:04:   integer, intent(in   ) :: np1 ! number of velocity grid points
:04:04:   real(longdouble_kind) :: pts(np1)
:-3:04:! Local variables
:04:04:   real(longdouble_kind) :: alpha, beta
:04:04:   real(longdouble_kind) :: xjac(0:np1-1)
:04:04:   real(longdouble_kind) :: jac(0:np1)
:04:04:   real(longdouble_kind) :: djac(0:np1)
:04:04:   integer prec ! number of mantissa bits
:04:04:   real(longdouble_kind) eps ! machine epsilon
:04:04:   real(longdouble_kind), parameter :: convthresh = 10 ! convergence threshold relative\
:-3:04:! to machine epsilon
:04:04:   integer, parameter :: kstop = 30 ! max iterations for polynomial deflation
:04:04:   real(longdouble_kind) :: poly
:04:04:   real(longdouble_kind) :: pder
:04:04:   real(longdouble_kind) :: recsum, thresh
:04:04:   real(longdouble_kind) :: dth
:04:04:   real(longdouble_kind) :: x
:04:04:   real(longdouble_kind) :: delx
:04:04:   real(longdouble_kind) :: c0, c1, c2, c10
:04:04:   integer i, j, k
:04:04:   integer n, nh
!
:06:06:   n = np1-1
:06:06:   c0 = 0.0_longdouble_kind
:06:06:   c1 = 1.0_longdouble_kind
:06:06:   c2 = 2.0_longdouble_kind
:06:06:   c10 = 10.0_longdouble_kind
:06:06:   alpha = c0
:06:06:   beta = c0
:-3:06:! =========================================================
:-3:06:! compute machine precision and set the convergence
:-3:06:! threshold thresh to 10 times that level
:-3:06:! =========================================================
:06:06:   prec = precision(c10)
:06:06:   eps = c10**(-prec)
:06:06:   thresh = convthresh*eps
:-3:06:! ============================================================
:-3:06:! Compute first half of the roots by "polynomial deflation".
:-3:06:! ============================================================
:06:06:   dth = qq_pi/(2*n+2)
:06:06:   nh =(n+1)/2
:05:05:   do j = 0,nh-1
:06:06:     x = cos((c2*j+1)*dth) ! first guess at root
:06:06:     k = 0
:06:06:     delx = c1
:05:05:     do while(k<kstop.and.abs(delx)>thresh)
:06:06:       call jacobi(n+1,x,alpha,beta,jac(0:n+1),djac(0:n+1))
:06:06:       poly = jac(n+1)
:06:06:       pder = djac(n+1)
:06:06:       recsum = c0
:05:05:       do i = 0,j-1
:06:06:         recsum = recsum+c1/(x-xjac(i))
:-1:06:       enddo
:06:06:       delx =-poly/(pder-recsum*poly)
:06:06:       x = x+delx
:06:06:       k = k+1
:-1:06:     enddo
:06:06:     xjac(j) = x
:-1:06:   enddo
:-3:06:! ================================================
:-3:06:! compute the second half of the roots by symmetry
:-3:06:! ================================================
:05:05:   do j = 0,nh
:06:06:     xjac(n-j) =-xjac(j)
:-1:06:   enddo
:06:06:   if (modulo(n,2)==0) xjac(nh) = c0
:-3:06:! ====================================================
:-3:06:! Reverse the sign of everything so that indexing
:-3:06:! increases with position
:-3:06:! ====================================================
:05:05:   do j = 0,n
:06:06:     pts(j+1) =-xjac(j)
:-1:06:   enddo
!
:-1:01:   end function gauss_pts
:-3:01:! ================================================
:-3:01:! gauss_wts:
:-3:01:!
:-3:01:! Gauss Legendre Weights
:-3:01:! ================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function gauss_wts(np1, gpts) result(wts)
!
:04:04:   integer, intent(in   ) :: np1
:04:04:   real(longdouble_kind), intent(in   ) :: gpts(np1) ! gauss-legendre points
:04:04:   real(longdouble_kind) :: wts(np1) ! gauss-legendre weights
:-3:04:! Local variables
:04:04:   real(longdouble_kind) :: c0, c1, c2
:04:04:   real(longdouble_kind) :: alpha
:04:04:   real(longdouble_kind) :: beta
:04:04:   real(longdouble_kind) :: djac(np1)
:04:04:   integer i, n
!
:06:06:   c0 = 0.0_longdouble_kind
:06:06:   c1 = 1.0_longdouble_kind
:06:06:   c2 = 2.0_longdouble_kind
:06:06:   alpha = c0
:06:06:   beta = c0
:06:06:   n = np1-1
:06:06:   djac = jacobi_derivatives(np1,alpha,beta,np1,gpts)
:05:05:   do i = 1,np1
:06:06:     wts(i) = c2/((c1-gpts(i)**2)*djac(i)*djac(i))
:-1:06:   enddo
!
:-1:01:   end function gauss_wts
:-3:01:#endif
:-3:01:! ==============================================================
:-3:01:! test_gauss:
:-3:01:!
:-3:01:! Unit Tester for Gaussian Points, Weights
:-3:01:! ==============================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   subroutine test_gauss(npts)
!
:04:04:   integer, intent(in   ) :: npts
:04:04:   type(quadrature_t) :: gs
:04:04:   integer i
:04:04:   real(real_kind) :: gssum
!
:06:06:   gs = gauss(npts)
:06:06:   print*
:06:06:   print*,"============================================"
:06:06:   print*," testing gaussian quadrature..."
:06:06:   print*
:06:06:   print*," points weights"
:06:06:   print*,"============================================"
:05:05:   do i = 1,npts
:06:06:     print*,i,gs%points(i),gs%weights(i)
:-1:06:   enddo
:06:06:   print*,"============================================"
:06:06:   gssum = sum(gs%weights(:))
:06:06:   print*,"sum of gaussian weights = ",gssum
:06:06:   print*,"============================================"
:06:06:   deallocate(gs%points)
:06:06:   deallocate(gs%weights)
!
:-1:01:   end subroutine test_gauss
:-3:01:! ==============================================================
:-3:01:! gausslobatto:
:-3:01:!
:-3:01:! Find the Gauss-Lobatto Legendre collocation points xgl(i) and the
:-3:01:! corresponding weights.
:-3:01:!
:-3:01:! ==============================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function gausslobatto(npts) result(gll)
!
:04:04:   integer, intent(in   ) :: npts
:04:04:   type(quadrature_t) :: gll
!
:06:06:   allocate(gll%points(npts))
:06:06:   allocate(gll%weights(npts))
:06:06:   gll%points = gausslobatto_pts(npts)
:06:06:   gll%weights = gausslobatto_wts(npts,gll%points)
!
:-1:01:   end function gausslobatto
:-3:01:! ==============================================================
:-3:01:! gausslobatto_pts:
:-3:01:!
:-3:01:! Compute the Gauss-Lobatto Collocation points
:-3:01:! for Jacobi Polynomials
:-3:01:!
:-3:01:! ==============================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function gausslobatto_pts(np1) result(pts)
!-------------------------------------------------------------------------------
:03:03:   use physicalconstants, only: qq_pi
!
:04:04:   integer, intent(in   ) :: np1 ! number of velocity grid points
:04:04:   real(longdouble_kind) :: pts(np1)
:-3:04:! Local variables
:04:04:   real(longdouble_kind) :: alpha, beta
:04:04:   real(longdouble_kind) :: xjac(0:np1-1)
:04:04:   real(longdouble_kind) :: jac(0:np1)
:04:04:   real(longdouble_kind) :: jacm1(0:np1)
:04:04:   real(longdouble_kind) :: djac(0:np1)
:04:04:   integer prec ! number of mantissa bits
:04:04:   real(longdouble_kind) eps ! machine epsilon
:04:04:   real(longdouble_kind), parameter :: convthresh = 10 ! convergence threshold relative
:-3:04:! to machine epsilon
:04:04:   integer, parameter :: kstop = 30 ! max iterations for polynomial deflation
:04:04:   real(longdouble_kind) :: a, b, det
:04:04:   real(longdouble_kind) :: poly
:04:04:   real(longdouble_kind) :: pder
:04:04:   real(longdouble_kind) :: recsum, thresh
:04:04:   real(longdouble_kind) :: dth, cd, sd, cs, ss, cstmp
:04:04:   real(longdouble_kind) :: x
:04:04:   real(longdouble_kind) :: delx
:04:04:   real(longdouble_kind) :: c0, c1, c2, c10
:04:04:   integer i, j, k
:04:04:   integer n, nh
!
:06:06:   n = np1-1
:06:06:   c0 = 0.0_longdouble_kind
:06:06:   c1 = 1.0_longdouble_kind
:06:06:   c2 = 2.0_longdouble_kind
:06:06:   c10 = 10.0_longdouble_kind
:06:06:   alpha = c0
:06:06:   beta = c0
:-3:06:! =========================================================
:-3:06:! compute machine precision and set the convergence
:-3:06:! threshold thresh to 10 times that level
:-3:06:! =========================================================
:06:06:   prec = precision(c10)
:06:06:   eps = c10**(-prec)
:06:06:   thresh = convthresh*eps
:-3:06:! =====================================================
:-3:06:! initialize the end points
:-3:06:! =====================================================
:06:06:   xjac(0) = c1
:06:06:   xjac(n) =-c1
:-3:06:! ============================================================
:-3:06:! Compute first half of the roots by "polynomial deflation".
:-3:06:! ============================================================
:-3:06:! ============================================================
:-3:06:! compute the parameters in the polynomial whose
:-3:06:! roots are desired...
:-3:06:! ============================================================
:06:06:   call jacobi(n+1,c1,alpha,beta,jac(0:n+1),djac(0:n+1))
:06:06:   call jacobi(n+1,-c1,alpha,beta,jacm1(0:n+1),djac(0:n+1))
:06:06:   det = jac(n)*jacm1(n-1)-jacm1(n)*jac(n-1)
:06:06:   a =-(jac(n+1)*jacm1(n-1)-jacm1(n+1)*jac(n-1))/det
:06:06:   b =-(jac(n)*jacm1(n+1)-jacm1(n)*jac(n+1))/det
:06:06:   dth = qq_pi/(2*n+1)
:06:06:   cd = cos(c2*dth)
:06:06:   sd = sin(c2*dth)
:06:06:   cs = cos(dth)
:06:06:   ss = sin(dth)
:06:06:   nh =(n+1)/2
:05:05:   do j = 1,nh-1
:06:06:     x = cs ! first guess at root
:06:06:     k = 0
:06:06:     delx = c1
:05:05:     do while(k<kstop.and.abs(delx)>thresh)
:06:06:       call jacobi(n+1,x,alpha,beta,jac(0:n+1),djac(0:n+1))
:06:06:       poly = jac(n+1)+a*jac(n)+b*jac(n-1)
:06:06:       pder = djac(n+1)+a*djac(n)+b*djac(n-1)
:06:06:       recsum = c0
:05:05:       do i = 0,j-1
:06:06:         recsum = recsum+c1/(x-xjac(i))
:-1:06:       enddo
:06:06:       delx =-poly/(pder-recsum*poly)
:06:06:       x = x+delx
:06:06:       k = k+1
:-1:06:     enddo
:06:06:     xjac(j) = x
:-3:06:! =====================================================
:-3:06:!  compute the guesses for the roots
:-3:06:!  for the next points,i.e :
:-3:06:!
:-3:06:!  ss = sn(theta) => sin(theta+2*dth)
:-3:06:!  cs = cs(theta) => cs(theta+2*dth)
:-3:06:! =====================================================
:06:06:     cstmp = cs*cd-ss*sd
:06:06:     ss = cs*sd+ss*cd
:06:06:     cs = cstmp
:-1:06:   enddo
:-3:06:! ================================================
:-3:06:! compute the second half of the roots by symmetry
:-3:06:! ================================================
:05:05:   do j = 1,nh
:06:06:     xjac(n-j) =-xjac(j)
:-1:06:   enddo
:06:06:   if (modulo(n,2)==0) xjac(nh) = c0
:-3:06:! ====================================================
:-3:06:! Reverse the sign of everything so that indexing
:-3:06:! increases with position
:-3:06:! ====================================================
:05:05:   do j = 0,n
:06:06:     pts(j+1) =-xjac(j)
:-1:06:   enddo
!
:-1:01:   end function gausslobatto_pts
:-3:01:! ================================================
:-3:01:! Gauss Lobatto Legendre Weights
:-3:01:! ================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function gausslobatto_wts(np1, glpts) result(wts)
!
:04:04:   integer, intent(in   ) :: np1
:04:04:   real(longdouble_kind), intent(in   ) :: glpts(np1)
:04:04:   real(longdouble_kind) :: wts(np1)
:-3:04:! Local variables
:04:04:   real(longdouble_kind) :: c0, c2
:04:04:   real(longdouble_kind) :: alpha
:04:04:   real(longdouble_kind) :: beta
:04:04:   real(longdouble_kind) :: jac(np1)
:04:04:   integer i, n
!
:06:06:   c0 = 0.0_longdouble_kind
:06:06:   c2 = 2.0_longdouble_kind
:06:06:   alpha = c0
:06:06:   beta = c0
:06:06:   n = np1-1
:06:06:   jac = jacobi_polynomials(n,alpha,beta,np1,glpts)
:05:05:   do i = 1,np1
:06:06:     wts(i) = c2/(n*(n+1)*jac(i)*jac(i))
:-1:06:   enddo
!
:-1:01:   end function gausslobatto_wts
:-3:01:! ==============================================================
:-3:01:! test_gausslobatto:
:-3:01:!
:-3:01:! Unit Tester for Gaussian Lobatto Quadrature...
:-3:01:! ==============================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   subroutine test_gausslobatto(npts)
!
:04:04:   integer, intent(in   ) :: npts
:04:04:   type(quadrature_t) :: gll
:04:04:   integer i
:04:04:   real(real_kind) :: gllsum
!
:06:06:   gll = gausslobatto(npts)
:06:06:   print*
:06:06:   print*,"============================================"
:06:06:   print*," testing gauss-lobatto quadrature..."
:06:06:   print*
:06:06:   print*," points weights"
:06:06:   print*,"============================================"
:05:05:   do i = 1,npts
:06:06:     print*,i,gll%points(i),gll%weights(i)
:-1:06:   enddo
:06:06:   print*,"============================================"
:06:06:   gllsum = sum(gll%weights(:))
:06:06:   print*,"sum of gauss-lobatto weights = ",gllsum
:06:06:   print*,"============================================"
:06:06:   deallocate(gll%points)
:06:06:   deallocate(gll%weights)
!
:-1:01:   end subroutine test_gausslobatto
:-3:01:! ================================================
:-3:01:!
:-3:01:! subroutine jacobi:
:-3:01:!
:-3:01:!  Computes the Jacobi Polynomials (jac) and their
:-3:01:!  first derivatives up to and including degree n
:-3:01:!  at point x on the interval (-1,1).
:-3:01:!
:-3:01:!    See for example the recurrence relations
:-3:01:!    in equation 2.5.4 (page 70) in
:-3:01:!
:-3:01:!    "Spectral Methods in Fluid Dynamics",
:-3:01:!    by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
:-3:01:!    Springer-Verlag, 1988.
:-3:01:! ================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   subroutine jacobi(n, x, alpha, beta, jac, djac)
!
:04:04:   integer, intent(in   ) :: n
:04:04:   real(longdouble_kind), intent(in   ) :: x
:04:04:   real(longdouble_kind), intent(in   ) :: alpha
:04:04:   real(longdouble_kind), intent(in   ) :: beta
:04:04:   real(longdouble_kind) :: jac(0:n)
:04:04:   real(longdouble_kind) :: djac(0:n)
:-3:04:! Local variables
:04:04:   real(longdouble_kind) :: a1k
:04:04:   real(longdouble_kind) :: a2k
:04:04:   real(longdouble_kind) :: a3k
:04:04:   real(longdouble_kind) :: da2kdx
:04:04:   real(longdouble_kind) :: c2, c1, c0
:04:04:   integer :: k
!
:06:06:   c0 = 0.0_longdouble_kind
:06:06:   c1 = 1.0_longdouble_kind
:06:06:   c2 = 2.0_longdouble_kind
:06:06:   jac(0) = c1
:06:06:   jac(1) =(c1+alpha)*x
:06:06:   djac(0) = c0
:06:06:   djac(1) =(c1+alpha)
:05:05:   do k = 1,n-1
:06:06:     a1k = c2*(k+c1)*(k+alpha+beta+c1)*(c2*k+alpha+beta)
:06:06:     da2kdx =(c2*(k+c1)+alpha+beta)*(c2*k+alpha+beta+c1)*(c2*k+alpha+beta)
:06:06:     a2k =(c2*k+alpha+beta+c1)*(alpha*alpha-beta*beta)+x*da2kdx
:06:06:     a3k = c2*(k+alpha)*(k+beta)*(c2*k+alpha+beta+c2)
:06:06:     jac(k+1) =(a2k*jac(k)-a3k*jac(k-1))/a1k
:06:06:     djac(k+1) =(a2k*djac(k)+da2kdx*jac(k)-a3k*djac(k-1))/a1k
:-1:06:   enddo
!
:-1:01:   end subroutine jacobi
:-3:01:! ==========================================================
:-3:01:! This routine computes the Nth order Jacobi Polynomials
:-3:01:! (jac) for a vector of positions x on the interval (-1,1),
:-3:01:! of length npoints.
:-3:01:!
:-3:01:!    See for example the recurrence relations
:-3:01:!    in equation 2.5.4 (page 70) in
:-3:01:!
:-3:01:!     "Spectral Methods in Fluid Dynamics",
:-3:01:!     by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
:-3:01:!     Springer-Verlag, 1988.
:-3:01:!
:-3:01:! ===========================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function jacobi_polynomials(n, alpha, beta, npoints, x) result(jac)
!
:04:04:   integer, intent(in   ) :: n ! order of the jacobi polynomial
:04:04:   real(longdouble_kind) :: alpha
:04:04:   real(longdouble_kind) :: beta
:04:04:   integer, intent(in   ) :: npoints
:04:04:   real(longdouble_kind) :: x(npoints)
:04:04:   real(longdouble_kind) :: jac(npoints)
:-3:04:! Local variables
:04:04:   real(longdouble_kind) :: a1k
:04:04:   real(longdouble_kind) :: a2k
:04:04:   real(longdouble_kind) :: a3k
:04:04:   real(longdouble_kind) :: da2kdx
:04:04:   real(longdouble_kind) :: jacp1
:04:04:   real(longdouble_kind) :: jacm1
:04:04:   real(longdouble_kind) :: jac0
:04:04:   real(longdouble_kind) :: xtmp
:04:04:   real(longdouble_kind) :: c2, c1, c0
:04:04:   integer j, k
!
:06:06:   c0 = 0.0_longdouble_kind
:06:06:   c1 = 1.0_longdouble_kind
:06:06:   c2 = 2.0_longdouble_kind
:05:05:   do j = 1,npoints
:06:06:     xtmp = x(j)
:06:06:     jacm1 = c1
:06:06:     jac0 =(c1+alpha)*xtmp
:05:05:     do k = 1,n-1
:06:06:       a1k = c2*(k+c1)*(k+alpha+beta+c1)*(c2*k+alpha+beta)
:06:06:       da2kdx =(c2*k+alpha+beta+c2)*(c2*k+alpha+beta+c1)*(c2*k+alpha+beta)
:06:06:       a2k =(c2*k+alpha+beta+c1)*(alpha*alpha-beta*beta)+xtmp*da2kdx
:06:06:       a3k = c2*(k+alpha)*(k+beta)*(c2*k+alpha+beta+c2)
:06:06:       jacp1 =(a2k*jac0-a3k*jacm1)/a1k
:06:06:       jacm1 = jac0
:06:06:       jac0 = jacp1
:-1:06:     enddo
:06:06:     if (n==0) jac0 = jacm1
:06:06:     jac(j) = jac0
:-1:06:   enddo
!
:-1:01:   end function jacobi_polynomials
:-3:01:! ================================================
:-3:01:! This routine computes the first derivatives of Nth
:-3:01:! order Jacobi Polynomials (djac) for a vector of
:-3:01:! positions x on the interval (-1,1), of length npoints.
:-3:01:!
:-3:01:! See for example the recurrence relations
:-3:01:! in equation 2.5.4 (page 70) in
:-3:01:!
:-3:01:! "Spectral Methods in Fluid Dynamics",
:-3:01:! by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
:-3:01:! Springer-Verlag, 1988.
:-3:01:!
:-3:01:! ================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function jacobi_derivatives(n, alpha, beta, npoints, x) result(djac)
!
:04:04:   integer, intent(in   ) :: n ! order of the jacobi polynomial
:04:04:   real(longdouble_kind), intent(in   ) :: alpha
:04:04:   real(longdouble_kind), intent(in   ) :: beta
:04:04:   integer, intent(in   ) :: npoints
:04:04:   real(longdouble_kind), intent(in   ) :: x(npoints)
:04:04:   real(longdouble_kind) :: djac(npoints)
:-3:04:! Local variables
:-3:04:! Local variables
:04:04:   real(longdouble_kind) :: a1k
:04:04:   real(longdouble_kind) :: a2k
:04:04:   real(longdouble_kind) :: a3k
:04:04:   real(longdouble_kind) :: da2kdx
:04:04:   real(longdouble_kind) :: jacp1
:04:04:   real(longdouble_kind) :: jacm1
:04:04:   real(longdouble_kind) :: jac0
:04:04:   real(longdouble_kind) :: djacp1
:04:04:   real(longdouble_kind) :: djacm1
:04:04:   real(longdouble_kind) :: djac0
:04:04:   real(longdouble_kind) :: xtmp
:04:04:   real(longdouble_kind) :: c2, c1, c0
:04:04:   integer j, k
!
:06:06:   c0 = 0.0_longdouble_kind
:06:06:   c1 = 1.0_longdouble_kind
:06:06:   c2 = 2.0_longdouble_kind
:05:05:   do j = 1,npoints
:06:06:     xtmp = x(j)
:06:06:     jacm1 = c1
:06:06:     jac0 =(c1+alpha)*xtmp
:06:06:     djacm1 = c0
:06:06:     djac0 =(c1+alpha)
:05:05:     do k = 1,n-1
:06:06:       a1k = c2*(k+c1)*(k+alpha+beta+c1)*(c2*k+alpha+beta)
:06:06:       da2kdx =(c2*k+alpha+beta+c2)*(c2*k+alpha+beta+c1)*(c2*k+alpha+beta)
:06:06:       a2k =(c2*k+alpha+beta+c1)*(alpha*alpha-beta*beta)+xtmp*da2kdx
:06:06:       a3k = c2*(k+alpha)*(k+beta)*(c2*k+alpha+beta+c2)
:06:06:       jacp1 =(a2k*jac0-a3k*jacm1)/a1k
:06:06:       djacp1 =(a2k*djac0+da2kdx*jac0-a3k*djacm1)/a1k
:06:06:       jacm1 = jac0
:06:06:       jac0 = jacp1
:06:06:       djacm1 = djac0
:06:06:       djac0 = djacp1
:-1:06:     enddo
:06:06:     if (n==0) djac0 = djacm1
:06:06:     djac(j) = djac0
:-1:06:   enddo
!
:-1:01:   end function jacobi_derivatives
:-3:01:! ===================================================
:-3:01:!
:-3:01:! legendre:
:-3:01:!
:-3:01:! Compute the legendre polynomials using
:-3:01:! the recurrence relationship.
:-3:01:! return leg(m+1) = P_N(x) for m=0..N
:-3:01:! p_3 = Legendre polynomial of degree N
:-3:01:! p_2 = Legendre polynomial of degree N-1 at x
:-3:01:! p_1 = Legendre polynomial of degree N-2 at x
:-3:01:!
:-3:01:! ===================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function legendre(x, n) result(leg)
!
:04:04:   integer :: n
:04:04:   real(longdouble_kind) :: x
:04:04:   real(longdouble_kind) :: leg(n+1)
:04:04:   real(longdouble_kind) :: p_1, p_2, p_3
:04:04:   integer :: k
!
:06:06:   p_3 = 1.0_longdouble_kind
:06:06:   leg(1) = p_3
:05:05:   if (n.ne.0) then
:06:06:     p_2 = p_3
:06:06:     p_3 = x
:06:06:     leg(2) = p_3
:05:05:     do k = 2,n
:06:06:       p_1 = p_2
:06:06:       p_2 = p_3
:06:06:       p_3 =((2*k-1)*x*p_2-(k-1)*p_1)/k
:06:06:       leg(k+1) = p_3
:-1:06:     enddo
:-1:06:   endif
!
:-1:01:   end function legendre
:-3:01:! ===========================================
:-3:01:! quad_norm:
:-3:01:!
:-3:01:! compute normalization constants
:-3:01:! for k=1,N order Legendre polynomials
:-3:01:!
:-3:01:! e.g. gamma(k) in Canuto, page 58.
:-3:01:!
:-3:01:! ===========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function quad_norm(gquad, n) result(gamma)
!
:04:04:   type(quadrature_t), intent(in   ) :: gquad
:04:04:   integer, intent(in   ) :: n
:04:04:   real(longdouble_kind) :: gamma(n)
:-3:04:! Local variables
:04:04:   real(longdouble_kind) :: leg(n)
:04:04:   integer :: i, k
!
:06:06:   gamma(:) = 0.0_longdouble_kind
:05:05:   do i = 1,n
:06:06:     leg = legendre(gquad%points(i),n-1)
:05:05:     do k = 1,n
:06:06:       gamma(k) = gamma(k)+leg(k)*leg(k)*gquad%weights(i)
:-1:06:     enddo
:-1:06:   enddo
!
:-1:01:   end function quad_norm
:-3:01:! =======================
:-3:01:! TrapN:
:-3:01:! Numerical recipes
:-3:01:! =======================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   subroutine trapn(f, a, b, n, it, s)
!
:05:05:   interface
:02:02:     function f(x) result(f_x) ! function to be integrated
!-------------------------------------------------------------------------------
:03:03:     use kiapsbase, only: real_kind=>kim_real8_kind
!
:04:04:     real(real_kind), intent(in   ) :: x
:04:04:     real(real_kind) :: f_x
:-1:02:     end function f
:-1:02:   end interface
!
:04:04:   real(real_kind), intent(in   ) :: a, b
:04:04:   integer, intent(in   ) :: n
:04:04:   integer, intent(inout) :: it
:04:04:   real(real_kind), intent(inout) :: s
:04:04:   real(real_kind) :: ssum
:04:04:   real(real_kind) :: del
:04:04:   real(real_kind) :: rtnm
:04:04:   real(real_kind) :: x
:04:04:   integer :: j
!
:05:05:   if (n==1) then
:06:06:     s = 0.5d0*(b-a)*(f(a)+f(b))
:06:06:     it = 1
:05:05:   else
:06:06:     ssum = 0.0d0
:06:06:     rtnm = 1.0d0/it
:06:06:     del =(b-a)*rtnm
:06:06:     x = a+0.5*del
:05:05:     do j = 1,it
:06:06:       ssum = ssum+f(x)
:06:06:       x = x+del
:-1:06:     enddo
:06:06:     s = 0.5d0*(s+del*ssum)
:06:06:     it = 2*it
:-1:04:   endif
!
:-1:01:   end subroutine trapn
:-3:01:! ==========================================
:-3:01:! Trapezoid Rule for integrating functions
:-3:01:! from a to b with residual error eps
:-3:01:! ==========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function trapezoid(f, a, b, eps) result(integral)
!
:04:04:   integer, parameter :: nmax = 25 ! at most 2^nmax+1 points in integral
!
:05:05:   interface
:02:02:     function f(x) result(f_x) ! function to be integrated
!-------------------------------------------------------------------------------
:03:03:     use kiapsbase, only: real_kind=>kim_real8_kind
!
:04:04:     real(real_kind), intent(in   ) :: x
:04:04:     real(real_kind) :: f_x
:-1:02:     end function f
!
:-1:04:   end interface
:04:04:   real(real_kind), intent(in   ) :: a, b ! the integral bounds
:04:04:   real(real_kind), intent(in   ) :: eps ! relative error bound for integral
:04:04:   real(real_kind) :: integral ! the integral result(within eps)
:04:04:   real(real_kind) :: s ! integral approximation
:04:04:   real(real_kind) :: sold ! previous integral approx
:04:04:   integer :: n
:04:04:   integer :: it
:-3:04:! ==============================================================
:-3:04:! Calculate I here using trapezoid rule using f and a DO loop...
:-3:04:! ==============================================================
!
:06:06:   s = 1.0d30
:06:06:   sold = 0.0d0
:06:06:   n = 1
:06:06:   it = 0
:05:05:   do while(n<=nmax.and.abs(s-sold)>eps*abs(sold))
:06:06:     sold = s
:06:06:     call trapn(f,a,b,n,it,s)
:-3:06:#ifdef _QUAD_DBG
:06:06:     print*,"n = ",n," abs(s-sold) ",abs(s-sold)," threshold = ",abs(sold)*eps
:-3:06:#endif
:06:06:     n = n+1
:-1:06:   enddo
:06:06:   integral = s
!
:-1:01:   end function trapezoid
:-3:01:! ==========================================
:-3:01:! Simpsons Rule for integrating functions
:-3:01:! from a to b with residual error eps
:-3:01:! ==========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function simpsons(f, a, b, eps) result(integral)
!
:04:04:   integer, parameter :: nmax = 25 ! at most 2^nmax+1 points in integral
!
:05:05:   interface
:02:02:     function f(x) result(f_x) ! function to be integrated
!-------------------------------------------------------------------------------
:03:03:     use kiapsbase, only: real_kind=>kim_real8_kind
!
:04:04:     real(real_kind), intent(in   ) :: x
:04:04:     real(real_kind) :: f_x
:-1:02:     end function f
!
:-1:04:   end interface
:04:04:   real(real_kind), intent(in   ) :: a, b ! the integral bounds
:04:04:   real(real_kind), intent(in   ) :: eps ! relative error bound for integral
:04:04:   real(real_kind) :: integral ! the integral result(within eps)
:04:04:   real(real_kind) :: s ! integral approximation
:04:04:   real(real_kind) :: os ! previous integral approx
:04:04:   real(real_kind) :: st ! integral approximation
:04:04:   real(real_kind) :: ost ! previous integral approx
:04:04:   integer :: n
:04:04:   integer :: it
:-3:04:! ==============================================================
:-3:04:! Calculate I here using trapezoid rule using f and a DO loop...
:-3:04:! ==============================================================
!
:06:06:   ost = 0.0d0
:06:06:   s = 1.0d30
:06:06:   os = 0.0d0
:06:06:   n = 1
:06:06:   it = 0
:05:05:   do while((n<=nmax.and.abs(s-os)>eps*abs(os)).or.n<=2)
:06:06:     os = s
:06:06:     call trapn(f,a,b,n,it,st)
:06:06:     s =(4.0d0*st-ost)/3.0d0
:-3:06:#ifdef _QUAD_DBG
:06:06:     print*,"n = ",n," abs(s-os) = ",abs(s-os)," threshold = ",abs(os)*eps
:-3:06:#endif
:06:06:     ost = st
:06:06:     n = n+1
:-1:06:   enddo
:06:06:   integral = s
!
:-1:01:   end function simpsons
:-3:01:! ==========================================
:-3:01:! gaussian_int:
:-3:01:!
:-3:01:! Gaussian Quadrature Rule for integrating
:-3:01:! function f from a to b  with gs weights and
:-3:01:! points with precomputed gaussian quadrature
:-3:01:! and weights.
:-3:01:! ==========================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:02:02:   function gaussian_int(f, a, b, gs) result(integral)
!
:04:04:   integer, parameter :: nmax = 10 ! at most 2^nmax+1 points in integral
!
:05:05:   interface
:02:02:     function f(x) result(f_x) ! function to be integrated
!-------------------------------------------------------------------------------
:03:03:     use kiapsbase, only: real_kind=>kim_real8_kind
!
:04:04:     real(real_kind), intent(in   ) :: x
:04:04:     real(real_kind) :: f_x
:-1:02:     end function f
!
:-1:04:   end interface
:04:04:   real(real_kind), intent(in   ) :: a, b ! the integral bounds
:04:04:   type(quadrature_t), intent(in   ) :: gs ! gaussian points/wts
:04:04:   real(real_kind) :: integral ! the integral result(within eps)
:04:04:   integer :: i
:04:04:   real(real_kind) :: s, x
:-3:04:! ==============================================================
:-3:04:! Calculate I = S f(x)dx here using gaussian quadrature
:-3:04:! ==============================================================
!
:06:06:   s = 0.0d0
:05:05:   do i = 1,size(gs%points)
:06:06:     x = 0.50d0*((b-a)*gs%points(i)+(b+a))
:06:06:     s = s+gs%weights(i)*f(x)
:-1:06:   enddo
:06:06:   integral = s*(0.5d0*(b-a))
!
:-1:01:   end function gaussian_int
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
:-1:00:   end module quadrature
!-------------------------------------------------------------------------------
