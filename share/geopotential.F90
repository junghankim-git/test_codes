!-------------------------------------------------------------------------------
   module geopotential
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    2017-09-14  junghan kim    initial setup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds       , only: i4, r8, r16
   use constant    , only: p0_, cp_, rd_, akapa_, g_, t0_, cvpm_
   use hybrid_coord, only: read_hybrid_coord
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   integer(i4) :: nlev, shyb, ehyb
   real(r8), dimension(:), allocatable :: a_i, b_i, eta_i, p_i, ex_i
   real(r8), dimension(:), allocatable :: a_m, b_m, eta_m, p_m, ex_m
   real(r8) :: ps = 102400.0
!   real(r8) :: ps =   50000.0
!   real(r8) :: ps =   120000.0
!
   public :: initialize, finalize, get_height_p, get_height_ex, get_height_et, p_i, p_m
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine initialize(nlevs, hvfile_i, hvfile_m)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4)     , intent(in   ) :: nlevs
   character(len=*), intent(in   ) :: hvfile_i, hvfile_m
! local variables
   real(r8) :: tmp(nlevs+1)
   integer(i4) :: k
!
   nlev = nlevs
   allocate(a_i(nlev+1),b_i(nlev+1),eta_i(nlev+1),p_i(nlev+1),ex_i(nlev+1))
   allocate(a_m(nlev),b_m(nlev),eta_m(nlev),p_m(nlev),ex_m(nlev))
   call read_hybrid_coord(nlev,hvfile_i,hvfile_m,shyb,ehyb,a_i,b_i,eta_i,a_m,b_m,eta_m)
!
#if 0
   do k = 1,nlev
     p_m(k)  = p0_*eta_m(k)
     p_i(k)  = p0_*eta_i(k)
   enddo
   p_i(nlev+1)  = p0_*eta_i(nlev+1)
#endif
   do k = 1,nlev
     p_m(k)  = p0_*a_m(k)+ps*b_m(k)
     p_i(k)  = p0_*a_i(k)+ps*b_i(k)
   enddo
   p_i(nlev+1)  = p0_*a_i(nlev+1)+ps*b_i(nlev+1)
! reverse
   tmp(:) = eta_i(:)
   do k = 1,nlev+1
     eta_i(k) = tmp(nlev-k+2)
   enddo
   tmp(1:nlev) = eta_m(:)
   do k = 1,nlev
     eta_m(k) = tmp(nlev-k+1)
   enddo
   !
   tmp(:) = p_i(:)
   do k = 1,nlev+1
     p_i(k) = tmp(nlev-k+2)
   enddo
   tmp(1:nlev) = p_m(:)
   do k = 1,nlev
     p_m(k) = tmp(nlev-k+1)
   enddo
!
   do k = 1,nlev
     ex_m(k) = (p_m(k)/p0_)**akapa_
     ex_i(k) = (p_i(k)/p0_)**akapa_
   enddo
   ex_i(nlev+1) = (p_i(nlev+1)/p0_)**akapa_
!
   return
   end subroutine initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine finalize()
!-------------------------------------------------------------------------------
   implicit none
! local variables
!
   deallocate(a_i,b_i,eta_i,p_i,ex_i)
   deallocate(a_m,b_m,eta_m,p_m,ex_m)
!
   return
   end subroutine finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_height_p(h_i,h_m)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(nlev+1) :: h_i
   real(r8), dimension(nlev)   :: h_m
! local variables
   integer(i4) :: k
   real(r8)    :: dh
!
   h_i(1) = 0.0
   do k = 1,nlev
      dh       = rd_*t0_/p_m(k)*(p_i(k)-p_i(k+1))/g_
      h_m(k)   = h_i(k)+0.5*dh
      h_i(k+1) = h_i(k)+    dh
   enddo
!
   return
   end subroutine get_height_p
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_height_ex(h_i,h_m)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(nlev+1) :: h_i
   real(r8), dimension(nlev)   :: h_m
! local variables
   integer(i4) :: k
   real(r8)    :: dh1, dh2, tmp
!
   h_i(1) = 0.0
   do k = 1,nlev
      tmp      = cp_*t0_/ex_m(k)
      dh1      = tmp*(ex_i(k)-ex_m(k))/g_
      dh2      = tmp*(ex_m(k)-ex_i(k+1))/g_
      h_m(k)   = h_i(k)+dh1
      h_i(k+1) = h_m(k)+dh2
   enddo
!
   return
   end subroutine get_height_ex
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_height_et(h_i,h_m)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(nlev+1) :: h_i
   real(r8), dimension(nlev)   :: h_m
! local variables
   integer(i4) :: k
   real(r8)    :: dh
   real(r8), dimension(nlev)   :: deta
   real(r8), dimension(nlev)   :: t, tb, t_2
   real(r8), dimension(nlev)   :: mu, mu_2, mub, al, alb, al_2
   real(r8), dimension(nlev+1) :: p_2, pb
   real(r8), dimension(nlev+1) :: hb_i, h_i_2
   real(r8), dimension(nlev+1) :: p_2_i, pb_i
!
   do k = 1,nlev+1
     pb_i(k) = p0_*eta_i(k)
     !pb_i(k) = 80000.0*eta_i(k)
   enddo
!
   do k = 1,nlev
     pb(k) = p0_*eta_m(k)
     !pb(k) = 80000.0*eta_m(k)
   enddo
!
   do k = 1,nlev
     deta(k) = eta_i(k+1)-eta_i(k)
     mub(k)  = (pb_i(k+1)-pb_i(k))/deta(k)
     mu(k)   = (p_i(k+1)-p_i(k))/deta(k)
     mu_2(k) = mu(k)-mub(k)
     p_2(k)  = p_m(k)-pb(k)
   enddo
!
   do k = 1,nlev
#if 0
     tb(k)  = t0_*((p0_/pb(k))**akapa_)
     t(k)   = t0_*((p0_/p_m(k))**akapa_)
     t_2(k) = t(k)-tb(k)
#endif
     t(k)   = t0_*((p0_/p_m(k))**akapa_)
     tb(k)  = t0_
     t_2(k) = t(k)-tb(k)
     alb(k) = t(k)/p0_*((pb(k)/p0_)**cvpm_)
     al(k)  = t(k)/p0_*((p_m(k)/p0_)**cvpm_)
     al_2(k) = al(k)-alb(k)
   enddo
!
   hb_i(1) = 0.0
   do k = 1,nlev
      dh        = rd_*alb(k)*mub(k)*deta(k)/g_
      hb_i(k+1) = hb_i(k)-dh
   enddo
!
   h_i_2(1) = 0.0
   do k = 1,nlev
      dh         = rd_*(mu(k)*al_2(k)+mu_2(k)*alb(k))*deta(k)/g_
      h_i_2(k+1) = h_i_2(k)-dh
   enddo
!
   h_i(1) = 0.0
   do k = 1,nlev
      h_i(k+1) = hb_i(k+1)+h_i_2(k+1)
   enddo
   do k = 1,nlev
     h_m(k) = 0.5*(h_i(k)+h_i(k+1))
   enddo
!
   return
   end subroutine get_height_et
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_height_et_org(h_i,h_m)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), dimension(nlev+1) :: h_i
   real(r8), dimension(nlev)   :: h_m
! local variables
   integer(i4) :: k
   real(r8)    :: dh, theta
!#if 0
   real(r8), dimension(nlev) :: mu, deta, p_2, pb, al_2, alb
!
   do k = 1,nlev
     deta(k) = eta_i(k+1)-eta_i(k)
     mu(k)   = (p_i(k+1)-p_i(k))/deta(k)
   enddo
!#endif
!
   h_i(1) = 0.0
   do k = 1,nlev
      theta    = t0_*((p0_/p_m(k))**akapa_)
      dh       = rd_*theta/p0_*((p_m(k)/p0_)**cvpm_)*mu(k)*deta(k)/g_
      h_m(k)   = h_i(k)-0.5*dh
      h_i(k+1) = h_i(k)-    dh
      !dh       = rd_*t0_/p_m(k)*mu(k)*deta(k)/g_ ! eta
   enddo
!
   return
   end subroutine get_height_et_org
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module geopotential
!-------------------------------------------------------------------------------
