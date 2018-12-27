!-------------------------------------------------------------------------------
!
!  module hybrid_coord
!
!> @brief
!>  - module for generating of hybrid coefficients.
!>
!> @date 15feb2014
!>  - junghan kim: first working version
!>                    * eta = f(x), where x is the normalized level
!>                      (0 <= eta, x <= 1)
!>                    * In hybrid region, the bottom boundary of the pure
!>                      pressure region is linearly connected with the top
!>                      boundary of the sigma region.
!> @date 19feb2014
!>  - junghan kim: some new features added.
!>                    * Particular eta profile (e.g., ECMWF-137) can be used
!>                      through cubic spline interpolation.
!>                    * For given top and surface eta values,
!>                      smooth eta profiles can be obtained.
!>                    * A given eta profile [eta = f(x)] can be shifted
!>                      with respect to x.
!>                    * Atmospheric states for generated eta are displayed.
!
!-------------------------------------------------------------------------------
#define USE_CS
   module hybrid_coord
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    2014-02-19  junghan kim    initial setup
!    2017-02-13  junghan kim    code clean-up
!    2017-06-20  junghan kim    code refactoring, remove shift
!    2017-06-26  junghan kim    add read_file, code refactoring
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds,               only: i4, r8, l4
   use constant,            only: rgas=>rd_, g=>g_, p0=>p0_
   use std_atmosphere,      only: std_atm_initialize, get_pressure, get_height, std_ps
   use read_hybrid_file,    only: read_coefficient
   use cubic_spline_interp, only: cs_interp_initialize, cs_interp_finalize, cs_interp_get
   use spline_interp,       only: spline_t, spline_initialize, spline_finalize,&
                          spline_set_x, spline_set, spline_get, spline_integral
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   !real(r8), public :: def_topeta  = 0.000020_r8  ! default top eta (half lev(   2), 0.02 hPa)
   real(r8), public :: def_topeta  = 0.00060_r8   ! default top eta (half lev(   2), 0.02 hPa)
   real(r8), public :: def_sfceta  = 0.99763106564372372_r8   ! default sfc eta (half lev(nlev),   20 m  )
   real(r8), public :: def_pt      = 0.0_r8       ! default pt (p = ax(p0-pt) + bx(ps-pt))
   real(r8), public :: def_tmean   = 222.56_r8    ! default mean temperature
   integer(i4)      :: order_si    = 2            ! continuous on 2rd derivation
!
! interface
   interface get_eta
     module procedure get_eta_ecmwf
     module procedure get_eta_polynomial
     module procedure get_eta_equidistance
     module procedure get_eta_file
     module procedure get_eta_interpolation
   end interface get_eta
   interface get_ab
     module procedure get_ab_ecmwf
     module procedure get_ab_line
   end interface get_ab
!
! hybrid coordinate
   type hyb_coord_t
     integer(i4) :: nlevs, shyb, ehyb ! (shyb:start of hybrid, ehyb: end of hybrid
     integer(i4) :: method
     real(r8)    :: pt, tmean
     logical(l4) :: use_std_atm
     logical(l4) :: adjust_top
     logical(l4) :: adjust_sfc
     real(r8)    :: topeta
     real(r8)    :: sfceta
     logical(l4) :: initialized = .false.
     !
     real(r8), dimension(:), allocatable :: a_i, b_i, eta_i, a_m, b_m, eta_m
     real(r8), dimension(:), allocatable :: p_i, h_i, p_m, h_m
   end type hyb_coord_t
!
   public :: hyb_coord_t, initialize, get_eta, get_ab, driver, finalize, read_hybrid_coord
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine initialize(hvc, nlevs, adjust_top, adjust_sfc, top_hpa, sfc_m,     &
                                                        use_std_atm, pt, tmean)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hyb_coord_t),     intent(inout) :: hvc
   integer(i4),           intent(in   ) :: nlevs
   logical(l4), optional, intent(in   ) :: adjust_top, adjust_sfc
   real(r8)   , optional, intent(in   ) :: top_hpa, sfc_m
   logical(l4), optional, intent(in   ) :: use_std_atm
   real(r8)   , optional, intent(in   ) :: pt, tmean
! local variables
   integer(i4) :: k
   real(r8)    :: prs_2int
!
   hvc%nlevs = nlevs    ! # of model levels
   hvc%ehyb  = nlevs-1  ! start # of sigma level at interface
!
   allocate(hvc%a_i(nlevs+1))
   allocate(hvc%b_i(nlevs+1))
   allocate(hvc%eta_i(nlevs+1))
   allocate(hvc%a_m(nlevs))
   allocate(hvc%b_m(nlevs))
   allocate(hvc%eta_m(nlevs))
   allocate(hvc%p_i(nlevs+1))
   allocate(hvc%h_i(nlevs+1))
   allocate(hvc%p_m(nlevs))
   allocate(hvc%h_m(nlevs))
!
   if (present(adjust_top)) then
     hvc%adjust_top  = adjust_top
   else
     hvc%adjust_top  = .false.
   endif
   if (present(adjust_sfc)) then
     hvc%adjust_sfc  = adjust_sfc
   else
     hvc%adjust_sfc  = .false.
   endif
!
   if (present(use_std_atm)) then
     hvc%use_std_atm = use_std_atm
   else
     hvc%use_std_atm = .true.
   endif
   if (hvc%use_std_atm) then
     call std_atm_initialize(2)
   endif
!
   if (hvc%adjust_top) then
     if (present(top_hpa)) then
       hvc%topeta    = 2.0_r8*top_hpa/1000.0_r8  ! hPa -> eta
     else
       hvc%topeta    = def_topeta
     endif
   endif
   if (hvc%adjust_sfc) then
     if (present(sfc_m)) then
       if (hvc%use_std_atm) then
         prs_2int    = get_pressure(2.0_r8*sfc_m)
         hvc%sfceta  = prs_2int/std_ps  ! Pa -> eta !sfc_m
       else
         prs_2int    = get_pressure_at_point(tmean,2.0_r8*sfc_m)
         hvc%sfceta  = prs_2int/p0      ! Pa -> eta !sfc_m
       endif
     else
       hvc%sfceta    = def_sfceta
     endif
   endif
!
   if (present(pt)) then
     hvc%pt          = pt       ! top pressure at interface
   else
     hvc%pt          = def_pt
   endif
   if (present(tmean)) then
     hvc%tmean       = tmean    ! mean temperature in the atmosphere
   else
     hvc%tmean       = def_tmean
   endif
!
   hvc%initialized = .true.
!
   return
   end subroutine initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_eta_ecmwf(hvc,nlevs_ecmwf)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hyb_coord_t), intent(inout) :: hvc
   integer(i4),       intent(in   ) :: nlevs_ecmwf
! local variables
   integer(i4) :: k, nlevs, shyb_ecmwf
   character(len=512)  :: filename
   real(r8), dimension(:), allocatable :: xlev, r_xlev, r_eta, r_a, r_b
!
   if (.not.hvc%initialized) then
     print *, 'the module need basic setting values...'
     stop
   endif
!
   hvc%method = 1
!
   allocate(r_xlev(nlevs_ecmwf+1))
   allocate(r_eta(nlevs_ecmwf+1))
   allocate(r_a(nlevs_ecmwf+1))
   allocate(r_b(nlevs_ecmwf+1))
!
! reference 
   nlevs         = hvc%nlevs
   write(filename,'(a,i4.4,a)') './coefficient/ecmwf/ecmwf_half_',nlevs_ecmwf,'nlevs.ascii'
   call read_coefficient(trim(filename),nlevs_ecmwf,r_a,r_b,r_eta,.true.)
   !
   call gen_nomalize_x(nlevs_ecmwf,r_xlev)
   do k = 1,nlevs_ecmwf+1
     if (r_b(k).ne.0.0_r8) then
       shyb_ecmwf = k
       exit
     endif
   enddo
   hvc%shyb = dble(shyb_ecmwf-1)/dble(nlevs_ecmwf+1)*dble(nlevs+1)+1
!
   allocate(xlev(nlevs+1))
   call gen_nomalize_x(nlevs,xlev)
   call apply_cubic_spline(nlevs_ecmwf,r_xlev,r_eta,nlevs,xlev,hvc%eta_i)
   deallocate(xlev)
!
   deallocate(r_xlev)
   deallocate(r_eta)
   deallocate(r_a)
   deallocate(r_b)
!
   call get_eta_post(nlevs,hvc%adjust_top,hvc%adjust_sfc,hvc%topeta,hvc%sfceta,hvc%eta_i,hvc%eta_m)
!
   return
   end subroutine get_eta_ecmwf
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_eta_polynomial(hvc,shyb,coefpoly)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hyb_coord_t), intent(inout) :: hvc
   integer(i4),       intent(in   ) :: shyb
   real(r8),          intent(in   ) :: coefpoly(15)
! local variables
   integer(i4) :: k, nlevs
   real(r8)    :: sk
!
   if (.not.hvc%initialized) then
     print *, 'the module need basic setting values...'
     stop
   endif
!
   hvc%method = 2
   hvc%shyb   = shyb
   nlevs      = hvc%nlevs
!
   hvc%eta_i(1)           = 0.0_r8
   hvc%eta_i(nlevs+1) = 1.0_r8
   do k = 2,nlevs
     sk           = (dble(k)-1.0_r8)/dble(nlevs)
     hvc%eta_i(k) = polynomial(15,coefpoly,sk)
   enddo
!
   call get_eta_post(nlevs,hvc%adjust_top,hvc%adjust_sfc,hvc%topeta,hvc%sfceta,hvc%eta_i,hvc%eta_m)
!
   return
   end subroutine get_eta_polynomial
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_eta_equidistance(hvc,shyb,pt)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hyb_coord_t), intent(inout) :: hvc
   integer(i4),       intent(in   ) :: shyb
   real(r8),          intent(in   ) :: pt
! local variables
   integer(i4) :: k, nlevs
   real(r8)    :: sk, zt
!
   if (.not.hvc%initialized) then
     print *, 'the module need basic setting values...'
     stop
   endif
!
   if (hvc%pt.ne.pt) then
     print *, 'method(equidistance): check pt'
     stop
   endif
!
   hvc%method = 3
   hvc%shyb   = shyb
   nlevs      = hvc%nlevs
!
   if (pt.eq.0.0_r8) then
     print *, 'method 3: pt must not be 0'
     stop
   endif
   hvc%eta_i(1)       = 0.0_r8
   hvc%eta_i(nlevs+1) = 1.0_r8
   zt = get_height_at_point(hvc%tmean,pt)  ! .true.(p->z)
   do k = 2,nlevs
     sk           = (dble(k)-1.0_r8)/dble(nlevs)
     hvc%eta_i(k) = (p0*dexp(-g/(rgas*hvc%tmean)*(zt*(1.0d0-sk)))-pt)/(p0-pt)
   enddo
!
   call get_eta_post(nlevs,hvc%adjust_top,hvc%adjust_sfc,hvc%topeta,hvc%sfceta,hvc%eta_i,hvc%eta_m)
!
   return
   end subroutine get_eta_equidistance
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_eta_file(hvc, hvfile_i, hvfile_m)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hyb_coord_t), intent(inout) :: hvc
   character(len=*),  intent(in   ) :: hvfile_i
   character(len=*),  intent(in   ) :: hvfile_m
! local variables
!
   if (.not.hvc%initialized) then
     print *, 'the module need basic setting values...'
     stop
   endif
!
   hvc%method = 4
!
   call read_hybrid_coord(hvc%nlevs,trim(hvfile_i),trim(hvfile_m),             &
          hvc%shyb,hvc%ehyb,hvc%a_i,hvc%b_i,hvc%eta_i,hvc%a_m,hvc%b_m,hvc%eta_m)
   return
   end subroutine get_eta_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_eta_interpolation(hvc, ratio)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hyb_coord_t),      intent(inout) :: hvc
   real(r8), dimension(4), intent(in   ) :: ratio ! 0.01, 1, 100, 1000 (hPa)
! local variables
   logical(l4) :: has_thermo
   integer(i4) :: nlevs, r_nlevs, k, total
   real(r8), dimension(:), allocatable :: r_lev, r_eta, xlev
   type(spline_t) :: si
!
   if (.not.hvc%initialized) then
     print *, 'the module need basic setting values...'
     stop
   endif
!
   nlevs = hvc%nlevs
!
   hvc%method = 5
!
   if (sum(ratio).gt.100.0_r8) then
     print *, 'some error...', sum(ratio)
     stop
   endif
!
!      0hPa      0.01hPa                       1hPa          100hPa                    1000hPa
!         <- ? ->       <-         ?        ->       <-  ?  ->     <-        ?        ->
!   (0.0,0.0), (?,topeta), (?,0.00001), ..., (?,0.001), ..., (?,0.1), ...,( n   ,eta), (n+1,1.0)
!   (0.0,0.0), (?,topeta), (?,0.00001), ..., (?,0.001), ..., (?,0.1), ...,(0.xxx,eta), (1.0,1.0)
!
   if (hvc%topeta.ge.0.00001_r8) then
     has_thermo = .true.
     r_nlevs    = 6
   else
     has_thermo = .false.
     r_nlevs    = 7
   endif
   allocate(r_lev(r_nlevs))
   allocate(r_eta(r_nlevs))
!
   r_lev(1)         = 0.000_r8
   r_eta(1)         = 0.000_r8

   r_lev(2)         = (dble(2)-1.0_r8)/dble(nlevs)
   r_eta(2)         = hvc%topeta

   r_lev(r_nlevs-3) = ratio(2)*0.01_r8
   r_eta(r_nlevs-3) = 0.001_r8

   r_lev(r_nlevs-2) = ratio(3)*0.01_r8
   r_eta(r_nlevs-2) = 0.100_r8

   r_lev(r_nlevs-1) = (dble(nlevs)-1.0_r8)/dble(nlevs)
   r_eta(r_nlevs-1) = hvc%sfceta

   r_lev(r_nlevs  ) = 1.000_r8
   r_eta(r_nlevs  ) = 1.000_r8
   if (hvc%topeta.lt.0.00001_r8) then
     r_lev(3)       = ratio(1)*0.01_r8
     r_eta(3)       = 0.00001_r8
   endif
!
#ifdef USE_CS
   call cs_interp_initialize(r_nlevs,r_lev,r_eta,0.0_r8,0.0_r8)
#else
   call spline_initialize(si,r_nlevs,order_si)
   call spline_set_x(si,r_lev)
   call spline_set(si,r_eta,0.0_r8,0.0_r8,0.0_r8)
#endif
!
   hvc%shyb = 41
!
   allocate(xlev(nlevs+1))
   call gen_nomalize_x(nlevs,xlev)
!
#ifdef USE_CS
   do k = 1,nlevs+1
     hvc%eta_i(k) = cs_interp_get(xlev(k))
   enddo
#else
   call spline_get(si,nlevs+1,xlev,hvc%eta_i,0)
#endif
!#if 0
 print '(a,i3)',' nlevs = ',r_nlevs
 print '(a,6(f7.5,x))',' lev = ',r_lev
 print '(a,6(f7.5,x))',' eta = ',r_eta
! print *, xlev(:)
! print *, hvc%eta_i(:)
!#endif
!
#ifdef USE_CS
   call cs_interp_finalize()
#else
   call spline_finalize(si)
#endif
!
   deallocate(r_lev,r_eta)
   deallocate(xlev)
!
   call get_eta_post(nlevs,hvc%adjust_top,hvc%adjust_sfc,hvc%topeta,hvc%sfceta,hvc%eta_i,hvc%eta_m)
!
   return
   end subroutine get_eta_interpolation
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_eta_post(nlevs, adjust_top, adjust_sfc, topeta, sfceta, eta_i, eta_m)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: nlevs
   logical(l4),                  intent(in   ) :: adjust_top, adjust_sfc
   real(r8),                     intent(in   ) :: topeta, sfceta
   real(r8), dimension(nlevs+1), intent(  out) :: eta_i
   real(r8), dimension(nlevs)  , intent(  out) :: eta_m
!
   integer(i4) :: k
!
   if (adjust_top) call adjust_top_bndry(nlevs,eta_i,topeta)
   if (adjust_sfc) call adjust_sfc_bndry(nlevs,eta_i,sfceta)
!
   call check_eta(nlevs,eta_i)
!
   do k = 1,nlevs
     eta_m(k) = 0.5_r8*(eta_i(k)+eta_i(k+1))
   enddo
!
   return
   end subroutine get_eta_post
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_ab_ecmwf(hvc, nlevs_ecmwf)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hyb_coord_t), intent(inout) :: hvc
   integer(i4),       intent(in   ) :: nlevs_ecmwf
! local variables
   integer(i4) :: nlevs, shyb, ehyb, k
!
   if (hvc%method.eq.4) then
     print *, 'could not create ab in method = 4'
     stop
   endif
!
   nlevs = hvc%nlevs
   shyb  = hvc%shyb
   ehyb  = hvc%ehyb
!
! calculates coefficients a and b for interface (half) levels
!
   ! pressure layer (A)
   do k = 1,shyb-1
     hvc%a_i(k) = hvc%eta_i(k)
     hvc%b_i(k) = 0.0_r8
   enddo
!
   ! hybrid layer (A, B)
   call get_hybrid_ab_ecmwf(nlevs,shyb,ehyb,hvc%eta_i,nlevs_ecmwf,hvc%a_i,hvc%b_i)
!
   ! sigma layer (A)
   do k = ehyb+1,nlevs+1
     hvc%a_i(k) = 0.0_r8
     hvc%b_i(k) = hvc%eta_i(k)
   enddo
!
! calculates coefficients a and b for middle (full) levels
   do k = 1,nlevs
     hvc%a_m(k) = 0.5_r8*(hvc%a_i(k+1)+hvc%a_i(k))
     hvc%b_m(k) = 0.5_r8*(hvc%b_i(k+1)+hvc%b_i(k))
   enddo
!
   return
   end subroutine get_ab_ecmwf
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_ab_line(hvc)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hyb_coord_t), intent(inout) :: hvc
! local variables
   integer(i4) :: nlevs, shyb, ehyb, k
!
   if (hvc%method.eq.4) then
     print *, 'could not create ab in method = 4'
     stop
   endif
!
   nlevs = hvc%nlevs
   shyb  = hvc%shyb
   ehyb  = hvc%ehyb
!
! calculates coefficients a and b for interface (half) levels
!
   ! pressure layer (A)
   do k = 1,shyb-1
     hvc%a_i(k) = hvc%eta_i(k)
     hvc%b_i(k) = 0.0_r8
   enddo
!
   ! hybrid layer (A, B)
   call get_hybrid_ab_line(nlevs,shyb,ehyb,hvc%eta_i,hvc%a_i,hvc%b_i)
!
   ! sigma layer (A)
   do k = ehyb+1,nlevs+1
     hvc%a_i(k) = 0.0_r8
     hvc%b_i(k) = hvc%eta_i(k)
   enddo
!
! calculates coefficients a and b for middle (full) levels
   do k = 1,nlevs
     hvc%a_m(k) = 0.5_r8*(hvc%a_i(k+1)+hvc%a_i(k))
     hvc%b_m(k) = 0.5_r8*(hvc%b_i(k+1)+hvc%b_i(k))
   enddo
!
   return
   end subroutine get_ab_line
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_options(hvc)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hyb_coord_t), intent(in   ) :: hvc
! local variables
!
   if (.not.hvc%initialized) then
     print *, 'the module not initialized...'
     stop
   endif
!
   print *, ' - nlevs    = ', hvc%nlevs
   print *, ' - shyb     = ', hvc%shyb
   print *, ' - ehyb     = ', hvc%ehyb
   print *, ' - pt       = ', hvc%pt
   print *, ' - top(eta) = ', hvc%eta_i(2)
   print *, ' - sfc(eta) = ', hvc%eta_i(hvc%nlevs)
!
   return
   end subroutine print_options
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine finalize(hvc)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hyb_coord_t), intent(inout) :: hvc
!
   deallocate(hvc%a_i)
   deallocate(hvc%b_i)
   deallocate(hvc%eta_i)
   deallocate(hvc%a_m)
   deallocate(hvc%b_m)
   deallocate(hvc%eta_m)
   deallocate(hvc%p_i)
   deallocate(hvc%h_i)
   deallocate(hvc%p_m)
   deallocate(hvc%h_m)
!
   return
   end subroutine finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine driver(hvc)
!-------------------------------------------------------------------------------
   implicit none
!
   type(hyb_coord_t), intent(inout) :: hvc
!
   if (.not.hvc%initialized) then
     print *, 'the module not initialized...'
     stop
   endif
!
   call print_options(hvc)
!
   call get_pressure_height(hvc%nlevs,hvc%pt,hvc%tmean,hvc%use_std_atm,        &
                            hvc%eta_i,hvc%a_i,hvc%b_i,hvc%eta_m,               &
                            hvc%p_i,hvc%p_m,hvc%h_i,hvc%h_m)
!
   if (hvc%method.ne.4) then
   call write_hybrid_coord(hvc%nlevs,hvc%shyb,hvc%ehyb,hvc%p_m(1),hvc%h_m(1),  &
                           hvc%tmean,hvc%a_i,hvc%a_m,hvc%b_i,hvc%b_m)
   endif
!
   call print_hybrid_coord(hvc%nlevs,hvc%p_i,hvc%p_m,hvc%h_i,hvc%h_m,          &
                           hvc%eta_i,hvc%eta_m,hvc%a_i,hvc%a_m,hvc%b_i,hvc%b_m)
!
   call analysis_hybrid_coord(hvc%nlevs,hvc%p_i,hvc%p_m,hvc%h_i,hvc%h_m)
!
   return
   end subroutine driver
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine gen_nomalize_x(n, lev)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),              intent(in   ) :: n
   real(r8), dimension(n+1), intent(  out) :: lev
!
   integer(i4) :: k
!
   lev(1)   = 0.0_r8
   lev(n+1) = 1.0_r8
   do k = 2,n
     lev(k) =(dble(k)-1.0_r8)/dble(n)
   enddo
!
   return
   end subroutine gen_nomalize_x
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_height_at_point(temperature, pressure) result(height)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: temperature, pressure
   real(r8)                :: height
!
   height = rgas*temperature/g*dlog(p0/pressure)
!
   end function get_height_at_point
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_pressure_at_point(temperature, height) result(pressure)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: temperature, height
   real(r8)                :: pressure
!
   pressure = p0*dexp(-g*height/(rgas*temperature))
!
   end function get_pressure_at_point
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function polynomial(norder, coefpoly, x)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                 intent(in   ) :: norder
   real(r8), dimension(norder), intent(in   ) :: coefpoly
   real(r8),                    intent(in   ) :: x
   real(r8)                                   :: polynomial
!
   integer(i4) :: i
   real(r8)    :: w
!
   w = 0.0_r8
   do i = 1,norder+1
     w = w+coefpoly(i)*(x**dble(i-1))
   enddo
!
   polynomial = w
!
   end function polynomial
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine adjust_top_bndry(nlevs, eta, etat)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: nlevs
   real(r8), dimension(nlevs+1), intent(inout) :: eta
   real(r8),                     intent(in   ) :: etat
!
   integer(i4) :: k
   real(r8)    :: pi, etact, alpha
!
   pi = dacos(-1.0_r8)
!
   alpha = 20.0_r8
   etact = eta(2)
!
   do k = 2,nlevs
!     eta(k) = eta(k) - (etact-etat)*  &
!                     dexp(-0.01_r8*dble(k-2)/dble(nlevs+1))
!     eta(k) = eta(k) - (etact-etat)*  &
!                     dcos(dble(k-2)/dble(nlevs+1)*0.5_r8*pi)**0.1_r8
     eta(k) = eta(k)*(1.0_r8-(1.0_r8-etat/etact)*dexp(-alpha*dble(k-2)/dble(nlevs+1)))
   enddo
!
   return
!
   return
   end subroutine adjust_top_bndry
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine adjust_sfc_bndry(nlevs, eta, etas)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: nlevs
   real(r8), dimension(nlevs+1), intent(inout) :: eta
   real(r8),                     intent(in   ) :: etas
!
   integer(i4) :: k
   real(r8)    :: pi, etact, beta
!
   pi = dacos(-1.0_r8)
!
   beta = 20.0_r8
   etact = eta(nlevs)
!
   do k = 2,nlevs
!     eta(k) = eta(k) - (etact-etas)*  &
!                     dexp(-0.01_r8*dble(k-2)/dble(nlevs+1))
!     eta(k) = eta(k) - (etact-etas)*  &
!                     dcos(dble(k-2)/dble(nlevs+1)*0.5_r8*pi)**0.1_r8
     eta(k) = eta(k)*(1.0_r8-(1.0_r8-etas/etact)*dexp(beta*dble(k-nlevs)/dble(nlevs+1)))
   enddo
!
   return
   end subroutine adjust_sfc_bndry
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine apply_cubic_spline(r_nlevs, r_xlev, r_eta, nlevs, xlev, eta)
!-------------------------------------------------------------------------------
   implicit none
!
!
   integer(i4),                    intent(in   ) :: r_nlevs
   real(r8), dimension(r_nlevs+1), intent(in   ) :: r_xlev
   real(r8), dimension(r_nlevs+1), intent(in   ) :: r_eta
   integer(i4),                    intent(in   ) :: nlevs
   real(r8), dimension(nlevs+1),   intent(in   ) :: xlev
   real(r8), dimension(nlevs+1),   intent(inout) :: eta
!
   integer(i4)                  :: k
   type(spline_t) :: si
!
#ifdef USE_CS
   call cs_interp_initialize(r_nlevs+1,r_xlev,r_eta,0.0_r8,0.0_r8)
!
   eta(1)       = 0.0_r8
   eta(nlevs+1) = 1.0_r8
   do k = 2,nlevs
     eta(k) = cs_interp_get(xlev(k))
   enddo
!
   call cs_interp_finalize()
#else
   call spline_initialize(si,r_nlevs+1,order_si)
   call spline_set_x(si,r_xlev)
   call spline_set(si,r_eta,0.0_r8,0.0_r8,0.0_r8)
   call spline_get(si,nlevs+1,xlev,eta,0)
   call spline_finalize(si)
#endif
!
   return
   end subroutine apply_cubic_spline
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine shift_xaxis(nlevs, xaxis_in, shift, xaxis_out)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: nlevs
   real(r8), dimension(nlevs+1), intent(in   ) :: xaxis_in
   real(r8),                     intent(in   ) :: shift
   real(r8), dimension(nlevs+1), intent(  out) :: xaxis_out
!
   real(r8)    :: pi
   integer(i4) :: k
!
   pi = dacos(-1.0_r8)
!
   xaxis_out(1)       = xaxis_in(1)
   xaxis_out(nlevs+1) = xaxis_in(nlevs+1)
!
   do k = 2,nlevs
     xaxis_out(k) = xaxis_in(k)-shift*dsin(dble(k)/dble(nlevs+1)*pi)**3.0_r8
   enddo
!
   return
   end subroutine shift_xaxis
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine check_eta(nlevs, eta_i)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: nlevs
   real(r8), dimension(nlevs+1), intent(in   ) :: eta_i
!
   integer(i4) :: k
!
   ! check eta < 0.0
   do k = 2,nlevs
     if (eta_i(k)<0.0_r8) then
       write(*,'(a,i3,a,f13.9,a)') ' - error: eta_i(',k,') = ',eta_i(k),' < 0'
     endif
   enddo
!
   ! check inverse
   do k = 2,nlevs+1
     if (eta_i(k)<eta_i(k-1)) then
       write(6,'(a,i3,a,i3,a,f13.9,a,f13.9)') ' - error: eta_i(',k-1,') > eta_i(',k,'): ',eta_i(k-1),' < ',eta_i(k)
!       call exit(-1)
     endif
   enddo
!
   return
   end subroutine check_eta
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_hybrid_ab_ecmwf(nlevs,shyb,ehyb,eta_i,nlevs_ecmwf,a_i,b_i)
!-------------------------------------------------------------------------------
   implicit none
!
!
   integer(i4),                  intent(in   ) :: nlevs, shyb, ehyb
   real(r8), dimension(nlevs+1), intent(in   ) :: eta_i
   integer(l4),                  intent(in   ) :: nlevs_ecmwf
   real(r8), dimension(nlevs+1), intent(  out) :: a_i, b_i
!
   integer(i4) :: k
   real(r8), dimension(:), allocatable :: r_eta, r_a, r_b
   character(len=512)                  :: filename
!
   allocate(r_eta(nlevs_ecmwf+1))
   allocate(r_a(nlevs_ecmwf+1))
   allocate(r_b(nlevs_ecmwf+1))
   !
   write(filename,'(a,i4.4,a)') './coefficient/ecmwf/ecmwf_half_',nlevs_ecmwf,'nlevs.ascii'
   call read_coefficient(trim(filename),nlevs_ecmwf,r_a,r_b,r_eta,.true.)
   !
   call apply_cubic_spline_for_b(nlevs_ecmwf,r_eta,r_b,nlevs,shyb,ehyb,eta_i,b_i)
   do k = shyb,ehyb
     a_i(k) = eta_i(k)-b_i(k)
   enddo
   !
   deallocate(r_eta)
   deallocate(r_a)
   deallocate(r_b)
!
   return
   end subroutine get_hybrid_ab_ecmwf
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_hybrid_ab_line(nlevs,shyb,ehyb,eta_i,a_i,b_i)
!-------------------------------------------------------------------------------
   implicit none
!
!
   integer(i4),                  intent(in   ) :: nlevs, shyb, ehyb 
   real(r8), dimension(nlevs+1), intent(in   ) :: eta_i
   real(r8), dimension(nlevs+1), intent(  out) :: a_i, b_i
!
   integer(i4) :: k
!
   do k = shyb,ehyb        ! hybrid layer
     !a_i(k) = eta_i(nprs)*((eta_i(nsig)-eta_i(k))/(eta_i(nsig)-eta_i(nprs)))
     a_i(k) = eta_i(shyb-1)*((eta_i(ehyb+1)-eta_i(k))/(eta_i(ehyb+1)-eta_i(shyb-1)))
     b_i(k) = eta_i(k)-a_i(k)
   enddo
!
   return
   end subroutine get_hybrid_ab_line
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine apply_cubic_spline_for_b(r_nlevs, r_eta, r_b_i, nlevs, shyb, ehyb, eta_i, coef_b)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: r_nlevs
   real(r8), dimension(:),       intent(in   ) :: r_eta, r_b_i
   integer(i4),                  intent(in   ) :: nlevs
   integer(i4),                  intent(in   ) :: shyb, ehyb
   real(r8), dimension(nlevs+1), intent(in   ) :: eta_i
   real(r8), dimension(nlevs+1), intent(inout) :: coef_b
   integer(i4) :: k
   type(spline_t) :: si
!
#ifdef USE_CS
   call cs_interp_initialize(r_nlevs+1,r_eta,r_b_i,0.0_r8,1.0_r8)
!
   coef_b(1) = 0.0_r8
   coef_b(nlevs+1) = 1.0_r8
   do k = shyb,ehyb
     coef_b(k) = cs_interp_get(eta_i(k))
     if (coef_b(k)<0.0_r8) coef_b(k) = 0.0_r8
   enddo
!
   call cs_interp_finalize()
#else
   call spline_initialize(si,r_nlevs+1,order_si)
   call spline_set_x(si,r_eta)
   call spline_set(si,r_b_i,0.0_r8,0.0_r8,0.0_r8)
   coef_b(1) = 0.0_r8
   coef_b(nlevs+1) = 1.0_r8
   call spline_get(si,ehyb-shyb+1,eta_i(shyb:ehyb),coef_b(shyb:ehyb))
   call spline_finalize(si)
#endif
!
   return
   end subroutine apply_cubic_spline_for_b
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_pressure_height(nlevs, pt, tmean, use_std_atm,          &
                                  eta_i, a_i, b_i, eta_m, p_i, p_m, h_i, h_m)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: nlevs
   real(r8),                     intent(in   ) :: pt, tmean
   real(r8), dimension(nlevs+1), intent(in   ) :: eta_i
   real(r8), dimension(nlevs+1), intent(in   ) :: a_i, b_i
   real(r8), dimension(nlevs),   intent(in   ) :: eta_m
   logical(l4),                  intent(in   ) :: use_std_atm
   real(r8), dimension(nlevs+1), intent(  out) :: p_i, h_i
   real(r8), dimension(nlevs),   intent(  out) :: p_m, h_m
!
   real(r8) :: hs, zt
   integer(i4) :: k
!
   ! calcuate values of pressure on the interface levels
   do k = 1,nlevs+1
     p_i(k) =(p0-pt)*a_i(k)+(std_ps-pt)*b_i(k)+pt ! std_ps = 1015.25 hpa
   enddo
!
   ! calculate values of pressure on the middle levels
   do k = 1,nlevs
     p_m(k) =(p_i(k)+p_i(k+1))*0.5_r8
   enddo
!
   ! calculate values of height on the interface levels
   if (use_std_atm) then
     do k = 1,nlevs
       h_i(k) = get_height(p_i(k))
       h_m(k) = get_height(p_m(k))
     enddo
     h_i(nlevs+1) = get_height(p_i(nlevs+1))
   else
     hs = rgas*tmean/g
     h_i(1) = get_height_at_point(tmean,pt)  ! .true.(p->z)
     do k = 1,nlevs
       h_i(k+1) = hs*dlog(p0/p_i(k+1))
       h_m(k)   = hs*dlog(p0/p_m(k))
     enddo
   endif
!
   return
   end subroutine get_pressure_height
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine print_hybrid_coord(nlevs, p_i_in, p_m_in, h_i_in, h_m_in,        &
                                              eta_i, eta_m, a_i, a_m, b_i, b_m)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: nlevs
   real(r8), dimension(nlevs+1), intent(in   ) :: p_i_in, h_i_in, eta_i, a_i, b_i
   real(r8), dimension(nlevs),   intent(in   ) :: p_m_in, h_m_in, eta_m, a_m, b_m
!
   integer(i4)                  :: k
   real(r8), dimension(nlevs+1) :: p_i, h_i
   real(r8), dimension(nlevs)   :: p_m, h_m
   character(len=60)            :: fmtlabel, fmtstrt, fmtstr
!
   p_i(:) = p_i_in(:)/100.0_r8
   p_m(:) = p_m_in(:)/100.0_r8
   h_i(:) = h_i_in(:)/1000.0_r8
   h_m(:) = h_m_in(:)/1000.0_r8
!
   fmtlabel = '(a7,x,a13,x,a15,x,a12,x,a12,x,a12) '
   fmtstrt  = '(a7,x,f13.9,x,f15.10,x,f12.10,x,f12.10,x,f12.10) '
   fmtstr   = '(i3,a4,x,f13.9,x,f15.10,x,f12.10,x,f12.10,x,f12.10) '
!
   write(6,*) ' '
   write(6,'(a)') ' All levels :'
   write(6,fmtlabel) 'Level','Height[km]','Pressure[hPa]','eta','A','B'
   write(6,fmtstrt) '1/2',h_i(1),p_i(1),eta_i(1),a_i(1),b_i(1)
   do k = 1,nlevs
     write(6,fmtstr) k,' ',h_m(k),p_m(k),eta_m(k),a_m(k),b_m(k)
     write(6,fmtstr) k,'+1/2',h_i(k+1),p_i(k+1),eta_i(k+1),a_i(k+1),b_i(k+1)
   enddo
   write(6,'(a)') ' '
!
   write(6,'(a)') ' Half level :'
   write(6,fmtlabel) 'Level','Height[km]','Pressure[hPa]','eta','A','B'
   write(6,fmtstrt) '1/2',h_i(1),p_i(1),eta_i(1),a_i(1),b_i(1)
   do k = 1,nlevs
     write(6,fmtstr) k,'+1/2',h_i(k+1),p_i(k+1),eta_i(k+1),a_i(k+1),b_i(k+1)
   enddo
   write(6,'(a)') ' '
!
   write(6,'(a)') ' Full level :'
   write(6,fmtlabel) 'Level','Height[km]','Pressure[hPa]','eta','A','B'
   do k = 1,nlevs
     write(6,fmtstr) k,' ',h_m(k),p_m(k),eta_m(k),a_m(k),b_m(k)
   enddo
   write(6,'(a)') ' '
!
   return
   end subroutine print_hybrid_coord
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine analysis_hybrid_coord(nlevs, p_i, p_m, h_i, h_m)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: nlevs
   real(r8), dimension(nlevs+1), intent(in   ) :: p_i, h_i
   real(r8), dimension(nlevs),   intent(in   ) :: p_m, h_m
!
   integer(i4) :: ntropo_i, nstrato_i, nmeso_i, nthermo_i
   integer(i4) :: n1km_i, n2km_i, n3km_i, nkm_i
   integer(i4) :: ntropo_m, nstrato_m, nmeso_m, nthermo_m
   integer(i4) :: n1km_m, n2km_m, n3km_m, nkm_m
   real(r8)    :: tropo, strato, meso, thermo
   integer(i4) :: k
   character(len=20) :: fmtchar, fmtcharr
!
   fmtcharr = '(a46,f8.3,x,f8.3)'
!
   write(6,'(a)') ' '
   write(6,'(a46,a8,x,a8,a4)') ' ### surface,top info. : ','half','full','###'
   write(6,fmtcharr) '* pt: p(1) : ',p_i(1)/100.0_r8,p_m(1)/100.0_r8
   write(6,fmtcharr) '* zt: z(1) : ',h_i(1)/1000.0_r8,h_m(1)/1000.0_r8
   write(6,fmtcharr) '* ps:p(-2) : ',p_i(nlevs)/100.0_r8,p_m(nlevs)/100.0_r8
   write(6,fmtcharr) '* hs:h(-2) : ',h_i(nlevs)/100.0_r8,h_m(nlevs)/100.0_r8
   write(6,'(a)') ' '
!
   tropo     = 10000.0_r8
   strato    = 100.0_r8
   meso      = 1.0_r8
   thermo    = 0.0_r8
   ntropo_i  = 0
   nstrato_i = 0
   nmeso_i   = 0
   nthermo_i = 0
   ntropo_m  = 0
   nstrato_m = 0
   nmeso_m   = 0
   nthermo_m = 0
!
   do k = 1,nlevs+1
     if ((p_i(k)>tropo).and.(p_i(k)<=std_ps)) then
       ntropo_i = ntropo_i+1
     elseif ((p_i(k)>strato).and.(p_i(k)<=tropo)) then
       nstrato_i = nstrato_i+1
     elseif ((p_i(k)>meso).and.(p_i(k)<=strato)) then
       nmeso_i = nmeso_i+1
     elseif ((p_i(k)>= thermo).and.(p_i(k)<=meso)) then
       nthermo_i = nthermo_i+1
     else
       write(6,*) 'Check values of pressure....',k,p_i(k)
       call exit(-1)
     endif
   enddo
!
   do k = 1,nlevs
     if ((p_m(k)>tropo).and.(p_m(k)<=std_ps)) then
       ntropo_m = ntropo_m+1
     elseif ((p_m(k)>strato).and.(p_m(k)<=tropo)) then
       nstrato_m = nstrato_m+1
     elseif ((p_m(k)>meso).and.(p_m(k)<=strato)) then
       nmeso_m = nmeso_m+1
     elseif ((p_m(k)>= thermo).and.(p_m(k)<=meso)) then
       nthermo_m = nthermo_m+1
     else
       write(6,*) 'Check values of pressure....',k,p_m(k)
       call exit(-1)
     endif
   enddo
!
   fmtchar = '(a46,i8,x,i8) '
   write(6,'(a)') ' '
   write(6,'(a46,a8,x,a8,a4)') '### pressure-level info.(#) : ','half','full','###'
   write(6,fmtchar) '* troposphere ( 100 hpa < p <= 1000 hpa) : ',ntropo_i,ntropo_m
   write(6,fmtchar) '* stratosphere(   1 hpa < p <=  100 hpa) : ',nstrato_i,nstrato_m
   write(6,fmtchar) '* mesosphere  (0.01 hpa < p <=    1 hpa) : ',nmeso_i,nmeso_m
   write(6,fmtchar) '* thermosphere(   0 hpa <=p <= 0.01 hpa) : ',nthermo_i,nthermo_m
   write(6,'(a) ') ' '
!
   fmtchar = '(a46,f8.2,x,f8.2) '
   write(6,'(a46,a8,x,a8,a4)') '### pressure-level info.(#) : ','half','full','###'
   write(6,fmtchar) '* troposphere ( 100 hpa < p <= 1000 hpa) : ',dble(ntropo_i)/dble(nlevs+1)*100.0_r8,dble(ntropo_m)/dble(nlevs)*100.0_r8
   write(6,fmtchar) '* stratosphere(   1 hpa < p <=  100 hpa) : ',dble(nstrato_i)/dble(nlevs+1)*100.0_r8,dble(nstrato_m)/dble(nlevs)*100.0_r8
   write(6,fmtchar) '* mesosphere  (0.01 hpa < p <=    1 hpa) : ',dble(nmeso_i)/dble(nlevs+1)*100.0_r8,dble(nmeso_m)/dble(nlevs)*100.0_r8
   write(6,fmtchar) '* thermosphere(   0 hpa <=p <= 0.01 hpa) : ',dble(nthermo_i)/dble(nlevs+1)*100.0_r8,dble(nthermo_m)/dble(nlevs)*100.0_r8
   write(6,'(a) '),' '
!
   n1km_i = 0
   n2km_i = 0
   n3km_i = 0
   nkm_i = 0
   n1km_m = 0
   n2km_m = 0
   n3km_m = 0
   nkm_m = 0
!
   do k = 1,nlevs+1
     if ((h_i(k)>= 0.0_r8).and.(h_i(k)<=1000.0_r8)) then
       n1km_i = n1km_i+1
     elseif ((h_i(k)>1000.0_r8).and.(h_i(k)<=2000.0_r8)) then
       n2km_i = n2km_i+1
     elseif ((h_i(k)>2000.0_r8).and.(h_i(k)<=3000.0_r8)) then
       n3km_i = n3km_i+1
     elseif ((h_i(k)>3000.0_r8).and.(h_i(k)<=huge(h_i(k)))) then
       nkm_i = nkm_i+1
     elseif (abs(h_i(k))>=huge(h_i(k))) then
       nkm_i = nkm_i+1
     else
       write(6,*) 'Check values of height....',k,h_i(k)
       call exit(-1)
     endif
   enddo
!
   do k = 1,nlevs
     if ((h_m(k)>= 0.0_r8).and.(h_m(k)<=1000.0_r8)) then
       n1km_m = n1km_m+1
     elseif ((h_m(k)>1000.0_r8).and.(h_m(k)<=2000.0_r8)) then
       n2km_m = n2km_m+1
     elseif ((h_m(k)>2000.0_r8).and.(h_m(k)<=3000.0_r8)) then
       n3km_m = n3km_m+1
     elseif ((h_m(k)>3000.0_r8).and.(h_m(k)<=huge(h_m(k)))) then
       nkm_m = nkm_m+1
     elseif (abs(h_m(k))>= huge(h_m(k))) then
       nkm_m = nkm_m+1
     else
       write(6,*) 'Check values of height....',k,h_m(k)
       call exit(-1)
     endif
   enddo
!
   fmtchar = '(a46,i8,x,i8) '
   write(6,'(a)') ' '
   write(6,'(a46,a8,x,a8,a4)') adjustl('### height-level info.(#) : '),'half','full','###'
   write(6,fmtchar) ' * (0 km <= h <= 1 km) : ',n1km_i,n1km_m
   write(6,fmtchar) ' * (0 km <= h <= 2 km) : ',n1km_i+n2km_i,n1km_m+n2km_m
   write(6,fmtchar) ' * (0 km <= h <= 3 km) : ',n1km_i+n2km_i+n3km_i,n1km_m+n2km_m+n3km_m
   write(6,'(a)') ' '
!
   fmtchar = '(a46,f8.2,x,f8.2) '
   write(6,'(a46,a8,x,a8,a4)') adjustl('### height-level info.(%) : '),'half','full','###'
   write(6,fmtchar) ' * (0 km <= h <= 1 km) : ',dble(n1km_i)/dble(nlevs+1)*100.0_r8,dble(n1km_m)/dble(nlevs)*100.0_r8
   write(6,fmtchar) ' * (0 km <= h <= 2 km) : ',dble(n1km_i+n2km_i)/dble(nlevs+1)*100.0_r8,dble(n1km_m+n2km_m)/dble(nlevs)*100.0_r8
   write(6,fmtchar) ' * (0 km <= h <= 3 km) : ',dble(n1km_i+n2km_i+n3km_i)/dble(nlevs+1)*100.0_r8,dble(n1km_m+n2km_m+n3km_m)/dble(nlevs)*100.0_r8
   write(6,'(a)') ' '
!
   return
   end subroutine analysis_hybrid_coord
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine write_hybrid_coord(nlevs, shyb, ehyb, pt_m, zt_m, tmean, a_i, a_m, b_i, b_m)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: nlevs, shyb, ehyb
   real(r8), dimension(nlevs+1), intent(in   ) :: a_i, b_i
   real(r8), dimension(nlevs),   intent(in   ) :: a_m, b_m
   real(r8),                     intent(in   ) :: pt_m, zt_m, tmean
!
   integer(i4) :: k
   character(len=4)  :: char_lev, char_prs, char_sig
   character(len=90) :: outfilename_a, outfilename_i, outfilename_m
   character(len=90) :: outfilenameasc_i, outfilenameasc_m
   character(len=12) :: kind_coord
   character(len=15) :: char_pt_m, char_zt_m, char_tmean
!
   kind_coord = './data/hybrid'
   write(char_lev,'(i4.4)') nlevs
   write(char_prs,'(i4.4)') shyb
   write(char_sig,'(i4.4)') ehyb
   write(char_tmean,'(i5.5)') int(tmean)
   write(char_pt_m,'(i5.5)') int(pt_m)
   if (pt_m==0.0_r8) then
     char_zt_m = '_Inf_'
   else
     write(char_zt_m,'(i5.5)') int(zt_m/1000.0_r8)
   endif
!
   outfilename_a = './output/hv_'
   outfilename_i = './output/half_'
   outfilename_m = './output/full_'
!
   outfilename_a = trim(outfilename_a)//trim(adjustr(char_lev))//'nlevs_'
   outfilename_i = trim(outfilename_i)//trim(adjustr(char_lev))//'nlevs_'
   outfilename_m = trim(outfilename_m)//trim(adjustr(char_lev))//'nlevs_'
!
   outfilename_a = trim(outfilename_a)//trim(adjustr(char_prs))//'s_'
   outfilename_i = trim(outfilename_i)//trim(adjustr(char_prs))//'s_'
   outfilename_m = trim(outfilename_m)//trim(adjustr(char_prs))//'s_'
!
   outfilename_a = trim(outfilename_a)//trim(adjustr(char_sig))//'e_'
   outfilename_i = trim(outfilename_i)//trim(adjustr(char_sig))//'e_'
   outfilename_m = trim(outfilename_m)//trim(adjustr(char_sig))//'e_'
!
   outfilename_a = trim(outfilename_a)//trim(char_tmean)//'K_'
   outfilename_i = trim(outfilename_i)//trim(char_tmean)//'K_'
   outfilename_m = trim(outfilename_m)//trim(char_tmean)//'K_'
!
   outfilename_a = trim(outfilename_a)//trim(char_pt_m)//'pa_'
   outfilename_i = trim(outfilename_i)//trim(char_pt_m)//'pa_'
   outfilename_m = trim(outfilename_m)//trim(char_pt_m)//'pa_'
!
   outfilename_a = trim(outfilename_a)//trim(char_zt_m)//"km.dat"
   outfilename_i = trim(outfilename_i)//trim(char_zt_m)//"km.dat"
   outfilename_m = trim(outfilename_m)//trim(char_zt_m)//"km.dat"
!
   write(6,*) '* Output file( all level): ',trim(outfilename_a)
   write(6,*) '* Output file(half level): ',trim(outfilename_i)
   write(6,*) '* Output file(full level): ',trim(outfilename_m)
!
   open(unit=70,file=outfilename_a,status='unknown',access='sequential',form='unformatted')
   open(unit=71,file=outfilename_i,status='unknown',access='sequential',form='unformatted')
   open(unit=72,file=outfilename_m,status='unknown',access='sequential',form='unformatted')
!
   outfilenameasc_i = './coefficient/kiaps/kiaps_half_'
   outfilenameasc_m = './coefficient/kiaps/kiaps_full_'
   outfilenameasc_i = trim(outfilenameasc_i)//trim(adjustr(char_lev))//'nlevs_'//trim(char_pt_m)//'pa.ascii'
   outfilenameasc_m = trim(outfilenameasc_m)//trim(adjustr(char_lev))//'nlevs_'//trim(char_pt_m)//'pa.ascii'
!
   write(6,*) '* Output file(half level,ascii): ',trim(outfilenameasc_i)
   write(6,*) '* Output file(full level,ascii): ',trim(outfilenameasc_m)
!
   open(unit=81,file=outfilenameasc_i,status='unknown',form='formatted')
   open(unit=82,file=outfilenameasc_m,status='unknown',form='formatted')
!
   ! all levels
   write(70) pt_m
   write(70) shyb
   write(70) ehyb
   do k = 1,nlevs+1
     write(70) a_i(k)
     write(70) b_i(k)
   enddo
   do k = 1,nlevs
   write(70) a_m(k)
   write(70) b_m(k)
   enddo
   close(70)
!
   ! half levels
   do k = 1,nlevs+1
   write(71) a_i(k)
   write(71) b_i(k)
   write(81,*) a_i(k),b_i(k)
   enddo
   close(71)
   close(81)
!
   ! full levels
   do k = 1,nlevs
     write(72) a_m(k)
     write(72) b_m(k)
     write(82,*) a_m(k),b_m(k)
   enddo
   close(72)
   close(82)
!
   call system('cp '//trim(outfilename_i)//' ./test.dat')
!
   return
   end subroutine write_hybrid_coord
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_hybrid_coord(nlevs, hvfile_i, hvfile_m, shyb, ehyb, a_i, b_i, eta_i, a_m, b_m, eta_m)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                  intent(in   ) :: nlevs
   character(len=*),             intent(in   ) :: hvfile_i, hvfile_m
   integer(i4),                  intent(inout) :: shyb, ehyb
   real(r8), dimension(nlevs+1), intent(inout) :: a_i, b_i, eta_i
   real(r8), dimension(nlevs)  , intent(inout) :: a_m, b_m, eta_m
!
   integer(i4) :: k, err
!
   open(unit=70,file=trim(hvfile_i),access='sequential',form='unformatted',iostat=err)
   if (err.ne.0) then
     print *, 'cannot open file: ', trim(hvfile_i)
   endif
   do k = 1,nlevs+1
     read(70) a_i(k)
     read(70) b_i(k)
     eta_i(k) = a_i(k)+b_i(k)
   enddo
   close(70)
!
   open(unit=71,file=trim(hvfile_m),access='sequential',form='unformatted',iostat=err)
   if (err.ne.0) then
     print *, 'cannot open file: ', trim(hvfile_m)
   endif
   do k = 1,nlevs
     read(71) a_m(k)
     read(71) b_m(k)
     eta_m(k) = a_m(k)+b_m(k)
   enddo
   close(71)
!
   do k = 1,nlevs+1
     if (b_i(k).ne.0.0_r8) then
       shyb = k
       exit
     endif
   enddo
   do k = 2,nlevs+1
     if (a_i(k).eq.0.0_r8) then
       ehyb = k-1
       exit
     endif
   enddo
!
   return
   end subroutine read_hybrid_coord
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module hybrid_coord
!-------------------------------------------------------------------------------
