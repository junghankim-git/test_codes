!-------------------------------------------------------------------------------
   program genvcoord
!-------------------------------------------------------------------------------
!
!  abstract : main program of the coupled system (atm, wav, ocn, ... + cpl)
!
!  history log :
!    2014-04-15   junghan kim     initial setup
!    2014-08-15   junghan kim     fully working version
!    2017-02-13   junghan kim     code clean-up
!    2017-06-20   junghan kim     code refactoring
!
!-------------------------------------------------------------------------------
   use kinds,        only: i4, r8, l4
   use constant,     only: p0=>p0_
   use hybrid_coord, only: hyb_coord_t, initialize, get_eta, get_ab, driver, finalize
!-------------------------------------------------------------------------------
   implicit none
!
!
! nlevs       :the number of vertical levels
! nprs        :the number of pure pressure levels
! nsig        :the number of sigma levels
! pt          :the top pressure
! zt          :the top height
! tmean       :constant reference atmospheric temperature
!              tmean = 222.56k is the temperature such that the top-most full level of
! coefpoly    :polynomial coefficients for the polynomial fit for eta
!             :eta = polynomial of normalized vertical index
! smooth      :integer flag for smoothing
!             :1 = straight line
!             :2 = a, b by polynomial fitting
!             :3 = a, b by differential equation
! method      :1 = logical flag if cubic-spline interpolation based on ecmwf eta is used
!             :2 = polynomial
!             :3 = 
! use_std_atm :logical flag if the standard atmosphere is used
!
   integer(i4)             :: nlevs, shyb, ehyb, nlevs_ecmwf, k
   real(r8)                :: pt, tmean, top_hpa, sfc_m
   real(r8), dimension(15) :: coefpoly
   logical(l4)             :: adjust_top, adjust_sfc
   logical(l4)             :: use_std_atm
   type(hyb_coord_t)       :: hvc
!
   nlevs       = 50
   shyb        = -1
!
   use_std_atm = .true.
   pt          = 0.0_r8
   tmean       = 222.56_r8
   adjust_top  = .true.
   adjust_sfc  = .true.
   top_hpa     = 0.3_r8  ! hPa
   sfc_m       = 10.0_r8 ! m
!
   nlevs_ecmwf = 137
!
   coefpoly    = 0.0_r8
!
#if 0
   !call initialize(hvc,nlevs,shyb,adjust_top,adjust_sfc,top_hpa,sfc_m,use_std_atm,pt,tmean)
   call initialize(hvc,nlevs,adjust_top,adjust_sfc)
   call get_eta(hvc,nlevs_ecmwf)
   call get_ab(hvc,nlevs_ecmwf)
#endif
#if 0
! test
   call initialize(hvc,nlevs,adjust_top,adjust_sfc)
   !                    0.01    1.00   100.0     1000.0
   call get_eta(hvc,(/0.99_r8,1.98_r8,40.59_r8,56.43_r8/))
   call get_ab(hvc)
#endif
!#if 0
! just check and print
   !call initialize(hvc,100,adjust_top,adjust_sfc)
   !call get_eta(hvc,'/home/jhkim/TestBed/KIM/3.0.micros/3.0.04/03.dev/Exp_10h/SW_G/vcoord/Half_0100nlevs_0030plev_0101slev_00222K_00030pa_00057km.dat','/home/jhkim/TestBed/KIM/3.0.micros/3.0.04/03.dev/Exp_10h/SW_G/vcoord/Full_0100nlevs_0030plev_0101slev_00222K_00030pa_00057km.dat')
   !call initialize(hvc,50,adjust_top,adjust_sfc)
   !call get_eta(hvc,'/home/jhkim/TestBed/KIM/3.0.micros/3.0.04/03.dev/Exp_10h/SW_G/vcoord/half_50_0.3.dat','/home/jhkim/TestBed/KIM/3.0.micros/3.0.04/03.dev/Exp_10h/SW_G/vcoord/full_50_0.3.dat')
   !call get_eta(hvc,'/home/yslee/half_50_0.3.dat','/home/yslee/full_50_0.3.dat')
   call initialize(hvc,91,adjust_top,adjust_sfc)
   call get_eta(hvc,'/home/jhkim/TestBed/KIM/3.0.micros/3.0.11/26.test/exp_10h/vcoord/half_ecmwf91.dat','/home/jhkim/TestBed/KIM/3.0.micros/3.0.11/26.test/exp_10h/vcoord/full_ecmwf91.dat')
!#endif
!
   call driver(hvc)
print *, '# half level (a, b, eta)'
do k = 1,50+1
print *, hvc%a_i(k), hvc%b_i(k), hvc%eta_i(k)
enddo
print *, ' '
print *, '# full level (a, b, eta) '
do k = 1,50
print *, hvc%a_m(k), hvc%b_m(k), hvc%eta_m(k)
enddo
print *, ' '
   call finalize(hvc)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_exp_case(jexp)
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: jexp
!
   return
   end subroutine set_exp_case
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end
!-------------------------------------------------------------------------------
