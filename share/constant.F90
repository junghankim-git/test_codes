!-------------------------------------------------------------------------------
   module constant
!-------------------------------------------------------------------------------
!
!  abstract : define various constants
!
!  history log :
!    2015-07-01  in-sun song    initial setup 
!    2016-04-01  suk-jin choi   unify constants of dyn and phy
!    2016-10-01  jung-eun kim   clean-up & apply unified constant
!
!  reference :
!    - bull. amer. meteor. soc. (aug. 1974) vol 55, no.8, pp.926-930 (bam)
!    - emanuel, k. a., atmospheric convection, isbn 0-19-506630-8, oxford (ema)
!
!  variable :  
!    pi_        the ratio of a circle's circumference to its diameter (bam)                
!    g_         gravitational acceleration at sea level (m/s^2) (bam)                   
!    rerth_     mean radius of the earth (m) (bam)                                        
!    rrerth_    inverse of rerth_ [=1.0e0/rerth_]
!    tsdrerth_  earth's sidereal day (bam)                                              
!    omega_     earth's angular speed of rotation (rad/s}) (bam) 
!               [=2.0e0*pi_/tsdrerth_]
!    karman_    von karman constant                                                       
!    sbc_       stefan-boltzmann constant (W/m^2/K^4) (bam)                        
!    avogad_    avogadro number (bam)                                                   
!    boltz_     boltzmann constant (J/K) (bam)                                     
!    mair_      molecular weight of dry air (kg) (bam)                                
!    mvap_      molecular weight of water vapor (kg) (bam)                           
!    r_         universal gas constant (J/K) (bam) [=avogad_*boltz_]
!    rd_        gas constant for dry air (J/K/kg) (bam) [=r_/mair_]
!    rv_        gas constant for water vapor (J/K/kg) (bam) [=r_/mvap_]             
!    cp_        specific heat of dry air at constant pressure (J/K/kg) (bam) 
!    cv_        specific heat of dry air at constant volume (J/K/kg) (bam)
!               [=cp_-rd_]  
!    cvpm_      [=-cv_/cp_]
!    cpocv_     [=cp_/cv_]
!    cvap_      specific heat of water vapor at constant pressure (J/K/kg) (ema)
!    cliq_      specific heat of liqid water at constant pressure (J/K/kg) (ema)
!    cice_      specific heat of ice         at constant pressure (J/K/kg) (ema)
!    hvap_      latent heat of vaporization at 0'C (J/kg) (ema)
!    hfus_      latent heat of fusion at 0'C       (J/kg) (ema)
!    hsub_      latent heat of sublimation         (J/kg) (ema)
!    psat_      saturation vapor pressure at 0'C (Pa)
!    akapa_     [=rd_/cp_]
!    rdorv_     ratio of dry air to water vapor gas constant [=rd_/rv_]
!    rdorvm1_   [=rd_/rv_-1.]
!    rdog_      [=rd_/g_]
!    rvord_     [=rv_/rd_]
!    rvordm1_   [=rv_/rd_-1.]
!    bp_        base state pressure (Pa)
!    bt_        base state sea level temperature (K)
!    blps_      base state temperature difference between 
!               base pres and 1/e of atm depth (K)
!    tiso_      isothermal temperature in stratosphere (K)
!    p0_        reference surface pressure (Pa)
!    t0_        base state potential temperature (K)
!    t0c_       ice/water mix temperature (K)
!    rhoh2o_    water denstiy (kg/m^3)
!    convrad_   [=cal_*1.e4/60.e0]
!    elocp_     [=hvap_/cp_]
!    degrad_    [=180.0e0/pi_]
!
!-------------------------------------------------------------------------------
   use kinds, only : r8, i4
!
   implicit none
!
   public
!
   real(kind=r8),    parameter :: pi_       = 3.141592653589793238462643383279_r8
   real(kind=r8),    parameter :: g_        = 9.80616e+0_r8        
   real(kind=r8),    parameter :: rerth_    = 6.37122e+6_r8        
   real(kind=r8),    parameter :: rrerth_   = 1.0e0_r8/rerth_  
!                                           = 1.5695581066106648e-7
   real(kind=r8),    parameter :: tsdrerth_ = 8.6164e+4_r8       
   real(kind=r8),    parameter :: omega_    = 2.0e0_r8*pi_/tsdrerth_
!                                           = 7.2921235169903748e-5
   real(kind=r8),    parameter :: karman_   = 4.0e-1_r8
   real(kind=r8),    parameter :: sbc_      = 5.67e-8_r8
   real(kind=r8),    parameter :: avogad_   = 6.02214e+26_r8
   real(kind=r8),    parameter :: boltz_    = 1.38065e-23_r8
   real(kind=r8),    parameter :: mair_     = 2.8966e+1_r8
   real(kind=r8),    parameter :: mvap_     = 1.8016e+1_r8
   real(kind=r8),    parameter :: r_        = avogad_*boltz_   
!                                           = 8.314467591e+3
   real(kind=r8),    parameter :: rd_       = r_/mair_         
!                                           = 2.870423113650487e+2
   real(kind=r8),    parameter :: rv_       = r_/mvap_         
!                                           = 4.6150463982015992e+2
   real(kind=r8),    parameter :: cp_       = 1.00464e+3_r8
   real(kind=r8),    parameter :: cv_       = cp_-rd_          
!                                           = 7.1759768863495128e+2
   real(kind=r8),    parameter :: cvpm_     = -cv_/cp_         
!                                           = -7.1428341359586645e-1
   real(kind=r8),    parameter :: cpocv_    = cp_/cv_          
!                                           = 1.4000045093666262e+0
   real(kind=r8),    parameter :: cvap_     = 1.8700e+3_r8
   real(kind=r8),    parameter :: cliq_     = 4.1900e+3_r8
   real(kind=r8),    parameter :: cice_     = 2.1060e+3_r8
   real(kind=r8),    parameter :: hvap_     = 2.5010e+6_r8
   real(kind=r8),    parameter :: hfus_     = 3.3370e+5_r8
   real(kind=r8),    parameter :: hsub_     = 2.8340e+6_r8
   real(kind=r8),    parameter :: psat_     = 6.1078e+2_r8
   real(kind=r8),    parameter :: akapa_    = rd_/cp_          
!                                           = 2.8571658640413355e-1
   real(kind=r8),    parameter :: rdorv_    = rd_/rv_          
!                                           = 6.2197058620451562e-1
   real(kind=r8),    parameter :: rdorvm1_  = rd_/rv_-1._r8    
!                                           = -3.7802941379548438e-1
   real(kind=r8),    parameter :: rdog_     = rd_/g_           
!                                           = 2.9271632460111675e+1
   real(kind=r8),    parameter :: rvord_    = rv_/rd_          
!                                           = 1.6077930728241563e+0
   real(kind=r8),    parameter :: rvrdm1_   = 0.6077338_r8
!
   real(kind=r8),    parameter :: rvordm1_  = rv_/rd_-1._r8    
!                                           = 6.0779307282415629e-1
   real(kind=r8),    parameter :: fv_       = rvordm1_
   real(kind=r8),    parameter :: bp_       = 1.0e+5_r8
   real(kind=r8),    parameter :: bt_       = 2.9e+2_r8
   real(kind=r8),    parameter :: blps_     = 5.0e+1_r8
   real(kind=r8),    parameter :: tiso_     = 2.0e+2_r8
   real(kind=r8),    parameter :: p0_       = 1.0e+5_r8
   real(kind=r8),    parameter :: t0_       = 3.0e+2_r8
   integer(kind=i4), parameter :: iqv_      = 1
!
   real(kind=r8),    parameter :: solr_     = 1.3533e+3_r8
   real(kind=r8),    parameter :: t0c_      = 2.7315e+2_r8
   real(kind=r8),    parameter :: ttp_      = 2.7316e+2_r8
   real(kind=r8),    parameter :: t0cs_     = 2.7135e+2_r8
   real(kind=r8),    parameter :: cal_      = 4.1855e+0_r8
   real(kind=r8),    parameter :: convrad_  = cal_*1.e4_r8/60.e0_r8  
!                                           = 6.9758333333333333e+2
   real(kind=r8),    parameter :: rhoh2o_   = 1.0e+3_r8
   real(kind=r8),    parameter :: qmin_     = 1.0e-30_r8
   real(kind=r8),    parameter :: qmin8_    = 1.0e-8_r8
   real(kind=r8),    parameter :: qmin9_    = 1.0e-9_r8
   real(kind=r8),    parameter :: qmin30_   = 1.0e-30_r8
   real(kind=r8),    parameter :: aday_     = 8.64e+4_r8
   real(kind=r8),    parameter :: elocp_    = hvap_/cp_        
!                                           = 2.4894489568402614e+3
   real(kind=r8),    parameter :: rhoair0_  = 1.28e0_r8
   real(kind=r8),    parameter :: rhosnow_  = 1.00e+2_r8
   real(kind=r8),    parameter :: cb2pa_    = 1.00e+3_r8
   real(kind=r8),    parameter :: degrad_   = 180.0e0_r8/pi_   
!                                           = 5.7295779513082323e+1
   real(kind=r8),    parameter :: daysec_   = 1.1574e-5_r8
!
   end module constant
!-------------------------------------------------------------------------------
