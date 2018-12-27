!-------------------------------------------------------------------------------
!>
!> @brief
!>  - element module.
!>
!> @date ?????2012
!>  - JH KIM  : First written from the HOMME and was modified for KIAPSGM framework.
!> @date 30JAN2015
!>  - JH KIM  : Added to "dp3d", "dpdiss_ave", "dpdiss_biharmonic", "tensorVisc", "vec_sphere2cart" variables. (CAM-SE 5.3)
!> @date 02JUN2015
!>  - SH CHOI : arrange for KIM_SW
!-------------------------------------------------------------------------------

#include "KIM.h"

MODULE Element

  !```````````````````````````````````````````````````````
  USE KiapsBase,         ONLY: int_kind  => KIM_INT_KIND   &
                              ,real_kind => KIM_REAL8_KIND &
                              ,r8 => KIM_REAL8_KIND &
                              ,long_kind => KIM_LONG_KIND  &
                              ,longdouble_kind => KIM_LONGDOUBLE_KIND
  USE CoordinateSystems, ONLY: spherical_polar_t  &
                              ,cartesian2D_t      &
                              ,cartesian3D_t      &
                              ,distance
  USE Dimensions,        ONLY: np    &
                              ,npsq  &
                              ,nlev  &
                              ,nlevp &
                              ,qsize_d  &
                              ,nc    &
                              ,ntrac_d &
                              ,n_moist &
                              ,n_trace 
  USE Edge,              ONLY: edgedescriptor_t  &
                              ,rotation_t
  USE Gridgraph,         ONLY: gridvertex_t

!#if defined(CORE_SW)
!  USE PhysicalConstants, ONLY: n_moist, n_trace
!#endif
  !```````````````````````````````````````````````````````


  IMPLICIT NONE
  !-------------------------------------------------------
  PRIVATE


#if defined(CORE_SW)
  TYPE, PUBLIC :: elem_state_t
     SEQUENCE

    !------------------------------------------- Previous step
     REAL(r8), DIMENSION(np,np)           :: pds_1     ! the hydrostatic surface pressure for dry air
     REAL(r8), DIMENSION(np,np,nlev+1)    :: mu_1      ! perturbation dry air mass in column 
     REAL(r8), DIMENSION(np,np,nlev+1)    :: mu_1_w    ! perturbation dry air mass in column : full(w) levels
     REAL(r8), DIMENSION(np,np,2,nlev+1)  :: v_1       ! x, y-wind component
     REAL(r8), DIMENSION(np,np,nlev+1)    :: w_1       ! z-wind component
     REAL(r8), DIMENSION(np,np,nlev+1)    :: ph_1      ! perturbation geopotential
     REAL(r8), DIMENSION(np,np,nlev+1)    :: t_1       ! perturbation theta
     REAL(r8), DIMENSION(np,np,nlev+1)    :: ww1       ! OMEGA(mu eta_dot)

    !------------------------------------------- Current step
     REAL(r8), DIMENSION(np,np)           :: pds_2     ! the hydrostatic surface pressure for dry air
     REAL(r8), DIMENSION(np,np,nlev+1)    :: mu_2      ! perturbation dry air mass in column
     REAL(r8), DIMENSION(np,np,nlev+1)    :: mu_2_w    ! perturbation dry air mass in column : full(w) levels
     REAL(r8), DIMENSION(np,np,2,nlev+1)  :: v_2       ! x, y-wind component
     REAL(r8), DIMENSION(np,np,nlev+1)    :: w_2       ! z-wind component
     REAL(r8), DIMENSION(np,np,nlev+1)    :: ph_2      ! perturbation geopotential
     REAL(r8), DIMENSION(np,np,nlev+1)    :: t_2       ! perturbation theta
     REAL(r8), DIMENSION(np,np,nlev+1)    :: t_2_w     ! perturbation theta
     REAL(r8), DIMENSION(np,np,nlev+1)    :: ww        ! omega(mu eta_dot) 
     REAL(r8), DIMENSION(np,np,nlev+1)    :: omega     ! omega(dp/dt)

    !-------------------------------------------
     REAL(r8), DIMENSION(np,np)           :: pds_tend  ! RHS
     REAL(r8), DIMENSION(np,np,2,nlev+1)  :: rv_tend   ! RHS
     REAL(r8), DIMENSION(np,np,nlev+1)    :: rw_tend   ! RHS
     REAL(r8), DIMENSION(np,np,nlev+1)    :: ph_tend   ! RHS
     REAL(r8), DIMENSION(np,np,nlev+1)    :: t_tend    ! RHS
 
    !-------------------------------------------
     REAL(r8), DIMENSION(np,np)           :: pds_tendf ! RHS for forcing
     REAL(r8), DIMENSION(np,np,2,nlev+1)  :: rv_tendf  ! RHS for forcing
     REAL(r8), DIMENSION(np,np,nlev+1)    :: rw_tendf  ! RHS for forcing
     REAL(r8), DIMENSION(np,np,nlev+1)    :: ph_tendf  ! RHS for forcing
     REAL(r8), DIMENSION(np,np,nlev+1)    :: t_tendf   ! RHS for forcing
 
    !-------------------------------------------for debub
     REAL(r8), DIMENSION(np,np,2,nlev+1)  :: v_tend_rk   ! RHS for RK    time ! last one
     REAL(r8), DIMENSION(np,np,2,nlev+1)  :: v_tend_tau  ! RHS for small time ! last one

    !-------------------------------------------
     REAL(r8), DIMENSION(np,np)           :: pds_save  !
     REAL(r8), DIMENSION(np,np,nlev+1)    :: mu_save   !   
     REAL(r8), DIMENSION(np,np,2,nlev+1)  :: v_save    !   
     REAL(r8), DIMENSION(np,np,nlev+1)    :: w_save    !   
     REAL(r8), DIMENSION(np,np,nlev+1)    :: ph_save   !   
     REAL(r8), DIMENSION(np,np,nlev+1)    :: t_save    !   

    !-------------------------------------------
     REAL(r8), DIMENSION(np,np,2,nlev+1)  :: v_init    ! initial x, y-wind component
     REAL(r8), DIMENSION(2,nlev+1)        :: v_base    ! base-state 
     REAL(r8), DIMENSION(nlev+1)          :: z_base    !   
     REAL(r8), DIMENSION(nlev+1)          :: t_base    !   

    !-------------------------------------------
     REAL(r8), DIMENSION(np,np,2,nlev+1)  :: rv        ! coupled
     REAL(r8), DIMENSION(np,np,nlev+1)    :: rw        !   

    !-------------------------------------------
     REAL(r8), DIMENSION(np,np,2,nlev+1)  :: rv_m      ! mean flux in small-step loop
     REAL(r8), DIMENSION(np,np,nlev+1)    :: ww_m

    !-------------------------------------------
     REAL(r8), DIMENSION(np,np)           :: psfc      ! surface pressure ... for dry+moist air
     REAL(r8), DIMENSION(np,np)           :: ptop      ! top     pressure ... for dry+moist air

     REAL(r8), DIMENSION(np,np,nlev+1)    :: pd        ! dry pressure
     REAL(r8), DIMENSION(np,np)           :: pdst      ! bar() : the hydrostatic surface pressure for dry air
     REAL(r8), DIMENSION(np,np)           :: pdsts     ! bar() in small step 
     REAL(r8), DIMENSION(np,np)           :: pdsb      ! base-state
 
     REAL(r8), DIMENSION(np,np,nlev+1)    :: phb       ! base-state geopotential Z staggered
     REAL(r8), DIMENSION(np,np,nlev+1)    :: php       ! full geopotential
     REAL(r8), DIMENSION(np,np,nlev+1)    :: ph0       ! initial geopotential

     REAL(r8), DIMENSION(np,np,nlev+1)    :: mu0       ! for initial dry pressure
     REAL(r8), DIMENSION(np,np,nlev+1)    :: mut       ! half(mass) levels : dp/deta ! total      mu base state dry air mass in column
     REAL(r8), DIMENSION(np,np,nlev+1)    :: mut_w     ! full(w) levels              ! total      mu base state dry air mass in column
     REAL(r8), DIMENSION(np,np,nlev+1)    :: muts      ! 
     REAL(r8), DIMENSION(np,np,nlev+1)    :: muts_w    ! 
     REAL(r8), DIMENSION(np,np,nlev+1)    :: mub       ! base-state mu base state dry air mass in column
     REAL(r8), DIMENSION(np,np,nlev+1)    :: mub_w     ! base-state mu base state dry air mass in column
     REAL(r8), DIMENSION(np,np,nlev+1)    :: muave
     REAL(r8), DIMENSION(np,np,nlev+1)    :: muave_w
     REAL(r8), DIMENSION(np,np       )    :: mudf      ! 

     REAL(r8), DIMENSION(np,np,nlev+1)    :: pb        ! base-state pressure
     REAL(r8), DIMENSION(np,np,nlev+1)    :: pb_w      ! 
     REAL(r8), DIMENSION(np,np,nlev+1)    :: p         ! perturbation pressure
     REAL(r8), DIMENSION(np,np,nlev+1)    :: p_w       ! perturbation pressure
     REAL(r8), DIMENSION(np,np,nlev+1)    :: p_hyd     !
     REAL(r8), DIMENSION(np,np,nlev+1)    :: p_hyd_w   !

    !-------------- 150107
     REAL(r8), DIMENSION(np,np,nlev+1)    :: p_phy     ! [dry+moist] middle    level pressure
     REAL(r8), DIMENSION(np,np,nlev+1)    :: p8w       ! [dry+moist] interface level pressure

     REAL(r8), DIMENSION(np,np,nlev+1)    :: rho       ! total density [dry+moist]
     REAL(r8), DIMENSION(np,np,nlev+1)    :: alt       ! total inverse density [dry]
     REAL(r8), DIMENSION(np,np,nlev+1)    :: alt_w     ! total inverse density [dry]
     REAL(r8), DIMENSION(np,np,nlev+1)    :: alb       ! base-state inverse density [dry]
     REAL(r8), DIMENSION(np,np,nlev+1)    :: alb_w     ! base-state inverse density [dry]
     REAL(r8), DIMENSION(np,np,nlev+1)    :: al        ! perturbation inverse density [dry]
     REAL(r8), DIMENSION(np,np,nlev+1)    :: al_w      ! perturbation inverse density [dry]

     REAL(r8), DIMENSION(np,np,nlev+1)    :: temp_init ! temperature for useInifile
     REAL(r8), DIMENSION(np,np,nlev+1)    :: t_init    ! theta 
     REAL(r8), DIMENSION(np,np,nlev+1)    :: t_dum     ! theta ! mean theta IGW
     REAL(r8), DIMENSION(np,np,nlev+1)    :: t_2save

     REAL(r8), DIMENSION(np,np)           :: ht        ! HGT "Terrain Height"
     REAL(r8), DIMENSION(np,np,2)         :: dht       ! gradient HGT "Terrain Height"

     REAL(r8), DIMENSION(np,np,nlev+1)    :: h_diabatic

    !-------------------------------------------
     REAL(r8), DIMENSION(np,np)           :: sina      !
     REAL(r8), DIMENSION(np,np)           :: cosa      !
     REAL(r8), DIMENSION(np,np)           :: e         !
     REAL(r8), DIMENSION(np,np)           :: f         !

    !------------ i1 variable dyn
     REAL(r8), DIMENSION(np,np,nlev+1)    :: c2a       !
     REAL(r8), DIMENSION(np,np,nlev+1)    :: pm1       !

     REAL(r8), DIMENSION(np,np,nlev+1)    :: alpha     !
     REAL(r8), DIMENSION(np,np,nlev+1)    :: a         !
     REAL(r8), DIMENSION(np,np,nlev+1)    :: gamma     !

     REAL(r8), DIMENSION(np,np,2,nlev+1)  :: cqv       !
     REAL(r8), DIMENSION(np,np,nlev+1)    :: cqw       !

     REAL(r8), DIMENSION(np,np,nlev+1)    :: divv      !
     REAL(r8), DIMENSION(np,np,2)         :: grad_pb   !
     REAL(r8), DIMENSION(np,np,2)         :: grad_phb  !
     REAL(r8), DIMENSION(np,np,2)         :: grad_p    !
     REAL(r8), DIMENSION(np,np,2)         :: grad_ph   !

   !------------ tracer 
     REAL(r8), DIMENSION(np,np,nlev+1,n_trace)  ::  Q_1     ! old
     REAL(r8), DIMENSION(np,np,nlev+1,n_trace)  ::  Q_2     ! current
     REAL(r8), DIMENSION(np,np,nlev+1,n_trace)  ::  Q_tend  ! 
     REAL(r8), DIMENSION(np,np,nlev+1,n_trace)  ::  Q_tendf ! 

   !------------ moist : hydrometeors
     REAL(r8), DIMENSION(np,np,nlev+1,n_moist)  :: moist_1           ! np_x, nlev+1, n_moist=3, currently : "kg kg-1"
     REAL(r8), DIMENSION(np,np,nlev+1,n_moist)  :: moist_2           ! np_x, nlev+1, n_moist=3, currently : "kg kg-1"
     REAL(r8), DIMENSION(np,np,nlev+1,n_moist)  :: moist_tend        ! np_x, nlev+1, n_moist=3, currently
     REAL(r8), DIMENSION(np,np,nlev+1,n_moist)  :: moist_tendf       ! np_x, nlev+1, n_moist=3, currently

   !------------ GRIMS physics
     REAL(r8), DIMENSION(np,np)           :: rn2d
     REAL(r8), DIMENSION(np,np)           :: rc2d
     REAL(r8), DIMENSION(np,np)           :: sn2d
 
  END TYPE elem_state_t


  TYPE, PUBLIC :: elem_diags_t
     SEQUENCE
     REAL(r8), DIMENSION(np,np)         :: olr          ! outgoing longwave rad flux
     REAL(r8), DIMENSION(np,np)         :: toa_sw_up    ! toa upward   shortwave flux
     REAL(r8), DIMENSION(np,np)         :: toa_sw_dn    ! toa downward shortwave flux
     REAL(r8), DIMENSION(np,np)         :: sfc_lw_dn    ! sfc downward longwave flux
     REAL(r8), DIMENSION(np,np)         :: sfc_lw_up    ! sfc upward longwave flux
     REAL(r8), DIMENSION(np,np)         :: sfc_sw_up    ! sfc upward shorwave flux
     REAL(r8), DIMENSION(np,np)         :: sfc_sw_dn    ! sfc downward shortwave flux
     REAL(r8), DIMENSION(np,np)         :: xmu          ! coeff for diurnal cycle
     REAL(r8), DIMENSION(np,np)         :: sfc_shflx    ! sfc sensible heat flux
     REAL(r8), DIMENSION(np,np)         :: sfc_lhflx    ! sfc latent heat flux

     REAL(r8), DIMENSION(np,np,nlev)    :: swhr   ! vertical sw heating rate
     REAL(r8), DIMENSION(np,np,nlev)    :: lwhr   ! vertical lw heating rate
     REAL(r8), DIMENSION(np,np,nlev)    :: ttenpbl! pbl tendency
     REAL(r8), DIMENSION(np,np,nlev)    :: o3     ! ozone vertical distribution

     REAL(r8), DIMENSION(np,np)         :: evap     ! evaporation from surface (lsm output)
     REAL(r8), DIMENSION(np,np)         :: pblh     ! pbl height 

     REAL(r8), DIMENSION(np,np,nlev+1)  :: rh        ! relative humidity
     REAL(r8), DIMENSION(np,np)         :: u10m
     REAL(r8), DIMENSION(np,np)         :: v10m
     REAL(r8), DIMENSION(np,np)         :: t2m
     REAL(r8), DIMENSION(np,np)         :: q2m
     REAL(r8), DIMENSION(np,np)         :: rh2m

     REAL(r8), DIMENSION(np,np)         :: landfrac
     REAL(r8), DIMENSION(np,np)         :: seaice
  END TYPE elem_diags_t

  TYPE, PUBLIC :: derived_state_t
     SEQUENCE
     ! Primitive equations forcings
     !REAL(r8), DIMENSION(np,np,nlev,n_moist+n_trace) :: FQ ! F-Tracers  
     REAL(r8), DIMENSION(np,np,nlev,n_trace)         :: FQ ! F-Tracers  
     REAL(r8), DIMENSION(np,np,nlev,n_moist)         :: FM ! F-Tracers  
     REAL(r8), DIMENSION(np,np,2,nlev)               :: FV       ! F-momentum       
     REAL(r8), DIMENSION(np,np,nlev)                 :: FT         ! F-Potential Temperature   
  END TYPE derived_state_t
#endif

#if defined(CORE_SH)
  integer, public, parameter :: timelevels=3

  type, public :: elem_state_t
     sequence
!
!    note: variables (and sequence) must match that in prim_restart_mod.F90
!
! prognostic variables
!
     real (kind=real_kind) :: v(np,np,2,nlev,timelevels)   ! velocity                                            1
     real (kind=real_kind) :: T(np,np,nlev,timelevels)     ! temperature                                         2
     real (kind=real_kind) :: lnps(np,np,timelevels)       ! log surface pressure                                3
     real (kind=real_kind) :: ps_v(np,np,timelevels)       ! surface pressure on v grid                          4
     real (kind=real_kind) :: phis(np,np)         ! surface geopotential (prescribed)                            5
     ! [ Add for V-Lag.
     ! vetically lagrangian code advects dp instead of ps_v
     real (kind=real_kind) :: dp3d(np,np,nlev,timelevels)  ! delta p on levels
     ! ]
     ! qsize = 1 is related to the mixing ratio
     ! everything else are passive tracers that can eventually
     ! be forced by the column model.
     real (kind=real_kind) :: Q(np,np,nlev,qsize_d,timelevels)  ! Tracer concentration
     real (kind=real_kind) :: Qdp(np,np,nlev,qsize_d,timelevels)  ! Tracer mass           

  end type elem_state_t

  ! JPE: This parameter must match the number of variables in the state
  ! structure all of which are assumed to be of kind=real_kind.  This is a
  ! requirement for restart I/O. 
  ! MT: this is now obsolete?
  integer(kind=int_kind),public,parameter :: StateComponents=7


  TYPE, PUBLIC :: derived_state_t
     sequence

     real (kind=real_kind) :: vn0(np,np,2,nlev)   ! velocity at n0 saved for use by tracers when sub-cycling
     real (kind=real_kind) :: vstar(np,np,2,nlev) ! velocity on Lagrangian surfaces

     real (kind=real_kind) :: phi(np,np,nlev)     ! geopotential                      
     real (kind=real_kind) :: omega_p(np,np,nlev) ! vertical tendency (derived)       
     real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)  ! vertical velocity to be used for tracers

     real (kind=real_kind) :: grad_lnps(np,np,2)  ! gradient of log surface pressure               
     real (kind=real_kind) :: zeta(np,np,nlev)    ! relative vorticity                             
     real (kind=real_kind) :: div(np,np,nlev,timelevels)     ! divergence                          

     real (kind=real_kind) :: dp(np,np,nlev)       
     real (kind=real_kind) :: divdp(np,np,nlev) 
     real (kind=real_kind) :: divdp_proj(np,np,nlev) 

     ! Primitive equations forcings
#ifdef CAM
     real (kind=real_kind) :: FQ(np,np,nlev,qsize_d, 1) ! F-Tracers   
     real (kind=real_kind) :: FM(np,np,2,nlev, 1)       ! F-momentum  
     real (kind=real_kind) :: FT(np,np,nlev, 1)         ! F-Temperature
#else
     ! when does forcing ever need more than 1 timelevel?  
     real (kind=real_kind) :: FQ(np,np,nlev,qsize_d, timelevels) ! F-Tracers  
     real (kind=real_kind) :: FQ_dyn(np,np,nlev,qsize_d)    ! F-Tracers for GPPACK
     real (kind=real_kind) :: FV(np,np,2,nlev, timelevels)  ! F-momentum       
     real (kind=real_kind) :: FT(np,np,nlev, timelevels)    ! F-Temperature   
#endif

     real (kind=real_kind) :: pecnd(np,np,nlev)         ! pressure pert. from condensate
     real (kind=real_kind) :: FQps(np,np,timelevels)         ! implied forcing of FQ on ps_v

! 20150107 jh.kim
     real (kind=real_kind) :: dpdiss_ave(np,np,nlev)
     real (kind=real_kind) :: dpdiss_biharmonic(np,np,nlev)

  END TYPE derived_state_t

  type, public :: elem_accum_t
     sequence
!     real (kind=real_kind) :: u(np,np,nlev)       ! zonal velocity on sphere
!     real (kind=real_kind) :: T(np,np,nlev)       ! temperature
!     real (kind=real_kind) :: Q(np,np,nlev) ! tracers
!     real (kind=real_kind) :: ke(np,np,nlev)      ! kinetic energy

!
! **** ENERGY DIAGNOSTICS ****    
!
! Energy equation:   
! KE_t  = T1 + T2  + D1   + Err   +  vertical & horizontal advection terms
! IE_t  = S1 + D2                 +  vertical & horizontal advection terms
! PE_t  = S2        
!
!  KEvert*  =  KE net vertical advection    (should be zero) 
!  KEhoriz* =  KE net horizonatl advection  (should be zero)
!  IEvert*  =  IE net vertical advection    (should be zero) 
!  IEhoriz* =  IE net horizonatl advection  (should be zero)
!
! With leapfrog, energy equations are all exact except KE 
! (has an Err term that goes to zero as dt**2)
!
! Transfer terms:
! T1   = -< dp/dn u, RT_v/p grad_p >     KE<->IE:   T1 + T2-T2_s = S1
! T2   = -< dp/dn u, grad_phi >          KE<->PE:   T2_s         = S2
! T2_s = -< dp/dn u, grad_phis >
! S1   = < Cp_star dp/dn , RT omega_p/Cp_star >  
! S2   = -< div (u dp/dn), phis >                
!
#ifdef ENERGY_DIAGNOSTICS
     real (kind=real_kind) :: KEvert1(np,np)      ! term from continuity equ
     real (kind=real_kind) :: KEvert2(np,np)      ! term from momentum equ
     real (kind=real_kind) :: IEvert1(np,np)      ! term from continuity equ
     real (kind=real_kind) :: IEvert2(np,np)      ! term from T equ
     real (kind=real_kind) :: IEvert1_wet(np,np)  ! wet term from continuity equ
     real (kind=real_kind) :: IEvert2_wet(np,np)  ! wet term from T equ

     real (kind=real_kind) :: KEhorz1(np,np)      ! at time t
     real (kind=real_kind) :: KEhorz2(np,np)      ! after calling time_advance, these will be at time t-1
     real (kind=real_kind) :: IEhorz1(np,np)
     real (kind=real_kind) :: IEhorz2(np,np)
     real (kind=real_kind) :: IEhorz1_wet(np,np)
     real (kind=real_kind) :: IEhorz2_wet(np,np)

     real (kind=real_kind) :: T1(np,np)
     real (kind=real_kind) :: T2(np,np)
     real (kind=real_kind) :: T2_s(np,np)
     real (kind=real_kind) :: S1(np,np)
     real (kind=real_kind) :: S1_wet(np,np)
     real (kind=real_kind) :: S2(np,np)

     ! the KE conversion term and diffusion term
     real (kind=real_kind) :: DIFF(np,np,2,nlev) ! net hypervis term
     real (kind=real_kind) :: DIFFT(np,np,nlev) ! net hypervis term
     real (kind=real_kind) :: CONV(np,np,2,nlev) ! dpdn u dot CONV = T1 + T2
#endif
!   ^
!   |
! endif  ENERGY_DIAGNOSTICS

     ! the "4" timelevels represents data computed at:
     !  1  t-.5   
     !  2  t+.5   after dynamics
     !  3  t+.5   after forcing
     !  4  t+.5   after Robert
     ! after calling TimeLevelUpdate, all time above decrease by 1.0
     real (kind=real_kind) :: KEner(np,np,4)
     real (kind=real_kind) :: PEner(np,np,4)
     real (kind=real_kind) :: IEner(np,np,4)
     real (kind=real_kind) :: IEner_wet(np,np,4)
     real (kind=real_kind) :: Qvar(np,np,qsize_d,4)  ! Q variance at half time levels   
     real (kind=real_kind) :: Qmass(np,np,qsize_d,4) ! Q mass at half time levels
     real (kind=real_kind) :: Q1mass(np,np,qsize_d)  ! Q mass at full time levels
     real (kind=real_kind) :: mass_added(qsize_d)    ! mass added by qneg fixer

  end type elem_accum_t
#endif

  TYPE, PUBLIC :: index_t
     SEQUENCE
     INTEGER(KIND=int_kind) :: ia(npsq),ja(npsq)
     INTEGER(KIND=int_kind) :: is,ie
     INTEGER(KIND=int_kind) :: NumUniquePts
     INTEGER(KIND=int_kind) :: UniquePtOffset
  END TYPE index_t


  TYPE, PUBLIC :: element_t
     SEQUENCE

     INTEGER(KIND=int_kind) :: LocalId
     INTEGER(KIND=int_kind) :: GlobalId
     
     TYPE (element_t), POINTER :: prev
     TYPE (element_t), POINTER :: next
!#if defined(CORE_SH)
!#endif

     ! Coordinate values of element points
     TYPE (spherical_polar_t) :: spherep(np,np)           ! Spherical coordinates of GLL points

     ! equ-angular gnomonic projection coordinates
     TYPE (cartesian2D_t)     :: cartp(np,np)  ! gnomonic coordinates of GLL points 
     TYPE (cartesian2D_t)     :: corners(4)    ! gnomonic coordinates of element corners
     REAL (KIND=real_kind)    :: u2qmap(4,2)   ! bilinear map from ref element to quad in cubedsphere coordinates

     ! 3D cartesian coordinates
     TYPE (cartesian3D_t)     :: corners3D(4)  

     ! element diagnostics
     REAL (KIND=real_kind)    :: area               ! Area of element
     REAL (KIND=real_kind)    :: max_eig        ! max singular value of metinv
     REAL (KIND=real_kind)    :: min_eig        ! min singular value of metinv
     REAL (KIND=real_kind)    :: dx_short       ! short length scale
     REAL (KIND=real_kind)    :: dx_long        ! long length scale
     REAL (KIND=real_kind)    :: variable_hyperviscosity(np,np)  ! hyperviscosity based on above
     REAL (KIND=real_kind)    :: tensorVisc(2,2,np,np)    !og, matrix V for tensor viscosity
     REAL (KIND=real_kind)    :: courant !  advective courant number
     REAL (KIND=real_kind)    :: hv_courant ! hyperviscosity courant number

     ! Edge connectivity information
     INTEGER(KIND=int_kind)   :: node_numbers(4)
     INTEGER(KIND=int_kind)   :: node_multiplicity(4)      ! number of elements sharing corner node

     TYPE (GridVertex_t)      :: vertex                 ! Element grid vertex information
     TYPE (EdgeDescriptor_t)  :: desc

     TYPE (elem_state_t)      :: state

#if defined(CORE_SW)
!--------------------------------------140930 start
!#if defined(GFS_SFC) /* sjchoi */
!     TYPE (elem_phys_state_t) :: phys_state
!#endif
!sjchoi
     TYPE (elem_diags_t)      :: diags
!--------------------------------------140930 end
#endif

     type (derived_state_t)    :: derived
#if defined(CORE_SH)
     type (elem_accum_t)       :: accum
#endif

     ! Metric terms 
     REAL (KIND=real_kind)    :: met(2,2,np,np)      ! metric tensor on velocity and pressure grid
     REAL (KIND=real_kind)    :: metinv(2,2,np,np)   ! metric tensor on velocity and pressure grid
     REAL (KIND=real_kind)    :: metdet(np,np)       ! g = SQRT(det(g_ij)) on velocity and pressure grid
     REAL (KIND=real_kind)    :: rmetdet(np,np)      ! 1/metdet on velocity pressure grid
     REAL (KIND=real_kind)    :: D(2,2,np,np)        ! Map covariant field on cube to vector field on the sphere
     REAL (KIND=real_kind)    :: Dinv(2,2,np,np)     ! Map vector field on the sphere to covariant v on cube

     ! Convert vector fields from spherical to rectangular components
     ! The transpose of this operation is its pseudoinverse.
     REAL (KIND=real_kind)    :: vec_sphere2cart(np,np,3,2)

     ! Mass matrix terms for an element on a cube face
     REAL (KIND=real_kind)    :: mp(np,np)          ! mass matrix on  velocity and pressure grid
     REAL (KIND=real_kind)    :: rmp(np,np)         ! inverse mass matrix on velocity and pressure grid

     ! Mass matrix terms for an element on the sphere
     ! This mass matrix is used when solving the equations in weak form
     ! with the natural (surface area of the sphere) inner product
     REAL (KIND=real_kind)    :: spheremp(np,np)          ! mass matrix on velocity and pressure grid
     REAL (KIND=real_kind)    :: rspheremp(np,np)         ! inverse mass matrix on velocity and pressure grid

     INTEGER(KIND=long_kind)  :: gdofP(np,np)        ! Global degree of freedom (P-grid)

     ! Coreolis term
     REAL (KIND=real_kind)    :: fcor(np,np)        ! coreolis term !sjchoi 2*omega*sin(phi)
#if defined(CORE_SW)
     REAL (KIND=real_kind)    :: ecor(np,np)        ! coreolis term !sjchoi 2*omega*cos(phi)
#endif

     ! Solver weights (used only for non-staggered grid
     REAL (KIND=real_kind)    :: solver_wts(np,np)

     TYPE (index_t) :: idxP
     TYPE (index_t), POINTER  :: idxV
     INTEGER :: FaceNum

     ! force element_t to be a multiple of 8 bytes.  
     ! on BGP, code will crash (signal 7, or signal 15) if 8 byte alignment is off
     ! check core file for:
     ! core.63:Generated by interrupt..(Alignment Exception DEAR=0xa1ef671c ESR=0x01800000 CCR0=0x4800a002)
     INTEGER :: dummy
  END TYPE element_t


  TYPE, PUBLIC :: eroot_t
     SEQUENCE
     TYPE(element_t), POINTER :: first
  END TYPE eroot_t

  INTEGER, PUBLIC :: NumEdges

  TYPE (eroot_t), PUBLIC :: eroot

  PUBLIC :: element_coordinates
  PUBLIC :: element_var_coordinates
  PUBLIC :: element_var_coordinates3D

  PUBLIC :: LLAddEdge,LLFindEdge, LLInsertEdge
  PUBLIC :: LLSetEdgeCount,LLGetEdgeCount
  PUBLIC :: LLFree


  PUBLIC :: GetColumnIdP,GetColumnIdV

CONTAINS

  ! =======================================
  !  GetColumnIdP:
  !  
  !  Gives a unique identifier for a Physics 
  !  column on the P-grid
  ! =======================================
  function GetColumnIdP(elem,i,j) result(col_id)
    type(element_t), intent(in) :: elem
    integer, intent(in) :: i,j
    integer :: col_id

    col_id = elem%gdofP(i,j)

  end function GetColumnIdP
  ! =======================================
  !  GetColumnIdV:
  !  
  !  Gives a unique identifier for a Physics 
  !  column on the V-grid
  ! =======================================
  function GetColumnIdV(elem,i,j) result(col_id)
    type(element_t), intent(in) :: elem
    integer, intent(in) :: i,j
    integer :: col_id

    col_id = elem%gdofP(i,j)

  end function GetColumnIdV

  ! =======================================
  ! element_coordinates:
  !
  ! Initialize 2D rectilinear element 
  ! colocation points
  !
  ! =======================================

  function element_coordinates(start,end,points) result(cart)
    type (cartesian2D_t), intent(in) :: start
    type (cartesian2D_t), intent(in) :: end
    real (kind=longdouble_kind), intent(in) :: points(:)

    type (cartesian2D_t) :: cart(SIZE(points),SIZE(points))
    type (cartesian2D_t) :: length, centroid
    real (kind=longdouble_kind) :: y 
    integer i,j

    length%x   = 0.50D0*(end%x-start%x)
    length%y   = 0.50D0*(end%y-start%y)
    centroid%x = 0.50D0*(end%x+start%x) 
    centroid%y = 0.50D0*(end%y+start%y) 
    do j=1,SIZE(points)
       y = centroid%y + length%y*points(j)
       do i=1,SIZE(points)
          cart(i,j)%x = centroid%x + length%x*points(i)
          cart(i,j)%y = y
       end do
    end do

  end function element_coordinates

  function element_var_coordinates(c,points) result(cart)
    type (cartesian2D_t), intent(in) :: c(4)
    real (kind=longdouble_kind), intent(in) :: points(:)
    type (cartesian2D_t) :: cart(SIZE(points),SIZE(points))

    real (kind=longdouble_kind) :: p(size(points))
    real (kind=longdouble_kind) :: q(size(points))
    integer i,j

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do j=1,SIZE(points)
       do i=1,SIZE(points)
          cart(i,j)%x = p(i)*p(j)*c(1)%x &
                      + q(i)*p(j)*c(2)%x &
                      + q(i)*q(j)*c(3)%x &
                      + p(i)*q(j)*c(4)%x 
          cart(i,j)%y = p(i)*p(j)*c(1)%y &
                      + q(i)*p(j)*c(2)%y &
                      + q(i)*q(j)*c(3)%y &
                      + p(i)*q(j)*c(4)%y 
       end do
    end do

  end function element_var_coordinates


  function element_var_coordinates3d(c,points) result(cart)
    type (cartesian3D_t), intent(in) :: c(4)
    real (kind=longdouble_kind), intent(in) :: points(:)
    type (cartesian3D_t) :: cart(SIZE(points),SIZE(points))

    real (kind=longdouble_kind) :: p(size(points))
    real (kind=longdouble_kind) :: q(size(points)),r
    integer i,j

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do j=1,SIZE(points)
       do i=1,SIZE(points)
          cart(i,j)%x = p(i)*p(j)*c(1)%x &
                      + q(i)*p(j)*c(2)%x &
                      + q(i)*q(j)*c(3)%x &
                      + p(i)*q(j)*c(4)%x 
          cart(i,j)%y = p(i)*p(j)*c(1)%y &
                      + q(i)*p(j)*c(2)%y &
                      + q(i)*q(j)*c(3)%y &
                      + p(i)*q(j)*c(4)%y 
          cart(i,j)%z = p(i)*p(j)*c(1)%z &
                      + q(i)*p(j)*c(2)%z &
                      + q(i)*q(j)*c(3)%z &
                      + p(i)*q(j)*c(4)%z 

          ! project back to sphere:
          r = distance(cart(i,j))
          cart(i,j)%x = cart(i,j)%x/r
          cart(i,j)%y = cart(i,j)%y/r
          cart(i,j)%z = cart(i,j)%z/r
       end do
    end do

  end function element_var_coordinates3d


  subroutine LLSetEdgeCount(value)
    implicit none
    integer,intent(in)   :: value
    NumEdges=value
  end subroutine LLSetEdgeCount


  subroutine LLGetEdgeCount(value)
    implicit none
    integer,intent(out)  :: value
    value=NumEdges
  end subroutine LLGetEdgeCount


  recursive subroutine copy_node(node2,node1)

    type (element_t), intent(out) :: node2
    type (element_t), intent(in)  :: node1

    node2%LocalId = node1%LocalId
    node2%GlobalId = node1%GlobalId
    node2%prev       = node1%prev
    node2%next       = node1%next

  end subroutine copy_node


  subroutine LLFree(List)

    implicit none
    type(eroot_t) :: List
    type(element_t), pointer :: temp_node
    integer :: nlist,i


    temp_node => List%first
    ! Find the end of the list
    do while(associated(temp_node%next))
       temp_node => temp_node%next
    enddo

    temp_node => temp_node%prev
    !Now step back and deallocate all entries  
    do while(associated(temp_node))
       deallocate(temp_node%next)
       temp_node => temp_node%prev
    enddo

  end subroutine LLFree


  subroutine LLInsertEdge(EdgeList,Gid,Lid)
    type (eroot_t), intent(inout) :: EdgeList
    integer, intent(in)  :: Gid
    integer, intent(inout) :: Lid
    logical :: found

    call LLFindEdge(EdgeList,Gid,found) 
    if(.not. found) then 
       call LLAddEdge(EdgeList,Gid,Lid) 
    endif

  end subroutine LLInsertEdge


  subroutine LLFindEdge(Edge,Gid,found)

    type (eroot_t), intent(in) :: Edge
    integer, intent(in)  :: Gid
    logical, intent(out) :: found

    type (element_t), pointer :: temp_node

    found =.FALSE.

    temp_node => Edge%first
    do while(associated(temp_node) .and. (.not. found))
       if(Gid .eq. temp_node%GlobalId) then 
          found = .TRUE. 
       else
          temp_node => temp_node%next
       endif
    enddo
  end subroutine LLFindEdge


  subroutine LLAddEdge(EdgeList,Gid,Lid)
    type (eroot_t), intent(inout) :: EdgeList
    integer, intent(in)  :: Gid
    integer, intent(out)  :: Lid

    type(element_t), pointer :: temp_node
    type(element_t), pointer  :: new_node
    type(element_t), pointer :: parent

    temp_node => EdgeList%first
    parent    => EdgeList%first

    do while(associated(temp_node))
       parent => temp_node
       temp_node => parent%next
    enddo
    allocate(new_node)
    NumEdges = NumEdges + 1

    new_node%GlobalId=Gid
    new_node%LocalId=NumEdges
    NULLIFY(new_node%next)
    new_node%prev => parent

    if(associated(EdgeList%first)) then
       parent%next => new_node 
    else
       EdgeList%first => new_node 
    endif
    Lid = NumEdges

  end subroutine LLAddEdge


END MODULE Element
