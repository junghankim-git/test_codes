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
!-------------------------------------------------------------------------------
   module element
!```````````````````````````````````````````````````````
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
   use kiapsbase, only: int_kind=>kim_int_kind &
                                                , real_kind=>kim_real8_kind &
                                                       , r8=>kim_real8_kind &
                                                 , long_kind=>kim_long_kind &
                                     , longdouble_kind=>kim_longdouble_kind
   use coordinatesystems, only: spherical_polar_t &
                                                               , cartesian2d_t &
                                                               , cartesian3d_t &
                                                                    , distance
   use dimensions, only: np &
                                                                        , npsq &
                                                                        , nlev &
                                                                       , nlevp &
                                                                     , qsize_d &
                                                                          , nc &
                                                                     , ntrac_d &
                                                                     , n_moist &
                                                                     , n_trace
   use edge, only: edgedescriptor_t &
                                                                  , rotation_t
   use gridgraph, only: gridvertex_t
!#if defined(CORE_SW)
!  USE PhysicalConstants, ONLY: n_moist, n_trace
!#endif
!```````````````````````````````````````````````````````
   implicit none
!-------------------------------------------------------
!
   private
#if defined(CORE_SW)
!
   type, public :: elem_state_t
     sequence
!------------------------------------------- Previous step
     real(r8), dimension(np, np) :: pds_1 ! the hydrostatic surface pressure for dry air
     real(r8), dimension(np, np, nlev+1) :: mu_1 ! perturbation dry air mass in column
     real(r8), dimension(np, np, nlev+1) :: mu_1_w ! perturbation dry air mass in column:full(w) levels
     real(r8), dimension(np, np, 2, nlev+1) :: v_1 ! x, y-wind component
     real(r8), dimension(np, np, nlev+1) :: w_1 ! z-wind component
     real(r8), dimension(np, np, nlev+1) :: ph_1 ! perturbation geopotential
     real(r8), dimension(np, np, nlev+1) :: t_1 ! perturbation theta
     real(r8), dimension(np, np, nlev+1) :: ww1 ! omega(mu eta_dot)
!------------------------------------------- Current step
     real(r8), dimension(np, np) :: pds_2 ! the hydrostatic surface pressure for dry air
     real(r8), dimension(np, np, nlev+1) :: mu_2 ! perturbation dry air mass in column
     real(r8), dimension(np, np, nlev+1) :: mu_2_w ! perturbation dry air mass in column:full(w) levels
     real(r8), dimension(np, np, 2, nlev+1) :: v_2 ! x, y-wind component
     real(r8), dimension(np, np, nlev+1) :: w_2 ! z-wind component
     real(r8), dimension(np, np, nlev+1) :: ph_2 ! perturbation geopotential
     real(r8), dimension(np, np, nlev+1) :: t_2 ! perturbation theta
     real(r8), dimension(np, np, nlev+1) :: t_2_w ! perturbation theta
     real(r8), dimension(np, np, nlev+1) :: ww ! omega(mu eta_dot)
     real(r8), dimension(np, np, nlev+1) :: omega ! omega(dp/dt)
!-------------------------------------------
     real(r8), dimension(np, np) :: pds_tend ! rhs
     real(r8), dimension(np, np, 2, nlev+1) :: rv_tend ! rhs
     real(r8), dimension(np, np, nlev+1) :: rw_tend ! rhs
     real(r8), dimension(np, np, nlev+1) :: ph_tend ! rhs
     real(r8), dimension(np, np, nlev+1) :: t_tend ! rhs
!-------------------------------------------
     real(r8), dimension(np, np) :: pds_tendf ! rhs for forcing
     real(r8), dimension(np, np, 2, nlev+1) :: rv_tendf ! rhs for forcing
     real(r8), dimension(np, np, nlev+1) :: rw_tendf ! rhs for forcing
     real(r8), dimension(np, np, nlev+1) :: ph_tendf ! rhs for forcing
     real(r8), dimension(np, np, nlev+1) :: t_tendf ! rhs for forcing
!-------------------------------------------for debub
     real(r8), dimension(np, np, 2, nlev+1) :: v_tend_rk ! rhs for rk time ! last one
     real(r8), dimension(np, np, 2, nlev+1) :: v_tend_tau ! rhs for small time ! last one
!-------------------------------------------
     real(r8), dimension(np, np) :: pds_save !
     real(r8), dimension(np, np, nlev+1) :: mu_save !
     real(r8), dimension(np, np, 2, nlev+1) :: v_save !
     real(r8), dimension(np, np, nlev+1) :: w_save !
     real(r8), dimension(np, np, nlev+1) :: ph_save !
     real(r8), dimension(np, np, nlev+1) :: t_save !
!-------------------------------------------
     real(r8), dimension(np, np, 2, nlev+1) :: v_init ! initial x, y-wind component
     real(r8), dimension(2, nlev+1) :: v_base ! base-state
     real(r8), dimension(nlev+1) :: z_base !
     real(r8), dimension(nlev+1) :: t_base !
!-------------------------------------------
     real(r8), dimension(np, np, 2, nlev+1) :: rv ! coupled
     real(r8), dimension(np, np, nlev+1) :: rw !
!-------------------------------------------
     real(r8), dimension(np, np, 2, nlev+1) :: rv_m ! mean flux in small-step loop
     real(r8), dimension(np, np, nlev+1) :: ww_m
!-------------------------------------------
     real(r8), dimension(np, np) :: psfc ! surface pressure ... for dry+moist air
     real(r8), dimension(np, np) :: ptop ! top pressure ... for dry+moist air
     real(r8), dimension(np, np, nlev+1) :: pd ! dry pressure
     real(r8), dimension(np, np) :: pdst ! bar():the hydrostatic surface pressure for dry air
     real(r8), dimension(np, np) :: pdsts ! bar() in small step
     real(r8), dimension(np, np) :: pdsb ! base-state
     real(r8), dimension(np, np, nlev+1) :: phb ! base-state geopotential z staggered
     real(r8), dimension(np, np, nlev+1) :: php ! full geopotential
     real(r8), dimension(np, np, nlev+1) :: ph0 ! initial geopotential
     real(r8), dimension(np, np, nlev+1) :: mu0 ! for initial dry pressure
     real(r8), dimension(np, np, nlev+1) :: mut ! half(mass) levels:dp/deta ! total mu base state dry air mass in column
     real(r8), dimension(np, np, nlev+1) :: mut_w ! full(w) levels ! total mu base state dry air mass in column
     real(r8), dimension(np, np, nlev+1) :: muts !
     real(r8), dimension(np, np, nlev+1) :: muts_w !
     real(r8), dimension(np, np, nlev+1) :: mub ! base-state mu base state dry air mass in column
     real(r8), dimension(np, np, nlev+1) :: mub_w ! base-state mu base state dry air mass in column
     real(r8), dimension(np, np, nlev+1) :: muave
     real(r8), dimension(np, np, nlev+1) :: muave_w
     real(r8), dimension(np, np) :: mudf !
     real(r8), dimension(np, np, nlev+1) :: pb ! base-state pressure
     real(r8), dimension(np, np, nlev+1) :: pb_w !
     real(r8), dimension(np, np, nlev+1) :: p ! perturbation pressure
     real(r8), dimension(np, np, nlev+1) :: p_w ! perturbation pressure
     real(r8), dimension(np, np, nlev+1) :: p_hyd !
     real(r8), dimension(np, np, nlev+1) :: p_hyd_w !
!-------------- 150107
     real(r8), dimension(np, np, nlev+1) :: p_phy ! [dry+moist] middle level pressure
     real(r8), dimension(np, np, nlev+1) :: p8w ! [dry+moist] interface level pressure
     real(r8), dimension(np, np, nlev+1) :: rho ! total density [dry+moist]
     real(r8), dimension(np, np, nlev+1) :: alt ! total inverse density [dry]
     real(r8), dimension(np, np, nlev+1) :: alt_w ! total inverse density [dry]
     real(r8), dimension(np, np, nlev+1) :: alb ! base-state inverse density [dry]
     real(r8), dimension(np, np, nlev+1) :: alb_w ! base-state inverse density [dry]
     real(r8), dimension(np, np, nlev+1) :: al ! perturbation inverse density [dry]
     real(r8), dimension(np, np, nlev+1) :: al_w ! perturbation inverse density [dry]
     real(r8), dimension(np, np, nlev+1) :: temp_init ! temperature for useinifile
     real(r8), dimension(np, np, nlev+1) :: t_init ! theta
     real(r8), dimension(np, np, nlev+1) :: t_dum ! theta ! mean theta igw
     real(r8), dimension(np, np, nlev+1) :: t_2save
     real(r8), dimension(np, np) :: ht ! hgt "terrain height"
     real(r8), dimension(np, np, 2) :: dht ! gradient hgt "terrain height"
     real(r8), dimension(np, np, nlev+1) :: h_diabatic
!-------------------------------------------
     real(r8), dimension(np, np) :: sina !
     real(r8), dimension(np, np) :: cosa !
     real(r8), dimension(np, np) :: e !
     real(r8), dimension(np, np) :: f !
!------------ i1 variable dyn
     real(r8), dimension(np, np, nlev+1) :: c2a !
     real(r8), dimension(np, np, nlev+1) :: pm1 !
     real(r8), dimension(np, np, nlev+1) :: alpha !
     real(r8), dimension(np, np, nlev+1) :: a !
     real(r8), dimension(np, np, nlev+1) :: gamma !
     real(r8), dimension(np, np, 2, nlev+1) :: cqv !
     real(r8), dimension(np, np, nlev+1) :: cqw !
     real(r8), dimension(np, np, nlev+1) :: divv !
     real(r8), dimension(np, np, 2) :: grad_pb !
     real(r8), dimension(np, np, 2) :: grad_phb !
     real(r8), dimension(np, np, 2) :: grad_p !
     real(r8), dimension(np, np, 2) :: grad_ph !
!------------ tracer
     real(r8), dimension(np, np, nlev+1, n_trace) :: q_1 ! old
     real(r8), dimension(np, np, nlev+1, n_trace) :: q_2 ! current
     real(r8), dimension(np, np, nlev+1, n_trace) :: q_tend !
     real(r8), dimension(np, np, nlev+1, n_trace) :: q_tendf !
!------------ moist : hydrometeors
     real(r8), dimension(np, np, nlev+1, n_moist) :: moist_1 ! np_x, nlev+1, n_moist = 3, currently:"kg kg-1"
     real(r8), dimension(np, np, nlev+1, n_moist) :: moist_2 ! np_x, nlev+1, n_moist = 3, currently:"kg kg-1"
     real(r8), dimension(np, np, nlev+1, n_moist) :: moist_tend ! np_x, nlev+1, n_moist = 3, currently
     real(r8), dimension(np, np, nlev+1, n_moist) :: moist_tendf ! np_x, nlev+1, n_moist = 3, currently
!------------ GRIMS physics
     real(r8), dimension(np, np) :: rn2d
     real(r8), dimension(np, np) :: rc2d
     real(r8), dimension(np, np) :: sn2d
   end type elem_state_t
!
   type, public :: elem_diags_t
     sequence
     real(r8), dimension(np, np) :: olr ! outgoing longwave rad flux
     real(r8), dimension(np, np) :: toa_sw_up ! toa upward shortwave flux
     real(r8), dimension(np, np) :: toa_sw_dn ! toa downward shortwave flux
     real(r8), dimension(np, np) :: sfc_lw_dn ! sfc downward longwave flux
     real(r8), dimension(np, np) :: sfc_lw_up ! sfc upward longwave flux
     real(r8), dimension(np, np) :: sfc_sw_up ! sfc upward shorwave flux
     real(r8), dimension(np, np) :: sfc_sw_dn ! sfc downward shortwave flux
     real(r8), dimension(np, np) :: xmu ! coeff for diurnal cycle
     real(r8), dimension(np, np) :: sfc_shflx ! sfc sensible heat flux
     real(r8), dimension(np, np) :: sfc_lhflx ! sfc latent heat flux
     real(r8), dimension(np, np, nlev) :: swhr ! vertical sw heating rate
     real(r8), dimension(np, np, nlev) :: lwhr ! vertical lw heating rate
     real(r8), dimension(np, np, nlev) :: ttenpbl! pbl tendency
     real(r8), dimension(np, np, nlev) :: o3 ! ozone vertical distribution
     real(r8), dimension(np, np) :: evap ! evaporation from surface(lsm output)
     real(r8), dimension(np, np) :: pblh ! pbl height
     real(r8), dimension(np, np, nlev+1) :: rh ! relative humidity
     real(r8), dimension(np, np) :: u10m
     real(r8), dimension(np, np) :: v10m
     real(r8), dimension(np, np) :: t2m
     real(r8), dimension(np, np) :: q2m
     real(r8), dimension(np, np) :: rh2m
     real(r8), dimension(np, np) :: landfrac
     real(r8), dimension(np, np) :: seaice
   end type elem_diags_t
!
   type, public :: derived_state_t
     sequence
! Primitive equations forcings
!REAL(r8), DIMENSION(np,np,nlev,n_moist+n_trace) :: FQ ! F-Tracers
     real(r8), dimension(np, np, nlev, n_trace) :: fq ! f-tracers
     real(r8), dimension(np, np, nlev, n_moist) :: fm ! f-tracers
     real(r8), dimension(np, np, 2, nlev) :: fv ! f-momentum
     real(r8), dimension(np, np, nlev) :: ft ! f-potential temperature
   end type derived_state_t
#endif
#if defined(CORE_SH)
   integer, public, parameter :: timelevels = 3
!
   type, public :: elem_state_t
     sequence
!
!    note: variables (and sequence) must match that in prim_restart_mod.F90
!
! prognostic variables
!
     real(real_kind) :: v(np, np, 2, nlev, timelevels) ! velocity 1
     real(real_kind) :: t(np, np, nlev, timelevels) ! temperature 2
     real(real_kind) :: lnps(np, np, timelevels) ! log surface pressure 3
     real(real_kind) :: ps_v(np, np, timelevels) ! surface pressure on v grid 4
     real(real_kind) :: phis(np, np) ! surface geopotential(prescribed) 5
! [ Add for V-Lag.
! vetically lagrangian code advects dp instead of ps_v
     real(real_kind) :: dp3d(np, np, nlev, timelevels) ! delta p on levels
! ]
! qsize = 1 is related to the mixing ratio
! everything else are passive tracers that can eventually
! be forced by the column model.
     real(real_kind) :: q(np, np, nlev, qsize_d, timelevels) ! tracer concentration
     real(real_kind) :: qdp(np, np, nlev, qsize_d, timelevels) ! tracer mass
   end type elem_state_t
! JPE: This parameter must match the number of variables in the state
! structure all of which are assumed to be of kind=real_kind.  This is a
! requirement for restart I/O.
! MT: this is now obsolete?
   integer(int_kind), public, parameter :: statecomponents = 7
!
   type, public :: derived_state_t
     sequence
     real(real_kind) :: vn0(np, np, 2, nlev) ! velocity at n0 saved for use by tracers when sub-cycling
     real(real_kind) :: vstar(np, np, 2, nlev) ! velocity on lagrangian surfaces
     real(real_kind) :: phi(np, np, nlev) ! geopotential
     real(real_kind) :: omega_p(np, np, nlev) ! vertical tendency(derived)
     real(real_kind) :: eta_dot_dpdn(np, np, nlevp) ! vertical velocity to be used for tracers
     real(real_kind) :: grad_lnps(np, np, 2) ! gradient of log surface pressure
     real(real_kind) :: zeta(np, np, nlev) ! relative vorticity
     real(real_kind) :: div(np, np, nlev, timelevels) ! divergence
     real(real_kind) :: dp(np, np, nlev)
     real(real_kind) :: divdp(np, np, nlev)
     real(real_kind) :: divdp_proj(np, np, nlev)
! Primitive equations forcings
#ifdef CAM
     real(real_kind) :: fq(np, np, nlev, qsize_d, 1) ! f-tracers
     real(real_kind) :: fm(np, np, 2, nlev, 1) ! f-momentum
     real(real_kind) :: ft(np, np, nlev, 1) ! f-temperature
#else
! when does forcing ever need more than 1 timelevel?
     real(real_kind) :: fq(np, np, nlev, qsize_d, timelevels) ! f-tracers
     real(real_kind) :: fq_dyn(np, np, nlev, qsize_d) ! f-tracers for gppack
     real(real_kind) :: fv(np, np, 2, nlev, timelevels) ! f-momentum
     real(real_kind) :: ft(np, np, nlev, timelevels) ! f-temperature
#endif
     real(real_kind) :: pecnd(np, np, nlev) ! pressure pert. from condensate
     real(real_kind) :: fqps(np, np, timelevels) ! implied forcing of fq on ps_v
! 20150107 jh.kim
     real(real_kind) :: dpdiss_ave(np, np, nlev)
     real(real_kind) :: dpdiss_biharmonic(np, np, nlev)
   end type derived_state_t
!
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
     real(real_kind) :: kevert1(np, np) ! term from continuity equ
     real(real_kind) :: kevert2(np, np) ! term from momentum equ
     real(real_kind) :: ievert1(np, np) ! term from continuity equ
     real(real_kind) :: ievert2(np, np) ! term from t equ
     real(real_kind) :: ievert1_wet(np, np) ! wet term from continuity equ
     real(real_kind) :: ievert2_wet(np, np) ! wet term from t equ
     real(real_kind) :: kehorz1(np, np) ! at time t
     real(real_kind) :: kehorz2(np, np) ! after calling time_advance, these will be at time t-1
     real(real_kind) :: iehorz1(np, np)
     real(real_kind) :: iehorz2(np, np)
     real(real_kind) :: iehorz1_wet(np, np)
     real(real_kind) :: iehorz2_wet(np, np)
     real(real_kind) :: t1(np, np)
     real(real_kind) :: t2(np, np)
     real(real_kind) :: t2_s(np, np)
     real(real_kind) :: s1(np, np)
     real(real_kind) :: s1_wet(np, np)
     real(real_kind) :: s2(np, np)
! the KE conversion term and diffusion term
     real(real_kind) :: diff(np, np, 2, nlev) ! net hypervis term
     real(real_kind) :: difft(np, np, nlev) ! net hypervis term
     real(real_kind) :: conv(np, np, 2, nlev) ! dpdn u dot conv = t1+t2
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
     real(real_kind) :: kener(np, np, 4)
     real(real_kind) :: pener(np, np, 4)
     real(real_kind) :: iener(np, np, 4)
     real(real_kind) :: iener_wet(np, np, 4)
     real(real_kind) :: qvar(np, np, qsize_d, 4) ! q variance at half time levels
     real(real_kind) :: qmass(np, np, qsize_d, 4) ! q mass at half time levels
     real(real_kind) :: q1mass(np, np, qsize_d) ! q mass at full time levels
     real(real_kind) :: mass_added(qsize_d) ! mass added by qneg fixer
   end type elem_accum_t
#endif
!
   type, public :: index_t
     sequence
     integer(int_kind) :: ia(npsq), ja(npsq)
     integer(int_kind) :: is, ie
     integer(int_kind) :: numuniquepts
     integer(int_kind) :: uniqueptoffset
   end type index_t
!
   type, public :: element_t
     sequence
     integer(int_kind) :: localid
     integer(int_kind) :: globalid
     type(element_t), pointer :: prev
     type(element_t), pointer :: next
!#if defined(CORE_SH)
!#endif
! Coordinate values of element points
     type(spherical_polar_t) :: spherep(np, np) ! spherical coordinates of gll points
! equ-angular gnomonic projection coordinates
     type(cartesian2d_t) :: cartp(np, np) ! gnomonic coordinates of gll points
     type(cartesian2d_t) :: corners(4) ! gnomonic coordinates of element corners
     real(real_kind) :: u2qmap(4, 2) ! bilinear map from ref element to quad in cubedsphere coordinates
! 3D cartesian coordinates
     type(cartesian3d_t) :: corners3d(4)
! element diagnostics
     real(real_kind) :: area ! area of element
     real(real_kind) :: max_eig ! max singular value of metinv
     real(real_kind) :: min_eig ! min singular value of metinv
     real(real_kind) :: dx_short ! short length scale
     real(real_kind) :: dx_long ! long length scale
     real(real_kind) :: variable_hyperviscosity(np, np) ! hyperviscosity based on above
     real(real_kind) :: tensorvisc(2, 2, np, np) !og, matrix v for tensor viscosity
     real(real_kind) :: courant ! advective courant number
     real(real_kind) :: hv_courant ! hyperviscosity courant number
! Edge connectivity information
     integer(int_kind) :: node_numbers(4)
     integer(int_kind) :: node_multiplicity(4) ! number of elements sharing corner node
     type(gridvertex_t) :: vertex ! element grid vertex information
     type(edgedescriptor_t) :: desc
     type(elem_state_t) :: state
#if defined(CORE_SW)
!--------------------------------------140930 start
!#if defined(GFS_SFC) /* sjchoi */
!     TYPE (elem_phys_state_t) :: phys_state
!#endif
!sjchoi
     type(elem_diags_t) :: diags
!--------------------------------------140930 end
#endif
     type(derived_state_t) :: derived
#if defined(CORE_SH)
     type(elem_accum_t) :: accum
#endif
! Metric terms
     real(real_kind) :: met(2, 2, np, np) ! metric tensor on velocity and pressure grid
     real(real_kind) :: metinv(2, 2, np, np) ! metric tensor on velocity and pressure grid
     real(real_kind) :: metdet(np, np) ! g = sqrt(det(g_ij)) on velocity and pressure grid
     real(real_kind) :: rmetdet(np, np) ! 1/metdet on velocity pressure grid
     real(real_kind) :: d(2, 2, np, np) ! map covariant field on cube to vector field on the sphere
     real(real_kind) :: dinv(2, 2, np, np) ! map vector field on the sphere to covariant v on cube
! Convert vector fields from spherical to rectangular components
! The transpose of this operation is its pseudoinverse.
     real(real_kind) :: vec_sphere2cart(np, np, 3, 2)
! Mass matrix terms for an element on a cube face
     real(real_kind) :: mp(np, np) ! mass matrix on velocity and pressure grid
     real(real_kind) :: rmp(np, np) ! inverse mass matrix on velocity and pressure grid
! Mass matrix terms for an element on the sphere
! This mass matrix is used when solving the equations in weak form
! with the natural (surface area of the sphere) inner product
     real(real_kind) :: spheremp(np, np) ! mass matrix on velocity and pressure grid
     real(real_kind) :: rspheremp(np, np) ! inverse mass matrix on velocity and pressure grid
     integer(long_kind) :: gdofp(np, np) ! global degree of freedom(p-grid)
! Coreolis term
     real(real_kind) :: fcor(np, np) ! coreolis term !sjchoi 2*omega*sin(phi)
#if defined(CORE_SW)
     real(real_kind) :: ecor(np, np) ! coreolis term !sjchoi 2*omega*cos(phi)
#endif
! Solver weights (used only for non-staggered grid
     real(real_kind) :: solver_wts(np, np)
     type(index_t) :: idxp
     type(index_t), pointer :: idxv
     integer :: facenum
! force element_t to be a multiple of 8 bytes.
! on BGP, code will crash (signal 7, or signal 15) if 8 byte alignment is off
! check core file for:
! core.63:Generated by interrupt..(Alignment Exception DEAR=0xa1ef671c ESR=0x01800000 CCR0=0x4800a002)
     integer :: dummy
   end type element_t
!
   type, public :: eroot_t
     sequence
     type(element_t), pointer :: first
   end type eroot_t
   integer, public :: numedges
   type(eroot_t), public :: eroot
   public :: element_coordinates
   public :: element_var_coordinates
   public :: element_var_coordinates3d
   public :: lladdedge, llfindedge, llinsertedge
   public :: llsetedgecount, llgetedgecount
   public :: llfree
   public :: getcolumnidp, getcolumnidv
!
   contains
! =======================================
!  GetColumnIdP:
!
!  Gives a unique identifier for a Physics
!  column on the P-grid
! =======================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function getcolumnidp(elem, i, j) result(col_id)
!
   type(element_t), intent(in   ) :: elem
   integer, intent(in   ) :: i, j
   integer :: col_id
!
   col_id = elem%gdofp(i,j)
!
   end function getcolumnidp
! =======================================
!  GetColumnIdV:
!
!  Gives a unique identifier for a Physics
!  column on the V-grid
! =======================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function getcolumnidv(elem, i, j) result(col_id)
!
   type(element_t), intent(in   ) :: elem
   integer, intent(in   ) :: i, j
   integer :: col_id
!
   col_id = elem%gdofp(i,j)
!
   end function getcolumnidv
! =======================================
! element_coordinates:
!
! Initialize 2D rectilinear element
! colocation points
!
! =======================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function element_coordinates(start, end, points) result(cart)
!
   type(cartesian2d_t), intent(in   ) :: start
   type(cartesian2d_t), intent(in   ) :: end
   real(longdouble_kind), intent(in   ) :: points(:)
   type(cartesian2d_t) :: cart(size(points), size(points))
   type(cartesian2d_t) :: length, centroid
   real(longdouble_kind) :: y
   integer i, j
!
   length%x = 0.50d0*(end%x-start%x)
   length%y = 0.50d0*(end%y-start%y)
   centroid%x = 0.50d0*(end%x+start%x)
   centroid%y = 0.50d0*(end%y+start%y)
   do j = 1,size(points)
     y = centroid%y+length%y*points(j)
     do i = 1,size(points)
       cart(i,j)%x = centroid%x+length%x*points(i)
       cart(i,j)%y = y
     enddo
   enddo
!
   end function element_coordinates
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function element_var_coordinates(c, points) result(cart)
!
   type(cartesian2d_t), intent(in   ) :: c(4)
   real(longdouble_kind), intent(in   ) :: points(:)
   type(cartesian2d_t) :: cart(size(points), size(points))
   real(longdouble_kind) :: p(size(points))
   real(longdouble_kind) :: q(size(points))
   integer i, j
!
   p(:) =(1.0d0-points(:))/2.0d0
   q(:) =(1.0d0+points(:))/2.0d0
   do j = 1,size(points)
     do i = 1, size(points)
       cart(i,j)%x = p(i)*p(j)*c(1)%x &
                                             +q(i)*p(j)*c(2)%x &
                                             +q(i)*q(j)*c(3)%x &
                                             +p(i)*q(j)*c(4)%x
       cart(i,j)%y = p(i)*p(j)*c(1)%y &
                                             +q(i)*p(j)*c(2)%y &
                                             +q(i)*q(j)*c(3)%y &
                                             +p(i)*q(j)*c(4)%y
     enddo
   enddo
!
   end function element_var_coordinates
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function element_var_coordinates3d(c, points) result(cart)
!
   type(cartesian3d_t), intent(in   ) :: c(4)
   real(longdouble_kind), intent(in   ) :: points(:)
   type(cartesian3d_t) :: cart(size(points), size(points))
   real(longdouble_kind) :: p(size(points))
   real(longdouble_kind) :: q(size(points)), r
   integer i, j
!
   p(:) =(1.0d0-points(:))/2.0d0
   q(:) =(1.0d0+points(:))/2.0d0
   do j = 1,size(points)
     do i = 1, size(points)
       cart(i,j)%x = p(i)*p(j)*c(1)%x &
                                             +q(i)*p(j)*c(2)%x &
                                             +q(i)*q(j)*c(3)%x &
                                             +p(i)*q(j)*c(4)%x
       cart(i,j)%y = p(i)*p(j)*c(1)%y &
                                             +q(i)*p(j)*c(2)%y &
                                             +q(i)*q(j)*c(3)%y &
                                             +p(i)*q(j)*c(4)%y
       cart(i,j)%z = p(i)*p(j)*c(1)%z &
                                             +q(i)*p(j)*c(2)%z &
                                             +q(i)*q(j)*c(3)%z &
                                             +p(i)*q(j)*c(4)%z
! project back to sphere:
       r = distance(cart(i,j))
       cart(i,j)%x = cart(i,j)%x/r
       cart(i,j)%y = cart(i,j)%y/r
       cart(i,j)%z = cart(i,j)%z/r
     enddo
   enddo
!
   end function element_var_coordinates3d
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine llsetedgecount(value)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: value
!
   numedges = value
!
   end subroutine llsetedgecount
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine llgetedgecount(value)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(  out) :: value
!
   value = numedges
!
   end subroutine llgetedgecount
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   recursive subroutine copy_node(node2, node1)
!
   type(element_t), intent(  out) :: node2
   type(element_t), intent(in   ) :: node1
!
   node2%localid = node1%localid
   node2%globalid = node1%globalid
   node2%prev = node1%prev
   node2%next = node1%next
!
   end subroutine copy_node
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine llfree(list)
!-------------------------------------------------------------------------------
   implicit none
!
   type(eroot_t) :: list
   type(element_t), pointer :: temp_node
   integer :: nlist, i
!
   temp_node=>list%first
! Find the end of the list
   do while(associated(temp_node%next))
     temp_node=>temp_node%next
   enddo
   temp_node=>temp_node%prev
!Now step back and deallocate all entries
   do while(associated(temp_node))
     deallocate(temp_node%next)
     temp_node=>temp_node%prev
   enddo
!
   end subroutine llfree
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine llinsertedge(edgelist, gid, lid)
!
   type(eroot_t), intent(inout) :: edgelist
   integer, intent(in   ) :: gid
   integer, intent(inout) :: lid
   logical :: found
!
   call llfindedge(edgelist,gid,found)
   if (.not.found) then
     call lladdedge(edgelist,gid,lid)
   endif
!
   end subroutine llinsertedge
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine llfindedge(edge, gid, found)
!
   type(eroot_t), intent(in   ) :: edge
   integer, intent(in   ) :: gid
   logical, intent(  out) :: found
   type(element_t), pointer :: temp_node
!
   found = .false.
   temp_node=>edge%first
   do while(associated(temp_node).and.(.not.found))
     if (gid.eq.temp_node%globalid) then
       found = .true.
     else
       temp_node=>temp_node%next
     endif
   enddo
!
   end subroutine llfindedge
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine lladdedge(edgelist, gid, lid)
!
   type(eroot_t), intent(inout) :: edgelist
   integer, intent(in   ) :: gid
   integer, intent(  out) :: lid
   type(element_t), pointer :: temp_node
   type(element_t), pointer :: new_node
   type(element_t), pointer :: parent
!
   temp_node=>edgelist%first
   parent=>edgelist%first
   do while(associated(temp_node))
     parent=>temp_node
     temp_node=>parent%next
   enddo
   allocate(new_node)
   numedges = numedges+1
   new_node%globalid = gid
   new_node%localid = numedges
   nullify(new_node%next)
   new_node%prev=>parent
   if (associated(edgelist%first)) then
     parent%next=>new_node
   else
     edgelist%first=>new_node
   endif
   lid = numedges
!
   end subroutine lladdedge
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module element
!-------------------------------------------------------------------------------
