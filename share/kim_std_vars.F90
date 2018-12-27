!-------------------------------------------------------------------------------
   module kim_std_vars
!-------------------------------------------------------------------------------
!
!  abstract :  bi-linear interpolation module
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-15  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only : i4, l4, r4, r8
   use statistics, only : get_lnorm, get_lerror
   use netcdf
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   integer(i4), parameter :: len_filename = 512
   integer(i4), parameter :: len_varname = 32
   integer(i4), parameter :: max_nfiles = 512
   ! standard outputs for kim
   integer(i4), parameter :: std_nvars = 111 !81 !79 !59
   character(len=len_varname), dimension(std_nvars) :: std_varnames
   integer(i4), dimension(std_nvars) :: std_varndims
   logical(l4), dimension(std_nvars) :: std_usedb
   ! value and coordinate
   type values_t
     real(r8) :: value, mu, minmax
     real(r8) :: width, nwidth
   end type values_t
   type coordinate_t
     integer(i4) :: iup
     integer(i4) :: iface, ie, je
     real(r8) :: lon, lat
   end type coordinate_t
   ! analysis
   integer(i4), parameter :: nanals = 7
   integer(i4), parameter :: nstats = 4 !(min, max, mu, sigma)
   type analysis_t
     integer(i4), dimension(:), allocatable :: npts ! nanals
     integer(i4) :: nmins, nmaxs
     type(values_t), dimension(:), allocatable :: value_min !
     type(values_t), dimension(:), allocatable :: value_max !
     type(coordinate_t), dimension(:), allocatable :: coord_min !
     type(coordinate_t), dimension(:), allocatable :: coord_max !
   end type analysis_t
   ! kim variables
   type kim_var_t ! nvars
     character(len=len_varname) :: varname
     integer(i4) :: varndim
     logical(l4) :: isvar
     real(r8), dimension(:,:), allocatable :: level1 ! nstats x 0 : nlevs ! normal
     real(r8), dimension(:,:,:), allocatable :: level2 ! nstats x nhoriz x nlevs ! for db
     type(analysis_t) :: anal
   end type kim_var_t
   ! kim files
   type kim_file_t
     character(len=len_filename) :: filename
     integer(i4) :: nhoriz, nlevs ! resolutions
     type(coordinate_t), dimension(:), allocatable :: coord
     integer(i4) :: nvars ! nhoriz
     type(kim_var_t), dimension(:), allocatable :: vars ! nvars
   end type kim_file_t
!
   interface putvarlevel1
     module procedure putvarlevel1_1d
     module procedure putvarlevel1_2d
   end interface putvarlevel1
   interface putvarlevels_nfiles
     module procedure putvarlevels_nfiles_1d
     module procedure putvarlevels_nfiles_2d
   end interface putvarlevels_nfiles
   interface analysisvars
     module procedure analysisvars_1d
     module procedure analysisvars_2d
   end interface analysisvars
!
   public :: len_filename, len_varname, max_nfiles
   public :: kim_var_t, kim_file_t, nc_check
   public :: inikimvars, setreskimfile, finkimvars
   public :: inikimfile, finkimfile, inikimfiles, finkimfiles, writekimfile, readkimfile
   public :: putvarlevel1, putvarlevels_nfiles, analysisvars, printfilestatus, printfileanalysis, printfileminmax
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine inistdout()
!-------------------------------------------------------------------------------
   implicit none
!
   ! set varnames
   std_varnames(001) = 'u';         std_varndims(001) = 2;   std_useDB(001) = .TRUE.
   std_varnames(002) = 'v';         std_varndims(002) = 2;   std_useDB(002) = .TRUE.
   std_varnames(003) = 'T';         std_varndims(003) = 2;   std_useDB(003) = .TRUE.
   std_varnames(004) = 'p';         std_varndims(004) = 2;   std_useDB(004) = .TRUE.
   std_varnames(005) = 'pint';      std_varndims(005) = 3;   std_useDB(005) = .TRUE.
   std_varnames(006) = 'ps';        std_varndims(006) = 1;   std_useDB(006) = .TRUE.
   std_varnames(007) = 'psl';       std_varndims(007) = 1;   std_useDB(007) = .TRUE.
   std_varnames(008) = 'omega';     std_varndims(008) = 2;   std_useDB(008) = .TRUE.
   std_varnames(009) = 'q';         std_varndims(009) = 2;   std_useDB(009) = .TRUE.
   std_varnames(010) = 'qc';        std_varndims(010) = 2;   std_useDB(010) = .FALSE.
   std_varnames(011) = 'qr';        std_varndims(011) = 2;   std_useDB(011) = .FALSE.
   std_varnames(012) = 'qs';        std_varndims(012) = 2;   std_useDB(012) = .FALSE.
   std_varnames(013) = 'qi';        std_varndims(013) = 2;   std_useDB(013) = .FALSE.
   std_varnames(014) = 'sst';       std_varndims(014) = 1;   std_useDB(014) = .TRUE.
   std_varnames(015) = 'topo';      std_varndims(015) = 1;   std_useDB(015) = .TRUE.
   std_varnames(016) = 'tsfc';      std_varndims(016) = 1;   std_useDB(016) = .TRUE.
   std_varnames(017) = 'seaice';    std_varndims(017) = 1;   std_useDB(017) = .TRUE.
   std_varnames(018) = 'landfrac';  std_varndims(018) = 1;   std_useDB(018) = .TRUE.
   std_varnames(019) = 'tg3';       std_varndims(019) = 1;   std_useDB(019) = .TRUE.
   std_varnames(020) = 'znt';       std_varndims(020) = 1;   std_useDB(020) = .TRUE.
   std_varnames(021) = 'weasd';     std_varndims(021) = 1;   std_useDB(021) = .TRUE.
   std_varnames(022) = 'ncrain';    std_varndims(022) = 1;   std_useDB(022) = .FALSE.
   std_varnames(023) = 'crain';     std_varndims(023) = 1;   std_useDB(023) = .FALSE.
   std_varnames(024) = 'soilm1';    std_varndims(024) = 1;   std_useDB(024) = .TRUE.
   std_varnames(025) = 'soilm2';    std_varndims(025) = 1;   std_useDB(025) = .TRUE.
   std_varnames(026) = 'soilm3';    std_varndims(026) = 1;   std_useDB(026) = .TRUE.
   std_varnames(027) = 'soilm4';    std_varndims(027) = 1;   std_useDB(027) = .TRUE.
   std_varnames(028) = 'soilt1';    std_varndims(028) = 1;   std_useDB(028) = .TRUE.
   std_varnames(029) = 'soilt2';    std_varndims(029) = 1;   std_useDB(029) = .TRUE.
   std_varnames(030) = 'soilt3';    std_varndims(030) = 1;   std_useDB(030) = .TRUE.
   std_varnames(031) = 'soilt4';    std_varndims(031) = 1;   std_useDB(031) = .TRUE.
   std_varnames(032) = 'canopy';    std_varndims(032) = 1;   std_useDB(032) = .TRUE.
   std_varnames(033) = 'rh';        std_varndims(033) = 2;   std_useDB(033) = .TRUE.
   std_varnames(034) = 'hgt';       std_varndims(034) = 2;   std_useDB(034) = .TRUE.
   std_varnames(035) = 'u10m';      std_varndims(035) = 1;   std_useDB(035) = .FALSE.
   std_varnames(036) = 'v10m';      std_varndims(036) = 1;   std_useDB(036) = .FALSE.
   std_varnames(037) = 't2m';       std_varndims(037) = 1;   std_useDB(037) = .FALSE.
   std_varnames(038) = 'q2m';       std_varndims(038) = 1;   std_useDB(038) = .FALSE.
   std_varnames(039) = 'rh2m';      std_varndims(039) = 1;   std_useDB(039) = .FALSE.
   std_varnames(040) = 'Tmsl';      std_varndims(040) = 1;   std_useDB(040) = .TRUE.
   std_varnames(041) = 'slmsk';     std_varndims(041) = 1;   std_useDB(041) = .TRUE.
   std_varnames(042) = 'tmax';      std_varndims(042) = 1;   std_useDB(042) = .FALSE.
   std_varnames(043) = 'tmin';      std_varndims(043) = 1;   std_useDB(043) = .FALSE.
   std_varnames(044) = 'tcld';      std_varndims(044) = 1;   std_useDB(044) = .FALSE.
   std_varnames(045) = 'hcld';      std_varndims(045) = 1;   std_useDB(045) = .FALSE.
   std_varnames(046) = 'mcld';      std_varndims(046) = 1;   std_useDB(046) = .FALSE.
   std_varnames(047) = 'lcld';      std_varndims(047) = 1;   std_useDB(047) = .FALSE.
   std_varnames(048) = 'shtfl';     std_varndims(048) = 1;   std_useDB(048) = .FALSE.
   std_varnames(049) = 'lhtfl';     std_varndims(049) = 1;   std_useDB(049) = .FALSE.
   std_varnames(050) = 'dlwrsfc';   std_varndims(050) = 1;   std_useDB(050) = .FALSE.
   std_varnames(051) = 'ulwrsfc';   std_varndims(051) = 1;   std_useDB(051) = .FALSE.
   std_varnames(052) = 'ulwrtoa';   std_varndims(052) = 1;   std_useDB(052) = .FALSE.
   std_varnames(053) = 'uswrtoa';   std_varndims(053) = 1;   std_useDB(053) = .FALSE.
   std_varnames(054) = 'uswrsfc';   std_varndims(054) = 1;   std_useDB(054) = .FALSE.
   std_varnames(055) = 'dswrsfc';   std_varndims(055) = 1;   std_useDB(055) = .FALSE.
   std_varnames(056) = 'dswrtoa';   std_varndims(056) = 1;   std_useDB(056) = .FALSE.
   std_varnames(057) = 'xmu';       std_varndims(057) = 1;   std_useDB(057) = .FALSE.
   std_varnames(058) = 'lats';      std_varndims(058) = 1;   std_useDB(058) = .FALSE.
   std_varnames(059) = 'lons';      std_varndims(059) = 1;   std_useDB(059) = .FALSE.
   std_varnames(060) = 'cld';       std_varndims(060) = 2;   std_useDB(060) = .FALSE.
   std_varnames(061) = 'oz';        std_varndims(061) = 2;   std_useDB(061) = .FALSE.
   std_varnames(062) = 'frontgf';   std_varndims(062) = 2;   std_useDB(062) = .FALSE.
! dynamics
   std_varnames(063) = 'psfc';      std_varndims(063) = 1;   std_useDB(063) = .FALSE.
   std_varnames(064) = 'gsfc';      std_varndims(064) = 1;   std_useDB(064) = .FALSE.
   std_varnames(065) = 'pdsb';      std_varndims(065) = 1;   std_useDB(065) = .FALSE.
   std_varnames(066) = 'pdsp';      std_varndims(066) = 1;   std_useDB(066) = .FALSE.
   std_varnames(067) = 'pds' ;      std_varndims(067) = 1;   std_useDB(067) = .FALSE.
   std_varnames(068) = 'ht'  ;      std_varndims(068) = 1;   std_useDB(068) = .FALSE.
   std_varnames(069) = 'dht_u';     std_varndims(069) = 1;   std_useDB(069) = .FALSE.
   std_varnames(070) = 'dht_v';     std_varndims(070) = 1;   std_useDB(070) = .FALSE.
   std_varnames(071) = 'u_m';       std_varndims(071) = 2;   std_useDB(071) = .FALSE.
   std_varnames(072) = 'v_m';       std_varndims(072) = 2;   std_useDB(072) = .FALSE.
   std_varnames(073) = 'theta_m';   std_varndims(073) = 2;   std_useDB(073) = .FALSE.
   std_varnames(074) = 'thetap_m';  std_varndims(074) = 2;   std_useDB(074) = .FALSE.
   std_varnames(075) = 'omega_m';   std_varndims(075) = 2;   std_useDB(075) = .FALSE.
   std_varnames(076) = 'rho_m';     std_varndims(076) = 2;   std_useDB(076) = .FALSE.
   std_varnames(077) = 'alb_m';     std_varndims(077) = 2;   std_useDB(077) = .FALSE.
   std_varnames(078) = 'alp_m';     std_varndims(078) = 2;   std_useDB(078) = .FALSE.
   std_varnames(079) = 'al_m';      std_varndims(079) = 2;   std_useDB(079) = .FALSE.
   std_varnames(080) = 'mub_m';     std_varndims(080) = 2;   std_useDB(080) = .FALSE.
   std_varnames(081) = 'mup_m';     std_varndims(081) = 2;   std_useDB(081) = .FALSE.
   std_varnames(082) = 'mu_m';      std_varndims(082) = 2;   std_useDB(082) = .FALSE.
   std_varnames(083) = 'pb_m';      std_varndims(083) = 2;   std_useDB(083) = .FALSE.
   std_varnames(084) = 'pp_m';      std_varndims(084) = 2;   std_useDB(084) = .FALSE.
   std_varnames(085) = 'p_m';       std_varndims(085) = 2;   std_useDB(085) = .FALSE.
   std_varnames(086) = 'p_hyd_m';   std_varndims(086) = 2;   std_useDB(086) = .FALSE.
   std_varnames(087) = 'ru_tend';   std_varndims(087) = 2;   std_useDB(087) = .FALSE.
   std_varnames(088) = 'rv_tend';   std_varndims(088) = 2;   std_useDB(088) = .FALSE.
   std_varnames(089) = 't_tend';    std_varndims(089) = 2;   std_useDB(089) = .FALSE.
   std_varnames(090) = 't_init';    std_varndims(090) = 2;   std_useDB(090) = .FALSE.
   std_varnames(091) = 't_m';       std_varndims(091) = 2;   std_useDB(091) = .FALSE.
   std_varnames(092) = 'h_m';       std_varndims(092) = 2;   std_useDB(092) = .FALSE.
   std_varnames(093) = 'divv';      std_varndims(093) = 3;   std_useDB(093) = .FALSE.
   std_varnames(094) = 'w_i';       std_varndims(094) = 3;   std_useDB(094) = .FALSE.
   std_varnames(095) = 'mub_i';     std_varndims(095) = 3;   std_useDB(095) = .FALSE.
   std_varnames(096) = 'mup_i';     std_varndims(096) = 3;   std_useDB(096) = .FALSE.
   std_varnames(097) = 'mu_i';      std_varndims(097) = 3;   std_useDB(097) = .FALSE.
   std_varnames(098) = 'p_hyd_i';   std_varndims(098) = 3;   std_useDB(098) = .FALSE.
   std_varnames(099) = 'ph0_i';     std_varndims(099) = 3;   std_useDB(099) = .FALSE.
   std_varnames(100) = 'phb_i';     std_varndims(100) = 3;   std_useDB(100) = .FALSE.
   std_varnames(101) = 'php_i';     std_varndims(101) = 3;   std_useDB(101) = .FALSE.
   std_varnames(102) = 'p8w_i';     std_varndims(102) = 3;   std_useDB(102) = .FALSE.
   std_varnames(103) = 'rw_tend';   std_varndims(103) = 3;   std_useDB(103) = .FALSE.
   std_varnames(104) = 'h_i';       std_varndims(104) = 3;   std_useDB(104) = .FALSE.
   std_varnames(105) = 'qv_dyn';    std_varndims(105) = 2;   std_useDB(105) = .FALSE.
   std_varnames(106) = 'qc_dyn';    std_varndims(106) = 2;   std_useDB(106) = .FALSE.
   std_varnames(107) = 'qr_dyn';    std_varndims(107) = 2;   std_useDB(107) = .FALSE.
   std_varnames(108) = 'qi_dyn';    std_varndims(108) = 2;   std_useDB(108) = .FALSE.
   std_varnames(109) = 'qs_dyn';    std_varndims(109) = 2;   std_useDB(109) = .FALSE.
   std_varnames(110) = 'oz_dyn';    std_varndims(110) = 2;   std_useDB(110) = .FALSE.
   std_varnames(111) = 'cld_dyn';   std_varndims(111) = 2;   std_useDB(111) = .FALSE.
#if 0
   std_varnames(110) = 'm1_tend';   std_varndims(110) = 2;   std_useDB(110) = .FALSE.
   std_varnames(111) = 'm2_tend';   std_varndims(111) = 2;   std_useDB(111) = .FALSE.
   std_varnames(112) = 'm3_tend';   std_varndims(112) = 2;   std_useDB(112) = .FALSE.
   std_varnames(113) = 'm4_tend';   std_varndims(113) = 2;   std_useDB(113) = .FALSE.
   std_varnames(114) = 'm5_tend';   std_varndims(114) = 2;   std_useDB(114) = .FALSE.
   std_varnames(115) = 'm1_tendf';  std_varndims(115) = 2;   std_useDB(115) = .FALSE.
   std_varnames(116) = 'm2_tendf';  std_varndims(116) = 2;   std_useDB(116) = .FALSE.
   std_varnames(117) = 'm3_tendf';  std_varndims(117) = 2;   std_useDB(117) = .FALSE.
   std_varnames(118) = 'm4_tendf';  std_varndims(118) = 2;   std_useDB(118) = .FALSE.
   std_varnames(119) = 'm5_tendf';  std_varndims(119) = 2;   std_useDB(119) = .FALSE.
#endif
!
   return
   end subroutine inistdout
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
! inikimvars (vars_t, nvars, varnames)  <= use 'varnames' vars
! inikimvars (vars_t, nvars)            <= use db output vars
   subroutine inikimvars(vars_t, nvars, varnames)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_var_t), dimension(:), allocatable, intent(inout) :: vars_t
   integer(i4),                                intent(  out) :: nvars
   character(len=*), dimension(:), optional,   intent(in   ) :: varnames
! local variables
   integer(i4) :: ivar, std_ivar
!
   ! standard variables
   call inistdout()
!
   if (present(varnames)) then
!
     nvars = size(varnames)
     allocate(vars_t(nvars))
     !
     do ivar = 1,nvars
     do std_ivar = 1,std_nvars
       if (trim(varnames(ivar)).eq.trim(std_varnames(std_ivar))) then
         vars_t(ivar)%varname = std_varnames(std_ivar)
         vars_t(ivar)%varndim=std_varndims(std_ivar)
         vars_t(ivar)%isvar = .false.
       endif
     enddo
     enddo
!
   else
!
     nvars = 0
     do std_ivar = 1,std_nvars
       if (std_usedb(std_ivar)) then
         nvars = nvars+1
       endif
     enddo
     allocate(vars_t(nvars))
     ivar = 0
     do std_ivar = 1,std_nvars
       if (std_usedb(std_ivar)) then
         ivar = ivar+1
         vars_t(ivar)%varname = std_varnames(std_ivar)
         vars_t(ivar)%varndim = std_varndims(std_ivar)
         vars_t(ivar)%isvar   = .true.
       endif
     enddo
!
   endif
!
   return
   end subroutine inikimvars
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine setreskimvars(vars_t, nhoriz, nlevs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_var_t), dimension(:), allocatable, intent(inout) :: vars_t
   integer(i4),                                intent(in   ) :: nhoriz, nlevs
! local variables
   integer(i4) :: nvars, ivar
!
  ! standard variables
   nvars = size(vars_t)
!
   do ivar = 1,nvars
!
     if (allocated(vars_t(ivar)%level1)) deallocate(vars_t(ivar)%level1)
     if (allocated(vars_t(ivar)%level2)) deallocate(vars_t(ivar)%level2)
     if (vars_t(ivar)%varndim==1) then
       allocate(vars_t(ivar)%level1(nstats,0:1))
       allocate(vars_t(ivar)%level2(nstats,nhoriz,1))
     elseif (vars_t(ivar)%varndim==2) then
       allocate(vars_t(ivar)%level1(nstats,0:nlevs))
       allocate(vars_t(ivar)%level2(nstats,nhoriz,nlevs))
     elseif (vars_t(ivar)%varndim==3) then
       allocate(vars_t(ivar)%level1(nstats,0:nlevs+1))
       allocate(vars_t(ivar)%level2(nstats,nhoriz,nlevs+1))
     else
       print*,'error(varndim) in setreskimvars'
       stop
     endif
   enddo
!
   return
   end subroutine setreskimvars
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine setreskimfile(file_t, nhoriz, nlevs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_file_t), intent(inout) :: file_t
   integer(i4),      intent(in   ) :: nhoriz, nlevs
!
   ! standard variables
   file_t%nhoriz = nhoriz
   file_t%nlevs = nlevs
   if (allocated(file_t%coord)) deallocate(file_t%coord)
   allocate(file_t%coord(nhoriz))
!
   call setreskimvars(file_t%vars,nhoriz,nlevs)
!
   return
   end subroutine setreskimfile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine finkimvars(vars_t)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_var_t), dimension(:), allocatable, intent(inout) :: vars_t
! local variables
   integer(i4) :: nvars, ivar
!
   nvars = size(vars_t)
!
   do ivar = 1,nvars
     if (allocated(vars_t(ivar)%level1)) deallocate(vars_t(ivar)%level1)
     if (allocated(vars_t(ivar)%level2)) deallocate(vars_t(ivar)%level2)
     call finanalysis(vars_t(ivar)%anal)
   enddo
   if (allocated(vars_t)) deallocate(vars_t)
!
   return
   end subroutine finkimvars
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine inianalysis(anal_t, nmins, nmaxs)
!-------------------------------------------------------------------------------
   implicit none
!
   type(analysis_t), intent(inout) :: anal_t
   integer(i4),      intent(in   ) :: nmins, nmaxs
!
   if (allocated(anal_t%npts))      deallocate(anal_t%npts)
   if (allocated(anal_t%value_min)) deallocate(anal_t%value_min)
   if (allocated(anal_t%value_max)) deallocate(anal_t%value_max)
   if (allocated(anal_t%coord_min)) deallocate(anal_t%coord_min)
   if (allocated(anal_t%coord_max)) deallocate(anal_t%coord_max)
   anal_t%nmins = nmins
   anal_t%nmaxs = nmaxs
   allocate(anal_t%npts(nanals))
   allocate(anal_t%value_min(nmins))
   allocate(anal_t%value_max(nmaxs))
   allocate(anal_t%coord_min(nmins))
   allocate(anal_t%coord_max(nmaxs))
   anal_t%coord_min(:)%iface =-1
   anal_t%coord_min(:)%ie    =-1
   anal_t%coord_min(:)%je    =-1
   anal_t%coord_max(:)%iface =-1
   anal_t%coord_max(:)%ie    =-1
   anal_t%coord_max(:)%je    =-1
!
   return
   end subroutine inianalysis
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine finanalysis(anal_t)
!-------------------------------------------------------------------------------
   implicit none
!
   type(analysis_t), intent(inout) :: anal_t
!
   if (allocated(anal_t%npts)) deallocate(anal_t%npts)
   if (allocated(anal_t%value_min)) deallocate(anal_t%value_min)
   if (allocated(anal_t%value_max)) deallocate(anal_t%value_max)
   if (allocated(anal_t%coord_min)) deallocate(anal_t%coord_min)
   if (allocated(anal_t%coord_max)) deallocate(anal_t%coord_max)
!
   return
   end subroutine finanalysis
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
! inikimfile(file_t, filename)                                      <- use standard output vars
! inikimfile(file_t, filename, nvars(<=0), varnames)                <- use standard output vars
! inikimfile(file_t, filename, nvars( >0), varnames)                <- use specific output vars
! inikimfile(file_t, filename, nvars( >0), varnames, usedb=.true.)  <- use db output vars
   subroutine inikimfile(file_t, filename, nvars, varnames, usedb)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_file_t),                         intent(inout) :: file_t
   character(len=*),                         intent(in   ) :: filename
   integer(i4), optional,                    intent(in   ) :: nvars
   character(len=*), dimension(:), optional, intent(in   ) :: varnames
   logical(l4), optional,                    intent(in   ) :: usedb
! local variables
   integer(i4) :: l_nvars
   logical(l4) :: l_usedb
!
   if (present(usedb)) then
     l_usedb = usedb
   else
     l_usedb = .false.
   endif
!
   file_t%filename = filename
   if (present(nvars).and.present(varnames)) then
     if (nvars.le.0) then
       call inikimvars(file_t%vars,l_nvars,std_varnames)
     else
       call inikimvars(file_t%vars,l_nvars,varnames)
     endif
   else
     if (l_usedb) then
       call inikimvars(file_t%vars,l_nvars) ! for read/write db
     else
       call inikimvars(file_t%vars,l_nvars,std_varnames)
     endif
   endif
   file_t%nvars = l_nvars
!
   return
   end subroutine inikimfile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine finkimfile(file_t)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_file_t), intent(inout) :: file_t
!
   if (allocated(file_t%coord)) deallocate(file_t%coord)
   call finkimvars(file_t%vars)
!
   return
   end subroutine finkimfile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine inikimfiles(files_t, filenames, nvars, varnames, usedb)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_file_t), dimension(:), allocatable, intent(inout) :: files_t
   character(len=*), dimension(:),              intent(in   ) :: filenames
   integer(i4), optional,                       intent(in   ) :: nvars
   character(len=*), dimension(:), optional,    intent(in   ) :: varnames
   logical(l4), optional,                       intent(in   ) :: usedb
! local variables
   integer(i4) :: nfiles, ifile, l_nvars
   logical(l4) :: l_usedb
!
   nfiles = size(filenames)
   allocate(files_t(nfiles))
!
   if (present(usedb)) then
     l_usedb = usedb
   else
     l_usedb = .false.
   endif
!
   do ifile = 1,nfiles
     call inikimfile(files_t(ifile),filenames(ifile),nvars,varnames,l_usedb)
   enddo
!
   return
   end subroutine inikimfiles
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine finkimfiles(files_t)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_file_t), dimension(:), allocatable, intent(inout) :: files_t
! local variables
   integer(i4) :: nfiles, ifile
!
   nfiles = size(files_t)
!
   do ifile = 1, nfiles
     call finkimfile(files_t(ifile))
   enddo
!
   return
   end subroutine finkimfiles
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine writekimfile(file_t)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_file_t), intent(in   ) :: file_t
! local variables
   integer(i4) :: ierr, ivar, nvars, nhoriz, nlevs, idx_dim
   integer(i4) :: wmod, fileid
   integer(i4) :: did_nhoriz, did_1lev, did_2lev, did_nlevs, did_nlevsp, did_nlevsp2, did_nvars, did_nstats
   integer(i4), dimension(2, 3) :: dimids_lev1
   integer(i4), dimension(3, 3) :: dimids_lev2
   integer(i4) :: varid_lons, varid_lats
   integer(i4), dimension(:), allocatable :: varids_lev1, varids_lev2
!
   nvars = file_t%nvars
   nhoriz = file_t%nhoriz
   nlevs = file_t%nlevs
!
   ! file create
   wmod = ior(nf90_clobber,nf90_64bit_offset)
   ierr = nf90_create(trim(file_t%filename),wmod,fileid)
   ! define demensions
   ierr = nf90_def_dim(fileid,'ncol',nhoriz,did_nhoriz)
   ierr = nf90_def_dim(fileid,'nhoriz',nhoriz,did_nhoriz)
   ierr = nf90_def_dim(fileid,'one',1,did_1lev)
   ierr = nf90_def_dim(fileid,'two',2,did_2lev)
   ierr = nf90_def_dim(fileid,'nlevs',nlevs,did_nlevs)
   ierr = nf90_def_dim(fileid,'nlevsp1',nlevs+1,did_nlevsp)
   ierr = nf90_def_dim(fileid,'nlevsp2',nlevs+2,did_nlevsp2)
   ierr = nf90_def_dim(fileid,'nvars',nvars,did_nvars)
   ierr = nf90_def_dim(fileid,'nstats',nstats,did_nstats)
   !
   dimids_lev1(:,1) =(/did_nstats,did_2lev/)
   dimids_lev1(:,2) =(/did_nstats,did_nlevsp/)
   dimids_lev1(:,3) =(/did_nstats,did_nlevsp2/)
   dimids_lev2(:,1) =(/did_nstats,did_nhoriz,did_1lev/)
   dimids_lev2(:,2) =(/did_nstats,did_nhoriz,did_nlevs/)
   dimids_lev2(:,3) =(/did_nstats,did_nhoriz,did_nlevsp/)
   ! define coordinates
   ierr = nf90_def_var(fileid,'lons',nf90_real8,did_nhoriz,varid_lons)
   ierr = nf90_def_var(fileid,'lats',nf90_real8,did_nhoriz,varid_lats)
   ! define varaibles level1 and lev2
   allocate(varids_lev1(nvars),varids_lev2(nvars))
   do ivar = 1,nvars
     idx_dim=file_t%vars(ivar)%varndim
     ierr = nf90_def_var(fileid,trim(file_t%vars(ivar)%varname)//'_lev1',nf90_real8,dimids_lev1(:,idx_dim),varids_lev1(ivar))
   enddo
   do ivar = 1,nvars
     idx_dim=file_t%vars(ivar)%varndim
     ierr = nf90_def_var(fileid,trim(file_t%vars(ivar)%varname)//'_lev2',nf90_real8,dimids_lev2(:,idx_dim),varids_lev2(ivar))
   enddo
   ierr = nf90_enddef(fileid)
   ! put coordinate
   ierr = nf90_put_var(fileid,varid_lons,file_t%coord(:)%lon)
   ierr = nf90_put_var(fileid,varid_lats,file_t%coord(:)%lat)
   ! put varaibles
   do ivar = 1,nvars
     ierr = nf90_put_var(fileid,varids_lev1(ivar),file_t%vars(ivar)%level1(:,:))
     ierr = nf90_put_var(fileid,varids_lev2(ivar),file_t%vars(ivar)%level2(:,:,:))
   enddo
   ! close file
   ierr = nf90_close(fileid)
!
   deallocate(varids_lev1,varids_lev2)
!
   return
   end subroutine writekimfile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine readkimfile(file_t)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_file_t), intent(inout) :: file_t
! local variables
   integer(i4) :: ierr, ivar, nvars, nhoriz, nlevs, idx_dim
   integer(i4) :: fileid
   integer(i4) :: did_nhoriz, did_1lev, did_2lev, did_nlevs, did_nlevsp, did_nlevsp2, did_nvars, did_nstats
   integer(i4), dimension(2, 3) :: dimids_lev1
   integer(i4), dimension(3, 3) :: dimids_lev2
   integer(i4) :: varid
   integer(i4), dimension(:), allocatable :: varids_lev1, varids_lev2
!
   ! file create
   ierr = nf90_open(trim(file_t%filename),nf90_nowrite,fileid)
   call nc_check(ierr)
   ! define demensions
   ierr = nf90_inq_dimid(fileid,'nhoriz',did_nhoriz)
   ierr = nf90_inquire_dimension(fileid,did_nhoriz,len=nhoriz)
   !ierr = nf90_inq_dimid(fileid,'one',    did_1lev)
   !ierr = nf90_inq_dimid(fileid,'two',    did_2lev)
   ierr = nf90_inq_dimid(fileid,'nlevs',did_nlevs)
   ierr = nf90_inquire_dimension(fileid,did_nlevs,len=nlevs)
   !ierr = nf90_inq_dimid(fileid,'nlevsp1',did_nlevsp)
   !ierr = nf90_inq_dimid(fileid,'nlevsp2',did_nlevsp2)
   ierr = nf90_inq_dimid(fileid,'nvars',did_nvars)
   ierr = nf90_inquire_dimension(fileid,did_nvars,len=nvars)
   !ierr = nf90_inq_dimid(fileid,'nstats', did_nstats)
!
   call setreskimfile(file_t,nhoriz,nlevs)
!
   ! get coordinates
   ierr = nf90_inq_varid(fileid,'lons',varid)
   ierr = nf90_get_var(fileid,varid,file_t%coord(:)%lon)
   ierr = nf90_inq_varid(fileid,'lats',varid)
   ierr = nf90_get_var(fileid,varid,file_t%coord(:)%lat)
!
   file_t%nvars = nvars
   ! get varaibles level1 and lev2
   allocate(varids_lev1(nvars),varids_lev2(nvars))
   do ivar = 1,nvars
     ierr = nf90_inq_varid(fileid,trim(file_t%vars(ivar)%varname)//'_lev1',varids_lev1(ivar))
     ierr = nf90_get_var(fileid,varids_lev1(ivar),file_t%vars(ivar)%level1(:,:))
   enddo
   do ivar = 1,nvars
     ierr = nf90_inq_varid(fileid,trim(file_t%vars(ivar)%varname)//'_lev2',varids_lev2(ivar))
     ierr = nf90_get_var(fileid,varids_lev2(ivar),file_t%vars(ivar)%level2(:,:,:))
   enddo
!
   ! close file
   ierr = nf90_close(fileid)
!
   deallocate(varids_lev1,varids_lev2)
!
   return
   end subroutine readkimfile
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function getvarindex(file_t, varname) result(idx)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_file_t), intent(in   ) :: file_t
   character(len=*), intent(in   ) :: varname
   integer(i4) :: idx
! local variables
   integer(i4) :: ivar, nvars
!
   idx =-1
!
   nvars = size(file_t%vars)
!
   do ivar = 1,nvars
     if (trim(varname).eq.(file_t%vars(ivar)%varname)) then
       idx = ivar
       exit
     endif
   enddo
!
   end function getvarindex
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine putvarlevel1_1d(var_t, variable)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_var_t),        intent(inout) :: var_t
   real(r8), dimension(:), intent(in   ) :: variable
! local variables
   integer(i4) :: nhoriz
!
   nhoriz = size(variable(:))
!
   if (var_t%varndim==1) then
     var_t%level1(1,1) = minval(variable(:))
     var_t%level1(2,1) = maxval(variable(:))
     var_t%level1(3,1) = sum(variable(:))/real(nhoriz,8)
     var_t%level1(4,1) = getdeviation(var_t%level1(3,1),nhoriz,variable)
     var_t%level1(1,0) = var_t%level1(1,1)
     var_t%level1(2,0) = var_t%level1(2,1)
     var_t%level1(3,0) = var_t%level1(3,1)
     var_t%level1(4,0) = var_t%level1(4,1)
   else
     print*,'check varndim in putvarlevel1_1d'
     stop
   endif
!
   return
   end subroutine putvarlevel1_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine putvarlevel1_2d(var_t, variable)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_var_t),          intent(inout) :: var_t
   real(r8), dimension(:,:), intent(in   ) :: variable
! local variables
   integer(i4) :: nhoriz, ilev, nlevs
!
   nhoriz = size(variable,dim=1)
   nlevs = size(variable,dim=2)
!
   if (var_t%varndim==2.or.var_t%varndim==3) then
     do ilev = 1,nlevs
       var_t%level1(1,ilev) = minval(variable(:,ilev))
       var_t%level1(2,ilev) = maxval(variable(:,ilev))
       var_t%level1(3,ilev) = sum(variable(:,ilev))/real(nhoriz,8)
       var_t%level1(4,ilev) = getdeviation(var_t%level1(3,ilev),nhoriz,variable(:,ilev))
     enddo
     var_t%level1(1,0) = minval(variable(:,:))
     var_t%level1(2,0) = maxval(variable(:,:))
     var_t%level1(3,0) = sum(variable(:,:))/real(nhoriz*nlevs,8)
     var_t%level1(4,0) = getdeviation(var_t%level1(3,0),nhoriz*nlevs,reshape(variable,shape =(/size(variable)/)))
   else
     print*,'check varndim in putvarlevel1_2d'
     stop
   endif
!
   return
   end subroutine putvarlevel1_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine putvarlevels_nfiles_1d(var_t, variable)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_var_t),          intent(inout) :: var_t
   real(r8), dimension(:,:), intent(in   ) :: variable
! local variables
   integer(i4) :: nhoriz, nfiles, ihoriz
!
   nhoriz = size(variable,dim=1)
   nfiles = size(variable,dim=2)
!
   if (var_t%varndim==1) then
     do ihoriz = 1,nhoriz
       var_t%level2(1,ihoriz,1) = minval(variable(ihoriz,:))
       var_t%level2(2,ihoriz,1) = maxval(variable(ihoriz,:))
       var_t%level2(3,ihoriz,1) = sum(variable(ihoriz,:))/real(nfiles,8)
       var_t%level2(4,ihoriz,1) = getdeviation(var_t%level2(3,ihoriz,1),nfiles,variable(ihoriz,:))
     enddo
     var_t%level1(1,1) = minval(variable(:,:))
     var_t%level1(2,1) = maxval(variable(:,:))
     var_t%level1(3,1) = sum(variable(:,:))/real(nhoriz*nfiles,8)
     var_t%level1(4,1) = getdeviation(var_t%level1(3,1),nhoriz*nfiles,reshape(variable,shape =(/size(variable)/)))
     var_t%level1(1,0) = var_t%level1(1,1)
     var_t%level1(2,0) = var_t%level1(2,1)
     var_t%level1(3,0) = var_t%level1(3,1)
     var_t%level1(4,0) = var_t%level1(4,1)
   else
     print*,'check varndim in putvarlevel1_nfiles_1d'
     stop
   endif
!
   return
   end subroutine putvarlevels_nfiles_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine putvarlevels_nfiles_2d(var_t, variable)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_var_t),            intent(inout) :: var_t
   real(r8), dimension(:,:,:), intent(in   ) :: variable
! local variables
   integer(i4) :: ihoriz, nhoriz, ilev, nlevs, nfiles
!
   nhoriz = size(variable,dim=1)
   nlevs = size(variable,dim=2)
   nfiles = size(variable,dim=3)
!
   if (var_t%varndim==2.or.var_t%varndim==3) then
     do ilev = 1,nlevs
     do ihoriz = 1,nhoriz
       var_t%level2(1,ihoriz,ilev) = minval(variable(ihoriz,ilev,:))
       var_t%level2(2,ihoriz,ilev) = maxval(variable(ihoriz,ilev,:))
       var_t%level2(3,ihoriz,ilev) = sum(variable(ihoriz,ilev,:))/real(nfiles,8)
       var_t%level2(4,ihoriz,ilev) = getdeviation(var_t%level2(3,ihoriz,ilev),nfiles,variable(ihoriz,ilev,:))
     enddo
     enddo
     !
     do ilev = 1,nlevs
       var_t%level1(1,ilev) = minval(variable(:,ilev,:))
       var_t%level1(2,ilev) = maxval(variable(:,ilev,:))
       var_t%level1(3,ilev) = sum(variable(:,ilev,:))/real(nhoriz*nfiles,8)
       var_t%level1(4,ilev) = getdeviation(var_t%level1(3,ilev),nhoriz*nfiles,reshape(variable,shape =(/size(variable(:,ilev,:))/)))
     enddo
     var_t%level1(1,0) = minval(variable(:,:,:))
     var_t%level1(2,0) = maxval(variable(:,:,:))
     var_t%level1(3,0) = sum(variable(:,:,:))/real(nhoriz*nlevs*nfiles,8)
     var_t%level1(4,0) = getdeviation(var_t%level1(3,0),nhoriz*nlevs*nfiles,reshape(variable,shape =(/size(variable)/)))
   else
     print*,'check varndim in putvarlevel1_nfiles_2d'
     stop
   endif
!
   return
   end subroutine putvarlevels_nfiles_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function getdeviation(mu, n, variable) result(res)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),               intent(in   ) :: mu
   integer(l4),            intent(in   ) :: n
   real(r8), dimension(:), intent(in   ) :: variable
   real(r8) :: res
! local variables
   integer(i4) :: i
   real(r8) :: dev
!
   dev = 0.0_r8
!
   if (size(variable).ne.n) then
     print*,'check size in getdeviation...'
     stop
   endif
!
   do i = 1,n
     dev = dev+(variable(i)-mu)**2.0_r8
   enddo
!
   res = sqrt(dev/real(n,8))
!
   end function getdeviation
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine analysisvars_1d(var_t, db_t, variable)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_var_t),        intent(inout) :: var_t
   type(kim_file_t),       intent(in   ) :: db_t
   real(r8), dimension(:), intent(in   ) :: variable
! local variables
   integer(i4) :: nhoriz, ihoriz, idx
   integer(i4) :: n_ltmin, n_gtmax, n_le1sig, n_le2sig, n_le3sig, n_gt3sig
   real(r8) :: mindb, maxdb, mudb, sigdb
!
   nhoriz = size(variable,dim=1)
!
   idx = getvarindex(db_t,trim(var_t%varname))
!
   if (idx.gt.0) then
   if (var_t%varndim==1) then
!
     n_ltmin = 0; n_gtmax = 0; n_le1sig = 0; n_le2sig = 0; n_le3sig = 0; n_gt3sig = 0
     do ihoriz = 1,nhoriz
       mindb = db_t%vars(idx)%level2(1,ihoriz,1)
       maxdb = db_t%vars(idx)%level2(2,ihoriz,1)
       mudb = db_t%vars(idx)%level2(3,ihoriz,1)
       sigdb = db_t%vars(idx)%level2(4,ihoriz,1)
       if (variable(ihoriz).lt.mindb) n_ltmin = n_ltmin+1
       if (variable(ihoriz).gt.maxdb) n_gtmax = n_gtmax+1
       if (variable(ihoriz).ge.mudb-1.0_r8*sigdb.and.variable(ihoriz).le.mudb+1.0_r8*sigdb) n_le1sig = n_le1sig+1
       if (variable(ihoriz).ge.mudb-2.0_r8*sigdb.and.variable(ihoriz).le.mudb+2.0_r8*sigdb) n_le2sig = n_le2sig+1
       if (variable(ihoriz).ge.mudb-3.0_r8*sigdb.and.variable(ihoriz).le.mudb+3.0_r8*sigdb) n_le3sig = n_le3sig+1
       if (variable(ihoriz).lt.mudb-3.0_r8*sigdb.or.variable(ihoriz).gt.mudb+3.0_r8*sigdb) n_gt3sig = n_gt3sig+1
     enddo
     !
     call inianalysis(var_t%anal,n_ltmin,n_gtmax)
     var_t%anal%npts(1) = nhoriz
     var_t%anal%npts(2) = n_ltmin
     var_t%anal%npts(3) = n_gtmax
     var_t%anal%npts(4) = n_le1sig
     var_t%anal%npts(6) = n_le3sig-n_le2sig
     var_t%anal%npts(5) = n_le2sig-n_le1sig
     var_t%anal%npts(7) = n_gt3sig
  
     ! min / max
     n_ltmin = 0; n_gtmax = 0
     do ihoriz = 1,nhoriz
       mindb = db_t%vars(idx)%level2(1,ihoriz,1)
       maxdb = db_t%vars(idx)%level2(2,ihoriz,1)
       mudb = db_t%vars(idx)%level2(3,ihoriz,1)
       sigdb = db_t%vars(idx)%level2(4,ihoriz,1)
       if (variable(ihoriz).lt.mindb) then
         n_ltmin = n_ltmin+1
         var_t%anal%value_min(n_ltmin)%value = variable(ihoriz)
         var_t%anal%value_min(n_ltmin)%mu = mudb
         var_t%anal%value_min(n_ltmin)%minmax = mindb
         var_t%anal%value_min(n_ltmin)%width = abs(variable(ihoriz)-mudb)
         var_t%anal%value_min(n_ltmin)%nwidth = abs(variable(ihoriz)-mudb)/abs(mindb-mudb)
         var_t%anal%coord_min(n_ltmin)%iup = ihoriz
         var_t%anal%coord_min(n_ltmin)%lon = db_t%coord(ihoriz)%lon
         var_t%anal%coord_min(n_ltmin)%lat = db_t%coord(ihoriz)%lat
       endif
       if (variable(ihoriz).gt.maxdb) then
         n_gtmax = n_gtmax+1
         var_t%anal%value_max(n_gtmax)%value = variable(ihoriz)
         var_t%anal%value_max(n_gtmax)%mu = mudb
         var_t%anal%value_max(n_gtmax)%minmax = maxdb
         var_t%anal%value_max(n_gtmax)%width = abs(variable(ihoriz)-mudb)
         var_t%anal%value_max(n_gtmax)%nwidth = abs(variable(ihoriz)-mudb)/abs(maxdb-mudb)
         var_t%anal%coord_max(n_gtmax)%iup = ihoriz
         var_t%anal%coord_max(n_gtmax)%lon = db_t%coord(ihoriz)%lon
         var_t%anal%coord_max(n_gtmax)%lat = db_t%coord(ihoriz)%lat
       endif
     enddo
!
   endif
   endif
!
   return
   end subroutine analysisvars_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine analysisvars_2d(var_t, db_t, variable)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_var_t),          intent(inout) :: var_t
   type(kim_file_t),         intent(in   ) :: db_t
   real(r8), dimension(:,:), intent(in   ) :: variable
! local variables
   integer(i4) :: nhoriz, ihoriz, nlevs, ilev, idx
   integer(i4) :: n_ltmin, n_gtmax, n_le1sig, n_le2sig, n_le3sig, n_gt3sig
   real(r8) :: mindb, maxdb, mudb, sigdb
!
   nhoriz = size(variable,dim=1)
   nlevs = size(variable,dim=2)
!
   idx = getvarindex(db_t,trim(var_t%varname))
!
   if (idx.gt.0) then
   if (var_t%varndim==2.or.var_t%varndim==3) then
!
     n_ltmin = 0; n_gtmax = 0; n_le1sig = 0; n_le2sig = 0; n_le3sig = 0; n_gt3sig = 0
     do ilev = 1,nlevs
     do ihoriz = 1,nhoriz
       mindb = db_t%vars(idx)%level2(1,ihoriz,ilev)
       maxdb = db_t%vars(idx)%level2(2,ihoriz,ilev)
       mudb = db_t%vars(idx)%level2(3,ihoriz,ilev)
       sigdb = db_t%vars(idx)%level2(4,ihoriz,ilev)
       if (variable(ihoriz,ilev).lt.mindb) n_ltmin = n_ltmin+1
       if (variable(ihoriz,ilev).gt.maxdb) n_gtmax = n_gtmax+1
       if (variable(ihoriz,ilev).ge.mudb-1.0_r8*sigdb.and.variable(ihoriz,ilev).le.mudb+1.0_r8*sigdb) n_le1sig = n_le1sig+1
       if (variable(ihoriz,ilev).ge.mudb-2.0_r8*sigdb.and.variable(ihoriz,ilev).le.mudb+2.0_r8*sigdb) n_le2sig = n_le2sig+1
       if (variable(ihoriz,ilev).ge.mudb-3.0_r8*sigdb.and.variable(ihoriz,ilev).le.mudb+3.0_r8*sigdb) n_le3sig = n_le3sig+1
       if (variable(ihoriz,ilev).lt.mudb-3.0_r8*sigdb.or.variable(ihoriz,ilev).gt.mudb+3.0_r8*sigdb) n_gt3sig = n_gt3sig+1
     enddo ! ihoriz
     enddo ! ilev
     call inianalysis(var_t%anal,n_ltmin,n_gtmax)
     var_t%anal%npts(1) = nhoriz*nlevs
     var_t%anal%npts(2) = n_ltmin
     var_t%anal%npts(3) = n_gtmax
     var_t%anal%npts(4) = n_le1sig
     var_t%anal%npts(6) = n_le3sig-n_le2sig
     var_t%anal%npts(5) = n_le2sig-n_le1sig
     var_t%anal%npts(7) = n_gt3sig
     ! min / max
     n_ltmin = 0; n_gtmax = 0
     do ilev = 1,nlevs
     do ihoriz = 1,nhoriz
       mindb = db_t%vars(idx)%level2(1,ihoriz,ilev)
       maxdb = db_t%vars(idx)%level2(2,ihoriz,ilev)
       mudb = db_t%vars(idx)%level2(3,ihoriz,ilev)
       sigdb = db_t%vars(idx)%level2(4,ihoriz,ilev)
       if (variable(ihoriz,ilev).lt.mindb) then
         n_ltmin = n_ltmin+1
         var_t%anal%value_min(n_ltmin)%value = variable(ihoriz,ilev)
         var_t%anal%value_min(n_ltmin)%mu = mudb
         var_t%anal%value_min(n_ltmin)%minmax = mindb
         var_t%anal%value_min(n_ltmin)%width = abs(variable(ihoriz,ilev)-mudb)
         var_t%anal%value_min(n_ltmin)%nwidth = abs(variable(ihoriz,ilev)-mudb)/abs(mindb-mudb)
         var_t%anal%coord_min(n_ltmin)%iup = ihoriz
         var_t%anal%coord_min(n_ltmin)%lon = db_t%coord(ihoriz)%lon
         var_t%anal%coord_min(n_ltmin)%lat = db_t%coord(ihoriz)%lat
       endif
       if (variable(ihoriz,ilev).gt.maxdb) then
         n_gtmax = n_gtmax+1
         var_t%anal%value_max(n_gtmax)%value = variable(ihoriz,ilev)
         var_t%anal%value_max(n_gtmax)%mu = mudb
         var_t%anal%value_max(n_gtmax)%minmax = maxdb
         var_t%anal%value_max(n_gtmax)%width = abs(variable(ihoriz,ilev)-mudb)
         var_t%anal%value_max(n_gtmax)%nwidth = abs(variable(ihoriz,ilev)-mudb)/abs(maxdb-mudb)
         var_t%anal%coord_max(n_gtmax)%iup = ihoriz
         var_t%anal%coord_max(n_gtmax)%lon = db_t%coord(ihoriz)%lon
         var_t%anal%coord_max(n_gtmax)%lat = db_t%coord(ihoriz)%lat
       endif
     enddo !ihoriz
     enddo !ilev
!
   endif ! vardim==2 or 3
   endif ! idx>0
!
   return
   end subroutine analysisvars_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine printfilestatus(files_t, ilev)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_file_t), dimension(:), allocatable, intent(in   ) :: files_t
   integer(i4),                                 intent(in   ) :: ilev
! local variables
   integer(i4) :: nfiles, ifile, nvars, ivar, nlevs, i
   character(len=256) :: fmt_t, fmt_c
   character(len=512) :: sline, bline
   character(len=512) :: string
   character(len=16), dimension(:,:), allocatable :: values
!
   nfiles = size(files_t)
   nvars = files_t(1)%nvars
   nlevs = files_t(1)%nlevs
!
   write(fmt_t, '(a) ') '(a1, a12, a1, a7, a1, a15, a1, a15, a1, a15, a1, a15, a1) '
   write(sline, '(a) ') '--------------------------------------------------------------------------------------'
   write(bline, '(a) ') '======================================================================================'
!
   allocate(values(nstats, 2+1))
!
   ! min / max / sum
   write(string, '(a, i4) ') ' statistics : min/max/average/deviation, ilev = ', ilev
   print*, trim(string)
   print*, trim(bline)
   write(string, trim(fmt_t)) '|', 'varname', '|', ' ', '|', 'MIN', '|', 'MAX', '|', 'AVERAGE', '|', 'DEVIATION', '|'
   print*, trim(string)
   print*, trim(bline)

   do ivar = 1, nvars

   if (nfiles==1) then

   do i = 1, nstats
   if (files_t(1)%vars(ivar)%isvar) then
   write(values(i, 1), '(e15.6) ') files_t(1)%vars(ivar)%level1(i, ilev)
   else
   write(values(i, 1), '(a15) ') '-'
   endif
   write(values(i, 2), '(a15) ') '-'
   write(values(i, 3), '(a15) ') '-'
   enddo
   write(string, trim(fmt_t)) '|', trim(files_t(1)%vars(ivar)%varname), '|', 'File 1', '|', values(1, 1), '|', values(2, 1), '|', values(3, 1), '|', values(4, 1), '|'
   print*, trim(string)

   else
   do i = 1, nstats
   if (files_t(1)%vars(ivar)%isvar) then
   write(values(i, 1), '(e15.6) ') files_t(1)%vars(ivar)%level1(i, ilev)
   else
   write(values(i, 1), '(a15) ') '-'
   endif

   if (files_t(2)%vars(ivar)%isvar) then
   write(values(i, 2), '(e15.6) ') files_t(2)%vars(ivar)%level1(i, ilev)
   else
   write(values(i, 2), '(a15) ') '-'
   endif

   if (files_t(1)%vars(ivar)%isvar.and.files_t(2)%vars(ivar)%isvar) then
   write(values(i, 3), '(e15.6) ') files_t(2)%vars(ivar)%level1(i, ilev)-files_t(1)%vars(ivar)%level1(i, ilev)
   else
   write(values(i, 3), '(a15) ') '-'
   endif
   enddo

   write(string, trim(fmt_t)) '|', '', '|', 'File 1', '|', values(1, 1), '|', values(2, 1), '|', values(3, 1), '|', values(4, 1), '|'
   print*, trim(string)
   write(string, trim(fmt_t)) '|', trim(files_t(1)%vars(ivar)%varname), '|', 'File 2', '|', values(1, 2), '|', values(2, 2), '|', values(3, 2), '|', values(4, 2), '|'
   print*, trim(string)
   write(string, trim(fmt_t)) '|', '', '|', 'Diff.', '|', values(1, 3), '|', values(2, 3), '|', values(3, 3), '|', values(4, 3), '|'
   print*, trim(string)

   endif

   print*, trim(sline)
   enddo

   write(string, '(a) ') ' '
   print*, trim(string)
   write(string, '(a) ') ' '
   print*, trim(string)

   deallocate(values)
!
   return
   end subroutine printfilestatus
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine printfileanalysis(file_t, file_db)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_file_t), intent(in   ) :: file_t
   type(kim_file_t), intent(in   ) :: file_db
! local variables
   integer(i4) :: nvars, ivar, i, idx
   integer(i4) :: n_npts, n_ltmin, n_gtmax, n_le1sig, n_le2sig, n_le3sig, n_gt3sig
   real(r8) :: r_npts, r_ltmin, r_gtmax, r_le1sig, r_le2sig, r_le3sig, r_gt3sig
   character(len=256) :: fmt_t, fmt_c
   character(len=512) :: sline, bline
   character(len=512) :: string
!
   nvars = file_t%nvars
!
   write(fmt_t, '(a) ') '(a1, a11, a1, a17, a1, a14, a1, a14, a1, a17, a1, a17, a1, a17, a1, a17, a1) '
   write(fmt_c, '(a) ') '(a1, a11, a1, i9, a2, f5.1, a2, i6, a2, f5.1, a2, i6, a2, f5.1, a2, i9, a2, f5.1, a2, i9, a2, f5.1, a2, i9, a2, f5.1, a2, i9, a2, f5.1, a2) '
   write(sline, '(a) ') '-------------------------------------------------------------------------------------------------------------------------------------'
   write(bline, '(a) ') '===================================================================================================================================== '
!
   ! min / max / sum
   write(string, '(a) ') ' statistics'
   print*, trim(string)
   print*, trim(bline)
   write(string, trim(fmt_t)) '|', 'varname', '|', '# of points(%) ', '|', '<min(%) ', '|', '>max(%) ', '|', '1 sigma(68.2%) ', '|', '2 sigma(27.2%) ', '|', '3 sigma(4.2%) ', '|', '>3sigma(0.2%) ', '|'
   print*, trim(string)
   print*, trim(bline)

   do ivar = 1, nvars

   idx = getvarindex(file_db, trim(file_t%vars(ivar)%varname))

   if (idx.gt.0) then

   n_npts = file_t%vars(ivar)%anal%npts(1)
   n_ltmin = file_t%vars(ivar)%anal%npts(2)
   n_gtmax = file_t%vars(ivar)%anal%npts(3)
   n_le1sig = file_t%vars(ivar)%anal%npts(4)
   n_le2sig = file_t%vars(ivar)%anal%npts(5)
   n_le3sig = file_t%vars(ivar)%anal%npts(6)
   n_gt3sig = file_t%vars(ivar)%anal%npts(7)

   r_npts = 100.0_r8
   r_ltmin = real(n_ltmin, 8)/real(n_npts, 8)*100.0_r8
   r_gtmax = real(n_gtmax, 8)/real(n_npts, 8)*100.0_r8
   r_le1sig = real(n_le1sig, 8)/real(n_npts, 8)*100.0_r8
   r_le2sig = real(n_le2sig, 8)/real(n_npts, 8)*100.0_r8
   r_le3sig = real(n_le3sig, 8)/real(n_npts, 8)*100.0_r8
   r_gt3sig = real(n_gt3sig, 8)/real(n_npts, 8)*100.0_r8


   write(string, trim(fmt_c)) '|', trim(file_t%vars(ivar)%varname), '|', n_npts, '(', r_npts, ') |', n_ltmin, '(', r_ltmin, ') |', n_gtmax, '(', r_gtmax, ') |', n_le1sig, '(', r_le1sig, ') |', n_le2sig, '(', r_le2sig, ') |', n_le3sig, '(', r_le3sig, ') |', n_gt3sig, '(', r_gt3sig, ') |'
   print*, trim(string)

   print*, trim(sline)

   endif

   enddo

   write(string, '(a) ') ' '
   print*, trim(string)
   write(string, '(a) ') ' '
   print*, trim(string)
!
   return
   end subroutine printfileanalysis
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine printfileminmax(file_t, file_db)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_file_t), intent(in   ) :: file_t
   type(kim_file_t), intent(in   ) :: file_db
! local variables
   integer(i4) :: nvars, ivar, cnts, idx
   integer(i4) :: imin, imax, nmins, nmaxs, iup, iface, ie, je
   real(r8) :: lon, lat
   real(r8) :: value, mu, sigma, minmax, width, nwidth
   character(len=256) :: fmt_t, fmt_c
   character(len=512) :: sline, bline
   character(len=512) :: string
!
   nvars = file_t%nvars
!
   write(fmt_t, '(a) ') '(a1, a11, a1, a9, a1, a11, a1, a14, a1, a7, a1, a7, a1, a15, a1, a15, a1, a15, a1, a15, a1, a6, a1) '
   write(fmt_c, '(a) ') '(a1, a11, a1, a9, a1, i11, a1, i14, a1, f7.2, a1, f7.2, a1, e15.6, a1, e15.6, a1, e15.6, a1, e15.6, a1, f6.2, a1) '
   write(sline, '(a) ') '-----------------------------------------------------------------------------------------------------------------------------------------'
   write(bline, '(a) ') '========================================================================================================================================= '

!-----------------------------------------
! MIN / MAX / SUM
   write(string, '(a) ') ' statistics(min/max), width =(mu-value) or(value-mu), nwidth =(mu-value)/(mu-min) or(value-mu)/(max-mu) '
   print*, trim(string)
   print*, trim(bline)
   write(string, trim(fmt_t)) '|', 'varname', '|', 'MIN/max', '|', 'imin/imax', '|', 'iUP', '|', 'lon', '|', 'lat', '|', 'average', '|', 'min or max', '|', 'value', '|', 'width', '|', 'nwidth', '|'
   print*, trim(string)
   print*, trim(bline)

   do ivar = 1, nvars

   idx = getvarindex(file_db, trim(file_t%vars(ivar)%varname))

   if (idx.gt.0) then

   cnts = 0
   nmins = file_t%vars(ivar)%anal%nmins
   nmaxs = file_t%vars(ivar)%anal%nmaxs

   if (nmins.gt.0) then
   do imin = 1, nmins
   iup = file_t%vars(ivar)%anal%coord_min(imin)%iup
!           iface = file_t%vars(ivar)%anal%coord_min(imin)%iface
!           ie    = file_t%vars(ivar)%anal%coord_min(imin)%ie
!           je    = file_t%vars(ivar)%anal%coord_min(imin)%je
   lon = file_t%vars(ivar)%anal%coord_min(imin)%lon
   lat = file_t%vars(ivar)%anal%coord_min(imin)%lat
   value = file_t%vars(ivar)%anal%value_min(imin)%value
   mu = file_t%vars(ivar)%anal%value_min(imin)%mu
   minmax = file_t%vars(ivar)%anal%value_min(imin)%minmax
   width = file_t%vars(ivar)%anal%value_min(imin)%width
   nwidth = file_t%vars(ivar)%anal%value_min(imin)%nwidth
   write(string, fmt_c) '|', trim(file_t%vars(ivar)%varname), '|', 'MIN', '|', imin, '|', iup, '|', lon, '|', lat, '|', mu, '|', minmax, '|', value, '|', width, '|', nwidth, '|'
   print*, trim(string)
   enddo
   print*, trim(sline)
   endif

   if (nmaxs.gt.0) then
   do imax = 1, nmaxs
   iup = file_t%vars(ivar)%anal%coord_max(imax)%iup
!           iface = file_t%vars(ivar)%anal%coord_max(imax)%iface
!           ie    = file_t%vars(ivar)%anal%coord_max(imax)%ie
!           je    = file_t%vars(ivar)%anal%coord_max(imax)%je
   lon = file_t%vars(ivar)%anal%coord_max(imax)%lon
   lat = file_t%vars(ivar)%anal%coord_max(imax)%lat
   value = file_t%vars(ivar)%anal%value_max(imax)%value
   mu = file_t%vars(ivar)%anal%value_max(imax)%mu
   minmax = file_t%vars(ivar)%anal%value_max(imax)%minmax
   width = file_t%vars(ivar)%anal%value_max(imax)%width
   nwidth = file_t%vars(ivar)%anal%value_max(imax)%nwidth
   write(string, fmt_c) '|', trim(file_t%vars(ivar)%varname), '|', 'MAX', '|', imax, '|', iup, '|', lon, '|', lat, '|', mu, '|', minmax, '|', value, '|', width, '|', nwidth, '|'
   print*, trim(string)
   enddo
   print*, trim(sline)
   endif

   endif

   enddo

   write(string, '(a) ') ' '
   print*, trim(string)
   write(string, '(a) ') ' '
   print*, trim(string)
!
   return
   end subroutine printfileminmax
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine nc_check(ncstat_in)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: ncstat_in
!
   if (ncstat_in.ne.nf90_noerr) then
     print*,trim(nf90_strerror(ncstat_in))
     stop
   endif
!
   return
   end subroutine nc_check
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module kim_std_vars
!-------------------------------------------------------------------------------
