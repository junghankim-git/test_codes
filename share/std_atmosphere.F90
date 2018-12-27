!-------------------------------------------------------------------------------
   module std_atmosphere
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    2015-08-02  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only:i4, r8, r16
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   integer(i4), parameter :: natmoslayers = 8
   real(r8), dimension(natmoslayers) :: basetemp, lapserate
   real(r8), dimension(natmoslayers+1) :: altitude, basepres
   real(r8), parameter, public :: std_g = 9.80665_r8
!
   !real(r8), parameter, public :: std_ps     = 100000.0_r8
   !real(r8), parameter, public :: std_ps     = 101325.0_r8
   !real(r8), parameter, public :: stdrgas   =  287.1_r8
   !real(r8), parameter, public :: stdrgas   =  287.058_r8
   !real(r8), parameter, public :: stdrgas   =  287.053_r8
   !real(r8), parameter, public :: stdrgas   =  287.04_r8
   !real(r8), parameter, public :: stdkelvin =   273.15_r8
   !real(r8), parameter, public :: stdkelvin =   273.0_r8
   !real(r8), parameter, public :: std_ps     = 101325.0_r8
   !real(r8), parameter, public :: stdrgas   =  287.04_r8
   !real(r8), parameter, public :: stdkelvin = 273.15_r8
!
   real(r8), parameter, public :: std_ps     = 101325.0_r8
   real(r8), parameter, public :: stdrgas   = 287.053_r8
   real(r8), parameter, public :: stdkelvin = 273.15_r8
!
   ! http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/atmos/atmos.html
   !real(r8), parameter, public         :: std_ps     = 101300.0_r8
   !real(r8), parameter, public         :: stdrgas   =  287.1_r8
   !real(r8), parameter, public         :: stdkelvin = 273.2_r8
!
   public :: std_atm_initialize
   public :: get_temperature ! height to temperature
   public :: get_height ! pressure to height
   public :: get_pressure ! height to pressure
!
   contains
!
! International Standard Atmosphere
!================================================================================================
! Atmosphere  :      height(km)    :Lapse Rate(K/km): Base Temperature:Base Pressure (Pa)
!================================================================================================
! Troposphere :0.0    < z <= 11.000:     -6.6       :+15.0  (288.15 K):  101325.0000
! Tropopause  :11.000 < z <= 20.000:     +0.0       :-56.5  (216.65 K):   22632.0000
! Stratosphere:20.000 < z <= 32.000:     +1.0       :-56.5  (216.65 K):    5474.9000
! Stratosphere:32.000 < z <= 47.000:     +2.8       :-44.5  (228.65 K):     868.0200
! Stratopause :47.000 < z <= 51.000:     +0.0       : -2.5  (270.65 K):     110.9100
! Mesosphere  :51.000 < z <= 71.000:     -2.8       : -2.5  (270.65 K):      66.9390
! Mesosphere  :71.000 < z <= 86.852:     -2.0       :-58.5  (214.65 K):       3.9564
! Mesopause   :84.852 < z          :         .      :-86.28 (186.87 K):       0.3734
!================================================================================================
! http://en.wikipedia.org/wiki/International_Standard_Atmosphere
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine std_atm_initialize(standard)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: standard
   ! standard = 1: ISA
   ! standard = 2: U.S.
! local variables
   integer(i4) :: k
!
   ! altitude (g: constant)
   altitude(1) = 0.0_r8
   altitude(2) = 11000.0_r8
   altitude(3) = 20000.0_r8
   altitude(4) = 32000.0_r8
   altitude(5) = 47000.0_r8
   altitude(6) = 51000.0_r8
   altitude(7) = 71000.0_r8
   altitude(8) = 84852.0_r8
   ! temp
   altitude(9) = 200000.0_r8
!
   ! base temperature
   basetemp(1) = stdkelvin+15.0_r8
   basetemp(2) = stdkelvin-56.5_r8
   basetemp(3) = stdkelvin-56.5_r8
   basetemp(4) = stdkelvin-44.5_r8
   basetemp(5) = stdkelvin-2.5_r8
   basetemp(6) = stdkelvin-2.5_r8
   basetemp(7) = stdkelvin-58.5_r8
   basetemp(8) = stdkelvin-86.28_r8
!
! lapse rate: method 1
!   do k = 1,natmoslayers-1
!     lapserate(k) = (basetemp(k+1)-basetemp(k))/(altitude(k+1)-altitude(k))
!   end do
!   lapserate(natmoslayers) = 0.0_r8
!
   ! lapse rate: method 2
   lapserate(1) =-6.5d-3
   lapserate(2) =+0.0d-3
   lapserate(3) =+1.0d-3
   lapserate(4) =+2.8d-3
   lapserate(5) =+0.0d-3
   lapserate(6) =-2.8d-3
   lapserate(7) =-2.0d-3
   lapserate(8) =+0.0d-3
!
   ! not fix pressure
   basepres(1) = std_ps
   do k = 2,natmoslayers
     if (lapserate(k-1).eq.0.0_r8) then
       basepres(k) = basepres(k-1)*dexp(-std_g/stdrgas/basetemp(k-1)*(altitude(k)-altitude(k-1)))
     else
       basepres(k) = basepres(k-1)*  (1.0_r8+lapserate(k-1)/basetemp(k-1)*(altitude(k)-altitude(k-1)))**(-std_g/stdrgas/lapserate(k-1))
     endif
   enddo
!
! fix base pressure
!   if (standard == 1) then
!     basepres(1) = std_ps
!     basepres(2) = 22632.0_r8
!     basepres(3) =  5474.9_r8
!     basepres(4) =  868.02_r8
!     basepres(5) =  110.91_r8
!     basepres(6) =  66.939_r8
!     basepres(7) =  3.9564_r8
!     basepres(8) =  0.3734_r8
!     basepres(9) =  0.0_r8
!   else
!     basepres(1) = std_ps
!     basepres(2) = 22632.1_r8
!     basepres(3) =  5474.89_r8
!     basepres(4) =  868.019_r8
!     basepres(5) =  110.906_r8
!     basepres(6) =  66.9389_r8
!     basepres(7) =  3.95642_r8
!     basepres(8) =  0.3734_r8
!     basepres(9) =  0.0_r8
!   end if
!
   return
   end subroutine std_atm_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_temperature(height)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: height
   real(r8) :: get_temperature
! local variables
   integer(i4) :: k
!   if (height > altitude(natmoslayers)) then
!     print *,'out of height ...'
!     stop
!   end if
!
   do k = 1,natmoslayers-1
     if (height==altitude(1)) then
       get_temperature = basetemp(1)
     endif
     if (height>altitude(k).and.height<=altitude(k+1)) then
       get_temperature = basetemp(k)+lapserate(k)*(height-altitude(k))
     endif
   enddo
!
   if (height>altitude(natmoslayers)) then
     get_temperature = basetemp(natmoslayers)+lapserate(natmoslayers)*(height-altitude(natmoslayers))
   endif
!
   end function get_temperature
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_height(pressure)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: pressure
   real(r8) :: get_height
! local variables
   integer(i4) :: k
!
   do k = 1,natmoslayers
     if (pressure<basepres(k).and.pressure>= basepres(k+1)) then
       if (lapserate(k).eq.0.0_r8) then
         get_height = altitude(k)+stdrgas*basetemp(k)/std_g*dlog(basepres(k)/pressure)
       else
         get_height = altitude(k)+1.0_r8/lapserate(k)*basetemp(k)*((basepres(k)/pressure)**(lapserate(k)*stdrgas/std_g)-1.0_r8)
       endif
     endif
   enddo
   if (pressure==basepres(1)) get_height = 0.0_r8
!
   end function get_height
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_pressure(height)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: height
   real(r8) :: get_pressure
! local variables
   integer(i4) :: k
!
   do k = 1,natmoslayers
     if (height>altitude(k).and.height<=altitude(k+1)) then
       if (lapserate(k).eq.0.0_r8) then
         get_pressure = basepres(k)*dexp(-std_g/stdrgas/basetemp(k)*(height-altitude(k)))
       else
         get_pressure = basepres(k)*(1.0_r8+lapserate(k)/basetemp(k)*(height-altitude(k)))**(-std_g/lapserate(k)/stdrgas)
       endif
     endif
   enddo
   if (height==altitude(1)) get_pressure = std_ps
!
   end function get_pressure
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_height_int(pressure)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: pressure
   real(r8) :: get_height_int
! local variables
   integer(i4), parameter :: n = 9000000
   real(r8), parameter :: htop = 90000.0_r8
   real(r8) :: dh = htop/n
   real(r8) :: left, right
   real(r8) :: h, t
   integer(i4) :: k
   logical :: lfind
!
   lfind = .false.
!
   right = stdrgas/std_g*dlog(std_ps/pressure)
!
   get_height_int = 0.0_r8
   h = 0.0_r8
   left = 0.0_r8
!
   do k = 1,n
     left = left+1.0_r8/get_temperature(h)*dh
     if (left>= right) then
       get_height_int = h
       lfind = .true.
       exit
     endif
     h = h+dh
   enddo
!
   if (pressure==0.0_r8) then
     get_height_int = 200000.0_r8
     lfind = .true.
   endif
!
   if (.not.lfind) then
     print*,'Could not find height...'
     stop
   endif
!
   end function get_height_int
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module std_atmosphere
!-------------------------------------------------------------------------------
