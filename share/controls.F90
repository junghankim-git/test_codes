!-------------------------------------------------------------------------------
   module controls
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    2015-09-30  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only : i4, l4, r4, r8
   use logger, only : printlog, tostring
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   integer(i4), public :: np = 4
   integer(i4), public :: ne = 1
   integer(i4), parameter, public :: ntimelevs = 2
   real(r8), public :: dt = 0.1_r8
   real(r8), public :: velocity = 1.0_r8
!
! experiment
   integer(i4), public :: iexp = 1 ! 1 : square, 2 : gaussian
!
   real(r8), public :: amp = 1.0_r8 ! ampltute(hight of gaussian)
! squre
! gaussian
   real(r8), public :: mu = 20.0_r8 ! mu(initial center position of gaussian)
   real(r8), public :: sigma = 4.0_r8 ! sigma(standard deviation of gaussian)
!
   end module controls
!-------------------------------------------------------------------------------
