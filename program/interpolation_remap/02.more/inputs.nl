&inputs
! spline (spline remap or not)
use_spline = .true.

! method
order    = 2        ! order: 0(constant), 1(linear), 2(parabolic), 3(cubic), 4(quartic)
bndry    = 2        ! bndry: 0(mirror), 1(mirror+slope), 2(periodic)
ismono   = .true.

! testcase (1:1st poly, 2:2rd poly, 3:3rd, 4:harmonic, 5:gaussian, 6:step, 7:some)
testcase = 5

! source grid size
n        = 12
!n        = 91
!n        = 400
!n        = 5

! target grid size
m        = 11
!m        = 90
!m        = 399
!m        = 17

! inner function
l        = 1000
/
