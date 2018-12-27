&inputs

iExp = 3
!iExp = 4064


! which_method = 0: polynomial, 1: use other coefficient
which_method = -1
! adjust top and surface eta
useAdjustTop = .FALSE.
useAdjustSfc = .FALSE.
topEta       = 0.0D0
sfcEta       = 0.0D0
! shifting
useShift     = .FALSE.
shift        = 0.0D0

! for which_method = 0
nPoly = 5
coefPoly(2:5) = 0.002, 0.0, 4.0, -3.002

! for which_method = 1
!  which_coef = 0: ECMWF, 1: GRIMs
which_coef  = -1
nlevs_other = -1
apply_B     = .FALSE.

! pressure, hybrid, sigma
autoDet     = .FALSE.
nprs        = -1
nhyb        = -1

! for output (pressure and height)
useStdAtm = .TRUE.

/
