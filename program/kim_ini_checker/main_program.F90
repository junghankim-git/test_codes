
program main
use kinds
use chk_kim_ini, only: kim_ini_t, initialize, finalize, analysis
implicit none
!
integer(i4)        :: nfiles
type(kim_ini_t)    :: kim1_t, kim2_t
character(len=512) :: filename1, filename2

filename1 = '/s3/data/kiaps/kiaps-val/khseol/KIM_INPUT_3.0.08/NE240L50/2017060100/20170601_00UTC_SW_ENDGameUM_50lev.nc'
filename2 = '/s3/scratch/kiaps/kiaps-sys/jhkim/final/UP-20170601000000-000000.nc'

nfiles = iargc()

if (nfiles.ge.2) then
  call getarg(1,filename1)
  call getarg(2,filename2)
endif

call initialize(kim1_t,trim(filename1))
call initialize(kim2_t,trim(filename2))
call analysis(kim1_t,kim2_t)
call finalize(kim1_t)
call finalize(kim2_t)

end program main
