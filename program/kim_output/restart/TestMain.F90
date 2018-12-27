PROGRAM Main
 USE CheckRestart, ONLY: Ini, Run, Fin
 IMPLICIT NONE
 CHARACTER(LEN=512) :: filename1, filename2

  filename1 = '/scratch/jhkim/TestBed/Data/2.4.04_01/10h/SW_G/restart/initialRun.r0000000040.nc'
  filename2 = '/scratch/jhkim/TestBed/Data/2.4.04_01/10h/SW_G/restart/restart-20110725130000-000001.nc'
  CALL Ini(filename1, filename2)

  CALL Run()

  CALL Fin()

END PROGRAM Main
