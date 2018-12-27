PROGRAM GenPost
 USE ResHis2Post
 IMPLICIT NONE

 CHARACTER(LEN=512) :: resfile, hisfile, outfile

 resfile = '/scratch/jhkim/TestBed/KIM/Output/restart/v2.4/initialRun.r0000000040.nc'
 hisfile = '/scratch/jhkim/TestBed/KIM/Output/2.4/ne30/gnu/10h/SW_G/101/UP-20110725130000-000002.nc'
 outfile = './out.nc'

 IF (IARGC() .LT. 3) THEN
   PRINT *, '# Need arguments... (./genpos restarfile.nc historyfile.nc outputfile.nc)'
   PRINT *, '# PROGRAM : STOP'
   STOP
 END IF

 CALL GETARG(1,resfile)
 CALL GETARG(2,hisfile)
 CALL GETARG(3,outfile)
 PRINT *, '# PROGRAM : START'
 PRINT *, ' * setting files'
 PRINT *, '  - restart file: ', TRIM(resfile)
 PRINT *, '  - history file: ', TRIM(hisfile)
 PRINT *, '  - output  file: ', TRIM(outfile)
 PRINT *, ' '

 PRINT *, ' * initialize '
 CALL Ini(TRIM(resfile), TRIM(hisfile), TRIM(outfile))

 PRINT *, ' * run'
 CALL Run()

 PRINT *, ' * finalize'
 CALL Fin()

 PRINT *, '# PROGRAM : END'

END PROGRAM GenPost
