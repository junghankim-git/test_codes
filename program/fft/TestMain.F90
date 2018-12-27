PROGRAM TestMain
  IMPLICIT NONE
#include <fftw3.f>
  INTEGER, PARAMETER                      :: i4 = 4, i8 = 8, r8 = 8
  REAL(r8), PARAMETER                     :: pi = 3.14159265358979323846264338_r8

  ! data
  INTEGER(i4)                             :: n, testcase
  REAL(r8)                                :: zmin, zmax, dz, z, freq
  REAL(r8)                                :: nwave

  ! output
  INTEGER(i4)                             :: unitnum

  COMPLEX(r8), DIMENSION(:), ALLOCATABLE  :: input, freqs, differential, integration, tmp
  INTEGER(i8)                             :: forward_p, backward_p

  INTEGER(i4)  :: i
  

  n        = 100
  testcase = 5 ! 1(sin(x)), 2(cos(x)), 3(y=(2pix)), 4(y=x^2), 5(y=x), 6(step function)
  nwave    = 1.0_r8
  IF(testcase == 1) THEN
    zmin     = 0.0_r8
    zmax     = 2.0_r8*pi*nwave
  ELSE IF(testcase == 2) THEN
    zmin     = 0.0_r8
    zmax     = 2.0_r8*pi*nwave
  ELSE IF(testcase == 3) THEN
    zmin     = 0.0_r8
    zmax     = 1.0_r8
  ELSE IF(testcase == 4) THEN
    zmin     = -1.0_r8
    zmax     = 1.0_r8
  ELSE IF(testcase == 5) THEN
    zmin     = 0.0_r8
    zmax     = 1.0_r8
  ELSE IF(testcase == 6) THEN
    zmin     = 0.0_r8
    zmax     = 1.0_r8
  END IF

  ! Allocate Variables
  ALLOCATE(input(n))
  ALLOCATE(freqs(n))
  ALLOCATE(differential(n))
  ALLOCATE(integration(n))

  CALL WriteMetaData('info.dat', testcase, n, zmin, zmax)
  CALL OpenFile('Result.dat', unitnum)


!########### Frequency 
!  CALL dfftw_plan_dft_r2c_1d(forward_p, N, input, freqs, FFTW_MEASURE)
  CALL dfftw_plan_dft_1d(forward_p, N, input, freqs, FFTW_FORWARD, FFTW_MEASURE)


  IF(testcase == 1) THEN
    DO i = 1, N
      z = (zmax-zmin)/DBLE(N)*DBLE(i-1)+zmin
      input(i) = DSIN(z)
    END DO
  ELSE IF(testcase == 2) THEN
    DO i = 1, N
      z = (zmax-zmin)/DBLE(N)*DBLE(i-1)+zmin
      input(i) = DCOS(z)
    END DO
  ELSE IF(testcase == 3) THEN
    DO i = 1, N
      z = (zmax-zmin)/DBLE(N)*DBLE(i-1)+zmin
      input(i) = DCOS(2.0_r8*pi*z)
    END DO
  ELSE IF(testcase == 4) THEN
    DO i = 1, N
      z = (zmax-zmin)/DBLE(N)*DBLE(i-1)+zmin
      input(i) = z*z
    END DO
  ELSE IF(testcase == 5) THEN
    DO i = 1, N
      z = (zmax-zmin)/DBLE(N)*DBLE(i-1)+zmin
      input(i) = z
    END DO
  ELSE IF(testcase == 6) THEN
    DO i = 1, N
      IF(i .LE. N/2+1) THEN
        !input(i) = 0.0_r8
        input(i) = COMPLEX(0.0_r8, 0.0_r8)
      ELSE
        !input(i) = 1.0_r8
        input(i) = COMPLEX(1.0_r8, 0.0_r8)
      END IF
    END DO
  END IF
  CALL WriteData_Cmplx(unitnum, n, input)



  freqs(:) = COMPLEX(0.0_r8, 0.0_r8)
  CALL dfftw_execute_dft(forward_p, input, freqs)

  CALL WriteData_Cmplx(unitnum, n, freqs)
  tmp = freqs



!########### Differential
  CALL dfftw_plan_dft_1d(backward_p, N, freqs, differential, FFTW_BACKWARD, FFTW_MEASURE)
  freqs = tmp

  DO i = 1, n
    IF(i .LE. n/2+1) THEN
      freq = DBLE(i-1)
    ELSE
      freq = -DBLE(n-i+1)
    ENDIF
! remove high frequency
    IF(i .GE. n/3 .AND. i .LE. 2*N/3) THEN
      freqs(i) = COMPLEX(0.0_r8, 0.0_r8)
    END IF
    freqs(i) = 2.0_r8*pi*freq/(DBLE(N)*(zmax-zmin))*COMPLEX(-DIMAG(freqs(i)), DBLE(freqs(i)))
  END DO

  CALL dfftw_execute_dft(backward_p, freqs, differential)

  CALL WriteData_Cmplx(unitnum, n, differential)



!########### integration
  !tmp = freqs
  CALL dfftw_plan_dft_1d(backward_p, N, freqs, integration, FFTW_BACKWARD, FFTW_MEASURE)
  freqs = tmp

  DO i = 1, n
    IF(i .LE. n/2+1) THEN
      freq = DBLE(i-1)
    ELSE
      freq = -DBLE(n-i+1)
    ENDIF
    IF(i==1) freq = 1.0_r8
! remove high frequency
    IF(i .GE. n/3 .AND. i .LE. 2*N/3) THEN
      freqs(i) = COMPLEX(0.0_r8, 0.0_r8)
    END IF
    freqs(i) = (zmax-zmin)/(2.0_r8*pi*DBLE(N)*freq)*COMPLEX(DIMAG(freqs(i)), -DBLE(freqs(i)))
    !freqs(i) = (2.0_r8*pi*DBLE(N)*freq)*COMPLEX(DIMAG(freqs(i)), -DBLE(freqs(i)))
  END DO

  CALL dfftw_execute_dft(backward_p, freqs, integration)

  CALL WriteData_Cmplx(unitnum, n, integration)









  CALL dfftw_destroy_plan(forward_p)
  CALL dfftw_destroy_plan(backward_p)


  CALL CloseFile(unitnum)

  ! Deallocate Variables
  DEALLOCATE(input)
  DEALLOCATE(freqs)
  DEALLOCATE(differential)
  DEALLOCATE(integration)


END PROGRAM TestMain





SUBROUTINE WriteMetaData(filename, testcase, n, zmin, zmax)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)           :: filename
  INTEGER(4), INTENT(IN)             :: testcase, n
  REAL(8), INTENT(IN)                :: zmin, zmax

  ! local
  INTEGER(4)     :: unitnum

  unitnum = 21

  OPEN(unitnum, FILE=filename, STATUS='UNKNOWN', FORM='FORMATTED')

  WRITE(unitnum, *) testcase
  WRITE(unitnum, *) n
  WRITE(unitnum, *) zmin
  WRITE(unitnum, *) zmax

  CLOSE(unitnum)

END SUBROUTINE WriteMetaData






SUBROUTINE OpenFile(filename, unitnum)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)           :: filename
  INTEGER(4), INTENT(IN)             :: unitnum

  OPEN(unitnum, FILE=filename, STATUS='UNKNOWN', FORM='FORMATTED')


END SUBROUTINE OpenFile






SUBROUTINE CloseFile(unitnum)
  IMPLICIT NONE
  INTEGER(4), INTENT(IN)             :: unitnum

  CLOSE(unitnum)

END SUBROUTINE CloseFile






SUBROUTINE WriteData_Real(unitnum, n, values)
  IMPLICIT NONE
  INTEGER(4), INTENT(IN)             :: unitnum
  INTEGER(4), INTENT(IN)             :: n
  REAL(8), DIMENSION(n), INTENT(IN)  :: values

  ! local
  INTEGER(4)     :: i

  DO i = 1, n
    WRITE(unitnum, *) values(i), 0.0D0
  END DO

END SUBROUTINE WriteData_Real





SUBROUTINE WriteData_Cmplx(unitnum, n, values)
  IMPLICIT NONE
  INTEGER(4), INTENT(IN)             :: unitnum
  INTEGER(4), INTENT(IN)                   :: n
  DOUBLE COMPLEX, DIMENSION(n), INTENT(IN) :: values

  ! local
  INTEGER(4)     :: i

  DO i = 1, n
    WRITE(unitnum, *) DBLE(values(i)), DIMAG(values(i))
  END DO

END SUBROUTINE WriteData_Cmplx
