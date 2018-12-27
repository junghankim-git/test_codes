



MODULE test
use mymodule1,     only: abc, def,fjdk,&
                                                fjdsakl, &
          fjdoia
use mymodule2, ONLY: N
implicit none

contains

subroutine my_sub1(a, b,c)
   IMPLICIT NONE
integer, dimension(:), intent(in) ::a




REAL, DIMension(:), intent(out) ::b


REAL, DIMension(:), intent(INout) ::c
integer :: i, j, k, var
    !***********************************************
    !
    ! my main algorithm 
    !
    ! calculate y = ax + b
    ! z = \frac {} {} \times psi + delta, ...
    !
    !***********************************************
do i = 1,N
b(i) = REAL(a(i))
             enddo
do k = 1,N
do j = 1,N
do i = 1,N





C(I) = C(I)+B(I)+SOME_CALCULATION ...SOME_CALCULATION... SOME_CALCULATION...SOME_CALCULATION...+A(I)
enddo
enddo
enddo
END SUBROUTINE

function my_fun(a, b,c)
implicit none
integer, dimension(:), intent(in) ::a
REAL, DIMension(:), intent(out) ::b
REAL, DIMension(:), intent(INout) ::c

! my main algorithm 
if (a==b) then
  if (b==c) then
    ! some code



endif
 endif
END SUBROUTINE
  
    end MODULE test
