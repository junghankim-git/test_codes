
program test
implicit none
integer, parameter :: i4 = 4, i8 = 8
integer :: tmp
integer(i8) :: v8

!tmp = 2717909376
!print *,tmp,huge(tmp)

v8 = 1023_i8*1023_i8*4104_i8
print *,v8
print *,4.*huge(tmp)/1000./1000./1000

end program test
