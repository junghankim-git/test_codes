PROGRAM test
 IMPLICIT NONE
 INTEGER, PARAMETER        :: i4 = 8, r8 = 8
 INTEGER(i4), PARAMETER    :: n = 3
 INTEGER(i4), DIMENSION(n) :: my_n
 INTEGER(i4)               :: my_gcd, my_lcm
!
 integer(i4) :: num, nn
 integer(i4), dimension(:), allocatable :: lnn


!  my_n(1) = 5
!  my_n(2) = 2
!  my_n(3) = 3

  !my_gcd = GCD(my_n(1), my_n(2))
!  my_gcd = GetGCD(n, my_n)
!  my_lcm = GetLCM(n, my_n)

!  PRINT *, 'nums = ', my_n(:)
!  PRINT *, 'GCD  = ', my_gcd
!  PRINT *, 'LCM  = ', my_lcm

  !num = 13195
  num = 600851475143_i4
  !print *, huge(num)
  !call get_somes(num,nn,lnn)
  call get_somes(num,nn,lnn)
 
  print *, nn
  print *, lnn
  deallocate(lnn)

 CONTAINS


 FUNCTION GetLCM(n, nums) RESULT(res)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN)   :: n
  INTEGER(i4), DIMENSION(n) :: nums
  INTEGER(i4)               :: res
  ! local
  INTEGER(i4) :: i, mylcm

   res = LCM(nums(1), nums(2))
   DO i = 3, n
     res = LCM(res, nums(i))
   END DO

 END FUNCTION GetLCM




 ! Least Common Multiple
 RECURSIVE FUNCTION LCM(a, b) RESULT(res)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN) :: a, b
  INTEGER(i4)             :: res
  ! local
  INTEGER(i4) :: mygcd

  mygcd = GCD(a, b)
  res   = a*b/mygcd

 END FUNCTION LCM




 FUNCTION GetGCD(n, nums) RESULT(res)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN)   :: n
  INTEGER(i4), DIMENSION(n) :: nums
  INTEGER(i4)               :: res
  ! local
  INTEGER(i4) :: i

   res = GCD(nums(1), nums(2))
   DO i = 3, n
     res = GCD(res, nums(i))
   END DO

 END FUNCTION GetGCD




 ! Greatest Commin Divisor
 RECURSIVE FUNCTION GCD(a, b) RESULT(res)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN) :: a, b
  INTEGER(i4)             :: res

  IF (b == 0) THEN
    res = a
  ELSE
    res = GCD(b, MOD(a, b))
  END IF

 END FUNCTION GCD


 ! Divisors
 SUBROUTINE GetDivisor(num, ndivs, divs)
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN)  :: num
  INTEGER(i4), INTENT(OUT) :: ndivs
  INTEGER(i4), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: divs
  ! local
  INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ldivs
  INTEGER(i4) :: div, nlds, nrds

   ALLOCATE(ldivs(num))
   ! left
   nlds       = 1
   ldivs(1)   = 1
   ! right
   nrds       = 1
   ldivs(num) = num

   div = 2
   DO WHILE(div <= SQRT(REAL(num)))

     IF (MOD(num,div) .EQ. 0) THEN
       nlds  = nlds + 1
       ldivs(nlds) = div
       IF (div .NE. (num/div)) THEN
         nrds = nrds + 1
         ldivs(num-nrds+1) = num/div
       END IF
     END IF

     div = div + 1

   END DO


   ndivs = nlds+nrds
   ALLOCATE(divs(ndivs))
   divs(1:nlds)     = ldivs(1:nlds)
   divs(nlds+1:ndivs) = ldivs(num-nrds+1:num)

   DEALLOCATE(ldivs)

 END SUBROUTINE GetDivisor


   subroutine get_somes(num, n, somes)
   implicit none
   integer(i4), intent(in   ) :: num
   integer(i4), intent(  out) :: n
   integer(i4), dimension(:), allocatable, intent(  out) :: somes
!
   integer(i4), parameter :: max_n = 10000
   integer(i4) :: i, n1, n2
   integer(i4), dimension(:), allocatable :: tmp1, tmp2
   integer(i4), dimension(max_n) :: res
!
print *, 'start'
   call getdivisor(num, n1, tmp1)
print *, 'n1 = ', n1
print *, 'tmp1 = ', tmp1
!
#if 0
   n = 0
   do i=1,n1
     call getdivisor(tmp1(i),n2,tmp2)
print *, 'n2 = ', n2
     if (n2.eq.2) then
       n = n + 1
       res(n) = tmp1(i)
     endif
     deallocate(tmp2)
   enddo
   allocate(somes(n))
   somes(1:n) = res(1:n)
#endif
!
   end subroutine get_somes
!
END PROGRAM test

