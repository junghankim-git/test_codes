PROGRAM test
 IMPLICIT NONE
 INTEGER, PARAMETER        :: i4 = 8, r8 = 8
 INTEGER(i4), PARAMETER    :: n = 3
 INTEGER(i4), DIMENSION(n) :: my_n
 INTEGER(i4)               :: my_gcd, my_lcm
!
 integer(i4) :: num, nn
 integer(i4), dimension(:), allocatable :: lnn
 integer(i4) :: i



!  !num = 13195_i4
!  num = 600851475143_i4
!  call get_somes(num,nn,lnn)
!  print *, nn
!  print *, lnn
!  deallocate(lnn)
!!
!  num = 6857_i4
!  call getdivisor(num,nn,lnn)
!  print *, nn
!  print *, lnn
!  deallocate(lnn)
!
!  num = 417_i4
!  call getdivisor(num,nn,lnn)
!  print *, nn
!  print *, lnn
!  deallocate(lnn)
!  num = 418_i4
!  call getdivisor(num,nn,lnn)
!  print *, nn
!  print *, lnn
!  deallocate(lnn)
!  num = 420_i4
!  call getdivisor(num,nn,lnn)
!  print *, nn
!  print *, lnn
!  deallocate(lnn)
  num = 237600_i4
  call getdivisor(num,nn,lnn)
  print *, nn
  do i=1,nn,2
  print *, lnn(i+1), lnn(i)
  enddo
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
 SUBROUTINE GetDivisor_old(num, ndivs, divs)
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

 END SUBROUTINE GetDivisor_old


 ! Divisors
 SUBROUTINE GetDivisor(num, ndivs, divs)
  USE buffer, only: buffer_t, buf_i8, buf_initialize, buf_finalize, buf_add, buf_add_i8
  IMPLICIT NONE
  INTEGER(i4), INTENT(IN)  :: num
  INTEGER(i4), INTENT(OUT) :: ndivs
  INTEGER(i4), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: divs
  ! local
  type(buffer_t) :: buf
  INTEGER(i4) :: div
!
   call buf_initialize(buf, buf_i8)
   call buf_add(buf, 1_i4)
   call buf_add(buf, num)
   ndivs = 2

   div = 2
   DO WHILE(div <= SQRT(REAL(num)))

     IF (MOD(num,div) .EQ. 0) THEN
       ndivs = ndivs+1
       call buf_add(buf, div)
       IF (div .NE. (num/div)) THEN
         ndivs = ndivs+1
         call buf_add(buf, num/div)
       END IF
     END IF

     div = div + 1

   END DO


   ALLOCATE(divs(ndivs))
   divs(:) = buf%i8(:)
!
   call buf_finalize(buf)

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
   call getdivisor(num, n1, tmp1)
!
   n = 0
   do i=1,n1
     call getdivisor(tmp1(i),n2,tmp2)
     if (n2.eq.2) then
       n = n + 1
       res(n) = tmp1(i)
     endif
     deallocate(tmp2)
   enddo
   allocate(somes(n))
   somes(1:n) = res(1:n)
!
   end subroutine get_somes
!
END PROGRAM test

