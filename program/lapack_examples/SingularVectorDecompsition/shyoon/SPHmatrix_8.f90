Program matrix

USE RealSPH

IMPLICIT NONE

integer                   :: l,m,i,j
integer, parameter        :: n=8
integer, dimension(3)     :: ll, mm
real(r8), dimension(n)    :: lat, long
real(r8), dimension(n)    :: cylm,sylm
real(r8), dimension(n,n)  :: mat
real(r8), dimension(2,n)  :: CY
real(r8), dimension(3,n)  :: SY
real(r8), dimension(n)    :: Y00, Y10, Y20
character                 :: JOBU, JOBVT
integer, parameter        :: LDA=8, LDU=8, LDVT=8
integer,parameter         :: LWORK=5*n
real(r8), dimension(n)    :: S
real(r8), dimension(n,n)  :: mat_a, U, W, VT, TMP1, TMP2
real(r8), dimension(5*n)  :: WORK
integer                   :: INFO

!set lat, long

 lat(1)  = asin(1._r8/(sqrt(3._r8)))
 long(1) = 0._r8 
 lat(2)  = -asin(1._r8/(sqrt(3._r8)))
 long(2) = 0._r8

 lat(3)  = asin(1._r8/(sqrt(3._r8)))
 long(3) = pi/2._r8

 lat(4)  = -asin(1._r8/(sqrt(3._r8)))
 long(4) = pi/2._r8

 lat(5)  = asin(1._r8/(sqrt(3._r8)))
 long(5) = pi

 lat(6)  = -asin(1._r8/(sqrt(3._r8)))
 long(6) = pi

 lat(7)  = asin(1._r8/(sqrt(3._r8)))
 long(7) = 3._r8*pi/2._r8

 lat(8)  = -asin(1._r8/(sqrt(3._r8)))
 long(8) = 3._r8*pi/2._r8

!compute Y00, Y10, Y20

do i=1,n
 Y00(i) = 1._r8/(2._r8*(sqrt(pi)))
 Y10(i) = (1._r8/2._r8)*(sqrt(3._r8/pi))*(sin(lat(i)))
 Y20(i) = (1._r8/4._r8)*(sqrt(5._r8/pi))*((3._r8*(sin(lat(i)))**2)-1._r8)
enddo

!read grid

open(10,file='l_m.txt')
do i=1,3
read(10,*) ll(i),mm(i)

!compute SPH

l=ll(i)
m=mm(i)
CALL IniRecurCoef
CALL ComputSPH (l, m, n, lat, long, cylm, sylm)

if (l==1.and.m==1) then
CY(1,:)= cylm(:)
SY(1,:)= sylm(:)
elseif(l==2.and.m==1) then
CY(2,:)= cylm(:)
SY(2,:)= sylm(:)
elseif(l==2.and.m==2) then
SY(3,:)= sylm(:)
endif

enddo

!write matrix

do j=1,n
 do i=1,n
 if (j==1)then
  mat(i,j)= Y00(i)
 elseif (j==2) then
  mat(i,j)= SY(1,i)
 elseif (j==3) then
  mat(i,j)= Y10(i)
 elseif (j==4) then
  mat(i,j)= CY(1,i)
 elseif (j==5) then
  mat(i,j)= SY(3,i)
 elseif (j==6) then
  mat(i,j)= SY(2,i)
 elseif (j==7) then
  mat(i,j)= Y20(i)
 elseif (j==8) then
  mat(i,j)= CY(2,i)
 endif
 enddo
enddo

!output
open(11, file='out.txt')
do i=1,n
write(11,900)(mat(i,j),j=1,n)
enddo
900 format (8(F23.15))

do i=1,n
write(*,900)(mat(i,j),j=1,n)
enddo

mat_a = mat
CALL dgesvd('A','A',n,n,mat_a,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO)
 
write(*,*) 'Singular Value'
write(*,900) S

write(*,*) 'U'
write(*,900) U

write(*,*) 'transposeV'
write(*,900) VT

write(*,*) 'EXIT'
write(*,900) INFO

  do j = 1, n
    do i = 1, n
      if (i == j) then
        W(i,j) = S(i)
      else
        W(i,j) = 0.0_r8
      end if
    end do
  end do

  call matmul(n, U, W, TMP1)
  call matmul(n, TMP1, VT, TMP2)

write(*,900) TMP2

contains

  subroutine matmul(n, a, b, c)

    integer, intent(in) :: n
    real(r8), dimension(n,n), intent(in) :: a, b
    real(r8), dimension(n,n), intent(out) :: c

    integer :: i, j, k

!   sum_k a_{i,k) b_{k,j} = c_{i,j}

    do j = 1, n
      do i = 1, n
        c(i,j) = 0.0_r8
        do k = 1, n
          c(i,j) = c(i,j) + a(i,k)*b(k,j) 
        end do
      end do
    end do

    return
  end subroutine matmul

End program matrix
