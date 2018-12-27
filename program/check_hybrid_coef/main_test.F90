!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
   use kinds
   use read_hybrid_file, only : read_coefficient
!-------------------------------------------------------------------------------
   implicit none
! local variables
   integer(i4), parameter :: nlevs = 137
   real(r8), dimension(nlevs+1) :: coefa, coefb, eta
   character(len=10) :: levstr
   character(len=100) :: basename
   character(len=100) :: coef_file
   integer(i4) :: k
!
   print*, '######## start read coefficient file ########'
   print*, ' '
!
   basename = '/home/jhkim/work/program/hybrid_coord/coefficient/ecmwf/ecmwf_half_'
   write(levstr, '(i4.4) ') nlevs
   print*, levstr
   coef_file = trim(basename)//trim(levstr)//'nlevs.ascii'
   print*, coef_file
!
   call read_coefficient(coef_file, nlevs, coefa, coefb, eta, .true.)
!
   do k = 1, nlevs+1
     print*, coefa(k), coefb(k), eta(k)
   enddo
!
   end program test
