!-------------------------------------------------------------------------------
   program test
!-------------------------------------------------------------------------------
   use kinds
   use space_filling_curve, only: do_factorize
!-------------------------------------------------------------------------------
   implicit none
   integer(i4) :: n, i
   logical(l4) :: isfac!, is3, is2
   integer(i4) :: nfac
   integer(i4), allocatable :: facs(:)
   character(len=1) :: is3, is2
   real(r8) :: a, res4, res3
!
   a = 10020.0_r8
   n = 1000
!
   do i = 1,n
     !is3 = .false.
     !is2 = .false.
     !if (mod(i,3).eq.0) is3 = .true.
     !if (mod(i,2).eq.0) is2 = .true.
     is3 = 'X'
     is2 = 'X'
     if (mod(i,3).eq.0) is3 = 'O'
     if (mod(i,2).eq.0) is2 = 'O'
     !
     isfac = do_factorize(i,nfac,facs)
     if (nfac.gt.0.and.isfac) deallocate(facs)
     if (isfac) then
       res4 = a/real(i*3,r8)
       res3 = a/real(i*2,r8)
       !print '(i4,a,2x,a,x,i5,x,f7.2,3x,a,x,i5,x,f7.2)',i,':',is2,i*3,res4,is3,i*2,res3
       print '(i4,3x,a,x,i5,3x,a,x,i5)',i,is2,i*3,is3,i*2
     endif
   enddo
!
   end program test
