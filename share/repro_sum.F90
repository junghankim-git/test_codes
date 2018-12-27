!-------------------------------------------------------------------------------
   module repro_sum
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds
!-------------------------------------------------------------------------------
   implicit none
!
   public :: scs, ddpdd, ddpddf
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine scs(add_a, err, add_b)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(inout) :: add_a, err
   real(r8), intent(in   ) :: add_b
! local variables
   real(r8) :: tsum, new_err, tmp
!
   tmp = err+add_b
   tsum = add_a+tmp
   new_err = tmp-(tsum-add_a)
!
   add_a = tsum
   err = new_err
!
   return
   end subroutine scs
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ddpdd(add_a, err, add_b)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(inout) :: add_a, err
   real(r8), intent(in   ) :: add_b
! local variables
   real(r8) :: tsum, new_err
!
   tsum = add_a+add_b
   new_err =(add_b-(tsum-add_a))+(add_a-(tsum-(tsum-add_a)))
   call scs(tsum,err,new_err)
!
   add_a = tsum
!
   return
   end subroutine ddpdd
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function ddpddf(add_a, err, add_b)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(inout) :: err
   real(r8), intent(in   ) :: add_a, add_b
   real(r8) :: ddpddf
! local variables
   real(r8) :: tsum, new_err
!
   tsum = add_a+add_b
   new_err =(add_b-(tsum-add_a))+(add_a-(tsum-(tsum-add_a)))
   call scs(tsum,err,new_err)
!
   ddpddf = tsum
!
   end function ddpddf
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module repro_sum
!-------------------------------------------------------------------------------
