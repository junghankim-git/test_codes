!-------------------------------------------------------------------------------
   module test
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  your   name    code comment
!
!  structure:
!
!-------------------------------------------------------------------------------
   use mymodule1, only: abc, def, fjdk, &
                                                                      fjdsakl, &
                                                                         fjdoia
   use mymodule2, only: n
   implicit none
!
   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine my_sub1(a, b, c)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, dimension(:), intent(in   ) :: a
   real, dimension(:), intent(  out) :: b
   real, dimension(:), intent(inout) :: c
   integer :: i, j, k, var
!***********************************************
!
! my main algorithm
!
! calculate y = ax + b
! z = \frac {} {} \times psi + delta, ...
!
!***********************************************
!
   do i = 1, n
     b(i) = real(a(i))
   enddo
!
   do k = 1, n
     do j = 1, n
       do i = 1, n
         c(i) = c(i)+b(i)+some_calculation ...some_calculation... some_calculation...some_calculation...+a(i)
       enddo
       enddo
   enddo
!
   end subroutine
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function my_fun(a, b, c)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, dimension(:), intent(in   ) :: a
   real, dimension(:), intent(  out) :: b
   real, dimension(:), intent(inout) :: c
! my main algorithm
!
   if (a==b) then
     if (b==c) then
! some code
     endif
   endif
!
   end subroutine
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module test
!-------------------------------------------------------------------------------
