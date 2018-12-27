!-------------------------------------------------------------------------------
   module statistics
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    2015-09-24  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds,  only: i4, l4, r4, r8
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   interface get_lnorm
     module procedure get_lnorm_1d
     module procedure get_lnorm_2d
   end interface get_lnorm
!
   interface get_lerror
     module procedure get_lerror_1d
     module procedure get_lerror_2d
   end interface get_lerror
!
   public :: get_lnorm, get_lerror
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_lnorm_1d(nn, var, pnorm) result(norm)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),             intent(in   ) :: nn
   real(r8), dimension(nn), intent(in   ) :: var
   character(len=*),        intent(in   ) :: pnorm
   real(r8) :: norm
! local variables
   integer(i4) :: i
!
   norm = 0.0_r8
!
   if (trim(pnorm)=='L1') then
     do i = 1,nn
       norm = norm+dabs(var(i))
     enddo
   elseif (trim(pnorm)=='L2') then
     do i = 1,nn
       norm = norm+var(i)**2.0_r8
     enddo
     norm = dsqrt(norm)
   elseif (trim(pnorm)=='Linf') then
     norm = maxval(dabs(var))
   else
     !call printlog('check pnorm... '//pnorm,-2,2)
     print *,'check pnorm...',pnorm
     stop
   endif
!
   end function get_lnorm_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_lnorm_2d(nn, mm, var, pnorm) result(norm)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                intent(in   ) :: nn, mm
   real(r8), dimension(nn,mm), intent(in   ) :: var
   character(len=*),           intent(in   ) :: pnorm
   real(r8) :: norm
! local variables
   integer(i4) :: i, j
!
   norm = 0.0_r8
!
   if (trim(pnorm)=='L1') then
     do j = 1,mm
     do i = 1,nn
       norm = norm+dabs(var(i,j))
     enddo
     enddo
   elseif (trim(pnorm)=='L2') then
     do j = 1,mm
     do i = 1,nn
       norm = norm+var(i,j)**2.0_r8
     enddo
     enddo
     norm = dsqrt(norm)
   elseif (trim(pnorm)=='Linf') then
     norm = maxval(dabs(var))
   else
     !call printlog('Check pnorm... '//pnorm,-2,2)
     print *,'check pnorm...',pnorm
     stop
   endif
!
   end function get_lnorm_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_lerror_1d(nn, var1, var2, pnorm) result(error)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),             intent(in   ) :: nn
   real(r8), dimension(nn), intent(in   ) :: var1, var2
   character(len=*),        intent(in   ) :: pnorm
   real(r8) :: error
! local variables
   integer(i4) :: i
   real(r8), dimension(nn) :: diff
   real(r8) :: basenorm, diffnorm
!
   error = 0.0_r8
!
   if (trim(pnorm)=='L1'.or.trim(pnorm)=='L2'.or.trim(pnorm)=='Linf') then
     do i = 1,nn
       diff(i) = var2(i)-var1(i)
     enddo
     basenorm = get_lnorm(nn,var1,pnorm)
     diffnorm = get_lnorm(nn,diff,pnorm)
     if (basenorm==0) then
       error = 0.0_r8
     else
       error = diffnorm/basenorm
     endif
   else
     !call printlog('Check pnorm... '//pnorm,-2,2)
     print *,'check pnorm...',pnorm
     stop
   endif
!
   end function get_lerror_1d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function get_lerror_2d(nn, mm, var1, var2, pnorm) result(error)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                intent(in   ) :: nn, mm
   real(r8), dimension(nn,mm), intent(in   ) :: var1, var2
   character(len=*),           intent(in   ) :: pnorm
   real(r8) :: error
! local variables
   integer(i4) :: i, j
   real(r8), dimension(nn, mm) :: diff
   real(r8) :: basenorm, diffnorm
!
   error = 0.0_r8
!
   if (trim(pnorm)=='L1'.or.trim(pnorm)=='L2'.or.trim(pnorm)=='Linf') then
     do j = 1,mm
     do i = 1,nn
       diff(i,j) = var2(i,j)-var1(i,j)
     enddo
     enddo
     basenorm = get_lnorm(nn,mm,var1,pnorm)
     diffnorm = get_lnorm(nn,mm,diff,pnorm)
     if (basenorm==0) then
       error = 0.0_r8
     else
       error = diffnorm/basenorm
     endif
   else
     !call printlog('Check pnorm... '//pnorm,-2,2)
     print *,'check pnorm...',pnorm
     stop
   endif
!
   end function get_lerror_2d
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module statistics
!-------------------------------------------------------------------------------
