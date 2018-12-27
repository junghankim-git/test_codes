!-------------------------------------------------------------------------------
   module read_hybrid_file
!-------------------------------------------------------------------------------
!
!  abstract :
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-15  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   public :: read_coefficient
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_coefficient(filename, nlevs, coef_a, coef_b, eta, lecmwf)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),             intent(in   ) :: filename
   integer(i4),                  intent(in   ) :: nlevs
   real(r8), dimension(nlevs+1), intent(inout) :: eta, coef_a, coef_b
   logical(l4),                  intent(in   ) :: lecmwf
! local variables
   integer(i4) :: k
!
   open(unit=77,file=filename,status='unknown',form='formatted')
!
   do k = 1,nlevs+1
     read(77,*) coef_a(k),coef_b(k)
   enddo
!
   if (lecmwf) then
     do k = 1,nlevs+1
       coef_a(k) = coef_a(k)/100000.0_r8 ! pascal to eta
     enddo
   endif
!
   do k = 1,nlevs+1
     eta(k) = coef_a(k)+coef_b(k)
   enddo
!
   close(77)
!
   return
   end subroutine read_coefficient
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module read_hybrid_file
!-------------------------------------------------------------------------------
