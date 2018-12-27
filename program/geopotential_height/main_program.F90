!-------------------------------------------------------------------------------
   program genvcoord
!-------------------------------------------------------------------------------
!
!  abstract : main program of the coupled system (atm, wav, ocn, ... + cpl)
!
!  history log :
!    2014-04-15   junghan kim     initial setup
!
!-------------------------------------------------------------------------------
   use kinds       , only: i4, r8, l4
   use geopotential, only: initialize, finalize, get_height_p, get_height_ex, get_height_et, p_i, p_m
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4) :: nlevs
   character(len=512) :: file_i, file_m
   real(r8), dimension(:), allocatable :: hi_p, hm_p, hi_ex, hm_ex, hi_et, hm_et
   integer(i4) :: k
!
   allocate(hi_p(nlevs+1),hi_ex(nlevs+1),hi_et(nlevs+1))
   allocate(hm_p(nlevs),hm_ex(nlevs),hm_et(nlevs))
!
   nlevs = 50
   file_i = '/home/jhkim/TestBed/KIM/3.0.micros/3.0.10/04.dev/tools/Inputs/vcoord/half_50_0.3.dat'
   file_m = '/home/jhkim/TestBed/KIM/3.0.micros/3.0.10/04.dev/tools/Inputs/vcoord/full_50_0.3.dat'
   call initialize(nlevs,trim(file_i),trim(file_m))
   call get_height_p(hi_p,hm_p)
   call get_height_ex(hi_ex,hm_ex)
   call get_height_et(hi_et,hm_et)
!
   do k = 1,nlevs+1
     print '(i2,x,f18.11,2x,f18.11,2x,f18.11,2x,f18.11)',k,p_i(k),hi_p(k),hi_ex(k),hi_p(k)-hi_ex(k) !hi_et(k)
     !print *,k,p_i(k),hi_p(k),hi_ex(k),hi_et(k)
   enddo
   print *, ' '
   do k = 1,nlevs
     print '(i2,x,f18.11,2x,f18.11,2x,f18.11,2x,f18.11)',k,p_m(k),hm_p(k),hm_ex(k),hm_p(k)-hm_ex(k)
     !print *,k,p_m(k),hm_p(k),hm_ex(k),hm_et(k)
   enddo
!
   call finalize()
   deallocate(hi_p,hi_ex,hi_et)
   deallocate(hm_p,hm_ex,hm_et)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set_exp_case(jexp)
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: jexp
!
   return
   end subroutine set_exp_case
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end
!-------------------------------------------------------------------------------
