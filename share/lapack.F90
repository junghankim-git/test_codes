!-------------------------------------------------------------------------------
   module lapack
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
   use kinds, only : i4, r4, r8
!-------------------------------------------------------------------------------
   implicit none
!
   public :: gesvd
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine gesvd(m_in, n_in, a_in, s_in, u_in, vt_in)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),                          intent(in   ) :: m_in, n_in
   real(r8), dimension(m_in, n_in),      intent(inout) :: a_in
   real(r8), dimension(m_in, m_in),      intent(inout) :: u_in
   real(r8), dimension(n_in, n_in),      intent(inout) :: vt_in
   real(r8), dimension(min(n_in, m_in)), intent(inout) :: s_in
! local variables
   integer(i4) :: lwork, info
   real(r8), dimension(5*m_in*n_in) :: work
!
   lwork = 5*m_in*n_in
!
   call dgesvd('A','A',m_in,n_in,a_in,m_in,s_in,u_in,m_in,vt_in,n_in,work,lwork,info)
   if (info.ne.0) then
     print*,'DGESVD : check arguments.... program will be stoped.'
     stop
   endif
!
   return
   end subroutine gesvd
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module lapack
