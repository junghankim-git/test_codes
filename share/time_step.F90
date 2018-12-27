!-------------------------------------------------------------------------------
   module time_step
!-------------------------------------------------------------------------------
!
!  abstract : control time step
!
!  history log :
!    2017-05-25   junghan kim     first written
!
!-------------------------------------------------------------------------------
   use kinds, only: i4, l4, r4, r8
!
   implicit none
!
   private
!
! time step
!
   integer(i4), parameter, public :: ntls = 3
!
! definition for time_step object
!
   type :: time_step_t
     integer(i4) :: nm, n0, np
     integer(i4) :: istep, nsteps, osteps
     logical(l4) :: is_max_step
   end type
!
! internal time step object
!
   public :: time_step_t, time_step_set, time_step_update
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine time_step_set(tl, nsteps) ! not used!!!
!-------------------------------------------------------------------------------
   implicit none
!
   type(time_step_t), intent(inout) :: tl
   integer(i4),       intent(in   ) :: nsteps
!
   tl%nm     = 1
   tl%n0     = 2
   tl%np     = 3
   tl%istep  = 0
   tl%nsteps = nsteps
   tl%is_max_step = .false.
!
   return
   end subroutine time_step_set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine time_step_update(tl, update_step)
!-------------------------------------------------------------------------------
   implicit none
!
   type(time_step_t),           intent(inout) :: tl
   logical(l4)      , optional, intent(in   ) :: update_step
! local variables
   integer(i4) :: tmp
   logical(l4) :: update
!
   if (present(update_step)) then
     update = update_step
   else
     update = .true.
   endif
!
   tmp   = tl%np
   tl%np = tl%nm
   tl%nm = tl%n0
   tl%n0 = tmp
!
   if (update) then
     tl%istep = tl%istep+1
     if (tl%istep.ge.tl%nsteps) tl%is_max_step = .true.
   endif
!
   return
   end subroutine time_step_update
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module time_step
!-------------------------------------------------------------------------------
