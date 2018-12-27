!-------------------------------------------------------------------------------
   module logger
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
   use kinds, only : i4, r4, r8, l4
!   use parallel, only : kim_par
!
   private
!
   ! maximium message size
   integer(i4), parameter, public :: max_str = 4096
   integer(i4), parameter, public :: max_tostr = 15
   integer(i4), parameter, public :: max_lev = 5
   ! error levels
   integer(i4), parameter, public :: err_normal = 0
   integer(i4), parameter, public :: err_warning =-1
   integer(i4), parameter, public :: err_fatal =-2
!
   public :: printlog, tostring
!
   interface tostring
     module procedure tostring_i4
     module procedure tostring_l4
     module procedure tostring_r4
     module procedure tostring_r8
   end interface tostring
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine printlog(string, elev, llev, isall)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*), intent(in   ) :: string
   integer(i4),      intent(in   ), optional :: elev, llev
   logical(l4),      intent(in   ), optional :: isall
! local variables
   integer(i4) :: l_l_lev, l_e_lev, ierr
   logical(l4) :: l_isall
   character(len=max_str) :: message
   character(len=8) :: header
!
   if (present(elev)) then
     l_e_lev = elev
   else
     l_e_lev =-99
   endif
!
   if (present(llev)) then
     l_l_lev = llev
   else
     l_l_lev = 5
   endif
!
   if (present(isall)) then
     l_isall = isall
   else
     l_isall = .false.
   endif
!
   if (l_e_lev==err_normal) then
     write(header,'(a8) ') '[ log] : '
   elseif (l_e_lev==err_warning) then
     write(header,'(a8) ') '[ warn] : '
   elseif (l_e_lev==err_fatal) then
     write(header,'(a8) ') '[fatal] : '
   else
     write(header,'(a8) ') ' '
   endif
!
   if (l_l_lev==1) then
     write(message,'(a8,a1,x,a) ') header,'#',string
   elseif (l_l_lev==2) then
     write(message,'(a8,a2,x,a) ') header,'*',string
   elseif (l_l_lev==3) then
     write(message,'(a8,a3,x,a) ') header,'-',string
   elseif (l_l_lev==4) then
     write(message,'(a8,a4,x,a) ') header,' .',string
   else
     write(message,'(a8,a5,x,a) ') header,' ',string
   endif
!
!   if (l_isall) then
!     write(*,*) trim(message),kim_par%iproc
!   else
!     if (kim_par%ismaster) then
       write(*,*) trim(message)
!     endif
!   endif
!
   if (l_e_lev==err_fatal) then
     stop
   endif
!
   return
   end subroutine printlog
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function tostring_i4(value) result(outstring)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: value
   character(len=16) :: outstring
!
   if (abs(value)>= 0.and.abs(value)<10) then
     write(outstring,'(i2)') value
   elseif (abs(value)>= 10.and.abs(value)<100) then
     write(outstring,'(i3)') value
   elseif (abs(value)>= 100.and.abs(value)<1000) then
     write(outstring,'(i4)') value
   elseif (abs(value)>= 1000.and.abs(value)<10000) then
     write(outstring,'(i5)') value
   elseif (abs(value)>= 10000.and.abs(value)<100000) then
     write(outstring,'(i6)') value
   elseif (abs(value)>= 100000.and.abs(value)<1000000) then
     write(outstring,'(i7)') value
   elseif (abs(value)>= 1000000.and.abs(value)<10000000) then
     write(outstring,'(i8)') value
   elseif (abs(value)>= 10000000.and.abs(value)<100000000) then
     write(outstring,'(i9)') value
   elseif (abs(value)>= 100000000.and.abs(value)<1000000000) then
     write(outstring,'(i10)') value
   elseif (abs(value)>= 1000000000.and.abs(value)<huge(value)) then
     write(outstring,'(i11)') value
   else
     write(outstring,'(i15)') value
   endif
!
   end function tostring_i4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function tostring_l4(value) result(outstring)
!-------------------------------------------------------------------------------
   implicit none
!
   logical(l4), intent(in   ) :: value
   character(len=8) :: outstring
!
   if (value) then
     write(outstring,'(a)') 'TRUE'
   else
     write(outstring,'(a)') 'FALSE'
   endif
!
   end function tostring_l4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function tostring_r4(value) result(outstring)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), intent(in   ) :: value
   character(len=8) :: outstring
!
   write(outstring,'(f8.5)') value
!
   end function tostring_r4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function tostring_r8(value) result(outstring)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: value
   character(len=14) :: outstring
!
   write(outstring,'(f14.11)') value
!
   end function tostring_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module logger
!-------------------------------------------------------------------------------
