!-------------------------------------------------------------------------------
   module calendar
!-------------------------------------------------------------------------------
!
!  abstract : a calendar module for NWP(Numerical Weather Prediction)
!
!  history log :
!    2013-01-02   junghan kim     first written
!    2014-04-11   sang-yoon jun   modify calendar_inc
!    2016-09-27   junghan kim     code refactoring
!
!-------------------------------------------------------------------------------
!
   implicit none
!
!
   private
!
! kinds
!
   integer, parameter :: i4 = 4, l4 = 4
!
! common resources
!
   integer(i4), parameter, public :: calkind_none      = 0
   integer(i4), parameter, public :: calkind_gregorian = 1
   integer(i4), parameter, public :: calkind_julian    = 2
   integer(i4), parameter, public :: calkind_noleap    = 3
   integer(i4), parameter, public :: calkind_360day    = 4
!
! week
!
   integer(i4), parameter, public :: sun = 0
   integer(i4), parameter, public :: mon = 1
   integer(i4), parameter, public :: tue = 2
   integer(i4), parameter, public :: wed = 3
   integer(i4), parameter, public :: thu = 4
   integer(i4), parameter, public :: fri = 5
   integer(i4), parameter, public :: sat = 6
!
   character(len=3), dimension(0:6), public :: nameweek
   data nameweek/'SUN','MON','TUE','WED','THU','FRI','SAT'/
!
! month
!
   integer(i4), parameter, public :: jan = 1
   integer(i4), parameter, public :: feb = 2
   integer(i4), parameter, public :: mar = 3
   integer(i4), parameter, public :: apr = 4
   integer(i4), parameter, public :: may = 5
   integer(i4), parameter, public :: jun = 6
   integer(i4), parameter, public :: jul = 7
   integer(i4), parameter, public :: aug = 8
   integer(i4), parameter, public :: sep = 9
   integer(i4), parameter, public :: oct = 10
   integer(i4), parameter, public :: nov = 11
   integer(i4), parameter, public :: dec = 12
!
   character(len=3), dimension(12), public :: namemonth
   data namemonth /'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG', 'SEP','OCT','NOV','DEC'/
!
! initailization
!
   logical(l4), public :: is_initialized = .false.
!
! types
!
   type calendar_t
     integer(i4) :: kindcal
     integer(i4) :: year
     integer(i4) :: month
     integer(i4) :: day
     integer(i4) :: monthsperyear = 12
     integer(i4), dimension(12) :: dayspermonth
     integer(i4) :: daysperyear
     integer(i4) :: secondsperday = 24*60*60
     integer(i4) :: secondsperyear
   end type
!
   public :: calendar_t
   public :: calendar_set, calendar_get, calendar_inc, isleapyear, daysinmonth
   public :: daysinyear, dayofweek, date2numdays, numdays2date
   public :: calendar_convert, calendar_sync
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine calendar_set(calendar, kindcal, year, month, day)
!-------------------------------------------------------------------------------
   implicit none
!
   type(calendar_t), intent(inout) :: calendar
   integer(i4),      intent(in   ) :: kindcal
   integer(i4),      intent(in   ) :: year, month, day
!
   if ((kindcal>calkind_360day).or.(kindcal<=calkind_none)) then
     print *, 'check kindcal in calendar_set...'
     stop
   endif
!
   calendar%kindcal = kindcal
   if (isleapyear(kindcal, year)) then
     calendar%daysperyear = 366
     calendar%dayspermonth(:) =   (/31,29,31,30,31,30,31,31,30,31,30,31/)
   else
     calendar%daysperyear = 365
     calendar%dayspermonth(:) =   (/31,28,31,30,31,30,31,31,30,31,30,31/)
     if (kindcal==calkind_360day) then
       calendar%daysperyear = 360
       calendar%dayspermonth(:) = (/30,30,30,30,30,30,30,30,30,30,30,30/)
     endif
   endif
!
   calendar%year  = year
   calendar%month = month
   calendar%day   = day
!
   calendar%secondsperyear = sum(calendar%dayspermonth(:))*calendar%secondsperday
!
   is_initialized = .true.
!
   return
   end subroutine calendar_set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine calendar_get(calendar, year, month, day)
!-------------------------------------------------------------------------------
   implicit none
!
   type(calendar_t), intent(in   ) :: calendar
   integer(i4),      intent(  out) :: year, month, day
!
   year  = calendar%year
   month = calendar%month
   day   = calendar%day
!
   return
   end subroutine calendar_get
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function calendar_create(kindcal, year, month, day)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: kindcal
   integer(i4), intent(in   ) :: year, month, day
   type(calendar_t)           :: calendar_create
!
   calendar_create%kindcal = kindcal
   calendar_create%year    = year
   calendar_create%month   = month
   calendar_create%day     = day
!
   end function calendar_create
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine calendar_inc(calendar, year, month, day)
!-------------------------------------------------------------------------------
   implicit none
!
   type(calendar_t), intent(inout) :: calendar
   integer(i4),      intent(in   ) :: year, month, day
!
! local variables
!
   integer(i4) :: l_daysinmonth
!
   calendar%year  = calendar%year+year
   calendar%month = calendar%month+month
   calendar%day   = calendar%day+day
!
   if (calendar%month>12) then
     calendar%year  = calendar%year+(calendar%month/12)
     calendar%month = mod(calendar%month,12)
     call calendar_set(calendar,calendar%kindcal,calendar%year,calendar%month,calendar%day)
   endif
!
   l_daysinmonth = daysinmonth(calendar%kindcal,calendar%year,calendar%month)
!
   do while(calendar%day>l_daysinmonth)
     calendar%day   = calendar%day-l_daysinmonth
     calendar%month = calendar%month+1
     if (calendar%month>12) then
       calendar%year  = calendar%year+(calendar%month/12)
       calendar%month = mod(calendar%month,12)
       call calendar_set(calendar,calendar%kindcal,calendar%year,calendar%month,calendar%day)
     endif
     l_daysinmonth = daysinmonth(calendar%kindcal,calendar%year,calendar%month)
   enddo
!
   return
   end subroutine calendar_inc
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function isleapyear(kindcal, year)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: kindcal
   integer(i4), intent(in   ) :: year
   logical(l4)                :: isleapyear
!
   isleapyear = .false.
!
   if (kindcal==calkind_gregorian) then
     if (mod(year,4)==0) then
       isleapyear = .true.
     elseif(mod(year,100)==0) then
       isleapyear = .false.
     elseif(mod(year,400)==0) then
       isleapyear = .true.
     endif
   elseif(kindcal==calkind_julian) then
     if (mod(year,4)==0) then
       isleapyear = .true.
     endif
   elseif(kindcal==calkind_noleap) then
     isleapyear = .false.
   elseif(kindcal==calkind_360day) then
     isleapyear = .false.
   else
     print *, 'check kindcal in isleapyear...'
     stop
   endif
!
   end function isleapyear
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function daysinmonth(kindcal, year, month)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: kindcal
   integer(i4), intent(in   ) :: year
   integer(i4), intent(in   ) :: month
!
! local variables
!
   integer(i4) :: daysinmonth
   integer(i4) :: l_imonth, l_days
   integer(i4) :: l_dayspermonth(12)
!
   l_dayspermonth(:) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
!
   if (kindcal==calkind_gregorian) then
     daysinmonth = l_dayspermonth(month)
     if ((month==2).and.isleapyear(kindcal, year)) then
       daysinmonth = daysinmonth+1
     endif
   elseif(kindcal==calkind_julian) then
     daysinmonth = l_dayspermonth(month)
     if ((month==2).and.isleapyear(kindcal, year)) then
       daysinmonth = daysinmonth+1
     endif
   elseif(kindcal==calkind_noleap) then
     daysinmonth = l_dayspermonth(month)
   elseif(kindcal==calkind_360day) then
     daysinmonth = 30
   else
     print *, 'check kindcal in daysinmonth...'
     stop
   endif
!
   return
   end function daysinmonth
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function daysinyear(kindcal, year, month, day)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: kindcal
   integer(i4),           intent(in   ) :: year
   integer(i4), optional, intent(in   ) :: month
   integer(i4), optional, intent(in   ) :: day
!
! local variables
!
   integer(i4) :: daysinyear
   integer(i4) :: l_kindcal
   integer(i4) :: l_imonth, l_tmonth
!
   daysinyear = 0
!
   if (present(month)) then
     l_tmonth = month
   else
     l_tmonth = 13
   endif
!
   do l_imonth = 1,l_tmonth-1
     daysinyear = daysinyear+daysinmonth(kindcal,year,l_imonth)
   enddo
!
   if (present(month).and.present(day)) then
     daysinyear = daysinyear+day
   else
     print *,'check kindcal in daysinyear...'
     stop
   endif
!
   return
   end function daysinyear
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine numdays2date(kindcal, numdays, year, month, day)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: kindcal
   integer(i4),           intent(in   ) :: numdays
   integer(i4),           intent(  out) :: year
   integer(i4), optional, intent(  out) :: month
   integer(i4), optional, intent(  out) :: day
!
! local variables
!
   integer(i4) :: l_year, l_numleapyear
   integer(i4) :: l_dayspermonth(12)
   integer(i4) :: l_month, l_day
!
   l_dayspermonth(:) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
!
   if (kindcal==calkind_gregorian) then
     year = numdays/366
     l_numleapyear = year/4
     l_numleapyear = l_numleapyear-year/100
     l_numleapyear = l_numleapyear+year/400
     l_day = numdays-(year-l_numleapyear)*365-l_numleapyear*366
     do while(l_day>daysinyear(kindcal,year+1))
       l_day = l_day-daysinyear(kindcal,year+1)
       year  = year+1
     enddo
     year = year+1
     l_month = 1
     do while(l_day>daysinmonth(kindcal,year,l_month))
       l_day = l_day-daysinmonth(kindcal,year,l_month)
       l_month = l_month+1
     enddo
     if (present(month)) month = l_month
     if (present(day)) day = l_day
   elseif(kindcal==calkind_julian) then
     year = numdays/366
     l_numleapyear =(year)/4
     l_day = numdays -(year-l_numleapyear)*365-l_numleapyear*366
     do while(l_day>daysinyear(kindcal,year+1))
       l_day = l_day-daysinyear(kindcal,year+1)
       year = year+1
     enddo
     year = year+1
     l_month = 1
     do while(l_day>daysinmonth(kindcal,year,l_month))
       l_day = l_day-daysinmonth(kindcal,year,l_month)
       l_month = l_month+1
     enddo
     if (present(month)) month = l_month
     if (present(day)) day = l_day
   elseif(kindcal==calkind_noleap) then
     year = numdays/365+1
     l_day = numdays-year*365
     l_month = 1
     do while(l_day>l_dayspermonth(l_month))
       l_day = l_day-l_dayspermonth(l_month)
       l_month = l_month+1
     enddo
     if (present(month)) month = l_month
     if (present(day)) day = l_day
   elseif(kindcal==calkind_360day) then
     year = numdays/360+1
     l_month = mod(numdays,360)/30+1
     l_day = mod(mod(numdays,360),30)+1
     if (present(month)) month = l_month
     if (present(day)) day = l_day
   else
     print *,'check kindcal in numdays2date...'
     stop
   endif
!
   return
   end subroutine numdays2date
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine date2numdays(kindcal, year, month, day, numdays)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: kindcal
   integer(i4),           intent(in   ) :: year
   integer(i4), optional, intent(in   ) :: month
   integer(i4), optional, intent(in   ) :: day
   integer(i4),           intent(  out) :: numdays
!
! local variables
!
   integer(i4) :: l_year, l_numleapyear
!
   numdays = 0
   l_year = year-1
!
   if (kindcal==calkind_gregorian) then
     l_numleapyear =(l_year)/4
     l_numleapyear = l_numleapyear-(l_year)/100
     l_numleapyear = l_numleapyear+(l_year)/400
     numdays = numdays+(l_year-l_numleapyear)*365 ! days of regular year
     numdays = numdays+(l_numleapyear)*366 ! days of leap year
   elseif(kindcal==calkind_julian) then
     l_numleapyear =(l_year)/4
     numdays = numdays+(l_year-l_numleapyear)*365 ! days of regular year
     numdays = numdays+(l_numleapyear)*366 ! days of leap year
   elseif(kindcal==calkind_noleap) then
     numdays = numdays+l_year*365 ! days of regular year
! daysinyear below means days of this year
     numdays = numdays+daysinyear(kindcal, year, month, day)
   elseif(kindcal==calkind_360day) then
     numdays = numdays+l_year*360 ! days of regular year
! daysinyear below means days of this year
     numdays = numdays+daysinyear(kindcal, year, month, day)
   else
     print *, 'check kindcal in date2numdays...'
     stop
   endif
!
   if (present(month).and.present(day)) then
! daysinyear below means days of this year
     numdays = numdays+daysinyear(kindcal, year, month, day)
   elseif(present(month).and..not.present(day)) then
! daysinyear below means days of this year
     numdays = numdays+daysinyear(kindcal, year, month)
   elseif(.not.present(month).and..not.present(day)) then
! daysinyear below means days of this year
     numdays = numdays+daysinyear(kindcal, year)
   else
     print *, 'check date2numdays...'
     stop
   endif
!
   return
   end subroutine date2numdays
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function dayofweek(kindcal, year, month, day)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: kindcal
   integer(i4), intent(in   ) :: year
   integer(i4), intent(in   ) :: month
   integer(i4), intent(in   ) :: day
!
! local variables
!
   integer(i4) :: dayofweek
   integer(i4) :: numdays
!
   numdays = 0
!
   call date2numdays(kindcal,year,month,day,numdays)
!
   if (kindcal==calkind_gregorian) then
     dayofweek = sun ! 01-jan-0001 is monday,so base is sunday
   elseif(kindcal==calkind_julian) then
     dayofweek = fri ! 01-jan-0001 is saturday,so base is friday
   elseif(kindcal==calkind_noleap) then
     dayofweek = sun ! 01-jan-0001 is monday,so base is sunday
   elseif(kindcal==calkind_360day) then
     dayofweek = sun ! 01-jan-0001 is monday,so base is sunday
   else
     print *,'check kindcal in dayofweek...'
     stop
   endif
!
   dayofweek = dayofweek+numdays
   dayofweek = mod(dayofweek,7)
!
   return
   end function dayofweek
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine calendar_convert(kindcal1, year1, month1, day1,                  &
                                                 kindcal2, year2, month2, day2)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: kindcal1
   integer(i4), intent(in   ) :: year1
   integer(i4), intent(in   ) :: month1
   integer(i4), intent(in   ) :: day1
   integer(i4), intent(in   ) :: kindcal2
   integer(i4), intent(  out) :: year2
   integer(i4), intent(  out) :: month2
   integer(i4), intent(  out) :: day2
!
! local variables
!
   integer(i4) :: l_numdays
!
   call date2numdays(kindcal1,year1,month1,day1,l_numdays)
!
   if ((kindcal1==calkind_julian).and.(kindcal2/=calkind_julian)) then
     l_numdays = l_numdays-2
   elseif((kindcal1/=calkind_julian).and.(kindcal2==calkind_julian)) then
     l_numdays = l_numdays+2
   endif
!
   call numdays2date(kindcal2,l_numdays,year2,month2,day2)
!
   return
   end subroutine calendar_convert
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine calendar_sync(calendar, datetime)
!-------------------------------------------------------------------------------
   implicit none
!
   type(calendar_t),                    intent(  out) :: calendar
   integer(i4), dimension(3), optional, intent(in   ) :: datetime
!
! local variables
!
   character(len=10), dimension(3) :: l_charval
   integer(i4), dimension(8)       :: l_datetime
!
   if (present(datetime)) then
     calendar%year  = datetime(1)
     calendar%month = datetime(2)
     calendar%day   = datetime(3)
   else
     call date_and_time(l_charval(1),l_charval(2),l_charval(3),l_datetime)
     calendar%year  = l_datetime(1)
     calendar%month = l_datetime(2)
     calendar%day   = l_datetime(3)
   endif
!
   return
   end subroutine calendar_sync
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module calendar
!-------------------------------------------------------------------------------
