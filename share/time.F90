!-------------------------------------------------------------------------------
   module time
!-------------------------------------------------------------------------------
!
!  abstract : main program of the coupled system (atm, wav, ocn, ...+cpl)
!
!  history log :
!    2013-01-07   junghan kim     first written
!    2016-07-12   junghan kim     code refactoring
!    2016-07-25   junghan kim     add string2second, time2string
!    2016-09-29   junghan kim     refactoring (remove optional calendar, interval)
!
!-------------------------------------------------------------------------------
   use calendar, only : calendar_t, calendar_set, calendar_inc
   use calendar, only : dayofweek, daysinmonth, daysinyear, nameweek, namemonth
!
   implicit none
!
   private
!
   integer, parameter :: i4 = 4, l4 = 4, r4 = 4, r8 = 8
!
! definition for time object
!
   interface getdigits
     module procedure getdigits_i4
     module procedure getdigits_r4
     module procedure getdigits_r8
   end interface getdigits
!
   interface get_julianday
     module procedure get_julianday_r4
     module procedure get_julianday_r8
   end interface get_julianday
!
   interface floatsecond2integer
     module procedure floatsecond2integer_r4
     module procedure floatsecond2integer_r8
   end interface floatsecond2integer
!
   type timeinterval_t
     integer(i4) :: year
     integer(i4) :: month
     integer(i4) :: day
     integer(i4) :: hour
     integer(i4) :: minute
     integer(i4) :: second
     integer(i4) :: millisecond
     integer(i4) :: microsecond
   end type
!
   type time_t
     type(calendar_t)     :: cal
     type(timeinterval_t) :: interval
     integer(i4)          :: hour
     integer(i4)          :: minute
     integer(i4)          :: second
     integer(i4)          :: millisecond
     integer(i4)          :: microsecond
     logical(l4)          :: isperpetual
   end type
!
! time step
!
   integer(i4), parameter, public :: ntimelevels = 3
!
   integer(i4), public :: nsplit
   integer(i4), public :: nmaxstep
   integer(i4), public :: nstep
   integer(i4), public :: ndays
   real(r8),    public :: timestep ! dynamics time step must be set
   real(r8),    public :: timescale_phys ! physics time scale must be set
!
! definition for timelevel object
!
   type :: timelevel_t
     integer(i4) :: nm1
     integer(i4) :: n
     integer(i4) :: np1
     integer(i4) :: istep
     integer(i4) :: istepleap
   end type
!
! internal timeLevel object
!
   public :: time_t, timeinterval_t
   public :: time_set, time_get, time_create, time_sync
   public :: time_interval_set, time_interval_create
   public :: time_inc, time_print, floatsecond2integer
   public :: time_level_set, time_level_update, current_timestep
   public :: get_julianday, time_getstring
   public :: time2string, string2second
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine time_set(time, caltype, year, month, day, hour, minute, second,  &
                                          millisecond, microsecond,isperpetual)
!-------------------------------------------------------------------------------
   implicit none
!
   type(time_t),          intent(inout) :: time
   integer(i4),           intent(in   ) :: caltype
   integer(i4), optional, intent(in   ) :: year, month, day
   integer(i4), optional, intent(in   ) :: hour, minute, second
   integer(i4), optional, intent(in   ) :: millisecond, microsecond
   logical(l4), optional, intent(in   ) :: isperpetual
!
! local variables
!
   integer(i4) :: lyear, lmonth, lday
!
! calendar
!
   lyear  = 1
   lmonth = 1
   lday   = 1
   if (present(year))  lyear  = year
   if (present(month)) lmonth = month
   if (present(day))   lday   = day
   call calendar_set(time%cal, caltype, lyear, lmonth, lday)
!
! interval
!
#if 0
   time%interval%year        = 0
   time%interval%month       = 0
   time%interval%day         = 0
   time%interval%hour        = 0
   time%interval%minute      = 0
   time%interval%second      = 0
   time%interval%millisecond = 0
   time%interval%microsecond = 0
#endif
!
! time
!
   time%hour        = 0
   time%minute      = 0
   time%second      = 0
   time%millisecond = 0
   time%microsecond = 0
   if (present(hour))        time%hour        = hour
   if (present(minute))      time%minute      = minute
   if (present(second))      time%second      = second
   if (present(millisecond)) time%millisecond = millisecond
   if (present(microsecond)) time%microsecond = microsecond
!
! perpetual
!
   if (present(isperpetual)) then
     time%isperpetual = isperpetual
   else
     time%isperpetual = .false.
   endif
!
   return
   end subroutine time_set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function time_create(caltype, year, month, day, hour, minute, second,       &
                            millisecond, microsecond, isperpetual) result(time)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),           intent(in   ) :: caltype
   integer(i4), optional, intent(in   ) :: year, month, day
   integer(i4), optional, intent(in   ) :: hour, minute, second
   integer(i4), optional, intent(in   ) :: millisecond, microsecond
   logical(l4), optional, intent(in   ) :: isperpetual
   type(time_t)                         :: time
!
! local variables
!
   integer(i4) :: lyear, lmonth, lday
!
! calendar
!
   lyear  = 1
   lmonth = 1
   lday   = 1
   if (present(year))  lyear  = year
   if (present(month)) lmonth = month
   if (present(day))   lday   = day
   call calendar_set(time%cal, caltype, lyear, lmonth, lday)
!
! interval
!
   time%interval%year        = 0
   time%interval%month       = 0
   time%interval%day         = 0
   time%interval%hour        = 0
   time%interval%minute      = 0
   time%interval%second      = 0
   time%interval%millisecond = 0
   time%interval%microsecond = 0
!
! time
!
   time%hour        = 0
   time%minute      = 0
   time%second      = 0
   time%millisecond = 0
   time%microsecond = 0
   if (present(hour))        time%hour        = hour
   if (present(minute))      time%minute      = minute
   if (present(second))      time%second      = second
   if (present(millisecond)) time%millisecond = millisecond
   if (present(microsecond)) time%microsecond = microsecond
!
! perpetual
!
   if (present(isperpetual)) then
     time%isperpetual = isperpetual
   else
     time%isperpetual = .false.
   endif
!
   end function time_create
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine time_get(time, year, month, day, hour, minute, second)
!-------------------------------------------------------------------------------
   implicit none
!
   type(time_t),          intent(in   ) :: time
   integer(i4), optional, intent(  out) :: year, month, day
   integer(i4),           intent(  out) :: hour, minute, second
!
! local variables
!
   if (present(year))  year  = time%cal%year
   if (present(month)) month = time%cal%month
   if (present(day))   day   = time%cal%day
!
   hour   = time%hour
   minute = time%minute
   second = time%second
!
   return
   end subroutine time_get
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_julianday_r4(time, julianday)
!-------------------------------------------------------------------------------
   implicit none
!
   type(time_t), intent(in   ) :: time
   real(r4),     intent(  out) :: julianday
!
! local variables
!
   real(r4)    :: l_secspertotalsecsinday
   integer(i4) :: l_dayinyear, l_secsinday
   integer(i4) :: l_year, l_month, l_day, l_hour, l_minute, l_second
!
   l_year  = time%cal%year
   l_month = time%cal%month
   l_day   = time%cal%day
!
   l_hour   = time%hour
   l_minute = time%minute
   l_second = time%second
!
   l_secsinday = l_hour*60*60+l_minute*60+l_second
   l_dayinyear = daysinyear(time%cal%kindcal,l_year,l_month,l_day)
   l_secspertotalsecsinday = real(l_secsinday)/real(time%cal%secondsperday)
   julianday = real(l_dayinyear)+l_secspertotalsecsinday
!
   return
   end subroutine get_julianday_r4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine get_julianday_r8(time, julianday)
!-------------------------------------------------------------------------------
   implicit none
!
   type(time_t), intent(in   ) :: time
   real(r8),     intent(  out) :: julianday
!
! local variables
!
   real(r8)    :: l_secspertotalsecsinday
   integer(i4) :: l_dayinyear, l_secsinday
   integer(i4) :: l_year, l_month, l_day, l_hour, l_minute, l_second
!
   l_year  = time%cal%year
   l_month = time%cal%month
   l_day   = time%cal%day
!
   l_hour   = time%hour
   l_minute = time%minute
   l_second = time%second
!
   l_secsinday = l_hour*60*60+l_minute*60+l_second
   l_dayinyear = daysinyear(time%cal%kindcal,l_year,l_month,l_day)
   l_secspertotalsecsinday = real(l_secsinday,r8)/real(time%cal%secondsperday,r8)
   julianday = real(l_dayinyear,r8)+l_secspertotalsecsinday
!
   return
   end subroutine get_julianday_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine time_sync(time)
!-------------------------------------------------------------------------------
   implicit none
!
   type(time_t), intent(inout) :: time
!
! local variables
!
   character(len=10), dimension(3) :: l_charval
   integer(i4), dimension(8) :: l_datetime
!
   call date_and_time(l_charval(1),l_charval(2),l_charval(3),l_datetime)
!
   time%cal%year    = l_datetime(1)
   time%cal%month   = l_datetime(2)
   time%cal%day     = l_datetime(3)
   time%hour        = l_datetime(5)
   time%minute      = l_datetime(6)
   time%second      = l_datetime(7)
   time%millisecond = l_datetime(8)
   time%microsecond = 0
!
   return
   end subroutine time_sync
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine time_interval_set(timeinterval, year, month, day, hour, minute,    &
                              second, millisec, microsec)
!-------------------------------------------------------------------------------
   implicit none
!
   type(timeinterval_t),  intent(inout) :: timeinterval
   integer(i4), optional, intent(in   ) :: year, month, day
   integer(i4), optional, intent(in   ) :: hour, minute, second
   integer(i4), optional, intent(in   ) :: millisec, microsec
!
   timeinterval%year        = 0
   timeinterval%month       = 0
   timeinterval%day         = 0
   timeinterval%hour        = 0
   timeinterval%minute      = 0
   timeinterval%second      = 0
   timeinterval%millisecond = 0
   timeinterval%microsecond = 0
!
   if (present(year))     timeinterval%year        = year
   if (present(month))    timeinterval%month       = month
   if (present(day))      timeinterval%day         = day
   if (present(hour))     timeinterval%hour        = hour
   if (present(minute))   timeinterval%minute      = minute
   if (present(second))   timeinterval%second      = second
   if (present(millisec)) timeinterval%millisecond = millisec
   if (present(microsec)) timeinterval%microsecond = microsec
!
   return
   end subroutine time_interval_set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function time_interval_create(year, month, day, hour, minute, second)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), optional, intent(in   ) :: year, month, day
   integer(i4),           intent(in   ) :: hour, minute, second
   type(timeinterval_t)                 :: time_interval_create
!
   time_interval_create%year  = 0
   time_interval_create%month = 0
   time_interval_create%day   = 0
   if (present(year))  time_interval_create%year  = year
   if (present(month)) time_interval_create%month = month
   if (present(day))   time_interval_create%day   = day
   time_interval_create%hour   = hour
   time_interval_create%minute = minute
   time_interval_create%second = second
!
   return
   end function time_interval_create
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine time_inc(time)
!-------------------------------------------------------------------------------
   implicit none
!
   type(time_t), intent(inout) :: time
!
! local variables
!
   integer(i4) :: l_day, l_hour, l_minute
   integer(i4) :: l_second, l_millisecond, l_microsecond
!
   l_microsecond = time%microsecond+time%interval%microsecond
   l_millisecond = time%millisecond+time%interval%millisecond
   l_second      = time%second+time%interval%second
   l_minute      = time%minute+time%interval%minute
   l_hour        = time%hour+time%interval%hour
!
   l_millisecond    = l_millisecond+l_microsecond/1000
   time%microsecond = mod(l_microsecond,1000)
!
   l_second         = l_second+l_millisecond/1000
   time%millisecond = mod(l_millisecond,1000)
!
   l_minute    = l_minute+l_second/60
   time%second = mod(l_second,60)
!
   l_hour      = l_hour+l_minute/60
   time%minute = mod(l_minute,60)
!
   l_day     = l_hour/24
   time%hour = mod(l_hour,24)
!
   if (.not.time%isperpetual) then
     call calendar_inc(time%cal,year=time%interval%year,                       &
                         month=time%interval%month,day=time%interval%day+l_day)
   endif
!
   return
   end subroutine time_inc
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine time_getstring(time, currtime, noblank)
!-------------------------------------------------------------------------------
   implicit none
!
   type(time_t),          intent(in   ) :: time
   character(len=*),      intent(  out) :: currtime
   logical(l4), optional, intent(in   ) :: noblank
!
! local variables
!
   integer(i4) :: l_dayweek
   logical(l4) :: l_noblank
!
   if (present(noblank)) then
     l_noblank = noblank
   else
     l_noblank = .false.
   endif
   l_dayweek = dayofweek(time%cal%kindcal,time%cal%year,time%cal%month,        &
                                                                   time%cal%day)
!
   if (l_noblank) then
     write(currtime,101),nameweek(l_dayweek),'.',namemonth(time%cal%month),'.',&
    time%cal%day,'.',time%hour,':',time%minute,':',time%second,'.',time%cal%year
   else
     write(currtime,100),nameweek(l_dayweek),namemonth(time%cal%month),        &
            time%cal%day,time%hour,':',time%minute,':',time%second,time%cal%year
   endif
!
   100 format(a3,x,a3,x,i2.2,x,i2.2,a,i2.2,a,i2.2,x,i4.4)
   101 format(a3,a1,a3,a1,i2.2,a1,i2.2,a,i2.2,a,i2.2,a1,i4.4)
!
   return
   end subroutine time_getstring
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine time_print(time)
!-------------------------------------------------------------------------------
   implicit none
!
   type(time_t), intent(in   ) :: time
!
! local variables
!
   integer(i4) :: l_dayweek
   character(len=32) :: currtime
!
   call time_getstring(time,currtime)
   write(*,'(a)') trim(currtime)
!
   return
   end subroutine time_print
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine time_level_set(timelevel, nm1, n, np1, istep) ! not used!!!
!-------------------------------------------------------------------------------
   implicit none
!
   type(timelevel_t),     intent(inout) :: timelevel
   integer(i4), optional, intent(in   ) :: nm1, n, np1, istep
!
   timelevel%nm1       = 1
   timelevel%n         = 2
   timelevel%np1       = 3
   timelevel%istep     = 0
   timelevel%istepleap = 2
!
   if (present(nm1).and.present(n).and.present(np1).and.present(istep)) then
     timelevel%nm1 = nm1
     timelevel%n = n
     timelevel%np1 = np1
     timelevel%istep = istep
   elseif(present(nm1).or.present(n).or.present(np1).or.present(istep)) then
     print *, 'check arguments in time_level_set...'
     stop
   endif
!
   return
   end subroutine time_level_set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine time_level_update(timelevel, updatetype) ! not used
!-------------------------------------------------------------------------------
   implicit none
!
   type(timelevel_t), intent(inout) :: timelevel
   character(len=*),  intent(in   ) :: updatetype
!
! local variables
!
   integer(i4) :: l_tmp
!
   if (trim(adjustl(updatetype))=='leapfrog') then
     l_tmp         = timelevel%np1
     timelevel%np1 = timelevel%nm1
     timelevel%nm1 = timelevel%n
     timelevel%n   = l_tmp
   elseif(trim(adjustl(updatetype))=='forward') then
     l_tmp         = timelevel%np1
     timelevel%np1 = timelevel%n
     timelevel%n   = l_tmp
   else
     print *, 'check updatetype in time_level_update...'
     stop
   endif
!
   timelevel%istep = timelevel%istep+1
!
   return
   end subroutine time_level_update
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function current_timestep(timelevel)
!-------------------------------------------------------------------------------
   implicit none
!
   type(timelevel_t), intent(in   ) :: timelevel
!
! local variables
!
   real(r8)    :: current_timestep
   integer(i4) :: l_tmp
!
   current_timestep = dble(timelevel%istep)*timestep
!
   return
   end function current_timestep
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine floatsecond2integer_r4(second, sec, millisec, microsec)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4),    intent(in   ) :: second
   integer(i4), intent(  out) :: sec, millisec, microsec
!
! local variables
!
   real(r4) :: tmp
!
   sec = int(second)
   tmp = second-float(sec)
!
   millisec = int(tmp*1000.0)
   tmp      = tmp-float(millisec)/1000.0
!
   microsec = int(tmp*1000.0*1000.0)
   tmp      = tmp-float(microsec)/1000.0/1000.0
!
   return
   end subroutine floatsecond2integer_r4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine floatsecond2integer_r8(second, sec, millisec, microsec)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8),    intent(in   ) :: second
   integer(i4), intent(  out) :: sec, millisec, microsec
!
! local variables
!
   real(r8) :: tmp
!
   sec = int(second)
   tmp = second-float(sec)
!
   millisec = int(tmp*1000.0_r8)
   tmp      = tmp-float(millisec)/1000.0_r8
!
   microsec = int(tmp*1000.0_r8*1000.0_r8)
   tmp      = tmp-float(microsec)/1000.0_r8/1000.0_r8
!
   return
   end subroutine floatsecond2integer_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function string2second(string) result(second)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*), intent(in   ) :: string
   integer(i4)                     :: second
!
! local variables
!
   integer(i4), parameter        :: n_nums = 10, n_units = 4
   character, dimension(n_nums)  :: num_string
   character, dimension(n_units) :: unit_string
   integer, dimension(n_units)   :: pos_units
   integer, dimension(n_units)   :: times, factors
   integer(i4)                   :: nchars, i, j, is, ie, nn
   character(len=32)             :: nfmt
!
   num_string  = (/'0','1','2','3','4','5','6','7','8','9'/)
   unit_string = (/'d','h','m','s'/)
   times(:)    = 0 !(/day,hour,minute,second/)
   factors(:)  = (/24*60*60,60*60,60,1/)
   second      = 0
!
! get times
!
   pos_units = getposunits(string,n_nums,num_string,n_units,unit_string)
   do i = 1,n_units
     is = -1
     if (pos_units(i).gt.0) then
       do j = 1, n_units
         if (pos_units(j).lt.pos_units(i).and.pos_units(j).ne.-1) then
           if (pos_units(j)+1.gt.is) is = pos_units(j)+1
         endif
       enddo ! i i
       if (is.eq.-1) is = 1 ! s e
       ie = pos_units(i)-1 ! |....|
       nn = ie-is+1 ! xxxd123456h ->nn = 6
       write(nfmt,'(a,i2,a)') '(i',nn,') '
       read(string(is:ie),trim(nfmt)) times(i)
     endif
   enddo
!
! times to second
!
   do i = 1,n_units
     if (times(i).ne.0) then
       second = second+times(i)*factors(i)
     endif
   enddo
!
   end function string2second
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function getposunits(string, n_nums, num_string, n_units, unit_string)      &
                                                                result(unit_pos)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*),              intent(in   ) :: string
   integer(i4),                   intent(in   ) :: n_nums, n_units
   character, dimension(n_nums),  intent(in   ) :: num_string
   character, dimension(n_units), intent(in   ) :: unit_string
   integer(i4), dimension(n_units)              :: unit_pos
!
! local variables
!
   integer(i4) :: nchars, i, j
   integer(i4), dimension(n_units) :: unit_cnts
   character, dimension(n_nums+n_units) :: chk_string
   logical(l4) :: isok
!
   nchars = len(trim(string))
!
   chk_string(1:n_nums) = num_string(:)
   chk_string(n_nums+1:n_nums+n_units) = unit_string
!
   do j = 1,nchars
     isok = .false.
     do i = 1,n_nums+n_units
       if (string(j:j).eq.chk_string(i)) then
         isok = .true.
       endif
     enddo
     if (.not. isok) then
       print *, 'error : not allow unit('//string(j:j)//') '
       stop
     endif
   enddo
!
! get positions of units
!
   unit_cnts(:) = 0
   unit_pos(:)  = -1
   do j = 1,nchars
     do i = 1,n_units
       if (string(j:j).eq.unit_string(i)) then
         unit_cnts(i) = unit_cnts(i)+1
         unit_pos(i)  = j
       endif
     enddo
   enddo
!
! check last unit
!
   if (maxval(unit_pos(:)).lt.nchars) then
     unit_cnts(n_units) = unit_cnts(n_units)+1
     unit_pos(n_units)  = nchars+1
   endif
!
! check double units
!
   do i = 1,n_units
     if (unit_cnts(i).gt.1) then
       print *, 'error : double unit('//unit_string(i)//')'
       stop
     endif
   enddo
!
! check null numbers
!
   do i = 2,n_units -1
     do j = 1,n_units
       if (i.ne.j.and.unit_pos(i).ne.-1.and.unit_pos(j).ne.-1) then
         if (abs(unit_pos(i)-unit_pos(j)).lt.2) then
           print *, 'error : check null number('//unit_string(i)//', '//       &
                                                             unit_string(j)//')'
           stop
         endif
       endif
     enddo
   enddo

   end function getposunits
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function i4tostring(var) result(outstring)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: var
   character(len=16) :: outstring
!
! local variables
!
   integer(i4)      :: ndigits
   character(len=2) :: fmts
!
   ndigits = getdigits(var)
!
   write(fmts, '(i2)') ndigits+1
   write(outstring, '(i'//trim(fmts)//')') var
!
   end function i4tostring
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function time2string(second) result(ftimes)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: second
   character(len=32) :: ftimes
!
! local variables
!
   integer(i4), parameter :: n_units = 4
   character, dimension(n_units) :: unit_string
   integer, dimension(n_units) :: times, factors
   integer(i4) :: ndigits, i, sec
!
   sec = second
   unit_string = (/'d','h','m','s'/)
   times(:)    = 0 !(/day,hour,minute,second/)
   factors(:)  = (/24*60*60,60*60,60,1/)
!
! get times
!
   do i = 1,n_units
     times(i) = sec/factors(i)
     sec      = mod(sec,factors(i))
   enddo

! write string
   ftimes = ''
   do i = 1,n_units
     if (times(i).gt.0) then
       ndigits = getdigits(times(i))
       write(ftimes,'(a,i'//trim(i4tostring(ndigits))//',a)') trim(ftimes),times(i),unit_string(i)
     endif
   enddo
!
   end function time2string
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function getdigits_i4(var) result(length)
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4), intent(in   ) :: var
   integer(i4) :: length
!
! local variables
!
   integer(i4) :: maxi, i
!
   maxi = 10
!
   if (var==0) then
     length = 1
     return
   endif
!
   if (abs(var)>=10**(maxi-1).and.abs(var)<=huge(var)) then
     length = maxi
     return
   endif
!
   do i = 1,maxi-1
     if (abs(var)>=10**(i-1).and.abs(var)<10**i) then
       length = i
       return
     endif
   enddo
!
   end function getdigits_i4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function getdigits_r4(var) result(length)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r4), intent(in   ) :: var
   integer(i4)             :: length
!
! local variables
!
   integer(i4) :: maxi, i
!
   maxi = 39
!
   if (var==0.0_r4) then
     length = 1
     return
   endif
!
   if (abs(var)>=10.0_r4**(maxi-1).and.abs(var)<=huge(var)) then
     length = maxi
     return
   endif
!
   do i = 1,maxi-1
     if (abs(var)>=10.0_r4**(i-1).and.abs(var)<10.0_r4**i) then
       length = i
       return
     endif
   enddo
!
   end function getdigits_r4
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function getdigits_r8(var) result(length)
!-------------------------------------------------------------------------------
   implicit none
!
   real(r8), intent(in   ) :: var
   integer(i4)             :: length
!
! local variables
!
   integer(i4) :: maxi, i
!
   maxi = 309
!
   if (var==0.0_r8) then
     length = 1
     return
   endif
!
   if (abs(var)>=10.0_r8**(maxi-1).and.abs(var)<=huge(var)) then
     length = maxi
     return
   endif
!
   do i = 1,maxi-1
     if (abs(var)>=10.0_r8**(i-1).and.abs(var)<10.0_r8**i) then
       length = i
       return
     endif
   enddo
!
   end function getdigits_r8
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module time
!-------------------------------------------------------------------------------
