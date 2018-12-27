!------------------------------------------------------------------------------
   module opt_parser
!-------------------------------------------------------------------------------
!
!  abstract :  fortran option parser
!
!  history log :
!    2016-12-21  junghan kim    initial setup
!    2017-02-15  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use list, only: list_t, list_l4, list_i4, list_i8, list_r4, list_r8, list_st, &
                   list_initialize, list_finalize, list_add, list_change
   private
! kinds
   integer, parameter :: i4 = 4, r4 = 4, r8 = 8, l4 = 4
! internal parameters
   integer(i4), parameter :: max_nopts_    = 20
   integer(i4), parameter :: max_len_opt_  = 2
   integer(i4), parameter :: max_len_lopt_ = 16
   integer(i4), parameter :: max_len_str_  = 512 ! external
   integer(i4), parameter :: max_len_help_ = 32
   integer(i4), parameter :: max_len_arrs_ = 32
! internal parameters (actions)
   integer(i4), parameter :: opt_bool_     = 0
   integer(i4), parameter :: opt_store_    = 1
   integer(i4), parameter :: opt_append_   = 2
! external parameters (types)
   integer(i4), parameter :: opt_logical_  = 1
   integer(i4), parameter :: opt_int_      = 2
   integer(i4), parameter :: opt_real_     = 3
   integer(i4), parameter :: opt_real8_    = 4
   integer(i4), parameter :: opt_string_   = 5
!
   type opt_args_t
     ! help
     logical(l4) :: onhelp
     ! arguments
     integer(i4) :: nargs
     character(len=max_len_str_ ), dimension(:), allocatable :: args
     ! options
     integer(i4)  :: nopts
     character(len=max_len_opt_ ), dimension(max_nopts_) :: opts
     character(len=max_len_lopt_), dimension(max_nopts_) :: lopts
     integer(i4),                  dimension(max_nopts_) :: acts
     integer(i4),                  dimension(max_nopts_) :: types
     character(len=max_len_help_), dimension(max_nopts_) :: help
     type(list_t),                 dimension(max_nopts_) :: defs
     type(list_t),                 dimension(max_nopts_) :: vals
   end type opt_args_t
!
   interface opt_add
     module procedure opt_add_l4
     module procedure opt_add_i4
     module procedure opt_add_r4
     module procedure opt_add_r8
     module procedure opt_add_st
   end interface opt_add
   interface opt_get
     module procedure opt_get_l4
     module procedure opt_get_i4
     module procedure opt_get_r4
     module procedure opt_get_r8
     module procedure opt_get_st
   end interface opt_get
   interface string2array
     module procedure string2array_l4
     module procedure string2array_i4
     module procedure string2array_r4
     module procedure string2array_r8
     module procedure string2array_st
   end interface string2array
!
   public :: opt_args_t, max_len_str_, max_len_arrs_
   public :: opt_bool_, opt_store_, opt_append_
   public :: opt_logical_, opt_int_, opt_real_, opt_real8_, opt_string_
   public :: opt_initialize, opt_finalize, opt_add, opt_print
   public :: opt_get, opt_get_arg, opt_show_help
   public :: string2array
!
   contains
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_initialize(op, nargs, args)
   implicit none
   type(opt_args_t), intent(inout) :: op
   integer(i4),                    optional, intent(in   ) :: nargs
   character(len=*), dimension(:), optional, intent(in   ) :: args
! local variables
   integer(i4)        :: i
   character(len=128) :: message
   integer(i4)        :: nargs_l
   character(len=max_len_str_), dimension(:), allocatable :: args_l
!
   op%onhelp = .false.
!
   if (present(nargs).and.present(args)) then
     nargs_l = nargs
     allocate(args_l(nargs_l))
     do i = 1, nargs_l
       args_l(i) = trim(args(i))
     enddo
   else
     nargs_l = iargc()
     allocate(args_l(nargs_l))
     do i = 1, nargs_l
       call getarg(i,args_l(i))
       if (trim(args_l(i)).eq.'-h'.or.trim(args_l(i)).eq.'--help') then
         op%onhelp = .true.
       endif
     enddo
   endif
!
   op%nargs = nargs_l
   if (len(args_l).lt.nargs_l) then
     write(message,'(a,x,i2,x,i2)') 'opt_initialize: check size of args and nargs :',nargs_l,len(args_l)
     call opt_die(trim(message))
   endif
   allocate(op%args(nargs_l))
   do i = 1,nargs_l
     op%args(i) = trim(args_l(i))
   enddo
!
   deallocate(args_l)
!
   return
   end subroutine opt_initialize
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_finalize(op)
   implicit none
   type(opt_args_t), intent(inout) :: op
! local variables
   integer(i4) :: i
!
   deallocate(op%args)
!
   return
   end subroutine opt_finalize
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_add_l4(op, opt, lopt, act, def, help)
   implicit none
!
   type(opt_args_t),            intent(inout) :: op
   character(len=max_len_opt_), intent(in   ) :: opt
   character(len=*),            intent(in   ) :: lopt
   integer(i4),                 intent(in   ) :: act
   logical(l4),                 intent(in   ) :: def
   character(len=*),            intent(in   ) :: help
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   if (op%nopts.ge.max_nopts_) then
     write(message,'(a,i3)') 'opt_add: exceeded the maximum number of options: ',max_nopts_
     call opt_die(trim(message))
   endif
   if (act.ne.opt_bool_) then
     write(message,'(a)') 'opt_add: act for logical must be opt_bool_ '
     call opt_die(trim(message))
   endif
!
   call opt_check(op,opt,lopt)
!
   op%nopts = op%nopts+1
   op%opts(op%nopts)  = opt
   op%lopts(op%nopts) = trim(lopt)
   op%help(op%nopts)  = trim(help)
   op%types(op%nopts) = opt_logical_
   op%acts(op%nopts)  = opt_bool_
!
   call list_initialize(op%defs(op%nopts),list_l4)
   call list_add(op%defs(op%nopts),def)
   call list_initialize(op%vals(op%nopts),list_l4)
   call list_add(op%vals(op%nopts),def)
!
   call opt_set(op,trim(opt),op%nopts)
!
   return
   end subroutine opt_add_l4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_add_i4(op, opt, lopt, act, def, help)
   implicit none
!
   type(opt_args_t),            intent(inout) :: op
   character(len=max_len_opt_), intent(in   ) :: opt
   character(len=*),            intent(in   ) :: lopt
   integer(i4),                 intent(in   ) :: act
   integer(i4),                 intent(in   ) :: def
   character(len=*),            intent(in   ) :: help
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   if (op%nopts.ge.max_nopts_) then
     write(message,'(a,i3)') 'opt_add: exceeded the maximum number of options: ',max_nopts_
     call opt_die(trim(message))
   endif
   if (act.ne.opt_store_.and.act.ne.opt_append_) then
     write(message,'(a)') 'opt_add: check the act argument for i4 '
     call opt_die(trim(message))
   endif
!
   call opt_check(op,opt,lopt)
!
   op%nopts = op%nopts+1
   op%opts(op%nopts)  = opt
   op%lopts(op%nopts) = trim(lopt)
   op%help(op%nopts)  = trim(help)
   op%types(op%nopts) = opt_int_
   op%acts(op%nopts)  = act
!
   call list_initialize(op%defs(op%nopts),list_i4)
   call list_add(op%defs(op%nopts),def)
   call list_initialize(op%vals(op%nopts),list_i4)
   call list_add(op%vals(op%nopts),def)
!
   call opt_set(op,trim(opt),op%nopts)
!
   return
   end subroutine opt_add_i4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_add_r4(op, opt, lopt, act, def, help)
   implicit none
!
   type(opt_args_t),            intent(inout) :: op
   character(len=max_len_opt_), intent(in   ) :: opt
   character(len=*),            intent(in   ) :: lopt
   integer(i4),                 intent(in   ) :: act
   real(r4),                    intent(in   ) :: def
   character(len=*),            intent(in   ) :: help
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   if (op%nopts.ge.max_nopts_) then
     write(message,'(a,i3)') 'opt_add: exceeded the maximum number of options: ',max_nopts_
     call opt_die(trim(message))
   endif
   if (act.ne.opt_store_.and.act.ne.opt_append_) then
     write(message,'(a)') 'opt_add: check the act argument for r4 '
     call opt_die(trim(message))
   endif
!
   call opt_check(op,opt,lopt)
!
   op%nopts = op%nopts+1
   op%opts(op%nopts)  = opt
   op%lopts(op%nopts) = trim(lopt)
   op%help(op%nopts)  = trim(help)
   op%types(op%nopts) = opt_real_
   op%acts(op%nopts)  = act
!
   call list_initialize(op%defs(op%nopts),list_r4)
   call list_add(op%defs(op%nopts),def)
   call list_initialize(op%vals(op%nopts),list_r4)
   call list_add(op%vals(op%nopts),def)
!
   call opt_set(op,trim(opt),op%nopts)
!
   return
   end subroutine opt_add_r4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_add_r8(op, opt, lopt, act, def, help)
   implicit none
!
   type(opt_args_t),            intent(inout) :: op
   character(len=max_len_opt_), intent(in   ) :: opt
   character(len=*),            intent(in   ) :: lopt
   integer(i4),                 intent(in   ) :: act
   real(r8),                    intent(in   ) :: def
   character(len=*),            intent(in   ) :: help
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   if (op%nopts.ge.max_nopts_) then
     write(message,'(a,i3)') 'opt_add: exceeded the maximum number of options: ',max_nopts_
     call opt_die(trim(message))
   endif
   if (act.ne.opt_store_.and.act.ne.opt_append_) then
     write(message,'(a)') 'opt_add: check the act argument for r8 '
     call opt_die(trim(message))
   endif
!
   call opt_check(op,opt,lopt)
!
   op%nopts = op%nopts+1
   op%opts(op%nopts)  = opt
   op%lopts(op%nopts) = trim(lopt)
   op%help(op%nopts)  = trim(help)
   op%types(op%nopts) = opt_real8_
   op%acts(op%nopts)  = act
!
   call list_initialize(op%defs(op%nopts),list_r8)
   call list_add(op%defs(op%nopts),def)
   call list_initialize(op%vals(op%nopts),list_r8)
   call list_add(op%vals(op%nopts),def)
!
   call opt_set(op,trim(opt),op%nopts)
!
   return
   end subroutine opt_add_r8
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_add_st(op, opt, lopt, act, def, help)
   implicit none
!
   type(opt_args_t),            intent(inout) :: op
   character(len=max_len_opt_), intent(in   ) :: opt
   character(len=*),            intent(in   ) :: lopt
   integer(i4),                 intent(in   ) :: act
   character(len=*),            intent(in   ) :: def
   character(len=*),            intent(in   ) :: help
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   if (op%nopts.ge.max_nopts_) then
     write(message,'(a,i3)') 'opt_add: exceeded the maximum number of options: ',max_nopts_
     call opt_die(trim(message))
   endif
   if (act.ne.opt_store_.and.act.ne.opt_append_) then
     write(message,'(a)') 'opt_add: check the act argument for i4 '
     call opt_die(trim(message))
   endif
!
   call opt_check(op,opt,lopt)
!
   op%nopts = op%nopts+1
   op%opts(op%nopts)  = opt
   op%lopts(op%nopts) = trim(lopt)
   op%help(op%nopts)  = trim(help)
   op%types(op%nopts) = opt_string_
   op%acts(op%nopts)  = act
!
   call list_initialize(op%defs(op%nopts),list_st)
   call list_add(op%defs(op%nopts),def)
   call list_initialize(op%vals(op%nopts),list_st)
   call list_add(op%vals(op%nopts),def)
!
   call opt_set(op,trim(opt),op%nopts)
!
   return
   end subroutine opt_add_st
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_set(op, opt_name, iopt)
   implicit none
   type(opt_args_t), intent(inout) :: op
   character(len=*), intent(in   ) :: opt_name
   integer(i4),      intent(in   ) :: iopt
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   do i = 1,op%nargs
     if (trim(opt_name).eq.trim(op%args(i))) then
!
       if     (op%acts(iopt).eq.opt_bool_ ) then
!
         if (op%vals(iopt)%l4(1)) then
           op%vals(iopt)%l4(1) = .false.
         else
           op%vals(iopt)%l4(1) = .true.
         endif
         call opt_remove_arg(op,i,1)
         exit
!
       elseif (op%acts(iopt).eq.opt_store_) then ! not logical
!
         if (op%nargs.ge.i+1) then
           if (op%args(i+1)(1:1).ne.'-') then
             if     (op%types(iopt).eq.opt_int_   ) then  ! need error processing
               read(op%args(i+1),*) op%vals(iopt)%i4(1)
             elseif (op%types(iopt).eq.opt_real_  ) then
               read(op%args(i+1),*) op%vals(iopt)%r4(1)
             elseif (op%types(iopt).eq.opt_real8_ ) then
               read(op%args(i+1),*) op%vals(iopt)%r8(1)
             elseif (op%types(iopt).eq.opt_string_) then
               op%vals(iopt)%st(1) = trim(op%args(i+1))
               !read(op%args(i+1),*) op%vals(iopt)%st(1)
             endif
           else
             write(message,'(a,a)') "opt_set: check VALUE for the option: ", trim(opt_name)
             call opt_die(trim(message))
           endif
         else ! if nargs.nt.i+1
           write(message,'(a,a)') "opt_set: need VALUE for the option: ", trim(opt_name)
           call opt_die(trim(message))
         endif
         call opt_remove_arg(op,i,2)
         exit
!
       else  ! ! if .not.opt_name
         write(message,'(a,a,x,i2)') 'opt_set: check action for the option: ', trim(opt_name), op%acts(iopt)
         call opt_die(trim(message))
       endif
!
     endif ! if opt_name
   enddo
!
   return
   end subroutine opt_set
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_print(op)
   implicit none
   type(opt_args_t), intent(inout) :: op
! local variables
   integer(i4)        :: i
!   character(len=128) :: message
!
   do i=1,op%nopts
     if     (op%types(i).eq.opt_logical_) then
       write(*,'(i2,a,a,a,l)')     i, 'th options   : ', trim(op%opts(i)), ' : ', op%vals(i)%l4(1)
     elseif (op%types(i).eq.opt_int_    ) then
       write(*,'(i2,a,a,a,i4)')    i, 'th options   : ', trim(op%opts(i)), ' : ', op%vals(i)%i4(1)
     elseif (op%types(i).eq.opt_real_   ) then
       write(*,'(i2,a,a,a,e10.5)') i, 'th options   : ', trim(op%opts(i)), ' : ', op%vals(i)%r4(1)
     elseif (op%types(i).eq.opt_real8_  ) then
       write(*,'(i2,a,a,a,e10.5)') i, 'th options   : ', trim(op%opts(i)), ' : ', op%vals(i)%r8(1)
     elseif (op%types(i).eq.opt_string_ ) then
       write(*,'(i2,a,a,a,a)')     i, 'th options   : ', trim(op%opts(i)), ' : ', trim(op%vals(i)%st(1))
     endif
   enddo
   do i=1,op%nargs
     if (op%args(i).ne.'-h') then
     write(*,'(i2,a,a)')     i, 'th arguments :', trim(op%args(i))
     endif
   enddo
!
   return
   end subroutine opt_print
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_show_help(op)
   implicit none
   type(opt_args_t), intent(inout) :: op
! local variables
   integer(i4)        :: i
   character(len=08)  :: str_opt
   character(len=20)  :: str_lopt
   character(len=10)  :: str_def
   character(len=128) :: message
!
   call get_command(message)
   write(*,'(a)') 'Usage: '//message(1:index(message,' ')-1) //' [options]'
   write(*,'(a)') ' '
   write(*,'(a)') 'Options:'
   str_opt = '-h'; str_lopt = '--help'; str_def = 'off'
   write(message,'(2x,a8,a2,a20,a10,x,a)') adjustl(str_opt),', ',adjustl(str_lopt),adjustl(str_def),'show this help message and exit'
   write(*,'(a)') trim(message)
!
   do i = 1,op%nopts
     if     (op%types(i).eq.opt_logical_) then
       str_opt = op%opts(i); str_lopt = trim(op%lopts(i))
       if (op%defs(i)%l4(1)) then
         str_def = 'on'
       else
         str_def = 'off'
       endif
     else
       str_opt = op%opts(i)//' VALUE'; str_lopt = trim(op%lopts(i))//' VALUE'
       if     (op%types(i).eq.opt_int_   ) then
         write(str_def,'(i10)') op%defs(i)%i4(1)
       elseif (op%types(i).eq.opt_real_  ) then
         write(str_def,'(e10.4)') op%defs(i)%r4(1)
       elseif (op%types(i).eq.opt_real8_ ) then
         write(str_def,'(e10.4)') op%defs(i)%r8(1)
       elseif (op%types(i).eq.opt_string_) then
         str_def = trim(op%defs(i)%st(1))
       endif
     endif
     write(message,'(2x,a8,a2,a20,a10,x,a)') adjustl(str_opt),', ',adjustl(str_lopt),adjustl(str_def),trim(op%help(i))
     write(*,'(a)') trim(message)
   enddo
   stop
!
   return
   end subroutine opt_show_help
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_remove_arg(op, ipos, num)
   implicit none
   type(opt_args_t), intent(inout) :: op
   integer(i4),      intent(in   ) :: ipos, num
! local variables
   integer(i4)        :: j
   character(len=max_len_str_), dimension(:), allocatable :: targs
!
   allocate(targs(op%nargs))
   targs(:) = op%args(:)
!
   op%nargs = op%nargs-num
   deallocate(op%args)
   allocate(op%args(op%nargs))
   do j = 1,op%nargs
     if     (j.lt.ipos) then
       op%args(j) = targs(j)
     elseif (j.ge.ipos) then
       op%args(j) = targs(j+num)
     endif
   enddo
!
   deallocate(targs)
!
   return
   end subroutine opt_remove_arg
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_check(op, opt, lopt)
   implicit none
   type(opt_args_t), intent(inout) :: op
   character(len=*), intent(in   ) :: opt
   character(len=*), intent(in   ) :: lopt
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   do i = 1,op%nopts
     if (trim(op%opts(i)).eq.trim(opt).or.trim(op%lopts(i)).eq.trim(lopt)) then
       write(message,'(a,a,x,a)') "opt_check: the option's name already exists : ", trim(opt), trim(lopt)
       call opt_die(trim(message))
     elseif ('-h'.eq.trim(opt).or.'--help'.eq.trim(lopt)) then
       write(message,'(a,a,x,a)') "opt_check: the option's name already exists : ", trim(opt), trim(lopt)
       call opt_die(trim(message))
     endif
   enddo
!
   return
   end subroutine opt_check
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_check_type(type)
   implicit none
   integer(i4), intent(in   ) :: type
! local variables
   character(len=128) :: message
!
   if (type.ne.opt_logical_.and.type.ne.opt_int_.and.type.ne.opt_real_.and.type.ne.opt_real8_.and.type.ne.opt_string_) then
     write(message,'(a,i2)') "opt_check_type: the option's type : ", type
     call opt_die(trim(message))
   endif
!
   return
   end subroutine opt_check_type
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_get_l4(op, opt, res)
   implicit none
   type(opt_args_t), intent(inout) :: op
   character(len=*), intent(in   ) :: opt
   logical(l4),      intent(  out) :: res
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   if (op%onhelp) then
     call opt_show_help(op)
   endif
!
   do i = 1,op%nopts
     if (trim(op%opts(i)).eq.trim(opt)) then
       if (op%types(i).eq.opt_logical_) then
         res = op%vals(i)%l4(1)
       else
         write(message,'(a,a)') 'opt_get: check the type of option: ', opt
         call opt_die(trim(message))
       endif
     endif
   enddo
!
   return
   end subroutine opt_get_l4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_get_i4(op, opt, res)
   implicit none
   type(opt_args_t), intent(inout) :: op
   character(len=*), intent(in   ) :: opt
   integer(i4),      intent(  out) :: res
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   if (op%onhelp) then
     call opt_show_help(op)
   endif
!
   do i = 1,op%nopts
     if (trim(op%opts(i)).eq.trim(opt)) then
       if (op%types(i).eq.opt_int_) then
         res = op%vals(i)%i4(1)
       else
         write(message,'(a,a)') 'opt_get: check the type of option: ', opt
         call opt_die(trim(message))
       endif
     endif
   enddo
!
   return
   end subroutine opt_get_i4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_get_r4(op, opt, res)
   implicit none
   type(opt_args_t), intent(inout) :: op
   character(len=*), intent(in   ) :: opt
   real(r4),         intent(  out) :: res
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   if (op%onhelp) then
     call opt_show_help(op)
   endif
!
   do i = 1,op%nopts
     if (trim(op%opts(i)).eq.trim(opt)) then
       if (op%types(i).eq.opt_real_) then
         res = op%vals(i)%r4(1)
       else
         write(message,'(a,a)') 'opt_get: check the type of option: ', opt
         call opt_die(trim(message))
       endif
     endif
   enddo
!
   return
   end subroutine opt_get_r4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_get_r8(op, opt, res)
   implicit none
   type(opt_args_t), intent(inout) :: op
   character(len=*), intent(in   ) :: opt
   real(r8),         intent(  out) :: res
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   if (op%onhelp) then
     call opt_show_help(op)
   endif
!
   do i = 1,op%nopts
     if (trim(op%opts(i)).eq.trim(opt)) then
       if (op%types(i).eq.opt_real8_) then
         res = op%vals(i)%r8(1)
       else
         write(message,'(a,a)') 'opt_get: check the type of option: ', opt
         call opt_die(trim(message))
       endif
     endif
   enddo
!
   return
   end subroutine opt_get_r8
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_get_st(op, opt, res)
   implicit none
   type(opt_args_t),            intent(inout) :: op
   character(len=*),            intent(in   ) :: opt
   character(len=max_len_str_), intent(  out) :: res
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   if (op%onhelp) then
     call opt_show_help(op)
   endif
!
   do i = 1,op%nopts
     if (trim(op%opts(i)).eq.trim(opt)) then
       if (op%types(i).eq.opt_string_) then
         res = trim(op%vals(i)%st(1))
       else
         write(message,'(a,a)') 'opt_get: check the type of option: ', opt
         call opt_die(trim(message))
       endif
     endif
   enddo
!
   return
   end subroutine opt_get_st
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_get_arg(op, pos, res)
   implicit none
   type(opt_args_t),            intent(inout) :: op
   integer(i4),                 intent(in   ) :: pos
   character(len=max_len_str_), intent(  out) :: res
! local variables
   integer(i4)        :: i
   character(len=128) :: message
!
   if (op%onhelp) then
     call opt_show_help(op)
   endif
!
   if (pos.gt.op%nargs) then
     write(message,'(a,i2,a,i2)') 'check argument position: ', pos, ', nargs = ', op%nargs
     call opt_die(trim(message))
   endif
!
   res = trim(op%args(pos))
!
   return
   end subroutine opt_get_arg
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine opt_die(message)
   implicit none
   character(len=*), intent(in   ) :: message
! local variables
!
   write(*,*) trim(message)
   stop
!
   return
   end subroutine opt_die
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   function string2format(string) result(res)
   implicit none
   character(len=*), intent(in   ) :: string
   character(len=8)                :: res
! local variables
   logical(l4) :: isok
   integer(i4) :: i, j, len_str, pos, npoints, len_chk
   character, dimension(11) :: chk_num
!
   len_chk = 11
   chk_num = (/'0','1','2','3','4','5','6','7','8','9','.'/)
!
   len_str = len(trim(string))
!
! check number
!
   do i = 1, len_str
     isok = .false.
     do j = 1, len_chk
       if (string(i:i).eq.chk_num(j)) then
         isok = .true.
       end if
     end do
     if (.not.isok) then
       print *, 'argument is not number...', string
       stop
     end if
   end do
!
! determine number
!
   pos = -1
   do i = 1, len_str
     if (string(i:i).eq.'.') then
       if (pos.eq.-1) then
         pos = i
       else
         print *, 'check argument (point)...'
         stop
       end if
     end if
   end do
!
   write(res,'(a2,i2.2,a,i2.2,a1)') '(f', (len_str+1), '.', (len_str-pos), ')'
!
   end function string2format
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine string2array_l4(string, array, divisor)
   implicit none
   character(len=*),                       intent(in   ) :: string
   logical(l4), dimension(:), allocatable, intent(inout) :: array
   character(len=1), optional,             intent(in   ) :: divisor
! local variables
   character(len=1) :: div
   integer(i4)      :: len_str, ndivs, narray
   integer(i4)      :: n, i, j, is, ie
   integer(i4), dimension(100) :: idivs
!
   if (present(divisor)) then
     div = divisor
   else
     div = ','
   endif
!
   len_str = len(string)
   ndivs = 0
   do i = 1,len_str
     if (string(i:i).eq.div) then
       ndivs = ndivs+1
       idivs(ndivs) = i
     endif
   enddo
!
   allocate(array(ndivs+1))
   if (ndivs.gt.0) then
     do i = 1, ndivs+1
       !
       if (i.eq.1) then
         is =            1; ie = idivs(i)-1
       elseif (i.eq.ndivs+1) then
         is = idivs(i-1)+1; ie = len_str
       else
         is = idivs(i-1)+1; ie = idivs(i)-1
       endif
       !
       if (string(is:ie).eq.'.true.') then
         array(i) = .true.
       else
         array(i) = .false.
       endif
       !
     enddo
   endif
!
   return
   end subroutine string2array_l4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine string2array_i4(string, array, divisor)
   implicit none
   integer(i4), dimension(:), allocatable, intent(inout) :: array
#include <string2array.inc>
!
   return
   end subroutine string2array_i4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine string2array_r4(string, array, divisor)
   implicit none
   real(r4), dimension(:), allocatable, intent(inout) :: array
#include <string2array.inc>
!
   return
   end subroutine string2array_r4
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine string2array_r8(string, array, divisor)
   implicit none
   real(r8), dimension(:), allocatable, intent(inout) :: array
#include <string2array.inc>
!
   return
   end subroutine string2array_r8
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   subroutine string2array_st(string, array, divisor)
   implicit none
   character(len=max_len_arrs_), dimension(:), allocatable, intent(inout) :: array
#include <string2array.inc>
!
   return
   end subroutine string2array_st
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
   end module opt_parser
!------------------------------------------------------------------------------
